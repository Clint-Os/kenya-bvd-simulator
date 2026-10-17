"""Hand-checked cases for the facility stock planner, plus a sync check with index.html."""
import json
import re
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

from supplies import daily_demand, load_defaults, plan  # noqa: E402


def test_ten_patients_hand_calculation():
    d = daily_demand(np.full(5, 10.0))
    # 4 visits x 10 patients / 4 per entry = 10 entries x 2 staff = 20; plus 0.5 x 10 = 5 support
    assert d["ppe_sets"][0] == pytest.approx(25)
    assert d["coverall"][0] == pytest.approx(25)
    assert d["faceshield"][0] == pytest.approx(2.5)            # reused 10 times
    # WHO 2025: two pairs of gloves, outer pair renewed between patients (3 extra patients per entry);
    # support staff (5 entries) wear the inner pair plus heavy-duty outer gloves
    assert d["gloves"][0] == pytest.approx(20 * (2 + 3) + 5 * 1)
    assert d["apron"][0] == pytest.approx(20 * (1 + 3))            # disposable apron changed between patients; none for support staff
    assert d["mask"][0] == pytest.approx(25)                       # medical mask, not N95
    assert d["respirator"][0] == 0                                 # N95 only for aerosol-generating procedures
    assert d["heavyapron"][0] == pytest.approx(5 / 30)             # reusable, support staff only
    assert d["handrub"][0] == pytest.approx(20 * (0.15 + 0.05 * 3) + 5 * 0.15)
    assert d["bodybag"][0] == pytest.approx(10 * 0.31 / 6)
    assert d["water"][0] == pytest.approx(1.5)
    assert d["chlorine"][0] == pytest.approx((10 * 15 * 5 + 25 * 2 * 0.5) / 0.675 / 1000)


def test_moh_listed_items_and_fixtures():
    d = daily_demand(np.full(2, 12.0))
    assert d["scrubs"][0] == pytest.approx(d["ppe_sets"][0] / 3)           # laundered, 3 entries per set
    assert d["heavygloves"][0] == pytest.approx(0.5 * 12 / 10)             # support entries only
    assert d["basins"][0] == pytest.approx(2.0, abs=1e-3)                  # MoH: one sink per six-bed cubicle
    assert d["sharps"][0] == pytest.approx(12 * 0.05)


def test_fixtures_are_peak_counts_not_window_sums():
    occ = np.concatenate([np.full(10, 6.0), np.full(10, 12.0)])
    rows, _ = plan(occ)
    basins = next(r for r in rows if r["id"] == "basins")
    assert basins["window_total"] == pytest.approx(2.0, abs=1e-3)
    assert basins["stock_needed"] == 3                                    # ceil(2.0 x 1.25)


def test_chlorine_batches_follow_moh_refresh_rule():
    cfg = load_defaults()
    base = daily_demand(np.full(2, 10.0), cfg)["chlorine"][0]
    cfg["chlorine"]["refreshes_per_day"] = 1
    assert daily_demand(np.full(2, 10.0), cfg)["chlorine"][0] < 0.5 * base


def test_empty_and_sub_half_bed_wards_use_nothing():
    d = daily_demand(np.array([0.0, 0.3]))
    assert d["ppe_sets"].sum() == 0 and d["gloves"].sum() == 0 and d["basins"].sum() == 0


def test_single_patient_still_needs_a_team_each_shift():
    d = daily_demand(np.array([1.0]))
    assert d["ppe_sets"][0] == pytest.approx(3 * 2 + 0.5)


def test_stock_is_worst_window_plus_buffer():
    occ = np.concatenate([np.full(20, 2.0), np.full(60, 10.0)])
    rows, _ = plan(occ)
    cov = next(r for r in rows if r["id"] == "coverall")
    assert cov["window_total"] == pytest.approx(25 * 30)
    assert cov["stock_needed"] == 938                           # ceil(750 x 1.25)


def test_doubling_staff_per_entry_roughly_doubles_clinical_sets():
    cfg = load_defaults()
    base = daily_demand(np.full(3, 20.0), cfg)["ppe_sets"][0]
    cfg["drivers"]["staff_per_entry"] = 4
    assert daily_demand(np.full(3, 20.0), cfg)["ppe_sets"][0] > 1.7 * base


def test_html_embeds_the_same_defaults():
    html = (ROOT / "index.html").read_text()
    m = re.search(r'<script type="application/json" id="supplies-defaults">(.*?)</script>', html, re.S)
    assert m, "defaults block missing from index.html"
    assert json.loads(m.group(1)) == load_defaults()
