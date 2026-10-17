"""Facility stock planner: turns expected beds-in-use into PPE and IPC supplies.

Demand is built bottom-up from a care model, not from a per-patient norm:

    occupied beds -> red-zone entries -> staff entries -> PPE sets (+ gloves, waste,
    chlorine) ; deaths -> body bags ; bed-days -> water and sharps ; beds -> fixtures
    (basins, buckets, bins), which are counted at peak occupancy rather than summed.

All defaults live in ../data/supplies_defaults.json, which index.html embeds, so the
simulator and this module compute identical numbers (tests/test_supplies.py checks
the copies match). Defaults are planning assumptions; replace them with the facility
SOP or Ministry of Health IPC guidance.

Usage:
    python supplies.py --city Mombasa --beds 26 --days 365
    python supplies.py --full-unit --beds 26            # worst case: every bed occupied
"""
import argparse
import json
import math
from pathlib import Path

import numpy as np

DEFAULTS_PATH = Path(__file__).resolve().parent.parent / "data" / "supplies_defaults.json"
CFR_HOSP = 0.31          # model.NATURAL_HISTORY
HOSP_TO_OUTCOME = 6.0    # days, model.NATURAL_HISTORY


def load_defaults(path=DEFAULTS_PATH):
    return json.loads(Path(path).read_text())


def daily_demand(occ, cfg=None):
    """Daily consumption for each item given a daily series of occupied beds.

    occ: array of expected occupied beds in the facility, one value per day.
    Returns {item_id: array}, plus 'staff_entries', 'ppe_sets' and 'deaths'.
    """
    cfg = cfg or load_defaults()
    d = cfg["drivers"]
    occ = np.asarray(occ, dtype=float)
    on = occ >= 0.5                                   # fractional patients below half a bed do not open a ward
    clin_entries = np.where(
        on, np.maximum(d["visits_per_patient_day"] * occ / d["patients_per_entry"],
                       d["min_entries_per_day_if_occupied"]), 0.0)
    clin_staff = clin_entries * d["staff_per_entry"]
    support = np.where(on, d["support_entries_per_bed_day"] * occ, 0.0)
    sets = clin_staff + support                        # one full PPE set per staff entry
    seen = np.clip(occ, 1.0, d["patients_per_entry"])  # patients seen on one entry
    deaths = occ * CFR_HOSP / HOSP_TO_OUTCOME

    ch = cfg["chlorine"]
    litres_05 = ch["litres_0_5pct_standing_per_bed"] * ch["refreshes_per_day"]   # MoH: fresh solution every 8 h
    chlorine_kg = (occ * litres_05 * 5.0 +
                   sets * ch["litres_0_05pct_per_entry"] * 0.5) / ch["available_fraction"] / 1000.0

    out = {"staff_entries": clin_staff, "ppe_sets": sets, "deaths": deaths, "occupied": occ}
    for it in cfg["items"]:
        b, q, reuse = it["basis"], it["qty"], it["reuse"]
        if b == "entry":
            # clinical entries use qty (+ per_change for each extra patient seen); support entries use support_qty
            sq = it.get("support_qty", q)
            v = (clin_staff * (q + it.get("per_change", 0.0) * (seen - 1.0)) + support * sq) / reuse
        elif b == "support":
            v = support * q / reuse
        elif b == "fixture":
            v = np.where(on, occ * q, 0.0)               # equipment per occupied bed, counted at peak (see plan)
        elif b == "death":
            v = deaths * q
        elif b == "chlorine":
            v = chlorine_kg
        elif b == "bedday":
            v = occ * q
        else:
            raise ValueError(f"unknown basis {b!r}")
        out[it["id"]] = v
    return out


def plan(occ, cfg=None):
    """Stock requirement per item: the worst rolling window, plus a safety buffer."""
    cfg = cfg or load_defaults()
    d = cfg["drivers"]
    dem = daily_demand(occ, cfg)
    w = int(d["window_days"])
    rows = []
    for it in cfg["items"]:
        x = dem[it["id"]]
        if it["basis"] == "fixture":                      # equipment is counted at peak occupancy, not summed over days
            peak_window = float(x.max())
        else:
            win = np.convolve(x, np.ones(min(w, len(x))), mode="valid")   # rolling window sums
            peak_window = float(win.max()) if len(win) else 0.0
        rows.append({
            "id": it["id"], "name": it["name"], "unit": it["unit"],
            "peak_per_day": float(x.max()),
            "window_total": peak_window,
            "stock_needed": math.ceil(peak_window * (1 + d["buffer"]) - 1e-9),
        })
    return rows, dem


def occupancy_from_model(city, beds, days=365, R0_ref=2.0, share=1.0, t_resp=21,
                         contact_cut=0.6, isolation_speedup=1.5, s_c_resp=0.9, **overrides):
    """Expected (mean-field) occupied beds from the deterministic model in model.py.

    The response starts at fixed day `t_resp` here; index.html instead starts it a set
    delay after the first isolation admission, so the two agree only approximately.
    """
    from model import build_params, simulate
    p = build_params(city, R0_ref=R0_ref, beds=beds)
    p.update(t_resp=t_resp, contact_cut=contact_cut, isolation_speedup=isolation_speedup,
             s_c_resp=s_c_resp, **overrides)
    h = simulate(p, days=days, seed_I=(0, 2), seed_E=(0, 0))["H_total"]
    return share * np.minimum(h, beds)      # beds are a hard cap; the soft cap in the ODE leaks slightly


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--city", default="Mombasa", choices=["Nairobi", "Mombasa"])
    ap.add_argument("--beds", type=int, default=26)
    ap.add_argument("--days", type=int, default=365)
    ap.add_argument("--share", type=float, default=1.0, help="share of city isolation patients at this facility")
    ap.add_argument("--full-unit", action="store_true", help="assume every bed occupied")
    a = ap.parse_args()
    cfg = load_defaults()
    cfg["drivers"]["facility_share"] = a.share
    occ = (np.full(a.days + 1, a.beds * a.share, float) if a.full_unit
           else occupancy_from_model(a.city, a.beds, a.days, share=a.share))
    rows, dem = plan(occ, cfg)
    print(f"peak occupied beds {occ.max():.1f}; peak PPE sets/day {dem['ppe_sets'].max():.0f}")
    for r in rows:
        print(f"{r['name'][:52]:<54}{r['peak_per_day']:>10.1f}/day{r['stock_needed']:>10,d} {r['unit']}")


if __name__ == "__main__":
    main()
