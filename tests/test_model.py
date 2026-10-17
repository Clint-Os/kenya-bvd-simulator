"""Sanity checks that pin the headline numbers quoted in the article."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from model import R0, build_params  # noqa: E402
from stoch import simulate  # noqa: E402


def test_r0_nairobi_with_isolation_beds():
    p = build_params("Nairobi", R0_ref=2.0, beds=57)
    assert abs(R0(p, True) - 1.95) < 0.01


def test_r0_mombasa_with_isolation_beds():
    p = build_params("Mombasa", R0_ref=2.0, beds=4)
    assert abs(R0(p, True) - 2.18) < 0.01


def test_r0_roughly_doubles_when_wards_are_full():
    for city, beds in (("Nairobi", 57), ("Mombasa", 4)):
        p = build_params(city, R0_ref=2.0, beds=beds)
        assert R0(p, False) > 1.9 * R0(p, True)


def test_more_beds_improve_containment():
    few = (simulate("Mombasa", 4, nsim=400, seed=1)["total"] < 20_000).mean()
    many = (simulate("Mombasa", 40, nsim=400, seed=1)["total"] < 20_000).mean()
    assert many > few
