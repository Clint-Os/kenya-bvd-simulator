"""Mombasa isolation-bed sweep: probability of containment vs total beds.

Reproduces the curves behind panel 2 of the simulator and the Utange sizing in
the article. Containment = cumulative cases below 20,000 within 400 days, with
two index cases and a response (60% contact cut, 1.5x isolation speed, 90% safe
burial) starting `tau` days after the first isolation admission.

Reported both as total Mombasa beds and as Utange beds (= total - 4 at Port Reitz).

Usage:  python containment_sweep.py [--nsim 1200] [--out ../data/containment_sweep.json]
"""
import argparse
import json
from pathlib import Path

import numpy as np

from stoch import simulate

PORT_REITZ = 4
CONTAINED = 20_000
BEDS = list(range(4, 80, 3))
TAUS = (3, 7, 14, 21)
TARGETS = (0.90, 0.95, 0.99)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nsim", type=int, default=1200)
    ap.add_argument("--out", default=str(Path(__file__).resolve().parent.parent / "data" / "containment_sweep.json"))
    a = ap.parse_args()

    curves = {}
    for tau in TAUS:
        curves[str(tau)] = [
            float((simulate("Mombasa", b, nsim=a.nsim, tau=tau, seed=b)["total"] < CONTAINED).mean())
            for b in BEDS
        ]
    # smallest total bed count whose containment probability reaches each target:
    # locate the crossing on the coarse grid, then refine bed-by-bed with 3x the runs
    def p_contain(b, tau):
        return float((simulate("Mombasa", b, nsim=3 * a.nsim, tau=tau, seed=1000 + b)["total"] < CONTAINED).mean())

    thresholds = {}
    for tau, c in curves.items():
        thresholds[tau] = {}
        for tg in TARGETS:
            j = next((i for i, p in enumerate(c) if p >= tg), None)
            if j is None:
                thresholds[tau][str(tg)] = None
                continue
            lo = BEDS[max(j - 1, 0)] + 1 if j else BEDS[0]
            total = next((b for b in range(lo, BEDS[j] + 1) if p_contain(b, int(tau)) >= tg), BEDS[j])
            thresholds[tau][str(tg)] = {"total_beds": total, "utange_beds": total - PORT_REITZ}
    out = {"beds_total": BEDS, "curves": curves, "beds_for_target": thresholds,
           "nsim": a.nsim, "contained_below": CONTAINED}
    Path(a.out).write_text(json.dumps(out, indent=1))
    print(json.dumps(thresholds, indent=1))


if __name__ == "__main__":
    main()
