import numpy as np, json
from model import *
from bifurcation import params as bparams, equilibria
from pathlib import Path
ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / 'data'; FIGS = ROOT / 'figures'; FIGS.mkdir(exist_ok=True)

RESP = dict(t_resp=21, isolation_speedup=1.5, s_c_resp=0.9)
CONTAINED = 20_000

def final_cases(city, R, B, cut, omega=0.0):
    p = build_params(city, R0_ref=R, beds=B); p.update(RESP, contact_cut=cut, omega=omega)
    return simulate(p, days=540)['C'][-1]

def critical_beds(city, R, cut, omega=0.0, lo=5, hi=20000):
    if final_cases(city, R, hi, cut, omega) > CONTAINED: return np.inf
    if final_cases(city, R, lo, cut, omega) < CONTAINED: return lo
    for _ in range(14):
        m = np.sqrt(lo*hi)
        if final_cases(city, R, m, cut, omega) < CONTAINED: hi = m
        else: lo = m
    return hi

out = {}
for city in CITIES:
    N = CITIES[city]['population']
    for R in (1.5, 2.0, 2.5):
        p = build_params(city, R0_ref=R)
        row = dict(R0=R0(p), R_no_beds=R0(p, False), parts=R0_by_setting(p), crit={})
        for cut in (0.45, 0.5, 0.55, 0.6, 0.65, 0.7):
            b = critical_beds(city, R, cut)
            row['crit'][cut] = dict(beds=b, per100k=b/N*1e5)
        row['crit_omega30'] = {cut: critical_beds(city, R, cut, omega=0.3) for cut in (0.5, 0.55, 0.6)}
        out[f'{city}|{R}'] = row
        print(city, R, 'R0=%.2f Rnob=%.2f'%(row['R0'], row['R_no_beds']),
              {k: (round(v['beds']) if np.isfinite(v['beds']) else 'inf') for k, v in row['crit'].items()},
              'omega30', {k: (round(v) if np.isfinite(v) else 'inf') for k, v in row['crit_omega30'].items()}, flush=True)
json.dump(out, open(DATA / 'critical_beds.json', 'w'), default=lambda o: None if not np.isfinite(o) else float(o))
