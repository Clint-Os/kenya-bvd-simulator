"""Stochastic tau-leap of the two-group SEIHFR for isolation-capacity planning.

Vectorised over `nsim` runs. Beds are a hard cap: admission happens only when a
bed is free (Kucharski et al. 2015). Response starts `tau` days after the first
admission (detection). Used to size Mombasa's isolation beds (Port Reitz fixed
at 4; Utange the free variable) so the outbreak is contained with target
probability. Containment = cumulative cases below `cap` at the horizon.
"""
import numpy as np
from model import build_params


def simulate(city, beds, nsim=1500, days=400, seed_cases=2, tau=7, cut=0.6,
             speedup=1.5, s_c_resp=0.9, dt=0.25, seed=0, draws=None):
    rng = np.random.default_rng(seed)
    p = build_params(city, R0_ref=2.0, beds=beds)
    Ng = np.round(p["Ng"]).astype(np.int64)
    n = nsim
    S = np.tile(Ng, (n, 1)); E = np.zeros((n, 2), np.int64)
    I = np.zeros((n, 2), np.int64); H = np.zeros((n, 2), np.int64)
    Fc = np.zeros((n, 2), np.int64); C = np.zeros((n, 2), np.int64)
    I[:, 1] = seed_cases; S[:, 1] -= seed_cases; C[:, 1] = seed_cases
    tdet = np.full(n, np.inf); peakH = np.zeros(n)
    c = p["c"]; Ctot = np.sum(c * Ng)
    g = {k: (np.full(n, p[k]) if draws is None or k not in draws else np.asarray(draws[k]))
         for k in ("betaI", "betaH", "betaF", "sigma", "gI", "gH", "gF", "cfr_c", "cfr_h")}
    rho = p["rho"]
    pois = rng.poisson
    for step in range(int(days / dt)):
        t = step * dt
        live = (E.sum(1) + I.sum(1) + H.sum(1) + Fc.sum(1)) > 0
        if not live.any():
            break
        on = t >= (tdet + tau)
        cutv = np.where(on, cut, 0.0)[:, None]
        sc = np.where(on, s_c_resp, p["s_c"])
        sp = np.where(on, speedup, 1.0)
        Htot = H.sum(1)
        lamI = (1 - cutv[:, 0]) * g["betaI"] * (c * I).sum(1) / Ctot
        lamF = g["betaF"] * (c * Fc).sum(1) / Ctot
        lamH = g["betaH"] * Htot / p["N"]
        foi = c[None, :] * (lamI + lamF)[:, None] + lamH[:, None]
        newE = pois(np.clip(foi * S, 0, None) * dt).astype(np.int64)
        newE = np.minimum(newE, S)
        prog = lambda x, rate: np.minimum(pois(x * rate * dt).astype(np.int64), x)
        EtoI = prog(E, g["sigma"][:, None])
        free = np.clip(beds - Htot, 0, None)
        wantadm = prog(I, (rho * sp[:, None]))
        # admit up to free beds, informal group first (higher need)
        adm = np.zeros_like(I)
        for grp in (1, 0):
            take = np.minimum(wantadm[:, grp], free)
            adm[:, grp] = take; free = free - take
        ItoO = prog(I - adm, g["gI"][:, None])
        HtoO = prog(H, g["gH"][:, None])
        FtoR = prog(Fc, g["gF"][:, None])
        dc = rng.binomial(ItoO, g["cfr_c"][:, None]); dh = rng.binomial(HtoO, g["cfr_h"][:, None])
        S -= newE; E += newE - EtoI
        I += EtoI - adm - ItoO; H += adm - HtoO
        # only unsafe burials enter F (safe & dignified burials skip it)
        unsafe_c = rng.binomial(dc, np.clip(1 - sc, 0, 1)[:, None])
        unsafe_h = rng.binomial(dh, np.clip(1 - p["s_h"], 0, 1))
        Fc += unsafe_c + unsafe_h - FtoR
        Fc = np.clip(Fc, 0, None)
        C += EtoI
        newdet = np.isinf(tdet) & (adm.sum(1) > 0)
        tdet[newdet] = t
        peakH = np.maximum(peakH, H.sum(1))
    total = C.sum(1)
    return dict(total=total, peakH=peakH, tdet=tdet)
