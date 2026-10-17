"""
Equilibrium and Hopf bifurcation analysis with bed capacity B as parameter.

A plain SEIR outbreak model has no endemic equilibrium, so it cannot undergo a
Hopf bifurcation. Following the hospital-bed literature (Shan & Zhu 2014;
Wan & Cui 2013; SIHR models with empty-bed-dependent admission), we add
population turnover (births and non-Ebola deaths, rate mu) and an admission
rate that saturates with occupancy: rho * B / (B + H). We then
  1. trace the endemic equilibrium as B varies (numerical continuation),
  2. compute the Jacobian eigenvalues along the branch,
  3. flag Hopf points where a complex pair crosses the imaginary axis,
     and saddle-node/backward-bifurcation folds where a real eigenvalue does.

Single homogeneous population; per-capita scaling (fractions), B in beds per
100,000 people so the two cities can be compared directly.
"""
import numpy as np
from scipy.optimize import fsolve
from model import NATURAL_HISTORY, R0_SHARES


def params(R0_ref=1.5, rho_days=5.0, mu_years=None, c=1.0, s_c=0.5, s_h=0.9,
           omega_wane_years=None):
    nh = NATURAL_HISTORY
    sigma, gI, gH, gF = 1/nh["incubation"], 1/nh["onset_to_outcome"], 1/nh["hosp_to_outcome"], 1/nh["death_to_burial"]
    rho_ref = 1/5.0
    t_I = 1/(rho_ref+gI); p_h = rho_ref/(rho_ref+gI)
    p_cd = gI*nh["cfr_comm"]/(rho_ref+gI)
    bI = R0_ref*R0_SHARES["community"]/t_I
    bH = R0_ref*R0_SHARES["hospital"]/(p_h/gH)
    bF = R0_ref*R0_SHARES["funeral"]/((p_cd+p_h*nh["cfr_hosp"])/gF)
    mu = 1/(365*(mu_years or 64.0))  # Kenya life expectancy ~64 y
    w = 0.0 if omega_wane_years is None else 1/(365*omega_wane_years)
    return dict(sigma=sigma, gI=gI, gH=gH, gF=gF, cc=nh["cfr_comm"], ch=nh["cfr_hosp"],
                bI=bI*c, bH=bH, bF=bF*c, rho=1/rho_days, mu=mu, s_c=s_c, s_h=s_h, w=w)


def f(x, p, B):
    """x = (S,E,I,H,F) as fractions; R = 1 - sum - (dead removed, births replace)."""
    S, E, I, H, F = x
    b = B / 1e5                     # beds as a fraction of population
    lam = p["bI"]*I + p["bH"]*H + p["bF"]*F
    adm = p["rho"] * b/(b + H) * I
    mu = p["mu"]
    R = 1 - S - E - I - H - F
    dS = mu - lam*S - mu*S + p["w"]*R
    dE = lam*S - (p["sigma"]+mu)*E
    dI = p["sigma"]*E - adm - (p["gI"]+mu)*I
    dH = adm - (p["gH"]+mu)*H
    dF = p["gI"]*p["cc"]*(1-p["s_c"])*I + p["gH"]*p["ch"]*(1-p["s_h"])*H - p["gF"]*F
    # Ebola deaths are replaced through the birth term mu (constant population
    # assumption); acceptable at equilibrium since prevalence is small.
    return np.array([dS, dE, dI, dH, dF])


def jac(x, p, B, h=1e-9):
    J = np.zeros((5, 5)); f0 = f(x, p, B)
    for k in range(5):
        dx = np.zeros(5); dx[k] = h*max(1, abs(x[k]))
        J[:, k] = (f(x+dx, p, B) - f0)/dx[k]
    return J


def continuation(p, B_grid, x0=None):
    """Trace endemic equilibrium over B_grid (descending from large B)."""
    out = []
    x = x0 if x0 is not None else np.array([0.6, 1e-4, 1e-4, 1e-5, 1e-6])
    for B in B_grid:
        sol, info, ier, msg = fsolve(f, x, args=(p, B), full_output=True, xtol=1e-14)
        if ier == 1 and sol[2] > 1e-12 and np.all(sol > -1e-12):
            x = sol
            ev = np.linalg.eigvals(jac(sol, p, B))
            out.append((B, sol.copy(), ev))
        else:
            out.append((B, None, None))
    return out


def find_hopf(branch):
    """Indices where the leading complex pair's real part changes sign."""
    pts = []
    prev = None
    for B, x, ev in branch:
        if ev is None:
            prev = None; continue
        cpx = ev[np.abs(ev.imag) > 1e-10]
        lead = cpx.real.max() if len(cpx) else np.nan
        if prev is not None and np.isfinite(lead) and np.isfinite(prev[1]) and np.sign(lead) != np.sign(prev[1]):
            pts.append((prev[0], B, prev[1], lead))
        prev = (B, lead)
    return pts


def equilibria(p, B, n=4000):
    """All endemic equilibria at bed level B (per 100k): 1-D root find in I."""
    b = B/1e5
    mu, s, g = p["mu"], p["sigma"], p
    def build(I):
        k = g["gH"]+mu
        # rho*b*I/(b+H) = k H  ->  k H^2 + k b H - rho b I = 0
        H = (-k*b + np.sqrt((k*b)**2 + 4*k*g["rho"]*b*I))/(2*k)
        adm = g["rho"]*b/(b+H)*I
        F = (g["gI"]*g["cc"]*(1-g["s_c"])*I + g["gH"]*g["ch"]*(1-g["s_h"])*H)/g["gF"]
        E = (adm + (g["gI"]+mu)*I)/s
        lam = g["bI"]*I + g["bH"]*H + g["bF"]*F
        S = (s+mu)*E/lam
        R = 1-S-E-I-H-F
        res = mu - lam*S - mu*S + g["w"]*R
        return np.array([S,E,I,H,F]), res
    Is = np.logspace(-10, -1, n)
    vals = np.array([build(I)[1] for I in Is])
    roots = []
    for k in range(n-1):
        if np.sign(vals[k]) != np.sign(vals[k+1]):
            lo, hi = Is[k], Is[k+1]
            for _ in range(80):
                m = np.sqrt(lo*hi)
                if np.sign(build(m)[1]) == np.sign(build(lo)[1]): lo = m
                else: hi = m
            x, _ = build(np.sqrt(lo*hi))
            if x[0] > 0 and x[0] < 1:
                ev = np.linalg.eigvals(jac(x, p, B))
                roots.append((x, ev))
    return roots
