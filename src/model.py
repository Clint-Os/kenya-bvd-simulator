"""
Bundibugyo virus disease (BVD) transmission model for Nairobi and Mombasa.

Structure: two-group (formal / informal settlement) SEIHFR model after
Legrand et al. (2007), with a hospitalisation (isolation) rate that falls as
isolation beds fill up. Bed capacity B is the control / bifurcation parameter.

Time unit: days.
"""
import numpy as np
from scipy.integrate import solve_ivp

# ---------------------------------------------------------------------------
# Natural history (BDBV where available, Legrand 2007 delays otherwise)
# ---------------------------------------------------------------------------
NATURAL_HISTORY = dict(
    incubation=7.0,        # days; median 7 (Wamala 2010), mean 6.3 (MacNeil 2010)
    onset_to_outcome=10.0, # days in community before death/recovery (Legrand: ~9.6-10)
    hosp_to_outcome=6.0,   # days in isolation before death/discharge (Legrand: 4.6 death)
    death_to_burial=2.0,   # days (Legrand 2007)
    cfr_hosp=0.31,         # BDBV adjusted CFR, Uganda 2007 (de Padua & Akhmetzhanov 2026)
    cfr_comm=0.47,         # crude CFR, DRC 2026 confirmed cases (WHO, Aug 2026)
    safe_burial_hosp=0.90, # share of isolation deaths given safe & dignified burial
)

# Legrand-style split of R0 by setting (Gomes et al. 2014: 0.8 / 0.4 / 0.6)
R0_SHARES = dict(community=0.8/1.8, hospital=0.4/1.8, funeral=0.6/1.8)

# ---------------------------------------------------------------------------
# City profiles (KNBS 2019 KPHC + cited secondary sources)
# ---------------------------------------------------------------------------
CITIES = {
    "Nairobi": dict(
        population=4_397_073,       # KNBS 2019 KPHC Vol I
        density=6_247,              # persons / km2, KNBS 2019
        household_size=2.9,         # KNBS 2019
        informal_share=0.50,        # 41% (World Bank) to 60% (UN-Habitat) -> midpoint
        isolation_beds_covid=1_085, # MoH county readiness report, Aug 2020
        onset_to_isolation_formal=4.0,    # days, assumption (on-site NPHL/KEMRI labs)
        onset_to_isolation_informal=6.0,  # days, assumption
        safe_burial_comm=0.5,
    ),
    "Mombasa": dict(
        population=1_208_333,
        density=5_495,
        household_size=3.1,
        informal_share=0.60,        # "nearly 60%" (Standard, 2026, citing studies)
        isolation_beds_covid=404,
        onset_to_isolation_formal=5.0,    # +1 day: no designated Ebola lab in county
        onset_to_isolation_informal=7.0,
        safe_burial_comm=0.5,
    ),
}

INFORMAL_CONTACT_MULTIPLIER = 1.6   # assumption: crowding, shared sanitation


def build_params(city, R0_ref=1.5, beds=None, contact_mult=INFORMAL_CONTACT_MULTIPLIER,
                 overrides=None):
    """Return a parameter dict for one city.

    Transmission coefficients are calibrated so that a *reference* homogeneous
    population (contact multiplier 1, unlimited beds, onset-to-isolation 5 days,
    no safe burial) has R0 = R0_ref split by the Legrand shares. City-specific
    R0 then follows from that city's structure.
    """
    nh = dict(NATURAL_HISTORY)
    cp = dict(CITIES[city])
    if overrides:
        for k, v in overrides.items():
            (cp if k in cp else nh)[k] = v

    sigma = 1 / nh["incubation"]
    gI = 1 / nh["onset_to_outcome"]
    gH = 1 / nh["hosp_to_outcome"]
    gF = 1 / nh["death_to_burial"]

    # reference calibration
    rho_ref = 1 / 5.0
    t_I = 1 / (rho_ref + gI)
    p_hosp = rho_ref / (rho_ref + gI)
    p_comm_death = gI * nh["cfr_comm"] / (rho_ref + gI)
    betaI = R0_ref * R0_SHARES["community"] / t_I
    betaH = R0_ref * R0_SHARES["hospital"] / (p_hosp / gH)
    funeral_inflow = p_comm_death + p_hosp * nh["cfr_hosp"]  # no safe burial at reference
    betaF = R0_ref * R0_SHARES["funeral"] / (funeral_inflow / gF)

    N = cp["population"]
    q = cp["informal_share"]
    return dict(
        city=city, N=N,
        Ng=np.array([N * (1 - q), N * q]),
        c=np.array([1.0, contact_mult]),
        rho=np.array([1 / cp["onset_to_isolation_formal"], 1 / cp["onset_to_isolation_informal"]]),
        sigma=sigma, gI=gI, gH=gH, gF=gF,
        cfr_c=nh["cfr_comm"], cfr_h=nh["cfr_hosp"],
        s_c=cp["safe_burial_comm"], s_h=nh["safe_burial_hosp"],
        betaI=betaI, betaH=betaH, betaF=betaF,
        B=cp["isolation_beds_covid"] if beds is None else beds,
    )


def bed_access(Htot, B, k=0.02):
    """Fraction of normal admission rate still possible. Smoothly -> 0 as H -> B."""
    x = (B - Htot) / max(B, 1e-9)
    # softplus-smoothed max(0, x), capped at 1
    return np.minimum(1.0, k * np.log1p(np.exp(x / k)))


def rhs(t, y, p):
    S, E, I, H, F, R, D, C = y.reshape(8, 2)
    Ng, c = p["Ng"], p["c"]
    Ctot = np.sum(c * Ng)

    # response phase (from day t_resp): fewer community contacts, faster
    # isolation, more safe and dignified burials
    on = t >= p.get("t_resp", 1e9)
    cut = p.get("contact_cut", 0.0) if on else 0.0
    rho = p["rho"] * (p.get("isolation_speedup", 1.0) if on else 1.0)
    s_c = p.get("s_c_resp", p["s_c"]) if on else p["s_c"]

    lam_I = (1 - cut) * p["betaI"] * np.sum(c * I) / Ctot
    lam_F = p["betaF"] * np.sum(c * F) / Ctot
    lam_H = p["betaH"] * np.sum(H) / p["N"]
    inf = (c * (lam_I + lam_F) + lam_H) * S

    # admissions fall as isolation beds fill; overflow can partly be absorbed by
    # community care centres / supported home isolation (share omega)
    phi = bed_access(np.sum(H), p["B"])
    adm = rho * (phi + (1 - phi) * p.get("omega", 0.0)) * I
    out_I = p["gI"] * I
    out_H = p["gH"] * H
    deaths_c = out_I * p["cfr_c"]
    deaths_h = out_H * p["cfr_h"]

    dS = -inf
    dE = inf - p["sigma"] * E
    dI = p["sigma"] * E - adm - out_I
    dH = adm - out_H
    dF = deaths_c * (1 - s_c) + deaths_h * (1 - p["s_h"]) - p["gF"] * F
    dR = out_I * (1 - p["cfr_c"]) + out_H * (1 - p["cfr_h"]) + p["gF"] * F \
        + deaths_c * s_c + deaths_h * p["s_h"]
    dD = deaths_c + deaths_h
    dC = p["sigma"] * E          # cumulative symptomatic cases
    return np.concatenate([dS, dE, dI, dH, dF, dR, dD, dC])


def simulate(p, days=365, seed_I=(0, 5), seed_E=(0, 5)):
    y0 = np.zeros((8, 2))
    y0[0] = p["Ng"]
    y0[1] = seed_E
    y0[2] = seed_I
    y0[0] -= y0[1] + y0[2]
    sol = solve_ivp(rhs, (0, days), y0.ravel(), args=(p,), t_eval=np.arange(days + 1),
                    method="LSODA", rtol=1e-7, atol=1e-6, max_step=0.5)
    Y = sol.y.reshape(8, 2, -1)
    names = "S E I H F R D C".split()
    out = {n: Y[i].sum(0) for i, n in enumerate(names)}
    out["t"] = sol.t
    out["I_informal"] = Y[2, 1]
    out["H_total"] = Y[3].sum(0)
    out["bed_access"] = bed_access(out["H_total"], p["B"])
    return out


def R0(p, unlimited_beds=True):
    """Next-generation-matrix R0 at the disease-free state."""
    Ng, c, rho = p["Ng"], p["c"], p["rho"]
    Ctot = np.sum(c * Ng)
    n = 2
    # state order per group: E, I, H, F  -> index 4*g + j
    Fm = np.zeros((8, 8))
    Vm = np.zeros((8, 8))
    for g in range(n):
        e = 4 * g
        for k in range(n):
            Fm[e, 4 * k + 1] = Ng[g] * c[g] * p["betaI"] * c[k] / Ctot
            Fm[e, 4 * k + 3] = Ng[g] * c[g] * p["betaF"] * c[k] / Ctot
            Fm[e, 4 * k + 2] = Ng[g] * p["betaH"] / p["N"]
        r = rho[g] if unlimited_beds else 0.0
        Vm[e, e] = p["sigma"]
        Vm[e + 1, e] = -p["sigma"]
        Vm[e + 1, e + 1] = r + p["gI"]
        Vm[e + 2, e + 1] = -r
        Vm[e + 2, e + 2] = p["gH"]
        Vm[e + 3, e + 1] = -p["gI"] * p["cfr_c"] * (1 - p["s_c"])
        Vm[e + 3, e + 2] = -p["gH"] * p["cfr_h"] * (1 - p["s_h"])
        Vm[e + 3, e + 3] = p["gF"]
    K = Fm @ np.linalg.inv(Vm)
    return float(np.max(np.abs(np.linalg.eigvals(K))))


def R0_by_setting(p):
    """Split R0 into community / hospital / funeral parts (each setting alone)."""
    base = dict(p)
    parts = {}
    for name, keys in [("community", ("betaH", "betaF")), ("hospital", ("betaI", "betaF")),
                       ("funeral", ("betaI", "betaH"))]:
        q = dict(base)
        for kk in keys:
            q[kk] = 0.0
        parts[name] = R0(q)
    return parts
