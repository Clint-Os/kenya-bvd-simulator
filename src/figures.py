import json, numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, ListedColormap
from model import *
from bifurcation import params as bparams, equilibria
from pathlib import Path
ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / 'data'; FIGS = ROOT / 'figures'; FIGS.mkdir(exist_ok=True)

COL = {"Nairobi": "#00939A", "Mombasa": "#D97757"}
INK, MUTED, GRID = "#0F3B40", "#5B6B6D", "#E6E9E8"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10, "axes.edgecolor": MUTED,
                     "axes.labelcolor": INK, "xtick.color": MUTED, "ytick.color": MUTED,
                     "axes.spines.top": False, "axes.spines.right": False,
                     "axes.titleweight": "bold", "axes.titlecolor": INK, "axes.grid": True,
                     "grid.color": GRID, "grid.linewidth": 0.8, "lines.linewidth": 2})
RESP = dict(t_resp=21, isolation_speedup=1.5, s_c_resp=0.9)
results = {}

# 1. Same beds, same response, different outcome -------------------------------
fig, ax = plt.subplots(1, 2, figsize=(11, 4))
for city in CITIES:
    p = build_params(city, R0_ref=2.0, beds=100); p.update(RESP, contact_cut=0.6)
    s = simulate(p, days=365)
    new = np.diff(s["C"], prepend=0)
    ax[0].plot(s["t"], np.maximum(new, 1e-2), color=COL[city], label=city)
    ax[1].plot(s["t"], s["H_total"], color=COL[city], label=city)
    results[f"curve_{city}"] = dict(cases=float(s["C"][-1]), deaths=float(s["D"][-1]),
                                    peak_new=float(new.max()), peak_day=int(np.argmax(new)))
ax[0].set_yscale("log"); ax[0].set_ylim(0.1, None)
ax[0].set_title("Daily new cases (log scale)"); ax[0].set_xlabel("Days after seeding")
ax[1].axhline(100, color=INK, lw=1, ls="--"); ax[1].text(5, 104, "100 isolation beds", color=INK, fontsize=9)
ax[1].set_title("Patients in isolation"); ax[1].set_xlabel("Days after seeding")
for a in ax: a.legend(frameon=False)
fig.suptitle("Same beds (100), same response (60% contact reduction from day 21), R_ref = 2.0",
             color=INK, fontsize=11)
fig.tight_layout(); fig.savefig(FIGS / 'fig1_curves.png', dpi=200); plt.close()

# 2. Tipping map -------------------------------------------------------------
d = json.load(open(DATA / 'grid_R2.0.json'))
cuts, beds = np.array(d["cuts"]), np.array(d["beds"])
fig, ax = plt.subplots(1, 2, figsize=(11, 4.2), sharey=True)
cmap = plt.get_cmap("Greys")
for a, city in zip(ax, CITIES):
    M = np.array(d["res"][city])
    im = a.pcolormesh(beds, cuts * 100, np.maximum(M, 1), norm=LogNorm(1e2, 5e6), cmap=cmap, shading="nearest")
    a.contour(beds, cuts * 100, M, levels=[2e4], colors=[COL[city]], linewidths=2.5)
    a.set_xscale("log"); a.set_xlabel("Ebola isolation beds"); a.set_title(city)
    a.grid(False)
ax[0].set_ylabel("Community contact reduction (%)")
cb = fig.colorbar(im, ax=ax, shrink=0.85); cb.set_label("Total cases at 18 months")
fig.suptitle("Tipping boundary (coloured line): below-left the outbreak escapes, above-right it is contained",
             color=INK, fontsize=11)
fig.savefig(FIGS / 'fig2_tipping_map.png', dpi=200, bbox_inches="tight"); plt.close()

# 3. Critical beds vs contact reduction -----------------------------------------
cb_ = json.load(open(DATA / 'critical_beds.json'))
fig, ax = plt.subplots(figsize=(6.5, 4))
for city in CITIES:
    row = cb_[f"{city}|2.0"]["crit"]
    xs = [float(k) * 100 for k in row]; ys = [v["beds"] if v["beds"] else np.nan for v in row.values()]
    ax.plot(xs, ys, color=COL[city], marker="o", ms=6, label=city)
    for x, y in zip(xs, ys):
        if np.isfinite(y): ax.annotate(f"{y:.0f}", (x, y), textcoords="offset points", xytext=(6, 4), fontsize=8, color=MUTED)
ax.set_yscale("log"); ax.set_xlabel("Community contact reduction (%)"); ax.set_ylabel("Minimum isolation beds")
ax.set_title("Beds needed to keep the outbreak contained (R_ref = 2.0)"); ax.legend(frameon=False)
fig.tight_layout(); fig.savefig(FIGS / 'fig3_critical_beds.png', dpi=200); plt.close()

# 4. Bifurcation diagram ------------------------------------------------------
def city_bparams(city, R):
    p = build_params(city, R0_ref=R)
    w = p["Ng"] / p["N"]
    c_bar = float(np.sum(w * p["c"]))
    rho_days = float(np.sum(w / p["rho"]))
    return bparams(R0_ref=R, rho_days=rho_days, c=c_bar, s_c=p["s_c"])

Bgrid = np.logspace(4, -1, 160)
fig, ax = plt.subplots(1, 2, figsize=(11, 4))
hopf = {}
for city in CITIES:
    for R, ls in ((1.5, "-"), (2.0, "--")):
        bp = city_bparams(city, R)
        lead_re, cpx = [], []
        for B in Bgrid:
            eq = equilibria(bp, B)
            if not eq: lead_re.append(np.nan); cpx.append(False); continue
            x, ev = max(eq, key=lambda e: e[0][2])
            k = np.argmax(ev.real); lead_re.append(ev[k].real * 365); cpx.append(abs(ev[k].imag) > 1e-12)
        lead_re = np.array(lead_re); cpx = np.array(cpx)
        ax[0].plot(Bgrid, lead_re, color=COL[city], ls=ls, label=f"{city}, R_ref {R}")
        # Hopf: sign change of real part while complex
        idx = [i for i in range(1, len(Bgrid)) if cpx[i] and cpx[i-1] and np.sign(lead_re[i]) != np.sign(lead_re[i-1])]
        if idx:
            i = idx[0]; Bh = np.sqrt(Bgrid[i] * Bgrid[i-1])
            hopf[f"{city}|{R}"] = Bh
            ax[0].plot(Bh, 0, "o", ms=8, color=COL[city], mec="white", mew=2)
        # bed-collapse transition: complex pair turns real (oscillatory -> runaway node)
        idx2 = [i for i in range(1, len(Bgrid)) if cpx[i-1] and not cpx[i]]
        if idx2: hopf[f"{city}|{R}|node"] = float(np.sqrt(Bgrid[idx2[0]] * Bgrid[idx2[0]-1]))
ax[0].axhline(0, color=INK, lw=1)
ax[0].set_xscale("log"); ax[0].set_yscale("symlog", linthresh=1e-2)
ax[0].set_xlabel("Isolation beds per 100,000 people"); ax[0].set_ylabel("Leading eigenvalue, real part (per year)")
ax[0].set_title("Stability of the endemic equilibrium"); ax[0].legend(frameon=False, fontsize=8)
ax[0].annotate("Hopf points\n(~2,800-3,500 beds/100k)", (hopf.get("Nairobi|1.5", 3000), 0), xytext=(-150, 40),
               textcoords="offset points", arrowprops=dict(arrowstyle="-", color=MUTED), color=MUTED, fontsize=9)
# time series either side of Hopf
from scipy.integrate import solve_ivp
from bifurcation import f as bf
bp = city_bparams("Nairobi", 1.5)
for B, col, lab in ((20000, MUTED, "20,000 beds/100k (stable side)"), (100, COL["Nairobi"], "100 beds/100k (unstable side)")):
    x0 = max(equilibria(bp, B), key=lambda e: e[0][2])[0] * np.array([1, 1.05, 1.05, 1.05, 1.05])
    sol = solve_ivp(lambda t, y: bf(y, bp, B), (0, 365 * 120), x0, t_eval=np.linspace(0, 365 * 120, 8000),
                    method="LSODA", rtol=1e-9, atol=1e-14)
    ax[1].plot(sol.t / 365, np.maximum(sol.y[2] * 1e5, 1e-3), color=col, label=lab, lw=1.5)
ax[1].axhspan(1e-3, 1e5 / CITIES["Nairobi"]["population"], color=GRID, zorder=0)
ax[1].text(2, 2e-3, "fewer than 1 infectious person in the city: chains die out", fontsize=8, color=MUTED)
ax[1].set_yscale("log"); ax[1].set_ylim(1e-3, 1e4); ax[1].set_xlabel("Years"); ax[1].set_ylabel("Infectious per 100,000")
ax[1].set_title("Either side of the Hopf point (Nairobi, R_ref 1.5)"); ax[1].legend(frameon=False, fontsize=8)
fig.tight_layout(); fig.savefig(FIGS / 'fig4_bifurcation.png', dpi=200); plt.close()
results["hopf"] = hopf

# 5. R0 by setting --------------------------------------------------------------
fig, ax = plt.subplots(figsize=(6.5, 3.6))
settings = ["community", "hospital", "funeral"]; w = 0.36; y = np.arange(3)
r0tab = {}
for k, city in enumerate(CITIES):
    p = build_params(city, R0_ref=2.0); parts = R0_by_setting(p)
    r0tab[city] = dict(total=R0(p), no_beds=R0(p, False), **parts)
    vals = [parts[s] for s in settings]
    ax.barh(y + (k - 0.5) * w, vals, height=w - 0.04, color=COL[city], label=f"{city} (R0 = {R0(p):.2f})")
    for yy, v in zip(y + (k - 0.5) * w, vals): ax.text(v + 0.02, yy, f"{v:.2f}", va="center", fontsize=8, color=MUTED)
ax.set_yticks(y, [s.title() for s in settings]); ax.invert_yaxis(); ax.grid(axis="y", visible=False)
ax.set_xlabel("Reproduction number from this setting alone"); ax.legend(frameon=False, loc="lower right")
ax.set_title("Where transmission happens (R_ref = 2.0)")
fig.tight_layout(); fig.savefig(FIGS / 'fig5_r0_settings.png', dpi=200); plt.close()
results["r0"] = r0tab
json.dump(results, open(DATA / 'results.json', "w"), indent=1, default=float)
print(json.dumps(results, indent=1, default=float))
