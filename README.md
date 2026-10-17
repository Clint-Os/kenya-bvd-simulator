# Kenya BVD Outbreak Simulator

An interactive model of **Bundibugyo virus disease (Ebola) spread in Nairobi and Mombasa**, built to support isolation-bed planning in public hospitals. It is a **planning tool, not a forecast**.

**Live simulator:** https://clint-os.github.io/kenya-bvd-simulator/
**Article (PDF):** [`docs/Bundibugyo-isolation-bed-capacity-article.pdf`](docs/Bundibugyo-isolation-bed-capacity-article.pdf)

## What the simulator shows

1. **Will two imported cases take off?** Branching-process probability that two index cases start a self-sustaining outbreak, live in the reproduction number and the overdispersion parameter *k*.
2. **How many isolation beds does Mombasa need?** Stochastic containment results for Mombasa (Port Reitz's 4 beds plus a county unit such as Utange) at different detection delays, with a 90 / 95 / 99% containment target.
3. **Watch the outbreak unfold.** A live deterministic two-group SEIHFR model of both cities with every assumption adjustable.

## Headline results (reference R₀ = 2.0, two index cases)

| | Nairobi | Mombasa |
|---|---|---|
| R₀ with isolation beds available | 1.95 | 2.18 |
| R₀ once wards are full | 4.24 | 4.41 |
| Beds assumed | 57 (public) | 4 at Port Reitz + Utange (the quantity being sized) |

With a response starting 7 days after the first isolation admission, Mombasa needs about 22 isolation beds in total (18 at Utange) for a 95% probability of containment, and about 36 for 99%. Detection speed matters most: delaying detection from 3 to 21 days roughly triples the beds required. See the article for definitions, parameters and limitations.

## Repository layout

```
index.html        The simulator (single self-contained file, no build step)
docs/             Article PDF and literature synthesis
src/              Model code (Python)
  model.py          Deterministic two-group SEIHFR model, R0 by next-generation matrix
  stoch.py          Stochastic tau-leap engine with a hard bed cap
  containment_sweep.py   Containment probability vs beds; Utange sizing
  bifurcation.py    Bed-capacity Hopf bifurcation and equilibria
  analysis.py       Critical-bed search
  figures.py        Article figures
data/             Result files used by the simulator and the article
tests/            Sanity checks on the headline numbers
```

## Running the simulator locally

Open `index.html` in any modern browser. Nothing to install. The only external resource is the Google Fonts stylesheet; the page falls back to system fonts if it is blocked.

## Running the model code

```bash
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
python -m pytest -q tests          # checks the headline numbers
cd src
python figures.py                  # regenerates article figures into ../figures
python containment_sweep.py        # ~4 min; writes ../data/containment_sweep.json
```

The sweep reproduces the Utange bed thresholds to within Monte Carlo and grid noise (a few beds; the 99% target is the least stable). The thresholds quoted in the article are in `data/utange_thresholds.json`.

## Method in one paragraph

Each city is split into formal and informal-settlement populations moving through Susceptible → Exposed → Infectious → Hospitalised (isolated) → Funeral → Removed, after Legrand et al. (2007). Transmission occurs in the community, in hospital and at unsafe burials. Admission to isolation is **capped by beds**: a case that cannot be admitted keeps transmitting, so R₀ roughly doubles once wards fill. The response (60% contact reduction, faster isolation, 90% safe burial) starts a set delay after the first admission. Outbreak probability uses a negative-binomial branching process (k = 0.18); containment is cumulative cases below 20,000 within 400 days in a binomial tau-leap simulation. Two Hopf bifurcations are analysed (bed capacity and surveillance delay). Parameters come from the Bundibugyo and wider Ebola literature and the 2019 Kenya census; all are listed in the "How the model works" table on the page.

## Publishing on GitHub Pages

1. Settings → Pages → Build and deployment → **Deploy from a branch**.
2. Branch **main**, folder **/ (root)**, Save.
3. After about a minute the page is live at the URL at the top of this file.

`.nojekyll` is included so Pages serves the files as they are.

## Status and caveats

First-draft working model. All numbers are planning scenarios, not forecasts. Intervention effect sizes are borrowed from Zaire and Sudan ebolavirus outbreaks, and several Bundibugyo parameters are substituted from the wider Ebola literature (labelled in the data table). The two cities are modelled separately, with no travel between them. Some literature citations in the article are flagged for verification before external circulation.

## Versions

- **v0.1** (October 2026): first public version. Simulator updated to match the article (two index cases, Mombasa bed sizing, bed-capacity and surveillance-delay analyses).

## License

MIT (see `LICENSE`). Attribution appreciated.
