# Kenya BVD Outbreak Simulator

An interactive model of **Bundibugyo virus disease (Ebola) spread in Nairobi and Mombasa**, built to support isolation-bed planning in public hospitals. It is a **planning tool, not a forecast**.

**Live simulator:** https://clint-os.github.io/kenya-bvd-simulator/
**Article (PDF):** [`docs/Bundibugyo-isolation-bed-capacity-article.pdf`](docs/Bundibugyo-isolation-bed-capacity-article.pdf)

## What the simulator shows

1. **Will two imported cases take off?** Branching-process probability that two index cases start a self-sustaining outbreak, live in the reproduction number and the overdispersion parameter *k*.
2. **How many isolation beds does Mombasa need?** Stochastic containment results for Mombasa (Port Reitz's 4 beds plus a county unit such as Utange) at different detection delays, with a 90 / 95 / 99% containment target.
3. **Watch the outbreak unfold.** A live deterministic two-group SEIHFR model of both cities with every assumption adjustable.
4. **What should an isolation unit stock?** A facility stock planner that turns the expected number of patients in isolation into PPE and infection-control quantities (coveralls, gloves, respirators, face shields, boots, waste bags, hand rub, body bags, chlorine, water). Adjust staff per entry, entries per patient, the stocking window and a safety buffer, edit any per-item quantity, enter your stock on hand to see days of cover, and download the list as CSV.

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
docs/             Article PDF, literature synthesis, DEPLOYMENT.md (hand to Claude Code in VS Code)
src/              Model code (Python)
  model.py          Deterministic two-group SEIHFR model, R0 by next-generation matrix
  stoch.py          Stochastic tau-leap engine with a hard bed cap
  containment_sweep.py   Containment probability vs beds; Utange sizing
  bifurcation.py    Bed-capacity Hopf bifurcation and equilibria
  analysis.py       Critical-bed search
  supplies.py       Facility stock planner (PPE and IPC consumables)
  figures.py        Article figures
data/             Result files used by the simulator and the article
  supplies_defaults.json   Defaults for the stock planner (also embedded in index.html)
tests/            Sanity checks on the headline numbers and the stock planner
CLAUDE.md         Project notes and rules for Claude sessions
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
python supplies.py --city Mombasa --beds 26   # stock list for a unit; add --full-unit for the ceiling
python containment_sweep.py        # ~4 min; writes ../data/containment_sweep.json
```

The sweep reproduces the Utange bed thresholds to within Monte Carlo and grid noise (a few beds; the 99% target is the least stable). The thresholds quoted in the article are in `data/utange_thresholds.json`.

## Stock planner: how the numbers are built

Occupied beds → red-zone entries (visits per patient per day ÷ patients seen per entry, with at least one team entry per shift while any patient is in) → staff entries (× staff per entry) plus support entries (cleaners, waste, burial) → one full PPE set per entry. Following the WHO 2025 IPC guideline for Ebola and Marburg diseases: coverall, head-and-neck covering, a medical mask with eye protection, two pairs of nitrile gloves with the outer pair renewed between patients, and an apron changed between patients; cleaners and burial teams wear heavy-duty outer gloves and aprons. N95 respirators are only for aerosol-generating procedures, so their default is zero. Body bags follow expected deaths (31% case fatality, 6-day stay). Chlorine follows the Kenya MoH rule of fresh 0.5% solution every 8 hours plus 0.05% for hands and gloves. Water is per occupied bed-day. Ward fixtures are counted at peak occupancy. Stock to hold is the busiest run of *N* days (your resupply lead time) plus a safety buffer. Beds are a hard cap, so "every bed full" is the most any scenario can demand.

**What is sourced and what is assumed.** Practice rules: WHO IPC guideline for Ebola and Marburg diseases (2025). Item lists: Kenya National IPC Guidelines (2010), WHO Disease Commodity Package (2018), KMA Ebola circular (2026). Water: UNICEF's 100–400 L per person per day and an IMC field figure of 150 L per patient-day. **None of these publishes per-patient quantities, so the per-entry quantities and care-model settings are planning assumptions.** Replace them with your facility SOP before ordering. The planner uses the mean trajectory, so individual outbreaks will vary around it.

`data/supplies_defaults.json` holds the defaults (including WHO-quoted unit prices for later costing). `index.html` embeds a copy, and `tests/test_supplies.py` fails if the two drift apart. `src/supplies.py` and the page compute identical numbers.

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

- **v0.3** (17 October 2026): stock planner aligned with the WHO 2025 IPC guideline (medical mask, glove and apron changes, heavy-duty gear for cleaners), MoH-listed items and ward fixtures, chlorine refresh rule.
- **v0.2** (17 October 2026): adds the facility stock planner (panel 4, `src/supplies.py`, `data/supplies_defaults.json`).
- **v0.1** (October 2026): first public version. Simulator updated to match the article (two index cases, Mombasa bed sizing, bed-capacity and surveillance-delay analyses).

## License

MIT (see `LICENSE`). Attribution appreciated.
