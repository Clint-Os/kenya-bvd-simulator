# Project notes for Claude

## What this is

An interactive planning tool and model code for **Bundibugyo virus disease (Ebola) spread in Nairobi and Mombasa**, built to size isolation beds and, now, facility stock. It is a planning tool, not a forecast. Owner: Clint (Kenyan health-data, regulatory and digital-health consultancy). Visual identity: deep teal `#0F3B40`, terracotta `#D97757`. He prefers concise, high-impact output and production-ready files.

## Layout

`index.html` single-file simulator (no build) · `src/` Python model · `data/` results and `supplies_defaults.json` · `docs/` article PDF, literature synthesis, `DEPLOYMENT.md` · `tests/`. See `README.md`.

## Commands

```bash
pip install -r requirements.txt
python -m pytest -q tests
cd src && python figures.py            # article figures into ../figures (git-ignored)
cd src && python containment_sweep.py  # ~4 min
cd src && python supplies.py --city Mombasa --beds 26
```

## Rules that matter

- **Two copies of the stock-planner defaults** exist: `data/supplies_defaults.json` and the `<script id="supplies-defaults">` block in `index.html`. Change both. `tests/test_supplies.py` fails if they differ.
- **The page (JS) and `src/supplies.py` must compute identical numbers.** Change both together and keep the hand-checked tests.
- **Label assumptions as assumptions.** No source gives per-patient PPE quantities. Practice rules come from the WHO IPC guideline for Ebola and Marburg diseases (2025); item lists from the Kenya National IPC Guidelines (2010), the WHO Disease Commodity Package (2018) and the KMA circular (2026). Quantities and care-model drivers are assumptions until a Kenyan SOP replaces them.
- **Several people or sessions may edit this folder.** Run `git status` first. If there are uncommitted changes you did not make, read the diff and merge. Never overwrite a file wholesale without reading it.
- Beds are a hard cap in the model, so "every bed full" is the ceiling for any stock estimate.
- Do not commit tokens, `.venv/`, `figures/` or `__pycache__/`.
- Commit messages end with the attribution lines given in the session.

## Decisions already made

Nairobi 57 public beds (Nairobi Hospital's 130 private beds excluded). Mombasa has 4 beds at Port Reitz; Utange is the free variable. Two index cases. R₀ 2.0 reference. Containment means fewer than 20,000 cumulative cases in 400 days. Both Hopf analyses are in the article (bed capacity; surveillance delay, read as overshoot, not literal cycles).

## Open items

1. Replace stock-planner assumptions with Kenyan SOP figures when found (Clint has the MoH 2010 policy, which has no quantities).
2. Next phase: specific pharmaceuticals and further non-pharmaceutical supplies.
3. Optional costing: `reference_prices_usd` in the defaults holds WHO-quoted unit costs, not yet used.
4. Technical report doc and slide deck still carry first-draft numbers.
5. Verify the literature citations flagged as unverified in the article.
