# Deployment brief

Hand this file to Claude Code in VS Code (or follow it yourself). Goal: publish this folder as the public GitHub repo `Clint-Os/kenya-bvd-simulator` and serve `index.html` through GitHub Pages at <https://clint-os.github.io/kenya-bvd-simulator/>.

## State of play

- A repo of that name already exists. Its first upload went through the web page and lost the folders (`src/`, `data/`, `docs/`, `tests/`). Treat the local folder as the source of truth.
- `index.html` is a single self-contained page at the repo root. There is no build step. The only external request is the Google Fonts stylesheet.

## Steps

1. **Check tools.** `git --version`, `gh auth status` (or sign in to GitHub from VS Code), `python --version`.
2. **Get the repo locally.** `gh repo clone Clint-Os/kenya-bvd-simulator` into a new folder. If it is empty or has only a few flat files, that is expected.
3. **Sync the bundle into the clone.** Copy everything from this project folder except `.git`, `.venv`, `figures/` and `__pycache__`. Preview first with `rsync -avn --delete --exclude .git ./ ../kenya-bvd-simulator/`, read the list of deletions (they should only be stray flat copies of files that now live in folders), then run it without `-n`.
4. **Test.** `python -m venv .venv`, activate it, `pip install -r requirements.txt`, `python -m pytest -q tests`. All tests must pass. One test fails if the defaults embedded in `index.html` drift from `data/supplies_defaults.json`.
5. **Look at the page.** Open `index.html` in a browser. Check the four panels render, the sliders respond, and there are no console errors. Check at phone width.
6. **Commit and push to `main`.** Use a clear message. Do not force-push. Never commit tokens, `.venv/` or `figures/`.
7. **Enable GitHub Pages.** `gh api -X POST repos/Clint-Os/kenya-bvd-simulator/pages -f "source[branch]=main" -f "source[path]=/"`. If it answers 409 (already enabled), run `gh api repos/Clint-Os/kenya-bvd-simulator/pages --jq .html_url` and confirm the source is `main` at `/`. The web route is Settings → Pages → Deploy from a branch.
8. **Verify.** Wait about two minutes, then `curl -sI https://clint-os.github.io/kenya-bvd-simulator/` should return `200`. Open the live URL and check the article link in the README resolves (`docs/Bundibugyo-isolation-bed-capacity-article.pdf`).
9. **Tidy the repo page.** `gh repo edit Clint-Os/kenya-bvd-simulator --description "Bundibugyo (Ebola) outbreak simulator for Nairobi and Mombasa: isolation beds and facility stock planning" --homepage https://clint-os.github.io/kenya-bvd-simulator/ --add-topic ebola --add-topic kenya --add-topic epidemiology --add-topic public-health --add-topic simulation`.
10. **Report back** with the live URL, the commit hash and anything that failed.

## Later updates

Edit, run the tests, commit, push. Pages rebuilds on every push to `main`. Bump the version line in the README when behaviour changes.

## If something goes wrong

- Pages shows 404: confirm `index.html` is at the repo root on `main` and the Pages source is `main` at `/`.
- Fonts look different: the Google Fonts request was blocked. The page falls back to system fonts, which is fine.
- Push rejected: pull first (`git pull --rebase`), re-run the tests, then push.
