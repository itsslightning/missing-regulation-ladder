# Deploying to Streamlit Community Cloud

Everything the deployed app needs is already committed, so this is a
configuration step rather than a build step. **It has to be done from your own
account** — authorising a Community Cloud app against your GitHub is yours to
click, not something to automate on your behalf.

## Why it should just work

- The entry point is `app/streamlit_app.py`.
- `requirements.txt` is the deploy floor (Community Cloud installs from it;
  `uv.lock` stays as the exact pin for reproducing an analysis run).
- The app **reads only committed parquet files** in `data/processed/` and does
  no network I/O and no computation at runtime, so there is nothing to
  configure: no secrets, no API keys, no database.
- Nothing under `data/raw/` or `data/restricted/` is needed at runtime, which
  is why the restricted sources never have to leave this machine.

## Steps

1. Go to https://share.streamlit.io and sign in with GitHub.
2. **New app** → **Deploy a public app from a repo**.
3. Fill in:
   - Repository: `itsslightning/missing-regulation-ladder`
   - Branch: `main`
   - Main file path: `app/streamlit_app.py`
4. Deploy. First build takes a few minutes while it installs pandas, pyarrow,
   plotly, scipy and statsmodels.

Optionally set a custom subdomain under **Advanced settings** — something like
`missing-regulation-ladder` — so the URL is quotable in an application.

## After it is live

Add the URL to the README badge line at the top and to the sibling
`scz-target-prioritization` README, so the two projects cross-link in both
directions.

## If the build fails

- **Memory.** Community Cloud gives ~1 GB. The committed tables total ~35 MB
  and the largest single one is `smr_exposures.parquet` at 17 MB, so this
  should be comfortable. If it ever becomes tight, `detection_long.parquet`
  (12 MB) is the one to drop — only the SCHEMA heatmap reads it, and it could
  be pre-aggregated to just the 32 SCHEMA genes.
- **A missing table.** The pages degrade with an explicit error rather than a
  traceback if a parquet is absent; run `uv run python scripts/run_all.py` and
  commit the regenerated `data/processed/`.
- **Package resolution.** `requirements.txt` pins floors, not exact versions,
  deliberately — Community Cloud resolves against its own Python. If a
  resolution conflict appears, pin the offending package to the version in
  `uv.lock`.
