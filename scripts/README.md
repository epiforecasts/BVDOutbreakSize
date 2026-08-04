# Scripts

Helper scripts for the analysis.
Each runs against a specific Julia project.
Run it as `julia --project=<project> scripts/<name>.jl` from the repository root, unless noted otherwise below.

## Data pipeline (`--project=scripts`)

These are the scripts to reach for when a new SitRep lands, run in this order.
See `data/README.md` for the full data-update procedure, including the manual transcription steps these scripts do not cover.

| Script | What it does |
| --- | --- |
| `check_new_sitreps.jl` | Lists INSP SitReps not yet in `data/insp_sitrep_scanned.csv`. Exits non-zero if any are missing. |
| `download_sitreps.jl` | Downloads the INSP SitRep PDFs into `data/sitrep_pdfs/` (git-ignored). Also `task download-sitreps`. |
| `confirm_insp_data.jl` | Regenerates the confirmed-case and confirmed-death TOML blocks for `data/observations.toml` from the INRB-UMIE mirror, and cross-checks them against the scanned totals. Also `task confirm-data`. |
| `refresh_releases.jl` | Pulls each tagged results release's headline estimate into `data/released_estimates.csv`. Also `task refresh-releases`. |

`check_new_sitreps.jl`, `download_sitreps.jl` and `confirm_insp_data.jl` need no Julia packages beyond `Downloads`.
`refresh_releases.jl` also needs the `gh` CLI, authenticated against the repo.

### Onset-curve digitiser

`digitize_onset_curve.jl` digitises the symptom-onset epidemic-curve figure from the analytique SitRep PDFs into `data/onset_curve_scanned.csv`.
The figure is a raster bar chart with no data table.
It needs no Julia packages beyond stdlib, so run it as `julia scripts/digitize_onset_curve.jl [pdf_dir] [out_csv]`.
`digitize_onset_curve.py` is a byte-identical Python port for the automated data-updater, which has Python but not Julia access.
Run it with `uv run scripts/digitize_onset_curve.py`, which fetches Pillow and numpy from its PEP 723 inline metadata.
Both need poppler's `pdfimages`, `pdftotext` and `pdfinfo` on `PATH` (`apt install poppler-utils`, or `brew install poppler`).
Download the PDFs first with `download_sitreps.jl`.

## Publishing the results (`--project=docs`)

| Script | What it does |
| --- | --- |
| `score_releases.jl` | Scores every past release's saved forecasts against the now-observed data and refreshes the scoring and per-release R_T/C_T/R0 overlay CSVs. |
| `standalone_report.jl` | Lifts the rendered Vitepress analysis page into one self-contained offline HTML file. |

Both also run under `--project=.`.
CI uses `--project=docs` because that environment is already instantiated at that point in the build.

## Entry points and reproduction

| Script | Project | What it does |
| --- | --- | --- |
| `run.jl` | `.` | Regenerates the published results by running the analysis and sensitivity pages. |
| `reproduce.jl` | none | Bootstraps a full reproduction from a fresh clone; run with `curl -fsSL https://raw.githubusercontent.com/epiforecasts/BVDOutbreakSize/main/scripts/reproduce.jl \| julia`. |
| `backfill_forecasts.jl` | none | Reconstructs the one-week-ahead forecast each past release made but never saved, for `score_releases.jl` to score. Needs only `Dates` itself: it checks out each release tag into its own worktree and runs there under that tag's own project. |
| `backfill_drivers/driver_v1.0.0.jl`, `driver_v1.1.0.jl` | that tag's own `docs` | Standalone drivers for the two release tags whose model predates the fit registry. Run by `backfill_forecasts.jl`, not directly. |

## Formatting

`run_formatter.sh` runs JuliaFormatter over `src/`, `test/`, `docs/` and `scripts/` from its own isolated `test/formatter/` sub-environment.
The pre-commit hook and `task lint` both call it directly.
It takes no project flag.

## Gradient-cost diagnostic (`--project=.`)

`bench_convolve.jl` times the Mooncake gradient of the renewal convolution, scalar loop against a vectorised rewrite.
Its finding is cited from a docstring in `src/renewal.jl`, which justifies keeping the simpler scalar loop.
Run it with `julia --project=. scripts/bench_convolve.jl` if you touch that code path.
