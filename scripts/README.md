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
| `confirm_insp_data.jl` | Cross-checks the scanned confirmed-case and confirmed-death totals against the INRB-UMIE mirror, and reports the dates each source carries alone. Also `task confirm-data`. |
| `scan_zone_tableau2.jl` | Scans Tableau 2 of the PDFs (per-health-zone confirmed cases and deaths within each province) into the `[zone_confirmed_history]` and `[zone_death_history]` blocks, admitting a vintage only when its zone rows partition the committed province cumulatives. Also `task zone-tableau2`. |
| `confirm_zone_data.jl` | Cross-checks the two zone blocks against the INRB-UMIE mirror's per-zone CSVs and lists every disagreement. Also `task confirm-zone-data`. |
| `build_health_zones.py` | Writes `data/health_zones.csv` and `data/health_zones.geojson` from the INRB-UMIE health-zone GeoJSON (Python standard library only; pass the GeoJSON path). Also `task health-zones`. |
| `refresh_releases.jl` | Pulls each tagged results release's headline estimate into `data/released_estimates.csv`. Also `task refresh-releases`. |

`check_new_sitreps.jl`, `download_sitreps.jl` and `confirm_insp_data.jl` need no Julia packages beyond `Downloads`.
`scan_zone_tableau2.jl` and `confirm_zone_data.jl` need only `TOML`, `Printf` and `Downloads`.
`refresh_releases.jl` also needs the `gh` CLI, authenticated against the repo.

### Onset-curve digitiser

`digitize_onset_curve.jl` digitises the symptom-onset epidemic-curve figure from the analytique SitRep PDFs into `data/onset_curve_scanned.csv`.
The figure is a raster bar chart with no data table.
It needs no Julia packages beyond stdlib, so run it as `julia scripts/digitize_onset_curve.jl [pdf_dir] [out_csv]`.
`digitize_onset_curve.py` is a byte-identical Python port for the automated data-updater, which has Python but not Julia access.
The Julia script is the reference: `data/onset_curve_scanned.csv` is its output, and a disagreement between the two is a fault in the port.
`test/test_onset_digitiser.jl` runs both against the committed file when the PDFs are present.
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

`run_formatter.sh` runs Runic over `src/`, `test/`, `docs/`, `scripts/`, `benchmark/` and `ext/` from its own isolated `test/formatter/` sub-environment.
`task format` and `task lint` both call it directly.
It takes no project flag.
The pre-commit hook does not: pre-commit builds the hook its own environment from the Runic version in `.pre-commit-config.yaml`, and `test/package/CodeFormatting.jl` checks that version against the pin in `test/formatter/Project.toml`.

## Developer diagnostics (`--project=.`)

These are not wired into CI, the Taskfile or any workflow.
Reach for them when a fit misbehaves or a hot path changes.
Each documents its own invocation in its header comment.

| Script | What it does |
| --- | --- |
| `bench_convolve.jl` | Times the Mooncake gradient of the renewal convolution, scalar loop against a vectorised rewrite. Its finding is cited from a docstring in `src/renewal.jl`, which justifies keeping the simpler scalar loop. Run it if you touch that code path. |
| `bench_discretise.jl` | Times the Mooncake gradient of the censored delay discretisation. Documents the reasoning behind the CDF-difference form in `src/renewal.jl`. |
| `summarise_chain.jl` | Prints posterior summaries and fit diagnostics for a saved chain. |
| `prior_vs_posterior.jl` | Samples the prior and sets it beside a saved posterior, so a parameter the data does not inform is visible. |
| `zone_fit_report.jl` | Builds a self-contained HTML fit report for the health-zone model from a joint parent chain: prior predictive check, simulation-based recovery, sampler diagnostics, posterior predictive checks, ranking and map, under `logs/zone_report/`. Runs under `--project=docs`; see its header for the flags. |
| `fit_joint_stream.jl` | Fits the full joint model outside the fit cache and CI, streaming live progress to `logs/joint_fit.log` so a long fit can be watched with `tail -f`. Writes the chain to `logs/joint_chain.jls`. The name refers to the streamed progress, not to fitting a single data stream. |
