# Phase 2 — per-dataset (individual fit) by-release scoring and presentation

Hand-off for the archie `add-scoring` agent. Branch `forecast-scoring` (pushed).
Do the compute-heavy work here (fits, docs build, backfill) on archie, checking
`~/.claude/hooks/compute-budget.sh` and running Julia-heavy steps serially
(parallel `Pkg.instantiate`/formatter contend on the shared depot and stall).

## What is already done (Phase 1, committed on this branch)

Five commits, `341a9139`..`1b962b66`:

- **C1** forecasts saved as release assets: `forecast_archive` (src/forecast.jl)
  writes `output/forecast.csv` at horizons 7/14/21/28 (analysis.jl saving
  block) + `forecast_validation.csv` (sensitivity.jl). Long schema
  `made_date,horizon,target_date,stream,draw,value`; streams `confirmed
  cases`/`confirmed deaths`/`recovered` (incident new-over-horizon) + `isolation
  beds` (occupancy level).
- **C2** `src/scoring.jl` on **ScoringRules.jl** (unregistered git dep in
  Project.toml + docs/Project.toml + test/Project.toml via `[sources]`):
  `crps_sample` (`crps(samples, obs)`), `log_crps_sample` (log1p then crps),
  `score_draws`, `forecast_score_summary`. Tests in test/test_scoring.jl.
- **C3** `scripts/score_releases.jl`: pulls each `results-v*` release's
  `forecast.csv`, scores vs current obs + a Poisson persistence baseline, writes
  `data/forecast_scores.csv`, `data/forecast_overlay.csv`,
  `data/rt_by_release.csv`. Excluded from the fit-cache hash
  (`FIT_DATA_EXCLUDE` in docs/fits/registry.jl); run in CI before render.
- **C4/C5** sensitivity.jl sections "Forecast scoring across releases" and
  "Reproduction number by release" (`plot_forecast_overlay`,
  `forecast_score_summary`, `plot_estimate_evolution` with new `refline`).
- **C6** `scripts/backfill_forecasts.jl` (not run): worktree per release tag,
  refit, forecast, publish to a `forecasts-backfill` release. Pre-registry
  releases v1.4.0-v1.7.0 flagged for manual handling.

## Phase 2 goal

Save and present the individual single-stream fits per release, the per-dataset
analogue of the existing "Rt by dataset" (`plot_rt_streams`) and "outbreak size
by each data stream" (`streams_table`) sections, and score them so we compare,
per dataset, the persistence baseline vs the individual fit vs the joint. Each
individual model forecasts only its own dataset (the deaths fit forecasts
deaths, not cases). Save R_T + outbreak size + thinned draws on release, and
present one panel per dataset (faceted).

## Key findings (grounded — do not re-derive)

Single-stream fits are loaded in docs/examples/_setup.jl:127-132 as
`chn_exports, chn_deaths, chn_cases, chn_confirmed, chn_confirmed_deaths,
chn_treatment`; their `C_T` posteriors as `posterior_C_*` (_setup.jl:144-149).

Because every single-stream model includes the latent submodel un-prefixed
(`latent ~ to_submodel(_latent(...), false)`), these are **top-level** on all
single-stream chains: `:C_T`, `:R_T`, `:r`, `:cumulative_infections`,
`rt_state.log_R0/sigma_rw/z/intervention_effect`. So per-stream R_T and size
need no new machinery (`chn[:R_T]`, `chn[:C_T]`; `reconstruct_rt` already works,
that is how `plot_rt_streams` renders).

Stream-specific forecast deterministics are **nested** under the submodel
binding, so `forecast_reported` (which reads top-level `expected_reports_T`,
`expected_deaths_T`, `k`, ...) cannot run on a single-stream chain. Bindings
(src/models/joint.jl):

| model | stream-specific binding | dispersion |
|---|---|---|
| cases_only | `cases_state.*` | `dispersion_state.k` |
| deaths_only | `deaths_state.*` | `dispersion_state.k` |
| confirmed_only | `confirmed_state.*` | `dispersion_state.k` |
| confirmed_deaths_only | `confirmed_deaths_state.*` | `dispersion_state.k` |
| treatment_only | `treatment_state.*` | `dispersion_state.k` |
| exports_joint_only | `exports_state.*` | `dispersion_state.k` |

## Tasks

### P2.0 — introspect exact nested field names (archie, cheap)
Run a 20-draw fit of each single-stream model (`nuts_sample(...; samples=20,
chains=1, check_model=false)`), print `keys(chn)`, and record the exact nested
names for each stream's cut-off expected count and cumulative trajectory (e.g.
`confirmed_state.expected_confirmed_T`, `cases_state.expected_reports_T`,
`treatment_state.expected_bed_demand_T`/`bed_capacity`,
`deaths_state.expected_deaths_T`, `exports_state.*`). The joint's field names
are the un-prefixed versions (see `forecast_reported`, src/forecast.jl:202+).

### P2.1 — `forecast_stream` (src/forecast.jl)
Add a per-stream forecaster that projects ONE stream from a chain, reusing the
existing helpers `_evolving_rates`, `_geometric_new`/`_evolving_new`,
`_nb_rand`, `_daily_at_cutoff` (all read the top-level renewal/rt-walk params,
which single-stream chains carry). Read the stream-specific expected count,
dispersion and (where relevant) cumulative trajectory by trying the nested name
first and the top-level name second, so the same function serves both the
individual fit and the joint. Signature:

```
forecast_stream(chn, stream::Symbol; horizon, obs_value) -> Vector  # draws
```

`stream ∈ (:reported_cases, :suspected_deaths, :confirmed_cases,
:confirmed_deaths, :isolation_beds, :exports)`. Incident streams return the new
count over the horizon (matching `forecast_archive`'s incident convention);
`:isolation_beds` returns the occupancy level. Add unit tests (synthetic chains
mirroring test/test_forecast.jl, testing both a nested-name and a top-level-name
chain).

### P2.2 — save per-stream on release (analysis.jl saving block)
For each single-stream chain and the joint, write:
- `output/stream_estimates.csv`: one row per fit with R_T (median + 30/60/90
  from `chn[:R_T]`) and size C_T (from `chn[:C_T]`), plus `fit` (the stream
  name / "joint").
- `output/stream_forecasts.csv`: each fit's forecast of ITS OWN dataset via
  `forecast_stream`, in the `forecast_archive` long schema plus a `fit` column;
  the joint's forecasts of the shared streams tagged `fit="joint"`.
- `output/stream_draws.csv`: thinned R_T and C_T draws per fit.
All attach via the existing `output/*` release glob (docs.yml).

### P2.3 — aggregate + score across releases (scripts/score_releases.jl)
- Pull each release's `stream_estimates.csv` → `data/rt_by_release_by_stream.csv`
  and `data/size_by_release_by_stream.csv` (columns `release,date,fit,median,
  lo30,hi30,lo60,hi60,lo90,hi90`).
- Pull each release's `stream_forecasts.csv` → score each fit's forecast of its
  own dataset vs the now-observed truth AND vs the persistence baseline; extend
  `data/forecast_scores.csv` and `data/forecast_overlay.csv` with a `fit` column
  (`baseline` / `individual:<stream>` / `joint`). Per (stream, horizon) we then
  have baseline, the individual fit, and (for the shared streams the joint
  forecasts) the joint. Add the three new CSVs to `FIT_DATA_EXCLUDE`.

### P2.4 — presentation (sensitivity.jl, faceted per dataset)
- "Reproduction number by release and dataset": one panel per dataset (reuse
  `plot_estimate_evolution` per stream with `refline=1.0`, or add a faceted
  `plot_evolution_by_group`), from `rt_by_release_by_stream.csv`.
- "Outbreak size by release and dataset": same, from
  `size_by_release_by_stream.csv`.
- Extend `forecast_score_summary` to group by (stream, horizon, fit) so the
  table shows baseline vs individual vs joint (CRPS, log-CRPS, coverage, bias,
  relative skill). Colour the forecasts-vs-now overlay by `fit`.
Terse prose only (plots and tables, matching the current sensitivity page).

### P2.5 — backfill extension (scripts/backfill_forecasts.jl)
Extend the per-tag driver to also run the single-stream fits (they are in
`build_fit_specs`) and write `stream_estimates.csv` + `stream_forecasts.csv`, so
the per-dataset history is backfilled alongside the joint.

## Verification (archie)
- P2.0 introspection confirms field names before P2.1.
- Unit tests: `task test` (TestItemRunner; note it scans sibling worktrees, use
  an explicit test path — see the local-test-tooling note).
- Run `scripts/score_releases.jl` and confirm the new CSVs; run a full docs
  build and eyeball the faceted sections render (host
  `python3 -m http.server --directory docs/build/1`).
- Format with `./scripts/run_formatter.sh` (pinned) before committing.

## Conventions
- UK English, ≤80 cols, no trailing whitespace, sentence-case comments (no
  ALL-CAPS), one sentence per line in prose.
- Commit per unit; do not push to main; keep the branch current with
  origin/main.
