# Contributing

Issues and pull requests are welcome at
[epiforecasts/BVDOutbreakSize](https://github.com/epiforecasts/BVDOutbreakSize).
This page covers how the project is laid out, how to run it, and the
conventions to follow when changing it.

## Repository layout

- `src/BVDOutbreakSize.jl` — the package: data loading
  (`load_observations`), NUTS sampling (`nuts_sample`), the shared
  Gauss-Legendre integrators (`integrate`, `delay_convolution`,
  `integrate_cumulative`, `integrate_exports_deaths`), summary and
  comparison tables, plotting, the no-onward-deaths projection
  (`predict_no_onward_deaths`) and forecast helpers
  (`forecast_reported`).
  The published Imperial point estimates live here as
  `REPORT_SCENARIOS`.
- `docs/examples/analysis.jl` — the Literate walkthrough that is the
  analysis.
  It defines the Turing submodels and composers, runs the fits, and
  writes every output.
  This is the main artifact.
- `docs/make.jl` — DocumenterVitepress build.
  Copies `README.md` to `index.md`, executes the literate to
  `analysis.md`, and builds the bibliography.
- `data/observations.toml` — single source of truth for observation
  data (case and death counts, traveller volumes, sources).
  Loaded via `load_observations()` and never hardcoded.
  Update this one file for a new situation report and the analysis
  picks it up.
  The literate re-binds its observation `const`s from the loaded TOML,
  so the package constants are defaults only.
- `scripts/run.jl` — regenerates published results by including the
  literate and writes CSVs to `output/`.
- `test/` — one file per feature, driven by `test/runtests.jl`.
- `external/bdbv-linelist-analysis` — git submodule, source of the
  onset-to-death delay priors.

## Running and testing

There is no Taskfile.
Use the `julia --project` commands:

```bash
# Instantiate the package environment
julia --project=. -e 'using Pkg; Pkg.instantiate()'

# Run the analysis
julia --project=. docs/examples/analysis.jl

# Regenerate the published output CSVs into output/
julia --project=. scripts/run.jl

# Run the full test suite
julia --project=. -e 'using Pkg; Pkg.test()'

# Build the docs (executes the literate, HTML in docs/build/)
julia --project=docs -e 'using Pkg; \
  Pkg.develop(PackageSpec(path=pwd())); Pkg.instantiate()'
julia --project=docs docs/make.jl
```

`test/runtests.jl` includes each `test/test_*.jl`.
To iterate on one file, run it inside a REPL after
`using BVDOutbreakSize`, or temporarily comment out the others in
`runtests.jl`.

CI runs the test suite (`.github/workflows/test.yml`) and builds the
docs, publishing `output/` as a GitHub Release on each push to `main`
(`.github/workflows/docs.yml`).

## Model architecture

The model is assembled from small, swappable Turing submodels rather
than one monolithic block (the build-up is drawn as a flowchart on the
[Analysis](analysis.md) page).
There are three layers.

**Building-block submodels**, one per parameter family, each owning
its own priors:

- `exponential_growth_model` samples the doubling time `τ` and the
  doubling-time multiplier `m = T/τ`, not `τ` and `T` directly, to
  break the `C(T) = exp(rT)` ridge.
- `delay_model` is the gamma onset-to-death delay.
- `cfr_model` is the case-fatality ratio.
- `detection_window_model` is the McCabe rectangular export detection
  window `w`; the default export mechanism instead reuses the DRC
  onset-to-report delay `f_rep` (`report_delay_model`) as the
  onset-to-detection delay.
- `surveillance_dispersion_model` samples on the `1/√k` scale.
- `pooled_ascertainment_model` partially pools the DRC and Uganda
  reporting fractions `p_drc` and `p_uganda` on the logit scale.

**Observation submodels**, one per data stream, each taking the growth
state, adding its forward integral and likelihood: `exports_model`
(Poisson), `deaths_model` (NegBinomial), `reported_cases_model`
(NegBinomial), `confirmed_cases_model` (NegBinomial), and
`exports_deaths_model` (Poisson). The Uganda export streams have
delay-convolution variants (`exports_delay_model`,
`exports_deaths_delay_model`, `exports_detection_timing_delay_model`)
that replace the rectangular detection window with an onset-to-detection
delay reusing the DRC `f_rep`; these are the `bvd_joint` defaults.

**Composers** stitch the blocks into full generative models:
`exports_only_model`, `deaths_only_model`, `cases_only_model`,
`exports_deaths_only_model`, `imperial_only_model` (exports and
deaths, the Imperial joint configuration), and `bvd_joint` (all four
streams).
Each composer conditionally includes only the likelihoods for the
streams it carries.
A single-stream composer never instantiates the other observation
submodels, so a discrete stream is never left sampled, which would
trip Turing's model check.
Pass a stream as `missing` to drop its likelihood; `bvd_joint` with
all streams missing is the generator used for the prior and posterior
predictive checks.

## Conventions

- Maximum 80 characters per line of code.
- One sentence per line in prose and markdown; do not wrap prose at 80
  characters.
- The abstract is single-sourced in `README.md`, wrapped in
  `<!-- ABSTRACT:START -->` / `<!-- ABSTRACT:END -->` markers.
  Edit the abstract in `README.md` only.
  `docs/examples/analysis.jl` loads it at build time via a Documenter
  `@eval` block that reads `README.md` and regex-extracts the text
  between those markers, so do not duplicate it into the analysis
  page.
- Table-construction and other setup code in `analysis.jl` is hidden
  inside `<details>` dropdowns via `#md # @raw html` blocks; the bare
  result object follows (with `#hide`) so only the output renders.
- The surveillance dispersion prior is a half-normal
  `truncated(Normal(0, 1); lower = 0)` on `inv_sqrt_k`.
- Docstrings use DocStringExtensions (`$(TYPEDSIGNATURES)`).
- The AD backend is Mooncake reverse-mode; integrals use Gauss-Legendre
  quadrature (`DEATH_INTEGRAL_ALG` with `n = 64`,
  `CUMULATIVE_INTEGRAL_ALG` with `n = 32`); models compose via
  `~ to_submodel(...)`.
  The deaths-among-exports CDF is written as an inner integral of the
  density because the reverse-mode AD backend does not support the
  gamma CDF shape-parameter derivative.
- NaN and Inf safe clamps (`safe_nbinomial`, `eps`-flooring of
  expected counts) guard against extreme NUTS warmup proposals; keep
  them when editing the likelihoods.

### Analysis report prose

These apply to the narrative prose in `docs/examples/analysis.jl`.
Use the existing report text as the template for tone.

- No code references in the narrative.
  Do not name functions, parameters, files, or `:symbols` in the prose.
  Describe each quantity in words, and define a derived quantity in
  words the first time it appears, near its figure or table.
- Concise and direct.
  Cut filler and adjectives.
  Avoid the LLM-indicator words: comprehensive, leverage, robust,
  framework (when vague), utilise, facilitate, novel, landscape, foster,
  harness, streamline, pivotal, nuanced, multifaceted, cornerstone,
  synergy, overarching.
- Report intervals as sentences, without a leading median.
  Write the credible interval as a phrase, not a
  "median (lower, upper)" construction.
- Minimise colons and dashes in prose; use them only when needed.
- UK English throughout.
- Section and subsection titles are just the title.
  No descriptive suffix after a title (not "Reproduction number —
  weekly random walk with intervention ramp", just "Reproduction
  number"), and no detail-dump in the first sentence after a heading.
- Order the methods generatively, infections through to observation
  endpoints: the infection process first, then the epidemiological
  processes (delays, case-fatality ratio), then the observation models
  (surveillance streams before exports), then the joint model.
- Define every quantity before it is used.
  Define the reproduction number before the seeding that relies on it;
  introduce the initial infection count before describing how it
  arises; define every symbol and operator (including convolution) the
  first time it appears. Never use a symbol the reader has not met.
- Do not repeat. State a convention once (the credible-interval levels,
  the delay discretisation) and do not restate it per bullet or
  subsection. Cut sentences that duplicate earlier content.
- Cite the source of each prior and carry the uncertainty the source
  reports. When a source gives a distribution with uncertainty (a shape
  and scale with intervals), propagate that, not a self-assigned
  weakly-informative spread. Do not write "with an assumed
  weakly-informative spread" repeatedly. If a prior is our own choice,
  say so plainly ("we use a prior of ...").
- State assumptions as assumptions ("we assume a single seed case",
  "we assume the response scale-up takes about three weeks"). Do not
  assert a false rationale for a modelling choice (not "a Poisson
  because the count is small").
- Do not editorialise or justify priors in the narrative (not "a
  diffuse prior would let the background absorb the whole stream").
  State what the model does.
- Methods belong in the methods. Do not leave model description (the
  intervention model, the counterfactual, the forecast, the evaluation)
  in the results; move it to the methods and keep the results to
  findings.
- Label quantities accurately. Do not call suspected cases onsets;
  prefer "current cumulative" over "final cumulative".
- For a latent quantity (infections, onsets, deaths) report the
  modelled estimate without overlaying observed data that sits
  downstream of unmodelled processes.
- Plots use the same credible-interval ribbons as the tables, not a
  bare median, and show only the period being estimated rather than
  greying out the rest.
- Model code shown in the report is clean; strip working comments
  before it is displayed.
- Flag a future improvement as a GitHub issue, not a buried caveat in
  the prose.

## Pull requests

- `main` is branch-protected; changes go through pull requests.
- Run the test suite before opening a pull request.
- Add a bullet to the [News](news.md) page under `Unreleased` for any
  user-visible change.
