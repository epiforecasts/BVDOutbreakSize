# Contributing

Issues and pull requests are welcome at [epiforecasts/BVDOutbreakSize](https://github.com/epiforecasts/BVDOutbreakSize).
This page covers how the project is laid out, how to run it, and the conventions to follow when changing it.

## Repository layout

- `src/BVDOutbreakSize.jl` — the module entry point: dependencies, the export list, and the `include` order for the rest of `src/`.
  The functionality itself is split across single-purpose files, listed next.
- `src/data.jl` — `load_observations` and `freeze_observations`, the observation-manifest loader.
- `src/constants.jl` — fixed constants, including the published Imperial point estimates (`REPORT_SCENARIOS`, `REPORT_SCENARIOS_CI`).
- `src/sampling.jl` — NUTS sampling (`nuts_sample`, `fit_parallel`) and the AD backend setup.
- `src/renewal.jl` — the shared renewal-process helpers: `renewal_infections`, the delay convolutions `convolve_delay`, `convolve_survival` and `convolve_pmf`, and `discretise_censored`.
- `src/models/priors.jl`, `src/models/observations.jl`, `src/models/joint.jl` — the building-block submodels, the observation submodels and the composers.
  See [Model architecture](#Model-architecture) below.
- `src/summaries.jl`, `src/scoring.jl` — summary and comparison tables, and forecast scoring.
- `src/counterfactual.jl` — the no-onward-deaths projection (`predict_no_onward_deaths`).
- `src/forecast.jl` — forecast helpers (`forecast_reported`).
- `src/confirmed_cfr.jl` — delay-corrected confirmed-case-fatality-ratio helpers.
- `src/plots.jl` — plotting.
- `docs/examples/analysis.jl` and `docs/examples/sensitivity.jl` — the Literate walkthroughs, split so the expensive fits and the render can fan out across CI runners.
  `analysis.jl` carries the methods, results and one-week-ahead forecast.
  `sensitivity.jl` carries the forecast validation and the comparison/sensitivity analyses.
  Both load their fits through the shared `docs/examples/_setup.jl`.
  `analysis.jl` is the main artifact.
- `docs/fits/` — the fit-cache machinery: `registry.jl` (the fit-id list), `cache.jl` (content-addressed fit caching under `logs/fit_cache`), `one.jl` (fit and cache a single id, `task fit`), `all.jl` (fit every model, `task fit-all`), and `list.jl` (print fit ids for the CI matrix).
- `docs/execute.jl` — runs one Literate page against the fit cache and writes its markdown, figures and half of `output/` (used by `task docs-main` and `task docs-sensitivity`).
- `docs/make.jl` — the Vitepress combine step: copies `README.md` to `index.md`, assembles the site from the already-rendered markdown, and builds the bibliography (used by `task docs`).
- `data/observations.toml` — single source of truth for observation data (case and death counts, traveller volumes, sources).
  Loaded via `load_observations()` and never hardcoded.
  Update this one file for a new situation report and the analysis picks it up.
  The literate re-binds its observation `const`s from the loaded TOML, so the package constants are defaults only.
- `scripts/run.jl` — regenerates published results by including the literate and writes CSVs to `output/`.
- `test/` — one file per feature, driven by `test/runtests.jl`, plus `enzyme/`, `formatter/`, `jet/` and `package/` for the quality checks.
- `external/bdbv-linelist-analysis` — git submodule, source of the onset-to-death delay priors.

## Running and testing

`Taskfile.yml` at the repository root wraps the common commands.
Run `task --list` for the full set with descriptions.
The ones used day to day:

```bash
# Instantiate the package environment (no task wraps this)
julia --project=. -e 'using Pkg; Pkg.instantiate()'

# Fit and cache every model, then render the docs from the cache
task fit-all
task docs

# Fit and cache a single model by id (list ids with
# `julia --project=docs docs/fits/list.jl`)
BVD_FIT_ID=deaths task fit

# Render just the main analysis page from the cache, for fast iteration
task docs-main
task docs-sensitivity

# Format (SciML style) over src/, test/, docs/, scripts/
task format

# Run the full test suite, or skip the quality checks (Aqua/JET/format)
# for a quicker one
task test
task test-quick

# Regenerate the published output CSVs into output/ (no task wraps this)
julia --project=. scripts/run.jl
```

Running `julia --project=. docs/examples/analysis.jl` directly instead steps through the full narrative without going through the fit cache, fitting every model inline.
This is the slow path.
Prefer `task fit-all && task docs` for anything beyond reading the source.

A build streams per-fit progress by default: every NUTS fit writes `logs/<fit>.log` (iteration, log-density, divergences) and a TensorBoard run under `logs/tensorboard/<fit>/`, controlled by `BVD_FIT_LOG` (`all` when unset, or `progress`, `tensorboard`, `none`).
CI release builds set `BVD_FIT_LOG=none`.
Tail a log for quick liveness, or run `task tensorboard` to view all fits in the worktree.
The logs live under the git-ignored `logs/`, so each worktree keeps its own.
A cached fit lives under `logs/fit_cache` and is keyed by a content hash of the fit-relevant code and data, so a Turing or dependency version bump invalidates it.
Refit rather than debugging a `KeyError` on a stale chain.

`test/runtests.jl` includes each `test/test_*.jl`.
To iterate on one file, run it inside a REPL after `using BVDOutbreakSize`, or temporarily comment out the others in `runtests.jl`.

CI runs the test suite (`.github/workflows/test.yml`) and builds the docs, publishing `output/` as a GitHub Release on each push to `main` (`.github/workflows/docs.yml`).

## Model architecture

The model is assembled from small, swappable Turing submodels rather than one monolithic block (the build-up is drawn as a flowchart on the [Analysis](analysis.md) page).
There are three layers.

**Building-block submodels**, one per parameter family, each owning its own priors:

- `exponential_growth_model` samples the doubling time `τ` and the doubling-time multiplier `m = T/τ`, not `τ` and `T` directly, to break the `C(T) = exp(rT)` ridge.
- `onset_to_death_model` is the gamma onset-to-death delay, built by convolving two atomic `gamma_delay_model` delays (onset→admission, admission→death) rather than fitting the onset→death delay directly.
- `cfr_model` is the case-fatality ratio.
- `surveillance_dispersion_model` samples on the `1/√k` scale.
- `pooled_ascertainment_model` partially pools the DRC and Uganda reporting fractions `p_drc` and `p_uganda` on the logit scale.

**Observation submodels**, one per data stream, each taking the growth state, adding its forward integral and likelihood: `exports_model` and `exports_deaths_model` (the Uganda imports and import deaths, each an inhomogeneous Poisson process over the dated detection days rather than a rectangular detection window), `deaths_model`, `reported_cases_model`, `confirmed_cases_model` and `confirmed_deaths_model` (NegBinomial DRC surveillance and laboratory streams), and `treatment_flow_model` and `recovered_model` (the isolation-occupancy and recovered-among-confirmed streams).

**Composers** stitch the blocks into full generative models: `exports_only_model`, `deaths_only_model`, `cases_only_model`, `confirmed_only_model`, `confirmed_deaths_only_model`, `treatment_only_model`, `exports_deaths_only_model`, `exports_joint_only_model` (the Uganda export cases and deaths fitted jointly over one shared at-risk prevalence), and `bvd_joint` (every stream together).
Each composer conditionally includes only the likelihoods for the streams it carries.
A single-stream composer never instantiates the other observation submodels, so a discrete stream is never left sampled, which would trip Turing's model check.
Pass a stream as `missing` to drop its likelihood.
`bvd_joint` with all streams missing is the generator used for the prior and posterior predictive checks.

## Conventions

- Maximum 80 characters per line of code.
- One sentence per line in write-up prose and markdown.
  Do not wrap prose at 80 characters.
  Code, code comments and docstrings keep the 80-character limit and wrap normally.
  This rule is for prose only.
- The shared front matter (title, authors, abstract, scope) is single-sourced in `README.md`, up to the `<!-- SHARED:END -->` marker.
  Edit it in `README.md` only.
  `docs/examples/analysis.jl` loads it at build time via a Documenter `@eval` block that reads `README.md` and extracts everything before that marker, so do not duplicate it into the analysis page.
- Table-construction and other setup code in `analysis.jl` is hidden inside `<details>` dropdowns via `#md # @raw html` blocks.
  The bare result object follows (with `#hide`) so only the output renders.
- The surveillance dispersion prior is a half-normal `truncated(Normal(0.6, 0.2); lower = 0)` on `inv_sqrt_k`, following the Stan prior-choice recommendations.
- Docstrings use DocStringExtensions (`$(TYPEDSIGNATURES)`).
- Two implementation details: the AD backend is Mooncake reverse-mode, and models compose via `~ to_submodel(...)`.
- NaN and Inf safe clamps (`safe_nbinomial`, `eps`-flooring of expected counts) guard against extreme NUTS warmup proposals.
  Keep them when editing the likelihoods.

### Analysis report prose

These apply to the narrative prose in `docs/examples/analysis.jl`, and to write-up prose generally.
Use the existing report text as the template for tone.
The measured sentence- and paragraph-level rules below were reverse-engineered from a manuscript the maintainers are happy with.
The repo-specific rules that follow take precedence where the two disagree.
`AGENTS.md` carries a short version of this section for agent sessions, pointing back here for the full rule.
Keep the two in agreement: if this section changes in a way that affects the summary, update `AGENTS.md` too.

- No code references in the narrative.
  Do not name functions, parameters, files, or `:symbols` in the prose.
  Describe each quantity in words, and define a derived quantity in words the first time it appears, near its figure or table.
- Concise and direct.
  Cut filler and adjectives.
  Avoid the LLM-indicator words: comprehensive, leverage, robust, framework (when vague), utilise, facilitate, novel, landscape, foster, harness, streamline, pivotal, nuanced, multifaceted, cornerstone, synergy, overarching.
  Avoid filler and stance markers: very, significant, crucial, essential, clearly, obviously, "it is important to note that".
  Delete the LLM connective openers `Furthermore`, `Moreover`, `Additionally`, `In conclusion`, `Overall`.
  Restrict genuine connectives to `However,` `Whilst`, `Yet`, `Unlike`, `Similarly`, `Alternatively`, `Finally,`, `Instead`.
- Target 15-30 words per sentence and treat 40 as the ceiling.
  Exceed it only for a flat, colon-led enumeration, never for nested subordination.
  Use at most two commas per sentence for internal structure, a third only in a serial list.
- Never use a dash as punctuation.
  Replace it with a full stop, or with parentheses when the aside is a gloss or example.
  Use a colon only to introduce a list or an expansion, never to join two independent statements.
  Use a semicolon only to separate list items that themselves contain commas.
- Instead of a trailing qualifier (", which means that...", ", thereby enabling..."), cut the comma and start a new sentence with `This`, `These` or `It also`.
  Keep any surviving trailing `, which` clause to about one sentence in ten, and only where the clause carries the point rather than padding it.
- Write paragraphs of 2-6 sentences.
  Open each with its claim, not with context, and close it with a consequence or a limitation, never a restatement of the opening.
- Hedge with `can`, `could`, `may` or `potential`, one hedge per claim.
  Never stack hedges ("may potentially", "could possibly suggest").
- Report intervals as sentences, without a leading median.
  Write the credible interval as a phrase, not a "median (lower, upper)" construction.
- Report numbers with their provenance and units.
  Delete a vague quantity ("a large number of", "significantly more") rather than leave it unsourced.
- UK English throughout.
- Section and subsection titles are just the title.
  No descriptive suffix after a title (not "Reproduction number — weekly random walk with intervention ramp", just "Reproduction number"), and no detail-dump in the first sentence after a heading.
- Order the methods generatively, infections through to observation endpoints: the infection process first, then the epidemiological processes (delays, case-fatality ratio), then the observation models (surveillance streams before exports), then the joint model.
- Define every quantity before it is used.
  Define the reproduction number before the seeding that relies on it; introduce the initial infection count before describing how it arises; define every symbol and operator (including convolution) the first time it appears.
  Never use a symbol the reader has not met.
  Give an equation once in display maths and refer to its symbols inline afterwards; never narrate a displayed equation in words.
- Do not repeat.
  State a convention once (the credible-interval levels, the delay discretisation) and do not restate it per bullet or subsection.
  Cut sentences that duplicate earlier content.
- Cite the source of each prior and carry the uncertainty the source reports.
  When a source gives a distribution with uncertainty (a shape and scale with intervals), propagate that, not a self-assigned weakly-informative spread.
  Do not write "with an assumed weakly-informative spread" repeatedly.
  If a prior is our own choice, say so plainly ("we use a prior of ...").
- State assumptions as assumptions ("we assume a single seed case", "we assume the response scale-up takes about three weeks").
  Do not assert a false rationale for a modelling choice (not "a Poisson because the count is small").
- Do not editorialise or justify priors in the narrative (not "a diffuse prior would let the background absorb the whole stream").
  State what the model does.
- Methods belong in the methods.
  Do not leave model description (the intervention model, the counterfactual, the forecast, the evaluation) in the results.
  Move it to the methods and keep the results to findings.
- No project or issue history in the narrative: no PR or issue numbers, no "we previously used X", "this was changed to Y", "after refactoring".
  State the current design only.
  A difference from work being replicated, or a limitation of an external dependency, is the exception, written as a present-tense fact.
  An issue number may also stay where the issue is the record for something still unresolved, such as a workaround pointing at the open root cause it does not fix, or a data decision whose reasoning lives in the issue.
  The test is whether the number still tells the reader something they cannot get from the code.
  "Tracked by #495" earns its place; "changed in #495" does not.
- Label quantities accurately.
  Do not call suspected cases onsets, and prefer "current cumulative" over "final cumulative".
- For a latent quantity (infections, onsets, deaths) report the modelled estimate without overlaying observed data that sits downstream of unmodelled processes.
- Plots use the same credible-interval ribbons as the tables, not a bare median, and show only the period being estimated rather than greying out the rest.
- Model code shown in the report is clean.
  Strip working comments before it is displayed.
- Flag a future improvement as a GitHub issue, not a buried caveat in the prose.
- Bullet lists, bold and italics are advisory rather than enforced here, since this is a technical report rather than a journal manuscript.
  Do not add new bullet lists or bold to the narrative prose, and do not restructure an existing list that carries genuinely parallel content, but do strip bold or italics used for mid-sentence emphasis.

## Pull requests

- `main` is branch-protected.
  Changes go through pull requests.
- Run the test suite before opening a pull request.
- Add a bullet to the [News](news.md) page under `Unreleased` for any user-visible change.
