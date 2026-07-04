# Model review — BVDOutbreakSize (July 2026)

A detailed review of the renewal joint model, its supporting code, and its
validation. The aim is to add to — not repeat — the existing issue backlog and
the queued work in `notes/review-queue.md`. Where a point is already tracked it
is cross-referenced rather than re-argued.

Reviewed at commit on branch `claude/model-review-improvements-mmkbz2`
(`src/models/{joint,priors,observations}.jl`, `src/{renewal,sampling,summaries,
forecast,counterfactual,confirmed_cfr,data,constants}.jl`, the `test/` suite and
`docs/examples/`).

## 1. Overall assessment

This is an unusually careful and mature real-time model. The generative story
is coherent end to end: a two-phase renewal (analytic cryptic exponential →
daily recursion) seeded from a molecular-clock growth rate, one shared onset
incidence routed into every surveillance stream through sampled, double-interval
-censored delays, and a genetic TMRCA bound that anchors the outbreak age. The
renewal primitives are AD-transparent, NaN/Inf-guarded, and genuinely
well-unit-tested against hand calculations and independent references. The
engineering around the model — content-addressed fit caching, freeze-and-refit
matched-in-time comparisons, per-stream calibration, package-quality gates
(Aqua/JET/ExplicitImports/doctests), and sampler-tuning scripts — is excellent.

The review therefore concentrates on where the *inference* is fragile rather
than where the code is wrong, plus a small number of concrete
coherence/consistency defects found in the code.

The single most important theme: **the abundant data (suspected and confirmed
cases) does not identify the headline quantity (outbreak size `C_T`)**. Size is
pinned by the deaths stream, the Uganda exports, the genetic bound, and the `m`
prior — all comparatively thin or prior-heavy — while the suspected/confirmed
streams are structurally degenerate with the non-BVD background and
ascertainment. Most of the findings below are facets of that fact.

## 2. Priority findings and recommendations

### P1 — The headline size rests on an undecided structural choice

`notes/review-queue.md` records an unresolved, mutually-exclusive headline
decision: the **composition-linked** lab positivity (the current default,
`confirmed_positivity_link = :composition`) makes the laboratory data identify
`λ_bg` and implies a large size (~5000), whereas the **free-positivity** link
leaves `λ_bg` weak and lets deaths/exports pin a much smaller size (~2280). A
~2× swing in the headline turns on a modelling switch that is currently a
default, not a conclusion.

Recommendation: treat this as the top-line caveat, not an internal note. Publish
**both** links side by side as a first-class sensitivity output (the machinery
already exists — `positivity_link` is a keyword on `confirmed_cases_model` and
`bvd_joint`), and state explicitly in the report which streams are doing the
identifying under each. Until the choice is defensible on external grounds, the
headline interval should arguably span both, or the report should lead with the
quantity that is robust to it. This subsumes and sharpens the identifiability
concern in issue #243 (why `R_t`/size look flat/free given linear data).

### P1 — No simulation-based calibration or parameter-recovery test

For a Bayesian model of this size (dozens of submodels, a length-6 pooled
dispersion, several weakly-identified nuisance blocks), the headline
correctness check — simulate data at known parameters, confirm the posterior
recovers them — is absent. The `:slow` end-to-end tests assert only
finiteness/positivity/shape, so a likelihood mis-specification that still
produces positive finite draws would pass silently.

Recommendation: add a simulate-and-recover test on `infection_model` and on a
reduced `bvd_joint` (a few streams) that checks `C_T`, `R0`, `T`, `CFR`, and the
key ascertainment/background parameters are recovered within their intervals at
a handful of seeds; longer term, rank-statistic SBC on the reduced joint. This
is the biggest single gap and it directly de-risks the P1 identifiability
concern above — recovery under the composition link vs the free link would show
which structural choice is actually recoverable from simulated ground truth.

### P1 — No convergence gating and no proper scoring rule

- Convergence (`fit_diagnostics`) is only smoke-tested on a trivial Gaussian;
  no R-hat/ESS/divergence threshold is asserted anywhere, and never against the
  real joint. Divergence-geometry fragility is documented in code comments
  (`src/models/observations.jl:266`, `src/models/priors.jl:1058`) but not
  guarded. With `target_accept` lowered to 0.85 for build speed
  (`src/sampling.jl`), a regression that reintroduces divergences would not fail
  CI.
- Forecast evaluation stops at bias + interval coverage + a single `within_90`
  flag on one ~1-week-back freeze. There is **no CRPS / weighted-interval /
  log score** and no rolling multi-horizon scoring. For a model whose stated
  purpose is a real-time forecast, that is thin.

Recommendation: (a) add a convergence assertion to the docs/CI path on the real
fit (max R-hat < ~1.01, min bulk ESS above a floor, divergences at/near zero),
so the published build fails loudly if geometry degrades; (b) add a proper
scoring rule (CRPS or WIS via a small helper, or `scoringutils`) and score the
one-week-ahead forecast over a rolling set of frozen cut-offs, not one.

### P2 — Overlapping streams scored as conditionally independent (tracked, quantify it)

Already issue #307. Suspected cases, analysed volume, confirmed cases, isolation
occupancy, recovered, and the treatment flows are largely deterministic
functions of the *same* latent `bvd_reports_daily`/`bg_daily`, yet each
contributes an independent NB/Binomial term. As the v1.6+ streams multiplied,
the number of terms sharing one case pool grew, so the `C_T` interval is likely
too narrow. Recommendation: run the leave-a-stream-out sensitivity #307
proposes and report the interval-width delta; this pairs naturally with the SBC
work (coverage under simulation would expose over-precision directly).

### P2 — Heavy dependence on fixed geometry constants and a fixed clock

`SEEDING_LEAD_DAYS=30`, `RENEWAL_START_LEAD=14`, `RT_WALK_LEAD=28`, the single
molecular clock rate, and `tmrca_days_sd=15` are fixed and materially shape the
outbreak age `T` and hence size. Clock-rate sensitivity exists
(`sens_fast_clock`) but is gated behind `RUN_SENSITIVITY`, which is **off on
main/PR builds** and only on for release-tag builds — so the published stable
docs may show no geometry/clock sensitivity at all. Recommendation: run at least
a lightweight `RENEWAL_START_LEAD` and clock-rate sensitivity on every docs
build (or clearly state in the report that the stable build omits them), and
consider propagating clock-rate uncertainty into `tmrca_days_sd` rather than
fixing it.

### P2 — `C_T = 2^m` size prior is heavy-tailed and unbounded by any test

`m ~ truncated(Normal(3, 3); lower=0)` gives `2^m` a very heavy right tail
(`m` at +2SD ≈ 9 → 512; the tail runs much higher before truncation bites).
The genetic bound and deaths/exports are relied on to tame it, but no
prior-predictive test bounds the *implied* `C_T` scale, and prior-predictive
coverage in the suite is uneven (only `λ_bg` gets a real prior-predictive
sanity bound). Recommendation: add a prior-predictive check asserting the `C_T`
prior stays within a plausible order of magnitude, mirroring the existing
`λ_bg` background-vs-observed bound in `test/test_test_positivity.jl`.

## 3. Concrete code-level defects (low-risk fixes)

These are small, verified, and safe to correct independently of the modelling
discussion.

1. **`m_prior_centre` is vestigial and its docstring is stale/misleading.**
   `src/data.jl:252` (and the cross-reference at `src/constants.jl:115`) still
   describe `m` via the *old* total-doublings interpretation — "18 May ⇒
   `m ≈ log2(501) ≈ 9`", "`C_T = 2^m` is the cumulative infection count". Under
   the current model `m` counts **only the cryptic-phase doublings**
   (`M_PRIOR_BASE = 3.0`), and the joint uses the fixed `Normal(3, 3)` prior —
   `m_prior_centre` is exported and documented but **never wired into any
   prior**. Worse, its advancing-centre logic (`m_base + elapsed/doubling_days`)
   would be *wrong* to use now: the observation-window growth is already added
   separately as `τ_obs` in `infection_model`, so advancing the cryptic centre
   with calendar time would double-count it. Recommendation: either delete
   `m_prior_centre` (and drop the exports/cross-refs), or rewrite its docstring
   to the cryptic-phase interpretation and add a warning that it must not be
   combined with the `τ_obs` term. At minimum fix the `m ≈ 9` / "cumulative
   infection count" wording.

2. **`forecast_vs_truth_trajectory` uses the projection method the headline
   forecast deliberately abandoned.** `src/forecast.jl:531-533` projects
   `cases_T * grow` with `grow = exp(r·horizon)` and `cases_T =
   expected_reports_T` (a *cumulative* expected total) — i.e. it scales the
   cumulative stock by the growth factor. `forecast_reported` explicitly rejects
   exactly this (`src/forecast.jl:315-319`: scaling the stock "would shrink it
   whenever the growth rate is negative") and instead adds projected *new*
   incidence to the observed cut-off cumulative. So this exported function will
   shrink the cumulative below the cut-off when `r < 0` and is inconsistent with
   the report's own forecast. It is exported but not called anywhere in
   `docs/examples/`. Recommendation: rewrite it to the additive `_new_h` method
   (or remove it), so the public API can't be used to produce a projection that
   contradicts the headline one.

3. **`reconstruct_rt` / `_evolving_rates` / `_gi_pmf` re-encode the model's
   walk and generation-interval discretisation.** The `R_t` knot/interpolation/
   ramp logic and the GI PMF are duplicated across `src/plots.jl:1107`
   (`reconstruct_rt`) and `src/forecast.jl` (`_evolving_rates`, `_gi_pmf`),
   separate from `rt_walk_model`/`generation_interval_model`. Any change to the
   model's walk (knot spacing, ramp, clamp behaviour) must be mirrored in three
   places or the plots/forecast silently drift from the fit. Recommendation:
   factor a single shared reconstruction helper the model, plots, and forecast
   all call, and add a test that its output matches the model's own `Rt` on a
   fixed draw.

4. **Doc drift in sampler defaults.** `src/sampling.jl` sets `target_accept =
   0.85`, but `scripts/trial_target_accept.jl` still calls 0.90 "the new
   default" in its header. Harmless, but reconcile so the recorded default is
   unambiguous (and note the divergence-gating point in P1 — 0.85 is on the
   aggressive side for a model with documented geometry fragility).

## 4. Smaller observations

- **Euler–Lotka refinement depth.** `euler_lotka_r` uses 2 Newton steps from the
  small-`r` seed (`src/renewal.jl:117`). It feeds only the *reported* current
  growth `r`/doubling time (not the likelihood), so accuracy is non-critical,
  but for `R_t` far from 1 (warmup, or a strong intervention) two steps can be
  short of convergence. If ever used somewhere load-bearing, raise `steps` or
  iterate to tolerance.
- **`check_model = true` default with a documented exception.** `nuts_sample`
  defaults to `check_model = true`; `exports_deaths_only_model` needs it off
  (sampled discrete `Poisson` in predictive mode). This is a caller
  responsibility that is easy to trip on a new single-stream composer — worth a
  one-line note at each composer that runs a stream in predictive mode.
- **Confirmed-CFR rescaling.** `delay_corrected_confirmed_cfr`
  (`src/confirmed_cfr.jl:110-116`) rescales the modelled daily confirmed
  incidence so its total matches the scored `expected_confirmed_T`, reconciling
  the modelled-analysed-volume series with the observed-denominator scoring.
  Defensible, but it is a reconciliation step that should be called out in the
  report's methods, since the corrected CFR denominator depends on it.
- **Background/case-finding scaling** (issue #374) and the **onset-to-sample
  grounding** (issues #359/#377) are the most valuable *modelling* follow-ups
  already scoped; they would replace a prior-dominated nuisance (`lab_delay`,
  issue #128) and the additive-background assumption with data-anchored
  structure, which is the right direction given the P1 identifiability theme.

## 5. Data and provenance

`load_observations`/`freeze_observations` (`src/data.jl`) are clean: dated TOML
in, grid indices derived at load, freezing recomputes scalars from truncated
histories so matched-in-time fits are genuine. Two provenance points the report
already discloses but that a reviewer should keep prominent: the primary
situation-report figures are **LLM-extracted from PDFs** (with a second re-read
pass), and the model code/priors/analysis were **LLM-drafted and not
independently replicated** against the source authors' code. Both are honestly
flagged in `analysis.jl`'s Limitations; given the headline's dependence on a few
thin streams (deaths, exports), an independent human spot-check of the deaths
and exports series against the source sitreps would be worth more than any
additional modelling.

## 6. Suggested sequencing

1. Land the low-risk code fixes in §3 (independent, no refit needed).
2. Add the SBC/recovery test and the convergence + prior-predictive `C_T`
   assertions (§2 P1/P2) — these gate everything else.
3. Elevate the composition-vs-free lab link to a published sensitivity and add a
   proper scoring rule over rolling freezes (§2 P1).
4. Then the tracked modelling follow-ups (#374, #359/#377, #307 quantification).

Nothing here calls the model's construction into question — it is a thoughtful,
well-built real-time analysis. The recommendations are about making the
headline's dependence on its thin, prior-heavy identifying streams *visible and
tested*, rather than implicit.
