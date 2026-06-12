# Suspected- and confirmed-death model redesign

This note records the redesign of the suspected-death (`deaths_model`) and
confirmed-death (`confirmed_deaths_model`) submodels, the reasoning behind each
change, and the posterior checks used to validate it.
It is the working record behind the pull request.

## Motivation

The previous death pathway had several gaps relative to the case pathway.

1. No death ascertainment.
   Suspected deaths were `CFR ·` (onset-to-death-convolved infections) with no
   capture fraction, i.e. every BVD death was assumed to enter the suspected-
   death count, while suspected cases carried `p_drc ≈ 0.75`.
2. No non-BVD background deaths.
   The suspected-death background (`λ_bg_death`) was switched off by default,
   so the suspected deaths were treated as pure BVD even though the suspect-
   death definition (symptomatic-then-deceased) admits non-BVD deaths.
3. Confirmed deaths grounded on the CASE composition.
   The confirmed-death confirmation probability was the suspected-CASE BVD
   composition `q_susp` enriched on the odds scale by a free scalar `m_death`.
   This drew the death-confirmation rate from the case pool and added a free
   enrichment that was hard to interpret.
4. The confirmed-death model was not parallel to the confirmed-case model.
   The confirmed cases fit a modelled analysed-specimen volume and score the
   positives as that volume times a composition-linked positivity; the deaths
   thinned the suspected-death total with no volume analogue.
5. The case analysed-volume modelling was under-documented.

## What changed

### Suspected deaths (`deaths_model`)

- Death ascertainment `p_death` (new `death_ascertainment_model`), an
  informative logit-Normal centred on 0.9 (a death is more reliably reported
  than a living suspect). BVD suspected deaths are now
  `p_death · CFR ·` (onset-to-death convolution).
- Non-BVD background deaths tied to the case background. The per-day
  background suspected deaths are `cfr_bg · λ_bg,t`, a background CFR
  (`background_cfr_model`, `Beta(2, 6)`) applied to the already-identified
  per-day non-BVD suspected-CASE background. This is the user's suggestion of
  "a scaled version of the case one, i.e. assuming some CFR for background
  cases." It removes the degeneracy that kept a free `λ_bg_death` switched off:
  the case background is pinned by the laboratory positivity link, so scaling
  it gives the death background a level and time profile without a second free,
  outbreak-size-degenerate rate.
- The legacy free `death_background` and per-vintage `background_re` death
  paths are kept for sensitivity analyses but are no longer the joint default.
- The joint now builds the case stream BEFORE the death stream so the case
  background is available to scale.

### Confirmed deaths (`confirmed_deaths_model`)

Rebuilt to mirror the confirmed-case laboratory pipeline:

- Death "analysed" volume = `τ_death ·` (suspected deaths carried to laboratory
  receipt by the same `receipt_pmf` the confirmed cases use). `τ_death`
  (`death_testing_fraction_model`, `Beta(2, 4)`) is the fraction of suspected
  deaths that reach the laboratory, the death analogue of the case testing
  fraction `τ_test`, centred lower because post-mortem swabbing is rarer.
- Death-pool composition `q_death = bvd_death / (bvd_death + bg_death)` from the
  death series' OWN BVD and background components, NOT the case composition.
  This is well-defined precisely because the death background is now on.
- Assay positivity `p = s·q_death + (1−spec)(1−q_death)` with the same PCR
  sensitivity and specificity priors as the confirmed cases (the death stream
  draws its own values, so it is self-contained and does not pull the case
  composition, `p_drc` or `λ_bg` directly).
- Confirmed-death increments scored as `NegBinomial(k)` counts of the modelled
  confirmed-death volume (the modelled-volume route the early and post-lab
  confirmed CASE windows also use), since no published death-analysed
  denominator exists.
- The `m_death` odds enrichment and the dependence on `q_susp` are removed.

### Deterministics

Removed `m_death`. Repurposed `death_composition` to `q_death` (death-pool
composition) and `death_confirmation` to the death-confirmation positivity.
Added `death_ascertainment` (`p_death`), `background_cfr` (`cfr_bg`) and
`tau_death`. `lambda_bg_death`, `background_death_total`,
`expected_confirmed_deaths_T` and `onset_to_death_confirmation_pmf` are kept,
so the forecast, confirmed-CFR and summary consumers are unchanged.

## Mapping to the requested changes

| Request | Status |
|---|---|
| Prior on death ascertainment (~0.9) | done — `death_ascertainment_model` |
| Report-to-receipt delay suspected→confirmed deaths | done (already present; reused in the death analysed volume) |
| Non-BVD background deaths, scaled from background cases | done — `cfr_bg · λ_bg` |
| Confirmed-death testing model not drawn from confirmed cases | done — death confirmation is self-contained (own composition, own assay draws) |
| Remove `m_death`; get `q` from death params | done — `q_death` from the death pool |
| Document the case analysed-volume modelling | done — methods + docstrings |
| Death volume as a scaling, confirmed deaths more like cases | done — `τ_death ·` suspected-death volume + composition positivity |

The one partial item is "the testing model draws a lot from confirmed cases."
The confirmed-DEATH model is now fully decoupled from the case pool. The
confirmed-CASE positivity still uses the suspect-pool composition by design —
that link is how the laboratory data identify the non-BVD background `λ_bg`,
and a free alternative already exists (`confirmed_positivity_link = :free`) for
sensitivity. Fully replacing the case composition link is a larger, separate
change and is out of scope here.

## Posterior checks

Exploratory streaming joint fits (`background_re = true`,
`positivity_link = :composition`, genetic bound on, `n = 108`). The short
300x2 fits are under-warmed (NUTS adapts only ~150 steps at `N = 300`), so the
convergence numbers are not publishable; they are reported here only to show
the redesign matches the baseline. A longer 1000x4 fit is used for the headline.

### 300x2 exploratory fit: baseline (old model) vs redesign

| Quantity | Baseline (old) | Redesign |
|---|---|---|
| max R-hat / min ESS | 2.13 / 1.5 | 2.13 / 1.6 |
| divergences (of 600) | 0 | 0 |
| C_T median [90%] | 6888 [3191, 24725] | 5300 [2953, 11423] |
| CFR | 0.37 | 0.39 |
| p_drc | 0.67 | 0.60 |
| expected_deaths_T | 828 | 1006 |
| expected_confirmed_T | 900 | 880 |
| expected_confirmed_deaths_T | — | 192 |

New death parameters (redesign, 300x2):

- `death_ascertainment` 0.90 [0.81, 0.95] — sits at its informative prior, as
  expected for a weakly-identified ascertainment that leans on the prior.
- `background_cfr` 0.087 [0.015, 0.32] — pulled DOWN from the prior mean 0.25:
  the data want few non-BVD background deaths.
- `death_composition` (q_death) 0.78 [0.47, 0.95] — ~78% of suspected deaths
  are BVD at the cut-off.
- `death_confirmation` 0.67 [0.41, 0.88] — consistent with
  `s·q_death + (1−spec)(1−q_death)` at q ≈ 0.78, s ≈ 0.85.
- `tau_death` 0.44 [0.28, 0.66] — a moderate share of suspected deaths reach
  the laboratory.

Reading:

- The new death parameters are all well-behaved and interpretable; the death
  model now mirrors the case model.
- Convergence is IDENTICAL to the baseline (R-hat ≈ 2.1 at the short fit), so
  the redesign neither improves nor degrades sampler geometry; the multimodal
  short-warmup behaviour is pre-existing.
- The redesign slightly TIGHTENS the C_T upper tail (24725 → 11423).
- The model over-predicts the fitted cut-off counts on EVERY stream (confirmed
  cases 880 vs observed 676; confirmed deaths 192 vs 136; the suspected-case
  expectation is similarly high). This is systemic, NOT a death-model artefact:
  it is the first-vintage / cryptic-phase background over-accrual (the first
  per-vintage bin sums the modelled series from grid day 1, so the non-BVD
  background accrues over the ~85-day pre-surveillance span). A separate piece
  of work fixes that root cause; because the death background is now scaled
  from the case background, it will benefit automatically once the case-side
  fix lands.

### 1000x4 headline fit

(appended after the longer fit completes.)
