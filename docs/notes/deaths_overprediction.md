# Confirmed-deaths over-prediction and the death-background prior

Investigation on `main`'s deaths model (branch `deaths-investigate`, off
commit 738bcaf, data cut-off 4 June / SitRep 021).
Fits use the headline joint built by
`scripts/deaths_overprediction.jl` at relaxed settings
(600 samples, 2 chains, `target_accept = 0.9`).

## Summary

The confirmed-deaths stream over-predicts.
Conditional posterior-predictive confirmed deaths sit near 146 against the
82 observed at the 4 June cut-off (the figure was 64 at the earlier 3 June
cut-off the task was written against; the data have since advanced).
The cause is that the non-BVD suspected-death background `λ_bg_death` is
pinned low by its tight prior, which forces the BVD share of the
suspect-death pool `q_death` to ~0.97 and the death positivity to ~0.77,
far above the ~0.33 confirmation the data imply.

Adjusting the death-background prior does not fix the level cleanly.
A wider half-normal (SD 1.0) leaves confirmed deaths essentially unchanged
(PP 147 vs 146) and degrades mixing (divergences 2 → 30, bulk-ESS 432 →
122).
A prior re-centred on a large background (Normal(1.24, 0.6)) does pull the
level down (PP 99) but only by splitting the posterior into the degenerate
bimodal regime (R-hat 2.125, bulk-ESS 3), the death analogue of the
case-side multimodality that broke the joint at R-hat ~2.1 when the case
background was forced up (issue #151).
The level cannot be moved by the background prior alone because the death
positivity has no severity-enrichment / random-effect freedom: it is tied
rigidly to the pool composition `q_death`, so the model resolves the
data-vs-model conflict elsewhere (in `p_deaths`, pinned at its ceiling)
rather than by lowering positivity.

## Diagnosis (baseline, current model)

Death-background prior `truncated(Normal(0, 0.25); lower = 0)` (the
`main` default).

| Quantity | Median [5%, 95%] |
|---|---|
| `CFR` | 0.295 [0.180, 0.450] (≈ prior mean) |
| `p_deaths` (deaths ascertainment) | 0.994 [0.915, 1.077] — pinned at ceiling |
| `λ_bg_death` | 0.146 [0.011, 0.440] — pinned low |
| `q_death` at cut-off | 0.965 [0.896, 0.998] |
| death positivity at cut-off | 0.765 [0.679, 0.851] |
| `τ_death` (forwarding) | 0.274 [0.154, 0.459] |
| expected confirmed deaths (det.) | 176 [110, 290] |
| PP confirmed deaths (cond.) | 146 [49, 501] (obs = 82) |

Mixing: max R-hat 1.010, min bulk-ESS 432, 2 divergences.

The observed confirmation rate is 82/246 ≈ 0.33, so the suspect-death pool
is far from all-BVD.
The tight `Normal(0, 0.25)` background prior (median ≈ 0.17/day, ≈ 9% of
the suspect-death total) forces `q_death` high and over-confirms.

## Wider death-background prior

Prior widened to `truncated(Normal(0, 1.0); lower = 0)`, matching the case
`λ_bg` default width (SD 1.0 vs the death side's 0.25).

| Quantity | Baseline | Wider Normal(0,1.0) |
|---|---|---|
| `λ_bg_death` | 0.146 [0.011, 0.440] | 0.435 [0.037, 1.568] |
| `q_death` at cut-off | 0.965 | 0.905 [0.726, 0.992] |
| death positivity at cut-off | 0.765 | 0.717 [0.572, 0.827] |
| `p_deaths` | 0.994 | 0.999 — still at ceiling |
| expected confirmed deaths (det.) | 176 | 181 |
| PP confirmed deaths (cond.) | 146 | 147 (obs = 82) |
| max R-hat | 1.010 | 1.013 |
| min bulk-ESS | 432 | 122 |
| divergences | 2 | 30 |

The wider prior lets `λ_bg_death` rise (median 0.146 → 0.435) but it stays
weakly identified (95% upper bound 1.57).
`q_death` and positivity drop only modestly and confirmed deaths do not
move toward 82 (PP 147 vs 146).
Mixing degrades sharply: divergences 2 → 30, bulk-ESS 432 → 122, and the
fit ran roughly twice as long (a tiny-step NUTS crawl), the signature of
the warned multimodality.
R-hat did not frankly split here, but the divergence and ESS collapse are
the same pathology in milder form.

## Centred-high death-background prior

Prior re-centred on a large background, `truncated(Normal(1.24, 0.6);
lower = 0)`: with ~246 suspect deaths over T ≈ 132 days, a non-BVD
background carrying ~2/3 of the pool is ≈ 1.24/day, so this prior pushes
toward the ~0.33 BVD share the confirmation rate implies.

| Quantity | Baseline | Centred Normal(1.24, 0.6) |
|---|---|---|
| `λ_bg_death` | 0.146 | 0.442 [0.349, 1.756] |
| `q_death` at cut-off | 0.965 | 0.905 [0.700, 0.926] |
| death positivity at cut-off | 0.765 | 0.678 [0.557, 0.752] |
| `τ_death` | 0.274 | 0.091 [0.089, 0.430] |
| expected confirmed deaths (det.) | 176 | 85 [84, 273] |
| PP confirmed deaths (cond.) | 146 | 99 [38, 405] (obs = 82) |
| max R-hat | 1.010 | 2.125 |
| min bulk-ESS | 432 | 3 |
| divergences | 2 | 19 |

The centred prior does pull the level toward the observed 82 (PP 99, and
one mode of `expected_confirmed_deaths_total` sits at 85).
But it does so only by splitting the posterior: max R-hat 2.125 with
bulk-ESS 3, the same frank multimodality the case side hit when its
background was forced up (issue #151).
The 85/99 figures are not a converged posterior, they are one mode of an
unmixed chain; `τ_death` collapses to ~0.09 in that mode.
This fit also ran the slowest of the three (a long tiny-step crawl).
So the centred prior is not a usable fix: the only way it lowers the level
is by entering the degenerate bimodal regime.


## Why the prior alone cannot fix the level

The over-prediction is not a free parameter the background can absorb.
The death positivity is

```
p_pos_death = s · q_death + (1 − spec) · (1 − q_death),
q_death = μ_BVD_death / N_death_susp,
```

with `s`, `spec` shared (fixed) from the case lab pipeline.
`q_death` is the raw pool composition with no upward or downward freedom.
On the case side the confirmed positivity is the composition `φ` UPSAMPLED
by a severity-enrichment `δ₀·e^{−c/decay}` plus a per-vintage random
effect (`confirmed_q_re_model`), so the tested BVD share can sit above the
pool composition while the data still identify a large non-BVD background
`λ_bg`.
The death stream has neither, so to match an observed positivity the model
must move `q_death`, which means moving the pool composition, which fights
the suspect-death count likelihood and `λ_bg_death`.
The model instead pins `p_deaths` at its ceiling and keeps `q_death` high,
over-confirming.

## Recommendation

Do not adopt a changed death-background prior as the fix.
A wider prior leaves the level unchanged and degrades mixing; a
centred-high prior only lowers the level by entering the degenerate
bimodal regime (R-hat 2.1, ESS 3).
Neither converges to a usable answer.

The level is better addressed by giving the death positivity the same
identifying freedom the case stream has, so the tested death share can
differ from the raw pool composition.
Two options, in order of preference:

1. Share the case stream's composition-link / severity-enrichment
   structure for the death positivity, deriving `q_death` from the pool
   composition upsampled by a severity-enrichment (and optionally a small
   per-vintage random effect), with `s`/`spec` still shared.
   This ties death positivity to `λ_bg_death` the way the case side ties
   positivity to `λ_bg`, so the confirmed-death data can identify a
   genuinely large non-BVD background without forcing `q_death` to ~1.
2. Failing that, relate the death testing/forwarding rate to the case rate
   with an enrichment factor (issue #206's proposal), which changes the
   level through `τ_death` rather than positivity; this is a weaker fix
   because it does not address the `q_death` pinning.

This is a structural change to the death positivity, not a prior tweak, so
it should be implemented and re-fitted as a separate piece of work rather
than landed here.

## Relationship to draft #210 (`deaths-testing`)

PR #210 removed `τ_death` and made deaths share the case capacity queue.
That over-drained the queue and over-predicted confirmed deaths more, not
less.
The informed-background route investigated here does not over-predict more
but also does not fix the level, and it breaks mixing.
Neither the shared-queue (#210) nor a wider background prior is the right
path; the composition-link/enrichment structure (recommendation 1) is the
candidate that addresses the actual mechanism (the rigid `q_death`).

## Reproduce

```
julia --project=. scripts/deaths_overprediction.jl baseline
julia --project=. scripts/deaths_overprediction.jl wide
julia --project=. scripts/deaths_overprediction.jl centred
```

Each invocation fits one joint and appends a result block to
`scripts/deaths_results.txt` (flushed per line).
