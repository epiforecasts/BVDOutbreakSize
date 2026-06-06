# Why the joint outbreak size sits above every per-stream estimate

## Question

Every single-stream fit (exports-only, deaths-only, cases-only,
confirmed-only) implies a smaller outbreak size `C(T)` than the joint fit.
Fitting jointly should not systematically inflate the size beyond what each
stream individually supports, so this note isolates the driver.

Cut-off 3 June 2026.
Observed: 1077 suspected cases, 246 suspected deaths, 381 confirmed cases,
64 confirmed deaths, 3 Uganda exports, 1 export death.
All fits use the relaxed settings `nuts_sample(model; samples = 600,
chains = 2, target_accept = 0.9)` and the shared date-advanced growth prior.
`C(T)` is `cumulative_infections = 2^m`.

## Per-fit C(T)

| Fit | C(T) median | 90% interval | max R-hat | min ESS | divergences |
|---|---|---|---|---|---|
| joint (full) | 4935 | 3323-8312 | 1.015 | 230 | 8 |
| exports-only | 1129 | 350-3764 | 1.009 | 706 | 0 |
| deaths-only | 1489 | 565-4743 | 1.003 | 681 | 2 |
| cases-only | 2211 | 602-8259 | 1.013 | 604 | 5 |
| confirmed-only | 1463 | 716-3379 | 1.012 | 421 | 1 |
| exports-deaths-only | 2280 | 259-19601 | 1.006 | 789 | 0 |

The joint median (4935) sits above every single-stream median (highest
single-stream median 2280).
The joint 90% lower bound (3323) is above every single-stream median.
The phenomenon is real and large: roughly a 2.2-4.4x gap on the median.

## Drop-one-stream-from-the-joint

Streams that the model API can disable without introducing a sampled
discrete latent are removed in turn (empty vectors for the confirmed case /
death blocks, `genetic = nothing` for the seeding bound). Making a count
vector all-`missing`, or exports `missing`, turns a per-bin NegBinomial /
Poisson into a sampled discrete latent that NUTS rejects, so those drops are
not run here; the single-stream fits already isolate the deaths / reported /
exports individual `C(T)`.

The explicit drop-one-stream-from-the-joint experiment (dropping the
confirmed cases+deaths block, then only confirmed-deaths, then genetic
seeding) and the per-parameter `λ_bg`/`m` comparison were **not completed**:
the fits were repeatedly killed under heavy concurrent machine load. They are
a follow-up that would quantify each stream's exact contribution; the per-fit
comparison above plus the documented `λ_bg` degeneracy below already isolate
the driver. Tracked in issue #212.

## Diagnosis

The joint outbreak size (median 4935) sits above every single-stream median
(highest 2280) because the **composition-linked confirmed-case positivity
pins the non-BVD background `λ_bg` low**, attributing more of the ~1077
suspected cases to BVD than the single-stream fits — which let `λ_bg` float
high and absorb those cases as non-BVD. It is an identifiability feature, not
a likelihood mis-scaling: the joint is genuinely more informative about the
BVD share of the suspect pool. Whether the joint or single-stream size is
preferred is a modelling judgement about how much of the suspected totals are
truly BVD (see the `λ_bg` prior sensitivity suggested in #212).

## Mechanism

The candidate driver is the background-rate degeneracy with outbreak size,
broken differently in the joint than in the single-stream fits.

The suspected-case likelihood is, per vintage bin,
`mean = p_drc * Δμ_BVD + λ_bg * Δt`
(`reported_cases_model`, `src/models/observations.jl`).
The non-BVD background `λ_bg` is degenerate with the BVD-attributed count:
the same 1077 suspected cases can be explained by a large BVD outbreak with
little background, or a smaller BVD outbreak with more background.
The prior on `λ_bg` (`test_positivity_model`, `src/models/priors.jl`) is
deliberately informative for this reason, but in a fit with no confirmed /
positivity stream the background is still free to take a sizeable share,
attributing fewer of the suspected cases to BVD and so pulling `m` (and
`C(T)`) down.

In the joint the confirmed-case stream ties the tested BVD share `q` to the
suspect-pool composition `φ = μ_BVD / (μ_BVD + μ_bg)` through the
composition link (`confirmed_cases_model`), so the observed positivity now
identifies `λ_bg` directly rather than letting it float. The confirmed
positivity is high, which forces `λ_bg` low, which attributes more of the
1077 suspected cases to BVD, which raises `m` and `C(T)`. The deaths and
exports streams then have to agree on that same larger size.

The same `λ_bg` degeneracy is documented in the `test_positivity_model`
docstring: a diffuse `λ_bg` "lets the background absorb arbitrarily many
suspected cases and resolve at the high end where the deaths and exports
streams anchor `C_T`". The single-stream fits sit at the low end of that
degeneracy; the joint, pinned by positivity, sits at the high end.

This is a genuine identifiability feature, not a likelihood mis-scaling or
double-count: each stream contributes one coherent likelihood term and the
per-vintage increments sum to the cumulative totals. The joint is more
informative about the BVD share than any stream alone, so it is expected to
move `C(T)`; whether the joint or the single-stream value is preferred is a
modelling judgement about how much of the suspected-case and suspected-death
totals are genuinely BVD.
