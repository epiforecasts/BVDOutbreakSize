# Patch (meta-population) model for the 2026 DRC BVD outbreak

## Motivation

The headline model treats the DRC outbreak as a single well-mixed population.
Transmission is spatially structured: the outbreak began in Ituri Province and
has since spread to Nord-Kivu and, marginally, Sud-Kivu.
The INSP sitreps publish per-province spatial tables (Tableau 1) giving
confirmed cases by province.
A patch model lets us estimate province-specific reproduction numbers and
per-province outbreak sizes while still fitting every national data stream.

## What the province data can and cannot support

Before choosing a structure, look at what the spatial tables actually contain
over the fitted window (18 June to 6 July, 17 vintages).

| | 18 Jun | 6 Jul | growth |
|---|---|---|---|
| Ituri | 853 | 1556 | 1.82x |
| Nord-Kivu | 77 | 149 | 1.94x |
| Sud-Kivu | 3 | 3 | 1.00x |
| Ituri share | 91.4% | 91.1% | — |

Three facts follow, and they drive every design decision below.

1. **The provincial shares are flat over this window.**
Ituri's share of confirmed cases sits between 91.0% and 91.4% at every one of
the 17 vintages, and Ituri and Nord-Kivu have grown at very nearly the same
rate.
This is a *finding*, not a licence to hard-code a constant Rt ratio: the window
is short, the response is concentrated on Ituri (so divergence is expected), and
a model that cannot represent divergence cannot detect it either.
The model therefore lets the provincial Rt trajectories separate and *estimates*
how much they do (`region_drift_sd`), rather than assuming they cannot.
See "Reproduction numbers" below.

2. **Sud-Kivu is static.**
A constant 3 confirmed cases, so zero increments in every interval.
It informs the model that transmission there is negligible, and nothing else.

3. **The province counts are an exact partition of the national counts.**
At every vintage, Ituri + Nord-Kivu + Sud-Kivu equals the national confirmed
total that the existing `confirmed_cases_model` already fits.
This is verified, not assumed.

## Latent process

### Per-patch renewal

For each patch `p` on day `t`:

```
I_{p,t} = R_{p,t} · Σ_s I_{p,t−s} · g_s  +  ε · Σ_q K_{p,q} · I_{q,t−1}
```

with `g_s` the shared generation-interval PMF (the biology of transmission does
not depend on province, so it is sampled once) and the incubation PMF likewise
shared.
`patch_infections()` in `renewal.jl` implements this and reduces exactly to
`renewal_infections()` per row when the kernel is zero.

### Reproduction numbers: common trend + correlated MVN deviations

`patch_rt_model` builds the per-patch Rt from the existing national
weekly-knot random walk (`rt_walk_model`, unchanged) plus per-patch
deviations that vary in space AND in time, drawn from a multivariate-normal
random walk with a learned cross-patch correlation:

```
log R_{p,t} = mu(t) + delta_p(t),          sum_p delta_p(t) = 0
Delta delta(t_k) ~ MVN(0, Sigma)
Sigma = diag(sigma_delta) . Omega . diag(sigma_delta),   Omega ~ LKJ(2)
```

Three choices here are deliberate.

**Rt varies in time as well as in space, and by how much is estimated.**
An earlier version of this branch fixed the per-patch modifier to be constant
in time, on the grounds that the confirmed-case shares are flat over the
observed window.
That was the wrong direction of inference: the flatness of the shares is a
fact to be *estimated*, not a constraint to impose.
A constant modifier makes "the provinces share one temporal Rt shape" a
structural assumption, so the model could never detect divergence even once it
appeared — and the window is short.
It also forbids the leading hypothesis outright: the response is concentrated
on the Ituri epicentre, so Ituri's Rt falling faster than Nord-Kivu's is
exactly what one would expect, and a constant modifier cannot represent it.
`sigma_delta -> 0` recovers the constant modifier exactly, so nothing is lost
if the provinces really are moving together.
`sigma_delta` is therefore the headline spatial diagnostic (surfaced as
`region_drift_sd`): a posterior pushed away from zero is direct evidence that
provincial Rt trajectories are separating.

**A common trend, not independent walks.**
The model keeps the national walk `mu(t)` as the pooling target.
An unstructured multivariate walk — one free log-Rt trajectory per province,
correlated only through an LKJ prior, with no common trend — is strictly more
flexible, but it shrinks the wrong way.
`LKJ(eta)` has density proportional to `det(Omega)^(eta-1)`, maximised at
`Omega = I`, so `eta = 2` mildly favours *independent* provincial walks.
The provinces are not equally observed: over the fitted window Nord-Kivu
contributes 74 laboratory positives and Sud-Kivu contributes none at all.
Shrinking toward independence estimates their Rt almost entirely from that;
shrinking toward a common trend lets them borrow strength from Ituri and
deviate only where the data insist.
For a meta-population under one national response, partial pooling toward a
shared trend is the right inductive bias.

This costs nothing in generality.
With `mu(t)` present the deviation covariance `Sigma` is still free, so the
cross-patch correlation is *learned* rather than assumed.
The model **is** a multivariate-normal random walk; it simply carries a common
factor rather than leaving the correlation structure to carry it.
(Writing the increments out: `Cov(dlog R_p, dlog R_q) = sigma_rw^2 + Sigma_pq`,
a one-factor-plus-free-residual covariance.)

**Sum-to-zero, not a reference patch.**
The deviations are centred at every knot, so no province is privileged.
Fixing `delta_1 = 0` (reference coding) also identifies the model, but it
forces the primary patch to have no idiosyncratic deviation at all: Ituri
would *be* the national trend by construction while the other provinces carry
their own noise.
That asymmetry is an artefact of the identifiability fix, not epidemiology.
Centring introduces one redundant coordinate per knot (the mean of the raw
innovations, which the likelihood never sees); it is drawn from its proper
prior and is Gaussian and well-conditioned, so it costs a few cheap sampled
dimensions rather than a posterior ridge.

The reported national `R_T` is the incidence-weighted aggregate implied by the
summed patch infections (`implied_national_Rt`), which inverts the renewal
equation on the total.
This is provably the force-of-infection-weighted mean of the patch Rts, so it
is the reproduction number that reproduces the national trajectory.

**What the data can identify.**
The per-province composition identifies the *contrast* between provinces.
With three patches that is essentially one number, the Ituri / Nord-Kivu
contrast, since Sud-Kivu carries no signal.
Expect `Omega` to be largely prior-driven and Sud-Kivu's Rt to be pinned by
the deviation prior rather than by data.
That is honest, and it is why `Sigma` gets a proper shrinkage prior rather
than a flat one.

### Importation: off by default, and honest about why

There is **no mobility or origin-destination data** for this outbreak.
Worse, importation is not identifiable here even in principle: with flat shares
(fact 1), any pair of (importation rate `ε`, secondary-patch seed) that
reproduces the observed Nord-Kivu level fits equally well.
The two are confounded.

So `ε` is sampled **only when a non-zero `importation_kernel` is supplied**.
The default kernel is all-zero, the patches are uncoupled, and no `ε` enters the
parameter space.
Sampling `ε` against a zero kernel — as an earlier version of this branch did —
adds a dimension the likelihood never touches, giving a parameter whose
posterior is exactly its prior and which only slows the sampler.

With the default, each secondary patch is explained by its own sampled seed
(standing in for the unobserved introductions from Ituri) and its own `R_t`.
The kernel machinery is retained and tested, so a mobility-informed kernel can
be passed as a sensitivity analysis if data ever appear.

## Data streams

### National (unchanged)

Every national stream is fitted exactly as in `bvd_joint`, against the summed
patch onsets: suspected cases, suspected deaths, confirmed cases, confirmed
deaths, the laboratory pipeline, treatment flows and recoveries.
Uganda exports are driven by the **primary patch alone**, since the border
crossings that stream describes are Ituri-to-Uganda.

### Province-level: a composition, not a second count likelihood

Because the province counts are an exact partition of the national counts
(fact 3), fitting them with their own count likelihood would put the same
observations into the joint density **twice**, double-weighting the confirmed
stream against every other stream in the model.
The first version of this branch did exactly that.

The information the spatial tables add over the national series is purely
*spatial*.
So the likelihood is factorised:

```
P(y_1, y_2, y_3) = P(N) · P(y_1, y_2, y_3 | N),    N = Σ_p y_p
```

The existing national `confirmed_cases_model` supplies the total term `P(N)`.
`province_composition_model` supplies **only** the conditional composition term,
scored by stick-breaking over the patches with an overdispersed Binomial
(`safe_betabinomial`), the sequential form of a Dirichlet-multinomial.
The vintage totals are conditioned on and never scored, so nothing is counted
twice.
A shared overdispersion `ρ` absorbs extra-Binomial variation in how cases are
attributed to provinces (reporting lag between the provincial and national
tables, reassignment between health zones).

The test suite pins this down directly: the composition's log-density is
invariant to rescaling the modelled confirmed level (it sees only normalised
shares) but responds to a change in the split.

## The province-invariant ascertainment assumption

The composition consumes the modelled per-patch confirmed cases only through
*normalised shares*, so every factor common to all provinces cancels: test
sensitivity, ascertainment, background, test positivity.
That cancellation is what stops the composition re-scoring the national total.
It also encodes an assumption, and it is the model's main soft spot.

**The probability that a true infection becomes a confirmed case is assumed
the same in every province.**

If Nord-Kivu ascertains a smaller fraction of its infections than Ituri, its
confirmed share understates its true share of infections, and the model cannot
tell that apart from a genuinely lower Nord-Kivu `Rt`.
The bias lands in `δ_2`, and the per-patch `C_T` split is skewed toward Ituri.
Given that Ituri is the epicentre with the response concentrated there, a
lower Nord-Kivu ascertainment is plausible, so this is not a hypothetical.

**The fix needs a per-province denominator**, and one may exist.
The sitreps appear to carry per-province laboratory throughput (samples
analysed by province, section 4.3 / Tableau 6), which is exactly the
denominator that would identify a province-specific test positivity and so
separate provincial ascertainment from provincial `Rt`.
It has not been scanned: `data/insp_sitrep_scanned.csv` currently holds only
national `samples_analysed`.
Extracting it is the single highest-value addition to this model.
Until then, read the per-patch `C_T` split as conditional on province-invariant
ascertainment, and the national headline (which does not depend on the split)
as unaffected.

## Identifiability summary

| Quantity | Identified by | Status |
|---|---|---|
| `δ_2` (Nord-Kivu log-Rt offset) | level of the confirmed-case split | identified, but **confounded with province ascertainment** |
| `δ_3` (Sud-Kivu) | zero increments throughout | weakly identified, pinned low |
| secondary-patch seeds | early provincial level | identified given `ε = 0` |
| `ε` (importation) | nothing | **not identified**; off by default |
| `ρ` (composition overdispersion) | vintage-to-vintage scatter in shares | identified |
| national `C_T`, `R_T` | all national streams, as before | unaffected by the split |

## Status

- [x] `patch_infections()`, `importation_from_kernel()`, `implied_national_Rt()`
  in `renewal.jl`
- [x] `patch_rt_model()` (reference-coded constant modifiers) in `priors.jl`
- [x] `patch_infection_model()` in `priors.jl`
- [x] `province_composition_model()` in `observations.jl`
- [x] `_patch_latent()` and `bvd_patch_joint()` in `joint.jl`
- [x] Per-province data pipeline: `province_confirmed_history` in
  `observations.toml` → `load_observations()` → `province_increment_matrix()`
- [x] Headline quantities (`C_T`, `R_T`, `r`, `T`, `CFR`, `R0`,
  `doubling_time`) surfaced under the same names as `bvd_joint`, so a patch
  chain drops into `summary_table` and the existing reporting unchanged
- [x] `patch_summary_table()` for per-patch outbreak sizes and Rts
- [x] Test suite (`test/test_patch_model.jl`)

## Next steps

1. **Scan the per-province laboratory throughput** from the sitreps. This is
   the highest-value item: it is the denominator that separates provincial
   ascertainment from provincial `Rt`, and without it `δ_2` carries both.
2. Per-location forecasts (`forecast_patch`): project each patch over the
   horizon using the national walk's terminal drift plus the constant `δ_p`.
3. Forecast-vs-truth against the per-province spatial tables.
4. Wire the patch model into `docs/examples/analysis.jl` as a spatial section
   alongside the national headline fit.
5. Further province-level streams (isolation occupancy by province, Tableau 6)
   need the same composition treatment, since they too are partitions of
   national totals the model already fits. Do not add them as raw counts.
