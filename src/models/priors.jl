# Prior and latent-state submodels: the building blocks shared across the
# observation submodels and the joint composer. Each `@model` is one piece
# of the generative process — a single prior, a delay, the reproduction
# number, the generating infection process, or the onset staging — so it
# can be reused across composers without duplication. Delays are sampled
# from priors and discretised with CensoredDistributions; nothing is fixed.

## --- Delay submodels (priors only; all delays sampled) ------------------

"""
Generic delay submodel parameterised by mean and SD, discretised to a
daily PMF over lags `0 … nmax` by double interval censoring of a
moment-matched LogNormal (see [`lognormal_meansd`](@ref) and
[`discretise_censored`](@ref)). The LogNormal CDF differentiates cleanly
under Mooncake, so this is the AD-safe discretisation route for every
delay in the renewal convolutions. The mean and SD carry
weakly-informative priors, so the delay is estimated rather than fixed.
Returns `(; pmf, dist, mean, sd)`.
"""
@model function censored_delay_model(nmax::Integer; mean_prior, sd_prior)
    delay_mean ~ mean_prior
    delay_sd ~ sd_prior
    dist = lognormal_meansd(delay_mean, delay_sd)
    return (; pmf = discretise_censored(dist, nmax), dist,
        mean = delay_mean, sd = delay_sd)
end

"""
Generation-interval submodel. The mean and SD are sampled from priors
centred on the Ebola virus disease serial interval as a generation-time
proxy (mean 15.3 d, SD 9.3 d; WHO Ebola Response Team 2014, NEJM), so the
generation time is estimated around the published value rather than
fixed. Both the mean and SD are sampled, so the published uncertainty
propagates into the generation interval. Discretised with
[`censored_delay_model`](@ref); the lag-0 bin is dropped and the remainder
renormalised, left-truncating the generation interval at one day so an
infectee is infected strictly after its infector. Returns
`(; g, gi_mean, gi_sd)`.
"""
@model function generation_interval_model(nmax::Integer;
        mean_prior = truncated(Normal(15.3, 3.0); lower = 1),
        sd_prior = truncated(Normal(9.3, 2.0); lower = 1))
    d ~ to_submodel(censored_delay_model(nmax; mean_prior, sd_prior))
    g = d.pmf[2:end] ./ sum(d.pmf[2:end])
    return (; g, gi_mean = d.mean, gi_sd = d.sd)
end

## --- Reproduction number ------------------------------------------------

"""
Weekly piecewise-linear log-scale reproduction number over `n` days, with
a smooth intervention ramp. Knots sit at weekly spacing
([`knot_days`](@ref)) and follow a Gaussian random walk in non-centred
cumulative-sum form: standard-normal innovations are scaled by `sigma_rw`
and accumulated, avoiding the funnel geometry of the centred recursion and
matching the non-centred ascertainment block. Daily log-`R_t` is the
linear interpolation between knots ([`interpolate_knots`](@ref)). An
intervention at `breakpoint` (e.g. the first WHO situation report) adds a
sampled effect `intervention_effect` shaped by a logistic ramp
([`sigmoid_ramp`](@ref)), so transmission changes gradually over the ramp
window rather than instantly; `breakpoint = missing` drops the term. The
ramp scale `ramp` defaults to 21 days, an about-three-week transition
reflecting that a response (case finding, isolation, vaccination) takes
weeks to bite rather than switching at a single date; pass `ramp` to widen
or narrow it.
`Rt = exp.(log_Rt)`. Returns
`(; Rt, log_R, days, sigma_rw, log_R0, intervention_effect)`.

The initial reproduction number prior is anchored on the molecular-clock
growth estimate for this outbreak. Cuomo-Dannenburg & Ghafari's
phylodynamic reanalysis of the first ten BDBV genomes puts the mean
epidemic doubling time at 15.2–24.5 d (centre 20 d). With the renewal
generation interval that 20-day doubling implies `R0 ≈ 1.6`, and the
15.2–24.5 d range maps to `R0 ≈ 1.47–1.84`, so the prior is
`Normal(log(1.6), 0.10)`: centred on the molecular-clock doubling with a
log-SD slightly wider than the genetic spread (≈0.07) to allow for the
wide per-assumption clock intervals. Sampling on `R0` (with `r0` derived
through Euler–Lotka) rather than on `r` directly keeps the renewal seeding
consistent; the implied doubling-time prior tracks the genetic estimate.
This molecular-clock prior anchors the ESTABLISHED reproduction number (the
walk base at the genetic bound). The pre-establishment seeding `R_t` is a
SEPARATE, low and wide prior (`seed_log_r_prior`, default `Normal(log(1),
0.5)`, centred on no growth with broad support over sub-critical to slow
growth), so the genetic growth estimate for the established outbreak does
not act on the introduction phase and push the seeding time earlier.

The random-walk step SD prior is a tight half-normal SD 0.05. The
observed sitrep window is only the final ≈9 days of a ≈90-day inferred
outbreak, so the weekly log-`R_t` knots over the unobserved stretch are
free to drift; with a wider step the walk climbed from `R0 ≈ 1.9` to a
terminal `R_T ≈ 3.1` (a faster terminal doubling, 5.5 d, than the observed
7.5–8.8 d), inflating `C_T` through a large pool of recent, not-yet-
observed infections. A tight step keeps `R_t` near the seeding value so
the outbreak size is not driven by an unsupported terminal acceleration.

The intervention effect is constrained to be non-positive
(`truncated(Normal(0, 0.4); upper = 0)`): a declared WHO response (case
finding, isolation, vaccination) can only reduce transmission or leave it
unchanged, never increase it, so the ramp's effect on log-`R_t` is bounded
at or below zero. This is stronger than the earlier symmetric
`Normal(0, 0.5)` (under which the effect was unidentified and `R_t`
drifted up unchecked) and the one-sided lean `Normal(-0.3, 0.4)`: the
response now has a definite non-increasing effect, with the half-normal
admitting anything from no effect (mode) to a substantial decline. Because
the breakpoint is only ≈11 days before the cut-off and the ramp is a
fortnight, the response damps `R_t` only partially by the cut-off.
"""
@model function rt_walk_model(n::Integer;
        week::Integer = 7,
        breakpoint::Union{Missing, Real} = missing,
        rt_start::Integer = 1,
        ramp::Real = 21.0,
        log_r0_prior = Normal(log(1.6), 0.10),
        seed_log_r_prior = Normal(log(1.0), 0.5),
        sigma_prior = truncated(Normal(0, 0.05); lower = 0),
        effect_prior = truncated(Normal(0, 0.4); upper = 0))
    days = knot_days(n; week, start = rt_start)
    nb = length(days)
    ## Two distinct reproduction numbers: a LOW seeding `R_t` over the
    ## pre-establishment window (before the genetic bound, days < rt_start),
    ## and the established `R0` at the bound that the random walk grows from.
    ## Keeping them separate stops the molecular-clock growth prior (which
    ## describes the *established* outbreak) from acting on the introduction
    ## phase and pushing the seeding time earlier.
    log_R0 ~ log_r0_prior
    log_R0_seed ~ seed_log_r_prior
    sigma_rw ~ sigma_prior
    z ~ product_distribution(fill(Normal(0, 1), max(nb - 1, 1)))
    intervention_effect ~ effect_prior
    steps = sigma_rw .* z[1:(nb - 1)]
    log_R = log_R0 .+ vcat(zero(log_R0), cumsum(steps))
    log_Rt_walk = interpolate_knots(log_R, days, n)
    log_Rt = [t < rt_start ? log_R0_seed : log_Rt_walk[t] for t in 1:n]
    log_Rt = log_Rt .+ intervention_effect .* sigmoid_ramp(n, breakpoint; ramp)
    Rt = exp.(log_Rt)
    return (; Rt, log_R, days, sigma_rw, log_R0, log_R0_seed,
        intervention_effect)
end

## --- Seeding and the generating infection process -----------------------

"""
Seed submodel: the latent infection count `I0` on the last day of the
seeding window, representing the zoonotic introduction. Default prior is a
truncated Normal centred on a single seed; the prior is injectable. The
seeding window is filled by exponential growth at the implied rate in
[`infection_model`](@ref).

The lower truncation is `0.02` rather than `0`. The implied outbreak start
`T` ([`seeding_age`](@ref)) is a monotone, saturating function of `I0`: as
`I0 → 0` the whole seeded trajectory shrinks toward zero and the cumulative
never crosses one over the unobserved seeding stretch, collapsing into a
degenerate tiny-outbreak / very-early-start mode (`C_T ≈ 3`, `T` pinned at
the grid edge). Bounding `I0` away from zero removes that mode without
biasing the start: in the infection-plus-genetic submodel IN ISOLATION the
implied `T` is then wide and unimodal (median ≈ 82 d, 95% ≈ 64–93 d on the
report grid, no grid-edge pile-up, no divergences). The exponential-forcing
renewal seed cannot make the start a free parameter decoupled from outbreak
size — `T` is a monotone, saturating function of `I0` and the implied growth
`r0` — so the start is controlled through this `I0` bound and the
molecular-clock `T` bound rather than a direct `T` prior (a likelihood term
on the derived crossing funnels the geometry badly and was rejected). In the
FULL joint, the streams that anchor outbreak size (deaths, exports) can
still pull `I0` / `R_t` up enough to drive `T` toward the grid edge: the
edge-pinning is then a symptom of the size streams demanding a large `C_T`,
not of the seed prior, and is addressed by the size-stream tensions (the
exports truncation, the background identification) rather than the seed
alone. Fitting `T` explicitly (e.g. via a Gibbs step on a seeding day) is
the alternative when a fully size-decoupled start is required.
"""
@model function seed_model(; i0_prior = truncated(Normal(0.1, 0.1); lower = 0.02))
    I0 ~ i0_prior
    return (; I0)
end

"""
Generating infection process: the latent submodel whose expected daily
infections every downstream stream consumes. Samples the reproduction
number trajectory, the generation interval and the seed via injected
submodels, derives the seeding-window growth rate from the initial
reproduction number through the Euler–Lotka relation
([`euler_lotka_r`](@ref)), seeds the first `length(g)` days as exponential
growth ([`seed_infections`](@ref)), then runs the discrete renewal
recursion ([`renewal_infections`](@ref)). Replaces the integral model's
exponential-growth trajectory. The `breakpoint` is forwarded to the
reproduction-number submodel so the intervention ramp lands on the right
day. Exposes the daily infection incidence, its cumulative sum, and the
reported growth summaries: the current growth rate `r` and
`doubling_time`, the implied initial growth rate `r0`, and the outbreak
age `T` (the seeding-to-cut-off time, derived as the smooth crossing where
cumulative infections reach one; see [`seeding_age`](@ref)). Requires a
grid of `n ≥ 2` days for the current growth rate (the cut-off grid always
spans the seeding window); a single-day grid falls back to zero growth.
Returns `(; infections, cumulative, Rt, g, I0, r0, r, T, C_T, doubling_time)`.
"""
@model function infection_model(n::Integer;
        breakpoint::Union{Missing, Real} = missing,
        rt_start::Integer = 1,
        rt = rt_walk_model,
        gi = generation_interval_model,
        seed = seed_model,
        gi_nmax::Integer = 40,
        seed_len::Union{Nothing, Integer} = nothing)
    rt_state ~ to_submodel(rt(n; breakpoint, rt_start))
    gi_state ~ to_submodel(gi(gi_nmax))
    seed_state ~ to_submodel(seed())
    Rt = rt_state.Rt
    g = gi_state.g
    r0 = euler_lotka_r(Rt[1], g)
    ## Seeding-window length. Defaults to the generation-interval support
    ## (`length(g)`), the natural lookback the renewal recursion needs. A
    ## shorter window concentrates the seed near the establishment day and
    ## reduces the flat low-count stretch over which the cumulative=1
    ## crossing (the implied outbreak start) wanders, which otherwise admits
    ## a second, very-early start regime.
    L = seed_len === nothing ? length(g) : min(Int(seed_len), length(g))
    seed_vec = seed_infections(seed_state.I0, r0, L)
    infections = renewal_infections(Rt, g, seed_vec)
    cumulative = cumsum(infections)
    ## Current growth rate from the last two days; n ≥ 2 (the cut-off grid
    ## always spans the seeding window, so this holds in practice), with a
    ## zero-growth fallback on a degenerate single-day grid.
    r = n >= 2 ?
        log(safe_rate(infections[n])) - log(safe_rate(infections[n - 1])) :
        zero(eltype(infections))
    return (; infections, cumulative, Rt, g, I0 = seed_state.I0, r0, r,
        T = seeding_age(cumulative, n), C_T = cumulative[n],
        doubling_time = doubling_time(r))
end

"""
Onset-incidence submodel: convolve the renewal infections with the
sampled incubation-period PMF to get daily symptom-onset incidence.
Computed once per draw and reused by every downstream observation stream,
so the staging infections → onsets → each observed event is explicit. The
incubation delay submodel is injected. The incubation period cannot be
fitted from the BDBV line list (no exposure dates), so the line-list
reanalysis recommends the MacNeil et al. (2010) Bundibugyo estimate from
the 2007 Uganda outbreak: mean 6.3 d (95% CI 5.2-7.3, n = 24). The mean
prior `Normal(6.3, 0.54)` reproduces MacNeil's reported 95% CI (SD = CI
half-width / 1.96); MacNeil give no interval on the spread, so the SD
prior is a weakly-informative modelling choice centred on the
CV-implied spread (≈ 3.5 d). Returns
`(; onsets, incubation_pmf, incubation_mean, incubation_sd)`.
"""
@model function onset_incidence_model(infections::AbstractVector;
        incubation = (nmax) -> censored_delay_model(nmax;
            mean_prior = truncated(Normal(6.3, 0.54); lower = 1),
            sd_prior = truncated(Normal(3.5, 0.8); lower = 1)),
        incubation_nmax::Integer = 30)
    inc_state ~ to_submodel(incubation(incubation_nmax))
    onsets = convolve_delay(infections, inc_state.pmf)
    return (; onsets, incubation_pmf = inc_state.pmf,
        incubation_mean = inc_state.mean, incubation_sd = inc_state.sd)
end

## --- Genetic seeding bound ----------------------------------------------

"""
One-sided molecular-clock seeding bound on the outbreak age `T` (see
[`infection_model`](@ref)). The TMRCA is treated as a right-censored,
noisy reading of the seeding time, so deeper or wider sampling only pushes
it older; the likelihood contributes `P(read ≥ tmrca_days)`. Passing
`tmrca_days = missing` makes the submodel a no-op.
"""
@model function genetic_seeding_model(T::Real,
        tmrca_days::Union{Missing, Real}; tmrca_days_sd::Real = 15.0)
    if !ismissing(tmrca_days)
        tmrca_days ~ censored(Normal(T, tmrca_days_sd); upper = tmrca_days)
    end
    return (; T, tmrca_days_sd)
end

## --- Shared nuisance priors ---------------------------------------------

"""
Case-fatality ratio prior. Default `Beta(6.6, 13.4)` has mean ≈ 0.33,
matching the CDC summary for past BVD outbreaks. Used by the deaths and
deaths-among-exports streams.
"""
@model function cfr_model(; cfr_prior = Beta(6.6, 13.4))
    CFR ~ cfr_prior
    return (; CFR)
end

"""
Prior on the mean daily traveller volume from the source area to Uganda.
Default centred on `ITURI_DAILY_TRAVEL` with SD `ITURI_DAILY_TRAVEL_SD`,
truncated at zero. Sets the per-capita travel rate for the exports stream.
"""
@model function traveller_volume_model(;
        mean::Real = ITURI_DAILY_TRAVEL,
        sd::Real = ITURI_DAILY_TRAVEL_SD)
    daily_travellers ~ truncated(Normal(mean, sd); lower = 0)
    return (; daily_travellers)
end

"""
Non-BVD background rate for the suspected-death stream
([`deaths_model`](@ref)), the death analogue of the suspected-case
background `λ_bg` ([`test_positivity_model`](@ref)). The DRC sitrep
suspected-death definition is symptomatic-then-deceased, so a death need
not be a true BVD death; this submodel samples the per-day non-BVD
background death rate `λ_bg_death`. Its cumulative contribution over the
grid is `λ_bg_death · n`.

The default `truncated(Normal(0, 0.25); lower = 0)` is deliberately
informative, mirroring the case background: deaths are far fewer than
suspected cases (≈ 246 suspected deaths at the cut-off vs ≈ 1077
suspected cases), so the background rate is scaled down accordingly. With
SD 0.25 the median background is ≈ 0.17/day (≈ 22 deaths over a ≈ 132-day
grid, a modest minority of the suspected-death total) while still
admitting a genuine non-BVD signal. The background is degenerate with
outbreak size, so a diffuse prior would let it absorb arbitrarily many
suspected deaths. Pass `lambda_prior` to override. Returns
`(; λ_bg_death)`.
"""
@model function death_background_model(;
        lambda_prior = truncated(Normal(0.0, 0.25); lower = 0))
    λ_bg_death ~ lambda_prior
    return (; λ_bg_death)
end

"""
Test-positivity machinery shared by the suspected- and confirmed-case
streams. Samples

- `λ_bg` — the per-day non-BVD background suspected-case rate, on a
  half-normal scale. Drives the suspected/confirmed contrast: suspected
  cases mix the BVD onset-to-report signal with this additive background,
  while the laboratory pipeline only confirms the BVD share.
- `τ_test` — the fraction of suspected cases that are sampled and routed
  to the laboratory pipeline.

The default `λ_bg` prior is `truncated(Normal(0.5, 0.3); lower = 0)`, a
half-normal with its CENTRE shifted above zero. Its total contribution to
the expected suspected-case count over the grid is `λ_bg · T`, with `T`
the seeding-to-cut-off span. The prior is deliberately informative and
non-zero-centred because `λ_bg` is degenerate with outbreak size (the
per-vintage reported mean mixes the `p_drc`-scaled BVD increment with
`λ_bg · Δt`): a diffuse prior lets the background absorb arbitrarily many
suspected cases and resolve at the high end where the deaths and exports
streams anchor `C_T`, while a prior with its mode AT zero leaves a second
mode that collapses the background to zero (the suspected-death definition
and the laboratory's large negative fraction imply a genuine non-BVD
signal exists). Centring at 0.5/day with SD 0.3 pins a sane non-zero
background — mode 0.5/day, 95% bound ≈ 1.1/day, ≈ 50 of the ≈ 1077
suspected cases over a ≈ 100-day grid — a modest minority that still
admits a real non-BVD signal but cannot explain the majority of suspected
cases. Earlier diffuse SD-1 / SD-5 priors left a second posterior mode in
which the background explained most suspected cases (positivity ≈ 0.2);
the tightened above-zero centre removes that mode. Pass `lambda_prior` to
override. `τ_test` defaults to `Beta(5, 2)` (mean ≈ 0.71).

The derived per-suspected positivity is exposed inside
[`reported_cases_model`](@ref); the per-test positivity inside
[`confirmed_cases_model`](@ref). Returns `(; λ_bg, τ_test)`.
"""
@model function test_positivity_model(;
        lambda_prior = truncated(Normal(0.5, 0.3); lower = 0),
        fraction_tested_prior = Beta(5.0, 2.0))
    λ_bg ~ lambda_prior
    τ_test ~ fraction_tested_prior
    return (; λ_bg, τ_test)
end

"""
Per-vintage non-BVD background rate as a partially-pooled, non-centred
random effect, the time-varying generalisation of the scalar `λ_bg` /
`λ_bg_death`. Used by the suspected-case ([`reported_cases_model`](@ref))
and suspected-death ([`deaths_model`](@ref)) streams when their
`background_re` switch is on. The same non-BVD reporting environment
plausibly drives both streams, so the two backgrounds can share this
submodel's hyperparameters.

The baseline `λ_mu` is the scalar background rate on its natural
half-normal scale, with the same informative above-zero-centred default as
the scalar `λ_bg` (`truncated(Normal(0.5, 0.3); lower = 0)` for cases;
pass a tighter `baseline_prior` for deaths). The per-vintage rate is a
multiplicative
log-normal deviation from this baseline,

```math
\\lambda_v = \\lambda_\\mu \\,
    \\exp(\\sigma_{bg}\\, z_v), \\qquad z_v \\sim \\mathcal N(0, 1),
```

with `σ_bg` the pooling SD, passed in rather than sampled here so the
suspected-case and suspected-death streams can SHARE one pooling SD (the
same non-BVD reporting environment drives both); see
[`background_pooling_model`](@ref), which samples it once at the composer
level. The deviation is multiplicative so the per-vintage rate stays
positive without a clamp and `σ_bg → 0` recovers the scalar baseline
exactly (every `λ_v = λ_mu`). Each stream still samples its own baseline
`λ_mu` and per-vintage deviations `z`. `nv` is the number of vintage
windows. Returns `(; λ, λ_mu, σ_bg, z)` with `λ` a length-`nv` vector of
per-vintage rates.
"""
@model function background_re_model(nv::Integer, σ_bg::Real;
        baseline_prior = truncated(Normal(0.5, 0.3); lower = 0))
    m = max(nv, 1)
    λ_mu ~ baseline_prior
    z ~ product_distribution(fill(Normal(0, 1), m))
    λ := λ_mu .* exp.(σ_bg .* z[1:nv])
    return (; λ, λ_mu, σ_bg, z = z[1:nv])
end

"""
Shared pooling SD `σ_bg` for the per-vintage background random effect
([`background_re_model`](@ref)). Sampled once at the composer level and
passed to both the suspected-case and suspected-death backgrounds, so the
two streams share one time-variation scale rather than each estimating its
own from few vintages. The prior is a tight half-normal
`truncated(Normal(0, 0.3); lower = 0)`, deliberately small: the background
is degenerate with outbreak size, so a wide random effect would let
individual windows absorb arbitrary suspected counts and re-open the
second posterior mode in which the background explains the majority of
suspected cases. Regularising `σ_bg` toward zero keeps the time variation
a perturbation of the informative scalar baselines. Returns `(; σ_bg)`.
"""
@model function background_pooling_model(;
        pooling_prior = truncated(Normal(0.0, 0.3); lower = 0))
    σ_bg ~ pooling_prior
    return (; σ_bg)
end

"""
PCR sensitivity prior for the GeneXpert Ebola assay. `Beta(30, 2)` has
mean 0.94 and 95% interval 0.84–0.99, sitting just below the field whole
blood clinical sensitivity reported in the Sierra Leone Zaire-ebolavirus
field evaluation, leaving room for early-infection low-viral-load
specimens and field handling. Scales the confirmed-case stream so the
confirmed counts reflect imperfect detection of true BVD infections.
Matches the integral `main` prior. Returns `(; s_test)`.
"""
@model function test_sensitivity_model(; sensitivity_prior = Beta(30.0, 2.0))
    s_test ~ sensitivity_prior
    return (; s_test)
end

"""
PCR specificity prior for the Ebola assay. `Beta(60, 2)` has mean ≈ 0.97
and 95% interval ≈ 0.91–0.998, a high-but-imperfect specificity reflecting
that a small fraction of non-BVD specimens test positive (cross-reaction,
contamination, low-level false positives). Used by the composition-linked
confirmed-case positivity so the tested-positive probability is
`p = s · q + (1 − spec)(1 − q)` with `q` the tested BVD share: the
false-positive term `(1 − spec)(1 − q)` makes the confirmed counts respond
to the non-BVD share `1 − q`, so the laboratory data identify the
background `λ_bg` rather than only the BVD signal. Returns `(; spec)`.
"""
@model function test_specificity_model(; specificity_prior = Beta(60.0, 2.0))
    spec ~ specificity_prior
    return (; spec)
end

"""
Report-to-laboratory-confirmation (lab-turnaround) delay submodel. The
delay from a suspected case being reported to its specimen being
laboratory confirmed, discretised to a daily PMF over lags `0 … nmax`
by [`censored_delay_model`](@ref) so it convolves cleanly onto the
renewal onsets. The mean and SD carry weakly-informative priors centred
on a short turnaround with a heavy right tail allowing for specimen
shipment to a confirmatory lab; no per-sample outbreak data anchors this
prior, matching integral `main`. Returns `(; pmf, dist, mean, sd)`.
"""
@model function lab_delay_model(nmax::Integer = 30;
        mean_prior = truncated(Normal(4.5, 2.0); lower = 1),
        sd_prior = truncated(Normal(4.0, 1.5); lower = 1))
    d ~ to_submodel(censored_delay_model(nmax; mean_prior, sd_prior))
    return (; pmf = d.pmf, dist = d.dist, mean = d.mean, sd = d.sd)
end

"""
Per-vintage laboratory positivity for the confirmed-case stream, in
partially-pooled non-centred form. The confirmed positives in each
laboratory window are scored as a `Binomial` of the observed
specimens-analysed denominator (see [`confirmed_cases_model`](@ref)), so
the positivity is the probability a tested specimen is confirmed. The
per-window logit positivity shares a baseline `q_mu` and is perturbed by
non-centred deviations `z_q` scaled by the pooling SD `σ_q`, so a window
with little data is shrunk toward the baseline while a window with a
strong signal can depart from it. `σ_q → 0` recovers a single shared
positivity. The baseline prior is centred on the cut-off cumulative
positivity (≈ 0.28, i.e. 210 / 755 on the 28 May data) on the logit
scale. Conditioning on the observed denominator and giving the positivity
its own random effect decouples the confirmed counts from the
multiplicative ascertainment ridge (`p_drc · s_test · τ_test`) that
basin-split the joint, so the outbreak size is pinned by the deaths and
exports streams rather than forced through the laboratory positivity.
Returns `(; p_pos, q_mu, σ_q)` with `p_pos` a length-`nv` vector.
"""
@model function confirmed_positivity_model(nv::Integer;
        baseline_prior = Normal(logit(0.28), 0.7),
        pooling_prior = truncated(Normal(0.0, 1.0); lower = 0))
    m = max(nv, 1)
    q_mu ~ baseline_prior
    σ_q ~ pooling_prior
    z_q ~ product_distribution(fill(Normal(0, 1), m))
    logit_p = q_mu .+ σ_q .* z_q[1:nv]
    p_pos := logistic.(logit_p)
    return (; p_pos, q_mu, σ_q)
end

"""
Severity-enrichment prior for the COMPOSITION-LINKED confirmed positivity
(`positivity_link = :composition` in [`confirmed_cases_model`](@ref)). In
that mode the per-window tested BVD share is not a free random effect; it
is the suspect-pool composition `φ_v = (p_drc · BVD)_v / ((p_drc · BVD)_v +
λ_bg_v)` over each laboratory window, UPSAMPLED by a severity enrichment
that decays as testing widens:

```math
\\mathrm{logit}(q_v) = \\mathrm{logit}(\\varphi_v)
    + \\delta_0\\, e^{-c_v / \\text{decay}},
```

with `c_v` the cumulative analysed volume at window `v` (the testing
clock). The lab over-tests BVD early (severe cases are triaged first and
are more likely BVD), the enrichment `δ₀·e^{−c/decay}` relaxing toward zero
as testing widens, at which point the tested share equals the pool
composition. This ties positivity to the background `λ_bg`, so the
confirmed/positivity data identify the non-BVD background rather than it
being absorbed by a free per-window random effect, correcting the model's
treatment of suspected cases as a large overestimate of true BVD.

`δ₀` is the early severity log-odds enrichment of BVD; lower-truncated at 0
because severity triage upsamples BVD, never down. The default
`truncated(Normal(1.5, 0.75); lower = 0)` is deliberately moderate /
bounded: even severity-triaged testing cannot be near-pure BVD (other
haemorrhagic / severe febrile illness is also triaged), so for a pool
composition `φ ≈ 0.4` the early tested share is `logistic(logit(0.4) +
1.5) ≈ 0.75`. `decay_scale` is the relaxation timescale on the analysed-
volume clock. Pass `logodds_prior` / `decay_prior` to override. Used by
[`confirmed_cases_model`](@ref) in composition mode. Returns
`(; δ0, decay_scale)`.
"""
@model function severity_enrichment_model(;
        logodds_prior = truncated(Normal(1.5, 0.75); lower = 0),
        decay_prior = truncated(Normal(0.0, 200.0); lower = 0.0))
    δ0 ~ logodds_prior
    decay_scale ~ decay_prior
    return (; δ0, decay_scale)
end

"""
Confirmed-death enrichment scalar `m_death` for the confirmed-death
stream ([`confirmed_deaths_model`](@ref)). Confirmed deaths are a thinning
of the suspected deaths whose per-window confirmation probability is the
suspected-case BVD composition `q_susp` enriched on the odds scale by
`m_death`:

```math
p = \\operatorname{logistic}(\\operatorname{logit}(q_\\text{susp}) +
    \\log m_\\text{death}),
```

a multiplicative effect on the confirmation odds that stays in `(0, 1)`
without a hard clamp. `m_death > 1` means a suspected death is more likely
to be laboratory-confirmed BVD than a suspected case; `m_death < 1` means
it is less likely, as when deaths are under-swabbed relative to living
suspected cases (post-mortem confirmation is rarer). `m_death = 1` ties
the death-confirmation rate to the case composition. The default
`LogNormal(0, 1.0)` is weakly informative and centred on no enrichment, so
the single informative confirmed-death total (17 of 246 suspected deaths
on the 28 May data) determines the differential rather than the prior:
with a suspected-case BVD composition near 0.4, that count needs
`m_death ≈ 0.1`, which a tight prior centred on 1 would fight. Returns
`(; m_death)`.
"""
@model function confirmed_death_enrichment_model(;
        enrichment_prior = LogNormal(0.0, 1.0))
    m_death ~ enrichment_prior
    return (; m_death)
end

"""
Shared negative-binomial dispersion `k` for the passive-surveillance
streams (suspected deaths, reported cases and confirmed cases). Sampled on
the `1/sqrt(k)` scale with a weakly-informative half-normal prior
following the Stan prior-choice recommendations. Returns
`(; k, inv_sqrt_k)`.
"""
@model function surveillance_dispersion_model(;
        inv_sqrt_k_prior = truncated(Normal(0.6, 0.2); lower = 0))
    inv_sqrt_k ~ inv_sqrt_k_prior
    k := 1.0 / (inv_sqrt_k^2 + eps(typeof(inv_sqrt_k)))
    return (; k, inv_sqrt_k)
end

"""
Independent ascertainment fractions for the DRC and Uganda surveillance
systems. The two countries run different surveillance systems — DRC
passive community surveillance and Uganda point-of-entry / hospital
detection — so each ascertainment fraction has its own logit-scale prior
with no shared parameter. An alternative to the composer-default
[`pooled_ascertainment_model`](@ref) for sensitivity analyses that do
not share strength between the two systems.

Both fractions default to a logit-Normal prior centred on a reporting
fraction of 0.75 with SD 0.6 (95% support roughly 0.48–0.91), reflecting
the active case-finding and contact tracing of a declared Ebola response
rather than baseline passive surveillance. Pass `drc_prior` /
`uganda_prior` to set the two systems' priors separately.
"""
@model function independent_ascertainment_model(;
        drc_prior = Normal(logit(0.75), 0.6),
        uganda_prior = Normal(logit(0.75), 0.6))
    logit_p_drc ~ drc_prior
    logit_p_uganda ~ uganda_prior
    p_drc := logistic(logit_p_drc)
    p_uganda := logistic(logit_p_uganda)
    return (; logit_p_drc, logit_p_uganda, p_drc, p_uganda)
end

"""
Partially pooled ascertainment fractions for the DRC and Uganda
surveillance systems, sampled in non-centred form to avoid the funnel
geometry. Both logit-scale fractions share a hyperprior with mean `μ`
and pooling strength `τ`. Used by [`reported_cases_model`](@ref),
[`exports_model`](@ref) and [`exports_deaths_model`](@ref); this is the
composer default. The shared mean defaults to a reporting fraction of
0.75 (logit scale), reflecting the active case-finding of a declared
Ebola response; a lower ascertainment inflates the inferred infections
(and so the outbreak size `C_T`) for the same observed counts.
"""
@model function pooled_ascertainment_model(;
        mu_prior = Normal(logit(0.75), 1.0),
        tau_prior = truncated(Normal(0, 0.5); lower = 1e-4))
    μ_logit ~ mu_prior
    τ_logit ~ tau_prior
    z_drc ~ Normal(0, 1)
    z_uganda ~ Normal(0, 1)
    logit_p_drc = μ_logit + τ_logit * z_drc
    logit_p_uganda = μ_logit + τ_logit * z_uganda
    p_drc := logistic(logit_p_drc)
    p_uganda := logistic(logit_p_uganda)
    return (; μ_logit, τ_logit, p_drc, p_uganda)
end
