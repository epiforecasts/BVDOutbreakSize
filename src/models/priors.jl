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
ramp scale `ramp` defaults to 14 days, a fortnight transition reflecting
that a response (case finding, isolation, vaccination) takes weeks to bite
rather than switching at a single date; pass `ramp` to widen or narrow it.
`Rt = exp.(log_Rt)`. Returns
`(; Rt, log_R, days, sigma_rw, log_R0, intervention_effect)`.
"""
@model function rt_walk_model(n::Integer;
        week::Integer = 7,
        breakpoint::Union{Missing, Real} = missing,
        ramp::Real = 14.0,
        log_r0_prior = Normal(log(1.3), 0.4),
        sigma_prior = truncated(Normal(0, 0.2); lower = 0),
        effect_prior = Normal(0, 0.5))
    days = knot_days(n; week)
    nb = length(days)
    log_R0 ~ log_r0_prior
    sigma_rw ~ sigma_prior
    z ~ product_distribution(fill(Normal(0, 1), max(nb - 1, 1)))
    intervention_effect ~ effect_prior
    steps = sigma_rw .* z[1:(nb - 1)]
    log_R = log_R0 .+ vcat(zero(log_R0), cumsum(steps))
    log_Rt = interpolate_knots(log_R, days, n)
    log_Rt = log_Rt .+ intervention_effect .* sigmoid_ramp(n, breakpoint; ramp)
    Rt = exp.(log_Rt)
    return (; Rt, log_R, days, sigma_rw, log_R0, intervention_effect)
end

## --- Seeding and the generating infection process -----------------------

"""
Seed submodel: the latent infection count `I0` on the last day of the
seeding window, representing the zoonotic introduction. Default prior is a
truncated Normal centred on a single seed; the prior is injectable. The
seeding window is filled by exponential growth at the implied rate in
[`infection_model`](@ref).
"""
@model function seed_model(; i0_prior = truncated(Normal(1.0, 1.0); lower = 0))
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
        rt = rt_walk_model,
        gi = generation_interval_model,
        seed = seed_model,
        gi_nmax::Integer = 40)
    rt_state ~ to_submodel(rt(n; breakpoint))
    gi_state ~ to_submodel(gi(gi_nmax))
    seed_state ~ to_submodel(seed())
    Rt = rt_state.Rt
    g = gi_state.g
    r0 = euler_lotka_r(Rt[1], g)
    seed_vec = seed_infections(seed_state.I0, r0, length(g))
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
incubation delay submodel is injected, defaulting to a prior centred on
the Ebola virus disease incubation period (mean 9.7 d, SD 5.4 d; WHO Ebola
Response Team 2014, NEJM). Returns
`(; onsets, incubation_pmf, incubation_mean, incubation_sd)`.
"""
@model function onset_incidence_model(infections::AbstractVector;
        incubation = (nmax) -> censored_delay_model(nmax;
            mean_prior = truncated(Normal(9.7, 2.0); lower = 1),
            sd_prior = truncated(Normal(5.4, 1.5); lower = 1)),
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
Test-positivity machinery shared by the suspected- and confirmed-case
streams. Samples

- `λ_bg` — the per-day non-BVD background suspected-case rate, on a
  half-normal scale. Drives the suspected/confirmed contrast: suspected
  cases mix the BVD onset-to-report signal with this additive background,
  while the laboratory pipeline only confirms the BVD share.
- `τ_test` — the fraction of suspected cases that are sampled and routed
  to the laboratory pipeline.

The derived per-suspected positivity is exposed inside
[`reported_cases_model`](@ref); the per-test positivity is exposed inside
[`confirmed_cases_model`](@ref). The prior centres follow integral
`main`: a wide half-normal on the background rate and a `Beta(5, 2)`
testing fraction (mean ≈ 0.71). Returns `(; λ_bg, τ_test)`.
"""
@model function test_positivity_model(;
        lambda_prior = truncated(Normal(0.0, 10.0); lower = 0),
        fraction_tested_prior = Beta(5.0, 2.0))
    λ_bg ~ lambda_prior
    τ_test ~ fraction_tested_prior
    return (; λ_bg, τ_test)
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
with no shared parameter. Used by [`reported_cases_model`](@ref),
[`exports_model`](@ref) and [`exports_deaths_model`](@ref); this is the
composer default.

Both fractions default to a logit-Normal prior centred on a reporting
fraction of 0.25 with SD 0.6 (95% support roughly 0.09–0.51), weakly
informative about the low-ascertainment regime typical of passive BVD /
Ebola surveillance in rural Ituri. Pass `drc_prior` / `uganda_prior` to
set the two systems' priors separately.
"""
@model function independent_ascertainment_model(;
        drc_prior = Normal(logit(0.25), 0.6),
        uganda_prior = Normal(logit(0.25), 0.6))
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
and pooling strength `τ`. An alternative to the composer-default
[`independent_ascertainment_model`](@ref) for sensitivity analyses that
share strength between the two systems.
"""
@model function pooled_ascertainment_model(;
        mu_prior = Normal(logit(0.25), 1.0),
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
