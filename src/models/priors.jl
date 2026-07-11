# Prior and latent-state submodels: the building blocks shared across the
# observation submodels and the joint composer. Each `@model` is one piece
# of the generative process — a single prior, a delay, the reproduction
# number, the generating infection process, or the onset staging — so it
# can be reused across composers without duplication. Delays are sampled
# from priors and discretised with CensoredDistributions. Nothing is fixed.

## --- Delay submodels (priors only, all delays sampled) -------------------

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
Generation-interval submodel. Samples a Gamma shape `α` and scale `θ`,
parameterised directly from the source's reported generation-time
distribution rather than moment-matching a LogNormal from mean/SD priors.
The source is the Ebola virus disease serial interval as a generation-time
proxy (mean 15.3 d, SD 9.3 d; WHO Ebola Response Team 2014, NEJM), which
maps once to a Gamma shape `α ≈ 2.71` and scale `θ ≈ 5.65`
(`α = (mean/sd)²`, `θ = sd²/mean`). The priors are centred on those values:
`α ~ Normal⁺(2.71, 0.7)` and `θ ~ Normal⁺(5.65, 1.5)`, lower-truncated to
keep the Gamma well defined. The SDs propagate the source's reported
uncertainty (the NEJM serial-interval mean carries a 95% CI of 13.0–17.6 d,
an SD on the mean of ≈1.17 d), not a self-assigned spread. Gamma
shape/scale differentiate cleanly under Mooncake, so this is AD-stable.

Discretised through the same double-interval-censoring route as the other
delays ([`discretise_censored`](@ref)). The lag-0 bin is dropped and the
remainder renormalised, left-truncating the generation interval at one day
so an infectee is infected strictly after its infector. Returns
`(; g, gi_mean, gi_sd, gi_alpha, gi_theta)`.
"""
@model function generation_interval_model(nmax::Integer;
        alpha_prior = truncated(Normal(2.71, 0.7); lower = 0.1),
        theta_prior = truncated(Normal(5.65, 1.5); lower = 0.1))
    α ~ alpha_prior
    θ ~ theta_prior
    dist = Gamma(α, θ)
    pmf = discretise_censored(dist, nmax)
    g = pmf[2:end] ./ sum(pmf[2:end])
    return (; g, gi_mean = α * θ, gi_sd = sqrt(α) * θ,
        gi_alpha = α, gi_theta = θ)
end

"""
Natural-parameter Gamma delay submodel. Samples a Gamma shape `α` and
scale `θ` directly from the priors, builds `Gamma(α, θ)`, and discretises
to a daily PMF over lags `0 … nmax` by double interval censoring
([`discretise_censored`](@ref)), keeping the lag-0 bin (an onset-to-event
delay can be same-day, unlike the generation interval). This carries a
line-list delay reanalysis through on its natural parameters with the
reported posterior uncertainty, rather than moment-matching a LogNormal
from mean/SD priors. Gamma shape/scale differentiate cleanly under
Mooncake. Returns `(; pmf, dist, mean, sd, alpha, theta)`.
"""
@model function gamma_delay_model(nmax::Integer; alpha_prior, theta_prior)
    α ~ alpha_prior
    θ ~ theta_prior
    dist = Gamma(α, θ)
    return (; pmf = discretise_censored(dist, nmax), dist,
        mean = α * θ, sd = sqrt(α) * θ, alpha = α, theta = θ)
end

"""
Onset-to-death delay as the convolution of two natural-parameter Gamma
atomic delays, onset→admission (`oa`) and admission→death (`ad`), each
sampled through [`gamma_delay_model`](@ref) and combined by convolving
their PMFs. This matches the companion line-list reanalysis, which fits
the atomic components and convolves them rather than fitting onset→death
directly, so no moment-matching is needed: each atomic delay keeps its
own Gamma shape and scale prior with the reanalysis's reported
uncertainty. The convolved PMF is truncated back to lags `0 … nmax` and
renormalised. Returns `(; pmf, mean, sd, oa_mean, ad_mean)`.
"""
@model function onset_to_death_model(nmax::Integer;
        oa_alpha_prior, oa_theta_prior, ad_alpha_prior, ad_theta_prior)
    oa ~ to_submodel(gamma_delay_model(nmax; alpha_prior = oa_alpha_prior,
        theta_prior = oa_theta_prior))
    ad ~ to_submodel(gamma_delay_model(nmax; alpha_prior = ad_alpha_prior,
        theta_prior = ad_theta_prior))
    full = convolve_pmf(oa.pmf, ad.pmf)
    trimmed = full[1:(nmax + 1)]
    pmf = trimmed ./ sum(trimmed)
    return (; pmf, mean = oa.mean + ad.mean,
        sd = sqrt(oa.sd^2 + ad.sd^2), oa_mean = oa.mean, ad_mean = ad.mean)
end

"""
Wilson–Hilferty approximation to the continuous median of a `Gamma` as a
function of its `mean` and `sd`: `median ≈ mean·(1 − sd²/(9·mean²))³`.
Smooth in the mean and SD with no quantile inversion, so it enters a
gradient-based likelihood cleanly, and accurate to a few percent for the
shapes here. Used to match the confirmed onset-to-sample convolution's
continuous median to the cohort's reported median.
"""
gamma_median_wh(mean::Real, sd::Real) = mean * (1 - sd^2 / (9 * mean^2))^3

"""
Onset-to-sample prior configuration from the NEJM DRC 2026 BVD cohort
(Akilimali et al. 2026, doi:10.1056/NEJMc2608070). The confirmed-positive
onset-to-sample interval (N = 129) was estimated as a continuous Gamma through
the `epidist` marginal model correcting for double interval censoring and right
truncation, chosen over lognormal and Weibull by LOOIC. The cohort reports a
continuous mean of 7.4 d (95% CrI 5.3–13.5) and median of 4.8 d (95% CrI
3.46–7.84).

The confirmed onset→report→receipt convolution is grounded on its mean and
median. The mean is the sum of the two legs' means and the variance the
sum of their variances. The median follows by [`gamma_median_wh`](@ref).
Each is fitted to the reported value as a Normal observation whose SD is
the reported 95% CrI half-width over 1.96 (`mean_se`, `median_se`), so
the cohort's uncertainty enters directly and the constraint is soft.
Returns a NamedTuple
`(; mean_obs, mean_se, median_obs, median_se)` for
the `onset_to_sample` argument of [`bvd_joint`](@ref).
"""
function nejm_onset_to_sample(; mean::Real = 7.4,
        mean_se::Real = (13.5 - 5.3) / 2 / 1.96, median::Real = 4.8,
        median_se::Real = (7.84 - 3.46) / 2 / 1.96)
    return (; mean_obs = mean, mean_se, median_obs = median, median_se)
end

"""
Log-density grounding the confirmed onset-to-sample convolution on the cohort:
soft Normal fits of the convolution's continuous mean (the sum of the report
and receipt leg means) and continuous median (from [`gamma_median_wh`](@ref) of
the summed leg variances) to the reported `mean_obs`/`median_obs` with SDs
`mean_se`/`median_se`. The fixed Gaussian normalising constants are dropped.
"""
function onset_to_sample_logweight(report_mean::Real, report_sd::Real,
        receipt_mean::Real, receipt_sd::Real, cfg)
    μ = report_mean + receipt_mean
    sd = sqrt(report_sd^2 + receipt_sd^2)
    med = gamma_median_wh(μ, sd)
    return -0.5 * ((cfg.mean_obs - μ) / cfg.mean_se)^2 -
           0.5 * ((cfg.median_obs - med) / cfg.median_se)^2
end

## --- Reproduction number ------------------------------------------------

"""
Weekly piecewise-linear log-scale reproduction number over `n` days, with
a smooth intervention ramp. Knots sit at weekly spacing
([`knot_days`](@ref)) and follow a Gaussian random walk in non-centred
cumulative-sum form: standard-normal innovations are scaled by `sigma_rw`
and accumulated, avoiding the funnel geometry of the centred recursion.
Daily log-`R_t` is the linear interpolation between knots
([`interpolate_knots`](@ref)). An intervention at `breakpoint` (e.g. the
first WHO situation report) adds a sampled effect `intervention_effect`
shaped by a logistic ramp ([`sigmoid_ramp`](@ref)) of scale `ramp`
(default 21 days, roughly the time a response takes to bite), so
transmission changes gradually rather than instantly. `breakpoint =
missing` drops the term. `Rt = exp.(log_Rt)`.

The walk base `log_R0` is not sampled here. It is passed in as a derived
quantity: the first reproduction number is derived forward from the sampled
growth rate `r` and the generation interval through Euler–Lotka
(`R0 = r_to_R0(r, g)` in [`infection_model`](@ref)). The prior therefore
sits on the growth rate (see [`exponential_growth_model`](@ref)), grounded
on the BEAST X molecular-clock doubling time (mbalaplacide2026). The
established reproduction number is whatever that growth implies under our
generation interval, rather than a separately asserted `R0` prior. The
genetic report gives `R0 ≈ 1.31–1.55` under its own generation interval.
Deriving `R0` forward from the shared growth rate under our generation
interval is the consistent thing to do. This single growth source pins both
the established reproduction number (the walk base at the genetic bound)
and, through the renewal seeding, the cryptic exponential phase, so the
outbreak has one growth source. The grid days before the renewal start are
filled by the analytic cryptic exponential and so are unused by the walk.
The walk clamps to `R0` before its first knot.

The random-walk step SD prior is a half-normal SD 0.1, so the weekly
log-`R_t` is unlikely to change by more than about 20% (two SD ≈ 0.2) from
one week to the next. The walk starts at the first situation report
(`rt_start`), so every knot sits in the observed window and `R_t` bends to
a slowdown or acceleration over the sitreps rather than drifting over the
unobserved pre-report stretch.

The intervention effect is constrained non-positive
(`truncated(Normal(0, 0.4); upper = 0)`): a declared WHO response (case
finding, isolation, vaccination) can only reduce transmission or leave it
unchanged, so the ramp's effect on log-`R_t` is bounded at or below zero.
The half-normal admits anything from no effect (mode) to a substantial
decline. The breakpoint sits only around 11 days before the cut-off and
the ramp runs three weeks, so the response damps `R_t` only partially by
the cut-off.

Returns `(; Rt, log_R, days, sigma_rw, log_R0, intervention_effect)`.
"""
@model function rt_walk_model(n::Integer, log_R0_base::Real;
        week::Integer = 7,
        breakpoint::Union{Missing, Real} = missing,
        rt_start::Integer = 1,
        ramp::Real = RT_INTERVENTION_RAMP,
        sigma_prior = truncated(Normal(0, 0.1); lower = 0),
        effect_prior = truncated(Normal(0, 0.4); upper = 0))
    days = knot_days(n; week, start = rt_start)
    nb = length(days)
    ## The established `R0` at the genetic bound is the base the random walk
    ## grows from. It is derived (forward Euler–Lotka from the sampled growth
    ## rate) and passed in, not sampled here. It is tracked as a deterministic
    ## so the walk base stays available on the chain. The days before the
    ## renewal start (`rt_start`) are filled by the analytic cryptic
    ## exponential in `infection_model`, so the walk values there are unused.
    ## The interpolation clamps to `log_R0` before the first knot, which is
    ## harmless.
    log_R0 := log_R0_base
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
Molecular-clock growth-and-size prior for the renewal cryptic phase.
Samples the exponential growth rate `r` directly (the primary
epidemiological assumption, placed on the genetic doubling time) and the
doubling count `m`, then exposes

```math
\\tau = \\log 2 / r,\\qquad T_\\text{cryptic} = m\\,\\tau,
\\qquad C_T = 2^m,
```

as deterministics. `T = m·τ` is the cryptic-phase duration (origin →
renewal start): `m` counts the doublings during the cryptic phase, so the
cryptic phase grows a single import to `2^m` infections at the renewal
start, independent of `r`. The composer ([`infection_model`](@ref)) adds
the observation span `τ_obs = n − renewal_start` to get the total outbreak
age `T_total = m·τ + τ_obs`, which carries the genetic seeding bound
([`genetic_seeding_model`](@ref)).

The growth rate carries the prior `r ~ LogNormal(log(log2 /
M_PRIOR_DOUBLING_DAYS), 0.40)`, with median doubling time (11.7 d)
matching the BEAST X estimate (mbalaplacide2026, Exponential growth
model, 95% HPD 6.8–17.5). The log-SD 0.40 is wider than the ≈0.24 that
HPD implies. The HPD is conditional on a single-rate coalescent, the
assumption the field epidemiology contradicts (kupferschmidt2026). An
independent reanalysis puts the doubling time at 15.2–24.5 d
(cuomodannenburg2026). The induced doubling-time prior is
`LogNormal(log 11.7, 0.40)`, 5.3–25.6 d at 95%. The first
reproduction number is then derived forward from this `r` and our generation
interval through Euler–Lotka (`R0 = r_to_R0(r, g)` in
[`infection_model`](@ref)), so the cryptic exponential phase and the
established renewal share one growth source, and the established
reproduction number is consistent with the genetic growth under our
generation interval rather than pinned by a separate `R0` prior.

`m ~ truncated(Normal(5, 4); lower = 0)` is deliberately wide. The centre
and spread are elicited from the field epidemiology and the genetics:
field work in Mongbwalu traced a sustained transmission chain back to a
death on 25 January 2026 and identified 500+ suspected cases between
mid-January and mid-May (kupferschmidt2026), and the genetic TMRCA
(mbalaplacide2026) is a lower bound on the outbreak age consistent with an
origin that early. `m` counts only the cryptic doublings, not the cut-off
case total. With the ≈11.7-day doubling, a centre of 5 spans
`5 · τ ≈ 58.5` cryptic days, placing the implied prior origin at the end
of January (the documented first deaths). The SD of 4 lets the data and
the genetic seeding bound pull the origin earlier (into December) or
later without over-committing.

In the renewal, `2^m` is the prior seed at the renewal start, which the
renewal recursion grows forward under `R_t`. Pass `m_prior` to override,
e.g. one whose centre advances via [`m_prior_centre`](@ref) for a later
cut-off. Returns `(; τ, r, m, T, C_T)`.
"""
@model function exponential_growth_model(;
        r_prior = LogNormal(log(log(2) / M_PRIOR_DOUBLING_DAYS), 0.40),
        m_prior = truncated(Normal(5.0, 4.0); lower = 0))
    r ~ r_prior
    m ~ m_prior
    τ := log(2) / r
    T := m * τ
    C_T := 2.0^m
    return (; τ, r, m, T, C_T)
end

"""
Seed submodel: the latent infection count `I0` on the last day of the
seeding window, representing the zoonotic introduction. Default prior is a
truncated Normal centred on a single seed. The prior is injectable. The
seeding window is filled by exponential growth at the implied rate in
[`infection_model`](@ref).
"""
@model function seed_model(; i0_prior = truncated(Normal(0.1, 0.1); lower = 0))
    I0 ~ i0_prior
    return (; I0)
end

"""
Generating infection process for the two-phase renewal seeding. Samples
the generation interval and the cryptic exponential growth rate `r` (the
prior sits on `r`, the molecular-clock growth, in
[`exponential_growth_model`](@ref)), then derives the established
reproduction number `R0` (= the first `R_t`) forward from that `r` and the
generation interval through Euler–Lotka (`R0 = r_to_R0(r, g)`) and passes
`log R0` as the walk base to the reproduction-number submodel. The
cryptic phase and the established renewal share one growth source, the
sampled growth rate `r`, rather than the cryptic phase carrying a
separate clock-rate prior or the walk asserting a separate `R0` prior.

The renewal runs only over the observation window
`[renewal_start, cut-off]`, where `renewal_start` is the genetic-TMRCA
grid day `rt_start` (the day the reproduction-number walk starts, before
which `R_t` is held flat). The cryptic exponential phase from the origin to
the renewal start is analytic and off the renewal grid except for the
days needed as recursion history. The seed at the renewal start is the
cryptic-phase realised size `2^m` ([`seed_at_renewal_start`](@ref)), where
`m` counts the doublings during the cryptic phase: the magnitude is
independent of `r`, so `r` (hence the derived `R0`) leaves the seed
magnitude alone and appears only in the renewal growth. A back-scaled
seed `2^m e^{−r·τ_obs}` would put `r` into both the seed and the renewal
growth, cancelling for a fixed realised size and opening a flat ridge
along which `R0` slides to 1. Keeping the seed `r`-independent avoids
this. The grid days `1 … renewal_start` before the renewal start are
filled smoothly by the cryptic exponential curve at rate `r` ending at
`2^m` ([`seed_infections`](@ref)), giving the recursion a full generation
interval of differentiable history. The renewal recursion
([`renewal_infections`](@ref)) then grows the trajectory over
`renewal_start+1 … n` under the time-varying `R_t`.

The total outbreak age is `T = m·τ + τ_obs` (cryptic duration plus the
observation span `τ_obs = n − renewal_start`). The genetic seeding bound
is applied to this total `T` at the composer. The renewal start sits a
small lead after the genetic TMRCA day, past the TMRCA uncertainty where
sustained transmission is confident, so `τ_obs = n − renewal_start <
tmrca_days` and the censored bound
`tmrca ~ censored(Normal(T, sd); upper = tmrca_days)` stays informative:
it pulls the origin to sit at or before the MRCA, so the cryptic duration
`m·τ` cannot be too short. The genetic bound therefore defines the
cryptic-phase length through `T`.

The realised cut-off size is `C_T = cumulative[n]`. The `breakpoint` is
forwarded to the reproduction-number submodel. Exposes the daily
infections and cumulative sum, the total prior outbreak age `T`, cryptic
doubling count `m`/`τ` and prior size scale `C_T_prior`, the realised
cut-off size `C_T`, the established reproduction number `R0` and its
implied cryptic rate `r0` (with `doubling_time_initial`), the current
growth `r`/`doubling_time` derived from the cut-off reproduction number
`Rt[n]` through forward Euler–Lotka (so `r` is sign-consistent with
`R_T := Rt[n]` by construction), and the diagnostic-only `seeding_age`.
Returns `(; infections, cumulative, Rt, g, seed_at_renewal_start, m, τ,
R0, r0, r, doubling_time_initial, T, C_T, C_T_prior, doubling_time,
seeding_age)`.
"""
@model function infection_model(n::Integer;
        breakpoint::Union{Missing, Real} = missing,
        rt_start::Integer = 1,
        rt_walk_start::Integer = rt_start,
        rt = rt_walk_model,
        gi = generation_interval_model,
        growth = exponential_growth_model,
        gi_nmax::Integer = cdf_nmax(Gamma(2.71, 5.65)))
    gi_state ~ to_submodel(gi(gi_nmax))
    g = gi_state.g
    ## One growth source: the prior is on the cryptic exponential growth rate
    ## `r` (sampled in `growth`), and the established reproduction number
    ## `R0` (= the walk base, the first `R_t`) is derived forward from that
    ## `r` and the generation interval through Euler–Lotka. The cryptic
    ## phase and the established renewal therefore share `r`.
    growth_state ~ to_submodel(growth())
    r_clock = growth_state.r
    R0 = r_to_R0(r_clock, g)
    ## The random walk's first knot sits at `rt_walk_start`, decoupled from
    ## the renewal start `rt_start`: the renewal seeds and grows from the
    ## genetic-TMRCA renewal start, but `R_t` is held flat at `R0` until
    ## `rt_walk_start` (the first situation report). Before any case or death
    ## surveillance the dynamics are unidentified, so a free walk there only
    ## adds unsupported drift. `rt_walk_start` defaults to `rt_start`, the
    ## walk-from-renewal-start case.
    rt_state ~ to_submodel(rt(n, log(R0); breakpoint, rt_start = rt_walk_start))
    Rt = rt_state.Rt
    ## renewal_start = genetic-TMRCA grid day (`rt_start`). The observation
    ## span is τ_obs = n − renewal_start. The renewal-start seed magnitude is
    ## `2^m` directly (the cryptic phase grows one import to `2^m` over `m`
    ## doublings, `r`-free). Fill grid days 1…renewal_start with the cryptic
    ## exponential curve at rate `r` ending at `2^m` (a full GI of history),
    ## then run the renewal forward from renewal_start+1.
    renewal_start = clamp(rt_start, 1, n)
    τ_obs = n - renewal_start
    seed0 = seed_at_renewal_start(growth_state.C_T)
    seed_vec = seed_infections(seed0, r_clock, renewal_start)
    infections = renewal_infections(Rt, g, seed_vec)
    cumulative = cumsum(infections)
    ## Total outbreak age: cryptic duration (m·τ) plus the observation span.
    T_total = growth_state.T + τ_obs
    ## Current growth rate at the cut-off, derived from the cut-off
    ## reproduction number `Rt[n]` and the generation interval through forward
    ## Euler–Lotka (the inverse of the `r_to_R0` that derives `R0` from the
    ## clock growth above). This makes the reported current growth rate
    ## consistent with `R_T := Rt[n]` by construction: `r < 0` iff `R_T < 1`.
    ## The realised last-two-days slope `log I[n] − log I[n-1]` is not used
    ## for this: the intervention ramp depresses the final renewal step
    ## (`I[n] < I[n-1]` while `Rt[n] ≥ 1`), so that slope can disagree in
    ## sign with `R_T` at the cut-off, an end-of-trajectory edge artifact
    ## rather than the instantaneous growth.
    r = euler_lotka_r(@inbounds(Rt[n]), g)
    return (; infections, cumulative, Rt, g, seed_at_renewal_start = seed0,
        m = growth_state.m, τ = growth_state.τ, R0, r0 = r_clock, r,
        doubling_time_initial = doubling_time(r_clock),
        T = T_total, C_T = cumulative[n],
        C_T_prior = growth_state.C_T, doubling_time = doubling_time(r),
        seeding_age = seeding_age(cumulative, n))
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
half-width / 1.96). MacNeil give no interval on the spread, so the SD
prior is a weakly-informative modelling choice centred on the
CV-implied spread (≈ 3.5 d). Returns
`(; onsets, incubation_pmf, incubation_mean, incubation_sd)`.
"""
@model function onset_incidence_model(infections::AbstractVector;
        incubation = (nmax) -> censored_delay_model(nmax;
            mean_prior = truncated(Normal(6.3, 0.54); lower = 1),
            sd_prior = truncated(Normal(3.5, 0.8); lower = 1)),
        incubation_nmax::Integer = cdf_nmax(lognormal_meansd(6.3, 3.5)))
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
it older. The likelihood contributes `P(read ≥ tmrca_days)`. Passing
`tmrca_days = missing` makes the submodel a no-op.
"""
@model function genetic_seeding_model(T::Real,
        tmrca_days::Union{Missing, Real}; tmrca_days_sd::Real = 16.0)
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
not be a true BVD death. This submodel samples the per-day non-BVD
background death rate `λ_bg_death`. Its cumulative contribution over the
grid is `λ_bg_death · n`.

The default `truncated(Normal(0, 0.25); lower = 0)` is deliberately
informative, mirroring the case background: deaths are far fewer than
suspected cases (≈ 246 suspected deaths vs ≈ 1077 suspected cases at the
last stable suspected vintages), so the background rate is scaled down
accordingly. With SD 0.25 the median background is ≈ 0.17/day (a modest
minority of the suspected-death total over the grid) while still
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
Ascertainment of the suspected-death stream ([`deaths_model`](@ref)): the
fraction `p_death` of true BVD deaths that enter the INSP suspected-death
count by the cut-off. The suspected-death definition is symptomatic-then-
deceased, so a fatal BVD infection only counts once the death is reported,
and not every BVD death is captured. The expected BVD suspected deaths are
therefore `p_death · CFR` of the onset-to-death-convolved infections, the
death analogue of the suspected-case ascertainment `p_drc`
([`pooled_ascertainment_model`](@ref)).

The default `Normal(logit(0.9), 0.5)` on the logit scale is deliberately
informative and centred on a high ascertainment: a death is a salient event
in an Ebola response and is reported more reliably than a living suspected
case (`p_drc ≈ 0.75`). The SD 0.5 gives a 90% prior interval of roughly
0.80–0.95, admitting moderate under-ascertainment without letting the death
stream slide to an implausibly low capture. `p_death` is weakly identified
on its own (it trades off with the CFR for the suspected-death level), so it
leans on this prior. The export-death stream and the CFR prior pin the CFR
separately. Pass `ascertainment_prior` to override. Returns
`(; p_death, logit_p_death)`.
"""
@model function death_ascertainment_model(;
        ascertainment_prior = Normal(logit(0.9), 0.5))
    logit_p_death ~ ascertainment_prior
    p_death := logistic(logit_p_death)
    return (; p_death, logit_p_death)
end

"""
Background case-fatality ratio `cfr_bg` for the non-BVD suspected-death
background ([`deaths_model`](@ref)). The renewal joint ties the non-BVD
suspected-death background to the (already identified) non-BVD
suspected-case background `λ_bg` ([`test_positivity_model`](@ref)) rather
than giving deaths a free, outbreak-size-degenerate background rate of
their own: the per-day background suspected deaths are `cfr_bg · λ_bg_v`,
a background CFR applied to the per-day non-BVD suspected-case rate. A
non-BVD suspected case (other severe febrile or haemorrhagic illness that
meets the suspect definition) carries its own fatality risk, and
`cfr_bg` is the share of that background pool that is reported as a
suspected death.

Tying the death background to the case background removes the degeneracy that
keeps a free `λ_bg_death` switched off: the case background is pinned by the
laboratory positivity link, so scaling it by `cfr_bg` gives the death
background a level and time profile without a second free rate competing with
outbreak size. The default `Beta(2, 6)` (mean ≈ 0.25, 90% ≈ 0.05–0.52) is
weakly informative and centred below the BVD CFR (non-BVD suspect illness is
on average less lethal than BVD). Pass `cfr_prior` to override. Returns
`(; cfr_bg)`.
"""
@model function background_cfr_model(; cfr_prior = Beta(2.0, 6.0))
    cfr_bg ~ cfr_prior
    return (; cfr_bg)
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

The default `λ_bg` prior is a half-normal
`truncated(Normal(0, 1.0); lower = 0)`. Its total contribution to the
expected suspected-case count over the grid is `λ_bg · T`, with `T` the
seeding-to-cut-off span. The prior is deliberately informative because
`λ_bg` is degenerate with outbreak size (the per-vintage reported mean
mixes the `p_drc`-scaled BVD increment with `λ_bg · Δt`), so a diffuse
prior lets the background absorb arbitrarily many suspected cases and
resolve at the high end where the deaths and exports streams pin `C_T`.
A background-noise process must not be able to explain more suspected
cases than were ever reported. With SD 1.0 the median background is
≈ 0.67/day and the 95% prior bound ≈ 2.0/day, a modest minority of the
≈ 1077 suspected cases observed by the last stable suspected-case vintage
while still admitting a genuine non-BVD signal. A wider SD (e.g. SD 5)
left a second posterior mode in which the background explains the
majority of suspected cases (positivity ≈ 0.2, background ≈ 2.3× the
observed total). Pass
`lambda_prior` to override. `τ_test` defaults to `Beta(5, 2)`
(mean ≈ 0.71).

The derived per-suspected positivity is exposed inside
[`reported_cases_model`](@ref). The per-test positivity inside
[`confirmed_cases_model`](@ref). Returns `(; λ_bg, τ_test)`.
"""
@model function test_positivity_model(;
        lambda_prior = truncated(Normal(0.0, 1.0); lower = 0),
        fraction_tested_prior = Beta(5.0, 2.0))
    λ_bg ~ lambda_prior
    τ_test ~ fraction_tested_prior
    return (; λ_bg, τ_test)
end

"""
Base treatment-admission probability for the isolation-occupancy stream
([`treatment_flow_model`](@ref)). Samples `p_iso`, the fraction of
ascertained suspected cases that are admitted to and retained in an
isolation/treatment bed at the base (non-BVD rule-out) intensity, so the
modelled bed occupancy is `p_iso` times the survival-convolution of the
admission inflow. BVD suspects are admitted at a higher rate skewed up from
this base by a severity log-odds ([`isolation_severity_model`](@ref)),
since triage admits the sicker patients and BVD presents more severely.

The default `Beta(2, 2)` is weakly informative on `(0, 1)` with mean ½ and
no mass piled at the bounds. `p_iso` is partially confounded with the
length-of-stay mean for the occupancy level (Little's law: mean occupancy
≈ `p_iso · admissions · (E[LOS] + 1)`), so the length-of-stay prior carries
the duration and `p_iso` absorbs the admission/retention fraction. The
length-of-stay also sets the lag and smoothing of occupancy relative to the
inflow, which the daily occupancy series identifies. Pass `p_prior` to
override. Returns `(; p_iso)`.
"""
@model function isolation_admission_model(; p_prior = Beta(2.0, 2.0))
    p_iso ~ p_prior
    return (; p_iso)
end

"""
Severity skew for isolation admission ([`treatment_flow_model`](@ref)).
Samples `δ_iso ≥ 0`, the log-odds by which a BVD suspect is more likely to be
admitted to and retained in an isolation bed than a non-BVD rule-out at the
same base intensity `p_iso` ([`isolation_admission_model`](@ref)), so the BVD
admission probability is `logistic(logit(p_iso) + δ_iso)`. Admission cannot
condition on the unobserved BVD status of a suspect. The skew instead
represents the net effect of severity-based triage, where the sicker
patients are isolated and BVD presents more severely, enriching BVD among
the admitted. The non-negative truncation keeps a BVD suspect at least as
likely to be admitted as a rule-out.

The isolation stream observes only total occupancy, so the skew is weakly
identified from it alone. Its effect is to enrich the long-stay BVD
component of demand, which the occupancy persistence informs only mildly,
so the half-normal `truncated(Normal(0, 0.75); lower = 0)` carries most of
the weight (`δ_iso = 0` recovers a shared admission rate). Pass
`logodds_prior` to override. Returns `(; δ_iso)`.
"""
@model function isolation_severity_model(;
        logodds_prior = truncated(Normal(0.0, 0.75); lower = 0))
    δ_iso ~ logodds_prior
    return (; δ_iso)
end

"""
Isolation/treatment-bed capacity for the supply-limited occupancy stream
([`treatment_flow_model`](@ref)). Samples the number of beds available,
`capacity`, the ceiling the latent bed demand saturates against. Bed
occupancy has been supply-driven (demand has outstripped supply, with
occupancy catching up as capacity is expanded), so the modelled occupancy is
the demand passed through a soft cap at `capacity` rather than tracking
demand directly.

The default `LogNormal(log 450, 0.42)` is weakly informative and positive by
construction (no truncation boundary), with median 450 and a ≈0.44
coefficient of variation, centred on the bed count implied by the reported
occupancy rates (the "Taux d'occupation" gives `capacity = occupancy /
rate ≈ 400–452` over 9–13 June). The capacity is identified by the
implied-capacity series the isolation submodel fits, so the prior only
has to bracket it. A single national capacity is a limitation: it
averages over a growing capacity and cannot represent local saturation
(one province full while another has slack, see
[`bed_capacity_walk_model`](@ref) for the time-varying form). Pass
`capacity_prior` to override. Returns `(; capacity)`.
"""
@model function bed_capacity_model(;
        capacity_prior = LogNormal(log(450.0), 0.42))
    capacity ~ capacity_prior
    return (; capacity)
end

"""
Time-varying isolation/treatment-bed capacity over the daily grid, a
multiplicative random walk: the supply-limited occupancy stream
([`treatment_flow_model`](@ref)) uses `C(t)` as the ceiling the latent
bed demand saturates against on each day. Capacity is not fixed (beds are
being added: SitRep 030 records mattress and bed deliveries and new
treatment centres opening), so the walk tracks the growth a single scalar
capacity ([`bed_capacity_model`](@ref)) cannot.

The walk is a non-centred cumulative log-deviation from a baseline bed
count `C0` on weekly knots, linearly interpolated to the daily grid (the
same parameterisation as the reproduction-number and background walks):
with knot values `\\log C` and knot days `d`,
`C(t) = C0 · exp(\\text{interp}(σ_cap · cumsum(z)))` with `z ~ Normal(0, 1)`
per knot and a tight innovation SD `σ_cap`, keeping capacity a gentle
drift rather than per-day jumps. Knots need far fewer innovations than a
daily walk, avoiding the high-dimensional funnel. The baseline carries
the same weakly-informative `LogNormal(log 450, 0.42)` prior as the
scalar model (median 450 beds, ≈0.44 CV), so `C0` is sampled on the log
scale and the whole capacity `log C(t) = log C0 + walk` is fully
log-scale. The implied-capacity series the isolation submodel fits pins
`C(t)` on the days a rate is published.

Knots run only from `start`, the first day with occupancy or capacity
data. Capacity is flat at `C0` before it. Off-window capacity carries no
likelihood, so walking it adds unidentified innovations that leave the
posterior poorly conditioned, and `start` keeps the knots to the days the
data speaks to. Pass `start = 1` for knots over the whole grid, or `week`
to change the knot spacing. A single national capacity remains a
limitation: it cannot represent local saturation, one province full
while another has slack. Pass `baseline_prior` / `innovation_prior` to
override. Returns `(; C, C0, σ_cap)` with `C` a length-`n` vector.
"""
@model function bed_capacity_walk_model(n::Integer; start::Integer = 1,
        week::Integer = 7,
        baseline_prior = LogNormal(log(450.0), 0.42),
        innovation_prior = truncated(Normal(0.0, 0.05); lower = 0))
    C0 ~ baseline_prior
    σ_cap ~ innovation_prior
    s = clamp(Int(start), 1, n)
    days = knot_days(n; week = week, start = s)
    nb = length(days)
    ## Non-negative innovations, so capacity is non-decreasing: beds are added
    ## over the response and are not taken away, so `C(t)` cannot drop below an
    ## already-reached level. This also keeps the effective ceiling from
    ## jittering down into the observed occupancy.
    z ~ product_distribution(fill(truncated(Normal(0, 1); lower = 0),
        max(nb - 1, 1)))
    steps = σ_cap .* z[1:max(nb - 1, 0)]
    log_knots = vcat(zero(σ_cap), cumsum(steps))
    walk = interpolate_knots(log_knots, days, n)
    C = C0 .* exp.(walk)
    return (; C, C0, σ_cap)
end

"""
Recovery probability for the recovered-among-confirmed stream
([`recovered_model`](@ref)). The fraction of confirmed cases whose outcome
is recovery rather than death is the confirmed-case survival fraction, the
complement of the case-fatality ratio, so it is grounded on the model's
CFR rather than estimated from scratch. The confirmed cases are a slightly
different population from the one the CFR is defined over (they have been
laboratory-confirmed and brought into care), so the survival fraction is the
complement of the CFR adjusted on the log-odds scale by a sampled offset,

```math
p_\\text{recover} = \\operatorname{logistic}\\!\\bigl(
    \\operatorname{logit}(1 - \\mathrm{CFR}) + \\delta_\\text{rec}\\bigr),
```

with `δ_rec ~ Normal(0, 0.5)` centred at zero, so the default recovery
fraction is exactly `1 − CFR` and the data move it only as far as they
support. The offset keeps `p_recover` in `(0, 1)` without a hard clamp and
lets the confirmed-population survival differ modestly from the CFR
complement. `p_recover` is partially confounded with the
confirmation-to-recovery delay for the count of recoveries observed by the
cut-off (a long delay right-censors recoveries that have not yet resolved),
so the delay carries the timing and `p_recover` the eventual survival
fraction. Pass `offset_prior` to override. Returns
`(; p_recover, recovery_offset)`.
"""
@model function recovery_probability_model(CFR::Real;
        offset_prior = Normal(0.0, 0.5))
    recovery_offset ~ offset_prior
    base = clamp(1 - CFR, eps(typeof(CFR)), one(CFR) - eps(typeof(CFR)))
    p_recover := logistic(logit(base) + recovery_offset)
    return (; p_recover, recovery_offset)
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
half-normal scale, with the same informative default as the scalar
`λ_bg` (`truncated(Normal(0, 1.0); lower = 0)` for cases. Pass a tighter
`baseline_prior` for deaths). The per-vintage rate is a multiplicative
log-normal deviation from this baseline,

```math
\\lambda_v = \\lambda_\\mu \\,
    \\exp(\\sigma_{bg}\\, z_v), \\qquad z_v \\sim \\mathcal N(0, 1),
```

with `σ_bg` the pooling SD, passed in rather than sampled here so the
suspected-case and suspected-death streams can share one pooling SD. See
[`background_pooling_model`](@ref), which samples it once at the composer
level. The deviation is multiplicative so the per-vintage rate stays
positive without a clamp and `σ_bg → 0` recovers the scalar baseline
exactly (every `λ_v = λ_mu`). Each stream still samples its own baseline
`λ_mu` and per-vintage deviations `z`. `nv` is the number of vintage
windows. Returns `(; λ, λ_mu, σ_bg, z)` with `λ` a length-`nv` vector of
per-vintage rates.
"""
@model function background_re_model(nv::Integer, σ_bg::Real;
        baseline_prior = truncated(Normal(0.0, 1.0); lower = 0))
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
        pooling_prior = truncated(Normal(0.0, 0.1); lower = 0))
    σ_bg ~ pooling_prior
    return (; σ_bg)
end

"""
Non-BVD background rate as a smooth weekly lognormal random walk over the
surveillance window, the alternative to the per-vintage step random
effect ([`background_re_model`](@ref)). The log-rate follows a
non-centred random walk on weekly knots and is linearly interpolated to
the daily grid, the same parameterisation as the reproduction-number walk
([`rt_walk_model`](@ref)): the background is a slow drift, so a knot per
`week` carries the time variation with far fewer innovations than a daily
walk, avoiding the high-dimensional funnel a daily walk over a long
window opens. The series is gated to zero before the surveillance
`onset` (the non-BVD background does not exist before surveillance
began) and ramps in over the first `onset_ramp` days of the window. With
knot values `\\log\\lambda` and knot days `d`,

```math
\\log \\lambda_d = \\log \\lambda_0 + \\sigma_{rw} \\sum_{s < d} z_s,
\\qquad z_s \\sim \\mathcal N(0, 1),
\\qquad \\lambda_t = 0 \\ \\text{for}\\ t < \\text{onset}.
```

`σ_rw` (passed in, shared across the suspected-case and suspected-death
streams via [`background_pooling_model`](@ref)) is the per-knot innovation
SD on the log scale. A tight prior keeps the background close to constant
(a gentle drift, not week-to-week jumps), which regularises the
background/outbreak-size degeneracy (closing the high-background second
posterior mode that breaks convergence) and keeps the series smooth (so a
death background scaled from it carries no steps). Knots run only over
the surveillance window `[onset, n]`, so the number of innovations is
small. `onset ≤ 1` runs it over the whole grid. Pass `week` to change the
knot spacing. Returns `(; λ, λ_mu, σ_bg)` with `λ` the length-`n` daily
series (zero before `onset`).
"""
@model function background_walk_model(n::Integer, σ_rw::Real;
        onset::Integer = 1, onset_ramp::Integer = 7, week::Integer = 7,
        baseline_prior = truncated(Normal(0.0, 8.0); lower = 0))
    t0 = clamp(Int(onset), 1, n)
    nw = n - t0 + 1
    ## Weekly knots over the window, linearly interpolated to the daily grid
    ## (see [`knot_days`](@ref) and [`interpolate_knots`](@ref)).
    days = knot_days(n; week = week, start = t0)
    nb = length(days)
    ## Half-normal baseline on the natural scale, the same informative prior
    ## as the scalar `λ_bg` ([`test_positivity_model`](@ref)). It bounds the
    ## background level tightly (a lognormal/log-scale level has a heavy right
    ## tail the background/outbreak-size degeneracy exploits to run away), so
    ## the background cannot blow up to explain the suspected stream.
    λ_mu ~ baseline_prior
    z ~ product_distribution(fill(Normal(0, 1), max(nb - 1, 1)))
    ## Smooth multiplicative deviation: a non-centred cumulative (random-walk)
    ## log-deviation from the baseline, anchored at the baseline on the first
    ## knot. A tight `σ_rw` keeps the walk a gentle drift around the bounded
    ## baseline. Interpolated to daily so a death background scaled from it is
    ## smooth.
    steps = σ_rw .* z[1:max(nb - 1, 0)]
    log_knots = vcat(zero(σ_rw), cumsum(steps))
    walk = interpolate_knots(log_knots, days, n)[t0:n]
    λ_window = λ_mu .* exp.(walk)
    ## Linear onset ramp `0 → 1` over the first `onset_ramp` days of the
    ## window, so the gated background grows in from zero instead of stepping
    ## straight to `λ_mu` at the surveillance boundary (which would put a
    ## one-day jump into the suspected-death trajectory scaled from it). The
    ## ramp reaches 1 within the window. `onset_ramp ≤ 1` gives a hard onset
    ## (a step to `λ_mu`).
    rr = clamp(Int(onset_ramp), 1, nw)
    ramp = [min(i, rr) / rr for i in 1:nw]
    λ_window = ramp .* λ_window
    T = eltype(λ_window)
    λ = vcat(zeros(T, t0 - 1), λ_window)
    return (; λ, λ_mu, σ_bg = σ_rw)
end

"""
Confirmation-process sensitivity prior. `Beta(38, 2)` centres near a mean
of 0.95 with a tight spread. Confirmation runs on the altona RealStar
Filovirus Screen RT-PCR, which detects Bundibugyo virus at 11–67 RNA
copies per reaction. The Zaire-specific GeneXpert Ebola assay does not
reliably detect Bundibugyo. A single assay draw is sensitive to about
0.85, but a suspect is confirmed or ruled out through repeat control
tests rather than one PCR, so the effective process sensitivity is
higher (two controls give about 0.98) and `Beta(38, 2)` credits that
process. Under the severe-first backlog the
first vintage's analysed batch is near-pure BVD (`q ≈ 1`), so the v1
positivity ≈ `s` identifies the sensitivity from the early data. Returns
`(; s_test)`.
"""
@model function test_sensitivity_model(;
        sensitivity_prior = Beta(38.0, 2.0))
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
Confirmed-positives overdispersion prior. The confirmed positives in each
laboratory window are scored as an overdispersed `BetaBinomial` of the
observed analysed denominator (see [`safe_betabinomial`](@ref) and
[`confirmed_cases_model`](@ref)) rather than a plain `Binomial`, because
the per-window positivity `p_pos` is a smooth pooled / composition-linked
curve that does not capture the day-to-day laboratory batching and
within-window positivity heterogeneity the confirmed counts carry. A
plain `Binomial` on denominators of several hundred specimens gives
predictive intervals far too tight, so the confirmed stream is
systematically under-covered. The intra-window correlation
`ρ ∈ (0, 1)` inflates the window variance to
`n·p·(1 − p)·(1 + (n − 1)·ρ)`, with `ρ → 0` recovering the `Binomial`.
The default `Beta(1, 24)` (mean ≈ 0.04, 90% ≈ 0.002–0.12) is weakly
informative: it favours a small overdispersion and shrinks toward the
`Binomial` when the data support it, while a single scalar identified
across the laboratory windows lets the confirmed positives themselves
set the spread. Returns `(; ρ)`.
"""
@model function confirmed_overdispersion_model(;
        overdispersion_prior = Beta(1.0, 24.0))
    ρ ~ overdispersion_prior
    return (; ρ)
end

"""
Report-to-laboratory-confirmation (lab-turnaround) delay submodel. The
delay from a suspected case being reported to its specimen being
laboratory confirmed, discretised to a daily PMF over lags `0 … nmax`
by [`censored_delay_model`](@ref) so it convolves cleanly onto the
renewal onsets. The mean and SD carry weakly-informative priors centred
on a short turnaround with a heavy right tail allowing for specimen
shipment to a confirmatory lab. No per-sample outbreak data grounds this
delay, so the likelihood does not identify the turnaround mean or SD. The
priors are therefore kept tight around the documented turnaround belief
(mean ≈ 4.5 d, SD ≈ 4 d) rather than wide, since a wide prior on an
unidentified nuisance delay only makes the sampler wander it (it was the
worst-mixing quantity in the joint, dragging the confirmation PMFs convolved
from it). Returns `(; pmf, dist, mean, sd)`.
"""
@model function lab_delay_model(
        nmax::Integer = cdf_nmax(lognormal_meansd(4.5, 4.0));
        mean_prior = truncated(Normal(4.5, 1.0); lower = 1),
        sd_prior = truncated(Normal(4.0, 0.75); lower = 1))
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
Severity-enrichment prior for the composition-linked confirmed positivity
(`positivity_link = :composition` in [`confirmed_cases_model`](@ref)). In
that mode the per-window tested BVD share is not a free random effect. It
is the suspect-pool composition `φ_v = (p_drc · BVD)_v / ((p_drc · BVD)_v +
λ_bg_v)` over each laboratory window, upsampled by a severity enrichment
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

`δ₀` is the early severity log-odds enrichment of BVD. Lower-truncated at 0
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
Death testing fraction `τ_death`, the fallback for the death-only composer
([`confirmed_deaths_only_model`](@ref)), which has no case stream to set the
death testing volume from. It thins the suspected deaths to a death "analysed"
volume at the case testing rate, drawing `τ_death` from the same prior as the
case testing fraction (`Beta(5, 2)`, mean ≈ 0.71). The full joint instead
scales the modelled case analysed volume (see [`confirmed_deaths_model`](@ref)
and [`death_testing_scaling_model`](@ref)) and does not draw this submodel.
Pass `fraction_prior` to override. Returns `(; τ_death)`.
"""
@model function death_testing_fraction_model(; fraction_prior = Beta(5.0, 2.0))
    τ_death ~ fraction_prior
    return (; τ_death)
end

"""
Death testing-intensity scaling for the confirmed-death volume in the joint
([`confirmed_deaths_model`](@ref)). The death analysed volume is the modelled
case analysed volume carried at the per-day suspected death-to-case ratio,
times this scaling. That ratio already carries the suspect-pool severity and
the suspected-death level, so the scaling is the per-suspect testing-intensity
difference between deaths and living suspects alone. No death-testing data
grounds it, so it is a tight log-normal centred on one (`LogNormal(0, 0.25)`,
median 1, 90% ≈ 0.66–1.51): deaths are tested at the case intensity unless the
confirmed-death counts pull the scaling off one. Pass `scaling_prior` to
override. Returns `(; scaling)`.
"""
@model function death_testing_scaling_model(;
        scaling_prior = LogNormal(0.0, 0.25))
    scaling ~ scaling_prior
    return (; scaling)
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
Partially-pooled negative-binomial dispersions for the `n_streams`
passive-surveillance count streams in the joint model (suspected cases,
suspected deaths, confirmed cases and confirmed deaths). Each stream draws
its own dispersion from a shared population, so heterogeneous streams (a
handful of deaths versus hundreds of suspects versus a daily laboratory
volume) no longer share one global `k` that the dominant stream pulls
around, while the sparse streams still borrow strength through the common
hyper-parameters rather than going noisy on a fully independent draw.

The pooling is on the `log(1/sqrt(k))` scale: a population mean `μ_log`, a
pooling SD `τ`, and per-stream deviations, with `k_s = 1 / inv_sqrt_k_s^2`. The
default is the **centred** form, `log(1/sqrt(k))_s ~ Normal(μ_log, τ)` drawn
directly: the passive-surveillance streams are data-rich (hundreds of
suspected/confirmed counts), so each stream's dispersion is strongly informed,
and the non-centred form `inv_sqrt_k_s = exp(μ_log + τ z_s)` then funnels
(`z_s = (log_isk_s − μ_log)/τ` diverges as `τ → 0`), stretching NUTS
trajectories. On the joint, centring removes that funnel — worst dispersion
bulk-ESS ≈ 102 → 156, divergences 5 → 2 and ~10% faster wall-clock at a 150×2
fit — for the same posterior over `k`. Pass `centred = false` for the
non-centred form, which is the better choice when the streams are data-poor and
prior-dominated. The population mean is centred on the shared `1/sqrt(k)` prior
of [`surveillance_dispersion_model`](@ref) (`exp(μ_log)` near 0.6), and the
half-normal `τ` keeps the per-stream dispersions close unless the data pull
them apart (`τ = 0` collapses every stream to the population value, the
shared-`k` model). Returns `(; k, inv_sqrt_k, k_pop, μ_log, τ)` with `k` a
length-`n_streams` vector.
"""
@model function pooled_dispersion_model(n_streams::Integer;
        mean_prior = Normal(log(0.6), 0.33),
        sd_prior = truncated(Normal(0, 0.6); lower = 0),
        centred::Bool = true)
    μ_log ~ mean_prior
    τ ~ sd_prior
    m = max(n_streams, 1)
    if centred
        ## Draw each stream's `log(1/sqrt(k))` directly from the population.
        ## `eps` floors the SD so a `τ ≈ 0` draw stays a proper distribution.
        log_isk ~ product_distribution(
            fill(Normal(μ_log, τ + eps(typeof(τ))), m))
        inv_sqrt_k = exp.(log_isk[1:n_streams])
    else
        z ~ product_distribution(fill(Normal(0, 1), m))
        inv_sqrt_k = exp.(μ_log .+ τ .* z[1:n_streams])
    end
    k = 1.0 ./ (inv_sqrt_k .^ 2 .+ eps(eltype(inv_sqrt_k)))
    k_pop = 1.0 / (exp(μ_log)^2 + eps(typeof(float(μ_log))))
    return (; k, inv_sqrt_k, k_pop, μ_log, τ)
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
[`exports_model`](@ref) and [`exports_deaths_model`](@ref). This is the
composer default. The shared mean defaults to a reporting fraction of
0.75 (logit scale), reflecting the active case-finding of a declared
Ebola response. A lower ascertainment inflates the inferred infections
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

## --- Patch (multi-population) models -----------------------------------

"""
Reproduction numbers for several spatial patches (Ituri, Nord-Kivu,
Sud-Kivu): a common national trend plus per-patch deviations that are free
to vary in space AND in time, drawn from a multivariate-normal random walk
with a learned cross-patch correlation.

```math
\\log R_{p,t} = \\mu(t) + \\delta_p(t), \\qquad
\\textstyle\\sum_p \\delta_p(t) = 0,
```

```math
\\Delta\\boldsymbol{\\delta}(t_k) \\sim
    \\mathrm{MVN}\\bigl(\\mathbf{0},\\, \\Sigma\\bigr), \\qquad
\\Sigma = \\mathrm{diag}(\\sigma_\\delta)\\, \\Omega\\,
          \\mathrm{diag}(\\sigma_\\delta), \\qquad
\\Omega \\sim \\mathrm{LKJ}(2),
```

on the same weekly knots as the national walk, linearly interpolated to the
daily grid.

### A common trend, not independent walks

`μ(t)` is the existing national weekly-knot walk ([`rt_walk_model`](@ref)),
kept intact: the national streams see exactly the `Rt` process the headline
model fits, and it is the target the provinces pool toward.

That pooling target is the point. An unstructured multivariate walk — one
free `log Rt` trajectory per province, correlated through an LKJ prior and
no common trend — is strictly more flexible, but it shrinks the wrong way.
`LKJ(η)` has density proportional to `det(Ω)^(η-1)`, maximised at `Ω = I`,
so `η = 2` mildly favours INDEPENDENT provincial walks. The provinces are
not equally observed: over the fitted window Nord-Kivu contributes 74
laboratory positives and Sud-Kivu contributes none at all. Shrinking toward
independence estimates their `Rt` almost entirely from that, while shrinking
toward a common trend lets them borrow strength from Ituri and deviate only
where the data insist. For a meta-population under one national response,
partial pooling toward a shared trend is the right inductive bias.

Writing the model this way costs nothing in generality: with `μ(t)` present,
the deviation covariance `Σ` is still free, so the cross-patch correlation is
LEARNED rather than assumed. This IS a multivariate-normal random walk — it
just carries a common factor rather than leaving the correlation structure to
carry it.

### Sum-to-zero, not a reference patch

The deviations are centred at every knot, so no province is privileged.
Fixing `δ_1 ≡ 0` instead (reference coding) would also identify the model,
but it forces the primary patch to have no idiosyncratic deviation at all —
Ituri would BE the national trend by construction while the other provinces
carry their own noise. That asymmetry is an artefact of the identifiability
fix, not epidemiology.

Centring introduces one redundant coordinate per knot (the mean of the raw
innovations, which the likelihood never sees). It is drawn from its proper
prior and is Gaussian and well-conditioned, so it costs a few cheap sampled
dimensions rather than a posterior ridge.

### What the data can and cannot identify here

The composition of the confirmed cases identifies the CONTRAST between
provinces. With three patches that is essentially one number, the Ituri /
Nord-Kivu contrast, since Sud-Kivu carries no signal. Expect `Ω` to be
largely prior-driven and Sud-Kivu's `Rt` to be pinned by the deviation prior
rather than by data; that is honest, and it is why `Σ` is given a proper
shrinkage prior rather than a flat one.

`σ_δ → 0` recovers a common `Rt` shape shared by every province (a fixed
ratio between them). It is a special case of this model, not an assumption
baked into it: whether the provinces are really moving together is
ESTIMATED. `σ_δ` is therefore the headline spatial diagnostic — a posterior
pushed away from zero is direct evidence that provincial `Rt` trajectories
are separating, which is exactly what a response concentrated on the Ituri
epicentre would produce.

Returns the Rt matrix `(n_patches × n)`, the national trend, the full
deviation trajectory `δ_patch` `(n_patches × n)`, the per-patch deviation
scales and the correlation matrix.
"""
@model function patch_rt_model(n::Integer, n_patches::Integer,
        log_R0_base::Real;
        breakpoint::Union{Missing, Real} = missing,
        week::Integer = 7,
        rt_start::Integer = 1,
        rt_walk_start::Integer = rt_start,
        rt = rt_walk_model,
        region_sd_prior = truncated(Normal(0, 0.15); lower = 0),
        region_drift_sd_prior = truncated(Normal(0, 0.05); lower = 0),
        lkj_prior = LKJCholesky(max(n_patches, 2), 2.0),
        region_offset_prior = Normal(0, 1))
    ## Common national trend: the existing single-patch walk, unchanged.
    ## `rt_walk_start` maps to `rt_start` in the inner model, matching the
    ## convention in [`infection_model`](@ref).
    rt_state ~ to_submodel(
        rt(n, log_R0_base; breakpoint, rt_start = rt_walk_start), false)
    Rt_national = rt_state.Rt
    log_Rt_national = log.(Rt_national)
    ## The deviations live on the SAME weekly knots as the national walk, so
    ## both processes are described at the same resolution.
    days = knot_days(n; week, start = rt_walk_start)
    nb = length(days)
    ## Deviation scales (one per patch) and their cross-patch correlation.
    ## `LKJCholesky` samples the Cholesky FACTOR directly, so the
    ## decomposition never lands on the AD tape.
    σ_level ~ region_sd_prior
    σ_δ ~ product_distribution(fill(region_drift_sd_prior, n_patches))
    Ω_L ~ lkj_prior
    L = Ω_L.L
    ## Standard-normal draws for the level and for each knot's innovation.
    z_level ~ product_distribution(fill(region_offset_prior, n_patches))
    z_drift ~ product_distribution(
        fill(region_offset_prior, max(n_patches * (nb - 1), 1)))
    Tp = promote_type(eltype(Rt_national), typeof(float(σ_level)),
        eltype(σ_δ), eltype(L), eltype(z_level), eltype(z_drift))
    ## Correlated deviations, centred at every knot so the patches sum to
    ## zero and no province is privileged.
    δ_knots = zeros(Tp, n_patches, nb)
    lvl = zeros(Tp, n_patches)
    @inbounds for i in 1:n_patches
        acc = zero(Tp)
        for j in 1:i
            acc += L[i, j] * z_level[j]
        end
        lvl[i] = σ_level * acc
    end
    lvl_bar = sum(lvl) / n_patches
    @inbounds for i in 1:n_patches
        δ_knots[i, 1] = lvl[i] - lvl_bar
    end
    innov = zeros(Tp, n_patches)
    @inbounds for k in 2:nb
        for i in 1:n_patches
            acc = zero(Tp)
            for j in 1:i
                acc += L[i, j] * z_drift[(k - 2) * n_patches + j]
            end
            innov[i] = σ_δ[i] * acc
        end
        innov_bar = sum(innov) / n_patches
        for i in 1:n_patches
            δ_knots[i, k] = δ_knots[i, k - 1] + (innov[i] - innov_bar)
        end
    end
    ## Interpolate each patch's deviation to the daily grid and build Rt.
    δ_patch = zeros(Tp, n_patches, n)
    Rt_matrix = zeros(Tp, n_patches, n)
    @inbounds for p in 1:n_patches
        δ_daily = interpolate_knots(δ_knots[p, :], days, n)
        for t in 1:n
            δ_patch[p, t] = δ_daily[t]
            Rt_matrix[p, t] = exp(log_Rt_national[t] + δ_daily[t])
        end
    end
    ## Cross-patch correlation matrix, reconstructed from its factor for
    ## reporting (Ω = L Lᵀ).
    Ω = zeros(Tp, n_patches, n_patches)
    @inbounds for i in 1:n_patches, j in 1:n_patches

        acc = zero(Tp)
        for k in 1:min(i, j)
            acc += L[i, k] * L[j, k]
        end
        Ω[i, j] = acc
    end
    return (; Rt_matrix, Rt_national, log_Rt_national, δ_patch, δ_knots,
        σ_level, σ_δ, Ω, sigma_rw = rt_state.sigma_rw,
        log_R0 = rt_state.log_R0,
        intervention_effect = rt_state.intervention_effect)
end

"""
Multi-patch latent infection process. Runs a renewal equation per spatial
patch (Ituri, Nord-Kivu, Sud-Kivu) with a shared generation interval, a
shared incubation period, and an optional between-patch importation
kernel.

### Structure

For patch `p` on day `t`:

```math
I_{p,t} = R_{p,t} \\cdot \\sum_{s \\ge 1} I_{p,t-s}\\, g_s
          + \\varepsilon \\sum_{q \\neq p} K_{p,q}\\, I_{q,t-1}
```

with `g_s` the shared generation-interval PMF (sampled once, the biology
of transmission does not depend on province), `R_{p,t}` from
[`patch_rt_model`](@ref), `K` the importation kernel, and `ε` the
importation intensity.

### Importation

`ε` is sampled only when a non-zero `importation_kernel` is supplied. The
default kernel is all-zero, so by default the patches are uncoupled and no
`ε` enters the parameter space. This is deliberate. There is no mobility
or origin-destination data for this outbreak, and the per-province case
composition is flat over the observed window, which leaves importation and
the secondary-patch seeds confounded: any `(ε, seed)` pair that reproduces
the observed Nord-Kivu level fits the data equally well. Sampling `ε`
against an all-zero kernel would add a dimension the likelihood never
touches, giving a parameter whose posterior is exactly its prior. With the
default, each secondary patch is explained by its own seed and its own
`R_t`. Passing a kernel turns the coupling on for a sensitivity analysis.

### Seeding

The primary patch (`p = 1`, Ituri) uses the existing cryptic-exponential
seed: growth at the sampled molecular-clock rate `r` over the cryptic
window, reaching `seed_at_renewal_start(C_T)` at the renewal start.

Each secondary patch is seeded as a FRACTION of the primary patch's seed
(`seed_fraction_prior`), not as an absolute count. The fraction stands in
for the unobserved early introductions from Ituri, and it is the natural
scale because it is what the data speak to: with importation off, the
per-province *level* of the case split is set by the relative seed, so the
observed Nord-Kivu share (~9% of confirmed, near-constant across the
window) maps almost directly onto a seed fraction of roughly the same size.

This matters more than it looks. An absolute seed prior on the secondary
patches — an earlier version of this model used
`N⁺(0.01, 0.01)` — puts them ~4 orders of magnitude below the primary
patch's `2^m ≈ 164`, a seed ratio of about 13,700 : 1. With importation off
by default, a secondary patch then has only two routes to infections: its
own seed and its own `R_t`. If the seed is pinned that far below what the
data need, the log-Rt deviation `δ_p` is forced to absorb the entire
shortfall: reaching a 9% Nord-Kivu share requires `δ ≈ 1.0` (an `Rt` ratio
of 2.7), which is a 3.4-sigma draw on the deviation prior, and the mapping
from `δ` to the share is a knife-edge (`δ = 0.5` gives 0.5%, `δ = 1.0`
gives 54%). The reported provincial `Rt` difference would then be an
artefact of the seed prior rather than an epidemiological finding, which is
precisely the quantity the patch model exists to estimate.

Parameterising the seed as a fraction decouples the two: the seed explains
the LEVEL of the provincial split and `δ_p` is identified by its TIME
TREND. That is the decomposition the data actually support.

The default `LogNormal(log(0.05), 1)` has a median of 5% of the primary
seed and a 90% interval of roughly 1% to 26%, so it spans the observed
share comfortably without asserting it.

### Returns

The per-patch state, plus the national aggregates the observation models
and the headline summaries consume: `infections_total`, `cumulative_total`,
`C_T` (the national cut-off cumulative), and `Rt_national_implied`, the
incidence-weighted aggregate reproduction number obtained by inverting the
renewal equation on the summed infections ([`implied_national_Rt`](@ref)).
`R_T`, `r`, `T` and `doubling_time` mirror [`infection_model`](@ref) so a
patch chain carries the same headline quantities as a single-patch one.
"""
@model function patch_infection_model(n::Integer, n_patches::Integer;
        breakpoint::Union{Missing, Real} = missing,
        rt_start::Integer = 1,
        rt_walk_start::Integer = rt_start,
        rt = patch_rt_model,
        gi = generation_interval_model,
        growth = exponential_growth_model,
        gi_nmax::Integer = cdf_nmax(Gamma(2.71, 5.65)),
        importation_kernel::AbstractMatrix = province_importation_kernel(
            PROVINCE_POPULATIONS[1:min(n_patches, end)]),
        importation_epsilon_prior = Beta(1, 100),
        seed_fraction_prior = LogNormal(log(0.05), 1.0),
        incubation = (nmax) -> censored_delay_model(nmax;
            mean_prior = truncated(Normal(6.3, 0.54); lower = 1),
            sd_prior = truncated(Normal(3.5, 0.8); lower = 1)),
        incubation_nmax::Integer = cdf_nmax(lognormal_meansd(6.3, 3.5)))
    ## 1. Shared generation interval.
    gi_state ~ to_submodel(gi(gi_nmax))
    g = gi_state.g
    ## 2. ONE growth source, as in [`infection_model`](@ref): the prior is on
    ##    the cryptic growth rate `r`, and the established `R0` (the walk
    ##    base) is derived forward from `r` and the generation interval
    ##    through Euler-Lotka.
    growth_state ~ to_submodel(growth())
    r_clock = growth_state.r
    R0 = r_to_R0(r_clock, g)
    ## 3. Per-patch Rt: national trend plus per-patch deviations.
    rt_state ~ to_submodel(
        rt(n, n_patches, log(R0); breakpoint, rt_start, rt_walk_start), false)
    Rt_matrix = rt_state.Rt_matrix
    δ_patch = rt_state.δ_patch
    ## 4. Per-patch seeds. The primary patch takes the cryptic exponential.
    ##    Each secondary patch takes a FRACTION of that seed, which is the
    ##    scale the data speak to: with importation off, the relative seed
    ##    sets the LEVEL of the provincial case split, leaving `δ_p` to be
    ##    identified by its TIME TREND. An absolute seed prior pinned far
    ##    below the primary's `2^m` would force `δ_p` to absorb the whole
    ##    level difference, making the reported provincial Rt gap an artefact
    ##    of the seed prior (see the docstring).
    renewal_start = clamp(rt_start, 1, n)
    τ_obs = n - renewal_start
    seed0_primary = seed_at_renewal_start(growth_state.C_T)
    seed_primary = seed_infections(seed0_primary, r_clock, renewal_start)
    seed_fraction ~ product_distribution(
        fill(seed_fraction_prior, max(n_patches - 1, 1)))
    Tp = promote_type(eltype(Rt_matrix), eltype(g), typeof(float(r_clock)),
        eltype(seed_fraction), typeof(float(seed0_primary)))
    seeds_matrix = zeros(Tp, n_patches, renewal_start)
    @inbounds for j in 1:renewal_start
        seeds_matrix[1, j] = seed_primary[j]
    end
    @inbounds for p in 2:n_patches
        ## A scaled copy of the primary patch's cryptic curve: the secondary
        ## province is a fraction of the same epidemic, so it grows at the
        ## same clock rate `r` over the cryptic window.
        s_p = seed_infections(
            seed_fraction[p - 1] * seed0_primary, r_clock, renewal_start)
        for j in 1:renewal_start
            seeds_matrix[p, j] = s_p[j]
        end
    end
    ## 5. Importation intensity, only when the kernel actually couples the
    ##    patches (see the docstring: unidentified against an all-zero
    ##    kernel, and confounded with the secondary seeds in any case).
    coupled = any(!iszero, importation_kernel)
    ε = zero(Tp)
    if coupled
        ε ~ importation_epsilon_prior
        importation_epsilon := ε
    end
    ## 6. Multi-patch renewal.
    infections_matrix = patch_infections(Rt_matrix, g, seeds_matrix,
        importation_kernel, ε)
    ## 7. Per-patch cumulatives and the national aggregate.
    cumulative_matrix = zeros(Tp, n_patches, n)
    @inbounds for p in 1:n_patches
        acc = zero(Tp)
        for t in 1:n
            acc += infections_matrix[p, t]
            cumulative_matrix[p, t] = acc
        end
    end
    C_T_patch = [@inbounds(cumulative_matrix[p, n]) for p in 1:n_patches]
    infections_total = vec(sum(infections_matrix; dims = 1))
    cumulative_total = cumsum(infections_total)
    ## 8. Per-patch onsets through the SHARED incubation PMF.
    inc_state ~ to_submodel(incubation(incubation_nmax))
    onsets_matrix = zeros(Tp, n_patches, n)
    @inbounds for p in 1:n_patches
        @views onsets_matrix[p, :] = convolve_delay(
            infections_matrix[p, :], inc_state.pmf)
    end
    ## 9. Aggregate reproduction number. Inverting the renewal equation on
    ##    the summed infections gives the incidence-weighted mean of the
    ##    patch `Rt`s: with `I_{p,t} = R_{p,t} · force_{p,t}`, summing over
    ##    patches gives `I_t / Σ_p force_{p,t} = Σ_p R_{p,t} force_{p,t} /
    ##    Σ_p force_{p,t}`. This is the `Rt` that reproduces the national
    ##    trajectory, so it is the one the headline `R_T` reports.
    Rt_national_implied = implied_national_Rt(infections_total, g)
    R_T = @inbounds(Rt_national_implied[n])
    ## 10. Headline quantities, mirroring [`infection_model`](@ref) so a
    ##     patch chain summarises exactly like a single-patch one.
    r = euler_lotka_r(R_T, g)
    T_total = growth_state.T + τ_obs
    return (; infections_matrix, cumulative_matrix, onsets_matrix,
        Rt_matrix, δ_patch, C_T_patch,
        σ_level = rt_state.σ_level,
        σ_δ = rt_state.σ_δ,
        Ω = rt_state.Ω,
        infections_total, cumulative_total,
        Rt_national = rt_state.Rt_national,
        Rt_national_implied, g, R0, r0 = r_clock, r, R_T,
        m = growth_state.m, τ = growth_state.τ,
        T = T_total, C_T = @inbounds(cumulative_total[n]),
        doubling_time = doubling_time(r),
        seed_at_renewal_start = seed0_primary, seed_fraction,
        seeding_age = seeding_age(cumulative_total, n),
        incubation_pmf = inc_state.pmf)
end
