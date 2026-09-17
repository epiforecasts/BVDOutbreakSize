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
under Mooncake, so this is the AD-safe discretisation route for every delay
in the renewal convolutions. The mean and SD carry weakly-informative
priors, so the delay is estimated rather than fixed. Returns
`(; pmf, dist, mean, sd)`.
"""
@model function censored_delay_model(nmax::Integer; mean_prior, sd_prior)
    delay_mean ~ mean_prior
    delay_sd ~ sd_prior
    dist = lognormal_meansd(delay_mean, delay_sd)
    return (; pmf = discretise_censored(dist, nmax), dist,
        mean = delay_mean, sd = delay_sd)
end

"""
Generation-interval submodel, on the Gamma shape `α` and scale `θ`. The
source is the Ebola virus disease serial interval as a generation-time
proxy (mean 15.3 d, SD 9.3 d; WHO Ebola Response Team 2014, NEJM), which
maps once to `α ≈ 2.71` and `θ ≈ 5.65` (`α = (mean/sd)²`,
`θ = sd²/mean`). The priors are centred there,
`α ~ Normal⁺(2.71, 0.7)` and `θ ~ Normal⁺(5.65, 1.5)`, lower-truncated to
keep the Gamma well defined. The SDs propagate the source's reported
uncertainty, the NEJM serial-interval mean carrying a 95% CI of
13.0–17.6 d, an SD on the mean of ≈1.17 d.

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
([`discretise_censored`](@ref)), keeping the lag-0 bin, since an
onset-to-event delay can be same-day unlike the generation interval. This
carries a line-list delay reanalysis through on its natural parameters with
the reported posterior uncertainty. Returns
`(; pmf, dist, mean, sd, alpha, theta)`.
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
their PMFs. This matches the companion line-list reanalysis, which fits the
atomic components rather than onset→death directly, so each atomic delay
keeps its own Gamma shape and scale prior with the reanalysis's reported
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
[wilson1931](@citet) approximation to the continuous
median of a `Gamma` as a function of its `mean` and `sd`:
`median ≈ mean·(1 − sd²/(9·mean²))³`.
Smooth in the mean and SD with no quantile inversion, so it enters a
gradient-based likelihood cleanly, and accurate to a few percent for the
shapes here. Used to match the confirmed onset-to-sample convolution's
continuous median to the cohort's reported median.
"""
gamma_median_wh(mean::Real, sd::Real) = mean * (1 - sd^2 / (9 * mean^2))^3

"""
Onset-to-sample prior configuration from the NEJM DRC 2026 BVD cohort
[akilimali2026](@cite). The confirmed-positive onset-to-sample interval
(N = 129) was estimated as a continuous Gamma through the `epidist` marginal
model [epidist](@cite), correcting for double interval censoring and right
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
cumulative-sum form, with standard-normal innovations scaled by `sigma_rw`
and accumulated, avoiding the funnel geometry of the centred recursion.
Daily log-`R_t` is the linear interpolation between knots
([`interpolate_knots`](@ref)). An intervention at `breakpoint` (e.g. the
first WHO situation report) adds a sampled effect `intervention_effect`
shaped by a logistic ramp ([`sigmoid_ramp`](@ref)) of scale `ramp` (default
21 days, roughly the time a response takes to bite), so transmission
changes gradually rather than instantly. `breakpoint = missing` drops the
term. `Rt = exp.(log_Rt)`.

The walk base `log_R0` is not sampled here. It is passed in, derived
forward from the sampled growth rate `r` and the generation interval
through Euler–Lotka (`R0 = r_to_R0(r, g)` in [`infection_model`](@ref)),
so the prior sits on the growth rate instead (see
[`exponential_growth_model`](@ref)). That single growth source pins both
the established reproduction number and, through the renewal seeding, the
cryptic exponential phase. The grid days before the renewal start are
filled by the analytic cryptic exponential and are unused by the walk,
which clamps to `R0` before its first knot.

The random-walk step SD prior is a half-normal SD 0.1, so the weekly
log-`R_t` is unlikely to change by more than about 20% (two SD ≈ 0.2) from
one week to the next. The walk starts at `rt_start`, so every knot sits in
the observed window rather than drifting over the unobserved pre-report
stretch.

The intervention effect is constrained non-positive
(`truncated(Normal(0, 0.4); upper = 0)`). A declared WHO response (case
finding, isolation, vaccination) can only reduce transmission or leave it
unchanged. The half-normal admits anything from no effect (mode) to a
substantial decline.

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
    ## The established `R0` at the genetic bound is the base the walk grows
    ## from. It is derived and passed in, not sampled here, and tracked as a
    ## deterministic so it stays available on the chain.
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
Samples the exponential growth rate `r` (the primary epidemiological
assumption, placed on the genetic doubling time) and the generation count
`m`, then exposes

```math
\\tau = \\log 2 / r,\\qquad T_\\text{cryptic} = m\\,G,
\\qquad C_T = e^{r T_\\text{cryptic}},
```

as deterministics, with `G` the mean generation interval. `T = m·G` is the
cryptic-phase duration (origin → renewal start). `m` counts the
transmission generations during the cryptic phase, and a single import
grows over them to a daily incidence `C_T` at the renewal start. The
composer ([`infection_model`](@ref)) adds the observation span
`τ_obs = n − renewal_start` to get the total outbreak age
`T_total = m·G + τ_obs`, which carries the genetic seeding bound
([`genetic_seeding_model`](@ref)).

The growth rate carries the prior
`r ~ LogNormal(log(log2 / M_PRIOR_DOUBLING_DAYS), 0.40)`, with median
doubling time (11.7 d) matching the BEAST X estimate (mbalaplacide2026,
exponential growth model, 95% HPD 6.8–17.5). The log-SD 0.40 is wider than
the ≈0.24 that HPD implies, because the HPD is conditional on a
single-rate coalescent, which the field epidemiology contradicts
(kupferschmidt2026), and an independent reanalysis puts the doubling time
at 15.2–24.5 d (cuomodannenburg2026). The induced doubling-time prior is
`LogNormal(log 11.7, 0.40)`, 5.3–25.6 d at 95%. The first reproduction
number is derived forward from this `r` and our generation interval through
Euler–Lotka (`R0 = r_to_R0(r, g)` in [`infection_model`](@ref)), so the
cryptic exponential phase and the established renewal share one growth
source.

`m ~ truncated(Normal(2.75, 1.2); lower = 0)` counts the transmission
generations between the index infection and the renewal start, so the
origin sits `T = m · G` days back and the cryptic phase grows one infection
per day there to `C_T = exp(r · T)` per day at the renewal start.

The centre puts the origin in mid-February 2026, 2.75 generation intervals
before the renewal start. The 90% prior origin runs mid-January to
mid-March, so the traced 25 January 2026 index death (kupferschmidt2026)
sits at about the 87th percentile rather than at the centre. It is the
earliest chain field work reached, which bounds the origin rather than
dating it. The genetic TMRCA (mbalaplacide2026) is a lower bound consistent
with an origin that early. The 99th percentile seed is about 150 infections
per day, against a fitted outbreak of order ten thousand in total.

In the renewal, `C_T` is the prior seed at the renewal start, which the
renewal recursion grows forward under `R_t`. Pass `m_prior` to override. Do
not pass a centre from [`m_prior_centre`](@ref), which is scaled for the
backfill parameterisation where `2^m` is the cut-off cumulative case total
and would give a renewal seed of order half a million per day. Returns
`(; τ, r, m, T, C_T, G)`.
"""
@model function exponential_growth_model(g::AbstractVector;
        r_prior = LogNormal(log(log(2) / M_PRIOR_DOUBLING_DAYS), 0.40),
        m_prior = truncated(Normal(2.75, 1.2); lower = 0))
    r ~ r_prior
    m ~ m_prior
    ## Mean generation interval, the unit `m` is counted in. `g` is indexed
    ## from one day, so the mean is `Σ i·g[i]`.
    G := sum(i * g[i] for i in eachindex(g))
    τ := log(2) / r
    ## Outbreak age is generations times the generation interval, so it does
    ## not depend on `r`.
    T := m * G
    ## Daily incidence at the renewal start, grown from one infection per day
    ## at the origin over `T` days at the cryptic rate.
    C_T := exp(r * T)
    return (; τ, r, m, T, C_T, G)
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
[`exponential_growth_model`](@ref)), derives the established reproduction
number `R0` (= the first `R_t`) forward from that `r` and the generation
interval through Euler–Lotka (`R0 = r_to_R0(r, g)`), and passes `log R0` as
the walk base to the reproduction-number submodel. The cryptic phase and
the established renewal therefore share one growth source.

The renewal runs only over the observation window
`[renewal_start, cut-off]`, where `renewal_start` is the genetic-TMRCA grid
day `rt_start` (the day the reproduction-number walk starts, before which
`R_t` is held flat). The cryptic exponential phase from the origin to the
renewal start is analytic and off the renewal grid except for the days
needed as recursion history. The seed at the renewal start is the daily
incidence the cryptic phase reaches, `C_T = e^{r·m·G}`
([`seed_at_renewal_start`](@ref)), where `m` counts the transmission
generations during the cryptic phase. The magnitude is referenced to the
origin, so a larger `r` raises both the seed and the derived `R0` and the
two compound. A cut-off-referenced seed `C_T e^{−r·τ_obs}` would instead put
`r` into the seed and the renewal growth in opposing directions, opening a
flat ridge along which `R0` slides to 1. The grid days
`1 … renewal_start` are filled smoothly by the cryptic exponential curve at
rate `r` ending at `C_T` ([`seed_infections`](@ref)), giving the recursion a
full generation interval of differentiable history. The renewal recursion
([`renewal_infections`](@ref)) then grows the trajectory over
`renewal_start+1 … n` under the time-varying `R_t`.

The total outbreak age is `T = m·τ + τ_obs` (cryptic duration plus the
observation span `τ_obs = n − renewal_start`). The genetic seeding bound is
applied to this total `T` at the composer. The renewal start sits a small
lead after the genetic TMRCA day, past the TMRCA uncertainty where
sustained transmission is confident, so `τ_obs < tmrca_days` and the
censored bound `tmrca ~ censored(Normal(T, sd); upper = tmrca_days)` stays
informative. It pulls the origin to sit at or before the MRCA, so the
cryptic duration `m·τ` cannot be too short.

The realised cut-off size is `C_T = cumulative[n]`. The `breakpoint` is
forwarded to the reproduction-number submodel. Returns
`(; infections, cumulative, Rt, g, seed_at_renewal_start, m, τ, R0, r0, r,
doubling_time_initial, T, C_T, C_T_prior, doubling_time, seeding_age)`,
where `r`/`doubling_time` are the current growth derived from the cut-off
reproduction number `Rt[n]` through forward Euler–Lotka (so `r` is
sign-consistent with `R_T := Rt[n]` by construction), `r0` the cryptic rate
implied by `R0`, and `seeding_age` is diagnostic only.
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
    ## One growth source. The prior is on the cryptic exponential growth rate
    ## `r`, and the established reproduction number `R0` (the walk base) is
    ## derived forward from it through Euler–Lotka.
    growth_state ~ to_submodel(growth(g))
    r_clock = growth_state.r
    R0 = r_to_R0(r_clock, g)
    ## The random walk's first knot sits at `rt_walk_start`, decoupled from
    ## the renewal start. The renewal seeds and grows from the genetic-TMRCA
    ## renewal start, but `R_t` is held flat at `R0` until the first
    ## situation report, because before any case or death surveillance the
    ## dynamics are unidentified and a free walk there only adds unsupported
    ## drift. `rt_walk_start` defaults to `rt_start`.
    rt_state ~ to_submodel(rt(n, log(R0); breakpoint, rt_start = rt_walk_start))
    Rt = rt_state.Rt
    ## The renewal-start seed is the daily incidence `C_T = exp(r·T)` reached
    ## after the cryptic phase's `m` generations. Grid days
    ## `1…renewal_start` are filled with the cryptic exponential curve at
    ## rate `r` ending at that seed, a full generation interval of history,
    ## and the renewal then runs forward from `renewal_start+1`.
    renewal_start = clamp(rt_start, 1, n)
    τ_obs = n - renewal_start
    seed0 = seed_at_renewal_start(growth_state.C_T)
    seed_vec = seed_infections(seed0, r_clock, renewal_start)
    infections = renewal_infections(Rt, g, seed_vec)
    cumulative = cumsum(infections)
    ## Total outbreak age: cryptic duration (m generations) plus the span.
    T_total = growth_state.T + τ_obs
    ## Current growth rate at the cut-off, derived from the cut-off
    ## reproduction number `Rt[n]` and the generation interval through
    ## forward Euler–Lotka, the inverse of the `r_to_R0` above. This makes
    ## the reported growth rate consistent with `R_T := Rt[n]` by
    ## construction, so `r < 0` iff `R_T < 1`. The realised last-two-days
    ## slope is not used: the intervention ramp depresses the final renewal
    ## step, so that slope can disagree in sign with `R_T`.
    r = euler_lotka_r(@inbounds(Rt[n]), g)
    return (; infections, cumulative, Rt, g, seed_at_renewal_start = seed0,
        m = growth_state.m, τ = growth_state.τ, R0, r0 = r_clock, r,
        doubling_time_initial = doubling_time(r_clock),
        T = T_total, C_T = cumulative[n],
        C_T_prior = growth_state.C_T, doubling_time = doubling_time(r),
        seeding_age = seeding_age(cumulative, n))
end

"""
Onset-incidence submodel. Convolves the renewal infections with the sampled
incubation-period PMF to get daily symptom-onset incidence, computed once
per draw and reused by every downstream observation stream. The incubation
delay submodel is injected. The incubation period cannot be fitted from the
BDBV line list (no exposure dates), so it follows the MacNeil et al. (2010)
Bundibugyo estimate from the 2007 Uganda outbreak, mean 6.3 d (95% CI
5.2-7.3, n = 24). The mean prior `Normal(6.3, 0.54)` reproduces that 95% CI
(SD = CI half-width / 1.96). MacNeil give no interval on the spread, so the
SD prior is a weakly-informative choice centred on the CV-implied spread
(≈ 3.5 d). Returns
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
[`infection_model`](@ref)). The TMRCA is treated as a right-censored, noisy
reading of the seeding time, so deeper or wider sampling only pushes it
older. The likelihood contributes `P(read ≥ tmrca_days)`.
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

The default `truncated(Normal(0, 0.25); lower = 0)` is informative,
mirroring the case background. Deaths are far fewer than suspected cases
(≈ 246 suspected deaths against ≈ 1077 suspected cases at the last stable
vintages), so the background rate is scaled down accordingly. With SD 0.25
the median background is ≈ 0.17/day, a modest minority of the
suspected-death total over the grid, while still admitting a genuine
non-BVD signal. The background is degenerate with outbreak size, so a
diffuse prior would let it absorb arbitrarily many suspected deaths. Pass
`lambda_prior` to override. Returns `(; λ_bg_death)`.
"""
@model function death_background_model(;
        lambda_prior = truncated(Normal(0.0, 0.25); lower = 0))
    λ_bg_death ~ lambda_prior
    return (; λ_bg_death)
end

"""
Ascertainment of the suspected-death stream ([`deaths_model`](@ref)), the
fraction `p_death` of true BVD deaths that enter the INSP suspected-death
count by the cut-off. The suspected-death definition is
symptomatic-then-deceased, so a fatal BVD infection only counts once the
death is reported, and not every BVD death is captured. The expected BVD
suspected deaths are `p_death · CFR` of the onset-to-death-convolved
infections, the death analogue of the suspected-case ascertainment `p_drc`
([`pooled_ascertainment_model`](@ref)).

The default `Normal(logit(0.9), 0.5)` on the logit scale is informative and
centred on a high ascertainment, since a death is a salient event in an
Ebola response and is reported more reliably than a living suspected case
(`p_drc ≈ 0.75`). The SD 0.5 gives a 90% prior interval of roughly
0.80–0.95. `p_death` is weakly identified on its own, trading off with the
CFR for the suspected-death level, so it leans on this prior. The
export-death stream and the CFR prior pin the CFR separately. Pass
`ascertainment_prior` to override. Returns `(; p_death, logit_p_death)`.
"""
@model function death_ascertainment_model(;
        ascertainment_prior = Normal(logit(0.9), 0.5))
    logit_p_death ~ ascertainment_prior
    p_death := logistic(logit_p_death)
    return (; p_death, logit_p_death)
end

"""
Background case-fatality ratio `cfr_bg` for the non-BVD suspected-death
background ([`deaths_model`](@ref)). The per-day background suspected
deaths are `cfr_bg · λ_bg_v`, a background CFR applied to the per-day
non-BVD suspected-case rate ([`test_positivity_model`](@ref)). A non-BVD
suspected case (other severe febrile or haemorrhagic illness that meets the
suspect definition) carries its own fatality risk, and `cfr_bg` is the
share of that background pool reported as a suspected death.

Tying the death background to the case background rather than giving deaths
a free rate of their own removes the degeneracy that keeps a free
`λ_bg_death` switched off. The case background is pinned by the laboratory
positivity link, so scaling it by `cfr_bg` gives the death background a
level and time profile without a second free rate competing with outbreak
size. The default `Beta(2, 6)` (mean ≈ 0.25, 90% ≈ 0.05–0.52) is weakly
informative and centred below the BVD CFR, since non-BVD suspect illness is
on average less lethal. Pass `cfr_prior` to override. Returns `(; cfr_bg)`.
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
seeding-to-cut-off span. The prior is informative because `λ_bg` is
degenerate with outbreak size (the per-vintage reported mean mixes the
`p_drc`-scaled BVD increment with `λ_bg · Δt`), and a diffuse prior lets the
background absorb arbitrarily many suspected cases and opens a second
posterior mode in which it explains the majority of them. With SD 1.0 the
median background is ≈ 0.67/day and the 95% prior bound ≈ 2.0/day, a modest
minority of the ≈ 1077 suspected cases observed by the last stable
suspected-case vintage while still admitting a genuine non-BVD signal. Pass
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
length-of-stay mean for the occupancy level (Little's law, mean occupancy
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
identified and the half-normal `truncated(Normal(0, 0.75); lower = 0)`
carries most of the weight. `δ_iso = 0` recovers a shared admission rate.
Pass `logodds_prior` to override. Returns `(; δ_iso)`.
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

The default `LogNormal(log 450, 0.42)` is weakly informative and positive
by construction, with median 450 and a ≈0.44 coefficient of variation,
centred on the bed count implied by the reported occupancy rates (the
"Taux d'occupation" gives `capacity = occupancy / rate ≈ 400–452` over
9–13 June). The capacity is identified by the implied-capacity series the
isolation submodel fits, so the prior only has to bracket it. A single
national capacity averages over a growing capacity and cannot represent
local saturation, one province full while another has slack. See
[`bed_capacity_walk_model`](@ref) for the time-varying form. Pass
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
count `C0` on weekly knots, linearly interpolated to the daily grid, the
same parameterisation as the reproduction-number and background walks. With
knot values `\\log C` and knot days `d`,
`C(t) = C0 · exp(\\text{interp}(σ_cap · cumsum(z)))` with `z ~ Normal(0, 1)`
per knot and a tight innovation SD `σ_cap`, keeping capacity a gentle
drift rather than per-day jumps. Knots need far fewer innovations than a
daily walk, avoiding the high-dimensional funnel. The baseline carries
the same weakly-informative `LogNormal(log 450, 0.42)` prior as the
scalar model (median 450 beds, ≈0.44 CV), so `C0` is sampled on the log
scale and the whole capacity `log C(t) = log C0 + walk` is fully
log-scale. The implied-capacity series the isolation submodel fits pins
`C(t)` on the days a rate is published.

Knots run only from `start`, the first day with occupancy or capacity data,
and capacity is flat at `C0` before it. Off-window capacity carries no
likelihood, so walking it there would add unidentified innovations. Pass
`start = 1` for knots over the whole grid, or `week` to change the knot
spacing. A single national capacity cannot represent local saturation, one
province full while another has slack. Pass
`baseline_prior` / `innovation_prior` to override. Returns
`(; C, C0, σ_cap)` with `C` a length-`n` vector.
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
    ## Non-negative innovations, so capacity is non-decreasing. Beds are
    ## added over the response and not taken away, so `C(t)` cannot drop
    ## below an already-reached level, and the effective ceiling cannot
    ## jitter down into the observed occupancy.
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
fraction is exactly `1 − CFR`. The offset keeps `p_recover` in `(0, 1)`
without a hard clamp. `p_recover` is partially confounded with the
confirmation-to-recovery delay for the count of recoveries observed by the
cut-off, since a long delay right-censors recoveries that have not yet
resolved, so the delay carries the timing and `p_recover` the eventual
survival fraction. Pass `offset_prior` to override. Returns
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
suspected-case and suspected-death streams share one. See
[`background_pooling_model`](@ref), which samples it once at the composer
level. The deviation is multiplicative, so the per-vintage rate stays
positive without a clamp and `σ_bg → 0` recovers the scalar baseline
exactly. Each stream still samples its own baseline `λ_mu` and per-vintage
deviations `z`. `nv` is the number of vintage windows. Returns
`(; λ, λ_mu, σ_bg, z)` with `λ` a length-`nv` vector of per-vintage rates.
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
own from few vintages. The prior is a half-normal of scale 0.3.

The background is degenerate with outbreak size, so the prior still
regularises rather than freeing the scale. A wide random effect would let
individual windows absorb arbitrary suspected counts and re-open the
posterior mode in which the background explains the majority of suspected
cases, and regularising `σ_bg` toward zero keeps the time variation a
perturbation of the informative scalar baselines.

The scale was 0.1, which the data have outgrown. With the daily
new-suspect series carried to the cut-off the posterior sits at 0.17 to
0.22, about twice that scale, and the walk mixes badly against a prior
pulling the other way. At 0.3 the same posterior sits below the scale, so
the data set the time variation rather than the prior. Returns
`(; σ_bg)`.
"""
@model function background_pooling_model(;
        pooling_prior = truncated(Normal(0.0, 0.3); lower = 0))
    σ_bg ~ pooling_prior
    return (; σ_bg)
end

"""
Non-BVD background rate as a smooth weekly lognormal random walk over the
surveillance window, the alternative to the per-vintage step random
effect ([`background_re_model`](@ref)). The log-rate follows a Gaussian
random walk on weekly knots and is linearly interpolated to the daily
grid, the same knot layout as the reproduction-number walk
([`rt_walk_model`](@ref)). The background is a slow drift, so a knot per
`week` carries the time variation with far fewer innovations than a daily
walk. The series is gated to zero before the surveillance `onset`, since
the non-BVD background does not exist before surveillance began, and ramps
in over the first `onset_ramp` days of the window. With knot values
`\\log\\lambda` and knot days `d`,

```math
\\log \\lambda_d = \\log \\lambda_0 + \\sum_{s < d} \\epsilon_s,
\\qquad \\epsilon_s \\sim \\mathcal N(0, \\sigma_{rw}),
\\qquad \\lambda_t = 0 \\ \\text{for}\\ t < \\text{onset}.
```

`σ_rw` is the per-knot innovation SD on the log scale, passed in and shared
across the suspected-case and suspected-death streams via
[`background_pooling_model`](@ref). Its regularising prior keeps the
background a slow drift, which holds down the background/outbreak-size
degeneracy and keeps the series smooth, so a death background scaled from
it carries no steps.

The default is the centred form, each knot step drawn directly at
`Normal(0, σ_rw)`. The daily new-suspect series runs to the cut-off, so
the walk is strongly informed, and the non-centred form
`steps = σ_rw .* z` funnels, with `z` diverging as `σ_rw → 0` and
stretching NUTS trajectories. That funnel is what broke the joint fit
when the series resumed. The worst R-hat went to 1.6 with 5 bulk
effective samples, against 1.05 and 65 on the vintage before, and the
last walk innovation was among the worst-mixing parameters. Pass
`centred = false` for the non-centred form, the better choice when the
walk is weakly informed and prior-dominated. Both forms carry the same
prior, a cumulative sum of `Normal(0, σ_rw)` steps, so only the sampled
coordinates differ.
[`pooled_dispersion_model`](@ref) carries the same switch for the same
reason.

Knots run only over the surveillance window `[onset, n]`, so the number of
innovations is small. `onset ≤ 1` runs it over the whole grid. Pass
`week` to change the knot spacing.

`λ_mu ~ truncated(Normal(0, 20); lower = 0)` is the scale the walk
multiplies, not the level over the window. The log-deviation is pinned to
zero on the first knot, so `λ_mu` anchors the window's start and the
innovations carry the series from there. The half-normal shrinks that
anchor toward zero, which stops the background out-explaining the outbreak
signal, and the scale is wide enough not to truncate the anchor the
suspected-case data support. Pass `baseline_prior` to override.

Returns `(; λ, λ_mu, σ_bg)` with `λ` the length-`n` daily series (zero
before `onset`).
"""
@model function background_walk_model(n::Integer, σ_rw::Real;
        onset::Integer = 1, onset_ramp::Integer = 7, week::Integer = 7,
        baseline_prior = truncated(Normal(0.0, 20.0); lower = 0),
        centred::Bool = true)
    t0 = clamp(Int(onset), 1, n)
    nw = n - t0 + 1
    ## Weekly knots over the window, linearly interpolated to the daily grid
    ## (see [`knot_days`](@ref) and [`interpolate_knots`](@ref)).
    days = knot_days(n; week = week, start = t0)
    nb = length(days)
    ## Half-normal rather than lognormal. A log-scale level has a heavy right
    ## tail the background/outbreak-size degeneracy exploits to run away.
    λ_mu ~ baseline_prior
    m = max(nb - 1, 1)
    if centred
        ## Draw each knot step on its own scale. `eps` floors the SD so a
        ## `σ_rw ≈ 0` draw stays a proper distribution.
        steps ~ product_distribution(
            fill(Normal(0, σ_rw + eps(typeof(float(σ_rw)))), m))
        walk_steps = steps[1:max(nb - 1, 0)]
    else
        z ~ product_distribution(fill(Normal(0, 1), m))
        walk_steps = σ_rw .* z[1:max(nb - 1, 0)]
    end
    ## Smooth multiplicative deviation, a cumulative log-deviation from the
    ## baseline, anchored there on the first knot. Interpolated to daily so a
    ## death background scaled from it is smooth.
    log_knots = vcat(zero(σ_rw), cumsum(walk_steps))
    walk = interpolate_knots(log_knots, days, n)[t0:n]
    λ_window = λ_mu .* exp.(walk)
    ## Linear onset ramp `0 → 1` over the first `onset_ramp` days of the
    ## window, so the gated background grows in from zero rather than stepping
    ## straight to `λ_mu` at the surveillance boundary and putting a one-day
    ## jump into the suspected-death trajectory scaled from it.
    ## `onset_ramp ≤ 1` gives a hard onset.
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
Filovirus Screen RT-PCR [rieger2016](@cite), which detects Bundibugyo
virus at 11–67 RNA copies per reaction. The Zaire-specific GeneXpert Ebola
assay does not reliably detect Bundibugyo
[cepheid_xpert_ebola_ifu, pinsky2015, semper2016](@cite). A single assay
draw is sensitive to about 0.85, but a suspect is confirmed or ruled out
through repeat control tests rather than one PCR, so the effective process
sensitivity is higher (two controls give about 0.98). Under the
severe-first backlog the first vintage's analysed batch is near-pure BVD
(`q ≈ 1`), so the v1 positivity ≈ `s` identifies the sensitivity from the
early data. Returns `(; s_test)`.
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
`p = s · q + (1 − spec)(1 − q)` with `q` the tested BVD share. The
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
[`confirmed_cases_model`](@ref)). The per-window positivity `p_pos` is a
smooth curve that does not capture the day-to-day laboratory batching and
within-window positivity heterogeneity the confirmed counts carry, and a
plain `Binomial` on denominators of several hundred specimens gives
predictive intervals far too tight. The intra-window correlation
`ρ ∈ (0, 1)` inflates the window variance to
`n·p·(1 − p)·(1 + (n − 1)·ρ)`, with `ρ → 0` recovering the `Binomial`.
The default `Beta(1, 24)` (mean ≈ 0.04, 90% ≈ 0.002–0.12) is weakly
informative and shrinks toward the `Binomial` when the data support it. One
scalar is identified across the laboratory windows, so the confirmed
positives themselves set the spread. Returns `(; ρ)`.
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
priors are kept tight around the documented turnaround belief (mean
≈ 4.5 d, SD ≈ 4 d), since a wide prior on an unidentified nuisance delay
only makes the sampler wander it, dragging the confirmation PMFs convolved
from it. Returns `(; pmf, dist, mean, sd)`.
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
positivity (≈ 0.28, that is 210 / 755 on the 28 May data) on the logit
scale. Conditioning on the observed denominator and giving the positivity
its own random effect decouples the confirmed counts from the
multiplicative ascertainment ridge `p_drc · s_test · τ_test`, so the
outbreak size is pinned by the deaths and exports streams rather than
forced through the laboratory positivity. Returns `(; p_pos, q_mu, σ_q)`
with `p_pos` a length-`nv` vector.
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

with `c_v` the cumulative analysed volume at window `v`, the testing clock.
The lab over-tests BVD early, since severe cases are triaged first and are
more likely BVD, and the enrichment `δ₀·e^{−c/decay}` relaxes toward zero as
testing widens, at which point the tested share equals the pool
composition. This ties positivity to the background `λ_bg`, so the
confirmed and positivity data identify the non-BVD background rather than
it being absorbed by a free per-window random effect.

`δ₀` is the early severity log-odds enrichment of BVD, lower-truncated at 0
because severity triage upsamples BVD, never down. The default
`truncated(Normal(1.5, 0.75); lower = 0)` is moderate and bounded, since
even severity-triaged testing cannot be near-pure BVD (other haemorrhagic
or severe febrile illness is also triaged). For a pool composition
`φ ≈ 0.4` the early tested share is `logistic(logit(0.4) + 1.5) ≈ 0.75`.
`decay_scale` is the relaxation timescale on the analysed-volume clock.
Pass `logodds_prior` / `decay_prior` to override. Used by
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
([`confirmed_deaths_only_model`](@ref)), which has no case stream to set
the death testing volume from. It thins the suspected deaths to a death
"analysed" volume at the case testing rate, drawing `τ_death` from the same
prior as the case testing fraction (`Beta(5, 2)`, mean ≈ 0.71). The full
joint instead scales the modelled case analysed volume (see
[`confirmed_deaths_model`](@ref) and [`death_testing_scaling_model`](@ref))
and does not draw this submodel. Pass `fraction_prior` to override. Returns
`(; τ_death)`.
"""
@model function death_testing_fraction_model(; fraction_prior = Beta(5.0, 2.0))
    τ_death ~ fraction_prior
    return (; τ_death)
end

"""
Death testing-intensity scaling for the confirmed-death volume in the joint
([`confirmed_deaths_model`](@ref)). The death analysed volume is the
modelled case analysed volume carried at the per-day suspected
death-to-case ratio, times this scaling. That ratio already carries the
suspect-pool severity and the suspected-death level, so the scaling is the
per-suspect testing-intensity difference between deaths and living suspects
alone. No death-testing data grounds it, so it is a tight log-normal
centred on one (`LogNormal(0, 0.25)`, median 1, 90% ≈ 0.66–1.51). Deaths
are tested at the case intensity unless the confirmed-death counts pull the
scaling off one. Pass `scaling_prior` to override. Returns `(; scaling)`.
"""
@model function death_testing_scaling_model(;
        scaling_prior = LogNormal(0.0, 0.25))
    scaling ~ scaling_prior
    return (; scaling)
end

"""
Specimens analysed per suspect sampled.

[`confirmed_cases_model`](@ref) scales the laboratory volume by this factor
as well as by `τ_test`. `τ_test` is a probability and the receipt kernel
conserves mass, so it alone bounds the modelled analysed volume below the
modelled suspect inflow. Specimens are not persons: a suspect can yield
several through repeat exclusion testing, and swabbed community deaths and
screened contacts enter the laboratory denominator without being counted
as suspects reported, so the ratio can exceed one.

`κ ~ LogNormal(0, 0.25)` has median 1 and a 90% range of about 0.66 to
1.51, matching [`death_testing_scaling_model`](@ref).
"""
@model function specimen_intensity_model(;
        intensity_prior = LogNormal(0.0, 0.25))
    κ ~ intensity_prior
    return (; κ)
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
handful of deaths against hundreds of suspects against a daily laboratory
volume) do not share one global `k` that the dominant stream pulls around,
while the sparse streams still borrow strength through the common
hyper-parameters rather than going noisy on a fully independent draw.

The pooling is on the `log(1/sqrt(k))` scale, with a population mean
`μ_log`, a pooling SD `τ`, and per-stream deviations, and
`k_s = 1 / inv_sqrt_k_s^2`. The default is the centred form,
`log(1/sqrt(k))_s ~ Normal(μ_log, τ)` drawn directly. The
passive-surveillance streams are data-rich, so each stream's dispersion is
strongly informed and the non-centred form
`inv_sqrt_k_s = exp(μ_log + τ z_s)` funnels, with
`z_s = (log_isk_s − μ_log)/τ` diverging as `τ → 0` and stretching NUTS
trajectories. On the joint, centring removes that funnel: worst dispersion
bulk-ESS 102 → 156, divergences 5 → 2, and about 10% faster wall-clock at a
150×2 fit. Pass `centred = false` for the non-centred form, the better
choice when the streams are data-poor and prior-dominated. The population
mean is centred on the shared `1/sqrt(k)` prior of
[`surveillance_dispersion_model`](@ref) (`exp(μ_log)` near 0.6), and the
half-normal `τ` keeps the per-stream dispersions close unless the data pull
them apart. `τ = 0` collapses every stream to the population value, the
shared-`k` model. Returns `(; k, inv_sqrt_k, k_pop, μ_log, τ)` with `k` a
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
systems. The two countries run different systems, DRC passive community
surveillance and Uganda point-of-entry or hospital detection, so each
ascertainment fraction has its own logit-scale prior with no shared
parameter. An alternative to the composer-default
[`pooled_ascertainment_model`](@ref) for sensitivity analyses that do not
share strength between the two systems.

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
to vary in space and in time, drawn from a multivariate-normal random walk
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

`μ(t)` is the national weekly-knot walk ([`rt_walk_model`](@ref)), kept
intact, so the national streams see exactly the `Rt` process the headline
model fits and it is the target the provinces pool toward. The provinces
are not equally observed: over the fitted window Nord-Kivu contributes 74
laboratory positives and Sud-Kivu contributes none at all. Shrinking toward
a common trend lets them borrow strength from Ituri and deviate only where
the data insist. With `μ(t)` present the deviation covariance `Σ` is still
free, so the cross-patch correlation is learned rather than assumed.

### Sum-to-zero, not a reference patch

The deviations are centred at every knot, so no province is privileged.
Fixing `δ_1 ≡ 0` would also identify the model, but it forces the primary
patch to have no idiosyncratic deviation at all. Centring introduces one
redundant coordinate per knot (the mean of the raw innovations, which the
likelihood never sees), drawn from its proper prior and well-conditioned,
so it costs a few cheap sampled dimensions rather than a posterior ridge.

### What the data can and cannot identify here

The composition of the confirmed cases identifies the contrast between
provinces. With three patches that is essentially one number, the Ituri /
Nord-Kivu contrast, since Sud-Kivu carries no signal. Expect `Ω` to be
largely prior-driven and Sud-Kivu's `Rt` to be pinned by the deviation
prior rather than by data, which is why `Σ` is given a proper shrinkage
prior rather than a flat one.

`σ_δ → 0` recovers a common `Rt` shape shared by every province, a fixed
ratio between them. It is a special case of this model rather than an
assumption baked into it. `σ_δ` is therefore the headline spatial
diagnostic, and a posterior pushed away from zero is evidence that
provincial `Rt` trajectories are separating, which is what a response
concentrated on the Ituri epicentre would produce.

### Mean reversion, not a random walk

The deviations mean-revert to zero rather than random-walk. A random walk
has no mean, so a province that happens to sit above the national trend at
the last vintage is projected to stay above it for ever with the gap as
likely to widen as to close. That matters here because the per-province
vintages stop well before the cut-off. The knots therefore follow

```math
\\boldsymbol{\\delta}(t_k) = \\phi\\, \\boldsymbol{\\delta}(t_{k-1})
    + \\boldsymbol{\\eta}_k, \\qquad
\\phi = 2^{-\\text{week} / h},
```

with `h` the sampled half-life of a provincial divergence in days, so a
province with persistent divergence can still show one. The innovations are
centred at every knot and `φ` is one shared scalar, so the sum-to-zero
constraint survives exactly. A half-life far longer than the window
recovers the random walk.

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
        region_halflife_prior = LogNormal(log(42), 0.6),
        lkj_prior = LKJCholesky(max(n_patches, 2), 2.0),
        region_offset_prior = Normal(0, 1))
    ## Common national trend, the single-patch walk unchanged.
    ## `rt_walk_start` maps to `rt_start` in the inner model, matching the
    ## convention in [`infection_model`](@ref). Attached prefixed (no
    ## `false`), so the walk's parameters reach the chain as
    ## `rt_state.sigma_rw`, `rt_state.log_R0`, `rt_state.z` and
    ## `rt_state.intervention_effect`, the names the analysis and sensitivity
    ## pages read. Attaching it unprefixed surfaces them bare and fails at
    ## render time on a KeyError.
    rt_state ~ to_submodel(
        rt(n, log_R0_base; breakpoint, rt_start = rt_walk_start))
    Rt_national = rt_state.Rt
    log_Rt_national = log.(Rt_national)
    ## The deviations live on the same weekly knots as the national walk, so
    ## both processes are described at the same resolution.
    days = knot_days(n; week, start = rt_walk_start)
    nb = length(days)
    ## Single patch. The deviations are sum-to-zero across the patches, so
    ## with one patch delta is identically zero and the patch Rt is the
    ## national walk. Sampling the deviation machinery would then add
    ## prior-only dimensions the likelihood never touches, so it is skipped
    ## entirely and `n_patches = 1` collapses this model exactly onto the
    ## single-population one.
    if n_patches == 1
        Tp1 = eltype(Rt_national)
        δ_patch1 = zeros(Tp1, 1, n)
        Rt_matrix1 = zeros(Tp1, 1, n)
        @inbounds for t in 1:n
            Rt_matrix1[1, t] = Rt_national[t]
        end
        return (; Rt_matrix = Rt_matrix1, Rt_national, log_Rt_national,
            δ_patch = δ_patch1, δ_knots = zeros(Tp1, 1, nb),
            σ_level = zero(Tp1), σ_δ = zeros(Tp1, 1),
            Ω = ones(Tp1, 1, 1), δ_halflife = zero(Tp1),
            sigma_rw = rt_state.sigma_rw,
            log_R0 = rt_state.log_R0,
            intervention_effect = rt_state.intervention_effect)
    end
    ## Deviation scales (one per patch) and their cross-patch correlation.
    ## `LKJCholesky` samples the Cholesky factor directly, so the
    ## decomposition never lands on the AD tape.
    σ_level ~ region_sd_prior
    σ_δ ~ product_distribution(fill(region_drift_sd_prior, n_patches))
    ## Mean reversion. The deviations are an AR(1) toward zero on the knots,
    ## parameterised by the half-life of a provincial divergence in days,
    ## which is the elicitable quantity. The per-knot retention is
    ## `phi = 2^(-week / halflife)`, so a half-life far longer than the
    ## window recovers the random walk and a short one pulls each province
    ## back to the national trend between knots. One half-life is shared
    ## across provinces, not one each: the retention multiplies the whole
    ## deviation vector, and a centred vector scaled by a scalar is still
    ## centred, so the sum-to-zero constraint survives exactly.
    δ_halflife ~ region_halflife_prior
    φ = exp2(-week / δ_halflife)
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
            δ_knots[i, k] = φ * δ_knots[i, k - 1] + (innov[i] - innov_bar)
        end
    end
    ## Interpolate each patch's deviation to the daily grid and build Rt.
    δ_patch = zeros(Tp, n_patches, n)
    Rt_matrix = zeros(Tp, n_patches, n)
    @inbounds for p in 1:n_patches
        ## A view, not a copy: `interpolate_knots` only reads its knots, and
        ## the copy put one `getindex` per knot on the gradient tape.
        δ_daily = interpolate_knots(view(δ_knots, p, :), days, n)
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
        σ_level, σ_δ, Ω, δ_halflife, sigma_rw = rt_state.sigma_rw,
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

The default `importation_kernel` is the gravity kernel of
[`province_importation_kernel`](@ref), a fixed weighting by destination
population, so the provinces are coupled and the intensity `ε` is sampled.
There is no mobility or origin-destination data for this outbreak, so the
kernel is a structural assumption, and `ε` is weakly identified against the
secondary-patch seeds, since both raise a secondary province's early
incidence. Read `ε` as the scale of coupling the data will tolerate rather
than as a measured flow.

Passing an all-zero kernel uncouples the provinces. `ε` is then not
sampled, since against a zero kernel it would be a dimension the likelihood
never touches, and each secondary patch is explained by its own seed and
its own `R_t`.

### Seeding

The outbreak is assumed to have begun in Ituri, so the primary patch carries
the whole cryptic seed, growing at the sampled molecular-clock rate `r` over
the cryptic window to reach `seed_at_renewal_start(C_T)` at the renewal
start. The
secondary patches start empty and are seeded by importation from it, so a
province's arrival is a consequence of the kernel and of `ε` rather than a
parameter. The first-appearance dates carry real information about how fast
the outbreak spread between provinces, and a free seed fraction per
province would absorb exactly that, since both a larger seed and a stronger
coupling raise a secondary province's early incidence.

An all-zero kernel leaves a secondary patch no route to infections at all,
so the uncoupled path keeps the sampled fractions (`seed_fraction_prior`, a
`LogNormal` on the fraction of the primary seed). They partition the
national cryptic seed rather than adding to it, so `2^m` stays the
country's cryptic size for any patch count and `C_T` stays comparable
across them.

### Returns

The per-patch state, plus the national aggregates the observation models
and the headline summaries consume: `infections_total`, `cumulative_total`,
and `C_T` (the national cut-off cumulative). `R_T`, `r`, `T` and
`doubling_time` mirror [`infection_model`](@ref) so a patch chain carries
the same headline quantities as a single-patch one. `R_T` is the
incidence-weighted aggregate reproduction number at the cut-off, obtained
by inverting the renewal equation on the summed infections
([`implied_national_Rt_at`](@ref)) on the cut-off day alone.
`importation_matrix` is the daily infections each province received from
the others, which is what the imports figure on the analysis page draws.
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
        importation_sd_prior = truncated(Normal(0, 0.5); lower = 0),
        importation_effect_prior = Normal(0, 0.5),
        seed_fraction_prior = LogNormal(log(0.05), 1.0),
        incubation = (nmax) -> censored_delay_model(nmax;
            mean_prior = truncated(Normal(6.3, 0.54); lower = 1),
            sd_prior = truncated(Normal(3.5, 0.8); lower = 1)),
        incubation_nmax::Integer = cdf_nmax(lognormal_meansd(6.3, 3.5)))
    ## 1. Shared generation interval.
    gi_state ~ to_submodel(gi(gi_nmax))
    g = gi_state.g
    ## 2. One growth source, as in [`infection_model`](@ref). The prior is on
    ##    the cryptic growth rate `r`, and the established `R0` (the walk
    ##    base) is derived forward from it through Euler-Lotka.
    growth_state ~ to_submodel(growth(g))
    r_clock = growth_state.r
    R0 = r_to_R0(r_clock, g)
    ## 3. Per-patch Rt: national trend plus per-patch deviations.
    rt_state ~ to_submodel(
        rt(n, n_patches, log(R0); breakpoint, rt_start, rt_walk_start), false)
    Rt_matrix = rt_state.Rt_matrix
    δ_patch = rt_state.δ_patch
    ## 4. Per-patch seeds. The primary patch takes the cryptic exponential.
    ##    Each secondary patch takes a fraction of that seed, which is the
    ##    scale the data speak to. With importation off the relative seed
    ##    sets the level of the provincial case split, leaving `δ_p` to be
    ##    identified by its time trend. An absolute seed prior pinned far
    ##    below the primary's `2^m` would force `δ_p` to absorb the whole
    ##    level difference, making the reported provincial Rt gap an artefact
    ##    of the seed prior.
    renewal_start = clamp(rt_start, 1, n)
    τ_obs = n - renewal_start
    seed0_total = seed_at_renewal_start(growth_state.C_T)
    ## The outbreak is assumed to have begun in Ituri, so the primary patch
    ## takes the whole cryptic seed and the others are seeded by importation
    ## from it. An
    ## all-zero kernel leaves a secondary patch no route to infections at
    ## all, so the uncoupled path keeps the sampled fractions. With one patch
    ## there is nothing to seed and the fraction would be a prior-only
    ## dimension either way.
    coupled = any(!iszero, importation_kernel)
    if n_patches > 1 && !coupled
        seed_fraction ~ product_distribution(
            fill(seed_fraction_prior, n_patches - 1))
    else
        seed_fraction = Float64[]
    end
    Tp = promote_type(eltype(Rt_matrix), eltype(g), typeof(float(r_clock)),
        eltype(seed_fraction), typeof(float(seed0_total)))
    ## The fractions partition the national cryptic seed, they do not add to
    ## it. `growth_state.C_T` is `2^m`, and the `m` prior is elicited as a
    ## national quantity, so it is the size of the whole cryptic phase.
    ## Dividing through by `(1 + Σf)` keeps the national seed at `2^m` for any
    ## number of patches, so `C_T` stays comparable across `n_patches` and the
    ## genetic prior keeps its meaning.
    seed_shares = zeros(Tp, n_patches)
    if isempty(seed_fraction)
        seed_shares[1] = one(Tp)
    else
        seed_denom = one(Tp) + sum(seed_fraction)
        seed_shares[1] = one(Tp) / seed_denom
        @inbounds for p in 2:n_patches
            seed_shares[p] = seed_fraction[p - 1] / seed_denom
        end
    end
    seeds_matrix = zeros(Tp, n_patches, renewal_start)
    @inbounds for p in 1:n_patches
        ## A scaled copy of the same cryptic curve: each province is a share
        ## of one epidemic, so it grows at the same clock rate `r` over the
        ## cryptic window.
        s_p = seed_infections(
            seed_shares[p] * seed0_total, r_clock, renewal_start)
        for j in 1:renewal_start
            seeds_matrix[p, j] = s_p[j]
        end
    end
    ## 5. Importation intensity: one level per origin, partially pooled, and
    ##    a common change at detection on the ramp the reproduction number
    ##    already uses. Only sampled when the kernel couples the patches,
    ##    since against an all-zero kernel it would be a prior-only dimension.
    ##    Per origin because the provinces do not export alike and the kernel
    ##    only carries size and distance. Pooled because the small provinces
    ##    export too little for their own level to be identified, so
    ##    `σ_ε → 0` recovers one shared intensity and a province the data say
    ##    nothing about sits at the pooled mean. The deviations are centred,
    ##    so `ε_bar` stays the overall level. Time-varying because the
    ##    outbreak being known changes movement, and the provinces that arrive
    ##    either side of the breakpoint are what separates `β_ε`.
    ε_matrix = zeros(Tp, n_patches, n)
    if coupled
        ε_bar ~ importation_epsilon_prior
        σ_ε ~ importation_sd_prior
        z_ε ~ product_distribution(fill(Normal(0, 1), n_patches))
        β_ε ~ importation_effect_prior
        z_bar = sum(z_ε) / n_patches
        ramp = sigmoid_ramp(n, breakpoint)
        @inbounds for q in 1:n_patches
            lvl = ε_bar * exp(σ_ε * (z_ε[q] - z_bar))
            for t in 1:n
                ## Capped at one: the origin cannot send away more than it
                ## generates. The prior sits four orders of magnitude below
                ## the cap, so this binds only in the far tail.
                ε_matrix[q, t] = min(lvl * exp(β_ε * ramp[t]), one(Tp))
            end
        end
        importation_epsilon := ε_bar
        importation_epsilon_sd := σ_ε
        importation_epsilon_effect := β_ε
        importation_epsilon_patch := [ε_matrix[q, n] for q in 1:n_patches]
    end
    ## 6. Multi-patch renewal. Each province runs its own renewal at its own
    ##    reproduction number and the national trajectory is their sum. There
    ##    is no separate national process and nothing rescales the patches to
    ##    match one. `mu(t)` is a central trend the provinces pool toward, and
    ##    the reproduction number the country actually ran at is read back off
    ##    the summed infections in step 9.
    renewal_state = patch_infections(Rt_matrix, g, seeds_matrix,
        importation_kernel, ε_matrix)
    infections_matrix = renewal_state.infections
    importation_matrix = renewal_state.importation
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
    ## 8. Per-patch onsets through the shared incubation PMF.
    inc_state ~ to_submodel(incubation(incubation_nmax))
    onsets_matrix = zeros(Tp, n_patches, n)
    @inbounds for p in 1:n_patches
        @views onsets_matrix[p, :] = convolve_delay(
            infections_matrix[p, :], inc_state.pmf)
    end
    ## 9. Aggregate reproduction number. Inverting the renewal equation on
    ##    the summed infections gives the incidence-weighted mean of the
    ##    patch `Rt`s. With `I_{p,t} = R_{p,t} · force_{p,t}`, summing over
    ##    patches gives `I_t / Σ_p force_{p,t} = Σ_p R_{p,t} force_{p,t} /
    ##    Σ_p force_{p,t}`. This is the `Rt` that reproduces the national
    ##    trajectory, so it is the one the headline `R_T` reports. Only the
    ##    cut-off day is inverted, since that is the only day reported.
    R_T = implied_national_Rt_at(infections_total, g, n)
    ## 10. Headline quantities, mirroring [`infection_model`](@ref) so a
    ##     patch chain summarises exactly like a single-patch one.
    r = euler_lotka_r(R_T, g)
    T_total = growth_state.T + τ_obs
    return (; infections_matrix, cumulative_matrix, onsets_matrix,
        Rt_matrix, importation_matrix,
        δ_patch, δ_knots = rt_state.δ_knots,
        C_T_patch,
        σ_level = rt_state.σ_level,
        σ_δ = rt_state.σ_δ,
        δ_halflife = rt_state.δ_halflife,
        Ω = rt_state.Ω,
        infections_total, cumulative_total,
        Rt_national = rt_state.Rt_national,
        g, R0, r0 = r_clock, r, R_T,
        m = growth_state.m, τ = growth_state.τ,
        T = T_total, C_T = @inbounds(cumulative_total[n]),
        doubling_time = doubling_time(r),
        seed_at_renewal_start = seed0_total, seed_fraction,
        seeding_age = seeding_age(cumulative_total, n),
        incubation_pmf = inc_state.pmf)
end

"""
Relative export propensity of each province into Uganda, partially pooled.

Ituri is the reference, at one, because the traveller volume and the source
population [`exports_model`](@ref) carries are Ituri's, the point-of-entry
counts having been collected there. Every other province is measured
against it, as the chance that one of its infections is detected crossing
into Uganda relative to one of Ituri's.

The secondary weights are drawn from a common log-normal whose location and
spread are both sampled, so they pool toward each other rather than toward
a fixed number, and `tau -> 0` makes them one shared weight. The default
location prior has a median of 15% of Ituri's propensity with a 90%
interval of roughly 3% to 80%.

### What the data can say

Very little, and that is the point of reporting it. Four events reach these
streams over the whole window, three exported cases and one exported death,
and the weights enter only through the summed exporting infections. Expect
the posterior to track the prior. Read it as what the model assumes rather
than as an estimate, and read the difference it makes to the provincial
incidence split rather than the weight itself.

Returns `(; weights, pooling_sd, location)`, with `weights[1] = 1`.
"""
@model function province_export_pressure_model(n_patches::Integer;
        location_prior = Normal(log(0.15), 1.0),
        pooling_sd_prior = truncated(Normal(0, 0.5); lower = 0),
        offset_prior = Normal(0, 1))
    ## One province has nothing to pool with and no secondary weight to
    ## sample, so only the reference weight is returned.
    if n_patches <= 1
        return (; weights = ones(Float64, max(n_patches, 1)),
            pooling_sd = 0.0, location = 0.0)
    end
    μ_w ~ location_prior
    τ_w ~ pooling_sd_prior
    z_w ~ product_distribution(fill(offset_prior, n_patches - 1))
    Tw = promote_type(typeof(float(μ_w)), typeof(float(τ_w)), eltype(z_w))
    weights = ones(Tw, n_patches)
    @inbounds for p in 2:n_patches
        weights[p] = exp(μ_w + τ_w * z_w[p - 1])
    end
    return (; weights, pooling_sd = τ_w, location = μ_w)
end
