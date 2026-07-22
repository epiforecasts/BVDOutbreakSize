#!/usr/bin/env julia
#
# Standalone integral-era (v1.0.0) backfill driver.
#
# v1.0.0 forecasts the reported/suspected streams only: its
# `forecast_reported` emits `cases_new` (reported cases), `deaths_new`
# (suspected deaths) and `exports_new` (exports). Confirmed, recovered
# and isolation did not exist at this tag, so three streams are archived
# ("reported cases", "suspected deaths", "exports").
#
# Run inside a worktree checked out at v1.0.0, against its `docs` project.
# At v1.0.0 the joint model, the single-stream models and every submodel
# they use are defined INLINE in `docs/examples/analysis.jl` (there is no
# `src/forecast.jl` and the package does not export `bvd_joint`). That model
# block (analysis.jl:399-1193) is reproduced verbatim below, unchanged except
# for dropping the interleaved figure/prose lines. `forecast_reported`,
# `nuts_sample`, `load_observations` and the integration helpers the models
# call are package-exported. The headline joint call (analysis.jl:1254) and
# `forecast_reported` call (analysis.jl:1637) are reproduced verbatim in the
# footer, then the archive is built with the schema the renewal backfill
# writes (made_date, horizon, target_date, stream, draw, value).
#
# Environment overrides:
#   BVD_BACKFILL_DEST     output CSV path (default forecast_v1.0.0.csv)
#   BVD_BACKFILL_SAMPLES  post-warmup draws per chain (default 1000)
#   BVD_BACKFILL_CHAINS   chains (default 2)

using Turing
using Turing: to_submodel
using Distributions
using StatsFuns: logit, logistic
using DataFrames: DataFrame, propertynames, nrow
using CSV: CSV
using Random
using Dates: Date, Day
using Serialization: serialize
using BVDOutbreakSize

const HORIZONS = (7, 14, 21, 28)
const THIN = 5

DEST = get(ENV, "BVD_BACKFILL_DEST",
    joinpath(pkgdir(BVDOutbreakSize), "output", "backfill",
        "forecast_v1.0.0.csv"))
SAMPLES = parse(Int, get(ENV, "BVD_BACKFILL_SAMPLES", "1000"))
CHAINS = parse(Int, get(ENV, "BVD_BACKFILL_CHAINS", "2"))

obs = load_observations()
cutoff = Date(obs.as_of_date)

## Constants the analysis.jl-local model block closes over (analysis.jl:263).
const ITURI_POPULATION = obs.source_population
const ITURI_DAILY_TRAVEL = obs.daily_outbound_travellers
const EXPORTED_CASES = obs.exported_cases
const TOTAL_DEATHS = obs.total_deaths
const REPORTED_CASES = obs.reported_cases

## ---- v1.0.0 analysis.jl model block (lines 399-1193), verbatim ----

@model function exponential_growth_model(;
        tau_prior = LogNormal(log(14), 0.4),
        m_prior = truncated(Normal(7.0, 2.5);
            lower = 0, upper = 13.0))
    τ ~ tau_prior
    m ~ m_prior
    r := log(2) / τ
    T := m * τ
    C_T := 2.0 ^ m
    cumulative = s -> exp(r * s)
    return (; τ, r, m, T, C_T, cumulative)
end

#md # ```@raw html
#md # </details>
#md # ```

# ##### Onset-to-death delay
#
# Following McCabe et al., we assume the symptom-onset-to-death delay is
# gamma distributed with shape $\alpha$ and scale $\theta$, with density
# $f$ and CDF $F_d$:
#
# ```math
# \text{delay} \sim \mathrm{Gamma}(\alpha,\ \theta). \tag{4}
# ```
#
# The McCabe et al. report uses the point estimate of
# [rosello2015](@cite). We instead use the companion Bayesian reanalysis
# of the same Isiro line list [bdbv_linelist_analysis_2026](@cite),
# which re-estimates the delay with uncertainty. We carry that
# uncertainty into the fit through truncated Normal priors centred on
# its estimates:
#
# ```math
# \alpha \sim \mathrm{Normal}^{+}(4.3,\ 1.22), \qquad
# \theta \sim \mathrm{Normal}^{+}(2.6,\ 0.82). \tag{5}
# ```
#
# The delay estimation in that reanalysis follows the recommendations
# of [charniga2024](@cite).

#md # ```@raw html
#md # <details><summary>Submodel: delay_model</summary>
#md # ```

@model function delay_model(;
        alpha_prior = truncated(Normal(4.3, 1.22); lower = 0),
        theta_prior = truncated(Normal(2.6, 0.82); lower = 0))
    α ~ alpha_prior
    θ ~ theta_prior
    return (; α, θ, dist = Gamma(α, θ))
end

#md # ```@raw html
#md # </details>
#md # ```

# ##### Case-fatality ratio
#
# The US Centers for Disease Control and Prevention (CDC) summary for
# the two previous BVD outbreaks is $55$ deaths in $169$ cases
# ($\approx 33\%$), with confidence bands spanning roughly
# $24$-$40\%$. The companion Bundibugyo virus (BDBV) reanalysis reports
# a baseline of $0.47$ ($95\%$ CrI $0.31$-$0.65$) for
# non-healthcare-worker (non-HCW) confirmed cases. The prior on the
# case-fatality ratio is
#
# ```math
# \mathrm{CFR} \sim \mathrm{Beta}(6,\ 14), \tag{6}
# ```
#
# with mean $0.30$ and $95\%$ interval roughly $0.13$-$0.51$.

#md # ```@raw html
#md # <details><summary>Submodel: cfr_model</summary>
#md # ```

@model function cfr_model(; cfr_prior = Beta(6.0, 14.0))
    CFR ~ cfr_prior
    return (; CFR)
end

#md # ```@raw html
#md # </details>
#md # ```

# ##### Detection window
#
# $w$ is the mean time during which a case is still infectious and
# detectable abroad (incubation + onset-to-detection). The prior is
# based on the detection windows McCabe et al. sweep in their Method 1
# scenarios (10, 15 and 20 days): it is centred on their central 15-day
# value with an SD wide enough to cover the 10–20 day range.
#
# ```math
# w \sim \mathrm{Normal}^{+}(15,\ 5)\ \text{days}. \tag{7}
# ```

#md # ```@raw html
#md # <details><summary>Submodel: detection_window_model</summary>
#md # ```

@model function detection_window_model(;
        window_prior = truncated(Normal(15.0, 5.0); lower = 0))
    w ~ window_prior
    return (; w)
end

#md # ```@raw html
#md # </details>
#md # ```

# ##### Daily traveller volume
#
# The number of people crossing from the source area to Uganda each day
# sets the travel rate in the exports likelihood. We treat it as an
# estimated quantity rather than a fixed input. McCabe et al. Table 3
# records mean weekly passenger counts across seven points of entry; the
# Ituri-side daily total of $1871$ is a sample mean across roughly
# $15$-$21$ point-of-entry-weeks. We use a Normal prior centred on
# $1871$ with SD $200$ ($\approx 10\%$ CV), truncated at zero, covering
# point-of-entry variation and the sitrep sampling uncertainty; the
# source population is kept fixed (census).

#md # ```@raw html
#md # <details><summary>Submodel: traveller volume</summary>
#md # ```

@model function traveller_volume_model(;
        mean::Real = ITURI_DAILY_TRAVEL,
        sd::Real = ITURI_DAILY_TRAVEL_SD)
    daily_travellers ~ truncated(Normal(mean, sd); lower = 0)
    return (; daily_travellers)
end

#md # ```@raw html
#md # </details>
#md # ```

# ##### Surveillance dispersion
#
# We assume the passive-surveillance counts are reported with negative
# binomial observation error around their expected value, using the
# same error model for both streams it applies to — suspected deaths
# and reported cases in the DRC — with a single shared dispersion $k$
# because they come from the same surveillance system. Under the
# mean-$\mu$ / dispersion-$k$ parameterisation a count $Y$ has
#
# ```math
# Y \sim \mathrm{NegBinomial}(\mu,\ k), \qquad
# \mathrm{Var}(Y) = \mu + \frac{\mu^2}{k}. \tag{8}
# ```
#
# The dispersion captures passive-surveillance noise (under-reporting
# that varies by district, weekend reporting effects, and batched
# updates), not transmission heterogeneity.
# We judge this noise to be substantial, so a priori we expect
# meaningful overdispersion rather than near-Poisson counts.
# Following the Stan prior-choice recommendations
# [stan_prior_choice](@cite), the dispersion is sampled on the
# $1/\sqrt{k}$ scale, which behaves like a standard deviation, with a
# weakly-informative prior centred on that expected overdispersion,
#
# ```math
# 1/\sqrt{k} \sim \mathrm{Normal}^{+}(0.6,\ 0.2), \tag{9}
# ```
#
# giving $k$ a prior median near $3$ with a 90% range of about $1$-$14$.
# Because each stream contributes essentially one aggregate count, $k$
# is only weakly identified, so this prior carries the inference and is
# set to reflect the overdispersion we expect from passive surveillance.
# This extends the McCabe et al. report, which uses a Poisson likelihood
# for the Method 2 deaths and does not model the reported case counts at
# all; the negative binomial adds overdispersion to absorb
# passive-surveillance noise.

#md # ```@raw html
#md # <details><summary>Submodel: surveillance_dispersion_model</summary>
#md # ```

@model function surveillance_dispersion_model(;
        inv_sqrt_k_prior = truncated(Normal(0.6, 0.2); lower = 0))
    inv_sqrt_k ~ inv_sqrt_k_prior
    k := 1.0 / (inv_sqrt_k^2 + eps(typeof(inv_sqrt_k)))
    return (; k, inv_sqrt_k)
end

#md # ```@raw html
#md # </details>
#md # ```

# ##### Ascertainment — partial pooling between DRC and Uganda
#
# Two surveillance systems detect cases: DRC passive surveillance (the
# reported suspected-case count) and Uganda's point-of-entry / hospital
# surveillance (the exported-case count). Each captures only a fraction
# of the true cases passing through it, and each fraction is informed
# by essentially a single aggregate data point — the one reported-case
# total and the one export count — so neither is well identified on its
# own. We therefore centre the prior on an assumed reporting fraction
# of 25% and partially pool the two fractions so they share strength:
# treating them as identical would conflate two different systems,
# while treating them as independent would leave the Uganda fraction
# almost wholly prior-driven.
#
# Both ascertainment fractions $p_{\text{drc}}$ and $p_{\text{uganda}}$
# share a logit-scale hyperprior with mean $\mu$ and SD $\tau$:
#
# ```math
# \mu \sim \mathrm{Normal}(\mathrm{logit}(0.25),\ 1),
# \qquad
# \tau \sim \mathrm{Normal}^{+}(0,\ 0.5), \tag{10}
# ```
#
# ```math
# \mathrm{logit}(p_{\text{drc}}) \sim
#     \mathrm{Normal}(\mu,\ \tau),
# \qquad
# \mathrm{logit}(p_{\text{uganda}}) \sim
#     \mathrm{Normal}(\mu,\ \tau), \tag{11}
# ```
#
# with $p = \mathrm{logistic}(\mathrm{logit}\,p)$. Here $\tau$ is the
# pooling strength: small $\tau$ pulls the two fractions together (the
# shared-fraction limit), large $\tau$ lets them move independently
# (the separate-fraction limit). The cases likelihood uses
# $p_{\text{drc}}$; the two Uganda-side likelihoods use
# $p_{\text{uganda}}$.
#
# We sample this prior in its non-centred form: draw offsets
# $z_{\text{drc}}, z_{\text{uganda}} \sim \mathrm{Normal}(0, 1)$ and set
# $\mathrm{logit}(p) = \mu + \tau z$. This is the same prior but avoids
# the funnel geometry of the centred form, which gave hundreds of
# divergent transitions.

#md # ```@raw html
#md # <details><summary>Submodel: pooled_ascertainment_model</summary>
#md # ```

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

#md # ```@raw html
#md # </details>
#md # ```

# We model both the deaths and the reported cases with a negative
# binomial, so we define one small constructor for it and share it.
# It is parameterised by mean $\mu$ and dispersion $k$ (so the variance
# is given by equation (8)), with NaN / Inf-safe clamping on the
# success probability so extreme NUTS proposals during warmup do not
# trip the distribution domain check. It is used by the deaths and
# cases observation submodels below.

#md # ```@raw html
#md # <details><summary>Function: safe_nbinomial</summary>
#md # ```

function safe_nbinomial(k, μ)
    p_raw = k / (k + max(μ, eps(typeof(μ))))
    p = isfinite(p_raw) ?
        clamp(p_raw, eps(typeof(k)), one(k) - eps(typeof(k))) :
        eps(typeof(k))
    return NegativeBinomial(k, p)
end

#md # ```@raw html
#md # </details>
#md # ```

# #### Observation submodels
#
# With the building blocks in place, each observation submodel takes
# the growth state as input, introduces the forward integral it needs,
# and ties one data stream to the latent $C(T)$. The forward integrals
# (the at-risk person-time integral for exports, the gamma convolution
# for deaths, and the deaths-among-exports convolution) are solved
# numerically, so they support any onset-to-death delay or growth curve
# without re-derivation. Each submodel introduces its likelihood by
# referring back to the parameters defined in equations (1)-(11).

# ##### Exports — Method 1 (geographic spread)
#
# Each case in the source population travels to Uganda on any given
# day with probability
# $q = \text{daily travellers}/\text{source population}$,
# treating cases as exchangeable with the general population. A case
# is *detection-eligible* for $w$ days from infection (equation (7)).
# For a case infected at outbreak age $s \leq T$, the accumulated
# probability of being detected in Uganda by $T$ is
#
# ```math
# P(\text{detected by } T \mid \text{infected at } s)
#     = q \cdot \min(T - s,\ w). \tag{12}
# ```
#
# Splitting at $s = T - w$ (full window elapsed before $T-w$, partial
# window after) and summing across incidence $i(s)$ gives the full
# export integral
#
# ```math
# \mathbb{E}[\text{exports by }T]
#     = q \cdot \Bigl[ w \cdot C(T-w)
#          + \int_{T-w}^{T} i(s) \, (T - s) \, ds \Bigr], \tag{13}
# ```
#
# which integration by parts collapses to the cleaner at-risk person-
# time form using the cumulative-incidence trajectory $C(s)$ of
# equation (1):
#
# ```math
# \mathbb{E}[\text{exports by }T] = q \cdot \int_{T-w}^{T} C(s)\, ds. \tag{14}
# ```
#
# For exponential growth this evaluates to
# $q\cdot(C(T) - C(T-w))/r$. We evaluate equation (14) numerically so
# the same form works for any growth parameterisation, scale by the
# Uganda ascertainment fraction $p_{\text{uganda}}$ (equation (11)),
# and apply a Poisson likelihood:
#
# ```math
# \mu_e = p_{\text{uganda}} \cdot q \cdot \int_{T-w}^{T} C(s)\, ds,
# \qquad
# Y_{\text{exports}} \sim \mathrm{Poisson}(\mu_e). \tag{15}
# ```
#
# !!! note "Comparison with McCabe et al. / Imai 2020"
#     McCabe et al. use the small-$rw$ simplification
#     $\mu_e \approx q\cdot w\cdot C(T)$, the limit of equation (14)
#     as $r \to 0$.
#     For BVD's prior range $rw \in 0.33 - 2.0$ the simplification
#     under-estimates $C(T)$ by roughly $15$-$57\%$. We use a Poisson
#     likelihood for the detected exports; at the small detection
#     probability here ($p \approx q\cdot w \approx 6\cdot 10^{-3}$) it
#     is indistinguishable from a binomial detection model.

#md # ```@raw html
#md # <details><summary>Submodel: exports_model</summary>
#md # ```

@model function exports_model(
        exported_cases::Union{Missing, Integer},
        growth_state, p_uganda::Real;
        source_population::Real = ITURI_POPULATION,
        window = detection_window_model(),
        traveller = traveller_volume_model())
    cumulative = growth_state.cumulative
    T = growth_state.T

    window_state ~ to_submodel(window, false)
    w = window_state.w

    travel_state ~ to_submodel(traveller, false)
    daily_travellers = travel_state.daily_travellers

    window_start = max(T - w, zero(T))
    cumulative_window_integral := integrate_cumulative(
        cumulative, window_start, T)
    expected_exports := max(
        p_uganda * (daily_travellers / source_population) *
        cumulative_window_integral,
        eps(typeof(daily_travellers * one(T) * p_uganda))
    )

    exported_cases ~ Poisson(expected_exports)

    return (; w, daily_travellers, p_uganda,
        cumulative_window_integral, expected_exports)
end

#md # ```@raw html
#md # </details>
#md # ```

# ##### Deaths — Method 2 (back-calculation from deaths)
#
# Expected cumulative deaths by $T$ from a single seeding case is the
# CFR-weighted convolution of the cumulative-incidence trajectory
# $C(s)$ (equation (1)) with the onset-to-death density $f$
# (equation (4)):
#
# ```math
# \mathbb{E}[D(T)] = \mathrm{CFR} \cdot
#     \int_0^T e^{r s}\, f(T - s)\, ds. \tag{16}
# ```
#
# For a gamma delay this integral has an exact closed form carrying a
# $\gamma(\alpha, (\beta + r)T)/\Gamma(\alpha)$ factor from the finite
# upper limit; McCabe et
# al. use the large-$T$ simplification
# $D(T) \approx \mathrm{CFR}\cdot C(T)\cdot(1 + r/\beta)^{-\alpha}$
# (valid for $T \gtrsim 12/(\beta+r)$), which
# drops that factor. We evaluate equation (16) numerically instead,
# which is exact and lets the delay family be swapped with no change to
# the quadrature. The
# observed deaths follow the NegBinomial likelihood of equation (8)
# with the dispersion $k$ of equation (9), supplied by the composer so
# it can be shared with the cases likelihood:
#
# ```math
# Y_{\text{deaths}} \sim \mathrm{NegBinomial}(\mathbb{E}[D(T)],\ k). \tag{17}
# ```

#md # ```@raw html
#md # <details><summary>Submodel: deaths_model</summary>
#md # ```

@model function deaths_model(
        total_deaths::Union{Missing, Integer},
        growth_state, k::Real;
        delay = delay_model(),
        cfr = cfr_model())
    C_T = growth_state.C_T
    r = growth_state.r
    T = growth_state.T

    delay_state ~ to_submodel(delay, false)
    cfr_state ~ to_submodel(cfr, false)

    CFR = cfr_state.CFR

    ## NaN-safe clamp: extreme NUTS proposals during warmup can push
    ## the expected count to NaN / Inf.
    raw_deaths = expected_deaths(CFR, r, T, delay_state.dist)
    expected_deaths_T := isfinite(raw_deaths) ?
                         max(raw_deaths, eps(typeof(raw_deaths))) :
                         eps(typeof(raw_deaths))

    total_deaths ~ safe_nbinomial(k, expected_deaths_T)

    return (; CFR, delay_dist = delay_state.dist, expected_deaths_T)
end

#md # ```@raw html
#md # </details>
#md # ```

# ##### Cases — ascertainment extension (no McCabe et al. counterpart)
#
# Methods 1 and 2 use exports and deaths only. Reported
# suspected cases from the same passive-surveillance system carry
# information about $C(T)$ once the DRC ascertainment fraction $p_{\text{drc}}$
# (equation (11)) is introduced:
#
# ```math
# \mu_c = p_{\text{drc}} \cdot C(T),
# \qquad
# Y_{\text{cases}} \sim \mathrm{NegBinomial}(\mu_c,\ k). \tag{18}
# ```
#
# The dispersion $k$ (equation (9)) is shared with the deaths
# likelihood; both are sampled once by the composer.

#md # ```@raw html
#md # <details><summary>Submodel: cases_model</summary>
#md # ```

@model function cases_model(
        reported_cases::Union{Missing, Integer},
        growth_state, k::Real, p_drc::Real)
    C_T = growth_state.C_T

    raw_reports = p_drc * C_T
    expected_reports := isfinite(raw_reports) ?
                        max(raw_reports, eps(typeof(raw_reports))) :
                        eps(typeof(raw_reports))

    reported_cases ~ safe_nbinomial(k, expected_reports)

    return (; p_drc, expected_reports)
end

#md # ```@raw html
#md # </details>
#md # ```

# ##### Deaths among exports — fourth observation likelihood
#
# Uganda reports a single death among its detected exports. That count
# is informative: if the exports happened long ago, more of them would
# have died by now under the onset-to-death gamma, so the observed
# deaths-among-exports bound how recently the exports occurred and help
# constrain $T$ (and $C(T)$). The expected count reuses the at-risk
# export integrand of equation (14) but weights each case by its
# probability of having died by $T$, the onset-to-death CDF
# $F_d(T - s)$ (equation (4)), then scales by the CFR, the travel rate
# $q$ and the Uganda ascertainment fraction $p_{\text{uganda}}$:
#
# ```math
# \mathbb{E}[D_{\text{uganda}}]
#     = \mathrm{CFR} \cdot p_{\text{uganda}} \cdot q
#       \cdot \int_{T-w}^{T} C(s)\, F_d(T - s)\, ds. \tag{19}
# ```
#
# Equation (19) is evaluated numerically, writing $F_d$ as the inner
# integral of the density $f$ so the whole expression differentiates
# through $f$ alone (the reverse-mode AD backend does not support the
# gamma CDF shape-parameter derivative). The detection window $w$ and
# daily traveller volume are shared with the exports likelihood so the
# two Uganda-side observations use the same person-time, and a Poisson
# likelihood ties the observed count to equation (19):
#
# ```math
# Y_{\text{exports-deaths}} \sim
#     \mathrm{Poisson}(\mathbb{E}[D_{\text{uganda}}]). \tag{20}
# ```
#
# !!! note "Selection-bias caveat"
#     This assumes Uganda's surveillance retains detected exports
#     through to any subsequent death. If the system instead loses
#     cases that progress to death, the observed deaths-among-exports
#     count is selection-biased downward and the constraint it places
#     on $T$ is overstated.

#md # ```@raw html
#md # <details><summary>Submodel: exports_deaths_model</summary>
#md # ```

@model function exports_deaths_model(
        exports_deaths::Union{Missing, Integer},
        growth_state, CFR::Real, delay_dist, p_uganda::Real;
        window::Real,
        daily_travellers::Real,
        source_population::Real = ITURI_POPULATION)
    cumulative = growth_state.cumulative
    T = growth_state.T

    window_start = max(T - window, zero(T))
    exports_deaths_integral := integrate_exports_deaths(
        cumulative, delay_dist, window_start, T, T)
    q = daily_travellers / source_population
    raw = CFR * p_uganda * q * exports_deaths_integral
    expected_exports_deaths := isfinite(raw) ?
                               max(raw, eps(typeof(raw))) : eps(typeof(raw))

    exports_deaths ~ Poisson(expected_exports_deaths)

    return (; expected_exports_deaths)
end

#md # ```@raw html
#md # </details>
#md # ```

# #### Composers
#
# These composers stitch the building blocks into the **full
# generative models** for each analysis. McCabe et al. invert a
# deterministic summary formula at fixed nuisance parameters; here we
# sample the entire generative process — growth, delay, CFR, detection
# window, traveller volume, dispersion, ascertainment — and condition
# on the observed counts. Each composer conditionally includes only the
# likelihoods
# for the streams it carries, so a single-stream composer never
# instantiates the other observation submodels (and so a discrete
# stream is never left sampled, which would trip Turing's model check).
#
# The joint composer samples a single `surveillance_dispersion_model`
# and passes that same $k$ to both deaths and cases likelihoods, so
# they share one passive-surveillance noise scale. It also samples a
# single `pooled_ascertainment_model`, threading $p_{\text{drc}}$ to the cases
# likelihood and $p_{\text{uganda}}$ to the two Uganda-side likelihoods. The
# window $w$ and daily traveller volume sampled by the exports
# likelihood are reused by the deaths-among-exports likelihood so the
# two share person-time.

# ##### Exports-only fit — Method 1 analogue

#md # ```@raw html
#md # <details><summary>Composer: exports-only fit</summary>
#md # ```

@model function exports_only_model(
        exported_cases::Union{Missing, Integer};
        growth = exponential_growth_model(),
        exports = exports_model,
        ascertainment = pooled_ascertainment_model())
    growth_state ~ to_submodel(growth, false)
    asc_state ~ to_submodel(ascertainment, false)

    exports_state ~ to_submodel(
        exports(exported_cases, growth_state, asc_state.p_uganda), false)

    cumulative_cases := growth_state.C_T
end

#md # ```@raw html
#md # </details>
#md # ```

# ##### Deaths-only fit — Method 2 analogue

#md # ```@raw html
#md # <details><summary>Composer: deaths-only fit</summary>
#md # ```

@model function deaths_only_model(
        total_deaths::Union{Missing, Integer};
        growth = exponential_growth_model(),
        deaths = deaths_model,
        dispersion = surveillance_dispersion_model())
    growth_state ~ to_submodel(growth, false)
    dispersion_state ~ to_submodel(dispersion, false)
    k = dispersion_state.k

    deaths_state ~ to_submodel(
        deaths(total_deaths, growth_state, k), false)

    cumulative_cases := growth_state.C_T
end

#md # ```@raw html
#md # </details>
#md # ```

# ##### Cases-only fit — ascertainment extension (no McCabe et al. counterpart)

#md # ```@raw html
#md # <details><summary>Composer: cases-only fit</summary>
#md # ```

@model function cases_only_model(
        reported_cases::Union{Missing, Integer};
        growth = exponential_growth_model(),
        cases = cases_model,
        dispersion = surveillance_dispersion_model(),
        ascertainment = pooled_ascertainment_model())
    growth_state ~ to_submodel(growth, false)
    dispersion_state ~ to_submodel(dispersion, false)
    asc_state ~ to_submodel(ascertainment, false)
    k = dispersion_state.k

    cases_state ~ to_submodel(
        cases(reported_cases, growth_state, k, asc_state.p_drc), false)

    cumulative_cases := growth_state.C_T
end

#md # ```@raw html
#md # </details>
#md # ```

# ##### Deaths-among-exports-only fit (no McCabe et al. counterpart)

#md # ```@raw html
#md # <details><summary>Composer: exports-deaths-only fit</summary>
#md # ```

@model function exports_deaths_only_model(
        exports_deaths::Union{Missing, Integer};
        growth = exponential_growth_model(),
        delay = delay_model(),
        cfr = cfr_model(),
        window = detection_window_model(),
        traveller = traveller_volume_model(),
        exports_deaths_model = exports_deaths_model,
        ascertainment = pooled_ascertainment_model(),
        source_population::Real = ITURI_POPULATION)
    growth_state ~ to_submodel(growth, false)
    delay_state ~ to_submodel(delay, false)
    cfr_state ~ to_submodel(cfr, false)
    window_state ~ to_submodel(window, false)
    asc_state ~ to_submodel(ascertainment, false)

    travel_state ~ to_submodel(traveller, false)
    daily_travellers = travel_state.daily_travellers

    exports_deaths_state ~ to_submodel(
        exports_deaths_model(exports_deaths, growth_state,
            cfr_state.CFR, delay_state.dist, asc_state.p_uganda;
            window = window_state.w,
            daily_travellers = daily_travellers,
            source_population = source_population),
        false)

    cumulative_cases := growth_state.C_T
end

#md # ```@raw html
#md # </details>
#md # ```

# ##### Joint fit — full posterior over $C(T)$ from all data streams
#
# The joint composer is the same generative process conditioned on all four
# observed streams simultaneously. Each stream argument may be passed
# as `missing` to drop it; the matching likelihood is then not
# instantiated, so the composer doubles as a generator (all streams
# missing) for the prior- and posterior-predictive checks.

#md # ```@raw html
#md # <details><summary>Composer: joint fit</summary>
#md # ```

@model function bvd_joint(
        exported_cases::Union{Missing, Integer},
        total_deaths::Union{Missing, Integer},
        reported_cases::Union{Missing, Integer} = missing,
        exports_deaths::Union{Missing, Integer} = missing;
        growth = exponential_growth_model(),
        exports = exports_model,
        deaths = deaths_model,
        cases = cases_model,
        exports_deaths_model = exports_deaths_model,
        dispersion = surveillance_dispersion_model(),
        ascertainment = pooled_ascertainment_model(),
        source_population::Real = ITURI_POPULATION)
    growth_state ~ to_submodel(growth, false)
    dispersion_state ~ to_submodel(dispersion, false)
    asc_state ~ to_submodel(ascertainment, false)
    k = dispersion_state.k
    p_drc = asc_state.p_drc
    p_uganda = asc_state.p_uganda

    exports_state ~ to_submodel(
        exports(exported_cases, growth_state, p_uganda), false)
    deaths_state ~ to_submodel(
        deaths(total_deaths, growth_state, k), false)
    cases_state ~ to_submodel(
        cases(reported_cases, growth_state, k, p_drc), false)
    exports_deaths_state ~ to_submodel(
        exports_deaths_model(exports_deaths, growth_state,
            deaths_state.CFR, deaths_state.delay_dist, p_uganda;
            window = exports_state.w,
            daily_travellers = exports_state.daily_travellers,
            source_population = source_population),
        false)

    cumulative_cases := growth_state.C_T
end

#md # ```@raw html
#md # </details>
#md # ```

# ##### McCabe et al. reimplementation — exports and deaths only
#
# McCabe et al.'s joint configuration uses exactly two data sources: the
# geographic-spread exports (Method 1) and the back-calculation from
# deaths (Method 2). It has no reported-cases ascertainment model and
# no deaths-among-exports likelihood. This composer wraps just
# those two observation submodels so the sense-check can fix the model
# down to exactly the McCabe et al. joint configuration. Either stream may
# be `missing`: passing `missing` for exports recovers a pure Method 2
# (deaths-only) fit without instantiating the exports likelihood.

#md # ```@raw html
#md # <details><summary>Composer: report reimplementation</summary>
#md # ```

@model function imperial_only_model(
        exported_cases::Union{Missing, Integer},
        total_deaths::Union{Missing, Integer};
        growth = exponential_growth_model(),
        exports = exports_model,
        deaths = deaths_model,
        dispersion = surveillance_dispersion_model(),
        ascertainment = pooled_ascertainment_model())
    growth_state ~ to_submodel(growth, false)
    dispersion_state ~ to_submodel(dispersion, false)
    asc_state ~ to_submodel(ascertainment, false)
    k = dispersion_state.k
    p_uganda = asc_state.p_uganda

    if !ismissing(exported_cases)
        exports_state ~ to_submodel(
            exports(exported_cases, growth_state, p_uganda), false)
    end
    deaths_state ~ to_submodel(
        deaths(total_deaths, growth_state, k), false)

    cumulative_cases := growth_state.C_T
end

## ---- fit + forecast + archive ----

@info "fitting v1.0.0 integral joint" cutoff samples=SAMPLES chains=CHAINS

## Headline joint call, verbatim from v1.0.0 analysis.jl:1254. v1.0.0's
## `bvd_joint` takes four positional streams and no genetic/detection kwargs.
chn = nuts_sample(
    bvd_joint(obs.exported_cases, obs.total_deaths,
        obs.reported_cases, obs.exports_deaths);
    samples = SAMPLES, chains = CHAINS)

## Keep the chain: this tag has no `fit_or_load`, so a failure in the
## forecast/archive step below would otherwise cost a whole refit.
chain_path = "$(DEST).chain"
mkpath(dirname(DEST))
try
    serialize(chain_path, chn)
    @info "chain saved" chain_path
catch e
    @warn "could not save chain" exception = e
end

## Forecast call, verbatim from v1.0.0 analysis.jl:1637.
runs = [(h,
            forecast_reported(chn;
                horizon = h,
                daily_travellers = ITURI_DAILY_TRAVEL,
                source_population = ITURI_POPULATION,
                obs_cases = REPORTED_CASES,
                obs_deaths = TOTAL_DEATHS,
                obs_exports = EXPORTED_CASES))
        for h in HORIZONS]

## Three scoreable reported/suspected streams; confirmed, recovered and
## isolation did not exist yet. Labels match the scorer's STREAM_HISTORY.
streams = (
    (:cases_new, "reported cases"),
    (:deaths_new, "suspected deaths"),
    (:exports_new, "exports"))
out = DataFrame(made_date = Date[], horizon = Int[], target_date = Date[],
    stream = String[], draw = Int[], value = Float64[])
for (horizon, fc) in runs
    h = Int(horizon)
    target = cutoff + Day(h)
    for (col, label) in streams
        col in propertynames(fc) || continue
        vals = fc[!, col]
        for (d, i) in enumerate(1:THIN:length(vals))
            push!(out, (cutoff, h, target, label, d, Float64(vals[i])))
        end
    end
end
isempty(out) && error("no archived streams in the forecast")
CSV.write(DEST, out)
@info "wrote archive" dest=DEST rows=nrow(out) cutoff streams=unique(out.stream)
