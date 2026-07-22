#!/usr/bin/env julia
#
# Standalone integral-era (v1.1.0) backfill driver.
#
# v1.1.0 forecasts the reported/suspected streams only: its
# `forecast_reported` emits `cases_new` (reported cases), `deaths_new`
# (suspected deaths) and `exports_new` (exports). Confirmed, recovered
# and isolation did not exist at this tag, so three streams are archived
# ("reported cases", "suspected deaths", "exports").
#
# Run inside a worktree checked out at v1.1.0, against its `docs` project.
# At v1.1.0 the joint model, the single-stream models, every submodel they
# use AND `genetic_seeding_model` are defined INLINE in
# `docs/examples/analysis.jl` (the package does not export them). That model
# block (analysis.jl:445-1415) is reproduced verbatim below, unchanged except
# for dropping the interleaved figure/prose lines (the `cfr_prior_fig` plot).
# `forecast_reported`, `nuts_sample`, `load_observations` and the integration
# helpers the models call are package-exported. The headline joint call
# (analysis.jl:1479) and `forecast_reported` call (analysis.jl:1934) are
# reproduced verbatim in the footer, then the archive is built with the
# schema the renewal backfill writes (made_date, horizon, target_date,
# stream, draw, value).
#
# Environment overrides:
#   BVD_BACKFILL_DEST     output CSV path (default forecast_v1.1.0.csv)
#   BVD_BACKFILL_SAMPLES  post-warmup draws per chain (default 1000)
#   BVD_BACKFILL_CHAINS   chains (default 2)

using Turing
using Turing: to_submodel, @varname
using Distributions
using StatsFuns: logit, logistic
using DataFrames: DataFrame, propertynames, nrow
using CSV: CSV
using Random
using Dates: Date, Day, value
using Serialization: serialize
using BVDOutbreakSize

const HORIZONS = (7, 14, 21, 28)
const THIN = 5

DEST = get(ENV, "BVD_BACKFILL_DEST",
    joinpath(pkgdir(BVDOutbreakSize), "output", "backfill",
        "forecast_v1.1.0.csv"))
SAMPLES = parse(Int, get(ENV, "BVD_BACKFILL_SAMPLES", "1000"))
CHAINS = parse(Int, get(ENV, "BVD_BACKFILL_CHAINS", "2"))

obs = load_observations()
cutoff = Date(obs.as_of_date)

## Constants the analysis.jl-local model block closes over (analysis.jl:300).
const ITURI_POPULATION = obs.source_population
const ITURI_DAILY_TRAVEL = obs.daily_outbound_travellers
const EXPORTED_CASES = obs.exported_cases
const EXPORTS_DEATHS = obs.exports_deaths

## ---- v1.1.0 analysis.jl model block (lines 445-1415), verbatim ----

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

# ##### Genetic seeding bound
#
# A BEAST time tree of the first ten sequenced genomes
# [virological2026](@cite) places the TMRCA, the age of the oldest
# internal node of the tree, at a mean of 25 March 2026, with a 95% HPD
# interval of about $\pm 30$ days.
# The temporal sampling range is too short to estimate the clock, so it
# is fixed. The source analysis considers two literature rates for the
# 2013–2016 West African Ebola epidemic [holmes2016](@cite): a
# $1.2\times10^{-3}$ substitutions/site/year rate across all public
# data, and a faster $1.9\times10^{-3}$ early-epidemic rate that dates
# the TMRCA more recently. We use the $1.2\times10^{-3}$ rate in the
# main analysis and the $1.9\times10^{-3}$ rate in the
# [clock-rate sensitivity](#Clock-rate-sensitivity).
# This is a lower bound on the seeding time $T$: adding sequences, or
# more geographically representative ones, can only push the TMRCA
# earlier, never later (the sampled tree is almost entirely from Bunia).
# Combining the genetic TMRCA with the other data streams as a seeding
# bound follows a suggestion of N. Ferguson [ferguson2026](@cite).
# We parameterise the bound as an uncertain threshold
# $B \sim \mathrm{Normal}(g, \sigma)$, where
# $g = t_{\mathrm{cut}} - t_{\mathrm{TMRCA}}$ is the data cut-off date
# minus the reported TMRCA date (so it tracks the cut-off rather than a
# fixed offset), and require $T \ge B$, leaving $T$ free above it.
# Marginalising over $B$ gives a soft one-sided likelihood,
#
# ```math
# p_\text{gen}(T) = \Pr[B \le T] = \Phi\!\left(\frac{T - g}{\sigma}\right),
# \qquad \sigma = 15\ \text{d}. \tag{3a}
# ```

#md # ```@raw html
#md # <details><summary>Submodel: genetic_seeding_model</summary>
#md # ```

@model function genetic_seeding_model(T, tmrca_days::Real;
        tmrca_days_sd::Real)
    ## The molecular-clock TMRCA is a right-censored, noisy reading of the
    ## true seeding time T: deeper or wider sampling only pushes it older,
    ## so we learn the reading is at least `tmrca_days`, i.e. P(read ≥ g).
    tmrca_days ~ censored(Normal(T, tmrca_days_sd); upper = tmrca_days)
    return (; tmrca_days, tmrca_days_sd)
end

# Observing $g$ (`tmrca_days`) at the upper censoring point of
# $\mathrm{Normal}(T, \sigma)$ contributes the log probability of the
# censored upper tail, which is the soft bound of Eq. (3a):
#
# ```math
# \Pr[\mathrm{Normal}(T, \sigma) \ge g]
#   = 1 - \Phi\!\left(\frac{g - T}{\sigma}\right)
#   = \Phi\!\left(\frac{T - g}{\sigma}\right).
# ```

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
# ($\approx 33\%$;
# [CDC outbreak history](https://www.cdc.gov/ebola/outbreaks/index.html)),
# with confidence bands spanning
# roughly $26$-$40\%$. The companion Bundibugyo virus (BDBV) reanalysis
# reports a baseline of $0.47$ ($95\%$ CrI $0.31$-$0.65$) for
# non-healthcare-worker (non-HCW) confirmed cases. The prior on the
# case-fatality ratio is
#
# ```math
# \mathrm{CFR} \sim \mathrm{Beta}(6.6,\ 13.4), \tag{6}
# ```
#
# with mean $0.33$ and $95\%$ interval roughly $0.15$-$0.54$. The mean
# matches the CDC $55/169 \approx 33\%$ figure and the corrected central
# CFR in the 20 May report [mccabe2026update](@cite).

#md # ```@raw html
#md # <details><summary>Submodel: cfr_model</summary>
#md # ```

@model function cfr_model(; cfr_prior = Beta(6.6, 13.4))
    CFR ~ cfr_prior
    return (; CFR)
end

#md # ```@raw html
#md # </details>
#md # ```

# The prior density, with the CDC $0.33$ figure marked, as a sense
# check.

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
# Because the prior allows near-Poisson counts, $k$ itself ranges over
# many orders of magnitude, so the pair plots and summary table show
# dispersion on both the sampled $1/\sqrt{k}$ scale, which is easier to
# read, and the more familiar $k$ scale.
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
# Both ascertainment fractions $p_{\text{DRC}}$ and $p_{\text{Uganda}}$
# share a logit-scale hyperprior with mean $\mu$ and SD $\tau$:
#
# ```math
# \mu \sim \mathrm{Normal}(\mathrm{logit}(0.25),\ 1),
# \qquad
# \tau \sim \mathrm{Normal}^{+}(0,\ 0.5), \tag{10}
# ```
#
# ```math
# \mathrm{logit}(p_{\text{DRC}}) \sim
#     \mathrm{Normal}(\mu,\ \tau),
# \qquad
# \mathrm{logit}(p_{\text{Uganda}}) \sim
#     \mathrm{Normal}(\mu,\ \tau), \tag{11}
# ```
#
# with $p = \mathrm{logistic}(\mathrm{logit}\,p)$. Here $\tau$ is the
# pooling strength: small $\tau$ pulls the two fractions together (the
# shared-fraction limit), large $\tau$ lets them move independently
# (the separate-fraction limit). The cases likelihood uses
# $p_{\text{DRC}}$; the two Uganda-side likelihoods use
# $p_{\text{Uganda}}$.
#
# We sample this prior in its non-centred form: draw offsets
# $z_{\text{DRC}}, z_{\text{Uganda}} \sim \mathrm{Normal}(0, 1)$ and set
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
# numerically. Each submodel introduces its likelihood by
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
# which integration by parts collapses to the at-risk person-time form
# using the cumulative-incidence trajectory $C(s)$ of
# equation (1):
#
# ```math
# \mathbb{E}[\text{exports by }T] = q \cdot \int_{T-w}^{T} C(s)\, ds. \tag{14}
# ```
#
# For exponential growth this evaluates to
# $q\cdot(C(T) - C(T-w))/r$. We evaluate equation (14) numerically so
# the same form works for any growth parameterisation, scale by the
# Uganda ascertainment fraction $p_{\text{Uganda}}$ (equation (11)),
# and apply a Poisson likelihood:
#
# ```math
# \mu_e = p_{\text{Uganda}} \cdot q \cdot \int_{T-w}^{T} C(s)\, ds,
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

    q = daily_travellers / source_population
    expected_exports_T := expected_exports(cumulative, p_uganda, q, T, w)

    exported_cases ~ Poisson(expected_exports_T)

    return (; w, daily_travellers, p_uganda,
        expected_exports = expected_exports_T)
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
# which is exact. The
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

# ##### Cases — ascertainment extension
# Methods 1 and 2 use exports and deaths only. Reported
# suspected cases from the same passive-surveillance system carry
# information about $C(T)$ once the DRC ascertainment fraction $p_{\text{DRC}}$
# (equation (11)) is introduced:
#
# ```math
# \mu_c = p_{\text{DRC}} \cdot C(T),
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

# ##### Deaths among exports
#
# Uganda reports a single death among its detected exports. That count
# is informative. If the exports happened long ago, we would expect
# more deaths among them by now under the onset-to-death gamma, so the
# observed deaths-among-exports bound how recently the exports occurred
# and help constrain $T$ (and $C(T)$). The expected count reuses the at-risk
# export integrand of equation (14) but weights each case by its
# probability of having died by $T$, the onset-to-death CDF
# $F_d(T - s)$ (equation (4)), then scales by the CFR, the travel rate
# $q$ and the Uganda ascertainment fraction $p_{\text{Uganda}}$:
#
# ```math
# \mathbb{E}[D_{\text{Uganda}}]
#     = \mathrm{CFR} \cdot p_{\text{Uganda}} \cdot q
#       \cdot \int_{T-w}^{T} C(s)\, F_d(T - s)\, ds. \tag{19}
# ```
#
# Equation (19) is evaluated numerically, writing $F_d$ as the inner
# integral of the density $f$ so the whole expression differentiates
# through $f$ alone (the reverse-mode AD backend does not support the
# gamma CDF shape-parameter derivative). The detection window $w$ and
# daily traveller volume are shared with the exports likelihood so the
# two Uganda-side observations use the same person-time.
#
# Uganda's export deaths are point-of-entry / hospital-detected, so
# their dates are recorded directly and carry information beyond the
# total count: a death seen early bounds how old the outbreak can be. We
# model the detected export deaths as an inhomogeneous Poisson process
# with cumulative intensity $\mathbb{E}[D_{\text{Uganda}}(t)]$
# (equation (19), at any elapsed time $t$) and use its time-resolved
# likelihood. Splitting $[0, T]$ at the earliest
# dated death (offset $\Delta_1$ before the cut-off, elapsed time
# $T-\Delta_1$): before it no export death was seen, contributing one
# continuous survival term over $[0, T-\Delta_1]$; from that day to the
# cut-off each day $d$ carries a Poisson count of the export deaths that
# day, with bin mean $\mu_d$,
#
# ```math
# \log L = \sum_d y_d \log \mu_d - \mathbb{E}[D_{\text{Uganda}}(T)],
# \quad
# \mu_d = \mathbb{E}[D_{\text{Uganda}}(t_d)]
#         - \mathbb{E}[D_{\text{Uganda}}(t_{d-1})]. \tag{20}
# ```
#
# The continuous term collapses the long pre-death zero stretch into a
# single weight, so there is no per-day vector of zeros, while the
# recent window is resolved per day. The number of daily bins is fixed
# by the earliest death's date, so the likelihood is well defined even
# though $T$ is latent; with one death the series is a single count, and
# it takes more dated deaths as they are reported.
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
        export_deaths_daily::AbstractVector,
        growth_state, CFR::Real, delay_dist, p_uganda::Real;
        pre_start_deaths::Union{Missing, Integer} = 0,
        window::Real,
        daily_travellers::Real,
        source_population::Real = ITURI_POPULATION)
    cumulative = growth_state.cumulative
    T = growth_state.T
    q = daily_travellers / source_population
    n = length(export_deaths_daily)   # days from earliest death to cut-off

    ## Precompute the onset-to-death CDF once and reuse it across every
    ## bin edge below (`T - s ≤ window` over the domain; see
    ## `ExportDeathDelay`).
    delay = ExportDeathDelay(delay_dist, window)
    Λ(t) = expected_exports_deaths(
        cumulative, delay, CFR, p_uganda, q, t, window)

    ## Pre-death zero stretch as one Poisson observed at 0; `missing`
    ## generates it for predictive checks (see equation (20)).
    pre = T - n > zero(T) ? Λ(T - n) : zero(T)
    pre_start_deaths ~ Poisson(max(pre, zero(pre)))

    ## Carry the upper edge forward so each Λ is evaluated once.
    λlo = pre
    for i in 1:n
        λhi = Λ(T - n + i)
        μ_day = max(λhi - λlo, eps(typeof(λhi)))
        export_deaths_daily[i] ~ Poisson(μ_day)
        λlo = λhi
    end

    return (;)
end

#md # ```@raw html
#md # </details>
#md # ```

# ##### First export detection — timing survival term
#
# The same logic applies to the *first detected export case* (Uganda's
# first hospital admission), using the at-risk export person-time
# intensity $\mathbb{E}[\text{exports}(t)]$ (equation (15)) in place of
# the export-death intensity. With $\Delta$ the offset from the first
# admission date to the cut-off and $t_1 = T - \Delta$,
#
# ```math
# \Pr(\text{no export detected before } t_1)
#     = \exp\!\bigl(-\mathbb{E}[\text{exports}(t_1)]\bigr). \tag{22}
# ```
#
# As with the export-death term, this is one-sided and only marginally
# constrains the posterior because the Uganda detections sit only days
# before the cut-off.
# Passing `delta = missing` makes the submodel a no-op.

#md # ```@raw html
#md # <details><summary>Submodel: exports_detection_timing_model</summary>
#md # ```

@model function exports_detection_timing_model(
        growth_state, p_uganda::Real;
        delta::Union{Missing, Real},
        pre_detection_exports::Union{Missing, Integer} = 0,
        window::Real,
        daily_travellers::Real,
        source_population::Real = ITURI_POPULATION)
    if !ismissing(delta)
        cumulative = growth_state.cumulative
        T = growth_state.T
        t1 = T - delta
        q = daily_travellers / source_population
        survived_exports := t1 <= zero(T) ? zero(T) :
                            expected_exports(cumulative, p_uganda, q, t1, window)
        ## No detection before t1 as a Poisson observed at 0; `missing`
        ## generates it for predictive checks (see equation (22)).
        pre_detection_exports ~ Poisson(max(survived_exports, zero(T)))
    end

    return (;)
end

#md # ```@raw html
#md # </details>
#md # ```

# #### Composers
#
# These composers combine the building blocks into the full model for
# each analysis. McCabe et al. invert a
# deterministic summary formula at fixed nuisance parameters; here we
# sample all of them — growth, delay, CFR, detection
# window, traveller volume, dispersion, ascertainment — and condition
# on the observed counts. Each composer conditionally includes only the
# likelihoods
# for the streams it carries, so a single-stream composer never
# instantiates the other observation submodels (and so a discrete
# stream is never left sampled).
# Of the observation streams, the geographic-spread exports reproduce
# McCabe et al.'s Method 1 and the back-calculation from deaths their
# Method 2; the reported-cases ascertainment, the deaths-among-exports,
# the export-detection-timing, and the genetic seeding terms are
# extensions with no counterpart in McCabe et al.
#
# The joint composer samples a single dispersion scale $k$ and passes it
# to both the deaths and cases likelihoods, so they share one passive-
# surveillance noise scale. It also samples a single pooled set of
# ascertainment fractions, threading $p_{\text{DRC}}$ to the cases
# likelihood and $p_{\text{Uganda}}$ to the two Uganda-side likelihoods. The
# window $w$ and daily traveller volume sampled by the exports
# likelihood are reused by the deaths-among-exports likelihood so the
# two share person-time.
#
# We write single-stream composers for the four count-based streams
# only.
# The export-detection-timing and genetic seeding terms constrain the
# outbreak start $T$ rather than the size directly and are weakly
# identified in isolation, so we do not fit them on their own; the
# joint composer still conditions on both.

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

# ##### Cases-only fit — ascertainment extension
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

# ##### Deaths-among-exports-only fit
#md # ```@raw html
#md # <details><summary>Composer: exports-deaths-only fit</summary>
#md # ```

@model function exports_deaths_only_model(
        export_deaths_daily::AbstractVector;
        growth = exponential_growth_model(),
        delay = delay_model(),
        cfr = cfr_model(),
        window = detection_window_model(),
        traveller = traveller_volume_model(),
        exports_deaths_model = exports_deaths_model,
        ascertainment = pooled_ascertainment_model(),
        source_population::Real = ITURI_POPULATION,
        pre_start_deaths::Union{Missing, Integer} = 0)
    growth_state ~ to_submodel(growth, false)
    delay_state ~ to_submodel(delay, false)
    cfr_state ~ to_submodel(cfr, false)
    window_state ~ to_submodel(window, false)
    asc_state ~ to_submodel(ascertainment, false)

    travel_state ~ to_submodel(traveller, false)
    daily_travellers = travel_state.daily_travellers

    exports_deaths_state ~ to_submodel(
        exports_deaths_model(export_deaths_daily, growth_state,
            cfr_state.CFR, delay_state.dist, asc_state.p_uganda;
            pre_start_deaths = pre_start_deaths,
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
        export_deaths_daily::AbstractVector = Int[];
        growth = exponential_growth_model(),
        exports = exports_model,
        deaths = deaths_model,
        cases = cases_model,
        exports_deaths_model = exports_deaths_model,
        exports_detection_timing = exports_detection_timing_model,
        dispersion = surveillance_dispersion_model(),
        ascertainment = pooled_ascertainment_model(),
        genetic = nothing,
        source_population::Real = ITURI_POPULATION,
        pre_start_deaths::Union{Missing, Integer} = 0,
        pre_detection_exports::Union{Missing, Integer} = 0,
        first_export_detection_delta::Union{Missing, Real} = missing)
    growth_state ~ to_submodel(growth, false)
    if genetic !== nothing
        genetic_state ~ to_submodel(genetic(growth_state.T), false)
    end
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
        exports_deaths_model(export_deaths_daily, growth_state,
            deaths_state.CFR, deaths_state.delay_dist, p_uganda;
            pre_start_deaths = pre_start_deaths,
            window = exports_state.w,
            daily_travellers = exports_state.daily_travellers,
            source_population = source_population),
        false)
    detection_timing_state ~ to_submodel(
        exports_detection_timing(growth_state, p_uganda;
            delta = first_export_detection_delta,
            pre_detection_exports = pre_detection_exports,
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

@info "fitting v1.1.0 integral joint" cutoff samples=SAMPLES chains=CHAINS

## Genetic-seeding closure, analysis.jl-local (analysis.jl:1476).
genetic_seeding = T -> genetic_seeding_model(T, obs.genetic_tmrca_days;
    tmrca_days_sd = obs.genetic_tmrca_days_sd)

## Headline joint call, verbatim from v1.1.0 analysis.jl:1479.
chn = nuts_sample(
    bvd_joint(obs.exported_cases, obs.total_deaths,
        obs.reported_cases, obs.export_deaths_daily;
        first_export_detection_delta = obs.first_export_detection_delta,
        genetic = genetic_seeding);
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

## Forecast call, verbatim from v1.1.0 analysis.jl:1934.
runs = [(h,
            forecast_reported(chn;
                horizon = h,
                daily_travellers = ITURI_DAILY_TRAVEL,
                source_population = ITURI_POPULATION,
                obs_cases = obs.reported_cases,
                obs_deaths = obs.total_deaths,
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
