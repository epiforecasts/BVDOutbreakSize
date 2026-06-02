# Prior submodels: building blocks shared across the observation
# submodels and joint composers. Each `@model` is a small piece of the
# generative process — a single prior or a small group of related
# priors — so it can be reused in different composers without code
# duplication.

"""
Prior on the cumulative infection count `C(T) = exp(r·T)` via a doublings
parameterisation. Samples the exponential growth rate `r` and a
doubling-count `m = T/τ`, then exposes `(τ, m, r, T, C_T, cumulative)`
as deterministics for downstream submodels. `C(T)` is the latent
*infection* count; the composers map it onto symptom onsets through the
incubation period (see [`incubation_model`](@ref)) before the per-stream
delays act.

McCabe et al.'s primary assumption is the doubling time (their 7/14/21-day
sweep); each doubling time implies a growth rate `r = log(2)/τ`, and the
prior is placed on that implied `r`, with `τ = log(2)/r` recovered as a
deterministic. The default `r ~ LogNormal(log(log(2)/14), 0.4)` is
exactly equivalent to a `LogNormal(log(14), 0.4)` prior on the doubling
time: because `r = log(2)/τ` is a reciprocal, the log-scale SD `0.4` is
preserved, so the implied doubling-time prior (and every derived
quantity) matches that prior on `τ`.

The default doubling-count prior `m ~ Normal(M_PRIOR_BASE, 3)` (truncated
at 0) is centred on `M_PRIOR_BASE = 9` (`C_T = 2^9 = 512`), the doubling
count implied by McCabe et al.'s first-report (18 May 2026) Method 2
central scenario of 501 cases (`log2(501) ≈ 9`). `C_T = 2^m` is now the
cumulative *infection* count; 9 is kept as a weakly-informative centre of
the same order rather than rescaled to a case-equivalent infection count.
For a fit at a later cut-off, pass an `m_prior` whose centre advances from
that base date via [`m_prior_centre`](@ref) so the prior tracks the
elapsed time. This is a weakly-informative centring choice (SD 3 gives 95%
support ≈ `m ∈ (3, 15)`, `C_T ∈ (8, 32000)`); the fit is dominated by the
likelihood, so it mainly sets where the joint sampler starts.
"""
@model function exponential_growth_model(;
        r_prior = LogNormal(log(log(2) / 14), 0.4),
        m_prior = truncated(Normal(M_PRIOR_BASE, 3.0); lower = 0))
    r ~ r_prior
    m ~ m_prior
    τ := log(2) / r
    T := m * τ
    C_T := 2.0 ^ m
    cumulative = s -> exp(r * s)
    return (; τ, r, m, T, C_T, cumulative)
end

"""
One-sided molecular-clock seeding bound: the TMRCA is treated as a
right-censored, noisy reading of the latent seeding time `T`, so the
likelihood contributes `P(read >= tmrca_days)`. See also
[`exponential_growth_model`](@ref).
"""
@model function genetic_seeding_model(T, tmrca_days::Real;
        tmrca_days_sd::Real)
    ## The molecular-clock TMRCA is a right-censored, noisy reading of the
    ## true seeding time T: deeper or wider sampling only pushes it older,
    ## so we learn the reading is at least `tmrca_days`, i.e. P(read ≥ g).
    tmrca_days ~ censored(Normal(T, tmrca_days_sd); upper = tmrca_days)
    return (; tmrca_days, tmrca_days_sd)
end

"""
Onset-to-death delay prior. Samples a gamma shape `α` and scale `θ`
from truncated-normal priors centred on the Bayesian BDBV line-list
reanalysis estimates, and returns the resulting `Gamma(α, θ)`
distribution for use in [`deaths_model`](@ref) and
[`exports_deaths_model`](@ref).
"""
@model function delay_model(;
        alpha_prior = truncated(Normal(4.3, 1.22); lower = 0),
        theta_prior = truncated(Normal(2.6, 0.82); lower = 0))
    α ~ alpha_prior
    θ ~ theta_prior
    return (; α, θ, dist = Gamma(α, θ))
end

"""
Onset-to-report delay prior. Samples a gamma shape `α_rep` and scale
`θ_rep` from truncated-normal priors centred on the BDBV linelist
posterior on the symptom-onset to suspected-case-notification delay,
loosened to allow for 2026-specific deviations. Used by
[`reported_cases_model`](@ref).
"""
@model function report_delay_model(;
        alpha_prior = truncated(Normal(2.5, 1.0); lower = 0.1),
        theta_prior = truncated(Normal(4.5, 1.5); lower = 0.1))
    α_rep ~ alpha_prior
    θ_rep ~ theta_prior
    return (; α = α_rep, θ = θ_rep, dist = Gamma(α_rep, θ_rep))
end

"""
Report-to-lab-confirmation delay prior. Samples a gamma shape `α_lab`
and scale `θ_lab` from truncated-normal priors with a heavy right tail
to allow for sample shipment to a confirmatory lab. No per-sample
outbreak data anchors this prior. Used by [`confirmed_cases_model`](@ref).
"""
@model function lab_delay_model(;
        alpha_prior = truncated(Normal(1.5, 1.0); lower = 0.1),
        theta_prior = truncated(Normal(3.0, 2.0); lower = 0.1))
    α_lab ~ alpha_prior
    θ_lab ~ theta_prior
    return (; α = α_lab, θ = θ_lab, dist = Gamma(α_lab, θ_lab))
end

"""
Incubation-period prior, the infection-to-symptom-onset delay that sits
at the head of the generative process. Samples the mean and coefficient
of variation, recovers the gamma shape `α_inc = 1/CV²` and scale
`θ_inc = mean·CV²` as deterministics, and returns the resulting
`Gamma(α_inc, θ_inc)` distribution.

The incubation period cannot be fitted from the BDBV line list (the
Rosello deposit has no exposure dates), so the line-list reanalysis
recommends the MacNeil et al. (2010) Bundibugyo estimate from the 2007
Uganda outbreak instead: a mean of 6.3 days (95% CI 5.2-7.3, n = 24).
The prior is placed on the mean and the CV so MacNeil's uncertainty is
carried directly: `mean ~ Normal(6.3, 0.54)` reproduces the reported 95%
CI 5.2-7.3 (SD = CI half-width / 1.96). MacNeil give no interval on the
spread, so the CV prior `Normal(0.55, 0.12)` is a weakly-informative
modelling choice spanning the dispersion implied by their 2-20 day
observed range; the CV prior needs human sign-off, the mean prior does
not.

Used by the joint composer and the onset-driven single-stream composers
to map the latent cumulative *infections* `C(T) = exp(r·T)` onto the
cumulative symptom onsets the downstream delays act on. Under exponential
growth the convolution is the exact constant rescale `mgf(incubation,
−r)`, applied as the `onset_fraction` of [`deaths_model`](@ref),
[`reported_cases_model`](@ref), [`confirmed_cases_model`](@ref) and
[`exports_deaths_model`](@ref).
"""
@model function incubation_model(;
        mean_prior = truncated(Normal(6.3, 0.54); lower = 0.1),
        cv_prior = truncated(Normal(0.55, 0.12); lower = 0.1))
    mean_inc ~ mean_prior
    cv_inc ~ cv_prior
    α_inc := inv(cv_inc^2)
    θ_inc := mean_inc * cv_inc^2
    return (; α = α_inc, θ = θ_inc, mean = mean_inc, cv = cv_inc,
        dist = Gamma(α_inc, θ_inc))
end

"""
PCR sensitivity prior. Beta(6, 2): mean 0.75, 95% interval 0.39-0.97.
Confirmation runs on the altona RealStar Filovirus Screen RT-PCR, which
detects Bundibugyo virus at 11-67 RNA copies per reaction; the rapid
Cepheid GeneXpert Ebola assay is Zaire-ebolavirus-specific and does not
reliably detect Bundibugyo. The prior keeps good analytical sensitivity
plausible while carrying substantial downside mass for early
low-viral-load specimens and field handling. Used by
[`confirmed_cases_model`](@ref).
"""
@model function test_sensitivity_model(;
        sensitivity_prior = Beta(6.0, 2.0))
    s_test ~ sensitivity_prior
    return (; s_test)
end

"""
Priority (triage) testing strength prior. The laboratory does not analyse
the suspect pool in arrival order; it triages, testing the
highest-pre-test-probability (most-likely-BVD) samples first. The
cumulative BVD samples *tested* therefore saturate toward the BVD
available in the pool as the analysed count grows, while later batches
drain a low-yield background backlog. This produces the observed
saturating confirmed positives (101, 105, 106, 121) against a near-
doubling analysed volume (211, 295, 295, 403) and the falling per-vintage
positivity (0.48 → 0.30) that a constant or ramped background cannot.

Parameterised by a single front-loading strength `κ ≥ 1` (see
[`priority_bvd_tested`](@ref)): with `B` BVD available, `N` the total
available pool and `A` analysed, the cumulative BVD tested is
`B · (1 − (1 − A/N)^κ)`. `κ = 1` is proportional (no-priority) sampling,
`bvd_tested = B·A/N`, which makes the per-test positivity the available
BVD share `s·B/N` — exactly the old constant-background model, so
`Δκ = 0` (the `priority_off` flag) is the backward-compatible reduction.
Larger `κ` front-loads BVD more strongly; `κ → ∞` creams all BVD into the
first analysed samples (positivity → `s` early), which over-predicts the
observed early 0.48, so the prior keeps `κ` soft. `κ = 1 + Δκ` with the
default `Δκ ~ Normal+(0, 2)` (κ median ≈ 2.3, 95% ≈ 1–5): with only four
lab vintages the strength is weakly identified, so the prior is
deliberately informative and the BVD pool is tied to the trajectory (not
a free pool-size parameter). Pass `delta_kappa_prior` to override.
"""
@model function test_priority_model(;
        delta_kappa_prior = truncated(Normal(0.0, 2.0); lower = 0),
        priority_off::Bool = false)
    if priority_off
        Δκ = 0.0
        κ_priority = 1.0
    else
        Δκ ~ delta_kappa_prior
        κ_priority = 1.0 + Δκ
    end
    return (; Δκ, κ_priority)
end

"""
Test-positivity machinery. Samples
- `λ_bg` — the per-day non-BVD background suspected-case rate, on a
  half-normal scale. Underlies the suspected/confirmed contrast.
- `τ` — the fraction of suspected cases that get sampled and routed
  to the laboratory pipeline; together with the lab-delay CDF this
  handles right-truncation of the per-test positivity observation.

The default `λ_bg` prior is a half-normal
`truncated(Normal(0, 1.0); lower = 0)`. Its total contribution to the
expected suspected-case count over the window is `λ_bg · T`, where `T`
is the latent seeding-to-cut-off time (≈ 132 days on current data).
The prior is deliberately informative because `λ_bg` is degenerate
with outbreak size (the per-bin reported mean is
`p_drc · ΔμBVD0 + λ_bg · Δt`), so a diffuse prior lets the background
absorb arbitrarily many suspected cases and resolve at the high end
where the deaths and exports streams anchor `C_T`. A background-noise
process must not be able to explain more suspected cases than were
ever reported. With SD 1.0 the median background is ≈ 0.67/day (≈ 89
cases, ≈ 8% of the 1077 observed at the 26 May cut-off) and the 95%
prior bound is ≈ 2.0/day (≈ 259 cases, ≈ 24% of observed), keeping the
background a modest minority of the suspected total while still
admitting a genuine non-BVD signal. A wider SD (≈ 1.5) was tried but
left a second posterior mode in which `λ_bg` runs to ≈ 8/day and the
background explains the majority of suspected cases; SD 1.0 keeps the
fit in the regime where the BVD trajectory, not the background, drives
the suspected total. Pass `lambda_prior` to override.

The derived per-suspected positivity `μ_BVD / μ_cases` is exposed
inside [`reported_cases_model`](@ref); the per-test positivity
`s · BVD_tested / (BVD_tested + bg_tested)` is exposed inside
[`confirmed_cases_model`](@ref).

## Saturating ramp background, anchored to reporting onset

The background is a saturating ramp in time, anchored to surveillance /
reporting onset rather than seeding (issue: lab-stream data-model
mismatch). With the ramp clock starting at the reporting-onset elapsed
time `t_report` (`Δt = t − t_report`), the per-day non-BVD rate is
`λ_bg(t) = 0` for `t ≤ t_report` and
`λ_bg(t) = λ0 + Δλ·(1 − e^(−Δt/scale))` afterwards, rising from the
baseline `λ0` towards `λ0 + Δλ` over the surveillance scale-up timescale
`scale` once case-finding has begun (see [`BackgroundRamp`](@ref)). A
constant-rate background plus exponential BVD forces the per-test
positivity to *rise* over the lab vintages, but the observed cumulative
positivity *falls* (0.48 → 0.30 over 23-26 May); a background that
broadens as the suspected-case definition widens, anchored so it is
still climbing across the late lab window, pulls positivity down there
(a seeding-anchored ramp has saturated long before the cut-off).

`λ0` is the baseline rate (default half-normal `Normal+(0, 1)`, the
previous `λ_bg` prior) and `Δλ` the ramp amplitude (default half-normal
`Normal+(0, 1)`). `scale` is a fixed keyword (default
[`BACKGROUND_RAMP_SCALE`](@ref)): the four lab vintages cannot identify
it, so it is not sampled. The reporting-onset offset is supplied by the
observation submodels via [`report_onset_offset`](@ref) and the latent
`T`, so `t_report = T − offset`; this submodel returns the sampled
`λ0` / `Δλ` and the fixed `scale` for the observation models to build the
anchored [`BackgroundRamp`](@ref). The `delta_zero` flag fixes `Δλ = 0`,
which with `t_report = 0` (the default seeding-anchored ramp returned
here) recovers the constant background `λ_bg(t) = λ0`, `μ_bg(t) = λ0·t`,
exactly.

Returns the legacy scalar `λ_bg = λ0` (the baseline rate, so existing
constant-background callers keep working), the sampled `λ0` / `Δλ`, the
fixed `scale`, a seeding-anchored [`BackgroundRamp`](@ref) `bg`
(`t_report = 0`) and a `μ_bg` cumulative closure for any caller that does
not supply a reporting offset.
"""
@model function test_positivity_model(;
        lambda_prior = truncated(Normal(0.0, 1.0); lower = 0),
        ramp_amplitude_prior = truncated(Normal(0.0, 1.0); lower = 0),
        fraction_tested_prior = Beta(5.0, 2.0),
        scale::Real = BACKGROUND_RAMP_SCALE,
        delta_zero::Bool = false)
    λ0 ~ lambda_prior
    if delta_zero
        Δλ = zero(λ0)
    else
        Δλ ~ ramp_amplitude_prior
    end
    τ_test ~ fraction_tested_prior
    bg = background_ramp(λ0, Δλ, scale)
    μ_bg = t -> bg_cumulative(bg, t)
    return (; λ_bg = λ0, λ0, Δλ, scale, bg, μ_bg, τ_test)
end

"""
Laboratory daily processing-capacity prior for the constant-capacity
batch model (issue #174). Samples a single constant daily analysis
capacity `κ` (samples a laboratory can process per day) on a log scale.
The FIFO backlog of received specimens is drained at `κ` per day, so
`κ` sets the lab turnaround in place of a fixed onset-to-lab Gamma delay:
when arrivals exceed `κ` a backlog builds and analysis lags receipt,
recovering the observed mid-period stall without a separate delay
distribution.

Default `LogNormal(log(80), 0.5)`: median ≈ 80 samples/day, 95% interval
≈ 30-210. Anchored on the observed daily analysed increments across the
four 23-26 May vintages (211 → 295 → 295 → 403, i.e. ≈ 84, 0, 108 per
day; the zero is the 25 May Ituri stoppage fixed externally), so the
prior centres on the order of the realised batch size while staying wide.
Used by [`lab_throughput_model`](@ref).
"""
@model function lab_capacity_model(;
        capacity_prior = LogNormal(log(80.0), 0.5))
    κ_lab ~ capacity_prior
    return (; κ_lab)
end

"""
Case-fatality ratio prior. Default `Beta(6.6, 13.4)` has mean ≈ 0.33,
matching the CDC summary for past BVD outbreaks. Used by
[`deaths_model`](@ref) and [`exports_deaths_model`](@ref).
"""
@model function cfr_model(; cfr_prior = Beta(6.6, 13.4))
    CFR ~ cfr_prior
    return (; CFR)
end

"""
Prior on the detection window `w` — mean days during which a case is
still infectious and detectable abroad. Default centred on 15 days
(the McCabe et al. central scenario) with SD 5. Used by
[`exports_model`](@ref).
"""
@model function detection_window_model(;
        window_prior = truncated(Normal(15.0, 5.0); lower = 0))
    w ~ window_prior
    return (; w)
end

"""
Prior on the mean daily traveller volume from the source area to
Uganda. Default centred on `ITURI_DAILY_TRAVEL` with SD
`ITURI_DAILY_TRAVEL_SD`, truncated at zero. Used by
[`exports_model`](@ref).
"""
@model function traveller_volume_model(;
        mean::Real = ITURI_DAILY_TRAVEL,
        sd::Real = ITURI_DAILY_TRAVEL_SD)
    daily_travellers ~ truncated(Normal(mean, sd); lower = 0)
    return (; daily_travellers)
end

"""
Shared negative-binomial dispersion `k` for both passive-surveillance
streams (suspected deaths and reported cases). Sampled on the
`1/sqrt(k)` scale with a weakly-informative half-normal prior
following the Stan prior-choice recommendations.
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
and pooling strength `τ`. Used by [`reported_cases_model`](@ref),
[`exports_model`](@ref) and [`exports_deaths_model`](@ref); this is the
composer default.
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

"""
Per-bin random-effect DRC ascertainment for the per-vintage reported
and confirmed streams. Given the pooled hyperparameters `μ_logit`
and `τ_logit` from [`pooled_ascertainment_model`](@ref), draws `n` IID
non-centred logit-scale offsets `z_drc_t ~ Normal(0, 1)` and exposes
`p_drc_t = logistic(μ_logit + τ_logit · z_drc_t)` as a length-`n`
vector. Used by [`reported_cases_model`](@ref) and
[`confirmed_cases_model`](@ref) so each vintage bin draws its own
ascertainment fraction from the same population distribution that the
pooled scalar `p_drc` is drawn from. With `n = 1` the draw matches that
pooled scalar.

No autocorrelation is imposed: this is an IID random effect over bins,
not a random walk. Identification leans on the pooling — `τ_logit`
shrinks the per-bin draws back toward the hyperprior mean when the
data are uninformative about a particular bin.
"""
@model function daily_ascertainment_model(n::Integer,
        μ_logit::Real, τ_logit::Real)
    z_drc_t ~ filldist(Normal(0, 1), n)
    p_drc_t := logistic.(μ_logit .+ τ_logit .* z_drc_t)
    return (; z_drc_t, p_drc_t)
end

"""
Deaths-reporting ascertainment factor, allowing the observed
*suspected* deaths to drift around the BVD-driven CFR-weighted
expectation. Default `Normal(1.0, 0.05)` truncated at zero: ~95%
prior mass within 10% of unity, but the prior allows both slight
under-reporting (`p_deaths < 1`, missed BVD deaths) and slight
over-reporting (`p_deaths > 1`, non-BVD deaths captured by the
suspected case definition). The prior is judgement-based — there is
no external surveillance-completeness study for this outbreak — and
is intentionally tight so it cannot absorb the bulk of a
data-vs-model conflict; widen `sd_prior` for sensitivity. Used by
[`deaths_model`](@ref) (multiplies the expected-deaths trajectory).
"""
@model function deaths_ascertainment_model(;
        ascertainment_prior = truncated(Normal(1.0, 0.05);
        lower = 0))
    p_deaths ~ ascertainment_prior
    return (; p_deaths)
end
