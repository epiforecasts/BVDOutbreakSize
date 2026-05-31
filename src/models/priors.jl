# Prior submodels: building blocks shared across the observation
# submodels and joint composers. Each `@model` is a small piece of the
# generative process — a single prior or a small group of related
# priors — so it can be reused in different composers without code
# duplication.

"""
Prior on the cumulative case count `C(T) = exp(r·T)` via a doublings
parameterisation. Samples the exponential growth rate `r` and a
doubling-count `m = T/τ`, then exposes `(τ, m, r, T, C_T, cumulative)`
as deterministics for downstream submodels.

McCabe et al.'s primary assumption is the doubling time (their 7/14/21-day
sweep); each doubling time implies a growth rate `r = log(2)/τ`, and the
prior is placed on that implied `r`, with `τ = log(2)/r` recovered as a
deterministic. The default `r ~ LogNormal(log(log(2)/14), 0.4)` is
exactly equivalent to a `LogNormal(log(14), 0.4)` prior on the doubling
time: because `r = log(2)/τ` is a reciprocal, the log-scale SD `0.4` is
preserved, so the implied doubling-time prior (and every derived
quantity) matches that prior on `τ`.

The default doubling-count prior `m ~ Normal(9, 3)` (truncated at 0) is
centred on `m = 9` (`C_T = 2^9 = 512`), the doubling count implied by
McCabe et al.'s first-report (18 May 2026) Method 2 central scenario of
501 cases (`log2(501) ≈ 9`). For a fit at a later cut-off, pass an
`m_prior` whose centre advances from that base date via
[`m_prior_centre`](@ref) so the prior tracks the elapsed time. This is a
weakly-informative centring choice (SD 3 gives 95% support ≈ `m ∈ (3, 15)`,
`C_T ∈ (8, 32000)`); the fit is dominated by the likelihood, so it mainly
sets where the joint sampler starts.
"""
@model function exponential_growth_model(;
        r_prior = LogNormal(log(log(2) / 14), 0.4),
        m_prior = truncated(Normal(9.0, 3.0); lower = 0))
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
PCR sensitivity prior for the GeneXpert Ebola assay. Beta(30, 2): mean
0.94, 95% interval 0.84-0.99. Sits just below the 100% (95% CI
84.6-100%, n = 22) clinical sensitivity on field whole blood reported
in the Sierra Leone Zaire-ebolavirus field evaluation, leaving room for
early-infection low-viral-load specimens, field handling, and the lack
of Bundibugyo-specific validations. Used by
[`confirmed_cases_model`](@ref).
"""
@model function test_sensitivity_model(;
        sensitivity_prior = Beta(30.0, 2.0))
    s_test ~ sensitivity_prior
    return (; s_test)
end

"""
Test-positivity machinery. Samples
- `λ_bg` — the per-day non-BVD background suspected-case rate, on a
  half-normal scale. Underlies the suspected/confirmed contrast.
- `τ` — the fraction of suspected cases that get sampled and routed
  to the laboratory pipeline; together with the lab-delay CDF this
  handles right-truncation of the per-test positivity observation.

The derived per-suspected positivity `μ_BVD / μ_cases` is exposed
inside [`reported_cases_model`](@ref); the per-test positivity
`s · BVD_tested / (BVD_tested + bg_tested)` is exposed inside
[`confirmed_cases_model`](@ref).
"""
@model function test_positivity_model(;
        lambda_prior = truncated(Normal(0.0, 10.0); lower = 0),
        fraction_tested_prior = Beta(5.0, 2.0))
    λ_bg ~ lambda_prior
    τ_test ~ fraction_tested_prior
    return (; λ_bg, τ_test)
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
Partially pooled ascertainment fractions for the DRC and Uganda
surveillance systems, sampled in non-centred form to avoid the funnel
geometry. Both logit-scale fractions share a hyperprior with mean `μ`
and pooling strength `τ`. Used by [`reported_cases_model`](@ref),
[`exports_model`](@ref) and [`exports_deaths_model`](@ref).
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
