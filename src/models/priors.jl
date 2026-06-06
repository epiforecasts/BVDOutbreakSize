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

The prior is placed on the doubling time (the primary epidemiological
assumption); each doubling time implies a growth rate `r = log(2)/τ`, so
the prior is sampled on that implied `r`, with `τ = log(2)/r` recovered as
a deterministic. The default centres on the molecular-clock doubling-time
estimate for this outbreak: Cuomo-Dannenburg & Ghafari's phylodynamic
reanalysis of the first ten BDBV genomes (cuomodannenburg2026) puts the
mean doubling time at 15.2-24.5 days across six substitution-rate
assumptions, so the centre is
set to 20 days (`M_PRIOR_DOUBLING_DAYS`), the middle of that range, slower
than McCabe et al.'s 14-day central scenario. The log-SD is set to `0.15`:
reading the 15.2-24.5 day spread of mean doubling times as roughly a 95%
interval implies a log-SD near `0.12`, and `0.15` inflates that a little to
allow for the wide per-assumption intervals without drifting off the
molecular-clock estimate. The default `r ~ LogNormal(log(log(2)/20), 0.15)`
is exactly equivalent to a `LogNormal(log(20), 0.15)` prior on the doubling
time, because `r = log(2)/τ` is a reciprocal that preserves the log-scale
SD. The implied 95% doubling-time interval is roughly `τ ∈ (14.9, 26.8)`
days, sitting on the reported range with a light inflation. The McCabe et
al. Method 2 reproduction instead pins the growth rate to their 14-day
central doubling time (see the report). Fits remain likelihood-dominated,
so the
situation-report trajectory still drives the posterior growth rate.

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
        r_prior = LogNormal(log(log(2) / M_PRIOR_DOUBLING_DAYS), 0.15),
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
PCR sensitivity prior. `Beta(6, 2)` (mean 0.75, 95% interval 0.39-0.97),
untruncated. Confirmation runs on the altona RealStar
Filovirus Screen RT-PCR, which detects Bundibugyo virus at 11-67 RNA
copies per reaction; the rapid Cepheid GeneXpert Ebola assay is
Zaire-ebolavirus-specific and does not reliably detect Bundibugyo. The
Beta keeps good analytical sensitivity plausible while carrying downside
mass for early low-viral-load specimens and field handling.

The Beta keeps good analytical sensitivity plausible while carrying
downside mass for early low-viral-load specimens. Under the severe-first
backlog model the first vintage's analysed batch is near-pure BVD
(`q ≈ 1` when selection is strong), so the v1 positivity ≈ `s` identifies
the sensitivity directly from the early data; no lower truncation is
imposed. Pass `sensitivity_prior` to override. Used by
[`confirmed_cases_model`](@ref).
"""
@model function test_sensitivity_model(;
        sensitivity_prior = Beta(6.0, 2.0))
    s_test ~ sensitivity_prior
    return (; s_test)
end

"""
Severity-enrichment prior for the composition-linked confirmed positivity
(see [`confirmed_cases_model`](@ref)). The tested BVD share is the
suspect-pool composition `φ = μ_BVD / (μ_BVD + μ_bg)` UPSAMPLED by a
severity-enrichment that decays as testing widens:

```math
\\mathrm{logit}(q_v) = \\mathrm{logit}(\\varphi_v)
    + \\delta_0\\, e^{-c_v / \\text{decay}},
```

so the lab over-tests BVD early (severe cases are triaged first and are
more likely BVD), the enrichment `δ₀·e^{−c/decay}` relaxing toward zero as
testing widens, at which point the tested share equals the pool composition.
This ties positivity to `μ_bg`, so the confirmed/positivity data identify
the non-BVD background `λ_bg` rather than it being absorbed by a free
selection curve.

`δ₀` is the early severity log-odds enrichment of BVD; lower-truncated at 0
because severity triage upsamples BVD, never down. The default
`truncated(Normal(1.5, 0.75); lower = 0)` is deliberately moderate /
bounded: even severity-triaged testing cannot be near-pure BVD (other
haemorrhagic / severe febrile illness is also triaged), so for a pool
composition `φ ≈ 0.4` the early tested share is `logistic(logit(0.4) +
1.5) ≈ 0.75`, not the near-1 of a radical free curve. `decay_scale`
is the relaxation timescale on the analysed-volume clock. Pass
`logodds_prior` / `decay_prior` to override. Used by
[`confirmed_cases_model`](@ref).
"""
@model function severity_enrichment_model(;
        logodds_prior = truncated(Normal(1.5, 0.75); lower = 0),
        decay_prior = truncated(Normal(0.0, 10.0); lower = 0.0))
    δ0 ~ logodds_prior
    decay_scale ~ decay_prior
    return (; δ0, decay_scale)
end

"""
PCR specificity prior. The per-test positivity of the confirmed stream is
`p_pos = s·q + (1−spec)·(1−q)`, where `q` is the BVD share of the
analysed batch (see [`confirmed_cases_model`](@ref)): true positives at
sensitivity `s` plus false positives at false-positive rate `1−spec` on
the non-BVD share. Modelling specificity lets the observed positivity
exceed `s·q` (it can include false positives), so the observed peak no
longer pins a hard lower bound on the sensitivity, replacing the earlier
`s`-floor.

The default `Beta(50, 1.5)` is informative and high — mean ≈ 0.97, with
most mass above 0.95 and a modest lower tail — reflecting the altona
RealStar RT-PCR confirmation assay's strong analytical specificity. The
prior is deliberately tight: a loose specificity prior lets the
false-positive term absorb the positivity signal and reintroduces the
low-sensitivity / small-outbreak multimodality (issue #151) the `s`-floor
had suppressed. Pass `specificity_prior` to override (e.g. a looser Beta
for a sensitivity analysis). Used by [`confirmed_cases_model`](@ref).
"""
@model function test_specificity_model(;
        specificity_prior = Beta(50.0, 1.5))
    spec_test ~ specificity_prior
    return (; spec_test)
end

"""
Per-vintage tested-BVD-share random effect for the confirmed stream. The
per-window positivity on the 28 May data is non-monotone (0.48, 0.05,
0.15, 0.02, 0.79): the smooth composition baseline cannot match both the
early high batch and the late 28 May surge. This submodel adds a
partially-pooled logit-scale offset to the smooth baseline share, so each
vintage's tested BVD fraction `q_v = logistic(logit(q_base,v) + σ_q·z_v)`
can fit its own positivity while sharing strength through the pooling
scale `σ_q`.

The sensitivity `s` stays fixed (an assay property); the positivity wobble
is absorbed entirely by the mixture share `q` (the fraction of the batch
that is truly BVD), honouring the principle that positivity moves through
who is tested, not through `s`. `n` IID standard-normal offsets `z_q` and
the pooling SD `σ_q ~ truncated(Normal(0, 1); lower = 0)` (logit-scale)
are sampled; `σ_q → 0` recovers the smooth baseline curve. Used by
[`confirmed_cases_model`](@ref) when `q_random_effect = true`.
"""
@model function confirmed_q_re_model(n::Integer;
        sigma_prior = truncated(Normal(0.0, 1.0); lower = 0))
    σ_q ~ sigma_prior
    z_q ~ filldist(Normal(0, 1), n)
    return (; σ_q, z_q)
end

"""
Epidemiological-exclusion fraction `e` for the laboratory-throughput
queue (see [`confirmed_cases_model`](@ref)).
`e` is the share of suspected cases ruled out by epidemiological
follow-up and never sampled, so the cumulative received backlog
asymptotes to `(1 − e)·N_susp` rather than the whole suspect total. It
gives the received-count likelihood the correct ceiling.

`e` is NOT identified by the six observed received/analysed points: the
sitrep `retrait des non-cas` mixes tested-negatives with the
epidemiologically excluded, so no published figure pins it. It is
therefore a weakly-informative prior-only scalar. The default
`Beta(2, 12)` has mean ≈ 0.14 and is used as the sensitivity arm; the
headline fit pins `e = 0` (forward fraction 1) by passing
`epi_exclusion = nothing`. Returns `(; e, forward = 1 − e)`.
"""
@model function epi_exclusion_model(; e_prior = Beta(2.0, 12.0))
    e ~ e_prior
    return (; e, forward = one(e) - e)
end

"""
Test-positivity machinery. Samples
- `λ_bg` — the per-day non-BVD background suspected-case rate, on a
  half-normal scale. Underlies the suspected/confirmed contrast; the
  cumulative background is the constant-rate `μ_bg(t) = λ_bg · t`.
- `τ_forward` — the fraction of suspected cases forwarded to the
  laboratory. It scales the cumulative suspect backlog into the expected
  received count (`received_v ~ NegBinomial(τ_forward · N_susp,v, k)`, see
  [`confirmed_cases_model`](@ref)), so the received-count stream pins it
  directly from received-versus-suspected.

The default `λ_bg` prior is a half-normal
`truncated(Normal(0, 1.0); lower = 0)`. Its total contribution to the
expected suspected-case count over the window is `λ_bg · T`, where `T`
is the latent seeding-to-cut-off time (≈ 132 days on current data).
The prior is deliberately informative because `λ_bg` is degenerate
with outbreak size (the per-bin reported mean is
`p_drc · ΔμBVD0 + λ_bg · Δt`), so a diffuse prior lets the background
absorb arbitrarily many suspected cases and resolve at the high end
where the deaths and exports streams anchor `C_T`. With SD 1.0 the median
background is ≈ 0.67/day (≈ 89 cases, ≈ 8% of the 1077 observed at the
26 May cut-off) and the 95% prior bound is ≈ 2.0/day (≈ 259 cases),
keeping the background a modest minority of the suspected total while
still admitting a genuine non-BVD signal. Pass `lambda_prior` to override.

`τ_forward` has a `Beta(5, 2)` default (mean ≈ 0.71). The derived
per-suspected positivity `μ_BVD / μ_cases` is exposed inside
[`reported_cases_model`](@ref); the per-test positivity and the tested
BVD share are exposed inside [`confirmed_cases_model`](@ref). Pass
`fraction_forwarded_prior` to override.

Set `sample_forward = false` to fix `τ_forward = 1` without sampling it.
The queue path (`confirmed_queue = true`) uses this because forwarding is
governed by `1 − e` from [`epi_exclusion_model`](@ref), leaving
`τ_forward` a dead dimension that would otherwise clutter the posterior.
"""
@model function test_positivity_model(;
        lambda_prior = truncated(Normal(0.0, 1.0); lower = 0),
        fraction_forwarded_prior = Beta(5.0, 2.0),
        sample_forward::Bool = true)
    λ_bg ~ lambda_prior
    if sample_forward
        τ_forward ~ fraction_forwarded_prior
    else
        ## Queue path: forwarding is governed by `1 − e` from the
        ## exclusion submodel, so `τ_forward` is a fixed constant and is
        ## not sampled (avoids a dead posterior dimension).
        τ_forward = one(λ_bg)
    end
    return (; λ_bg, τ_forward)
end

"""
Report-to-lab-receipt delay prior. Samples a Gamma shape `α` and scale
`θ` for the time from a suspected-case report to specimen receipt at the
Bunia/INRB laboratory (transport and intake), centred near 3 days. The
confirmed model convolves the suspect backlog with this delay to get the
received queue. Prior needs human sign-off against the sitrep turnaround.
Used by [`confirmed_cases_model`](@ref).
"""
@model function lab_receipt_delay_model(;
        alpha_prior = truncated(Normal(2.0, 1.0); lower = 0.1),
        theta_prior = truncated(Normal(1.5, 0.75); lower = 0.1))
    α_recv ~ alpha_prior
    θ_recv ~ theta_prior
    return (; α = α_recv, θ = θ_recv, dist = Gamma(α_recv, θ_recv))
end

"""
Laboratory analysis-capacity prior over the `n` confirmed vintages. The
per-window daily capacity follows a log random walk centred on
`capacity_centre` samples analysed per day (external DRC sitrep / NYT
figure, ~150): `log κ_1 ~ Normal(log centre, level_sd)` and
`log κ_v = log κ_{v−1} + rw_sd·z_v` with non-centred innovations `z`. The
throughput likelihood in [`confirmed_cases_model`](@ref) scales the
available received backlog by `1 − exp(−κ_v·Δt_v / backlog)`, so capacity
limits analysis when the backlog is large and the whole backlog is
processed when capacity is ample. Returns the length-`n` `capacity`.
"""
@model function lab_capacity_model(n::Integer;
        capacity_centre::Real = 150.0,
        level_sd::Real = 0.5,
        rw_sd_prior = truncated(Normal(0.0, 0.5); lower = 0))
    z_cap ~ filldist(Normal(0, 1), n)
    rw_sd ~ rw_sd_prior
    ## Non-centred log random walk: window 1 at the centred level, later
    ## windows accumulate rw_sd-scaled innovations.
    steps = cumsum(vcat(zero(eltype(z_cap)), z_cap[2:end]))
    logκ = log(capacity_centre) .+ level_sd * z_cap[1] .+ rw_sd .* steps
    capacity := exp.(logκ)
    return (; capacity, rw_sd)
end

"""
Case-fatality ratio prior. Default `Beta(6.6, 13.4)` has mean ≈ 0.33,
matching the CDC summary for past BVD outbreaks. Used by
[`deaths_model`](@ref) and [`exports_deaths_model`](@ref).

The CFR multiplies the latent *infection* trajectory: with no asymptomatic
fraction and no case-ascertainment on the death denominator it is applied
as deaths per infection, coinciding with the infection-fatality ratio (IFR)
under those assumptions. The conventional CFR label is kept.
"""
@model function cfr_model(; cfr_prior = Beta(6.6, 13.4))
    CFR ~ cfr_prior
    return (; CFR)
end

"""
Non-BVD background rate for the suspected-death stream. The DRC sitrep
suspected deaths are a noisy passive-surveillance count: some are true
BVD deaths, some are non-BVD deaths swept in by the broad suspect case
definition. This submodel samples the per-day non-BVD background death
rate `λ_bg_death`, the death analogue of the case background `λ_bg`
([`test_positivity_model`](@ref)): the cumulative background is the
constant-rate `μ_bg_death(t) = λ_bg_death · t`, added to the BVD-driven
CFR-weighted expectation inside [`deaths_model`](@ref).

The default prior is a half-normal `truncated(Normal(0, 0.25); lower = 0)`.
Its total contribution to the expected suspected-death count over the
window is `λ_bg_death · T` (`T` ≈ 132 days on current data). The prior is
deliberately informative, mirroring `λ_bg` for cases: a diffuse background
is degenerate with outbreak size (the per-bin suspected-death mean is
`Δμ_BVD_death + λ_bg_death · Δt`), so the prior keeps the background a
modest minority of the suspected-death total. With SD 0.25 the median
background is ≈ 0.17/day (≈ 22 deaths, ≈ 9% of the 246 suspected deaths at
the 28 May cut-off) and the 95% prior bound is ≈ 0.49/day (≈ 65 deaths),
admitting a genuine non-BVD signal while leaving the bulk to BVD. A
time-varying background is a wanted follow-up (issue #194); this
constant-rate version is identifiable and is the death analogue of the
case `λ_bg`. Pass `lambda_prior` to override. Used by
[`deaths_model`](@ref) and [`confirmed_deaths_model`](@ref).
"""
@model function death_background_model(;
        lambda_prior = truncated(Normal(0.0, 0.25); lower = 0))
    λ_bg_death ~ lambda_prior
    return (; λ_bg_death)
end

"""
Death-specimen forwarding fraction `τ_death`: the fraction of the
suspect-death backlog whose post-mortem specimen reaches the laboratory
and is analysed. It is the death analogue of the case forwarding fraction
`τ_forward` ([`test_positivity_model`](@ref)). The confirmed-death
increment is a genuine lab/positivity process on those analysed death
specimens (see [`confirmed_deaths_model`](@ref)):

```math
\\Delta D_{conf,v} \\sim
    \\mathrm{NegBinomial}(
        \\tau_{death}\\cdot p_{pos,death,v}\\cdot\\Delta N_{death,v},\\ k),
\\qquad
p_{pos,death,v} = s\\, q_{death,v} + (1 - \\text{spec})(1 - q_{death,v}),
```

with `ΔN_death,v` the suspect-death backlog increment (BVD plus
background), `q_death,v = μ_BVD_death / N_death_susp` the BVD share of that
pool, and `s`, `spec` the PCR sensitivity / specificity *shared* with the
confirmed-case lab pipeline. `τ_death` carries the forwarding rate; the
BVD-share signal lives in the positivity, not in `τ_death`, so it is not
forced to a boundary the way the old `coverage_death` thinning was.

Default `Beta(2, 8)` is weakly-informative favouring low forwarding (mean
0.20, 95% interval ≈ 0.03-0.48): post-mortem specimen submission is
sparse. It is identified by the 17 confirmed deaths against the
positivity-weighted suspect-death backlog. Pass `fraction_prior` to
override. Used by [`confirmed_deaths_model`](@ref).
"""
@model function death_forward_model(;
        fraction_prior = Beta(2.0, 8.0))
    τ_death ~ fraction_prior
    return (; τ_death)
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
Separate negative-binomial dispersion for the suspected (reported) stream
only, so that stream can be down-weighted (loosened) without changing the
deaths / received dispersion `k`. Identical structure to
[`surveillance_dispersion_model`](@ref) but with distinctly-named sampled
and deterministic variables (`inv_sqrt_k_rep`, `k_rep`), so the two
dispersions can coexist in the joint composer without a `:=` name clash.
A wider default prior (`truncated(Normal(1.2, 0.4); lower = 0)`, smaller
`k`, heavier-tailed NegBinomial) down-weights the suspected counts so the
confirmed stream drives the shared parameters; pass `inv_sqrt_k_prior` to
override. Used by [`bvd_joint`](@ref) via its `reported_dispersion`
keyword.
"""
@model function reported_dispersion_model(;
        inv_sqrt_k_prior = truncated(Normal(1.2, 0.4); lower = 0))
    inv_sqrt_k_rep ~ inv_sqrt_k_prior
    k_rep := 1.0 / (inv_sqrt_k_rep^2 + eps(typeof(inv_sqrt_k_rep)))
    return (; k = k_rep, inv_sqrt_k = inv_sqrt_k_rep)
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
fraction of 0.75 with SD 0.6, matching the pooled composer default. This
de-pooled variant is a sensitivity alternative to
[`pooled_ascertainment_model`](@ref) (the headline default) and is not
used by the joint composer. Pass `drc_prior` / `uganda_prior` to set the
two systems' priors separately.
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
composer default.
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
