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
Severe-first testing-selection prior for the laboratory confirmed stream.
The lab tests the most-likely-BVD (severest / most obvious) suspects
first, so the BVD share of the tested pool starts high and decays to a
broad-pool baseline as testing widens (see [`confirmed_cases_model`](@ref)
and [`severe_first_share`](@ref)):

```math
q_v = q_\\infty + (q_0 - q_\\infty)\\, e^{-c_v / \\text{scale}},
```

with `c_v` the elapsed time since testing (reporting) onset. This submodel
samples the two shape parameters of that curve:

- `q0` — the early severe-cluster BVD fraction (the share of the
  first-tested obvious cases that are true BVD). Its prior is near 1
  (`Beta(20, 1.5)`, mean ≈ 0.93, most mass above 0.85): the first batch
  the lab runs is the obvious-BVD cluster. With `q0 ≈ 1` the first-vintage
  positivity is `≈ s` (true positives on a near-pure-BVD batch), so the
  sensitivity `s` is identified directly from the early data and needs no
  lower floor.
- `decay_scale` — the timescale (days) over which the tested BVD share
  relaxes from `q0` toward the baseline `q∞`. Identified by how fast the
  observed positivity falls across the four vintages; a half-normal
  `Normal+(0, 10)` (median ≈ 6.7 days) spans the lab window.
- `qinf` — the baseline BVD fraction the broad suspect pool settles at once
  the severe cluster is exhausted. Its prior `Beta(6, 6)` (mean 0.5,
  central) sits near the cut-off test-positivity-implied share: with the
  late positivity ≈ 0.30 and `s ≈ 0.48`, `qinf ≈ 0.6`. The plateau
  positivity is `s·qinf + (1−spec)(1−qinf)`. Outbreak size stays pinned by
  the deaths / exports streams and by the received-count likelihood (the
  forwarded fraction of the suspect backlog), not by `qinf`. The
  count-implied composition `μ_BVD / (μ_BVD + μ_bg)` is still exposed as a
  diagnostic inside [`confirmed_cases_model`](@ref).

Pass `q0_prior` / `decay_prior` / `qinf_prior` to override.
"""
@model function test_selection_model(;
        q0_prior = Beta(20.0, 1.5),
        decay_prior = truncated(Normal(0.0, 10.0); lower = 0.0),
        qinf_prior = Beta(6.0, 6.0))
    q0 ~ q0_prior
    decay_scale ~ decay_prior
    qinf ~ qinf_prior
    return (; q0, decay_scale, qinf)
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
Beta-Binomial overdispersion prior for the per-vintage confirmed
likelihood. The cumulative-to-vintage analysed denominators carry
laboratory reporting noise (backfill, batching, lab stalls) that a plain
Binomial cannot absorb, so the per-window positivity is wildly
non-monotone (0.48, 0.05, 0.15, 0.02, 0.79 on the 28 May data). A
Beta-Binomial with concentration `φ` lets each window's positive count
disperse around its mean `ΔA_v · p_pos,v` while keeping the same expected
positivity, so the smooth positivity curve no longer has to thread every
noisy window exactly. The Binomial is recovered as `φ → ∞`.

The prior is on `φ` directly (samples per concentration); the default
`Gamma(2, 25)` (mean 50, mass on 5–150) admits both near-Binomial
(`φ` large) and strongly-overdispersed (`φ` small) windows. Pass
`concentration_prior` to override. Used by [`confirmed_cases_model`](@ref)
when `overdispersed = true`.
"""
@model function confirmed_overdispersion_model(;
        concentration_prior = Gamma(2.0, 25.0))
    φ_conf ~ concentration_prior
    return (; φ_conf)
end

"""
Per-vintage tested-BVD-share random effect for the confirmed stream. The
per-window positivity on the 28 May data is non-monotone (0.48, 0.05,
0.15, 0.02, 0.79): a monotone severe-first / volume q-curve cannot match
both the early high batch and the late 28 May surge. This submodel adds a
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
Imputed analysed-denominator latent for confirmed vintages with no
observed analysed count (the 18-22 May early and 29-31 May late lab
windows, where the national `Nbre d'échantillons analysés` is missing).
The denominator is imputed as a TIGHT log-scale latent anchored to the
observed analysed increments, so it is a smooth extrapolation of the
known 23-28 May series rather than a free dimension.

For `n` missing vintages the latent log increments are a log-random-walk
anchored at `log_anchor` (the log geometric-mean of the observed analysed
increments), `log ΔA_j = log_anchor + σ_A · cumsum(z_A)_j`, with `n` IID
standard-normal steps `z_A` and a small walk SD `σ_A` (default
`truncated(Normal(0, 0.3); lower = 0)`). The walk keeps the imputed
denominators close to the observed scale. `σ_A → 0` pins every imputed
denominator at the anchor. Used by [`confirmed_cases_model`](@ref) when
`analysed_impute` is set and some `samples_analysed` entries are missing.

NEGATIVE RESULT — DEFAULT OFF. Imputing the no-denominator confirmed
vintages this way still funnels: a four-chain joint fit over the early
18-22 May extension does not converge (R-hat ≈ 2.1, bulk ESS ≈ 1, every
NUTS transition saturates the maximum tree depth with zero divergences),
and tightening the pooling from `σ_A = 0.3` to `0.05` does not help. The
ridge is the per-vintage `q` random effect times the imputed `ΔA`: each
no-denominator vintage carries both a free positivity offset and a latent
denominator, and `ΔC ≈ ΔA · p_pos(q)` leaves them jointly unidentified.
The free analysis-capacity (`μ_A`) imputation funnels the same way. The
confirmed stream is therefore fitted only over the 23-28 May vintages
that carry an observed analysed denominator; this submodel is retained as
an opt-in experiment, off by default.
"""
@model function analysed_impute_model(n::Integer, log_anchor::Real;
        sigma_prior = truncated(Normal(0.0, 0.3); lower = 0))
    σ_A ~ sigma_prior
    z_A ~ filldist(Normal(0, 1), n)
    ## Log-random-walk around the observed-increment anchor; the cumulative
    ## sum gives a smooth drift while the tight σ_A keeps the imputed
    ## denominators on the observed scale.
    log_ΔA = log_anchor .+ σ_A .* cumsum(z_A)
    return (; σ_A, z_A, log_ΔA)
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
"""
@model function test_positivity_model(;
        lambda_prior = truncated(Normal(0.0, 1.0); lower = 0),
        fraction_forwarded_prior = Beta(5.0, 2.0))
    λ_bg ~ lambda_prior
    τ_forward ~ fraction_forwarded_prior
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
"""
@model function cfr_model(; cfr_prior = Beta(6.6, 13.4))
    CFR ~ cfr_prior
    return (; CFR)
end

"""
Lab-confirmation coverage for BVD deaths: the fraction of true BVD deaths
whose specimen reaches the laboratory and is tested. The confirmed-death
increment thins the *modelled* BVD-death trajectory by `coverage_death · s`,
where `s` is the shared confirmed-case PCR sensitivity
([`test_sensitivity_model`](@ref)). A death specimen is BVD (`q = 1`), so
its positivity is `s`, not `s · q`; `coverage_death` carries the
death-specimen submission rate. It is identified by the single point (17
confirmed against the modelled BVD-death total ≈ 246 suspected deaths)
given `s` from the cases, so no degeneracy. Default `Beta(2, 18)` is
weakly-informative favouring low coverage (mean 0.10, 95% interval ≈
0.01-0.26), reflecting that post-mortem specimen submission is sparse;
the data-implied `coverage · s ≈ 17/246` pins it near this centre rather
than the prior dominating. Used by [`confirmed_deaths_model`](@ref).
"""
@model function death_coverage_model(;
        coverage_prior = Beta(2.0, 18.0))
    coverage_death ~ coverage_prior
    return (; coverage_death)
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
