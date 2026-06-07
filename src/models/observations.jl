# Observation submodels: each ties one data stream to the latent
# cumulative case count `C(T)`. They consume the growth state and any
# shared nuisance parameters (dispersion, ascertainment) from the
# composer above, introduce only the priors they need on top, and add
# a single likelihood term.

"""
NaN / Inf-safe `NegativeBinomial` constructor parameterised by mean
`μ` and dispersion `k`, with clamping on the success probability so
extreme NUTS proposals during warmup do not trip the distribution
domain check. Shared by [`deaths_model`](@ref),
[`reported_cases_model`](@ref) and [`confirmed_cases_model`](@ref).
"""
function safe_nbinomial(k, μ)
    p_raw = k / (k + max(μ, eps(typeof(μ))))
    p = isfinite(p_raw) ?
        clamp(p_raw, eps(typeof(k)), one(k) - eps(typeof(k))) :
        eps(typeof(k))
    return NegativeBinomial(k, p)
end

"""
Exports likelihood (Method 1, geographic spread). Couples the
observed exported-case count to `C(T)` through the at-risk
person-time integral over the detection window, with a Poisson
likelihood. Samples the detection window and traveller volume
submodels internally.
"""
@model function exports_model(
        exported_cases::Union{Missing, Integer},
        growth_state, p_uganda::Real;
        source_population::Real = ITURI_POPULATION,
        window = detection_window_model(),
        traveller = traveller_volume_model(),
        onset_fraction::Real = 1.0)
    cumulative = growth_state.cumulative
    T = growth_state.T

    window_state ~ to_submodel(window, false)
    w = window_state.w

    travel_state ~ to_submodel(traveller, false)
    daily_travellers = travel_state.daily_travellers

    q = daily_travellers / source_population
    ## `onset_fraction` defaults to 1.0 so the McCabe-style window path is
    ## unchanged; pass the incubation mgf to onset-rescale the at-risk
    ## person-time if the window is interpreted as onset→detection.
    expected_exports_T := onset_fraction *
                          expected_exports(cumulative, p_uganda, q, T, w;
        r = growth_state.r)

    exported_cases ~ Poisson(expected_exports_T)

    return (; w, daily_travellers, p_uganda,
        expected_exports = expected_exports_T)
end

"""
Exports likelihood with an explicit infection→detection delay. The
delay-convolution counterpart of [`exports_model`](@ref): the at-risk
export person-time is the infection→detection delay-survival integral
[`expected_exports_delay`](@ref) rather than the McCabe et al.
rectangular detection window. The exports stream is travel-gated, so the
at-risk clock starts at infection (the traveller moves during
incubation, pre-symptomatic). The infection→detection delay `f_det` is
passed in rather than sampled here: the export count is a single datum
and cannot identify its own delay, and entering surveillance is taken to
be the same process abroad as in the DRC, so the composers pass the
incubation ⊕ onset-to-report delay (see [`combined_delay`](@ref),
[`report_delay_model`](@ref) and [`reported_cases_model`](@ref)), letting
the reported and incubation streams pin it. There is no separate
incubation rescale here: the incubation period already sits inside
`f_det`. Samples the traveller volume submodel internally and exposes
`f_det` and `daily_travellers` for the export-deaths and detection-timing
delay submodels. Couples to `C(T)` through a Poisson likelihood.
"""
@model function exports_delay_model(
        exported_cases::Union{Missing, Integer},
        growth_state, p_uganda::Real, f_det;
        source_population::Real = ITURI_POPULATION,
        traveller = traveller_volume_model(),
        last_offset::Real = 0)
    r = growth_state.r
    T = growth_state.T

    travel_state ~ to_submodel(traveller, false)
    daily_travellers = travel_state.daily_travellers

    q = daily_travellers / source_population
    ## The exports stopped accruing `last_offset` days before the cut-off
    ## (the last import), so the at-risk export integral runs only to
    ## `t_last = T - last_offset`. The latent size `growth_state.C_T` is
    ## still reported at the cut-off `T` by the caller.
    t_last = T - last_offset
    expected_exports_T := expected_exports_delay(r, p_uganda, q, t_last,
        f_det)

    exported_cases ~ Poisson(expected_exports_T)

    return (; f_det, daily_travellers, p_uganda,
        expected_exports = expected_exports_T)
end

"""
Dated exports likelihood with an explicit infection→detection delay. The
time-resolved counterpart of [`exports_delay_model`](@ref): instead of a
single cumulative count at the cut-off, the observed Uganda exports are a
daily series `exported_cases_daily` (index 1 = earliest detection day,
last index = the cut-off day), modelled as an inhomogeneous Poisson
process. The cumulative export expectation
[`expected_exports_delay`](@ref) is placed onto the daily grid through the
shared between-edge differencing [`daily_increment_kernel`](@ref), the
same construction as [`exports_deaths_delay_model`](@ref): a continuous
survival weight over the pre-detection stretch followed by per-day Poisson
counts.

The pre-detection survival term `pre_detection_exports ~ Poisson(Λ(T-n))`
observed at zero subsumes the first-export-detection timing bound of
[`exports_detection_timing_delay_model`](@ref): it asserts that no export
was expected before the earliest observed detection. The total expected
(pre-detection weight ⊕ daily increments) equals the cumulative
expectation `Λ(T)` the scalar single-total likelihood uses, so the daily
model conditions on the same total while also using the detection timing.
A single-day series (`length 1`) reduces to that cumulative split.

The exports stream is travel-gated, so the at-risk clock starts at
infection; `f_det` (the incubation ⊕ onset-to-report delay, see
[`combined_delay`](@ref)) and the traveller volume submodel carry the
delay and per-day per-capita travel rate, exactly as in
[`exports_delay_model`](@ref). Exposes `f_det` and `daily_travellers` for
the export-deaths delay submodel.

`last_offset` (default 0) places the last series day `last_offset` days
before the cut-off rather than at the cut-off. Because the export stream
is travel-gated and cross-border movement patterns likely shift over the
outbreak (and the most recent days are right-truncated by reporting lag),
the stream is run only up to the most recent reported import to Uganda;
the at-risk clock, per-day grid and reported `expected_exports` all anchor
on that day (`t_last = T - last_offset`) instead of the cut-off.
"""
@model function exports_daily_delay_model(
        exported_cases_daily::AbstractVector,
        growth_state, p_uganda::Real, f_det;
        pre_detection_exports::Union{Missing, Integer} = 0,
        last_offset::Real = 0,
        source_population::Real = ITURI_POPULATION,
        traveller = traveller_volume_model())
    r = growth_state.r
    T = growth_state.T

    travel_state ~ to_submodel(traveller, false)
    daily_travellers = travel_state.daily_travellers

    q = daily_travellers / source_population
    n = length(exported_cases_daily)   # earliest detection to last series day

    ## The series ends `last_offset` days before the cut-off (the date of
    ## the most recent reported import to Uganda). The travel-gated export
    ## stream is only run up to that day: cross-border movement patterns
    ## likely shift over the outbreak and the days after the last import
    ## are right-truncated by reporting lag, so they carry no informative
    ## zero. The at-risk clock and the per-day grid anchor on this date.
    t_last = T - last_offset

    ## One at-risk trajectory over [0, t_last], reused across every edge
    ## (the pre-detection edge at i = 0 plus the n daily edges), the
    ## travel-gated analogue of the DRC `DailyBVDTrajectory` convolution.
    traj = ExportRiskTrajectory(t_last, r)
    edges = [t_last - n + i for i in 0:n]
    Λ = p_uganda .* q .* export_at_risk(traj, edges, f_det)

    ## Pre-detection zero stretch as one Poisson observed at 0; `missing`
    ## generates it for predictive checks. This is the first-detection
    ## timing bound: no export expected before the earliest detection.
    pre = max(Λ[1], zero(eltype(Λ)))
    pre_detection_exports ~ Poisson(pre)

    ## Per-day means via the shared between-edge differencing, starting
    ## from the pre-detection survival weight `pre`.
    μ_day = daily_increment_kernel(Λ[2:end], pre)
    for i in 1:n
        exported_cases_daily[i] ~ Poisson(μ_day[i])
    end

    expected_exports_T := max(Λ[end], eps(eltype(Λ)))
    return (; f_det, daily_travellers, p_uganda,
        expected_exports = expected_exports_T)
end

"""
Deaths likelihood (Method 2, back-calculation from deaths). Couples the
observed DRC suspected deaths to `C(T)` through the CFR-weighted gamma
convolution. The observation is a per-vintage vector `total_deaths`
aligned with `t_edges` (elapsed time since seeding to each sitrep date,
ascending, latest = `T`); the model fits the between-vintage
cumulative-count increments as independent NegBinomial terms sharing the
dispersion `k` of [`surveillance_dispersion_model`](@ref). A single
observation (`length 1`, `t_edges = [T]`) reduces to the cumulative
single-total NegBinomial likelihood, matching the McCabe et al. Method 2
configuration. Samples the [`delay_model`](@ref) and [`cfr_model`](@ref)
submodels internally.

`p_deaths` multiplies the expected-deaths trajectory to allow the
observed *suspected* deaths to drift around the BVD-driven CFR-weighted
expectation; pass it from [`deaths_ascertainment_model`](@ref) at the
joint-composer level. Defaults to `1.0` so the single-stream paths
reduce to the original likelihood.

`λ_bg_death` adds a constant-rate non-BVD background to the suspected
deaths, the death analogue of the case background `λ_bg`
([`death_background_model`](@ref), [`reported_cases_model`](@ref)): the
cumulative background is `μ_bg_death(s) = λ_bg_death · s`, so each bin's
background increment is `λ_bg_death · Δt`. The per-bin mean is then
`Δμ_BVD_death + Δμ_bg`. It defaults to `0.0`, so the single-stream and
McCabe Method 2 paths reduce to the pure BVD likelihood. The BVD-only and
total (BVD plus background) cumulative trajectories are returned so
[`confirmed_deaths_model`](@ref) can build the death-specimen positivity
from the suspect-death composition.
"""
@model function deaths_model(
        total_deaths::AbstractVector,
        growth_state, k::Real, t_edges::AbstractVector;
        delay = delay_model(),
        cfr = cfr_model(),
        p_deaths::Real = 1.0,
        λ_bg_death::Real = 0.0,
        onset_fraction::Real = 1.0)
    r = growth_state.r

    delay_state ~ to_submodel(delay, false)
    cfr_state ~ to_submodel(cfr, false)
    f_death = delay_state.dist
    CFR = cfr_state.CFR

    n = length(t_edges)
    length(total_deaths) == n ||
        error("total_deaths length must match t_edges (got " *
              "$(length(total_deaths)) vs $n)")

    ## CFR-weighted cumulative expected BVD deaths at each bin edge. The
    ## latent trajectory is cumulative *infections*, so `onset_fraction`
    ## (the incubation mgf) maps it onto onsets before the onset-to-death
    ## convolution; the drift factor `p_deaths` scales the whole
    ## trajectory.
    μ_BVD_at_edges = [s <= zero(s) ? zero(s) :
                      p_deaths * onset_fraction *
                      delay_convolution(CFR, r, s, f_death)
                      for s in t_edges]

    ## Constant-rate non-BVD background, the death analogue of the case
    ## `λ_bg`: cumulative μ_bg_death(s) = λ_bg_death·s, so the suspect-death
    ## cumulative is the BVD trajectory plus this linear background. The
    ## per-bin NegBinomial mean is the between-edge increment of the total
    ## (with a NaN / Inf-safe positive clamp via `daily_increment_kernel`).
    Tt = typeof(μ_BVD_at_edges[1] * float(λ_bg_death))
    Λ_at_edges = Vector{Tt}(undef, n)
    for i in 1:n
        s_k = t_edges[i]
        μ_bg_k = max(convert(Tt, λ_bg_death) * s_k, zero(Tt))
        Λ_at_edges[i] = convert(Tt, μ_BVD_at_edges[i]) + μ_bg_k
    end
    bin_means = daily_increment_kernel(Λ_at_edges)

    for i in 1:n
        total_deaths[i] ~ safe_nbinomial(k, bin_means[i])
    end

    expected_deaths_T := Λ_at_edges[end]

    return (; CFR, p_deaths, λ_bg_death, delay_dist = f_death,
        expected_deaths_T, Λ_at_edges, μ_BVD_at_edges)
end

"""
Reported (suspected) cases likelihood. Couples the observed DRC
suspected-case counts to `C(T)` as the sum of (i) a BVD-driven
contribution `p_drc · onset_fraction · ∫₀^s exp(r·u) · f_rep(s-u) du` using
the onset-to-report delay [`report_delay_model`](@ref) and (ii) a non-BVD
background at the constant rate `λ_bg` (see [`test_positivity_model`](@ref)),
so each bin's background increment is the cumulative difference
`λ_bg·s_k − λ_bg·s_{k−1} = λ_bg · Δt`. The latent trajectory is cumulative
*infections*, so `onset_fraction` (the incubation mgf) maps it onto
onsets before the onset-to-report convolution; the background term is
non-BVD and unscaled.

The observation is a per-vintage vector `reported_cases` aligned with
`t_edges` (elapsed time since seeding to each sitrep date, ascending,
latest = `T`); the model fits the between-vintage cumulative-count
increments as independent NegBinomial terms sharing `k` with
[`deaths_model`](@ref) and [`confirmed_cases_model`](@ref). Each bin
carries its own ascertainment `p_drc_per_bin[v]` (a per-bin random
effect from [`daily_ascertainment_model`](@ref)), applied to the
between-edge BVD increment. The first bin's "increment" is the
cumulative at the first edge, so a single observation (`length 1`,
`t_edges = [T]`) reduces to the cumulative single-total likelihood used
by McCabe et al.

The unit-ascertainment BVD cumulative `μ_BVD,0(s)` is evaluated at every
edge through the Gamma closed-form [`delay_convolution`](@ref). Returns
`report_delay_dist = f_rep` so [`confirmed_cases_model`](@ref) can reuse
the same onset-to-report kernel, and the derived per-suspected
positivity `μ_BVD / μ_cases` at the cut-off as a diagnostic.
"""
@model function reported_cases_model(
        reported_cases::AbstractVector,
        growth_state, k::Real,
        p_drc_per_bin::AbstractVector, t_edges::AbstractVector;
        report_delay = report_delay_model(),
        test_positivity = test_positivity_model(),
        report_onset_offset::Union{Nothing, Real} = nothing,
        onset_fraction::Real = 1.0)
    report_state ~ to_submodel(report_delay, false)
    test_positivity_state ~ to_submodel(test_positivity, false)
    λ_bg = test_positivity_state.λ_bg
    τ_forward = test_positivity_state.τ_forward
    f_rep = report_state.dist
    r = growth_state.r

    n = length(t_edges)
    n == length(p_drc_per_bin) ||
        error("p_drc_per_bin length must match t_edges (got " *
              "$(length(p_drc_per_bin)) vs $n)")
    n == length(reported_cases) ||
        error("reported_cases length must match t_edges (got " *
              "$(length(reported_cases)) vs $n)")

    ## Unit-ascertainment BVD cumulative onsets at each bin edge (Gamma
    ## closed form), the infection trajectory mapped onto onsets by
    ## `onset_fraction` (the incubation mgf). Per-bin ascertainment is
    ## applied below on the between-edge increment so each bin's mean
    ## tracks its own random-effect draw.
    μ_BVD0_at_edges = [s <= zero(s) ? zero(s) :
                       onset_fraction *
                       delay_convolution(one(eltype(p_drc_per_bin)),
                           r, s, f_rep) for s in t_edges]
    ΔμBVD0 = daily_increment_kernel(μ_BVD0_at_edges)

    Tt = eltype(ΔμBVD0)
    Λ_at_edges = Vector{Tt}(undef, n)
    bin_means = Vector{Tt}(undef, n)
    Λ_prev = zero(Tt)
    μ_bg_prev = zero(Tt)
    for i in 1:n
        s_k = t_edges[i]
        ## Constant-rate non-BVD background: cumulative μ_bg(s_k) = λ_bg·s_k,
        ## so the between-edge increment is λ_bg·(s_k − s_{k−1}).
        μ_bg_k = max(λ_bg * s_k, zero(Tt))
        Δμ_bg = max(μ_bg_k - μ_bg_prev, zero(Tt))
        μ_bg_prev = μ_bg_k
        raw = p_drc_per_bin[i] * ΔμBVD0[i] + Δμ_bg
        μ_i = isfinite(raw) ? max(raw, eps(typeof(raw))) :
              eps(typeof(raw))
        bin_means[i] = μ_i
        Λ_prev += μ_i
        Λ_at_edges[i] = Λ_prev
    end

    for i in 1:n
        reported_cases[i] ~ safe_nbinomial(k, bin_means[i])
    end

    ## Implied per-suspected positivity at the cut-off (BVD share of the
    ## expected suspected total). Exposed for comparison with the sitrep.
    μ_BVD_cum = sum(p_drc_per_bin[i] * ΔμBVD0[i] for i in 1:n)
    positivity := μ_BVD_cum / Λ_at_edges[end]

    expected_reports_total := Λ_at_edges[end]

    return (; p_drc_per_bin, λ_bg, τ_forward,
        expected_reports_total, positivity,
        report_delay_dist = f_rep,
        μ_BVD0_at_edges, Λ_at_edges)
end

"""
Laboratory pipeline likelihood. Models the lab-confirmed cases over time
and, where available, the testing volume that gates them.

`confirmed_cases` is a per-vintage vector of confirmed-case *increments*
(between-vintage differences of `Cumul positifs`) aligned with `t_edges`.
Each window's new positives are observed on the samples newly analysed in
that window, the same increment construction as the reported and deaths
streams (equations (23)-(24)):

```math
\\Delta C_v \\sim \\mathrm{Binomial}(\\Delta A_v,\\ p_{pos,v}),
\\qquad
p_{pos,v} = s_{test}\\, q_v + (1 - \\text{spec})\\,(1 - q_v),
```

where `ΔA_v` is the analysed count newly added in window `v`
(`samples_analysed`, a known denominator), `q_v` the BVD share of that
slice, `s_test` the PCR sensitivity and `1 − spec` the false-positive rate
(see [`test_specificity_model`](@ref)). Conditioning on the slice rather
than the cumulative total avoids re-using the early analysed samples in
every later denominator.

`q_v` is the tested BVD share set by the **composition link**: the lab
over-tests BVD early (severe cases are triaged first and are more likely
BVD), the severity enrichment relaxing from the early surge `δ0` toward a
persistent selection floor `δ∞` over the suspect-pool composition
`φ_v = μ_BVD,v / (μ_BVD,v + μ_bg,v)` as testing widens (see
[`severity_enrichment_model`](@ref)):

```math
\\mathrm{logit}(q_v) = \\mathrm{logit}(\\varphi_v)
    + \\delta_\\infty + (\\delta_0 - \\delta_\\infty)\\, e^{-c_v / \\text{decay}},
\\qquad
c_v = \\frac{\\sum_{j \\le v} \\Delta A_j}{\\text{volume\\_scale}},
```

with `c_v` the cumulative analysed VOLUME (in units of `volume_scale`
samples), so a lab stall pauses the relaxation. `δ0` is the early severity
log-odds enrichment, `decay` its timescale, and `δ∞` the persistent
selection floor. With `δ∞ = 0` the enrichment relaxes to zero, tying `q`
to `φ` so the confirmed/positivity data identify the non-BVD background
`λ_bg`; a positive `δ∞` keeps the tested share above `φ` (persistent
testing selection), which the weakly-informative floor prior lets the data
decide.

`q_random_effect` adds a per-vintage partially-pooled logit-scale offset
to this baseline share (`q_v = logistic(logit(q_base,v) + σ_q·z_v)`, see
[`confirmed_q_re_model`](@ref)) so each window's positivity can fit the
non-monotone wobble while `s` stays fixed. It is on by default because the
observed positivity is non-monotone; pass `nothing` to recover the smooth
composition baseline.

The cumulative suspect backlog at each edge is split BVD / background:
`μ_BVD(s_v) = p_{DRC,v} · onset_fraction · ∫₀^{s_v} e^{r u} f_rep(s_v−u) du`
(the onset→report Gamma closed form, the same quantity the reported stream
produces) and `μ_bg(s_v) = λ_bg · s_v` (constant-rate non-BVD). Their sum
`N_susp,v` is the report-level backlog; their ratio `μ_BVD / N_susp` is the
count-implied BVD composition, exposed as a diagnostic
(`q_baseline_count`).

Specimens reach the lab after a transport delay, so the cumulative
*received* backlog is `N_susp` convolved with the receipt-delay kernel
`f_receipt` (see [`lab_receipt_delay_model`](@ref)).

`samples_received` is the per-vintage received-count *increment* vector
(differenced `Cumul échantillons reçus`). When supplied each window is
observed as `ΔR_v ~ NegBinomial(τ_forward · ΔN_recv,v, k)`, pinning the
forwarded fraction `τ_forward`. A vector of `missing` entries drops the
received likelihood (and is generated for predictive checks).

`samples_analysed` is the per-vintage analysed-increment denominator
vector (`ΔA_v`) aligned with `t_edges`. It both denominates the confirmed
Binomial and conditions the analysis-capacity random walk (see
[`lab_capacity_model`](@ref)): for windows 2..n the lab processes
`μ_A = backlog·(1 − exp(−κ_v·Δt_v/backlog))` of the available received
backlog, observed as `ΔA_v ~ Poisson(μ_A)`. Window 1's cumulative is the
initial condition. Pass `missing` entries (with a non-missing
`tests_analysed` fallback) to generate predictive draws.

`tests_analysed` is the cut-off cumulative analysed count, retained for
the predictive denominator and the lab-cut-off diagnostic. The per-test
positivity `p_pos`, tested BVD share `q_cutoff` and baseline `q_baseline`
are exposed as derived quantities.

The confirmed stream is a coherent lab-throughput queue that fits ALL
vintages, including the dark windows that lack a published analysed count.
The capacity-limited drain
`μ_A_v = backlog·(1 − exp(−κ_v·Δt_v / backlog))` is computed for every
window (the cumulative analysed advances by the observed increment where
present, else by `μ_A_v`). Observed-denominator windows condition on the
real count: `ΔC_v ~ Binomial(ΔA_obs_v, p_pos_v)` and
`ΔA_obs_v ~ Poisson(μ_A_v)`. Dark windows use the exact marginal of
`Binomial(Poisson(μ_A_v), p_pos_v)`, namely `ΔC_v ~ Poisson(μ_A_v·p_pos_v)`,
so the unobserved denominator is integrated out (Poisson thinning) rather
than carried as a free per-vintage latent. This removes the funnel. The
predicted dark-window denominators are exposed as `μ_A_pred` /
`dark_analysed_total`.

`epi_exclusion` (default `nothing`) supplies the queue's forwarded fraction
`(1 − e)`: with `nothing` the headline fit pins `e = 0` (forward 1); pass
[`epi_exclusion_model`](@ref) for the opt-in `e ~ Beta(2, 12)` sensitivity.
In the queue path `τ_forward` is unused (the received asymptote is
`(1 − e)·N_susp`).

A single observation (`length 1`, `t_edges = [T]`) reduces to the
cumulative confirmed Binomial.
"""
@model function confirmed_cases_model(
        confirmed_cases::AbstractVector,
        samples_analysed::AbstractVector,
        samples_received::AbstractVector,
        tests_analysed::Union{Missing, Integer},
        growth_state, k::Real,
        p_drc_per_bin::AbstractVector, λ_bg::Real, τ_forward::Real, f_rep,
        t_edges::AbstractVector, tests_edge::Real;
        test_sensitivity = test_sensitivity_model(),
        test_specificity = test_specificity_model(),
        severity_enrichment = severity_enrichment_model(),
        receipt_delay = lab_receipt_delay_model(),
        capacity_model = lab_capacity_model,
        capacity_centre::Real = 150.0,
        report_onset_offset::Union{Nothing, Real} = nothing,
        q_random_effect = confirmed_q_re_model,
        epi_exclusion = nothing,
        volume_scale::Real = 200.0,
        onset_fraction::Real = 1.0)
    sensitivity_state ~ to_submodel(test_sensitivity, false)
    specificity_state ~ to_submodel(test_specificity, false)
    ## Composition-linked positivity. The tested BVD share is tied to the
    ## suspect-pool composition `φ = μ_BVD/(μ_BVD+μ_bg)` upsampled by a
    ## decaying severity enrichment `δ0` (see [`severity_enrichment_model`]
    ## (@ref)), so the positivity data identify the background `λ_bg` rather
    ## than a free curve.
    enrich_state ~ to_submodel(severity_enrichment, false)
    δ0 = enrich_state.δ0
    decay_scale = enrich_state.decay_scale
    δ∞ = enrich_state.δ∞
    receipt_state ~ to_submodel(receipt_delay, false)
    ## Epi-exclusion fraction `e` for the throughput queue: the share of
    ## suspects ruled out by epi follow-up and never sampled, so the
    ## received backlog asymptotes to `(1 − e)·N_susp`. Default OFF (e = 0,
    ## forward = 1) for the headline fit; the opt-in Beta(2, 12) prior is a
    ## sensitivity arm.
    if epi_exclusion !== nothing
        epi_state ~ to_submodel(epi_exclusion, false)
        forward_frac = epi_state.forward
    else
        forward_frac = nothing
    end
    n_edges_re = length(t_edges)
    if q_random_effect !== nothing
        q_re_state ~ to_submodel(q_random_effect(n_edges_re), false)
        σ_q = q_re_state.σ_q
        z_q = q_re_state.z_q
    else
        σ_q = nothing
        z_q = nothing
    end
    n_edges = length(t_edges)
    capacity_state ~ to_submodel(
        capacity_model(n_edges; capacity_centre = capacity_centre), false)
    s_test = sensitivity_state.s_test
    spec_test = specificity_state.spec_test
    ## `δ0`, `decay_scale` are set above by the severity-enrichment submodel.
    f_receipt = receipt_state.dist
    κ = capacity_state.capacity
    r = growth_state.r
    T = growth_state.T

    n = length(t_edges)
    n == length(p_drc_per_bin) ||
        error("p_drc_per_bin length must match t_edges (got " *
              "$(length(p_drc_per_bin)) vs $n)")
    n == length(confirmed_cases) ||
        error("confirmed_cases length must match t_edges (got " *
              "$(length(confirmed_cases)) vs $n)")
    n == length(samples_analysed) ||
        error("samples_analysed length must match t_edges (got " *
              "$(length(samples_analysed)) vs $n)")
    n == length(samples_received) ||
        error("samples_received length must match t_edges (got " *
              "$(length(samples_received)) vs $n)")

    ## Cumulative suspect backlog at each edge, split BVD / background. The
    ## BVD-suspected cumulative is `μ_BVD(s) = p_DRC · onset_fraction ·
    ## ∫₀^s e^{r u} f_rep(s−u) du` (the onset→report Gamma closed form, the
    ## same quantity the reported stream produces); the non-BVD background
    ## cumulative is `μ_bg(s) = λ_bg · s` (constant-rate). Their sum
    ## `N_susp,v` is the backlog the received-count likelihood conditions on
    ## (below). Their ratio `μ_BVD / N_susp` is the count-implied BVD
    ## composition `φ`, the plateau the composition link relaxes to.
    t_report = report_onset_offset === nothing ? zero(T) :
               max(T - report_onset_offset, zero(T))

    Tt = typeof(float(r) * float(λ_bg) * onset_fraction *
                float(δ0) * float(decay_scale))
    ## Per-edge analysed denominator: the observed analysed increment where
    ## present, else zero (dark windows are owned by the Poisson-thinned
    ## marginal below). Both the volume clock and the confirmed Binomial use
    ## this.
    A_imp = Vector{Tt}(undef, n)
    for i in 1:n
        A_imp[i] = samples_analysed[i] === missing ? zero(Tt) :
                   convert(Tt, samples_analysed[i])
    end
    p_pos = Vector{Tt}(undef, n)
    q_at = Vector{Tt}(undef, n)
    qinf_count_at = Vector{Tt}(undef, n)
    Nsusp_at = Vector{Tt}(undef, n)
    Λ_at_edges = Vector{Tt}(undef, n)
    ## Composition-linked positivity with modelled specificity. The tested
    ## BVD share is the suspect-pool composition `φ = μ_BVD/(μ_BVD+μ_bg)`
    ## upsampled by a severity log-odds enrichment `δ0·exp(−c/decay)` that
    ## decays as testing widens, so `q → φ` at the plateau. The per-test
    ## positivity combines true and false positives,
    ##   p_pos = s·q + (1 − spec)·(1 − q),
    ## with `1 − spec` the false-positive rate. The binomial conditions on
    ## the observed analysed denominator.
    ## Cumulative analysed volume at each edge, for the volume-indexed
    ## clock. The severe-first share then relaxes with how many samples have
    ## been processed rather than calendar time, so a lab stall (no new
    ## analysed) does not advance the curve.
    analysed_cum_at = Vector{Tt}(undef, n)
    acc_a = zero(Tt)
    for i in 1:n
        acc_a += A_imp[i]
        analysed_cum_at[i] = acc_a
    end
    for i in 1:n
        s_i = oftype(T, t_edges[i])
        μ_bvd = s_i <= zero(s_i) ? zero(Tt) :
                convert(Tt, p_drc_per_bin[i]) * onset_fraction *
                delay_convolution(one(r), r, s_i, f_rep)
        μ_bg = max(convert(Tt, λ_bg) * s_i, zero(Tt))
        denom = μ_bvd + μ_bg
        Nsusp_at[i] = denom
        qinf_count_at[i] = denom > zero(Tt) ?
                           clamp(μ_bvd / denom, zero(Tt), one(Tt)) :
                           zero(Tt)
        ## Volume clock for the severity enrichment: cumulative analysed
        ## volume in units of `volume_scale` samples, so a lab stall (no new
        ## analysed) pauses the relaxation.
        c_i = analysed_cum_at[i] / convert(Tt, volume_scale)
        ## Composition link with persistent testing selection: the suspect-
        ## pool composition φ = μ_BVD/(μ_BVD+μ_bg) upsampled by a severity
        ## log-odds enrichment δ∞ + (δ0−δ∞)·exp(−c/decay) that relaxes from
        ## the early severity surge δ0 to a persistent selection floor δ∞.
        ## δ∞ = 0 recovers the pure composition link (tested share → φ, so
        ## positivity pins λ_bg); δ∞ > 0 keeps the tested share above φ, so
        ## positivity no longer pins the background.
        φ = clamp(qinf_count_at[i], eps(Tt), one(Tt) - eps(Tt))
        δ∞_c = convert(Tt, δ∞)
        δ_i = δ∞_c + (convert(Tt, δ0) - δ∞_c) *
                     exp(-c_i / convert(Tt, decay_scale))
        q_base = logistic(logit(φ) + δ_i)
        ## Optional partially-pooled per-window offset on the tested BVD
        ## share, logit-scale, so each vintage's positivity can fit the
        ## non-monotone wobble while `s` stays fixed. `σ_q → 0` recovers the
        ## smooth baseline curve.
        if z_q === nothing
            q = q_base
        else
            qb = clamp(q_base, eps(Tt), one(Tt) - eps(Tt))
            q = clamp(
                logistic(logit(qb) + convert(Tt, σ_q) * convert(Tt, z_q[i])),
                zero(Tt), one(Tt))
        end
        q_at[i] = q
        p_raw = s_test * q + (one(Tt) - spec_test) * (one(Tt) - q)
        p_pos[i] = isfinite(p_raw) ?
                   clamp(p_raw, eps(Tt), one(Tt) - eps(Tt)) : eps(Tt)
        ## Expected confirmed increment: newly-analysed slice ΔA_v times its
        ## positivity (a diagnostic). Uses the imputed denominator where the
        ## analysed count is missing.
        Λ_at_edges[i] = p_pos[i] * A_imp[i]
    end

    ## Received queue: the suspect backlog reaches the lab after a transport
    ## delay, so the cumulative received is the suspect backlog convolved with
    ## the receipt-delay kernel f_receipt. Confirmed uses a pooled (constant)
    ## p_drc, so a single representative value drives the continuous backlog.
    ## Built before the confirmed likelihood because the queue path's
    ## per-window analysed mean μ_A depends on this received backlog.
    p_drc_c = convert(Tt, p_drc_per_bin[1])
    N_susp_fn = let r = r, f_rep = f_rep, p_drc_c = p_drc_c,
        onset_fraction = onset_fraction, λ_bg = λ_bg

        u -> u <= zero(u) ? zero(Tt) :
             p_drc_c * onset_fraction *
             delay_convolution(one(r), r, u, f_rep) +
             max(convert(Tt, λ_bg) * u, zero(Tt))
    end
    N_recv_at = Vector{Tt}(undef, n)
    for i in 1:n
        s_i = oftype(T, t_edges[i])
        raw = s_i <= zero(s_i) ? zero(Tt) :
              convert(Tt, delay_convolution(N_susp_fn, s_i, f_receipt))
        N_recv_at[i] = max(raw, zero(Tt))
    end

    ## Forwarded fraction of the received backlog `(1 − e)` from the
    ## epi-exclusion submodel (default `e = 0`, forward = 1).
    fwd = forward_frac === nothing ? one(Tt) : convert(Tt, forward_frac)

    ## Per-window predicted analysed mean μ_A from the capacity-limited drain
    ## of the received backlog. Computed for ALL windows (the queue path needs
    ## μ_A for the dark windows that lack an observed denominator). Window 1's
    ## cumulative analysed (observed or predicted) is the initial condition.
    μ_A_at = Vector{Tt}(undef, n)
    let analysed_cum = zero(Tt)
        for i in 1:n
            ## Window 1 has no preceding edge. Use the spacing to the next
            ## vintage as the representative window length, NOT the full
            ## seeding-to-edge-1 elapsed time: the latter (≈ t_edges[1],
            ## ~120 days) makes cap ≫ backlog and drains the entire
            ## pre-window received backlog in one step, blowing up the
            ## first dark window's μ_A. A single-vintage fit falls back to
            ## the cumulative.
            Δt = if i == 1
                n >= 2 ?
                max(oftype(T, t_edges[2]) - oftype(T, t_edges[1]),
                    zero(T)) :
                max(oftype(T, t_edges[1]), zero(T))
            else
                max(oftype(T, t_edges[i]) - oftype(T, t_edges[i - 1]),
                    zero(T))
            end
            recv_cum = fwd * N_recv_at[i]
            backlog = max(recv_cum - analysed_cum, eps(Tt))
            cap = max(convert(Tt, κ[i]) * Δt, zero(Tt))
            μ_A_raw = backlog * (one(Tt) - exp(-cap / backlog))
            μ_A_at[i] = isfinite(μ_A_raw) ? max(μ_A_raw, eps(Tt)) : eps(Tt)
            ## Advance the analysed cumulative by the observed increment where
            ## present, else the predicted mean (carries the queue through the
            ## dark windows).
            analysed_cum += samples_analysed[i] === missing ? μ_A_at[i] :
                            convert(Tt, samples_analysed[i])
        end
    end

    ## Coherent lab-throughput queue with a Poisson-thinned denominator.
    ## OBSERVED-denominator windows condition on the real analysed
    ## increment: confirmed_v ~ Binomial(ΔA_obs, p_pos) AND
    ## ΔA_obs ~ Poisson(μ_A) (pins positivity and capacity). DARK windows
    ## (no published analysed) use the exact marginal of
    ## Binomial(Poisson(μ_A), p_pos), namely
    ## confirmed_v ~ Poisson(μ_A · p_pos): the denominator is integrated
    ## out, NOT a free latent, so there is no per-vintage funnel.
    for i in 1:n
        if samples_analysed[i] === missing
            μ_c = isfinite(μ_A_at[i] * p_pos[i]) ?
                  max(μ_A_at[i] * p_pos[i], eps(Tt)) : eps(Tt)
            confirmed_cases[i] ~ Poisson(μ_c)
        else
            samples_analysed[i] ~ Poisson(μ_A_at[i])
            A_i = Int(samples_analysed[i])
            confirmed_cases[i] ~ Binomial(A_i, p_pos[i])
        end
    end

    ## Received increments observed where present, mean (1 − e)·ΔN_recv.
    ΔN_recv = daily_increment_kernel(N_recv_at)
    recv_means = Vector{Tt}(undef, n)
    for i in 1:n
        raw = fwd * ΔN_recv[i]
        recv_means[i] = isfinite(raw) ? max(raw, eps(Tt)) : eps(Tt)
        samples_received[i] === missing && continue
        samples_received[i] ~ safe_nbinomial(k, recv_means[i])
    end

    ## Per-test positivity at the lab cut-off (last edge ≤ tests_edge):
    ## exposed for comparison with the sitrep `Taux de positivité`.
    te_idx = something(
        findlast(
            i -> t_edges[i] <= tests_edge + sqrt(eps(typeof(tests_edge))),
            1:n), n)
    p_positive := p_pos[te_idx]
    q_cutoff := q_at[te_idx]
    ## Plateau tested share: the severity enrichment relaxes to the pool
    ## composition, so the plateau is the count-implied composition itself.
    q_baseline := qinf_count_at[te_idx]
    q_baseline_count := qinf_count_at[te_idx]
    δ0_out := δ0
    τ_forward_out := τ_forward
    ## Daily analysis capacity at the cut-off vintage (samples/day).
    capacity_cutoff := κ[te_idx]

    ## Cumulative expected totals at the cut-off: sum the per-window
    ## increment means (confirmed positives, received samples).
    expected_confirmed_total := sum(Λ_at_edges)
    expected_received_total := sum(recv_means)

    ## Queue-path diagnostic: the predicted analysed denominator μ_A summed
    ## over the dark windows (those without a published analysed count). Lets
    ## the dark-window denominators be compared against the observed 23-28 May
    ## anchors.
    dark_analysed_total := sum(samples_analysed[i] === missing ? μ_A_at[i] :
                               zero(eltype(μ_A_at)) for i in 1:n)
    ## Per-window predicted analysed mean and received mean, exposed for the
    ## dark-window-versus-anchor comparison and the received posterior check.
    μ_A_pred := copy(μ_A_at)
    recv_pred := copy(recv_means)

    return (; p_positive, p_pos, q_at, qinf_count_at, Nsusp_at,
        N_recv_at, recv_means, μ_A_at, dark_analysed_total, capacity = κ,
        s_test, spec_test, δ0, decay_scale, τ_forward, p_drc_per_bin,
        expected_confirmed_total, expected_received_total, Λ_at_edges)
end

"""
Count-implied BVD composition among suspects at each elapsed-time edge:
`q(s) = μ_BVD(s) / N_susp(s)`, with `μ_BVD = p_drc · onset_fraction ·
∫₀^s e^{r u} f_rep(s−u) du` the cumulative BVD-suspect backlog and
`N_susp = μ_BVD + λ_bg·s` adding the constant-rate non-BVD background. The
same `q_baseline_count` quantity [`confirmed_cases_model`](@ref) exposes;
factored out so the confirmed-death stream can evaluate it at the
suspected-death edges from the shared report delay, ascertainment,
background and growth state. Returns a vector aligned with `t_edges`.
"""
function bvd_count_composition(r, p_drc, λ_bg, f_rep,
        t_edges::AbstractVector; onset_fraction::Real = 1.0)
    Tt = typeof(float(r) * float(λ_bg) * float(p_drc) * onset_fraction)
    q = Vector{Tt}(undef, length(t_edges))
    for i in eachindex(t_edges)
        s_i = convert(Tt, t_edges[i])
        μ_bvd = s_i <= zero(s_i) ? zero(Tt) :
                convert(Tt, p_drc) * onset_fraction *
                delay_convolution(one(r), r, s_i, f_rep)
        μ_bg = max(convert(Tt, λ_bg) * s_i, zero(Tt))
        denom = μ_bvd + μ_bg
        ## Guard the ratio: at extreme growth `μ_bvd` can overflow to Inf,
        ## giving Inf/Inf = NaN. Fall back to zero share when non-finite.
        ratio = denom > zero(Tt) ? μ_bvd / denom : zero(Tt)
        q[i] = isfinite(ratio) ? clamp(ratio, zero(Tt), one(Tt)) : zero(Tt)
    end
    return q
end

"""
Modelled cumulative BVD-death trajectory at each elapsed-time edge:
`μ_death(s) = CFR · p_deaths · onset_fraction · ∫₀^s e^{r u} f_death(s−u) du`,
the same CFR-weighted convolution [`deaths_model`](@ref) integrates for the
suspected-death expectation. Factored out so the confirmed-death stream can
evaluate the modelled BVD-death expectation at the confirmed-death edges
from the shared growth, CFR and onset-to-death delay. Returns a vector
aligned with `t_edges`.
"""
function bvd_death_trajectory(r, CFR, f_death, t_edges::AbstractVector;
        p_deaths::Real = 1.0, onset_fraction::Real = 1.0)
    return [s <= zero(s) ? zero(float(r) * float(CFR)) :
            p_deaths * onset_fraction *
            delay_convolution(CFR, r, s, f_death)
            for s in t_edges]
end

"""
Laboratory-confirmed-deaths likelihood. The DRC sitrep front page reports
`Cumul décès parmi les confirmés`, the cumulative deaths that have been
laboratory-confirmed (17 at the 28 May cut-off). This is a genuine
lab/positivity process on the post-mortem death specimens that reach the
laboratory, mirroring the confirmed-case pipeline
([`confirmed_cases_model`](@ref)) rather than the earlier
`coverage_death · s` thinning of the modelled BVD-death trajectory (issue
#193): that thinning conditioned directly on the latent BVD-death
trajectory with no observation/background process, so all the slack sat in
`coverage_death` at its boundary.

A fraction `τ_death` of the suspect-death backlog is forwarded to and
analysed by the laboratory; the BVD share of that pool sets the positivity:

```math
\\Delta D_{conf,v} \\sim
    \\mathrm{NegBinomial}(
        \\tau_{death}\\cdot p_{pos,death,v}\\cdot\\Delta N_{death,v},\\ k),
\\qquad
p_{pos,death,v} = s\\, q_{death,v} + (1 - \\text{spec})(1 - q_{death,v}),
```

with `ΔN_death,v` the between-edge increment of `nsusp_death_at_edges` (the
cumulative suspect-death backlog, BVD plus non-BVD background, that
[`deaths_model`](@ref) builds), `q_death,v = μ_BVD_death / N_death_susp` the
BVD share of the suspect-death pool at edge `v` (from `bvd_death_at_edges`
and `nsusp_death_at_edges`), and `k` the shared count dispersion. The PCR
sensitivity `s` and specificity `spec` are *imported shared* from the
confirmed-case lab submodel (assay properties, not re-sampled), so the
death stream inherits the case lab calibration. `τ_death` is the
death-specimen forwarding fraction from [`death_forward_model`](@ref); the
BVD-share signal lives in the composition-driven positivity, not in
`τ_death`, so it is identified by the 17 confirmed deaths without sitting
at a boundary. `confirmed_deaths` accepts `missing` entries for
posterior-predictive generation.
"""
@model function confirmed_deaths_model(
        confirmed_deaths::AbstractVector,
        bvd_death_at_edges::AbstractVector,
        nsusp_death_at_edges::AbstractVector,
        s::Real, spec::Real, k::Real;
        death_forward = death_forward_model())
    n = length(nsusp_death_at_edges)
    n == length(confirmed_deaths) ||
        error("confirmed_deaths length must match nsusp_death_at_edges " *
              "(got $(length(confirmed_deaths)) vs $n)")
    n == length(bvd_death_at_edges) ||
        error("bvd_death_at_edges length must match nsusp_death_at_edges " *
              "(got $(length(bvd_death_at_edges)) vs $n)")

    forward_state ~ to_submodel(death_forward, false)
    τ_death = forward_state.τ_death

    ## Per-edge BVD share of the suspect-death pool (composition), then the
    ## death-specimen positivity combining true positives at sensitivity `s`
    ## and false positives at false-positive rate `1 − spec`, the same form
    ## the confirmed-case stream uses. The expected confirmed increment is
    ## the forwarded fraction `τ_death` of the suspect-death backlog
    ## increment `ΔN_death`, weighted by that positivity. Clamped positive
    ## inside `safe_nbinomial`.
    Tt = typeof(float(s) * float(spec) * float(τ_death) *
                float(nsusp_death_at_edges[1]))
    ΔN = daily_increment_kernel(convert(Vector{Tt}, nsusp_death_at_edges))
    conf_means = Vector{Tt}(undef, n)
    q_death_at = Vector{Tt}(undef, n)
    p_pos_death_at = Vector{Tt}(undef, n)
    for i in 1:n
        denom = convert(Tt, nsusp_death_at_edges[i])
        q_d = denom > zero(Tt) ?
              clamp(convert(Tt, bvd_death_at_edges[i]) / denom,
            zero(Tt), one(Tt)) : zero(Tt)
        q_death_at[i] = q_d
        p_raw = convert(Tt, s) * q_d +
                (one(Tt) - convert(Tt, spec)) * (one(Tt) - q_d)
        p_pos_death_at[i] = isfinite(p_raw) ?
                            clamp(p_raw, eps(Tt), one(Tt) - eps(Tt)) :
                            eps(Tt)
        raw = convert(Tt, τ_death) * p_pos_death_at[i] * ΔN[i]
        conf_means[i] = isfinite(raw) ? max(raw, eps(Tt)) : eps(Tt)
    end

    for i in 1:n
        confirmed_deaths[i] ~ safe_nbinomial(k, conf_means[i])
    end

    τ_death_out := τ_death
    q_death_cutoff := q_death_at[end]
    p_pos_death_cutoff := p_pos_death_at[end]
    expected_confirmed_deaths_total := sum(conf_means)

    return (; τ_death, s, spec, q_death_at, p_pos_death_at, conf_means,
        expected_confirmed_deaths_total)
end

"""
Bin-mean kernel for the per-bin count likelihoods. For `n`
cumulative-intensity values ``\\Lambda(t_k)`` at the bin edges, returns
the `n` between-edge increments `[Λ(t_1) - Λ_0, Λ(t_2)-Λ(t_1), ...]`
clamped to be strictly positive and finite. Shared by
[`reported_cases_model`](@ref), [`confirmed_cases_model`](@ref),
[`deaths_model`](@ref) and [`exports_deaths_model`](@ref) so the
bin-difference logic lives in one place. `init` is the cumulative value
the first bin is measured from (`Λ_0`); it defaults to zero, so the
first element is the cumulative at the first edge and a single edge
reduces to the cumulative single-total mean. The exported-deaths stream
passes its pre-death survival weight as `init`.
"""
function daily_increment_kernel(Λ_at_edges::AbstractVector, init = nothing)
    n = length(Λ_at_edges)
    Tt = eltype(Λ_at_edges)
    means = Vector{Tt}(undef, n)
    ## `init` is positional, not a keyword: a keyword bundles it into a
    ## `(; init)` NamedTuple, and when the caller's `init` is a constant
    ## (e.g. the pre-death stretch collapses to `zero(T)` for extreme NUTS
    ## proposals) that single-field struct is wholly inactive, which
    ## Enzyme's reverse-mode struct rule rejects. A positional argument is
    ## Const-annotated directly and avoids the mixed-activity struct.
    Λ_prev = init === nothing ? zero(Tt) : convert(Tt, init)
    @inbounds for i in 1:n
        raw = Λ_at_edges[i] - Λ_prev
        means[i] = isfinite(raw) ? max(raw, eps(typeof(raw))) :
                   eps(typeof(raw))
        Λ_prev = Λ_at_edges[i]
    end
    return means
end

"""
Time-resolved deaths-among-exports likelihood. Models the dated
Uganda export deaths as an inhomogeneous Poisson process: a
continuous survival term over the pre-death stretch followed by
per-day Poisson counts. The detection window and traveller volume
are supplied by [`exports_model`](@ref) so the two Uganda-side
likelihoods share person-time.
"""
@model function exports_deaths_model(
        export_deaths_daily::AbstractVector,
        growth_state, CFR::Real, delay_dist, p_uganda::Real;
        pre_start_deaths::Union{Missing, Integer} = 0,
        window::Real,
        daily_travellers::Real,
        source_population::Real = ITURI_POPULATION,
        onset_fraction::Real = 1.0)
    cumulative = growth_state.cumulative
    T = growth_state.T
    q = daily_travellers / source_population
    n = length(export_deaths_daily)   # days from earliest death to cut-off

    ## Precompute the onset-to-death CDF once and reuse it across every
    ## bin edge below (`T - s ≤ window` over the domain; see
    ## `ExportDeathDelay`). The latent trajectory is cumulative
    ## *infections*, so `onset_fraction` (the incubation mgf) maps it onto
    ## onsets before the onset-to-death convolution.
    delay = ExportDeathDelay(delay_dist, window)
    Λ(t) = onset_fraction * expected_exports_deaths(
        cumulative, delay, CFR, p_uganda, q, t, window)

    ## Pre-death zero stretch as one Poisson observed at 0; `missing`
    ## generates it for predictive checks (see equation (20)).
    pre = T - n > zero(T) ? Λ(T - n) : zero(T)
    pre_start_deaths ~ Poisson(max(pre, zero(pre)))

    ## Per-day means via the shared between-edge differencing
    ## (`daily_increment_kernel`), the same construction as the DRC
    ## streams but starting from the pre-death survival weight `pre`
    ## rather than zero.
    Λ_at_edges = [Λ(T - n + i) for i in 1:n]
    μ_day = daily_increment_kernel(Λ_at_edges, pre)
    for i in 1:n
        export_deaths_daily[i] ~ Poisson(μ_day[i])
    end

    return (;)
end

"""
Time-resolved deaths-among-exports likelihood with an explicit
infection→detection delay. The delay-convolution counterpart of
[`exports_deaths_model`](@ref): the per-edge intensity is the
infection→detection delay-survival expectation
[`expected_exports_deaths_delay`](@ref), where the rectangular detection
window of the McCabe path is replaced by the infection→detection survival
`f_det`. Because the at-risk clock starts at infection, the detection
survival is 1 at age 0 (a just-infected traveller is certainly not yet
detected). The death timing is likewise keyed to infection: `delay_dist`
is the infection→death delay (incubation ⊕ onset-to-death, see
[`combined_delay`](@ref)), not the bare onset-to-death delay, so deaths
are not timed one incubation period too early. The pre-death survival
stretch and per-day [`daily_increment_kernel`](@ref) structure are
unchanged. `f_det` (the incubation ⊕ onset-to-report delay) and
`daily_travellers` are supplied by [`exports_delay_model`](@ref) so the
two Uganda-side likelihoods share person-time. Incubation enters both
`f_det` and `delay_dist`, a slight accepted double-count of the shared
incubation period (better than omitting it on death entirely).

`last_offset` (default 0) anchors the last series day `last_offset` days
before the cut-off, matching [`exports_daily_delay_model`](@ref): both
travel-gated streams are run only up to the most recent reported import to
Uganda, since movement patterns likely shift and the most recent days are
right-truncated by reporting lag.
"""
@model function exports_deaths_delay_model(
        export_deaths_daily::AbstractVector,
        growth_state, CFR::Real, delay_dist, p_uganda::Real;
        pre_start_deaths::Union{Missing, Integer} = 0,
        last_offset::Real = 0,
        f_det,
        daily_travellers::Real,
        source_population::Real = ITURI_POPULATION)
    r = growth_state.r
    T = growth_state.T
    q = daily_travellers / source_population
    n = length(export_deaths_daily)   # earliest death to last series day

    ## The series ends `last_offset` days before the cut-off (the date of
    ## the most recent reported import to Uganda). Like the export-case
    ## stream, this travel-gated death stream is only run up to that day:
    ## movement patterns likely shift and the days after the last import
    ## are right-truncated by reporting lag. The at-risk clock anchors on
    ## this date.
    t_last = T - last_offset

    ## One at-risk trajectory over [0, t_last] shared across every edge,
    ## the same engine as the export-case stream. The per-edge intensity
    ## is the infection-keyed delay-survival weight: the infection→death
    ## CDF weighted by the infection→detection survival (1 at age 0).
    traj = ExportRiskTrajectory(t_last, r)
    edges = [t_last - n + i for i in 0:n]
    Λ = (CFR * p_uganda * q) .*
        export_death_at_risk(traj, edges, f_det, delay_dist)

    ## Pre-death zero stretch as one Poisson observed at 0; `missing`
    ## generates it for predictive checks.
    pre = max(Λ[1], zero(eltype(Λ)))
    pre_start_deaths ~ Poisson(pre)

    ## Per-day means via the shared between-edge differencing, starting
    ## from the pre-death survival weight `pre`.
    μ_day = daily_increment_kernel(Λ[2:end], pre)
    for i in 1:n
        export_deaths_daily[i] ~ Poisson(μ_day[i])
    end

    return (;)
end

"""
First-export-detection timing survival term. Adds a one-sided
`Pr(no export detected before t1)` Poisson observation at zero,
matching the at-risk export person-time intensity. Passing
`delta = missing` makes the submodel a no-op.
"""
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
                            expected_exports(cumulative, p_uganda, q, t1,
            window; r = growth_state.r)
        ## No detection before t1 as a Poisson observed at 0; `missing`
        ## generates it for predictive checks (see equation (22)).
        pre_detection_exports ~ Poisson(max(survived_exports, zero(T)))
    end

    return (;)
end

"""
First-export-detection timing survival term with an explicit
infection→detection delay. The delay-convolution counterpart of
[`exports_detection_timing_model`](@ref): the at-risk export
person-time uses [`expected_exports_delay`](@ref) with the
infection→detection survival `f_det` instead of the rectangular window.
Adds the same one-sided `Pr(no export detected before t1)` Poisson
observation at zero. `f_det` (the incubation ⊕ onset-to-report delay, see
[`combined_delay`](@ref)) and `daily_travellers` are supplied by
[`exports_delay_model`](@ref). Passing `delta = missing` makes the
submodel a no-op.
"""
@model function exports_detection_timing_delay_model(
        growth_state, p_uganda::Real;
        delta::Union{Missing, Real},
        pre_detection_exports::Union{Missing, Integer} = 0,
        f_det,
        daily_travellers::Real,
        source_population::Real = ITURI_POPULATION)
    if !ismissing(delta)
        r = growth_state.r
        T = growth_state.T
        t1 = T - delta
        q = daily_travellers / source_population
        survived_exports := t1 <= zero(T) ? zero(T) :
                            expected_exports_delay(r, p_uganda, q, t1,
            f_det)
        ## No detection before t1 as a Poisson observed at 0; `missing`
        ## generates it for predictive checks.
        pre_detection_exports ~ Poisson(max(survived_exports, zero(T)))
    end

    return (;)
end
