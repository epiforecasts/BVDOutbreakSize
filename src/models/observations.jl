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
Exports likelihood with an explicit onset-to-detection delay. The
delay-convolution counterpart of [`exports_model`](@ref): the at-risk
export person-time is the onset-to-detection delay-survival integral
[`expected_exports_delay`](@ref) rather than the McCabe et al.
rectangular detection window. The onset-to-detection delay `f_det` is
passed in rather than sampled here: the export count is a single datum
and cannot identify its own delay, and "symptom onset → entering
surveillance as a suspected case" is taken to be the same process abroad
as the DRC onset-to-report delay, so the composers pass the `f_rep` of
[`report_delay_model`](@ref) (the same draw the reported-cases stream
uses, see [`reported_cases_model`](@ref)) and the reported stream pins
it. Samples the traveller volume submodel internally and exposes `f_det`
and `daily_travellers` for the export-deaths and detection-timing delay
submodels. `onset_fraction` is the incubation mgf, mapping the latent
infection trajectory onto onsets before the onset-to-detection delay
acts. Couples to `C(T)` through a Poisson likelihood.
"""
@model function exports_delay_model(
        exported_cases::Union{Missing, Integer},
        growth_state, p_uganda::Real, f_det;
        source_population::Real = ITURI_POPULATION,
        traveller = traveller_volume_model(),
        onset_fraction::Real = 1.0)
    r = growth_state.r
    T = growth_state.T

    travel_state ~ to_submodel(traveller, false)
    daily_travellers = travel_state.daily_travellers

    q = daily_travellers / source_population
    expected_exports_T := expected_exports_delay(r, p_uganda, q, T, f_det;
        onset_fraction = onset_fraction)

    exported_cases ~ Poisson(expected_exports_T)

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
"""
@model function deaths_model(
        total_deaths::AbstractVector,
        growth_state, k::Real, t_edges::AbstractVector;
        delay = delay_model(),
        cfr = cfr_model(),
        p_deaths::Real = 1.0,
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

    ## CFR-weighted cumulative expected deaths at each bin edge. The
    ## latent trajectory is cumulative *infections*, so `onset_fraction`
    ## (the incubation mgf) maps it onto onsets before the onset-to-death
    ## convolution; the drift factor `p_deaths` scales the whole
    ## trajectory. The between-edge increment is the per-bin NegBinomial
    ## mean (with a NaN / Inf-safe positive clamp via
    ## `daily_increment_kernel`).
    Λ_at_edges = [s <= zero(s) ? zero(s) :
                  p_deaths * onset_fraction *
                  delay_convolution(CFR, r, s, f_death)
                  for s in t_edges]
    bin_means = daily_increment_kernel(Λ_at_edges)

    for i in 1:n
        total_deaths[i] ~ safe_nbinomial(k, bin_means[i])
    end

    expected_deaths_T := Λ_at_edges[end]

    return (; CFR, p_deaths, delay_dist = f_death,
        expected_deaths_T, Λ_at_edges)
end

"""
Reported (suspected) cases likelihood. Couples the observed DRC
suspected-case counts to `C(T)` as the sum of (i) a BVD-driven
contribution `p_drc · onset_fraction · ∫₀^s exp(r·u) · f_rep(s-u) du` using
the onset-to-report delay [`report_delay_model`](@ref) and (ii) a non-BVD
background `λ_bg · s`, with `λ_bg` and the testing fraction `τ_test`
supplied by [`test_positivity_model`](@ref). The latent trajectory is
cumulative *infections*, so `onset_fraction` (the incubation mgf) maps it
onto onsets before the onset-to-report convolution; the background term
is non-BVD and unscaled.

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
        onset_fraction::Real = 1.0)
    report_state ~ to_submodel(report_delay, false)
    test_positivity_state ~ to_submodel(test_positivity, false)
    λ_bg = test_positivity_state.λ_bg
    τ_test = test_positivity_state.τ_test
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
    for i in 1:n
        s_k = t_edges[i]
        Δt = i == 1 ? s_k : (s_k - t_edges[i - 1])
        raw = p_drc_per_bin[i] * ΔμBVD0[i] + λ_bg * max(Δt, zero(Δt))
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

    return (; p_drc_per_bin, λ_bg, τ_test,
        expected_reports_total, positivity,
        report_delay_dist = f_rep,
        μ_BVD0_at_edges, Λ_at_edges)
end

"""
Laboratory pipeline likelihood. Models the lab-confirmed cases over time
and, where available, the testing volume that gates them.

`confirmed_cases` (`Cumul positifs`) is a per-vintage vector aligned with
`t_edges`; the model fits the between-vintage cumulative-count increments
as independent NegBinomial terms with per-bin mean

```math
\\mu_v^{conf} = p_{DRC,v} \\cdot s_{test} \\cdot \\tau_{test} \\cdot
               \\bigl(I_{lab,0}(s_v) - I_{lab,0}(s_{v-1})\\bigr),
\\qquad
I_{lab,0}(s) = \\int_0^{s} \\mu_{BVD,0}(u)\\, f_{lab}(s - u)\\,du,
```

where `μ_BVD,0` is the unit-ascertainment onset-to-report BVD cumulative
(kernel `f_rep` passed in from [`reported_cases_model`](@ref)) and
`s_test` is the PCR sensitivity. A [`DailyBVDTrajectory`](@ref)
precomputes `μ_BVD,0` once on a shared Gauss-Legendre node set so the
outer quadrature of the lab convolution is reused across every edge.

`tests_analysed` (`Cumul échantillons analysés`) is a single cumulative
count observed at its own elapsed time `tests_edge` (which may lag the
case cut-off `T` if lab reporting stops earlier). When present it enters
as one NegBinomial term with mean
``\\tau\\,(\\text{BVD}_\\text{tested} + \\text{bg}_\\text{tested})``
accumulated to `tests_edge`; the background tested volume uses the
lab-delay CDF ``F_\\text{lab}`` (via `_gamma_cdf`, finite α gradient
under Mooncake, see #138) so right-truncation is handled by the
integration limits. The per-test positivity
`s · BVD_tested / (BVD_tested + bg_tested)` is exposed as a derived
quantity for comparison with the sitrep figure; confirmed counts are not
re-observed conditional on tests, so the two streams are not
double-counted. Pass `tests_analysed = missing` to drop the testing
stream.

A single confirmed observation (`length 1`, `t_edges = [T]`) reduces to
the cumulative confirmed likelihood.
"""
@model function confirmed_cases_model(
        confirmed_cases::AbstractVector,
        tests_analysed::Union{Missing, Integer},
        growth_state, k::Real,
        p_drc_per_bin::AbstractVector, λ_bg::Real, τ_test::Real, f_rep,
        t_edges::AbstractVector, tests_edge::Real;
        lab_delay = lab_delay_model(),
        test_sensitivity = test_sensitivity_model(),
        onset_fraction::Real = 1.0)
    lab_state ~ to_submodel(lab_delay, false)
    sensitivity_state ~ to_submodel(test_sensitivity, false)
    f_lab = lab_state.dist
    s_test = sensitivity_state.s_test
    r = growth_state.r
    T = growth_state.T
    α_lab = lab_state.α
    θ_lab = lab_state.θ

    n = length(t_edges)
    n == length(p_drc_per_bin) ||
        error("p_drc_per_bin length must match t_edges (got " *
              "$(length(p_drc_per_bin)) vs $n)")
    n == length(confirmed_cases) ||
        error("confirmed_cases length must match t_edges (got " *
              "$(length(confirmed_cases)) vs $n)")

    ## Unit-ascertainment μ_BVD,0 at the shared Gauss-Legendre nodes over
    ## [0, T]; the outer quadrature of I_lab,0 is reused across every bin
    ## edge. Per-bin ascertainment is then applied on the between-edge
    ## I_lab,0 increment so each bin's mean tracks its own draw.
    ## The latent trajectory is cumulative *infections*, so `onset_fraction`
    ## (the incubation mgf) maps it onto onsets before the onset-to-report
    ## and report-to-lab convolutions; it scales the whole BVD lab
    ## cumulative, hence both the confirmed counts and the BVD tested
    ## volume below.
    trajectory = DailyBVDTrajectory(T, r, f_rep)
    I_lab0_edges = onset_fraction .*
                   delay_convolution(trajectory, t_edges, f_lab)
    ΔIlab0 = daily_increment_kernel(I_lab0_edges)

    Tt = eltype(I_lab0_edges)
    conf_means = Vector{Tt}(undef, n)
    Λ_at_edges = Vector{Tt}(undef, n)
    ## Cumulative unit-ascertainment-weighted BVD tested volume (before
    ## the τ gate) accumulated up to `tests_edge`, summing the same
    ## per-bin ascertainment increments the confirmed likelihood uses.
    bvd_tested_unit = zero(Tt)
    Λ_prev = zero(Tt)
    for i in 1:n
        raw_bvd = p_drc_per_bin[i] * ΔIlab0[i]
        bvd_inc = isfinite(raw_bvd) ? max(raw_bvd, eps(typeof(raw_bvd))) :
                  eps(typeof(raw_bvd))
        μ_i = s_test * τ_test * bvd_inc
        conf_means[i] = μ_i
        Λ_prev += μ_i
        Λ_at_edges[i] = Λ_prev
        if t_edges[i] <= tests_edge + sqrt(eps(typeof(tests_edge)))
            bvd_tested_unit += bvd_inc
        end
    end

    for i in 1:n
        confirmed_cases[i] ~ safe_nbinomial(k, conf_means[i])
    end

    ## Tested-volume mean at `tests_edge`: τ · (BVD tested + background
    ## tested). The background arrives at constant rate λ_bg and is
    ## convolved against F_lab so only samples processed by `tests_edge`
    ## count.
    raw_bvd_tested = τ_test * bvd_tested_unit
    BVD_tested := isfinite(raw_bvd_tested) ?
                  max(raw_bvd_tested, eps(typeof(raw_bvd_tested))) :
                  eps(typeof(raw_bvd_tested))
    te = oftype(T, tests_edge)
    ## Background tested volume at `te`: constant-rate λ_bg arrivals
    ## convolved against the lab-delay CDF and right-truncated at `te`,
    ## which is ∫₀^te F_lab(te - u) du = ∫₀^te F_lab(v) dv, the closed
    ## form `_gamma_cdf_integral`.
    bg_integral = _gamma_cdf_integral(α_lab, θ_lab, te)
    raw_bg_tested = τ_test * λ_bg * bg_integral
    bg_tested := isfinite(raw_bg_tested) ?
                 max(raw_bg_tested, eps(typeof(raw_bg_tested))) :
                 eps(typeof(raw_bg_tested))

    expected_tested := BVD_tested + bg_tested
    ## Per-test positivity: s × BVD fraction in the tested pool. Exposed
    ## for comparison against the sitrep `Taux de positivité`.
    p_pos_raw = s_test * BVD_tested / expected_tested
    p_positive := isfinite(p_pos_raw) ?
                  clamp(p_pos_raw,
        eps(typeof(p_pos_raw)),
        one(p_pos_raw) - eps(typeof(p_pos_raw))) :
                  eps(typeof(p_pos_raw))

    ## Single tested-volume NegBinomial at the laboratory stream's own
    ## cut-off. Sampled unconditionally so it conditions on the data when
    ## present and is generated for posterior-predictive checks when
    ## `missing`. Confirmed counts are not re-observed conditional on it,
    ## so the two lab streams are not double-counted.
    tests_analysed ~ safe_nbinomial(k, expected_tested)

    expected_confirmed_total := Λ_at_edges[end]

    return (; expected_tested, p_positive,
        BVD_tested, bg_tested, s_test, τ_test, p_drc_per_bin,
        expected_confirmed_total,
        lab_delay_dist = f_lab,
        I_lab0_edges, Λ_at_edges)
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
onset-to-detection delay. The delay-convolution counterpart of
[`exports_deaths_model`](@ref): the per-edge intensity is the
onset-to-detection delay-survival expectation
[`expected_exports_deaths_delay`](@ref), where the rectangular detection
window of the McCabe path is replaced by the onset-to-detection survival
`f_det`. The pre-death survival stretch and per-day
[`daily_increment_kernel`](@ref) structure are unchanged; `onset_fraction`
(the incubation mgf) is folded into the expectation. `f_det` and
`daily_travellers` are supplied by [`exports_delay_model`](@ref) so the
two Uganda-side likelihoods share person-time.
"""
@model function exports_deaths_delay_model(
        export_deaths_daily::AbstractVector,
        growth_state, CFR::Real, delay_dist, p_uganda::Real;
        pre_start_deaths::Union{Missing, Integer} = 0,
        f_det,
        daily_travellers::Real,
        source_population::Real = ITURI_POPULATION,
        onset_fraction::Real = 1.0)
    cumulative = growth_state.cumulative
    T = growth_state.T
    q = daily_travellers / source_population
    n = length(export_deaths_daily)   # days from earliest death to cut-off

    ## Per-edge intensity is the onset-to-detection delay-survival
    ## expectation: the onset-to-death CDF weighted by the onset-to-
    ## detection survival, with `onset_fraction` (the incubation mgf)
    ## mapping the infection trajectory onto onsets.
    Λ(t) = expected_exports_deaths_delay(
        cumulative, f_det, delay_dist, CFR, p_uganda, q, t;
        onset_fraction = onset_fraction)

    ## Pre-death zero stretch as one Poisson observed at 0; `missing`
    ## generates it for predictive checks.
    pre = T - n > zero(T) ? Λ(T - n) : zero(T)
    pre_start_deaths ~ Poisson(max(pre, zero(pre)))

    ## Per-day means via the shared between-edge differencing, starting
    ## from the pre-death survival weight `pre`.
    Λ_at_edges = [Λ(T - n + i) for i in 1:n]
    μ_day = daily_increment_kernel(Λ_at_edges, pre)
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
onset-to-detection delay. The delay-convolution counterpart of
[`exports_detection_timing_model`](@ref): the at-risk export
person-time uses [`expected_exports_delay`](@ref) with the
onset-to-detection survival `f_det` instead of the rectangular window.
Adds the same one-sided `Pr(no export detected before t1)` Poisson
observation at zero. `f_det`, `daily_travellers` and `onset_fraction`
are supplied by [`exports_delay_model`](@ref). Passing `delta = missing`
makes the submodel a no-op.
"""
@model function exports_detection_timing_delay_model(
        growth_state, p_uganda::Real;
        delta::Union{Missing, Real},
        pre_detection_exports::Union{Missing, Integer} = 0,
        f_det,
        daily_travellers::Real,
        source_population::Real = ITURI_POPULATION,
        onset_fraction::Real = 1.0)
    if !ismissing(delta)
        r = growth_state.r
        T = growth_state.T
        t1 = T - delta
        q = daily_travellers / source_population
        survived_exports := t1 <= zero(T) ? zero(T) :
                            expected_exports_delay(r, p_uganda, q, t1,
            f_det; onset_fraction = onset_fraction)
        ## No detection before t1 as a Poisson observed at 0; `missing`
        ## generates it for predictive checks.
        pre_detection_exports ~ Poisson(max(survived_exports, zero(T)))
    end

    return (;)
end
