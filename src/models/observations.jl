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
        traveller = traveller_volume_model())
    cumulative = growth_state.cumulative
    T = growth_state.T

    window_state ~ to_submodel(window, false)
    w = window_state.w

    travel_state ~ to_submodel(traveller, false)
    daily_travellers = travel_state.daily_travellers

    q = daily_travellers / source_population
    expected_exports_T := expected_exports(cumulative, p_uganda, q, T, w;
        r = growth_state.r)

    exported_cases ~ Poisson(expected_exports_T)

    return (; w, daily_travellers, p_uganda,
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

`q_v` is the **severe-first** tested BVD share: the lab tests the
most-likely-BVD (severest / most obvious) suspects first, so `q` starts at
the early severe-cluster fraction `q0` (near 1) and relaxes to a broad-pool
baseline `qinf` as testing widens (see [`severe_first_share`](@ref) and
[`test_selection_model`](@ref)):

```math
q_v = q_\\infty + (q_0 - q_\\infty)\\, e^{-c_v / \\text{decay}},
\\qquad c_v = \\max(s_v - t_{report},\\ 0),
```

with `c_v` the elapsed time since testing (reporting) onset
`t_report = T − report_onset_offset`. `q` FALLS from `q0` to `qinf` and
LEVELS OFF there (no overshoot to zero). With `q0 ≈ 1` the first vintage
reads positivity `≈ s`, so the sensitivity is identified from the early
data and needs no lower floor; the plateau positivity
`s·qinf + (1−spec)(1−qinf)` holds the late vintages. `q0`, `qinf` and the
decay timescale are sampled by [`test_selection_model`](@ref).

`q_random_effect` adds a per-vintage partially-pooled logit-scale offset
to this baseline share (`q_v = logistic(logit(q_base,v) + σ_q·z_v)`, see
[`confirmed_q_re_model`](@ref)) so each window's positivity can fit the
non-monotone wobble while `s` stays fixed. It is on by default because the
observed positivity is non-monotone; pass `nothing` to recover the smooth
severe-first baseline.

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

`tests_analysed` is the fallback Binomial denominator for any `missing`
`samples_analysed` entry under predictive generation. The per-test
positivity `p_pos`, tested BVD share `q_cutoff` and baseline `q_baseline`
are exposed as derived quantities.

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
        test_selection = test_selection_model(),
        receipt_delay = lab_receipt_delay_model(),
        capacity_model = lab_capacity_model,
        capacity_centre::Real = 150.0,
        report_onset_offset::Union{Nothing, Real} = nothing,
        overdispersion = nothing,
        q_random_effect = confirmed_q_re_model,
        selection_clock::Symbol = :time,
        volume_scale::Real = 200.0,
        onset_fraction::Real = 1.0)
    sensitivity_state ~ to_submodel(test_sensitivity, false)
    specificity_state ~ to_submodel(test_specificity, false)
    selection_state ~ to_submodel(test_selection, false)
    receipt_state ~ to_submodel(receipt_delay, false)
    if overdispersion !== nothing
        overdispersion_state ~ to_submodel(overdispersion, false)
        φ_conf = overdispersion_state.φ_conf
    else
        φ_conf = nothing
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
    q0 = selection_state.q0
    decay_scale = selection_state.decay_scale
    qinf = selection_state.qinf
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
    ## composition, exposed as a diagnostic; the q-curve baseline `qinf` is a
    ## free parameter near the cut-off positivity, so the outbreak size stays
    ## pinned by the deaths / exports streams and the received counts rather
    ## than forced through this ratio.
    t_report = report_onset_offset === nothing ? zero(T) :
               max(T - report_onset_offset, zero(T))

    Tt = typeof(float(r) * float(λ_bg) * onset_fraction * float(q0) *
                float(decay_scale) * float(qinf))
    p_pos = Vector{Tt}(undef, n)
    q_at = Vector{Tt}(undef, n)
    qinf_count_at = Vector{Tt}(undef, n)
    Nsusp_at = Vector{Tt}(undef, n)
    Λ_at_edges = Vector{Tt}(undef, n)
    qinf_c = convert(Tt, qinf)
    ## Severe-first testing with modelled specificity, on cumulative counts.
    ## The lab tests the most-likely-BVD (severest / obvious) suspects first,
    ## so the BVD share of the analysed batch starts at the early
    ## severe-cluster fraction `q0` (near 1) and relaxes to the broad-pool
    ## baseline `qinf` as testing widens (`severe_first_share`):
    ##   q_v = qinf + (q0 − qinf)·exp(−c_v / decay_scale),
    ## with `c_v` the elapsed time since testing (reporting) onset. q FALLS
    ## from q0 to qinf and LEVELS OFF there (no overshoot to zero). The
    ## per-test positivity combines true and false positives,
    ##   p_pos = s·q + (1 − spec)·(1 − q),
    ## with `1 − spec` the false-positive rate. With q0 ≈ 1 the first vintage
    ## reads positivity ≈ s, identifying the sensitivity from the early data
    ## (no floor); the plateau s·qinf + (1−spec)(1−qinf) holds the late
    ## positivity. The binomial conditions on the observed analysed
    ## denominator.
    ## Cumulative analysed volume at each edge, for the volume-indexed
    ## clock. The severe-first share then relaxes with how many samples have
    ## been processed rather than calendar time, so a lab stall (no new
    ## analysed) does not advance the curve.
    analysed_cum_at = Vector{Tt}(undef, n)
    acc_a = zero(Tt)
    for i in 1:n
        acc_a += samples_analysed[i] === missing ? zero(Tt) :
                 convert(Tt, samples_analysed[i])
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
        ## Clock for the severe-first share: elapsed time since testing
        ## onset (`:time`) or cumulative analysed volume in units of
        ## `volume_scale` samples (`:volume`).
        c_i = selection_clock === :volume ?
              analysed_cum_at[i] / convert(Tt, volume_scale) :
              max(s_i - oftype(s_i, t_report), zero(s_i))
        q_base = severe_first_share(q0, qinf_c, c_i, decay_scale)
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
        ## positivity (a diagnostic).
        Λ_at_edges[i] = p_pos[i] *
                        (samples_analysed[i] === missing ? zero(Tt) :
                         convert(Tt, samples_analysed[i]))
    end

    ## New positives observed on the newly-analysed slice ΔA_v via a
    ## Binomial, the same between-vintage increment form as the reported and
    ## deaths streams (first window's increment is the cumulative, so a
    ## single vintage reduces to the cumulative Binomial). A `missing`
    ## denominator falls back to `tests_analysed` for predictive draws.
    for i in 1:n
        A_i = samples_analysed[i] === missing ?
              (tests_analysed === missing ? 0 : tests_analysed) :
              samples_analysed[i]
        if φ_conf === nothing
            confirmed_cases[i] ~ Binomial(Int(A_i), p_pos[i])
        else
            ## Beta-Binomial: positives disperse around ΔA·p with
            ## concentration φ (α = φ·p, β = φ·(1−p)), absorbing the
            ## laboratory reporting noise the plain Binomial cannot.
            α_bb = max(φ_conf * p_pos[i], eps(Tt))
            β_bb = max(φ_conf * (one(Tt) - p_pos[i]), eps(Tt))
            confirmed_cases[i] ~ BetaBinomial(Int(A_i), α_bb, β_bb)
        end
    end

    ## Received queue: the suspect backlog reaches the lab after a transport
    ## delay, so the cumulative received is the suspect backlog convolved with
    ## the receipt-delay kernel f_receipt. Confirmed uses a pooled (constant)
    ## p_drc, so a single representative value drives the continuous backlog.
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

    ## Received increments: the lab forwards a fraction τ_forward of the
    ## received backlog, so the per-window mean is τ_forward·ΔN_recv, observed
    ## through a NegBinomial sharing `k`. A `missing` entry is generated for
    ## predictive checks.
    ΔN_recv = daily_increment_kernel(N_recv_at)
    recv_means = Vector{Tt}(undef, n)
    for i in 1:n
        raw = convert(Tt, τ_forward) * ΔN_recv[i]
        recv_means[i] = isfinite(raw) ? max(raw, eps(Tt)) : eps(Tt)
        samples_received[i] ~ safe_nbinomial(k, recv_means[i])
    end

    ## Analysis throughput: the lab processes a fraction of the available
    ## received backlog set by its daily capacity κ_v over the window length
    ## Δt_v: μ_A = backlog·(1 − exp(−κ_v·Δt_v/backlog)) — capacity-limited when
    ## the backlog is large, the whole backlog when capacity is ample. Window
    ## 1's cumulative analysed is the initial condition (it accumulated over an
    ## unknown pre-window period); windows 2..n condition the capacity random
    ## walk on the genuine between-sitrep gaps. A `missing` entry is generated.
    analysed_cum = (n >= 1 && samples_analysed[1] !== missing) ?
                   convert(Tt, samples_analysed[1]) : zero(Tt)
    for i in 2:n
        Δt = max(oftype(T, t_edges[i]) - oftype(T, t_edges[i - 1]), zero(T))
        recv_cum = convert(Tt, τ_forward) * N_recv_at[i]
        backlog = max(recv_cum - analysed_cum, eps(Tt))
        cap = max(convert(Tt, κ[i]) * Δt, zero(Tt))
        μ_A_raw = backlog * (one(Tt) - exp(-cap / backlog))
        μ_A = isfinite(μ_A_raw) ? max(μ_A_raw, eps(Tt)) : eps(Tt)
        samples_analysed[i] ~ Poisson(μ_A)
        analysed_cum += samples_analysed[i] === missing ? μ_A :
                        convert(Tt, samples_analysed[i])
    end

    ## Per-test positivity at the lab cut-off (last edge ≤ tests_edge):
    ## exposed for comparison with the sitrep `Taux de positivité`.
    te_idx = something(
        findlast(
            i -> t_edges[i] <= tests_edge + sqrt(eps(typeof(tests_edge))),
            1:n), n)
    p_positive := p_pos[te_idx]
    q_cutoff := q_at[te_idx]
    q_baseline := qinf
    q_baseline_count := qinf_count_at[te_idx]
    τ_forward_out := τ_forward
    ## Daily analysis capacity at the cut-off vintage (samples/day).
    capacity_cutoff := κ[te_idx]

    ## Cumulative expected totals at the cut-off: sum the per-window
    ## increment means (confirmed positives, received samples).
    expected_confirmed_total := sum(Λ_at_edges)
    expected_received_total := sum(recv_means)

    return (; p_positive, p_pos, q_at, qinf, qinf_count_at, Nsusp_at,
        N_recv_at, recv_means, capacity = κ, s_test, spec_test, q0,
        decay_scale, τ_forward, p_drc_per_bin, expected_confirmed_total,
        expected_received_total, Λ_at_edges)
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
        ## giving Inf/Inf = NaN that would propagate (via `m_death·q`) into
        ## a NaN Binomial gradient in the confirmed-death stream and break
        ## the joint AD init. Fall back to zero share when non-finite.
        ratio = denom > zero(Tt) ? μ_bvd / denom : zero(Tt)
        q[i] = isfinite(ratio) ? clamp(ratio, zero(Tt), one(Tt)) : zero(Tt)
    end
    return q
end

"""
Laboratory-confirmed-deaths likelihood. The DRC sitrep front page reports
`Cumul décès parmi les confirmés`, the cumulative deaths that have been
laboratory-confirmed (17 at the 28 May cut-off). In a filovirus outbreak
deaths are tested aggressively (post-mortem / safe-burial swabs), so a
death is confirmed *because* it occurred; the natural denominator is the
suspected-death pool, not the confirmed-case cohort. The confirmed-death
increment is a Binomial thinning of the suspected-death increment:

```math
\\Delta D_{conf,v} \\sim
    \\mathrm{Binomial}(\\Delta D_{susp,v},\\ p_{death,conf}),
```

with `ΔD_susp,v` the between-vintage suspected-death increment (the
`death_history` increments [`deaths_model`](@ref) consumes) aligned with
`t_edges`.

There is no death-specific laboratory denominator (no analysed/tested
count for deaths), so `p_death_conf` cannot be identified as a free
parameter. It is *derived* from the confirmed-case BVD composition
`q_at_edges = μ_BVD / N_susp` (the count-implied BVD share among suspects
the confirmed-cases stream estimates from the shared lab machinery) times
a single tight-prior enrichment scalar `m_death` from
[`confirmed_death_enrichment_model`](@ref):

```math
p_{death,conf}(v) = \\mathrm{clamp}(m_{death} \\cdot q(v),\\ 0,\\ 1).
```

The stream is therefore almost entirely an *input* from the confirmed-case
model: it adds at most the one scalar `m_death` (centred on 1, no
enrichment), the most the single informative point (17/246) supports. Set
the enrichment prior SD to zero to tie deaths directly to the case
composition with no free parameter. `confirmed_deaths` accepts `missing`
entries for posterior-predictive generation.
"""
@model function confirmed_deaths_model(
        confirmed_deaths::AbstractVector,
        susp_death_increments::AbstractVector,
        q_at_edges::AbstractVector;
        confirmed_death_enrichment = confirmed_death_enrichment_model())
    n = length(q_at_edges)
    n == length(confirmed_deaths) ||
        error("confirmed_deaths length must match q_at_edges (got " *
              "$(length(confirmed_deaths)) vs $n)")
    n == length(susp_death_increments) ||
        error("susp_death_increments length must match q_at_edges " *
              "(got $(length(susp_death_increments)) vs $n)")

    enrich_state ~ to_submodel(confirmed_death_enrichment, false)
    m_death = enrich_state.m_death

    Tt = typeof(float(m_death) * float(first(q_at_edges)))
    p_death_conf = Vector{Tt}(undef, n)
    for i in 1:n
        ## Keep the confirmation probability strictly inside (0, 1): at the
        ## saturated boundary the Binomial logpdf for the partial 17/246
        ## point is log(0) and its reverse-mode gradient is NaN, which
        ## breaks the joint AD init. The eps margins leave the composition
        ## link unchanged away from saturation.
        p_death_conf[i] = clamp(convert(Tt, m_death) *
                                convert(Tt, q_at_edges[i]),
            eps(Tt), one(Tt) - eps(Tt))
    end

    for i in 1:n
        ΔD = susp_death_increments[i] === missing ? 0 :
             Int(susp_death_increments[i])
        confirmed_deaths[i] ~ Binomial(ΔD, p_death_conf[i])
    end

    p_death_conf_cutoff := p_death_conf[end]
    expected_confirmed_deaths_total := sum(p_death_conf[i] *
                                           (susp_death_increments[i] ===
                                            missing ? zero(Tt) :
                                            convert(Tt,
                                               susp_death_increments[i]))
    for i in 1:n)

    return (; m_death, p_death_conf, p_death_conf_cutoff,
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
