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
background whose cumulative is the saturating ramp `μ_bg(s)` (see
[`BackgroundRamp`](@ref) and [`test_positivity_model`](@ref)), so each
bin's background increment is `μ_bg(s_k) − μ_bg(s_{k−1})`. With `Δλ = 0`
this is the constant `λ0 · Δt`. The latent trajectory is cumulative
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
    τ_test = test_positivity_state.τ_test
    f_rep = report_state.dist
    r = growth_state.r
    ## Anchor the background ramp to reporting onset: the background clock
    ## starts at `t_report = T − report_onset_offset` (clamped at 0), so
    ## non-BVD suspects accrue only once case-finding has begun, and the
    ## ramp is still climbing across the late lab window. `nothing`
    ## (the default) keeps the seeding-anchored ramp (`t_report = 0`),
    ## which with `Δλ = 0` recovers the constant background exactly.
    t_report = report_onset_offset === nothing ? zero(growth_state.T) :
               max(growth_state.T - report_onset_offset,
        zero(growth_state.T))
    bg = background_ramp(test_positivity_state.λ0,
        test_positivity_state.Δλ, test_positivity_state.scale, t_report)

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
        ## Background increment is the cumulative ramp difference
        ## μ_bg(s_k) − μ_bg(s_{k−1}); with Δλ = 0 this is λ0·(s_k − s_{k−1}),
        ## the old constant-rate `λ_bg · Δt`.
        μ_bg_k = bg_cumulative(bg, s_k)
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

    return (; p_drc_per_bin, λ_bg, bg, τ_test,
        expected_reports_total, positivity,
        report_delay_dist = f_rep,
        μ_BVD0_at_edges, Λ_at_edges)
end

"""
Laboratory pipeline likelihood. Models the lab-confirmed cases over time
and, where available, the testing volume that gates them.

`confirmed_cases` (`Cumul positifs`) is a per-vintage *cumulative* vector
aligned with `t_edges`. Each vintage's cumulative confirmed count is
observed conditional on its cumulative analysed denominator through a
Binomial (equations (23)-(24) of the walkthrough):

```math
C_v \\sim \\mathrm{Binomial}(A_v,\\ p_{pos,v}),
\\qquad
p_{pos,v} = \\frac{s_{test}\\, \\text{BVD}_{tested}(A_v)}{A_v},
```

where `A_v` is the cumulative analysed count to edge `v`
(`samples_analysed`, a *known* denominator, not modelled) and
`BVD_tested(A_v)` is the cumulative BVD among those `A_v` analysed
samples under **priority (triage) testing**: the lab analyses the suspect
pool in order of pre-test BVD probability, creaming the most-likely-BVD
samples first. With `B_v` BVD available, `N_v` the total available pool
and strength `κ` (see [`test_priority_model`](@ref),
[`priority_bvd_tested`](@ref)),

```math
\\text{BVD}_{tested}(A_v) = B_v\\,\\bigl(1 - (1 - A_v/N_v)^{\\kappa}\\bigr),
\\qquad \\text{bg}_{tested}(A_v) = A_v - \\text{BVD}_{tested}(A_v).
```

The cumulative BVD tested saturates toward `B_v` as `A_v` grows, so the
confirmed positives plateau while the analysed volume climbs and the
per-test positivity falls — the observed saturating-positives /
falling-positivity pattern that a constant or ramped background cannot
produce. `κ = 1` is proportional sampling, `BVD_tested = B·A/N`, giving
`p_pos = s·B/N` (the old available-BVD-share model); `Δκ = 0` is the
backward-compatible reduction. The cumulative form is robust to vintages
with zero newly-analysed samples (e.g. the 25 May Ituri stoppage).

The available pool is the FULL cumulative lab-eligible backlog to each
vintage (the high-suspicion samples the lab triages first are those
accumulated since seeding, not just the last few days): `B_v` is
`p_{DRC,v}` times
``I_{lab,0}(s_v) = \\int_0^{s_v} \\mu_{BVD,0}(u)\\, f_{lab}(s_v - u)\\,du``,
`μ_BVD,0` the unit-ascertainment onset-to-report BVD cumulative (kernel
`f_rep` from [`reported_cases_model`](@ref)) and `s_test` the PCR
sensitivity. A [`DailyBVDTrajectory`](@ref) precomputes `μ_BVD,0` once on
a shared Gauss-Legendre node set. The background available is the
constant-rate non-BVD backlog
``\\int_0^{s_v} \\lambda_{bg}\\, F_{lab}(s_v - u)\\,du`` (via
[`bg_tested_integral`](@ref); the ramp is dormant, `Δλ = 0`), so the
time-variation lives in the testing order / exhaustion, not the
background rate. Anchoring the pool to the short reporting window instead
starves it (the analysed counts far exceed an 8-day BVD slice), so
`report_onset_offset` is retained only for API symmetry with the
reported/ramp path and does not shrink the priority pool.

`samples_analysed` is the per-vintage analysed denominator vector aligned
with `t_edges`. Pass a vector of `missing` entries (with a non-missing
fallback `tests_analysed` for the denominators) to generate
posterior-predictive confirmed draws.

`tests_analysed` (`Cumul échantillons analysés`) is a single cumulative
count observed at its own elapsed time `tests_edge` (which may lag the
case cut-off `T` if lab reporting stops earlier). When present it enters
as one NegBinomial term with mean
``\\tau\\,(\\text{BVD}_\\text{tested} + \\text{bg}_\\text{tested})``
accumulated to `tests_edge`. It also supplies the binomial denominator
for any `missing` `samples_analysed` entry under predictive generation.
The per-test positivity `s · BVD_tested / (BVD_tested + bg_tested)` is
exposed as a derived quantity for comparison with the sitrep figure. Pass
`tests_analysed = missing` to drop the tested-volume NegBinomial.

A single confirmed observation (`length 1`, `t_edges = [T]`, single
`samples_analysed`) reduces to the cumulative confirmed Binomial.
"""
@model function confirmed_cases_model(
        confirmed_cases::AbstractVector,
        samples_analysed::AbstractVector,
        tests_analysed::Union{Missing, Integer},
        growth_state, k::Real,
        p_drc_per_bin::AbstractVector, bg::BackgroundRamp, τ_test::Real,
        f_rep,
        t_edges::AbstractVector, tests_edge::Real;
        lab_delay = lab_delay_model(),
        test_sensitivity = test_sensitivity_model(),
        test_priority = test_priority_model(),
        report_onset_offset::Union{Nothing, Real} = nothing,
        onset_fraction::Real = 1.0)
    lab_state ~ to_submodel(lab_delay, false)
    sensitivity_state ~ to_submodel(test_sensitivity, false)
    priority_state ~ to_submodel(test_priority, false)
    f_lab = lab_state.dist
    s_test = sensitivity_state.s_test
    κ_priority = priority_state.κ_priority
    r = growth_state.r
    T = growth_state.T
    α_lab = lab_state.α
    θ_lab = lab_state.θ
    ## The available suspect backlog the lab triages is the FULL cumulative
    ## lab-eligible pool to each vintage (BVD accumulated since seeding plus
    ## the non-BVD background), not just the testing-window increment:
    ## anchoring the pool to the short reporting window starves it (the
    ## analysed counts 211-403 far exceed an 8-day BVD slice), collapsing
    ## positivity. `report_onset_offset` is retained for API symmetry with
    ## the reported/ramp path but does not shrink the priority pool.
    background_pool = background_ramp(bg.λ0, bg.Δλ, bg.scale, zero(T))

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

    ## Unit-ascertainment μ_BVD,0 at the shared Gauss-Legendre nodes over
    ## [0, T]; the outer quadrature of I_lab,0 is reused across every bin
    ## edge. Per-bin ascertainment is then applied on the between-edge
    ## I_lab,0 increment so each bin's mean tracks its own draw.
    ## The latent trajectory is cumulative *infections*, so `onset_fraction`
    ## (the incubation mgf) maps it onto onsets before the onset-to-report
    ## and report-to-lab convolutions; it scales the whole BVD lab
    ## cumulative, hence both the confirmed counts and the BVD tested
    ## volume below.
    ## Full cumulative lab-eligible BVD at each bin edge (from seeding).
    trajectory = DailyBVDTrajectory(T, r, f_rep)
    I_lab0_edges = onset_fraction .*
                   delay_convolution(trajectory, t_edges, f_lab)

    Tt = eltype(I_lab0_edges)
    p_pos = Vector{Tt}(undef, n)
    Λ_at_edges = Vector{Tt}(undef, n)
    ## Priority/exhaustion testing on cumulative counts. The lab analyses
    ## the suspect backlog in order of pre-test BVD probability, creaming
    ## the most-likely-BVD samples first, so the BVD share of the analysed
    ## batch is highest early and falls as the backlog drains. The available
    ## pool is the FULL cumulative lab-eligible backlog to each vintage:
    ##  - BVD available  B_i = p_DRC · I_lab,0(s_i)  (all BVD eligible to
    ##    date — the high-suspicion samples the lab triages first);
    ##  - background available bg_i = ∫_0^{s_i} λ_bg F_lab du (constant-rate
    ##    non-BVD backlog, lab-convolved, via `bg_tested_integral`);
    ##  - total pool N_i = B_i + bg_i.
    ## The observed cumulative analysed `A_i` (data) draws from this pool
    ## with BVD front-loaded by strength κ: BVD_tested = B·(1−(1−A/N)^κ)
    ## (`priority_bvd_tested`), background-tested = A − BVD_tested, and the
    ## per-test positivity = s · BVD_tested / A. When the BVD backlog is
    ## exhausted (A ≥ N) all of B is tested regardless of κ, so positives
    ## plateau at s·B while A climbs and positivity = s·B/A falls — the
    ## exhaustion limit that reproduces the saturating positives. κ = 1
    ## (priority off) gives BVD_tested = B·A/N, i.e. positivity = s·B/N,
    ## the proportional (no-priority) reduction. The binomial conditions on
    ## the observed analysed denominator. The pool is the lab-eligible
    ## suspect count (not thinned by `τ_test`, which stays the
    ## suspected-stream sampling fraction in `reported_cases_model`).
    bvd_tested_unit = zero(Tt)
    bvd_avail_te = zero(Tt)
    bg_avail_te = zero(Tt)
    N_te = zero(Tt)
    for i in 1:n
        raw_bvd = p_drc_per_bin[i] * I_lab0_edges[i]
        bvd_avail = isfinite(raw_bvd) ? max(raw_bvd, zero(Tt)) : zero(Tt)
        bg_avail = max(
            bg_tested_integral(background_pool, α_lab, θ_lab,
                oftype(T, t_edges[i])),
            zero(Tt))
        N_i = bvd_avail + bg_avail
        A_i = samples_analysed[i] === missing ? N_i :
              convert(Tt, samples_analysed[i])
        bvd_tested = priority_bvd_tested(bvd_avail, N_i, A_i, κ_priority)
        bvd_tested = max(bvd_tested, eps(Tt))
        p_raw = A_i > zero(Tt) ? s_test * bvd_tested / A_i : zero(Tt)
        p_pos[i] = isfinite(p_raw) ?
                   clamp(p_raw, eps(Tt), one(Tt) - eps(Tt)) : eps(Tt)
        ## Expected cumulative confirmed: analysed denominator × cumulative
        ## positivity, a derived diagnostic.
        Λ_at_edges[i] = p_pos[i] *
                        (samples_analysed[i] === missing ? zero(Tt) :
                         convert(Tt, samples_analysed[i]))
        if t_edges[i] <= tests_edge + sqrt(eps(typeof(tests_edge)))
            bvd_tested_unit = bvd_tested
            bvd_avail_te = bvd_avail
            bg_avail_te = bg_avail
            N_te = N_i
        end
    end

    ## Confirmed counts observed conditional on the analysed denominator
    ## via a cumulative Binomial (equation (24)). `A_v` is data, not
    ## modelled, which removes the multiplicative ridge of the
    ## NegBinomial-increment form. A `missing` denominator falls back to
    ## the cumulative `tests_analysed` so predictive draws still have an
    ## integer denominator.
    for i in 1:n
        A_i = samples_analysed[i] === missing ?
              (tests_analysed === missing ? 0 : tests_analysed) :
              samples_analysed[i]
        confirmed_cases[i] ~ Binomial(Int(A_i), p_pos[i])
    end

    ## Tested-volume mean at `tests_edge`: the total available pool
    ## N = BVD available + background available (lab-eligible counts
    ## accrued over the testing window from `t_report`); this is the
    ## expected processed volume the analysed-count NegBinomial conditions
    ## on. `BVD_tested` is the priority-creamed BVD among the analysed
    ## samples at `tests_edge`.
    BVD_tested := isfinite(bvd_tested_unit) ?
                  max(bvd_tested_unit, eps(typeof(bvd_tested_unit))) :
                  eps(typeof(bvd_tested_unit))
    bg_tested := isfinite(bg_avail_te) ?
                 max(bg_avail_te, eps(typeof(bg_avail_te))) :
                 eps(typeof(bg_avail_te))

    expected_tested := isfinite(N_te) ? max(N_te, eps(typeof(N_te))) :
                       eps(typeof(N_te))
    ## Per-test positivity at `tests_edge`: s × priority BVD share of the
    ## analysed batch. Exposed for comparison with the sitrep
    ## `Taux de positivité`. Uses the observed analysed count at the
    ## tested-volume edge when present, else the modelled pool N.
    A_te = tests_analysed === missing ? N_te : convert(Tt, tests_analysed)
    p_pos_raw = A_te > zero(Tt) ? s_test * bvd_tested_unit / A_te :
                zero(Tt)
    p_positive := isfinite(p_pos_raw) ?
                  clamp(p_pos_raw,
        eps(typeof(p_pos_raw)),
        one(p_pos_raw) - eps(typeof(p_pos_raw))) :
                  eps(typeof(p_pos_raw))

    ## Single tested-volume NegBinomial at the laboratory stream's own
    ## cut-off. Sampled unconditionally so it conditions on the data when
    ## present and is generated for posterior-predictive checks when
    ## `missing`. The confirmed binomial conditions on the *observed*
    ## analysed denominator, not this modelled count, so the two streams
    ## are not double-counted.
    tests_analysed ~ safe_nbinomial(k, expected_tested)

    expected_confirmed_total := Λ_at_edges[end]

    return (; expected_tested, p_positive, p_pos,
        BVD_tested, bg_tested, s_test, τ_test, p_drc_per_bin,
        expected_confirmed_total,
        lab_delay_dist = f_lab,
        I_lab0_edges, Λ_at_edges)
end

"""
Constant-capacity laboratory-throughput likelihood (issue #174). Models
the three lab streams — samples received, samples analysed and confirmed
positives — through a single daily FIFO backlog drained at a constant
capacity `κ`, rather than a fixed onset-to-lab Gamma delay. This replaces
the lab-delay convolution of [`confirmed_cases_model`](@ref) for the
confirmed/analysed path: the backlog provides the lag, so the
report-to-lab delay and the `p_DRC · s · τ` scale no longer compete for
the same turnaround (the multiplicative ridge that breaks the per-vintage
NegBinomial / binomial confirmed fit, #163).

The daily generative process over `1:D` (`D = round(T)`):

- *Arrivals.* The suspected stream feeds the lab. Daily BVD-suspected
  onsets-to-report are the increments of the Gamma closed-form
  ``p_{DRC}\\,\\mu_{BVD,0}(d)`` (kernel `f_rep`, shared with
  [`reported_cases_model`](@ref)); daily non-BVD background arrives at the
  saturating-ramp rate `λ_bg(d)` (the ramp cumulative increment over day
  `d`; constant `λ0` when `Δλ = 0`). A receipt fraction `τ_recv` thins
  both into
  samples *received* by the lab on day `d`.
- *Backlog drain.* Received samples join a FIFO backlog drained at `κ`
  per day, `processed_d = min(κ_d, backlog_{d-1} + arrivals_d)`, with the
  initial backlog fixed to zero and `κ_d = 0` on the externally-supplied
  `stoppage_days` (the 25 May Ituri stoppage, documented in
  `data/insp_sitrep_scanned.csv`). BVD and background are drained in
  proportion to their current backlog shares (a well-mixed approximation
  to strict FIFO: smooth for AD and adequate for cumulative totals).
- *Observations.* At each vintage edge the cumulative received and
  cumulative processed are NegBinomial; cumulative positives are observed
  conditional on the *observed* cumulative analysed denominator through a
  Binomial with success `s_test · q_v`, where `q_v` is the BVD share of
  the cumulative processed pool.

```math
R_v \\sim \\mathrm{NegBinomial}(\\textstyle\\sum_{d \\le v} \\text{recv}_d, k),
\\quad
A_v \\sim \\mathrm{NegBinomial}(\\textstyle\\sum_{d \\le v} \\text{proc}_d, k),
\\quad
C_v \\sim \\mathrm{Binomial}(A_v^{obs},\\ s_{test}\\, q_v).
```

`samples_received`, `samples_analysed` and `confirmed_cases` are
per-vintage *cumulative* vectors aligned with `t_edges`; `missing`
entries are generated for predictive checks (the binomial denominator
then falls back to the modelled cumulative processed, rounded). A single
edge reduces each stream to its cumulative single-total likelihood.
"""
@model function lab_throughput_model(
        samples_received::AbstractVector,
        samples_analysed::AbstractVector,
        confirmed_cases::AbstractVector,
        growth_state, k::Real,
        p_drc_per_bin::AbstractVector, bg::BackgroundRamp, τ_recv::Real,
        f_rep,
        κ_lab::Real, t_edges::AbstractVector;
        test_sensitivity = test_sensitivity_model(),
        stoppage_days::AbstractVector = Int[],
        onset_fraction::Real = 1.0)
    sensitivity_state ~ to_submodel(test_sensitivity, false)
    s_test = sensitivity_state.s_test
    r = growth_state.r
    T = growth_state.T

    n = length(t_edges)
    (n == length(p_drc_per_bin) == length(samples_received) ==
     length(samples_analysed) == length(confirmed_cases)) ||
        error("lab_throughput_model: stream lengths must match t_edges")

    Tt = typeof(float(r) * float(τ_recv) * float(κ_lab) * onset_fraction)
    D = max(1, Int(round(T)))
    ## Map each vintage edge to its daily index (1..D). Edges are elapsed
    ## times; the bin a vintage falls in is `round(edge)`.
    edge_day = [clamp(Int(round(t_edges[i])), 1, D) for i in 1:n]
    ## Per-vintage ascertainment broadcast onto its daily window; days up
    ## to a vintage's edge use that vintage's `p_drc`.
    p_for_day = Vector{Tt}(undef, D)
    let vi = 1
        for d in 1:D
            while vi < n && d > edge_day[vi]
                vi += 1
            end
            p_for_day[d] = convert(Tt, p_drc_per_bin[vi])
        end
    end

    ## Daily BVD-suspected onset-to-report cumulative (Gamma closed form),
    ## mapped onto onsets by `onset_fraction`. Differenced to daily
    ## arrivals below.
    μ_bvd_cum = Tt[d <= 0 ? zero(Tt) :
                   onset_fraction *
                   delay_convolution(one(Tt), r, oftype(T, d), f_rep)
                   for d in 1:D]

    ## FIFO backlog drained at constant capacity κ (zero on stoppage
    ## days). BVD and background tracked separately; drained in proportion
    ## to current backlog shares (well-mixed approximation to strict FIFO).
    stopset = Set(Int.(stoppage_days))
    bl_bvd = zero(Tt)
    bl_bg = zero(Tt)
    recv_cum = zero(Tt)
    proc_cum = zero(Tt)
    proc_bvd_cum = zero(Tt)
    recv_at = Vector{Tt}(undef, n)
    proc_at = Vector{Tt}(undef, n)
    q_at = Vector{Tt}(undef, n)
    μ_prev = zero(Tt)
    μ_bg_prev = zero(Tt)
    vi = 1
    for d in 1:D
        Δμ = max(μ_bvd_cum[d] - μ_prev, zero(Tt))
        μ_prev = μ_bvd_cum[d]
        a_bvd = τ_recv * p_for_day[d] * Δμ
        ## Daily background arrivals = ramp cumulative increment over day d,
        ## μ_bg(d) − μ_bg(d−1). With Δλ = 0 this is the constant λ0 per day.
        μ_bg_d = bg_cumulative(bg, oftype(T, d))
        Δμ_bg = max(μ_bg_d - μ_bg_prev, zero(Tt))
        μ_bg_prev = μ_bg_d
        a_bg = τ_recv * Δμ_bg
        recv_cum += a_bvd + a_bg
        bl_bvd += a_bvd
        bl_bg += a_bg
        total = bl_bvd + bl_bg
        κ_d = (d in stopset) ? zero(Tt) : convert(Tt, κ_lab)
        proc = min(κ_d, total)
        frac = total > zero(Tt) ? proc / total : zero(Tt)
        dr_bvd = frac * bl_bvd
        bl_bvd -= dr_bvd
        bl_bg -= frac * bl_bg
        proc_cum += proc
        proc_bvd_cum += dr_bvd
        ## Record cumulative pools at every vintage edge landing on day d.
        while vi <= n && edge_day[vi] == d
            recv_at[vi] = max(recv_cum, eps(Tt))
            proc_at[vi] = max(proc_cum, eps(Tt))
            q_raw = proc_cum > zero(Tt) ? proc_bvd_cum / proc_cum : zero(Tt)
            q_at[vi] = clamp(q_raw, eps(Tt), one(Tt) - eps(Tt))
            vi += 1
        end
    end

    ## Received and analysed as NegBinomial on the cumulative pools.
    for i in 1:n
        samples_received[i] ~ safe_nbinomial(k, recv_at[i])
        samples_analysed[i] ~ safe_nbinomial(k, proc_at[i])
    end

    ## Positives conditional on the *observed* cumulative analysed
    ## denominator via a Binomial with success s · q (BVD share of the
    ## processed pool). A `missing` denominator falls back to the modelled
    ## cumulative processed so predictive draws still have a denominator.
    p_pos = Vector{Tt}(undef, n)
    for i in 1:n
        p_raw = s_test * q_at[i]
        p_pos[i] = clamp(p_raw, eps(Tt), one(Tt) - eps(Tt))
        A_i = samples_analysed[i] === missing ?
              Int(round(proc_at[i])) : Int(samples_analysed[i])
        confirmed_cases[i] ~ Binomial(A_i, p_pos[i])
    end

    expected_received_total := recv_at[end]
    expected_analysed_total := proc_at[end]
    p_positive := p_pos[end]

    return (; recv_at, proc_at, q_at, p_pos, s_test, κ_lab,
        expected_received_total, expected_analysed_total, p_positive)
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
