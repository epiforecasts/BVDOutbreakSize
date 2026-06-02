## Tests for the explicit infection→detection delay export expectations
## (`expected_exports_delay`, `expected_exports_deaths_delay`,
## `combined_delay`) and the delay-mechanism `bvd_joint` defaults. The
## exports stream is travel-gated, so its at-risk clock runs from
## infection: the detection delay is incubation ⊕ onset-to-report,
## moment-matched to one Gamma. The delay forms generalise the McCabe
## rectangular detection window and reduce to it exactly as the
## infection→detection delay collapses to a point mass.

@testitem "expected_exports_delay matches the manual reconstruction" begin
    using Distributions: Gamma
    using BVDOutbreakSize: expected_exports_delay, delay_convolution,
                           integrate, _exp_cumulative_integral
    r = 0.05
    p = 0.25
    q = 1871 / 4_392_200
    T = 90.0
    f_det = Gamma(2.0, 4.5)

    got = expected_exports_delay(r, p, q, T, f_det)

    ## Independent rebuild of p·q·(∫₀ᵀ C − ∫₀ᵀ conv(f_det)).
    conv_integral = integrate(
        s -> delay_convolution(1.0, r, s, f_det), 0.0, T)
    want = p * q * (_exp_cumulative_integral(r, 0.0, T) - conv_integral)

    @test got ≈ want rtol = 1e-8
    @test got > 0
end

@testitem "combined_delay moment-matches the sum of two Gammas" begin
    using Distributions: Gamma, mean, var
    using BVDOutbreakSize: combined_delay

    a = Gamma(1.1, 5.7)   # incubation-like
    b = Gamma(2.5, 4.5)   # onset-to-report-like
    c = combined_delay(a, b)

    @test c isa Gamma
    ## Mean and variance add; the matched Gamma reproduces both.
    @test mean(c) ≈ mean(a) + mean(b) rtol = 1e-12
    @test var(c) ≈ var(a) + var(b) rtol = 1e-12
end

@testitem "expected_exports_delay effective window is the delay mean" begin
    using Distributions: Gamma, mean
    using BVDOutbreakSize: expected_exports_delay, combined_delay,
                           delay_convolution, integrate,
                           _exp_cumulative_integral

    ## The effective export window (∫C − ∫conv)/C(T) is the
    ## growth-weighted at-risk person-time per current infection. As
    ## r → 0 it equals the mean infection→detection delay (the at-risk
    ## pool is everyone infected within one mean delay of the cut-off).
    ## This is the ~17.5-day infection→detection delay (incubation ⊕
    ## onset-to-report), restoring McCabe's infection→detection window
    ## meaning rather than the ~8-day onset-to-report×os the flat-`os`
    ## bug produced.
    incubation = Gamma(1.1, 5.7)         # mean ≈ 6.27 d
    report = Gamma(2.5, 4.5)             # mean ≈ 11.25 d
    f_det = combined_delay(incubation, report)
    expected_mean = mean(incubation) + mean(report)   # ≈ 17.5 d

    eff_window(r, T) =
        let
            cum = _exp_cumulative_integral(r, 0.0, T)
            conv = integrate(s -> delay_convolution(1.0, r, s, f_det), 0.0, T)
            (cum - conv) / exp(r * T)
        end

    ## r → 0 limit recovers the full mean delay.
    @test eff_window(1e-4, 2_000.0) ≈ expected_mean rtol = 5e-2

    ## At the BVD growth rate the growth weighting shortens the window,
    ## but it must still substantially exceed the ~8-day window the old
    ## flat-os onset-to-report mechanism implied (onset-to-report mean ×
    ## os ≈ 11.25 × 0.72 ≈ 8 d).
    @test eff_window(0.05, 400.0) > 10.0
end

@testitem "expected_exports_delay reduces to the McCabe window" begin
    using Distributions: Gamma
    using BVDOutbreakSize: expected_exports_delay, expected_exports
    r = 0.05
    p = 0.25
    q = 1871 / 4_392_200
    T = 90.0
    w = 9.0
    cumulative = s -> exp(r * s)

    mccabe = expected_exports(cumulative, p, q, T, w; r = r)

    ## As the onset-to-detection Gamma collapses toward a point mass at
    ## `w` (CV = 1/√α → 0) the delay survival becomes the top-hat
    ## 1{t-s<w}, so the delay form converges to the window form. The
    ## tolerance tightens as the CV shrinks.
    rels = Float64[]
    for α in (20.0, 80.0, 320.0)
        f_det = Gamma(α, w / α)
        got = expected_exports_delay(r, p, q, T, f_det)
        push!(rels, abs(got - mccabe) / mccabe)
    end
    @test rels[end] < 1e-2
    @test rels[end] < rels[1]
end

@testitem "expected_exports_delay grows with elapsed time" begin
    using Distributions: Gamma
    using BVDOutbreakSize: expected_exports_delay
    r = 0.05
    f_det = Gamma(2.0, 4.5)
    f(t) = expected_exports_delay(r, 0.25, 1871 / 4_392_200, t, f_det)
    @test f(60.0) < f(90.0) < f(120.0)
    ## Always strictly positive (clamped), even at t = 0.
    @test f(0.0) > 0
end

@testitem "expected_exports_delay scales linearly with p and q" begin
    using Distributions: Gamma
    using BVDOutbreakSize: expected_exports_delay
    r = 0.05
    p = 0.25
    q = 1871 / 4_392_200
    T = 90.0
    f_det = Gamma(2.0, 4.5)

    base = expected_exports_delay(r, p, q, T, f_det)
    @test expected_exports_delay(r, 0.5 * p, q, T, f_det) ≈ 0.5 * base rtol = 1e-8
    @test expected_exports_delay(r, p, 0.5 * q, T, f_det) ≈ 0.5 * base rtol = 1e-8
end

@testitem "expected_exports_deaths_delay matches the manual convolution" begin
    using Distributions: Gamma
    using BVDOutbreakSize: expected_exports_deaths_delay, integrate,
                           _gamma_cdf, _delay_scale
    r = 0.05
    cumulative = s -> exp(r * s)
    f_det = Gamma(2.0, 4.5)
    f_death = Gamma(4.3, 2.6)
    CFR = 0.30
    p = 0.25
    q = 1871 / 4_392_200
    T = 90.0

    got = expected_exports_deaths_delay(
        cumulative, f_det, f_death, CFR, p, q, T)

    ## Independent reconstruction of the documented integrand:
    ## CFR·p·q·∫₀ᵀ C(s)·ccdf(f_det, T-s)·F_death(T-s) ds.
    g = s -> begin
        a = T - s
        a <= 0 ? 0.0 :
        cumulative(s) * (1 - _gamma_cdf(2.0, 4.5, a)) *
        _gamma_cdf(4.3, 2.6, a)
    end
    scale = _delay_scale(f_det) + _delay_scale(f_death)
    want = CFR * p * q * integrate(g, 0.0, T, scale)

    @test got ≈ want rtol = 1e-8
    @test got > 0
end

@testitem "expected_exports_deaths_delay grows with elapsed time" begin
    using Distributions: Gamma
    using BVDOutbreakSize: expected_exports_deaths_delay
    r = 0.05
    cumulative = s -> exp(r * s)
    f_det = Gamma(2.0, 4.5)
    f_death = Gamma(4.3, 2.6)
    f(t) = expected_exports_deaths_delay(
        cumulative, f_det, f_death, 0.30, 0.25, 1871 / 4_392_200, t)
    @test f(60.0) < f(90.0) < f(120.0)
    @test f(0.0) > 0
end

@testitem "expected_exports_deaths_delay reduces to the McCabe window" begin
    using Distributions: Gamma
    using BVDOutbreakSize: expected_exports_deaths_delay,
                           expected_exports_deaths
    r = 0.05
    cumulative = s -> exp(r * s)
    f_death = Gamma(4.3, 2.6)
    CFR = 0.30
    p = 0.25
    q = 1871 / 4_392_200
    T = 90.0
    w = 9.0

    mccabe = expected_exports_deaths(
        cumulative, f_death, CFR, p, q, T, w)

    ## As the infection→detection Gamma collapses toward a point mass at
    ## `w` the detection survival becomes the top-hat 1{t-s<w}, so the
    ## delay form converges to the window form (with the same death CDF
    ## `f_death` in both).
    rels = Float64[]
    for α in (20.0, 80.0, 320.0)
        f_det = Gamma(α, w / α)
        got = expected_exports_deaths_delay(
            cumulative, f_det, f_death, CFR, p, q, T)
        push!(rels, abs(got - mccabe) / mccabe)
    end
    @test rels[end] < 1e-2
    @test rels[end] < rels[1]
end

@testitem "export deaths timed from infection shift later than onset" begin
    using Distributions: Gamma, mean
    using BVDOutbreakSize: expected_exports_deaths_delay, combined_delay
    r = 0.05
    cumulative = s -> exp(r * s)
    f_det = Gamma(2.0, 4.5)
    incubation = Gamma(1.1, 5.7)          # mean ≈ 6.27 d
    onset_to_death = Gamma(4.3, 2.6)      # bare onset-to-death

    ## Timing death from infection (incubation ⊕ onset-to-death) pushes
    ## the death CDF later by one incubation period, so fewer export
    ## deaths have accrued by any given cut-off than under the bare
    ## onset-to-death delay (which omits incubation and fires too early).
    infection_to_death = combined_delay(incubation, onset_to_death)
    @test mean(infection_to_death) > mean(onset_to_death)

    args = (cumulative, f_det)
    tail = (0.30, 0.25, 1871 / 4_392_200)
    for T in (40.0, 70.0, 100.0)
        with_inc = expected_exports_deaths_delay(
            args..., infection_to_death, tail..., T)
        without_inc = expected_exports_deaths_delay(
            args..., onset_to_death, tail..., T)
        ## Later death timing => strictly fewer accrued export deaths.
        @test with_inc < without_inc
    end
end

## The closed-form Gamma-CDF path must differentiate: gradients of both
## delay expectations w.r.t. the delay (α, θ) and growth r must be finite
## and match central finite differences. Mirrors the Mooncake /
## FiniteDifferences gradient test in test_exports_death_timing.jl.

@testitem "expected_exports_delay has correct gradients" tags=[:ad] begin
    using Distributions: Gamma
    using FiniteDifferences: central_fdm, grad
    using Mooncake: Mooncake
    using BVDOutbreakSize: expected_exports_delay, combined_delay
    p = 0.25
    q = 1871 / 4_392_200
    T = 90.0

    ## x = [r, α_det, θ_det].
    fast(x) = expected_exports_delay(x[1], p, q, T, Gamma(x[2], x[3]))

    x = [0.05, 2.0, 4.5]
    cache = Mooncake.prepare_gradient_cache(fast, x)
    gf = Mooncake.value_and_gradient!!(cache, fast, x)[2][2]
    gd = grad(central_fdm(5, 1), fast, x)[1]
    @test all(isfinite, gf)
    @test gf ≈ gd rtol = 1e-4

    ## The composer path differentiates through combined_delay: x =
    ## [r, α_inc, θ_inc, α_rep, θ_rep].
    via(x) = expected_exports_delay(x[1], p, q, T,
        combined_delay(Gamma(x[2], x[3]), Gamma(x[4], x[5])))
    xv = [0.05, 1.1, 5.7, 2.5, 4.5]
    cache_v = Mooncake.prepare_gradient_cache(via, xv)
    gfv = Mooncake.value_and_gradient!!(cache_v, via, xv)[2][2]
    gdv = grad(central_fdm(5, 1), via, xv)[1]
    @test all(isfinite, gfv)
    @test gfv ≈ gdv rtol = 1e-4
end

@testitem "expected_exports_deaths_delay has correct gradients" tags=[:ad] begin
    using Distributions: Gamma
    using FiniteDifferences: central_fdm, grad
    using Mooncake: Mooncake
    using BVDOutbreakSize: expected_exports_deaths_delay, combined_delay
    CFR = 0.30
    p = 0.25
    q = 1871 / 4_392_200
    T = 90.0

    ## x = [r, α_det, θ_det, α_death, θ_death]; cumulative closes over r.
    fast(x) = expected_exports_deaths_delay(
        s -> exp(x[1] * s), Gamma(x[2], x[3]), Gamma(x[4], x[5]),
        CFR, p, q, T)

    x = [0.05, 2.0, 4.5, 4.3, 2.6]
    cache = Mooncake.prepare_gradient_cache(fast, x)
    gf = Mooncake.value_and_gradient!!(cache, fast, x)[2][2]
    gd = grad(central_fdm(5, 1), fast, x)[1]
    @test all(isfinite, gf)
    @test gf ≈ gd rtol = 1e-4

    ## The composer path differentiates through combined_delay on BOTH
    ## the detection (incubation ⊕ onset-to-report) and death (incubation
    ## ⊕ onset-to-death) delays: x =
    ## [r, α_inc, θ_inc, α_rep, θ_rep, α_death, θ_death].
    via(x) = expected_exports_deaths_delay(
        s -> exp(x[1] * s),
        combined_delay(Gamma(x[2], x[3]), Gamma(x[4], x[5])),
        combined_delay(Gamma(x[2], x[3]), Gamma(x[6], x[7])),
        CFR, p, q, T)
    xv = [0.05, 1.1, 5.7, 2.5, 4.5, 4.3, 2.6]
    cache_v = Mooncake.prepare_gradient_cache(via, xv)
    gfv = Mooncake.value_and_gradient!!(cache_v, via, xv)[2][2]
    gdv = grad(central_fdm(5, 1), via, xv)[1]
    @test all(isfinite, gfv)
    @test gfv ≈ gdv rtol = 1e-4
end

@testitem "bvd_joint delay defaults sample from the prior" tags=[:slow] begin
    using Turing: sample, Prior
    using BVDOutbreakSize: bvd_joint
    import FlexiChains

    ## Small synthetic data: a single cumulative vintage per stream plus
    ## a short dated export-death series. Defaults use the delay
    ## mechanism, whose infection→detection delay is built from the shared
    ## incubation and DRC onset-to-report draws, so the trace carries no
    ## separate detection delay and no rectangular window `w`.
    model = bvd_joint(2, [120], [500], [0, 0, 1];
        reported_offsets = [0],
        death_offsets = [0])

    chn = sample(model, Prior(), 50;
        chain_type = FlexiChains.VNChain, progress = false)

    _has(name) =
        try
            chn[name]
            true
        catch
            false
        end

    ## The detection delay is assembled from the shared incubation
    ## (α_inc / θ_inc) and report (α_rep / θ_rep) draws; there is no
    ## bespoke onset-to-detection delay (α_det / θ_det) nor a McCabe
    ## window `w`.
    @test _has(:α_rep)
    @test _has(:θ_rep)
    @test _has(:α_inc)
    @test _has(:θ_inc)
    @test !_has(:α_det)
    @test !_has(:θ_det)
    @test !_has(:w)
    infections = vec(Array(chn[:cumulative_infections]))
    cases = vec(Array(chn[:cumulative_cases]))
    @test all(isfinite, infections) && all(>(0), infections)
    @test all(isfinite, cases) && all(>(0), cases)
end

@testitem "bvd_joint exports use the combined infection→detection delay" tags=[:slow] begin
    using Turing: sample, Prior
    using BVDOutbreakSize: bvd_joint, exports_delay_model, combined_delay
    using Distributions: Gamma, mean
    import FlexiChains

    ## Swap in a spy export submodel that records the `f_det` it is
    ## handed (the 4th positional argument of exports_delay_model) on each
    ## draw, then sample the joint model. The export detection delay must
    ## be the incubation ⊕ onset-to-report combined delay, NOT the report
    ## delay alone: its mean equals the incubation mean plus the report
    ## mean from that draw.
    captured = Ref{Any}(nothing)

    spy_exports(args...; kwargs...) = begin
        captured[] = args[4]
        exports_delay_model(args...; kwargs...)
    end

    model = bvd_joint(2, [120], [500];
        reported_offsets = [0],
        death_offsets = [0],
        exports = spy_exports)

    chn = sample(model, Prior(), 1;
        chain_type = FlexiChains.VNChain, progress = false)

    f_det = captured[]
    @test f_det isa Gamma

    ## Reconstruct the expected combined delay from that draw's incubation
    ## and report parameters; the captured f_det matches it exactly.
    α_inc = only(vec(Array(chn[:α_inc])))
    θ_inc = only(vec(Array(chn[:θ_inc])))
    α_rep = only(vec(Array(chn[:α_rep])))
    θ_rep = only(vec(Array(chn[:θ_rep])))
    want = combined_delay(Gamma(α_inc, θ_inc), Gamma(α_rep, θ_rep))
    @test f_det.α ≈ want.α
    @test f_det.θ ≈ want.θ
    ## Combined mean exceeds the report mean alone (incubation is added).
    @test mean(f_det) > mean(Gamma(α_rep, θ_rep))
end

@testitem "window detection-timing survival matches expected_exports" begin
    using Distributions: Gamma
    using Turing: sample, Prior
    using FlexiChains: VNChain
    using BVDOutbreakSize: exports_detection_timing_model, expected_exports

    ## Exercises the window-based detection-timing submodel (the McCabe
    ## revert-path component, off the default joint path): its survival
    ## term must equal the at-risk person-time `expected_exports` at the
    ## earlier elapsed time `t1 = T - delta`.
    r = 0.05
    growth_state = (; cumulative = s -> exp(r * s), T = 90.0, r = r)
    p_uganda = 0.25
    daily_travellers = 1871.0
    source_population = 4_392_200.0
    window = 15.0
    delta = 10.0

    model = exports_detection_timing_model(growth_state, p_uganda;
        delta = delta, window = window,
        daily_travellers = daily_travellers,
        source_population = source_population)
    chn = sample(model, Prior(), 1; chain_type = VNChain, progress = false)

    got = only(vec(Array(chn[:survived_exports])))
    q = daily_travellers / source_population
    want = expected_exports(growth_state.cumulative, p_uganda, q,
        growth_state.T - delta, window; r = r)
    @test got ≈ want rtol = 1e-8
    @test got > 0
end
