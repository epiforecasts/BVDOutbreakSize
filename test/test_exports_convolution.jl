## Tests for the explicit onset-to-detection delay export expectations
## (`expected_exports_delay`, `expected_exports_deaths_delay`) and the
## delay-mechanism `bvd_joint` defaults. The delay forms generalise the
## McCabe rectangular detection window and reduce to it exactly as the
## onset-to-detection delay collapses to a point mass.

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

    ## Independent rebuild of p·q·os·(∫₀ᵀ C − ∫₀ᵀ conv).
    conv_integral = integrate(
        s -> delay_convolution(1.0, r, s, f_det), 0.0, T)
    want = p * q * (_exp_cumulative_integral(r, 0.0, T) - conv_integral)

    @test got ≈ want rtol = 1e-8
    @test got > 0
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

@testitem "expected_exports_delay scales linearly with onset_fraction" begin
    using Distributions: Gamma
    using BVDOutbreakSize: expected_exports_delay
    r = 0.05
    p = 0.25
    q = 1871 / 4_392_200
    T = 90.0
    f_det = Gamma(2.0, 4.5)

    base = expected_exports_delay(r, p, q, T, f_det; onset_fraction = 1.0)
    half = expected_exports_delay(r, p, q, T, f_det; onset_fraction = 0.5)
    @test half ≈ 0.5 * base rtol = 1e-8
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

    ## As the onset-to-detection Gamma collapses toward a point mass at
    ## `w` the detection survival becomes the top-hat 1{t-s<w}, so the
    ## delay form converges to the window form.
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

## The closed-form Gamma-CDF path must differentiate: gradients of both
## delay expectations w.r.t. the delay (α, θ) and growth r must be finite
## and match central finite differences. Mirrors the Mooncake /
## FiniteDifferences gradient test in test_exports_death_timing.jl.

@testitem "expected_exports_delay has correct gradients" tags=[:ad] begin
    using Distributions: Gamma
    using FiniteDifferences: central_fdm, grad
    using Mooncake: Mooncake
    using BVDOutbreakSize: expected_exports_delay
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
end

@testitem "expected_exports_deaths_delay has correct gradients" tags=[:ad] begin
    using Distributions: Gamma
    using FiniteDifferences: central_fdm, grad
    using Mooncake: Mooncake
    using BVDOutbreakSize: expected_exports_deaths_delay
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
end

@testitem "bvd_joint delay defaults sample from the prior" tags=[:slow] begin
    using Turing: sample, Prior
    using BVDOutbreakSize: bvd_joint
    import FlexiChains

    ## Small synthetic data: a single cumulative vintage per stream plus
    ## a short dated export-death series. Defaults use the delay
    ## mechanism, which reuses the DRC onset-to-report delay
    ## (α_rep, θ_rep) as the onset-to-detection delay, so the trace
    ## carries no separate detection delay and no rectangular window `w`.
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

    ## The detection delay is the shared report delay: α_rep / θ_rep are
    ## present, and there is no bespoke onset-to-detection delay
    ## (α_det / θ_det) nor a McCabe window `w`.
    @test _has(:α_rep)
    @test _has(:θ_rep)
    @test !_has(:α_det)
    @test !_has(:θ_det)
    @test !_has(:w)
    infections = vec(Array(chn[:cumulative_infections]))
    cases = vec(Array(chn[:cumulative_cases]))
    @test all(isfinite, infections) && all(>(0), infections)
    @test all(isfinite, cases) && all(>(0), cases)
end

@testitem "bvd_joint exports reuse the reported onset-to-report delay" tags=[:slow] begin
    using Turing: sample, Prior
    using BVDOutbreakSize: bvd_joint, exports_delay_model
    using Distributions: Gamma
    import FlexiChains

    ## Swap in a spy export submodel that records the `f_det` it is
    ## handed (the 4th positional argument of exports_delay_model) on each
    ## draw, then sample the joint model. The export detection delay must
    ## be the reported stream's onset-to-report delay Gamma(α_rep, θ_rep),
    ## whose parameters are exposed in the trace.
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

    ## The captured detection delay is the last draw's onset-to-report
    ## Gamma; it matches that draw's α_rep / θ_rep.
    α_rep = only(vec(Array(chn[:α_rep])))
    θ_rep = only(vec(Array(chn[:θ_rep])))
    @test f_det.α ≈ α_rep
    @test f_det.θ ≈ θ_rep
end
