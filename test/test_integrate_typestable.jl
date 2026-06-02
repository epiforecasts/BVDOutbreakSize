## Type stability and numerical-identity guards for the hand-rolled
## Gauss-Legendre `integrate`. The previous implementation routed the
## integrand closure through a SciML `IntegralProblem`/`solve`, which lost
## inference across the parameter boundary and returned `Any` from `.u`.
## That `Any` propagated into `expected_exports_deaths`, boxing the
## exports-deaths likelihood inside the NUTS inner loop. The fix evaluates
## the same fixed Gauss-Legendre rule directly so Julia specialises on the
## integrand's concrete return type. The numerical test pins agreement
## with the Integrals.jl reference computed live in the test environment.

@testsnippet IntegrateRef begin
    using Distributions: Gamma, pdf
    import Integrals

    f1(x) = exp(0.013 * x) + 0.5
    f2(x) = sin(0.1 * x) + 2.0
    f3(x) = pdf(Gamma(4.3, 2.6), x)

    ## Reference oracle: the previous Integrals.jl implementation of
    ## `integrate` — the same fixed Gauss-Legendre rule (n = 32 uniform,
    ## n = 64 clustered) routed through `IntegralProblem` / `solve`,
    ## computed live here so the hand-rolled rule is checked against the
    ## actual library, not pinned constants. Integrals.jl is a test-only
    ## dependency; it is dropped from the package itself.
    function ref_uniform(f, lo, hi; n = 32)
        hi <= lo && return zero(hi - lo)
        hw = (hi - lo) / 2
        prob = Integrals.IntegralProblem(
            (u, p) -> p.f(p.hw * (u + 1) + p.lo), (-1.0, 1.0),
            (; f, hw, lo))
        return hw * Integrals.solve(prob, Integrals.GaussLegendre(; n)).u
    end
    function ref_clustered(f, lo, hi, scale; n = 64)
        hi <= lo && return zero(hi - lo)
        (isfinite(scale) && scale > zero(scale)) ||
            return ref_uniform(f, lo, hi; n)
        span = hi - lo
        expo = max(one(span), log(span / scale) / log(oftype(span, 2)))
        prob = Integrals.IntegralProblem(
            (u, p) -> begin
                v = (u + one(u)) / 2
                d = p.span * v^p.expo
                p.f(p.hi - d) *
                (p.span * p.expo * v^(p.expo - one(p.expo)) / 2)
            end, (-1.0, 1.0), (; f, hi, span, expo))
        return Integrals.solve(prob, Integrals.GaussLegendre(; n)).u
    end

    UNIFORM_CASES = [(f1, 0.0, 10.0), (f1, 0.0, 1.0), (f2, 2.0, 50.0),
        (f3, 0.0, 40.0), (f1, 5.0, 5.0), (f1, 10.0, 5.0)]
    CLUSTERED_CASES = [(f1, 0.0, 10.0, 2.0), (f2, 2.0, 50.0, 5.0),
        (f3, 0.0, 360.0, 0.5), (f1, 0.0, 40.0, 1.0e6),
        (f1, 0.0, 10.0, 0.0), (f1, 0.0, 10.0, -1.0),
        (f1, 0.0, 10.0, Inf), (f1, 5.0, 5.0, 2.0)]
end

@testitem "integrate: matches Integrals.jl reference" setup=[IntegrateRef] begin
    using BVDOutbreakSize: integrate
    for (f, lo, hi) in UNIFORM_CASES
        @test isapprox(integrate(f, lo, hi), ref_uniform(f, lo, hi);
            rtol = 1e-10, atol = 1e-12)
    end
    for (f, lo, hi, sc) in CLUSTERED_CASES
        @test isapprox(integrate(f, lo, hi, sc),
            ref_clustered(f, lo, hi, sc); rtol = 1e-10, atol = 1e-12)
    end
end

@testitem "integrate: type stable" setup=[IntegrateRef] begin
    using BVDOutbreakSize: integrate
    using Test: @inferred
    @test (@inferred integrate(f1, 0.0, 10.0)) isa Float64
    @test (@inferred integrate(f1, 0.0, 10.0, 2.0)) isa Float64
    ## degenerate-scale fallback path stays stable too
    @test (@inferred integrate(f1, 0.0, 10.0, 0.0)) isa Float64
end

@testitem "expected_exports_deaths: type stable" begin
    using BVDOutbreakSize: expected_exports_deaths, ExportDeathDelay
    using Distributions: Gamma
    using Test: @inferred
    C(s) = exp(0.05 * s)
    ed = ExportDeathDelay(Gamma(4.3, 2.6), 30.0)
    val = @inferred expected_exports_deaths(
        C, ed, 0.3, 0.6, 0.001, 50.0, 10.0)
    @test val isa Float64
end

@testitem "integrate: allocation-free" setup=[IntegrateRef] begin
    using BVDOutbreakSize: integrate, expected_exports_deaths,
                           ExportDeathDelay
    using Distributions: Gamma
    integrate(f1, 0.0, 10.0)
    @test (@allocated integrate(f1, 0.0, 10.0)) == 0
    @test (@allocated integrate(f1, 0.0, 10.0, 2.0)) == 0
    ## `integrate` itself is allocation-free; the small residual on the
    ## full `expected_exports_deaths` is the integrand closure capturing
    ## the `ExportDeathDelay` grid, which predates this fix (was ~7 kB,
    ## dominated by the boxed `integrate` result, now a few hundred bytes).
    C(s) = exp(0.05 * s)
    ed = ExportDeathDelay(Gamma(4.3, 2.6), 30.0)
    expected_exports_deaths(C, ed, 0.3, 0.6, 0.001, 50.0, 10.0)
    @test (@allocated expected_exports_deaths(
        C, ed, 0.3, 0.6, 0.001, 50.0, 10.0)) < 1024
end
