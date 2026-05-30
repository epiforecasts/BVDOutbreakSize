## Type stability and numerical-identity guards for the hand-rolled
## Gauss-Legendre `integrate`. The previous implementation routed the
## integrand closure through a SciML `IntegralProblem`/`solve`, which lost
## inference across the parameter boundary and returned `Any` from `.u`.
## That `Any` propagated into `expected_exports_deaths`, boxing the
## exports-deaths likelihood inside the NUTS inner loop. The fix evaluates
## the same fixed Gauss-Legendre rule directly so Julia specialises on the
## integrand's concrete return type. These tests pin both the restored
## type stability and bit-for-bit-close agreement with the reference rule.

@testsnippet IntegrateGolden begin
    using Distributions: Gamma, pdf

    f1(x) = exp(0.013 * x) + 0.5
    f2(x) = sin(0.1 * x) + 2.0
    f3(x) = pdf(Gamma(4.3, 2.6), x)

    ## Reference values captured from the pre-fix `solve`-based
    ## implementation; the hand-rolled rule must reproduce them.
    UNIFORM_GOLDEN = [
        ((f1, 0.0, 10.0), 15.679106409586298),
        ((f1, 0.0, 1.0), 1.5065282584468584),
        ((f2, 2.0, 50.0), 102.96404392378018),
        ((f3, 0.0, 40.0), 0.9997571129065723),
        ((f1, 5.0, 5.0), 0.0),
        ((f1, 10.0, 5.0), 0.0)
    ]
    CLUSTERED_GOLDEN = [
        ((f1, 0.0, 10.0, 2.0), 15.679106399481002),
        ((f2, 2.0, 50.0, 5.0), 102.96404392380417),
        ((f3, 0.0, 360.0, 0.5), 1.0021325794366642),
        ((f1, 0.0, 40.0, 1.0e6), 72.46366536145278),
        ((f1, 0.0, 10.0, 0.0), 15.679106409586293),
        ((f1, 0.0, 10.0, -1.0), 15.679106409586293),
        ((f1, 0.0, 10.0, Inf), 15.679106409586293),
        ((f1, 5.0, 5.0, 2.0), 0.0)
    ]
end

@testitem "integrate: identical to reference" setup=[IntegrateGolden] begin
    using BVDOutbreakSize: integrate
    for ((f, lo, hi), want) in UNIFORM_GOLDEN
        got = integrate(f, lo, hi)
        @test isapprox(got, want; rtol = 1e-10, atol = 0.0)
    end
    for ((f, lo, hi, sc), want) in CLUSTERED_GOLDEN
        got = integrate(f, lo, hi, sc)
        @test isapprox(got, want; rtol = 1e-10, atol = 0.0)
    end
end

@testitem "integrate: type stable" setup=[IntegrateGolden] begin
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

@testitem "integrate: allocation-free" setup=[IntegrateGolden] begin
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
