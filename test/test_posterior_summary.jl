## Tests for posterior_summary: the six equal-tailed quantile fields
## are produced in the documented order and match `Statistics.quantile`
## on a known vector.

@testitem "posterior_summary returns six fields in documented order" begin
    using Statistics: quantile
    using BVDOutbreakSize: posterior_summary
    xs = collect(0.0:0.01:1.0)
    s = posterior_summary(xs)

    # Field names and order match the docstring.
    @test propertynames(s) ==
        (:lo90, :lo60, :lo30, :hi30, :hi60, :hi90)

    # Values match Statistics.quantile on the same vector.
    @test s.lo90 ≈ quantile(xs, 0.05)
    @test s.lo60 ≈ quantile(xs, 0.2)
    @test s.lo30 ≈ quantile(xs, 0.35)
    @test s.hi30 ≈ quantile(xs, 0.65)
    @test s.hi60 ≈ quantile(xs, 0.8)
    @test s.hi90 ≈ quantile(xs, 0.95)

    # Monotone in the credible level.
    @test s.lo90 < s.lo60 < s.lo30 < s.hi30 < s.hi60 < s.hi90
end
