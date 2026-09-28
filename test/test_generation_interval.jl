@testitem "generation_interval_model: GI priors match the NEJM sampling error" begin
    using BVDOutbreakSize: generation_interval_model, cdf_nmax
    using Distributions: Gamma
    using Random: seed!
    using Statistics: cor, quantile

    ## The source's fitted Gamma is mean 15.3 d, SD 9.3 d from 92 pairs. The
    ## 95% prior intervals are the sampling error of each from 92 pairs.
    seed!(1)
    m = generation_interval_model(cdf_nmax(Gamma(2.71, 5.65)))
    draws = [m() for _ in 1:4000]
    gi_mean = getfield.(draws, :gi_mean)
    gi_sd = getfield.(draws, :gi_sd)
    @test all(isapprox(sum(d.g), 1) for d in draws[1:10])
    lo, hi = quantile(gi_mean, [0.025, 0.975])
    @test isapprox(lo, 13.4; atol = 0.2)
    @test isapprox(hi, 17.2; atol = 0.2)
    lo, hi = quantile(gi_sd, [0.025, 0.975])
    @test isapprox(lo, 7.3; atol = 0.2)
    @test isapprox(hi, 11.3; atol = 0.2)
    ## The mean and SD carry separate sourced widths, so they are drawn
    ## independently.
    @test abs(cor(gi_mean, gi_sd)) < 0.05
end

@testitem "generation_interval_model: an overflowed mean is rejected, not thrown" begin
    ## The step-size search at the start of warm-up can push `gi_mean` to
    ## `Inf`, so the Gamma scale `sd^2 / mean` is exactly zero. The density
    ## is then `-Inf`, which the sampler rejects, rather than a DomainError.
    using BVDOutbreakSize: generation_interval_model, cdf_nmax
    using Distributions: Gamma
    using Turing: logjoint

    m = generation_interval_model(cdf_nmax(Gamma(2.71, 5.65)))
    @test logjoint(m, (gi_mean = Inf, gi_sd = 9.3)) == -Inf
    @test isfinite(logjoint(m, (gi_mean = 15.3, gi_sd = 9.3)))
end
