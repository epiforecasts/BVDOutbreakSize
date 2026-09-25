@testitem "generation_interval_model: GI mean prior matches the NEJM interval" begin
    using BVDOutbreakSize: generation_interval_model, cdf_nmax
    using Distributions: Gamma
    using Random: seed!
    using Statistics: quantile

    ## The source's 95% CI on the serial-interval mean is 13.0-17.6 d.
    seed!(1)
    m = generation_interval_model(cdf_nmax(Gamma(2.71, 5.65)))
    gi_mean = [m().gi_mean for _ in 1:4000]
    lo, hi = quantile(gi_mean, [0.025, 0.975])
    @test isapprox(lo, 13.0; atol = 0.3)
    @test isapprox(hi, 17.6; atol = 0.3)
end
