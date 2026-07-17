## Tests for the forecast-scoring primitives in src/scoring.jl:
## crps_sample, log_crps_sample and score_draws. Scored against
## ScoringRules directly where sensible.

@testitem "crps_sample matches ScoringRules.crps(samples, obs)" begin
    using ScoringRules: crps
    using BVDOutbreakSize: crps_sample

    samples = [1.0, 2.0, 3.0, 4.0, 5.0]
    @test crps_sample(2.5, samples) == crps(samples, 2.5)
end

@testitem "crps_sample of a point-mass ensemble is the absolute error" begin
    using BVDOutbreakSize: crps_sample

    @test crps_sample(5.0, fill(3.0, 100)) ≈ 2.0
    @test crps_sample(3.0, fill(5.0, 100)) ≈ 2.0
end

@testitem "crps_sample is zero at a point-mass forecast equal to obs" begin
    using BVDOutbreakSize: crps_sample

    @test crps_sample(5.0, fill(5.0, 100)) ≈ 0.0
end

@testitem "log_crps_sample scores on the log scale, not log of the score" begin
    using BVDOutbreakSize: crps_sample, log_crps_sample

    obs = 120.0
    samples = [10.0, 50.0, 100.0, 150.0, 400.0, 900.0]
    @test log_crps_sample(obs, samples) ==
          crps_sample(log1p(obs), log1p.(samples))
    @test log_crps_sample(obs, samples) != log(crps_sample(obs, samples))
end

@testitem "score_draws returns the documented NamedTuple" begin
    using Random: MersenneTwister
    using BVDOutbreakSize: score_draws, crps_sample, log_crps_sample,
                           bias_sample

    rng = MersenneTwister(1)
    samples = 100.0 .+ 10.0 .* randn(rng, 500)
    obs = 105.0

    s = score_draws(obs, samples)
    @test s.crps ≈ crps_sample(obs, samples)
    @test s.log_crps ≈ log_crps_sample(obs, samples)
    @test s.bias ≈ bias_sample(obs, samples)
    @test s.n == length(samples)
    @test s.coverage_50 isa Bool
    @test s.coverage_90 isa Bool
    @test -1 <= s.bias <= 1
end

@testitem "score_draws coverage flags respond to how far obs sits" begin
    using Random: MersenneTwister
    using BVDOutbreakSize: score_draws

    rng = MersenneTwister(2)
    samples = 100.0 .+ 10.0 .* randn(rng, 2_000)

    inside = score_draws(102.0, samples)
    @test inside.coverage_90 == true

    outside = score_draws(1000.0, samples)
    @test outside.coverage_90 == false
end

@testitem "score_draws handles empty samples gracefully" begin
    using BVDOutbreakSize: score_draws

    s = score_draws(5.0, Float64[])
    @test isnan(s.crps)
    @test isnan(s.log_crps)
    @test isnan(s.bias)
    @test s.coverage_50 == false
    @test s.coverage_90 == false
    @test s.n == 0
end
