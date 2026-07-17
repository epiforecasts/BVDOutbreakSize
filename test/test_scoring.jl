## Tests for the forecast-scoring primitives in src/scoring.jl:
## crps_sample, log_crps_sample and score_draws. Scored against
## ScoringRules directly where sensible. Also covers
## select_daily_releases, which picks the releases scripts/score_releases.jl
## scores.

@testitem "select_daily_releases keeps tagged and main-build releases" begin
    using Dates: DateTime
    using BVDOutbreakSize: select_daily_releases

    entries = [("results-v1.9.0", DateTime(2026, 7, 11, 10, 24, 23)),
        ("results-1223", DateTime(2026, 7, 14, 12, 16, 30)),
        ("v1.9.0", DateTime(2026, 7, 11, 10, 0, 0)),
        ("some-other-tag", DateTime(2026, 7, 12, 9, 0, 0))]

    @test select_daily_releases(entries) ==
          ["results-1223", "results-v1.9.0"]
end

@testitem "select_daily_releases keeps the newest main build of a day" begin
    using Dates: DateTime
    using BVDOutbreakSize: select_daily_releases

    ## Two main builds on 2026-07-06 collapse to the later one; the build
    ## on the next day is kept alongside it.
    entries = [("results-1160", DateTime(2026, 7, 6, 13, 36, 42)),
        ("results-1169", DateTime(2026, 7, 6, 22, 2, 50)),
        ("results-1172", DateTime(2026, 7, 7, 8, 0, 46))]

    @test select_daily_releases(entries) == ["results-1172", "results-1169"]
end

@testitem "select_daily_releases prefers the tagged release of a day" begin
    using Dates: DateTime
    using BVDOutbreakSize: select_daily_releases

    ## A tag build and a main build of the same commit publish identical
    ## forecasts under one timestamp; scoring both double-counts them.
    entries = [("results-1204", DateTime(2026, 7, 11, 10, 24, 23)),
        ("results-v1.9.0", DateTime(2026, 7, 11, 10, 24, 23))]
    @test select_daily_releases(entries) == ["results-v1.9.0"]

    ## The tag wins even when a later main build lands the same day.
    later = [("results-v1.9.0", DateTime(2026, 7, 11, 10, 24, 23)),
        ("results-1210", DateTime(2026, 7, 11, 23, 0, 0))]
    @test select_daily_releases(later) == ["results-v1.9.0"]
end

@testitem "select_daily_releases breaks timestamp ties deterministically" begin
    using Dates: DateTime
    using BVDOutbreakSize: select_daily_releases

    ## Two version tags and a main build share one timestamp: the higher
    ## version wins, by version order rather than string order.
    entries = [("results-v1.3.0", DateTime(2026, 6, 9, 22, 58, 13)),
        ("results-706", DateTime(2026, 6, 9, 22, 58, 13)),
        ("results-v1.4.0", DateTime(2026, 6, 9, 22, 58, 13))]
    @test select_daily_releases(entries) == ["results-v1.4.0"]

    ## Main builds sharing a timestamp fall back to the higher run number.
    mains = [("results-43", DateTime(2026, 5, 20, 9, 6, 47)),
        ("results-4", DateTime(2026, 5, 20, 9, 6, 47))]
    @test select_daily_releases(mains) == ["results-43"]
end

@testitem "select_daily_releases returns no tags for no releases" begin
    using Dates: DateTime
    using BVDOutbreakSize: select_daily_releases

    @test select_daily_releases(Tuple{String, DateTime}[]) == String[]
    @test select_daily_releases([("v1.0.0", DateTime(2026, 5, 1))]) ==
          String[]
end

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
