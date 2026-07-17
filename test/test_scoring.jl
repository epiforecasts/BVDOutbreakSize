## Tests for the forecast-scoring primitives in src/scoring.jl:
## crps_sample, log_crps_sample and score_draws. Scored against
## ScoringRules directly where sensible. Also covers
## select_daily_releases, which picks the releases scripts/score_releases.jl
## scores.

@testitem "select_daily_releases keeps tagged and main-build releases" begin
    using Dates: Date, DateTime
    using BVDOutbreakSize: select_daily_releases

    entries = [
        ("results-v1.9.0", DateTime(2026, 7, 11, 10, 24, 23),
            Date(2026, 7, 8)),
        ("results-1223", DateTime(2026, 7, 14, 12, 16, 30),
            Date(2026, 7, 11)),
        ("v1.9.0", DateTime(2026, 7, 11, 10, 0, 0), Date(2026, 7, 8)),
        ("some-other-tag", DateTime(2026, 7, 12, 9, 0, 0),
            Date(2026, 7, 9))]

    @test select_daily_releases(entries) ==
          ["results-1223", "results-v1.9.0"]
end

@testitem "select_daily_releases keeps the newest build of a data day" begin
    using Dates: Date, DateTime
    using BVDOutbreakSize: select_daily_releases

    ## Two main builds sharing cut-off 2026-07-03 collapse to the later
    ## one; the build on the next data day is kept alongside it.
    entries = [
        ("results-1160", DateTime(2026, 7, 6, 13, 36, 42),
            Date(2026, 7, 3)),
        ("results-1169", DateTime(2026, 7, 6, 22, 2, 50),
            Date(2026, 7, 3)),
        ("results-1172", DateTime(2026, 7, 7, 8, 0, 46),
            Date(2026, 7, 4))]

    @test select_daily_releases(entries) == ["results-1172", "results-1169"]
end

@testitem "select_daily_releases drops same-cutoff republications" begin
    using Dates: Date, DateTime
    using BVDOutbreakSize: select_daily_releases

    ## The report re-renders whenever main moves, so a main build and a
    ## later tag build can publish the same data days apart. Keying on the
    ## creation day would keep both and score one forecast twice; the
    ## shared cut-off collapses them, and the tag wins.
    entries = [
        ("results-1187", DateTime(2026, 7, 9, 8, 34, 2), Date(2026, 7, 6)),
        ("results-v1.8.0", DateTime(2026, 7, 10, 10, 1, 25),
            Date(2026, 7, 6))]
    @test select_daily_releases(entries) == ["results-v1.8.0"]

    ## Same mechanism between two main builds, with no tag to prefer: the
    ## later build of the same data survives.
    mains = [
        ("results-750", DateTime(2026, 6, 11, 21, 19, 50),
            Date(2026, 6, 10)),
        ("results-762", DateTime(2026, 6, 12, 9, 51, 19),
            Date(2026, 6, 10))]
    @test select_daily_releases(mains) == ["results-762"]
end

@testitem "select_daily_releases prefers the tagged release of a day" begin
    using Dates: Date, DateTime
    using BVDOutbreakSize: select_daily_releases

    ## A tag build and a main build of the same commit publish identical
    ## forecasts under one timestamp; scoring both double-counts them.
    entries = [
        ("results-1204", DateTime(2026, 7, 11, 10, 24, 23),
            Date(2026, 7, 8)),
        ("results-v1.9.0", DateTime(2026, 7, 11, 10, 24, 23),
            Date(2026, 7, 8))]
    @test select_daily_releases(entries) == ["results-v1.9.0"]

    ## The tag wins even when a later main build shares its cut-off.
    later = [
        ("results-v1.9.0", DateTime(2026, 7, 11, 10, 24, 23),
            Date(2026, 7, 8)),
        ("results-1210", DateTime(2026, 7, 11, 23, 0, 0),
            Date(2026, 7, 8))]
    @test select_daily_releases(later) == ["results-v1.9.0"]
end

@testitem "select_daily_releases breaks timestamp ties deterministically" begin
    using Dates: Date, DateTime
    using BVDOutbreakSize: select_daily_releases

    ## Two version tags and a main build share one timestamp: the higher
    ## version wins. v1.10.0 sorts BELOW v1.9.0 as a string and above it as
    ## a version, so this fails if the version is ever compared as text.
    entries = [
        ("results-v1.9.0", DateTime(2026, 6, 9, 22, 58, 13),
            Date(2026, 6, 7)),
        ("results-706", DateTime(2026, 6, 9, 22, 58, 13),
            Date(2026, 6, 7)),
        ("results-v1.10.0", DateTime(2026, 6, 9, 22, 58, 13),
            Date(2026, 6, 7))]
    @test select_daily_releases(entries) == ["results-v1.10.0"]

    ## Main builds sharing a timestamp fall back to the higher run number.
    ## "10" sorts below "9" as a string, so this fails on a text compare.
    mains = [("results-9", DateTime(2026, 5, 20, 9, 6, 47),
            Date(2026, 5, 18)),
        ("results-10", DateTime(2026, 5, 20, 9, 6, 47), Date(2026, 5, 18))]
    @test select_daily_releases(mains) == ["results-10"]
end

@testitem "select_daily_releases ignores non-results tags" begin
    using Dates: Date, DateTime
    using BVDOutbreakSize: is_results_release, select_daily_releases

    ## The backfill release stores reconstructed forecasts as assets and is
    ## never itself a candidate. The repo also publishes a release per code
    ## tag, which carries no results assets.
    @test !is_results_release("forecasts-backfill")
    @test !is_results_release("v1.9.0")
    @test is_results_release("results-v1.9.0")
    @test is_results_release("results-1243")

    entries = [
        ("forecasts-backfill", DateTime(2026, 7, 12, 9, 0, 0),
            Date(2026, 7, 8)),
        ("v1.9.0", DateTime(2026, 7, 11, 10, 0, 0), Date(2026, 7, 8))]
    @test select_daily_releases(entries) == String[]

    ## A tagged release must not be read as a main build, which would rank
    ## it below one and compare its run number as a version.
    tagged = [("results-v1.9.0", DateTime(2026, 7, 11, 10, 24, 23),
        Date(2026, 7, 8))]
    @test select_daily_releases(tagged) == ["results-v1.9.0"]
end

@testitem "select_daily_releases returns no tags for no releases" begin
    using Dates: Date, DateTime
    using BVDOutbreakSize: select_daily_releases

    @test select_daily_releases(Tuple{String, DateTime, Date}[]) == String[]
    @test select_daily_releases([("v1.0.0", DateTime(2026, 5, 1),
        Date(2026, 4, 28))]) == String[]
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
