## Tests for onsets_over_time: the per-date posterior summary of the latent
## symptom-onset trajectory used to share the "symptomatic cases" curve.

@testitem "onsets_over_time returns one dated row per grid day" begin
    using DataFrames: DataFrame, nrow
    using Dates: Date, Day
    using Random: MersenneTwister
    import FlexiChains
    using BVDOutbreakSize: onsets_over_time
    rng = MersenneTwister(7)
    ndraws = 80
    n = 30
    ## Each cumulative-onset deterministic is a draws×chains matrix of
    ## per-draw monotone vectors.
    traj = reshape(
        [cumsum(abs.(randn(rng, n))) for _ in 1:ndraws], ndraws, 1)
    chn = FlexiChains.FlexiChain{Symbol}(ndraws, 1,
        Dict(FlexiChains.Parameter(:cumulative_onsets) => traj))

    seeding = Date("2026-03-01")
    df = onsets_over_time(chn; n = n, seeding = seeding)

    @test df isa DataFrame
    @test nrow(df) == n
    @test names(df) == [
        "date",
        "new_onsets_lower_90", "new_onsets_lower_60", "new_onsets_lower_30",
        "new_onsets_upper_30", "new_onsets_upper_60", "new_onsets_upper_90",
        "cumulative_onsets_lower_90", "cumulative_onsets_lower_60",
        "cumulative_onsets_lower_30", "cumulative_onsets_upper_30",
        "cumulative_onsets_upper_60", "cumulative_onsets_upper_90"]
    @test df.date[1] == seeding
    @test df.date[end] == seeding + Day(n - 1)
end

@testitem "onsets_over_time intervals are ordered and cumulative grows" begin
    using Dates: Date
    using Random: MersenneTwister
    import FlexiChains
    using BVDOutbreakSize: onsets_over_time
    rng = MersenneTwister(11)
    ndraws = 120
    n = 25
    traj = reshape(
        [cumsum(abs.(randn(rng, n))) for _ in 1:ndraws], ndraws, 1)
    chn = FlexiChains.FlexiChain{Symbol}(ndraws, 1,
        Dict(FlexiChains.Parameter(:cumulative_onsets) => traj))

    df = onsets_over_time(chn; n = n, seeding = Date("2026-03-01"))

    ## Nested 30/60/90 endpoints are ordered, every day, both quantities.
    for q in ("new_onsets", "cumulative_onsets")
        @test all(df[!, "$(q)_lower_90"] .<= df[!, "$(q)_lower_60"] .<=
                  df[!, "$(q)_lower_30"] .<= df[!, "$(q)_upper_30"] .<=
                  df[!, "$(q)_upper_60"] .<= df[!, "$(q)_upper_90"])
    end
    ## The cumulative onsets are non-decreasing over time (monotone per draw).
    @test issorted(df.cumulative_onsets_upper_30)
    ## Daily new onsets are non-negative.
    @test all(df.new_onsets_lower_90 .>= 0)
end
