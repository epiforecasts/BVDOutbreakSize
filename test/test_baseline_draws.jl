## Tests for the persistence baseline in scripts/score_releases.jl
## (baseline_draws, _history_diffs): the Forecast Hub style spread taken
## from a stream's own past first differences, the widening with horizon,
## zero truncation, the short-history fallback to a plain Poisson draw,
## and reproducibility under a fixed seed.

@testitem "_history_diffs collects vintage gaps up to made_date" begin
    using Dates: Date, Day

    include(joinpath(@__DIR__, "..", "scripts", "score_releases.jl"))

    ## Four weekly vintages, ten days apart in the grid; the fourth is
    ## after made_date and must not contribute a difference.
    hist = (; days = [10, 17, 24, 31], counts = [100.0, 110.0, 125.0, 145.0])
    grid_date(day) = Date(2026, 1, 1) + Day(day)
    made_date = grid_date(24)

    diffs = _history_diffs(hist, grid_date, made_date)
    @test diffs == [(10.0, 7), (15.0, 7)]
end

@testitem "_history_diffs needs at least two vintages" begin
    using Dates: Date, Day

    include(joinpath(@__DIR__, "..", "scripts", "score_releases.jl"))

    grid_date(day) = Date(2026, 1, 1) + Day(day)
    @test _history_diffs(
        (; days = Int[], counts = Float64[]), grid_date, Date(2026, 1, 1)) ==
          Tuple{Float64, Int}[]
    @test _history_diffs(
        (; days = [5], counts = [10.0]), grid_date, Date(2026, 2, 1)) ==
          Tuple{Float64, Int}[]
end

@testitem "baseline_draws keeps the centre unchanged" begin
    using Dates: Date, Day
    using Random: MersenneTwister
    using Statistics: mean

    include(joinpath(@__DIR__, "..", "scripts", "score_releases.jl"))

    n = 60
    cutoff = Date(2026, 7, 15)
    grid_date(day) = cutoff - Day(n - day)
    ## Weekly vintages with a steady rise, giving four same-length (7-day)
    ## past differences by the made date.
    hist = (; days = [10, 17, 24, 31, 38],
        counts = [100.0, 110.0, 125.0, 145.0, 170.0])
    obs = (; confirmed_history = hist)
    made_date = grid_date(38)
    horizon = 7

    rng = MersenneTwister(1)
    draws = baseline_draws(
        obs, grid_date, "confirmed cases", made_date, horizon, 4_000, rng)
    ## The centre is the observed increment over the horizon-length window
    ## ending at made_date: 170 - 145 = 25. The spread pool is symmetric
    ## about zero, so the mean of many draws stays close to the centre.
    @test isapprox(mean(draws), 25.0; atol = 1.0)
end

@testitem "baseline_draws widens with horizon" begin
    using Dates: Date, Day
    using Random: MersenneTwister
    using Statistics: std

    include(joinpath(@__DIR__, "..", "scripts", "score_releases.jl"))

    n = 60
    cutoff = Date(2026, 7, 15)
    grid_date(day) = cutoff - Day(n - day)
    hist = (; days = [10, 17, 24, 31, 38],
        counts = [100.0, 110.0, 125.0, 145.0, 170.0])
    obs = (; confirmed_history = hist)
    made_date = grid_date(38)

    ## The past differences (their own window) do not depend on horizon, so
    ## comparing spreads at two horizons isolates the sqrt(horizon/window)
    ## scaling: a 28-day horizon is scaled twice as far as a 7-day one.
    short = baseline_draws(
        obs, grid_date, "confirmed cases", made_date, 7, 4_000,
        MersenneTwister(2))
    long = baseline_draws(
        obs, grid_date, "confirmed cases", made_date, 28, 4_000,
        MersenneTwister(2))
    @test std(long) > std(short)
    @test isapprox(std(long) / std(short), sqrt(28 / 7); atol = 0.15)
end

@testitem "baseline_draws truncates at zero" begin
    using Dates: Date, Day
    using Random: MersenneTwister

    include(joinpath(@__DIR__, "..", "scripts", "score_releases.jl"))

    n = 60
    cutoff = Date(2026, 7, 15)
    grid_date(day) = cutoff - Day(n - day)
    ## A small, falling level stream: the centre is low and the past
    ## differences include large negative swings, so an untruncated draw
    ## would often fall below zero.
    hist = (; days = [10, 17, 24, 31, 38],
        counts = [40.0, 5.0, 30.0, 2.0, 3.0])
    obs = (; isolation_history = hist)
    made_date = grid_date(38)

    draws = baseline_draws(
        obs, grid_date, "isolation beds", made_date, 7, 2_000,
        MersenneTwister(3))
    @test minimum(draws) >= 0.0
end

@testitem "baseline_draws falls back to Poisson under three differences" begin
    using Dates: Date, Day
    using Random: MersenneTwister
    using Distributions: Poisson

    include(joinpath(@__DIR__, "..", "scripts", "score_releases.jl"))

    n = 40
    cutoff = Date(2026, 7, 15)
    grid_date(day) = cutoff - Day(n - day)
    ## Two vintages give exactly one past difference, short of the
    ## three-difference floor, so the plain Poisson draw is used instead.
    hist = (; days = [10, 17], counts = [100.0, 110.0])
    obs = (; confirmed_history = hist)
    made_date = grid_date(17)
    horizon = 7

    centre = 10.0  # 110 - 100
    seed = 7
    draws = baseline_draws(
        obs, grid_date, "confirmed cases", made_date, horizon, 20,
        MersenneTwister(seed))
    expected = Float64.(rand(MersenneTwister(seed), Poisson(centre), 20))
    @test draws == expected
end

@testitem "baseline_draws is reproducible under a fixed seed" begin
    using Dates: Date, Day
    using Random: MersenneTwister

    include(joinpath(@__DIR__, "..", "scripts", "score_releases.jl"))

    n = 60
    cutoff = Date(2026, 7, 15)
    grid_date(day) = cutoff - Day(n - day)
    hist = (; days = [10, 17, 24, 31, 38],
        counts = [100.0, 110.0, 125.0, 145.0, 170.0])
    obs = (; confirmed_history = hist)
    made_date = grid_date(38)

    a = baseline_draws(
        obs, grid_date, "confirmed cases", made_date, 14, 500,
        MersenneTwister(42))
    b = baseline_draws(
        obs, grid_date, "confirmed cases", made_date, 14, 500,
        MersenneTwister(42))
    @test a == b
end
