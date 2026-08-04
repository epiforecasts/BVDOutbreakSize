## Tests for the persistence baseline in scripts/score_releases.jl
## (baseline_draws, _history_diffs, vintage_observations): the Forecast Hub
## style spread taken from a stream's own past first differences via an
## iterated daily random walk, the widening with horizon, zero truncation,
## the short-history fallback to a plain Poisson draw, reproducibility under
## a fixed seed, and the made_date vintage the baseline is built from.

@testitem "_history_diffs collects vintage gaps up to made_date" begin
    using Dates: Date, Day

    include(joinpath(@__DIR__, "..", "scripts", "score_releases.jl"))

    ## Four weekly vintages, ten days apart in the grid; the fourth is
    ## after made_date and must not contribute a difference.
    hist = (; days = [10, 17, 24, 31], counts = [100.0, 110.0, 125.0, 145.0])
    grid_date(day) = Date(2026, 1, 1) + Day(day)
    made_date = grid_date(24)
    ## No break-day fields, so `break_correction` is a no-op: this obs
    ## carries no `confirmed_break_days` at all.
    obs = (; confirmed_history = hist)

    diffs = _history_diffs(obs, grid_date, "confirmed cases", hist, made_date)
    @test diffs == [(10.0, 7), (15.0, 7)]
end

@testitem "_history_diffs needs at least two vintages" begin
    using Dates: Date, Day

    include(joinpath(@__DIR__, "..", "scripts", "score_releases.jl"))

    grid_date(day) = Date(2026, 1, 1) + Day(day)
    obs = (;)
    @test _history_diffs(obs, grid_date, "confirmed cases",
        (; days = Int[], counts = Float64[]), Date(2026, 1, 1)) ==
          Tuple{Float64, Int}[]
    @test _history_diffs(obs, grid_date, "confirmed cases",
        (; days = [5], counts = [10.0]), Date(2026, 2, 1)) ==
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

@testitem "baseline_draws' spread is an iterated daily random walk" begin
    using Dates: Date, Day
    using Random: MersenneTwister
    using Statistics: std

    include(joinpath(@__DIR__, "..", "scripts", "score_releases.jl"))

    ## A single 7-day-apart vintage pair gives one past difference of 70
    ## over a 7-day window; three copies clear the three-difference floor
    ## without changing the step pool. The one-day step this history
    ## produces is 70 / sqrt(7), not 70 / 7: a `window`-day change has
    ## variance `window * sigma^2` under a zero-drift walk, so dividing by
    ## `sqrt(window)` is what recovers the one-day step size `sigma`.
    n = 60
    cutoff = Date(2026, 7, 15)
    grid_date(day) = cutoff - Day(n - day)
    hist = (; days = [3, 10, 17, 24, 31, 38],
        counts = [30.0, 100.0, 170.0, 240.0, 310.0, 380.0])
    obs = (; confirmed_history = hist)
    made_date = grid_date(38)

    ## At horizon 1 the walk takes exactly one daily step, so the draws are
    ## the centre plus one draw from the symmetrised per-day-step pool
    ## (±70/sqrt(7), the only step this history produces): every draw lands
    ## on centre - 70/sqrt(7) or centre + 70/sqrt(7) (or is floored at
    ## zero).
    centre = 70.0  # 380 - 310, the observed increment over the last 7 days
    step = 70.0 / sqrt(7.0)
    draws1 = baseline_draws(
        obs, grid_date, "confirmed cases", made_date, 1, 4_000,
        MersenneTwister(11))
    @test all(d -> isapprox(d, centre - step) || isapprox(d, centre + step),
        draws1)

    ## Summing `horizon` iid steps gives a variance proportional to the step
    ## count: the population variance of a ±70/sqrt(7) Rademacher-like step
    ## is `step^2`, so `horizon` steps have variance `step^2 * horizon`,
    ## i.e. std growing with sqrt(horizon) — the same magnitude a single
    ## draw rescaled by `70 * sqrt(horizon / 7)` would give (the
    ## historically-correct target, since a `window`-day history diff
    ## rescaled to a `horizon`-day one is `d * sqrt(horizon / window)`).
    horizon = 16
    draws16 = baseline_draws(
        obs, grid_date, "confirmed cases", made_date, horizon, 20_000,
        MersenneTwister(12))
    expected_sd = 70.0 * sqrt(horizon / 7.0)
    @test isapprox(std(draws16), expected_sd; rtol = 0.1)
end

@testitem "vintage_observations falls back to the given obs with no path" begin
    using Dates: Date, Day

    include(joinpath(@__DIR__, "..", "scripts", "score_releases.jl"))

    obs = (; confirmed_history = (; days = [1], counts = [1.0]))
    grid_date(day) = Date(2026, 1, 1) + Day(day)

    ov, ogd = vintage_observations(nothing, Date(2026, 1, 5), obs, grid_date)
    @test ov === obs
    @test ogd === grid_date
end

## A minimal `observations.toml` fixture: just enough for `load_observations`
## (`as_of_date`, `genetic_tmrca.date`, one history block) so
## `vintage_observations` can load and freeze it without pulling in a real
## release manifest. A `@testsnippet`, not a plain top-level function: each
## `@testitem` below only sees its own code range, not sibling definitions
## in the same file, so the fixture must be shared through the snippet
## mechanism to be visible in more than one test item.
@testsnippet FixtureObs begin
    function _write_fixture_manifest(path; as_of_date, confirmed_dates,
            confirmed_values)
        dates_toml = join(("\"$d\"" for d in confirmed_dates), ", ")
        values_toml = join(confirmed_values, ", ")
        write(path, """
            as_of_date = "$as_of_date"

            [genetic_tmrca]
            date = "2026-03-01"

            [confirmed_case_history]
            dates = [$dates_toml]
            values = [$values_toml]
            """)
        return path
    end
end

@testitem "vintage_observations freezes the manifest" setup=[FixtureObs] begin
    using Dates: Date

    include(joinpath(@__DIR__, "..", "scripts", "score_releases.jl"))

    path = _write_fixture_manifest(joinpath(mktempdir(), "observations.toml");
        as_of_date = "2026-07-20",
        confirmed_dates = ["2026-07-01", "2026-07-08", "2026-07-15"],
        confirmed_values = [50, 80, 100])

    ## Placeholder current obs/grid_date, unused when a path is given.
    placeholder_obs = (;)
    placeholder_grid_date(day) = Date(2026, 1, 1)

    ov, ogd = vintage_observations(
        path, Date(2026, 7, 15), placeholder_obs, placeholder_grid_date)
    ## Frozen to made_date (15 July), not the manifest's own as_of_date
    ## (20 July): the vintage a forecast made on the 15th could have seen.
    @test ov.cutoff == Date(2026, 7, 15)
    ## The manifest's own history, converted back to calendar dates through
    ## `ogd` (not the placeholder), reproduces the vintages as written.
    @test [string(ogd(d)) for d in ov.confirmed_history.days] ==
          ["2026-07-01", "2026-07-08", "2026-07-15"]
    @test ov.confirmed_history.counts == [50, 80, 100]
end

@testitem "vintage_observations shields against a leak" setup=[FixtureObs] begin
    using Dates: Date, Day
    using Random: MersenneTwister
    using Statistics: mean

    include(joinpath(@__DIR__, "..", "scripts", "score_releases.jl"))

    ## The release's own snapshot at made_date: confirmed cases stood at
    ## 100 on 2026-07-15, an increment of 20 over the preceding week.
    made_date = Date(2026, 7, 15)
    release_path = _write_fixture_manifest(
        joinpath(mktempdir(), "observations.toml");
        as_of_date = "2026-07-15",
        confirmed_dates = ["2026-07-01", "2026-07-08", "2026-07-15"],
        confirmed_values = [50, 80, 100])

    ## The current manifest, as it stands today: INSP later revised the
    ## 15 July vintage upward to 150 (backfilled cases attributed to that
    ## report after the fact). A baseline built from this directly would
    ## see an increment of 70 the forecast itself never had.
    current_grid_date(day) = Date(2026, 7, 1) + Day(day)
    current_hist = (; days = [0, 7, 14], counts = [50.0, 80.0, 150.0])
    current_obs = (; confirmed_history = current_hist)

    vobs, vgrid_date = vintage_observations(
        release_path, made_date, current_obs, current_grid_date)

    rng = MersenneTwister(21)
    vintage_draws = baseline_draws(
        vobs, vgrid_date, "confirmed cases", made_date, 7, 2_000, rng)
    ## The vintage-shielded baseline centres on the pre-revision increment
    ## (20), not the leaked post-revision one (70).
    @test isapprox(mean(vintage_draws), 20.0; atol = 3.0)
    @test !isapprox(mean(vintage_draws), 70.0; atol = 5.0)

    ## Using the current manifest directly (the pre-fix behaviour) would
    ## instead see the leaked revision.
    leaked_draws = baseline_draws(
        current_obs, current_grid_date, "confirmed cases", made_date, 7,
        2_000, MersenneTwister(21))
    @test isapprox(mean(leaked_draws), 70.0; atol = 3.0)
end
