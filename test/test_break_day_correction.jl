## Tests for the harmonisation-break-day correction in
## scripts/score_releases.jl (break_correction, break_day_corrections,
## truth_at, baseline_draws, carry_break_days): a listed confirmed-stream
## break day's retrospective step (net vintage increment minus the printed
## 24h gross count) is subtracted from any window it falls inside, for both
## the forecast truth and the persistence baseline centre, and leaves every
## other stream untouched. See issue #511.

@testitem "truth_at subtracts a break day inside the window" begin
    using Dates: Date, Day

    include(joinpath(@__DIR__, "..", "scripts", "score_releases.jl"))

    grid_date(day) = Date(2026, 1, 1) + Day(day)
    ## Vintage at day 24 steps 110 -> 379, a net of 269 against a printed
    ## 24h gross of 97: a retrospective correction of 172.
    hist = (; days = [10, 17, 24, 31], counts = [100.0, 110.0, 379.0, 400.0])
    obs = (; cutoff = grid_date(40), confirmed_history = hist,
        confirmed_break_days = [24], confirmed_break_gross_cases = [97])

    made_date = grid_date(17)
    target_date = grid_date(31)  # (17, 31] contains the break day at 24
    ## Raw increment 400 - 110 = 290, less the 172 correction = 118.
    @test truth_at(obs, grid_date, "confirmed cases", made_date, target_date) ==
          118.0
end

@testitem "truth_at leaves a window without a break day unchanged" begin
    using Dates: Date, Day

    include(joinpath(@__DIR__, "..", "scripts", "score_releases.jl"))

    grid_date(day) = Date(2026, 1, 1) + Day(day)
    hist = (; days = [10, 17, 24, 31], counts = [100.0, 110.0, 379.0, 400.0])
    obs = (; cutoff = grid_date(40), confirmed_history = hist,
        confirmed_break_days = [24], confirmed_break_gross_cases = [97])

    ## made_date sits ON the break day, so (24, 31] does not contain it.
    made_date = grid_date(24)
    target_date = grid_date(31)
    @test truth_at(obs, grid_date, "confirmed cases", made_date, target_date) ==
          21.0  # 400 - 379, the raw increment, untouched
end

@testitem "truth_at leaves a stream with no break days untouched" begin
    using Dates: Date, Day

    include(joinpath(@__DIR__, "..", "scripts", "score_releases.jl"))

    grid_date(day) = Date(2026, 1, 1) + Day(day)
    ## Same shape and the same listed break day as the confirmed-cases
    ## fixture above, but stored under "reported cases": break days are
    ## only ever listed for the two confirmed streams, so this stream's
    ## truth must ignore them even though `obs` carries the fields.
    hist = (; days = [10, 17, 24, 31], counts = [100.0, 110.0, 379.0, 400.0])
    obs = (; cutoff = grid_date(40), reported_history = hist,
        confirmed_break_days = [24], confirmed_break_gross_cases = [97])

    made_date = grid_date(17)
    target_date = grid_date(31)
    @test truth_at(obs, grid_date, "reported cases", made_date, target_date) ==
          290.0  # 400 - 110, no correction applied
end

@testitem "baseline_draws' centre matches truth_at's correction" begin
    using Dates: Date, Day
    using Random: MersenneTwister
    using Statistics: mean

    include(joinpath(@__DIR__, "..", "scripts", "score_releases.jl"))

    grid_date(day) = Date(2026, 1, 1) + Day(day)
    ## Three vintages, so only two first differences are available and
    ## `baseline_draws` takes its Poisson fallback. That makes the sample
    ## mean a tight estimate of the centre (a standard error of 0.2 at
    ## these counts), which is the quantity under test; the random-walk
    ## spread is driven by the step pool `_history_diffs` builds and is
    ## covered by its own test below.
    hist = (; days = [17, 24, 31], counts = [110.0, 379.0, 400.0])
    obs = (; cutoff = grid_date(40), confirmed_history = hist,
        confirmed_break_days = [24], confirmed_break_gross_cases = [97])
    made_date = grid_date(31)

    ## The 7-day window (24, 31] opens on the break day, so it holds none
    ## and the centre is the raw 21: the control for the case below.
    draws = baseline_draws(
        obs, grid_date, "confirmed cases", made_date, 7, 4_000,
        MersenneTwister(1))
    @test isapprox(mean(draws), 21.0; atol = 1.0)

    ## The 14-day window (17, 31] does hold it, so the centre is the raw
    ## 290 less the 172 correction, the same corrected value `truth_at`
    ## gives for that window.
    draws14 = baseline_draws(
        obs, grid_date, "confirmed cases", made_date, 14, 4_000,
        MersenneTwister(1))
    @test isapprox(mean(draws14), 118.0; atol = 1.0)
end

@testitem "carry_break_days puts the declaration on a snapshot's grid" begin
    using Dates: Date, Day

    include(joinpath(@__DIR__, "..", "scripts", "score_releases.jl"))

    grid_date(day) = Date(2026, 1, 1) + Day(day)
    ## The current manifest declares two break days, at 25 and 60 days
    ## after the grid origin.
    obs = (; confirmed_break_days = [25, 60],
        confirmed_break_gross_cases = [97, 40],
        confirmed_break_gross_deaths = [62, 20])
    ## A snapshot cut off at day 30 on a grid of its own, carrying no
    ## declaration: day 25 is its day 35, and day 60 postdates it.
    ov = (; cutoff = grid_date(30), n = 40)

    carried = carry_break_days(ov, obs, grid_date)
    @test carried.confirmed_break_days == [35]
    @test carried.confirmed_break_gross_cases == [97]
    @test carried.confirmed_break_gross_deaths == [62]
    @test carried.cutoff == ov.cutoff

    ## Nothing to carry leaves the snapshot as it was.
    @test carry_break_days(ov, (;), grid_date) === ov
end

@testitem "break_correction sums every break day inside the window" begin
    using Dates: Date, Day

    include(joinpath(@__DIR__, "..", "scripts", "score_releases.jl"))

    grid_date(day) = Date(2026, 1, 1) + Day(day)
    ## Two listed break days: day 24 (net 269, gross 97, correction 172)
    ## and day 31 (net 90, gross 40, correction 50).
    hist = (; days = [10, 17, 24, 31], counts = [100.0, 110.0, 379.0, 469.0])
    obs = (; confirmed_history = hist,
        confirmed_break_days = [24, 31],
        confirmed_break_gross_cases = [97, 40])

    @test break_correction(
        obs, grid_date, "confirmed cases", grid_date(17), grid_date(31)) ==
          222.0  # 172 + 50, both days fall inside (17, 31]
    @test break_correction(
        obs, grid_date, "confirmed cases", grid_date(24), grid_date(31)) ==
          50.0  # only day 31 falls inside (24, 31]
end

@testitem "break_correction ignores a break day not yet in the history" begin
    using Dates: Date, Day

    include(joinpath(@__DIR__, "..", "scripts", "score_releases.jl"))

    grid_date(day) = Date(2026, 1, 1) + Day(day)
    ## A frozen history that stops before the listed break day, the state a
    ## baseline built from a snapshot older than the harmonisation sees.
    hist = (; days = [10, 17], counts = [100.0, 110.0])
    obs = (; confirmed_history = hist,
        confirmed_break_days = [24], confirmed_break_gross_cases = [97])

    @test break_correction(
        obs, grid_date, "confirmed cases", grid_date(10), grid_date(31)) ==
          0.0
end

@testitem "_history_diffs subtracts a break-day correction from its step" begin
    using Dates: Date, Day

    include(joinpath(@__DIR__, "..", "scripts", "score_releases.jl"))

    grid_date(day) = Date(2026, 1, 1) + Day(day)
    hist = (; days = [10, 17, 24, 31], counts = [100.0, 110.0, 379.0, 400.0])
    obs = (; confirmed_history = hist,
        confirmed_break_days = [24], confirmed_break_gross_cases = [97])
    made_date = grid_date(31)

    diffs = _history_diffs(obs, grid_date, "confirmed cases", hist, made_date)
    ## The 110 -> 379 step (day 17 -> 24) is corrected to 269 - 172 = 97;
    ## the other two steps are untouched.
    @test diffs == [(10.0, 7), (97.0, 7), (21.0, 7)]
end

@testitem "break_correction floors a gross above its own net step" begin
    using Dates: Date, Day

    include(joinpath(@__DIR__, "..", "scripts", "score_releases.jl"))

    grid_date(day) = Date(2026, 1, 1) + Day(day)
    ## A snapshot whose own vintage of the break day steps by only 40, less
    ## than the 97 printed alongside the day in the current manifest. The
    ## manifest loader refuses that pairing, so it can only arise on a
    ## snapshot the declaration was carried onto, and there the correction
    ## floors at zero rather than adding 57 back to the baseline.
    hist = (; days = [10, 17, 24], counts = [100.0, 110.0, 150.0])
    obs = (; confirmed_history = hist,
        confirmed_break_days = [24], confirmed_break_gross_cases = [97])

    @test break_correction(
        obs, grid_date, "confirmed cases", grid_date(17), grid_date(24)) ==
          0.0
end

@testitem "break_correction on a first vintage counts from zero" begin
    using Dates: Date, Day

    include(joinpath(@__DIR__, "..", "scripts", "score_releases.jl"))

    grid_date(day) = Date(2026, 1, 1) + Day(day)
    ## The break day is the stream's first vintage, so its step is the whole
    ## cumulative, matching what `cum_at` reads before a series has started.
    hist = (; days = [24, 31], counts = [379.0, 400.0])
    obs = (; confirmed_history = hist,
        confirmed_break_days = [24], confirmed_break_gross_cases = [97])

    @test break_correction(
        obs, grid_date, "confirmed cases", grid_date(17), grid_date(31)) ==
          282.0  # 379 - 97
end
