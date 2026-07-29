## Tests for the stream-reporting-coverage rule in scripts/score_releases.jl
## (stream_coverage_end, truth_at): a forecast is scorable only when its
## target falls inside the period its own truth source was still being
## reported, distinct from a target that is simply not yet observed.

@testitem "stream_coverage_end is the last vintage of a history stream" begin
    using Dates: Date, Day

    include(joinpath(@__DIR__, "..", "scripts", "score_releases.jl"))

    n = 40
    cutoff = Date(2026, 7, 15)
    grid_date(day) = cutoff - Day(n - day)
    obs = (; cutoff = cutoff, n = n,
        reported_history = (; days = [5, 12], counts = [50.0, 90.0]))
    @test stream_coverage_end(obs, grid_date, "reported cases") ==
          grid_date(12)
end

@testitem "stream_coverage_end is the last detection for exports" begin
    using Dates: Date, Day

    include(joinpath(@__DIR__, "..", "scripts", "score_releases.jl"))

    n = 40
    cutoff = Date(2026, 7, 15)
    grid_date(day) = cutoff - Day(n - day)
    obs = (; cutoff = cutoff, n = n, export_case_days = [8, 15, 20])
    @test stream_coverage_end(obs, grid_date, "exports") == grid_date(20)
end

@testitem "truth_at scores a target inside its stream's coverage" begin
    using Dates: Date, Day

    include(joinpath(@__DIR__, "..", "scripts", "score_releases.jl"))

    n = 40
    cutoff = Date(2026, 7, 15)
    grid_date(day) = cutoff - Day(n - day)
    obs = (; cutoff = cutoff, n = n,
        confirmed_history = (; days = [5, 12, 30],
            counts = [50.0, 90.0, 200.0]))
    made_date = grid_date(5)
    target_date = grid_date(12)  # inside coverage (last vintage at day 30)
    @test truth_at(obs, grid_date, "confirmed cases", made_date, target_date) ==
          40.0  # 90 - 50
end

@testitem "truth_at drops a target past its stream's coverage" begin
    using Dates: Date, Day

    include(joinpath(@__DIR__, "..", "scripts", "score_releases.jl"))

    n = 40
    cutoff = Date(2026, 7, 15)
    grid_date(day) = cutoff - Day(n - day)
    ## The reported-case history stops at day 12; a target on day 20 falls
    ## after it, so its unmoved cumulative total is not an observed zero.
    obs = (; cutoff = cutoff, n = n,
        reported_history = (; days = [5, 12], counts = [50.0, 90.0]))
    made_date = grid_date(5)
    target_date = grid_date(20)
    @test truth_at(obs, grid_date, "reported cases", made_date, target_date) ==
          :stopped_reporting
end

@testitem "truth_at drops a window straddling the coverage date" begin
    using Dates: Date, Day

    include(joinpath(@__DIR__, "..", "scripts", "score_releases.jl"))

    n = 40
    cutoff = Date(2026, 7, 15)
    grid_date(day) = cutoff - Day(n - day)
    obs = (; cutoff = cutoff, n = n,
        reported_history = (; days = [5, 12], counts = [50.0, 90.0]))
    ## made_date sits inside the reported period, target_date after it.
    made_date = grid_date(10)
    target_date = grid_date(18)
    @test truth_at(obs, grid_date, "reported cases", made_date, target_date) ==
          :stopped_reporting
end

@testitem "truth_at drops a window entirely past the coverage date" begin
    using Dates: Date, Day

    include(joinpath(@__DIR__, "..", "scripts", "score_releases.jl"))

    n = 40
    cutoff = Date(2026, 7, 15)
    grid_date(day) = cutoff - Day(n - day)
    obs = (; cutoff = cutoff, n = n,
        reported_history = (; days = [5, 12], counts = [50.0, 90.0]))
    made_date = grid_date(20)
    target_date = grid_date(27)
    @test truth_at(obs, grid_date, "reported cases", made_date, target_date) ==
          :stopped_reporting
end

@testitem "truth_at applies the coverage rule to a level stream" begin
    using Dates: Date, Day

    include(joinpath(@__DIR__, "..", "scripts", "score_releases.jl"))

    n = 40
    cutoff = Date(2026, 7, 15)
    grid_date(day) = cutoff - Day(n - day)
    obs = (; cutoff = cutoff, n = n,
        isolation_history = (; days = [5, 15], counts = [8.0, 20.0]))
    made_date = grid_date(5)

    inside = truth_at(
        obs, grid_date, "isolation beds", made_date, grid_date(15))
    @test inside == 20.0

    past = truth_at(obs, grid_date, "isolation beds", made_date, grid_date(25))
    @test past == :stopped_reporting
end

@testitem "truth_at applies the coverage rule to assembled exports" begin
    using Dates: Date, Day

    include(joinpath(@__DIR__, "..", "scripts", "score_releases.jl"))

    n = 40
    cutoff = Date(2026, 7, 15)
    grid_date(day) = cutoff - Day(n - day)
    obs = (; cutoff = cutoff, n = n, export_case_days = [5, 10])
    made_date = grid_date(5)

    inside = truth_at(obs, grid_date, "exports", made_date, grid_date(10))
    @test inside == 1.0  # one detection between day 5 and day 10

    past = truth_at(obs, grid_date, "exports", made_date, grid_date(20))
    @test past == :stopped_reporting
end

@testitem "truth_at tells not-yet-observed apart from stopped reporting" begin
    using Dates: Date, Day

    include(joinpath(@__DIR__, "..", "scripts", "score_releases.jl"))

    n = 40
    cutoff = Date(2026, 7, 15)
    grid_date(day) = cutoff - Day(n - day)
    ## Coverage ends well before the cut-off, but a target past the
    ## cut-off itself is reported as not-yet-observed, not stopped.
    obs = (; cutoff = cutoff, n = n,
        reported_history = (; days = [5, 12], counts = [50.0, 90.0]))
    made_date = grid_date(5)
    future = truth_at(
        obs, grid_date, "reported cases", made_date, cutoff + Day(7))
    @test future == :not_yet_observed
end

@testitem "truth_at keeps a genuine zero increment inside coverage" begin
    using Dates: Date, Day

    include(joinpath(@__DIR__, "..", "scripts", "score_releases.jl"))

    n = 40
    cutoff = Date(2026, 7, 15)
    grid_date(day) = cutoff - Day(n - day)
    ## The cumulative total is flat from day 5 to day 12 (a real week of no
    ## new cases), then rises again at day 20: coverage runs to day 20, so
    ## the flat week is a genuine zero, not a coverage artefact.
    obs = (; cutoff = cutoff, n = n,
        reported_history = (; days = [5, 12, 20], counts = [50.0, 50.0, 90.0]))
    made_date = grid_date(5)
    target_date = grid_date(12)
    @test truth_at(obs, grid_date, "reported cases", made_date, target_date) ==
          0.0

    ## A target past the last vintage is still dropped as a coverage
    ## artefact, so the rule tells the two zero-like cases apart rather
    ## than either scoring every zero or dropping every zero.
    past = truth_at(obs, grid_date, "reported cases", made_date, grid_date(30))
    @test past == :stopped_reporting
end

@testitem "score_release counts stopped-reporting groups apart" begin
    using Dates: Date, Day

    include(joinpath(@__DIR__, "..", "scripts", "score_releases.jl"))

    n = 40
    cutoff = Date(2026, 7, 15)
    grid_date(day) = cutoff - Day(n - day)
    ## "confirmed cases" is actively reported through the cut-off;
    ## "reported cases" stopped at day 12, well before the targets below.
    obs = (; cutoff = cutoff, n = n,
        confirmed_history = (; days = [5, 33, 40],
            counts = [10.0, 100.0, 150.0]),
        reported_history = (; days = [5, 12], counts = [50.0, 90.0]))
    made = string(grid_date(33))
    ## A confirmed-cases target beyond the cut-off (not yet observed) and a
    ## reported-cases target inside the window but past its own coverage.
    unobserved_target = string(cutoff + Day(3))
    stopped_target = string(grid_date(40))

    path = joinpath(mktempdir(), "forecast.csv")
    open(path, "w") do io
        println(io, "made_date,horizon,target_date,stream,draw,value")
        for d in 1:5
            println(io,
                join((made, 10, unobserved_target, "confirmed cases", d,
                        40 + d), ','))
            println(io,
                join((made, 7, stopped_target, "reported cases", d, 40 + d),
                    ','))
        end
    end

    result = score_release("results-vT.E.S", path, obs, grid_date)
    ## Each dropped group counts once per fit it carries (one score row per
    ## fit, the baseline never drawn for a group that is not scored): here
    ## one fit ("joint", the archive carries no `fit` column) per group.
    @test result.skipped == 1  # the confirmed-cases group
    @test result.stopped == 1  # the reported-cases group
    @test isempty(result.rows)
end

@testitem "stream_coverage_start is the first vintage of a history stream" begin
    using Dates: Date, Day

    include(joinpath(@__DIR__, "..", "scripts", "score_releases.jl"))

    n = 40
    cutoff = Date(2026, 7, 15)
    grid_date(day) = cutoff - Day(n - day)
    obs = (; cutoff = cutoff, n = n,
        isolation_history = (; days = [18, 25], counts = [200.0, 260.0]))
    @test stream_coverage_start(obs, grid_date, "isolation beds") ==
          grid_date(18)
end

@testitem "truth_at drops a window opening before reporting began" begin
    using Dates: Date, Day

    include(joinpath(@__DIR__, "..", "scripts", "score_releases.jl"))

    n = 40
    cutoff = Date(2026, 7, 15)
    grid_date(day) = cutoff - Day(n - day)
    ## The occupancy series begins at day 18, so a window opening at day 10
    ## reads its baseline and its increment off a series that has not
    ## started, which is the absence of a series rather than an occupancy of
    ## zero.
    obs = (; cutoff = cutoff, n = n,
        isolation_history = (; days = [18, 25], counts = [200.0, 260.0]))
    @test truth_at(obs, grid_date, "isolation beds",
        grid_date(10), grid_date(20)) === :not_yet_reporting
    ## A window opening on the first vintage is covered and scores.
    @test truth_at(obs, grid_date, "isolation beds",
        grid_date(18), grid_date(25)) == 260.0
end

@testitem "truth_at drops an incident window opening before reporting" begin
    using Dates: Date, Day

    include(joinpath(@__DIR__, "..", "scripts", "score_releases.jl"))

    n = 40
    cutoff = Date(2026, 7, 15)
    grid_date(day) = cutoff - Day(n - day)
    ## Measuring the increment from before the first vintage would count the
    ## whole cumulative total as new cases in the window.
    obs = (; cutoff = cutoff, n = n,
        recovered_history = (; days = [20, 30], counts = [100.0, 175.0]))
    @test truth_at(obs, grid_date, "recovered",
        grid_date(12), grid_date(30)) === :not_yet_reporting
    @test truth_at(obs, grid_date, "recovered",
        grid_date(20), grid_date(30)) == 75.0
end
