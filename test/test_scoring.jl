## Tests for the forecast-scoring primitives in src/scoring.jl:
## crps_sample, log_crps_sample and score_draws, checked against the
## closed-form ensemble CRPS. Also covers select_daily_releases, which
## picks the releases scripts/score_releases.jl scores.

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
    ## version wins. v1.10.0 sorts below v1.9.0 as a string and above it as
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

@testitem "forecast_score_overview aggregates across horizon and release" begin
    using Dates: Date
    using DataFrames: DataFrame, nrow, names
    using BVDOutbreakSize: forecast_score_overview

    ## Two releases at the same horizon; "confirmed cases" carries a joint,
    ## an individual ("confirmed") and a baseline fit in both, "recovered"
    ## carries only the joint and the baseline (no individual model fits
    ## it).
    scores = DataFrame(
        release = ["r1", "r1", "r1", "r2", "r2", "r2", "r1", "r1", "r2", "r2"],
        made_date = fill(Date(2026, 6, 1), 10), horizon = fill(7, 10),
        stream = vcat(fill("confirmed cases", 6), fill("recovered", 4)),
        fit = vcat(["joint", "confirmed", "baseline"],
            ["joint", "confirmed", "baseline"],
            ["joint", "baseline"], ["joint", "baseline"]),
        crps = [2.0, 4.0, 8.0, 4.0, 6.0, 12.0, 3.0, 6.0, 3.0, 6.0],
        log_crps = [0.2, 0.4, 0.8, 0.4, 0.6, 1.2, 0.3, 0.6, 0.3, 0.6],
        dispersion = fill(0.1, 10), overprediction = fill(0.05, 10),
        underprediction = fill(0.05, 10),
        coverage_50 = fill(1.0, 10), coverage_90 = fill(1.0, 10),
        bias = fill(0.0, 10))

    out = forecast_score_overview(scores)
    @test "fit" in names(out)
    ## One row per non-baseline fit; no baseline row survives.
    @test sort(out.fit) == ["confirmed", "joint", "joint"]
    @test !("baseline" in out.fit)

    row(stream, fit) = only(
        out[(out.stream .== stream) .& (out.fit .== fit), :])
    ## rel_to_baseline is the ratio of the two fits' mean CRPS over the
    ## matched forecasts, pooled across both releases and horizons: joint
    ## mean (2+4)/2 = 3, baseline mean (8+12)/2 = 10, ratio 0.3. The
    ## absolute baseline and individual CRPS are not published columns.
    cc_joint = row("confirmed cases", "joint")
    @test cc_joint.n == 2
    @test cc_joint.crps == 3.0
    @test cc_joint.rel_to_baseline == 0.3
    ## Against the individual fit: joint mean 3.0, individual mean
    ## (4+6)/2 = 5.0, ratio 0.6.
    @test cc_joint.rel_to_individual == 0.6
    ## Same ratio on the log scale: joint log_crps mean (0.2+0.4)/2 = 0.3,
    ## individual log_crps mean (0.4+0.6)/2 = 0.5, ratio 0.6.
    @test cc_joint.log_rel_to_individual == 0.6

    cc_indiv = row("confirmed cases", "confirmed")
    @test cc_indiv.crps == 5.0
    @test cc_indiv.rel_to_baseline == 0.5  # 5.0 / 10.0
    ## The individual fit's own row carries no individual-fit comparison.
    @test ismissing(cc_indiv.rel_to_individual)
    @test ismissing(cc_indiv.log_rel_to_individual)

    ## "recovered" has no individual model: those columns stay missing on
    ## its joint row, and its skill is against its own baseline only.
    rec_joint = row("recovered", "joint")
    @test rec_joint.rel_to_baseline == 0.5  # 3.0 / 6.0
    @test ismissing(rec_joint.rel_to_individual)
    @test ismissing(rec_joint.log_rel_to_individual)
end

@testitem "forecast_score_overview drops a group with no matched baseline" begin
    using Dates: Date
    using DataFrames: DataFrame
    using BVDOutbreakSize: forecast_score_overview

    ## The baseline is scored under a different horizon to the joint fit,
    ## so the matched set is empty: the row is dropped rather than
    ## published with a fabricated ratio.
    scores = DataFrame(
        release = ["r1", "r1"], made_date = fill(Date(2026, 6, 1), 2),
        horizon = [7, 14], stream = fill("confirmed cases", 2),
        fit = ["joint", "baseline"], crps = [2.0, 8.0],
        log_crps = [0.2, 0.8], dispersion = [0.1, 0.1],
        overprediction = [0.05, 0.05], underprediction = [0.05, 0.05],
        coverage_50 = [1.0, 1.0], coverage_90 = [1.0, 1.0],
        bias = [0.0, 0.0])

    out = forecast_score_overview(scores)
    @test isempty(out)
end

@testitem "forecast_score_overview guards a zero baseline mean" begin
    using Dates: Date
    using DataFrames: DataFrame
    using BVDOutbreakSize: forecast_score_overview

    ## A baseline scored to an exact zero CRPS (a coincidental perfect
    ## persistence match) would divide by zero; the guard reports the row
    ## with rel_to_baseline missing instead of Inf.
    scores = DataFrame(
        release = ["r1", "r1"], made_date = fill(Date(2026, 6, 1), 2),
        horizon = [7, 7], stream = fill("confirmed cases", 2),
        fit = ["joint", "baseline"], crps = [2.0, 0.0],
        log_crps = [0.2, 0.0], dispersion = [0.1, 0.0],
        overprediction = [0.05, 0.0], underprediction = [0.05, 0.0],
        coverage_50 = [1.0, 1.0], coverage_90 = [1.0, 1.0],
        bias = [0.0, 0.0])

    out = forecast_score_overview(scores)
    @test ismissing(only(out.rel_to_baseline))
    @test !any(isinf, skipmissing(out.rel_to_baseline))
end

@testitem "forecast_score_overview returns a typed empty frame" begin
    using Dates: Date
    using DataFrames: DataFrame, nrow, names
    using BVDOutbreakSize: forecast_score_overview

    empty_scores = DataFrame(
        release = String[], made_date = Date[], horizon = Int[],
        stream = String[], fit = String[], crps = Float64[],
        log_crps = Float64[], dispersion = Float64[],
        overprediction = Float64[], underprediction = Float64[],
        coverage_50 = Float64[], coverage_90 = Float64[], bias = Float64[])

    out = forecast_score_overview(empty_scores)
    @test nrow(out) == 0
    @test "stream" in names(out)
    @test "fit" in names(out)
    @test "rel_to_baseline" in names(out)
    @test "rel_to_individual" in names(out)
    @test !("crps_baseline" in names(out))
    @test !("crps_individual" in names(out))
end

@testitem "drop_individual_fit_columns omits the individual comparison" begin
    using Dates: Date
    using DataFrames: DataFrame, names
    using BVDOutbreakSize: forecast_score_overview, drop_individual_fit_columns

    ## A frozen-style table: only a "frozen" fit and its baseline, so the
    ## individual-fit columns are always missing throughout.
    scores = DataFrame(
        release = ["r1", "r1"], made_date = fill(Date(2026, 6, 1), 2),
        horizon = [7, 7], stream = fill("confirmed cases", 2),
        fit = ["frozen", "baseline"], crps = [2.0, 8.0],
        log_crps = [0.2, 0.8], dispersion = [0.1, 0.1],
        overprediction = [0.05, 0.05], underprediction = [0.05, 0.05],
        coverage_50 = [1.0, 1.0], coverage_90 = [1.0, 1.0],
        bias = [0.0, 0.0])

    table = forecast_score_overview(scores)
    @test "rel_to_individual" in names(table)

    out = drop_individual_fit_columns(table)
    @test !("rel_to_individual" in names(out))
    @test !("log_rel_to_individual" in names(out))
    ## Every other column, and the row itself, survives untouched.
    @test "rel_to_baseline" in names(out)
    @test only(out.fit) == "frozen"
end

@testitem "drop_degenerate_fit_column drops a single-valued fit column" begin
    using Dates: Date
    using DataFrames: DataFrame, names, nrow
    using BVDOutbreakSize: forecast_score_overview, drop_degenerate_fit_column

    ## The frozen-evaluation shape: only ever "frozen" and its baseline, so
    ## `fit` carries exactly one value ("frozen") throughout the table.
    scores = DataFrame(
        release = ["r1", "r1"], made_date = fill(Date(2026, 6, 1), 2),
        horizon = [7, 7], stream = fill("confirmed cases", 2),
        fit = ["frozen", "baseline"], crps = [2.0, 8.0],
        log_crps = [0.2, 0.8], dispersion = [0.1, 0.1],
        overprediction = [0.05, 0.05], underprediction = [0.05, 0.05],
        coverage_50 = [1.0, 1.0], coverage_90 = [1.0, 1.0],
        bias = [0.0, 0.0])
    table = forecast_score_overview(scores)

    out = drop_degenerate_fit_column(table)
    @test !("fit" in names(out))
    ## Every other column survives untouched, in place.
    @test "bias" in names(out)
    @test "dispersion" in names(out)
    @test "overprediction" in names(out)
    @test "underprediction" in names(out)
    @test "coverage_50" in names(out)
    @test "coverage_90" in names(out)
    @test "crps" in names(out)
    @test "rel_to_baseline" in names(out)
    @test nrow(out) == nrow(table)
end

@testitem "drop_degenerate_fit_column errors on a multi-fit table" begin
    using Dates: Date
    using DataFrames: DataFrame
    using BVDOutbreakSize: forecast_score_overview, drop_degenerate_fit_column

    ## The ordinary joint-vs-baseline shape: joint, confirmed and baseline,
    ## so `fit` carries more than one value and must not be dropped.
    scores = DataFrame(
        release = fill("r1", 3), made_date = fill(Date(2026, 6, 1), 3),
        horizon = fill(7, 3), stream = fill("confirmed cases", 3),
        fit = ["joint", "confirmed", "baseline"], crps = [2.0, 2.5, 8.0],
        log_crps = [0.2, 0.25, 0.8], dispersion = fill(0.1, 3),
        overprediction = fill(0.05, 3), underprediction = fill(0.05, 3),
        coverage_50 = fill(1.0, 3), coverage_90 = fill(1.0, 3),
        bias = fill(0.0, 3))
    table = forecast_score_overview(scores)

    @test_throws ErrorException drop_degenerate_fit_column(table)
end

@testitem "drop_degenerate_fit_column leaves an empty table's fit column" begin
    using BVDOutbreakSize: forecast_score_overview, drop_degenerate_fit_column
    using DataFrames: DataFrame, names, nrow

    empty = forecast_score_overview(DataFrame())
    out = drop_degenerate_fit_column(empty)
    @test "fit" in names(out)
    @test nrow(out) == 0
end

@testitem "forecast_score_by_horizon keeps a separate row per horizon" begin
    using Dates: Date
    using DataFrames: DataFrame
    using BVDOutbreakSize: forecast_score_by_horizon

    scores = DataFrame(
        release = ["r1", "r1", "r2", "r2"],
        made_date = fill(Date(2026, 6, 1), 4),
        horizon = [7, 7, 14, 14], stream = fill("confirmed cases", 4),
        fit = ["joint", "baseline", "joint", "baseline"],
        crps = [2.0, 8.0, 6.0, 12.0], log_crps = [0.2, 0.8, 0.6, 1.2],
        dispersion = fill(0.1, 4), overprediction = fill(0.05, 4),
        underprediction = fill(0.05, 4), coverage_50 = fill(1.0, 4),
        coverage_90 = fill(1.0, 4), bias = fill(0.0, 4))

    out = forecast_score_by_horizon(scores)
    @test sort(out.horizon) == [7, 14]
    h7 = only(out[out.horizon .== 7, :])
    h14 = only(out[out.horizon .== 14, :])
    @test h7.rel_to_baseline == 0.25   # 2.0 / 8.0
    @test h14.rel_to_baseline == 0.5   # 6.0 / 12.0
end

@testitem "forecast_score_by_release averages across horizon" begin
    using Dates: Date
    using DataFrames: DataFrame
    using BVDOutbreakSize: forecast_score_by_release

    ## One made date, two horizons: the by-release row pools both, unlike
    ## the by-horizon table which would keep them apart.
    scores = DataFrame(
        release = ["r1", "r1", "r1", "r1"],
        made_date = fill(Date(2026, 6, 1), 4),
        horizon = [7, 7, 14, 14], stream = fill("confirmed cases", 4),
        fit = ["joint", "baseline", "joint", "baseline"],
        crps = [2.0, 8.0, 6.0, 12.0], log_crps = [0.2, 0.8, 0.6, 1.2],
        dispersion = fill(0.1, 4), overprediction = fill(0.05, 4),
        underprediction = fill(0.05, 4), coverage_50 = fill(1.0, 4),
        coverage_90 = fill(1.0, 4), bias = fill(0.0, 4))

    out = forecast_score_by_release(scores)
    row = only(out)
    @test row.made_date == Date(2026, 6, 1)
    @test row.n == 2
    @test row.crps == 4.0          # mean(2.0, 6.0)
    @test row.rel_to_baseline == 0.4  # mean(2.0, 6.0) / mean(8.0, 12.0)
end

@testitem "score_release scores every mapped stream, skips the rest" begin
    using Dates: Date, Day
    using DataFrames: DataFrame

    ## Load the scorer's functions without its driver (guarded on
    ## PROGRAM_FILE), so score_release can be exercised directly.
    include(joinpath(@__DIR__, "..", "scripts", "score_releases.jl"))

    ## A synthetic now-observed manifest: a 40-day grid, each incident
    ## history stepping from 100 to 150 over the last week, isolation beds
    ## at 20, and three dated export detections in the final week. Each
    ## incident history opens a fortnight earlier so the baseline's own
    ## lookback window is covered (see `baseline_window_covered`).
    n = 40
    cutoff = Date(2026, 7, 15)
    inc = (; days = [19, 33, 40], counts = [60, 100, 150])
    obs = (; cutoff = cutoff, n = n,
        reported_history = inc, deaths_history = inc,
        confirmed_history = inc, confirmed_deaths_history = inc,
        recovered_history = inc,
        isolation_history = (; days = [33, 40], counts = [18, 20]),
        export_case_days = [35, 38, 40])
    grid_date(day) = obs.cutoff - Day(obs.n - day)

    ## A stream_forecasts.csv carrying all six real labels (with two fits of
    ## confirmed cases) plus one unmapped label that must not abort scoring.
    made = string(grid_date(33))
    target = string(grid_date(40))
    labels = [("reported cases", "cases"), ("suspected deaths", "deaths"),
        ("confirmed cases", "joint"), ("confirmed cases", "confirmed"),
        ("exports", "exports"), ("isolation beds", "treatment"),
        ("nonsense stream", "joint")]
    path = joinpath(mktempdir(), "stream_forecasts.csv")
    open(path, "w") do io
        println(io, "made_date,horizon,target_date,stream,draw,value,fit")
        for (stream, fit) in labels, d in 1:5

            println(io,
                join((made, 7, target, stream, d, 40 + d, fit), ','))
        end
    end

    result = score_release("results-vT.E.S", path, obs, grid_date)
    scored = DataFrame(result.rows)
    scored_streams = Set(scored.stream)

    ## Every mapped stream is scored, including exports (its truth built
    ## from export_case_days) and the two clean additions.
    for s in ["reported cases", "suspected deaths", "confirmed cases",
        "exports", "isolation beds"]
        @test s in scored_streams
    end
    ## The unmapped label is dropped, not scored, and did not abort: the
    ## mapped streams still produced rows.
    @test !("nonsense stream" in scored_streams)

    ## Each scored stream carries its fit rows plus one baseline; confirmed
    ## cases carries both its fits against the single baseline.
    cc = scored[scored.stream .== "confirmed cases", :]
    @test Set(cc.fit) == Set(["joint", "confirmed", "baseline"])
    ex = scored[scored.stream .== "exports", :]
    @test Set(ex.fit) == Set(["exports", "baseline"])
    ## The exports fit's five draws are all scored (n_samples), confirming
    ## its truth was assembled rather than the stream skipped.
    @test only(ex[ex.fit .== "exports", :].n_samples) == 5
end

@testitem "score_release tags a fit-less archive with default_fit" begin
    using Dates: Date, Day
    using DataFrames: DataFrame

    include(joinpath(@__DIR__, "..", "scripts", "score_releases.jl"))

    ## The frozen-fit archive (forecast_frozen.csv) carries the same schema as
    ## forecast.csv but no `fit` column, and is stamped with a past made date.
    n = 40
    cutoff = Date(2026, 7, 15)
    ## The confirmed history opens well before the made date so the
    ## baseline's own lookback window is covered (see
    ## `baseline_window_covered`).
    inc = (; days = [15, 26, 33], counts = [60, 100, 150])
    obs = (; cutoff = cutoff, n = n,
        reported_history = inc, deaths_history = inc,
        confirmed_history = inc, confirmed_deaths_history = inc,
        recovered_history = inc,
        isolation_history = (; days = [26, 33], counts = [18, 20]),
        export_case_days = [30, 32])
    grid_date(day) = obs.cutoff - Day(obs.n - day)

    ## A made date one week before the cut-off, so its 7-day target is already
    ## observed and the group is scored rather than skipped.
    made = string(grid_date(26))
    target = string(grid_date(33))
    path = joinpath(mktempdir(), "forecast_frozen.csv")
    open(path, "w") do io
        println(io, "made_date,horizon,target_date,stream,draw,value")
        for d in 1:5
            println(io,
                join((made, 7, target, "confirmed cases", d, 40 + d), ','))
        end
    end

    result = score_release("results-vT.E.S", path, obs, grid_date;
        default_fit = "frozen")
    scored = DataFrame(result.rows)
    ## The fit-less rows are tagged with default_fit, and the persistence
    ## baseline is still produced alongside.
    @test Set(scored.fit) == Set(["frozen", "baseline"])
    @test result.skipped == 0
end

@testitem "crps_sample matches the closed-form ensemble CRPS" begin
    using BVDOutbreakSize: crps_sample

    ## CRPS = mean|xᵢ - obs| - ½ mean|xᵢ - xⱼ|. For [1,2,3,4,5] at 2.5:
    ## mean|xᵢ - 2.5| = 6.5/5 = 1.3, and the pairwise term is 0.8, so 0.5.
    @test crps_sample(2.5, [1.0, 2.0, 3.0, 4.0, 5.0]) ≈ 0.5
    ## For [2,4] at 3: mean abs error 1, pairwise ½·1 = 0.5, so 0.5.
    @test crps_sample(3.0, [2.0, 4.0]) ≈ 0.5
    ## Symmetric ensemble about the observation.
    @test crps_sample(0.0, [-1.0, 1.0]) ≈ 0.5
    ## Order independence.
    @test crps_sample(2.5, [5.0, 1.0, 3.0, 2.0, 4.0]) ≈
          crps_sample(2.5, [1.0, 2.0, 3.0, 4.0, 5.0])
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
    @test isnan(s.dispersion)
    @test isnan(s.overprediction)
    @test isnan(s.underprediction)
    @test isnan(s.bias)
    @test s.coverage_50 == false
    @test s.coverage_90 == false
    @test s.n == 0
end

@testitem "crps_decomposition sums to the CRPS and stays non-negative" begin
    using Random: MersenneTwister
    using BVDOutbreakSize: crps_sample, crps_decomposition

    ## A mixed ensemble straddling the observation, so both the dispersion
    ## and the two directional components are all in play at once.
    rng = MersenneTwister(3)
    samples = 100.0 .+ 15.0 .* randn(rng, 41)
    obs = 105.0
    d = crps_decomposition(obs, samples)
    @test d.dispersion >= 0
    @test d.overprediction >= 0
    @test d.underprediction >= 0
    @test d.dispersion + d.overprediction + d.underprediction ≈
          crps_sample(obs, samples)
end

@testitem "crps_decomposition of a point-mass ensemble has no dispersion" begin
    using BVDOutbreakSize: crps_sample, crps_decomposition

    ## A point-mass forecast has no spread of its own, so its CRPS is all
    ## overprediction or underprediction depending on which side `obs` sits.
    above = crps_decomposition(5.0, fill(20.0, 50))
    @test above.dispersion ≈ 0.0
    @test above.underprediction ≈ 0.0
    @test above.overprediction ≈ crps_sample(5.0, fill(20.0, 50))

    below = crps_decomposition(20.0, fill(5.0, 50))
    @test below.dispersion ≈ 0.0
    @test below.overprediction ≈ 0.0
    @test below.underprediction ≈ crps_sample(20.0, fill(5.0, 50))

    ## A point-mass forecast equal to the observation scores (and
    ## decomposes to) zero throughout.
    exact = crps_decomposition(5.0, fill(5.0, 50))
    @test exact.dispersion ≈ 0.0
    @test exact.overprediction ≈ 0.0
    @test exact.underprediction ≈ 0.0
end

@testitem "crps_decomposition splits an ensemble entirely above obs" begin
    using Random: MersenneTwister
    using BVDOutbreakSize: crps_sample, crps_decomposition

    ## Every draw sits above the observation: the ensemble reads too high,
    ## so the extra cost beyond its own spread is all overprediction.
    rng = MersenneTwister(4)
    samples = 200.0 .+ 5.0 .* rand(rng, 60)
    obs = 100.0
    d = crps_decomposition(obs, samples)
    @test d.underprediction ≈ 0.0
    @test d.dispersion >= 0
    @test d.overprediction >= 0
    @test d.dispersion + d.overprediction ≈ crps_sample(obs, samples)
end

@testitem "crps_decomposition splits an ensemble entirely below obs" begin
    using Random: MersenneTwister
    using BVDOutbreakSize: crps_sample, crps_decomposition

    ## Every draw sits below the observation: the ensemble reads too low,
    ## so the extra cost beyond its own spread is all underprediction.
    rng = MersenneTwister(5)
    samples = 50.0 .+ 5.0 .* rand(rng, 60)
    obs = 200.0
    d = crps_decomposition(obs, samples)
    @test d.overprediction ≈ 0.0
    @test d.dispersion >= 0
    @test d.underprediction >= 0
    @test d.dispersion + d.underprediction ≈ crps_sample(obs, samples)
end

@testitem "select_fit_role splits a summary into joint and individual" begin
    using Dates: Date
    using DataFrames: DataFrame, nrow, names
    using BVDOutbreakSize: forecast_score_overview, forecast_score_by_horizon,
                           forecast_score_by_release, select_fit_role

    ## "confirmed cases" carries a joint, an individual ("confirmed") and a
    ## baseline fit; "reported cases" carries a joint, a differently named
    ## individual ("cases") and a baseline.
    scores = DataFrame(
        release = fill("r1", 6), made_date = fill(Date(2026, 6, 1), 6),
        horizon = fill(7, 6),
        stream = vcat(fill("confirmed cases", 3), fill("reported cases", 3)),
        fit = ["joint", "confirmed", "baseline",
            "joint", "cases", "baseline"],
        crps = [2.0, 4.0, 8.0, 3.0, 6.0, 12.0],
        log_crps = [0.2, 0.4, 0.8, 0.3, 0.6, 1.2],
        dispersion = fill(0.1, 6), overprediction = fill(0.05, 6),
        underprediction = fill(0.05, 6),
        coverage_50 = fill(1.0, 6), coverage_90 = fill(1.0, 6),
        bias = fill(0.0, 6))

    ## Each builder keeps every non-baseline fit, so the unfiltered table
    ## carries both roles and the relative-skill figure can compare them.
    for table in (forecast_score_overview(scores),
        forecast_score_by_horizon(scores), forecast_score_by_release(scores))
        @test sort(table.fit) == ["cases", "confirmed", "joint", "joint"]

        ## The joint role is the joint model's rows alone: no individual
        ## fit appears in a table headed as the joint model's.
        joint = select_fit_role(table, "joint")
        @test unique(joint.fit) == ["joint"]
        @test nrow(joint) == 2
        ## The individual role is each stream's own fit, whatever its id.
        indiv = select_fit_role(table, "individual")
        @test sort(indiv.fit) == ["cases", "confirmed"]
        ## The two roles partition the table, and the columns are unchanged.
        @test nrow(joint) + nrow(indiv) == nrow(table)
        @test names(joint) == names(table)
        @test names(indiv) == names(table)
        ## Baseline rows are excluded by the builders, so that role is empty.
        @test nrow(select_fit_role(table, "baseline")) == 0
    end
end

@testitem "select_fit_role reads a frozen fit in the joint role" begin
    using Dates: Date
    using DataFrames: DataFrame, nrow
    using BVDOutbreakSize: forecast_score_overview, select_fit_role

    ## The frozen evaluation is the joint model re-fit at a past cut-off,
    ## so its rows belong to the joint role rather than the individual one.
    scores = DataFrame(
        release = fill("r1", 2), made_date = fill(Date(2026, 6, 1), 2),
        horizon = [7, 7], stream = fill("confirmed cases", 2),
        fit = ["frozen", "baseline"], crps = [2.0, 8.0],
        log_crps = [0.2, 0.8], dispersion = fill(0.1, 2),
        overprediction = fill(0.05, 2), underprediction = fill(0.05, 2),
        coverage_50 = fill(1.0, 2), coverage_90 = fill(1.0, 2),
        bias = fill(0.0, 2))

    table = forecast_score_overview(scores)
    @test only(select_fit_role(table, "joint").fit) == "frozen"
    @test nrow(select_fit_role(table, "individual")) == 0
end

@testitem "select_fit_role errors on an unknown role" begin
    using DataFrames: DataFrame
    using BVDOutbreakSize: forecast_score_overview, select_fit_role

    empty = forecast_score_overview(DataFrame())
    ## A misspelt role must not read as a role with nothing scored.
    @test_throws ErrorException select_fit_role(empty, "individuals")
    @test_throws ErrorException select_fit_role(empty, "Joint")
end

@testitem "truth_at scores an increment for a cumulative stream" begin
    using Dates: Date, Day

    include(joinpath(@__DIR__, "..", "scripts", "score_releases.jl"))

    ## Every cumulative-count stream is scored on what it added over the
    ## window, matching what `forecast_archive` stores (`*_new` columns),
    ## never on the cumulative total standing at the target. The two are far
    ## apart here: the histories run from 1000 to 1150 while the window adds
    ## 50, so a level score cannot be mistaken for an increment score.
    n = 40
    cutoff = Date(2026, 7, 15)
    grid_date(day) = cutoff - Day(n - day)
    cum = (; days = [26, 33, 40], counts = [1000, 1100, 1150])
    obs = (; cutoff = cutoff, n = n,
        reported_history = cum, deaths_history = cum,
        confirmed_history = cum, confirmed_deaths_history = cum,
        recovered_history = cum, onset_report_history = cum,
        isolation_history = (; days = [26, 33, 40], counts = [18, 19, 20]),
        treatment_confirmed_incare_history =
        (; days = [26, 33, 40], counts = [8, 9, 12]),
        treatment_suspect_incare_history =
        (; days = [26, 33, 40], counts = [10, 10, 8]),
        export_case_days = [27, 35, 38, 40])
    made_date = grid_date(33)
    target_date = grid_date(40)

    for stream in ("reported cases", "suspected deaths", "confirmed cases",
        "confirmed deaths", "recovered", "onset reports")
        @test truth_at(obs, grid_date, stream, made_date, target_date) == 50.0
    end
    ## The assembled export stream counts detections the same way: three of
    ## the four landed inside the window.
    @test truth_at(obs, grid_date, "exports", made_date, target_date) == 3.0
end

@testitem "truth_at scores a level for a stock stream" begin
    using Dates: Date, Day

    include(joinpath(@__DIR__, "..", "scripts", "score_releases.jl"))

    ## Bed occupancy is a stock, not a cumulative count, so its truth is the
    ## occupancy standing at the target rather than a change over the window.
    ## Each ward sub-stock is read the same way, including the one that fell.
    n = 40
    cutoff = Date(2026, 7, 15)
    grid_date(day) = cutoff - Day(n - day)
    obs = (; cutoff = cutoff, n = n,
        isolation_history = (; days = [26, 33, 40], counts = [18, 19, 20]),
        treatment_confirmed_incare_history =
        (; days = [26, 33, 40], counts = [8, 9, 12]),
        treatment_suspect_incare_history =
        (; days = [26, 33, 40], counts = [10, 10, 8]))
    made_date = grid_date(33)
    target_date = grid_date(40)

    @test truth_at(obs, grid_date, "isolation beds", made_date,
        target_date) == 20.0
    @test truth_at(obs, grid_date, "treatment beds", made_date,
        target_date) == 12.0
    @test truth_at(obs, grid_date, "isolation beds (suspected)", made_date,
        target_date) == 8.0
end

@testitem "every scored stream is scored on its declared basis" begin
    using Dates: Date, Day

    include(joinpath(@__DIR__, "..", "scripts", "score_releases.jl"))

    ## The basis each stream is scored on has to match what the archive
    ## stores for it: the cumulative-count streams are archived as new counts
    ## over the horizon and the occupancy stocks as levels (see
    ## `forecast_archive` in `src/forecast.jl`).
    ##
    ## Driving this from the stream maps rather than a list means a stream
    ## added on the wrong basis fails here, and running each one through
    ## `truth_at` means the declared basis has to be the one it is actually
    ## scored on. An increment depends on where the window opens; a level
    ## does not. That separates the two without pinning a value per stream.
    n = 40
    cutoff = Date(2026, 7, 15)
    grid_date(day) = cutoff - Day(n - day)
    cum = (; days = [26, 33, 40], counts = [1000, 1100, 1150])
    obs = (; cutoff = cutoff, n = n,
        reported_history = cum, deaths_history = cum,
        confirmed_history = cum, confirmed_deaths_history = cum,
        recovered_history = cum, onset_report_history = cum,
        isolation_history = (; days = [26, 33, 40], counts = [18, 19, 20]),
        treatment_confirmed_incare_history =
        (; days = [26, 33, 40], counts = [8, 9, 12]),
        treatment_suspect_incare_history =
        (; days = [26, 33, 40], counts = [10, 10, 8]),
        export_case_days = [27, 35, 38, 40])
    target_date = grid_date(40)
    early, late = grid_date(26), grid_date(33)

    kinds = merge(Dict(k => v[2] for (k, v) in STREAM_HISTORY),
        Dict(STREAM_ASSEMBLED))
    @test !isempty(kinds)

    for (stream, kind) in kinds
        from_early = truth_at(obs, grid_date, stream, early, target_date)
        from_late = truth_at(obs, grid_date, stream, late, target_date)
        if kind === :incident
            @test from_early > from_late
        else
            @test kind === :level
            @test from_early == from_late
        end
    end
end
