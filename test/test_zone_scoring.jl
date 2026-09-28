## Tests for the per-zone release scoring in scripts/score_releases.jl
## (zone_stream_label, parse_zone_stream, zone_stream_history,
## spans_zone_reattribution, score_zone_release): each zone is scored as its
## own stream against its own history, a week without a zone vintage on its
## last day or holding a reattribution is not scored, and a release without
## the zone archive is a quiet skip.

## Two zones of `ituri` and one of `nord_kivu` on shared vintage days, with
## each province's `unallocated` row. `ituri`'s unallocated count falls at
## day 24, a reattribution into its named zones. No vintage lands on day 21.
@testsnippet ZoneScoringFixture begin
    function _zone_fixture()
        days = [10, 17, 24, 31]
        z(counts) = (; days, counts)
        return Dict(
            "ituri" => Dict(
                "bunia" => z([10, 20, 36, 40]),
                "rwampara" => z([2, 3, 9, 11]),
                "unallocated" => z([5, 6, 1, 1])
            ),
            "nord_kivu" => Dict(
                "beni" => z([4, 6, 7, 12]),
                "unallocated" => z([0, 0, 0, 0])
            )
        )
    end
    function _zone_archive(path, rows; method = "predict")
        open(path, "w") do io
            println(
                io, "made_date,horizon,target_date,province,zone,stream," *
                    "draw,value" * (isnothing(method) ? "" : ",method")
            )
            for (made, target, prov, zone, vals) in rows
                for (d, v) in enumerate(vals)
                    row = (
                        made, 7, target, prov, zone, "confirmed cases", d, v,
                    )
                    m = isnothing(method) ? () : (method,)
                    println(io, join((row..., m...), ','))
                end
            end
        end
        return path
    end
end

@testitem "zone labels are recognised and cannot collide" begin
    include(joinpath(@__DIR__, "..", "scripts", "score_releases.jl"))

    @test zone_stream_label("confirmed cases", "ituri.bunia") ==
        "confirmed cases {ituri.bunia}"
    @test parse_zone_stream("confirmed cases {ituri.bunia}") ==
        (; stream = "confirmed cases", province = "ituri", zone = "bunia")
    @test has_stream_truth("confirmed cases {ituri.bunia}")

    ## National and province labels are left alone.
    @test isnothing(parse_zone_stream("confirmed cases"))
    @test isnothing(parse_zone_stream("confirmed cases [ituri]"))
    ## A stream the zone archive does not split, or a key with no province,
    ## has no truth source.
    @test !has_stream_truth("confirmed deaths {ituri.bunia}")
    @test !has_stream_truth("confirmed cases {bunia}")
end

@testitem "zone truth is the zone's own increment" setup = [
    ZoneScoringFixture,
] begin
    using Dates: Date, Day

    include(joinpath(@__DIR__, "..", "scripts", "score_releases.jl"))

    grid_date(day) = Date(2026, 1, 1) + Day(day)
    obs = (; cutoff = grid_date(35), zone_confirmed_history = _zone_fixture())

    hist, kind = stream_history(obs, "confirmed cases {nord_kivu.beni}")
    @test kind == :incident
    @test collect(hist.counts) == [4, 6, 7, 12]
    @test truth_at(
        obs, grid_date, "confirmed cases {nord_kivu.beni}",
        grid_date(24), grid_date(31)
    ) == 5.0

    ## A zone the manifest does not carry has no history and is not scored.
    @test isempty(
        first(stream_history(obs, "confirmed cases {ituri.mambasa}")).days
    )
    @test truth_at(
        obs, grid_date, "confirmed cases {ituri.mambasa}",
        grid_date(24), grid_date(31)
    ) isa Symbol
    bare = (; cutoff = grid_date(35))
    @test isempty(
        first(stream_history(bare, "confirmed cases {ituri.bunia}")).days
    )
end

@testitem "a zone week without a vintage on its last day is not scored" setup = [
    ZoneScoringFixture,
] begin
    using Dates: Date, Day

    include(joinpath(@__DIR__, "..", "scripts", "score_releases.jl"))

    grid_date(day) = Date(2026, 1, 1) + Day(day)
    obs = (; cutoff = grid_date(35), zone_confirmed_history = _zone_fixture())

    ## No zone table is dated day 21, so the nearest vintage either side
    ## would count a different week.
    @test truth_at(
        obs, grid_date, "confirmed cases {nord_kivu.beni}",
        grid_date(14), grid_date(21)
    ) === :no_target_vintage
end

@testitem "a zone week holding a reattribution is not scored" setup = [
    ZoneScoringFixture,
] begin
    using Dates: Date, Day

    include(joinpath(@__DIR__, "..", "scripts", "score_releases.jl"))

    grid_date(day) = Date(2026, 1, 1) + Day(day)
    obs = (; cutoff = grid_date(35), zone_confirmed_history = _zone_fixture())

    ## `ituri`'s unallocated row falls at day 24, so (17, 24] carries counts
    ## moved into the named zones, not new cases.
    @test truth_at(
        obs, grid_date, "confirmed cases {ituri.bunia}",
        grid_date(17), grid_date(24)
    ) === :spans_reattribution
    @test truth_at(
        obs, grid_date, "confirmed cases {ituri.bunia}",
        grid_date(24), grid_date(31)
    ) == 4.0
    ## Another province's zones are unaffected.
    @test truth_at(
        obs, grid_date, "confirmed cases {nord_kivu.beni}",
        grid_date(17), grid_date(24)
    ) == 1.0

    ## The national harmonisation backfill is not split by zone either.
    broken = merge(
        obs,
        (; confirmed_break_days = [30], confirmed_break_gross_cases = [1])
    )
    @test truth_at(
        broken, grid_date, "confirmed cases {nord_kivu.beni}",
        grid_date(24), grid_date(31)
    ) === :spans_break
end

@testitem "score_zone_release scores each zone as its own stream" setup = [
    ZoneScoringFixture,
] begin
    using Dates: Date, Day

    include(joinpath(@__DIR__, "..", "scripts", "score_releases.jl"))

    grid_date(day) = Date(2026, 1, 1) + Day(day)
    obs = (; cutoff = grid_date(35), zone_confirmed_history = _zone_fixture())
    made, target = grid_date(24), grid_date(31)
    path = _zone_archive(
        joinpath(mktempdir(), "zone_forecast.csv"),
        [
            (made, target, "ituri", "ituri.bunia", 2.0:6.0),
            (made, target, "ituri", "ituri.rwampara", 0.0:4.0),
            (made, target, "nord_kivu", "nord_kivu.beni", 3.0:7.0),
            (grid_date(17), grid_date(24), "ituri", "ituri.bunia", 2.0:6.0),
            (
                grid_date(14), grid_date(21), "nord_kivu",
                "nord_kivu.beni", 3.0:7.0,
            ),
        ]
    )
    result = score_zone_release("results-1", path, obs, grid_date)

    ## Two zones of one patch are two streams, not one pooled sample.
    joint = filter(r -> r.fit == JOINT_FIT, result.rows)
    @test sort([r.stream for r in joint]) == [
        "confirmed cases {ituri.bunia}",
        "confirmed cases {ituri.rwampara}",
        "confirmed cases {nord_kivu.beni}",
    ]
    @test all(r -> r.n_samples == 5, joint)
    @test sort(unique(r.fit for r in result.rows)) == ["baseline", "joint"]
    @test Set(r.stream for r in result.overlay) ==
        Set(r.stream for r in joint)
    ## The week into the reattribution and the week ending on a day with no
    ## zone table are each counted for the log.
    @test result.spans_reattribution == 1
    @test result.no_target_vintage == 1

    ## Only rows of the zone forecast's own method are scored.
    other = _zone_archive(
        joinpath(mktempdir(), "zone_forecast.csv"),
        [(made, target, "ituri", "ituri.bunia", 2.0:6.0)];
        method = "share"
    )
    @test isempty(score_zone_release("results-1", other, obs, grid_date).rows)
end

@testitem "a release without the zone archive is skipped" setup = [
    ZoneScoringFixture,
] begin
    using Dates: Date, Day

    include(joinpath(@__DIR__, "..", "scripts", "score_releases.jl"))

    grid_date(day) = Date(2026, 1, 1) + Day(day)
    obs = (; cutoff = grid_date(35), zone_confirmed_history = _zone_fixture())

    ## Every release published before the zone archive carries no asset.
    @test isnothing(score_zone_release("results-1", nothing, obs, grid_date))
    ## An archive with no method column is not scored.
    bare = _zone_archive(
        joinpath(mktempdir(), "zone_forecast.csv"),
        [(grid_date(24), grid_date(31), "ituri", "ituri.bunia", 1.0:3.0)];
        method = nothing
    )
    @test score_zone_release("results-1", bare, obs, grid_date) === :no_method
end

@testitem "the zone baseline skips a window holding a reattribution" setup = [
    ZoneScoringFixture,
] begin
    using Dates: Date, Day

    include(joinpath(@__DIR__, "..", "scripts", "score_releases.jl"))

    grid_date(day) = Date(2026, 1, 1) + Day(day)
    obs = (; cutoff = grid_date(35), zone_confirmed_history = _zone_fixture())
    bunia = "confirmed cases {ituri.bunia}"
    beni = "confirmed cases {nord_kivu.beni}"

    ## A baseline centred on (17, 24] would count `ituri`'s reattribution.
    @test !baseline_window_covered(obs, grid_date, bunia, grid_date(24), 7)
    @test baseline_window_covered(obs, grid_date, beni, grid_date(24), 7)

    ## The step pool leaves out the vintage whose window is (17, 24], so the
    ## one step runs from the day-17 total to the day-31 total.
    h = first(stream_history(obs, bunia))
    steps = _window_total_steps(obs, grid_date, bunia, h, grid_date(31), 7)
    @test steps ≈ [(4 - 10) / sqrt(14)]
    hb = first(stream_history(obs, beni))
    @test length(
        _window_total_steps(obs, grid_date, beni, hb, grid_date(31), 7)
    ) == 2
end
