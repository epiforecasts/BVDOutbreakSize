## Tests for the per-province release scoring in scripts/score_releases.jl
## (province_stream_label, parse_province_stream, province_stream_history,
## spans_confirmed_break, score_province_release): a patch's truth is its
## member provinces' cumulative counts pooled vintage by vintage, a window
## holding a harmonisation-break day is not scored, and a release without the
## per-province archive is a quiet skip.

## The four source provinces the `other` patch pools, on shared vintage days,
## alongside a single-province patch. `sud_kivu` falls between the second and
## third vintage, so pooling before differencing is visible in the totals. A
## `@testsnippet`, not a top-level function: each `@testitem` below sees only
## its own code range, so a shared fixture has to come through the snippet.
@testsnippet ProvinceFixture begin
    _province_fixture() = Dict(
        "ituri" => (; days = [10, 17, 24], counts = [40, 70, 95]),
        "sud_kivu" => (; days = [10, 17, 24], counts = [5, 9, 7]),
        "tshopo" => (; days = [10, 17, 24], counts = [1, 3, 6]),
        "bas_uele" => (; days = [10, 17, 24], counts = [0, 2, 5]),
        "sud_ubangi" => (; days = [10, 17, 24], counts = [2, 2, 8])
    )
end

@testitem "province truth pools a patch's member provinces" setup = [
    ProvinceFixture,
] begin
    using Dates: Date, Day

    include(joinpath(@__DIR__, "..", "scripts", "score_releases.jl"))

    grid_date(day) = Date(2026, 1, 1) + Day(day)
    obs = (;
        cutoff = grid_date(30),
        province_confirmed_history = _province_fixture(),
    )

    ## `other` pools four provinces: 8, 16 and 26 cumulative by each vintage.
    hist, kind = stream_history(obs, "confirmed cases [other]")
    @test kind == :incident
    @test hist.days == [10, 17, 24]
    @test collect(hist.counts) == [8, 16, 26]

    ## A patch of one province is that province's own history.
    ituri, _ = stream_history(obs, "confirmed cases [ituri]")
    @test collect(ituri.counts) == [40, 70, 95]

    ## The scored window is the pooled increment over it, so `sud_kivu`'s
    ## fall nets against its neighbours' rises rather than being clamped away.
    @test truth_at(
        obs, grid_date, "confirmed cases [other]",
        grid_date(17), grid_date(24)
    ) == 10.0
    @test truth_at(
        obs, grid_date, "confirmed cases [ituri]",
        grid_date(17), grid_date(24)
    ) == 25.0
end

@testitem "province labels are recognised and cannot collide" begin
    include(joinpath(@__DIR__, "..", "scripts", "score_releases.jl"))

    @test province_stream_label("confirmed deaths", "nord_kivu") ==
        "confirmed deaths [nord_kivu]"
    @test has_stream_truth("confirmed deaths [nord_kivu]")
    @test parse_province_stream("confirmed cases [haut_uele]") ==
        (; stream = "confirmed cases", province = "haut_uele")

    ## Every national label is left alone, including the one that already
    ## carries a parenthesised qualifier.
    @test isnothing(parse_province_stream("confirmed cases"))
    @test isnothing(parse_province_stream("isolation beds (suspected)"))

    ## A stream the archive does not split, or a patch the model does not
    ## carry, has no truth source and is dropped like any unmapped label.
    @test !has_stream_truth("isolation beds [ituri]")
    @test !has_stream_truth("confirmed cases [kinshasa]")
end

@testitem "a patch with a member missing has no history" setup = [
    ProvinceFixture,
] begin
    using Dates: Date, Day

    include(joinpath(@__DIR__, "..", "scripts", "score_releases.jl"))

    grid_date(day) = Date(2026, 1, 1) + Day(day)
    partial = filter(p -> p.first != "sud_ubangi", _province_fixture())
    obs = (; cutoff = grid_date(30), province_confirmed_history = partial)

    ## A pool short one member would understate the patch, so the patch has
    ## no history at all and its groups are never scored.
    hist, _ = stream_history(obs, "confirmed cases [other]")
    @test isempty(hist.days)
    @test truth_at(
        obs, grid_date, "confirmed cases [other]",
        grid_date(17), grid_date(24)
    ) isa Symbol

    ## A manifest carrying no per-province block at all reads the same way.
    bare = (; cutoff = grid_date(30))
    @test isempty(first(stream_history(bare, "confirmed cases [ituri]")).days)
end

@testitem "a province window holding a break day is not scored" setup = [
    ProvinceFixture,
] begin
    using Dates: Date, Day

    include(joinpath(@__DIR__, "..", "scripts", "score_releases.jl"))

    grid_date(day) = Date(2026, 1, 1) + Day(day)
    obs = (;
        cutoff = grid_date(30),
        province_confirmed_history = _province_fixture(),
        confirmed_break_days = [20], confirmed_break_gross_cases = [10],
    )

    ## (17, 24] holds the break day at 20. No per-province split of that
    ## day's backfill is published, so the window is dropped rather than
    ## scored against a base integration.
    @test truth_at(
        obs, grid_date, "confirmed cases [other]",
        grid_date(17), grid_date(24)
    ) === :spans_break
    ## (10, 17] does not hold it, so it is scored as usual.
    @test truth_at(
        obs, grid_date, "confirmed cases [other]",
        grid_date(10), grid_date(17)
    ) == 8.0

    ## The rule is the province path's alone: a national stream corrects the
    ## window instead of dropping it.
    @test !spans_confirmed_break(
        obs, grid_date, "confirmed cases",
        grid_date(17), grid_date(24)
    )
end

@testitem "score_release scores the province archive per patch" setup = [
    ProvinceFixture,
] begin
    using Dates: Date, Day

    include(joinpath(@__DIR__, "..", "scripts", "score_releases.jl"))

    grid_date(day) = Date(2026, 1, 1) + Day(day)
    made, target = grid_date(17), grid_date(24)
    path = joinpath(mktempdir(), "province_forecast.csv")
    open(path, "w") do io
        println(io, "made_date,horizon,target_date,province,stream,draw,value")
        for province in ("ituri", "other")
            for (d, v) in enumerate(8.0:12.0)
                row = (made, 7, target, province, "confirmed cases", d, v)
                println(io, join(row, ','))
            end
        end
    end

    obs = (;
        cutoff = grid_date(30),
        province_confirmed_history = _province_fixture(),
    )
    result = score_release("results-1", path, obs, grid_date)

    ## Each patch's rows are scored as their own stream, against that patch's
    ## own truth and its own persistence baseline.
    @test sort(unique(r.stream for r in result.rows)) ==
        ["confirmed cases [ituri]", "confirmed cases [other]"]
    @test sort(unique(r.fit for r in result.rows)) == ["baseline", "joint"]
    @test result.spans_break == 0

    ## With a break day inside the window, both patches' groups are dropped
    ## and counted for the release log.
    broken = merge(
        obs,
        (; confirmed_break_days = [20], confirmed_break_gross_cases = [10])
    )
    skipped = score_release("results-1", path, broken, grid_date)
    @test isempty(skipped.rows)
    @test skipped.spans_break == 2
end

@testitem "a release without the province archive is skipped" setup = [
    ProvinceFixture,
] begin
    using Dates: Date, Day

    include(joinpath(@__DIR__, "..", "scripts", "score_releases.jl"))

    grid_date(day) = Date(2026, 1, 1) + Day(day)
    obs = (;
        cutoff = grid_date(30),
        province_confirmed_history = _province_fixture(),
    )

    ## `fetch_asset` returns `nothing` for a release carrying no such asset,
    ## which every release published before the archive existed does.
    @test isnothing(
        score_province_release("results-1", nothing, obs, grid_date)
    )
end

@testitem "the province scores read projection rows only" setup = [
    ProvinceFixture,
] begin
    using Dates: Date, Day

    include(joinpath(@__DIR__, "..", "scripts", "score_releases.jl"))

    grid_date(day) = Date(2026, 1, 1) + Day(day)
    made, target = grid_date(17), grid_date(24)
    obs = (;
        cutoff = grid_date(30),
        province_confirmed_history = _province_fixture(),
    )
    function archive(; method)
        path = joinpath(mktempdir(), "province_forecast.csv")
        open(path, "w") do io
            println(
                io, "made_date,horizon,target_date,province,stream,draw,value" *
                    (isnothing(method) ? "" : ",method")
            )
            for province in ("ituri", "other"), (d, v) in enumerate(8.0:12.0)
                row = (made, 7, target, province, "confirmed cases", d, v)
                m = isnothing(method) ? () :
                    (province == "ituri" ? method : "share",)
                println(io, join((row..., m...), ','))
            end
        end
        return path
    end

    ## Only rows the per-province projection made are scored.
    result = score_province_release(
        "results-1", archive(; method = "projection"), obs, grid_date
    )
    @test unique(r.stream for r in result.rows) == ["confirmed cases [ituri]"]
    @test unique(r.stream for r in result.overlay) ==
        ["confirmed cases [ituri]"]

    ## An archive with no method column predates the projection, so the
    ## whole release is skipped rather than scored as a share split.
    @test score_province_release(
        "results-1", archive(; method = nothing), obs, grid_date
    ) === :no_projection
end
