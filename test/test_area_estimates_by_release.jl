## Tests for the per-province and per-zone estimate collection in
## scripts/score_releases.jl (area_estimate_rows, write_area_estimates): each
## release's area table becomes long rows keyed by release, area and
## quantity, a release without the table is skipped, and a renamed column
## still reads.

## Two fake releases. The older carries no province table and writes its
## zone reproduction number under the `rt_` stem. The newer carries both in
## the current schema, with one province pooled into a shared patch and one
## zone that has no reproduction number.
@testsnippet AreaReleaseFixture begin
    using Dates: Date
    function _area_fixture(dir)
        old = mkpath(joinpath(dir, "results-100"))
        write(
            joinpath(old, "zone_estimates.csv"),
            """
            zone,label,province,patch,cases,deaths,rt_median,rt_lower,rt_upper,as_of
            mambasa,Mambasa,ituri,Ituri,40,12,1.25,0.8,1.9,2026-09-26
            """
        )
        new = mkpath(joinpath(dir, "results-200"))
        write(
            joinpath(new, "province_estimates.csv"),
            """
            province,patch,pooled,cases,deaths,R_T_median,R_T_lower,R_T_upper,as_of,cfr_median,cfr_lower,cfr_upper
            ituri,Ituri,0,7040,3271,1.0668,0.7330,1.5293,2026-10-20,0.48,0.35,0.61
            bas_uele,Other provinces,1,13,4,1.1031,0.6695,1.7218,2026-10-20,0.49,0.37,0.61
            """
        )
        write(
            joinpath(new, "zone_estimates.csv"),
            """
            zone,label,province,patch,cases,deaths,R_T_median,p_rt_above_one,R_T_lower,R_T_upper,as_of
            mambasa,Mambasa,ituri,Ituri,55,15,0.97,0.45,0.58,1.46,2026-10-20
            adja,Adja,ituri,Ituri,11,1,,,,,2026-10-20
            """
        )
        return [
            ("results-100", old, Date(2026, 9, 26)),
            ("results-200", new, Date(2026, 10, 20)),
        ]
    end
    ## Every release's rows of `level`, skipping one without the table.
    function _collect(releases, level, asset)
        rows = NamedTuple[]
        used = String[]
        for (tag, d, cutoff) in releases
            p = joinpath(d, asset)
            r = area_estimate_rows(tag, isfile(p) ? p : nothing, cutoff, level)
            isnothing(r) && continue
            append!(rows, r)
            push!(used, tag)
        end
        return (; rows, used)
    end
end

@testitem "area estimates skip a release without the table" setup = [
    AreaReleaseFixture,
] begin
    using Dates: Date
    include(joinpath(@__DIR__, "..", "scripts", "score_releases.jl"))

    mktempdir() do dir
        releases = _area_fixture(dir)
        prov = _collect(releases, :province, PROVINCE_ESTIMATES_ASSET)
        @test prov.used == ["results-200"]
        ## Two provinces, four quantities each.
        @test length(prov.rows) == 8
        rt = only(
            r for r in prov.rows if r.area == "ituri" && r.quantity == "R_T"
        )
        @test rt.date == "2026-10-20"
        @test (rt.median, rt.lo90, rt.hi90) == (1.0668, 0.733, 1.5293)
        ## A count carries no interval.
        cases = only(
            r for r in prov.rows if r.area == "ituri" && r.quantity == "cases"
        )
        @test cases.median == 7040
        @test ismissing(cases.lo90) && ismissing(cases.hi90)
        ## The pooled flag and the patch ride along on every row.
        pooled = [r for r in prov.rows if r.area == "bas_uele"]
        @test all(r -> r.pooled == "1" && r.patch == "Other provinces", pooled)
        @test any(r -> r.quantity == "CFR", pooled)
    end
end

@testitem "area estimates read a renamed column and key zones by province" setup = [
    AreaReleaseFixture,
] begin
    using Dates: Date
    include(joinpath(@__DIR__, "..", "scripts", "score_releases.jl"))

    mktempdir() do dir
        releases = _area_fixture(dir)
        zone = _collect(releases, :zone, ZONE_ESTIMATES_ASSET)
        @test zone.used == ["results-100", "results-200"]
        ## The older release's `rt_` stem reads as the reproduction number.
        old_rt = only(
            r for r in zone.rows
                if r.release == "results-100" && r.quantity == "R_T"
        )
        @test old_rt.area == "ituri.mambasa"
        @test (old_rt.median, old_rt.lo90, old_rt.hi90) == (1.25, 0.8, 1.9)
        ## A zone with no reproduction number gives no R_T row, but still
        ## gives its counts.
        adja = [r for r in zone.rows if r.area == "ituri.adja"]
        @test sort([r.quantity for r in adja]) == ["cases", "deaths"]
        ## A table with no area column is not read.
        p = joinpath(dir, "bad.csv")
        write(p, "cases,deaths\n1,2\n")
        @test isnothing(area_estimate_rows("x", p, Date(2026, 1, 1), :zone))
    end
end

@testitem "area estimates write a long table, header only when empty" setup = [
    AreaReleaseFixture,
] begin
    using Dates: Date
    include(joinpath(@__DIR__, "..", "scripts", "score_releases.jl"))

    mktempdir() do dir
        releases = _area_fixture(dir)
        dest = joinpath(dir, "out", "province", "released_estimates.csv")
        write_area_estimates(dest, NamedTuple[], :province)
        @test readlines(dest) ==
            ["release,date,area,patch,pooled,quantity,median,lo90,hi90"]

        prov = _collect(releases, :province, PROVINCE_ESTIMATES_ASSET)
        write_area_estimates(dest, prov.rows, :province)
        lines = readlines(dest)
        @test length(lines) == 1 + length(prov.rows)
        @test "results-200,2026-10-20,ituri,Ituri,0,R_T,1.067,0.733,1.529" in
            lines
        @test "results-200,2026-10-20,ituri,Ituri,0,cases,7040,," in lines

        zone = _collect(releases, :zone, ZONE_ESTIMATES_ASSET)
        zdest = joinpath(dir, "out", "zone", "released_estimates.csv")
        write_area_estimates(zdest, zone.rows, :zone)
        zlines = readlines(zdest)
        @test zlines[1] ==
            "release,date,area,label,province,patch,quantity,median,lo90,hi90"
        ## Sorted by date, so the older release comes first.
        @test startswith(zlines[2], "results-100,2026-09-26,ituri.mambasa")
    end
end
