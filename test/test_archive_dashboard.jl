## Tests for the network-free parts of scripts/archive_dashboard.jl. The
## script's `main` is guarded by a PROGRAM_FILE check, so including it here
## defines the helpers without calling `gh`.

@testitem "archive_dashboard keeps the newest release per cut-off" begin
    using BVDOutbreakSize: BVDOutbreakSize
    using Dates
    include(
        joinpath(pkgdir(BVDOutbreakSize), "scripts", "archive_dashboard.jl")
    )

    t(s) = DateTime(s)
    entries = [
        ("results-2968", t("2026-10-08T14:24:22"), Date(2026, 10, 7)),
        ("results-2978", t("2026-10-08T14:24:22"), Date(2026, 10, 7)),
        ("results-3000", t("2026-10-12T01:33:38"), Date(2026, 10, 10)),
        ("results-2935", t("2026-10-06T23:40:52"), Date(2026, 10, 5)),
        ("results-v3.0.0", t("2026-09-30T21:20:48"), Date(2026, 9, 30)),
        ("results-2817", t("2026-09-30T21:20:48"), Date(2026, 9, 30)),
    ]
    @test pick_vintages(entries) == [
        (Date(2026, 10, 10), "results-3000"),
        (Date(2026, 10, 7), "results-2978"),
        (Date(2026, 10, 5), "results-2935"),
        (Date(2026, 9, 30), "results-v3.0.0"),
    ]
    @test pick_vintages(entries; limit = 2) == [
        (Date(2026, 10, 10), "results-3000"),
        (Date(2026, 10, 7), "results-2978"),
    ]
    @test isempty(pick_vintages(empty(entries)))

    assets = [collect(DASHBOARD_CSVS); "observations.toml"; "site.zip"]
    @test carries_dashboard(assets)
    @test !carries_dashboard(filter(!=("zone_timeseries.csv"), assets))
    @test !carries_dashboard(filter(!=("observations.toml"), assets))
end

@testitem "archive_dashboard writes one folder per cut-off and index.json" begin
    using BVDOutbreakSize: BVDOutbreakSize
    using Dates
    include(
        joinpath(pkgdir(BVDOutbreakSize), "scripts", "archive_dashboard.jl")
    )

    mktempdir() do dir
        source = joinpath(dir, "releases")
        for tag in ("results-3000", "results-2978")
            mkpath(joinpath(source, tag))
            for f in DASHBOARD_CSVS
                write(joinpath(source, tag, f), "as_of\n$tag\n")
            end
        end
        ## A release whose download lost a table is left out.
        mkpath(joinpath(source, "results-2935"))
        write(joinpath(source, "results-2935", "zone_estimates.csv"), "x\n")

        dest = joinpath(dir, "archive")
        mkpath(joinpath(dest, "2026-01-01"))
        vintages = [
            (Date(2026, 10, 10), "results-3000"),
            (Date(2026, 10, 7), "results-2978"),
            (Date(2026, 10, 5), "results-2935"),
        ]
        written = @test_logs (:warn, r"results-2935") write_archive(
            dest, vintages, source
        )
        @test written == vintages[1:2]
        @test sort(readdir(dest)) == ["2026-10-07", "2026-10-10", "index.json"]
        @test sort(readdir(joinpath(dest, "2026-10-10"))) ==
            sort(collect(DASHBOARD_CSVS))
        @test read(joinpath(dest, "2026-10-07", "national_estimates.csv"), String) ==
            "as_of\nresults-2978\n"
        @test read(joinpath(dest, "index.json"), String) == """
            [
              {"as_of": "2026-10-10", "tag": "results-3000"},
              {"as_of": "2026-10-07", "tag": "results-2978"}
            ]
            """
        @test index_json(empty(vintages)) == "[]\n"
    end
end
