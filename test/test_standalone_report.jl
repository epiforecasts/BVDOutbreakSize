## Tests for the page lookup in scripts/standalone_report.jl. The script's
## entry point is guarded by a PROGRAM_FILE check, so including it here
## defines the helpers without building a report.

@testitem "standalone report resolves a page by path, not basename" begin
    using BVDOutbreakSize: BVDOutbreakSize
    include(
        joinpath(
            pkgdir(BVDOutbreakSize), "scripts",
            "standalone_report.jl"
        )
    )

    ## The rendered site as it actually is: two pages render to
    ## `national.html`, under `estimates/` and under `forecasts/`.
    root = mktempdir()
    for (dir, name) in (
            ("estimates", "national.html"), ("estimates", "province.html"),
            ("forecasts", "national.html"), (".", "methods.html"),
            ("assets", "vp-icons.css"),
        )
        mkpath(joinpath(root, dir))
        write(joinpath(root, dir, name), "$dir/$name")
    end

    ## A qualified suffix picks the page it names, not whichever the walk
    ## reached first.
    estimates = find_one(root, joinpath("estimates", "national.html"))
    @test read(estimates, String) == "estimates/national.html"

    forecasts = find_one(root, joinpath("forecasts", "national.html"))
    @test read(forecasts, String) == "forecasts/national.html"

    ## A basename matching two pages is refused rather than resolved, so a
    ## new page colliding with an existing one fails the build instead of
    ## silently publishing the wrong one.
    @test_throws ErrorException find_one(root, "national.html")

    ## An unambiguous basename still resolves, wherever it sits.
    @test read(find_one(root, "methods.html"), String) == "./methods.html"
    @test endswith(find_one(root, "vp-icons.css"), "vp-icons.css")

    ## A missing file is still an error, and a partial component is not a
    ## match: `ational.html` must not resolve `national.html`.
    @test_throws ErrorException find_one(root, "absent.html")
    @test_throws ErrorException find_one(root, "ational.html")
end
