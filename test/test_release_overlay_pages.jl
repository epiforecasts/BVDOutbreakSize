## The estimates across past releases are checks on the national fit, so they
## sit on the national in-sample page (#1058). The page reads the release
## overlays, so its CI render job must collect the rescored ones rather than
## fall back to the committed copies.

@testitem "report pages: past-release checks on the in-sample page" tags = [
    :quality,
] begin
    using BVDOutbreakSize
    pages = joinpath(pkgdir(BVDOutbreakSize), "docs", "pages")
    headings = [
        "# ## Estimate evolution across releases",
        "# ### Reproduction number by release",
        "# ### Basic reproduction number by release",
    ]
    hosts = Dict(h => String[] for h in headings)
    for (dir, _, files) in walkdir(pages), f in files
        endswith(f, ".jl") || continue
        path = joinpath(dir, f)
        lines = rstrip.(readlines(path))
        for h in headings
            h in lines && push!(hosts[h], relpath(path, pages))
        end
    end
    insample = joinpath("evaluation", "insample", "national.jl")
    for h in headings
        @test hosts[h] == [insample]
    end
end

@testitem "report pages: release-overlay readers collect the overlays" tags = [
    :quality,
] begin
    using BVDOutbreakSize
    root = pkgdir(BVDOutbreakSize)
    pages = joinpath(root, "docs", "pages")
    ## A code line that reads a release overlay table: through the typed
    ## fallback, by one of the per-release file names, or the released
    ## estimates loaded in the shared setup.
    reads(l) = !occursin(r"^\h*#", l) && (
        occursin("_release_data(", l) ||
            occursin(r"by_release[a-z_]*\.csv", l) ||
            occursin(r"\breleased_df\b", l)
    )
    readers = String[]
    for (dir, _, files) in walkdir(pages), f in files
        endswith(f, ".jl") && f != "_setup.jl" || continue
        path = joinpath(dir, f)
        any(reads, eachline(path)) &&
            push!(readers, splitext(relpath(path, pages))[1])
    end
    @test joinpath("evaluation", "insample", "national") in readers
    ## Each render matrix entry's `id` and its `overlays` flag, in order.
    workflow = read(joinpath(root, ".github", "workflows", "docs.yml"), String)
    flags = Dict(
        m[1] => m[2] == "true" for m in eachmatch(
                r"- id: (\S+)\n(?:\h+\w+:.*\n)*?\h+overlays: (true|false)",
                workflow
            )
    )
    for page in readers
        id = replace(page, '\\' => '/')
        @test get(flags, id, false)
    end
end

## Past forecasts against what was observed are scored on the forecast
## evaluation pages only, so no forecast page reads the archived forecast
## overlays or carries its own copy of the figure.
## A source-text lint: it reads the page files, it does not render them.
@testitem "report pages: past forecasts sit on the evaluation pages" tags = [
    :quality,
] begin
    using BVDOutbreakSize
    pages = joinpath(pkgdir(BVDOutbreakSize), "docs", "pages")
    code(path) = filter(l -> !occursin(r"^\h*#", l), readlines(path))
    readers = String[]
    plotters = String[]
    for (dir, _, files) in walkdir(pages), f in files
        endswith(f, ".jl") || continue
        path = joinpath(dir, f)
        page = replace(relpath(path, pages), '\\' => '/')
        lines = code(path)
        any(contains("forecast_overlay.csv"), lines) && push!(readers, page)
        any(contains("plot_forecast_overlay("), lines) &&
            push!(plotters, page)
    end
    @test !isempty(readers)
    @test all(startswith("evaluation/forecast/"), readers)
    @test !isempty(plotters)
    @test all(startswith("evaluation/forecast/"), plotters)
    for f in ("province.jl", "zone.jl")
        @test !any(
            contains("plot_forecast_overlay("),
            code(joinpath(pages, "forecasts", f))
        )
    end
end
