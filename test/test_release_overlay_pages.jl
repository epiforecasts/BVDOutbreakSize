## The estimates across past releases are checks on the fit, so they sit on
## the in-sample pages (#1058): the national page, and the province and zone
## pages for their own areas. These pages read the release overlays, so their
## CI render jobs must collect the rescored ones rather than fall back to the
## committed copies.

@testitem "report pages: past-release checks on the in-sample page" tags = [
    :quality,
] begin
    using BVDOutbreakSize
    pages = joinpath(pkgdir(BVDOutbreakSize), "docs", "pages")
    insample = joinpath("evaluation", "insample")
    national = joinpath(insample, "national.jl")
    areas = [joinpath(insample, "province.jl"), joinpath(insample, "zone.jl")]
    ## Each heading and the in-sample pages that carry it.
    headings = Dict(
        "# ## Estimate evolution across releases" => [national; areas],
        "# ### Reproduction number by release" => [national; areas],
        "# ### Confirmed cases by release" => areas,
        "# ### Basic reproduction number by release" => [national],
    )
    hosts = Dict(h => String[] for h in keys(headings))
    for (dir, _, files) in walkdir(pages), f in files
        endswith(f, ".jl") || continue
        path = joinpath(dir, f)
        lines = rstrip.(readlines(path))
        for h in keys(headings)
            h in lines && push!(hosts[h], relpath(path, pages))
        end
    end
    for (h, expected) in headings
        @test sort(hosts[h]) == sort(expected)
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
    for page in ("national", "province", "zone")
        @test joinpath("evaluation", "insample", page) in readers
    end
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
