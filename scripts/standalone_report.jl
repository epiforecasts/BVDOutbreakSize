# Build a self-contained single-file HTML copy of the rendered report: the
# front matter from `README.md`, then the methods and national estimates
# pages.
#
# The Vitepress build statically pre-renders the full analysis content
# (tables, math, images as base64 data URIs) into the page HTML. This
# script lifts that content out of
# the multi-file Vitepress site, inlines the stylesheet and the Inter
# web fonts, and drops the SPA JavaScript so the result is one HTML
# file that opens offline. It is published as a release asset by the
# docs workflow.
#
# Usage:
#   julia --project=docs scripts/standalone_report.jl [build_dir] [out_file]
#   build_dir defaults to docs/build, out_file to output/analysis.html.

using Base64: base64encode
using Markdown: Markdown
include(joinpath(@__DIR__, "..", "docs", "front_matter.jl"))

# Walk forward from the marker to the matching </div>, returning the
# whole balanced <div>…</div> that contains it.
function extract_balanced_div(html::AbstractString, marker::AbstractString)
    hit = findfirst(marker, html)
    hit === nothing && error("content marker not found: $marker")
    op = findprev("<div", html, first(hit))
    op === nothing && error("no opening <div before marker")
    i = first(op)
    depth = 0
    j = i
    n = lastindex(html)
    while j <= n
        if startswith(SubString(html, j), "<div")
            depth += 1
            j = nextind(html, j, 4)
        elseif startswith(SubString(html, j), "</div>")
            depth -= 1
            j = nextind(html, j, 6)
            depth == 0 && return SubString(html, i, prevind(html, j))
        else
            j = nextind(html, j)
        end
    end
    error("unbalanced <div> while extracting content")
end

# Replace woff2 url(...) references with base64 data URIs, reading the
# font files by basename from the same assets directory as the CSS. A
# referenced font that is not found is left untouched.
function embed_fonts(css::AbstractString, assets_dir::AbstractString)
    return replace(
        css,
        r"url\(([^)]*?([^/)]+\.woff2))\)" => function (m)
            name = match(r"url\([^)]*?([^/)]+\.woff2)\)", m).captures[1]
            path = joinpath(assets_dir, name)
            isfile(path) || return m
            data = base64encode(read(path))
            return "url(data:font/woff2;base64,$data)"
        end
    )
end

## Locate the one file under `root` whose trailing path components are
## `suffix`. The suffix is matched component by component rather than as a
## basename, and more than one match is an error rather than a choice,
## because two pages render to `national.html` (the national estimates and
## the forecasts). Taking whichever `walkdir` reached first would publish
## the wrong page as `analysis.html`, and would do it silently, since both
## are valid pages that lift cleanly.
function find_one(root::AbstractString, suffix::AbstractString)
    want = splitpath(suffix)
    hits = String[]
    for (dir, _, files) in walkdir(root)
        for f in files
            parts = splitpath(joinpath(dir, f))
            length(parts) >= length(want) &&
                parts[(end - length(want) + 1):end] == want &&
                push!(hits, joinpath(dir, f))
        end
    end
    isempty(hits) && error("$suffix not found under $root")
    length(hits) == 1 || error(
        "$suffix matches more than one file under $root: " *
            join(sort(hits), ", ")
    )
    return only(hits)
end

## The Vitepress `assets/` directory. It sits once at the site root, beside
## `vp-icons.css`, rather than beside each page.
function site_assets(build_dir::AbstractString)
    return joinpath(dirname(find_one(build_dir, "vp-icons.css")), "assets")
end

function build_standalone(build_dir::AbstractString, out_file::AbstractString)
    ## The methods and the national estimates, in that order, so the offline
    ## copy carries the model as well as the results. The published asset
    ## keeps the name `analysis.html` so the release download link does not
    ## move.
    page = find_one(build_dir, joinpath("estimates", "national.html"))
    methods_page = find_one(build_dir, "methods.html")
    assets = site_assets(build_dir)
    html = read(page, String)

    ## The front matter (title, authors, dates, abstract, scope) is not on
    ## any rendered page other than the home page, so render it from the
    ## README here. Its first line is the report title.
    front = front_matter()
    heading = match(r"^# (.*)$"m, front)
    heading === nothing &&
        error("front matter has no top-level heading to use as the title")
    title = heading.captures[1]

    content = join(
        (
            Markdown.html(Markdown.parse(front)),
            extract_balanced_div(read(methods_page, String), "vp-doc _"),
            extract_balanced_div(html, "vp-doc _"),
        ), "\n"
    )
    # Rewrite cross-references to either of the two pages in this file as
    # bare anchors, so the in-page jump links work in the standalone file
    # rather than navigating to the hosted site. The front matter links with
    # absolute URLs, so the host is optional.
    in_file = r"(?:https://epiforecasts\.io)?/BVDOutbreakSize/[^\"#]*" *
        r"(?:estimates/national|methods)#"
    content = replace(content, in_file => "#")
    # Any remaining root-relative site links point at other doc pages
    # (citations resolve to the references page, etc.) that are not part
    # of this single file. Absolutise them against the hosted site so
    # they resolve instead of breaking against the local filesystem.
    content = replace(content, "href=\"/BVDOutbreakSize/" => "href=\"https://epiforecasts.io/BVDOutbreakSize/")

    css = ""
    for f in readdir(assets)
        startswith(f, "style.") && endswith(f, ".css") || continue
        css *= embed_fonts(read(joinpath(assets, f), String), assets)
    end
    icons_css = read(find_one(build_dir, "vp-icons.css"), String)

    doc = """
    <!doctype html>
    <html lang="en">
    <head>
    <meta charset="utf-8">
    <meta name="viewport" content="width=device-width, initial-scale=1">
    <title>$title</title>
    <style>
    $css
    $icons_css
    body{max-width:900px;margin:2rem auto;padding:0 1rem;}
    </style>
    </head>
    <body>
    <div class="vp-doc">
    $content
    </div>
    </body>
    </html>
    """

    mkpath(dirname(out_file))
    write(out_file, doc)
    return out_file
end

if abspath(PROGRAM_FILE) == @__FILE__
    root = joinpath(@__DIR__, "..")
    build_dir = length(ARGS) >= 1 ? ARGS[1] :
        joinpath(root, "docs", "build")
    out_file = length(ARGS) >= 2 ? ARGS[2] :
        joinpath(root, "output", "analysis.html")
    out = build_standalone(build_dir, out_file)
    println(
        "wrote self-contained report: ", out,
        " (", filesize(out), " bytes)"
    )
end
