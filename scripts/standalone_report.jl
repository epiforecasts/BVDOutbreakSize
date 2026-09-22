# Build a self-contained single-file HTML copy of the rendered analysis
# page.
#
# The Vitepress build statically pre-renders the full analysis content
# (tables, math, images as base64 data URIs) into the page HTML. This
# script lifts that content out of
# the multi-file Vitepress site, inlines the stylesheet and the Inter
# web fonts, and drops the SPA JavaScript so the result is one HTML
# file that opens offline. It is published as a release asset by the
# docs workflow.
#
# Usage: julia scripts/standalone_report.jl [build_dir] [out_file]
#   build_dir defaults to docs/build, out_file to output/analysis.html.

using Base64: base64encode

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

function find_one(root::AbstractString, name::AbstractString)
    for (dir, _, files) in walkdir(root)
        name in files && return joinpath(dir, name)
    end
    error("$name not found under $root")
end

function build_standalone(build_dir::AbstractString, out_file::AbstractString)
    ## The national estimates page. The published asset keeps the name
    ## `analysis.html` so the release download link does not move.
    page = find_one(build_dir, "national.html")
    assets = joinpath(dirname(page), "assets")
    html = read(page, String)

    title = let m = match(r"<title>(.*?)\s*\|", html)
        m === nothing ? "Analysis" : m.captures[1]
    end

    content = extract_balanced_div(html, "vp-doc _")
    # Rewrite same-page cross-references to bare anchors so the in-page
    # jump links work in the standalone file rather than navigating to
    # the hosted site.
    content = replace(content, r"/BVDOutbreakSize/[^\"#]*national#" => "#")
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
