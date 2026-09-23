using Pkg: Pkg
Pkg.instantiate()

using Documenter
using DocumenterCitations
using DocumenterVitepress
using Literate
using BVDOutbreakSize

const REPO_ROOT = dirname(@__DIR__)
const PAGES_DIR = joinpath(@__DIR__, "pages")
const LITERATE_OUT = joinpath(@__DIR__, "src")

## The report is split across literate pages so the expensive fits and the
## render can fan out across CI runners. `methods` carries the data, the
## model and how it is fitted, `estimates/national` the national results,
## `estimates/province` the per-province estimates, `forecasts/national` the
## one-week-ahead projections, `evaluation/national` and `evaluation/province`
## the in-sample checks and the scoring against what arrived at each level,
## and `sensitivity` the comparison and sensitivity analyses. All load the
## same cached fits through the shared `docs/pages/_setup.jl`.
const PAGES = [
    "methods",
    "estimates/national", "estimates/province",
    "forecasts/national",
    "evaluation/national", "evaluation/province",
    "sensitivity",
]

## Build stage, so fitting and rendering can be split across jobs:
##   render-methods      → methods.jl → src/methods.md
##   render-main         → estimates/national.jl → src/estimates/national.md
##   render-province     → estimates/province.jl
##   render-forecast     → forecasts/national.jl
##   render-evaluation   → evaluation/national.jl
##   render-evaluation-province → evaluation/province.jl
##   render-sensitivity  → sensitivity.jl
##   combine             → assemble the Vitepress site from the pre-rendered
##                         markdown (no execution) and deploy
##   all (default)       → render both pages then combine, for local builds
const STAGE = get(ENV, "BVD_DOCS_STAGE", "all")

isdir(LITERATE_OUT) || mkpath(LITERATE_OUT)

## Literate-execute one page. `execute = true` inlines every output, so the
## combine step assembles the site without re-running any code.
function render_page(page)
    @info "Literate render" page
    ## A page id is its path under `docs/pages`, so a grouped page renders
    ## into the matching folder under `docs/src` and keeps its own file name.
    out = joinpath(LITERATE_OUT, dirname(page))
    isdir(out) || mkpath(out)
    return Literate.markdown(
        joinpath(PAGES_DIR, "$page.jl"), out;
        name = basename(page),
        flavor = Literate.DocumenterFlavor(),
        execute = true,
        credit = false
    )
end

include(joinpath(@__DIR__, "front_matter.jl"))

## Copy the README to the home page, with its live dates filled in and the
## SHARED-block marker comment stripped: it must not appear on the rendered
## home page (the Vitepress typographer mangles the `--` and shows it as
## text).
##
## The README links into the hosted report with absolute URLs so they work
## when read on GitHub. On the rendered home page those would pin to a fixed
## version (/stable/); rewrite them so they instead resolve within whichever
## version is being viewed. A link to a section becomes a Documenter `@ref`,
## which resolves by section title across every page and so survives a
## section moving page; the anchor's Documenter slug is the title with spaces
## replaced by dashes, so reversing that recovers the title. A link to a whole
## page becomes a relative link to that page's markdown. Only links to the
## documentation host are rewritten, so the badges and the repository links
## are left alone.
function write_index()
    readme = replace(readme_with_dates(), r"^<!-- SHARED:END -->\n"m => "")
    docs_url = r"\(https?://epiforecasts\.io/BVDOutbreakSize/[^)/]+/"
    readme = replace(
        readme,
        Regex(docs_url.pattern * "[a-z_/-]+#([^)]+)\\)") =>
            m -> begin
            slug = match(r"#([^)]+)\)$", m).captures[1]
            "(@ref \"" * replace(slug, '-' => ' ') * "\")"
        end,
        Regex(docs_url.pattern * "((?:[a-z_-]+/)*[a-z_-]+)\\)") =>
            m -> "(" * match(
            r"/BVDOutbreakSize/[^)/]+/((?:[a-z_-]+/)*[a-z_-]+)\)$", m
        ).captures[1] * ".md)"
    )
    return write(joinpath(LITERATE_OUT, "index.md"), readme)
end

## References page sourced from refs.bib through `@bibliography`.
function write_references()
    return open(joinpath(LITERATE_OUT, "references.md"), "w") do io
        println(io, "# References")
        println(io)
        println(io, "```@bibliography")
        println(io, "```")
    end
end

## Assemble and deploy the Vitepress site from the pre-rendered markdown. The
## two report pages are already executed (Literate `execute = true`), so
## makedocs does not re-run them; it resolves `@ref`/`@cite`/`@bibliography`
## across all pages in a single pass, so cross-page links resolve.
function combine()
    bib = CitationBibliography(
        joinpath(@__DIR__, "src", "refs.bib"); style = :authoryear
    )
    write_index()
    write_references()
    makedocs(;
        sitename = "BVDOutbreakSize",
        authors = "Sam Abbott and contributors",
        repo = "github.com/epiforecasts/BVDOutbreakSize",
        clean = true,
        doctest = false,
        warnonly = [:missing_docs, :linkcheck, :citations],
        plugins = [bib],
        pages = [
            "Home" => "index.md",
            "Estimates" => [
                "Summary" => "estimates/summary.md",
                "National" => "estimates/national.md",
                "Provinces" => "estimates/province.md",
            ],
            "Forecasts" => "forecasts/national.md",
            "Evaluation" => [
                "National" => "evaluation/national.md",
                "Provinces" => "evaluation/province.md",
            ],
            "Details" => [
                "Aim and origins" => "aim.md",
                "Methods" => "methods.md",
                "Limitations" => "limitations.md",
                "Sensitivity" => "sensitivity.md",
            ],
            "API" => [
                "Overview" => "lib/api.md",
                "Data and constants" => "lib/data.md",
                "Renewal and delays" => "lib/renewal.md",
                "Priors and latent submodels" => "lib/priors.md",
                "Observation models" => "lib/observations.md",
                "Joint and single-stream models" => "lib/joint.md",
                "Fitting" => "lib/fitting.md",
                "Summaries and diagnostics" => "lib/summaries.md",
                "Forecasts and scoring" => "lib/forecasts.md",
                "Plotting" => "lib/plotting.md",
                "Internals" => "lib/internals.md",
            ],
            "About" => [
                "Authors and funding" => "about.md",
                "Contributing" => "contributing.md",
                "News" => "news.md",
                "References" => "references.md",
            ],
        ],
        format = DocumenterVitepress.MarkdownVitepress(;
            repo = "github.com/epiforecasts/BVDOutbreakSize",
            devbranch = "main",
            devurl = "dev",
            ## Keep a docs version per minor release (v1.2, v1.3, …) in the
            ## version dropdown rather than only the major alias (v1), which
            ## is DocumenterVitepress's `:breaking` default.
            keep = :minor
        )
    )

    ## Use DocumenterVitepress.deploydocs, not the bare Documenter one:
    ## DocumenterVitepress 0.2 builds into numbered subfolders
    ## (docs/build/1/, …) and its deploydocs flattens each build/i/ to
    ## gh-pages/<base>/. Plain deploydocs leaves the numbered subdir, so
    ## the deployed site's asset URLs 404. Ref LuxDL/DocumenterVitepress.jl#280.
    return DocumenterVitepress.deploydocs(;
        repo = "github.com/epiforecasts/BVDOutbreakSize",
        target = "build",
        branch = "gh-pages",
        devbranch = "main",
        push_preview = true
    )
end

if STAGE == "render-methods"
    render_page("methods")
elseif STAGE == "render-main"
    render_page("estimates/national")
elseif STAGE == "render-province"
    render_page("estimates/province")
elseif STAGE == "render-forecast"
    render_page("forecasts/national")
elseif STAGE == "render-evaluation"
    render_page("evaluation/national")
elseif STAGE == "render-evaluation-province"
    render_page("evaluation/province")
elseif STAGE == "render-sensitivity"
    render_page("sensitivity")
elseif STAGE == "combine"
    combine()
elseif STAGE == "all"
    for page in PAGES
        render_page(page)
    end
    combine()
else
    error(
        "unknown BVD_DOCS_STAGE=$STAGE; expected one of render-methods, " *
            "render-main, " *
            "render-province, render-forecast, " *
            "render-evaluation, render-evaluation-province, " *
            "render-sensitivity, combine, all"
    )
end
