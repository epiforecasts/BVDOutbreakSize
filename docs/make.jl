using Pkg: Pkg
Pkg.instantiate()

using Documenter
using DocumenterCitations
using DocumenterVitepress
using Literate
using BVDOutbreakSize
import Dates

const REPO_ROOT = dirname(@__DIR__)
const EXAMPLES = joinpath(@__DIR__, "examples")
const LITERATE_OUT = joinpath(@__DIR__, "src")

## The report is split across two literate pages so the expensive fits and the
## render can fan out across CI runners: `analysis` (methods, results and the
## one-week-ahead forecast) and `sensitivity` (forecast validation and the
## comparison/sensitivity analyses). Both load the same cached fits through the
## shared `docs/examples/_setup.jl`.
const PAGES = ["analysis", "sensitivity"]

## Build stage, so fitting and rendering can be split across jobs:
##   render-main         → Literate-execute analysis.jl → docs/src/analysis.md
##   render-sensitivity  → Literate-execute sensitivity.jl → sensitivity.md
##   combine             → assemble the Vitepress site from the pre-rendered
##                         markdown (no execution) and deploy
##   all (default)       → render both pages then combine, for local builds
const STAGE = get(ENV, "BVD_DOCS_STAGE", "all")

isdir(LITERATE_OUT) || mkpath(LITERATE_OUT)

## Literate-execute one page. `execute = true` inlines every output, so the
## combine step assembles the site without re-running any code.
function render_page(page)
    @info "Literate render" page
    Literate.markdown(
        joinpath(EXAMPLES, "$page.jl"), LITERATE_OUT;
        name = page,
        flavor = Literate.DocumenterFlavor(),
        execute = true,
        credit = false)
end

## Copy the README to the home page, stripping the SHARED-block marker
## comments. The analysis page reads them from the source README to load
## the shared prose, but they must not appear on the rendered home page
## (the Vitepress typographer mangles the `--` and shows them as text).
##
## The README links to analysis-page sections with absolute hosted URLs
## so they work when read on GitHub. On the rendered home page those would
## pin to a fixed version (/stable/); rewrite them to Documenter `@ref`
## cross-references so they instead resolve within whichever version is
## being viewed. `@ref` resolves a section by title across all pages, so a
## link still resolves after a section moves from the analysis page to the
## sensitivity page. The link target is the section anchor, whose Documenter
## slug is the header title with spaces replaced by dashes, so reversing
## that recovers the title for `@ref`.
function write_index()
    readme = read(joinpath(REPO_ROOT, "README.md"), String)
    readme = replace(readme, r"^<!-- SHARED:END -->\n"m => "")
    ## Fill the live dates on the home page the same way the analysis page
    ## does: "Last updated" is the build date and "Data as of" is the loaded
    ## data cut-off, so a rebuild refreshes them without editing README.md.
    built = Dates.format(Dates.today(), "d U yyyy")
    asof = Dates.format(load_observations().cutoff, "d U yyyy")
    readme = replace(readme,
        r"\*\*Last updated:\*\* [^.]*\." => "**Last updated:** $built.",
        r"\*\*Data as of:\*\* [^.]*\." => "**Data as of:** $asof.")
    readme = replace(
        readme,
        r"\(https?://[^)]*?/(analysis|sensitivity)#([^)]+)\)" =>
            m -> begin
                slug = match(r"#([^)]+)\)$", m).captures[1]
                "(@ref \"" * replace(slug, '-' => ' ') * "\")"
            end
    )
    write(joinpath(LITERATE_OUT, "index.md"), readme)
end

## References page sourced from refs.bib through `@bibliography`.
function write_references()
    open(joinpath(LITERATE_OUT, "references.md"), "w") do io
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
        joinpath(@__DIR__, "src", "refs.bib"); style = :authoryear)
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
            "Summary" => "summary.md",
            "Analysis" => "analysis.md",
            "Sensitivity" => "sensitivity.md",
            "API" => "api.md",
            "Contributing" => "contributing.md",
            "News" => "news.md",
            "References" => "references.md"
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
    DocumenterVitepress.deploydocs(;
        repo = "github.com/epiforecasts/BVDOutbreakSize",
        target = "build",
        branch = "gh-pages",
        devbranch = "main",
        push_preview = true
    )
end

if STAGE == "render-main"
    render_page("analysis")
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
    error("unknown BVD_DOCS_STAGE=$STAGE; expected one of render-main, " *
          "render-sensitivity, combine, all")
end
