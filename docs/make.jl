using Pkg: Pkg
Pkg.instantiate()

using Documenter
using DocumenterCitations
using DocumenterVitepress
using Literate
using BVDOutbreakSize
import Dates

const bib = CitationBibliography(
    joinpath(@__DIR__, "src", "refs.bib");
    style = :authoryear
)

const REPO_ROOT = dirname(@__DIR__)
const LITERATE_SRC = joinpath(@__DIR__, "examples", "analysis.jl")
const LITERATE_OUT = joinpath(@__DIR__, "src")

isdir(LITERATE_OUT) || mkpath(LITERATE_OUT)

# Copy the README to the home page, stripping the SHARED-block marker
# comments. The analysis page reads them from the source README to load
# the shared prose, but they must not appear on the rendered home page
# (the Vitepress typographer mangles the `--` and shows them as text).
#
# The README links to analysis-page sections with absolute hosted URLs
# so they work when read on GitHub. On the rendered home page those would
# pin to a fixed version (/stable/); rewrite them to Documenter `@ref`
# cross-references so they instead resolve within whichever version is
# being viewed. The link target is the section anchor, whose Documenter
# slug is the header title with spaces replaced by dashes, so reversing
# that recovers the title for `@ref`.
let readme = read(joinpath(REPO_ROOT, "README.md"), String)
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
        r"\(https?://[^)]*?/analysis#([^)]+)\)" =>
            m -> begin
                slug = match(r"#([^)]+)\)$", m).captures[1]
                "(@ref \"" * replace(slug, '-' => ' ') * "\")"
            end
    )
    write(joinpath(LITERATE_OUT, "index.md"), readme)
end

Literate.markdown(
    LITERATE_SRC,
    LITERATE_OUT;
    name = "analysis",
    flavor = Literate.DocumenterFlavor(),
    execute = true,
    credit = false
)

# References page sourced from refs.bib through `@bibliography`.
open(joinpath(LITERATE_OUT, "references.md"), "w") do io
    println(io, "# References")
    println(io)
    println(io, "```@bibliography")
    println(io, "```")
end

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
        "Analysis" => "analysis.md",
        "API" => "api.md",
        "Contributing" => "contributing.md",
        "News" => "news.md",
        "References" => "references.md"
    ],
    format = DocumenterVitepress.MarkdownVitepress(;
        repo = "github.com/epiforecasts/BVDOutbreakSize",
        devbranch = "main",
        devurl = "dev",
        # Keep a docs version per minor release (v1.2, v1.3, …) in the
        # version dropdown rather than only the major alias (v1), which
        # is DocumenterVitepress's `:breaking` default.
        keep = :minor
    )
)

# Use DocumenterVitepress.deploydocs, not the bare Documenter one:
# DocumenterVitepress 0.2 builds into numbered subfolders
# (docs/build/1/, …) and its deploydocs flattens each build/i/ to
# gh-pages/<base>/. Plain deploydocs leaves the numbered subdir, so
# the deployed site's asset URLs 404. Ref LuxDL/DocumenterVitepress.jl#280.
DocumenterVitepress.deploydocs(;
    repo = "github.com/epiforecasts/BVDOutbreakSize",
    target = "build",
    branch = "gh-pages",
    devbranch = "main",
    push_preview = true
)
