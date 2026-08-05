# Execute one report page via Literate to produce docs/src/<page>.md. Select
# the page with the `BVD_DOC_PAGE` environment variable (`analysis` or
# `sensitivity`; default `analysis`). Every model fit is loaded from the
# content-addressed cache (`BVD_FIT_CACHE`) rather than refitted, so this is
# fast once the per-fit matrix (or an earlier run) has populated the cache.
#
#   BVD_DOC_PAGE=sensitivity julia --project=docs docs/execute.jl
#
# This is the execute step only: it runs the literate page (fits loaded from
# cache) and writes the markdown plus its figures and its half of the shared
# `output/` and `docs/src/summary_assets/`. The Vitepress render and deploy
# (docs/make.jl, which needs Node) is a separate combine step.
using Pkg: Pkg
Pkg.instantiate()

using Literate
using BVDOutbreakSize

const PAGE = String(strip(get(ENV, "BVD_DOC_PAGE", "analysis")))
PAGE in ("analysis", "sensitivity") ||
    error("BVD_DOC_PAGE must be \"analysis\" or \"sensitivity\", got \"$PAGE\"")

const LITERATE_SRC = joinpath(@__DIR__, "examples", "$PAGE.jl")
const LITERATE_OUT = joinpath(@__DIR__, "src")
isdir(LITERATE_OUT) || mkpath(LITERATE_OUT)

@info "Executing $PAGE.jl (fits are loaded from BVD_FIT_CACHE)…" cache = get(
    ENV, "BVD_FIT_CACHE", "logs/fit_cache")
Literate.markdown(
    LITERATE_SRC, LITERATE_OUT;
    name = PAGE,
    flavor = Literate.DocumenterFlavor(),
    execute = true,
    credit = false)
@info "Done" output = joinpath(LITERATE_OUT, "$PAGE.md")
