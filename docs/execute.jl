# Execute one report page via Literate to produce docs/src/<page>.md. Select
# the page with the `BVD_DOC_PAGE` environment variable (one of `methods`,
# `estimates/national`, `estimates/province`, `forecasts/national`,
# `evaluation/national`, `evaluation/province`, `sensitivity`). Every model
# fit is loaded from the content-addressed cache (`BVD_FIT_CACHE`) rather
# than refitted, so this is fast once the per-fit matrix (or an earlier run)
# has populated the cache.
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
PAGE in (
    "methods",
    "estimates/national", "estimates/province", "forecasts/national",
    "evaluation/national", "evaluation/province", "sensitivity",
) ||
    error(
    "BVD_DOC_PAGE must be one of methods, estimates/national, " *
        "estimates/province, forecasts/national, evaluation/national, " *
        "evaluation/province, sensitivity; got \"$PAGE\""
)

const LITERATE_SRC = joinpath(@__DIR__, "pages", "$PAGE.jl")
const LITERATE_OUT = joinpath(@__DIR__, "src")
isdir(LITERATE_OUT) || mkpath(LITERATE_OUT)

@info "Executing $PAGE.jl (fits are loaded from BVD_FIT_CACHE)…" cache = get(
    ENV, "BVD_FIT_CACHE", "logs/fit_cache"
)
const PAGE_OUT = joinpath(LITERATE_OUT, dirname(PAGE))
isdir(PAGE_OUT) || mkpath(PAGE_OUT)
Literate.markdown(
    LITERATE_SRC, PAGE_OUT;
    name = basename(PAGE),
    flavor = Literate.DocumenterFlavor(),
    execute = true,
    credit = false
)
@info "Done" output = joinpath(LITERATE_OUT, "$PAGE.md")
