# Execute the analysis example via Literate to produce docs/src/analysis.md.
# Every model fit in the report is loaded from the content-addressed cache
# (BVD_FIT_CACHE) rather than refitted, so this is fast once the per-fit matrix
# (or an earlier run) has populated the cache. This is the EXECUTE step only;
# the Vitepress render (docs/make.jl, which needs Node) is separate, and
# nothing is deployed here.
using Pkg: Pkg
Pkg.instantiate()

using Literate
using BVDOutbreakSize

const LITERATE_SRC = joinpath(@__DIR__, "examples", "analysis.jl")
const LITERATE_OUT = joinpath(@__DIR__, "src")
isdir(LITERATE_OUT) || mkpath(LITERATE_OUT)

@info "Executing analysis.jl (fits are loaded from BVD_FIT_CACHE)…" cache =
    get(ENV, "BVD_FIT_CACHE", "logs/fit_cache")
Literate.markdown(
    LITERATE_SRC, LITERATE_OUT;
    name = "analysis",
    flavor = Literate.DocumenterFlavor(),
    execute = true,
    credit = false)
@info "Done" output = joinpath(LITERATE_OUT, "analysis.md")
