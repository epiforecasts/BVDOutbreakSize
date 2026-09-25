# Run and cache a single fit from the registry, for the per-fit CI matrix (or
# an HPC task). Select the fit with the `BVD_FIT_ID` environment variable; the
# result is serialised into `BVD_FIT_CACHE` under its content-addressed key, so
# the docs build can load it instead of refitting.
#
#   BVD_FIT_ID=confirmed julia --project=docs docs/fits/one.jl
#
# Set `BVD_FIT_CACHE` to choose the cache directory (default `logs/fit_cache`),
# `BVD_REFIT=all` to ignore an existing cache entry, and `BVD_FIT_DRYRUN=1` to
# print the id and key without fitting (a cheap check that the id resolves).
#
# A dependent fit (the health-zone `local` and `local_frozen_validation`)
# loads its parent strictly from the same cache: fit `joint` or
# `frozen_validation` first.
#
# The fit's convergence diagnostics and headline posteriors are printed, and
# appended to the GitHub Actions job summary when `GITHUB_STEP_SUMMARY` is set.
# A diagnostics bundle (`<key>_diagnostics/*.csv`) and, for a joint chain, the
# parent extract the zone stage reads (`<key>.parent.jls`) are written next to
# the cached fit.
using Pkg: Pkg
Pkg.instantiate()

using BVDOutbreakSize
include(joinpath(@__DIR__, "registry.jl"))
include(joinpath(@__DIR__, "summary.jl"))

const ID = strip(get(ENV, "BVD_FIT_ID", ""))
const CACHE = fit_cache_dir()
const REFIT = lowercase(strip(get(ENV, "BVD_REFIT", ""))) in ("all", "true", "1")
const DRYRUN = lowercase(strip(get(ENV, "BVD_FIT_DRYRUN", ""))) in ("1", "true", "yes")

obs = load_observations()
# Build the full spec set (including the sensitivity re-fits and both stages)
# so this job can resolve any id either matrix enumerates — the matrices
# decide which fits run; here we only fit the single requested `id`, so
# listing all is cheap.
specs = build_fit_specs(obs; run_sensitivity = true, cache_dir = CACHE)
isempty(ID) && error(
    "set BVD_FIT_ID to one of: " *
        join((s.id for s in specs), ", ")
)
i = findfirst(s -> s.id == ID, specs)
i === nothing && error(
    "unknown BVD_FIT_ID=$ID; known: " *
        join((s.id for s in specs), ", ")
)

key = fit_key(ID)
@info "fit_one" id = ID key = key cache = CACHE dryrun = DRYRUN refit = REFIT
if DRYRUN
    println(key)
else
    result = fit_or_load(key, specs[i].thunk; cache_dir = CACHE, refit = REFIT)
    @info "cached" id = ID key = key
    ## Keys are rebuilt before the diagnostics read them: a chain that came
    ## back from the cache was serialised by another FlexiChains version.
    repaired = repair_chain_keys(result)
    write_fit_summary(ID, repaired)
    write_fit_extras(ID, key, repaired, CACHE)
end
