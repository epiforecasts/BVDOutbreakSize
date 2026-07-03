# Run and cache a SINGLE fit from the registry, for the per-fit CI matrix (or
# an HPC task). Select the fit with the `BVD_FIT_ID` environment variable; the
# result is serialised into `BVD_FIT_CACHE` under its content-addressed key, so
# the docs build can load it instead of refitting.
#
#   BVD_FIT_ID=confirmed julia --project=docs docs/fits/one.jl
#
# Set `BVD_FIT_CACHE` to choose the cache directory (default `logs/fit_cache`),
# `BVD_REFIT=all` to ignore an existing cache entry, and `BVD_FIT_DRYRUN=1` to
# print the id and key without fitting (a cheap check that the id resolves).
using Pkg: Pkg
Pkg.instantiate()

using BVDOutbreakSize
include(joinpath(@__DIR__, "registry.jl"))

const ID = strip(get(ENV, "BVD_FIT_ID", ""))
const CACHE = get(ENV, "BVD_FIT_CACHE",
    joinpath(pkgdir(BVDOutbreakSize), "logs", "fit_cache"))
const REFIT = lowercase(strip(get(ENV, "BVD_REFIT", ""))) in ("all", "true", "1")
const DRYRUN = lowercase(strip(get(ENV, "BVD_FIT_DRYRUN", ""))) in ("1", "true", "yes")

obs = load_observations()
specs = build_fit_specs(obs)
isempty(ID) && error("set BVD_FIT_ID to one of: " *
                     join((s.id for s in specs), ", "))
i = findfirst(s -> s.id == ID, specs)
i === nothing && error("unknown BVD_FIT_ID=$ID; known: " *
      join((s.id for s in specs), ", "))

key = fit_key(ID)
@info "fit_one" id=ID key=key cache=CACHE dryrun=DRYRUN refit=REFIT
if DRYRUN
    println(key)
else
    fit_or_load(key, specs[i].thunk; cache_dir = CACHE, refit = REFIT)
    @info "cached" id=ID key=key
end
