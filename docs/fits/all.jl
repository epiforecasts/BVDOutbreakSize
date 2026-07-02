# Fit and cache every model from the registry in one parallel pass, for a
# local full build. Mirrors the per-fit CI matrix (`docs/fits/one.jl`) but runs
# the whole registry through `fit_parallel`, writing each chain into the
# content-addressed cache (`BVD_FIT_CACHE`). Run this before
# `docs/make.jl`/`task docs`, which then render from the cache instead of
# refitting. Uses the same `docs/` environment as the render, so the chains are
# serialized and deserialized under identical package versions.
#
#   julia --project=docs docs/fits/all.jl
#
# Set `BVD_FIT_CACHE` to choose the cache directory (default `logs/fit_cache`),
# `BVD_REFIT=all` to ignore existing cache entries, and `BVD_RUN_SENSITIVITY`
# to include the sensitivity re-fits.
using Pkg: Pkg
Pkg.instantiate()

using BVDOutbreakSize
include(joinpath(@__DIR__, "registry.jl"))

const CACHE = get(ENV, "BVD_FIT_CACHE",
    joinpath(pkgdir(BVDOutbreakSize), "logs", "fit_cache"))
const REFIT = lowercase(strip(get(ENV, "BVD_REFIT", ""))) in ("all", "true", "1")

obs = load_observations()
specs = build_fit_specs(obs)
@info "Fitting $(length(specs)) models into the cache" cache=CACHE refit=REFIT
fit_parallel([() -> fit_or_load(fit_key(s.id), s.thunk;
                  cache_dir = CACHE, refit = REFIT) for s in specs])
@info "All fits cached" cache=CACHE
