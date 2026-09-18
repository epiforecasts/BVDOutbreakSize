# Print the fit ids as a JSON array, for the `fit` and `fit_dependent`
# matrices in `.github/workflows/docs.yml`. Set `BVD_RUN_SENSITIVITY=true` to
# include the sensitivity re-fits. `BVD_FIT_STAGE` picks the stage: `base`
# (the default) lists the fits with no parent, `dependent` the fits melded
# from a cached parent, and `all` both in registry order.
#
#   julia --project=docs docs/fits/list.jl   # ["joint","exports",...]
#   BVD_FIT_STAGE=dependent julia --project=docs docs/fits/list.jl
using Pkg: Pkg
Pkg.instantiate()

using BVDOutbreakSize
include(joinpath(@__DIR__, "registry.jl"))

ids = fit_ids(; stage = fit_stage_env(:base))
println("[", join(("\"" * id * "\"" for id in ids), ","), "]")
