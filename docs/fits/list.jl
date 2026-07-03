# Print the fit ids as a JSON array, for the `fit` matrix in
# `.github/workflows/fit-matrix.yml`. Set `BVD_RUN_SENSITIVITY=true` to include
# the sensitivity re-fits.
#
#   julia --project=docs docs/fits/list.jl   # ["joint","exports",...]
using Pkg: Pkg
Pkg.instantiate()

using BVDOutbreakSize
include(joinpath(@__DIR__, "registry.jl"))

ids = fit_ids()
println("[", join(("\"" * id * "\"" for id in ids), ","), "]")
