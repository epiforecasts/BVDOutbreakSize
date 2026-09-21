# Run this checkout's benchmark suite and save the results to JSON.
#
# Usage (from the repo root):
#   julia --project=benchmark benchmark/run.jl [out.json]
#
# `BVD_BENCH_JOINT=true` adds the full `bvd_joint` component (see
# `benchmarks.jl`); it is off by default.
using BenchmarkTools

out_file = get(ARGS, 1, "results.json")

include(joinpath(@__DIR__, "benchmarks.jl"))  # defines `SUITE`

# The budget comes from `DEFAULT_PARAMETERS` in `benchmarks.jl`, so a local
# run and a CI arm sample the same way.
results = run(SUITE; verbose = true)
BenchmarkTools.save(out_file, results)
println("Saved benchmark results to ", out_file)
