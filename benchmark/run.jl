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

# A short per-benchmark budget keeps the run affordable. The minimum-time
# estimator the comparison uses is stable well below the default 5 s, and
# the slowest components here are milliseconds per call.
results = run(SUITE; verbose = true, seconds = 1)
BenchmarkTools.save(out_file, results)
println("Saved benchmark results to ", out_file)
