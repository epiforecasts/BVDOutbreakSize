# Cold AD-compile cost per component, one fresh process each.
#
# Usage (from the repo root):
#   julia --project=benchmark benchmark/compile.jl [out.json]
#
# Why this is not part of `SUITE`: BenchmarkTools measures a steady-state
# call, and this cost is paid once per process. Mooncake builds the reverse
# rule when the `LogDensityFunction` is constructed, before any gradient is
# taken, so a suite that times `logdensity_and_gradient` in a warm session
# reports microseconds and never sees it. Measured on a 40-day grid, rule
# construction is roughly 87% of a component's cold build, and the full
# `bvd_joint` spends 969 s of its 1095 s cold build there.
#
# Each component runs in its own process because the first build in a
# session pays a fixed floor — about 38 s for a two-parameter model — that
# every later component in the same session avoids. Sharing one process
# would charge the whole floor to whichever component happened to run first.
#
# `BVD_BENCH_JOINT=true` adds the full joint, which costs ~18 min on its own.
using Printf

const OUT_FILE = get(ARGS, 1, "compile.json")
const BENCH_JOINT = lowercase(get(ENV, "BVD_BENCH_JOINT", "false")) == "true"
const REPO = dirname(@__DIR__)

include(joinpath(REPO, "test", "ad_fixtures.jl"))

names = [s.name for s in ADFixtures.scenarios(; joint = BENCH_JOINT)]

results = String[]
for name in names
    tmp = tempname() * ".json"
    cmd = `$(Base.julia_cmd()) --project=$(@__DIR__) $(joinpath(@__DIR__, "compile_one.jl")) $name $tmp`
    println(stderr, "[compile] ", name)
    ## A component that cannot build a rule is reported and skipped rather
    ## than aborting the sweep, matching the smoke test in `src/ad_gradients.jl`.
    if !success(pipeline(cmd; stdout = stderr, stderr = stderr))
        println(stderr, "[compile] skipping $name: build failed")
        continue
    end
    push!(results, read(tmp, String))
    rm(tmp; force = true)
end

open(OUT_FILE, "w") do io
    println(io, "[", join(results, ",\n"), "]")
end
println("Saved compile timings to ", OUT_FILE)
