# Benchmark suite for BVDOutbreakSize. Defines a BenchmarkTools
# `BenchmarkGroup` named `SUITE`, which `run.jl` executes and `compare.jl`
# turns into a PR comment.
#
# Groups:
#   "Log density"  — one unconstrained log-density evaluation per component,
#       the denominator the AD numbers are read against.
#   "AD gradients" — the gradient of that same log-density, per component
#       per backend, keyed `["AD gradients"][group][component][backend]`.
#
# The components come from `test/ad_fixtures.jl`, which
# `test/test_ad_gradients.jl` includes too, so the surface the tests assert
# is differentiable is the surface timed here.
#
# Enzyme is the second backend and is off by default. Loading it is what
# makes `BVDOutbreakSize.enzyme_adtype()` resolve, so leaving it out is
# what leaves `ADFixtures.backends()` reporting Mooncake alone.
# `BVD_BENCH_ENZYME=true` turns it on. Both backends together do not fit a
# CI run: the sweep was cancelled at the 90 minute cap, where Mooncake
# alone finished in 39 minutes on the same cold runner. Which pairs then
# register is decided by the smoke test in `src/ad_gradients.jl`.

using BenchmarkTools

include(joinpath(@__DIR__, "..", "test", "ad_fixtures.jl"))

const BENCH_ENZYME = lowercase(get(ENV, "BVD_BENCH_ENZYME", "false")) == "true"
if BENCH_ENZYME
    @eval using Enzyme
end

# The full `bvd_joint`: one gradient is ~14 ms over 76 parameters behind a
# cold compile of ~18 min under Mooncake. Off by default so a local run of
# the component suite stays quick; the benchmark workflow sets
# `BVD_BENCH_JOINT=true`, and that comparison is the only place the joint's
# gradient is timed. The test suite asserts the components differentiate and
# leaves the joint to this and to the docs build's fits.
const BENCH_JOINT = lowercase(get(ENV, "BVD_BENCH_JOINT", "false")) == "true"

const SCENARIOS = ADFixtures.scenarios(; joint = BENCH_JOINT)
const BACKENDS = ADFixtures.backends()

const SUITE = BenchmarkGroup()

include("src/log_density.jl")
include("src/ad_gradients.jl")
