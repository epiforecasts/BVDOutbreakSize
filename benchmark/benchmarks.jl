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
# Enzyme is loaded so `BVDOutbreakSize.enzyme_adtype()` resolves and
# `ADFixtures.backends()` offers it alongside Mooncake. Which pairs
# actually register is decided by the smoke test in `src/ad_gradients.jl`.

using BenchmarkTools
using Enzyme

include(joinpath(@__DIR__, "..", "test", "ad_fixtures.jl"))

# The full `bvd_joint` is off by default: one gradient is ~14 ms over 76
# parameters behind a cold compile of ~18 min under Mooncake, which no CI
# run can afford. `BVD_BENCH_JOINT=true` adds it for a local investigation.
const BENCH_JOINT = lowercase(get(ENV, "BVD_BENCH_JOINT", "false")) == "true"

const SCENARIOS = ADFixtures.scenarios(; joint = BENCH_JOINT)
const BACKENDS = ADFixtures.backends()

const SUITE = BenchmarkGroup()

include("src/log_density.jl")
include("src/ad_gradients.jl")
