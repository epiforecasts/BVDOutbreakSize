# Cold AD-compile cost of one component, in a fresh process.
#
# Usage (from the repo root):
#   julia --project=benchmark benchmark/compile_one.jl <component> [out.json]
#
# `benchmark/compile.jl` drives this once per component. It is a separate
# process per component because the cost being measured is paid once per
# process: Mooncake builds the reverse rule when the `LogDensityFunction` is
# constructed, and every later component in the same session reuses whatever
# that build already inferred.
using Printf
using BVDOutbreakSize: BVDOutbreakSize, default_adtype
using LogDensityProblems: LogDensityProblems, logdensity_and_gradient
using Turing: DynamicPPL

include(joinpath(@__DIR__, "..", "test", "ad_fixtures.jl"))

const COMPONENT = get(ARGS, 1, "")
const OUT_FILE = get(ARGS, 2, "")

scens = ADFixtures.scenarios(; joint = (COMPONENT == "bvd_joint"))
idx = findfirst(s -> s.name == COMPONENT, scens)
idx === nothing && error(
    "no component named $(repr(COMPONENT)). Have: " *
        join((s.name for s in scens), ", ")
)
scen = scens[idx]

vi, x = ADFixtures.linked_point(scen.model; seed = scen.seed)

## Forward-only first: the primal log density compiles the model body with
## no AD rule attached, so the difference is the rule rather than the model.
t0 = time_ns()
ldf_primal = DynamicPPL.LogDensityFunction(
    scen.model, DynamicPPL.getlogjoint, vi
)
LogDensityProblems.logdensity(ldf_primal, x)
primal_s = (time_ns() - t0) / 1.0e9

t0 = time_ns()
ldf = DynamicPPL.LogDensityFunction(
    scen.model, DynamicPPL.getlogjoint, vi; adtype = default_adtype()
)
rule_s = (time_ns() - t0) / 1.0e9

logp, grad = logdensity_and_gradient(ldf, x)

@printf("%-40s primal %7.2f s  rule %8.2f s\n", scen.name, primal_s, rule_s)

if !isempty(OUT_FILE)
    open(OUT_FILE, "w") do io
        println(io, "{")
        @printf(io, "  \"component\": %s,\n", repr(scen.name))
        @printf(io, "  \"group\": %s,\n", repr(scen.group))
        @printf(io, "  \"nparams\": %d,\n", length(x))
        @printf(io, "  \"primal_s\": %.4f,\n", primal_s)
        @printf(io, "  \"rule_s\": %.4f,\n", rule_s)
        @printf(io, "  \"logp\": %.17g,\n", logp)
        println(io, "  \"grad\": [", join((repr(g) for g in grad), ","), "]")
        println(io, "}")
    end
end
