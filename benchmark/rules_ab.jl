# Times the same components with `src/mooncake_rules.jl` loaded and not loaded,
# so the hand-written rules can be shown to still earn their place. The
# rules are chosen by a package preference read at load time, so each arm
# needs its own process; this script is one arm, and `task benchmark-rules`
# runs both.
#
# Run it when the AD backend is upgraded. If an arm's advantage has gone,
# the rules are dead weight and should be deleted rather than maintained.

using BVDOutbreakSize
using BVDOutbreakSize: default_adtype
using LogDensityProblems: logdensity_and_gradient
using Preferences: load_preference
using Turing: DynamicPPL
using Random: seed!

include(joinpath(@__DIR__, "..", "test", "ad_fixtures.jl"))

const RULES_ON = load_preference(BVDOutbreakSize, "mooncake_rules", true)

function timed(scen)
    seed!(scen.seed)
    vi = DynamicPPL.link(DynamicPPL.VarInfo(scen.model), scen.model)
    x = collect(vi[:])
    ldf = DynamicPPL.LogDensityFunction(
        scen.model, DynamicPPL.getlogjoint, vi; adtype = default_adtype()
    )
    logdensity_and_gradient(ldf, x)
    best = Inf
    for _ in 1:5
        t0 = time()
        for _ in 1:50
            logdensity_and_gradient(ldf, x)
        end
        best = min(best, (time() - t0) / 50)
    end
    return best
end

println("mooncake_rules = ", RULES_ON)
for scen in ADFixtures.scenarios()
    scen.name in ADFixtures.enzyme_skip_scenarios() && continue
    println(rpad(scen.name, 36), round(timed(scen) * 1.0e3; sigdigits = 4))
    flush(stdout)
end
