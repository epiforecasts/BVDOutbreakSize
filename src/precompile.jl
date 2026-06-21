# Precompile the expensive first-call work so a fresh process does not pay it
# on its first fit. The dominant cost is Mooncake building the reverse rule for
# `bvd_joint` the first time it differentiates the model (minutes). Compiling a
# single log-density gradient of a tiny synthetic joint here bakes those rules
# into the package precompile cache; CI persists that cache (julia-actions/
# cache), so both local and CI builds skip the compile. The Mooncake rules are
# keyed by method signature, not data size, so a 40-day synthetic fit compiles
# what the full-size fits reuse. The gradient path (not a NUTS fit) is the exact
# expensive operation and avoids sampler-adaptation fragility at tiny sizes.

using PrecompileTools: @setup_workload, @compile_workload
using LogDensityProblems: logdensity_and_gradient
import Turing: DynamicPPL

@setup_workload begin
    dh = (; days = [13, 18, 40], counts = [10, 14, 18])
    rh = (; days = [13, 18, 40], counts = [340, 516, 905])
    ch = (; days = [13, 18, 40], counts = [9, 17, 27])
    lh = (; days = [18, 40], counts = [30, 50])
    @compile_workload begin
        ## The full joint (the hot path) and a single-stream composer, each
        ## differentiated once under the default (Mooncake) backend.
        models = (
            bvd_joint(40, 2, 18, 905, 0, 27, 50;
                confirmed_deaths = 5, deaths_history = dh,
                reported_history = rh, confirmed_history = ch,
                lab_history = lh, breakpoint = 30),
            exports_only_model(40, 2))
        for m in models
            vi = DynamicPPL.link(DynamicPPL.VarInfo(m), m)
            ldf = DynamicPPL.LogDensityFunction(
                m, DynamicPPL.getlogjoint, vi; adtype = default_adtype())
            logdensity_and_gradient(ldf, collect(vi[:]))
        end
    end
end
