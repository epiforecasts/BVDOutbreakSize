## The two arms of the rules-on against rules-off comparison in
## `test/test_mooncake_rules.jl`. The test runs these in its own process,
## with the rules loaded, and in a child process that loads the package with
## the `mooncake_rules` preference off (`rules_off`).

using BVDOutbreakSize: BVDOutbreakSize, load_observations, production_joint,
    default_breakpoint, default_adtype
using Turing: DynamicPPL
using LogDensityProblems: logdensity_and_gradient
using Mooncake: Mooncake
using Random: Xoshiro

## Whether this process loaded the rules.
rules_loaded() = isdefined(BVDOutbreakSize, :_tangent_array)

## The production joint's log density, as NUTS samples it.
function production_ldf()
    obs = load_observations()
    model = production_joint(obs; breakpoint = default_breakpoint(obs))
    vi = DynamicPPL.link(DynamicPPL.VarInfo(model), model)
    ldf = DynamicPPL.LogDensityFunction(
        model, DynamicPPL.getlogjoint_internal, vi; adtype = default_adtype()
    )
    return model, ldf
end

## Unconstrained points at the `npoints` highest-density of `ndraws` seeded
## prior draws, so no point sits where a likelihood is clamped at its floor.
function prior_points(model; ndraws = 20, npoints = 3)
    draws = map(1:ndraws) do s
        vi = DynamicPPL.VarInfo(Xoshiro(s), model)
        (; lj = DynamicPPL.getlogjoint(vi), vi)
    end
    best = sort(filter(d -> isfinite(d.lj), draws); by = d -> -d.lj)
    return [collect(DynamicPPL.link(d.vi, model)[:]) for d in best[1:npoints]]
end

## Log density and gradient at each point.
function joint_values(ldf, xs)
    return map(xs) do x
        lp, g = logdensity_and_gradient(ldf, x)
        (; lp, g = copy(g))
    end
end

## Fastest time per call over batches of about 20 µs.
function fastest(run; seconds = 0.3)
    run()
    evals = max(1, round(Int, 2.0e-5 / max(@elapsed(run()), 1.0e-9)))
    best = Inf
    stop = time() + seconds
    while time() < stop
        best = min(best, @elapsed(foreach(_ -> run(), 1:evals)) / evals)
    end
    return best
end

## Each case's value and pullback, and the fastest time of one call.
function pullback_times(cases)
    return map(cases) do c
        ȳ = Mooncake.randn_tangent(Xoshiro(1), c.f(c.args...))
        cache = Mooncake.prepare_pullback_cache(c.f, c.args...)
        run() = Mooncake.value_and_pullback!!(cache, ȳ, c.f, c.args...)
        result = deepcopy(run())
        (; c.name, result, time = fastest(run))
    end
end
