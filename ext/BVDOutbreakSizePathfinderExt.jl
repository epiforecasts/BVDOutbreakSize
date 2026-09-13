module BVDOutbreakSizePathfinderExt

import BVDOutbreakSize
using Pathfinder: multipathfinder, pathfinder
using Turing.DynamicPPL: InitFromPrior, InitFromParams
using FlexiChains: FlexiChains
using Random: AbstractRNG, MersenneTwister, seed!

## One draw of a Pathfinder result, as an `InitFromParams` strategy keyed by
## the model's own variable names.
function _point(draws, i)
    vns = collect(FlexiChains.parameters(draws))
    return InitFromParams(Dict(k => draws[k][i] for k in vns))
end

# Pathfinder initialisation for NUTS: a quasi-Newton variational approximation
# of `model`, turned into one starting point per chain. Returns a length-
# `chains` vector of init strategies for `nuts_sample` to pass as
# `initial_params`.
#
# `multipath = true` (the original behaviour) runs `nruns` paths, importance-
# resamples `ndraws` constrained draws across them, and spreads the chains over
# that shared pool. The intent was for distinct modes to seed distinct chains.
# On the joint model it does the opposite: the importance weights are degenerate
# (Pareto shape well above 1, because the posterior is far from the Gaussian
# approximation), the pool collapses, and every chain receives what is
# effectively the same point. Two chains from one point is a poor basis for
# R-hat, which reads between-chain against within-chain variance.
#
# `multipath = false` instead gives each chain its own single-path run from its
# own seed and its own prior draw, so the starting points stay independent by
# construction and no importance sampling is involved. It costs `chains` paths
# rather than `nruns`.
function BVDOutbreakSize.pathfinder_init(
        model, chains::Integer;
        nruns::Integer, ndraws::Integer, adtype, rng::AbstractRNG,
        multipath::Bool = true)
    if !multipath
        ## A seed per chain, drawn up front, so each path is independent of the
        ## others and the whole set stays reproducible in the caller's `rng`.
        seeds = rand(rng, UInt, chains)
        return map(seeds) do s
            r = MersenneTwister()
            seed!(r, s)
            result = pathfinder(model; ndraws = 1,
                init_sampler = InitFromPrior(), adtype = adtype, rng = r)
            return _point(result.draws_transformed, 1)
        end
    end
    result = multipathfinder(model, ndraws;
        nruns = nruns, init_sampler = InitFromPrior(),
        adtype = adtype, rng = rng)
    ## Spread the `chains` initial points evenly across the resampled draw pool.
    draws = result.draws_transformed
    idx = ndraws <= chains ? collect(1:chains) :
          round.(Int, range(1, ndraws; length = chains))
    return [_point(draws, i) for i in idx]
end

end
