module BVDOutbreakSizePathfinderExt

import BVDOutbreakSize
using Pathfinder: multipathfinder
using Turing.DynamicPPL: InitFromPrior, InitFromParams
using FlexiChains: FlexiChains
using Random: AbstractRNG

# Multi-path Pathfinder initialisation for NUTS. Fit a quasi-Newton variational
# approximation of `model` from `nruns` independent runs (each seeded from the
# prior), importance-resample `ndraws` constrained draws across the runs, and
# turn `chains` of them into `InitFromParams` strategies, one per chain. Running
# several paths and resampling lets the initial points spread across the joint
# posterior's modes rather than collapse onto one, which matters for the
# small-outbreak ascertainment/background multimodality. Returns a length-
# `chains` vector of init strategies for `nuts_sample` to pass as
# `initial_params`.
function BVDOutbreakSize.pathfinder_init(
        model, chains::Integer;
        nruns::Integer, ndraws::Integer, adtype, rng::AbstractRNG)
    result = multipathfinder(model, ndraws;
        nruns = nruns, init_sampler = InitFromPrior(),
        adtype = adtype, rng = rng)
    ## Constrained per-draw parameters, keyed by the model's variable names.
    draws = result.draws_transformed
    vns = collect(FlexiChains.parameters(draws))
    ## Spread the `chains` initial points evenly across the resampled draw pool
    ## so distinct modes can seed distinct chains.
    idx = ndraws <= chains ? collect(1:chains) :
          round.(Int, range(1, ndraws; length = chains))
    point(i) = InitFromParams(Dict(k => draws[k][i] for k in vns))
    return [point(i) for i in idx]
end

end
