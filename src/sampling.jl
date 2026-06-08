# Sampler glue: the package-default AD type and the NUTS driver used
# to fit every Turing model.

"""
Mooncake reverse-mode AD with default `Mooncake.Config()`. Used as
the NUTS `adtype` keyword.
"""
default_adtype() = AutoMooncake(; config = Mooncake.Config())

"""
Enzyme reverse-mode AD with runtime activity and `Duplicated` function
annotation, the configuration the joint model's analytical gamma-CDF
rule needs (see `ext/BVDOutbreakSizeEnzymeExt.jl`). Returns an
`ADTypes.AutoEnzyme`; pass to `nuts_sample(model; adtype = ...)`.

`enzyme_adtype` is a stub; loading Enzyme (`using Enzyme`) activates
the method via `BVDOutbreakSizeEnzymeExt`. Calling `enzyme_adtype`
without Enzyme loaded raises a `MethodError`.
"""
function enzyme_adtype end

"""
TensorBoard streaming callback for `nuts_sample`, mirroring the
optional Enzyme backend. `tensorboard_callback(logdir; every = 20,
histograms = true)` opens a `TensorBoardLogger.TBLogger(logdir)` and, on
every kept post-warmup draw, logs through the sampler-agnostic
`AbstractMCMC.ParamsWithStats` interface to two grouped tag prefixes so
the dashboard stays navigable:

  * `params/<name>` — every sampled parameter
  * `diagnostics/<name>` — log-density (`logjoint`), divergence flag
    (`numerical_error`), step size, tree depth, acceptance rate, ...

Each scalar streams every step as a `.../value` time series. With
`histograms = true` (the default) a running histogram of the draws so
far is also logged every `every` steps as `.../distribution`,
populating the TensorBoard HISTOGRAMS and DISTRIBUTIONS dashboards. Set
`histograms = false` for scalar traces only, or widen `every` to log
histograms less often.

`tensorboard_callback` is a stub; loading TensorBoardLogger
(`using TensorBoardLogger`) activates the method via
`BVDOutbreakSizeTensorBoardLoggerExt`. Calling it without
TensorBoardLogger loaded raises an informative `ErrorException`.

Pass the result to `nuts_sample`:

```julia
using TensorBoardLogger
nuts_sample(model; callback = tensorboard_callback("logs/run"))
```

then view the run with `tensorboard --logdir logs/run`. Use `chains = 1`
for clean live traces; parallel chains share one logger and interleave.

See also: [`progress_callback`](@ref), [`nuts_sample`](@ref).
"""
function tensorboard_callback(logdir; kwargs...)
    throw(ErrorException(
        "tensorboard_callback requires TensorBoardLogger. Load it with " *
        "`using TensorBoardLogger` (it is an optional dependency) and " *
        "retry. logdir = $(repr(logdir))"))
end

@doc raw"""
A lightweight, dependency-free streaming progress callback for
`nuts_sample`. Returns a closure matching the AbstractMCMC callback
signature

```julia
callback(rng, model, sampler, transition, state, iteration; kwargs...)
```

which, every `every` iterations on each chain, appends one line to
the file at `path` recording the iteration number, the log joint
density, and a running count of divergent transitions.

The same closure is invoked from every thread `MCMCThreads()` spawns,
so the divergence tally and file writes are shared across chains and
guarded by a `ReentrantLock`. Tail the file live during a fit with
`tail -f <path>`.

Step-level statistics are read through the sampler-agnostic
`AbstractMCMC.ParamsWithStats(model, sampler, transition, state;
stats = true)` interface rather than by reaching into transition
fields, so the callback tracks Turing's transition format instead of a
fixed field layout. The log density is taken from the `logjoint`
statistic and the divergence flag from `numerical_error`. The whole
body is wrapped in `try`/`catch`, so a transition that does not expose
these statistics yields `missing` log-density / no divergence increment
rather than crashing the fit. The running divergence count is a single
total over all chains, reset implicitly each time a fresh callback is
constructed.

See also: [`tensorboard_callback`](@ref), [`nuts_sample`](@ref).
"""
function progress_callback(; path::AbstractString, every::Integer = 10)
    lock = ReentrantLock()
    ndivergent = Ref(0)
    return function (rng, model, sampler, transition, state, iteration;
            kwargs...)
        try
            stats = try
                AbstractMCMC.ParamsWithStats(
                    model, sampler, transition, state; stats = true).stats
            catch
                (;)
            end
            if get(stats, :numerical_error, false) === true
                Base.@lock lock (ndivergent[] += 1)
            end
            if iteration % every == 0
                lp = get(stats, :logjoint, missing)
                Base.@lock lock begin
                    open(path, "a") do io
                        println(io,
                            "iteration=$(iteration) lp=$(lp) " *
                            "divergences=$(ndivergent[])")
                    end
                end
            end
        catch
            # A streaming progress callback must never abort a fit.
        end
        return nothing
    end
end

"""
NUTS on `model`, parallel chains via `MCMCThreads`. Chains
initialise from the prior (`InitFromPrior()`) to keep the sampler
in regions with reasonable physical interpretation. Pass `init =
Turing.DynamicPPL.InitFromUniform()` to fall back to unconstrained
uniform initialisation.

Pass `callback` to stream live fit progress (iteration, log-density,
divergences) instead of waiting for the whole fit. Use
[`progress_callback`](@ref) for a dependency-free file/stdout stream,
or [`tensorboard_callback`](@ref) for a TensorBoard backend (requires
`using TensorBoardLogger`). The callback is forwarded to `sample` only
when non-`nothing`. Any additional `kwargs` are passed through to
`sample`.

A callback fires only on the samples that are kept, and NUTS discards
its adaptation phase by default, so warmup is silent. Set
`warmup = true` to keep the adaptation steps (`discard_adapt = false`),
which streams them to the callback so step-size adaptation and early
divergences are visible live. Those warmup draws are then also retained
in the returned chain, so the first `min(1000, samples ÷ 2)` draws are
adaptation steps rather than posterior samples; raise `samples`
accordingly or drop them before summarising.
"""
function nuts_sample(model;
        samples::Integer = 1_000,
        chains::Integer = 2,
        target_accept::Real = 0.9,
        seed::Integer = 20260518,
        progress::Bool = false,
        adtype = default_adtype(),
        init = InitFromPrior(),
        callback = nothing,
        warmup::Bool = false,
        kwargs...)
    rng = MersenneTwister(seed)
    cb_kwargs = callback === nothing ? (;) : (; callback = callback)
    warmup_kwargs = warmup ? (; discard_adapt = false) : (;)
    return sample(
        rng,
        model,
        NUTS(target_accept; adtype),
        MCMCThreads(),
        samples, chains;
        initial_params = fill(init, chains),
        progress = progress,
        cb_kwargs...,
        warmup_kwargs...,
        kwargs...
    )
end
