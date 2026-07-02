# Sampler glue: the package-default AD type and the NUTS driver used
# to fit every Turing model.

"""
Mooncake reverse-mode AD with default `Mooncake.Config()`. Used as
the NUTS `adtype` keyword.
"""
default_adtype() = AutoMooncake(; config = Mooncake.Config())

"""
    enzyme_adtype()

Enzyme reverse-mode AD type, an opt-in alternative to the default
[`default_adtype`](@ref) (Mooncake). Defined by the package's Enzyme
weak-dependency extension (`ext/BVDOutbreakSizeEnzymeExt.jl`); calling it
without `Enzyme` loaded raises a `MethodError`. Load `Enzyme` to activate
the extension. The `SpecialFunctions.gamma` `EnzymeRule` that the Beta and
NegativeBinomial normalising constants reach is supplied by
CensoredDistributions' own Enzyme extension. Enzyme differentiates the
single-stream composers and matches Mooncake; differentiating the full
joint is platform-dependent (it can hit an upstream Enzyme/LLVM compile
failure on some systems, see `test/test_enzyme.jl`), so Mooncake remains
the package default for fitting.
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
Compose several `nuts_sample` step callbacks into one. Each argument is
either a callback with the AbstractMCMC step signature or `nothing`;
`nothing` entries are dropped. The composite invokes the surviving
callbacks in order on every step. Returns the single callback unchanged
when only one survives, and `nothing` when none do (so `nuts_sample` sees
no callback at all rather than a no-op wrapper).

See also: [`fit_callback`](@ref), [`progress_callback`](@ref),
[`tensorboard_callback`](@ref).
"""
function combined_callback(callbacks...)
    cbs = filter(!isnothing, collect(callbacks))
    isempty(cbs) && return nothing
    length(cbs) == 1 && return only(cbs)
    return function (args...; kwargs...)
        for cb in cbs
            cb(args...; kwargs...)
        end
        return nothing
    end
end

"""
Build the logging callback for a named model fit, selected by the
`BVD_FIT_LOG` environment variable (or an explicit `spec`). This is the
wiring the report build uses so every fit streams its progress without
each call site repeating the callback construction.

Recognised `spec` values (case-insensitive), defaulting to `"all"` when
`BVD_FIT_LOG` is unset:

- `"all"` — both the dependency-free [`progress_callback`](@ref) (a
  `<name>.log` file under `logdir`) and the [`tensorboard_callback`](@ref)
  (a `tensorboard/<name>` run directory under `logdir`).
- `"progress"` — the file progress stream only.
- `"tensorboard"` (or `"tb"`) — the TensorBoard stream only.
- `"none"` — no logging; returns `nothing`. CI sets this to keep release
  builds quiet (`BVD_FIT_LOG=none`).

TensorBoard logging needs `TensorBoardLogger` loaded (it activates the
`tensorboard_callback` method through the package extension). When it is
requested but not loaded, the TensorBoard stream is skipped with a warning
rather than erroring, so a build without `using TensorBoardLogger` still
gets the file progress stream.

See also: [`combined_callback`](@ref), [`nuts_sample`](@ref).
"""
function fit_callback(name::AbstractString;
        logdir::AbstractString = "logs",
        spec::AbstractString = get(ENV, "BVD_FIT_LOG", "all"))
    mode = lowercase(strip(spec))
    mode == "none" && return nothing
    want_progress = !(mode in ("tensorboard", "tb"))
    want_tb = !(mode in ("progress", "file"))
    progress = if want_progress
        mkpath(logdir)
        progress_callback(; path = joinpath(logdir, "$(name).log"))
    end
    tb = want_tb ? _tensorboard_if_loaded(joinpath(logdir, "tensorboard",
        name)) : nothing
    return combined_callback(progress, tb)
end

## TensorBoard callback if the extension is active, otherwise `nothing`
## with a warning. Keeps `fit_callback`'s default `"all"` from erroring when
## TensorBoardLogger is not loaded.
function _tensorboard_if_loaded(logdir)
    ext = Base.get_extension(@__MODULE__,
        :BVDOutbreakSizeTensorBoardLoggerExt)
    if isnothing(ext)
        @warn "BVD_FIT_LOG requested TensorBoard logging but " *
              "TensorBoardLogger is not loaded; logging file progress only. " *
              "Add `using TensorBoardLogger` to enable it."
        return nothing
    end
    return tensorboard_callback(logdir)
end

"""
NUTS on `model`, parallel chains via `MCMCThreads`. Chains
initialise from the prior (`InitFromPrior()`) to keep the sampler
in regions with reasonable physical interpretation. Pass `init =
Turing.DynamicPPL.InitFromUniform()` to fall back to unconstrained
uniform initialisation.

`target_accept` defaults to 0.85. The earlier integral model needed 0.95
to keep the multimodal small-outbreak geometry from diverging, but the
renewal joint conditions the confirmed counts on the observed analysed
denominator (removing the multiplicative ascertainment ridge) and samples
the random-walk and ascertainment blocks in non-centred form, so the
posterior geometry is benign (the sanity fit converges with ≈1 divergence).
A lower target acceptance shortens the average NUTS trajectory, cutting
leapfrog steps (and so gradient evaluations) per iteration, so 0.85 trims
the per-iteration gradient cost while staying above the conventional 0.8
floor; raise it back toward 0.9–0.99 if a model variant reintroduces
divergences. The default
is two longer chains (1000 post-warmup draws each) rather than four shorter
ones, mirroring the integral model (#211), which roughly halves the docs
build at a similar total draw count.

`check_model = false` disables Turing's pre-sampling model check, which
rejects any model with a sampled discrete variable even when its value
feeds nothing downstream. The per-vintage DRC streams are now scored as
observed `~` data, so a composer conditioning on them passes the check
with the default `check_model = true`. The escape is needed only by
[`exports_deaths_only_model`](@ref), which runs the exports submodel in
predictive mode (`exported_cases ~ Poisson` with a `missing` count) purely
for the export onsets, leaving a sampled discrete `Poisson` draw. The
continuous parameters are unaffected.

Pass `callback` to stream live fit progress (iteration, log-density,
divergences) instead of waiting for the whole fit. Use
[`progress_callback`](@ref) for a dependency-free file/stdout stream,
or [`tensorboard_callback`](@ref) for a TensorBoard backend (requires
`using TensorBoardLogger`). The callback is forwarded to `sample` only
when non-`nothing`. Any additional `kwargs` are passed through to
`sample`.

`n_adapts` sets the NUTS warmup length (step-size and mass-matrix
adaptation), run in addition to `samples` and discarded by default. It
defaults to `min(250, samples ÷ 2)`, trimming the per-fit warmup from
Turing's default of `min(1000, samples ÷ 2)` (500 at the standard 1000
draws) to speed the report build.

A callback fires only on the samples that are kept, and NUTS discards
its adaptation phase by default, so warmup is silent. Set
`warmup = true` to keep the adaptation steps (`discard_adapt = false`),
which streams them to the callback so step-size adaptation and early
divergences are visible live. Those warmup draws are then also retained
in the returned chain, so the first `n_adapts` draws are adaptation
steps rather than posterior samples; raise `samples` accordingly or drop
them before summarising.
"""
function nuts_sample(model;
        samples::Integer = 1_000,
        chains::Integer = 2,
        target_accept::Real = 0.85,
        n_adapts::Integer = min(250, samples ÷ 2),
        seed::Integer = 20260518,
        progress::Bool = false,
        adtype = default_adtype(),
        init = InitFromPrior(),
        check_model::Bool = true,
        callback = nothing,
        warmup::Bool = false,
        kwargs...)
    rng = MersenneTwister(seed)
    cb_kwargs = callback === nothing ? (;) : (; callback = callback)
    warmup_kwargs = warmup ? (; discard_adapt = false) : (;)
    return sample(
        rng,
        model,
        NUTS(n_adapts, target_accept; adtype),
        MCMCThreads(),
        samples, chains;
        initial_params = fill(init, chains),
        progress = progress,
        check_model = check_model,
        cb_kwargs...,
        warmup_kwargs...,
        kwargs...
    )
end

"""
    fit_parallel(thunks; chains = 2)

Run independent model fits — each a zero-argument `thunk` returning a chain —
with model-level parallelism bounded by the available threads. At most
`Threads.nthreads() ÷ chains` fits run at once (so each fit keeps `chains`
threads for its own chains), clamped to the number of fits. This is
self-limiting and CI-safe: with two threads (e.g. CI's
`JULIA_NUM_THREADS=2`, the default `chains`) it runs the fits SEQUENTIALLY,
identical to a plain loop and with the same peak memory; on a many-core
machine with more threads it fans the fits out (eight threads → four fits at
once). Each fit seeds its own RNG, so the results do not depend on the
schedule. Returns the chains in input order.
"""
function fit_parallel(thunks::AbstractVector; chains::Integer = 2)
    n = length(thunks)
    maxconc = clamp(Threads.nthreads() ÷ max(chains, 1), 1, n)
    results = Vector{Any}(undef, n)
    if maxconc == 1
        for i in 1:n
            results[i] = thunks[i]()
        end
        return results
    end
    next = Threads.Atomic{Int}(1)
    @sync for _ in 1:maxconc
        Threads.@spawn begin
            while true
                i = Threads.atomic_add!(next, 1)
                i > n && break
                results[i] = thunks[i]()
            end
        end
    end
    return results
end
