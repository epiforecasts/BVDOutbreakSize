module BVDOutbreakSizeTensorBoardLoggerExt

import BVDOutbreakSize
import AbstractMCMC
using TensorBoardLogger: TBLogger, log_value, log_histogram

# Real method for the `tensorboard_callback` stub in `src/sampling.jl`.
# Loading TensorBoardLogger (`using TensorBoardLogger`) activates it.
#
# Parameters and step statistics are pulled through the sampler-agnostic
# `AbstractMCMC.ParamsWithStats(model, sampler, transition, state; params,
# stats)` interface (so the callback tracks Turing's transition format rather
# than a fixed field layout), then logged to the `TBLogger` at `logdir` under
# two grouped tag prefixes so the TensorBoard dashboard stays navigable:
#
#   * `params/<name>`      — every sampled parameter
#   * `diagnostics/<name>` — log-density (`logjoint`), divergence flag
#                            (`numerical_error`), step size, tree depth, ...
#
# Each scalar is logged every step as a `.../value` time series via
# `log_value` (with an explicit per-chain `step`, so parallel chains do not
# fight over a shared step counter). When `histograms = true` (the default) a
# running histogram of the draws so far is also logged every `every` steps as
# `.../distribution` via `log_histogram`, populating the TensorBoard HISTOGRAMS
# / DISTRIBUTIONS dashboards. The returned closure matches the AbstractMCMC
# callback signature and is passed straight to
# `nuts_sample(model; callback = ...)`. A single `TBLogger` is shared across
# the threads `MCMCThreads()` spawns; writes are guarded by a `ReentrantLock`
# and the whole body is wrapped in `try`/`catch` so a fit is never aborted by
# the callback. View the run with `tensorboard --logdir <logdir>`.
function BVDOutbreakSize.tensorboard_callback(
        logdir::AbstractString;
        every::Integer = 20,
        histograms::Bool = true,
        params_prefix::AbstractString = "params/",
        diagnostics_prefix::AbstractString = "diagnostics/")
    logger = TBLogger(logdir)
    lk = ReentrantLock()
    history = Dict{String, Vector{Float64}}()
    return function (rng, model, sampler, transition, state, iteration;
            kwargs...)
        try
            pws = AbstractMCMC.ParamsWithStats(
                model, sampler, transition, state; params = true, stats = true)
            Base.@lock lk begin
                _emit!(logger, history, params_prefix, pairs(pws.params),
                    iteration, every, histograms)
                _emit!(logger, history, diagnostics_prefix, pairs(pws.stats),
                    iteration, every, histograms)
            end
        catch
            # A streaming callback must never abort a fit.
        end
        return nothing
    end
end

# Log each real-valued entry of `nt_pairs` as a `<prefix><name>/value` scalar
# every step, and, when `histograms`, a `<prefix><name>/distribution` histogram
# of the accumulated draws every `every` steps.
function _emit!(logger, history, prefix, nt_pairs, iteration, every, histograms)
    for (name, value) in nt_pairs
        (value isa Real && isfinite(value)) || continue
        v = Float64(value)
        tag = string(prefix, name)
        log_value(logger, string(tag, "/value"), v; step = iteration)
        if histograms
            draws = get!(() -> Float64[], history, tag)
            push!(draws, v)
            if iteration % every == 0 && length(draws) > 1
                log_histogram(logger, string(tag, "/distribution"), draws;
                    step = iteration)
            end
        end
    end
    return nothing
end

end
