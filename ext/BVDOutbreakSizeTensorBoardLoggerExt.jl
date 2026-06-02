module BVDOutbreakSizeTensorBoardLoggerExt

import BVDOutbreakSize
using TensorBoardLogger: TBLogger, log_value

# Real method for the `tensorboard_callback` stub in `src/sampling.jl`.
# Loading TensorBoardLogger (`using TensorBoardLogger`) activates this
# method, so `nuts_sample(model; callback = tensorboard_callback(dir))`
# streams per-step NUTS diagnostics to the TensorBoard event files under
# `logdir`. View the run with `tensorboard --logdir <logdir>`.
#
# The returned closure matches the AbstractMCMC callback signature
# `callback(rng, model, sampler, transition, state, iteration; kwargs...)`
# and, every `every` steps on each chain, logs scalar summaries: the
# transition log-density (`transition.lp`), a divergence flag and running
# divergence count (`transition.stat.numerical_error`), and step size /
# tree depth when the sampler exposes them (`transition.stat.step_size`,
# `transition.stat.tree_depth`).
#
# A single `TBLogger` is shared across the threads spawned by
# `MCMCThreads()`; writes are guarded by a `ReentrantLock`. Field access
# is wrapped in `try`/`catch` so a transition missing any field logs what
# it can rather than aborting the fit. `kwargs` are forwarded to
# `TBLogger`.
function BVDOutbreakSize.tensorboard_callback(
        logdir::AbstractString; every::Integer = 1, kwargs...)
    logger = TBLogger(logdir; kwargs...)
    lock = ReentrantLock()
    ndivergent = Ref(0)
    return function (rng, model, sampler, transition, state, iteration;
            cbkwargs...)
        try
            divergent = false
            try
                divergent = transition.stat.numerical_error === true
            catch
                divergent = false
            end
            if divergent
                Base.@lock lock (ndivergent[] += 1)
            end
            if iteration % every == 0
                Base.@lock lock begin
                    _log_scalar(logger, "lp", _field(transition, :lp),
                        iteration)
                    log_value(logger, "divergent", divergent;
                        step = iteration)
                    log_value(logger, "divergences", ndivergent[];
                        step = iteration)
                    stat = _field(transition, :stat)
                    if stat !== missing
                        _log_scalar(logger, "step_size",
                            _field(stat, :step_size), iteration)
                        _log_scalar(logger, "tree_depth",
                            _field(stat, :tree_depth), iteration)
                    end
                end
            end
        catch
            # A streaming callback must never abort a fit.
        end
        return nothing
    end
end

# Read field `name` from `x`, returning `missing` if it is absent.
function _field(x, name::Symbol)
    return try
        getproperty(x, name)
    catch
        missing
    end
end

# Log `value` as a scalar only when it is a present, finite real.
function _log_scalar(logger, name, value, step)
    if value !== missing && value isa Real && isfinite(value)
        log_value(logger, name, float(value); step = step)
    end
    return nothing
end

end
