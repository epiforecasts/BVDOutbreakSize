# Warm-start adaptation for the registry fits. A fit that finishes writes a
# warm-start state: its draws in the linked parameter space, its step size,
# its convergence verdict and the model-source hash it ran on. The next fit
# of the same id reads that state, builds a dense metric from the draws,
# holds it fixed and adapts only the step size over a short warm-up
# (`FixedMetricNUTS`). With no usable state it runs the full `NUTS` warm-up.
#
# The state lives outside the fit cache, in `BVD_WARMSTART_OUT`, and the
# parent is read from `BVD_WARMSTART_PARENT`. In CI each is its own small
# Actions cache entry per fit id, so a parent survives the eviction of the
# large fit caches.
#
# The parent is not part of the fit-cache key. The list job computes one key
# before any parent is known, and a key that named the parent would change
# on every run, so the cache would never hit. The parent changes how fast a
# fit converges, not the distribution it samples, so a cached fit stays valid
# for its key. The key does carry the warm-start settings, and the state each
# fit writes records the parent it started from.

using Base.ScopedValues: ScopedValue, with
using Serialization: serialize, deserialize
using LinearAlgebra: Symmetric, cholesky, issuccess, diag, I
using Statistics: mean, median
using Dates: now

## Bump when the state layout changes, so older states are not read.
const WARM_START_SCHEMA = "v1"

warm_start_enabled() = lowercase(strip(get(ENV, "BVD_WARMSTART", "true"))) in
    ("1", "true", "yes")

## Warm-up iterations of a seeded fit. The full warm-up is the joint's.
warm_start_adapts() = parse(Int, get(ENV, "BVD_WARMSTART_ADAPTS", "300"))

## Weight on the identity when the parent's correlation matrix is shrunk.
warm_start_shrink() = parse(Float64, get(ENV, "BVD_WARMSTART_SHRINK", "0.75"))

## The least share of the metric that must map from the parent to seed.
## A model change since the parent does not stop a seed on its own: the
## metric affects only how fast the sampler moves, and the step size still
## adapts.
const WARM_START_MIN_COVERAGE = 0.9

warm_start_parent_dir() = _repo_dir(
    get(ENV, "BVD_WARMSTART_PARENT", joinpath("logs", "warmstart_parent"))
)
warm_start_out_dir() = _repo_dir(
    get(ENV, "BVD_WARMSTART_OUT", joinpath("logs", "warmstart"))
)
_repo_dir(d) = isabspath(d) ? d : joinpath(_PKG, d)

## The warm-start settings, for the fit-cache key.
function warm_start_key()
    warm_start_enabled() || return "warm:off"
    return string(
        "warm:", WARM_START_SCHEMA, ":", warm_start_adapts(), ":",
        warm_start_shrink(), ":", WARM_START_MIN_COVERAGE
    )
end

## Digest of the model sources, logged when it differs from the parent's.
warm_start_models_hash() = tree_sha256(joinpath(_PKG, "src", "models"))

warm_start_file(dir, key) = joinpath(dir, key * ".warmstart.jls")

"""
The warm start of one fit: its id and key, the parent state it may seed
from (`nothing` for none), and why. `path` records which warm-up ran.
"""
mutable struct FitWarmStart
    id::String
    key::String
    parent::Union{Nothing, NamedTuple}
    reason::String
    path::String
    coverage::Float64
end

## The warm start the current fit's sampler reads through `fit_sampler_args`.
const CURRENT_WARM_START = ScopedValue{Any}(nothing)

"""
The registry's sampler settings, with the current fit's warm start. Use in
place of `joint_sampler_args()` where a fit calls `nuts_sample`.
"""
fit_sampler_args() =
    (; joint_sampler_args()..., warm_start = CURRENT_WARM_START[])

"""
    fit_warm_start(id, key; parent_dir = warm_start_parent_dir())

The warm start for fit `id` at `key`: the newest state for `id` in
`parent_dir` that passes the checks, or none with the reason.
"""
function fit_warm_start(id, key; parent_dir = warm_start_parent_dir())
    none(reason) = FitWarmStart(id, key, nothing, reason, "full", 0.0)
    warm_start_enabled() || return none("warm start is off")
    pattern = Regex("^" * id * "__[0-9a-f]+\\.warmstart\\.jls\$")
    files = isdir(parent_dir) ?
        filter(f -> occursin(pattern, f), readdir(parent_dir)) : String[]
    isempty(files) && return none("no parent state for $id")
    newest = last(sort(files; by = f -> mtime(joinpath(parent_dir, f))))
    path = joinpath(parent_dir, newest)
    state = try
        deserialize(path)
    catch e
        return none("parent state unreadable: $(sprint(showerror, e))")
    end
    get(state, :schema, "") == WARM_START_SCHEMA ||
        return none("parent state schema $(get(state, :schema, "?"))")
    state.key == key && return none("parent is this fit")
    state.gate === :fail &&
        return none("parent $(state.key) failed the convergence gate")
    reason = "parent $(state.key)"
    if state.models_hash != warm_start_models_hash()
        reason *= ", with src/models changed since"
    end
    return FitWarmStart(id, key, state, reason, "full", 0.0)
end

"""
    seeded_inverse_metric(state, layout; shrink = warm_start_shrink())

The inverse metric for a model with linked `layout`, from the parent
`state`'s draws: their within-chain covariance with the correlations
shrunk towards the identity by `shrink`. A variable in both maps on the
leading entries the two share. A variable or entry new since the parent
takes the median parent variance of its variable, or of all parameters
when its variable is new, with no correlation. Returns
`(; inverse_metric, coverage)`, `coverage` the share of entries mapped.
"""
function seeded_inverse_metric(state, layout; shrink = warm_start_shrink())
    X = copy(state.draws)
    for c in unique(state.chain_index)
        rows = state.chain_index .== c
        X[rows, :] .-= mean(X[rows, :]; dims = 1)
    end
    S = (X' * X) ./ (size(X, 1) - 1)
    parent_ranges = Dict(state.layout)
    dim = isempty(layout) ? 0 : maximum(last(r) for (_, r) in layout)
    child = Int[]
    source = Int[]
    fill_var = zeros(dim)
    overall = median(diag(S))
    for (name, r) in layout
        pr = get(parent_ranges, name, nothing)
        n = pr === nothing ? 0 : min(length(r), length(pr))
        append!(child, r[1:n])
        append!(source, pr === nothing ? Int[] : pr[1:n])
        fill_var[r] .= pr === nothing ? overall : median(diag(S)[pr])
    end
    Σ = zeros(dim, dim)
    for j in 1:dim
        Σ[j, j] = fill_var[j]
    end
    Σ[child, child] .= S[source, source]
    sd = sqrt.(max.(diag(Σ), 1.0e-12))
    R = Σ ./ (sd * sd')
    Rs = (1 - shrink) .* R .+ shrink .* Matrix{Float64}(I, dim, dim)
    return (; inverse_metric = Rs .* (sd * sd'), coverage = length(child) / dim)
end

function BVDOutbreakSize.warm_start_settings(ws::FitWarmStart, model)
    if ws.parent === nothing
        @info "warm start: full warm-up" id = ws.id reason = ws.reason
        return nothing
    end
    layout = linked_layout(model)
    seeded = seeded_inverse_metric(ws.parent, layout)
    ws.coverage = seeded.coverage
    fallback(reason) = (
        ws.reason = reason;
        @info("warm start: full warm-up", id = ws.id, reason);
        nothing
    )
    if seeded.coverage < WARM_START_MIN_COVERAGE
        pct = round(100 * seeded.coverage; digits = 1)
        return fallback("only $pct% of the metric maps from the parent")
    end
    issuccess(cholesky(Symmetric(seeded.inverse_metric); check = false)) ||
        return fallback("seeded metric is not positive definite")
    ## A parent sampled with a fixed dense metric adapted its step size for
    ## a metric like this one. A parent's diagonal metric gives no guide.
    step_size = ws.parent.metric_kind === :fixed_dense ?
        ws.parent.step_size : 0.0
    ws.path = "seeded"
    @info "warm start: seeded" id = ws.id reason = ws.reason coverage =
        seeded.coverage step_size n_adapts = warm_start_adapts()
    return (;
        inverse_metric = seeded.inverse_metric, step_size,
        n_adapts = warm_start_adapts(),
    )
end

function BVDOutbreakSize.warm_start_record(ws::FitWarmStart, model, chn, alg)
    try
        ld = linked_draws(model, chn)
        steps = chn[BVDOutbreakSize.FlexiChains.Extra(:step_size)]
        step = filter(isfinite, Float64.(vec(collect(steps))))
        gate = if isdefined(Main, :fit_convergence)
            max_depth = fit_max_depth(ws.id)
            Main.fit_convergence(ws.id, chn; max_depth).status
        else
            :unknown
        end
        fixed = alg isa FixedMetricNUTS
        state = (;
            schema = WARM_START_SCHEMA, id = ws.id, key = ws.key,
            models_hash = warm_start_models_hash(), gate,
            metric_kind = fixed ? :fixed_dense : :adapted_diag,
            step_size = isempty(step) ? 0.0 : mean(step),
            layout = ld.layout, draws = ld.draws,
            chain_index = ld.chain_index, warm_up = ws.path,
            parent = ws.parent === nothing ? "" : ws.parent.key,
            reason = ws.reason, coverage = ws.coverage,
            traces = fixed ? alg.step_size_traces :
                Vector{NTuple{2, Float64}}[],
            created = string(now()),
        )
        dir = warm_start_out_dir()
        mkpath(dir)
        serialize(warm_start_file(dir, ws.key), state)
        write_warm_start_summary(state)
    catch e
        @warn "warm start: could not record the state" id = ws.id exception =
            (e, catch_backtrace())
    end
    return nothing
end

## The first warm-up iteration after which a chain's averaged step size
## stays within `tol` of its final value, or `nothing` for an empty trace.
function warm_start_settled(trace; tol = 0.1)
    isempty(trace) && return nothing
    final = last(trace)[2]
    k = length(trace)
    while k > 1 && abs(exp(trace[k - 1][2] - final) - 1) <= tol
        k -= 1
    end
    return k
end

## The dual-averaging trace in brief: each chain's averaged step size at a
## few warm-up iterations, as a share of its final value, and the iteration
## from which it stays within 10% of that value.
function warm_start_trace_table(traces; at = (25, 50, 100, 150, 200, 300))
    isempty(traces) && return ""
    io = IOBuffer()
    marks = [k for k in at if k <= maximum(length, traces)]
    println(
        io, "| chain | ", join(string.(marks), " | "),
        " | settled by | final step size |"
    )
    println(io, "|---", repeat(" | ---", length(marks) + 2), " |")
    for (c, tr) in enumerate(traces)
        isempty(tr) && continue
        final = last(tr)[2]
        share(k) = string(round(exp(tr[k][2] - final); digits = 2))
        cells = [k <= length(tr) ? share(k) : "" for k in marks]
        println(
            io, "| ", c, " | ", join(cells, " | "), " | ",
            warm_start_settled(tr), " | ", round(exp(final); sigdigits = 3),
            " |"
        )
    end
    return String(take!(io))
end

function write_warm_start_summary(state)
    io = IOBuffer()
    println(io, "### Warm start of `", state.id, "`\n")
    println(io, "- Warm-up: ", state.warm_up, ", ", state.reason, ".")
    state.warm_up == "seeded" && println(
        io, "- Metric entries mapped from the parent: ",
        round(100 * state.coverage; digits = 1), "%."
    )
    println(
        io, "- Convergence verdict recorded for the next fit: `",
        state.gate, "`."
    )
    t = warm_start_trace_table(state.traces)
    if !isempty(t)
        println(
            io, "\nAveraged step size over warm-up, as a share of its ",
            "final value:\n"
        )
        print(io, t)
    end
    md = String(take!(io))
    print(stdout, md)
    path = get(ENV, "GITHUB_STEP_SUMMARY", "")
    isempty(path) || open(f -> print(f, md), path, "a")
    return nothing
end
