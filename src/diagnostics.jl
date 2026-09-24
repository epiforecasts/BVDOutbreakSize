# Per-parameter fit diagnostics. `summaries.jl` reduces a fit to one worst
# R-hat, one smallest effective sample size and a divergence count. That is
# enough to see that a fit went wrong and never enough to see where. The
# functions here break those single numbers down: one row per scalar
# parameter element, the groups those elements fall into, where the
# divergent transitions sit in parameter space, and the same parameter
# compared between a reference fit and the fits it is assembled from.

## --- Per-parameter diagnostics ------------------------------------------

# A diagnostic entry is usable only when it is a finite number. A fixed or
# degenerate quantity carries an undefined R-hat or effective sample size,
# which would otherwise rank ahead of every genuinely badly mixing
# parameter.
_diag_value(x) = (x isa Number && isfinite(x)) ? Float64(x) : NaN

# Smallest finite entry, or `NaN` when every entry is undefined.
function _min_finite(v)
    f = filter(!isnan, v)
    return isempty(f) ? NaN : minimum(f)
end

# Largest finite entry, or `NaN` when every entry is undefined.
function _max_finite(v)
    f = filter(!isnan, v)
    return isempty(f) ? NaN : maximum(f)
end

# `i`th entry of a flattened diagnostic vector, or `NaN` past its end. The
# R-hat, bulk and tail summaries are built independently, so one can carry
# fewer entries than another for the same parameter.
_at(v, i) = i <= length(v) ? _diag_value(v[i]) : NaN

# Split a diagnostic name into its parameter and its element index. A
# vector-valued parameter is summarised as one indexed scalar per element,
# so the index arrives inside the name (`rt_state.log_R[13]`). A multi-index
# name keeps its first index, which is the axis the elements run along.
function _split_index(name::AbstractString)
    m = match(r"^(.*)\[([0-9]+(?:\s*,\s*[0-9]+)*)\]$", name)
    isnothing(m) && return (String(name), 0)
    return (
        String(m.captures[1]),
        parse(Int, strip(first(split(m.captures[2], ",")))),
    )
end

"""
Per-parameter fit diagnostics for one chain, one row per scalar parameter
element. A vector-valued parameter contributes one row per element, numbered
in `index` from one, and `parameter` carries the name without the index; a
scalar parameter carries index zero. Trajectories named in `exclude` are
skipped, as they are in the headline [`fit_diagnostics`](@ref) summary.

Columns are `:parameter`, `:index`, `:rhat`, `:ess_bulk` and `:ess_tail`.
Rows whose R-hat is undefined are dropped, so a fixed or degenerate quantity
cannot rank ahead of a parameter that genuinely mixes badly.
"""
function parameter_diagnostics(chn; exclude = _DIAGNOSTIC_EXCLUDE)
    rh = FlexiChains.rhat(chn)
    eb = FlexiChains.ess(chn; kind = :bulk)
    et = FlexiChains.ess(chn; kind = :tail)
    parameter = String[]
    index = Int[]
    rhat = Float64[]
    ess_bulk = Float64[]
    ess_tail = Float64[]
    for p in FlexiChains.parameters(rh)
        name = string(p)
        any(b -> startswith(name, b), exclude) && continue
        r, b, t = rh[p], eb[p], et[p]
        if r isa Number
            base, idx = _split_index(name)
            push!(parameter, base)
            push!(index, idx)
            push!(rhat, _diag_value(r))
            push!(ess_bulk, _diag_value(b))
            push!(ess_tail, _diag_value(t))
        else
            rv = vec(collect(r))
            bv = vec(collect(b))
            tv = vec(collect(t))
            for i in eachindex(rv)
                push!(parameter, name)
                push!(index, i)
                push!(rhat, _diag_value(rv[i]))
                push!(ess_bulk, _at(bv, i))
                push!(ess_tail, _at(tv, i))
            end
        end
    end
    df = DataFrame(; parameter, index, rhat, ess_bulk, ess_tail)
    return df[.!isnan.(df.rhat), :]
end

# Reader-facing name for one diagnostics row: the display label where the
# parameter has one, carrying the element index for a vector-valued
# parameter.
function _diag_label(parameter, index; labels = Dict{Symbol, String}())
    base = get(labels, Symbol(parameter), parameter)
    return index == 0 ? String(base) : string(base, "[", index, "]")
end

# Mask of the rows of `keys` that are the first to carry their combination of
# values. A model that names the same quantity twice, once sampled and once
# as a deterministic copy, gives both names identical diagnostics, and a
# ranked table then spends half its rows saying the same thing.
function _first_occurrence(keys)
    seen = Set{eltype(keys)}()
    return [k in seen ? false : (push!(seen, k); true) for k in keys]
end

# Diagnostics frame for `x`, which is either a chain or a frame
# `parameter_diagnostics` has already produced. R-hat and both effective
# sample sizes over several thousand parameters are not free to compute, so
# a page showing several views of one fit builds the frame once and hands it
# to each of them.
_as_diagnostics(x) = x isa DataFrame ? x : parameter_diagnostics(x)

"""
The `n` worst-mixing parameter elements of a fit, ranked by bulk effective
sample size from the lowest up. `fit` is a chain or a frame
[`parameter_diagnostics`](@ref) has already produced. `labels` maps a
parameter name to the display name the report uses for it. With
`collapse_aliases` a parameter whose diagnostics repeat those of one already
listed is dropped, since it is the same quantity under a second name.

Columns are `:parameter`, `:rhat`, `:ess_bulk` and `:ess_tail`.
"""
function worst_parameters_table(
        fit; n::Integer = 15,
        collapse_aliases::Bool = true, labels = Dict{Symbol, String}()
    )
    df = _as_diagnostics(fit)
    isempty(df) && return DataFrame(
        parameter = String[], rhat = Float64[],
        ess_bulk = Float64[], ess_tail = Float64[]
    )
    ## A stable sort so two parameters that tie keep model order and
    ## the published table does not reshuffle between builds.
    ranked = df[sortperm(df.ess_bulk; alg = MergeSort), :]
    if collapse_aliases
        ranked = ranked[
            _first_occurrence(
                collect(
                    zip(
                        ranked.rhat,
                        ranked.ess_bulk, ranked.ess_tail
                    )
                )
            ), :,
        ]
    end
    keep = first(ranked, n)
    return DataFrame(
        parameter = [
            _diag_label(r.parameter, r.index; labels = labels)
                for r in eachrow(keep)
        ],
        rhat = round.(keep.rhat; digits = 3),
        ess_bulk = round.(keep.ess_bulk; digits = 0),
        ess_tail = round.(keep.ess_tail; digits = 0)
    )
end

"""
Fit diagnostics grouped by parameter, one row per parameter name rather than
per element, ranked by bulk effective sample size from the lowest up. A
vector-valued parameter collapses to a single row, so a group whose elements
all mix badly is visible as one line rather than as hundreds.

`rhat_threshold` sets the R-hat an element has to exceed to be counted in
the last column, and `collapse_aliases` drops a parameter whose diagnostics
repeat those of one already listed. Columns are `:parameter`, `:elements`,
`:max_rhat`, `:min_ess_bulk` and the count above the threshold.
"""
function family_diagnostics_table(
        fit; n::Integer = 10,
        rhat_threshold::Real = 1.1, collapse_aliases::Bool = true,
        labels = Dict{Symbol, String}()
    )
    df = _as_diagnostics(fit)
    groups = unique(df.parameter)
    rows = [df.parameter .== g for g in groups]
    tbl = DataFrame(
        "parameter" => [String(get(labels, Symbol(g), g)) for g in groups],
        "elements" => [count(r) for r in rows],
        "max_rhat" => [
            round(_max_finite(df.rhat[r]); digits = 3)
                for r in rows
        ],
        "min_ess_bulk" => [
            round(_min_finite(df.ess_bulk[r]); digits = 0)
                for r in rows
        ],
        "above_$(rhat_threshold)" => [
            count(>(rhat_threshold), df.rhat[r])
                for r in rows
        ]
    )
    sort!(tbl, :min_ess_bulk; alg = MergeSort)
    if collapse_aliases
        tbl = tbl[
            _first_occurrence(
                collect(
                    zip(
                        tbl.elements, tbl.max_rhat,
                        tbl.min_ess_bulk
                    )
                )
            ), :,
        ]
    end
    return first(tbl, n)
end

"""
How far R-hat and the effective sample size spread across the parameters of
each fit, one row per fit. Pass each fit as `"label" => chain`, or as
`"label" => frame` where the frame is one [`parameter_diagnostics`](@ref)
has already produced.

A fit can carry a bad worst-case R-hat because one parameter is stuck, or
because most of the model is. The counts here separate those two cases. The
last column names the parameter whose elements reach the lowest bulk
effective sample size.
"""
function diagnostic_spread_table(
        fits::Pair{String}...;
        rhat_warn::Real = 1.01, rhat_bad::Real = 1.1,
        ess_low::Integer = 100, labels = Dict{Symbol, String}()
    )
    dfs = [_as_diagnostics(f.second) for f in fits]
    worst = map(dfs) do df
        isempty(df) && return "none"
        groups = unique(df.parameter)
        mins = [_min_finite(df.ess_bulk[df.parameter .== g]) for g in groups]
        g = groups[argmin(replace(mins, NaN => Inf))]
        return String(get(labels, Symbol(g), g))
    end
    return DataFrame(
        "fit" => [String(f.first) for f in fits],
        "parameters" => [nrow(df) for df in dfs],
        "rhat_above_$(rhat_warn)" => [
            count(>(rhat_warn), df.rhat)
                for df in dfs
        ],
        "rhat_above_$(rhat_bad)" => [
            count(>(rhat_bad), df.rhat)
                for df in dfs
        ],
        "percent_above_$(rhat_bad)" =>
            [
            round(
                100 * count(>(rhat_bad), df.rhat) / max(nrow(df), 1);
                digits = 1
            ) for df in dfs
        ],
        "ess_bulk_below_$(ess_low)" =>
            [count(x -> !isnan(x) && x < ess_low, df.ess_bulk) for df in dfs],
        "lowest_ess_parameter" => worst
    )
end

## --- Divergent transitions ----------------------------------------------

# Values of the sampler statistic `name`, or `nothing` when the chain does
# not carry it. Chains fitted by other samplers, and prior draws, carry no
# NUTS statistics at all.
function _extra_values(chn, name::Symbol)
    for e in FlexiChains.extras(chn)
        e.name === name || continue
        return chn[e]
    end
    return nothing
end

# Divergence flag per draw, flattened the same way parameter draws are, so
# the two line up entry for entry.
function _divergent_flags(chn)
    d = _extra_values(chn, :numerical_error)
    isnothing(d) && return Bool[]
    return [x === true for x in vec(collect(d))]
end

"""
Sampler behaviour per chain: how many draws each chain took, how many of
them were divergent, the step size it adapted to and the deepest tree it
built.

One chain with a far smaller step size and far deeper trees than the others
is a chain stuck somewhere the rest of the posterior never visits, which is
a different problem from divergences spread evenly across chains.
"""
function sampler_by_chain_table(chn)
    nd, nc = size(chn)
    ## Flatten each statistic the way parameter draws flatten, then read one
    ## chain's block out of it, so a chain's draws stay together.
    flat(name) =
    let x = _extra_values(chn, name)
        isnothing(x) ? nothing : vec(collect(x))
    end
    block(v, c) = v[((c - 1) * nd + 1):(c * nd)]
    div = flat(:numerical_error)
    step = flat(:step_size)
    depth = flat(:tree_depth)
    divergences = [
        isnothing(div) ? 0 :
            count(x -> x === true, block(div, c)) for c in 1:nc
    ]
    steps = [
        isnothing(step) ? NaN :
            _min_finite(Float64.(block(step, c))) for c in 1:nc
    ]
    depths = [
        isnothing(depth) ? 0 :
            round(Int, _max_finite(Float64.(block(depth, c))))
            for c in 1:nc
    ]
    return DataFrame(
        chain = collect(1:nc),
        draws = fill(nd, nc),
        divergences = divergences,
        percent_divergent = round.(100 .* divergences ./ nd; digits = 1),
        step_size = round.(steps; sigdigits = 3),
        deepest_tree = depths
    )
end

# Middle `width` interval of `x` as a printable range.
function _interval_string(x; width::Real = 0.9, digits::Integer = 3)
    lo = quantile(x, (1 - width) / 2)
    hi = quantile(x, 1 - (1 - width) / 2)
    return string(
        round(lo; sigdigits = digits), "–",
        round(hi; sigdigits = digits)
    )
end

"""
Where the divergent transitions of `chn` sit in parameter space, one row per
scalar parameter, ranked by separation from the largest down.

`separation` is how far the divergent draws sit from the posterior as a
whole, in standard deviations of the full posterior, signed by direction. A
separation near zero means divergences scattered across the posterior; a
large one means they concentrate in one region of that parameter.

The two interval columns give the middle `width` of all draws and of the
divergent draws alone. Vector-valued parameters are skipped, as are
parameters that never move, and `collapse_aliases` drops a parameter whose
intervals and separation repeat those of one already listed.
"""
function divergence_location_table(
        chn; n::Integer = 10, width::Real = 0.9,
        collapse_aliases::Bool = true, exclude = _DIAGNOSTIC_EXCLUDE,
        labels = Dict{Symbol, String}()
    )
    flag = _divergent_flags(chn)
    empty = DataFrame(
        parameter = String[], all_draws = String[],
        divergent_draws = String[], separation = Float64[]
    )
    any(flag) || return empty
    names_ = String[]
    alls = String[]
    divs = String[]
    seps = Float64[]
    for p in FlexiChains.parameters(chn)
        name = string(p)
        any(b -> startswith(name, b), exclude) && continue
        v = chn[p]
        x = vec(collect(v))
        eltype(x) <: Number || continue
        length(x) == length(flag) || continue
        xs = Float64.(x)
        all(isfinite, xs) || continue
        s = std(xs)
        (isfinite(s) && s > 0) || continue
        push!(names_, get(labels, Symbol(name), name))
        push!(alls, _interval_string(xs; width = width))
        push!(divs, _interval_string(xs[flag]; width = width))
        push!(seps, round((median(xs[flag]) - median(xs)) / s; digits = 2))
    end
    isempty(names_) && return empty
    tbl = DataFrame(
        parameter = names_, all_draws = alls,
        divergent_draws = divs, separation = seps
    )
    sort!(tbl, :separation; by = abs, rev = true, alg = MergeSort)
    if collapse_aliases
        tbl = tbl[
            _first_occurrence(
                collect(
                    zip(
                        tbl.all_draws,
                        tbl.divergent_draws, tbl.separation
                    )
                )
            ), :,
        ]
    end
    return first(tbl, n)
end

## --- Contrast between fits ----------------------------------------------

"""
Bulk effective sample size for every parameter element shared between a
reference fit and each of the fits it is compared with. Pass the reference
as `"label" => chain` first, then each comparison fit the same way. A frame
[`parameter_diagnostics`](@ref) has already produced stands in for a chain
on either side.

`ess_ratio` is the reference fit's bulk effective sample size over the
comparison fit's. A ratio well below one is a parameter that mixes in the
comparison fit and stops mixing in the reference, which points at what the
reference adds rather than at the parameter itself.

Columns are `:fit`, `:parameter`, `:index`, `:rhat`, `:ess_bulk`,
`:rhat_reference`, `:ess_bulk_reference` and `:ess_ratio`.
"""
function diagnostic_contrast(reference::Pair{String}, fits::Pair{String}...)
    ref = _as_diagnostics(reference.second)
    lookup = Dict((r.parameter, r.index) => r for r in eachrow(ref))
    fit = String[]
    parameter = String[]
    index = Int[]
    rhat = Float64[]
    ess_bulk = Float64[]
    rhat_reference = Float64[]
    ess_bulk_reference = Float64[]
    ess_ratio = Float64[]
    for (label, chn) in fits
        for r in eachrow(_as_diagnostics(chn))
            q = get(lookup, (r.parameter, r.index), nothing)
            isnothing(q) && continue
            (isnan(r.ess_bulk) || isnan(q.ess_bulk)) && continue
            r.ess_bulk > 0 || continue
            push!(fit, String(label))
            push!(parameter, r.parameter)
            push!(index, r.index)
            push!(rhat, r.rhat)
            push!(ess_bulk, r.ess_bulk)
            push!(rhat_reference, q.rhat)
            push!(ess_bulk_reference, q.ess_bulk)
            push!(ess_ratio, q.ess_bulk / r.ess_bulk)
        end
    end
    return DataFrame(;
        fit, parameter, index, rhat, ess_bulk,
        rhat_reference, ess_bulk_reference, ess_ratio
    )
end

"""
The `n` parameter elements whose mixing degrades most from a comparison fit
to the reference fit, taken from [`diagnostic_contrast`](@ref) and ranked by
the ratio of the two bulk effective sample sizes from the smallest up.
`collapse_aliases` drops a parameter whose numbers repeat those of one
already listed for the same fit.

Columns are `:fit`, `:parameter`, `:ess_bulk`, `:ess_bulk_reference` and
`:ess_ratio`.
"""
function diagnostic_contrast_table(
        df::DataFrame; n::Integer = 15,
        collapse_aliases::Bool = true, labels = Dict{Symbol, String}()
    )
    isempty(df) && return DataFrame(
        fit = String[], parameter = String[],
        ess_bulk = Float64[], ess_bulk_reference = Float64[],
        ess_ratio = Float64[]
    )
    ranked = sort(df, :ess_ratio; alg = MergeSort)
    if collapse_aliases
        ranked = ranked[
            _first_occurrence(
                collect(
                    zip(
                        ranked.fit,
                        ranked.ess_bulk, ranked.ess_bulk_reference
                    )
                )
            ), :,
        ]
    end
    keep = first(ranked, n)
    return DataFrame(
        fit = keep.fit,
        parameter = [
            _diag_label(r.parameter, r.index; labels = labels)
                for r in eachrow(keep)
        ],
        ess_bulk = round.(keep.ess_bulk; digits = 0),
        ess_bulk_reference = round.(keep.ess_bulk_reference; digits = 0),
        ess_ratio = round.(keep.ess_ratio; digits = 2)
    )
end
