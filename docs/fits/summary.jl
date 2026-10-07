# Short fit report for the per-fit CI matrix: the convergence diagnostics and
# a couple of headline posteriors for one fit, written to the GitHub Actions
# job summary (the file named by `GITHUB_STEP_SUMMARY`) and to stdout. A fit
# job then says whether its chain converged without anyone opening the log.

using BVDOutbreakSize: fit_diagnostics, posterior_summary,
    sampler_by_chain_table
using Statistics: median

include(joinpath(@__DIR__, "shared.jl"))

## Headline quantities reported when the chain carries them: the cumulative
## infections to the cut-off and the reproduction number there.
const SUMMARY_PARAMETERS = (:C_T, :R_T)


## Median and 90% credible interval of one quantity, or `nothing` when the
## chain does not carry it (a single-stream fit lacks some of them).
function _headline(chn, p::Symbol)
    draws = try
        vec(collect(chn[p]))
    catch
        return nothing
    end
    s = posterior_summary(draws)
    return (median(draws), s.lo90, s.hi90)
end

## One row per chain: the adapted step size, the share of draws whose tree
## reached the depth cap, the mean acceptance statistic and the E-BFMI. Left
## out for a chain that carries no NUTS statistics.
function _chain_rows(io, chn, max_depth)
    t = sampler_by_chain_table(chn; max_depth)
    all(isnan, t.mean_acceptance) && return nothing
    println(io)
    println(
        io,
        "| chain | step size | % at max depth $max_depth | mean acceptance | E-BFMI |"
    )
    println(io, "| --- | --- | --- | --- | --- |")
    for r in eachrow(t)
        println(
            io, "| ", r.chain, " | ", fmt_value(r.step_size),
            " | ", fmt_value(r.percent_at_max_depth),
            " | ", fmt_value(r.mean_acceptance),
            " | ", fmt_value(r.ebfmi), " |"
        )
    end
    return nothing
end

"""
    fit_summary_markdown(id, chn; max_depth = 10) -> String

Markdown block for fit `id`: the `fit_diagnostics` R-hat, effective sample
sizes and divergence count, the per-chain step size, share of draws at the
tree-depth cap `max_depth`, mean acceptance and E-BFMI, then the median and
90% credible interval of each quantity in `SUMMARY_PARAMETERS` the chain
carries.
"""
function fit_summary_markdown(id, chn; max_depth::Integer = 10)
    d = fit_diagnostics(chn)
    io = IOBuffer()
    println(io, "### Fit `", id, "`")
    println(io)
    println(io, "| max R-hat | min ESS bulk | min ESS tail | divergences |")
    println(io, "| --- | --- | --- | --- |")
    println(
        io, "| ", fmt_value(d.max_rhat),
        " | ", fmt_count(d.min_ess_bulk),
        " | ", fmt_count(d.min_ess_tail),
        " | ", d.n_divergent, " |"
    )
    _chain_rows(io, chn, max_depth)
    rows = Any[]
    for p in SUMMARY_PARAMETERS
        h = _headline(chn, p)
        h === nothing || push!(rows, (p, h))
    end
    if !isempty(rows)
        println(io)
        println(io, "| quantity | median | 90% credible interval |")
        println(io, "| --- | --- | --- |")
        for (p, (m, lo, hi)) in rows
            println(
                io, "| `", p, "` | ", fmt_value(m),
                " | ", fmt_value(lo), " to ", fmt_value(hi), " |"
            )
        end
    end
    println(io)
    return String(take!(io))
end

"""
    write_fit_summary(id, result; max_depth = 10) -> Nothing

Print the markdown summary of fit `id` to stdout and append it to the GitHub
Actions job summary when `GITHUB_STEP_SUMMARY` is set. A reporting failure is
warned about rather than thrown: the fit is already cached by this point and
is too expensive to lose to its summary.
"""
function write_fit_summary(id, result; max_depth::Integer = 10)
    try
        md = fit_summary_markdown(id, fit_chain(result); max_depth)
        print(stdout, md)
        path = get(ENV, "GITHUB_STEP_SUMMARY", "")
        isempty(path) || open(io -> print(io, md), path, "a")
    catch e
        @warn "could not summarise fit" id = id exception = e
    end
    return nothing
end

## --- Diagnostics bundle ---------------------------------------------------

using BVDOutbreakSize: parameter_diagnostics, sampler_by_chain_table
using DataFrames: DataFrame, nrow, names, eachrow
using Serialization: serialize

## A plain CSV writer, so the bundle needs no package the test environment
## lacks: header then one row per data row, strings quoted when they hold
## a comma or a quote, `missing` empty.
function _write_csv(path, df::DataFrame)
    cell(x) = ismissing(x) ? "" :
        x isa AbstractString && (occursin(',', x) || occursin('"', x)) ?
        '"' * replace(x, '"' => "\"\"") * '"' : string(x)
    open(path, "w") do io
        println(io, join(names(df), ","))
        for row in eachrow(df)
            println(io, join((cell(row[c]) for c in names(df)), ","))
        end
    end
    return path
end

"""
    write_fit_bundle(id, chn, dir; max_depth = 10) -> String

Write the diagnostics an agent needs to judge fit `id` without loading a
chain: `parameters.csv` (one row per scalar parameter element with its R-hat
and effective sample sizes, [`parameter_diagnostics`](@ref)),
`sampler_by_chain.csv` (draws, divergences, step size, deepest tree, share
of draws at the depth cap `max_depth`, mean acceptance and E-BFMI per chain,
[`sampler_by_chain_table`](@ref)) and `summary.csv` (one row: the
headline diagnostics of [`fit_diagnostics`](@ref) with the draw, chain and
parameter counts). Returns `dir`.
"""
function write_fit_bundle(id, chn, dir; max_depth::Integer = 10)
    mkpath(dir)
    params = parameter_diagnostics(chn)
    _write_csv(joinpath(dir, "parameters.csv"), params)
    _write_csv(joinpath(dir, "sampler_by_chain.csv"), sampler_by_chain_table(chn; max_depth))
    d = fit_diagnostics(chn)
    nd, nc = size(chn)
    _write_csv(
        joinpath(dir, "summary.csv"),
        DataFrame(
            fit = [string(id)], max_rhat = [d.max_rhat],
            min_ess_bulk = [d.min_ess_bulk], min_ess_tail = [d.min_ess_tail],
            divergences = [d.n_divergent], draws = [nd], chains = [nc],
            parameters = [nrow(params)]
        )
    )
    return dir
end

"""
    write_fit_extras(id, key, result, cache_dir; max_depth = 10) -> Nothing

Write, next to the cached fit `<cache_dir>/<key>.jls`, the diagnostics
bundle under `<key>_diagnostics/` and, for a chain the zone stage can meld
from (one carrying `infections_patch`), the parent extract
`<key>.parent.jls` ([`zone_parent_extract`](@ref)). A failure is warned
about rather than thrown: the fit is already cached.
"""
function write_fit_extras(id, key, result, cache_dir; max_depth::Integer = 10)
    chn = fit_chain(result)
    try
        write_fit_bundle(
            id, chn, joinpath(cache_dir, key * "_diagnostics"); max_depth
        )
        if BVDOutbreakSize._has_key(chn, :infections_patch)
            serialize(
                joinpath(cache_dir, key * ".parent.jls"),
                zone_parent_extract(chn; source = key)
            )
        end
    catch e
        @warn "could not write the fit extras" id = id exception = e
    end
    return nothing
end
