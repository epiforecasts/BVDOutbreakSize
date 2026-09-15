# Short fit report for the per-fit CI matrix: the convergence diagnostics and
# a couple of headline posteriors for one fit, written to the GitHub Actions
# job summary (the file named by `GITHUB_STEP_SUMMARY`) and to stdout. A fit
# job then says whether its chain converged without anyone opening the log.

using BVDOutbreakSize: fit_diagnostics, posterior_summary
using Statistics: median

## Headline quantities reported when the chain carries them: the cumulative
## infections to the cut-off and the reproduction number there.
const SUMMARY_PARAMETERS = (:C_T, :R_T)

## The frozen fits return `(; cutoff, o, chn)` rather than the chain itself.
_summary_chain(x) = x isa NamedTuple && haskey(x, :chn) ? x.chn : x

_fmt_count(x) = isfinite(x) ? string(round(Int, x)) : "n/a"
_fmt_value(x) = string(round(x; sigdigits = 3))

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

"""
    fit_summary_markdown(id, chn) -> String

Markdown block for fit `id`: the `fit_diagnostics` R-hat, effective sample
sizes and divergence count, then the median and 90% credible interval of each
quantity in `SUMMARY_PARAMETERS` the chain carries.
"""
function fit_summary_markdown(id, chn)
    d = fit_diagnostics(chn)
    io = IOBuffer()
    println(io, "### Fit `", id, "`")
    println(io)
    println(io, "| max R-hat | min ESS bulk | min ESS tail | divergences |")
    println(io, "| --- | --- | --- | --- |")
    println(io, "| ", _fmt_value(d.max_rhat),
        " | ", _fmt_count(d.min_ess_bulk),
        " | ", _fmt_count(d.min_ess_tail),
        " | ", d.n_divergent, " |")
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
            println(io, "| `", p, "` | ", _fmt_value(m),
                " | ", _fmt_value(lo), " to ", _fmt_value(hi), " |")
        end
    end
    println(io)
    return String(take!(io))
end

"""
    write_fit_summary(id, result) -> Nothing

Print the markdown summary of fit `id` to stdout and append it to the GitHub
Actions job summary when `GITHUB_STEP_SUMMARY` is set. A reporting failure is
warned about rather than thrown: the fit is already cached by this point and
is too expensive to lose to its summary.
"""
function write_fit_summary(id, result)
    try
        md = fit_summary_markdown(id, _summary_chain(result))
        print(stdout, md)
        path = get(ENV, "GITHUB_STEP_SUMMARY", "")
        isempty(path) || open(io -> print(io, md), path, "a")
    catch e
        @warn "could not summarise fit" id=id exception=e
    end
    return nothing
end
