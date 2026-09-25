# Convergence gate for the report's fits. `summary.jl` writes one fit's
# diagnostics to the job summary; this decides whether those diagnostics are
# good enough to publish a report from, and names the parameters responsible
# when they are not. The joint fit is the one the headline outbreak size comes
# from, so a build that publishes a joint fit which has not converged publishes
# a number nobody should read.
#
# Two thresholds per diagnostic rather than one. The fail threshold marks a
# fit that has not converged at all (an R-hat of 2.6 is a chain stuck at its
# initialisation, not a slow one). The warn threshold marks a fit that has
# converged but mixes poorly, which is where the joint fit currently sits and
# is tracked work rather than a build failure. Both are environment-settable,
# so tightening the gate as the fit improves is a variable, not a code change.

using BVDOutbreakSize: fit_diagnostics, parameter_diagnostics,
    worst_parameters_table, family_diagnostics_table,
    divergence_location_table, markdown_table

include(joinpath(@__DIR__, "shared.jl"))

## Fails a build. An R-hat this far from one, an effective sample size this
## small or a divergence rate this high is a sampler that has not explored the
## posterior, whatever the numbers it reports.
const CONVERGENCE_FAIL = (
    rhat = 1.1, ess_bulk = 25.0, ess_tail = 25.0, divergent_fraction = 0.05,
)

## Reported, but does not fail a build: converged and mixing badly.
const CONVERGENCE_WARN = (
    rhat = 1.05, ess_bulk = 100.0, ess_tail = 100.0,
    divergent_fraction = 0.01,
)

## The fits the gate covers by default. The joint fit carries the headline
## estimate and the health-zone fit the split of it that the zone maps and
## the zone forecast report. The single-stream fits are diagnostics of the
## joint and are allowed to be worse.
const CONVERGENCE_IDS = ["joint", "local"]

## A threshold read from the environment. A malformed value names the
## variable it came from: the parse failure alone says only that some number
## could not be read, and the report is never written to say which.
##
## `nan` and `inf` parse, and either one turns the threshold off without
## saying so, because every comparison against them is false. A gate that a
## typo can silently disable is the failure this exists to prevent, so a
## threshold has to be finite.
function _env_number(name, default)
    v = strip(get(ENV, name, ""))
    isempty(v) && return default
    n = tryparse(Float64, v)
    n === nothing && error("$name is \"$v\", which is not a number")
    isfinite(n) ||
        error("$name is \"$v\"; a threshold has to be a finite number")
    return n
end

"""
    convergence_thresholds() -> (; fail, warn)

The fail and warn thresholds the gate applies, each overridable by an
environment variable: `BVD_CONVERGENCE_FAIL_RHAT`,
`BVD_CONVERGENCE_FAIL_ESS_BULK`, `BVD_CONVERGENCE_FAIL_ESS_TAIL`,
`BVD_CONVERGENCE_FAIL_DIVERGENT` and the matching `BVD_CONVERGENCE_WARN_*`.
"""
function convergence_thresholds()
    each(prefix, defaults) = (;
        rhat = _env_number("BVD_CONVERGENCE_$(prefix)_RHAT", defaults.rhat),
        ess_bulk = _env_number(
            "BVD_CONVERGENCE_$(prefix)_ESS_BULK", defaults.ess_bulk
        ),
        ess_tail = _env_number(
            "BVD_CONVERGENCE_$(prefix)_ESS_TAIL", defaults.ess_tail
        ),
        divergent_fraction = _env_number(
            "BVD_CONVERGENCE_$(prefix)_DIVERGENT", defaults.divergent_fraction
        ),
    )
    return (;
        fail = each("FAIL", CONVERGENCE_FAIL),
        warn = each("WARN", CONVERGENCE_WARN),
    )
end


## A diagnostic that could not be computed (every parameter degenerate, or a
## chain carrying no sampler extras) is not evidence of convergence, but it is
## not evidence against it either, so it neither fails nor warns. The report
## prints it as `n/a`, which is visible.
_breaches(value, limit, worse) = isfinite(value) && worse(value, limit)

"""
    convergence_verdict(d; thresholds = convergence_thresholds())

Decide whether the headline diagnostics `d` (as [`fit_diagnostics`](@ref)
returns them) clear the gate. Returns `(; status, failures, warnings,
divergent_fraction)`, where `status` is `:pass`, `:warn` or `:fail` and the
two message vectors say which diagnostic breached which threshold.
"""
function convergence_verdict(d; thresholds = convergence_thresholds())
    fraction = d.n_draws > 0 ? d.n_divergent / d.n_draws : NaN
    checks = (
        (
            "max R-hat", d.max_rhat, >,
            thresholds.fail.rhat, thresholds.warn.rhat, fmt_value,
        ),
        (
            "min bulk ESS", d.min_ess_bulk, <,
            thresholds.fail.ess_bulk, thresholds.warn.ess_bulk, fmt_count,
        ),
        (
            "min tail ESS", d.min_ess_tail, <,
            thresholds.fail.ess_tail, thresholds.warn.ess_tail, fmt_count,
        ),
        (
            "divergent fraction", fraction, >,
            thresholds.fail.divergent_fraction,
            thresholds.warn.divergent_fraction, fmt_value,
        ),
    )
    failures = String[]
    warnings = String[]
    for (name, value, worse, fail, warn, fmt) in checks
        if _breaches(value, fail, worse)
            push!(
                failures,
                "$name is $(fmt(value)), past the failure threshold " *
                    "$(fmt(fail))"
            )
        elseif _breaches(value, warn, worse)
            push!(
                warnings,
                "$name is $(fmt(value)), past the warning threshold " *
                    "$(fmt(warn))"
            )
        end
    end
    status = !isempty(failures) ? :fail : (isempty(warnings) ? :pass : :warn)
    return (; status, failures, warnings, divergent_fraction = fraction)
end

const _STATUS_PHRASE = Dict(
    :pass => "passes the convergence checks",
    :warn => "passes the convergence checks with warnings",
    :fail => "fails the convergence checks",
)

## A fraction as a percentage, with no trailing `.0` on a whole number.
function _percent(x)
    p = round(100 * x; sigdigits = 3)
    return isinteger(p) ? string(Int(p)) : string(p)
end

"""
    convergence_summary(label, d; thresholds = convergence_thresholds())

A report reader's version of [`convergence_verdict`](@ref) for the headline
diagnostics `d` of the fit named `label`, as plain prose. `verdict` is the
one-sentence verdict. `detail` gives the values behind it and what a fit
needs to pass without warnings.
"""
function convergence_summary(label, d; thresholds = convergence_thresholds())
    v = convergence_verdict(d; thresholds = thresholds)
    w = thresholds.warn
    verdict = "The $label fit $(_STATUS_PHRASE[v.status])."
    detail = "Its worst R-hat is $(fmt_value(d.max_rhat)) and its lowest " *
        "bulk effective sample size is $(fmt_count(d.min_ess_bulk)), with " *
        "$(d.n_divergent) divergent transitions in $(d.n_draws) draws. " *
        "A fit passes without warnings when its worst R-hat is at most " *
        "$(fmt_value(w.rhat)), its effective sample sizes are at least " *
        "$(fmt_count(w.ess_bulk)) and at most " *
        "$(_percent(w.divergent_fraction))% of its draws diverge."
    return (; verdict, detail)
end

"""
    fit_convergence(id, chn; thresholds = convergence_thresholds())

Gate one fit. Returns the verdict merged with `(; id, diagnostics)`, and,
unless the fit passed, the per-parameter frame the report's diagnostics
tables are built from, so the caller can say which parameters are at fault
without recomputing several thousand R-hats.
"""
function fit_convergence(id, chn; thresholds = convergence_thresholds())
    d = fit_diagnostics(chn)
    v = convergence_verdict(d; thresholds = thresholds)
    per_parameter = v.status === :pass ? nothing : parameter_diagnostics(chn)
    return (; id = String(id), diagnostics = d, per_parameter, chain = chn, v...)
end

const _STATUS_MARK = Dict(
    :pass => "✅ passed", :warn => "⚠️ passed with warnings",
    :fail => "❌ failed",
)

## The worst status across the gated fits: one failure fails the build.
function _overall(checks)
    any(c -> c.status === :fail, checks) && return :fail
    any(c -> c.status === :warn, checks) && return :warn
    return :pass
end

## Headline row per fit, the same four diagnostics the report's fit-diagnostics
## table carries plus the draw count the divergences are out of.
function _headline_table(io, checks)
    println(
        io,
        "| fit | max R-hat | min ESS bulk | min ESS tail | divergences | draws |"
    )
    println(io, "| --- | --- | --- | --- | --- | --- |")
    for c in checks
        d = c.diagnostics
        println(
            io, "| `", c.id, "` | ", fmt_value(d.max_rhat),
            " | ", fmt_count(d.min_ess_bulk),
            " | ", fmt_count(d.min_ess_tail),
            " | ", d.n_divergent, " | ", d.n_draws, " |"
        )
    end
    return println(io)
end

## Which parameters are responsible, the same three views the sensitivity
## report's "Fit diagnostics by parameter" section shows, cut to what fits in
## a comment. Ranked by bulk effective sample size from the lowest up, then
## grouped so a whole badly mixing vector reads as one line, then where the
## divergent transitions concentrate.
function _culprit_tables(io, c, thresholds; n::Integer = 10)
    c.per_parameter === nothing && return nothing
    isempty(c.per_parameter) && return nothing
    println(
        io, "<details><summary>Worst-mixing parameters of `", c.id,
        "`</summary>\n"
    )
    println(io, markdown_table(worst_parameters_table(c.per_parameter; n = n)))
    println(io, "\nThe same parameters grouped, one row per parameter:\n")
    println(
        io, markdown_table(
            family_diagnostics_table(
                c.per_parameter; n = n,
                rhat_threshold = thresholds.warn.rhat
            )
        )
    )
    if c.diagnostics.n_divergent > 0
        divergences = divergence_location_table(c.chain; n = n)
        if !isempty(divergences)
            println(
                io,
                "\nWhere the divergent transitions sit, in standard " *
                    "deviations of the full posterior:\n"
            )
            println(io, markdown_table(divergences))
        end
    end
    return println(io, "\n</details>\n")
end

"""
    convergence_markdown(checks; thresholds = convergence_thresholds(),
                         marker = "", context = "") -> String

The gate's report: one verdict line, the headline diagnostics table, the
reasons each fit failed or warned and the parameters responsible. `marker` is
prepended verbatim, so a caller posting this as a sticky comment can find its
own comment again. `context` is a line of provenance (the commit and a link to
the run) printed under the verdict.
"""
function convergence_markdown(
        checks; thresholds = convergence_thresholds(),
        marker::AbstractString = "", context::AbstractString = ""
    )
    status = _overall(checks)
    io = IOBuffer()
    isempty(marker) || println(io, marker)
    println(io, "## Fit convergence: ", _STATUS_MARK[status])
    println(io)
    if status === :fail
        println(
            io,
            "The headline fit below has not converged, so the estimates this ",
            "build publishes should not be read. The documentation preview ",
            "is still built and deployed."
        )
    elseif status === :warn
        println(
            io,
            "Every gated fit converged. The diagnostics below are past the ",
            "warning thresholds, so the posterior is thinner than it should ",
            "be, but the build is not failed on it."
        )
    else
        println(io, "Every gated fit cleared the convergence thresholds.")
    end
    println(io)
    isempty(context) || (println(io, context); println(io))
    _headline_table(io, checks)
    for c in checks
        (isempty(c.failures) && isempty(c.warnings)) && continue
        println(io, "### `", c.id, "`: ", _STATUS_MARK[c.status])
        println(io)
        for m in c.failures
            println(io, "- **fail** — ", m)
        end
        for m in c.warnings
            println(io, "- warn — ", m)
        end
        println(io)
        _culprit_tables(io, c, thresholds)
    end
    println(
        io,
        "Thresholds — fail: R-hat above ", fmt_value(thresholds.fail.rhat),
        ", bulk ESS below ", fmt_count(thresholds.fail.ess_bulk),
        ", tail ESS below ", fmt_count(thresholds.fail.ess_tail),
        ", divergences above ", fmt_value(thresholds.fail.divergent_fraction),
        " of draws; warn: ", fmt_value(thresholds.warn.rhat), ", ",
        fmt_count(thresholds.warn.ess_bulk), ", ",
        fmt_count(thresholds.warn.ess_tail), ", ",
        fmt_value(thresholds.warn.divergent_fraction), "."
    )
    println(io)
    println(
        io,
        "The full per-parameter breakdown is on the sensitivity page, under ",
        "[Fit diagnostics by parameter]",
        "(https://epiforecasts.io/BVDOutbreakSize/dev/sensitivity",
        "#Fit-diagnostics-by-parameter)."
    )
    return String(take!(io))
end
