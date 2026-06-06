# Observation loading from the dated TOML manifest. The manifest stores
# calendar dates and cumulative counts (never grid day-indices), so the
# data stay aligned by date as vintages are added, revised or arrive
# sparsely. The renewal model works on a daily grid, so this loader
# derives the grid length and the per-vintage day-indices from the dates
# at load time: the cut-off is the last grid day, and the seeding day (the
# first grid day) is placed a fixed lead before the genetic TMRCA so the
# grid always contains the inferred seeding time.

## Days the grid extends before the genetic TMRCA date, so the seeding
## crossing has room to be inferred below the molecular-clock bound.
const SEEDING_LEAD_DAYS = 30

"""
Load the BVD observation manifest from `path` (a dated TOML file) and
return a named tuple for the renewal model. Calendar dates are converted
to 1-based grid day-indices (day 1 is the seeding day, day `n` the
cut-off) using the cut-off (`as_of_date`) and a seeding day placed
`seeding_lead` days before the genetic TMRCA date.

Returns the grid length `n`, the `cutoff` and `seeding` dates, the
per-stream cumulative totals at the cut-off (`reported_cases`,
`total_deaths`, `confirmed_cases`, `confirmed_deaths`, `tests_analysed`,
`exported_cases`, `exports_deaths`), the per-vintage histories as
`(; days, counts)` with `days` the grid day-indices (`reported_history`,
`confirmed_history`, `confirmed_deaths_history`, `deaths_history`,
`lab_history` from the analysed-specimen series, `tests_received_history`),
the genetic TMRCA bound `tmrca_days` (days before the cut-off), and
`who_first_sitrep_days` (days from the first situation report, the
earliest reported-case vintage, to the cut-off). The intervention
breakpoint grid day is `n - who_first_sitrep_days`. A cut-off scalar with
no explicit TOML block is derived from the final vintage of the matching
history.
"""
function load_observations(
        path::AbstractString = joinpath(@__DIR__, "..", "data",
            "observations.toml");
        seeding_lead::Integer = SEEDING_LEAD_DAYS,
        cutoff_date::Union{Nothing, Date, AbstractString} = nothing)
    raw = TOML.parsefile(path)
    _val(k) = raw[k]["value"]
    ## The cut-off is the manifest `as_of_date` unless an earlier
    ## `cutoff_date` is supplied (used to FREEZE the data to a past
    ## date, see `freeze_observations`). Freezing only ever moves the
    ## cut-off earlier and never invents vintages, so the grid stays
    ## date-aligned with the full-data fit.
    cutoff = isnothing(cutoff_date) ? Date(String(raw["as_of_date"])) :
             (cutoff_date isa Date ? cutoff_date :
              Date(String(cutoff_date)))
    ## Grid day-index (1-based) of a calendar date: seeding day is day 1.
    _gap(d) = Int(date2epochdays(cutoff) - date2epochdays(Date(String(d))))
    tmrca_date = Date(String(raw["genetic_tmrca"]["date"]))
    seeding = tmrca_date - Day(seeding_lead)
    n = Int(date2epochdays(cutoff) - date2epochdays(seeding)) + 1
    _index(d) = n - _gap(d)

    ## A dated cumulative history → grid day-indices and counts, sorted
    ## oldest-first so the model differences consecutive vintages into
    ## daily increments. Empty when the block is absent. Vintages dated
    ## after the cut-off are dropped, so freezing to an earlier date
    ## keeps only the data that was available by then.
    function history(key)
        haskey(raw, key) || return (; days = Int[], counts = Int[])
        block = raw[key]
        keep = [Date(String(d)) <= cutoff for d in block["dates"]]
        idx = Int[_index(d) for d in block["dates"][keep]]
        vals = Int.(block["values"][keep])
        ord = sortperm(idx)
        return (; days = idx[ord], counts = vals[ord])
    end

    reported_history = history("reported_case_history")
    confirmed_history = history("confirmed_case_history")
    confirmed_deaths_history = history("confirmed_death_history")
    deaths_history = history("death_history")
    ## The analysed-specimen series is the laboratory denominator; the
    ## received series is recorded for the pipeline view but not fitted.
    lab_history = history("tests_analysed_history")
    tests_received_history = history("tests_received_history")
    ## Cut-off scalar from an explicit TOML block, else the final
    ## (most recent) vintage of the matching history. When a `cutoff_date`
    ## freeze is active the explicit TOML scalars (which hold the final,
    ## full-data total) no longer match the truncated history, so the
    ## frozen final vintage is used instead.
    _hist_end(h) = isempty(h.counts) ? missing : h.counts[end]
    frozen = !isnothing(cutoff_date)
    _scalar(k, h) = (frozen || !haskey(raw, k)) ? _hist_end(h) : Int(_val(k))
    ## The first WHO joint situation report is the earliest reported-case
    ## vintage; days from it to the cut-off set the intervention breakpoint.
    who_first_sitrep_days = isempty(reported_history.days) ? n :
                            n - reported_history.days[1] + 1

    ## Days from the LAST observed Uganda import detection to the cut-off.
    ## The exported-case scalar counts imports detected by this last import
    ## day; no export was detected in the trailing window to the cut-off, so
    ## the exports likelihood truncates the at-risk person-time sum here
    ## rather than over-expecting exports in that window. Zero when no
    ## dated import series is present (sum runs to the cut-off).
    export_last_offset = if haskey(raw, "export_case_dates")
        dts = [Date(String(d)) for d in raw["export_case_dates"]["value"]]
        kept = filter(d -> d <= cutoff, dts)
        isempty(kept) ? 0 :
        Int(date2epochdays(cutoff) - date2epochdays(maximum(kept)))
    else
        0
    end

    return (; n, cutoff, seeding,
        exported_cases = Int(_val("exported_cases")),
        exports_deaths = Int(_val("exports_deaths")),
        total_deaths = frozen ?
                       _hist_end(deaths_history) : Int(_val("total_deaths")),
        reported_cases = frozen ?
                         _hist_end(reported_history) :
                         Int(_val("reported_cases")),
        confirmed_cases = _scalar("confirmed_cases", confirmed_history),
        confirmed_deaths = _scalar("confirmed_deaths",
            confirmed_deaths_history),
        tests_analysed = _scalar("cumulative_tests_analysed", lab_history),
        reported_history = reported_history,
        confirmed_history = confirmed_history,
        confirmed_deaths_history = confirmed_deaths_history,
        deaths_history = deaths_history,
        lab_history = lab_history,
        tests_received_history = tests_received_history,
        tmrca_days = _gap(raw["genetic_tmrca"]["date"]),
        who_first_sitrep_days, export_last_offset)
end

"""
    freeze_observations(cutoff_date; path = default manifest)

Load the observation manifest FROZEN to `cutoff_date`: the cut-off is
moved to `cutoff_date` and every dated history is truncated to the
vintages available by then, so the returned named tuple is what the
renewal model would have seen on that date. The cut-off scalar totals
(`reported_cases`, `total_deaths`, `confirmed_cases`, ...) are taken
from the truncated histories rather than the manifest's full-data
scalars. Use it to re-evaluate the renewal estimate at a past report
date (for example a McCabe et al. situation-report cut-off) for a
like-for-like, matched-in-time comparison.

`cutoff_date` accepts a `Date` or an ISO date string. It must be on or
after the earliest history vintage in the manifest (the renewal DRC
series begins 18 May 2026); an earlier date leaves the suspected
streams empty and is not a meaningful renewal fit.
"""
function freeze_observations(
        cutoff_date::Union{Date, AbstractString};
        path::AbstractString = joinpath(@__DIR__, "..", "data",
            "observations.toml"),
        seeding_lead::Integer = SEEDING_LEAD_DAYS)
    return load_observations(path; seeding_lead, cutoff_date)
end
