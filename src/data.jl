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
per-stream cumulative totals at the cut-off, the optional
`tests_analysed` scalar, the per-vintage histories as `(; days, counts)`
with `days` the grid day-indices, the genetic TMRCA bound `tmrca_days`
(days before the cut-off), and `who_first_sitrep_days` (days from the
first situation report, the earliest reported-case vintage, to the
cut-off). The intervention breakpoint grid day is `n - who_first_sitrep_days`.
"""
function load_observations(
        path::AbstractString = joinpath(@__DIR__, "..", "data",
            "observations.toml");
        seeding_lead::Integer = SEEDING_LEAD_DAYS)
    raw = TOML.parsefile(path)
    _val(k) = raw[k]["value"]
    cutoff = Date(String(raw["as_of_date"]))
    ## Grid day-index (1-based) of a calendar date: seeding day is day 1.
    _gap(d) = Int(date2epochdays(cutoff) - date2epochdays(Date(String(d))))
    tmrca_date = Date(String(raw["genetic_tmrca"]["date"]))
    seeding = tmrca_date - Day(seeding_lead)
    n = Int(date2epochdays(cutoff) - date2epochdays(seeding)) + 1
    _index(d) = n - _gap(d)

    ## A dated cumulative history → grid day-indices and counts, sorted
    ## oldest-first so the model differences consecutive vintages into
    ## daily increments. Empty when the block is absent.
    function history(key)
        haskey(raw, key) || return (; days = Int[], counts = Int[])
        block = raw[key]
        idx = Int[_index(d) for d in block["dates"]]
        vals = Int.(block["values"])
        ord = sortperm(idx)
        return (; days = idx[ord], counts = vals[ord])
    end

    reported_history = history("reported_case_history")
    ## The first WHO joint situation report is the earliest reported-case
    ## vintage; days from it to the cut-off set the intervention breakpoint.
    who_first_sitrep_days = isempty(reported_history.days) ? n :
                            n - reported_history.days[1] + 1

    return (; n, cutoff, seeding,
        exported_cases = Int(_val("exported_cases")),
        exports_deaths = Int(_val("exports_deaths")),
        total_deaths = Int(_val("total_deaths")),
        reported_cases = Int(_val("reported_cases")),
        confirmed_cases = haskey(raw, "confirmed_cases") ?
                          Int(_val("confirmed_cases")) : missing,
        tests_analysed = haskey(raw, "cumulative_tests_analysed") ?
                         Int(_val("cumulative_tests_analysed")) : missing,
        reported_history = reported_history,
        confirmed_history = history("confirmed_case_history"),
        deaths_history = history("death_history"),
        lab_history = history("lab_history"),
        tmrca_days = _gap(raw["genetic_tmrca"]["date"]),
        who_first_sitrep_days)
end
