# Observation data loading from a TOML manifest. The renewal model takes
# the cut-off grid length, the per-stream cumulative totals, and the
# per-vintage histories with their day grids; this loader returns them in
# one named tuple together with the first WHO situation-report offset used
# as the intervention breakpoint.

"""
Load the BVD observation manifest from `path` (a TOML file). Returns a
named tuple with the cut-off grid length `n`, the cut-off and seeding
dates, the per-stream cumulative totals, the per-vintage histories with
their day grids, the genetic TMRCA bound, and the first WHO
situation-report offset `who_first_sitrep_days` (days from that report to
the cut-off, inclusive). The intervention breakpoint grid day is
`n - who_first_sitrep_days`.
"""
function load_observations(
        path::AbstractString = joinpath(@__DIR__, "..", "data",
        "observations.toml"))
    raw = TOML.parsefile(path)
    cutoff = Date(raw["cutoff_date"])
    seeding = Date(raw["seeding_date"])
    n = Int(date2epochdays(cutoff) - date2epochdays(seeding)) + 1
    who_first = Date(raw["who_first_sitrep_date"])
    who_first_sitrep_days = Int(date2epochdays(cutoff) - date2epochdays(who_first)) + 1
    streams = raw["streams"]
    histories = get(raw, "histories", Dict{String, Any}())
    function history(key)
        h = get(histories, key, nothing)
        h === nothing && return (; days = Int[], counts = Int[])
        return (; days = Int.(h["days"]), counts = Int.(h["counts"]))
    end
    return (; n, cutoff, seeding,
        exported_cases = Int(streams["exported_cases"]),
        total_deaths = Int(streams["total_deaths"]),
        reported_cases = Int(streams["reported_cases"]),
        confirmed_cases = Int(streams["confirmed_cases"]),
        exports_deaths = Int(streams["exports_deaths"]),
        deaths_history = history("deaths"),
        reported_history = history("reported"),
        confirmed_history = history("confirmed"),
        lab_history = history("lab"),
        tmrca_days = get(raw, "tmrca_days", missing),
        who_first_sitrep_days)
end

## Doubling-count prior base: McCabe et al.'s first report (18 May
## 2026), whose Method 2 central scenario of 501 cases implies a
## doubling count `m = log2(501) ≈ 9`. The prior centre advances from
## this base by one doubling per `M_PRIOR_DOUBLING_DAYS` (the central
## 14-day doubling time) of elapsed time to the cut-off.
const M_PRIOR_BASE_DATE = "2026-05-18"
const M_PRIOR_BASE = 9.0
const M_PRIOR_DOUBLING_DAYS = 14.0

"""
    m_prior_centre(as_of_date; base_date, m_base, doubling_days)

Centre for the doubling-count prior `m`, based on `m_base` doublings at
`base_date` and advancing by one doubling per `doubling_days` of elapsed
time to `as_of_date`:

```math
m_0 = m_\\text{base} +
    \\frac{\\text{as\\_of} - \\text{base}}{\\text{doubling\\_days}}.
```

The base is McCabe et al.'s first report (18 May 2026; Method 2 central
501 cases ⇒ `m ≈ 9`), advancing at the central 14-day doubling time, so
the prior stays centred on the plausible outbreak size as the cut-off
moves — it tracks data refreshes without manual edits, and a
McCabe-date fit recovers the base value.
"""
function m_prior_centre(as_of_date::AbstractString;
        base_date::AbstractString = M_PRIOR_BASE_DATE,
        m_base::Real = M_PRIOR_BASE,
        doubling_days::Real = M_PRIOR_DOUBLING_DAYS)
    elapsed = date2epochdays(Date(as_of_date)) -
              date2epochdays(Date(base_date))
    return m_base + elapsed / doubling_days
end
