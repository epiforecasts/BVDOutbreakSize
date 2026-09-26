# Symptom-onset reporting-triangle loader. Parses the digitised onset-curve
# CSV (`data/onset_curve_scanned.csv`, one block per SitRep vintage),
# collapses byte-identical reprinted blocks to their earliest vintage,
# filters to the grid cut-off (the same convention `load_observations` uses
# for every other stream, see `src/data.jl`), and builds the between-vintage
# increment cells the reporting-delay hazard model
# (`onset_reporting_model`, `src/models/observations.jl`) fits.

"""
    ONSET_REPORT_MAX_DELAY

Maximum symptom-onset-to-report delay (days) the hazard model tracks,
`d = 0 … D-1`. The digitised triangle's own between-vintage increments
settle to noise by about three weeks. Roughly 90% of a bar's eventual
total is in by delay 14-17 and 95-97% by 20-25, and the tail beyond that
is scan noise. `28` sits where the reporting signal has decayed into that
noise floor.
"""
const ONSET_REPORT_MAX_DELAY = 28

"""
    _read_onset_curve_blocks(path)

Parse the digitised onset-curve CSV at `path` into one block per SitRep
vintage, in file order. Each row is `sitrep,report_date,onset_date,
confirmed_alive,confirmed_dead,confirmed_total`. Only `confirmed_total` is
used, since the alive/dead split is not modelled here. A manual line-split
parser rather than CSV.jl, as the file has a fixed six-column schema and
no quoting or embedded commas. Returns a `Vector` of
`(; sitrep::String, report_date::Date, onsets::Dict{Date,Int})`, one entry
per distinct `sitrep` id in first-seen order.
"""
function _read_onset_curve_blocks(path::AbstractString)
    blocks = NamedTuple{
        (:sitrep, :report_date, :onsets),
        Tuple{String, Date, Dict{Date, Int}},
    }[]
    order = Dict{String, Int}()
    open(path) do io
        first_line = true
        for line in eachline(io)
            if first_line
                first_line = false
                continue
            end
            isempty(strip(line)) && continue
            fields = split(line, ',')
            length(fields) == 6 || continue
            sitrep = String(strip(fields[1]))
            report_date = Date(strip(fields[2]))
            onset_date = Date(strip(fields[3]))
            total = parse(Int, strip(fields[6]))
            if haskey(order, sitrep)
                blocks[order[sitrep]].onsets[onset_date] = total
            else
                push!(
                    blocks, (;
                        sitrep, report_date,
                        onsets = Dict(onset_date => total),
                    )
                )
                order[sitrep] = length(blocks)
            end
        end
    end
    return blocks
end

"""
    _dedup_onset_blocks(blocks)

Collapse byte-identical reprinted onset-curve blocks. Some SitRep vintages
reprint the same figure as an earlier one (the digitisation reads an
identical onset-date -> total map), which would otherwise fabricate
increments of exactly zero across the reprint and bias the fitted delay
towards implausibly fast reporting.

Blocks are canonicalised as their sorted `(onset_date, total)` pairs, and
blocks sharing a canonical form collapse to the one with the earliest
report date. This is exact-value equality over the digitised content
rather than a hardcoded vintage-id list, so a future reprint is caught
automatically. Returns the surviving blocks sorted by ascending
`report_date`.
"""
function _dedup_onset_blocks(blocks)
    kept = Dict{Vector{Pair{Date, Int}}, Int}()
    out = @NamedTuple{
        sitrep::String, report_date::Date,
        onsets::Dict{Date, Int},
    }[]
    for b in blocks
        key = sort(collect(pairs(b.onsets)); by = first)
        if haskey(kept, key)
            j = kept[key]
            if b.report_date < out[j].report_date
                out[j] = (;
                    sitrep = out[j].sitrep,
                    report_date = b.report_date, onsets = out[j].onsets,
                )
            end
        else
            push!(
                out, (;
                    sitrep = b.sitrep, report_date = b.report_date,
                    onsets = b.onsets,
                )
            )
            kept[key] = length(out)
        end
    end
    sort!(out; by = x -> x.report_date)
    return out
end

"""
    load_onset_curve(path; cutoff, seeding)

Load the digitised symptom-onset reporting triangle at `path` and build the
cells [`onset_reporting_model`](@ref) fits: each onset date's level at its
first print, then its corrections while the reporting delay still moves it.

Distinct SitRep vintages are recovered by exact-value dedup
([`_dedup_onset_blocks`](@ref)), so reprinted figures collapse to their
earliest report date. Vintages reported after `cutoff` are dropped, the
same convention `load_observations` uses for every other stream, so
advancing the manifest `as_of_date` past a newly-digitised vintage's
report date picks that vintage up with no code change.

Every onset date is scored once as a level and then only through
increments, so no printed count enters the likelihood twice:

  - Level. The first vintage that prints onset date `u` gives one cell,
    differenced against an implicit empty predecessor (the sentinel
    `prev_report_days[i] = 0`), right-truncated at that vintage's report
    day. This is the only cell for dates first printed past the delay
    support, which is the complete curve back to the start of the
    digitised window.
  - Corrections. Each later vintage `s` that prints `u` while
    `report_day(s) - u < ONSET_REPORT_MAX_DELAY` gives one cell against
    the last earlier vintage that printed `u`:

```
y = confirmed_total(s, u) - confirmed_total(prev, u)
```

A level plus its corrections telescopes to the latest print inside the
delay support. Past the support the modelled increment is zero, so later
reprints of a settled date are not scored and late reclassification is not
modelled.

Each block has its own printed extent, its earliest to latest digitised
onset date, and the published figures stop their x axis short of the
report date by anything from zero to eight days. The axis simply ends, with
substantial counts on the last printed bar, rather than running on with
zero-height bars. A date outside a vintage's extent is unobserved in that
vintage rather than zero, so it gets no cell there, and its next
correction is taken against the last vintage that did print it. The choice
of axis limit does not depend on the counts it hides, so this treats them
as missing at random. Inside a block's own extent a missing row is a
digitisation omission of a zero-height bar and does read as a true zero.

[`onset_report_cdf`](@ref) returns `0` for any negative delay and grid day
`0` postdates no valid onset day, so the sentinel recovers the level as a
difference from nothing, with no extra branch downstream. Level cells are
what anchor `alpha`. Corrections only pin differences of `F`, so without a
level somewhere `alpha` would float.

Returns `(; onset_days, report_days, prev_report_days, increments,
total_days, total_counts, last_total)`. The first four are length-matched
`Vector{Int}`s (1-based grid day-indices for the first three, the observed
cell for the fourth) ready for [`onset_reporting_model`](@ref).
`total_days` and `total_counts` are the cumulative confirmed total printed
by each surviving vintage, keyed on its report day, in the same
`(days, counts)` shape every other stream's history carries. They are
built from every printed bar of a vintage, not from the scored cells.
`last_total` is the final entry of `total_counts`, or `missing` when no
vintage survives.

The per-vintage totals are not monotone across vintages. The roughly 4%
per-scan level error means a later scan can read a smaller total than an
earlier one even though late reporting only ever adds cases, and this
happens more than once in the current data. Consumers that need a
non-decreasing series must say what they do with a fall rather than assume
it cannot happen, and anything scored against this series should be scored
on its increments rather than its level.

A missing `path`, or a manifest with no in-cutoff vintage, returns the
same empty, `missing`-total shape, so the stream degrades to a no-op
rather than throwing.
"""
function load_onset_curve(
        path::AbstractString; cutoff::Date, seeding::Date
    )
    noop = (;
        onset_days = Int[], report_days = Int[],
        prev_report_days = Int[], increments = Int[],
        total_days = Int[], total_counts = Int[], last_total = missing,
    )
    isfile(path) || return noop

    snaps = filter(
        b -> b.report_date <= cutoff,
        _dedup_onset_blocks(_read_onset_curve_blocks(path))
    )
    isempty(snaps) && return noop

    ## Grid day-index of a calendar date: seeding day is day 1, matching
    ## `load_observations`'s `_index` for the same `(cutoff, seeding)` pair.
    _idx(d::Date) = Int(date2epochdays(d) - date2epochdays(seeding)) + 1
    _date(u::Integer) = seeding + Day(u - 1)

    onset_days = Int[]
    report_days = Int[]
    prev_report_days = Int[]
    increments = Int[]
    ## The last vintage to print each onset date so far, by index.
    last_print = Dict{Int, Int}()
    for s in eachindex(snaps)
        R = _idx(snaps[s].report_date)
        ## This vintage's printed extent. A date inside it with no row is a
        ## zero-height bar; a date outside it is not covered at all.
        lo = max(1, _idx(minimum(keys(snaps[s].onsets))))
        hi = min(R, _idx(maximum(keys(snaps[s].onsets))))
        for u in lo:hi
            d = _date(u)
            cur = get(snaps[s].onsets, d, 0)
            p = get(last_print, u, 0)
            if p == 0
                push!(onset_days, u)
                push!(report_days, R)
                push!(prev_report_days, 0)
                push!(increments, cur)
            elseif R - u < ONSET_REPORT_MAX_DELAY
                push!(onset_days, u)
                push!(report_days, R)
                push!(prev_report_days, _idx(snaps[p].report_date))
                push!(increments, cur - get(snaps[p].onsets, d, 0))
            end
            last_print[u] = s
        end
    end
    total_days = [_idx(snap.report_date) for snap in snaps]
    total_counts = [sum(values(snap.onsets)) for snap in snaps]
    return (;
        onset_days, report_days, prev_report_days, increments,
        total_days, total_counts, last_total = total_counts[end],
    )
end

"""
    onset_hazard_grid_start(onset_days, report_days; D = ONSET_REPORT_MAX_DELAY)

Report-date grid day the reporting-delay calendar walk `γ`
([`onset_report_hazard_model`](@ref)) starts from: bounded below by the
earliest scored onset date, otherwise pulled forward to one delay
support's width before the earliest report day,
`max(minimum(onset_days), minimum(report_days) - D + 1, 1)`. Returns `1`
for an empty `onset_days`.

[`onset_reporting_model`](@ref) uses this to build its own `γ`, and
returns it as `grid_start` beside `alpha_grid_start`. Every caller that
evaluates the fitted hazard outside the model (the report pages, through
[`fitted_onset_hazard`](@ref)) indexes `γ` from it.
`alpha`'s own grid ([`onset_ascertainment_model`](@ref)) always starts at
`minimum(onset_days)` and is unaffected. Pure, top-level.
"""
function onset_hazard_grid_start(
        onset_days::AbstractVector{<:Integer},
        report_days::AbstractVector{<:Integer};
        D::Integer = ONSET_REPORT_MAX_DELAY
    )
    isempty(onset_days) && return 1
    u_lo = minimum(onset_days)
    return max(u_lo, minimum(report_days) - Int(D) + 1, 1)
end
