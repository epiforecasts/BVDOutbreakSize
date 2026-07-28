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
settle to noise (mean increment near zero, frequent sign flips from
digitisation error rather than genuine late reporting) by about three
weeks: roughly 90% of a bar's eventual total is in by delay 14-17, 95-97%
by 20-25, and the tail beyond that is scan noise, not signal. `28` sits at
the point the true reporting signal has decayed into that noise floor,
matching the "~95% by 2.5 weeks, 98-99% by 3 weeks" completeness the
triangle supports.
"""
const ONSET_REPORT_MAX_DELAY = 28

"""
    _read_onset_curve_blocks(path)

Parse the digitised onset-curve CSV at `path` into one block per SitRep
vintage, in file order. Each row is `sitrep,report_date,onset_date,
confirmed_alive,confirmed_dead,confirmed_total`; only `confirmed_total` is
used (the alive/dead split is not modelled here, see
[`onset_reporting_model`](@ref)). A manual line-split parser, not CSV.jl:
the file has a fixed six-column schema and no quoting or embedded commas.
Returns a `Vector` of `(; sitrep::String, report_date::Date,
onsets::Dict{Date,Int})`, one entry per distinct `sitrep` id in first-seen
order.
"""
function _read_onset_curve_blocks(path::AbstractString)
    blocks = NamedTuple{(:sitrep, :report_date, :onsets),
        Tuple{String, Date, Dict{Date, Int}}}[]
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
                push!(blocks, (; sitrep, report_date,
                    onsets = Dict(onset_date => total)))
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

Blocks are canonicalised as their sorted `(onset_date, total)` pairs; blocks
sharing a canonical form collapse to the one with the earliest report date.
This is exact-value equality over the digitised content, not a hardcoded
vintage-id list, so a future reprint (of a different figure) is caught
automatically. Returns the surviving blocks sorted by ascending
`report_date`.
"""
function _dedup_onset_blocks(blocks)
    kept = Dict{Vector{Pair{Date, Int}}, Int}()
    out = @NamedTuple{sitrep::String, report_date::Date,
        onsets::Dict{Date, Int}}[]
    for b in blocks
        key = sort(collect(pairs(b.onsets)); by = first)
        if haskey(kept, key)
            j = kept[key]
            if b.report_date < out[j].report_date
                out[j] = (; sitrep = out[j].sitrep,
                    report_date = b.report_date, onsets = out[j].onsets)
            end
        else
            push!(out, (; sitrep = b.sitrep, report_date = b.report_date,
                onsets = b.onsets))
            kept[key] = length(out)
        end
    end
    sort!(out; by = x -> x.report_date)
    return out
end

"""
    load_onset_curve(path; cutoff, seeding,
        max_delay = ONSET_REPORT_MAX_DELAY, horizon = max_delay)

Load the digitised symptom-onset reporting triangle at `path` and build the
between-vintage increment cells [`onset_reporting_model`](@ref) fits.

Distinct SitRep vintages are recovered by exact-value dedup
([`_dedup_onset_blocks`](@ref)): reprinted figures collapse to their
earliest report date. Vintages reported after `cutoff` are dropped, the
same convention `load_observations`'s `history()`/`event_days()` closures
use for every other stream (see `src/data.jl`), so advancing the manifest
`as_of_date` past a newly-digitised vintage's report date picks that
vintage up automatically with no code change.

For each pair of consecutive (post-dedup, in-cutoff) vintages `s-1, s` and
each onset date `u` in the trailing `horizon`-day window
`[report_day(s) - horizon + 1, report_day(s)]` (clamped to grid day 1), one
increment cell is built:

```
y = confirmed_total(s, u) - confirmed_total(s-1, u)
```

The window is further clipped to the onset dates both vintages' figures
actually print. Each block has its own printed extent (its earliest to
latest digitised onset date), and the published figures stop their x axis
short of the report date by anything from zero to eight days: the axis
simply ends, with substantial counts on the last printed bar, rather than
running on with zero-height bars. Reading an uncovered date as zero would
therefore assert "nothing had been reported at this onset date yet" when the
figure says nothing at all about it, and because the axis gap (about five
days for most vintages) is close to the reporting delay itself, that
assertion lands squarely on the fastest part of the delay distribution: it
forces the hazard to about zero for delays 0-4 and then puts a spike at the
delay where the axis first covers the date. Cells outside either vintage's
extent are dropped instead, as unobserved rather than zero. The publisher's
choice of axis limit does not depend on the counts it hides, so dropping
them treats them as missing at random. Inside a block's own extent a missing
row is a digitisation omission of a zero-height bar and does read as a true
zero.

The cost of dropping them is that the shortest delays are observed rarely
or not at all: a correction cell needs both vintages to print the date, so
the smallest delay any consecutive pair reaches is `report_day(s) -
min(extent(s), extent(s-1))`, which for the current data is two days. Delay
zero is never observed, and the hazard there rests on partial pooling
towards the baseline (`σ_h0`, [`onset_report_hazard_model`](@ref)) rather
than on data.

The very first in-cutoff vintage is differenced against an implicit empty
predecessor via the sentinel `prev_report_days[i] = 0`: since
[`onset_report_cdf`](@ref) returns `0` for any negative delay and grid day
`0` postdates no valid onset day (`prev_report_days[i] - onset_days[i] < 0`
for every onset day `>= 1`), this recovers signal from the first digitised
vintage as a "difference from nothing" rather than discarding it, with no
extra branch anywhere downstream. Those cells score a level rather than a
correction, which is what anchors the hazard's asymptote: corrections only
ever pin differences of `F`, so without a level somewhere the asymptote
would float. [`onset_report_scales`](@ref) gives them the counting variation
a level carries and a correction largely does not.

`horizon` defaults to `max_delay` (a deliberate tie, not a second free
hyperparameter): it is exactly the hazard's own support, so scoring covers
every onset date the model expects a non-negligible increment for while
dropping settled dates that carry only digitisation noise. The settled
dates' levels are given up along with their noise, so the stream informs the
onset curve over the trailing four weeks only, not over the whole epidemic.

Returns `(; onset_days, report_days, prev_report_days, increments,
total_days, total_counts, last_total)`. The first four are length-matched
`Vector{Int}`s (1-based grid day-indices for the first three, the observed
increment for the fourth) ready for [`onset_reporting_model`](@ref).
`total_days`/`total_counts` are the per-vintage cumulative confirmed total
printed by each surviving vintage, keyed on its report day: the same
`(days, counts)` shape every other stream's history carries, so the
reported-onset total can be scored like any other observed series (see
`forecast_onsets`). It is built from every printed bar of a vintage, not
from the scored cells, which cover only the trailing `horizon` window.
`last_total` is the final entry of `total_counts` (a convenience scalar for
stream-comparison plots), or `missing` when no vintage survives.

The per-vintage totals are NOT monotone across vintages: the ≈4%
per-scan level error means a later scan can read a smaller total than an
earlier one even though late reporting only ever adds cases. Consumers
that need a non-decreasing series must say what they do with a fall
rather than assume it cannot happen.

A missing `path`, or a manifest with no in-cutoff vintage, returns the
same empty, `missing`-total shape, so the stream degrades to a no-op
rather than throwing.
"""
function load_onset_curve(path::AbstractString;
        cutoff::Date, seeding::Date,
        max_delay::Integer = ONSET_REPORT_MAX_DELAY,
        horizon::Integer = max_delay)
    noop = (; onset_days = Int[], report_days = Int[],
        prev_report_days = Int[], increments = Int[],
        total_days = Int[], total_counts = Int[], last_total = missing)
    isfile(path) || return noop

    snaps = filter(b -> b.report_date <= cutoff,
        _dedup_onset_blocks(_read_onset_curve_blocks(path)))
    isempty(snaps) && return noop

    ## Grid day-index of a calendar date: seeding day is day 1, matching
    ## `load_observations`'s `_index` for the same `(cutoff, seeding)` pair.
    _idx(d::Date) = Int(date2epochdays(d) - date2epochdays(seeding)) + 1
    _date(u::Integer) = seeding + Day(u - 1)

    ## Each block's own printed extent (earliest/latest digitised onset
    ## date), as grid day-indices. A date inside the extent that has no row
    ## is a zero-height bar; a date outside it is not covered by that
    ## figure at all and carries no observation. See the docstring above.
    extents = [(_idx(minimum(keys(snap.onsets))),
                   _idx(maximum(keys(snap.onsets)))) for snap in snaps]

    onset_days = Int[]
    report_days = Int[]
    prev_report_days = Int[]
    increments = Int[]
    H = max(Int(horizon), 1)
    for s in eachindex(snaps)
        R = _idx(snaps[s].report_date)
        Rprev = s == 1 ? 0 : _idx(snaps[s - 1].report_date)
        ## Onset dates both this vintage and its predecessor print. The
        ## first vintage has no predecessor to intersect with, so its cells
        ## span its own extent alone and score levels.
        cov_lo, cov_hi = extents[s]
        if s > 1
            cov_lo = max(cov_lo, extents[s - 1][1])
            cov_hi = min(cov_hi, extents[s - 1][2])
        end
        lo = max(R - H + 1, 1, cov_lo)
        hi = min(R, cov_hi)
        for u in lo:hi
            d = _date(u)
            prev = s == 1 ? 0 : get(snaps[s - 1].onsets, d, 0)
            cur = get(snaps[s].onsets, d, 0)
            push!(onset_days, u)
            push!(report_days, R)
            push!(prev_report_days, Rprev)
            push!(increments, cur - prev)
        end
    end
    ## Per-vintage cumulative confirmed total, over every printed bar rather
    ## than the scored cells: the scored window covers only the trailing
    ## `horizon` days, so summing the cells would give a rolling partial sum
    ## instead of the total the figure reports.
    total_days = [_idx(snap.report_date) for snap in snaps]
    total_counts = [sum(values(snap.onsets)) for snap in snaps]
    return (; onset_days, report_days, prev_report_days, increments,
        total_days, total_counts, last_total = total_counts[end])
end
