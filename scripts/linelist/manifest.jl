# Build a BVDOutbreakSize observation manifest with the case streams replaced.
#
# Shared by scripts/linelist/fit_joint.jl and
# scripts/linelist/fit_single.jl, so the two cannot drift: a fit of one
# stream and a fit of all of them must be reading the same substitution or the
# comparison between them means nothing.
#
# `include` this, do not run it.

using CSV
using DataFrames
using TOML

## The blocks the line list can supply. Everything else in the manifest stays
## as the situation reports gave it.
const LINELIST_BLOCKS = ("confirmed_case_history", "reported_case_history",
    "suspected_daily_history")

## Build the manifest: the released one with the case streams swapped out.
## Written to a file rather than passed in memory because `load_observations`
## reads a path, and because the manifest is then a recordable artefact of what
## was fitted.
##
## `source` records how the replacement streams were indexed, because that is
## the difference between a series with a ragged edge and one without, and it is
## not recoverable from the numbers alone.
function write_manifest(; released, streams, out,
        source = "DHIS2 case line list (bvd-data-mirror); built by " *
                 "bvd-analysis's stream builder")
    raw = TOML.parsefile(released)
    df = CSV.read(streams, DataFrame)

    for block in LINELIST_BLOCKS
        rows = sort(df[df.stream .== block, :], :date)
        isempty(rows) && error("no rows for $block in $streams")
        raw[block] = Dict{String, Any}(
            "dates" => [string(d) for d in rows.date],
            "values" => Int.(rows.value),
            "source" => source
        )
    end

    ## The cut-off scalar has to move with the history it belongs to.
    ##
    ## The released manifest freezes `reported_case_history` at 26 May, with
    ## `[reported_cases] value = 1077` matching its last vintage, because INSP
    ## began revising the suspected count downward from SitRep 013 and the series
    ## stopped being trustworthy. `load_observations` reads that scalar as the
    ## observed cumulative reported total at the cut-off and the models condition
    ## on it directly.
    ##
    ## Replacing the history without replacing the scalar leaves the model told
    ## that the cut-off total is 1077 while the history it fits reaches thirteen
    ## thousand at the same date, and it can only reconcile the two by distorting
    ## ascertainment and the growth path. So the scalar is rewritten from the
    ## replacement history's own last value.
    ##
    ## `confirmed_cases` needs no equivalent: the released manifest carries no
    ## such scalar, so `load_observations` already derives it from whichever
    ## confirmed history it is given.
    let rows = sort(df[df.stream .== "reported_case_history", :], :date)
        raw["reported_cases"] = Dict{String, Any}(
            "value" => Int(rows.value[end]),
            "source" => "Last vintage of the replacement " *
                        "reported_case_history above, so the cut-off total and " *
                        "the history it comes from cannot disagree."
        )
    end

    ## `confirmed_break_dates` marks days where INSP retrospectively harmonised
    ## its own confirmed series, with the gross counts printed in the situation
    ## reports. Those are artefacts of the situation-report series, not events
    ## in the outbreak: the line list carries its own revisions inside the case
    ## records and has no such steps. Kept, they would ask the model to absorb a
    ## harmonisation the data no longer contains, which is also why
    ## `load_observations` rejects them here. `occupancy_break_dates` stays,
    ## since the isolation stream is still the situation reports'.
    delete!(raw, "confirmed_break_dates")

    ## The cut-off moves to the last line-list day. `load_observations` drops
    ## vintages after the cut-off, so the situation-report streams are read to
    ## the same date and the two manifests stay comparable.
    last_day = maximum(CSV.read(streams, DataFrame).date)
    raw["as_of_date"] = string(last_day)

    open(out, "w") do io
        TOML.print(io, raw; sorted = true)
    end
    @info "manifest written" out as_of = string(last_day)
    return out
end
