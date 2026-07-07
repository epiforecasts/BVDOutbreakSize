#!/usr/bin/env julia
#
# Regenerate the cumulative confirmed-case and confirmed-death series for
# data/observations.toml from the INRB-UMIE national cumulative CSVs (the
# clean national transcription of the INSP situation reports,
# https://github.com/INRB-UMIE/BDBV2026-Data, data/insp_sitrep/processed),
# and confirm them against our own directly scanned headline figures in
# data/insp_sitrep_scanned.csv.
#
# These two streams are upstream-primary: the upstream national_* series
# are the source of truth and this script prints the ready-to-paste TOML
# blocks. The scan is the cross-check. Every other stream in the manifest
# (suspected, laboratory cumulatives, 24h analysed volume, daily new
# suspects, isolation, beds, recoveries) is scanned directly and is not
# handled here.
#
# Usage:
#
#   julia --project=scripts scripts/confirm_insp_data.jl
#
# Prints the regenerated [confirmed_case_history] and
# [confirmed_death_history] blocks, a per-date scan-vs-upstream
# reconciliation, and the report dates each source carries that the other
# does not. Exits non-zero if any overlapping date disagrees.

using CSV
using Chain
using DataFrames
using DataFramesMeta
using Dates
using Downloads

const BASE_URL = "https://raw.githubusercontent.com/INRB-UMIE/" *
                 "BDBV2026-Data/main/data/insp_sitrep/processed"
const SCANNED = joinpath(@__DIR__, "..", "data", "insp_sitrep_scanned.csv")

# (upstream file stem, scanned column, TOML key, human label) per series.
const SERIES = [
    ("insp_sitrep__national_cumulative_confirmed_cases__daily.csv",
        :confirmed_cases, "confirmed_case_history", "confirmed cases"),
    ("insp_sitrep__national_cumulative_confirmed_deaths__daily.csv",
        :confirmed_deaths, "confirmed_death_history", "confirmed deaths")
]

# Vintages the loader excludes from the fit (`drop_dates`): a cumulative
# whose between-vintage increment (the daily confirmed deaths) is a repeated
# round value rather than a genuine daily count. The upstream CSV carries
# these totals, so the block is regenerated with the dates present and this
# list re-emits the `drop_dates` key that marks them. Keyed by TOML block.
# The confirmed-death daily increment is +39 on 27, 29 and 30 June; the
# 29 and 30 June repeats after the first are dropped.
const DROP_DATES = Dict(
    "confirmed_death_history" => ["2026-06-29", "2026-06-30"])

# Latest non-missing scanned value per report date. The scanned file can
# carry two rows for one date (e.g. SitRep 012 and its revised re-issue
# 012_v2); take the last vintage that reports the column.
function scanned_series(scan, col)
    @chain scan begin
        @rsubset !ismissing($col)
        @select :date=:report_date :scanned=$col
        groupby(:date)
        @combine :scanned = last(:scanned)
    end
end

# National cumulative series from the upstream processed CSV, ND dropped,
# sorted oldest-first.
function upstream_series(file)
    df = CSV.read(Downloads.download("$BASE_URL/$file"),
        DataFrame; missingstring = ["ND"])
    value_col = names(df)[3]
    @chain df begin
        @rename :upstream = $value_col
        @transform :date = Date.(:date)
        @rsubset !ismissing(:upstream)
        @select :date :upstream
        @orderby :date
    end
end

# Print a TOML array wrapped to the observations.toml house style, with
# the continuation lines indented to align under the first item.
function print_wrapped(label, items, per_line)
    prefix = rpad(label, 6) * " = ["
    pad = " "^length(prefix)
    print(prefix)
    for (i, item) in enumerate(items)
        i > 1 && print(i % per_line == 1 ? ",\n$pad" : ", ")
        print(item)
    end
    println("]")
end

function print_block(key, df)
    println("[$key]")
    print_wrapped("dates", ("\"$(d)\"" for d in df.date), 4)
    print_wrapped("values", df.upstream, 14)
    if haskey(DROP_DATES, key)
        print_wrapped("drop_dates",
            ("\"$(d)\"" for d in DROP_DATES[key]), 4)
    end
    println()
end

scan = CSV.read(SCANNED, DataFrame; missingstring = [""])

println("# Regenerated from the INRB-UMIE national cumulative CSVs.")
println("# Paste into data/observations.toml (keep the surrounding ",
    "comments).\n")

any_mismatch = false
for (file, col, key, label) in SERIES
    upstream = upstream_series(file)
    print_block(key, upstream)
end

for (file, col, key, label) in SERIES
    println("=== $label: scan vs upstream ===")
    scanned = scanned_series(scan, col)
    upstream = upstream_series(file)

    shared = @chain innerjoin(scanned, upstream; on = :date) @orderby :date
    mismatched = @rsubset(shared, :scanned != :upstream)
    for row in eachrow(mismatched)
        println("  MISMATCH $(row.date): scanned $(row.scanned) vs " *
                "upstream $(row.upstream)")
    end
    nrow(mismatched) > 0 && (global any_mismatch = true)
    println("  $(nrow(shared)) dates compared, $(nrow(mismatched)) " *
            "mismatch(es)")

    new_dates = @chain antijoin(upstream, scanned; on = :date) @orderby :date
    if nrow(new_dates) > 0
        println("  upstream dates not scanned: " *
                join(("$(r.date)=$(r.upstream)" for r in eachrow(new_dates)),
            ", "))
    end
    println()
end

if any_mismatch
    println("Reconciliation FAILED: a scanned value disagrees with upstream.")
    exit(1)
else
    println("Reconciliation OK: all overlapping dates agree.")
end
