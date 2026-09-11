#!/usr/bin/env julia
#
# Cross-check the cumulative confirmed-case and confirmed-death series in
# data/observations.toml against the INRB-UMIE national cumulative CSVs
# (a national transcription of the INSP situation reports,
# https://github.com/INRB-UMIE/BDBV2026-Data, data/insp_sitrep/processed).
#
# insp.cd is the source of truth for both streams, as it is for every other
# stream in the manifest (data/README.md section 2, and the header of
# scripts/download_sitreps.jl). The fitted values are our own directly
# scanned headline and province-table figures in
# data/insp_sitrep_scanned.csv; the mirror is a second pair of eyes on the
# same PDFs, not an authority over them. It is also incomplete, carrying 86
# of the 102 report dates the manifest holds, so its series is never pasted
# into the manifest wholesale.
#
# Usage:
#
#   julia --project=scripts scripts/confirm_insp_data.jl
#
# Prints a per-date scan-vs-mirror reconciliation and the report dates each
# source carries that the other does not. Exits non-zero if an overlapping
# date disagrees, unless the disagreement is one of the documented mirror
# transcription errors in KNOWN_MIRROR_ERRORS below.

using CSV
using Chain
using DataFrames
using DataFramesMeta
using Dates
using Downloads

const BASE_URL = "https://raw.githubusercontent.com/INRB-UMIE/" *
                 "BDBV2026-Data/main/data/insp_sitrep/processed"
const SCANNED = joinpath(@__DIR__, "..", "data", "insp_sitrep_scanned.csv")

# Mirror transcription errors the primary PDFs have already settled, keyed
# by (report date, series label). Each entry records both values so a
# silently changed mirror is caught rather than waved through, and states
# what the SitRep's own pages say. A listed date is reported and does not
# fail the run; an unlisted disagreement still does, and so does an entry
# that no longer reproduces, so the table cannot quietly go stale. The
# reasoning behind each is written up in data/README.md; see issue #624.
const KNOWN_MIRROR_ERRORS = Dict(
    (Date("2026-08-08"), "confirmed cases") => (4294, 4209,
        "the mirror filed SitRep 085 (reporting date 07 August) under its " *
        "08 August publication date; SitRep 086's own page-1 headline and " *
        "province table both give 4294 for 08 August"),
    (Date("2026-08-08"), "confirmed deaths") => (1960, 1916,
        "same misdated SitRep 085 row; SitRep 086's headline and province " *
        "sum both give 1960"),
    (Date("2026-08-11"), "confirmed cases") => (4567, 4566,
        "SitRep 089 prints a Total of 4566 but its own province rows sum " *
        "to 4567 (3912+528+115+9+3), and SitRep 088's 4449 plus the " *
        "printed 118 new confirmed also gives 4567; the auditable table " *
        "sum wins, as for SitReps 009, 061 and 083"),
    (Date("2026-08-25"), "confirmed deaths") => (2744, 2755,
        "SitRep 103 gives 2744 in its headline, its Total row and its " *
        "province sum (2133+512+89+8+1+1); 2755 appears nowhere in " *
        "SitReps 101-104")
)

# (upstream file stem, scanned column, human label) per series. The label
# is what KNOWN_MIRROR_ERRORS is keyed on, so it is part of the contract.
const SERIES = [
    ("insp_sitrep__national_cumulative_confirmed_cases__daily.csv",
        :confirmed_cases, "confirmed cases"),
    ("insp_sitrep__national_cumulative_confirmed_deaths__daily.csv",
        :confirmed_deaths, "confirmed deaths")
]

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

scan = CSV.read(SCANNED, DataFrame; missingstring = [""])

any_mismatch = false
# Entries confirmed against this run's mirror, so a table that has outlived
# the disagreement it documents can be reported at the end.
seen_known = Set{Tuple{Date, String}}()

for (file, col, label) in SERIES
    println("=== $label: scan vs mirror ===")
    scanned = scanned_series(scan, col)
    upstream = upstream_series(file)

    shared = @chain innerjoin(scanned, upstream; on = :date) @orderby :date
    mismatched = @rsubset(shared, :scanned != :upstream)
    unexplained = 0
    for row in eachrow(mismatched)
        known = get(KNOWN_MIRROR_ERRORS, (row.date, label), nothing)
        if known !== nothing && known[1] == row.scanned &&
           known[2] == row.upstream
            push!(seen_known, (row.date, label))
            println("  known mirror error $(row.date): scan " *
                    "$(row.scanned) vs mirror $(row.upstream) - $(known[3])")
        else
            unexplained += 1
            println("  MISMATCH $(row.date): scanned $(row.scanned) vs " *
                    "mirror $(row.upstream)")
        end
    end
    unexplained > 0 && (global any_mismatch = true)
    println("  $(nrow(shared)) dates compared, $(nrow(mismatched)) " *
            "mismatch(es), $unexplained unexplained")

    new_dates = @chain antijoin(upstream, scanned; on = :date) @orderby :date
    if nrow(new_dates) > 0
        println("  mirror dates not scanned: " *
                join(("$(r.date)=$(r.upstream)" for r in eachrow(new_dates)),
            ", "))
    end
    only_scanned = nrow(antijoin(scanned, upstream; on = :date))
    if only_scanned > 0
        println("  $only_scanned scanned date(s) the mirror does " *
                "not carry, so its series is not a substitute for the scan")
    end
    println()
end

# An entry whose disagreement has gone (the mirror was corrected upstream,
# or the scan changed) is stale and must be removed rather than left to
# suppress a future real mismatch on the same date.
stale = setdiff(keys(KNOWN_MIRROR_ERRORS), seen_known)
if !isempty(stale)
    println("=== stale KNOWN_MIRROR_ERRORS entries ===")
    for (date, label) in sort(collect(stale))
        println("  $date $label no longer disagrees as recorded; remove it")
    end
    println()
    any_mismatch = true
end

if any_mismatch
    println("Reconciliation FAILED: an undocumented disagreement, or a " *
            "stale entry.")
    exit(1)
else
    println("Reconciliation OK: every disagreement is a documented " *
            "mirror error.")
end
