#!/usr/bin/env julia
#
# Confirm the hand-scanned DRC confirmed-case and confirmed-death totals
# in `data/insp_sitrep_scanned.csv` against the INRB-UMIE national
# cumulative series, the independent transcription of the same INSP
# situation reports (https://github.com/INRB-UMIE/BDBV2026-Data,
# data/insp_sitrep/processed). This is a cross-check, not a generator:
# we scan the situation-report PDFs directly (the headline national
# totals) and use the upstream `national_*` CSVs only to confirm those
# numbers and to flag report dates we have not yet scanned.
#
# Usage:
#
#   julia --project=scripts scripts/confirm_insp_data.jl
#
# Prints a per-date reconciliation for cumulative confirmed cases and
# confirmed deaths over the dates present in both sources, lists any
# upstream report dates missing from the scan (new data to add), and
# exits non-zero if any overlapping date disagrees.
#
# Only the cumulative confirmed-case and confirmed-death series are
# cross-checked: upstream publishes those as clean national daily CSVs
# that match our scanned headlines. The suspected, laboratory cumulative
# and 24h-analysed streams have no comparable clean national column
# upstream (suspected is reclassified downward; the 24h analysed volume
# is a per-province sum read from the report's laboratory block), so they
# stay scan-only and are not reconciled here.

using CSV
using Chain
using DataFrames
using DataFramesMeta
using Dates
using Downloads

const BASE_URL = "https://raw.githubusercontent.com/INRB-UMIE/" *
                 "BDBV2026-Data/main/data/insp_sitrep/processed"
const SCANNED = joinpath(@__DIR__, "..", "data", "insp_sitrep_scanned.csv")

# (upstream file stem, scanned column, human label) for each series we
# reconcile. The upstream value column is the third column of the CSV.
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

# National cumulative series from the upstream processed CSV, ND dropped.
function upstream_series(file)
    df = CSV.read(Downloads.download("$BASE_URL/$file"),
        DataFrame; missingstring = ["ND"])
    value_col = names(df)[3]
    @chain df begin
        @rename :upstream = $value_col
        @transform :date = Date.(:date)
        @rsubset !ismissing(:upstream)
        @select :date :upstream
    end
end

scan = CSV.read(SCANNED, DataFrame; missingstring = [""])

any_mismatch = false
for (file, col, label) in SERIES
    println("=== $label (scanned vs INRB-UMIE national) ===")
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
        println("  upstream dates not yet scanned:")
        for row in eachrow(new_dates)
            println("    $(row.date) -> $(row.upstream)")
        end
    end
    println()
end

if any_mismatch
    println("Reconciliation FAILED: scanned values disagree with upstream.")
    exit(1)
else
    println("Reconciliation OK: all overlapping dates agree.")
end
