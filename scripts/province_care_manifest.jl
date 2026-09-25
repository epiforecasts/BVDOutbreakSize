#!/usr/bin/env julia
#
# Emit the `[province_isolation_history]` and `[province_bed_capacity_history]`
# blocks for data/observations.toml from the two reads of the per-province
# isolation occupancy and bed figures:
#
#   data/province_care_scanned.csv   scripts/scan_province_care.jl, deterministic
#   data/province_care_read.csv      an independent blind read of the same PDFs
#
# A cell enters the manifest when both reads carry it and agree, or when one
# read carries it and the other has no figure for that (SitRep, province).
# A cell the two reads disagree on stops the script, since it has to be
# settled against the PDF by hand and the disagreement recorded in the CSV.
# Table-era bed counts (to SitRep 080) come from the scan alone: the scan
# reads the printed `Nombre de lits` row, while the blind read's table-era bed
# figures are implied from the printed occupancy rate and the national
# `bed_capacity_history` already carries that derivation. SitReps with no
# row in data/insp_sitrep_scanned.csv have no report date and are skipped,
# as the laboratory scan skips them.
#
# Each block is sparse: one `[block.province]` sub-table per province with its
# own `dates` and `values`, since the reports print the figures for different
# provinces on different days. A province with no printed figure on a day has
# no entry. Both blocks are printed to stdout for pasting into the manifest,
# as scan_province_lab.jl does for the laboratory block.
#
# Usage:
#
#   julia --project=scripts scripts/province_care_manifest.jl [scan.csv] [read.csv]

using CSV
using DataFrames
using Dates

const ROOT = normpath(joinpath(@__DIR__, ".."))
const SCAN_CSV = length(ARGS) >= 1 ? ARGS[1] :
    joinpath(ROOT, "data", "province_care_scanned.csv")
const READ_CSV = length(ARGS) >= 2 ? ARGS[2] :
    joinpath(ROOT, "data", "province_care_read.csv")
const SITREP_CSV = joinpath(ROOT, "data", "insp_sitrep_scanned.csv")
## Last SitRep of the occupation-table era; the prose era starts at 081.
const TABLE_ERA_END = 80

const KEYS = Dict(
    "Ituri" => "ituri", "Nord-Kivu" => "nord_kivu", "Sud-Kivu" => "sud_kivu",
    "Haut-Uele" => "haut_uele", "Haut-Uélé" => "haut_uele", "Tshopo" => "tshopo",
    "Bas-Uele" => "bas_uele", "Bas-Uélé" => "bas_uele",
    "Sud Ubangi" => "sud_ubangi", "Sud-Ubangi" => "sud_ubangi",
)
const ORDER = [
    "ituri", "nord_kivu", "sud_kivu", "haut_uele", "tshopo", "bas_uele",
    "sud_ubangi",
]

## SitRep numbers that have a report date in the hand-scanned headline CSV.
function dated_sitreps()
    out = Dict{Int, Date}()
    for (i, line) in enumerate(eachline(SITREP_CSV))
        i == 1 && continue
        fields = split(line, ',')
        length(fields) >= 2 || continue
        sr = tryparse(Int, fields[1])
        sr === nothing && continue
        d = tryparse(Date, fields[2])
        d === nothing || (out[sr] = d)
    end
    return out
end

function read_cells(path, column)
    df = CSV.read(path, DataFrame; types = Dict(column => Union{Missing, Int}, :sitrep => Int))
    out = Dict{Tuple{Int, String}, Int}()
    for row in eachrow(df)
        v = row[column]
        ismissing(v) && continue
        key = get(KEYS, String(row.province), nothing)
        key === nothing && error("$(path): unknown province `$(row.province)`")
        out[(Int(row.sitrep), key)] = Int(v)
    end
    return out
end

"""
Union of the two reads for one field, with the agreement rules above.
`table_era_scan_only` drops the blind read's table-era cells (the beds).
"""
function reconcile(scan, read, dates; table_era_scan_only::Bool = false)
    out = Dict{String, Vector{Tuple{Date, Int}}}()
    conflicts = String[]
    for key in union(keys(scan), keys(read))
        sr, prov = key
        haskey(dates, sr) || continue
        a = get(scan, key, nothing)
        b = get(read, key, nothing)
        if table_era_scan_only && sr <= TABLE_ERA_END
            b = nothing
        end
        v = if a !== nothing && b !== nothing
            a == b || push!(conflicts, "SitRep $(sr) $(prov): scan $(a), read $(b)")
            a
        else
            a === nothing ? b : a
        end
        v === nothing && continue
        push!(get!(out, prov, Tuple{Date, Int}[]), (dates[sr], v))
    end
    isempty(conflicts) || error(
        "the two reads disagree; settle these against the PDF and record " *
            "the decision in the CSV:\n  " * join(conflicts, "\n  ")
    )
    for pairs in values(out)
        sort!(pairs; by = first)
        ds = first.(pairs)
        allunique(ds) || error("more than one figure on one date: $(ds)")
    end
    return out
end

function emit(io, block, comment, source, data)
    println(io, comment)
    println(io, "[", block, "]")
    println(io, "source = \"", source, "\"")
    for key in ORDER
        haskey(data, key) || continue
        pairs = data[key]
        println(io)
        println(io, "[", block, ".", key, "]")
        println(io, "dates = [", join(("\"$(d)\"" for (d, _) in pairs), ", "), "]")
        println(io, "values = [", join((v for (_, v) in pairs), ", "), "]")
    end
    return println(io)
end

function main()
    dates = dated_sitreps()
    isolation = reconcile(
        read_cells(SCAN_CSV, :patients_isolated),
        read_cells(READ_CSV, :patients_isolated), dates
    )
    beds = reconcile(
        read_cells(SCAN_CSV, :beds), read_cells(READ_CSV, :beds), dates;
        table_era_scan_only = true
    )
    println("===== paste into data/observations.toml =====\n")
    emit(
        stdout, "province_isolation_history",
        """
        # Per-province patients in isolation at the end of the day, from the
        # occupation table of the INSP situation reports (Tableau IV, V, 5, 6
        # or 7 by vintage, to SitRep 080) and the per-province care-continuity
        # prose from SitRep 081. Generated by
        # `julia --project=scripts scripts/province_care_manifest.jl` from the
        # deterministic scan (data/province_care_scanned.csv) and an
        # independent blind read (data/province_care_read.csv), which agree on
        # every cell both carry. One sub-table per province with its own
        # dates, since coverage differs by province and by day; a province
        # that prints nothing on a day has no entry. Where a report separates
        # patients in normed structures from the total hospitalised (Nord-Kivu
        # from SitRep 124), the value is the total, which is what the national
        # tile counts. Fitted as a split of the printed sum of the provinces
        # present each day, alongside the national isolation_history.""",
        "INSP situation reports, occupation table (Patients en isolement, Fin J) to SitRep 080 and Continuité des soins / Prise en charge holistique prose from SitRep 081, per province; scripts/scan_province_care.jl reconciled with an independent blind read by scripts/province_care_manifest.jl.",
        isolation
    )
    emit(
        stdout, "province_bed_capacity_history",
        """
        # Per-province isolation and treatment beds (`lits`), from the
        # `Nombre de lits` row of the occupation table (from SitRep 062) and
        # the bed counts in the care-continuity prose from SitRep 081.
        # Generated by `julia --project=scripts scripts/province_care_manifest.jl`
        # from the same two reads as the occupancy block; table-era counts are
        # the scan's, since the blind read's table-era figures are implied
        # from the printed rate. The national bed_capacity_history is the sum
        # over the provinces that print a figure on a day; this block keeps
        # the split. Nord-Kivu's count from SitRep 124 is the beds in normed
        # structures, the figure its printed rate refers to. Fitted as a split
        # of the printed sum of the provinces present each day, over
        # per-patch capacity walks.""",
        "INSP situation reports, per-province bed counts (lits) from the occupation table and the care-continuity prose; scripts/scan_province_care.jl reconciled with an independent blind read by scripts/province_care_manifest.jl.",
        beds
    )
    for (name, data) in (("isolation", isolation), ("beds", beds))
        n = sum(length(v) for v in values(data); init = 0)
        println(stderr, "$(name): $(n) entries over $(length(data)) provinces")
    end
    return nothing
end

main()
