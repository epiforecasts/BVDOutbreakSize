#!/usr/bin/env julia
#
# Cross-check the per-health-zone cumulative confirmed-case and
# confirmed-death blocks in data/observations.toml against the INRB-UMIE
# per-zone daily CSVs (a second transcription of the same INSP situation
# reports, https://github.com/INRB-UMIE/BDBV2026-Data,
# data/insp_sitrep/processed).
#
# The PDFs are the source of truth and scripts/scan_zone_tableau2.jl reads
# them directly; the mirror is a second pair of eyes, never a replacement
# (data/README.md, and the header of scripts/confirm_insp_data.jl for the
# national series). The mirror carries its own transcription artefacts: on
# some dates its zone rows do not sum to the province total the report
# prints, and it records the unallocated row under the zone name NA
# without a province.
#
# Usage:
#
#   julia --project=scripts scripts/confirm_zone_data.jl
#   julia --project=scripts scripts/confirm_zone_data.jl path/to/mirror/dir
#
# With no argument the two CSVs are downloaded; with one they are read from
# that directory (the files are insp_sitrep__cumulative_confirmed_cases__daily.csv
# and insp_sitrep__cumulative_confirmed_deaths__daily.csv, with or without
# the `__daily` suffix). Prints, per block, the number of zone-dates the
# two transcriptions share and agree on, every disagreement with both
# values, the mirror zones that resolve to no zone in the manifest, and
# the dates each source carries alone. Always exits zero: the mirror is
# not a gate on the manifest, and the disagreements are characterised in
# data/README.md.

using Downloads
using Printf
using TOML

## Shares the name folding and the zone alias table with the scanner, so a
## spelling resolves the same way in both. The scanner's main() only runs
## when it is the program file.
include(joinpath(@__DIR__, "scan_zone_tableau2.jl"))

const BASE_URL = "https://raw.githubusercontent.com/INRB-UMIE/" *
    "BDBV2026-Data/main/data/insp_sitrep/processed"
const FILES = (
    ("zone_confirmed_history", "cumulative_confirmed_cases"),
    ("zone_death_history", "cumulative_confirmed_deaths"),
)

## Mirror zone names that are not SitRep spellings of a zone in the
## manifest, keyed on their folded form. The mirror disambiguates
## Lubunga with a province suffix; the others are zones it tracks from the
## alert tables of the first reports, which never carried a confirmed case.
const MIRROR_NAMES = Dict("lubunga tshopo" => "lubunga")

function mirror_rows(stem, dir)
    path = if dir === nothing
        Downloads.download("$(BASE_URL)/insp_sitrep__$(stem)__daily.csv")
    else
        candidates = [
            joinpath(dir, "insp_sitrep__$(stem)__daily.csv"),
            joinpath(dir, "insp_sitrep__$(stem).csv"),
        ]
        i = findfirst(isfile, candidates)
        i === nothing && error("no $(stem) CSV under $(dir)")
        candidates[i]
    end
    rows = Tuple{String, String, Union{Int, Nothing}}[]
    for (i, line) in enumerate(eachline(path))
        l = replace(strip(line), "\ufeff" => "")
        isempty(l) && continue
        f = split(l, ',')
        length(f) >= 3 || continue
        i == 1 && lowercase(f[1]) == "nom" && continue
        v = f[3] == "ND" ? nothing : tryparse(Int, f[3])
        push!(rows, (String(f[1]), String(f[2]), v))
    end
    return rows
end

function main()
    dir = length(ARGS) >= 1 ? ARGS[1] : nothing
    manifest = TOML.parsefile(joinpath(ROOT, "data", "observations.toml"))
    for (block, stem) in FILES
        blk = manifest[block]
        dates = String.(blk["dates"])
        ## Zone key → (province, series), and the summed unallocated row.
        series = Dict{String, Tuple{String, Vector{Int}}}()
        unallocated = zeros(Int, length(dates))
        for (prov, zones) in blk
            zones isa AbstractDict || continue
            for (z, v) in zones
                if z == "unallocated"
                    unallocated .+= Int.(v)
                else
                    haskey(series, z) && error("zone $(z) under two provinces")
                    series[z] = (prov, Int.(v))
                end
            end
        end
        date_idx = Dict(d => i for (i, d) in enumerate(dates))

        rows = mirror_rows(stem, dir)
        agree = 0
        disagree = Tuple{String, String, Int, Int}[]
        unmatched = Dict{String, Int}()
        nd = 0
        mirror_dates = Set{String}()
        ## A mirror date that is not an ISO day would otherwise be listed
        ## as a vintage the manifest is missing, which is the line the
        ## update routine reads to decide the manifest has fallen behind.
        malformed = Dict{String, Int}()
        for (nom, date, v) in rows
            if !occursin(r"^\d{4}-\d{2}-\d{2}$", date)
                malformed[date] = get(malformed, date, 0) + 1
                continue
            end
            push!(mirror_dates, date)
            v === nothing && (nd += 1; continue)
            c = canon(nom)
            key = c == "na" ? "unallocated" :
                replace(
                    get(MIRROR_NAMES, c, get(ZONE_ALIASES, c, c)),
                    " " => "_"
                )
            haskey(date_idx, date) || continue
            i = date_idx[date]
            ours = if key == "unallocated"
                unallocated[i]
            elseif haskey(series, key)
                series[key][2][i]
            else
                unmatched[nom] = max(get(unmatched, nom, 0), v)
                continue
            end
            ours == v ? (agree += 1) : push!(disagree, (date, key, ours, v))
        end

        println("=== $(block) against the mirror's $(stem) ===")
        println(
            "$(agree + length(disagree)) zone-dates compared: " *
                "$(agree) agree, $(length(disagree)) disagree " *
                "($(nd) mirror ND cells skipped)"
        )
        if !isempty(disagree)
            println("\nDisagreements (date, zone, ours, mirror):")
            for (d, k, o, m) in sort(disagree)
                @printf("  %s  %-18s %6d %6d  (%+d)\n", d, k, o, m, m - o)
            end
            bydate = Dict{String, Int}()
            for (d, _, _, _) in disagree
                bydate[d] = get(bydate, d, 0) + 1
            end
            println(
                "\nDisagreements per date: ",
                join(("$(d): $(n)" for (d, n) in sort(collect(bydate))), ", ")
            )
        end
        if !isempty(unmatched)
            println(
                "\nMirror zones with no series in the manifest " *
                    "(name: largest value on a shared date):"
            )
            for (n, v) in sort(collect(unmatched))
                println("  $(n): $(v)")
            end
        end
        isempty(malformed) || println(
            "\nMirror rows whose date is not an ISO day, skipped: ",
            join(("$(d) x$(n)" for (d, n) in sort(collect(malformed))), ", ")
        )
        only_mirror = sort(collect(setdiff(mirror_dates, dates)))
        only_ours = sort(setdiff(dates, mirror_dates))
        isempty(only_mirror) ||
            println(
            "\nDates the mirror carries and the manifest does not: ",
            join(only_mirror, ", ")
        )
        isempty(only_ours) ||
            println(
            "\nDates the manifest carries and the mirror does not: ",
            join(only_ours, ", ")
        )
        println()
    end
    return nothing
end

main()
