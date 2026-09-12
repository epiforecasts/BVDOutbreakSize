#!/usr/bin/env julia
#
# Scan Tableau 1 of the INSP situation reports, the per-province split of
# confirmed cases and confirmed deaths, and emit the
# `[province_confirmed_history]` and `[province_death_history]` blocks for
# data/observations.toml.
#
#   Répartition des cas et décès confirmés par province touchée
#
#    Province     Nouveaux cas(24h)   Cas confirmés   Décès   Létalité   Zones
#    Ituri                       47            5 508   2 501     45,4%   28/36
#    Nord-Kivu                   46            1 139     725     63,7%   16/34
#    Haut-Uélé                    6              264     110     41,7%    6/13
#    Tshopo                       0               24       9     37,5%    7/23
#    Sud-Kivu                     0                3       1     33,3%    1/34
#    Bas Uélé                     0                4       3     75,0%    3/11
#    Total                       99            6 942   3 349     48,2%  61/151
#
# The table grew from three provinces in four columns to six provinces in
# six, and the column order has changed at least once (the new-cases column
# moved from last to first). So the parse is layout-independent: on a
# province row it finds the first field after the province name holding a
# `%`, which is the Létalité, and reads the cumulative cases and deaths from
# the two fields immediately before it. The zones-touchées field also holds
# a `%`, hence the first one. A province missing from a vintage's table has
# not yet reported a confirmed case and is recorded as zero, not missing.
#
# Why the deaths column matters
#
# The per-province case split alone cannot separate "more infections" from
# "better case-finding": the confirmed count in a province is the product of
# its ascertainment and its incidence, and only that product is observed.
#
# The deaths column breaks the tie. Deaths are far harder to miss than cases,
# and the virus's case-fatality does not change at a provincial border. So
# under a shared true CFR, the per-province death split identifies the
# per-province incidence split, and the gap between the death split and the
# case split is what identifies the relative case ascertainment.
#
# The signal is large and sustained: Nord-Kivu holds a steady ~8.5-9% of
# confirmed cases but ~15-19% of confirmed deaths across every vintage, with a
# confirmed CFR of 54-59% against Ituri's 20-33%.
#
# Part of that gap is not ascertainment: in a fast-growing epidemic the
# observed CFR is biased down, because recent cases have not yet died. Ituri
# grows faster than Nord-Kivu, so some of its lower CFR is right-censoring.
# The model separates the two by applying the shared CFR and onset-to-death
# delay to each province's own incidence curve; this script only supplies the
# data.
#
# Validation. The per-province cases must sum to the national
# `confirmed_case_history` and the per-province deaths to the national
# `confirmed_death_history`, on every date. That reconciliation is what
# admits a vintage: a table that does not partition the national totals has
# been mis-parsed, and the script exits non-zero so it cannot reach the
# manifest. The one exception is a table that contradicts itself, its rows
# missing the national totals that the table's own printed Total row hits.
# There the report's arithmetic is at fault rather than the parse, and the
# vintage is reported and left out; SitRep 073 of 26 July is the only one so
# far, its Haut-Uélé row one case short of its own Total.
#
# Requires `pdftotext` (poppler-utils) and the sitrep PDFs
# (`task download-sitreps`).
#
# Usage:
#
#   julia --project=scripts scripts/scan_province_tableau1.jl
#   julia --project=scripts scripts/scan_province_tableau1.jl path/to/pdfdir

using TOML
using Printf

const ROOT = normpath(joinpath(@__DIR__, ".."))
const PDF_DIR = length(ARGS) >= 1 ? ARGS[1] :
                joinpath(ROOT, "data", "sitrep_pdfs")
const MANIFEST = joinpath(ROOT, "data", "observations.toml")
const SITREP_CSV = joinpath(ROOT, "data", "insp_sitrep_scanned.csv")

## Provinces in patch order: the first is the primary (origin) patch.
const PROVINCES = ["ituri", "nord_kivu", "sud_kivu", "haut_uele", "tshopo",
    "bas_uele"]
## Keyed on the folded name with every separator collapsed to one space, so
## "Bas Uélé", "Bas-Uélé" and "BAS UELE" all land on the same province.
const NAMES = Dict("ituri" => "ituri", "nord kivu" => "nord_kivu",
    "sud kivu" => "sud_kivu", "haut uele" => "haut_uele",
    "tshopo" => "tshopo", "bas uele" => "bas_uele")

fold(s) = lowercase(Base.Unicode.normalize(String(s); stripmark = true))

"""
Fold a table cell to the form the province names are keyed on: accents
stripped, every run of non-letters collapsed to one space. This absorbs the
hyphen/space and footnote-marker variation in the province column.
"""
function canon(s::AbstractString)
    t = replace(fold(s), r"[^a-z]+" => " ")
    return strip(t)
end

"""
Digits only: strips the thousands spaces ("1 792") and the footnote markers
("535*") the table carries.
"""
function digits_only(s::AbstractString)
    d = replace(s, r"[^0-9]" => "")
    return isempty(d) ? nothing : parse(Int, d)
end

"""
The lines of Tableau 1, from its caption to its Total row.

The caption has been printed three ways ("Tableau 1. Répartition ...",
"TABLEAU 1 — REPARTITION ...", and the bare heading), so it is matched on
the wording rather than the numbering. The table always closes with a Total
row; stopping there keeps the commentary that follows out of the parse.
"""
function tableau1_lines(text::AbstractString)
    lines = split(text, '\n')
    head = findfirst(
        l -> occursin(r"repartition des cas et deces confirmes par province",
            fold(l)), lines)
    head === nothing && return nothing
    out = String[]
    for l in lines[(head + 1):min(head + 40, length(lines))]
        ## Tableau 2 breaks the same province names down by health zone, so
        ## it must not be read as a continuation of Tableau 1.
        occursin("par province et zone de sante", fold(l)) && break
        push!(out, String(l))
        canon(first(split(strip(l), r"\s{2,}"))) == "total" && break
    end
    return out
end

"""
Cumulative `(cases, deaths)` from one Tableau 1 row, or `nothing`.

`start` is the first field that can hold a number. The Létalité percentage
is found rather than indexed, because the column order has changed between
vintages; the cases and deaths are the two fields before it.
"""
function parse_row(parts::Vector{<:AbstractString}, start::Int)
    pct = findfirst(i -> occursin('%', parts[i]), start:length(parts))
    pct === nothing && return nothing
    j = start + pct - 1
    j - 2 < start && return nothing
    c = digits_only(parts[j - 2])
    d = digits_only(parts[j - 1])
    (c === nothing || d === nothing) && return nothing
    return (c, d)
end

"""
Per-province `(cases, deaths)` for one sitrep, keyed on province, together
with the table's own printed Total row.

A footnote marker can displace the province name onto the line below its own
row (the harmonisation asterisk does this in the August vintages), so a row
whose first field carries no letters takes its name from the next line.
"""
function scan_tableau1(text::AbstractString)
    lines = tableau1_lines(text)
    lines === nothing && return nothing
    got = Dict{String, Tuple{Int, Int}}()
    total = nothing
    for (i, line) in enumerate(lines)
        ## Columns are separated by runs of 2+ spaces; a thousands
        ## separator is a single space, so it stays inside its column.
        parts = split(strip(line), r"\s{2,}")
        length(parts) < 3 && continue
        if canon(parts[1]) == "total"
            total === nothing && (total = parse_row(parts, 2))
            continue
        end
        name = get(NAMES, canon(parts[1]), nothing)
        start = 2
        if name === nothing && !occursin(r"[a-z]", fold(parts[1]))
            ## The row itself is all numbers: its province name sits on the
            ## next non-empty line, and every field is data.
            nxt = findfirst(l -> !isempty(strip(l)), lines[(i + 1):end])
            nxt === nothing && continue
            name = get(NAMES, canon(lines[i + nxt]), nothing)
            start = 1
        end
        name === nothing && continue
        haskey(got, name) && continue
        r = parse_row(parts, start)
        r === nothing || (got[name] = r)
    end
    return (got, total)
end

function sitrep_dates()
    dates = Dict{Int, String}()
    for (i, line) in enumerate(eachline(SITREP_CSV))
        i == 1 && continue
        f = split(line, ',')
        length(f) < 2 && continue
        n = tryparse(Int, f[1])
        n === nothing && continue
        dates[n] = f[2]
    end
    return dates
end

function main()
    isdir(PDF_DIR) || error("no sitrep PDFs at $(PDF_DIR); " *
          "run `task download-sitreps` first.")
    Sys.which("pdftotext") === nothing &&
        error("pdftotext not found; install poppler-utils.")

    dates = sitrep_dates()
    scanned = Dict{String, Dict{String, Tuple{Int, Int}}}()
    totals = Dict{String, Union{Nothing, Tuple{Int, Int}}}()
    srs = Dict{String, Int}()
    notable = Int[]

    for path in sort(filter(f -> endswith(f, ".pdf"),
        readdir(PDF_DIR; join = true)))
        m = match(r"(\d+)[_-]2026", basename(path))
        m === nothing && continue
        sr = parse(Int, m[1])
        haskey(dates, sr) || continue

        text = read(`pdftotext -layout $path -`, String)
        r = scan_tableau1(text)
        if r === nothing || isempty(r[1])
            push!(notable, sr)
            continue
        end
        scanned[dates[sr]] = r[1]
        totals[dates[sr]] = r[2]
        srs[dates[sr]] = sr
    end

    isempty(scanned) && error("no sitrep yielded a Tableau 1.")

    raw = TOML.parsefile(MANIFEST)
    natc = Dict(String(d) => v
    for (d, v) in zip(
        raw["confirmed_case_history"]["dates"],
        raw["confirmed_case_history"]["values"]))
    natd = Dict(String(d) => v
    for (d, v) in zip(
        raw["confirmed_death_history"]["dates"],
        raw["confirmed_death_history"]["values"]))

    cases(g) = sum(haskey(g, p) ? g[p][1] : 0 for p in PROVINCES)
    deaths(g) = sum(haskey(g, p) ? g[p][2] : 0 for p in PROVINCES)

    ## A date with no national totals cannot be reconciled, so it cannot be
    ## admitted; it is reported rather than dropped silently.
    all_dates = sort(collect(keys(scanned)))
    checkable = [d for d in all_dates if haskey(natc, d) && haskey(natd, d)]
    unchecked = setdiff(all_dates, checkable)

    ## Reconciliation against the national totals is what admits a vintage.
    ## The one exception is a table that contradicts itself: where the rows
    ## miss the national totals but the table's own printed Total hits them,
    ## the report's arithmetic is at fault, not the parse, and the split
    ## cannot be trusted. Those are reported and left out; everything else
    ## that misses is a mis-parse and stops the script.
    keep = String[]
    internal = String[]
    bad = String[]
    for d in checkable
        nat = (natc[d], natd[d])
        rows = (cases(scanned[d]), deaths(scanned[d]))
        if rows == nat
            push!(keep, d)
        elseif totals[d] == nat
            push!(internal, d)
        else
            push!(bad, d)
        end
    end

    ## A printed Total that misses while the rows reconcile is a typo in the
    ## Total cell alone; the split is still sound, so the date is kept.
    off_total = [d
                 for d in keep
                 if totals[d] !== nothing &&
        totals[d] != (cases(scanned[d]), deaths(scanned[d]))]

    println("Per-province confirmed cases and deaths (Tableau 1)\n")
    @printf("%11s %4s | %6s %6s | %6s %6s | %3s | %6s %6s\n", "date", "sr",
        "cases", "nat", "deaths", "nat", "np", "NK c%", "NK d%")
    println("-"^72)
    for d in sort(vcat(keep, bad))
        g = scanned[d]
        cs, ds = cases(g), deaths(g)
        nc, nd = natc[d], natd[d]
        okc = cs == nc ? string(nc) : "!$(nc)"
        okd = ds == nd ? string(nd) : "!$(nd)"
        nk = get(g, "nord_kivu", (0, 0))
        @printf("%11s %4d | %6d %6s | %6d %6s | %3d | %5.1f%% %5.1f%%\n",
            d, srs[d], cs, okc, ds, okd, length(g),
            100 * nk[1] / cs, 100 * nk[2] / ds)
    end

    if !isempty(notable)
        println("\nSitreps with no readable Tableau 1: ",
            join(notable, ", "))
    end
    if !isempty(internal)
        println("\nDates whose Tableau 1 contradicts itself, the rows " *
                "missing the national totals that its own printed Total " *
                "matches (not emitted):")
        for d in internal
            g, t = scanned[d], totals[d]
            @printf("  %11s sitrep %3d rows %d/%d, printed Total %d/%d\n",
                d, srs[d], cases(g), deaths(g), t[1], t[2])
        end
    end
    if !isempty(off_total)
        println("\nDates where the rows reconcile but the table's own " *
                "printed Total does not (kept):")
        for d in off_total
            g, t = scanned[d], totals[d]
            @printf("  %11s sitrep %3d rows %d/%d, printed Total %d/%d\n",
                d, srs[d], cases(g), deaths(g), t[1], t[2])
        end
    end
    if !isempty(unchecked)
        println(
            "\nDates with no national totals to reconcile against " *
            "(not emitted): ",
            join(unchecked, ", "))
    end

    if !isempty(bad)
        println()
        error("$(length(bad)) province/national disagreement(s) " *
              "($(join(bad, ", "))). The per-province " *
              "figures are an exact partition of the national totals, so a " *
              "mismatch means a mis-parse. Not emitting the blocks.")
    end
    println("\nAll $(length(keep)) dates reconcile with the national " *
            "confirmed case and death totals.")

    ## Nord-Kivu's death share sits well above its case share at every
    ## vintage; that gap is the identifying signal, so report it.
    ncs = [100 * get(scanned[d], "nord_kivu", (0, 0))[1] / cases(scanned[d])
           for d in keep]
    nds = [100 * get(scanned[d], "nord_kivu", (0, 0))[2] / deaths(scanned[d])
           for d in keep]
    @printf("\nNord-Kivu: case share %.1f-%.1f%%, death share %.1f-%.1f%%\n",
        minimum(ncs), maximum(ncs), minimum(nds), maximum(nds))
    @printf("Nord-Kivu death-to-case share ratio %.2fx at the cut-off\n",
        nds[end] / ncs[end])

    fmt(v) = join(v, ", ")
    println("\n\n===== paste into data/observations.toml =====\n")
    for (blk, idx, what) in (("province_confirmed_history", 1, "cases"),
        ("province_death_history", 2, "deaths"))
        println("[$(blk)]")
        println("dates = [", join(["\"$d\"" for d in keep], ", "), "]")
        for p in PROVINCES
            v = [haskey(scanned[d], p) ? scanned[d][p][idx] : 0 for d in keep]
            println("$(p) = [", fmt(v), "]")
        end
        println("source = \"INSP situation reports, Tableau 1 " *
                "(Répartition des cas et décès confirmés par province): " *
                "per-province cumulative confirmed $(what). A province " *
                "absent from a vintage's table has no confirmed $(what) " *
                "yet and is recorded as 0. Scanned by " *
                "scripts/scan_province_tableau1.jl, which requires the " *
                "per-province figures to sum exactly to the national " *
                "totals on every date.\"")
        println()
    end
    return nothing
end

main()
