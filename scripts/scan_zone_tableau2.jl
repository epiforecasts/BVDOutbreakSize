#!/usr/bin/env julia
#
# Scan Tableau 2 of the INSP situation reports, the per-health-zone split of
# confirmed cases and confirmed deaths within each province, and emit the
# `[zone_confirmed_history]` and `[zone_death_history]` blocks for
# data/observations.toml.
#
#   Répartition des cas et décès confirmés par province et zone de santé
#
#    Province / Zone de santé   Cas (n)   Décès (n)   Létalité   Nouv. cas ...
#    Ituri                        5 659       2 583      45,6%          44
#      Bunia                      1 583         505      31,9%          19
#      ...
#      A ventiler *                  NA         442         NA
#    Nord-Kivu                    1 281         797      62,2%          12
#      ...
#    Total                        7 258       3 510      48,4%          58
#
# The table has been printed in four layouts since it first appeared at
# SitRep 018 (1 June 2026):
#
# - 018-032: caption "... par zone de santé et par province", a bare
#   province heading, zone rows, an "Autres ZS (données non ventilées)" row
#   and a "Sous-total <province>" row carrying the province totals. The
#   Létalité column has no % sign before SitRep 021.
# - 034-058: caption "... par province et zone de santé", the province row
#   itself carries the totals, "Autres zones non encore identifiées" is the
#   unallocated row. SitRep 033 gives the split as prose only.
# - 059-083: adds the "Situation du jour (24 h)" columns to the right of the
#   Létalité; the unallocated row is "Non identifiées". Absent from 084-086
#   (the brief format). SitRep 087 prints the 24h flow alone with no
#   cumulative columns, so it has nothing this script can read.
# - 088-122: as 059-083, with an "A ventiler *" row of confirmed deaths in
#   the province's treatment centres not yet attributed to a zone (its
#   cases cell prints NA). Tshopo has a zone named Tshopo from SitRep 091.
#
# The parse is layout-independent in the same way as
# scan_province_tableau1.jl: on every row it finds the first field holding
# a `%` (the Létalité) and reads the cumulative cases and deaths from the
# two fields before it, falling back to the first two of three plain
# numeric fields for the vintages that print the Létalité without a sign.
# Which province a zone row belongs to is the most recent province row
# above it. A zone absent from a vintage's table has no confirmed case yet
# and is recorded as zero, as for the provinces.
#
# Zone names vary in spelling between vintages (Gety/Gethy, Rimba/Rumba,
# Tchomia/Tchomai, Mongbwalu/Mongbalu, ...). `ZONE_ALIASES` maps every
# variant seen to one key, and the script prints every zone key with the
# spellings that fed it so a new variant is visible rather than becoming a
# second series.
#
# Validation. Within each province the zone rows plus the unallocated row
# must sum exactly to that province's cumulative in the already-committed
# `[province_confirmed_history]` and `[province_death_history]` blocks, on
# every date, and those in turn to the national totals. A vintage whose
# zone rows do not partition the printed province row has a table that
# contradicts itself (SitRep 059 says so in its own footnote: the Ituri
# zone split is a stale consolidated distribution while the province row
# is current) and is reported and left out. A vintage whose zone rows
# match the printed province row but not the committed value is a
# mis-parse or a manifest discrepancy, and the script stops rather than
# emit it. Dates the province blocks do not carry (SitReps 018-031, before
# Tableau 1 existed) are reconciled against the printed province rows and
# the national totals instead, and reported as such.
#
# Requires `pdftotext` (poppler-utils) and the sitrep PDFs
# (`task download-sitreps`).
#
# Usage:
#
#   julia --project=scripts scripts/scan_zone_tableau2.jl
#   julia --project=scripts scripts/scan_zone_tableau2.jl path/to/pdfdir

using TOML
using Printf

const ROOT = normpath(joinpath(@__DIR__, ".."))
const PDF_DIR = length(ARGS) >= 1 ? ARGS[1] :
    joinpath(ROOT, "data", "sitrep_pdfs")
const MANIFEST = joinpath(ROOT, "data", "observations.toml")
const SITREP_CSV = joinpath(ROOT, "data", "insp_sitrep_scanned.csv")

## Provinces in patch order, matching PROVINCE_SOURCE_NAMES.
const PROVINCES = [
    "ituri", "nord_kivu", "sud_kivu", "haut_uele", "tshopo",
    "bas_uele", "sud_ubangi",
]
const NAMES = Dict(
    "ituri" => "ituri", "nord kivu" => "nord_kivu",
    "sud kivu" => "sud_kivu", "haut uele" => "haut_uele",
    "tshopo" => "tshopo", "bas uele" => "bas_uele",
    "sud ubangi" => "sud_ubangi"
)

## Spelling variants of a zone name, keyed on the canonical (folded) form
## of the variant and valued with the canonical form of the name used for
## the key. Every zone the table has ever printed resolves through this
## table or is its own key; the script lists the spellings behind each key
## and every key per province, so a new variant is visible rather than
## becoming a second series. Hyphen, space, case and accent differences
## ("Nia-Nia" / "Nia Nia", "Miti-Murhesa" / "Miti Murhesa") need no entry.
const ZONE_ALIASES = Dict(
    ## Seen in the archive's tables.
    "gethy" => "gety",                       # 059 onward
    "nai nia" => "nia nia",                  # 030
    "makisokisangani" => "makiso kisangani", # 064
    "makiso" => "makiso kisangani",          # wrapped name, 080 onward
    "boma" => "boma mangbetu",               # wrapped name, 080 onward
    "bambu mine" => "bambu",                 # 087 (24h-only table)
    "mungbwalu" => "mongbwalu",              # 005 (daily table)
    ## Variants the INRB-UMIE mirror's aliases.csv records from vintages
    ## this archive lacks (SitRep 057 of 10 July and the pre-018 reports).
    "nyakunde" => "nyankunde",
    "rumba" => "rimba",
    "mongbalu" => "mongbwalu",
    "tchomai" => "tchomia",
    "wanierukula" => "wanie rukula",
    "manguripa" => "manguredjipa"
)

## Rows of cases and deaths the report could not attribute to a zone. Each
## era words the row differently; all become `<province>.unallocated`.
const UNALLOCATED = (
    r"^autres zs", r"^autres zones", r"^non identifiees",
    r"^a ventiler", r"^non ventilees",
)

fold(s) = lowercase(Base.Unicode.normalize(String(s); stripmark = true))

"""
Fold a table cell to the form the names are keyed on: accents stripped,
every run of non-letters collapsed to one space. This absorbs the
hyphen/space, list-numbering and footnote-marker variation in the name
column ("1- Bunia", "Nord-Kivu*", "Nia Nia" / "Nia-Nia").
"""
function canon(s::AbstractString)
    t = replace(fold(s), r"[^a-z]+" => " ")
    return strip(t)
end

"""
Digits only: strips the thousands spaces ("1 792") and the footnote markers
("535*") the table carries. `nothing` for a cell with no digit, which covers
the NA and dash cells.
"""
function digits_only(s::AbstractString)
    d = replace(s, r"[^0-9]" => "")
    return isempty(d) ? nothing : parse(Int, d)
end

is_unallocated(c::AbstractString) = any(p -> occursin(p, c), UNALLOCATED)

zone_key(c::AbstractString) = replace(get(ZONE_ALIASES, c, c), " " => "_")

"""
The lines of Tableau 2, from its caption to its Total row, or `nothing`
with a reason when the report has no cumulative per-zone table.

The caption has been worded five ways, so it is matched on "cas", "décès"
and the province-and-zone phrase rather than on the numbering. The phrase
stops at "zone" because from SitRep 124 the caption reads "par province et
zone du 15 septembre 2026", dropping "de santé" for the date. The table
always closes with a Total row (TOTAL, TOTAL GÉNÉRAL, TOTAL NATIONAL) and
can run across a page break, so the search is long and stops at that row.
The header must name a cumulative column: SitRep 087 prints the same
caption over a 24h-flow-only table.
"""
function tableau2_lines(text::AbstractString)
    lines = split(text, '\n')
    head = findfirst(lines) do l
        f = fold(l)
        occursin("cas", f) && occursin("deces", f) &&
            (
            occursin("par province et zone", f) ||
                occursin("par zone de sante et", f) ||
                endswith(rstrip(f), "par zone de sante")
        )
    end
    if head === nothing
        ## A caption the match misses reads the same as a report that
        ## never printed the table, and the second is what gets written
        ## down. The numbering cannot tell them apart, since Tableau 2 is
        ## the alert table before SitRep 034. These two are the zone
        ## table's own: its footnote on deaths awaiting attribution, and
        ## any caption naming both a split and a zone.
        looks_present = any(lines) do l
            f = fold(l)
            occursin("ventilation pour etre classifies par zone", f) ||
                (occursin("repartition des cas", f) && occursin("zone", f))
        end
        return (
            nothing,
            looks_present ?
                "Tableau 2 is present but its caption did not match" :
                "no Tableau 2 caption",
        )
    end
    out = String[]
    for l in lines[(head + 1):min(head + 300, length(lines))]
        push!(out, String(l))
        ## The Total row carries numbers; a "Total décès" column header
        ## split onto its own line (SitRep 080) does not.
        startswith(canon(first(split(strip(l), r"\s{2,}"))), "total") &&
            occursin(r"[0-9]", l) && break
    end
    header = join(fold.(out[1:min(12, length(out))]), " ")
    (occursin("cumul", header) || occursin("cas (n)", header)) ||
        return (nothing, "Tableau 2 has no cumulative columns")
    return (out, "")
end

"""
Cumulative `(cases, deaths)` from one row, or `nothing`.

`start` is the first field that can hold a number. The Létalité percentage
is found rather than indexed, because the number of columns to its right
has changed between vintages; the cases and deaths are the two fields
before it. A row with no `%` is read from its first two plain integer
fields instead: SitReps 018-020 print the Létalité without a sign, SitRep
080 displaces the Létalité of a wrapped-name row onto the next line, and
SitRep 116 prints "1" in Buta's Létalité cell.
"""
function parse_row(parts::Vector{<:AbstractString}, start::Int)
    pct = findfirst(i -> occursin('%', parts[i]), start:length(parts))
    if pct !== nothing
        j = start + pct - 1
        j - 2 < start && return nothing
        c = digits_only(parts[j - 2])
        d = digits_only(parts[j - 1])
        (c === nothing || d === nothing) && return nothing
        return (c, d)
    end
    nums = [
        p for p in parts[start:end]
            if occursin(r"^[0-9][0-9 ]*(,[0-9]+)?$", p)
    ]
    length(nums) >= 2 || return nothing
    (occursin(',', nums[1]) || occursin(',', nums[2])) && return nothing
    return (digits_only(nums[1]), digits_only(nums[2]))
end

"""
Cumulative `(cases, deaths)` from an unallocated row, or `nothing`.

The "A ventiler" row prints NA in its cases cell (the row holds deaths in
the treatment centres not yet attributed to a zone) and a footnote marker
where the deaths figure can be displaced onto the following line, which
`next` carries. NA reads as zero. The earlier "Autres ZS" and "Non
identifiées" rows are ordinary rows.
"""
function parse_unallocated(
        parts::Vector{<:AbstractString},
        next::AbstractString
    )
    any(p -> occursin('%', p), parts) && return parse_row(parts, 2)
    cells = String[]
    for p in parts[2:end]
        if occursin(r"^NA$", strip(p))
            push!(cells, "0")
        elseif occursin(r"^[0-9][0-9 ]*$", strip(p))
            push!(cells, p)
        elseif strip(p) == "*"
            m = match(r"^\s*([0-9][0-9 ]*)\s*$", next)
            m === nothing && return nothing
            push!(cells, m[1])
        end
        length(cells) == 2 && break
    end
    length(cells) == 2 || return nothing
    return (digits_only(cells[1]), digits_only(cells[2]))
end

"""
Whether a table line is a stray name fragment: a single field of letters
with no digit, which is what a wrapped zone name ("Boma" / "Mangbetu",
"Makiso-" / "Kisangani") leaves above and below its all-number row.
"""
function is_fragment(line::AbstractString)
    s = strip(line)
    isempty(s) && return false
    occursin(r"[0-9%]", s) && return false
    return length(split(s, r"\s{2,}")) == 1
end

"""
Per-zone `(cases, deaths)` for one sitrep, keyed province → zone key, with
the printed per-province rows, the printed Total and a reason string when
nothing could be read.

A zone row belongs to the most recent province row above it. A province
name met a second time inside its own section is a zone of that name
(Tshopo has a ZS Tshopo). An all-number row takes its name from the
fragment lines around it; a bare heading with no numbers opens a province
section. The unallocated row is stored under `unallocated`.
"""
function scan_tableau2(text::AbstractString)
    lines, why = tableau2_lines(text)
    lines === nothing && return (nothing, nothing, nothing, why)
    rows = Dict{String, Dict{String, Tuple{Int, Int}}}()
    prov = Dict{String, Tuple{Int, Int}}()
    total = nothing
    cur = nothing
    consumed = falses(length(lines))
    for (i, line) in enumerate(lines)
        consumed[i] && continue
        isempty(strip(line)) && continue
        parts = split(strip(line), r"\s{2,}")
        name = parts[1]
        start = 2
        if !occursin(r"[a-z]", fold(name))
            ## No name on this row: it is a page number, a displaced cell,
            ## or a wrapped zone name whose fragments sit above and below.
            length(parts) >= 3 || continue
            frag = String[]
            if i > 1 && !consumed[i - 1] && is_fragment(lines[i - 1])
                push!(frag, strip(lines[i - 1]))
            end
            if i < length(lines) && is_fragment(lines[i + 1])
                push!(frag, strip(lines[i + 1]))
                consumed[i + 1] = true
            end
            isempty(frag) && continue
            name = join(frag, " ")
            start = 1
        end
        c = canon(name)
        if startswith(c, "sous total")
            p = get(NAMES, strip(c[11:end]), nothing)
            p === nothing && continue
            r = parse_row(parts, start)
            r === nothing || (prov[p] = r)
        elseif startswith(c, "total")
            occursin(r"[0-9]", line) || continue
            total = parse_row(parts, start)
            break
        elseif haskey(NAMES, c) && !(cur == NAMES[c] && haskey(prov, cur))
            cur = NAMES[c]
            get!(rows, cur, Dict{String, Tuple{Int, Int}}())
            r = parse_row(parts, start)
            r === nothing || (prov[cur] = r)
        elseif cur === nothing
            continue
        elseif is_unallocated(c)
            nxt = i < length(lines) ? lines[i + 1] : ""
            r = parse_unallocated(parts, nxt)
            r === nothing && continue
            occursin(r"^\s*[0-9][0-9 ]*\s*$", nxt) && (consumed[i + 1] = true)
            rows[cur]["unallocated"] = r
        else
            r = parse_row(parts, start)
            r === nothing && continue
            z = zone_key(c)
            haskey(rows[cur], z) && continue
            rows[cur][z] = r
            push!(get!(SPELLINGS, (cur, z), Set{String}()), c)
        end
    end
    isempty(rows) && return (nothing, nothing, nothing, "no zone rows parsed")
    return (rows, prov, total, "")
end

## Every spelling that fed each (province, zone key), for the report.
const SPELLINGS = Dict{Tuple{String, String}, Set{String}}()

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

zsum(z, idx) = sum(v[idx] for v in values(z); init = 0)

function main()
    isdir(PDF_DIR) || error(
        "no sitrep PDFs at $(PDF_DIR); " *
            "run `task download-sitreps` first."
    )
    Sys.which("pdftotext") === nothing &&
        error("pdftotext not found; install poppler-utils.")

    dates = sitrep_dates()
    scanned = Dict{String, Dict{String, Dict{String, Tuple{Int, Int}}}}()
    printed = Dict{String, Dict{String, Tuple{Int, Int}}}()
    totals = Dict{String, Union{Nothing, Tuple{Int, Int}}}()
    srs = Dict{String, Int}()
    unreadable = Tuple{Int, String}[]
    duplicates = String[]
    with_pdf = Set{Int}()

    for path in sort(
            filter(
                f -> endswith(f, ".pdf"),
                readdir(PDF_DIR; join = true)
            )
        )
        m = match(r"(\d+)[_-]2026", basename(path))
        m === nothing && continue
        sr = parse(Int, m[1])
        haskey(dates, sr) || continue
        push!(with_pdf, sr)

        text = read(`pdftotext -layout $path -`, String)
        rows, prov, total, why = scan_tableau2(text)
        if rows === nothing
            push!(unreadable, (sr, why))
            continue
        end
        d = dates[sr]
        ## Two files can carry one sitrep (the NNN-2026 and NNN_2026
        ## spellings, or a v2 re-issue); they must agree.
        if haskey(scanned, d) && scanned[d] != rows
            push!(duplicates, "$(basename(path)) (sitrep $(sr))")
        end
        scanned[d] = rows
        printed[d] = prov
        totals[d] = total
        srs[d] = sr
    end

    isempty(scanned) && error("no sitrep yielded a Tableau 2.")

    raw = TOML.parsefile(MANIFEST)
    hist(k) = Dict(
        String(d) => Int(v)
            for (d, v) in zip(raw[k]["dates"], raw[k]["values"])
    )
    natc, natd = hist("confirmed_case_history"), hist("confirmed_death_history")
    provblock(k) = Dict(
        String(d) => Dict(
            p => Int(raw[k][p][i])
                for p in PROVINCES if haskey(raw[k], p)
        )
            for (i, d) in enumerate(raw[k]["dates"])
    )
    provc = provblock("province_confirmed_history")
    provd = provblock("province_death_history")

    all_dates = sort(collect(keys(scanned)))
    checkable = [d for d in all_dates if haskey(natc, d) && haskey(natd, d)]
    unchecked = setdiff(all_dates, checkable)

    ## Reconciliation. For each date and province the zone rows plus the
    ## unallocated row must equal the committed province cumulative (cases
    ## and deaths), and the provinces must sum to the national totals. A
    ## date the province blocks do not carry is reconciled against the
    ## printed province rows and the national totals instead.
    keep = String[]
    national_only = String[]
    internal = Tuple{String, String}[]
    bad = Tuple{String, String}[]
    off_printed = Tuple{String, String}[]
    for d in checkable
        rows = scanned[d]
        ok = true
        via_province = haskey(provc, d) && haskey(provd, d)
        for p in union(
                keys(rows),
                via_province ? keys(provc[d]) : Set{String}()
            )
            z = get(rows, p, Dict{String, Tuple{Int, Int}}())
            s = (zsum(z, 1), zsum(z, 2))
            pr = get(printed[d], p, nothing)
            if via_province
                committed = (get(provc[d], p, 0), get(provd[d], p, 0))
                if s == committed
                    pr !== nothing && pr != s &&
                        push!(off_printed, (d, "$(p): zones $(s) printed $(pr)"))
                elseif pr !== nothing && pr == committed
                    push!(internal, (d, "$(p): zones $(s), printed row $(pr)"))
                    ok = false
                else
                    push!(
                        bad, (
                            d, "$(p): zones $(s), printed $(pr), " *
                                "committed $(committed)",
                        )
                    )
                    ok = false
                end
            elseif pr === nothing
                push!(bad, (d, "$(p): zones $(s), no printed province row"))
                ok = false
            elseif s != pr
                push!(internal, (d, "$(p): zones $(s), printed row $(pr)"))
                ok = false
            end
        end
        nat = (
            sum(zsum(z, 1) for z in values(rows)),
            sum(zsum(z, 2) for z in values(rows)),
        )
        if ok && nat != (natc[d], natd[d])
            push!(
                bad, (
                    d, "provinces sum to $(nat), national " *
                        "$((natc[d], natd[d]))",
                )
            )
            ok = false
        end
        ok || continue
        push!(keep, d)
        via_province || push!(national_only, d)
    end

    println("Per-zone confirmed cases and deaths (Tableau 2)\n")
    @printf(
        "%11s %4s | %6s %6s | %6s %6s | %5s %5s | %s\n", "date", "sr",
        "cases", "nat", "deaths", "nat", "zones", "unall", "note"
    )
    println("-"^78)
    for d in checkable
        rows = scanned[d]
        cs = sum(zsum(z, 1) for z in values(rows))
        ds = sum(zsum(z, 2) for z in values(rows))
        nz = sum(count(k -> k != "unallocated", keys(z)) for z in values(rows))
        un = sum(get(z, "unallocated", (0, 0))[1] for z in values(rows))
        note = d in keep ? (d in national_only ? "national only" : "") :
            any(t -> t[1] == d, internal) ? "table contradicts itself" :
            "MISMATCH"
        @printf(
            "%11s %4d | %6d %6s | %6d %6s | %5d %5d | %s\n",
            d, srs[d], cs, cs == natc[d] ? string(natc[d]) : "!$(natc[d])",
            ds, ds == natd[d] ? string(natd[d]) : "!$(natd[d])", nz, un, note
        )
    end

    early = sort(unique(sr for (sr, _) in unreadable if sr < 18))
    isempty(early) ||
        println("\nSitreps before 018 predate the table: ", join(early, ", "))
    later = sort(unique(t for t in unreadable if t[1] >= 18))
    if !isempty(later)
        println("\nSitreps with no readable Tableau 2:")
        for (sr, why) in later
            @printf("  %3d  %s\n", sr, why)
        end
    end
    ## A sitrep the CSV lists but the archive has no PDF for is absent from
    ## the blocks without any other trace, so it is named here.
    no_pdf = sort(setdiff(collect(keys(dates)), with_pdf))
    if !isempty(no_pdf)
        println(
            "\nSitreps in $(basename(SITREP_CSV)) with no PDF under " *
                "$(PDF_DIR) (not scanned): ",
            join((@sprintf("%03d", x) for x in no_pdf), ", ")
        )
    end
    if !isempty(duplicates)
        println(
            "\nDuplicate files that disagree with the file scanned " *
                "before them for the same sitrep: ",
            join(duplicates, ", ")
        )
    end
    if !isempty(internal)
        println(
            "\nDates whose Tableau 2 contradicts itself, the zone rows " *
                "not summing to the printed province row (not emitted):"
        )
        for (d, msg) in internal
            @printf("  %11s sitrep %3d  %s\n", d, srs[d], msg)
        end
    end
    if !isempty(off_printed)
        println(
            "\nDates where the zone rows reconcile with the committed " *
                "province value but the printed province row does not " *
                "(kept):"
        )
        for (d, msg) in off_printed
            @printf("  %11s sitrep %3d  %s\n", d, srs[d], msg)
        end
    end
    if !isempty(national_only)
        println(
            "\nDates before the province blocks begin, reconciled " *
                "against the printed province rows and the national " *
                "totals only (kept): ",
            join(national_only, ", ")
        )
    end
    if !isempty(unchecked)
        println(
            "\nDates with no national totals to reconcile against " *
                "(not emitted): ",
            join(unchecked, ", ")
        )
    end

    ## Every zone key with the spellings that fed it, so a new variant that
    ## has become a second series is visible.
    zones = Dict(p => Set{String}() for p in PROVINCES)
    for d in keep, (p, z) in scanned[d], k in keys(z)
        k == "unallocated" || push!(zones[p], k)
    end
    println(
        "\nZones per province (kept dates), with the spellings behind " *
            "any key fed by more than one:"
    )
    for p in PROVINCES
        isempty(zones[p]) && continue
        ks = sort(collect(zones[p]))
        println("  $(p) ($(length(ks))): ", join(ks, ", "))
        for k in ks
            sp = sort(collect(get(SPELLINGS, (p, k), Set{String}())))
            length(sp) > 1 && println("    $(k) <= ", join(sp, " / "))
        end
    end
    multi = [
        k for k in union(values(zones)...)
            if count(p -> k in zones[p], PROVINCES) > 1
    ]
    isempty(multi) ||
        println(
        "\nZone keys that appear under more than one province: ",
        join(multi, ", ")
    )

    if !isempty(bad)
        println("\nUnexplained disagreements:")
        for (d, msg) in bad
            @printf("  %11s sitrep %3d  %s\n", d, srs[d], msg)
        end
        println()
        error(
            "$(length(bad)) zone/province disagreement(s). The per-zone " *
                "rows plus the unallocated row partition each province, so " *
                "a mismatch that the printed province row does not explain " *
                "is a mis-parse. Not emitting the blocks."
        )
    end
    println(
        "\nAll $(length(keep)) kept dates reconcile with the province " *
            "and national confirmed case and death totals."
    )

    fmt(v) = join(v, ", ")
    dropped = sort(unique(srs[d] for (d, _) in internal))
    fmtsr(v) = join((@sprintf("%03d", x) for x in v), ", ")
    natsr = [srs[d] for d in national_only]
    println("\n\n===== paste into data/observations.toml =====\n")
    for (blk, idx, what) in (
            ("zone_confirmed_history", 1, "cases"),
            ("zone_death_history", 2, "deaths"),
        )
        println("[$(blk)]")
        println("dates = [", join(["\"$d\"" for d in keep], ", "), "]")
        for p in PROVINCES
            isempty(zones[p]) && continue
            for k in vcat(sort(collect(zones[p])), "unallocated")
                v = [
                    get(get(scanned[d], p, Dict()), k, (0, 0))[idx]
                        for d in keep
                ]
                println("$(p).$(k) = [", fmt(v), "]")
            end
        end
        println(
            "source = \"INSP situation reports, Tableau 2 " *
                "(Répartition des cas et décès confirmés par province et " *
                "zone de santé): per-health-zone cumulative confirmed " *
                "$(what) within each province. A zone absent from a " *
                "vintage's table has no confirmed $(what) yet and is " *
                "recorded as 0. `unallocated` is the report's own row of " *
                "$(what) not yet attributed to a zone (Autres ZS, Non " *
                "identifiées, A ventiler; NA is recorded as 0). Scanned " *
                "by scripts/scan_zone_tableau2.jl, which requires the " *
                "zone rows plus the unallocated row to sum exactly to " *
                "the province_$(what == "cases" ? "confirmed" : "death")" *
                "_history value on every date that block carries, and " *
                "to the printed province row and the national total on " *
                "the $(length(natsr)) dates it does not (SitReps " *
                "$(fmtsr(natsr))). SitReps $(fmtsr(dropped)) are left " *
                "out because their zone rows do not sum to their own " *
                "printed province row.\""
        )
        println()
    end
    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
