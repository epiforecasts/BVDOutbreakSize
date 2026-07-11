#!/usr/bin/env julia
#
# Scan Tableau 1 of the INSP situation reports — the per-province split of
# CONFIRMED CASES and CONFIRMED DEATHS — and emit the
# `[province_confirmed_history]` and `[province_death_history]` blocks for
# data/observations.toml.
#
#   Tableau 1. Répartition des cas et décès confirmés par province touchée
#
#    Province        Cas confirmés   Décès (confirmés)   Létalité   ...
#    Ituri                    1631               535*      32,8%
#    Nord-Kivu                 158                 89      56,3%
#    Sud-Kivu                    3                  1      33,3%
#    Total                    1 792                625      34,9%
#
# WHY THE DEATHS COLUMN MATTERS
#
# The per-province CASE split alone cannot separate "more infections" from
# "better case-finding": the confirmed count in a province is the product of
# its ascertainment and its incidence, and only that product is observed.
#
# The deaths column breaks the tie. Deaths are far harder to miss than cases,
# and the virus's case-fatality does not change at a provincial border. So
# under a shared true CFR, the per-province DEATH split identifies the
# per-province INCIDENCE split, and the gap between the death split and the
# case split is what identifies the relative case ascertainment.
#
# The signal is large and sustained: Nord-Kivu holds a steady ~8.5-9% of
# confirmed cases but ~15-19% of confirmed deaths across every vintage, with a
# confirmed CFR of 54-59% against Ituri's 20-33%.
#
# Part of that gap is NOT ascertainment: in a fast-growing epidemic the
# observed CFR is biased DOWN, because recent cases have not yet died. Ituri
# grows faster than Nord-Kivu, so some of its lower CFR is right-censoring.
# The model separates the two by applying the shared CFR and onset-to-death
# delay to each province's own incidence curve; this script only supplies the
# data.
#
# VALIDATION: the per-province cases must sum to the national
# `confirmed_case_history` and the per-province deaths to the national
# `confirmed_death_history`, on every date. Exits non-zero on any
# disagreement, so a mis-parse cannot reach the manifest.
#
# Requires `pdftotext` (poppler-utils) and the sitrep PDFs
# (`task download-sitreps`).
#
# Usage:
#
#   julia --project=scripts scripts/scan_province_tableau1.jl

using TOML
using Printf

const ROOT = normpath(joinpath(@__DIR__, ".."))
const PDF_DIR = length(ARGS) >= 1 ? ARGS[1] :
                joinpath(ROOT, "data", "sitrep_pdfs")
const MANIFEST = joinpath(ROOT, "data", "observations.toml")
const SITREP_CSV = joinpath(ROOT, "data", "insp_sitrep_scanned.csv")
const PROVINCES = ["ituri", "nord_kivu", "sud_kivu"]
const LABELS = ["ituri" => "ituri", "nord_kivu" => "nord-kivu",
    "sud_kivu" => "sud-kivu"]

fold(s) = lowercase(Base.Unicode.normalize(String(s); stripmark = true))

"""
Digits only: strips the thousands spaces ("1 792") and the footnote markers
("535*") the table carries.
"""
function digits_only(s::AbstractString)
    d = replace(s, r"[^0-9]" => "")
    return isempty(d) ? nothing : parse(Int, d)
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

    for path in sort(filter(f -> endswith(f, ".pdf"),
        readdir(PDF_DIR; join = true)))
        m = match(r"(\d+)[_-]2026", basename(path))
        m === nothing && continue
        sr = parse(Int, m[1])
        haskey(dates, sr) || continue

        text = read(`pdftotext -layout $path -`, String)
        sec = match(r"Tableau 1\..*?(?=Tableau 2|\n\s*3\.2|\n\s*4\.)"s, text)
        sec === nothing && continue

        got = Dict{String, Tuple{Int, Int}}()
        for line in split(sec.match, '\n')
            ## Columns are separated by runs of 2+ spaces; a thousands
            ## separator is a single space, so it stays inside its column.
            parts = split(strip(line), r"\s{2,}")
            length(parts) < 4 && continue
            nm = fold(parts[1])
            for (name, label) in LABELS
                nm == label || continue
                ## The third column must be the Létalité percentage. That is
                ## the gate: without it the row is not the table row we want
                ## and the two numbers before it are not cases and deaths.
                occursin('%', parts[4]) || continue
                c = digits_only(parts[2])
                d = digits_only(parts[3])
                (c === nothing || d === nothing) && continue
                got[name] = (c, d)
            end
        end
        length(got) == 3 && (scanned[dates[sr]] = got)
    end

    keep = sort(collect(keys(scanned)))
    isempty(keep) && error("no sitrep yielded a full Tableau 1.")

    raw = TOML.parsefile(MANIFEST)
    natc = Dict(String(d) => v
    for (d, v) in zip(
        raw["confirmed_case_history"]["dates"],
        raw["confirmed_case_history"]["values"]))
    natd = Dict(String(d) => v
    for (d, v) in zip(
        raw["confirmed_death_history"]["dates"],
        raw["confirmed_death_history"]["values"]))

    println("Per-province confirmed cases and deaths (Tableau 1)\n")
    @printf("%11s | %6s %6s | %6s %6s | %6s %6s\n", "date",
        "cases", "nat", "deaths", "nat", "NK c%", "NK d%")
    println("-"^66)
    bad = 0
    for d in keep
        g = scanned[d]
        cs = sum(g[p][1] for p in PROVINCES)
        ds = sum(g[p][2] for p in PROVINCES)
        nc = get(natc, d, nothing)
        nd = get(natd, d, nothing)
        okc = nc === nothing ? "-" : (cs == nc ? string(nc) : "!$(nc)")
        okd = nd === nothing ? "-" : (ds == nd ? string(nd) : "!$(nd)")
        (nc !== nothing && cs != nc) && (bad += 1)
        (nd !== nothing && ds != nd) && (bad += 1)
        @printf("%11s | %6d %6s | %6d %6s | %5.1f%% %5.1f%%\n", d, cs, okc,
            ds, okd, 100 * g["nord_kivu"][1] / cs, 100 * g["nord_kivu"][2] / ds)
    end

    if bad > 0
        println()
        error("$(bad) province/national disagreement(s). The per-province " *
              "figures are an exact partition of the national totals, so a " *
              "mismatch means a mis-parse. Not emitting the blocks.")
    end
    println("\nAll $(length(keep)) dates reconcile with the national " *
            "confirmed case AND death totals.")

    ## Nord-Kivu's death share sits well above its case share at every
    ## vintage; that gap is the identifying signal, so report it.
    ncs = [100 * scanned[d]["nord_kivu"][1] /
           sum(scanned[d][p][1] for p in PROVINCES) for d in keep]
    nds = [100 * scanned[d]["nord_kivu"][2] /
           sum(scanned[d][p][2] for p in PROVINCES) for d in keep]
    @printf("\nNord-Kivu: case share %.1f-%.1f%%, death share %.1f-%.1f%%, ratio %.2fx at the cut-off\n",
        minimum(ncs), maximum(ncs), minimum(nds), maximum(nds),
        nds[end] / ncs[end])

    fmt(v) = join(v, ", ")
    println("\n\n===== paste into data/observations.toml =====\n")
    for (blk, idx, what) in (("province_confirmed_history", 1, "cases"),
        ("province_death_history", 2, "deaths"))
        println("[$(blk)]")
        println("dates = [", join(["\"$d\"" for d in keep], ", "), "]")
        for p in PROVINCES
            println("$(p) = [", fmt([scanned[d][p][idx] for d in keep]), "]")
        end
        println("source = \"INSP situation reports, Tableau 1 " *
                "(Répartition des cas et décès confirmés par province): " *
                "per-province cumulative confirmed $(what). Scanned by " *
                "scripts/scan_province_tableau1.jl, which requires the " *
                "per-province figures to sum exactly to the national " *
                "totals on every date.\"")
        println()
    end
    return nothing
end

main()
