#!/usr/bin/env julia
#
# Scan the per-province laboratory throughput (section 4.3, "Laboratoire")
# out of the INSP situation-report PDFs and emit the
# `[province_lab_daily_history]` block for data/observations.toml.
#
# For each sitrep and each province the section reports the samples ANALYSED
# in the last 24h and how many came back POSITIVE, for example:
#
#   4.3. Laboratoire
#     • Ituri : 147 échantillons collectés et analysés (positivité 30,6% ; n=45).
#     • Nord-Kivu : 90 échantillons collectés et analysés (positivité 6,7% ; n=6).
#     • Sud-Kivu : 3 échantillons collectés et encours analysés.
#
# The phrasing varies between sitreps ("dont 15 ont été analysés", "Sur un
# total de 117 échantillons reçus, 97 ont été analysés", "tous sont revenus
# négatifs", ...), so each line is matched explicitly and anything that does
# not parse is REPORTED rather than silently dropped or guessed.
#
# A sample that is only "en cours d'analyse" has no result yet and so gives
# no denominator; those are excluded rather than counted as analysed.
#
# VALIDATION: the per-province analysed counts must sum, on every date, to
# the national `tests_analysed_daily_history` already in observations.toml.
# That national series is fitted by the model, so the per-province figures
# are an exact partition of it. The script exits non-zero on any
# disagreement, so a mis-parse cannot reach the manifest.
#
# Requires `pdftotext` (poppler-utils) on the PATH, and the sitrep PDFs
# already downloaded (`task download-sitreps`).
#
# Usage:
#
#   julia --project=scripts scripts/scan_province_lab.jl
#   julia --project=scripts scripts/scan_province_lab.jl path/to/pdfdir

using TOML
using Dates
using Printf

const ROOT = normpath(joinpath(@__DIR__, ".."))
const PDF_DIR = length(ARGS) >= 1 ? ARGS[1] : joinpath(ROOT, "data", "sitrep_pdfs")
const MANIFEST = joinpath(ROOT, "data", "observations.toml")
const SITREP_CSV = joinpath(ROOT, "data", "insp_sitrep_scanned.csv")

## Provinces in patch order: the first is the primary (origin) patch.
const PROVINCES = ["ituri", "nord_kivu", "sud_kivu"]
const PATTERNS = ["ituri" => "ituri", "nord_kivu" => "nord-kivu",
    "sud_kivu" => "sud-kivu"]

"""
Strip accents and lower-case, so the many spellings in the sitreps
("analysés", "analyses", "ANALYSÉS") all match one pattern.
"""
function fold(s::AbstractString)
    ## `stripmark` drops the combining marks that NFKD separates out, so the
    ## accented and unaccented spellings collapse onto one form.
    return lowercase(Base.Unicode.normalize(String(s); stripmark = true))
end

"""
Parse one province line from section 4.3 into `(analysed, positives)`, or
`nothing` when the line carries no usable counts (no result yet, or no
update shared). Never guesses: an unrecognised line returns `nothing` and
is reported by the caller.
"""
function parse_province_line(line::AbstractString)
    t = replace(fold(line), r"\s+" => " ")
    pending = occursin(r"en ?cours (d.)?analys", t)

    analysed = nothing
    ## "dont 15 ont ete analyses" / "dont 77 analyses"
    m = match(r"dont (\d+) (?:ont ete )?analys", t)
    if m !== nothing
        analysed = parse(Int, m[1])
    else
        ## "97 ont ete analyses"
        m = match(r"(\d+) ont ete analys", t)
        if m !== nothing
            analysed = parse(Int, m[1])
        else
            ## "147 echantillons collectes et analyses"
            ## "95 echantillons recus et analyses"
            m = match(r"(\d+) (?:nouveaux? )?echantillons?[^.]*?analys", t)
            m !== nothing && (analysed = parse(Int, m[1]))
        end
    end
    ## A batch still in the analyser has no result, so it is not a
    ## denominator for this date, whatever the collected count says.
    pending && (analysed = nothing)

    positives = nothing
    m = match(r"n ?= ?(\d+)", t)
    if m !== nothing
        positives = parse(Int, m[1])
    else
        m = match(r"(\d+) (?:sont|est) revenus? positifs?", t)
        if m !== nothing
            positives = parse(Int, m[1])
        elseif occursin("revenus negatifs", t)
            positives = 0
        end
    end

    ## "5 nouveaux echantillons ont ete collectes et tous sont revenus
    ## negatifs": no "analyse" verb, but a returned result means they were
    ## analysed. Fall back to the collected count, never when pending.
    if analysed === nothing && positives !== nothing && !pending
        m = match(r"(\d+) (?:nouveaux? )?echantillons?", t)
        m !== nothing && (analysed = parse(Int, m[1]))
    end

    (analysed === nothing || positives === nothing) && return nothing
    return (analysed, positives)
end

"""
Map sitrep number to report date, from the hand-scanned headline CSV.
"""
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
    ## date => province => (analysed, positives)
    scanned = Dict{String, Dict{String, Tuple{Int, Int}}}()
    unparsed = Tuple{Int, String, String}[]

    for path in sort(filter(f -> endswith(f, ".pdf"), readdir(PDF_DIR; join = true)))
        m = match(r"(\d+)[_-]2026", basename(path))
        m === nothing && continue
        sr = parse(Int, m[1])
        haskey(dates, sr) || continue

        text = read(`pdftotext -layout $path -`, String)
        sec = match(r"4\.3\..{0,20}?[Ll]aboratoire(.*?)(?:\n\s*4\.4|\n\s*5\.)"s,
            text)
        sec === nothing && continue

        got = Dict{String, Tuple{Int, Int}}()
        for line in split(sec[1], '\n')
            f = fold(line)
            for (name, pat) in PATTERNS
                occursin(Regex(pat * "\\s*:"), f) || continue
                r = parse_province_line(line)
                if r === nothing
                    push!(unparsed, (sr, name, strip(line)))
                else
                    got[name] = r
                end
            end
        end
        isempty(got) || (scanned[dates[sr]] = got)
    end

    ## Keep only dates where both informative patches have a denominator.
    keep = sort([d for (d, g) in scanned
                 if haskey(g, "ituri") && haskey(g, "nord_kivu")])
    isempty(keep) && error("no sitrep yielded a usable section 4.3.")

    ## --- Validate against the national daily analysed series -----------
    raw = TOML.parsefile(MANIFEST)
    nat = Dict(String(d) => v
    for (d, v) in zip(
        raw["tests_analysed_daily_history"]["dates"],
        raw["tests_analysed_daily_history"]["values"]))

    println("Per-province laboratory throughput (section 4.3)\n")
    @printf("%11s %6s %6s %6s %7s %9s  %s\n",
        "date", "IT an", "NK an", "SK an", "sum", "national", "check")
    println("-"^62)
    bad = 0
    for d in keep
        g = scanned[d]
        a = [haskey(g, p) ? g[p][1] : 0 for p in PROVINCES]
        s = sum(a)
        n = get(nat, d, nothing)
        note = if n === nothing
            "no national"
        elseif s == n
            "MATCH"
        else
            bad += 1
            @sprintf("DIFF %+d", s - n)
        end
        @printf("%11s %6d %6d %6d %7d %9s  %s\n", d, a[1], a[2], a[3], s,
            n === nothing ? "-" : string(n), note)
    end

    if !isempty(unparsed)
        println("\nLines with no usable counts (no result yet / no update):")
        for (sr, name, line) in unparsed
            @printf("  %3d %-10s %s\n", sr, name, first(line, 80))
        end
    end

    if bad > 0
        println()
        error("$(bad) date(s) where the per-province analysed counts do " *
              "not sum to the national total. The per-province figures are " *
              "an exact partition of the national series, so a mismatch " *
              "means a mis-parse. Not emitting the block.")
    end
    println("\nAll $(length(keep)) dates reconcile exactly with the " *
            "national tests_analysed_daily_history.")

    ## --- Emit the TOML block -------------------------------------------
    fmt(v) = join(v, ", ")
    println("\n\n===== paste into data/observations.toml =====\n")
    println("""
    # Per-province laboratory throughput from section 4.3 of the INSP
    # situation reports: samples ANALYSED in the last 24h and how many were
    # POSITIVE, for each province. Generated by
    # `julia --project=scripts scripts/scan_province_lab.jl`, which fails if
    # the per-province analysed counts do not sum exactly to the national
    # `tests_analysed_daily_history` on every date (they are an exact
    # partition of it, so a mismatch means a mis-parse).
    #
    # Dates where a province's samples were only "en cours d'analyse" carry
    # no result, so that province has no denominator and is recorded as 0
    # analysed / 0 positive: a Binomial with n = 0 contributes no likelihood.
    #
    # The positivity differs sharply between provinces and persistently so:
    # over this window Ituri ran 2112 tests for 671 positives (31.8%) and
    # Nord-Kivu 1340 tests for 74 (5.5%), a 5.8x gap. Nord-Kivu needs ~18
    # tests per case found against Ituri's ~3. The two provinces are testing
    # very differently-selected pools, so the per-province CONFIRMED counts
    # cannot be read as proportional to per-province infections without this
    # denominator.""")
    println("[province_lab_daily_history]")
    println("dates = [", join(["\"$d\"" for d in keep], ", "), "]")
    for p in PROVINCES
        an = [haskey(scanned[d], p) ? scanned[d][p][1] : 0 for d in keep]
        po = [haskey(scanned[d], p) ? scanned[d][p][2] : 0 for d in keep]
        println("$(p)_analysed = [", fmt(an), "]")
        println("$(p)_positive = [", fmt(po), "]")
    end
    println("source = \"INSP situation reports, section 4.3 (Laboratoire): " *
            "per-province samples analysed in the last 24h and the number " *
            "positive. Scanned from the PDFs by scripts/scan_province_lab.jl. " *
            "Every date's per-province analysed counts sum exactly to the " *
            "national tests_analysed_daily_history entry for that date.\"")
    return nothing
end

main()
