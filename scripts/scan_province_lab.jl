#!/usr/bin/env julia
#
# Scan the per-province laboratory throughput out of the INSP situation-report
# PDFs and emit the `[province_lab_daily_history]` block for
# data/observations.toml.
#
# For each sitrep and each province the laboratory section reports the samples
# analysed in the last 24h and how many came back positive, for example:
#
#   1.3. Laboratoire
#     • Ituri : 61 nouveaux résultats positifs (33 vivants et 28 décès) sur
#       342 nouveaux échantillons reçus et analysés (positivité : 17,8 %).
#     • Nord-Kivu : 28 nouveaux résultats positifs sur 157 échantillons
#       analysés (positivité de 17,8%).
#     • Bas-Uélé : 1 échantillon reçu et testé, le résultat est revenu négatif.
#
# Nothing about that section is stable. It has been numbered 4.3, 3.2 and 1.3,
# so it is found by its title rather than its number; it has grown from three
# provinces to six; its bullets wrap across lines, so an entry runs from a
# province name to the next bullet; and the phrasing of the counts changes
# from vintage to vintage ("dont 15 ont été analysés", "Sur un total de 117
# échantillons reçus, 97 ont été analysés", "tous sont revenus négatifs",
# ...). Each phrasing is matched explicitly and anything that does not parse
# is reported rather than silently dropped or guessed.
#
# A sample that is only "en cours d'analyse", or whose results have not been
# rendered, has no result yet and so gives no denominator; those are excluded
# rather than counted as analysed. Cumulative asides ("le cumul provincial
# s'établit à 3 positifs sur 296 échantillons analysés") are dropped before
# parsing, so a provincial cumulative cannot be read as a daily count.
#
# From SitRep 084 (6 August) the section gives only a national 24h total in
# prose, with no per-province bullets at all. Those vintages yield no
# per-province entry and are skipped, until the bullets return at SitRep 096
# (18 August).
#
# Validation. The per-province analysed counts must sum, on every date, to
# the national `tests_analysed_daily_history` already in observations.toml.
# That national series is fitted by the model, so the per-province figures
# are an exact partition of it, and the reconciliation is what admits a
# vintage. The script exits non-zero on any disagreement, so a mis-parse
# cannot reach the manifest.
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
const PDF_DIR = length(ARGS) >= 1 ? ARGS[1] :
                joinpath(ROOT, "data", "sitrep_pdfs")
const MANIFEST = joinpath(ROOT, "data", "observations.toml")
const SITREP_CSV = joinpath(ROOT, "data", "insp_sitrep_scanned.csv")

## Provinces in patch order: the first is the primary (origin) patch.
const PROVINCES = ["ituri", "nord_kivu", "sud_kivu", "haut_uele", "tshopo",
    "bas_uele"]
## A bullet names its province then a colon. The names are matched on the
## folded text, so the hyphen/space and accent variants all collapse here.
const ENTRY = Regex("^[^0-9a-z]*(ituri|nord[ -]?kivu|sud[ -]?kivu|" *
                    "haut[ -]?uele|tshopo|bas[ -]?uele)\\s*:+")
const KEYS = Dict("ituri" => "ituri", "nordkivu" => "nord_kivu",
    "sudkivu" => "sud_kivu", "hautuele" => "haut_uele",
    "tshopo" => "tshopo", "basuele" => "bas_uele")
## Folded spellings of each province, to spot a sentence that has moved on
## to another province.
const NAMES = Dict("ituri" => "ituri", "nord-kivu" => "nord_kivu",
    "nord kivu" => "nord_kivu", "sud-kivu" => "sud_kivu",
    "sud kivu" => "sud_kivu", "haut-uele" => "haut_uele",
    "haut uele" => "haut_uele", "tshopo" => "tshopo",
    "bas-uele" => "bas_uele", "bas uele" => "bas_uele")

## A batch of samples goes by several names.
const SAMPLES = "(?:echantillons?|swabs?|prelevements?)"

## The phrasings a positive count comes in, specific readings first.
const POSITIVES = Regex.(["n ?= ?(\\d+)",
    "(\\d+) ?nouveaux? resultats? positifs?",
    "dont (\\d+) ?(?:sont )?(?:revenus? )?positifs?",
    "parmi lesquels (\\d+) ?(?:sont )?(?:revenus? )?positifs?",
    "soit (\\d+) ?positifs?",
    "(\\d+) ?(?:sont|est) revenus? positifs?",
    "(\\d+) ?echantillons? testes? positifs?",
    "(\\d+) ?resultats? positifs?",
    "(\\d+) ?positifs?"])

## A batch split between two laboratories is analysed in full, so both
## halves count: "132 reçus, dont 47 analysés à Beni et 85 à Butembo".
const SPLIT_LABS = Regex("dont (\\d+) (?:ont ete )?(?:analyses?|testes?)" *
                         " a [^;]*? et (\\d+) a ")

## The phrasings an analysed count comes in. Each is tempered so that it
## cannot read across a pending phrase: a count is a denominator only if
## its own analysis is done.
const ANALYSED = Regex.([
    "dont (\\d+) (?:[a-z]+ ){0,2}?(?:ont ete |sont )?(?:analys|test)",
    "sur (?:les )?(\\d+) (?:nouveaux? )?(?:echantillons?|swabs?)",
    "et (\\d+) (?:ont ete )?(?:analyses?|testes?)",
    "(\\d+) ont ete (?:analyses?|testes?)",
    "(\\d+) (?:nouveaux? )?" * SAMPLES *
    "(?:(?!en ?cours|en attente)[^;])*?(?:analys|test)",
    "(\\d+) ?(?:echantillons? )?testes? positifs?",
    "(\\d+) (?:ont ete )?(?:analyses?|testes?)\\b(?! (?:restent )?en cours)"])

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
The lines of the laboratory section.

The section number has changed twice, so the heading is matched on the word
"Laboratoire" alone and the section runs to the next heading that is not one
of its own subsections.
"""
function lab_lines(text::AbstractString)
    lines = split(text, '\n')
    ## A heading is a section number and a capitalised title. The capital
    ## matters: without it a wrapped sentence that happens to start with a
    ## number ("13 laboratoires planifiés ...") ends the section early.
    head = r"^[^0-9\p{L}]*(\d+(?:\.\d+)*)\.?\s*[-—–]?\s*(\p{Lu}.*)$"
    start, num = nothing, ""
    for (i, l) in enumerate(lines)
        m = match(head, l)
        m === nothing && continue
        startswith(fold(m[2]), "laboratoire") || continue
        start, num = i, m[1]
        break
    end
    start === nothing && return String[]
    out = String[]
    for l in lines[(start + 1):min(start + 40, length(lines))]
        f = fold(l)
        (startswith(strip(f), "tableau ") || startswith(strip(f), "figure ")) &&
            break
        m = match(head, l)
        if m !== nothing && !startswith(m[1], num)
            break
        end
        push!(out, String(l))
    end
    return out
end

"""
The laboratory section split into one block of text per province.

A bullet wraps across as many lines as it needs, so an entry absorbs every
following line that does not itself start a new bullet.
"""
function province_entries(lines::Vector{String})
    got = Dict{String, String}()
    order = String[]
    current = nothing
    for line in lines
        isempty(strip(line)) && continue
        f = fold(line)
        m = match(ENTRY, f)
        if m !== nothing
            name = KEYS[replace(m[1], r"[ -]" => "")]
            current = haskey(got, name) ? nothing : name
            if current !== nothing
                got[name] = f[(m.offset + ncodeunits(m.match)):end]
                push!(order, name)
            end
        elseif occursin(r"^\s*[^0-9a-z\s]", f)
            ## A bullet that names no province is general laboratory news,
            ## not a continuation of the province above it.
            current = nothing
        elseif current !== nothing
            got[current] = got[current] * " " * f
        end
    end
    return got, order
end

"""
`(analysed, positives)` for one province entry, or `:unparsed` when the
phrasing is not recognised.

A province that states no completed analysis, whether because its samples
are still in the analyser, because no result was rendered, or because it
reports positives without a denominator, contributes `(0, 0)`: a Binomial
with n = 0 adds no likelihood, and the national reconciliation would catch
a denominator wrongly dropped this way. Never guesses the other direction:
a denominator with no readable numerator is reported by the caller and its
whole vintage is left out.
"""
function parse_province_entry(entry::AbstractString, name::AbstractString)
    t = replace(entry, r"\s+" => " ")
    clauses = split(t, r"[;.] ")
    ## Provincial cumulatives and the closing narrative, which recaps the
    ## other provinces, sit in the same block as the daily counts. Neither
    ## is about this province today, so both go before any number is read.
    others = [k for (k, v) in NAMES if v != name]
    keepc(c) = !occursin("cumul", c) &&
               !any(o -> occursin(o, c), others)
    clauses = [c for c in clauses if keepc(c)]
    text = join(clauses, " ; ")

    ## A batch still in the analyser, or one whose results were not
    ## rendered, is no denominator for this date whatever the collected
    ## count says.
    pending(c) = occursin(r"en ?cours (?:d.)?analys", c) ||
                 occursin(r"analyses? (?:restent )?en cours", c) ||
                 occursin(r"en attente", c) ||
                 occursin(r"resultats? attendus", c) ||
                 occursin(r"n.a ete rapporte", c) ||
                 occursin(r"n.a ete rendu", c)

    ## "Aucun nouveau résultat positif" reports that nothing came back, not
    ## a completed batch with no positives, so it blocks the fallback below
    ## while still fixing the positives at zero.
    norendu = occursin(r"aucun nouveau resultat", text)
    ## A batch reported as returned negative has been analysed, however
    ## loosely the collection itself is worded.
    completed = !norendu &&
                (occursin(r"revenus? negatifs?", text) ||
                 occursin(r"tous negatifs?", text) ||
                 occursin(r"[,(] ?negatif", text) ||
                 occursin(r"aucun[^;]{0,30}positif", text))

    positives = nothing
    for re in POSITIVES
        m = match(re, text)
        m === nothing && continue
        positives = parse(Int, m[1])
        break
    end
    (norendu || completed) && positives === nothing && (positives = 0)
    ## "le résultat est revenu positif" on a single sample.
    positives === nothing && occursin(r"resultat est revenu positif", text) &&
        (positives = 1)

    analysed = nothing
    for c in clauses
        m = match(SPLIT_LABS, c)
        if m !== nothing
            analysed = parse(Int, m[1]) + parse(Int, m[2])
            break
        end
        for re in ANALYSED
            m = match(re, c)
            m === nothing && continue
            analysed = parse(Int, m[1])
            break
        end
        analysed === nothing || break
    end

    ## "5 nouveaux echantillons ont ete collectes et tous sont revenus
    ## negatifs": no analysis verb, but a returned result means they were
    ## analysed, as does a stated positive count. Fall back to the collected
    ## count, never from a clause whose results have not come back.
    if analysed === nothing && !norendu &&
       (completed || (positives !== nothing && positives > 0))
        for c in clauses
            pending(c) && continue
            m = match(Regex("(\\d+) (?:nouveaux? )?" * SAMPLES), c)
            m === nothing && continue
            analysed = parse(Int, m[1])
            break
        end
    end

    ## No completed analysis stated: the province contributes no denominator,
    ## and so no numerator either.
    analysed === nothing && return (0, 0)
    positives === nothing && return :unparsed
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
    srs = Dict{String, Int}()
    unparsed = Tuple{Int, String, String}[]
    noentries = Int[]

    for path in sort(filter(f -> endswith(f, ".pdf"),
        readdir(PDF_DIR; join = true)))
        m = match(r"(\d+)[_-]2026", basename(path))
        m === nothing && continue
        sr = parse(Int, m[1])
        haskey(dates, sr) || continue

        text = read(`pdftotext -layout $path -`, String)
        entries, order = province_entries(lab_lines(text))
        if isempty(order)
            push!(noentries, sr)
            continue
        end

        got = Dict{String, Tuple{Int, Int}}()
        ok = true
        for name in order
            r = parse_province_entry(entries[name], name)
            if r === :unparsed
                push!(unparsed, (sr, name, first(entries[name], 90)))
                ok = false
            else
                got[name] = r
            end
        end
        ## A vintage with an unreadable bullet cannot be reconciled, because
        ## that province's contribution to the national total is unknown.
        ok || continue
        scanned[dates[sr]] = got
        srs[dates[sr]] = sr
    end

    isempty(scanned) && error("no sitrep yielded a usable laboratory " *
                              "section.")

    ## --- Validate against the national daily analysed series -----------
    raw = TOML.parsefile(MANIFEST)
    nat = Dict(String(d) => v
    for (d, v) in zip(
        raw["tests_analysed_daily_history"]["dates"],
        raw["tests_analysed_daily_history"]["values"]))

    analysed(g) = sum(haskey(g, p) ? g[p][1] : 0 for p in PROVINCES)
    all_dates = sort(collect(keys(scanned)))
    checkable = [d for d in all_dates if haskey(nat, d)]
    unchecked = setdiff(all_dates, checkable)
    keep = [d for d in checkable if analysed(scanned[d]) == nat[d]]
    bad = setdiff(checkable, keep)

    println("Per-province laboratory throughput (Laboratoire section)\n")
    @printf("%11s %4s %6s %6s %6s %6s %6s %6s %7s %9s  %s\n", "date", "sr",
        "IT", "NK", "SK", "HU", "TS", "BU", "sum", "national", "check")
    println("-"^88)
    for d in checkable
        g = scanned[d]
        a = [haskey(g, p) ? g[p][1] : 0 for p in PROVINCES]
        s = analysed(g)
        n = nat[d]
        note = s == n ? "match" : @sprintf("DIFF %+d", s - n)
        @printf("%11s %4d %6d %6d %6d %6d %6d %6d %7d %9d  %s\n", d, srs[d],
            a[1], a[2], a[3], a[4], a[5], a[6], s, n, note)
    end

    if !isempty(noentries)
        println(
            "\nSitreps whose laboratory section has no per-province " *
            "bullet: ", join(noentries, ", "))
    end
    if !isempty(unparsed)
        println("\nBullets with no usable counts (whole vintage left out):")
        for (sr, name, line) in unparsed
            @printf("  %3d %-10s %s\n", sr, name, strip(line))
        end
    end
    if !isempty(unchecked)
        println(
            "\nDates with no national analysed total to reconcile " *
            "against (not emitted): ",
            join(unchecked, ", "))
    end

    if !isempty(bad)
        println()
        error("$(length(bad)) date(s) where the per-province analysed " *
              "counts do not sum to the national total " *
              "($(join(bad, ", "))). The per-province figures are " *
              "an exact partition of the national series, so a mismatch " *
              "means a mis-parse. Not emitting the block.")
    end
    println("\nAll $(length(keep)) dates reconcile exactly with the " *
            "national tests_analysed_daily_history.")

    ## --- Emit the TOML block -------------------------------------------
    fmt(v) = join(v, ", ")
    println("\n\n===== paste into data/observations.toml =====\n")
    println("""
    # Per-province laboratory throughput from the laboratory section of the
    # INSP situation reports: samples analysed in the last 24h and how many
    # were positive, for each province. Generated by
    # `julia --project=scripts scripts/scan_province_lab.jl`, which fails if
    # the per-province analysed counts do not sum exactly to the national
    # `tests_analysed_daily_history` on every date (they are an exact
    # partition of it, so a mismatch means a mis-parse).
    #
    # Dates where a province's samples were only "en cours d'analyse", and
    # provinces the section does not mention, carry no result, so they have
    # no denominator and are recorded as 0 analysed / 0 positive: a Binomial
    # with n = 0 contributes no likelihood.
    #
    # The positivity differs sharply between provinces and persistently so.
    # Nord-Kivu needs several times as many tests per case found as Ituri.
    # The provinces are testing very differently-selected pools, so the
    # per-province confirmed counts cannot be read as proportional to
    # per-province infections without this denominator.""")
    println("[province_lab_daily_history]")
    println("dates = [", join(["\"$d\"" for d in keep], ", "), "]")
    for p in PROVINCES
        an = [haskey(scanned[d], p) ? scanned[d][p][1] : 0 for d in keep]
        po = [haskey(scanned[d], p) ? scanned[d][p][2] : 0 for d in keep]
        println("$(p)_analysed = [", fmt(an), "]")
        println("$(p)_positive = [", fmt(po), "]")
    end
    println("source = \"INSP situation reports, laboratory section " *
            "(4.3, then 3.2, then 1.3 Laboratoire): per-province samples " *
            "analysed in the last 24h and the number positive. Scanned " *
            "from the PDFs by scripts/scan_province_lab.jl. Every date's " *
            "per-province analysed counts sum exactly to the national " *
            "tests_analysed_daily_history entry for that date.\"")
    return nothing
end

main()
