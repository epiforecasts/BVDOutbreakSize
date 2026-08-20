#!/usr/bin/env julia
#
# Report how many INSP situation reports we have not yet recorded.
#
# The DRC INSP publishes each MVE-17 SitRep at
# https://insp.cd/ebola-17eme-epidemie/ (posts under the `sitrep`
# category). This script lists the published SitReps straight from INSP
# via the WordPress REST API, compares them with the latest report already
# transcribed in data/insp_sitrep_scanned.csv, and prints the gap plus the
# INSP post page for each missing report.
#
# Run it before an analysis refresh (and it is what the maintenance bot
# calls) so the manifest never drifts weeks behind again.
#
# Usage:
#
#   julia --project=scripts scripts/check_new_sitreps.jl
#
# Notes:
#  - INSP blocks some default user agents (HTTP 403); this script sends a
#    browser User-Agent, which returns HTTP 200. INSP is the primary source
#    both for currency checks and for the PDFs themselves
#    (download_sitreps.jl): the INRB-UMIE mirror usually lags INSP by a
#    report or two, has dropped individual vintages, and does not carry the
#    richer `analytique` PDFs, so it is a cross-check only
#    (confirm_insp_data.jl).
#  - Each printed post page embeds the real PDF in a `pdfemb-data` base64
#    blob (a base64 of {"url": ...}); data/README.md documents the one-line
#    decode to get a direct, fetchable PDF link.
#  - Exit code is 0 when up to date, 1 when one or more reports are missing
#    (so a caller / cron can branch on it). Exit code is unaffected by the
#    mirror-lead warning below - a mirror lead never licenses treating a
#    SitRep as missing when insp.cd itself has published nothing new.
#
# Mirror-lead check: insp.cd being "up to date" only means no new
# SitRep-numbered post/PDF exists yet. It does not mean the outbreak has
# stopped moving - INSP's own SitRep pipeline can fall behind INRB-UMIE's
# internal transcription pipeline for days at a time, and a silent "up to
# date" on a day like that reads as nothing happening when the response has
# in fact continued, just without a public PDF yet. So after the insp.cd
# comparison, this script also fetches the mirror's national cumulative
# confirmed-case CSV (one extra request) and compares its latest date
# against the manifest's `as_of_date`. If the mirror is ahead, it prints a
# loud warning - it does NOT change the exit code and does NOT license
# updating data/observations.toml from the mirror alone (see
# data/README.md's insp.cd-unreachable clause): that still requires either
# an insp.cd PDF or an explicit decision to record a mirror-only point. The
# point of the check is visibility, not automation. The fetch is best-effort
# and never allowed to crash the script - see check_mirror_lead().

using Dates
using Downloads

const UA = "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) " *
           "AppleWebKit/605.1.15"
const LIST_URL = "https://insp.cd/wp-json/wp/v2/posts?" *
                 "search=sitrep&per_page=100&_fields=id,slug,date"
const SCANNED = joinpath(@__DIR__, "..", "data", "insp_sitrep_scanned.csv")
const OBSERVATIONS = joinpath(@__DIR__, "..", "data", "observations.toml")
const MIRROR_CASES_URL = "https://raw.githubusercontent.com/INRB-UMIE/" *
                         "BDBV2026-Data/main/data/insp_sitrep/processed/" *
                         "insp_sitrep__national_cumulative_confirmed_cases__daily.csv"

fetch(url) = sprint() do io
    Downloads.download(url, io; headers = ["User-Agent" => UA])
end

"Highest SitRep number already transcribed in the scanned CSV."
function latest_recorded()
    isfile(SCANNED) || error("missing $SCANNED")
    best = 0
    for line in Iterators.drop(eachline(SCANNED), 1)  # skip header
        m = match(r"^\"?(\d+)", line)                 # numeric prefix
        m === nothing && continue
        best = max(best, parse(Int, m.captures[1]))
    end
    return best
end

"Published SitReps at INSP as (number, slug), newest first."
function published_sitreps()
    body = fetch(LIST_URL)
    out = Tuple{Int, String}[]
    for m in eachmatch(
        r"\"slug\":\"(sitrep[^\"]*)\"", body)
        slug = m.captures[1]
        num = match(r"-n0*(\d+)", slug)
        num === nothing && continue
        push!(out, (parse(Int, num.captures[1]), slug))
    end
    return sort!(unique(out); by = first, rev = true)
end

"Manifest `as_of_date` as a Date, read without a TOML dependency."
function manifest_as_of_date()
    for line in eachline(OBSERVATIONS)
        m = match(r"^as_of_date\s*=\s*\"(\d{4}-\d{2}-\d{2})\"", line)
        m === nothing || return Date(m.captures[1])
    end
    error("as_of_date not found in $OBSERVATIONS")
end

"Latest (date, value) row of the mirror's national cumulative confirmed-case CSV."
function mirror_latest_cases()
    body = fetch(MIRROR_CASES_URL)
    latest = nothing
    for line in Iterators.drop(eachline(IOBuffer(body)), 1)
        cols = split(line, ',')
        length(cols) < 3 && continue
        cols[3] == "ND" && continue
        d = tryparse(Date, cols[2])
        v = tryparse(Int, cols[3])
        (d === nothing || v === nothing) && continue
        (latest === nothing || d > latest[1]) && (latest = (d, v))
    end
    return latest
end

"Print a loud but non-fatal warning if the mirror is ahead of the manifest
even though insp.cd itself has nothing new. See the module-level comment
above for why this exists and what it deliberately does not do."
function check_mirror_lead()
    as_of = manifest_as_of_date()
    mirror = try
        mirror_latest_cases()
    catch e
        println("\n(Mirror-lead check skipped: could not fetch the ",
            "INRB-UMIE mirror CSV - ", sprint(showerror, e), ".)")
        return
    end
    if mirror === nothing
        println("\n(Mirror-lead check skipped: could not parse the ",
            "INRB-UMIE mirror CSV.)")
        return
    end
    mirror_date, mirror_value = mirror
    if mirror_date > as_of
        println("\n", "!"^70)
        println("WARNING: the INRB-UMIE mirror is AHEAD of the manifest ",
            "even though insp.cd has no new SitRep.")
        println("  Manifest as_of_date:       ", as_of)
        println("  Mirror latest cumulative:  ", mirror_date, " = ",
            mirror_value, " confirmed cases")
        println("This does NOT mean a SitRep was missed - insp.cd was ",
            "checked directly above and has nothing newer. It means the ",
            "response may be continuing without a public PDF yet (INRB-UMIE ",
            "is INSP's own modelling unit and may have line-list access ",
            "ahead of publication). Investigate before treating today as ",
            "fully quiet; do not update data/observations.toml from this ",
            "value alone without recording it as mirror-only provenance.")
        println("!"^70)
    else
        println("\nMirror check: no lead (mirror latest ", mirror_date,
            " <= manifest as_of_date ", as_of, ").")
    end
end

function main()
    recorded = latest_recorded()
    pubs = published_sitreps()
    isempty(pubs) &&
        error("no SitReps parsed from INSP (rate limited or markup change?)")
    latest = first(pubs)[1]
    missing_reports = sort([p for p in pubs if first(p) > recorded];
        by = first)

    println("Latest recorded (data/insp_sitrep_scanned.csv): SitRep ",
        lpad(recorded, 3, '0'))
    println("Latest published (insp.cd):                     SitRep ",
        lpad(latest, 3, '0'))

    if isempty(missing_reports)
        println("\nUp to date - no new SitReps to record.")
        check_mirror_lead()
        return 0
    end

    println("\n", length(missing_reports),
        " SitRep(s) not yet recorded:\n")
    for (num, slug) in missing_reports
        println("  SitRep ", lpad(num, 3, '0'),
            "  https://insp.cd/", slug, "/")
    end
    println("\nFetch, scan, and extend data/observations.toml + ",
        "data/insp_sitrep_scanned.csv per data/README.md.")
    return 1
end

exit(main())
