#!/usr/bin/env julia
#
# Report how many INSP situation reports we have NOT yet recorded.
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
#    (so a caller / cron can branch on it).

using Downloads

const UA = "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) " *
           "AppleWebKit/605.1.15"
const LIST_URL = "https://insp.cd/wp-json/wp/v2/posts?" *
                 "search=sitrep&per_page=100&_fields=id,slug,date"
const SCANNED = joinpath(@__DIR__, "..", "data", "insp_sitrep_scanned.csv")

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
