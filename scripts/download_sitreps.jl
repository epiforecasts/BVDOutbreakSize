#!/usr/bin/env julia
#
# Download the INSP situation-report PDFs from the official INSP site
# (https://insp.cd), the primary source of truth, so they can be scanned
# directly. The PDFs are WordPress media uploads; this script queries the
# WordPress REST media API, picks out the MVE SitRep PDFs, normalises each
# to SitRep_MVE_NNN_2026.pdf and downloads any not already present.
#
# The INSP site leads the INRB-UMIE GitHub mirror
# (https://github.com/INRB-UMIE/BDBV2026-Data), which lags by days and has
# dropped individual vintages, so it is no longer used as the source. The
# mirror's processed national CSVs remain a cross-check only
# (scripts/confirm_insp_data.jl).
#
# Usage:
#
#   julia --project=scripts scripts/download_sitreps.jl
#   julia --project=scripts scripts/download_sitreps.jl path/to/outdir
#
# With no argument the PDFs land in `data/sitrep_pdfs/` (git-ignored).
#
# Notes:
#  - INSP blocks default user agents with an HTTP 403, so every request
#    sends a browser User-Agent (the same one scripts/check_new_sitreps.jl
#    uses). Without it the media API returns nothing and the script cannot
#    tell an empty upstream from a refused one.
#  - The site answers slowly (15-20 s per call is normal), hence the
#    generous timeout and the retries.
#  - INSP publishes SitReps for several diseases through the same media
#    library, and their numbering collides with the MVE series
#    (`Sitrep-SGI-Rougeole-Rubeole-N°-29.pdf`, `VF_SITREP_SGI-GPM_N°10...`),
#    so the filename must name MVE to be mirrored.
#  - A corrected re-issue (`..._v2.pdf`) is invisible here: one URL is kept
#    per SitRep number. data/insp_sitrep_scanned.csv already carries a
#    `012_v2` row, so a re-issue has to be fetched by hand.

using Downloads

const MEDIA_API = "https://insp.cd/wp-json/wp/v2/media"
const UA = "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) " *
           "AppleWebKit/605.1.15"
const REQUEST_TIMEOUT = 120.0
const ATTEMPTS = 3

outdir = length(ARGS) >= 1 ? ARGS[1] :
         joinpath(@__DIR__, "..", "data", "sitrep_pdfs")
mkpath(outdir)

# Decode the JSON string escapes WordPress returns in source_url (`\/` and
# `\uXXXX`, e.g. `°` for the degree sign in `N°60`). The replacement
# receives the whole match, so the four hex digits start at index 3, and
# they are hexadecimal: `parse` defaults to base 10 and rejects `00b0`.
function json_unescape(s)
    s = replace(s, "\\/" => "/")
    replace(s,
        r"\\u([0-9a-fA-F]{4})" => m -> string(Char(parse(UInt16, m[3:6];
            base = 16))))
end

# Pull the SitRep number out of a filename, tolerating the inconsistent
# `N°60`, `N_61`, `No59`, `N°055` spellings, and zero-pad it to three digits.
function sitrep_number(name)
    m = match(r"(?i)sitrep.*?n[°o._\- ]*0*(\d{2,3})", name)
    m === nothing ? nothing : lpad(m.captures[1], 3, '0')
end

# One media-API page. Returns the body with its HTTP status so the caller
# can tell the end of the listing (the API 400s past the last page) from a
# refusal or a timeout, which are retried and then reported.
function api_page(url)
    for attempt in 1:ATTEMPTS
        io = IOBuffer()
        res = try
            Downloads.request(url; output = io, throw = false,
                headers = ["User-Agent" => UA], timeout = REQUEST_TIMEOUT)
        catch
            nothing
        end
        if res isa Downloads.Response
            res.status == 200 && return (; status = 200,
                body = String(take!(io)))
            res.status == 400 && return (; status = 400, body = "")
        end
        attempt < ATTEMPTS && sleep(2^attempt)
    end
    return (; status = 0, body = "")
end

# Walk the paginated media API and collect (number => source_url) for every
# MVE SitRep PDF, keeping the first URL seen per number.
function collect_sitrep_urls()
    urls = Dict{String, String}()
    page = 1
    while true
        res = api_page("$MEDIA_API?search=sitrep&per_page=100&page=$page")
        res.status == 400 && break  # past the last page
        res.status == 200 || error(
            "media API request failed after $ATTEMPTS attempts at " *
            "$MEDIA_API (page $page): site down, refusing this user " *
            "agent, or the API changed?")
        hits = collect(eachmatch(r"\"source_url\":\"([^\"]*?\.pdf)\"",
            res.body))
        isempty(hits) && break
        for h in hits
            url = json_unescape(h.captures[1])
            name = basename(url)
            occursin(r"(?i)sitrep", name) || continue
            occursin(r"(?i)mve", name) || continue
            num = sitrep_number(name)
            num === nothing && continue
            get!(urls, num, url)
        end
        page += 1
    end
    return urls
end

function fetch_pdf(url, dest)
    for attempt in 1:ATTEMPTS
        try
            Downloads.download(url, dest; headers = ["User-Agent" => UA],
                timeout = REQUEST_TIMEOUT)
            return true
        catch err
            attempt == ATTEMPTS && return false
            sleep(2^attempt)
        end
    end
    return false
end

urls = collect_sitrep_urls()
isempty(urls) &&
    error("no MVE SitRep PDFs found at $MEDIA_API (site down or API " *
          "changed?)")

downloaded = 0
for num in sort(collect(keys(urls)))
    dest = joinpath(outdir, "SitRep_MVE_$(num)_2026.pdf")
    if isfile(dest)
        println("skip  SitRep $num (already present)")
        continue
    end
    print("fetch SitRep $num ... ")
    if fetch_pdf(urls[num], dest)
        println("$(round(filesize(dest) / 1024; digits = 1)) KiB")
        global downloaded += 1
    else
        isfile(dest) && rm(dest)
        println("FAILED after $ATTEMPTS attempts")
    end
end

println("\n$downloaded new sitrep(s) into $outdir ($(length(urls)) upstream).")
