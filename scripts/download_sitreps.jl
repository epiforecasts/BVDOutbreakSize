#!/usr/bin/env julia
#
# Download the INSP situation-report PDFs from the official INSP site
# (https://insp.cd), the primary source of truth, so they can be scanned
# directly. The PDFs are WordPress media uploads; this script queries the
# WordPress REST media API, picks out the SitRep PDFs, normalises each to
# SitRep_MVE_NNN_2026.pdf and downloads any not already present.
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

using Downloads

const MEDIA_API = "https://insp.cd/wp-json/wp/v2/media"

outdir = length(ARGS) >= 1 ? ARGS[1] :
         joinpath(@__DIR__, "..", "data", "sitrep_pdfs")
mkpath(outdir)

# Decode the JSON string escapes WordPress returns in source_url (`\/` and
# `\uXXXX`, e.g. `°` for the degree sign in `N°60`).
function json_unescape(s)
    s = replace(s, "\\/" => "/")
    replace(s, r"\\u([0-9a-fA-F]{4})" => m -> string(Char(parse(UInt16, m[3:6]))))
end

# Pull the SitRep number out of a filename, tolerating the inconsistent
# `N°60`, `N_61`, `No59`, `N°055` spellings, and zero-pad it to three digits.
function sitrep_number(name)
    m = match(r"(?i)sitrep.*?n[°o._\- ]*0*(\d{2,3})", name)
    m === nothing ? nothing : lpad(m.captures[1], 3, '0')
end

# Walk the paginated media API and collect (number => source_url) for every
# SitRep PDF, keeping the first URL seen per number.
urls = Dict{String, String}()
page = 1
while true
    body = try
        sprint() do io
            Downloads.download("$MEDIA_API?search=sitrep&per_page=100&page=$page",
                io)
        end
    catch
        break  # past the last page the API returns a 400
    end
    hits = collect(eachmatch(r"\"source_url\":\"([^\"]*?\.pdf)\"", body))
    isempty(hits) && break
    for h in hits
        url = json_unescape(h.captures[1])
        name = basename(url)
        occursin(r"(?i)sitrep", name) || continue
        num = sitrep_number(name)
        num === nothing && continue
        get!(urls, num, url)
    end
    page += 1
end

isempty(urls) &&
    error("no SitRep PDFs found at $MEDIA_API (site down or API changed?)")

downloaded = 0
for num in sort(collect(keys(urls)))
    dest = joinpath(outdir, "SitRep_MVE_$(num)_2026.pdf")
    if isfile(dest)
        println("skip  SitRep $num (already present)")
        continue
    end
    print("fetch SitRep $num ... ")
    try
        Downloads.download(urls[num], dest)
        println("$(round(filesize(dest) / 1024; digits = 1)) KiB")
        global downloaded += 1
    catch err
        println("FAILED ($err)")
    end
end

println("\n$downloaded new sitrep(s) into $outdir ($(length(urls)) upstream).")
