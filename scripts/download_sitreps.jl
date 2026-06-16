#!/usr/bin/env julia
#
# Download the INSP situation-report PDFs from the INRB-UMIE mirror
# (https://github.com/INRB-UMIE/BDBV2026-Data, data/insp_sitrep/raw) so
# they can be scanned directly. The raw PDFs are stored with Git LFS; the
# `github.com/.../raw/...` URL serves the smudged binary, so no LFS client
# is needed. Any report already present in the output directory is
# skipped, so re-running fetches only the new sitreps.
#
# Usage:
#
#   julia --project=scripts scripts/download_sitreps.jl
#   julia --project=scripts scripts/download_sitreps.jl path/to/outdir
#
# With no argument the PDFs land in `data/sitrep_pdfs/` (git-ignored).
# Set GITHUB_TOKEN in the environment to lift the unauthenticated API
# rate limit if you fetch often.

using Downloads

const REPO = "INRB-UMIE/BDBV2026-Data"
const RAW_DIR = "data/insp_sitrep/raw"
const CONTENTS_API = "https://api.github.com/repos/$REPO/contents/$RAW_DIR"
const RAW_BASE = "https://github.com/$REPO/raw/main/$RAW_DIR"

outdir = length(ARGS) >= 1 ? ARGS[1] :
         joinpath(@__DIR__, "..", "data", "sitrep_pdfs")
mkpath(outdir)

api_headers = haskey(ENV, "GITHUB_TOKEN") ?
              ["Authorization" => "Bearer $(ENV["GITHUB_TOKEN"])"] :
              Pair{String, String}[]

# List the raw directory through the contents API and pull the PDF file
# names out of the JSON. A regex avoids a JSON dependency; the contents
# API returns one `"name": "..."` field per entry.
listing = sprint() do io
    Downloads.download(CONTENTS_API, io; headers = api_headers)
end
names = sort(unique(
    m.captures[1] for m in eachmatch(r"\"name\":\s*\"(SitRep[^\"]+\.pdf)\"",
    listing)))

isempty(names) &&
    error("no SitRep PDFs found at $CONTENTS_API (rate limited?)")

downloaded = 0
for name in names
    dest = joinpath(outdir, name)
    if isfile(dest)
        println("skip  $name (already present)")
        continue
    end
    print("fetch $name ... ")
    Downloads.download("$RAW_BASE/$name", dest)
    println("$(round(filesize(dest) / 1024; digits = 1)) KiB")
    global downloaded += 1
end

println("\n$downloaded new sitrep(s) into $outdir ($(length(names)) upstream).")
