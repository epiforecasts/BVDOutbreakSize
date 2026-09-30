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
# dropped individual vintages, so it is not the primary source. Its copies
# of the PDFs are byte-identical to INSP's, so it serves as the fallback
# below, and its processed national CSVs remain a cross-check
# (scripts/confirm_insp_data.jl).
#
# The media listing has dropped reports and on some days names no MVE PDF
# at all. So after it the script lists the published MVE posts and
# fills every report still missing on disk, rather than stopping, so a
# fresh checkout can rebuild the whole cache. It takes each from the mirror
# when the mirror holds it, which costs insp.cd nothing, and the rest
# through the per-post path below, pausing between requests.
#
# Usage:
#
#   julia --project=scripts scripts/download_sitreps.jl
#   julia --project=scripts scripts/download_sitreps.jl path/to/outdir
#   julia --project=scripts scripts/download_sitreps.jl --only 074[,073,...] [path/to/outdir]
#   julia --project=scripts scripts/download_sitreps.jl --limit 3 [path/to/outdir]
#
# With no argument the PDFs land in `data/sitrep_pdfs/` (git-ignored).
#
# `--only` fetches just the listed report numbers via a direct per-post
# lookup (the manual method documented in data/README.md, scripted here)
# instead of walking the paginated media-listing API. It costs two requests
# per report (find the post, decode its embedded PDF URL) rather than the
# ~50-request full-archive walk, so it is the considerate option when only
# a single new report is wanted, or when the media API is struggling but
# the posts API (used by check_new_sitreps.jl) still answers.
#
# `--limit N` fetches at most N missing reports, newest first, in either
# mode.
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
using Base64

const MEDIA_API = "https://insp.cd/wp-json/wp/v2/media"
const POSTS_API = "https://insp.cd/wp-json/wp/v2/posts"
const UA = "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) " *
    "AppleWebKit/605.1.15"
const REQUEST_TIMEOUT = 180.0
# insp.cd runs on modest infrastructure mid-outbreak: give it more time to
# answer and more space between retries before giving up, rather than
# treating a slow response as a reason to come back sooner.
const ATTEMPTS = 5
const BACKOFF_SECONDS = 5
# The walk stops itself: WordPress marks the end of the listing, but an API
# or cache that keeps answering page 1 would otherwise spin forever, and a
# nightly hang is worse than an error.
const MAX_PAGES = 50
# Pause between consecutive insp.cd requests on the fallback path, which can
# make a few dozen of them on a fresh checkout.
const PAUSE_SECONDS = 5
const MIRROR_TREE = "https://api.github.com/repos/INRB-UMIE/BDBV2026-Data/" *
    "git/trees/main?recursive=1"
const MIRROR_RAW = "https://raw.githubusercontent.com/INRB-UMIE/" *
    "BDBV2026-Data/main/"
# The mirror keeps its raw PDFs in Git LFS, so the raw URL serves a pointer
# for most of them and the media URL serves the file.
const MIRROR_LFS = "https://media.githubusercontent.com/media/INRB-UMIE/" *
    "BDBV2026-Data/main/"

# Parse `--only N1,N2,...` and `--limit N` out of ARGS, leaving any outdir
# argument in place.
function parse_args(args)
    rest = copy(args)
    function pop_flag!(flag)
        idx = findfirst(==(flag), rest)
        idx === nothing && return nothing
        idx == length(rest) && error("$flag requires a value")
        value = rest[idx + 1]
        deleteat!(rest, idx:(idx + 1))
        return value
    end
    only = pop_flag!("--only")
    numbers = only === nothing ? nothing :
        [lpad(strip(s), 3, '0') for s in split(only, ",")]
    limit = pop_flag!("--limit")
    return (;
        only_numbers = numbers,
        limit = limit === nothing ? typemax(Int) : parse(Int, limit), rest,
    )
end

parsed = parse_args(ARGS)
only_numbers = parsed.only_numbers
limit = parsed.limit
outdir = length(parsed.rest) >= 1 ? parsed.rest[1] :
    joinpath(@__DIR__, "..", "data", "sitrep_pdfs")
mkpath(outdir)

# Decode the JSON string escapes WordPress returns in source_url (`\/` and
# `\uXXXX`, e.g. `°` for the degree sign in `N°60`). The replacement
# receives the whole match, so the four hex digits start at index 3, and
# they are hexadecimal: `parse` defaults to base 10 and rejects `00b0`.
function json_unescape(s)
    s = replace(s, "\\/" => "/")
    return replace(
        s,
        r"\\u([0-9a-fA-F]{4})" => m -> string(
            Char(
                parse(
                    UInt16, m[3:6];
                    base = 16
                )
            )
        )
    )
end

# Pull the SitRep number out of a filename, tolerating the inconsistent
# `N°60`, `N_61`, `No59`, `N°055` spellings, and zero-pad it to three digits.
# From SitRep 084 the "MVEBDB brief format" renamed files to
# `SitRep_MVEBDB_NNN_...` with no `N`/`No`/`N°` before the digits at all, so
# the `n[...]` alternative alone silently matched nothing for 084 onward -
# not even landing in the "rejected" list, since `sitrep_number` returning
# `nothing` skips the file before the MVE-vs-other-disease check runs. The
# `mvebdb[...]` alternative below recovers that convention too.
function sitrep_number(name)
    m = match(r"(?i)sitrep.*?(?:n[°o._\- ]*|mvebdb[_\- ]*)0*(\d{2,3})", name)
    return m === nothing ? nothing : lpad(m.captures[1], 3, '0')
end

# One media-API page. Returns the body with the HTTP status and, when the
# request never got a reply, the last exception, so the caller can tell the
# end of the listing from a refusal or a timeout and can say which it was.
# A 400 body is returned too: WordPress marks the real end of the listing
# with `rest_post_invalid_page_number`, and a 400 from anything else (a rate
# limiter, say) must not be mistaken for it.
function api_page(url)
    last_status = 0
    last_err = nothing
    for attempt in 1:ATTEMPTS
        io = IOBuffer()
        res = try
            Downloads.request(
                url; output = io, throw = false,
                headers = ["User-Agent" => UA], timeout = REQUEST_TIMEOUT
            )
        catch err
            last_err = err
            nothing
        end
        if res isa Downloads.Response
            body = String(take!(io))
            (res.status == 200 || res.status == 400) &&
                return (; status = res.status, body, err = nothing)
            last_status = res.status
        end
        attempt < ATTEMPTS && sleep(BACKOFF_SECONDS * 2^attempt)
    end
    return (; status = last_status, body = "", err = last_err)
end

# What went wrong, in the words of whatever went wrong, rather than a list
# of guesses for the next person to work through.
function page_failure(res)
    res.status != 0 && return "HTTP $(res.status)"
    res.err === nothing && return "no response"
    return string(res.err)
end

# Walk the paginated media API and collect (number => source_url) for every
# MVE SitRep PDF, keeping the first URL seen per number. Also returns the
# SitRep-numbered files rejected for not naming MVE, so a renamed MVE report
# is visible rather than just missing.
function collect_sitrep_urls()
    urls = Dict{String, String}()
    rejected = String[]
    for page in 1:MAX_PAGES
        res = api_page("$MEDIA_API?search=sitrep&per_page=100&page=$page")
        ## Past the last page WordPress 400s with this code. Any other 400
        ## is a real failure and must not silently truncate the listing.
        if res.status == 400
            occursin("rest_post_invalid_page_number", res.body) && break
            error(
                "media API returned HTTP 400 at $MEDIA_API (page " *
                    "$page) without the end-of-listing code: $(res.body)"
            )
        end
        res.status == 200 || error(
            "media API request failed after $ATTEMPTS attempts at " *
                "$MEDIA_API (page $page): $(page_failure(res))"
        )
        hits = collect(
            eachmatch(
                r"\"source_url\":\"([^\"]*?\.pdf)\"",
                res.body
            )
        )
        isempty(hits) && break
        for h in hits
            url = json_unescape(h.captures[1])
            name = basename(url)
            occursin(r"(?i)sitrep", name) || continue
            num = sitrep_number(name)
            num === nothing && continue
            if !occursin(r"(?i)mve", name)
                push!(rejected, name)
                continue
            end
            get!(urls, num, url)
        end
        page == MAX_PAGES && error(
            "media API still returning results after $MAX_PAGES pages at " *
                "$MEDIA_API: pagination is not advancing?"
        )
    end
    return (; urls, rejected = sort!(unique(rejected)))
end

function fetch_pdf(url, dest)
    last_err = nothing
    for attempt in 1:ATTEMPTS
        try
            Downloads.download(
                url, dest; headers = ["User-Agent" => UA],
                timeout = REQUEST_TIMEOUT
            )
            return (; ok = true, err = nothing)
        catch err
            last_err = err
            attempt < ATTEMPTS && sleep(BACKOFF_SECONDS * 2^attempt)
        end
    end
    return (; ok = false, err = last_err)
end

# Published (number, post id, slug) triples straight from the posts API -
# the same endpoint and query check_new_sitreps.jl uses, so a report that
# is visible there is findable here even when the media-listing API isn't
# cooperating. The posts come newest first, 100 to a page, so the walk
# stops as soon as every number in `want` has been seen; with no `want` it
# reads to the end of the listing, which WordPress marks as the media walk
# does.
function published_posts(; want = nothing)
    out = Tuple{String, Int, String}[]
    for page in 1:MAX_PAGES
        page > 1 && sleep(PAUSE_SECONDS)
        res = api_page(
            "$POSTS_API?search=sitrep&per_page=100&_fields=id,slug,date" *
                "&page=$page"
        )
        res.status == 400 &&
            occursin("rest_post_invalid_page_number", res.body) && break
        res.status == 200 || error(
            "posts API request failed after $ATTEMPTS attempts at " *
                "$POSTS_API (page $page): $(page_failure(res))"
        )
        hits = collect(
            eachmatch(
                r"\{\"id\":(\d+),\"date\":\"[^\"]*\",\"slug\":\"(sitrep[^\"]*)\"\}",
                res.body
            )
        )
        isempty(hits) && break
        for m in hits
            slug = m.captures[2]
            num = match(r"-n0*(\d+)", slug)
            num === nothing && continue
            push!(
                out, (
                    lpad(num.captures[1], 3, '0'), parse(Int, m.captures[1]),
                    slug,
                )
            )
        end
        want !== nothing && want ⊆ Set(p[1] for p in out) && break
    end
    return out
end

# The MVE SitRep PDFs the INRB-UMIE mirror holds, as number => path, from
# one GitHub tree listing. Re-issues (`_v2`) are left out, as they are on
# the INSP path.
function mirror_urls()
    res = api_page(MIRROR_TREE)
    res.status == 200 || error(
        "mirror tree listing failed after $ATTEMPTS attempts at " *
            "$MIRROR_TREE: $(page_failure(res))"
    )
    urls = Dict{String, String}()
    for m in eachmatch(
            r"\"path\":\"(data/insp_sitrep/raw/SitRep_MVE_0*(\d+)[-_]2026\.pdf)\"",
            res.body
        )
        urls[lpad(m.captures[2], 3, '0')] = m.captures[1]
    end
    return urls
end

# Decode the `pdfemb-data` base64 blob embedded in a rendered post (the
# same mechanism data/README.md's manual fetch recipe documents) to recover
# the direct, fetchable PDF URL.
function embedded_pdf_url(content)
    m = match(r"pdfemb-data=([A-Za-z0-9_-]+)", content)
    m === nothing && return nothing
    b64 = replace(m.captures[1], '-' => '+', '_' => '/')
    b64 *= "="^mod(-length(b64), 4)
    decoded = String(base64decode(b64))
    um = match(r"\"url\":\"([^\"]*)\"", decoded)
    um === nothing && return nothing
    return json_unescape(um.captures[1])
end

function fetch_post_pdf_url(id)
    res = api_page("$POSTS_API/$id?_fields=content")
    res.status == 200 || error(
        "post fetch failed after $ATTEMPTS attempts at $POSTS_API/" *
            "$id: $(page_failure(res))"
    )
    return embedded_pdf_url(res.body)
end

mve_posts(posts) = Dict(
    p[1] => p for p in reverse(posts) if occursin(r"(?i)mve", p[3])
)

# Download `url` to `dest`, reporting the size or the failure. Returns
# whether the file landed.
function fetch_to(url, dest, label)
    print("fetch  $label ... ")
    res = fetch_pdf(url, dest)
    if res.ok
        println("$(round(filesize(dest) / 1024; digits = 1)) KiB")
        return true
    end
    isfile(dest) && rm(dest)
    println("FAILED after $ATTEMPTS attempts ($(res.err))")
    return false
end

# Download a mirror file, following its LFS pointer when the raw URL
# serves one, and keep it only when it is a PDF.
function fetch_mirror(path, dest, label)
    fetch_to(MIRROR_RAW * path, dest, "$label (mirror)") || return false
    if startswith(read(dest, String), "version https://git-lfs")
        fetch_to(MIRROR_LFS * path, dest, "$label (mirror LFS)") ||
            return false
    end
    startswith(read(dest, String), "%PDF") && return true
    rm(dest)
    println("  not a PDF, discarded")
    return false
end

# Look report `num` up through its post and download the PDF its content
# embeds.
function fetch_from_post(post, dest)
    num, id, slug = post
    print("lookup SitRep $num ($slug) ... ")
    pdf_url = fetch_post_pdf_url(id)
    if pdf_url === nothing
        println("no embedded PDF URL found in post content")
        return false
    end
    println("found")
    return fetch_to(pdf_url, dest, "SitRep $num")
end

dest_for(num) = joinpath(outdir, "SitRep_MVE_$(num)_2026.pdf")

# The newest `limit` of `numbers` not already on disk, in ascending order.
function missing_numbers(numbers, limit)
    todo = sort([n for n in numbers if !isfile(dest_for(n))]; rev = true)
    return sort(first(todo, min(limit, length(todo))))
end

# Selective mode: look each number up directly via the posts API and
# download just its PDF, never touching the paginated media listing.
function fetch_selected(numbers, limit)
    for n in numbers
        isfile(dest_for(n)) && println("skip   SitRep $n (already present)")
    end
    todo = missing_numbers(numbers, limit)
    posts = isempty(todo) ? Dict() :
        mve_posts(published_posts(; want = Set(todo)))
    downloaded = 0
    for num in todo
        if !haskey(posts, num)
            println(
                "SKIP   SitRep $num: no MVE post found among " *
                    "published sitreps"
            )
            continue
        end
        fetch_from_post(posts[num], dest_for(num)) && (downloaded += 1)
    end
    return println(
        "\n$downloaded new sitrep(s) into $outdir (selective mode)."
    )
end

# Fill the reports the media listing did not give: all of them when it
# names no MVE PDF at all, and the ones it has dropped otherwise.
# Every published MVE report still missing on disk is taken from the mirror
# when the mirror holds it, which costs insp.cd nothing, and through its
# post otherwise, pausing between insp.cd requests.
function fetch_gaps(limit)
    posts = mve_posts(published_posts())
    todo = missing_numbers(keys(posts), limit)
    isempty(todo) && return nothing
    println(
        "\n$(length(todo)) published report(s) missing after the media " *
            "listing: taking them from the INRB-UMIE mirror and the posts API"
    )
    mirror = mirror_urls()
    downloaded = 0
    from_posts = 0
    for num in todo
        if haskey(mirror, num)
            ok = fetch_mirror(mirror[num], dest_for(num), "SitRep $num")
        else
            from_posts > 0 && sleep(PAUSE_SECONDS)
            ok = fetch_from_post(posts[num], dest_for(num))
            from_posts += 1
        end
        ok && (downloaded += 1)
    end
    return println(
        "$downloaded of $(length(todo)) filled, $from_posts through the " *
            "posts API."
    )
end

const count_new = Ref(0)

function fetch_listed(listing, limit)
    ## Measles and SGI-GPM SitReps share this media library and their
    ## numbering collides with the MVE series, so these rejections are
    ## expected. They are printed anyway: if an MVE report is ever
    ## published without MVE in the filename it lands here, and a gap in
    ## the series would otherwise be indistinguishable from a rename.
    if !isempty(listing.rejected)
        println(
            "skipped $(length(listing.rejected)) sitrep-numbered " *
                "non-MVE file(s):"
        )
        for name in listing.rejected
            println("  $name")
        end
        println()
    end
    urls = listing.urls
    downloaded = 0
    for num in missing_numbers(keys(urls), limit)
        fetch_to(urls[num], dest_for(num), "SitRep $num") && (downloaded += 1)
    end
    count_new[] = downloaded
    return println(
        "\n$downloaded new sitrep(s) into $outdir ($(length(urls)) " *
            "upstream)."
    )
end

if only_numbers !== nothing
    fetch_selected(only_numbers, limit)
else
    listing = collect_sitrep_urls()
    isempty(listing.urls) && println(
        "no MVE SitRep PDFs in the media listing at $MEDIA_API"
    )
    fetch_listed(listing, limit)
    fetch_gaps(max(0, limit - count_new[]))
end
