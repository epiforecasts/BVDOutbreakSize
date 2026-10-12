#!/usr/bin/env julia
#
# Collect the dashboard map's inputs from past results releases, one
# vintage per data cut-off, so the map can show the estimates as they stood
# on an earlier date.
#
# Usage:
#
#   julia --project=docs scripts/archive_dashboard.jl
#   julia --project=docs scripts/archive_dashboard.jl DEST [owner/repo]
#
# DEST defaults to `logs/dashboard_archive`, which `stage_zone_map` in
# `docs/make.jl` copies to `zone_map/archive/` when it is present. The
# archive is written as
#
#   DEST/index.json               [{"as_of": "2026-10-10", "tag": "results-3000"}, ...]
#   DEST/<as_of>/<name>.csv       the five dashboard tables of that release
#
# with `index.json` listing the vintages newest first. A release is a
# candidate only when it carries all five tables (`DASHBOARD_CSVS`) and its
# `observations.toml`; releases from before the tables were published are
# skipped. Among the candidates, `select_daily_releases` keeps one per
# cut-off (`as_of_date`), as `scripts/score_releases.jl` does, and the
# newest `MAX_VINTAGES` cut-offs are kept.
#
# Requires the `gh` CLI to be authenticated against the repo. Releases are
# listed through the REST API, which also gives each one's asset names, so
# no table is downloaded for a release that would be skipped.

using BVDOutbreakSize: is_results_release, select_daily_releases
using Dates
using TOML

const DEFAULT_REPO = "epiforecasts/BVDOutbreakSize"
const DEFAULT_DEST = joinpath(dirname(@__DIR__), "logs", "dashboard_archive")

## The tables the dashboard map fetches, all but the health-zone geojson,
## which does not change between releases.
const DASHBOARD_CSVS = (
    "zone_estimates.csv", "province_estimates.csv",
    "national_estimates.csv", "zone_timeseries.csv",
    "province_timeseries.csv",
)

## Four weeks of daily cut-offs, at about 370 KB each.
const MAX_VINTAGES = 28

## Every published results release as `(tag, created, assets)`, drafts
## left out. Pages are requested by number rather than through
## `--paginate`, whose follow-up links some proxies refuse.
function release_entries(repo; per_page = 100)
    jq = ".[] | [.tag_name, .created_at, .draft,
        ([.assets[].name] | join(\",\"))] | @tsv"
    entries = Tuple{String, DateTime, Vector{String}}[]
    for page in Iterators.countfrom(1)
        url = "repos/$repo/releases?per_page=$per_page&page=$page"
        lines = filter(!isempty, split(read(`gh api $url --jq $jq`, String), '\n'))
        for line in lines
            tag, created, draft, assets = split(line, '\t')
            draft == "true" && continue
            is_results_release(tag) || continue
            push!(
                entries,
                (
                    String(tag),
                    DateTime(created, dateformat"yyyy-mm-ddTHH:MM:SSZ"),
                    String.(split(assets, ',')),
                )
            )
        end
        length(lines) < per_page && break
    end
    return entries
end

## Whether a release's assets include everything a vintage needs.
carries_dashboard(assets) =
    "observations.toml" in assets && all(in(assets), DASHBOARD_CSVS)

## Download a single asset from a release into `dir`, returning its path
## (or `nothing` when the download keeps failing). Mirrors
## `scripts/refresh_releases.jl`.
function fetch_asset(repo, tag, file, dir; attempts = 3)
    dest = joinpath(dir, file)
    for _ in 1:attempts
        try
            run(
                pipeline(
                    `gh release download $tag -R $repo -p $file
       -O $dest --clobber`; stdout = devnull, stderr = devnull
                )
            )
        catch
            continue
        end
        isfile(dest) && return dest
    end
    return nothing
end

## The data cut-off a release was built from, or `nothing`.
function release_cutoff(repo, tag, dir)
    obs = fetch_asset(repo, tag, "observations.toml", dir)
    isnothing(obs) && return nothing
    parsed = TOML.parsefile(obs)
    haskey(parsed, "as_of_date") || return nothing
    return Date(string(parsed["as_of_date"]))
end

"""
    pick_vintages(entries; limit = MAX_VINTAGES)

One release per data cut-off from `entries`, `(tag, created, cutoff)`
triples, as `(as_of, tag)` pairs newest first and at most `limit` long.
Which release a cut-off keeps is decided by `select_daily_releases`.
"""
function pick_vintages(entries; limit = MAX_VINTAGES)
    cutoff = Dict(String(t) => Date(d) for (t, _, d) in entries)
    tags = select_daily_releases(entries)
    return [(cutoff[t], t) for t in first(tags, limit)]
end

## `index.json` for `vintages`, `(as_of, tag)` pairs in the order given.
## Tags and dates carry no characters JSON would need escaped.
function index_json(vintages)
    rows = ["  {\"as_of\": \"$d\", \"tag\": \"$t\"}" for (d, t) in vintages]
    return isempty(rows) ? "[]\n" : "[\n" * join(rows, ",\n") * "\n]\n"
end

"""
    write_archive(dest, vintages, source)

Write the archive for `vintages`, `(as_of, tag)` pairs newest first, into a
fresh `dest`, copying each release's tables from `source/<tag>/`. A vintage
missing any of `DASHBOARD_CSVS` there is left out with a warning. Returns
the vintages written, which `dest/index.json` lists in the same order.
"""
function write_archive(dest, vintages, source)
    rm(dest; force = true, recursive = true)
    mkpath(dest)
    written = empty(vintages)
    for (d, t) in vintages
        files = [joinpath(source, t, f) for f in DASHBOARD_CSVS]
        if !all(isfile, files)
            @warn "skipping $t (dashboard tables missing)" as_of = d
            continue
        end
        out = mkpath(joinpath(dest, string(d)))
        for f in files
            cp(f, joinpath(out, basename(f)))
        end
        push!(written, (d, t))
    end
    write(joinpath(dest, "index.json"), index_json(written))
    return written
end

function main(args)
    dest = length(args) >= 1 ? args[1] : DEFAULT_DEST
    repo = length(args) >= 2 ? args[2] : DEFAULT_REPO
    all_releases = release_entries(repo)
    candidates = filter(e -> carries_dashboard(e[3]), all_releases)
    @info "results releases carrying the dashboard tables" carrying =
        length(candidates) total = length(all_releases)
    return mktempdir() do dir
        entries = Tuple{String, DateTime, Date}[]
        for (tag, created, _) in candidates
            d = release_cutoff(repo, tag, mkpath(joinpath(dir, tag)))
            if isnothing(d)
                @warn "skipping $tag (no parseable as_of_date)"
                continue
            end
            push!(entries, (tag, created, d))
        end
        vintages = pick_vintages(entries)
        for (_, tag) in vintages, f in DASHBOARD_CSVS
            fetch_asset(repo, tag, f, joinpath(dir, tag))
        end
        written = write_archive(dest, vintages, dir)
        println("Wrote $(length(written)) dashboard vintages to $dest")
        for (d, t) in written
            println("  $d  $t")
        end
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end
