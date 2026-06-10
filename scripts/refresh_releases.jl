#!/usr/bin/env julia
#
# Refresh the released outbreak-size estimates plotted in the
# "Estimate evolution across releases" figure by pulling each tagged
# results release straight from GitHub
# (https://github.com/epiforecasts/BVDOutbreakSize/releases) and
# writing them to `data/released_estimates.csv`.
#
# Usage:
#
#   julia --project=scripts scripts/refresh_releases.jl
#   julia --project=scripts scripts/refresh_releases.jl owner/repo
#
# Each `results-v*` release bundles the published summary tables and a
# copy of the input data. For every such release this reads:
#   * `observations.toml`            -> the data cut-off (`as_of_date`)
#   * `cumulative_cases_by_stream.csv` -> the published `joint` headline
#     outbreak size with its 30/60/90% bounds
#   * `posterior_draws.csv`          -> the posterior median
#
# `joint` is the headline outbreak size every release publishes: the
# integral model's `cumulative_cases` for the pre-renewal vintages and
# the renewal model's `C_T` from the first renewal release on. The model
# is detected from the draws columns (`C_T` for renewal, otherwise the
# integral `cumulative_cases`) and recorded so the figure can tell the
# two eras apart: once a release is renewal-model, its published estimate
# already is the renewal fit.
#
# Requires the `gh` CLI to be authenticated against the repo.

using CSV
using DataFrames
using Dates
using Statistics
using TOML

const DEFAULT_REPO = "epiforecasts/BVDOutbreakSize"

## The headline outbreak-size row in `cumulative_cases_by_stream.csv`,
## and the draws column carrying its posterior, per model.
const STREAM_ROW = "joint"
const RENEWAL_DRAW = "C_T"
const INTEGRAL_DRAW = "cumulative_cases"

repo = length(ARGS) >= 1 ? ARGS[1] : DEFAULT_REPO

## All tagged results releases, newest first, as reported by `gh`.
function results_release_tags(repo)
    out = read(`gh release list -R $repo -L 200 --json tagName
        --jq ".[].tagName"`, String)
    tags = filter(t -> occursin(r"^results-v[0-9]", t),
        split(strip(out), '\n'))
    return tags
end

## Download a single asset from a release into `dir`, returning its path
## (or `nothing` when the release does not carry that asset). Retries a
## couple of times so a transient `gh` failure does not drop a release.
function fetch_asset(repo, tag, file, dir; attempts = 3)
    dest = joinpath(dir, file)
    for _ in 1:attempts
        try
            run(pipeline(
                `gh release download $tag -R $repo -p $file
       -O $dest --clobber`; stdout = devnull, stderr = devnull))
        catch
            continue
        end
        isfile(dest) && return dest
    end
    return nothing
end

## Pull one release into a `(tag, date, model, median, lo30, hi30, lo60,
## hi60, lo90, hi90)` row, or `nothing` when it cannot be parsed.
function release_row(repo, tag, dir)
    obs = fetch_asset(repo, tag, "observations.toml", dir)
    streams = fetch_asset(repo, tag, "cumulative_cases_by_stream.csv", dir)
    draws_file = fetch_asset(repo, tag, "posterior_draws.csv", dir)
    any(isnothing, (obs, streams, draws_file)) && return nothing

    parsed = TOML.parsefile(obs)
    haskey(parsed, "as_of_date") || return nothing
    date = string(parsed["as_of_date"])

    df = CSV.read(streams, DataFrame)
    cols = names(df)              # [Stream, Lower 90%, ..., Upper 90%]
    row = findfirst(==(STREAM_ROW), df[!, 1])
    isnothing(row) && return nothing
    lo30, hi30 = df[row, cols[4]], df[row, cols[5]]
    lo60, hi60 = df[row, cols[3]], df[row, cols[6]]
    lo90, hi90 = df[row, cols[2]], df[row, cols[7]]

    ## The model, and the draws column carrying the headline posterior.
    draws = CSV.read(draws_file, DataFrame)
    model, drawcol = if RENEWAL_DRAW in names(draws)
        ("renewal", RENEWAL_DRAW)
    elseif INTEGRAL_DRAW in names(draws)
        ("integral", INTEGRAL_DRAW)
    else
        return nothing
    end
    med = median(draws[!, drawcol])

    return (tag, date, model, med, lo30, hi30, lo60, hi60, lo90, hi90)
end

tags = results_release_tags(repo)
isempty(tags) && error("no results-v* releases found for $repo")

rows = NamedTuple[]
mktempdir() do dir
    for tag in tags
        r = release_row(repo, tag, mkpath(joinpath(dir, tag)))
        if isnothing(r)
            @warn "skipping $tag (no parseable outbreak-size summary)"
            continue
        end
        push!(rows,
            (; tag = r[1], date = r[2], model = r[3],
                median = round(Int, r[4]),
                lo30 = round(Int, r[5]), hi30 = round(Int, r[6]),
                lo60 = round(Int, r[7]), hi60 = round(Int, r[8]),
                lo90 = round(Int, r[9]), hi90 = round(Int, r[10])))
    end
end

out = sort(DataFrame(rows), [:date, :model])
dest = joinpath(@__DIR__, "..", "data", "released_estimates.csv")
CSV.write(dest, out)
println("Wrote $(nrow(out)) releases to data/released_estimates.csv")
show(out; allrows = true)
println()
