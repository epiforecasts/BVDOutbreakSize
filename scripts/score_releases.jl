#!/usr/bin/env julia
#
# Score every past release's saved one-week-ahead forecasts against the
# now-observed data, plus a persistence baseline, and write tidy score
# tables to `data/forecast_scores.csv`. Also refreshes
# `data/rt_by_release.csv`, the per-release posterior R_T summary, in the
# same style as `data/released_estimates.csv` (see
# `scripts/refresh_releases.jl`).
#
# Every `results-v*` release that carries a `forecast.csv` asset (see
# `forecast_archive` in `src/forecast.jl`) has one row per (made_date,
# horizon, target_date, stream, draw): the forecast draws for the
# incident "confirmed cases", "confirmed deaths" and "recovered" streams
# and the LEVEL "isolation beds" stream. This script re-scores every one
# of those forecasts against the data that has since arrived, using
# `score_draws` (CRPS, log-CRPS, 50/90% coverage and bias), and compares
# it against a simple persistence baseline: for the incident streams, the
# observed count over the same-length window ending at the forecast's
# cut-off, replicated through a Poisson (a parameter-free, defensible
# spread; autoregression is explicitly out of scope); for isolation beds,
# the last observed occupancy carried forward through the same Poisson
# replication.
#
# MUST BE RE-RUN whenever a new `results-v*` release is cut: it always
# recomputes both tables from every release currently on GitHub, so
# re-running is idempotent and just refreshes the two output CSVs in
# place.
#
# Requires the `gh` CLI, authenticated against the repo. Runs under the
# package's own environment with no extra dependencies: `CSV.jl` is not a
# `[deps]` of the main project (only `scripts/Project.toml` carries it,
# for `refresh_releases.jl`), so this reads and writes the plain,
# unquoted CSV assets by hand instead of adding a dependency.
#
# Usage:
#
#   julia --project=. scripts/score_releases.jl
#   julia --project=. scripts/score_releases.jl owner/repo

using BVDOutbreakSize
using DataFrames
using Dates
using Distributions
using Random
using Statistics: median
using TOML

const DEFAULT_REPO = "epiforecasts/BVDOutbreakSize"
const FORECAST_ASSET = "forecast.csv"
const DRAWS_ASSET = "posterior_draws.csv"
const OBS_ASSET = "observations.toml"

repo = length(ARGS) >= 1 ? ARGS[1] : DEFAULT_REPO

# ----------------------------------------------------------------------
# gh CLI plumbing (mirrors `scripts/refresh_releases.jl`)
# ----------------------------------------------------------------------

## All tagged results releases, newest first, as reported by `gh`.
function results_release_tags(repo)
    out = read(`gh release list -R $repo -L 200 --json tagName
        --jq ".[].tagName"`, String)
    tags = filter(t -> occursin(r"^results-v[0-9]", t),
        split(strip(out), '\n'))
    return tags
end

## Download a single asset from a release into `dir`, returning its path
## (or `nothing` when the release does not carry that asset, or the
## download keeps failing). Retries a couple of times so a transient `gh`
## failure does not drop a release.
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

# ----------------------------------------------------------------------
# Minimal CSV I/O. `CSV.jl` is not available under `--project=.` (see
# above), and every asset this script reads or writes is a plain
# numeric/short-string table with no embedded commas or quoting, so a
# hand-rolled split is safe and avoids adding a dependency.
# ----------------------------------------------------------------------

## Read a simple CSV into `(header, rows)`: `header` is the first line
## split on commas, `rows` a `Vector` of `Vector{String}`, one per
## subsequent non-empty line, each split and trimmed the same way.
function read_simple_csv(path::AbstractString)
    lines = readlines(path)
    isempty(lines) && return (String[], Vector{String}[])
    header = String.(strip.(split(lines[1], ',')))
    rows = [String.(strip.(split(l, ',')))
            for l in lines[2:end] if !isempty(strip(l))]
    return (header, rows)
end

## Column index of `name` in `header`, erroring with a clear message
## when the asset does not carry the expected column.
function col(header, name)
    i = findfirst(==(name), header)
    isnothing(i) && error("column '$name' not found (have $header)")
    return i
end

## Write `cols` (a `Vector` of `name => values` pairs, all the same
## length) to `path` as a plain CSV, overwriting atomically: written to a
## sibling temp file, then renamed into place, so a concurrent reader
## never sees a half-written table.
function write_simple_csv(path::AbstractString, cols)
    n = length(last(cols[1]))
    tmp = path * ".tmp"
    open(tmp, "w") do io
        println(io, join(first.(cols), ','))
        for i in 1:n
            println(io, join((string(last(c)[i]) for c in cols), ','))
        end
    end
    mv(tmp, path; force = true)
end

# ----------------------------------------------------------------------
# Truth lookup against the CURRENT observations
# ----------------------------------------------------------------------

## The four forecast streams `forecast_archive` writes, mapped to the
## observed history each is scored against and whether it is an INCIDENT
## quantity (the new count over the horizon window ending at
## `target_date`) or a LEVEL quantity (isolation beds occupancy AT
## `target_date`).
const STREAM_HISTORY = Dict(
    "confirmed cases" => (:confirmed_history, :incident),
    "confirmed deaths" => (:confirmed_deaths_history, :incident),
    "recovered" => (:recovered_history, :incident),
    "isolation beds" => (:isolation_history, :level))

## Observed cumulative count at or before `date` from a `(days, counts)`
## history, held at the last report (cumulative counts are step-valued
## between vintages): `0` when no vintage predates `date`, i.e. the
## series had not started yet. `grid_date` converts a history's grid
## day-index back to a calendar `Date` (mirrors the local `grid_date`
## helper in `docs/examples/_setup.jl`).
function cum_at(hist, date::Date, grid_date)
    isempty(hist.days) && return 0
    idx = findlast(d -> grid_date(d) <= date, hist.days)
    isnothing(idx) && return 0
    return hist.counts[idx]
end

## Observed truth for one (stream, made_date, target_date) forecast group
## against the CURRENT `obs`, or `missing` when `target_date` is beyond
## the current data cut-off (not yet observed, so the group cannot be
## scored). Incident streams floor the new-count difference at zero, the
## same convention `forecast.jl` uses for `*_new` columns.
function truth_at(obs, grid_date, stream, made_date, target_date)
    target_date > obs.cutoff && return missing
    hist, kind = STREAM_HISTORY[stream]
    h = getproperty(obs, hist)
    kind == :level && return Float64(cum_at(h, target_date, grid_date))
    new = cum_at(h, target_date, grid_date) - cum_at(h, made_date, grid_date)
    return Float64(max(new, 0))
end

# ----------------------------------------------------------------------
# Persistence baseline
# ----------------------------------------------------------------------

## Persistence baseline draws for one forecast group: `n` replicates from
## a Poisson centred on the last-observed value, chosen as the simplest
## defensible spread (no free dispersion parameter to justify, and no
## autoregression, which is explicitly out of scope). For an INCIDENT
## stream the centre is the observed count over the horizon-length window
## ENDING at `made_date` (the same window length as the forecast, one
## horizon earlier); for the LEVEL stream ("isolation beds") it is the
## last observed occupancy at or before `made_date`, carried forward
## unchanged.
function baseline_draws(obs, grid_date, stream, made_date, horizon, n, rng)
    hist, kind = STREAM_HISTORY[stream]
    h = getproperty(obs, hist)
    centre = if kind == :level
        Float64(cum_at(h, made_date, grid_date))
    else
        prior = made_date - Day(horizon)
        Float64(max(
            cum_at(h, made_date, grid_date) - cum_at(h, prior, grid_date), 0))
    end
    return Float64.(rand(rng, Poisson(centre), n))
end

# ----------------------------------------------------------------------
# Score one release's forecast archive
# ----------------------------------------------------------------------

## Score every (made_date, horizon, target_date, stream) group in a
## release's `forecast.csv` against `obs` (the CURRENT observations),
## returning one row per (group, model) with `model` "ours" or
## "baseline". Groups whose `target_date` is not yet observed are
## skipped (counted in `.skipped` for the caller to log).
function score_release(tag, forecast_path, obs, grid_date)
    header, rows = read_simple_csv(forecast_path)
    made_i = col(header, "made_date")
    hor_i = col(header, "horizon")
    tgt_i = col(header, "target_date")
    str_i = col(header, "stream")
    val_i = col(header, "value")

    ## Group rows by (made_date, horizon, target_date, stream); a plain
    ## `Dict` keeps this independent of DataFrame group-by quirks.
    groups = Dict{Tuple{Date, Int, Date, String}, Vector{Float64}}()
    for r in rows
        key = (Date(r[made_i]), parse(Int, r[hor_i]),
            Date(r[tgt_i]), r[str_i])
        push!(get!(groups, key, Float64[]), parse(Float64, r[val_i]))
    end

    out = NamedTuple[]
    overlay = NamedTuple[]
    skipped = 0
    for (key, draws) in groups
        made_date, horizon, target_date, stream = key
        truth = truth_at(obs, grid_date, stream, made_date, target_date)
        if truth === missing
            skipped += 1
            continue
        end
        rng = MersenneTwister(hash((tag, stream, horizon, made_date)))
        base = baseline_draws(
            obs, grid_date, stream, made_date, horizon, length(draws), rng)
        for (model, samples) in (("ours", draws), ("baseline", base))
            s = score_draws(truth, samples)
            push!(out,
                (; release = tag, made_date, stream, horizon, target_date,
                    model, crps = s.crps, log_crps = s.log_crps,
                    coverage_50 = s.coverage_50, coverage_90 = s.coverage_90,
                    bias = s.bias, n_samples = s.n))
            ## Quantile summary of the same draws for the forecasts-versus-now
            ## overlay plot, alongside the observed truth so the docs build
            ## reads a plain table rather than re-pulling the release assets.
            q = posterior_summary(samples)
            r2(x) = round(x; digits = 2)
            push!(overlay,
                (; release = tag, made_date, stream, horizon, target_date,
                    model, observed = r2(truth), median = r2(median(samples)),
                    lo30 = r2(q.lo30), hi30 = r2(q.hi30),
                    lo60 = r2(q.lo60), hi60 = r2(q.hi60),
                    lo90 = r2(q.lo90), hi90 = r2(q.hi90)))
        end
    end
    return (; rows = out, overlay, skipped)
end

# ----------------------------------------------------------------------
# Per-release posterior R_T summary (mirrors `refresh_releases.jl`)
# ----------------------------------------------------------------------

## One `(tag, date, median, lo30, hi30, lo60, hi60, lo90, hi90)` row from
## a release's `posterior_draws.csv` R_T column and its own
## `observations.toml` cut-off, or `nothing` when either asset is
## missing or unparseable.
function rt_row(tag, draws_path, obs_path)
    (isnothing(draws_path) || isnothing(obs_path)) && return nothing
    parsed = TOML.parsefile(obs_path)
    haskey(parsed, "as_of_date") || return nothing
    date = string(parsed["as_of_date"])

    header, rows = read_simple_csv(draws_path)
    "R_T" in header || return nothing
    i = col(header, "R_T")
    rt = [parse(Float64, r[i]) for r in rows]
    isempty(rt) && return nothing

    s = posterior_summary(rt)
    r3(x) = round(x; digits = 3)
    return (tag, date, r3(median(rt)), r3(s.lo30), r3(s.hi30),
        r3(s.lo60), r3(s.hi60), r3(s.lo90), r3(s.hi90))
end

# ----------------------------------------------------------------------
# Driver
# ----------------------------------------------------------------------

tags = results_release_tags(repo)
isempty(tags) && error("no results-v* releases found for $repo")

## The CURRENT observations manifest is the now-observed truth every
## release's forecast is scored against.
obs = load_observations()
grid_date(day) = obs.cutoff - Day(obs.n - day)

score_rows = NamedTuple[]
overlay_rows = NamedTuple[]
rt_rows = NamedTuple[]
n_scored = 0
n_no_forecast = 0
n_no_rt = 0

mktempdir() do dir
    ## The `do` block is a function body (hard scope), so the running
    ## counters must be declared global to update the top-level bindings.
    global n_scored, n_no_forecast, n_no_rt
    for tag in tags
        tagdir = mkpath(joinpath(dir, tag))

        forecast_path = fetch_asset(repo, tag, FORECAST_ASSET, tagdir)
        if isnothing(forecast_path)
            n_no_forecast += 1
            @info "skipping $tag (no forecast.csv asset)"
        else
            result = try
                score_release(tag, forecast_path, obs, grid_date)
            catch e
                @warn "skipping $tag forecast scoring" exception=e
                nothing
            end
            if !isnothing(result)
                append!(score_rows, result.rows)
                append!(overlay_rows, result.overlay)
                n_scored += 1
                result.skipped > 0 && @info string(
                    tag, ": skipped ", result.skipped,
                    " not-yet-observed group(s)")
            end
        end

        draws_path = fetch_asset(repo, tag, DRAWS_ASSET, tagdir)
        obs_path = fetch_asset(repo, tag, OBS_ASSET, tagdir)
        r = try
            rt_row(tag, draws_path, obs_path)
        catch e
            @warn "skipping $tag R_T summary" exception=e
            nothing
        end
        if isnothing(r)
            n_no_rt += 1
        else
            push!(rt_rows,
                (; release = r[1], date = r[2], median = r[3],
                    lo30 = r[4], hi30 = r[5], lo60 = r[6], hi60 = r[7],
                    lo90 = r[8], hi90 = r[9]))
        end
    end
end

println("Scored $n_scored/$(length(tags)) releases " *
        "($n_no_forecast without a forecast.csv asset).")
println("R_T summary for $(length(rt_rows))/$(length(tags)) releases " *
        "($n_no_rt without posterior draws/observations).")

## `data/forecast_scores.csv`: one row per scored (release, made_date,
## stream, horizon, model) group.
scores = if isempty(score_rows)
    DataFrame(release = String[], made_date = Date[], stream = String[],
        horizon = Int[], target_date = Date[], model = String[],
        crps = Float64[], log_crps = Float64[], coverage_50 = Float64[],
        coverage_90 = Float64[], bias = Float64[], n_samples = Int[])
else
    sort(DataFrame(score_rows),
        [:release, :made_date, :stream, :horizon, :model])
end
scores_dest = joinpath(@__DIR__, "..", "data", "forecast_scores.csv")
write_simple_csv(scores_dest,
    [:release => scores.release,
        :made_date => string.(scores.made_date),
        :stream => scores.stream,
        :horizon => scores.horizon,
        :target_date => string.(scores.target_date),
        :model => scores.model,
        :crps => scores.crps,
        :log_crps => scores.log_crps,
        :coverage_50 => scores.coverage_50,
        :coverage_90 => scores.coverage_90,
        :bias => scores.bias,
        :n_samples => scores.n_samples])
println("Wrote $(nrow(scores)) scored forecasts to " *
        "data/forecast_scores.csv")

## `data/forecast_overlay.csv`: median and 30/60/90% bounds of each forecast
## group with the observed truth, for the forecasts-versus-now overlay plots.
overlay = if isempty(overlay_rows)
    DataFrame(release = String[], made_date = Date[], stream = String[],
        horizon = Int[], target_date = Date[], model = String[],
        observed = Float64[], median = Float64[], lo30 = Float64[],
        hi30 = Float64[], lo60 = Float64[], hi60 = Float64[],
        lo90 = Float64[], hi90 = Float64[])
else
    sort(DataFrame(overlay_rows),
        [:stream, :made_date, :horizon, :model])
end
overlay_dest = joinpath(@__DIR__, "..", "data", "forecast_overlay.csv")
write_simple_csv(overlay_dest,
    [:release => overlay.release,
        :made_date => string.(overlay.made_date),
        :stream => overlay.stream,
        :horizon => overlay.horizon,
        :target_date => string.(overlay.target_date),
        :model => overlay.model,
        :observed => overlay.observed,
        :median => overlay.median,
        :lo30 => overlay.lo30, :hi30 => overlay.hi30,
        :lo60 => overlay.lo60, :hi60 => overlay.hi60,
        :lo90 => overlay.lo90, :hi90 => overlay.hi90])
println("Wrote $(nrow(overlay)) forecast-overlay rows to " *
        "data/forecast_overlay.csv")

## `data/rt_by_release.csv`: mirrors `data/released_estimates.csv`.
rt = if isempty(rt_rows)
    DataFrame(release = String[], date = String[], median = Float64[],
        lo30 = Float64[], hi30 = Float64[], lo60 = Float64[],
        hi60 = Float64[], lo90 = Float64[], hi90 = Float64[])
else
    sort(DataFrame(rt_rows), [:date])
end
rt_dest = joinpath(@__DIR__, "..", "data", "rt_by_release.csv")
write_simple_csv(rt_dest,
    [:release => rt.release, :date => rt.date, :median => rt.median,
        :lo30 => rt.lo30, :hi30 => rt.hi30, :lo60 => rt.lo60,
        :hi60 => rt.hi60, :lo90 => rt.lo90, :hi90 => rt.hi90])
println("Wrote $(nrow(rt)) release R_T summaries to " *
        "data/rt_by_release.csv")
