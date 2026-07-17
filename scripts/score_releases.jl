#!/usr/bin/env julia
#
# Score every past release's saved one-week-ahead forecasts against the
# now-observed data, plus a persistence baseline, and write tidy score
# tables to `data/forecast_scores.csv`. Also refreshes
# `data/rt_by_release.csv`, the per-release posterior R_T summary, in the
# same style as `data/released_estimates.csv` (see
# `scripts/refresh_releases.jl`).
#
# Both kinds of results release are scored: `results-vX.Y.Z` from a
# version tag and `results-<run number>` from a push to `main`. At most one
# release per data day is kept, keyed on the cut-off each release was built
# from rather than when it was published, since releases sharing a cut-off
# carry the same forecast (see `select_daily_releases` in
# `src/scoring.jl`). Reading the cut-off means fetching every candidate's
# `observations.toml` before the selection, so the driver fetches those
# first and the rest of each winner's assets after.
#
# A release with no `forecast.csv` of its own falls back to the
# reconstruction `scripts/backfill_forecasts.jl` publishes for its code tag
# on the optional `forecasts-backfill` release (`forecast_v1.9.0.csv` for
# `results-v1.9.0`), scored under `<tag> (backfill)`. A saved forecast is
# always preferred over a reconstructed one. That release is absent by
# default, in which case releases without a stored forecast are skipped.
#
# Every release that carries a `forecast.csv` asset (see
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
# Re-run whenever a new release is published: it always recomputes both
# tables from every release currently on GitHub, so re-running is
# idempotent and just refreshes the output CSVs in place.
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

using BVDOutbreakSize: is_results_release, select_daily_releases

const DEFAULT_REPO = "epiforecasts/BVDOutbreakSize"
const FORECAST_ASSET = "forecast.csv"
const DRAWS_ASSET = "posterior_draws.csv"
const OBS_ASSET = "observations.toml"
const BACKFILL_TAG = "forecasts-backfill"

repo = length(ARGS) >= 1 ? ARGS[1] : DEFAULT_REPO

# ----------------------------------------------------------------------
# gh CLI plumbing (mirrors `scripts/refresh_releases.jl`)
# ----------------------------------------------------------------------

## Every results release `gh` reports, as `(tag, created)` pairs. The repo
## publishes a release per code tag too, which carries no results assets,
## so only results releases are kept. Which of them are scored is decided
## by `select_daily_releases` once each one's data cut-off is known.
function results_release_entries(repo)
    out = read(`gh release list -R $repo -L 200 --json tagName,createdAt
        --jq ".[] | [.tagName, .createdAt] | @tsv"`, String)
    entries = Tuple{String, DateTime}[]
    for line in split(strip(out), '\n')
        isempty(strip(line)) && continue
        tag, created = split(line, '\t')
        is_results_release(tag) || continue
        push!(entries,
            (String(tag),
                DateTime(created, dateformat"yyyy-mm-ddTHH:MM:SSZ")))
    end
    return entries
end

## Reconstructed-forecast asset for a tagged release on the optional
## `forecasts-backfill` release: `results-v1.9.0` was built from code tag
## `v1.9.0`, whose archive `scripts/backfill_forecasts.jl` writes as
## `forecast_v1.9.0.csv`. Main-build releases have no code tag, so no
## reconstruction exists for them.
function backfill_asset(tag)
    m = match(r"^results-(v[0-9][0-9A-Za-z.+-]*)$", tag)
    return isnothing(m) ? nothing : "forecast_$(m[1]).csv"
end

## Asset names carried by the `forecasts-backfill` release, empty when it
## has not been published. The release is optional: the backfill is a
## manual, compute-heavy job, so its absence is the normal case and is not
## an error.
function backfill_asset_names(repo)
    out = try
        ## `gh` reports the absent release on stderr, which is the expected
        ## case here, so it is silenced rather than shown as a failure.
        read(
            pipeline(
                `gh release view $BACKFILL_TAG -R $repo --json assets
          --jq ".assets[].name"`; stderr = devnull), String)
    catch
        return String[]
    end
    return [String(strip(l)) for l in split(strip(out), '\n')
            if !isempty(strip(l))]
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
# Release metadata
# ----------------------------------------------------------------------

## The data cut-off a release was built from, read from its
## `observations.toml` (`as_of_date`), or `nothing` when the asset is
## absent, unparseable or carries no cut-off.
function release_cutoff(obs_path)
    isnothing(obs_path) && return nothing
    parsed = try
        TOML.parsefile(obs_path)
    catch
        return nothing
    end
    haskey(parsed, "as_of_date") || return nothing
    return tryparse(Date, string(parsed["as_of_date"]))
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
## a release's `posterior_draws.csv` R_T column, dated by the `cutoff` the
## release was built from. `nothing` when the draws are missing or carry no
## R_T column, as pre-renewal releases do.
function rt_row(tag, draws_path, cutoff)
    isnothing(draws_path) && return nothing
    date = string(cutoff)

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

entries = results_release_entries(repo)
isempty(entries) && error("no releases found for $repo")

## The CURRENT observations manifest is the now-observed truth every
## release's forecast is scored against.
obs = load_observations()
grid_date(day) = obs.cutoff - Day(obs.n - day)

score_rows = NamedTuple[]
overlay_rows = NamedTuple[]
rt_rows = NamedTuple[]
n_scored = 0
n_no_forecast = 0
n_backfilled = 0
n_no_rt = 0

## Assets live under one temp tree for the whole run, so a release's
## `observations.toml` fetched for the selection is still there when its
## forecast is scored. `mktempdir` removes the tree when the process exits.
dir = mktempdir()
tagdir(tag) = mkpath(joinpath(dir, tag))

## Which forecasts a reconstruction is available for, when the optional
## backfill release exists at all.
backfill_names = backfill_asset_names(repo)
isempty(backfill_names) ||
    @info "$BACKFILL_TAG carries $(length(backfill_names)) asset(s)"

## Selection pass: a release is placed on a data day by its own cut-off, so
## every candidate's `observations.toml` is read before any is scored. A
## release whose cut-off cannot be read is dropped rather than scored on an
## assumed day, since its forecast could silently duplicate a kept one; it
## is listed rather than dropped quietly.
dated = Tuple{String, DateTime, Date}[]
undated = String[]
for (tag, created) in entries
    cutoff = release_cutoff(fetch_asset(repo, tag, OBS_ASSET, tagdir(tag)))
    isnothing(cutoff) ? push!(undated, tag) :
    push!(dated, (tag, created, cutoff))
end
isempty(undated) || @info string(
    "ignoring ", length(undated), " release(s) with no readable data ",
    "cut-off: ", join(undated, ", "))

tags = select_daily_releases(dated)
isempty(tags) && error("no results releases with a data cut-off for $repo")
cutoffs = Dict(t => c for (t, _, c) in dated)
println("Scoring $(length(tags)) release(s) of $(length(entries)), " *
        "one per data day.")

for tag in tags
    ## A top-level `for` is a soft scope, so the running counters must be
    ## declared global to update the bindings above rather than shadow them.
    global n_scored, n_no_forecast, n_backfilled, n_no_rt

    ## A saved forecast wins; a reconstruction stands in only where the
    ## release never stored one.
    label = tag
    forecast_path = fetch_asset(repo, tag, FORECAST_ASSET, tagdir(tag))
    if isnothing(forecast_path)
        asset = backfill_asset(tag)
        if !isnothing(asset) && asset in backfill_names
            forecast_path = fetch_asset(
                repo, BACKFILL_TAG, asset, tagdir(tag))
            isnothing(forecast_path) || (label = "$tag (backfill)")
        end
    end

    if isnothing(forecast_path)
        n_no_forecast += 1
        @info "skipping $tag (no stored or reconstructed forecast)"
    else
        result = try
            score_release(label, forecast_path, obs, grid_date)
        catch e
            @warn "skipping $label forecast scoring" exception=e
            nothing
        end
        if !isnothing(result)
            append!(score_rows, result.rows)
            append!(overlay_rows, result.overlay)
            n_scored += 1
            label == tag || (n_backfilled += 1)
            result.skipped > 0 && @info string(
                label, ": skipped ", result.skipped,
                " not-yet-observed group(s)")
        end
    end

    draws_path = fetch_asset(repo, tag, DRAWS_ASSET, tagdir(tag))
    r = try
        rt_row(tag, draws_path, cutoffs[tag])
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

println("Scored $n_scored/$(length(tags)) releases " *
        "($n_no_forecast without a stored or reconstructed forecast, " *
        "$n_backfilled from $BACKFILL_TAG).")
println("R_T summary for $(length(rt_rows))/$(length(tags)) releases " *
        "($n_no_rt without an R_T posterior).")

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
