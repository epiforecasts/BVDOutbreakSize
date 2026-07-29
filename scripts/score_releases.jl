#!/usr/bin/env julia
#
# Score every past release's saved one-week-ahead forecasts against the
# now-observed data, plus a persistence baseline, and write tidy score
# tables to `data/forecast_scores.csv`. Each score row also carries a
# log-CRPS relative skill to the persistence baseline,
# `log_rel_to_baseline`; the comparison against a stream's individual
# fit is computed at aggregation time, not stored per row. Also refreshes
# `data/rt_by_release.csv` and `data/r0_by_release.csv`, the per-release
# posterior R_T and initial-R0 summaries, in the same style as
# `data/released_estimates.csv` (see `scripts/refresh_releases.jl`).
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
# The forecast a release is scored from comes from the first of these it
# carries: `stream_forecasts.csv`, which holds every model's forecast in
# the archive schema plus a `fit` column; `forecast.csv`, the joint model's
# forecast alone, read as `fit = "joint"`; or the reconstruction
# `scripts/backfill_forecasts.jl` publishes for the release's code tag on
# the optional `forecasts-backfill` release (`forecast_v1.9.0.csv` for
# `results-v1.9.0`), scored under `<tag> (backfill)`. A saved forecast is
# always preferred over a reconstructed one; the backfill release is absent
# by default, in which case releases without a stored forecast are skipped.
#
# A forecast archive (see `forecast_archive` in `src/forecast.jl`) has one
# row per (made_date, horizon, target_date, stream, draw, fit): the
# forecast draws for the incident "confirmed cases", "confirmed deaths" and
# "recovered" streams and the LEVEL "isolation beds" stream. This script
# scores each fit of each stream against the data that has since arrived,
# using `score_draws` (CRPS, log-CRPS, its dispersion/overprediction/
# underprediction decomposition, 50/90% coverage and bias), against one
# persistence baseline per stream (`fit = "baseline"`, `baseline_draws`):
# the centre carried forward as before (the observed count over the
# same-length window ending at the forecast's cut-off for the incident
# streams, the last observed occupancy for isolation beds), with the
# spread simulated as an iterated day-by-day random walk from the stream's
# own past first differences, in the style of the COVID-19 Forecast Hub
# baseline, falling back to a Poisson spread when a stream has too little
# history yet. The baseline is built from the release's own
# `observations.toml` vintage, frozen no later than each group's own
# made_date (`vintage_observations`), not from the current, possibly
# later-revised manifest, so a correction or backfill that landed after
# the RELEASE's own cut-off cannot leak into its baseline (a residual gap
# remains for the frozen archive, where made_date can sit weeks before the
# release's own cut-off; see the driver's `tag_obs_path` comment). A
# stream no individual
# model fits, "recovered", is scored for the joint and the baseline only. A
# forecast group whose target runs past the stream's own reporting
# coverage (`truth_at`, `stream_coverage_end`) is not scored at all, the
# same way a not-yet-observed target is not.
#
# Each release's `stream_estimates.csv`, when present, also feeds the
# per-fit `data/rt_by_release_by_stream.csv` and
# `data/size_by_release_by_stream.csv` R_T and C_T overlays.
#
# Re-run whenever a new release is published: it always recomputes the
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
const FROZEN_FORECAST_ASSET = "forecast_frozen.csv"
const STREAM_FORECAST_ASSET = "stream_forecasts.csv"
const STREAM_ESTIMATES_ASSET = "stream_estimates.csv"
const DRAWS_ASSET = "posterior_draws.csv"
const OBS_ASSET = "observations.toml"
const BACKFILL_TAG = "forecasts-backfill"

## Values of the `fit` column, which names the model a row's forecast or
## estimate came from: a spec id, or `baseline` for the persistence
## baseline, which is no model's output. `joint` is also the fit a plain
## `forecast.csv` carries, that asset predating the per-stream fits.
const JOINT_FIT = "joint"
const BASELINE_FIT = "baseline"
## The fit label the frozen-fit forecasts (`forecast_frozen.csv`) are scored
## under. That asset carries no `fit` column (it is the joint model at each
## frozen cut-off), so the label is applied at scoring time.
const FROZEN_FIT = "frozen"

## Release tags whose forecast is excluded from scoring because the
## RECONSTRUCTION failed, not because the model performed badly: a chain
## that did not sample properly rather than a genuine forecast. Keyed on
## the release tag itself (before the "(backfill)" label is applied), each
## value names the signature that identifies it as a failed reconstruction,
## kept here rather than in a comment so it prints in the run log. Only the
## FORECAST is skipped for a listed tag; its `posterior_draws.csv` (the R_T
## and R0 summaries) is a separate, unaffected asset from the real release
## and is still scored normally.
const _FAILED_RECONSTRUCTIONS = Dict(
    "results-v1.6.0" =>
    "the reconstructed joint chain forecasts a " *
    "nearly-zero median at every horizon and every stream, with the " *
    "upper predictive tail occasionally exploding to five- and " *
    "six-digit values (e.g. a confirmed-cases CRPS above 100000 at " *
    "the 28-day horizon); that combination is the signature of a " *
    "chain that failed to sample properly, not a real forecast.")

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
## `observations.toml` (`as_of_date`), or `nothing` when the file is
## unparseable or carries no cut-off. The caller handles a missing file
## separately, since a failed fetch and a cut-off-less file are not the
## same event; passing `nothing` here still returns `nothing`.
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

## The forecast streams the archive writes, mapped to the observed history
## each is scored against (a `(; days, counts)` field on `obs`) and whether
## it is an INCIDENT quantity (the new count over the horizon window ending
## at `target_date`) or a LEVEL quantity (isolation beds occupancy AT
## `target_date`). The single-stream fits forecast six labels: the four
## here plus "reported cases" (the cases fit) and "suspected deaths" (the
## deaths fit), scored on the same incident basis as "confirmed cases".
## "exports" is scored too, but its truth is built by `stream_history` from
## the dated detection series rather than named here.
##
## The confirmed/suspect ward split ("treatment beds" and "isolation beds
## (suspected)") is a level quantity like the total: each ward occupancy is
## scored against its Tableau 6 occupancy sub-stock (`dont confirmes` /
## `dont suspects`, which sum to the total `isolation_history`). These are
## present only once a forecast carries the split, so before then no archive
## row bears the label and it is simply never scored.
const STREAM_HISTORY = Dict(
    "reported cases" => (:reported_history, :incident),
    "suspected deaths" => (:deaths_history, :incident),
    "confirmed cases" => (:confirmed_history, :incident),
    "confirmed deaths" => (:confirmed_deaths_history, :incident),
    "recovered" => (:recovered_history, :incident),
    "isolation beds" => (:isolation_history, :level),
    "treatment beds" => (:treatment_confirmed_incare_history, :level),
    "isolation beds (suspected)" =>
        (:treatment_suspect_incare_history, :level))

## Streams the archive can carry that have no `(; days, counts)` field, so
## `stream_history` assembles their truth. "exports" is the dated Uganda
## import series (`export_case_days`, sorted grid day-indices of detected
## imports): the same detections `exports_model` fits `expected_exports_T`
## to, so scoring the export forecast against their cumulative count uses
## the model's own truth. Each detection is one cumulative step, matching
## `dated_event_bins` in `src/models/observations.jl`.
const STREAM_ASSEMBLED = Dict("exports" => :incident)

## Whether this script has a truth source for `stream`, i.e. it is either a
## named history or an assembled one. An archived label that is neither
## (schema drift, a future stream) has no truth and is skipped, never
## scored against a wrong one.
function has_stream_truth(stream)
    haskey(STREAM_HISTORY, stream) || haskey(STREAM_ASSEMBLED, stream)
end

## The `(; days, counts)` cumulative history a stream is scored against and
## its incident/level kind, assembling the truth for streams stored in
## another shape. Errors on a stream with no truth source; callers guard
## with `has_stream_truth` first.
function stream_history(obs, stream)
    if haskey(STREAM_ASSEMBLED, stream)
        stream == "exports" || error("no assembler for stream '$stream'")
        d = obs.export_case_days
        return ((; days = d, counts = collect(1:length(d))),
            STREAM_ASSEMBLED[stream])
    end
    field, kind = STREAM_HISTORY[stream]
    return (getproperty(obs, field), kind)
end

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

## The last date a stream's own truth source was still being updated: the
## date of its last vintage, since beyond it the cumulative total is only
## ever repeated at its last reported value rather than genuinely observed.
## For the exports stream, assembled by `stream_history` from the dated
## detection series, this is the date of the last detection, since the
## manifest documents a later known import deliberately excluded as later
## than the cut-off that series is anchored to; the same generic rule
## reads that correctly with no special case for the stream. A history
## with no vintages at all has no coverage.
function stream_coverage_end(obs, grid_date, stream)
    h, _ = stream_history(obs, stream)
    isempty(h.days) && return typemin(Date)
    return grid_date(maximum(h.days))
end

## The first date a stream's own truth source carries an observation: the
## date of its first vintage. Before it `cum_at` reads zero because no
## vintage predates the date, which is the absence of a series rather than a
## count of zero. A window opening before this date therefore measures its
## increment from a non-observation, and a persistence baseline centred
## there carries zero forward against a real value, so such a window is not
## scorable. A history with no vintages at all has no coverage.
##
## An assembled stream is exempt. "exports" is a dated list of detections
## rather than a series of reporting vintages, so a date before its first
## detection carries the true statement that no export had been detected
## yet, not the absence of a source. Only its end anchor constrains it.
function stream_coverage_start(obs, grid_date, stream)
    haskey(STREAM_ASSEMBLED, stream) && return typemin(Date)
    h, _ = stream_history(obs, stream)
    isempty(h.days) && return typemax(Date)
    return grid_date(minimum(h.days))
end

## Observed truth for one (stream, made_date, target_date) forecast group
## against the CURRENT `obs`, or a `Symbol` naming why the group is not
## scorable: `:not_yet_observed` when `target_date` is beyond the current
## data cut-off, `:stopped_reporting` when it is beyond the stream's own
## coverage (see `stream_coverage_end`), so an unmoved cumulative total is
## not read as an observed zero, or `:not_yet_reporting` when the window
## opens before the stream began being reported (see
## `stream_coverage_start`), so a series that had not started is not read as
## a zero either. The window is scored only where the stream's own reporting
## covers it at both ends. Incident streams floor the new-count difference at
## zero, the same convention `forecast.jl` uses for `*_new` columns.
function truth_at(obs, grid_date, stream, made_date, target_date)
    target_date > obs.cutoff && return :not_yet_observed
    target_date > stream_coverage_end(obs, grid_date, stream) &&
        return :stopped_reporting
    made_date < stream_coverage_start(obs, grid_date, stream) &&
        return :not_yet_reporting
    h, kind = stream_history(obs, stream)
    kind == :level && return Float64(cum_at(h, target_date, grid_date))
    new = cum_at(h, target_date, grid_date) - cum_at(h, made_date, grid_date)
    return Float64(max(new, 0))
end

# ----------------------------------------------------------------------
# Persistence baseline
# ----------------------------------------------------------------------

## The `(value, window)` first differences between consecutive vintages of
## a `(; days, counts)` history at or before `made_date`: `value` the
## signed change in the count between the two vintages and `window` the
## number of days between them, kept apart because vintages do not arrive
## on a fixed cadence. Fewer than two vintages by `made_date` yields no
## differences.
function _history_diffs(hist, grid_date, made_date)
    out = Tuple{Float64, Int}[]
    n = length(hist.days)
    n < 2 && return out
    for i in 2:n
        d_prev = grid_date(hist.days[i - 1])
        d_this = grid_date(hist.days[i])
        d_this > made_date && break
        window = Dates.value(d_this - d_prev)
        window <= 0 && continue
        push!(out, (Float64(hist.counts[i] - hist.counts[i - 1]), window))
    end
    return out
end

## Persistence baseline draws for one forecast group, in the style of the
## COVID-19 Forecast Hub baseline (`COVIDhub-baseline`): what today looks
## like versus yesterday, where yesterday is itself a forecast from the day
## before it, iterated day by day out to the forecast horizon — a
## zero-drift random walk whose step distribution comes from the stream's
## own reporting history, rather than a Poisson's few-percent spread at a
## count of several hundred (which makes the baseline a near-point forecast
## and any skill ratio against it explode).
##
## The centre is unchanged: for an INCIDENT stream, the observed count over
## the horizon-length window ENDING at `made_date`; for the LEVEL stream
## ("isolation beds"), the last observed occupancy at or before
## `made_date`. `obs` must reflect only data available AT `made_date` — the
## caller is responsible for passing a `made_date`-vintage manifest (see
## `vintage_observations`), not the current, possibly later-revised one, so
## the baseline never sees a correction that landed after the forecast was
## made.
##
## The spread simulates the walk explicitly: each of the stream's past
## first differences (see `_history_diffs`), a `window`-day change, is
## converted to a ONE-DAY step by `value / sqrt(window)` rather than
## `value / window`. Under a zero-drift random walk a `window`-day change
## has variance `window * sigma^2` for the walk's own one-day variance
## `sigma^2`, so dividing by `sqrt(window)` (not `window`) is what recovers
## an estimate of `sigma` itself; dividing by `window` instead would still
## shrink with history length but by the wrong power, understating the
## spread by a further factor of `window`. The per-day steps are
## symmetrised about zero, so a run of only-rising or only-falling history
## does not bias the walk one way; one predictive draw then sums `horizon`
## independent daily steps sampled with replacement from that step pool,
## the same "iterate day by day to the horizon" construction the Hub
## baseline uses. Summing `horizon` iid steps of variance `sigma^2` gives a
## total variance `horizon * sigma^2`, so the spread grows with the square
## root of the horizon, matching a single draw rescaled by
## `sqrt(horizon / window)` (the construction this replaced) in expectation
## while genuinely iterating day by day as the comment above describes,
## rather than only matching its first two moments.
##
## Centre-versus-Hub check: `COVIDhub-baseline` centres each target on the
## single most recent observation; this baseline instead centres an
## INCIDENT stream on the observed count over the whole horizon-length
## window ending at `made_date`, since the streams here are noisy, sparse
## cumulative counts (a handful of vintages, days apart) rather than a
## dense daily series, so "yesterday's single value" would be a much
## noisier centre than the last `horizon` days pooled. The LEVEL stream
## (isolation beds) already centres on the single last observed occupancy,
## matching the Hub convention directly. Both centres were left as
## pre-existing; only the spread changed here.
##
## Falls back to the plain Poisson draw when fewer than three past
## differences are available, so an early made-date still scores.
function baseline_draws(obs, grid_date, stream, made_date, horizon, n, rng)
    h, kind = stream_history(obs, stream)
    centre = if kind == :level
        Float64(cum_at(h, made_date, grid_date))
    else
        prior = made_date - Day(horizon)
        Float64(max(
            cum_at(h, made_date, grid_date) - cum_at(h, prior, grid_date), 0))
    end

    diffs = _history_diffs(h, grid_date, made_date)
    length(diffs) < 3 && return Float64.(rand(rng, Poisson(centre), n))

    steps = [d / sqrt(window) for (d, window) in diffs]
    pool = vcat(steps, -steps)
    draws = Vector{Float64}(undef, n)
    for i in 1:n
        step = zero(Float64)
        for _ in 1:horizon
            step += rand(rng, pool)
        end
        draws[i] = centre + step
    end
    return max.(draws, 0.0)
end

# ----------------------------------------------------------------------
# Score one release's forecast archive
# ----------------------------------------------------------------------

## One score row and one overlay row for a set of predictive `samples`
## against the observed `truth`, pushed onto `out` and `overlay`.
function push_scored!(out, overlay, tag, key, fit, samples, truth)
    made_date, horizon, target_date, stream = key
    s = score_draws(truth, samples)
    push!(out,
        (; release = tag, made_date, stream, horizon, target_date,
            fit, crps = s.crps, log_crps = s.log_crps,
            dispersion = s.dispersion, overprediction = s.overprediction,
            underprediction = s.underprediction,
            coverage_50 = s.coverage_50, coverage_90 = s.coverage_90,
            bias = s.bias, n_samples = s.n))
    ## Quantile summary of the same draws for the forecasts-versus-now
    ## overlay plot, alongside the observed truth so the docs build reads a
    ## plain table rather than re-pulling the release assets.
    q = posterior_summary(samples)
    r2(x) = round(x; digits = 2)
    push!(overlay,
        (; release = tag, made_date, stream, horizon, target_date,
            fit, observed = r2(truth), median = r2(median(samples)),
            lo30 = r2(q.lo30), hi30 = r2(q.hi30),
            lo60 = r2(q.lo60), hi60 = r2(q.hi60),
            lo90 = r2(q.lo90), hi90 = r2(q.hi90)))
end

## In-process cache of vintage manifest loads, keyed by `(path, made_date)`:
## a release forecasts from one made date (or a handful, for the frozen
## archive), not once per scored group, so the same freeze is reused rather
## than re-parsed per (stream, horizon) row.
const _VINTAGE_CACHE = Dict{Tuple{String, Date}, Any}()

## The observation manifest and matching `grid_date` as they stood at
## `made_date`, for the persistence baseline: loaded from `obs_path` (a
## release's own `observations.toml` snapshot, already fetched by the
## driver to read its cut-off) and truncated no later than `made_date` with
## `load_observations`' `cutoff_date`, the same freezing `freeze_observations`
## does. This is what makes the baseline honest — it sees only the
## revisions and vintages that had landed by `made_date`, not whatever the
## manifest says today after any later correction or backfill, which would
## otherwise leak information the forecast itself never had.
##
## `obs_path === nothing` falls back to the given `obs`/`grid_date` (the
## CURRENT manifest) unchanged, e.g. from a test's synthetic `obs` NamedTuple
## that carries no manifest file to load, or a caller that has no per-release
## snapshot available.
##
## A `made_date` before the manifest's earliest vintage is valid: `history`
## and `event_days` (in `load_observations`) simply return empty series for
## it, the same "no vintage yet" state `cum_at`/`_history_diffs` already
## handle, so `baseline_draws` falls back to its Poisson floor rather than
## erroring.
function vintage_observations(obs_path, made_date, obs, grid_date)
    isnothing(obs_path) && return obs, grid_date
    key = (obs_path, made_date)
    ov = get!(_VINTAGE_CACHE, key) do
        load_observations(obs_path; cutoff_date = made_date)
    end
    vintage_grid_date(day) = ov.cutoff - Day(ov.n - day)
    return ov, vintage_grid_date
end

## Score every (made_date, horizon, target_date, stream, fit) group in a
## release's forecast archive against `obs` (the CURRENT observations),
## returning one row per group plus one persistence-baseline row per
## (made_date, horizon, target_date, stream), tagged `fit = "baseline"`.
##
## `fit` names the model a forecast came from. A `stream_forecasts.csv`
## carries it per row; a plain `forecast.csv` predates the per-stream fits
## and holds the joint model's forecast alone, so its rows are read as
## `default_fit` (`joint`). A stream no individual model fits, "recovered",
## therefore yields `joint` and `baseline` rows and no others: the absent fit
## is simply a fit that contributed no rows. The frozen-fit archive
## (`forecast_frozen.csv`) carries no `fit` column either, so it is scored
## with `default_fit = "frozen"`.
##
## `vintage_obs_path`, when given, is the release's own `observations.toml`
## snapshot: the persistence baseline is built from THAT manifest, frozen no
## later than each group's own `made_date` (see `vintage_observations`),
## rather than from `obs`, so a revision or backfill that landed after
## `made_date` cannot leak into the baseline. `obs`/`grid_date` are still used
## for the TRUTH every fit (including the baseline) is scored against, which
## is correctly the now-observed data regardless.
##
## Groups whose `target_date` is not yet observed are skipped (counted in
## `.skipped` for the caller to log), and so are groups whose `target_date`
## runs past the stream's own reporting coverage (counted in `.stopped`),
## kept apart because the two are not the same finding: one is a target the
## world has not reached yet, the other is a stream that stopped being
## updated before the target, whose unmoved cumulative total is not an
## observed zero (see `truth_at`, `stream_coverage_end`).
function score_release(tag, forecast_path, obs, grid_date;
        default_fit = JOINT_FIT, vintage_obs_path = nothing)
    header, rows = read_simple_csv(forecast_path)
    made_i = col(header, "made_date")
    hor_i = col(header, "horizon")
    tgt_i = col(header, "target_date")
    str_i = col(header, "stream")
    val_i = col(header, "value")
    fit_i = findfirst(==("fit"), header)

    ## Group rows by (made_date, horizon, target_date, stream) and fit; a
    ## plain `Dict` keeps this independent of DataFrame group-by quirks. A
    ## row whose stream label has no truth source costs itself, not the
    ## release: it is dropped here so the scoring loop only sees scorable
    ## streams, and the dropped labels are logged once per release rather
    ## than aborting on the first (an unmapped label used to KeyError
    ## through the whole release).
    Key = Tuple{Date, Int, Date, String}
    groups = Dict{Key, Dict{String, Vector{Float64}}}()
    unknown = 0
    seen_unknown = Set{String}()
    for r in rows
        stream = r[str_i]
        if !has_stream_truth(stream)
            unknown += 1
            push!(seen_unknown, stream)
            continue
        end
        key = (Date(r[made_i]), parse(Int, r[hor_i]), Date(r[tgt_i]), stream)
        fit = isnothing(fit_i) ? default_fit : r[fit_i]
        byfit = get!(() -> Dict{String, Vector{Float64}}(), groups, key)
        push!(get!(byfit, fit, Float64[]), parse(Float64, r[val_i]))
    end
    unknown > 0 && @warn string(
        tag, ": dropped ", unknown, " row(s) with no truth source for ",
        "stream label(s): ", join(sort(collect(seen_unknown)), ", "))

    out = NamedTuple[]
    overlay = NamedTuple[]
    skipped = 0
    stopped = 0
    unstarted = 0
    for (key, byfit) in groups
        made_date, horizon, target_date, stream = key
        truth = truth_at(obs, grid_date, stream, made_date, target_date)
        if truth === :not_yet_observed
            skipped += length(byfit)
            continue
        elseif truth === :stopped_reporting
            stopped += length(byfit)
            continue
        elseif truth === :not_yet_reporting
            unstarted += length(byfit)
            continue
        end
        for fit in sort(collect(keys(byfit)))
            push_scored!(out, overlay, tag, key, fit, byfit[fit], truth)
        end

        ## The baseline follows from the stream's own history, not from any
        ## fit, so every fit of one stream is scored against the same one
        ## rather than against a copy per fit. It is drawn at the width of
        ## the widest fit so its resolution never limits the comparison.
        ## Built from the made_date vintage, not the current manifest (see
        ## `vintage_observations`), so it cannot see a revision that landed
        ## after the forecast was made.
        n = maximum(length, values(byfit))
        rng = MersenneTwister(hash((tag, stream, horizon, made_date)))
        vobs, vgrid_date = vintage_observations(
            vintage_obs_path, made_date, obs, grid_date)
        base = baseline_draws(
            vobs, vgrid_date, stream, made_date, horizon, n, rng)
        push_scored!(out, overlay, tag, key, BASELINE_FIT, base, truth)
    end
    return (; rows = out, overlay, skipped, stopped, unstarted)
end

# ----------------------------------------------------------------------
# Relative forecast skill on the log-CRPS scale
# ----------------------------------------------------------------------

## The per-row relative skill against the persistence baseline, on the
## log-CRPS scale (`log_crps`, CRPS of `log1p`-transformed draws), for every
## row of a scored table, relative to other fits of the SAME (release,
## made_date, stream, horizon) group: a fit's `log_crps` over its stream's
## persistence-baseline `log_crps` in the group. Below 1 beats the baseline.
## The baseline row itself is `1.0`. Every scorable stream carries a
## baseline, so this is defined for every row of a group whose baseline was
## scored, and `missing` otherwise. A `missing` also results when the
## baseline's own `log_crps` is zero or non-finite, or the ratio itself
## comes out non-finite, the same degenerate-denominator guard `_safe_ratio`
## applies to the aggregate ratios in `src/scoring.jl`, so a chance
## zero-CRPS baseline row cannot write `Inf`/`NaN` into the CSV.
##
## The comparison against a stream's individual fit is not a per-row ratio:
## it is computed at aggregation time from matched group means, alongside
## the baseline comparison, by the table builders in `src/scoring.jl`.
##
## Returned as a `Vector{Union{Missing, Float64}}` aligned to `df`'s rows.
function rel_to_baseline_columns(df; baseline_fit = BASELINE_FIT)
    n = nrow(df)
    to_base = Vector{Union{Missing, Float64}}(missing, n)
    key(i) = (df.release[i], df.made_date[i], df.stream[i], df.horizon[i])
    groups = Dict{Any, Vector{Int}}()
    for i in 1:n
        push!(get!(() -> Int[], groups, key(i)), i)
    end
    for idx in values(groups)
        base_i = findfirst(i -> df.fit[i] == baseline_fit, idx)
        if !isnothing(base_i)
            base_lc = df.log_crps[idx[base_i]]
            valid_base = isfinite(base_lc) && base_lc != 0
            for i in idx
                to_base[i] = if df.fit[i] == baseline_fit
                    1.0
                elseif valid_base
                    r = df.log_crps[i] / base_lc
                    isfinite(r) ? r : missing
                else
                    missing
                end
            end
        end
    end
    return to_base
end

## Format one relative-skill value for the plain CSV: an empty cell for
## `missing` (an undefined ratio, e.g. a stream with no individual fit),
## else the ratio rounded to four places.
rel_to_baseline_cell(x) = ismissing(x) ? "" : string(round(x; digits = 4))

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

## The chain key for the renewal walk's base, `exp` of which is the
## established initial reproduction number R0 (see `reconstruct_rt` in
## `src/plots.jl`). Its draws ride along in `posterior_draws.csv` under this
## column name.
const LOG_R0_COL = "rt_state.log_R0"

## One `(tag, date, median, lo30, hi30, lo60, hi60, lo90, hi90)` row of the
## established initial reproduction number R0 = `exp(rt_state.log_R0)` from a
## release's `posterior_draws.csv`, the renewal walk's base: the first knot
## value, before the time-varying walk. Summarised the same way `rt_row`
## summarises the cut-off R_T, dated by the `cutoff` the release was built
## from. `nothing` when the draws are missing or carry no `log_R0` column, as
## pre-renewal releases and releases whose draws omit the walk base do.
function r0_row(tag, draws_path, cutoff)
    isnothing(draws_path) && return nothing
    date = string(cutoff)

    header, rows = read_simple_csv(draws_path)
    LOG_R0_COL in header || return nothing
    i = col(header, LOG_R0_COL)
    r0 = [exp(parse(Float64, r[i])) for r in rows]
    isempty(r0) && return nothing

    s = posterior_summary(r0)
    r3(x) = round(x; digits = 3)
    return (tag, date, r3(median(r0)), r3(s.lo30), r3(s.hi30),
        r3(s.lo60), r3(s.hi60), r3(s.lo90), r3(s.hi90))
end

# ----------------------------------------------------------------------
# Per-stream estimate summary
# ----------------------------------------------------------------------

## The two quantities a release's `stream_estimates.csv` carries, mapped to
## the per-release-by-stream table each is collected into.
const STREAM_QUANTITY_DEST = Dict(
    "R_T" => "rt_by_release_by_stream.csv",
    "C_T" => "size_by_release_by_stream.csv")

## The `(release, date, fit, median, lo30, ...)` rows from a release's
## `stream_estimates.csv`, one per `(fit, quantity)`, keyed by quantity so
## the caller writes each quantity to its own table. `date` is the release
## `cutoff`. Rows whose quantity is not one this script tabulates are
## ignored. `nothing` when the asset is absent or unreadable.
function stream_estimate_rows(tag, estimates_path, cutoff)
    isnothing(estimates_path) && return nothing
    header, rows = read_simple_csv(estimates_path)
    ("fit" in header && "quantity" in header) || return nothing
    fit_i = col(header, "fit")
    qty_i = col(header, "quantity")
    cols = Dict(c => col(header, c)
    for c in ("median", "lo30", "hi30", "lo60", "hi60", "lo90", "hi90"))

    byqty = Dict{String, Vector{NamedTuple}}()
    for r in rows
        qty = r[qty_i]
        haskey(STREAM_QUANTITY_DEST, qty) || continue
        vals = Dict(c => parse(Float64, r[i]) for (c, i) in cols)
        push!(get!(() -> NamedTuple[], byqty, qty),
            (; release = tag, date = string(cutoff), fit = r[fit_i],
                median = vals["median"], lo30 = vals["lo30"],
                hi30 = vals["hi30"], lo60 = vals["lo60"],
                hi60 = vals["hi60"], lo90 = vals["lo90"],
                hi90 = vals["hi90"]))
    end
    return byqty
end

# ----------------------------------------------------------------------
# Driver
# ----------------------------------------------------------------------

## Run only when invoked as a script, so a test can `include` this file for
## its functions (`score_release`, `truth_at`, …) without the driver's
## `gh` calls or file writes firing.
if abspath(PROGRAM_FILE) == @__FILE__
    entries = results_release_entries(repo)
    isempty(entries) && error("no releases found for $repo")

    ## The CURRENT observations manifest is the now-observed truth every
    ## release's forecast is scored against.
    obs = load_observations()
    grid_date(day) = obs.cutoff - Day(obs.n - day)

    score_rows = NamedTuple[]
    overlay_rows = NamedTuple[]
    ## Frozen-fit scores are kept in their own tables so the frozen made-dates
    ## never enter the cross-release baseline pool the main tables summarise.
    frozen_score_rows = NamedTuple[]
    frozen_overlay_rows = NamedTuple[]
    rt_rows = NamedTuple[]
    r0_rows = NamedTuple[]
    ## One accumulator per per-stream estimate table, keyed by its filename.
    stream_est_rows = Dict(
        f => NamedTuple[] for f in values(STREAM_QUANTITY_DEST))
    n_scored = 0
    n_no_forecast = 0
    n_backfilled = 0
    n_no_rt = 0
    n_frozen_scored = 0
    n_stopped = 0
    n_frozen_stopped = 0
    n_unstarted = 0
    n_frozen_unstarted = 0
    n_failed_reconstruction = 0

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
    ## assumed day, since its forecast could silently duplicate a kept one.
    ##
    ## The two reasons a cut-off is unreadable are logged apart, because they
    ## are not the same event. A release that carries no `as_of_date` is
    ## expected and permanent. One whose `observations.toml` could not be
    ## fetched (after `fetch_asset`'s retries) may be a transient `gh` failure,
    ## so a real distinct forecast could be dropped for this run; a re-run
    ## recovers it, but the log must let the two be told apart.
    dated = Tuple{String, DateTime, Date}[]
    no_cutoff = String[]
    unfetched = String[]
    for (tag, created) in entries
        obs_path = fetch_asset(repo, tag, OBS_ASSET, tagdir(tag))
        if isnothing(obs_path)
            push!(unfetched, tag)
            continue
        end
        cutoff = release_cutoff(obs_path)
        isnothing(cutoff) ? push!(no_cutoff, tag) :
        push!(dated, (tag, created, cutoff))
    end
    isempty(no_cutoff) || @info string(
        "ignoring ", length(no_cutoff), " release(s) that carry no data ",
        "cut-off: ", join(no_cutoff, ", "))
    isempty(unfetched) || @warn string(
        "could not fetch observations.toml for ", length(unfetched),
        " release(s); a distinct forecast may be dropped this run, a re-run ",
        "recovers it: ", join(unfetched, ", "))

    tags = select_daily_releases(dated)
    isempty(tags) && error("no results releases with a data cut-off for $repo")
    cutoffs = Dict(t => c for (t, _, c) in dated)
    println("Scoring $(length(tags)) release(s) of $(length(entries)), " *
            "one per data day.")

    for tag in tags
        ## A top-level `for` is a soft scope, so the running counters must be
        ## declared global to update the bindings above rather than shadow them.
        global n_scored, n_no_forecast, n_backfilled, n_no_rt, n_frozen_scored,
        n_stopped, n_frozen_stopped, n_unstarted, n_frozen_unstarted,
        n_failed_reconstruction

        ## The release's own `observations.toml` snapshot, already on disk
        ## from the selection pass above (`fetch_asset` is idempotent and
        ## `tagdir(tag)` names the same directory both times), so this is a
        ## local read, not a second network fetch. Used as the vintage
        ## manifest the persistence baseline is built from (see
        ## `vintage_observations`); `nothing` when the earlier fetch failed,
        ## in which case `score_release` falls back to the current `obs` and
        ## this release's baseline can leak a later revision — logged below
        ## rather than left a silent fallback, since this release was
        ## already selected to be scored.
        ##
        ## This snapshot is the release's OWN cut-off, not `made_date`'s: for
        ## the ordinary (non-frozen) scoring below `made_date` sits at or
        ## near the release's own cut-off, so the two are close together, but
        ## for the FROZEN scoring further down `made_date` is a fixed
        ## historical cut-off (e.g. 20 May) reused across many later release
        ## tags, so it can sit weeks before the snapshot used. `cutoff_date`
        ## in `vintage_observations`/`load_observations` only excludes
        ## vintages dated after `made_date`; it cannot un-revise a value that
        ## was already baked into the release-cutoff snapshot before
        ## `made_date`, if the correction landed between `made_date` and the
        ## release's own cut-off. This residual risk is not closed here — it
        ## would need a `made_date`-specific manifest snapshot per frozen
        ## cut-off, which is not archived — so the frozen baseline is honest
        ## about revisions that land after the RELEASE's cut-off but can
        ## still see one that lands between `made_date` and that cut-off.
        tag_obs_path = let p = joinpath(tagdir(tag), OBS_ASSET)
            isfile(p) ? p : nothing
        end
        isnothing(tag_obs_path) && @warn string(
            tag, ": no observations.toml on disk for this release; its ",
            "baseline falls back to the CURRENT manifest and can leak a ",
            "later revision or backfill into its persistence baseline")

        ## A saved forecast wins; a reconstruction stands in only where the
        ## release never stored one. The per-stream archive carries every fit's
        ## forecast, so it is preferred over the joint-only `forecast.csv`. A
        ## tag on `_FAILED_RECONSTRUCTIONS` is skipped before any of that: its
        ## forecast, wherever it comes from, is a failed reconstruction rather
        ## than a real one, so it is never worth fetching.
        label = tag
        forecast_path = if haskey(_FAILED_RECONSTRUCTIONS, tag)
            nothing
        else
            fp = fetch_asset(repo, tag, STREAM_FORECAST_ASSET, tagdir(tag))
            isnothing(fp) && (fp = fetch_asset(
                repo, tag, FORECAST_ASSET, tagdir(tag)))
            if isnothing(fp)
                asset = backfill_asset(tag)
                if !isnothing(asset) && asset in backfill_names
                    fp = fetch_asset(repo, BACKFILL_TAG, asset, tagdir(tag))
                    isnothing(fp) || (label = "$tag (backfill)")
                end
            end
            fp
        end

        if haskey(_FAILED_RECONSTRUCTIONS, tag)
            n_failed_reconstruction += 1
            @info string(tag, ": skipping known-failed reconstruction, not ",
                "scored as a model finding — ", _FAILED_RECONSTRUCTIONS[tag])
        elseif isnothing(forecast_path)
            n_no_forecast += 1
            @info "skipping $tag (no stored or reconstructed forecast)"
        else
            result = try
                score_release(label, forecast_path, obs, grid_date;
                    vintage_obs_path = tag_obs_path)
            catch e
                @warn "skipping $label forecast scoring" exception=e
                nothing
            end
            if !isnothing(result)
                append!(score_rows, result.rows)
                append!(overlay_rows, result.overlay)
                n_scored += 1
                label == tag || (n_backfilled += 1)
                n_stopped += result.stopped
                n_unstarted += result.unstarted
                result.skipped > 0 && @info string(
                    label, ": skipped ", result.skipped,
                    " not-yet-observed group(s)")
                result.stopped > 0 && @info string(
                    label, ": skipped ", result.stopped,
                    " group(s) past their stream's reporting coverage")
                result.unstarted > 0 && @info string(
                    label, ": skipped ", result.unstarted,
                    " group(s) before their stream began being reported")
            end
        end

        ## Frozen-fit forecasts: the current model re-fit at past cut-offs and
        ## forecast forward, each row stamped with its own frozen cut-off as
        ## the made date (see `forecast_frozen.csv` in the analysis page). It
        ## is scored the same way, tagged `frozen`, but into its own tables so
        ## the frozen made-dates never enter the cross-release baseline pool.
        ## Absent on releases built before the asset existed, a clean skip.
        frozen_path = fetch_asset(repo, tag, FROZEN_FORECAST_ASSET, tagdir(tag))
        if !isnothing(frozen_path)
            fresult = try
                score_release(tag, frozen_path, obs, grid_date;
                    default_fit = FROZEN_FIT, vintage_obs_path = tag_obs_path)
            catch e
                @warn "skipping $tag frozen forecast scoring" exception=e
                nothing
            end
            if !isnothing(fresult)
                append!(frozen_score_rows, fresult.rows)
                append!(frozen_overlay_rows, fresult.overlay)
                n_frozen_scored += 1
                n_frozen_stopped += fresult.stopped
                n_frozen_unstarted += fresult.unstarted
                fresult.skipped > 0 && @info string(
                    tag, " (frozen): skipped ", fresult.skipped,
                    " not-yet-observed group(s)")
                fresult.stopped > 0 && @info string(
                    tag, " (frozen): skipped ", fresult.stopped,
                    " group(s) past their stream's reporting coverage")
                fresult.unstarted > 0 && @info string(
                    tag, " (frozen): skipped ", fresult.unstarted,
                    " group(s) before their stream began being reported")
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

        ## R0 = exp(rt_state.log_R0) from the same draws asset, the renewal
        ## walk's base. Absent on pre-renewal releases and any whose draws
        ## omit the walk base, so a missing row is a clean skip.
        r0r = try
            r0_row(tag, draws_path, cutoffs[tag])
        catch e
            @warn "skipping $tag R0 summary" exception=e
            nothing
        end
        isnothing(r0r) || push!(r0_rows,
            (; release = r0r[1], date = r0r[2], median = r0r[3],
                lo30 = r0r[4], hi30 = r0r[5], lo60 = r0r[6], hi60 = r0r[7],
                lo90 = r0r[8], hi90 = r0r[9]))

        ## Per-fit R_T and C_T estimates, split into their own tables. Absent on
        ## every release until P2.2 lands, so a missing asset is a clean skip.
        est_path = fetch_asset(repo, tag, STREAM_ESTIMATES_ASSET, tagdir(tag))
        byqty = try
            stream_estimate_rows(tag, est_path, cutoffs[tag])
        catch e
            @warn "skipping $tag stream estimates" exception=e
            nothing
        end
        isnothing(byqty) || for (qty, rows) in byqty
            append!(stream_est_rows[STREAM_QUANTITY_DEST[qty]], rows)
        end
    end

    println("Scored $n_scored/$(length(tags)) releases " *
            "($n_no_forecast without a stored or reconstructed forecast, " *
            "$n_backfilled from $BACKFILL_TAG, $n_failed_reconstruction " *
            "excluded as known-failed reconstructions).")
    println("Dropped $n_stopped group(s) whose target ran past their " *
            "stream's own reporting coverage ($n_frozen_stopped in the " *
            "frozen tables).")
    println("Dropped $n_unstarted group(s) whose window opened before their " *
            "stream began being reported ($n_frozen_unstarted in the " *
            "frozen tables).")
    println("R_T summary for $(length(rt_rows))/$(length(tags)) releases " *
            "($n_no_rt without an R_T posterior).")

    ## `data/forecast_scores.csv`: one row per scored (release, made_date,
    ## stream, horizon, fit) group, `fit` the model or `baseline`. The two
    ## trailing column is the log-CRPS relative skill to the persistence
    ## baseline (baseline row `1.0`). The comparison against a stream's
    ## individual fit is computed at aggregation time by the table builders
    ## in `src/scoring.jl`, not stored per row here.
    scores = if isempty(score_rows)
        DataFrame(release = String[], made_date = Date[], stream = String[],
            horizon = Int[], target_date = Date[], fit = String[],
            crps = Float64[], log_crps = Float64[], dispersion = Float64[],
            overprediction = Float64[], underprediction = Float64[],
            coverage_50 = Float64[],
            coverage_90 = Float64[], bias = Float64[], n_samples = Int[])
    else
        sort(DataFrame(score_rows),
            [:release, :made_date, :stream, :horizon, :fit])
    end
    scores_base = rel_to_baseline_columns(scores)
    scores_dest = joinpath(@__DIR__, "..", "data", "forecast_scores.csv")
    write_simple_csv(scores_dest,
        [:release => scores.release,
            :made_date => string.(scores.made_date),
            :stream => scores.stream,
            :horizon => scores.horizon,
            :target_date => string.(scores.target_date),
            :fit => scores.fit,
            :crps => scores.crps,
            :log_crps => scores.log_crps,
            :dispersion => scores.dispersion,
            :overprediction => scores.overprediction,
            :underprediction => scores.underprediction,
            :coverage_50 => scores.coverage_50,
            :coverage_90 => scores.coverage_90,
            :bias => scores.bias,
            :n_samples => scores.n_samples,
            :log_rel_to_baseline => rel_to_baseline_cell.(scores_base)])
    println("Wrote $(nrow(scores)) scored forecasts to " *
            "data/forecast_scores.csv")

    ## `data/forecast_overlay.csv`: median and 30/60/90% bounds of each
    ## forecast group with the observed truth, for the overlay plots.
    overlay = if isempty(overlay_rows)
        DataFrame(release = String[], made_date = Date[], stream = String[],
            horizon = Int[], target_date = Date[], fit = String[],
            observed = Float64[], median = Float64[], lo30 = Float64[],
            hi30 = Float64[], lo60 = Float64[], hi60 = Float64[],
            lo90 = Float64[], hi90 = Float64[])
    else
        sort(DataFrame(overlay_rows),
            [:stream, :made_date, :horizon, :fit])
    end
    overlay_dest = joinpath(@__DIR__, "..", "data", "forecast_overlay.csv")
    write_simple_csv(overlay_dest,
        [:release => overlay.release,
            :made_date => string.(overlay.made_date),
            :stream => overlay.stream,
            :horizon => overlay.horizon,
            :target_date => string.(overlay.target_date),
            :fit => overlay.fit,
            :observed => overlay.observed,
            :median => overlay.median,
            :lo30 => overlay.lo30, :hi30 => overlay.hi30,
            :lo60 => overlay.lo60, :hi60 => overlay.hi60,
            :lo90 => overlay.lo90, :hi90 => overlay.hi90])
    println("Wrote $(nrow(overlay)) forecast-overlay rows to " *
            "data/forecast_overlay.csv")

    ## `data/forecast_scores_frozen.csv` and `data/forecast_overlay_frozen.csv`:
    ## the frozen-fit forecasts scored against the now-observed data, kept apart
    ## from the cross-release tables so their past made-dates never pool into
    ## the release baseline. Same schema as the release tables, with the model
    ## rows tagged `frozen`. Written even when empty (header only), since the
    ## docs build reads them and a missing file would throw.
    function _score_frame(rows)
        isempty(rows) && return DataFrame(release = String[],
            made_date = Date[], stream = String[], horizon = Int[],
            target_date = Date[], fit = String[], crps = Float64[],
            log_crps = Float64[], dispersion = Float64[],
            overprediction = Float64[], underprediction = Float64[],
            coverage_50 = Float64[],
            coverage_90 = Float64[], bias = Float64[], n_samples = Int[])
        return sort(DataFrame(rows),
            [:release, :made_date, :stream, :horizon, :fit])
    end
    function _overlay_frame(rows)
        isempty(rows) && return DataFrame(release = String[],
            made_date = Date[], stream = String[], horizon = Int[],
            target_date = Date[], fit = String[], observed = Float64[],
            median = Float64[], lo30 = Float64[], hi30 = Float64[],
            lo60 = Float64[], hi60 = Float64[], lo90 = Float64[],
            hi90 = Float64[])
        return sort(DataFrame(rows), [:stream, :made_date, :horizon, :fit])
    end

    ## The frozen scores get the same relative-skill column, computed
    ## against their OWN baseline within the frozen made-dates: the frozen
    ## fit's `log_crps` over its persistence baseline's.
    frozen_scores = _score_frame(frozen_score_rows)
    frozen_base = rel_to_baseline_columns(frozen_scores)
    write_simple_csv(
        joinpath(@__DIR__, "..", "data", "forecast_scores_frozen.csv"),
        [:release => frozen_scores.release,
            :made_date => string.(frozen_scores.made_date),
            :stream => frozen_scores.stream,
            :horizon => frozen_scores.horizon,
            :target_date => string.(frozen_scores.target_date),
            :fit => frozen_scores.fit,
            :crps => frozen_scores.crps,
            :log_crps => frozen_scores.log_crps,
            :dispersion => frozen_scores.dispersion,
            :overprediction => frozen_scores.overprediction,
            :underprediction => frozen_scores.underprediction,
            :coverage_50 => frozen_scores.coverage_50,
            :coverage_90 => frozen_scores.coverage_90,
            :bias => frozen_scores.bias,
            :n_samples => frozen_scores.n_samples,
            :log_rel_to_baseline => rel_to_baseline_cell.(frozen_base)])
    frozen_overlay = _overlay_frame(frozen_overlay_rows)
    write_simple_csv(
        joinpath(@__DIR__, "..", "data", "forecast_overlay_frozen.csv"),
        [:release => frozen_overlay.release,
            :made_date => string.(frozen_overlay.made_date),
            :stream => frozen_overlay.stream,
            :horizon => frozen_overlay.horizon,
            :target_date => string.(frozen_overlay.target_date),
            :fit => frozen_overlay.fit,
            :observed => frozen_overlay.observed,
            :median => frozen_overlay.median,
            :lo30 => frozen_overlay.lo30, :hi30 => frozen_overlay.hi30,
            :lo60 => frozen_overlay.lo60, :hi60 => frozen_overlay.hi60,
            :lo90 => frozen_overlay.lo90, :hi90 => frozen_overlay.hi90])
    println("Wrote $(nrow(frozen_scores)) frozen-fit scores over " *
            "$n_frozen_scored release(s) to data/forecast_scores_frozen.csv")

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

    ## `data/r0_by_release.csv`: the established initial reproduction number
    ## R0 = exp(rt_state.log_R0) per release, identical schema to
    ## `data/rt_by_release.csv`. Written even when empty, since the docs build
    ## reads it and a missing file throws.
    r0 = if isempty(r0_rows)
        DataFrame(release = String[], date = String[], median = Float64[],
            lo30 = Float64[], hi30 = Float64[], lo60 = Float64[],
            hi60 = Float64[], lo90 = Float64[], hi90 = Float64[])
    else
        sort(DataFrame(r0_rows), [:date])
    end
    r0_dest = joinpath(@__DIR__, "..", "data", "r0_by_release.csv")
    write_simple_csv(r0_dest,
        [:release => r0.release, :date => r0.date, :median => r0.median,
            :lo30 => r0.lo30, :hi30 => r0.hi30, :lo60 => r0.lo60,
            :hi60 => r0.hi60, :lo90 => r0.lo90, :hi90 => r0.hi90])
    println("Wrote $(nrow(r0)) release R0 summaries to " *
            "data/r0_by_release.csv")

    ## `data/rt_by_release_by_stream.csv` and
    ## `data/size_by_release_by_stream.csv`: the R_T and C_T posteriors of
    ## every fit, per release. Written even when empty, since the docs build
    ## reads them and a missing file throws rather than reading back empty.
    for file in sort(collect(values(STREAM_QUANTITY_DEST)))
        rows = stream_est_rows[file]
        df = if isempty(rows)
            DataFrame(release = String[], date = String[], fit = String[],
                median = Float64[], lo30 = Float64[], hi30 = Float64[],
                lo60 = Float64[], hi60 = Float64[], lo90 = Float64[],
                hi90 = Float64[])
        else
            sort(DataFrame(rows), [:date, :fit])
        end
        dest = joinpath(@__DIR__, "..", "data", file)
        write_simple_csv(dest,
            [:release => df.release, :date => df.date, :fit => df.fit,
                :median => df.median, :lo30 => df.lo30, :hi30 => df.hi30,
                :lo60 => df.lo60, :hi60 => df.hi60, :lo90 => df.lo90,
                :hi90 => df.hi90])
        println("Wrote $(nrow(df)) per-stream rows to data/$file")
    end
end
