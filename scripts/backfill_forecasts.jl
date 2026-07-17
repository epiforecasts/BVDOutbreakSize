#!/usr/bin/env julia
#
# Backfill the forecasts that past releases made but never stored.
#
# The report has always shown a one-week-ahead forecast at each release, but
# the forecast was never saved (only the posterior summaries and thinned draws
# were, and those do not carry the per-stream and reproduction-number-walk
# deterministics the forecast needs). To score past forecasts we reconstruct
# them honestly: check out the code as it was at each release, re-run its fit
# on that release's data, and export the forecast in the same archive schema
# the live build now writes (`forecast_archive` in `src/forecast.jl`).
#
# The reconstructions are collected under `output/backfill/` and are meant to
# be published as a dedicated `forecasts-backfill` release, which
# `scripts/score_releases.jl` then scores alongside the live release assets.
#
# This is compute-heavy: one full joint fit per release, and the fit cache is
# keyed on each tag's own source and data, so no cached fit carries over from
# the live build. Fits run through the cache under `logs/fit_cache`, so an
# interrupted run resumes, and a release whose archive already exists is
# skipped. Release tags are fit concurrently (`--concurrency`); the compute
# budget is checked before each fit starts.
#
# `forecast_archive` is newer than every release tag, so the driver cannot call
# it inside a release's worktree. It builds the archive columns itself, from
# the same stream/thinning conventions, and each tag's own `analysis.jl` call
# is reproduced by passing only the `obs_*` keywords that tag's
# `forecast_reported` declares.
#
# Coverage: only renewal-era releases (v1.4.0 on) carry the streams and the
# reproduction-number walk the forecast needs; integral-era releases are
# skipped. Releases whose code predates the fit registry
# (`docs/fits/registry.jl`) do not expose a single reusable fit entry point, so
# they are reported as needing a manual backfill rather than fit with a guessed
# call. In this repo that is v1.4.0-v1.6.0; extend the per-tag driver below to
# cover them if their forecasts are needed.
#
# The published reports only ever showed the 7-day horizon. The longer horizons
# reconstruct forecasts that tag's code would have made but never displayed.
#
# Usage:
#   julia scripts/backfill_forecasts.jl              # all renewal releases
#   julia scripts/backfill_forecasts.jl --only v1.8.0
#   julia scripts/backfill_forecasts.jl --keep       # keep the worktrees
#   julia scripts/backfill_forecasts.jl --concurrency 3

using Dates: Date

const REPO = "epiforecasts/BVDOutbreakSize"
const ROOT = normpath(joinpath(@__DIR__, ".."))
const OUT_DIR = joinpath(ROOT, "output", "backfill")
const WORKTREE_DIR = joinpath(ROOT, "worktrees")
## Shared across tags and kept outside the per-tag worktrees, so a fit survives
## `--keep` being off. Entries are content-addressed on each tag's own source
## and data, so sharing the directory cannot cross-contaminate tags.
const CACHE_DIR = joinpath(ROOT, "logs", "fit_cache")
const HORIZONS = (7, 14, 21, 28)
## Keep every fifth draw, matching `forecast_archive`'s use in the live build.
const THIN = 5
## v1.4.0 is the first renewal-model release (the `model` column of
## `data/released_estimates.csv`); everything earlier is integral-era.
const RENEWAL_FROM = v"1.4.0"
## Each fit samples with this many threads, so the default concurrency stays
## well inside a 16-core host. `--concurrency` overrides it.
const FIT_THREADS = 4
const DEFAULT_CONCURRENCY = 2

## Simple flag parsing: --keep, --only <tag>, --concurrency <n>.
function parse_args(args)
    keep = "--keep" in args
    only = nothing
    concurrency = DEFAULT_CONCURRENCY
    i = 1
    while i <= length(args)
        if args[i] == "--only"
            i < length(args) || error("--only needs a release tag")
            only = args[i + 1]
            i += 2
        elseif args[i] == "--concurrency"
            i < length(args) || error("--concurrency needs an integer")
            c = tryparse(Int, args[i + 1])
            isnothing(c) &&
                error("--concurrency must be an integer, got $(args[i + 1])")
            concurrency = c
            i += 2
        else
            i += 1
        end
    end
    concurrency >= 1 || error("--concurrency must be at least 1")
    return (; keep, only, concurrency)
end

## Renewal-era release tags to reconstruct, newest first, read from the tags
## themselves rather than from `data/released_estimates.csv` (that overlay
## lags the releases and would silently skip the newest ones).
function renewal_release_tags()
    out = readchomp(`git -C $ROOT tag --list "v*"`)
    tags = Tuple{VersionNumber, String}[]
    for line in split(out, '\n')
        tag = strip(line)
        isempty(tag) && continue
        v = tryparse(VersionNumber, lstrip(tag, 'v'))
        isnothing(v) && continue
        v >= RENEWAL_FROM && push!(tags, (v, String(tag)))
    end
    sort!(tags; by = first, rev = true)
    return [t for (_, t) in tags]
end

## Whether a tag's code carries the fit registry, the single reusable fit
## entry point the driver needs. Julia runs commands without a shell, so the
## `tag:path` revision is passed through literally.
function has_registry(code_tag)
    return success(pipeline(
        `git -C $ROOT cat-file -e $code_tag:docs/fits/registry.jl`;
        stderr = devnull))
end

## Whether the compute budget is clear to start another fit. Best-effort: when
## the guard script is absent or errors, assume it is fine (it fails open).
function compute_clear()
    guard = expanduser("~/.claude/hooks/compute-budget.sh")
    isfile(guard) || return true
    try
        return strip(read(`bash $guard`, String)) != "red"
    catch
        return true
    end
end

## Serialises the `Pkg.instantiate` step: every tag resolves into the same
## shared depot, and concurrent instantiate/precompile calls contend on the
## depot's package locks. Sampling itself still runs concurrently.
const INSTANTIATE_LOCK = ReentrantLock()

## Reconstruct one release's forecast in a worktree checked out at its code
## tag. Returns the path to the written forecast CSV, or `nothing` when the
## release cannot be reconstructed automatically (pre-registry code).
function backfill_one(code_tag; keep)
    dest = joinpath(OUT_DIR, "forecast_$(code_tag).csv")
    isfile(dest) && (@info "already done, skipping" code_tag; return dest)

    ## Pre-registry releases have no single reusable fit entry point.
    if !has_registry(code_tag)
        @warn "pre-registry release needs a manual backfill" code_tag
        return nothing
    end

    wt = joinpath(WORKTREE_DIR, "backfill-$(code_tag)")
    lock(INSTANTIATE_LOCK) do
        isdir(wt) || run(`git -C $ROOT worktree add --detach $wt $code_tag`)
        ## `registry.jl` includes `cache.jl`, which needs Serialization and
        ## SHA, and the archive is written with CSV; all three are docs-project
        ## dependencies, so the driver runs against `docs`, not the root
        ## project. Precompiling here, under the lock, keeps tags off each
        ## other's depot precompile locks and leaves the driver nothing to
        ## build, so it can load modules normally. Resolving with a stale image
        ## instead (`--compiled-modules=existing`) surfaces as an
        ## `UndefVarError` for a binding the loaded package really does define.
        run(`julia --project=$(joinpath(wt, "docs"))
             -e "using Pkg; Pkg.instantiate(); Pkg.precompile()"`)
    end
    try
        driver_path = joinpath(wt, "_backfill_driver.jl")
        write(driver_path, backfill_driver(dest))
        run(`julia --project=$(joinpath(wt, "docs"))
             --threads=$FIT_THREADS $driver_path`)
        rm(driver_path; force = true)
    finally
        keep || run(`git -C $ROOT worktree remove --force $wt`)
    end
    return isfile(dest) ? dest : nothing
end

## Source of the inline driver run inside a release's worktree. It mirrors the
## live build: build the headline joint fit spec from that tag's registry, run
## it through the content-addressed cache, forecast horizons 7-28 days and
## write the archive.
##
## `forecast_archive` postdates every release tag, so the driver reproduces it
## here: the incident streams as new counts over the horizon, isolation beds as
## an occupancy level, one row per stream/horizon/draw. Keep it in step with
## `forecast_archive` in `src/forecast.jl`.
function backfill_driver(dest)
    return """
    using Dates: Date, Day
    using DataFrames: DataFrame, propertynames, nrow
    using CSV: CSV
    using BVDOutbreakSize
    include(joinpath(pkgdir(BVDOutbreakSize), "docs", "fits", "registry.jl"))

    obs = load_observations()
    specs = build_fit_specs(obs)
    i = findfirst(s -> s.id == "joint", specs)
    isnothing(i) && error("this tag's registry has no \\"joint\\" fit; " *
                          "known ids: " * join((s.id for s in specs), ", "))
    chn = fit_or_load(fit_key("joint"), specs[i].thunk;
        cache_dir = "$(CACHE_DIR)")

    ## Pass only the observed-stream keywords this tag declares: the earlier
    ## renewal releases predate `obs_recovered`, so passing it unconditionally
    ## would fail on them.
    accepted = Base.kwarg_decl(first(methods(forecast_reported)))
    kw = Dict{Symbol, Any}(
        :obs_cases => obs.reported_cases, :obs_deaths => obs.total_deaths)
    for (kwname, field) in (
        :obs_confirmed => :confirmed_cases,
        :obs_confirmed_deaths => :confirmed_deaths,
        :obs_recovered => :recovered_cases)
        kwname in accepted && hasproperty(obs, field) &&
            (kw[kwname] = getproperty(obs, field))
    end
    missed = setdiff(
        (:obs_confirmed, :obs_confirmed_deaths, :obs_recovered), keys(kw))
    isempty(missed) ||
        @warn "tag does not take some observed-stream keywords" missed

    runs = [(h, forecast_reported(chn; horizon = h, kw...))
            for h in $(HORIZONS)]

    ## Mirror of `forecast_archive`: incident (new-over-horizon) and level
    ## quantities only, never cumulative.
    streams = (
        (:confirmed_new, "confirmed cases"),
        (:confirmed_deaths_new, "confirmed deaths"),
        (:recovered_new, "recovered"),
        (:isolation_level, "isolation beds"))
    out = DataFrame(made_date = Date[], horizon = Int[], target_date = Date[],
        stream = String[], draw = Int[], value = Float64[])
    for (horizon, fc) in runs
        h = Int(horizon)
        target = obs.cutoff + Day(h)
        for (col, label) in streams
            col in propertynames(fc) || continue
            vals = fc[!, col]
            for (d, i) in enumerate(1:$(THIN):length(vals))
                push!(out,
                    (obs.cutoff, h, target, label, d, Float64(vals[i])))
            end
        end
    end
    isempty(out) && error("no archived streams in the forecast")
    mkpath(dirname("$(dest)"))
    CSV.write("$(dest)", out)
    @info "wrote archive" dest="$(dest)" rows=nrow(out) cutoff=obs.cutoff
    """
end

# ----------------------------------------------------------------------
# Driver
# ----------------------------------------------------------------------

## Backfill every tag in `tags`, honouring `concurrency` fits in flight, and
## return the tags split into automatically reconstructed, still-manual and
## failed. Each fit samples on `FIT_THREADS`, so the default concurrency stays
## inside the host.
function run_backfill(tags; keep, concurrency)
    done = String[]
    manual = String[]
    failed = String[]
    results = ReentrantLock()
    gate = Base.Semaphore(concurrency)
    @sync for code_tag in tags
        @async begin
            Base.acquire(gate)
            try
                while !compute_clear()
                    @info "compute budget red; waiting 60s" code_tag
                    sleep(60)
                end
                path = try
                    backfill_one(code_tag; keep)
                catch e
                    @warn "backfill failed" code_tag exception=e
                    :failed
                end
                lock(results) do
                    path === :failed ? push!(failed, code_tag) :
                    isnothing(path) ? push!(manual, code_tag) :
                    push!(done, code_tag)
                end
            finally
                Base.release(gate)
            end
        end
    end
    sort!(done; rev = true)
    sort!(manual; rev = true)
    sort!(failed; rev = true)
    return (; done, manual, failed)
end

## Concatenate the per-release archives named in `done` into one combined file.
function combine_archives(done)
    combined = joinpath(OUT_DIR, "forecasts_backfill.csv")
    open(combined, "w") do io
        wrote_header = false
        for code_tag in done
            f = joinpath(OUT_DIR, "forecast_$(code_tag).csv")
            isfile(f) || continue
            lines = readlines(f)
            isempty(lines) && continue
            wrote_header || (println(io, lines[1]); wrote_header = true)
            for l in lines[2:end]
                isempty(strip(l)) || println(io, l)
            end
        end
    end
    return combined
end

function main(args = ARGS)
    opts = parse_args(args)
    mkpath(OUT_DIR)
    mkpath(CACHE_DIR)

    all_tags = renewal_release_tags()
    tags = isnothing(opts.only) ? all_tags : filter(==(opts.only), all_tags)
    isempty(tags) && error("no matching renewal release tags (have $all_tags)")

    @info "backfilling forecasts" tags concurrency=opts.concurrency
    r = run_backfill(tags; keep = opts.keep, concurrency = opts.concurrency)
    combined = combine_archives(r.done)

    println()
    println("Backfilled $(length(r.done)) release(s): ", join(r.done, ", "))
    isempty(r.manual) ||
        println("Manual backfill still needed for: ", join(r.manual, ", "))
    isempty(r.failed) ||
        println("Failed (see the log above): ", join(r.failed, ", "))
    println("Combined archive: $combined")
    println()
    println("Review the archives, then publish them as a dedicated release " *
            "with:")
    println()
    println("  gh release create forecasts-backfill -R $REPO \\")
    println("    --title 'Backfilled historical forecasts' \\")
    println("    --notes 'Forecasts reconstructed by re-running each past " *
            "release'\\''s own code on its data, because forecasts were " *
            "shown in the report but never stored. Scored by " *
            "scripts/score_releases.jl.' \\")
    println("    $OUT_DIR/forecasts_backfill.csv $OUT_DIR/forecast_v*.csv")
    return r
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
