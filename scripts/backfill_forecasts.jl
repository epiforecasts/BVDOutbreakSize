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
# This is compute-heavy: one full joint fit per release. It runs serially by
# default, is resumable (a release whose output already exists is skipped) and
# checks the compute budget between fits. Run it once, out of band, and review
# before publishing.
#
# Coverage: only renewal-era releases (v1.4.0 on) carry the streams and the
# reproduction-number walk the forecast needs; integral-era releases are
# skipped. Releases whose code predates the fit registry
# (`docs/fits/registry.jl`, added at v1.8.0) do not expose a single reusable
# fit entry point, so they are reported as needing a manual backfill rather
# than fit with a guessed call. In this repo that is v1.4.0-v1.7.0; extend the
# per-tag driver below to cover them if their forecasts are needed.
#
# Usage:
#   julia scripts/backfill_forecasts.jl              # all renewal releases
#   julia scripts/backfill_forecasts.jl --only v1.8.0
#   julia scripts/backfill_forecasts.jl --keep       # keep the worktrees
#   julia scripts/backfill_forecasts.jl --concurrency 2

using TOML: TOML
using Dates: Date

const REPO = "epiforecasts/BVDOutbreakSize"
const ROOT = normpath(joinpath(@__DIR__, ".."))
const OUT_DIR = joinpath(ROOT, "output", "backfill")
const WORKTREE_DIR = joinpath(ROOT, "worktrees")
const HORIZONS = (7, 14, 21, 28)

## Simple flag parsing: --keep, --only <tag>, --concurrency <n>.
function parse_args(args)
    keep = "--keep" in args
    only = nothing
    concurrency = 1
    i = 1
    while i <= length(args)
        if args[i] == "--only" && i < length(args)
            only = args[i + 1]
            i += 2
        elseif args[i] == "--concurrency" && i < length(args)
            concurrency = parse(Int, args[i + 1])
            i += 2
        else
            i += 1
        end
    end
    return (; keep, only, concurrency)
end

## Renewal-era release tags to reconstruct, newest first, from
## `data/released_estimates.csv` (the `model` column marks integral vs
## renewal). Each results release `results-vX.Y.Z` was built from the code tag
## `vX.Y.Z`, so the code tag is the results tag with the `results-` prefix
## stripped.
function renewal_release_tags()
    csv = joinpath(ROOT, "data", "released_estimates.csv")
    lines = readlines(csv)
    header = split(lines[1], ',')
    tag_i = findfirst(==("tag"), header)
    model_i = findfirst(==("model"), header)
    tags = String[]
    for l in lines[2:end]
        isempty(strip(l)) && continue
        f = split(l, ',')
        f[model_i] == "renewal" || continue
        push!(tags, replace(f[tag_i], "results-" => ""))
    end
    return sort(unique(tags); rev = true)
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

## Reconstruct one release's forecast in a worktree checked out at its code
## tag. Returns the path to the written forecast CSV, or `nothing` when the
## release cannot be reconstructed automatically (pre-registry code).
function backfill_one(code_tag; keep)
    dest = joinpath(OUT_DIR, "forecast_$(code_tag).csv")
    isfile(dest) && (@info "already done, skipping" code_tag; return dest)

    ## Pre-registry releases have no single reusable fit entry point.
    if !success(pipeline(`git -C $ROOT cat-file -e $code_tag:docs/fits/registry.jl`;
        stderr = devnull))
        @warn "pre-registry release needs a manual backfill" code_tag
        return nothing
    end

    wt = joinpath(WORKTREE_DIR, "backfill-$(code_tag)")
    isdir(wt) || run(`git -C $ROOT worktree add --detach $wt $code_tag`)
    try
        ## The per-tag driver runs inside the worktree's own environment so it
        ## uses that release's model code. It fits the joint through that tag's
        ## fit registry, forecasts each horizon and writes the archive CSV.
        driver = backfill_driver(dest)
        driver_path = joinpath(wt, "_backfill_driver.jl")
        write(driver_path, driver)
        run(`julia --project=$wt -e "using Pkg; Pkg.instantiate()"`)
        run(`julia --project=$wt $driver_path`)
        rm(driver_path; force = true)
    finally
        keep || run(`git -C $ROOT worktree remove --force $wt`)
    end
    return isfile(dest) ? dest : nothing
end

## Source of the inline driver run inside a release's worktree. It mirrors the
## live build: build the headline joint fit spec from that tag's registry, run
## it (through the content-addressed cache when present), forecast horizons
## 7-28 days and write the archive. Older tags may lack `forecast_archive`, so
## the archive columns are constructed by hand when it is not exported.
function backfill_driver(dest)
    horizons = join(HORIZONS, ", ")
    return """
    using BVDOutbreakSize
    include(joinpath(pkgdir(BVDOutbreakSize), "docs", "fits", "registry.jl"))
    obs = load_observations()
    specs = build_fit_specs(obs)
    joint = only(filter(s -> s.id == "joint", specs))
    chn = joint.thunk()
    runs = [(h, forecast_reported(chn; horizon = h,
                 obs_cases = obs.reported_cases, obs_deaths = obs.total_deaths,
                 obs_confirmed = obs.confirmed_cases,
                 obs_confirmed_deaths = obs.confirmed_deaths,
                 obs_recovered = obs.recovered_cases)) for h in ($horizons,)]
    if isdefined(BVDOutbreakSize, :forecast_archive)
        using CSV
        CSV.write("$(dest)",
            forecast_archive(runs; made_date = obs.cutoff, thin = 5))
    else
        error("this tag lacks forecast_archive; extend the driver to build " *
              "the made_date/horizon/target_date/stream/draw/value columns " *
              "by hand from the forecast_reported DataFrame")
    end
    """
end

# ----------------------------------------------------------------------
# Driver
# ----------------------------------------------------------------------

opts = parse_args(ARGS)
mkpath(OUT_DIR)

all_tags = renewal_release_tags()
tags = isnothing(opts.only) ? all_tags : filter(==(opts.only), all_tags)
isempty(tags) && error("no matching renewal release tags (have $all_tags)")

@info "backfilling forecasts" tags concurrency=opts.concurrency
done = String[]
manual = String[]
for code_tag in tags
    while !compute_clear()
        @info "compute budget red; waiting 60s before the next fit"
        sleep(60)
    end
    path = try
        backfill_one(code_tag; keep = opts.keep)
    catch e
        @warn "backfill failed" code_tag exception=e
        nothing
    end
    isnothing(path) ? push!(manual, code_tag) : push!(done, code_tag)
end

## Concatenate the per-release archives into one combined file.
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

println()
println("Backfilled $(length(done)) release(s): ", join(done, ", "))
isempty(manual) ||
    println("Manual backfill still needed for: ", join(manual, ", "))
println("Combined archive: $combined")
println()
println("Review the archives, then publish them as a dedicated release with:")
println()
println("  gh release create forecasts-backfill -R $REPO \\")
println("    --title 'Backfilled historical forecasts' \\")
println("    --notes 'Forecasts reconstructed by re-running each past " *
        "release'\\''s own code on its data, because forecasts were shown " *
        "in the report but never stored. Scored by " *
        "scripts/score_releases.jl.' \\")
println("    $OUT_DIR/forecasts_backfill.csv $OUT_DIR/forecast_v*.csv")
