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
# Coverage: renewal-era releases (v1.4.0 on) carry the streams and the
# reproduction-number walk the forecast needs. Registry-era releases (v1.7.0
# on) reach the fit through the fit registry (`docs/fits/registry.jl`). The
# earlier renewal releases (v1.4.0-v1.6.0) predate the registry, so a per-tag
# pre-registry driver reproduces each tag's own headline joint call verbatim
# instead. A renewal tag that is neither registry-era nor one of those recorded
# pre-registry tags is reported as needing a manual backfill.
#
# The integral-era releases (before v1.4.0) predate the renewal model and the
# reproduction-number walk, so each needs a driver reproducing its own joint
# call and forecast signature. The tags listed in `INTEGRAL_TAGS` are added to
# the run and dispatched to `integral_driver`; only v1.3.0 is listed so far (its
# `forecast_reported` is the only integral one that emits confirmed streams).
# The earlier integral releases (v1.0.0-v1.2.0) forecast the reported/suspected
# streams only and extend `INTEGRAL_TAGS` as their drivers land.
#
# The published reports only ever showed the 7-day horizon. The longer horizons
# reconstruct forecasts that tag's code would have made but never displayed.
#
# Usage:
#   julia scripts/backfill_forecasts.jl              # all release tags
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
## Renewal releases whose code predates the fit registry. `preregistry_driver`
## records each tag's own headline joint call, so it backfills exactly these.
const PREREGISTRY_TAGS = ("v1.4.0", "v1.5.0", "v1.6.0")
## Integral-era releases backfilled by `integral_driver`. These predate the
## renewal model, so `renewal_release_tags` never reaches them; they are added
## to the run explicitly. v1.3.0 is the only integral release whose
## `forecast_reported` emits confirmed streams (v1.0.0-v1.2.0 forecast the
## reported/suspected streams only and are added here as their drivers land).
## Extending this tuple is the single place to grow integral coverage.
const INTEGRAL_TAGS = ("v1.3.0",)
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
## tag. Returns the path to the written forecast CSV, or `nothing` when no
## driver is recorded for the tag (a renewal tag that is neither registry-era
## nor a recorded pre-registry tag).
function backfill_one(code_tag; keep)
    dest = joinpath(OUT_DIR, "forecast_$(code_tag).csv")
    isfile(dest) && (@info "already done, skipping" code_tag; return dest)

    ## Registry-era tags reach the fit through the registry; the earlier
    ## renewal tags predate it and use the per-tag pre-registry driver; the
    ## integral-era tags predate the renewal model and use the integral driver.
    driver_src = if has_registry(code_tag)
        backfill_driver(dest)
    elseif code_tag in PREREGISTRY_TAGS
        preregistry_driver(dest, code_tag)
    elseif code_tag in INTEGRAL_TAGS
        integral_driver(dest, code_tag)
    else
        @warn "no driver recorded; needs a manual backfill" code_tag
        return nothing
    end

    wt = joinpath(WORKTREE_DIR, "backfill-$(code_tag)")
    lock(INSTANTIATE_LOCK) do
        isdir(wt) || run(`git -C $ROOT worktree add --detach $wt $code_tag`)
        ## Both drivers run against the `docs` project, not the root project:
        ## the registry-era one includes `registry.jl` (which needs
        ## Serialization and SHA through `cache.jl`), the pre-registry one saves
        ## its chain with Serialization, and both write the archive with CSV,
        ## all docs-project dependencies. Precompiling here, under the lock,
        ## keeps tags off each other's depot precompile locks and leaves the
        ## driver nothing to build, so it can load modules normally. Resolving
        ## with a stale image instead (`--compiled-modules=existing`) surfaces
        ## as an `UndefVarError` for a binding the loaded package really does
        ## define.
        run(`julia --project=$(joinpath(wt, "docs"))
             -e "using Pkg; Pkg.instantiate(); Pkg.precompile()"`)
    end
    try
        driver_path = joinpath(wt, "_backfill_driver.jl")
        write(driver_path, driver_src)
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
    using Serialization: serialize
    using BVDOutbreakSize
    include(joinpath(pkgdir(BVDOutbreakSize), "docs", "fits", "registry.jl"))

    obs = load_observations()
    specs = build_fit_specs(obs)
    i = findfirst(s -> s.id == "joint", specs)
    isnothing(i) && error("this tag's registry has no \\"joint\\" fit; " *
                          "known ids: " * join((s.id for s in specs), ", "))
    chn = fit_or_load(fit_key("joint"), specs[i].thunk;
        cache_dir = "$(CACHE_DIR)")

    ## Keep the chain alongside the archive: `fit_or_load` already caches the
    ## fit under the content-addressed cache, but serialising it next to the
    ## archive lets a later re-archive (e.g. to add a stream) reuse the exact
    ## chain without re-resolving the cache or refitting.
    chain_path = "$(dest).chain"
    try
        serialize(chain_path, chn)
        @info "chain saved" chain_path
    catch e
        @warn "could not save chain" exception = e
    end

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

    ## Mirror of `forecast_archive`: the confirmed/recovered/isolation streams
    ## plus the reported-case and suspected-death incident streams the scorer
    ## maps (`STREAM_HISTORY` in `scripts/score_releases.jl`). Incident
    ## (new-over-horizon) and level quantities only, never the vintage-revised
    ## cumulative totals. Absent columns are skipped by the `propertynames`
    ## guard below, so a tag that lacks a stream carries fewer rows.
    streams = (
        (:confirmed_new, "confirmed cases"),
        (:confirmed_deaths_new, "confirmed deaths"),
        (:recovered_new, "recovered"),
        (:isolation_level, "isolation beds"),
        (:cases_new, "reported cases"),
        (:deaths_new, "suspected deaths"))
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

## The headline joint `nuts_sample` call as `code_tag`'s own `analysis.jl`
## writes it, ready to embed in the pre-registry driver.
##
## v1.4.0/v1.5.0 take neither the isolation/bed streams nor recovered: those
## columns do not exist in `src/data.jl` at those tags, so the model cannot
## ingest them and the forecast cannot emit them. v1.6.0 adds the
## daily-suspected, isolation, bed-capacity and recovered histories, which is
## what puts the bed and recovered deterministics in the chain that
## `forecast_reported` probes before emitting `isolation_level` and
## `recovered_new`.
##
## `target_accept` is not restated, so each tag's own `nuts_sample` default
## applies (v1.4.0/v1.5.0 default to 0.9, v1.6.0 to 0.85); restating one value
## would misfit the others. `samples`/`chains` are passed explicitly and
## default to the production setting.
function preregistry_joint_call(code_tag, samples, chains)
    if code_tag == "v1.4.0"
        return """
        nuts_sample(
            bvd_joint(
                obs.n, obs.exported_cases, obs.total_deaths,
                obs.reported_cases, obs.exports_deaths, obs.confirmed_cases,
                obs.tests_analysed;
                confirmed_deaths = obs.confirmed_deaths,
                deaths_history = obs.deaths_history,
                reported_history = obs.reported_history,
                confirmed_history = obs.confirmed_history,
                confirmed_deaths_history = obs.confirmed_deaths_history,
                lab_history = obs.lab_history,
                tests_received_history = obs.tests_received_history,
                export_case_days = obs.export_case_days,
                export_death_days = obs.export_death_days,
                breakpoint = bp,
                background_re = true,
                confirmed_positivity_link = :composition,
                genetic = genetic_seeding_model,
                tmrca_days = obs.tmrca_days);
            samples = $(samples), chains = $(chains))"""
    elseif code_tag == "v1.5.0"
        ## v1.5.0 dropped `tests_received_history` for the daily lab and
        ## suspected-case histories; its `bvd_joint` rejects the v1.4.0 kwargs.
        return """
        nuts_sample(
            bvd_joint(
                obs.n, obs.exported_cases, obs.total_deaths,
                obs.reported_cases, obs.exports_deaths, obs.confirmed_cases,
                obs.tests_analysed;
                confirmed_deaths = obs.confirmed_deaths,
                deaths_history = obs.deaths_history,
                reported_history = obs.reported_history,
                confirmed_history = obs.confirmed_history,
                confirmed_deaths_history = obs.confirmed_deaths_history,
                lab_history = obs.lab_history,
                lab_daily_history = obs.lab_daily_history,
                suspected_daily_history = obs.suspected_daily_history,
                export_case_days = obs.export_case_days,
                export_death_days = obs.export_death_days,
                breakpoint = bp,
                background_re = true,
                confirmed_positivity_link = :composition,
                genetic = genetic_seeding_model,
                tmrca_days = obs.tmrca_days);
            samples = $(samples), chains = $(chains))"""
    elseif code_tag == "v1.6.0"
        return """
        nuts_sample(
            bvd_joint(
                obs.n, obs.exported_cases, obs.total_deaths,
                obs.reported_cases, obs.exports_deaths, obs.confirmed_cases,
                obs.tests_analysed;
                confirmed_deaths = obs.confirmed_deaths,
                recovered_cases = obs.recovered_cases,
                deaths_history = obs.deaths_history,
                reported_history = obs.reported_history,
                confirmed_history = obs.confirmed_history,
                confirmed_deaths_history = obs.confirmed_deaths_history,
                lab_history = obs.lab_history,
                lab_daily_history = obs.lab_daily_history,
                suspected_daily_history = obs.suspected_daily_history,
                suspected_daily_deaths_history =
                    obs.suspected_daily_deaths_history,
                isolation_history = obs.isolation_history,
                bed_capacity_history = obs.bed_capacity_history,
                recovered_history = obs.recovered_history,
                export_case_days = obs.export_case_days,
                export_death_days = obs.export_death_days,
                breakpoint = bp,
                background_re = true,
                confirmed_positivity_link = :composition,
                genetic = genetic_seeding_model,
                tmrca_days = obs.tmrca_days);
            samples = $(samples), chains = $(chains))"""
    end
    error("no joint call recorded for $(code_tag)")
end

## Source of the inline driver run inside a pre-registry release's worktree
## (v1.4.0-v1.6.0). These tags predate the fit registry, so the driver
## reproduces each tag's own headline joint call verbatim (see
## `preregistry_joint_call`) rather than reaching the fit through the registry.
## The chain is serialised alongside the archive because these tags have no
## `fit_or_load`, so a failure in the archive step would otherwise cost a whole
## refit. Everything downstream of the fit matches `backfill_driver`, so the
## archive is schema-identical (made_date, horizon, target_date, stream, draw,
## value).
function preregistry_driver(dest, code_tag; samples = 1000, chains = 2)
    joint = preregistry_joint_call(code_tag, samples, chains)
    return """
    using Dates: Date, Day
    using DataFrames: DataFrame, propertynames, nrow
    using CSV: CSV
    using Serialization: serialize
    using BVDOutbreakSize

    obs = load_observations()
    @info "fitting" tag="$(code_tag)" cutoff=obs.cutoff
    bp = obs.n - obs.who_first_sitrep_days
    chn = $(joint)

    ## Keep the chain: these tags have no `fit_or_load`, so a failure in the
    ## archive step below would otherwise cost a whole refit.
    chain_path = "$(dest).chain"
    try
        serialize(chain_path, chn)
        @info "chain saved" chain_path
    catch e
        @warn "could not save chain" exception = e
    end

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

    ## Mirror of `forecast_archive`: the confirmed/recovered/isolation streams
    ## plus the reported-case and suspected-death incident streams the scorer
    ## maps (`STREAM_HISTORY` in `scripts/score_releases.jl`). Incident
    ## (new-over-horizon) and level quantities only, never the vintage-revised
    ## cumulative totals. Absent columns are skipped by the `propertynames`
    ## guard below, so a tag that lacks a stream carries fewer rows.
    streams = (
        (:confirmed_new, "confirmed cases"),
        (:confirmed_deaths_new, "confirmed deaths"),
        (:recovered_new, "recovered"),
        (:isolation_level, "isolation beds"),
        (:cases_new, "reported cases"),
        (:deaths_new, "suspected deaths"))
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

## Source of the inline driver run inside an integral-era release's worktree.
## Currently only v1.3.0, the one integral release whose `forecast_reported`
## emits confirmed streams, so it archives the two confirmed streams
## ("confirmed cases", "confirmed deaths"); recovered and isolation did not
## exist yet. `joint_obs`, `_increments` and `growth_for` live in v1.3.0's
## analysis.jl (not the package), so they are reproduced here verbatim;
## everything else is package-exported. The joint call (analysis.jl:2018) and
## forecast call (analysis.jl:2741) are reproduced exactly, including the
## integral `forecast_reported` signature (daily_travellers/source_population/
## obs_confirmed/obs_confirmed_deaths/obs_analysed). The chain is serialised
## alongside the archive because this tag has no `fit_or_load`. Downstream of
## the forecast the archive matches `backfill_driver` and `preregistry_driver`
## (made_date, horizon, target_date, stream, draw, value), so
## `score_releases.jl` scores it unchanged. v1.3.0 obs has no `cutoff` field,
## so the cut-off date is `Date(obs.as_of_date)`. The earlier integral tags
## (v1.0.0-v1.2.0) reproduce a different joint call and forecast signature and
## are added as their drivers land; the guarded stream superset below already
## covers the reported/suspected/export columns they emit.
function integral_driver(dest, code_tag; samples = 1000, chains = 2)
    return """
    using Dates: Date, Day
    using DataFrames: DataFrame, propertynames, nrow
    using CSV: CSV
    using Serialization: serialize
    using Turing
    using Distributions
    using BVDOutbreakSize

    ## Between-vintage increments of a cumulative sitrep history (the first
    ## increment is the cumulative at the first vintage).
    function _increments(v)
        d = similar(v, Int)
        prev = 0
        for i in eachindex(v)
            d[i] = v[i] - prev
            prev = v[i]
        end
        return d
    end

    ## Per-vintage observation arguments for the joint model, reproduced
    ## from v1.3.0 analysis.jl. Each DRC stream is fitted as between-vintage
    ## increments of its cumulative sitrep history.
    function joint_obs(o; observe = true)
        _stream(h,
            s) = h === missing ?
                 (Union{Missing, Int}[observe ? s : missing], [0]) :
                 (observe ? _increments(h.values) :
                  fill(missing, length(h.values)), h.offsets)
        rep, rep_off = _stream(o.reported_case_history, o.reported_cases)
        dth, dth_off = _stream(o.death_history, o.total_deaths)
        have_conf = o.confirmed_case_history !== missing ||
                    o.confirmed_cases !== missing
        have_pervintage = have_conf &&
                          o.confirmed_case_history !== missing &&
                          o.tests_analysed_history !== missing
        if have_pervintage
            sa = o.tests_analysed_history
            sr = o.tests_received_history
            ch = o.confirmed_case_history
            conf_off = collect(ch.offsets)
            conf = observe ? Union{Missing, Int}[_increments(ch.values)...] :
                   fill(missing, length(conf_off))
            keep = [i == 1 || sa.values[i] > sa.values[i - 1]
                    for i in eachindex(sa.values)]
            aoff = collect(sa.offsets)[keep]
            a_inc = _increments(sa.values[keep])
            ridx = [findfirst(==(off), sr.offsets) for off in aoff]
            r_inc = _increments([sr.values[i] for i in ridx])
            aobs = Dict(off => a_inc[k] for (k, off) in enumerate(aoff))
            robs = Dict(off => r_inc[k] for (k, off) in enumerate(aoff))
            analysed = Union{Missing, Int}[get(aobs, off, missing)
                                           for off in conf_off]
            received = Union{Missing, Int}[get(robs, off, missing)
                                           for off in conf_off]
        else
            conf_total = o.confirmed_cases !== missing ? o.confirmed_cases :
                         o.confirmed_case_history === missing ? missing :
                         o.confirmed_case_history.values[end]
            conf,
            conf_off = have_conf ?
                       (Union{Missing, Int}[observe ? conf_total : missing],
                [0]) : (Union{Missing, Int}[], Int[])
            analysed = Union{Missing, Int}[]
            received = Union{Missing, Int}[]
        end
        ec_full = o.exported_cases_daily
        last_import = isempty(ec_full) ? nothing : findlast(!=(0), ec_full)
        export_last_offset = last_import === nothing ? 0 :
                             length(ec_full) - last_import
        _truncate(v) = v[1:max(length(v) - export_last_offset, 0)]
        ecases = isempty(ec_full) ? ec_full :
                 (observe ? _truncate(ec_full) :
                  fill(missing, length(_truncate(ec_full))))
        ed_trunc = _truncate(o.export_deaths_daily)
        edaily = observe ? ed_trunc : fill(missing, length(ed_trunc))
        if o.confirmed_death_history !== missing
            cdh = o.confirmed_death_history
            cdeath = observe ?
                     Union{Missing, Int}[_increments(cdh.values)...] :
                     fill(missing, length(cdh.values))
            cdeath_off = collect(cdh.offsets)
        else
            cdeath = Union{Missing, Int}[]
            cdeath_off = Int[]
        end
        return (deaths = dth, reported = rep, export_deaths = edaily,
            kw = (; reported_offsets = rep_off, death_offsets = dth_off,
                confirmed_cases = conf, confirmed_offsets = conf_off,
                samples_analysed = analysed,
                samples_received = received,
                confirmed_deaths = cdeath,
                confirmed_death_offsets = cdeath_off,
                exported_cases_daily = ecases,
                export_last_offset = export_last_offset,
                confirmed_epi_exclusion = nothing,
                tests_analysed = observe ? o.cumulative_tests_analysed :
                                 missing, tests_offset = 0))
    end

    ## Growth submodel whose doubling-count prior centre advances with the
    ## cut-off date (analysis.jl:2004).
    function growth_for(o)
        exponential_growth_model(
            m_prior = truncated(Normal(m_prior_centre(o.as_of_date), 3.0);
                lower = 0))
    end

    obs = load_observations()
    cutoff = Date(obs.as_of_date)
    @info "fitting" tag="$(code_tag)" cutoff

    fit_args = joint_obs(obs)
    growth_now = growth_for(obs)
    genetic_seeding = T -> genetic_seeding_model(T, obs.genetic_tmrca_days;
        tmrca_days_sd = obs.genetic_tmrca_days_sd)

    ## Headline joint call, verbatim from v1.3.0 analysis.jl:2018.
    chn = nuts_sample(
        bvd_joint(obs.exported_cases, fit_args.deaths, fit_args.reported,
            fit_args.export_deaths; fit_args.kw...,
            growth = growth_now,
            first_export_detection_delta = obs.first_export_detection_delta,
            report_onset_offset = report_onset_offset(obs.as_of_date),
            confirmed_q_random_effect = confirmed_q_re_model,
            genetic = genetic_seeding);
        samples = $(samples), chains = $(chains), target_accept = 0.9)

    ## Keep the chain: this tag has no `fit_or_load`, so a failure in the
    ## forecast/archive step below would otherwise cost a whole refit.
    chain_path = "$(dest).chain"
    mkpath(dirname("$(dest)"))
    try
        serialize(chain_path, chn)
        @info "chain saved" chain_path
    catch e
        @warn "could not save chain" exception = e
    end

    ## Forecast call, verbatim from v1.3.0 analysis.jl:2741. The integral
    ## `forecast_reported` needs the traveller/population constants and the
    ## confirmed-queue anchors even though exports are dropped.
    runs = [(h, forecast_reported(chn;
                 horizon = h,
                 daily_travellers = obs.daily_outbound_travellers,
                 source_population = obs.source_population,
                 obs_confirmed = obs.confirmed_cases,
                 obs_confirmed_deaths = obs.confirmed_deaths,
                 obs_analysed = obs.cumulative_tests_analysed,
                 forecast_exports = false,
                 report_onset_offset = report_onset_offset(obs.as_of_date)))
            for h in $(HORIZONS)]

    ## Scoreable integral streams as a guarded superset: v1.3.0 emits the two
    ## confirmed streams (recovered and isolation did not exist yet); the
    ## earlier integral tags emit the reported/suspected/export streams. The
    ## `propertynames` guard below keeps each tag to the columns its
    ## `forecast_reported` actually produces. Labels match the scorer's
    ## `STREAM_HISTORY`/`STREAM_ASSEMBLED` in `scripts/score_releases.jl`.
    streams = (
        (:confirmed_new, "confirmed cases"),
        (:confirmed_deaths_new, "confirmed deaths"),
        (:cases_new, "reported cases"),
        (:deaths_new, "suspected deaths"),
        (:exports_new, "exports"))
    out = DataFrame(made_date = Date[], horizon = Int[],
        target_date = Date[], stream = String[], draw = Int[],
        value = Float64[])
    for (horizon, fc) in runs
        h = Int(horizon)
        target = cutoff + Day(h)
        for (col, label) in streams
            col in propertynames(fc) || continue
            vals = fc[!, col]
            for (d, i) in enumerate(1:$(THIN):length(vals))
                push!(out, (cutoff, h, target, label, d, Float64(vals[i])))
            end
        end
    end
    isempty(out) && error("no archived streams in the forecast")
    mkpath(dirname("$(dest)"))
    CSV.write("$(dest)", out)
    @info "wrote archive" dest="$(dest)" rows=nrow(out) cutoff
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

    ## `renewal_release_tags` only reaches v1.4.0 on; the integral tags predate
    ## it and are appended so `backfill_one` can dispatch them to the integral
    ## driver (and `--only v1.3.0` resolves).
    all_tags = vcat(renewal_release_tags(), collect(INTEGRAL_TAGS))
    tags = isnothing(opts.only) ? all_tags : filter(==(opts.only), all_tags)
    isempty(tags) && error("no matching release tags (have $all_tags)")

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
