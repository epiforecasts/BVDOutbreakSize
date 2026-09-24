# Shared setup for the report pages. This is plain Julia (not a Literate
# page): every page under `docs/pages/` includes it so each page can render
# on its own from the same fitted chains. It loads
# the packages, the observations and the fit registry (`docs/fits/registry.jl`),
# and defines the accessors the pages read fits and prior draws through
# (`load_fit`, `joint_prior_draws`). Each fit is loaded through the
# content-addressed cache (`fit_or_load`) the first time a page asks for it.
# In CI the fits are pre-populated by the per-fit matrix; locally a missing
# fit is computed and cached on first use.

_setup_t0 = time()
using Turing
using Distributions
using StatsFuns: logistic
using DataFrames: DataFrame, eachrow
import CSV
using Random
using Markdown
using Logging: SimpleLogger, with_logger
using Dates: Date, Day, value
using BVDOutbreakSize
import CairoMakie
## Loading TensorBoardLogger activates the `tensorboard_callback` extension
## so `fit_callback` can stream each fit to TensorBoard as well as a progress
## log. Set `BVD_FIT_LOG=none` to disable all fit logging (CI release builds).
using TensorBoardLogger

## Render figures at higher resolution so they stay crisp in the docs.
CairoMakie.activate!(type = "png", px_per_unit = 3)

## `activate!(type = "png")` leaves Figures still showable as `MIME"text/html"`,
## which Literate prefers over `image/png`, and Documenter's raw-block regex
## hits PCRE's ~64KB limit once a block grows past it. Disabling the
## html-family mimes forces every Figure display onto the image/png path.
CairoMakie.disable_mime!(
    "text/html", "application/vnd.webio.application+html",
    "application/prs.juno.plotpane+html", "juliavscode/html",
    "svg", "pdf"
)

Random.seed!(20260518)

## Guard the stateful setup against running twice in one module: the Literate
## render includes each page in its own sandbox module (always runs), but
## `scripts/run.jl` includes both pages into one session, where the second
## include reuses the state the first one loaded rather than refitting. The
## imports above are idempotent, so they stay outside the guard.
if !@isdefined(_BVD_SETUP_LOADED)
    _BVD_SETUP_LOADED = true

    ## Timing lines for the render. Literate captures a page's stderr and
    ## logs and drops them from the job log, so when `docs/execute.jl` points
    ## `BVD_RENDER_LOG` at a file the lines go there and are printed after
    ## the page. Otherwise (`scripts/run.jl`) they are logged.
    _render_log_path() = strip(get(ENV, "BVD_RENDER_LOG", ""))
    function _render_log(line)
        path = _render_log_path()
        isempty(path) && return @info line
        open(io -> println(io, line), path, "a")
        return nothing
    end
    _since(t0) = "$(round(time() - t0; digits = 1)) s"
    ## Run `f()` and log how long `what` took. Under `BVD_RENDER_LOG`,
    ## anything `f` writes to stderr (the fit cache hits) goes to the file
    ## too.
    function _timed(f, what)
        t0 = time()
        path = _render_log_path()
        result = if isempty(path)
            f()
        else
            open(path, "a") do io
                with_logger(() -> redirect_stderr(f, io), SimpleLogger(io))
            end
        end
        _render_log("$what: $(_since(t0))")
        return result
    end
    _render_log("setup packages loaded: $(_since(_setup_t0))")

    ## Observations and grid-date helpers shared by both pages.
    obs = load_observations()
    ## Grid day-index (day n is the cut-off) back to a calendar date.
    grid_date(day) = obs.cutoff - Day(obs.n - day)
    ## Date a cumulative history last reports (the cut-off for streams that
    ## run to it, or the freeze date for streams that stop earlier).
    hist_last_date(h) = history_last_date(grid_date, h)

    ## The forecast's count streams, split by whether the situation reports
    ## still update each one (`stream_reporting`). A stream that has stopped
    ## carries a cumulative total that only repeats its last reported value,
    ## so it can be projected but not validated against an observation.
    ## `scripts/score_releases.jl` withholds the same streams, though by its
    ## own per-target rule rather than this one. Both pages read the split
    ## from here so they agree.
    forecast_cum_cols = (
        :cases_cum, :deaths_cum, :confirmed_cum,
        :confirmed_deaths_cum, :recovered_cum,
    )
    reporting_cum_cols = Tuple(
        c for c in forecast_cum_cols
            if stream_reporting(obs, c)
    )
    stopped_cum_cols = Tuple(
        c for c in forecast_cum_cols
            if !stream_reporting(obs, c)
    )
    ## The matching new-count columns, for a figure that takes the forecast
    ## frame column by column rather than a keyed NamedTuple.
    new_cols(cols) = [stream_forecast_columns(c).new for c in cols]
    ## Keep the entries of a stream-keyed NamedTuple (`observed`, `baseline`,
    ## `individual`) belonging to `cols`, whichever of a stream's cumulative
    ## or new-count column each side is keyed by.
    function keep_streams(nt, cols)
        ids = [stream_id(c) for c in cols]
        return NamedTuple(
            k => v
                for (k, v) in pairs(nt) if stream_id(k) in ids
        )
    end

    ## The fits are defined once in `docs/fits/registry.jl` as a registry, so
    ## each can be run and cached independently — one per CI matrix job, or
    ## an HPC task — and loaded here through the content-addressed cache
    ## instead of being refitted inline. `_BREAKPOINT`, `validation_cutoff`
    ## and `frozen_cutoffs` come from the same registry so the report and the
    ## standalone fits agree.
    include(joinpath(pkgdir(BVDOutbreakSize), "docs", "fits", "registry.jl"))
    _BREAKPOINT = default_breakpoint(obs)

    ## The exports-deaths composer keeps the deaths and exports submodels only
    ## for their CFR, onset-to-death PMF and export onsets, leaving their own
    ## counts missing, which leaves two redundant sampled discrete draws, so its
    ## model check is disabled (see `nuts_sample`).
    validation_cutoff = default_validation_cutoff(obs)

    ## Published per-release estimates, pulled from the tagged results
    ## releases by `scripts/refresh_releases.jl` into
    ## `data/released_estimates.csv`. Columns: tag, date, model (integral or
    ## renewal), median and the 30/60/90% bounds.
    released_df = CSV.read(
        joinpath(pkgdir(BVDOutbreakSize), "data", "released_estimates.csv"),
        DataFrame
    )

    ## Frozen joint re-fits at the cut-offs McCabe et al. used (27 May matches
    ## the Lancet publication's cut-off), for the matched-in-time comparison
    ## further down. The estimate-evolution overlay relies on the published
    ## per-release estimates in `released_df`, so no per-release-date current-
    ## model re-fits are run here.
    frozen_cutoffs = default_frozen_cutoffs()
    ## Chamla et al.'s 8 June confirmed-case anchor, kept out of the
    ## McCabe-matched `frozen_cutoffs` but reused for the Chamla comparison and
    ## the estimate-evolution overlay.
    chamla_cutoff = default_chamla_cutoff()

    ## The frozen-joint, sensitivity-variant and delay/clock helpers used by
    ## the re-fits live in `docs/fits/registry.jl` (`build_fit_specs`), so
    ## they can be run from the standalone per-fit entry point too.

    ## Sensitivity refits (onset-to-death delay, molecular clock) are slow extra
    ## joint fits, gated on the `BVD_RUN_SENSITIVITY` env var. They run on
    ## release builds only (tag pushes, which deploy the versioned docs) and are
    ## skipped on main pushes and PR previews to keep those docs builds fast:
    ## `.github/workflows/docs.yml` sets the var to
    ## `startsWith(github.ref, 'refs/tags/')`. Set `BVD_RUN_SENSITIVITY=true`
    ## to run them locally.
    RUN_SENSITIVITY = lowercase(
        strip(
            get(
                ENV, "BVD_RUN_SENSITIVITY",
                "false"
            )
        )
    ) in ("true", "1", "yes", "on")

    ## Every fit is loaded through the content-addressed cache (`fit_or_load`):
    ## reused when a fit with the same model source, data and settings already
    ## exists — produced once by the `fit_*` jobs in
    ## `.github/workflows/docs.yml` or on the HPC — and refitted otherwise. Set
    ## `BVD_REFIT=all` to force a full refit. Each fit is loaded on first use
    ## (`load_fit`), so a cold cache fits one at a time; `task fit-all` runs
    ## the whole registry in parallel first.
    ## Resolve the cache dir against the package root, never the working
    ## directory: Literate executes the page with the cwd changed to docs/src,
    ## so a relative `BVD_FIT_CACHE` (as CI passes) would point at
    ## docs/src/logs/fit_cache, miss every cached chain and refit the whole
    ## report. Relative overrides are resolved against `pkgdir`; absolute ones
    ## are used as-is.
    _fit_cache_dir = let c = strip(get(ENV, "BVD_FIT_CACHE", ""))
        if isempty(c)
            joinpath(pkgdir(BVDOutbreakSize), "logs", "fit_cache")
        elseif isabspath(c)
            String(c)
        else
            joinpath(pkgdir(BVDOutbreakSize), c)
        end
    end
    _refit_all = lowercase(strip(get(ENV, "BVD_REFIT", ""))) in
        ("all", "true", "1")
    ## Fits are produced by the per-fit CI matrix (or `task fit-all` /
    ## `task fetch-fits` locally) before the render, so a render cache miss is
    ## a bug, usually a wrong `BVD_FIT_CACHE` or a fit the matrix never made.
    ## Strict mode fails in seconds naming the key rather than silently
    ## refitting the whole report into the render-job timeout. On by default;
    ## `BVD_FIT_STRICT=false` restores inline fitting for a cold local build.
    _strict = lowercase(strip(get(ENV, "BVD_FIT_STRICT", ""))) ∉
        ("false", "0", "no", "off")
    _fit_specs = build_fit_specs(
        obs;
        breakpoint = _BREAKPOINT, frozen_cutoffs = frozen_cutoffs,
        chamla_cutoff = chamla_cutoff,
        validation_cutoff = validation_cutoff,
        run_sensitivity = RUN_SENSITIVITY
    )
    _fit_spec_by_id = Dict(s.id => s for s in _fit_specs)
    _loaded_fits = Dict{String, Any}()
    ## Fit `id` from the registry, loaded on first use and kept for the
    ## session, so a page loads only the fits it reads and `scripts/run.jl`
    ## loads each one once across every page.
    function load_fit(id::AbstractString)
        haskey(_loaded_fits, id) && return _loaded_fits[id]
        spec = get(_fit_spec_by_id, id, nothing)
        spec === nothing && error(
            "no fit \"$id\" in the registry for this build; known ids: " *
                join(sort(collect(keys(_fit_spec_by_id))), ", ")
        )
        result = _timed("fit $id") do
            fit_or_load(
                fit_key(id), spec.thunk;
                cache_dir = _fit_cache_dir,
                refit = _refit_all,
                strict = _strict
            )
        end
        _loaded_fits[id] = result
        return result
    end

    ## The model each fit sampled, rebuilt to forecast from it.
    _fit_models = Dict(s.id => s.model for s in _fit_specs if haskey(s, :model))
    fit_model(id::AbstractString) = _fit_models[id]()

    ## Days past the cut-off every forecast is drawn to. The pages read the
    ## one-week forecast and the release archive the four weekly horizons.
    FORECAST_HORIZON = 28
    _forecast_cache = Dict{String, Any}()
    ## Posterior-predictive draws past the cut-off from fit `id`: its model
    ## run on with `predict` over its chain (see `forecast_draws`). Drawn once
    ## per page and shared by every forecast the page reads from that fit.
    function fit_forecast(id::AbstractString)
        return get!(_forecast_cache, id) do
            r = load_fit(id)
            chn = r isa NamedTuple ? r.chn : r
            forecast_draws(fit_model(id), chn; horizon = FORECAST_HORIZON)
        end
    end

    ## Draws from the joint prior, with every observation withheld. The
    ## in-sample page shows them as the prior predictive check; the national
    ## page overlays them on each posterior. The draw is kept after the first
    ## call and takes its own seeded generator, so both pages, and
    ## `scripts/run.jl`, overlay the same draws.
    _joint_prior_cache = Ref{Any}(nothing)
    function joint_prior_draws()
        cached = _joint_prior_cache[]
        cached !== nothing && return cached
        m = bvd_joint(
            obs.n, missing, missing, missing, missing, missing;
            deaths_history = (; days = Int[], counts = Int[]),
            reported_history = (; days = Int[], counts = Int[]),
            confirmed_history = (; days = Int[], counts = Int[]),
            export_case_days = obs.export_case_days,
            export_death_days = obs.export_death_days,
            breakpoint = obs.n - obs.who_first_sitrep_days,
            background_pooling = background_pooling_model,
            genetic = genetic_seeding_model,
            tmrca_days = obs.tmrca_days
        )
        chn = _timed("joint prior draws") do
            sample(Xoshiro(20260518), m, Prior(), 2_000; progress = false)
        end
        _joint_prior_cache[] = chn
        return chn
    end

    ## One frozen individual fit per stream at the validation cut-off, so the
    ## forecast validation can show each stream's own model alongside the
    ## frozen joint (`frozen_validation`). Keyed by the same ids the
    ## current-data individual fits use. Only the still-reported streams are
    ## fitted, so a stream that has stopped is absent rather than filtered
    ## out later (see `validation_stream_ids`).
    frozen_validation_stream_fits() = Dict(
        sid => load_fit("frozen_validation_$sid")
            for sid in validation_stream_ids(obs)
    )
    ## The frozen joint fits keyed by cut-off: the McCabe-matched
    ## `frozen_cutoffs` and Chamla's anchor. A new Dict on each call, since
    ## a page may add its own entries.
    function frozen_fits_by_cutoff()
        fits = Dict(c => load_fit("frozen_$c") for c in frozen_cutoffs)
        fits[chamla_cutoff] = load_fit("frozen_$chamla_cutoff")
        return fits
    end
    ## Basic reproduction number draws from a chain that walks its own
    ## renewal process, `exp` of the walk's log base `rt_state.log_R0`, the
    ## walk's starting value and a distinct quantity from the growth-clock
    ## rate `r0`. `nothing` for a chain carrying no walk base: an absent key
    ## throws rather than reading back empty, so the lookup is probed, the
    ## same way `_has_key` in src/forecast.jl probes a chain key.
    function r0_walk_draws(chn)
        return try
            exp.(vec(Array(chn[Symbol("rt_state.log_R0")])))
        catch
            nothing
        end
    end

    ## Per-province spatial-table data, reshaped once (a Dict{String} lookup
    ## inside a model body puts a memcmp foreigncall on the AD tape).
    ## `N_PATCHES` is the report's single source for the patch count, so a
    ## change to `PROVINCE_NAMES` reaches every table and figure at once.
    N_PATCHES = length(PROVINCE_NAMES)
    province_cases = province_increment_matrix(
        obs.province_confirmed_history, PROVINCE_NAMES,
        length(PROVINCE_NAMES)
    )
    province_deaths = province_increment_matrix(
        obs.province_death_history, PROVINCE_NAMES,
        length(PROVINCE_NAMES)
    )
    province_testing = province_testing_covariate(
        obs.province_lab_daily_history
    )

    ## Draws from the four-patch prior for the province pages. The
    ## `joint_prior_draws` are single-population and carry no province
    ## quantities. Every observation is withheld as there, but the province
    ## compositions are passed their observed counts: `Prior()` leaves them out
    ## of the density, so they only fix each vintage's total, and with
    ## `missing` an extreme prior draw can push that total past the integer
    ## range. A function rather than an eager draw, since every page includes
    ## this file and only the province pages need it. The draw is kept after
    ## the first call and takes its own seeded generator, so both pages, and
    ## `scripts/run.jl`, overlay the same draws.
    _patch_prior_cache = Ref{Any}(nothing)
    function patch_prior_draws(obs)
        cached = _patch_prior_cache[]
        cached !== nothing && first(cached) === obs && return last(cached)
        m = bvd_joint(
            obs.n, missing, missing, missing, missing, missing;
            deaths_history = (; days = Int[], counts = Int[]),
            reported_history = (; days = Int[], counts = Int[]),
            confirmed_history = (; days = Int[], counts = Int[]),
            export_case_days = obs.export_case_days,
            export_death_days = obs.export_death_days,
            breakpoint = obs.n - obs.who_first_sitrep_days,
            background_pooling = background_pooling_model,
            genetic = genetic_seeding_model,
            tmrca_days = obs.tmrca_days,
            n_patches = N_PATCHES,
            province_increments = province_cases.increments,
            province_days = province_cases.days,
            province_testing_covariate = province_testing,
            province_death_increments = province_deaths.increments,
            province_death_days = province_deaths.days
        )
        chn = _timed("patch prior draws") do
            sample(Xoshiro(20260518), m, Prior(), 1_000; progress = false)
        end
        _patch_prior_cache[] = (obs, chn)
        return chn
    end

    ## Joint posterior predictive for the in-sample pages. A function rather
    ## than an eager draw, since every page includes this file and only the
    ## in-sample pages need it. The draw is kept after the first call, so
    ## `scripts/run.jl` runs `predict` once for both pages.
    _joint_pp_cache = Ref{Any}(nothing)
    function joint_posterior_predictive()
        cached = _joint_pp_cache[]
        cached !== nothing && return cached
        ## The headline patch joint with each stream's counts dropped but its
        ## observation grid kept (`generator_joint`), so `predict` resamples
        ## every stream over the real cells from each fitted draw.
        chn_joint = load_fit("joint")
        t0 = time()
        pp = predict(
            generator_joint(obs; breakpoint = _BREAKPOINT), chn_joint
        )
        _render_log("joint posterior predictive: $(_since(t0))")
        _joint_pp_cache[] = pp
        return pp
    end
    ## The parameter-recovery results the `recovery` CI job writes
    ## (`scripts/recovery.jl`): one table of recovered quantities and one of
    ## forecasts, over every seed found, or empty tables when no run is
    ## available.
    function recovery_results()
        dir = joinpath(
            get(ENV, "BVD_OUTPUT_DIR", joinpath(pkgdir(BVDOutbreakSize), "output")),
            "recovery"
        )
        read_all(prefix) = begin
            files = isdir(dir) ? sort(
                    filter(
                        f -> startswith(f, prefix) && endswith(f, ".csv"),
                        readdir(dir)
                    )
                ) : String[]
            isempty(files) ? DataFrame() :
                reduce(vcat, [CSV.read(joinpath(dir, f), DataFrame) for f in files])
        end
        return (;
            params = read_all("recovery_"),
            forecasts = read_all("forecast_recovery_"),
        )
    end
    ## Clean display names for the summary tables and pair plots. The submodel
    ## prefixes (`rt_state.`, `gi_state.`, ...) are kept in the model so the
    ## nested submodels stay distinct; this map only relabels them for display.
    display_names = Dict{Symbol, String}(
        Symbol("rt_state.sigma_rw") => "Rt step size",
        Symbol("rt_state.intervention_effect") => "intervention effect",
        Symbol("gi_state.α") => "generation interval shape",
        Symbol("gi_state.θ") => "generation interval scale",
        Symbol("inc_state.delay_mean") => "incubation period mean",
        Symbol("inc_state.delay_sd") => "incubation period SD",
        Symbol("cases_state.report_state.α") => "onset-to-report shape",
        Symbol("cases_state.report_state.θ") => "onset-to-report scale",
        Symbol("deaths_state.od_state.oa.α") => "onset-to-admission shape",
        Symbol("deaths_state.od_state.oa.θ") => "onset-to-admission scale",
        Symbol("deaths_state.od_state.ad.α") => "admission-to-death shape",
        Symbol("deaths_state.od_state.ad.θ") => "admission-to-death scale",
        Symbol("exports_state.detect_state.α") => "onset-to-detection shape",
        Symbol("exports_state.detect_state.θ") => "onset-to-detection scale",
        Symbol("confirmed_state.receipt_state.d.delay_mean") => "report-to-receipt mean",
        Symbol("confirmed_state.receipt_state.d.delay_sd") => "report-to-receipt SD",
        :isolation_bvd_los_mean => "in-care BVD length-of-stay mean (mixture)",
        :isolation_death_los_mean => "in-care admission-to-death stay mean",
        :isolation_recovery_los_mean => "in-care admission-to-recovery stay mean",
        :isolation_admission_delay_mean => "suspected-to-admission delay mean",
        :isolation_ruleout_los_mean => "isolation non-BVD rule-out stay mean",
        :incare_cfr => "in-care fatality (CFR_iso)",
        :incare_cfr_modifier => "in-care fatality log-odds modifier",
        :incare_confirm_modifier => "in-care confirmation-rate modifier",
        :abscond_fraction => "daily abscond fraction",
        :recovery_delay_mean => "confirmation-to-recovery mean",
        Symbol("exports_state.travel_state.daily_travellers") => "daily travellers"
    )

    ## Renewal-start day used to align the reconstructed R(t) knot grid
    ## with the model, shared by the main and sensitivity R(t) plots.
    _rt_start_plot = clamp(
        obs.n - round(Int, obs.tmrca_days) + RENEWAL_START_LEAD, 1, obs.n
    )

    ## The symptom-onset triangle's own grid, derived from the observations
    ## exactly as `onset_reporting_model` derives it, since it is data rather
    ## than chain contents. The national page's reporting-delay section, the
    ## in-sample page's snapshot nowcasts and the forecast page's nowcast
    ## read it.
    _onset_grid_start = isempty(obs.onset_curve_history.onset_days) ? 1 :
        minimum(obs.onset_curve_history.onset_days)
    _onset_grid_end = isempty(obs.onset_curve_history.report_days) ?
        _onset_grid_start :
        max(
            maximum(obs.onset_curve_history.report_days),
            _onset_grid_start
        )

    ## The digitised onset snapshots up to the cut-off and, for each onset
    ## date they cover, the latest printed count and the report day it came
    ## from. The fitted stream holds only the corrections between snapshots,
    ## so the printed levels are read from the source blocks. A date inside a
    ## snapshot's printed extent but with no row is a zero-height bar. A date
    ## outside it is not covered by that snapshot.
    function onset_snapshot_readings()
        path = joinpath(
            pkgdir(BVDOutbreakSize), "data", "onset_curve_scanned.csv"
        )
        snaps = filter(
            b -> b.report_date <= obs.cutoff,
            BVDOutbreakSize._dedup_onset_blocks(
                BVDOutbreakSize._read_onset_curve_blocks(path)
            )
        )
        last_printed = Dict{Int, Float64}()
        last_report_day = Dict{Int, Int}()
        for snap in snaps
            lo, hi = extrema(keys(snap.onsets))
            R = obs.n - value(obs.cutoff - snap.report_date)
            for d in lo:Day(1):hi
                u = obs.n - value(obs.cutoff - d)
                (1 <= u <= obs.n) || continue
                last_printed[u] = Float64(get(snap.onsets, d, 0))
                last_report_day[u] = R
            end
        end
        return (; snaps, last_printed, last_report_day)
    end

    ## Daily modelled onsets per posterior draw, differenced from the
    ## cumulative onsets the chain stores.
    onset_daily_draws(chn) = [
        vcat(v[1], diff(v)) for v in vec(collect(chn[:cumulative_onsets]))
    ]

    ## A function giving posterior predictive draws of one digitised bar
    ## from per-draw modelled counts, through the onset stream's
    ## measurement error for a single read (`onset_report_scale`'s level
    ## case). Four replicates per draw keep the 90% band edge from being
    ## ragged with Monte Carlo error.
    function onset_bar_replicator(chn, rng)
        σ_mult = vec(collect(chn[Symbol("onset_report_state.σ_mult")]))
        σ_scan = vec(collect(chn[Symbol("onset_report_state.σ_scan")]))
        return draws -> [
            begin
                μ = draws[i]
                σ = σ_mult[i] *
                    onset_report_scale(μ, μ, 0.0, 1; scan_sd = σ_scan[i])
                μ + σ * rand(rng, TDist(4.0))
            end
                for _ in 1:4 for i in eachindex(draws)
        ]
    end

    _render_log("setup done: $(_since(_setup_t0))")

    ## Cross-release score tables written by `scripts/score_releases.jl`.
    ## The committed files are header-only until a release carries the asset,
    ## so the common path reads a real file to a zero-row frame; the typed
    ## `schema` is the fallback for a file that is absent entirely, since
    ## CSV.read throws on a missing path and would take the build with it.
    function _release_data(name, schema::NamedTuple)
        path = joinpath(pkgdir(BVDOutbreakSize), "data", name)
        isfile(path) && return CSV.read(path, DataFrame)
        return DataFrame([k => T[] for (k, T) in pairs(schema)])
    end

    ## One-week-ahead forecast from a frozen fit, the one the evaluation
    ## pages score against what has since been observed. A function rather
    ## than a value, so only the pages that validate pay for the forecast.
    ## `obs_recovered` is passed so the forecast carries a `recovered_new`
    ## column (materialised only when the recovered origin is given), letting
    ## the recovered stream be scored like the other streams. The frozen
    ## model is rebuilt from its own frozen observations, so the onset
    ## forecast runs on the triangle the frozen fit was fitted to.
    function validation_forecast_from(id::AbstractString)
        o = load_fit(id).o
        return forecast_reported(
            fit_forecast(id);
            horizon = 7,
            obs_cases = o.reported_cases,
            obs_deaths = o.total_deaths,
            obs_confirmed = o.confirmed_cases,
            obs_confirmed_deaths = o.confirmed_deaths,
            obs_recovered = o.recovered_cases
        )
    end
end # _BVD_SETUP_LOADED guard
