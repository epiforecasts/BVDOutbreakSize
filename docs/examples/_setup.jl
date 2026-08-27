# Shared setup for the analysis and sensitivity report pages. This is plain
# Julia (not a Literate page): both `analysis.jl` and `sensitivity.jl` include
# it so each page can render on its own from the same fitted chains. It loads
# the packages, the observations, the fit registry (`docs/fits/registry.jl`)
# and every model fit through the content-addressed cache (`fit_or_load`),
# then unpacks the named chains, cumulative-infection draws and display
# labels the pages share. In CI the fits are pre-populated by the per-fit
# matrix and loaded here; locally a missing fit is computed and cached on
# first use.

using Turing
using Distributions
using StatsFuns: logistic
using DataFrames: DataFrame, eachrow
import CSV
using Random
using Markdown
using Dates: Date, Day, value
using BVDOutbreakSize
import CairoMakie
## Loading TensorBoardLogger activates the `tensorboard_callback` extension
## so `fit_callback` can stream each fit to TensorBoard as well as a progress
## log. Set `BVD_FIT_LOG=none` to disable all fit logging (CI release builds).
using TensorBoardLogger

## Render figures at higher resolution so they stay crisp in the docs.
CairoMakie.activate!(type = "png", px_per_unit = 3)

Random.seed!(20260518)

## Guard the stateful setup against running twice in one module: the Literate
## render includes each page in its own sandbox module (always runs), but
## `scripts/run.jl` includes both pages into one session, where the second
## include reuses the state the first one loaded rather than refitting. The
## imports above are idempotent, so they stay outside the guard.
if !@isdefined(_BVD_SETUP_LOADED)
    _BVD_SETUP_LOADED = true

    ## Observations and grid-date helpers shared by both pages.
    obs = load_observations()
    ## Grid day-index (day n is the cut-off) back to a calendar date.
    grid_date(day) = obs.cutoff - Day(obs.n - day)
    ## Date a cumulative history last reports (the cut-off for streams that
    ## run to it, or the freeze date for streams that stop earlier).
    hist_last_date(h) = isempty(h.days) ? missing : grid_date(maximum(h.days))

    ## The forecast's count streams, split by whether the situation reports
    ## still update each one (`stream_reporting`). A stream that has stopped
    ## carries a cumulative total that only repeats its last reported value,
    ## so it can be projected but not validated against an observation.
    ## `scripts/score_releases.jl` withholds the same streams, though by its
    ## own per-target rule rather than this one. Both pages read the split
    ## from here so they agree.
    forecast_cum_cols = (:cases_cum, :deaths_cum, :confirmed_cum,
        :confirmed_deaths_cum, :recovered_cum)
    reporting_cum_cols = Tuple(c for c in forecast_cum_cols
    if stream_reporting(obs, c))
    stopped_cum_cols = Tuple(c for c in forecast_cum_cols
    if !stream_reporting(obs, c))
    ## The matching new-count columns, for a figure that takes the forecast
    ## frame column by column rather than a keyed NamedTuple.
    new_cols(cols) = [stream_forecast_columns(c).new for c in cols]
    ## Keep the entries of a stream-keyed NamedTuple (`observed`, `baseline`,
    ## `individual`) belonging to `cols`, whichever of a stream's cumulative
    ## or new-count column each side is keyed by.
    function keep_streams(nt, cols)
        ids = [stream_id(c) for c in cols]
        return NamedTuple(k => v
        for (k, v) in pairs(nt) if stream_id(k) in ids)
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
        DataFrame)

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
    RUN_SENSITIVITY = lowercase(strip(get(ENV, "BVD_RUN_SENSITIVITY",
        "false"))) in ("true", "1", "yes", "on")

    ## Every fit is loaded through the content-addressed cache (`fit_or_load`):
    ## reused when a fit with the same model source, data and settings already
    ## exists — produced once by the per-fit CI matrix (`.github/workflows/
    ## fit-matrix.yml`) or on the HPC — and refitted otherwise. Set
    ## `BVD_REFIT=all` to force a full refit. The loads still run through
    ## `fit_parallel`, so on a cold cache the joint overlaps the per-stream,
    ## frozen and (gated) sensitivity re-fits and keeps all cores busy; on a
    ## warm cache they deserialise in parallel.
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
    ## In CI the per-fit matrix produces every fit before the render, so a
    ## render cache miss is a bug (usually a wrong `BVD_FIT_CACHE`).
    ## `BVD_FIT_STRICT` makes such a miss fail in seconds naming the key,
    ## rather than silently refitting the whole report and hitting the
    ## render-job timeout. Off by default so a local cold build still fits.
    _strict = lowercase(strip(get(ENV, "BVD_FIT_STRICT", ""))) in
              ("all", "true", "1", "yes", "on")
    _fit_specs = build_fit_specs(obs;
        breakpoint = _BREAKPOINT, frozen_cutoffs = frozen_cutoffs,
        chamla_cutoff = chamla_cutoff,
        validation_cutoff = validation_cutoff,
        run_sensitivity = RUN_SENSITIVITY)
    _fit_results = fit_parallel([() -> fit_or_load(fit_key(s.id), s.thunk;
                                     cache_dir = _fit_cache_dir,
                                     refit = _refit_all,
                                     strict = _strict)
                                 for s in _fit_specs])
    _fits = Dict(s.id => r for (s, r) in zip(_fit_specs, _fit_results))

    chn_joint = _fits["joint"]
    chn_exports = _fits["exports"]
    chn_deaths = _fits["deaths"]
    chn_cases = _fits["cases"]
    chn_confirmed = _fits["confirmed"]
    chn_confirmed_deaths = _fits["confirmed_deaths"]
    chn_treatment = _fits["treatment"]
    chn_onsets = _fits["onsets"]
    frozen_lastweek = _fits["frozen_validation"]
    ## One frozen individual fit per stream at the same cut-off as
    ## `frozen_lastweek`, so the forecast validation section can show each
    ## stream's own model alongside the frozen joint. Keyed by the same
    ## `fit` ids the current-data individual fits use (`chn_cases`, …), so
    ## the two dicts read the same way.
    frozen_lastweek_streams = Dict(
        sid => _fits["frozen_validation_$sid"] for sid in VALIDATION_STREAM_IDS)
    frozen_results = [_fits["frozen_$c"] for c in frozen_cutoffs]
    frozen_by_cutoff = Dict(zip(frozen_cutoffs, frozen_results))
    frozen_by_cutoff[chamla_cutoff] = _fits["frozen_$chamla_cutoff"]
    frozen_C(c) = vec(Array(frozen_by_cutoff[c].chn[:C_T]))
    ## Basic reproduction number draws from a chain that walks its own
    ## renewal process, `exp` of the walk's log base `rt_state.log_R0`, the
    ## walk's starting value and a distinct quantity from the growth-clock
    ## rate `r0`. `nothing` for a chain carrying no walk base: an absent key
    ## throws rather than reading back empty, so the lookup is probed, the
    ## same way `_has_key` in src/forecast.jl probes a chain key.
    function r0_walk_draws(chn)
        try
            exp.(vec(Array(chn[Symbol("rt_state.log_R0")])))
        catch
            nothing
        end
    end
    ## Every frozen fit is a full joint fit, so it carries the same walk base
    ## chn_joint does.
    frozen_R0(c) = r0_walk_draws(frozen_by_cutoff[c].chn)
    if RUN_SENSITIVITY
        chn_joint_community_delay = _fits["sens_community_delay"]
        chn_joint_exp_growth_clock = _fits["sens_exp_growth_clock"]
    end

    posterior_C_joint = vec(Array(chn_joint[:C_T]))
    posterior_C_exports = vec(Array(chn_exports[:C_T]))
    posterior_C_deaths = vec(Array(chn_deaths[:C_T]))
    posterior_C_cases = vec(Array(chn_cases[:C_T]))
    posterior_C_confirmed = vec(Array(chn_confirmed[:C_T]))
    posterior_C_confirmed_deaths = vec(Array(chn_confirmed_deaths[:C_T]))
    posterior_C_treatment = vec(Array(chn_treatment[:C_T]))
    posterior_C_onsets = vec(Array(chn_onsets[:C_T]))

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
        Symbol("exports_state.travel_state.daily_travellers") => "daily travellers")

    ## Renewal-start day used to align the reconstructed R(t) knot grid
    ## with the model, shared by the main and sensitivity R(t) plots.
    _rt_start_plot = clamp(
        obs.n - round(Int, obs.tmrca_days) + RENEWAL_START_LEAD, 1, obs.n)
end # _BVD_SETUP_LOADED guard
