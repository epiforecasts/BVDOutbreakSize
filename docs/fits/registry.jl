# Registry of the expensive model fits in the analysis report. Each fit is
# defined once here as an `(id, kind, needs, thunk)` entry, so it can be run
# and cached independently — one per CI matrix job, or an HPC task — and the
# docs build then loads the chains through the content-addressed cache
# instead of refitting them inline. `build_fit_specs` mirrors the model
# calls in `docs/examples/analysis.jl`; keep the two in step.
#
# `needs` lists the ids a fit is initialised from: empty for a base fit,
# the parent id for a dependent fit, which runs in a second stage once its
# parent is cached (`base_fit_specs`, `dependent_fit_specs`).

include(joinpath(@__DIR__, "cache.jl"))

using BVDOutbreakSize
using Dates: Date, Day, value
using Distributions: truncated, Normal

const _PKG = pkgdir(BVDOutbreakSize)

## Source files whose contents define the fits: the model and its submodels,
## the renewal maths, the sampler, the data pipeline and this registry. A
## change to any of them invalidates every cached fit; plotting and reporting
## code (plots.jl, summaries.jl, ...) deliberately does not. `cache.jl` is
## here because it defines what the key covers: a change to the hashing rule
## that left the key alone would give two different rules the same key.
const FIT_SOURCE_FILES = [
    joinpath(_PKG, "src", "models", "priors.jl"),
    joinpath(_PKG, "src", "models", "observations.jl"),
    joinpath(_PKG, "src", "models", "joint.jl"),
    joinpath(_PKG, "src", "models", "zone.jl"),
    joinpath(_PKG, "src", "renewal.jl"),
    joinpath(_PKG, "src", "sampling.jl"),
    joinpath(_PKG, "src", "constants.jl"),
    joinpath(_PKG, "src", "data.jl"),
    joinpath(_PKG, "src", "onset_curve.jl"),
    joinpath(@__DIR__, "cache.jl"),
    @__FILE__
]

## Bump when the cache layout changes in a way that must invalidate old files.
const FIT_CACHE_SCHEMA = "v1"

## Data files under `data/` that are not fit inputs and so must not bust the
## cache. `released_estimates.csv` is a published-estimate overlay used only by
## the sensitivity page's evolution figure; the render job rewrites it before
## rendering, which would otherwise change every fit key and force a refit.
## `rt_by_release.csv`, `r0_by_release.csv`, `forecast_scores.csv`,
## `forecast_scores_frozen.csv`, `forecast_overlay.csv` and
## `forecast_overlay_frozen.csv` are the per-release R_T, R0 and
## forecast-scoring overlays that `scripts/score_releases.jl` rewrites before
## rendering from the published releases, so they are excluded for the same
## reason. `rt_by_release_by_stream.csv`, `size_by_release_by_stream.csv`
## and `r0_by_release_by_stream.csv` are the per-fit versions of the same
## overlay, rewritten the same way.
## Every file score_releases.jl writes into data/ must be listed here, or the
## render's data hash diverges from the fit matrix's and every fit misses.
##
## The digest covers every file under `data/` whatever its format, so this
## list is the only thing that keeps a file out of the key. `README.md`
## documents the directory and `sitrep_pdfs` holds the situation report PDFs
## the scans are taken from. Neither is read by the model, and hashing 230 MB
## of PDFs would refit the whole report each time one is downloaded. An entry
## naming a directory drops everything under it.
const FIT_DATA_EXCLUDE = (
    "released_estimates.csv",
    "rt_by_release.csv", "r0_by_release.csv",
    "forecast_scores.csv", "forecast_scores_frozen.csv",
    "forecast_overlay.csv", "forecast_overlay_frozen.csv",
    "rt_by_release_by_stream.csv", "size_by_release_by_stream.csv",
    "r0_by_release_by_stream.csv", "province_forecast_scores.csv",
    "README.md", "sitrep_pdfs",
)

"Content hash of the fit-relevant source, data and sampler settings."
function fit_content_hash(; samples::Integer = 500, chains::Integer = 2)
    return content_hash(
        FIT_SOURCE_FILES;
        data_dir = joinpath(_PKG, "data"),
        data_exclude = FIT_DATA_EXCLUDE,
        extra = string(FIT_CACHE_SCHEMA, ":", samples, "x", chains)
    )
end

"Content-addressed cache key for fit `id` at the given sampler settings."
function fit_key(id; samples::Integer = 500, chains::Integer = 2)
    return string(id, "__", fit_content_hash(; samples, chains))
end

"""
    fit_cache_dir() -> String

The fit cache directory: `BVD_FIT_CACHE` when set, else `logs/fit_cache`
under the package root. A relative override is resolved against the package
root, not the working directory, because Literate runs a page from
`docs/src`.
"""
function fit_cache_dir()
    c = strip(get(ENV, "BVD_FIT_CACHE", ""))
    return if isempty(c)
        joinpath(_PKG, "logs", "fit_cache")
    elseif isabspath(c)
        String(c)
    else
        joinpath(_PKG, c)
    end
end

## Canonical fit-setup values, so `analysis.jl`, `fit_one.jl` and this registry
## agree on the breakpoint, the validation cut-off and the frozen cut-offs.
default_breakpoint(obs) = obs.n - obs.who_first_sitrep_days
default_validation_cutoff(obs) = string(obs.cutoff - Day(7))
default_frozen_cutoffs() = ["2026-05-20", "2026-05-23", "2026-05-27"]
## Chamla et al.'s confirmed-case calibration anchor (598 confirmed by 8 June),
## fit separately from the McCabe frozen set so the confirmed-case projection
## rides a vintage with testing data rather than the near-empty 27 May stream.
default_chamla_cutoff() = "2026-06-08"
## The individual single-stream models that could be frozen at
## `validation_cutoff` for the "last week versus now" forecast validation
## (see `fit_frozen_stream` and `docs/examples/sensitivity.jl`'s "Forecast
## validation" section), paired with the observation stream each fits.
## "exports" is excluded (not forecast at all) and there is no individual
## model for "recovered".
const VALIDATION_STREAM_FITS = (
    ("cases", :suspected_cases),
    ("deaths", :suspected_deaths),
    ("confirmed", :confirmed_cases),
    ("confirmed_deaths", :confirmed_deaths),
    ("treatment", :isolation_beds),
)

"""
Fit ids of the single-stream models the forecast validation draws, for the
data in `obs`: the streams the situation reports are still updating.

The validation panels take the still-reported streams, so a fit of a stream
that has stopped feeds nothing drawn and is two NUTS fits a docs build does
not need. Deriving the list from the reporting status rather than holding a
second hand-maintained one means a stream that starts being reported again
comes back into validation on its own.
"""
function validation_stream_ids(obs)
    return Tuple(
        id
            for (id, stream) in VALIDATION_STREAM_FITS
            if stream_reporting(obs, stream)
    )
end
function run_sensitivity_env()
    return lowercase(strip(get(ENV, "BVD_RUN_SENSITIVITY", "false"))) in
        ("true", "1", "yes", "on")
end

"""
    build_fit_specs(obs; breakpoint, frozen_cutoffs, validation_cutoff,
                    run_sensitivity, samples = 500, chains = 2,
                    zone_fitter = nothing, cache_dir = fit_cache_dir())

Ordered list of the report's fits as `(; id, kind, needs, thunk)` named
tuples. `kind` is `:chain` for the headline joint and single-stream fits or
`:frozen` for the frozen/validation joints (whose thunk returns
`(; cutoff, o, chn)`). `needs` names the fits a thunk loads from `cache_dir`
before running: empty for a base fit, the parent id for the health-zone fits
`local` (from `joint`) and `local_frozen_validation` (from
`frozen_validation`). A dependent thunk loads its parent strictly; a missing
parent is an error, not a refit. `zone_fitter` fits a zone model from a
parent chain, `BVDOutbreakSize.fit_zone` when `nothing`, resolved when the
thunk runs. The sensitivity re-fits, including the zone variants
`local_mixing`, `local_parent_low` and
`local_parent_high`, are appended only when `run_sensitivity` is true.
"""
## Sampler settings for the headline and its spatial control.
##
## Both fits must use the same draw count, adaptation and acceptance target.
## They are the two halves of the spatial sensitivity, so a difference
## between their `C_T` posteriors reads as evidence about the spatial
## structure only if nothing else differs.
##
## 750 draws and 500 adaptation steps at a target acceptance of 0.80. The
## effective sample size is limited by adaptation rather than by the draw
## count, so the budget goes there: at 900 draws and 400 adaptation steps the
## headline returned 48 bulk and 39 tail against the control's 256 and 203,
## at a worst R-hat of 1.07 against 1.02.
##
## The binding constraint is the fit job's `timeout-minutes: 350`, under a
## hard six-hour ceiling on a GitHub-hosted job. At 900 draws and 400
## adaptation steps the headline took 320 minutes and the control 200. Total
## iterations are held at or below that, so adaptation is bought with draws
## rather than with wall-clock.
##
## NUTS terminates at the tree-depth cap on every iteration of both fits, at
## 1023 leapfrog steps, rather than at a U-turn. Raising `max_depth` is the
## direct fix and each extra level doubles the leapfrog steps, which does not
## fit the budget. A lower acceptance target lengthens the step and shortens
## the trajectories instead.
##
## `BVD_JOINT_SAMPLES`, `BVD_JOINT_WARMUP` and `BVD_JOINT_TARGET_ACCEPT`
## override all three without editing this file.
joint_target_accept() = parse(
    Float64,
    get(ENV, "BVD_JOINT_TARGET_ACCEPT", "0.80")
)
joint_samples(default::Integer) = parse(
    Int,
    get(ENV, "BVD_JOINT_SAMPLES", string(default))
)
## `nuts_sample` caps its own default at 200 adaptation steps, which is where
## this fit is short.
joint_warmup(default::Integer) = parse(
    Int,
    get(ENV, "BVD_JOINT_WARMUP", string(default))
)

## Sampler settings for the health-zone fits. `BVD_ZONE_SAMPLES` and
## `BVD_ZONE_WARMUP` override the draw and adaptation counts. Like
## `BVD_JOINT_*` they sit outside the content hash, so a short local run
## writes its chain under the production key.
zone_samples(default::Integer) = parse(
    Int,
    get(ENV, "BVD_ZONE_SAMPLES", string(default))
)
zone_warmup(default::Integer) = parse(
    Int,
    get(ENV, "BVD_ZONE_WARMUP", string(default))
)
zone_target_accept() = 0.8
zone_max_depth() = 8

## Looked up when a dependent thunk runs, so the registry builds without
## `fit_zone` and tests can inject a double through `zone_fitter`.
function default_zone_fitter()
    isdefined(BVDOutbreakSize, :fit_zone) ||
        error(
        "BVDOutbreakSize.fit_zone is not defined: the health-zone " *
            "model (src/models/zone.jl) is needed to run a dependent fit"
    )
    return BVDOutbreakSize.fit_zone
end

function build_fit_specs(
        obs;
        breakpoint = default_breakpoint(obs),
        frozen_cutoffs = default_frozen_cutoffs(),
        chamla_cutoff = default_chamla_cutoff(),
        validation_cutoff = default_validation_cutoff(obs),
        run_sensitivity = run_sensitivity_env(),
        samples::Integer = 500,
        chains::Integer = 2,
        zone_fitter = nothing,
        cache_dir::AbstractString = fit_cache_dir()
    )

    ## A joint fit at the headline settings to the data frozen at `cutoff_date`.
    ## `patches` turns the spatial structure on for this frozen fit. Only the
    ## one-week-back validation uses it, because it is the only frozen fit
    ## whose forecast is scored by province. The McCabe and Chamla
    ## comparisons stay single-population: they are set against external
    ## national estimates, where patches add cost and nothing else.
    function fit_frozen_joint(cutoff_date; patches::Bool = false)
        o = freeze_observations(cutoff_date)
        bp = o.n - o.who_first_sitrep_days
        pp = province_increment_matrix(
            o.province_confirmed_history,
            PROVINCE_NAMES, length(PROVINCE_NAMES)
        )
        pd = province_increment_matrix(
            o.province_death_history,
            PROVINCE_NAMES, length(PROVINCE_NAMES)
        )
        patch_args = patches ?
            (;
                n_patches = length(PROVINCE_NAMES),
                province_increments = pp.increments,
                province_days = pp.days,
                province_death_increments = pd.increments,
                province_death_days = pd.days,
                province_testing_covariate = province_testing_covariate(
                    o.province_lab_daily_history
                ),
            ) : (;)
        chn = nuts_sample(
            bvd_joint(
                o.n, o.exported_cases, o.total_deaths,
                o.reported_cases, o.exports_deaths, o.confirmed_cases,
                o.tests_analysed;
                confirmed_deaths = o.confirmed_deaths,
                deaths_history = o.deaths_history,
                reported_history = o.reported_history,
                confirmed_history = o.confirmed_history,
                confirmed_deaths_history = o.confirmed_deaths_history,
                lab_history = o.lab_history,
                lab_daily_history = o.lab_daily_history,
                isolation_history = o.isolation_history,
                bed_capacity_history = o.bed_capacity_history,
                occupancy_break_days = o.occupancy_break_days,
                confirmed_break_days = o.confirmed_break_days,
                confirmed_break_gross_cases = o.confirmed_break_gross_cases,
                confirmed_break_gross_deaths = o.confirmed_break_gross_deaths,
                export_case_days = o.export_case_days,
                export_death_days = o.export_death_days,
                onset_curve_history = o.onset_curve_history,
                breakpoint = bp,
                background_re = true,
                confirmed_positivity_link = :composition,
                genetic = genetic_seeding_model,
                tmrca_days = o.tmrca_days, patch_args...
            );
            samples = samples, chains = chains, target_accept = 0.9,
            callback = fit_callback("frozen_$(cutoff_date)")
        )
        return (; cutoff = o.cutoff, o, chn)
    end

    ## A single-stream fit at the headline settings to the data frozen at
    ## `cutoff_date`, mirroring `fit_frozen_joint` but for one of the
    ## individual per-stream models below (`model_id` one of "cases",
    ## "deaths", "confirmed", "confirmed_deaths", "treatment"). Used to show
    ## each stream's own individual model alongside the joint in the
    ## one-week-back forecast validation, not just the joint alone.
    function fit_frozen_stream(model_id::AbstractString, cutoff_date)
        o = freeze_observations(cutoff_date)
        bp = o.n - o.who_first_sitrep_days
        chn = if model_id == "cases"
            nuts_sample(
                cases_only_model(
                    o.n, o.reported_cases;
                    reported_history = o.reported_history,
                    suspected_daily_history = o.suspected_daily_history,
                    breakpoint = bp
                );
                samples = samples, chains = chains,
                callback = fit_callback("frozen_$(cutoff_date)_cases")
            )
        elseif model_id == "deaths"
            nuts_sample(
                deaths_only_model(
                    o.n, o.total_deaths;
                    deaths_history = o.deaths_history,
                    suspected_daily_deaths_history =
                        o.suspected_daily_deaths_history,
                    breakpoint = bp
                );
                samples = samples, chains = chains,
                callback = fit_callback("frozen_$(cutoff_date)_deaths")
            )
        elseif model_id == "confirmed"
            nuts_sample(
                confirmed_only_model(
                    o.n, o.confirmed_cases;
                    confirmed_history = o.confirmed_history,
                    lab_history = o.lab_history,
                    lab_daily_history = o.lab_daily_history,
                    confirmed_break_days = o.confirmed_break_days,
                    confirmed_break_gross_cases = o.confirmed_break_gross_cases,
                    breakpoint = bp
                );
                samples = samples, chains = chains,
                callback = fit_callback("frozen_$(cutoff_date)_confirmed")
            )
        elseif model_id == "confirmed_deaths"
            nuts_sample(
                confirmed_deaths_only_model(
                    o.n, o.confirmed_deaths,
                    o.total_deaths;
                    deaths_history = o.deaths_history,
                    confirmed_deaths_history = o.confirmed_deaths_history,
                    confirmed_break_days = o.confirmed_break_days,
                    confirmed_break_gross_deaths =
                        o.confirmed_break_gross_deaths,
                    breakpoint = bp
                );
                samples = samples, chains = chains,
                callback = fit_callback(
                    "frozen_$(cutoff_date)_confirmed_deaths"
                )
            )
        elseif model_id == "treatment"
            nuts_sample(
                treatment_only_model(
                    o.n;
                    isolation_history = o.isolation_history,
                    bed_capacity_history = o.bed_capacity_history,
                    treatment_admissions_history =
                        o.treatment_admissions_history,
                    treatment_deaths_history = o.treatment_deaths_history,
                    treatment_ruleout_history = o.treatment_ruleout_history,
                    treatment_absconded_history =
                        o.treatment_absconded_history,
                    occupancy_break_days = o.occupancy_break_days,
                    confirmed_break_days = o.confirmed_break_days,
                    confirmed_break_gross_cases = o.confirmed_break_gross_cases,
                    breakpoint = bp
                );
                samples = samples, chains = chains,
                callback = fit_callback("frozen_$(cutoff_date)_treatment")
            )
        else
            error(
                "fit_frozen_stream: no frozen single-stream model for " *
                    "id '$model_id'"
            )
        end
        return (; cutoff = o.cutoff, o, chn)
    end

    ## One joint re-fit on the live data, with hooks to override the deaths
    ## submodel and the molecular-clock bound for the sensitivity analyses.
    function refit_joint_variant(;
            deaths = deaths_model,
            confirmed = confirmed_cases_model,
            tmrca_days = obs.tmrca_days, tmrca_days_sd = 16.0
        )
        return nuts_sample(
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
                treatment_admissions_history = obs.treatment_admissions_history,
                treatment_deaths_history = obs.treatment_deaths_history,
                treatment_ruleout_history = obs.treatment_ruleout_history,
                treatment_absconded_history = obs.treatment_absconded_history,
                occupancy_break_days = obs.occupancy_break_days,
                confirmed_break_days = obs.confirmed_break_days,
                confirmed_break_gross_cases = obs.confirmed_break_gross_cases,
                confirmed_break_gross_deaths = obs.confirmed_break_gross_deaths,
                export_case_days = obs.export_case_days,
                export_death_days = obs.export_death_days,
                onset_curve_history = obs.onset_curve_history,
                breakpoint = breakpoint,
                background_re = true,
                confirmed_positivity_link = :composition,
                deaths = deaths,
                confirmed = confirmed,
                genetic = genetic_seeding_model,
                tmrca_days = tmrca_days,
                tmrca_days_sd = tmrca_days_sd
            );
            samples = samples, chains = chains, target_accept = 0.9,
            callback = fit_callback("variant")
        )
    end

    ## Community-pathway onset-to-death delay (Isiro 2012 line-list reanalysis).
    deaths_community_delay = (history, total, onsets, k; kwargs...) -> deaths_model(
        history, total, onsets, k;
        onset_to_death = gamma_delay_model(
            40;
            alpha_prior = truncated(Normal(5.48, 2.0); lower = 0.01),
            theta_prior = truncated(Normal(1.49, 0.5); lower = 0.1)
        ),
        kwargs...
    )

    ## Exponential growth tree prior: common ancestor ~7 days earlier than
    ## the Skygrid baseline (2026-03-08 vs 2026-03-15).
    clock_alt_offset = value(Date("2026-03-08") - Date("2026-03-15"))
    tmrca_days_alt = obs.tmrca_days - clock_alt_offset

    ## Per-province spatial-table data for the patch fit, reshaped once here.
    ## It must not be built inside the model body: it looks provinces up by
    ## name in a `Dict{String}`, and a string compare on the AD tape is a
    ## `memcmp` foreigncall Mooncake has no rule for, which aborts the
    ## gradient of the whole joint.
    patch_prov = province_increment_matrix(
        obs.province_confirmed_history, PROVINCE_NAMES,
        length(PROVINCE_NAMES)
    )
    patch_prov_deaths = province_increment_matrix(
        obs.province_death_history, PROVINCE_NAMES,
        length(PROVINCE_NAMES)
    )
    patch_testing = province_testing_covariate(obs.province_lab_daily_history)

    ## The headline fit and its spatial control must differ only in the patch
    ## structure. They are the two halves of the spatial sensitivity: a gap
    ## between their C_T posteriors is read as evidence about the spatial
    ## structure, which is only meaningful if nothing else differs. Splatting
    ## one shared NamedTuple into both is what keeps them from drifting apart
    ## a keyword at a time.
    joint_common = (;
        confirmed_deaths = obs.confirmed_deaths,
        recovered_cases = obs.recovered_cases,
        deaths_history = obs.deaths_history,
        reported_history = obs.reported_history,
        confirmed_history = obs.confirmed_history,
        confirmed_deaths_history = obs.confirmed_deaths_history,
        lab_history = obs.lab_history,
        lab_daily_history = obs.lab_daily_history,
        suspected_daily_history = obs.suspected_daily_history,
        suspected_daily_deaths_history = obs.suspected_daily_deaths_history,
        isolation_history = obs.isolation_history,
        bed_capacity_history = obs.bed_capacity_history,
        recovered_history = obs.recovered_history,
        treatment_admissions_history = obs.treatment_admissions_history,
        treatment_deaths_history = obs.treatment_deaths_history,
        treatment_ruleout_history = obs.treatment_ruleout_history,
        treatment_absconded_history = obs.treatment_absconded_history,
        treatment_confirmed_incare_history =
            obs.treatment_confirmed_incare_history,
        treatment_suspect_incare_history =
            obs.treatment_suspect_incare_history,
        occupancy_break_days = obs.occupancy_break_days,
        confirmed_break_days = obs.confirmed_break_days,
        confirmed_break_gross_cases = obs.confirmed_break_gross_cases,
        confirmed_break_gross_deaths = obs.confirmed_break_gross_deaths,
        export_case_days = obs.export_case_days,
        export_death_days = obs.export_death_days,
        onset_curve_history = obs.onset_curve_history,
        breakpoint = breakpoint,
        background_re = true,
        confirmed_positivity_link = :composition,
        genetic = genetic_seeding_model,
        tmrca_days = obs.tmrca_days,
    )

    ## The only difference between the headline and the control.
    patch_only = (;
        n_patches = length(PROVINCE_NAMES),
        province_increments = patch_prov.increments,
        province_days = patch_prov.days,
        province_death_increments = patch_prov_deaths.increments,
        province_death_days = patch_prov_deaths.days,
        province_testing_covariate = patch_testing,
    )

    specs = Any[
        ## Headline fit. The patch (meta-population) model is the joint. With
        ## `n_patches = 1` it collapses exactly onto the single-population
        ## model (the sum-to-zero deviations vanish, no importation, no
        ## composition terms), so there is one model rather than two. The
        ## headline runs it over the three affected provinces.
        (;
            id = "joint",
            kind = :chain, needs = String[],
            thunk = () -> nuts_sample(
                bvd_joint(
                    obs.n, obs.exported_cases, obs.total_deaths,
                    obs.reported_cases, obs.exports_deaths,
                    obs.confirmed_cases, obs.tests_analysed;
                    joint_common..., patch_only...
                );
                samples = joint_samples(800), chains = chains,
                n_adapts = joint_warmup(500),
                target_accept = joint_target_accept(),
                callback = fit_callback("joint")
            ),
        ),
        ## Sensitivity: the same model with the spatial structure turned off
        ## (`n_patches` defaults to 1). Splitting the country into provinces
        ## adds no national data, so the two C_T posteriors should agree; a
        ## gap is a defect in the spatial structure, not a finding about it.
        ## The provinces run free and the national trajectory is their sum,
        ## so the two are not identical by construction. The seed is
        ## partitioned across patches and importation is a transfer rather
        ## than a source, so with no deviations the two match exactly. What is
        ## left between them is the country running at the force-weighted mean
        ## of the provincial Rts rather than at the trend they are centred on;
        ## test/test_patch_model.jl pins the size of that. It runs at the headline's draw
        ## count, not the matrix one, so the comparison is like for like.
        (;
            id = "sens_no_patches",
            kind = :chain, needs = String[],
            thunk = () -> nuts_sample(
                bvd_joint(
                    obs.n, obs.exported_cases, obs.total_deaths,
                    obs.reported_cases, obs.exports_deaths,
                    obs.confirmed_cases, obs.tests_analysed;
                    joint_common...
                );
                samples = joint_samples(800), chains = chains,
                n_adapts = joint_warmup(500),
                target_accept = joint_target_accept(),
                callback = fit_callback("sens_no_patches")
            ),
        ),
        (;
            id = "exports",
            kind = :chain, needs = String[],
            thunk = () -> nuts_sample(
                exports_joint_only_model(
                    obs.n, obs.exported_cases,
                    obs.exports_deaths;
                    export_case_days = obs.export_case_days,
                    export_death_days = obs.export_death_days,
                    breakpoint = breakpoint
                );
                samples = samples, chains = chains,
                check_model = false, callback = fit_callback("exports")
            ),
        ),
        (;
            id = "deaths",
            kind = :chain, needs = String[],
            thunk = () -> nuts_sample(
                deaths_only_model(
                    obs.n, obs.total_deaths;
                    deaths_history = obs.deaths_history,
                    suspected_daily_deaths_history =
                        obs.suspected_daily_deaths_history,
                    breakpoint = breakpoint
                );
                samples = samples, chains = chains,
                callback = fit_callback("deaths")
            ),
        ),
        (;
            id = "cases",
            kind = :chain, needs = String[],
            thunk = () -> nuts_sample(
                cases_only_model(
                    obs.n, obs.reported_cases;
                    reported_history = obs.reported_history,
                    suspected_daily_history = obs.suspected_daily_history,
                    breakpoint = breakpoint
                );
                samples = samples, chains = chains,
                callback = fit_callback("cases")
            ),
        ),
        (;
            id = "confirmed",
            kind = :chain, needs = String[],
            thunk = () -> nuts_sample(
                confirmed_only_model(
                    obs.n, obs.confirmed_cases;
                    confirmed_history = obs.confirmed_history,
                    lab_history = obs.lab_history,
                    lab_daily_history = obs.lab_daily_history,
                    confirmed_break_days = obs.confirmed_break_days,
                    confirmed_break_gross_cases =
                        obs.confirmed_break_gross_cases,
                    breakpoint = breakpoint
                );
                samples = samples, chains = chains,
                callback = fit_callback("confirmed")
            ),
        ),
        (;
            id = "confirmed_deaths",
            kind = :chain, needs = String[],
            thunk = () -> nuts_sample(
                confirmed_deaths_only_model(
                    obs.n, obs.confirmed_deaths,
                    obs.total_deaths;
                    deaths_history = obs.deaths_history,
                    confirmed_deaths_history = obs.confirmed_deaths_history,
                    confirmed_break_days = obs.confirmed_break_days,
                    confirmed_break_gross_deaths =
                        obs.confirmed_break_gross_deaths,
                    breakpoint = breakpoint
                );
                samples = samples, chains = chains,
                callback = fit_callback("confirmed_deaths")
            ),
        ),
        (;
            id = "treatment",
            kind = :chain, needs = String[],
            thunk = () -> nuts_sample(
                treatment_only_model(
                    obs.n;
                    isolation_history = obs.isolation_history,
                    bed_capacity_history = obs.bed_capacity_history,
                    treatment_admissions_history =
                        obs.treatment_admissions_history,
                    treatment_deaths_history = obs.treatment_deaths_history,
                    treatment_ruleout_history = obs.treatment_ruleout_history,
                    treatment_absconded_history =
                        obs.treatment_absconded_history,
                    occupancy_break_days = obs.occupancy_break_days,
                    confirmed_break_days = obs.confirmed_break_days,
                    confirmed_break_gross_cases =
                        obs.confirmed_break_gross_cases,
                    breakpoint = breakpoint
                );
                samples = samples, chains = chains,
                callback = fit_callback("treatment")
            ),
        ),
        (;
            id = "onsets",
            kind = :chain, needs = String[],
            thunk = () -> nuts_sample(
                onsets_only_model(
                    obs.n;
                    onset_curve_history = obs.onset_curve_history,
                    breakpoint = breakpoint
                );
                samples = samples, chains = chains,
                callback = fit_callback("onsets")
            ),
        ),
        (;
            id = "frozen_validation", kind = :frozen, needs = String[],
            thunk = () -> fit_frozen_joint(
                validation_cutoff;
                patches = true
            ),
        ),
    ]
    ## One frozen individual fit per still-reported stream at the validation
    ## cut-off, so the "last week versus now" forecast validation can show
    ## each stream's own model alongside the frozen joint above, not the
    ## joint alone. `validation_stream_ids` names which single-stream models
    ## get one; "exports" is excluded (not forecast at all, see
    ## `forecast_reported`) and there is no individual model for
    ## "recovered".
    for sid in validation_stream_ids(obs)
        push!(
            specs,
            (;
                id = "frozen_validation_$sid", kind = :frozen, needs = String[],
                thunk = () -> fit_frozen_stream(sid, validation_cutoff),
            )
        )
    end
    for c in frozen_cutoffs
        push!(
            specs, (;
                id = "frozen_$c", kind = :frozen, needs = String[],
                thunk = () -> fit_frozen_joint(c),
            )
        )
    end
    ## Chamla's 8 June anchor, kept out of the McCabe-matched `frozen_cutoffs`;
    ## the estimate-evolution overlay pulls it in explicitly.
    push!(
        specs,
        (;
            id = "frozen_$chamla_cutoff", kind = :frozen, needs = String[],
            thunk = () -> fit_frozen_joint(chamla_cutoff),
        )
    )
    if run_sensitivity
        push!(
            specs,
            (;
                id = "sens_community_delay", kind = :chain, needs = String[],
                thunk = () -> refit_joint_variant(deaths = deaths_community_delay),
            ),
            (;
                id = "sens_exp_growth_clock", kind = :chain, needs = String[],
                thunk = () -> refit_joint_variant(
                    tmrca_days = tmrca_days_alt, tmrca_days_sd = 16.0
                ),
            )
        )
    end

    ## The health-zone fits meld from a cached joint chain. The parent is
    ## loaded strictly: a missing parent is an error, not a refit.
    function load_parent(parent_id)
        parent = specs[findfirst(s -> s.id == parent_id, specs)]
        return fit_or_load(
            fit_key(parent_id), parent.thunk;
            cache_dir = cache_dir, strict = true
        )
    end
    ## `variant` carries the zone model's own switches (`mixing`, `deaths`,
    ## `parent_summary`) for the sensitivity re-fits below.
    function fit_zone_from(parent_chn, o, name; variant...)
        fitter = zone_fitter === nothing ? default_zone_fitter() : zone_fitter
        return fitter(
            parent_chn, o;
            samples = zone_samples(600), chains = chains,
            n_adapts = zone_warmup(400),
            target_accept = zone_target_accept(),
            max_depth = zone_max_depth(),
            callback = fit_callback(name), variant...
        )
    end
    ## A dependent key shares the parent's content hash but says nothing
    ## about the parent chain's bytes, so a parent refit under an unchanged
    ## key pairs with the zone chain already cached against the old parent.
    ## The fixed per-chain seed keeps such a refit near-identical, so the
    ## pairing holds.
    push!(
        specs,
        ## The health-zone fit on the current data, melded from the headline.
        (;
            id = "local", kind = :chain, needs = ["joint"],
            thunk = () -> fit_zone_from(load_parent("joint"), obs, "local"),
        ),
        ## The same fit melded from the validation joint, on the observations
        ## that joint was fitted to, returned in the frozen shape.
        (;
            id = "local_frozen_validation", kind = :frozen,
            needs = ["frozen_validation"],
            thunk = () -> begin
                parent = load_parent("frozen_validation")
                chn = fit_zone_from(
                    parent.chn, parent.o,
                    "local_frozen_validation"
                )
                (; cutoff = parent.cutoff, o = parent.o, chn)
            end,
        )
    )
    ## Zone-model sensitivity re-fits, each one switch away from `local`:
    ## between-zone mixing on, the zone deaths stream on, and the parent
    ## summary taken from a low or high draw rather than the mean.
    if run_sensitivity
        for (id, variant) in (
                ("local_mixing", (; mixing = true)),
                ("local_parent_low", (; parent_summary = :draw_low)),
                ("local_parent_high", (; parent_summary = :draw_high)),
            )
            push!(
                specs,
                (;
                    id, kind = :chain, needs = ["joint"],
                    thunk = () -> fit_zone_from(
                        load_parent("joint"), obs, id;
                        variant...
                    ),
                )
            )
        end
    end
    validate_fit_specs(specs)
    return specs
end

"""
    validate_fit_specs(specs)

Check that every id in a spec's `needs` names a spec listed earlier. Throws
otherwise.
"""
function validate_fit_specs(specs)
    seen = Set{String}()
    for s in specs
        for parent in s.needs
            parent in seen || error(
                "fit '$(s.id)' needs '$parent', which " *
                    "is not listed before it in the registry"
            )
        end
        push!(seen, s.id)
    end
    return specs
end

"The specs with no `needs`: the first stage, run before anything else."
base_fit_specs(specs) = [s for s in specs if isempty(s.needs)]
"The specs that need a cached parent: the second stage."
dependent_fit_specs(specs) = [s for s in specs if !isempty(s.needs)]

const FIT_STAGES = (:all, :base, :dependent)

"""
    fit_stage_env(default) -> Symbol

The fit stage named by `BVD_FIT_STAGE` (`all`, `base` or `dependent`), or
`default` when unset.
"""
function fit_stage_env(default::Symbol)
    raw = lowercase(strip(get(ENV, "BVD_FIT_STAGE", "")))
    isempty(raw) && return default
    stage = Symbol(raw)
    stage in FIT_STAGES ||
        error(
        "BVD_FIT_STAGE=$raw; expected one of " *
            join(FIT_STAGES, ", ")
    )
    return stage
end

"The specs in `stage` (`:all`, `:base` or `:dependent`), in registry order."
function stage_fit_specs(specs, stage::Symbol)
    stage === :all && return specs
    stage === :base && return base_fit_specs(specs)
    stage === :dependent && return dependent_fit_specs(specs)
    error(
        "unknown fit stage $stage; expected one of " *
            join(FIT_STAGES, ", ")
    )
end

"Ordered fit ids for the current data, sensitivity setting and `stage`."
function fit_ids(
        obs = load_observations();
        run_sensitivity = run_sensitivity_env(), stage::Symbol = :all
    )
    specs = build_fit_specs(obs; run_sensitivity)
    return [s.id for s in stage_fit_specs(specs, stage)]
end
"Ids of the base fits, those with no parent."
function base_fit_ids(
        obs = load_observations();
        run_sensitivity = run_sensitivity_env()
    )
    return fit_ids(obs; run_sensitivity, stage = :base)
end
"Ids of the dependent fits, those melded from a cached parent."
function dependent_fit_ids(
        obs = load_observations();
        run_sensitivity = run_sensitivity_env()
    )
    return fit_ids(obs; run_sensitivity, stage = :dependent)
end
