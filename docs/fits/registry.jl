# Registry of the expensive model fits in the analysis report. Each fit is
# defined once here as an `(id, kind, thunk)` entry, so it can be run and
# cached independently — one per CI matrix job, or an HPC task — and the
# docs build then loads the chains through the content-addressed cache
# instead of refitting them inline. `build_fit_specs` mirrors the model
# calls in `docs/examples/analysis.jl`; keep the two in step.

include(joinpath(@__DIR__, "cache.jl"))

using BVDOutbreakSize
using Dates: Date, Day, value
using Distributions: truncated, Normal

const _PKG = pkgdir(BVDOutbreakSize)

## Source files whose contents define the fits: the model and its submodels,
## the renewal maths, the sampler, the data pipeline and this registry. A
## change to any of them invalidates every cached fit; plotting and reporting
## code (plots.jl, summaries.jl, ...) deliberately does not.
const FIT_SOURCE_FILES = [
    joinpath(_PKG, "src", "models", "priors.jl"),
    joinpath(_PKG, "src", "models", "observations.jl"),
    joinpath(_PKG, "src", "models", "joint.jl"),
    joinpath(_PKG, "src", "renewal.jl"),
    joinpath(_PKG, "src", "sampling.jl"),
    joinpath(_PKG, "src", "constants.jl"),
    joinpath(_PKG, "src", "data.jl"),
    joinpath(_PKG, "src", "onset_curve.jl"),
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
const FIT_DATA_EXCLUDE = ("released_estimates.csv",
    "rt_by_release.csv", "r0_by_release.csv",
    "forecast_scores.csv", "forecast_scores_frozen.csv",
    "forecast_overlay.csv", "forecast_overlay_frozen.csv",
    "rt_by_release_by_stream.csv", "size_by_release_by_stream.csv",
    "r0_by_release_by_stream.csv")

"Content hash of the fit-relevant source, data and sampler settings."
function fit_content_hash(; samples::Integer = 500, chains::Integer = 2)
    return content_hash(FIT_SOURCE_FILES;
        data_dir = joinpath(_PKG, "data"),
        data_exclude = FIT_DATA_EXCLUDE,
        extra = string(FIT_CACHE_SCHEMA, ":", samples, "x", chains))
end

"Content-addressed cache key for fit `id` at the given sampler settings."
function fit_key(id; samples::Integer = 500, chains::Integer = 2)
    string(id, "__", fit_content_hash(; samples, chains))
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
## The individual single-stream models frozen at `validation_cutoff` for
## the "last week versus now" forecast validation (see `fit_frozen_stream`
## and `docs/examples/sensitivity.jl`'s "Forecast validation" section).
## "exports" is excluded (not forecast at all) and there is no
## individual model for "recovered".
const VALIDATION_STREAM_IDS = ("cases", "deaths", "confirmed",
    "confirmed_deaths", "treatment")
function run_sensitivity_env()
    lowercase(strip(get(ENV, "BVD_RUN_SENSITIVITY", "false"))) in
    ("true", "1", "yes", "on")
end

"""
    build_fit_specs(obs; breakpoint, frozen_cutoffs, validation_cutoff,
                    run_sensitivity, samples = 500, chains = 2)

Ordered list of the report's fits as `(; id, kind, thunk)` named tuples. `kind`
is `:chain` for the headline joint and single-stream fits or `:frozen` for the
frozen/validation joints (whose thunk returns `(; cutoff, o, chn)`). The
sensitivity re-fits are appended only when `run_sensitivity` is true.
"""
function build_fit_specs(obs;
        breakpoint = default_breakpoint(obs),
        frozen_cutoffs = default_frozen_cutoffs(),
        chamla_cutoff = default_chamla_cutoff(),
        validation_cutoff = default_validation_cutoff(obs),
        run_sensitivity = run_sensitivity_env(),
        samples::Integer = 500, chains::Integer = 2)

    ## A joint fit at the headline settings to the data frozen at `cutoff_date`.
    function fit_frozen_joint(cutoff_date)
        o = freeze_observations(cutoff_date)
        bp = o.n - o.who_first_sitrep_days
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
                tmrca_days = o.tmrca_days);
            samples = samples, chains = chains, target_accept = 0.95,
            callback = fit_callback("frozen_$(cutoff_date)"))
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
                cases_only_model(o.n, o.reported_cases;
                    reported_history = o.reported_history,
                    suspected_daily_history = o.suspected_daily_history,
                    breakpoint = bp);
                samples = samples, chains = chains,
                callback = fit_callback("frozen_$(cutoff_date)_cases"))
        elseif model_id == "deaths"
            nuts_sample(
                deaths_only_model(o.n, o.total_deaths;
                    deaths_history = o.deaths_history,
                    suspected_daily_deaths_history =
                    o.suspected_daily_deaths_history,
                    breakpoint = bp);
                samples = samples, chains = chains,
                callback = fit_callback("frozen_$(cutoff_date)_deaths"))
        elseif model_id == "confirmed"
            nuts_sample(
                confirmed_only_model(o.n, o.confirmed_cases;
                    confirmed_history = o.confirmed_history,
                    lab_history = o.lab_history,
                    lab_daily_history = o.lab_daily_history,
                    confirmed_break_days = o.confirmed_break_days,
                    confirmed_break_gross_cases = o.confirmed_break_gross_cases,
                    breakpoint = bp);
                samples = samples, chains = chains,
                callback = fit_callback("frozen_$(cutoff_date)_confirmed"))
        elseif model_id == "confirmed_deaths"
            nuts_sample(
                confirmed_deaths_only_model(o.n, o.confirmed_deaths,
                    o.total_deaths;
                    deaths_history = o.deaths_history,
                    confirmed_deaths_history = o.confirmed_deaths_history,
                    confirmed_break_days = o.confirmed_break_days,
                    confirmed_break_gross_deaths =
                    o.confirmed_break_gross_deaths,
                    breakpoint = bp);
                samples = samples, chains = chains,
                callback = fit_callback(
                    "frozen_$(cutoff_date)_confirmed_deaths"))
        elseif model_id == "treatment"
            nuts_sample(
                treatment_only_model(o.n;
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
                    breakpoint = bp);
                samples = samples, chains = chains,
                callback = fit_callback("frozen_$(cutoff_date)_treatment"))
        else
            error("fit_frozen_stream: no frozen single-stream model for " *
                  "id '$model_id'")
        end
        return (; cutoff = o.cutoff, o, chn)
    end

    ## One joint re-fit on the live data, with hooks to override the deaths
    ## submodel and the molecular-clock bound for the sensitivity analyses.
    function refit_joint_variant(; deaths = deaths_model,
            confirmed = confirmed_cases_model,
            tmrca_days = obs.tmrca_days, tmrca_days_sd = 16.0)
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
                tmrca_days_sd = tmrca_days_sd);
            samples = samples, chains = chains, target_accept = 0.95,
            callback = fit_callback("variant"))
    end

    ## Community-pathway onset-to-death delay (Isiro 2012 line-list reanalysis).
    deaths_community_delay = (history, total, onsets, k; kwargs...) -> deaths_model(
        history, total, onsets, k;
        onset_to_death = gamma_delay_model(40;
            alpha_prior = truncated(Normal(5.48, 2.0); lower = 0.01),
            theta_prior = truncated(Normal(1.49, 0.5); lower = 0.1)),
        kwargs...)

    ## Exponential growth tree prior: common ancestor ~7 days earlier than
    ## the Skygrid baseline (2026-03-08 vs 2026-03-15).
    clock_alt_offset = value(Date("2026-03-08") - Date("2026-03-15"))
    tmrca_days_alt = obs.tmrca_days - clock_alt_offset

    ## Per-province spatial-table data for the patch fit, reshaped ONCE here.
    ## It must not be built inside the model body: it looks provinces up by
    ## name in a `Dict{String}`, and a string compare on the AD tape is a
    ## `memcmp` foreigncall Mooncake has no rule for, which aborts the
    ## gradient of the whole joint.
    patch_prov = province_increment_matrix(
        obs.province_confirmed_history, PROVINCE_NAMES, 3)
    patch_prov_deaths = province_increment_matrix(
        obs.province_death_history, PROVINCE_NAMES, 3)

    specs = Any[
        ## HEADLINE FIT. The patch (meta-population) model IS the joint: with
        ## `n_patches = 1` it collapses exactly onto the single-population
        ## model (the sum-to-zero deviations vanish, no importation, no
        ## composition terms), so there is one model rather than two. The
        ## headline runs it over the three affected provinces.
        (; id = "joint",
            kind = :chain,
            thunk = () -> nuts_sample(
                bvd_joint(obs.n, obs.exported_cases, obs.total_deaths,
                    obs.reported_cases, obs.exports_deaths,
                    obs.confirmed_cases, obs.tests_analysed;
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
                    treatment_admissions_history =
                    obs.treatment_admissions_history,
                    treatment_deaths_history = obs.treatment_deaths_history,
                    treatment_ruleout_history = obs.treatment_ruleout_history,
                    treatment_absconded_history =
                    obs.treatment_absconded_history,
                    occupancy_break_days = obs.occupancy_break_days,
                    export_case_days = obs.export_case_days,
                    export_death_days = obs.export_death_days,
                    breakpoint = breakpoint,
                    n_patches = 3,
                    province_increments = patch_prov.increments,
                    province_days = patch_prov.days,
                    province_death_increments = patch_prov_deaths.increments,
                    province_death_days = patch_prov_deaths.days,
                    tmrca_days = obs.tmrca_days);
                samples = samples, chains = chains, target_accept = 0.95,
                callback = fit_callback("joint"))),
        ## SENSITIVITY: the same model with the spatial structure turned OFF
        ## (`n_patches` defaults to 1). This is the check on the headline --
        ## if the patch structure were distorting the national fit, the two
        ## C_T posteriors would part company.
        (; id = "sens_no_patches",
            kind = :chain,
            thunk = () -> nuts_sample(
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
                    treatment_admissions_history =
                    obs.treatment_admissions_history,
                    treatment_deaths_history = obs.treatment_deaths_history,
                    treatment_ruleout_history = obs.treatment_ruleout_history,
                    treatment_absconded_history =
                    obs.treatment_absconded_history,
                    treatment_confirmed_incare_history =
                    obs.treatment_confirmed_incare_history,
                    treatment_suspect_incare_history =
                    obs.treatment_suspect_incare_history,
                    occupancy_break_days = obs.occupancy_break_days,
                    confirmed_break_days = obs.confirmed_break_days,
                    confirmed_break_gross_cases =
                    obs.confirmed_break_gross_cases,
                    confirmed_break_gross_deaths =
                    obs.confirmed_break_gross_deaths,
                    export_case_days = obs.export_case_days,
                    export_death_days = obs.export_death_days,
                    onset_curve_history = obs.onset_curve_history,
                    breakpoint = breakpoint,
                    background_re = true,
                    confirmed_positivity_link = :composition,
                    genetic = genetic_seeding_model,
                    tmrca_days = obs.tmrca_days);
<<<<<<< HEAD
                samples = samples, chains = chains, target_accept = 0.90,
                callback = fit_callback("sens_no_patches"))),
        ## The patch (meta-population) fit. Registered so that it is fitted
        ## end-to-end in CI like every other stream: the two defects that the
        ## patch model shipped with (a seed prior that forced the provincial
        ## Rt to absorb the case-split level, and an AD-breaking Dict lookup
        ## inside the model body) were invisible to the unit tests and only
        ## surfaced on a real fit. A standing fit makes the
        ## posterior-predictive check a gate rather than a manual step.
        (; id = "exports",
            kind = :chain,
            thunk = () -> nuts_sample(
                exports_joint_only_model(obs.n, obs.exported_cases,
                    obs.exports_deaths;
                    export_case_days = obs.export_case_days,
                    export_death_days = obs.export_death_days,
                    breakpoint = breakpoint);
                samples = samples, chains = chains,
                check_model = false, callback = fit_callback("exports"))),
        (; id = "deaths",
            kind = :chain,
            thunk = () -> nuts_sample(
                deaths_only_model(obs.n, obs.total_deaths;
                    deaths_history = obs.deaths_history,
                    suspected_daily_deaths_history =
                    obs.suspected_daily_deaths_history,
                    breakpoint = breakpoint);
                samples = samples, chains = chains,
                callback = fit_callback("deaths"))),
        (; id = "cases",
            kind = :chain,
            thunk = () -> nuts_sample(
                cases_only_model(obs.n, obs.reported_cases;
                    reported_history = obs.reported_history,
                    suspected_daily_history = obs.suspected_daily_history,
                    breakpoint = breakpoint);
                samples = samples, chains = chains,
                callback = fit_callback("cases"))),
        (; id = "confirmed",
            kind = :chain,
            thunk = () -> nuts_sample(
                confirmed_only_model(obs.n, obs.confirmed_cases;
                    confirmed_history = obs.confirmed_history,
                    lab_history = obs.lab_history,
                    lab_daily_history = obs.lab_daily_history,
                    confirmed_break_days = obs.confirmed_break_days,
                    confirmed_break_gross_cases =
                    obs.confirmed_break_gross_cases,
                    breakpoint = breakpoint);
                samples = samples, chains = chains,
                callback = fit_callback("confirmed"))),
        (; id = "confirmed_deaths",
            kind = :chain,
            thunk = () -> nuts_sample(
                confirmed_deaths_only_model(obs.n, obs.confirmed_deaths,
                    obs.total_deaths;
                    deaths_history = obs.deaths_history,
                    confirmed_deaths_history = obs.confirmed_deaths_history,
                    confirmed_break_days = obs.confirmed_break_days,
                    confirmed_break_gross_deaths =
                    obs.confirmed_break_gross_deaths,
                    breakpoint = breakpoint);
                samples = samples, chains = chains,
                callback = fit_callback("confirmed_deaths"))),
        (; id = "treatment",
            kind = :chain,
            thunk = () -> nuts_sample(
                treatment_only_model(obs.n;
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
                    breakpoint = breakpoint);
                samples = samples, chains = chains,
                callback = fit_callback("treatment"))),
        (; id = "onsets",
            kind = :chain,
            thunk = () -> nuts_sample(
                onsets_only_model(obs.n;
                    onset_curve_history = obs.onset_curve_history,
                    breakpoint = breakpoint);
                samples = samples, chains = chains,
                callback = fit_callback("onsets"))),
        (; id = "frozen_validation", kind = :frozen,
            thunk = () -> fit_frozen_joint(validation_cutoff))
    ]
    ## One frozen individual fit per stream at the validation cut-off, so
    ## the "last week versus now" forecast validation can show each
    ## stream's own model alongside the frozen joint above, not the joint
    ## alone. `VALIDATION_STREAM_IDS` names which single-stream models get
    ## one; "exports" is excluded (not forecast at all, see
    ## `forecast_reported`) and there is no individual model for
    ## "recovered".
    for sid in VALIDATION_STREAM_IDS
        push!(specs,
            (; id = "frozen_validation_$sid", kind = :frozen,
                thunk = () -> fit_frozen_stream(sid, validation_cutoff)))
    end
    for c in frozen_cutoffs
        push!(specs, (; id = "frozen_$c", kind = :frozen,
            thunk = () -> fit_frozen_joint(c)))
    end
    ## Chamla's 8 June anchor, kept out of the McCabe-matched `frozen_cutoffs`;
    ## the estimate-evolution overlay pulls it in explicitly.
    push!(specs,
        (; id = "frozen_$chamla_cutoff", kind = :frozen,
            thunk = () -> fit_frozen_joint(chamla_cutoff)))
    if run_sensitivity
        push!(specs,
            (; id = "sens_community_delay", kind = :chain,
                thunk = () -> refit_joint_variant(deaths = deaths_community_delay)),
            (; id = "sens_exp_growth_clock", kind = :chain,
                thunk = () -> refit_joint_variant(
                    tmrca_days = tmrca_days_alt, tmrca_days_sd = 16.0)))
    end
    return specs
end

"Ordered fit ids for the current data and sensitivity setting."
function fit_ids(
        obs = load_observations(); run_sensitivity = run_sensitivity_env())
    [s.id for s in build_fit_specs(obs; run_sensitivity)]
end
