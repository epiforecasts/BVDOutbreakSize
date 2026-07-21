# Registry of the expensive model fits in the analysis report. Each fit is
# defined once here as an `(id, kind, thunk)` entry, so it can be run and
# cached independently — one per CI matrix job, or an HPC task — and the docs
# build then loads the chains through the content-addressed cache instead of
# refitting them inline. `build_fit_specs` mirrors the model calls in
# `docs/examples/analysis.jl`; keep the two in step.

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
## reason. `rt_by_release_by_stream.csv` and `size_by_release_by_stream.csv`
## are the per-fit versions of the same overlay, rewritten the same way.
## Every file score_releases.jl writes into data/ MUST be listed here, or the
## render's data hash diverges from the fit matrix's and every fit misses.
const FIT_DATA_EXCLUDE = ("released_estimates.csv",
    "rt_by_release.csv", "r0_by_release.csv",
    "forecast_scores.csv", "forecast_scores_frozen.csv",
    "forecast_overlay.csv", "forecast_overlay_frozen.csv",
    "rt_by_release_by_stream.csv", "size_by_release_by_stream.csv")

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
                export_case_days = o.export_case_days,
                export_death_days = o.export_death_days,
                breakpoint = bp,
                background_re = true,
                confirmed_positivity_link = :composition,
                genetic = genetic_seeding_model,
                tmrca_days = o.tmrca_days);
            samples = samples, chains = chains, target_accept = 0.95,
            callback = fit_callback("frozen_$(cutoff_date)"))
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
                export_case_days = obs.export_case_days,
                export_death_days = obs.export_death_days,
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

    specs = Any[
        (; id = "joint",
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
                    treatment_admissions_history = obs.treatment_admissions_history,
                    treatment_deaths_history = obs.treatment_deaths_history,
                    treatment_ruleout_history = obs.treatment_ruleout_history,
                    treatment_absconded_history = obs.treatment_absconded_history,
                    treatment_confirmed_incare_history =
                    obs.treatment_confirmed_incare_history,
                    treatment_suspect_incare_history =
                    obs.treatment_suspect_incare_history,
                    occupancy_break_days = obs.occupancy_break_days,
                    export_case_days = obs.export_case_days,
                    export_death_days = obs.export_death_days,
                    breakpoint = breakpoint,
                    background_re = true,
                    confirmed_positivity_link = :composition,
                    genetic = genetic_seeding_model,
                    tmrca_days = obs.tmrca_days);
                samples = samples, chains = chains, target_accept = 0.90,
                callback = fit_callback("joint"))),
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
                    breakpoint = breakpoint);
                samples = samples, chains = chains,
                callback = fit_callback("confirmed_deaths"))),
        (; id = "treatment",
            kind = :chain,
            thunk = () -> nuts_sample(
                treatment_only_model(obs.n;
                    isolation_history = obs.isolation_history,
                    bed_capacity_history = obs.bed_capacity_history,
                    treatment_admissions_history = obs.treatment_admissions_history,
                    treatment_deaths_history = obs.treatment_deaths_history,
                    treatment_ruleout_history = obs.treatment_ruleout_history,
                    treatment_absconded_history = obs.treatment_absconded_history,
                    occupancy_break_days = obs.occupancy_break_days,
                    breakpoint = breakpoint);
                samples = samples, chains = chains,
                callback = fit_callback("treatment"))),
        (; id = "frozen_validation", kind = :frozen,
            thunk = () -> fit_frozen_joint(validation_cutoff))
    ]
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
function fit_ids(obs = load_observations(); run_sensitivity = run_sensitivity_env())
    [s.id for s in build_fit_specs(obs; run_sensitivity)]
end
