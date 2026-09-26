# Registry of the expensive model fits in the analysis report. Each fit is
# defined once here as an `(id, kind, thunk)` entry, so it can be run and
# cached independently — one per CI matrix job, or an HPC task — and the
# docs build then loads the chains through the content-addressed cache
# instead of refitting them inline. `build_fit_specs` mirrors the model
# calls in `docs/examples/analysis.jl`; keep the two in step.
#
# A health-zone fit adds `needs`, the ids it is melded from, and runs in a
# second stage once those are cached (`base_fit_specs`,
# `dependent_fit_specs`, and the `zone` group `fit_group` puts it in).

include(joinpath(@__DIR__, "cache.jl"))

using BVDOutbreakSize
using Dates: Date, Day, value
using Distributions: truncated, Normal

const _PKG = pkgdir(BVDOutbreakSize)

## Source files whose contents define the fits: the model and its submodels,
## the renewal maths, the hand-written Mooncake rules, the sampler, the data
## pipeline and this registry. A change to any of them invalidates every
## cached fit; plotting and reporting code (plots.jl, summaries.jl, ...)
## deliberately does not. The Mooncake rules are here because they change the
## floating-point gradients, and so the sampled chain. `cache.jl` is
## here because it defines what the key covers: a change to the hashing rule
## that left the key alone would give two different rules the same key.
const FIT_SOURCE_FILES = [
    joinpath(_PKG, "src", "models", "priors.jl"),
    joinpath(_PKG, "src", "models", "observations.jl"),
    joinpath(_PKG, "src", "models", "observation_distributions.jl"),
    joinpath(_PKG, "src", "models", "joint.jl"),
    joinpath(_PKG, "src", "models", "zone.jl"),
    joinpath(_PKG, "src", "models", "fit_args.jl"),
    joinpath(_PKG, "src", "renewal.jl"),
    joinpath(_PKG, "src", "sum_to_zero.jl"),
    joinpath(_PKG, "src", "mooncake_rules.jl"),
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

"""
Content hash of the fit-relevant source, data and sampler settings.
`sampler` is appended to the settings string when it is not empty.
"""
function fit_content_hash(;
        samples::Integer = 500, chains::Integer = 2,
        sampler::NamedTuple = (;)
    )
    extra = string(FIT_CACHE_SCHEMA, ":", samples, "x", chains)
    if !isempty(sampler)
        extra *= string(":", sampler)
    end
    return content_hash(
        FIT_SOURCE_FILES;
        data_dir = joinpath(_PKG, "data"),
        data_exclude = FIT_DATA_EXCLUDE,
        extra = extra
    )
end

"""
Content-addressed cache key for fit `id` at the given sampler settings.
The fits in `JOINT_SAMPLER_FITS` are also keyed on `joint_sampler_args()`
and those in `ZONE_SAMPLER_FITS` on `zone_sampler_args()`, so a run with the
`BVD_JOINT_*` overrides set writes its own cache entry.

The health-zone fits are keyed this way because they sample at the joint's
settings: without it an override moves the joint's key while leaving theirs
alone, and a zone chain melded from one joint is served for another.
"""
function fit_key(id; samples::Integer = 500, chains::Integer = 2)
    sampler = if id in JOINT_SAMPLER_FITS
        joint_sampler_args()
    elseif id in ZONE_SAMPLER_FITS
        zone_sampler_args()
    else
        (;)
    end
    return string(id, "__", fit_content_hash(; samples, chains, sampler))
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
## agree on the validation cut-off and the frozen cut-offs. The breakpoint
## comes from the package (`default_breakpoint`), alongside the fit keywords
## `src/precompile.jl` shares, so a copy here cannot shadow it and drift.
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
    fit_needs(spec) -> Vector{String}

The ids `spec` loads from the cache before it runs. Only a dependent fit
carries a `needs` field, so a base spec answers with an empty list and needs
no entry of its own.
"""
fit_needs(spec) = get(spec, :needs, String[])

"True for a fit melded from a cached parent chain."
is_dependent_fit(spec) = !isempty(fit_needs(spec))

"""
    fit_spec(id, model, sample)

A `:chain` fit spec whose thunk samples the model `model()` builds, with
`sample` taking that model to a chain. The spec carries `model` too, so a
forecast from the fit rebuilds exactly the model the fit sampled.
"""
fit_spec(id, model, sample) = (;
    id, kind = :chain, model, thunk = () -> sample(model()),
)

## Sampler settings for the fits in `JOINT_SAMPLER_FITS`.
##
## `BVD_JOINT_SAMPLES`, `BVD_JOINT_WARMUP`, `BVD_JOINT_TARGET_ACCEPT` and
## `BVD_JOINT_MAX_DEPTH` override all four without editing this file.
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

joint_max_depth() = parse(
    Int,
    get(ENV, "BVD_JOINT_MAX_DEPTH", "10")
)

## The fits that splat `joint_sampler_args()`: the headline, its spatial
## control, the one-week-back validation joint and the two sensitivity
## re-fits of the joint.
const JOINT_SAMPLER_FITS = (
    "joint", "sens_no_patches", "frozen_validation",
    "sens_community_delay", "sens_exp_growth_clock",
)

## The sampler budget every fit in `JOINT_SAMPLER_FITS` splats.
joint_sampler_args() = (;
    samples = joint_samples(1000), n_adapts = joint_warmup(500),
    target_accept = joint_target_accept(), max_depth = joint_max_depth(),
)

"""
Keyword arguments of the headline joint fit for `obs`, with `overrides`
merged over them. The headline and every sensitivity re-fit of it are built
from this, so a re-fit differs from the headline only in its overrides.
"""
function headline_joint_args(
        obs; breakpoint = default_breakpoint(obs), overrides...
    )
    return merge(
        joint_fit_args(obs; breakpoint), patch_fit_args(obs),
        values(overrides)
    )
end

"""
The keywords each sensitivity re-fit overrides on the headline joint, as a
`NamedTuple` keyed by fit id.
"""
function sensitivity_overrides(obs)
    ## Community-pathway onset-to-death delay (Isiro 2012 line-list
    ## reanalysis).
    function deaths_community_delay(history, total, onsets, k; kwargs...)
        return deaths_model(
            history, total, onsets, k;
            onset_to_death = gamma_delay_model(
                40;
                alpha_prior = truncated(Normal(5.48, 2.0); lower = 0.01),
                theta_prior = truncated(Normal(1.49, 0.5); lower = 0.1)
            ),
            kwargs...
        )
    end
    ## Exponential growth tree prior: common ancestor ~7 days earlier than
    ## the Skygrid baseline (2026-03-08 vs 2026-03-15).
    clock_alt_offset = value(Date("2026-03-08") - Date("2026-03-15"))
    return (;
        sens_community_delay = (; deaths = deaths_community_delay),
        sens_exp_growth_clock = (;
            tmrca_days = obs.tmrca_days - clock_alt_offset,
        ),
    )
end

## The fits that splat `zone_sampler_args()`.
const ZONE_SAMPLER_FITS = ("local", "local_frozen_validation")

## The sampler budget the health-zone stage splats. It follows the joint's
## overrides, so it is keyed like the joint's.
zone_sampler_args() = (;
    samples = joint_samples(800), joint_sampler_args().n_adapts,
    joint_sampler_args().target_accept,
)

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

"""
    build_fit_specs(obs; breakpoint, frozen_cutoffs, validation_cutoff,
                    run_sensitivity, samples = 500, chains = 2,
                    zone_fitter = nothing, cache_dir = fit_cache_dir())

Ordered list of the report's fits as `(; id, kind, model, thunk)` named
tuples. `kind` is `:chain` for the headline joint and single-stream fits or
`:frozen` for the frozen/validation joints (whose thunk returns
`(; cutoff, o, chn)`). `model` builds the model the thunk samples, which the
forecasts run past the cut-off. The two re-fits appended only when
`run_sensitivity` is true are not forecast and carry no `model`.

The health-zone fits `local` and `local_frozen_validation` add a `needs`
field naming the fit they meld from, which they load from `cache_dir`
strictly: a missing parent is an error, not a refit. They are the dependent
stage, and carry no `model`. `zone_fitter` fits a zone model from a parent
chain, `BVDOutbreakSize.fit_zone` when `nothing`, resolved when the thunk
runs.
"""
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
    function frozen_joint(cutoff_date; patches::Bool = false)
        o = freeze_observations(cutoff_date)
        bp = o.n - o.who_first_sitrep_days
        patch_args = patches ? patch_fit_args(o) : (;)
        model = bvd_joint(
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
            background_pooling = background_pooling_model,
            genetic = genetic_seeding_model,
            tmrca_days = o.tmrca_days, patch_args...
        )
        return (; o, model)
    end
    ## The patched validation joint takes the headline's sampler budget. The
    ## single-population frozen fits keep the smaller default.
    function fit_frozen_joint(cutoff_date; patches::Bool = false)
        f = frozen_joint(cutoff_date; patches)
        budget = patches ? joint_sampler_args() :
            (; samples = samples, target_accept = 0.9)
        chn = nuts_sample(
            f.model;
            budget..., chains = chains,
            callback = fit_callback("frozen_$(cutoff_date)")
        )
        return (; cutoff = f.o.cutoff, f.o, chn)
    end

    ## A single-stream fit at the headline settings to the data frozen at
    ## `cutoff_date`, mirroring `fit_frozen_joint` but for one of the
    ## individual per-stream models below (`model_id` one of "cases",
    ## "deaths", "confirmed", "confirmed_deaths", "treatment"). Used to show
    ## each stream's own individual model alongside the joint in the
    ## one-week-back forecast validation, not just the joint alone.
    function frozen_stream(model_id::AbstractString, cutoff_date)
        o = freeze_observations(cutoff_date)
        bp = o.n - o.who_first_sitrep_days
        model = if model_id == "cases"
            cases_only_model(
                o.n, o.reported_cases;
                reported_history = o.reported_history,
                suspected_daily_history = o.suspected_daily_history,
                breakpoint = bp
            )
        elseif model_id == "deaths"
            deaths_only_model(
                o.n, o.total_deaths;
                deaths_history = o.deaths_history,
                suspected_daily_deaths_history =
                    o.suspected_daily_deaths_history,
                breakpoint = bp
            )
        elseif model_id == "confirmed"
            confirmed_only_model(
                o.n, o.confirmed_cases;
                confirmed_history = o.confirmed_history,
                lab_history = o.lab_history,
                lab_daily_history = o.lab_daily_history,
                confirmed_break_days = o.confirmed_break_days,
                confirmed_break_gross_cases = o.confirmed_break_gross_cases,
                breakpoint = bp
            )
        elseif model_id == "confirmed_deaths"
            confirmed_deaths_only_model(
                o.n, o.confirmed_deaths,
                o.total_deaths;
                deaths_history = o.deaths_history,
                confirmed_deaths_history = o.confirmed_deaths_history,
                confirmed_break_days = o.confirmed_break_days,
                confirmed_break_gross_deaths =
                    o.confirmed_break_gross_deaths,
                breakpoint = bp
            )
        elseif model_id == "treatment"
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
            )
        else
            error(
                "fit_frozen_stream: no frozen single-stream model for " *
                    "id '$model_id'"
            )
        end
        return (; o, model)
    end
    function fit_frozen_stream(model_id::AbstractString, cutoff_date)
        f = frozen_stream(model_id, cutoff_date)
        chn = nuts_sample(
            f.model;
            samples = samples, chains = chains,
            callback = fit_callback("frozen_$(cutoff_date)_$(model_id)")
        )
        return (; cutoff = f.o.cutoff, f.o, chn)
    end

    ## One re-fit of the headline joint on the live data with the keywords
    ## `sensitivity_overrides` gives fit `id` replaced.
    function refit_joint_variant(id)
        return nuts_sample(
            bvd_joint(
                obs.n, obs.exported_cases, obs.total_deaths,
                obs.reported_cases, obs.exports_deaths, obs.confirmed_cases,
                obs.tests_analysed;
                headline_joint_args(
                    obs; breakpoint, sensitivity_overrides(obs)[Symbol(id)]...
                )...
            );
            joint_sampler_args()..., chains = chains,
            callback = fit_callback(id)
        )
    end

    ## The headline fit and its spatial control must differ only in the patch
    ## structure. They are the two halves of the spatial sensitivity: a gap
    ## between their C_T posteriors is read as evidence about the spatial
    ## structure, which is only meaningful if nothing else differs. Splatting
    ## one shared NamedTuple into both is what keeps them from drifting apart
    ## a keyword at a time.
    ##
    ## The list lives in the package (`src/models/fit_args.jl`) rather than
    ## here because `src/precompile.jl` builds its workload from the same
    ## call. Mooncake caches a reverse rule against a method signature, and a
    ## model's type carries the types of its arguments, so a workload that
    ## differs from this call in any argument type compiles a rule no fit
    ## reaches.
    joint_common = joint_fit_args(obs; breakpoint = breakpoint)

    specs = Any[
        ## Headline fit. The patch (meta-population) model is the joint. With
        ## `n_patches = 1` it collapses exactly onto the single-population
        ## model (the sum-to-zero deviations vanish, no importation, no
        ## composition terms), so there is one model rather than two. The
        ## headline runs it over the three affected provinces.
        fit_spec(
            "joint",
            () -> bvd_joint(
                obs.n, obs.exported_cases, obs.total_deaths,
                obs.reported_cases, obs.exports_deaths,
                obs.confirmed_cases, obs.tests_analysed;
                headline_joint_args(obs; breakpoint)...
            ),
            m -> nuts_sample(
                m;
                joint_sampler_args()..., chains = chains,
                callback = fit_callback("joint")
            )
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
        fit_spec(
            "sens_no_patches",
            () -> bvd_joint(
                obs.n, obs.exported_cases, obs.total_deaths,
                obs.reported_cases, obs.exports_deaths,
                obs.confirmed_cases, obs.tests_analysed;
                joint_common...
            ),
            m -> nuts_sample(
                m;
                joint_sampler_args()..., chains = chains,
                callback = fit_callback("sens_no_patches")
            )
        ),
        fit_spec(
            "exports",
            () -> exports_joint_only_model(
                obs.n, obs.exported_cases,
                obs.exports_deaths;
                export_case_days = obs.export_case_days,
                export_death_days = obs.export_death_days,
                breakpoint = breakpoint
            ),
            m -> nuts_sample(
                m;
                samples = samples, chains = chains,
                check_model = false, callback = fit_callback("exports")
            )
        ),
        fit_spec(
            "deaths",
            () -> deaths_only_model(
                obs.n, obs.total_deaths;
                deaths_history = obs.deaths_history,
                suspected_daily_deaths_history =
                    obs.suspected_daily_deaths_history,
                breakpoint = breakpoint
            ),
            m -> nuts_sample(
                m;
                samples = samples, chains = chains,
                callback = fit_callback("deaths")
            )
        ),
        fit_spec(
            "cases",
            () -> cases_only_model(
                obs.n, obs.reported_cases;
                reported_history = obs.reported_history,
                suspected_daily_history = obs.suspected_daily_history,
                breakpoint = breakpoint
            ),
            m -> nuts_sample(
                m;
                samples = samples, chains = chains,
                callback = fit_callback("cases")
            )
        ),
        fit_spec(
            "confirmed",
            () -> confirmed_only_model(
                obs.n, obs.confirmed_cases;
                confirmed_history = obs.confirmed_history,
                lab_history = obs.lab_history,
                lab_daily_history = obs.lab_daily_history,
                confirmed_break_days = obs.confirmed_break_days,
                confirmed_break_gross_cases =
                    obs.confirmed_break_gross_cases,
                breakpoint = breakpoint
            ),
            m -> nuts_sample(
                m;
                samples = samples, chains = chains,
                callback = fit_callback("confirmed")
            )
        ),
        fit_spec(
            "confirmed_deaths",
            () -> confirmed_deaths_only_model(
                obs.n, obs.confirmed_deaths,
                obs.total_deaths;
                deaths_history = obs.deaths_history,
                confirmed_deaths_history = obs.confirmed_deaths_history,
                confirmed_break_days = obs.confirmed_break_days,
                confirmed_break_gross_deaths =
                    obs.confirmed_break_gross_deaths,
                breakpoint = breakpoint
            ),
            m -> nuts_sample(
                m;
                samples = samples, chains = chains,
                callback = fit_callback("confirmed_deaths")
            )
        ),
        fit_spec(
            "treatment",
            () -> treatment_only_model(
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
            ),
            m -> nuts_sample(
                m;
                samples = samples, chains = chains,
                callback = fit_callback("treatment")
            )
        ),
        fit_spec(
            "onsets",
            () -> onsets_only_model(
                obs.n;
                onset_curve_history = obs.onset_curve_history,
                breakpoint = breakpoint
            ),
            m -> nuts_sample(
                m;
                samples = samples, chains = chains,
                callback = fit_callback("onsets")
            )
        ),
        (;
            id = "frozen_validation", kind = :frozen,
            model = () -> frozen_joint(validation_cutoff; patches = true).model,
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
                id = "frozen_validation_$sid", kind = :frozen,
                model = () -> frozen_stream(sid, validation_cutoff).model,
                thunk = () -> fit_frozen_stream(sid, validation_cutoff),
            )
        )
    end
    for c in frozen_cutoffs
        push!(
            specs, (;
                id = "frozen_$c", kind = :frozen,
                model = () -> frozen_joint(c).model,
                thunk = () -> fit_frozen_joint(c),
            )
        )
    end
    ## Chamla's 8 June anchor, kept out of the McCabe-matched `frozen_cutoffs`;
    ## the estimate-evolution overlay pulls it in explicitly.
    push!(
        specs,
        (;
            id = "frozen_$chamla_cutoff", kind = :frozen,
            model = () -> frozen_joint(chamla_cutoff).model,
            thunk = () -> fit_frozen_joint(chamla_cutoff),
        )
    )
    if run_sensitivity
        push!(
            specs,
            (
                (;
                    id, kind = :chain,
                    thunk = () -> refit_joint_variant(id),
                )
                    for id in string.(keys(sensitivity_overrides(obs)))
            )...
        )
    end

    ## The health-zone fits meld from a cached parent chain. The parent is
    ## loaded strictly: a missing parent is an error, not a refit.
    function load_parent(parent_id)
        parent = specs[findfirst(s -> s.id == parent_id, specs)]
        return fit_or_load(
            fit_key(parent_id), parent.thunk;
            cache_dir = cache_dir, strict = true
        )
    end
    ## The zone stage samples at the headline joint's settings, so a
    ## difference between the two levels is the model rather than the
    ## sampler.
    function fit_zone_from(parent_chn, o, name)
        fitter = zone_fitter === nothing ? default_zone_fitter() : zone_fitter
        return fitter(
            parent_chn, o;
            zone_sampler_args()..., chains = chains,
            callback = fit_callback(name)
        )
    end
    ## A dependent key shares the parent's content hash and moves with the
    ## joint's sampler overrides (`ZONE_SAMPLER_FITS`), but says nothing
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
    validate_fit_specs(specs)
    return specs
end

"""
    validate_fit_specs(specs)

Check that every id a spec needs names a spec listed earlier. Throws
otherwise.
"""
function validate_fit_specs(specs)
    seen = Set{String}()
    for s in specs
        for parent in fit_needs(s)
            parent in seen || error(
                "fit '$(s.id)' needs '$parent', which " *
                    "is not listed before it in the registry"
            )
        end
        push!(seen, s.id)
    end
    return specs
end

"The specs with no parent: the first stage, run before anything else."
base_fit_specs(specs) = [s for s in specs if !is_dependent_fit(s)]
"The specs that need a cached parent: the second stage."
dependent_fit_specs(specs) = [s for s in specs if is_dependent_fit(s)]

const FIT_STAGES = (:all, :base, :dependent)

## The CI jobs the fits are split across, in the order `list.jl --groups`
## prints them.
const FIT_GROUPS = ("joint", "streams", "frozen", "sensitivity", "zone")

"""
    fit_group(spec) -> String

The CI job that runs `spec`, one of `FIT_GROUPS`. `zone` is exactly
`dependent_fit_specs`, so a fit's parent decides both the local stage and
the CI job: a dependent fit runs after the group that cached its parent. The
other four partition the base fits by what a render job has to wait for.
"""
function fit_group(spec)
    is_dependent_fit(spec) && return "zone"
    spec.id in ("joint", "sens_no_patches") && return "joint"
    spec.kind === :frozen && return "frozen"
    startswith(spec.id, "sens_") && return "sensitivity"
    return "streams"
end

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
