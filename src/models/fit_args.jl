# The keyword shape every production `bvd_joint` fit is built from.
#
# A `DynamicPPL.Model`'s type carries the types of its arguments, and that
# type is the key Mooncake caches a reverse rule against. Two calls that
# differ in the type of any keyword are different method instances and share
# no compiled rule. The precompile workload therefore only saves a fit the
# cold compile when it passes the same keyword set, with the same types, as
# the fit does.
#
# Keeping that list in one place is what makes the two impossible to drift
# apart. `docs/fits/registry.jl` splats these into every joint fit, and
# `src/precompile.jl` splats them into the workload. A keyword added here
# reaches both; a keyword added to one call site alone reaches neither.

"""
$(TYPEDSIGNATURES)

Intervention day index for a `load_observations()` result: the day the WHO
published its first situation report, counted back from the grid end.
"""
default_breakpoint(obs) = obs.n - obs.who_first_sitrep_days

"""
$(TYPEDSIGNATURES)

Keyword arguments shared by every production [`bvd_joint`](@ref) fit, built
from a `load_observations()` result.

`breakpoint` is the intervention day index.

See also [`patch_fit_args`](@ref), which adds the spatial structure.
"""
function joint_fit_args(obs; breakpoint)
    return (;
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
        background_pooling = background_pooling_model,
        genetic = genetic_seeding_model,
        tmrca_days = obs.tmrca_days,
    )
end

"""
$(TYPEDSIGNATURES)

Spatial keyword arguments for the headline patch fit, built from a
`load_observations()` result. The single-population control passes
[`joint_fit_args`](@ref) without these.

The per-province tables are reshaped here rather than in the model body: the
lookup goes by province name through a `Dict{String}`, and a string compare
on the AD tape is a `memcmp` foreigncall Mooncake has no rule for, which
aborts the gradient of the whole joint.
"""
function patch_fit_args(obs)
    prov = province_increment_matrix(
        obs.province_confirmed_history, PROVINCE_NAMES,
        length(PROVINCE_NAMES)
    )
    prov_deaths = province_increment_matrix(
        obs.province_death_history, PROVINCE_NAMES,
        length(PROVINCE_NAMES)
    )
    prov_lab = province_lab_increment_matrix(
        obs.province_lab_daily_history, PROVINCE_NAMES,
        length(PROVINCE_NAMES)
    )
    return (;
        n_patches = length(PROVINCE_NAMES),
        province_increments = prov.increments,
        province_days = prov.days,
        province_death_increments = prov_deaths.increments,
        province_death_days = prov_deaths.days,
        province_testing_covariate =
            province_testing_covariate(obs.province_lab_daily_history),
        province_lab_increments = prov_lab.increments,
        province_lab_days = prov_lab.days,
        province_isolation = province_care_observations(
            obs.province_isolation_history, PROVINCE_NAMES
        ),
        province_capacity = province_care_observations(
            obs.province_bed_capacity_history, PROVINCE_NAMES
        ),
    )
end

"""
$(TYPEDSIGNATURES)

The headline [`bvd_joint`](@ref) model, for a `load_observations()` result.

This is the call the precompile workload compiles and the headline fit
samples, so both reach Mooncake as the same method instance.
"""
function production_joint(obs; breakpoint)
    return bvd_joint(
        obs.n, obs.exported_cases, obs.total_deaths,
        obs.reported_cases, obs.exports_deaths,
        obs.confirmed_cases, obs.tests_analysed;
        joint_fit_args(obs; breakpoint)...,
        patch_fit_args(obs)...
    )
end
