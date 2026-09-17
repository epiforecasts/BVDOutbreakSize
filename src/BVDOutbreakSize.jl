module BVDOutbreakSize

using Statistics: quantile, mean, cor, median, cov, std
using LinearAlgebra: cholesky, Symmetric, tr, I, Diagonal, diag,
    issuccess
using Statistics: quantile, mean, cor, median, std
using TOML: TOML
using DataFrames: DataFrame, rename, select, Not, nrow
using Chain: @chain
using Random: AbstractRNG, MersenneTwister
using Dates: Date, Day, date2epochdays, epochdays2date
using ADTypes: AutoMooncake
using Mooncake: Mooncake
using ChainRulesCore: ChainRulesCore
using Turing: @model, @addlogprob!, MCMCThreads, NUTS, sample, to_submodel
using Turing.DynamicPPL: InitFromPrior, InitFromVector, InitFromParams,
    LogDensityFunction, VarInfo, getlogjoint, link
using LogDensityProblems: LogDensityProblems
import AbstractMCMC
import FlexiChains
using DocStringExtensions: @template, DOCSTRING, EXPORTS, IMPORTS, TYPEDEF,
    TYPEDFIELDS, TYPEDSIGNATURES
using Distributions: Distribution, pdf, cdf, Poisson,
    NegativeBinomial, BetaBinomial, Normal,
    LogNormal, Beta, LKJCholesky,
    Gamma, TDist, truncated, censored, product_distribution,
    DirichletMultinomial
using CensoredDistributions: double_interval_censored
using SpecialFunctions: loggamma
using StatsFuns: logit, logistic
import CairoMakie
import AlgebraOfGraphics as AoG
import PairPlots
import JSON
using CairoMakie: Figure, Axis, hist!, density!, vlines!, hlines!, vspan!,
    lines!, scatter!, band!, linesegments!, scatterlines!

export JOINT_FIT, BASELINE_FIT, FROZEN_FIT,
    REPORT_SCENARIOS, REPORT_SCENARIOS_CI,
    CHAMLA_CONFIRMED_CENTRAL, CHAMLA_CONFIRMED_W12,
    ITURI_POPULATION, ITURI_DAILY_TRAVEL,
    ITURI_DAILY_TRAVEL_SD, RENEWAL_START_LEAD, RT_WALK_LEAD,
    RT_INTERVENTION_RAMP, ONSET_REPORT_MAX_DELAY,
    load_observations, freeze_observations, m_prior_centre,
    load_onset_curve,
    OBSERVATION_STREAMS, STREAM_REPORTING_GRACE_DAYS,
    stream_id, stream_forecast_columns, stream_last_date,
    stream_reporting,
    stream_report_status,
    summary_table, posterior_summary, markdown_table, MarkdownTable,
    patch_summary_table, patch_overview_table,
    province_cfr_table, province_forecast_table,
    province_forecast_vs_truth, plot_province_forecast,
    fit_diagnostics, diagnostics_table,
    parameter_diagnostics, worst_parameters_table,
    family_diagnostics_table, diagnostic_spread_table,
    sampler_by_chain_table, divergence_location_table,
    diagnostic_contrast, diagnostic_contrast_table,
    streams_table, comparison_table,
    bias_sample, stream_calibration, onsets_over_time,
    crps_sample, log_crps_sample, crps_decomposition, score_draws,
    forecast_score_overview, forecast_score_by_horizon,
    forecast_score_by_release, forecast_score_by_vintage,
    drop_individual_fit_columns,
    drop_degenerate_fit_column, select_fit_role, scored_overlay,
    nuts_sample, fit_parallel, default_adtype, enzyme_adtype,
    ViablePrior, viable_prior_init,
    progress_callback, tensorboard_callback,
    combined_callback, fit_callback,
    plot_cumulative_cases, plot_cumulative_trajectories,
    plot_stream_trajectories,
    plot_density_overlay, plot_prior_predictive,
    plot_posterior_predictive, plot_posterior_predictive_grid,
    plot_pair, plot_start_date_pair, plot_estimate_comparison,
    plot_correlation_heatmap, plot_stream_pairs,
    plot_estimate_evolution, plot_evolution_by_group,
    plot_forecast_overlay, plot_forecast_relative_skill,
    plot_forecast_skill_by_vintage,
    plot_projection_comparison,
    plot_scenario_comparison,
    plot_cfr_prior, plot_vintage_conditional_ppc,
    plot_vintage_incidence_ppc, plot_stream_calibration,
    plot_rt, plot_rt_streams, plot_rt_patches,
    plot_infections_patches, plot_imports_patches,
    plot_patch_summary,
    plot_province_composition_ppc,
    ZONE_MAP_PROVINCES, zone_key, zone_geojson_path,
    load_health_zones_geojson,
    plot_zone_map, plot_zone_map_panels, plot_rt_zones,
    plot_zone_shares, plot_zone_forecast, plot_zone_ranking,
    plot_zone_comparison, zone_summary_table,
    plot_rhat_spread, plot_parameter_index_diagnostics,
    plot_divergence_locations, plot_diagnostic_contrast,
    reconstruct_rt, reconstruct_patch_rt, reconstruct_onset_hazard,
    onset_nowcast_draws, plot_onset_nowcast_grid,
    predict_no_onward_deaths, plot_no_onward_deaths,
    forecast_reported, forecast_stream, forecast_table, forecast_archive,
    province_forecast_archive,
    forecast_onsets, onset_forecast_table,
    plot_forecast,
    plot_forecast_latent, plot_forecast_beds, plot_forecast_flows,
    forecast_vs_truth,
    plot_forecast_vs_truth, plot_forecast_vs_truth_latent,
    plot_forecast_beds_vs_truth,
    delay_corrected_cfr, delay_corrected_confirmed_cfr,
    confirmed_cfr_table, plot_confirmed_cfr,
    # renewal helpers
    renewal_infections, convolve_delay, convolve_survival, convolve_pmf,
    discretise_censored,
    euler_lotka_r, r_to_R0, doubling_time, seed_infections,
    confirmed_break_correction,
    seed_at_renewal_start,
    knot_days,
    interpolate_knots, sigmoid_ramp, seeding_age, lognormal_meansd,
    safe_rate,
    # prior / latent submodels
    censored_delay_model, gamma_delay_model, onset_to_death_model,
    nejm_onset_to_sample,
    generation_interval_model, rt_walk_model,
    seed_model, exponential_growth_model, infection_model,
    onset_incidence_model,
    genetic_seeding_model,
    cfr_model, traveller_volume_model, test_positivity_model,
    isolation_admission_model, isolation_severity_model, bed_capacity_model,
    bed_capacity_walk_model, recovery_probability_model,
    death_background_model, death_ascertainment_model, background_cfr_model,
    background_re_model, background_pooling_model,
    background_walk_model,
    expand_vintage_rate,
    test_sensitivity_model, test_specificity_model, lab_delay_model,
    confirmed_positivity_model, confirmed_overdispersion_model,
    severity_enrichment_model,
    death_testing_fraction_model, death_testing_scaling_model,
    specimen_intensity_model,
    surveillance_dispersion_model, pooled_dispersion_model,
    independent_ascertainment_model, pooled_ascertainment_model,
    # observation models
    deaths_model, reported_cases_model, confirmed_cases_model,
    confirmed_positivity_windows, confirmed_break_offset,
    break_step_centres,
    confirmed_deaths_model,
    treatment_flow_model, recovered_model,
    exports_model, exports_deaths_model,
    safe_studentt, onset_report_cdf, onset_report_cdf_extrapolated,
    onset_report_cdf_table,
    onset_report_G, onset_report_F, onset_nowcast,
    onset_report_anchor, onset_report_anchor_series,
    onset_report_moments, onset_report_scales, onset_report_scale,
    onset_vintage_indices, onset_scan_adjust,
    onset_report_expected_total,
    onset_report_ascertainment, onset_report_hazard_model,
    onset_ascertainment_model, onset_reporting_model,
    # joint composers
    exports_only_model, deaths_only_model, cases_only_model,
    confirmed_only_model, confirmed_deaths_only_model,
    treatment_only_model,
    exports_deaths_only_model, exports_joint_only_model, bvd_joint,
    onsets_only_model,
    PROVINCE_NAMES, PROVINCE_LABELS, PROVINCE_POPULATIONS,
    PROVINCE_CAPITALS, PROVINCE_MEMBERS,
    PROVINCE_SOURCE_NAMES, PROVINCE_SOURCE_POPULATIONS,
    PROVINCE_SOURCE_CAPITALS,
    PROVINCE_DISTANCE_DECAY, haversine_km,
    province_distance_matrix, province_importation_kernel,
    province_increment_matrix, province_testing_covariate,
    zone_increment_matrix, zone_reattribution_days,
    zone_cumulative_falls, load_health_zones,
    patch_infections, importation_from_kernel,
    implied_national_Rt, implied_national_Rt_at,
    patch_rt_model, patch_infection_model,
    province_export_pressure_model,
    province_composition_model,
    # health-zone model
    bvd_zone, fit_zone, zone_fit_inputs, zone_parent_inputs,
    zone_share_renewal, zone_meld_block, zone_week_midpoints,
    zone_parent_scale, zone_deformation, zone_meld_check,
    reconstruct_zone_shares, reconstruct_zone_rt, zone_infections,
    zone_forecast_shares, zone_forecast_draws, zone_forecast_archive,
    zone_overview_table,
    zone_forecast_table, zone_forecast_truth, zone_forecast_vs_truth,
    zone_forecast_scores, zone_composition_ppc, zone_composition_draws,
    zone_composition_calibration, plot_zone_composition_ppc,
    zone_diagnostics_table, zone_sampler_diagnostics

include("docstrings.jl")
include("constants.jl")
include("data.jl")
include("onset_curve.jl")
include("sampling.jl")
include("renewal.jl")
include("summaries.jl")
include("diagnostics.jl")
include("scoring.jl")
include("counterfactual.jl")
include("forecast.jl")
include("confirmed_cfr.jl")
include("plots.jl")
include("models/priors.jl")
include("models/observations.jl")
include("models/joint.jl")
include("models/zone.jl")
include("zone.jl")
include("precompile.jl")

end # module
