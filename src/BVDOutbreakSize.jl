module BVDOutbreakSize

using Statistics: quantile, mean, cor
using TOML: TOML
using DataFrames: DataFrame, rename
using Chain: @chain
using Random: MersenneTwister
using Dates: Date, Day, date2epochdays, epochdays2date
using ADTypes: AutoMooncake
using Mooncake: Mooncake
using ChainRulesCore: ChainRulesCore
using Turing: @model, @addlogprob!, MCMCThreads, NUTS, sample, to_submodel
using Turing.DynamicPPL: InitFromPrior
import AbstractMCMC
import FlexiChains
using DocStringExtensions: @template, DOCSTRING, EXPORTS, IMPORTS, TYPEDEF,
                           TYPEDFIELDS, TYPEDSIGNATURES
using Distributions: Distribution, pdf, cdf, Poisson,
                     NegativeBinomial, BetaBinomial, Normal,
                     LogNormal, Beta,
                     Gamma, truncated, censored, product_distribution
using CensoredDistributions: double_interval_censored
using StatsFuns: logit, logistic
using ScoringRules: crps
import CairoMakie
import AlgebraOfGraphics as AoG
import PairPlots
using CairoMakie: Figure, Axis, hist!, density!, vlines!, hlines!, vspan!,
                  lines!, scatter!, band!, linesegments!

export REPORT_SCENARIOS, REPORT_SCENARIOS_CI,
       CHAMLA_CONFIRMED_CENTRAL, CHAMLA_CONFIRMED_W12,
       ITURI_POPULATION, ITURI_DAILY_TRAVEL,
       ITURI_DAILY_TRAVEL_SD, RENEWAL_START_LEAD, RT_WALK_LEAD,
       RT_INTERVENTION_RAMP,
       load_observations, freeze_observations, m_prior_centre,
       summary_table, posterior_summary, markdown_table,
       fit_diagnostics, diagnostics_table,
       streams_table, comparison_table,
       bias_sample, stream_calibration, onsets_over_time,
       crps_sample, log_crps_sample, score_draws, forecast_score_summary,
       nuts_sample, fit_parallel, default_adtype, enzyme_adtype,
       progress_callback, tensorboard_callback,
       combined_callback, fit_callback,
       plot_cumulative_cases, plot_cumulative_trajectories,
       plot_stream_trajectories,
       plot_density_overlay, plot_prior_predictive,
       plot_posterior_predictive, plot_posterior_predictive_grid,
       plot_pair, plot_start_date_pair, plot_estimate_comparison,
       plot_correlation_heatmap, plot_stream_pairs,
       plot_estimate_evolution, plot_evolution_by_group,
       plot_forecast_overlay,
       plot_projection_comparison,
       plot_scenario_comparison,
       plot_cfr_prior, plot_vintage_conditional_ppc,
       plot_vintage_incidence_ppc, plot_stream_calibration,
       plot_rt, plot_rt_streams,
       reconstruct_rt,
       predict_no_onward_deaths, plot_no_onward_deaths,
       forecast_reported, forecast_stream, forecast_table, forecast_archive,
       plot_forecast,
       plot_forecast_latent, plot_forecast_beds, plot_forecast_flows,
       forecast_vs_truth, forecast_vs_truth_trajectory,
       plot_forecast_vs_truth, plot_forecast_vs_truth_latent,
       plot_forecast_beds_vs_truth,
       delay_corrected_cfr, delay_corrected_confirmed_cfr,
       confirmed_cfr_table, plot_confirmed_cfr,
# renewal helpers
       renewal_infections, convolve_delay, convolve_survival, convolve_pmf,
       discretise_censored,
       euler_lotka_r, r_to_R0, doubling_time, seed_infections,
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
       surveillance_dispersion_model, pooled_dispersion_model,
       independent_ascertainment_model, pooled_ascertainment_model,
# observation models
       deaths_model, reported_cases_model, confirmed_cases_model,
       confirmed_positivity_windows, confirmed_deaths_model,
       treatment_flow_model, recovered_model,
       exports_model, exports_deaths_model,
# joint composers
       exports_only_model, deaths_only_model, cases_only_model,
       confirmed_only_model, confirmed_deaths_only_model,
       treatment_only_model,
       exports_deaths_only_model, exports_joint_only_model, bvd_joint

include("docstrings.jl")
include("constants.jl")
include("data.jl")
include("sampling.jl")
include("renewal.jl")
include("summaries.jl")
include("scoring.jl")
include("counterfactual.jl")
include("forecast.jl")
include("confirmed_cfr.jl")
include("plots.jl")
include("models/priors.jl")
include("models/observations.jl")
include("models/joint.jl")
include("precompile.jl")

end # module
