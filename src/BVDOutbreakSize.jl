module BVDOutbreakSize

using Statistics: quantile
using TOML: TOML
using DataFrames: DataFrame, rename
using Chain: @chain
using Random: MersenneTwister
using Dates: Date, Day, date2epochdays, epochdays2date
using ADTypes: AutoMooncake
using Mooncake: Mooncake
using ChainRulesCore: ChainRulesCore
using Turing: @model, MCMCThreads, NUTS, sample, to_submodel
using Turing.DynamicPPL: InitFromPrior
import AbstractMCMC
import FlexiChains
using DocStringExtensions: @template, DOCSTRING, EXPORTS, IMPORTS, TYPEDEF,
                           TYPEDFIELDS, TYPEDSIGNATURES
using Distributions: Distribution, pdf, Poisson,
                     NegativeBinomial, Binomial, Normal, LogNormal, Beta,
                     Gamma, truncated, censored, product_distribution
using CensoredDistributions: double_interval_censored
using StatsFuns: logit, logistic
import CairoMakie
import AlgebraOfGraphics as AoG
import PairPlots
using CairoMakie: Figure, Axis, hist!, density!, vlines!, hlines!, vspan!,
                  lines!, scatter!, band!, linesegments!

export REPORT_SCENARIOS, REPORT_SCENARIOS_CI,
       ITURI_POPULATION, ITURI_DAILY_TRAVEL,
       ITURI_DAILY_TRAVEL_SD, RENEWAL_START_LEAD,
       load_observations, freeze_observations, m_prior_centre,
       summary_table, posterior_summary,
       fit_diagnostics, diagnostics_table,
       streams_table, comparison_table,
       nuts_sample, fit_parallel, default_adtype, enzyme_adtype,
       progress_callback, tensorboard_callback,
       plot_cumulative_cases, plot_cumulative_trajectories,
       plot_stream_trajectories,
       plot_density_overlay, plot_prior_predictive,
       plot_posterior_predictive, plot_posterior_predictive_grid,
       plot_pair, plot_start_date_pair, plot_estimate_comparison,
       plot_estimate_evolution,
       plot_cfr_prior, plot_vintage_conditional_ppc, plot_rt,
       predict_no_onward_deaths, plot_no_onward_deaths,
       forecast_reported, forecast_table, plot_forecast,
       plot_forecast_latent,
       forecast_vs_truth, forecast_vs_truth_trajectory,
       plot_forecast_vs_truth, plot_forecast_vs_truth_latent,
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
       generation_interval_model, rt_walk_model,
       seed_model, exponential_growth_model, infection_model,
       onset_incidence_model,
       genetic_seeding_model,
       cfr_model, traveller_volume_model, test_positivity_model,
       isolation_admission_model, recovery_probability_model,
       death_background_model, death_ascertainment_model, background_cfr_model,
       background_re_model, background_pooling_model,
       background_walk_model,
       expand_vintage_rate,
       test_sensitivity_model, test_specificity_model, lab_delay_model,
       confirmed_positivity_model, severity_enrichment_model,
       death_testing_fraction_model,
       surveillance_dispersion_model,
       independent_ascertainment_model, pooled_ascertainment_model,
# observation models
       deaths_model, reported_cases_model, confirmed_cases_model,
       confirmed_positivity_windows, confirmed_deaths_model,
       treatment_admission_model, recovered_model,
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
include("counterfactual.jl")
include("forecast.jl")
include("confirmed_cfr.jl")
include("plots.jl")
include("models/priors.jl")
include("models/observations.jl")
include("models/joint.jl")

end # module
