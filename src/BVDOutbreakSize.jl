module BVDOutbreakSize

using Statistics: quantile, mean, std
using TOML: TOML
using DataFrames: DataFrame, rename
using Chain: @chain
using Random: MersenneTwister
using Dates: Date, date2epochdays, epochdays2date
using ADTypes: AutoMooncake
using Mooncake: Mooncake
using ChainRulesCore: ChainRulesCore, NoTangent
using SpecialFunctions: digamma, loggamma
import SpecialFunctions
using Turing: @model, MCMCThreads, NUTS, sample, to_submodel, filldist
using Turing.DynamicPPL: InitFromPrior
import FlexiChains
using DocStringExtensions: @template, DOCSTRING, EXPORTS, IMPORTS, TYPEDEF,
                           TYPEDFIELDS, TYPEDSIGNATURES
using Distributions: Distribution, Gamma, cdf, ccdf, mgf, pdf, Poisson,
                     NegativeBinomial, Normal, LogNormal, Beta,
                     truncated, censored
using StatsFuns: logit, logistic
using Integrals: IntegralProblem, GaussLegendre, solve
import FastGaussQuadrature
import CairoMakie
import AlgebraOfGraphics as AoG
import PairPlots
using CairoMakie: Figure, Axis, hist!, density!, vlines!, vspan!,
                  lines!, scatter!, band!

export REPORT_SCENARIOS,
       ITURI_POPULATION, ITURI_DAILY_TRAVEL,
       ITURI_DAILY_TRAVEL_SD,
       load_observations, m_prior_centre,
       summary_table, posterior_summary,
       fit_diagnostics, diagnostics_table,
       streams_table, comparison_table,
       nuts_sample, default_adtype, enzyme_adtype,
       DEATH_INTEGRAL_ALG, CUMULATIVE_INTEGRAL_ALG,
       integrate, delay_convolution,
       integrate_cumulative, integrate_exports_deaths,
       expected_exports, expected_exports_deaths,
       ExportDeathDelay, EXPORT_DELAY_GRID_POINTS,
       DailyBVDTrajectory, daily_increment_kernel,
       plot_cumulative_cases, plot_density_overlay, plot_prior_predictive,
       plot_posterior_predictive, plot_posterior_predictive_grid,
       plot_pair, plot_start_date_pair, plot_estimate_comparison,
       plot_cfr_prior, plot_vintage_ppc,
       predict_no_onward_deaths, plot_no_onward_deaths,
       forecast_reported, forecast_table, plot_forecast,
       forecast_vs_truth, forecast_vs_truth_trajectory,
       plot_forecast_vs_truth,
# prior submodels
       exponential_growth_model, genetic_seeding_model, delay_model,
       report_delay_model, lab_delay_model, test_sensitivity_model,
       test_positivity_model,
       cfr_model, detection_window_model, traveller_volume_model,
       surveillance_dispersion_model, pooled_ascertainment_model,
       daily_ascertainment_model, deaths_ascertainment_model,
# observation models
       exports_model, deaths_model,
       reported_cases_model, confirmed_cases_model,
       exports_deaths_model,
       exports_detection_timing_model,
# joint composers
       exports_only_model, deaths_only_model, cases_only_model,
       confirmed_only_model,
       exports_deaths_only_model, bvd_joint,
       imperial_only_model

include("docstrings.jl")
include("constants.jl")
include("data.jl")
include("sampling.jl")
include("gamma_cdf.jl")
include("integrate.jl")
include("expectations.jl")
include("summaries.jl")
include("counterfactual.jl")
include("forecast.jl")
include("plots.jl")
include("models/priors.jl")
include("models/observations.jl")
include("models/joint.jl")

end # module
