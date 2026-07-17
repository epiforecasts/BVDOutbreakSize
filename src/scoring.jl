# Forecast-scoring primitives, following scoringutils conventions
# (https://epiforecasts.io/scoringutils/): score a single observation
# against a predictive sample using proper scoring rules, backed by
# ScoringRules.jl. Reuses the `bias_sample` and `_covered` helpers
# defined in summaries.jl.

"""
Continuous ranked probability score of a predictive `samples` ensemble
at a single observation `obs`, via `ScoringRules.crps`. Lower is
better; zero for a perfect point forecast. For a point-mass (constant)
ensemble the CRPS reduces to the absolute error `abs(obs - samples[1])`.
"""
function crps_sample(obs::Real, samples::AbstractVector{<:Real})
    return Float64(crps(samples, obs))
end

"""
CRPS on the log scale: both `obs` and every element of `samples` are
`log1p`-transformed before scoring. This is scoring *on the log
scale*, i.e. `crps_sample(log1p(obs), log1p.(samples))` — it is not
the logarithm of [`crps_sample`](@ref). Log-scale scoring downweights
the influence of large counts, appropriate when forecast accuracy
matters proportionally rather than in absolute terms (e.g. case
counts spanning orders of magnitude).
"""
function log_crps_sample(obs::Real, samples::AbstractVector{<:Real})
    return crps_sample(log1p(obs), log1p.(samples))
end

"""
Score a predictive `samples` ensemble against a single observation
`obs`, returning `(; crps, log_crps, coverage_50, coverage_90, bias, n)`:

  - `crps`: [`crps_sample`](@ref), the ensemble CRPS.
  - `log_crps`: [`log_crps_sample`](@ref), CRPS scored on the log scale.
  - `coverage_50`, `coverage_90`: whether `obs` falls inside the
    central 50% / 90% predictive interval (see `_covered`).
  - `bias`: [`bias_sample`](@ref), signed forecast bias in `[-1, 1]`.
  - `n`: number of predictive draws.

With no draws (`isempty(samples)`), the scores are `NaN`, the
coverage flags are `false`, and `n = 0`.
"""
function score_draws(obs::Real, samples::AbstractVector{<:Real})
    n = length(samples)
    if n == 0
        return (;
            crps = NaN, log_crps = NaN,
            coverage_50 = false, coverage_90 = false,
            bias = NaN, n = 0)
    end
    return (;
        crps = crps_sample(obs, samples),
        log_crps = log_crps_sample(obs, samples),
        coverage_50 = _covered(obs, samples, 0.5),
        coverage_90 = _covered(obs, samples, 0.9),
        bias = bias_sample(obs, samples),
        n = n)
end
