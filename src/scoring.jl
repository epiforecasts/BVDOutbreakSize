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

## The two kinds of results release `.github/workflows/docs.yml` publishes:
## `results-vX.Y.Z` from a version-tag push, and `results-<run number>`
## from a push to `main`.
const _VERSION_RELEASE = r"^results-v([0-9][0-9A-Za-z.+-]*)$"
const _MAIN_RELEASE = r"^results-([0-9]+)$"

## Rank one results release within its day: tagged releases above main
## builds, then the later timestamp, then the higher version or run number
## so releases sharing a timestamp still order deterministically. Returns
## `nothing` for a tag that is not a results release. The version and run
## number are only ever compared against their own kind, since the leading
## flag separates the two.
function _release_key(tag::AbstractString, created)
    m = match(_VERSION_RELEASE, tag)
    isnothing(m) || return (1, created, VersionNumber(m[1]))
    m = match(_MAIN_RELEASE, tag)
    isnothing(m) && return nothing
    return (0, created, parse(Int, m[1]))
end

"""
Select at most one results release per calendar day from `entries`, a
collection of `(tag, created)` pairs as reported by `gh release list`
(`created` a UTC `DateTime`). Both tagged releases (`results-vX.Y.Z`) and
main-build releases (`results-<run number>`, published by every push to
`main`) are eligible; any other tag is ignored. Days are the UTC calendar
days of the release timestamps.

Within a day a tagged release is preferred: a version tag and the `main`
push of the same commit publish two releases carrying identical forecasts,
so scoring both would count one forecast twice. Failing that the newest
main build of the day is kept, so the day is represented by the forecast
made from the most data. Releases sharing a timestamp are separated by
version or run number, keeping the selection deterministic.

Returns the selected tags, newest first.
"""
function select_daily_releases(entries)
    rels = [(String(t), c) for (t, c) in entries
            if !isnothing(_release_key(t, c))]
    isempty(rels) && return String[]

    days = Dict{Date, typeof(rels)}()
    for r in rels
        push!(get!(() -> empty(rels), days, Date(r[2])), r)
    end

    picks = [argmax(r -> _release_key(r[1], r[2]), days[d])
             for d in sort(collect(keys(days)))]
    sort!(picks; by = r -> r[2], rev = true)
    return first.(picks)
end

"""
Summarise a `data/forecast_scores.csv` table (as written by
`scripts/score_releases.jl`) into one row per `(stream, horizon)`: the number
of scored forecasts, the mean CRPS and log-scale CRPS of our forecasts, the
relative skill against the persistence baseline (`rel_skill`, the mean CRPS
ratio, below one when we beat the baseline), and the empirical 50% and 90%
coverage and mean bias. Returns an empty frame with the same columns when no
forecasts have been scored yet, so the report renders before any release
carries a stored forecast.
"""
function forecast_score_summary(scores::DataFrame)
    rows = NamedTuple[]
    if !isempty(scores)
        for s in sort(unique(scores.stream))
            hs = sort(unique(scores.horizon[scores.stream .== s]))
            for h in hs
                grp = scores[(scores.stream .== s) .& (scores.horizon .== h), :]
                ours = grp[grp.model .== "ours", :]
                base = grp[grp.model .== "baseline", :]
                isempty(ours) && continue
                mc = mean(ours.crps)
                mb = isempty(base) ? NaN : mean(base.crps)
                push!(rows,
                    (; stream = s, horizon = h, n = size(ours, 1),
                        crps = round(mc; digits = 2),
                        crps_baseline = round(mb; digits = 2),
                        rel_skill = round(mc / mb; digits = 2),
                        log_crps = round(mean(ours.log_crps); digits = 3),
                        coverage_50 = round(mean(ours.coverage_50); digits = 2),
                        coverage_90 = round(mean(ours.coverage_90); digits = 2),
                        bias = round(mean(ours.bias); digits = 2)))
            end
        end
    end
    isempty(rows) && return DataFrame(
        stream = String[], horizon = Int[], n = Int[], crps = Float64[],
        crps_baseline = Float64[], rel_skill = Float64[], log_crps = Float64[],
        coverage_50 = Float64[], coverage_90 = Float64[], bias = Float64[])
    return DataFrame(rows)
end
