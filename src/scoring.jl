# Forecast-scoring primitives, following scoringutils conventions
# (https://epiforecasts.io/scoringutils/): score a single observation
# against a predictive sample using proper scoring rules. Reuses the
# `bias_sample` and `_covered` helpers defined in summaries.jl.

## Continuous ranked probability score of an ensemble `samples` at a point
## observation `obs`, from the energy form `CRPS = E|X - obs| - ½ E|X - X'|`
## evaluated on the empirical ensemble. The pairwise term uses the sorted
## closed form `Σ_ij |x_i - x_j| = 2 Σ_i (2i - n - 1) x_(i)`, so the cost is
## the sort rather than the O(n²) double loop.
function _crps_ensemble(samples::AbstractVector{<:Real}, obs::Real)
    n = length(samples)
    x = sort!(Float64.(collect(samples)))
    mae = zero(Float64)
    for xi in x
        mae += abs(xi - obs)
    end
    spread = zero(Float64)
    for (i, xi) in enumerate(x)
        spread += (2i - n - 1) * xi
    end
    return mae / n - spread / n^2
end

"""
Continuous ranked probability score of a predictive `samples` ensemble
at a single observation `obs`. Lower is better; zero for a perfect point
forecast. For a point-mass (constant) ensemble the CRPS reduces to the
absolute error `abs(obs - samples[1])`.
"""
function crps_sample(obs::Real, samples::AbstractVector{<:Real})
    return _crps_ensemble(samples, obs)
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

## Whether `tag` names a results release of either kind. The repo also
## publishes a release per code tag (`v1.9.0`), and the reconstructed
## forecasts live under `forecasts-backfill`; neither is a results release.
function is_results_release(tag::AbstractString)
    return !isnothing(match(_VERSION_RELEASE, tag)) ||
           !isnothing(match(_MAIN_RELEASE, tag))
end

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
Select at most one results release per data day from `entries`, a
collection of `(tag, created, cutoff)` triples: the release tag, its
creation timestamp (a UTC `DateTime`, as reported by `gh release list`)
and the data cut-off it was built from (`as_of_date` in the release's
`observations.toml`). Both tagged releases (`results-vX.Y.Z`) and
main-build releases (`results-<run number>`, published by every push to
`main`) are eligible; any other tag is ignored.

Days are cut-off days, not creation days, because a release's forecast is
a function of the data it saw. Two releases sharing a cut-off carry the
same forecast however far apart they were published: the report re-renders
whenever `main` moves, so the same data is republished under a new tag.
Grouping on the creation timestamp would keep both and score that forecast
twice.

Within a day a tagged release is preferred, which also settles the case of
a version tag and the `main` push of one commit publishing at the same
instant. Failing that the newest build of the day is kept. Releases
sharing a timestamp are separated by version or run number, keeping the
selection deterministic.

Returns the selected tags, newest cut-off first.
"""
function select_daily_releases(entries)
    rels = [(String(t), c, Date(d))
            for (t, c, d) in entries
            if !isnothing(_release_key(t, c))]
    isempty(rels) && return String[]

    days = Dict{Date, typeof(rels)}()
    for r in rels
        push!(get!(() -> empty(rels), days, r[3]), r)
    end

    picks = [argmax(r -> _release_key(r[1], r[2]), days[d])
             for d in sort(collect(keys(days)))]
    sort!(picks; by = r -> r[3], rev = true)
    return first.(picks)
end

## Score columns summarised for one group of forecasts, against the mean
## baseline CRPS `mb` and mean baseline log-CRPS `mlb` of the same
## (stream, horizon).
function _score_stats(grp::DataFrame, mb::Real, mlb::Real)
    mc = mean(grp.crps)
    mlc = mean(grp.log_crps)
    rs = mc / mb
    lrs = mlc / mlb
    return (; n = size(grp, 1), crps = round(mc; digits = 2),
        crps_baseline = round(mb; digits = 2),
        rel_skill = round(rs; digits = 2),
        skill = round(1.0 - rs; digits = 3),
        log_crps = round(mlc; digits = 3),
        log_crps_baseline = round(mlb; digits = 3),
        log_rel_skill = round(lrs; digits = 3),
        log_skill = round(1.0 - lrs; digits = 3),
        coverage_50 = round(mean(grp.coverage_50); digits = 2),
        coverage_90 = round(mean(grp.coverage_90); digits = 2),
        bias = round(mean(grp.bias); digits = 2))
end

"""
Summarise a `data/forecast_scores.csv` table (as written by
`scripts/score_releases.jl`): the number of scored forecasts, the mean CRPS
and log-scale CRPS, the relative skill against the persistence baseline
(`rel_skill`, the mean CRPS ratio, below one when the forecast beats the
baseline), and the empirical 50% and 90% coverage and mean bias.

A table carrying a `fit` column, naming the model each row's forecast came
from (a spec id such as `joint` or `confirmed`, or `baseline` for the
persistence baseline), is summarised into one row per
`(stream, horizon, fit)`, one per non-baseline fit, with `rel_skill` taken
against that `(stream, horizon)`'s baseline rows. A stream no individual
model fits, such as "recovered", simply has no row for one.

An older table carrying `model` ∈ {ours, baseline} instead is summarised
into one row per `(stream, horizon)` with no `fit` column.

Returns an empty frame with the matching columns when nothing has been
scored yet, so the report renders before any release carries a forecast.
"""
function forecast_score_summary(scores::DataFrame)
    has_fit = "fit" in names(scores)
    rows = NamedTuple[]
    if !isempty(scores)
        for s in sort(unique(scores.stream))
            hs = sort(unique(scores.horizon[scores.stream .== s]))
            for h in hs
                grp = scores[(scores.stream .== s) .& (scores.horizon .== h), :]
                base = has_fit ? grp[grp.fit .== "baseline", :] :
                       grp[grp.model .== "baseline", :]
                mb = isempty(base) ? NaN : mean(base.crps)
                ## Every non-baseline fit of this stream is scored against
                ## the one baseline the stream's observations imply.
                fits = has_fit ? sort(unique(grp.fit[grp.fit .!= "baseline"])) :
                       ["ours"]
                for f in fits
                    ours = has_fit ? grp[grp.fit .== f, :] :
                           grp[grp.model .== "ours", :]
                    isempty(ours) && continue
                    mlb = isempty(base) ? NaN : mean(base.log_crps)
                    st = _score_stats(ours, mb, mlb)
                    push!(rows,
                        has_fit ? (; stream = s, horizon = h, fit = f, st...) :
                        (; stream = s, horizon = h, st...))
                end
            end
        end
    end
    if isempty(rows)
        empty = (; n = Int[], crps = Float64[], crps_baseline = Float64[],
            rel_skill = Float64[], skill = Float64[],
            log_crps = Float64[], log_crps_baseline = Float64[],
            log_rel_skill = Float64[], log_skill = Float64[],
            coverage_50 = Float64[], coverage_90 = Float64[],
            bias = Float64[])
        return has_fit ?
               DataFrame(; stream = String[], horizon = Int[],
            fit = String[], empty...) :
               DataFrame(; stream = String[], horizon = Int[], empty...)
    end
    return DataFrame(rows)
end
