# Forecast-scoring primitives, following scoringutils conventions
# (https://epiforecasts.io/scoringutils/): score a single observation
# against a predictive sample using proper scoring rules. Reuses the
# `bias_sample` and `_covered` helpers defined in summaries.jl.
#
# The overall, by-horizon and by-release summaries near the end of the
# file aggregate a scored table into report-ready rows, one per stream and
# fit (pooled, by horizon, or by release), each relative skill a ratio of
# aggregate mean scores over a matched set of forecasts rather than a mean
# of per-forecast ratios.

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
at a single observation `obs`. Lower is better, and it is zero for a
perfect point forecast. For a point-mass (constant) ensemble the CRPS
reduces to the absolute error `abs(obs - samples[1])`.
"""
function crps_sample(obs::Real, samples::AbstractVector{<:Real})
    return _crps_ensemble(samples, obs)
end

"""
CRPS on the log scale: both `obs` and every element of `samples` are
`log1p`-transformed before scoring, i.e.
`crps_sample(log1p(obs), log1p.(samples))`. It is not the logarithm of
[`crps_sample`](@ref). Log-scale scoring downweights
the influence of large counts, appropriate when forecast accuracy
matters proportionally rather than in absolute terms (e.g. case
counts spanning orders of magnitude).
"""
function log_crps_sample(obs::Real, samples::AbstractVector{<:Real})
    return crps_sample(log1p(obs), log1p.(samples))
end

"""
Sample-based CRPS decomposition of a predictive `samples` ensemble at a
single observation `obs`, into `(; dispersion, overprediction,
underprediction)`, following the convention scoringutils uses for the
weighted interval score, applied here to an ensemble's own order
statistics rather than a fixed set of quantiles.

`dispersion` is the ensemble's spread on its own terms: pair the `i`-th
smallest and `i`-th largest draw as a central interval and sum each
pair's width, weighted so the sum matches the CRPS the ensemble would
score against its own median, the best location it could be scored
against. `overprediction` and `underprediction` are the extra cost of
scoring against `obs` rather than that best case, split by which side of
each pair `obs` falls on: `overprediction` when a pair sits entirely
above `obs` (the ensemble reads too high), `underprediction` when a pair
sits entirely below it (the ensemble reads too low). All three are
non-negative and sum exactly to [`crps_sample`](@ref) at the same `obs`
and `samples`, since the construction is an exact re-partition of the
same order-statistic sum `crps_sample` evaluates, not an approximation.

With no draws (`isempty(samples)`), every component is `NaN`.
"""
function crps_decomposition(obs::Real, samples::AbstractVector{<:Real})
    n = length(samples)
    if n == 0
        return (; dispersion = NaN, overprediction = NaN,
            underprediction = NaN)
    end
    x = sort!(Float64.(collect(samples)))
    y = Float64(obs)
    m = n ÷ 2
    dispersion = zero(Float64)
    overprediction = zero(Float64)
    underprediction = zero(Float64)
    for i in 1:m
        lo = x[i]
        hi = x[n - i + 1]
        w = (2i - 1) / n^2
        dispersion += w * (hi - lo)
        overprediction += (2 / n) * max(lo - y, zero(Float64))
        underprediction += (2 / n) * max(y - hi, zero(Float64))
    end
    if isodd(n)
        med = x[m + 1]
        overprediction += max(med - y, zero(Float64)) / n
        underprediction += max(y - med, zero(Float64)) / n
    end
    return (; dispersion, overprediction, underprediction)
end

"""
Score a predictive `samples` ensemble against a single observation
`obs`, returning `(; crps, log_crps, dispersion, overprediction,
underprediction, coverage_50, coverage_90, bias, n)`:

  - `crps`: [`crps_sample`](@ref), the ensemble CRPS.
  - `log_crps`: [`log_crps_sample`](@ref), CRPS scored on the log scale.
  - `dispersion`, `overprediction`, `underprediction`:
    [`crps_decomposition`](@ref), summing to `crps`.
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
            dispersion = NaN, overprediction = NaN, underprediction = NaN,
            coverage_50 = false, coverage_90 = false,
            bias = NaN, n = 0)
    end
    decomp = crps_decomposition(obs, samples)
    return (;
        crps = crps_sample(obs, samples),
        log_crps = log_crps_sample(obs, samples),
        dispersion = decomp.dispersion,
        overprediction = decomp.overprediction,
        underprediction = decomp.underprediction,
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
## forecasts live under `forecasts-backfill`. Neither is a results release.
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
`main`) are eligible. Any other tag is ignored.

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

## Match `dfa` and `dfb` on `(release, horizon)`, the identity of one scored
## forecast target: the same release, at the same horizon within it, is
## always the same target window. Returns the two frames cut down to the
## keys they share, aligned row for row, so a mean taken over either side
## is a mean over the identical set of forecasts. Assumes at most one row
## per key in each frame, true of any single-fit subset of a scored table.
function _matched_scores(dfa::DataFrame, dfb::DataFrame)
    ia = Dict{Tuple{String, Int}, Int}()
    for i in 1:size(dfa, 1)
        ia[(dfa.release[i], dfa.horizon[i])] = i
    end
    idx_a = Int[]
    idx_b = Int[]
    for j in 1:size(dfb, 1)
        k = (dfb.release[j], dfb.horizon[j])
        haskey(ia, k) || continue
        push!(idx_a, ia[k])
        push!(idx_b, j)
    end
    return dfa[idx_a, :], dfb[idx_b, :]
end

## A ratio of two aggregate mean scores, `missing` rather than `Inf` or
## `NaN` when the denominator is not a finite non-zero number or the ratio
## itself is not finite (the numerical guard every skill aggregate needs:
## no aggregate may ever publish a non-finite value).
function _safe_ratio(num::Real, den::Real; digits::Int = 2)
    (isfinite(den) && den != 0) || return missing
    r = num / den
    return isfinite(r) ? round(r; digits = digits) : missing
end

## The id of `stream`'s individual single-stream fit: the one fit in
## `scores`, besides `joint_fit` and `baseline_fit`, scored anywhere for
## that stream. `missing` when the stream has none (e.g. "recovered"),
## since a stream is fit by at most one individual model.
function _individual_fit_id(scores::DataFrame, stream, joint_fit,
        baseline_fit)
    not_base_or_joint = (scores.fit .!= joint_fit) .&
                        (scores.fit .!= baseline_fit)
    ids = unique(scores.fit[(scores.stream .== stream) .& not_base_or_joint])
    return length(ids) == 1 ? ids[1] : missing
end

## Aggregate score stats for one `fit` of `stream`, restricted to `mask` (a
## `BitVector` over `scores`' rows already narrowed to the stream and,
## depending on which table is being built, a single horizon or made
## date). `nothing` when `fit` has no forecast in `mask` matched against
## its stream's baseline, so the caller can drop the row entirely.
##
## Every relative skill is the ratio of the two aggregate mean scores over
## the matched set of forecasts both the fit and its comparator scored
## (R25): a mean taken over one set of forecasts is never divided by a
## mean taken over another, and `_safe_ratio` keeps the result finite. The
## individual-fit columns use their own matched set against the stream's
## individual fit and are only ever populated on the joint fit's row of a
## stream that has one.
function _stream_fit_stats(scores::DataFrame, stream, fit, mask;
        joint_fit, baseline_fit)
    fit_grp = scores[mask .& (scores.fit .== fit), :]
    isempty(fit_grp) && return nothing
    baseline_grp = scores[mask .& (scores.fit .== baseline_fit), :]
    mfit, mbase = _matched_scores(fit_grp, baseline_grp)
    n = size(mfit, 1)
    n == 0 && return nothing

    crps = mean(mfit.crps)
    log_crps = mean(mfit.log_crps)
    ## The baseline's and the individual fit's own mean CRPS are only ever
    ## intermediate: they feed the ratios below and are not published as
    ## columns themselves (R25 keeps the table to the ratios a reader
    ## compares against one, not the absolute scores behind them).
    crps_baseline = mean(mbase.crps)
    log_crps_baseline = mean(mbase.log_crps)

    rel_to_individual = missing
    log_rel_to_individual = missing
    if fit == joint_fit
        indiv = _individual_fit_id(scores, stream, joint_fit, baseline_fit)
        if !ismissing(indiv)
            indiv_grp = scores[mask .& (scores.fit .== indiv), :]
            mfit2, mind = _matched_scores(fit_grp, indiv_grp)
            if size(mfit2, 1) > 0
                fit_c, ind_c = mean(mfit2.crps), mean(mind.crps)
                fit_lc, ind_lc = mean(mfit2.log_crps), mean(mind.log_crps)
                rel_to_individual = _safe_ratio(fit_c, ind_c)
                log_rel_to_individual = _safe_ratio(fit_lc, ind_lc)
            end
        end
    end

    return (; n,
        crps = round(crps; digits = 2),
        rel_to_baseline = _safe_ratio(crps, crps_baseline),
        log_crps = round(log_crps; digits = 3),
        log_rel_to_baseline = _safe_ratio(log_crps, log_crps_baseline),
        rel_to_individual, log_rel_to_individual,
        dispersion = round(mean(mfit.dispersion); digits = 2),
        overprediction = round(mean(mfit.overprediction); digits = 2),
        underprediction = round(mean(mfit.underprediction); digits = 2),
        coverage_50 = round(mean(mfit.coverage_50); digits = 2),
        coverage_90 = round(mean(mfit.coverage_90); digits = 2),
        bias = round(mean(mfit.bias); digits = 2))
end

## Column schema shared by the three score summaries below, `key_cols`
## first, the metric columns in the order the report shows them. An empty
## and a populated frame from the same builder always carry identical
## column names, order and types, so the report renders before any
## release carries a forecast.
function _score_summary_schema(key_cols::NamedTuple)
    metrics = (;
        n = Int[], crps = Float64[],
        rel_to_baseline = Union{Missing, Float64}[],
        log_crps = Float64[],
        log_rel_to_baseline = Union{Missing, Float64}[],
        rel_to_individual = Union{Missing, Float64}[],
        log_rel_to_individual = Union{Missing, Float64}[],
        dispersion = Float64[], overprediction = Float64[],
        underprediction = Float64[],
        coverage_50 = Float64[], coverage_90 = Float64[], bias = Float64[])
    return DataFrame(; key_cols..., metrics...)
end

"""
One row per `(stream, fit)`, pooled over every horizon and release in
`scores` (a `data/forecast_scores.csv`-shaped table, as written by
`scripts/score_releases.jl`), baseline rows excluded. This is the headline
scoring table: `n`, the mean CRPS and log-scale CRPS, the CRPS
decomposition (dispersion, overprediction, underprediction), the 50% and
90% coverage and bias, the relative skill against the stream's persistence
baseline on both scales, and, on the joint row of a stream that has an
individual single-stream fit, the relative skill against that fit on both
scales.

Every relative skill is the ratio of the two aggregate mean scores over the
matched set of forecasts both sides scored, not a mean of per-forecast
ratios, so one forecast whose comparator scored zero cannot make the ratio
infinite. A ratio is `missing`, never `Inf` or `NaN`, when the comparator's
mean is zero or the ratio itself is not finite.

Returns a typed zero-row frame when `scores` is empty, so the report
renders before any release carries a forecast.
"""
function forecast_score_overview(scores::DataFrame;
        joint_fit = JOINT_FIT, baseline_fit = BASELINE_FIT)
    empty = _score_summary_schema((; stream = String[], fit = String[]))
    isempty(scores) && return empty
    rows = NamedTuple[]
    for s in sort(unique(scores.stream))
        mask = scores.stream .== s
        fits = sort(unique(scores.fit[mask .& (scores.fit .!= baseline_fit)]))
        for f in fits
            st = _stream_fit_stats(scores, s, f, mask; joint_fit, baseline_fit)
            isnothing(st) && continue
            push!(rows, (; stream = s, fit = f, st...))
        end
    end
    isempty(rows) && return empty
    return DataFrame(rows)
end

"""
The same columns as [`forecast_score_overview`](@ref), one row per
`(stream, horizon, fit)` rather than pooled across every horizon, so a fit
that beats the baseline on average but not at every cut-off is visible.
"""
function forecast_score_by_horizon(scores::DataFrame;
        joint_fit = JOINT_FIT, baseline_fit = BASELINE_FIT)
    empty = _score_summary_schema((; stream = String[], horizon = Int[],
        fit = String[]))
    isempty(scores) && return empty
    rows = NamedTuple[]
    for s in sort(unique(scores.stream))
        smask = scores.stream .== s
        for h in sort(unique(scores.horizon[smask]))
            mask = smask .& (scores.horizon .== h)
            fits = sort(unique(
                scores.fit[mask .& (scores.fit .!= baseline_fit)]))
            for f in fits
                st = _stream_fit_stats(scores, s, f, mask; joint_fit,
                    baseline_fit)
                isnothing(st) && continue
                push!(rows, (; stream = s, horizon = h, fit = f, st...))
            end
        end
    end
    isempty(rows) && return empty
    return DataFrame(rows)
end

"""
The same columns as [`forecast_score_overview`](@ref), one row per
`(stream, fit, made_date)`, averaged across horizons rather than pooled
over releases too, baseline rows excluded, so the by-release detail table
stays compact.
"""
function forecast_score_by_release(scores::DataFrame;
        joint_fit = JOINT_FIT, baseline_fit = BASELINE_FIT)
    empty = _score_summary_schema((; made_date = Date[], stream = String[],
        fit = String[]))
    isempty(scores) && return empty
    rows = NamedTuple[]
    for s in sort(unique(scores.stream))
        smask = scores.stream .== s
        for d in sort(unique(scores.made_date[smask]))
            mask = smask .& (scores.made_date .== d)
            fits = sort(unique(
                scores.fit[mask .& (scores.fit .!= baseline_fit)]))
            for f in fits
                st = _stream_fit_stats(scores, s, f, mask; joint_fit,
                    baseline_fit)
                isnothing(st) && continue
                push!(rows, (; made_date = d, stream = s, fit = f, st...))
            end
        end
    end
    isempty(rows) && return empty
    return DataFrame(rows)
end

const _INDIVIDUAL_FIT_COLUMNS = ("rel_to_individual", "log_rel_to_individual")

"""
`table` (a table from [`forecast_score_overview`](@ref),
[`forecast_score_by_horizon`](@ref) or [`forecast_score_by_release`](@ref))
with the individual-fit comparison columns dropped rather than kept as an
all-`missing` column. Use this to display a table built from an evaluation
that never carries an individual single-stream fit, such as a frozen-fit
evaluation, which scores only the joint model at past cut-offs.
"""
function drop_individual_fit_columns(table::DataFrame)
    keep = [n for n in names(table) if !(n in _INDIVIDUAL_FIT_COLUMNS)]
    return table[:, keep]
end

"""
`table` (from [`forecast_score_overview`](@ref),
[`forecast_score_by_horizon`](@ref) or [`forecast_score_by_release`](@ref))
with its `fit` column dropped, for a table whose `fit` column is
single-valued by construction rather than merely by the data currently on
hand: the frozen-fit evaluation, which scores the current joint model
alone at past cut-offs, is the only such table in this report.

Every other column (`crps`, `dispersion`, `overprediction`,
`underprediction`, `coverage_50`, `coverage_90`, …) is left untouched. This
drops one column and nothing else. `table` is returned with `fit` still
present when it is empty, since an empty table carries no evidence either
way. State which model a de-columned table refers to in the surrounding
prose, so the table stays self-describing once the column is gone.

Errors when `table`'s `fit` column carries more than one distinct value,
rather than silently dropping it: a column that vanishes only when it
happens to be constant and reappears once a second model is scored is
worse than a constant column that never varies, so this is only safe to
call on a table whose single-fit shape is structural (as the frozen
evaluation's is), never on a table that could legitimately grow a second
model later.
"""
function drop_degenerate_fit_column(table::DataFrame)
    isempty(table.fit) && return table
    length(unique(table.fit)) == 1 ||
        error("drop_degenerate_fit_column: table's fit column carries " *
              "more than one value; only call this on a table whose " *
              "single-fit shape is structural, not incidental")
    keep = [n for n in names(table) if n != "fit"]
    return table[:, keep]
end
