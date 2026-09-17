# Forecast-scoring primitives, following scoringutils conventions
# (https://epiforecasts.io/scoringutils/). Score a single observation
# against a predictive sample using proper scoring rules. Reuses the
# `bias_sample` and `_covered` helpers defined in summaries.jl.
#
# The summaries near the end of the file aggregate a scored table into
# report-ready rows, one per stream and fit. Each relative skill is a ratio
# of aggregate mean scores over a matched set of forecasts rather than a
# mean of per-forecast ratios.

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
CRPS on the log scale. Both `obs` and every element of `samples` are
`log1p`-transformed before scoring, so this is
`crps_sample(log1p(obs), log1p.(samples))` and not the logarithm of
[`crps_sample`](@ref). Log-scale scoring downweights large counts, which
suits accuracy that matters proportionally rather than in absolute terms.
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

`dispersion` is the ensemble's spread on its own terms. Pair the `i`-th
smallest and `i`-th largest draw as a central interval and sum each pair's
width, weighted so the sum matches the CRPS the ensemble would score
against its own median. `overprediction` and `underprediction` are the
extra cost of scoring against `obs` rather than that best case, split by
which side of each pair `obs` falls on. All three are non-negative and sum
exactly to [`crps_sample`](@ref) at the same `obs` and `samples`, being an
exact re-partition of the same order-statistic sum rather than an
approximation.

With no draws (`isempty(samples)`), every component is `NaN`.
"""
function crps_decomposition(obs::Real, samples::AbstractVector{<:Real})
    n = length(samples)
    if n == 0
        return (;
            dispersion = NaN, overprediction = NaN,
            underprediction = NaN,
        )
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
            bias = NaN, n = 0,
        )
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
        n = n,
    )
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

## Rank one results release within its day. Tagged releases sort above
## main builds, then the later timestamp, then the higher version or run
## number so releases sharing a timestamp still order deterministically.
## Returns `nothing` for a tag that is not a results release. The leading
## flag keeps version and run number from being compared with each other.
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
same forecast however far apart they were published, so grouping on the
creation timestamp would score that forecast twice.

Within a day a tagged release is preferred, then the newest build of the
day. Releases sharing a timestamp are separated by version or run number,
keeping the selection deterministic.

Returns the selected tags, newest cut-off first.
"""
function select_daily_releases(entries)
    rels = [
        (String(t), c, Date(d))
            for (t, c, d) in entries
            if !isnothing(_release_key(t, c))
    ]
    isempty(rels) && return String[]

    days = Dict{Date, typeof(rels)}()
    for r in rels
        push!(get!(() -> empty(rels), days, r[3]), r)
    end

    picks = [
        argmax(r -> _release_key(r[1], r[2]), days[d])
            for d in sort(collect(keys(days)))
    ]
    sort!(picks; by = r -> r[3], rev = true)
    return first.(picks)
end

## Match `dfa` and `dfb` on `(release, horizon)`, which identifies one
## scored forecast target. Returns the two frames cut down to the keys they
## share, aligned row for row, so a mean taken over either side is a mean
## over the identical set of forecasts. Assumes at most one row per key in
## each frame, true of any single-fit subset of a scored table.
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
## itself is not finite. No aggregate may publish a non-finite value.
function _safe_ratio(num::Real, den::Real; digits::Int = 2)
    (isfinite(den) && den != 0) || return missing
    r = num / den
    return isfinite(r) ? round(r; digits = digits) : missing
end

## The id of `stream`'s individual single-stream fit: the one fit in
## `scores`, besides `joint_fit` and `baseline_fit`, scored anywhere for
## that stream. `missing` when the stream has none (e.g. "recovered"),
## since a stream is fit by at most one individual model.
function _individual_fit_id(
        scores::DataFrame, stream, joint_fit,
        baseline_fit
    )
    not_base_or_joint = (scores.fit .!= joint_fit) .&
        (scores.fit .!= baseline_fit)
    ids = unique(scores.fit[(scores.stream .== stream) .& not_base_or_joint])
    return length(ids) == 1 ? ids[1] : missing
end

## Aggregate score stats for one `fit` of `stream`, restricted to `mask`, a
## `BitVector` over `scores`' rows already narrowed to the stream and,
## depending on which table is being built, a single horizon or made date.
## `nothing` when `fit` has no forecast in `mask` matched against its
## stream's baseline, so the caller can drop the row entirely.
##
## Every relative skill is the ratio of the two aggregate mean scores over
## the matched set of forecasts both the fit and its comparator scored, so
## a mean over one set of forecasts is never divided by a mean over
## another. The individual-fit columns use their own matched set and are
## only populated on the joint fit's row of a stream that has one.
function _stream_fit_stats(
        scores::DataFrame, stream, fit, mask;
        joint_fit, baseline_fit
    )
    fit_grp = scores[mask .& (scores.fit .== fit), :]
    isempty(fit_grp) && return nothing
    baseline_grp = scores[mask .& (scores.fit .== baseline_fit), :]
    mfit, mbase = _matched_scores(fit_grp, baseline_grp)
    n = size(mfit, 1)
    n == 0 && return nothing

    crps = mean(mfit.crps)
    log_crps = mean(mfit.log_crps)
    ## The baseline's and the individual fit's own mean CRPS are
    ## intermediate. They feed the ratios below and are not published as
    ## columns themselves.
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

    return (;
        n,
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
        bias = round(mean(mfit.bias); digits = 2),
    )
end

## Column schema shared by the three score summaries below, `key_cols`
## first, then the metric columns in the order the report shows them. An
## empty and a populated frame from the same builder carry identical column
## names, order and types, so the report renders before any release carries
## a forecast.
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
        coverage_50 = Float64[], coverage_90 = Float64[], bias = Float64[],
    )
    return DataFrame(; key_cols..., metrics...)
end

"""
One row per `(stream, fit)`, pooled over every horizon and release in
`scores`, a `data/forecast_scores.csv`-shaped table, baseline rows
excluded. The headline scoring table. Carries `n`, the mean CRPS and
log-scale CRPS, the CRPS decomposition, the 50% and 90% coverage and bias,
the relative skill against the stream's persistence baseline on both
scales, and, on the joint row of a stream with an individual single-stream
fit, the relative skill against that fit.

Every relative skill is the ratio of the two aggregate mean scores over the
matched set of forecasts both sides scored, not a mean of per-forecast
ratios, so one forecast whose comparator scored zero cannot make the ratio
infinite. A ratio is `missing`, never `Inf` or `NaN`, when the comparator's
mean is zero or the ratio itself is not finite.

Returns a typed zero-row frame when `scores` is empty, so the report
renders before any release carries a forecast.
"""
function forecast_score_overview(
        scores::DataFrame;
        joint_fit = JOINT_FIT, baseline_fit = BASELINE_FIT
    )
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
function forecast_score_by_horizon(
        scores::DataFrame;
        joint_fit = JOINT_FIT, baseline_fit = BASELINE_FIT
    )
    empty = _score_summary_schema(
        (;
            stream = String[], horizon = Int[],
            fit = String[],
        )
    )
    isempty(scores) && return empty
    rows = NamedTuple[]
    for s in sort(unique(scores.stream))
        smask = scores.stream .== s
        for h in sort(unique(scores.horizon[smask]))
            mask = smask .& (scores.horizon .== h)
            fits = sort(
                unique(
                    scores.fit[mask .& (scores.fit .!= baseline_fit)]
                )
            )
            for f in fits
                st = _stream_fit_stats(
                    scores, s, f, mask; joint_fit,
                    baseline_fit
                )
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
over releases too, baseline rows excluded.
"""
function forecast_score_by_release(
        scores::DataFrame;
        joint_fit = JOINT_FIT, baseline_fit = BASELINE_FIT
    )
    empty = _score_summary_schema(
        (;
            made_date = Date[], stream = String[],
            fit = String[],
        )
    )
    isempty(scores) && return empty
    rows = NamedTuple[]
    for s in sort(unique(scores.stream))
        smask = scores.stream .== s
        for d in sort(unique(scores.made_date[smask]))
            mask = smask .& (scores.made_date .== d)
            fits = sort(
                unique(
                    scores.fit[mask .& (scores.fit .!= baseline_fit)]
                )
            )
            for f in fits
                st = _stream_fit_stats(
                    scores, s, f, mask; joint_fit,
                    baseline_fit
                )
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
with its `fit` column dropped. Only for a table whose `fit` column is
single-valued by construction rather than by the data currently on hand.
The frozen-fit evaluation, which scores the current joint model alone at
past cut-offs, is the only such table in this report.

Every other column is left untouched. `table` is returned with `fit` still
present when it is empty, since an empty table carries no evidence either
way. State which model a de-columned table refers to in the surrounding
prose, so the table stays self-describing once the column is gone.

Errors when `table`'s `fit` column carries more than one distinct value. A
column that vanishes only when it happens to be constant, and reappears
once a second model is scored, is worse than one that never varies.
"""
function drop_degenerate_fit_column(table::DataFrame)
    isempty(table.fit) && return table
    length(unique(table.fit)) == 1 ||
        error(
        "drop_degenerate_fit_column: table's fit column carries " *
            "more than one value; only call this on a table whose " *
            "single-fit shape is structural, not incidental"
    )
    keep = [n for n in names(table) if n != "fit"]
    return table[:, keep]
end

## Role of a forecast within a stream: the persistence baseline, the joint
## model or the stream's own single-stream fit. Deriving the role from the
## fit id keeps a figure's colour stable across streams whichever
## single-stream fit produced the individual row. A frozen row takes the
## joint role, being the joint model re-fit at a past cut-off. Every other
## id is a per-stream fit spec, an open-ended set, so it falls through to
## the individual role.
function _fit_role(fit)
    fit == BASELINE_FIT && return "baseline"
    (fit == JOINT_FIT || fit == FROZEN_FIT) && return "joint"
    return "individual"
end

const _FIT_ROLES = ("baseline", "individual", "joint")

"""
`table` (a table from [`forecast_score_overview`](@ref),
[`forecast_score_by_horizon`](@ref) or [`forecast_score_by_release`](@ref))
restricted to the rows of one `role`.

`"joint"` selects the joint model, including a frozen row, which is the
same model re-fit at a past cut-off. `"individual"` selects each stream's
own single-stream fit, whatever its id. `"baseline"` selects the
persistence baseline, which the three summaries above exclude, so it is
only ever non-empty on a table that kept it.

Use this to render one table per role. A figure comparing the roles
against each other reads the unfiltered table instead. Column names, order
and types are unchanged, as is the order of the rows that survive.

Errors on an unknown `role` rather than returning a zero-row table, since
a misspelt role and a role with nothing scored are otherwise
indistinguishable in the rendered report.
"""
function select_fit_role(table::DataFrame, role::AbstractString)
    role in _FIT_ROLES || error(
        "select_fit_role: unknown role " *
            "$(repr(role)); expected one of " *
            join(repr.(_FIT_ROLES), ", ")
    )
    return table[_fit_role.(table.fit) .== role, :]
end

"""
`overlay` (a `data/forecast_overlay.csv`-shaped table) restricted to the
streams that carry a persistence baseline, dropping every row of a stream
that carries none.

This is the rule the score summaries already apply. `_stream_fit_stats`
drops a fit whose matched set against the baseline is empty, so a stream
with no baseline row never reaches [`forecast_score_overview`](@ref) or
its by-horizon and by-release counterparts. Applying it to the overlay
figure too puts the tables and the figure on one definition of what has
been scored. A stream the situation reports have stopped publishing then
leaves both together, rather than holding a row of near-empty panels open
long after its last scored window.

The rule is stated on the baseline rather than on a count of made dates
because the baseline is what makes a window scoreable at all. A stream
with no baseline has nothing to be scored against, whatever its history.

Nothing is dropped from the scored data itself. The rows stay in
`data/forecast_scores.csv` and `data/forecast_overlay.csv`, which record
what was scored; this selects what is drawn.

Returns `overlay` unchanged when it is empty, so the report renders before
any release carries a forecast.
"""
function scored_overlay(overlay::DataFrame; baseline_fit = BASELINE_FIT)
    isempty(overlay) && return overlay
    kept = Set(overlay.stream[overlay.fit .== baseline_fit])
    return overlay[in.(overlay.stream, Ref(kept)), :]
end
