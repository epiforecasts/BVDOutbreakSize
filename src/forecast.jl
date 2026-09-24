# Posterior-predictive forecasts drawn from the fitted models. Each fitted
# model is rebuilt to run past its cut-off (`with_horizon`), and `predict`
# on it with the fitted chain keeps every draw's parameters and draws the
# future: the walks' future innovations, the renewal, every delay and
# ascertainment, and each stream's future counts through its own
# likelihood. The functions below only read those draws.

"""
    with_horizon(model, horizon) -> DynamicPPL.Model

The fitted `model` rebuilt to run `horizon` days past its cut-off: the same
arguments and keywords with `forecast = ForecastHorizon(horizon)`, or the
fitted model itself (`forecast = nothing`) for a horizon of zero. The
extended model keeps every fitted variable and adds only the future ones,
so `predict` on it with the fitted chain keeps each draw's parameters and
draws the future (see [`ForecastHorizon`](@ref)). Any conditioning or fixing
on `model` is kept.
"""
function with_horizon(model::DynamicPPL.Model, horizon::Integer)
    forecast = horizon > 0 ? ForecastHorizon(horizon) : nothing
    extended = model.f(values(model.args)...; model.defaults..., forecast)
    return DynamicPPL.contextualize(extended, model.context)
end

"""
    forecast_draws(model, chn; horizon = 28, seed = 20260520)

Posterior-predictive draws past the cut-off from a fitted `model` and its
chain `chn`: `predict` on [`with_horizon`](@ref)`(model, horizon)`. Each
draw keeps its fitted parameters. The future walk innovations, the future
latent series and every stream's future counts are drawn from the model.
The result is a chain with one draw per fitted draw, which
[`forecast_reported`](@ref), [`forecast_stream`](@ref),
[`forecast_onsets`](@ref) and [`forecast_provinces`](@ref) read at any
horizon up to `horizon`. `seed` fixes the draws, so the national and the
province forecasts read one set of draws.
"""
function forecast_draws(
        model::DynamicPPL.Model, chn; horizon::Integer = 28,
        seed::Integer = 20260520
    )
    return predict(MersenneTwister(seed), with_horizon(model, horizon), chn)
end

## Per-draw vectors of a forecast key, or `nothing` when the draws do not
## carry it (the fit has no such stream).
function _forecast_vectors(pp, key::AbstractString)
    k = Symbol(key)
    _has_key(pp, k) || return nothing
    return _draw_vectors(pp, k)
end

## Days the draws run past the cut-off, read off the future infections every
## composer tracks.
function _forecast_horizon(pp)
    v = _forecast_vectors(pp, "forecast_infections")
    isnothing(v) && throw(
        ArgumentError(
            "these draws carry no forecast; draw them with `forecast_draws`."
        )
    )
    return length(first(v))
end

function _check_horizon(pp, h::Integer)
    H = _forecast_horizon(pp)
    1 <= h <= H || throw(
        ArgumentError("horizon $h is outside the $H days the draws cover.")
    )
    return H
end

## Sum of each draw's first `h` future values, or `nothing`.
function _forecast_new(pp, key, h)
    v = _forecast_vectors(pp, key)
    isnothing(v) && return nothing
    return [sum(x[1:h]) for x in v]
end

## Each draw's value on future day `h`, or `nothing`.
function _forecast_at(pp, key, h)
    v = _forecast_vectors(pp, key)
    isnothing(v) && return nothing
    return [x[h] for x in v]
end

## Index of the future vintage `h` days past the cut-off among the weekly
## vintages of a forecast `H` days long, or `nothing` when `h` is not one.
_vintage_index(h, H) = findfirst(==(h), future_knot_days(0, H))

"""
    forecast_reported(pp; horizon = 7, obs_cases, obs_deaths, ...) -> DataFrame

The national forecast `horizon` days past the cut-off, one row per draw,
read from the posterior-predictive draws `pp` of the joint
([`forecast_draws`](@ref)). Every count is a draw of the fitted model's own
future observations: the renewal run on, each stream's delays and
ascertainment, and each stream's own likelihood. Columns:

- `:cases_cum`, `:deaths_cum`: cumulative suspected reported cases and
  deaths by the cut-off plus the horizon, the observed cut-off cumulative
  plus the drawn new counts.
- `:confirmed_cum`, `:confirmed_deaths_cum`: the laboratory-confirmed
  counterparts, present when `obs_confirmed` and `obs_confirmed_deaths` are
  supplied.
- `:cases_new`, … `:confirmed_deaths_new`: new counts over the horizon.
- `:isolation_level`: the reported isolation occupancy on the last day, the
  censored negative binomial around the modelled demand plus the
  reclassification offset. `:bed_demand` is the modelled demand that day
  and `:bed_shortfall` its excess over the modelled bed capacity.
- `:admissions_fc`, `:incare_deaths_fc`, `:ruleouts_fc`: the daily
  isolation flows on the last day.
- `:recovered_cum`, `:recovered_new`: recovered among confirmed. The
  cut-off cumulative is `obs_recovered` when supplied, otherwise the fitted
  `expected_recovered_T`, and `recovered_new` needs `obs_recovered`.
- `:infections_new`, `:onsets_new`, `:deaths_latent_new`: new latent
  infections, symptom onsets and deaths over the horizon, which carry
  parameter uncertainty only.
- `:rt_forecast`: the reproduction number on the last day.
- the [`forecast_onsets`](@ref) columns, when the fit carries the onset
  triangle and `horizon` is a whole number of weeks.
"""
function forecast_reported(
        pp;
        horizon::Integer = 7,
        obs_cases::Real,
        obs_deaths::Real,
        obs_confirmed::Union{Real, Missing} = missing,
        obs_confirmed_deaths::Union{Real, Missing} = missing,
        obs_recovered::Union{Real, Missing} = missing
    )
    h = Int(horizon)
    _check_horizon(pp, h)
    new(key) = _forecast_new(pp, key, h)
    at(key) = _forecast_at(pp, key, h)
    cases_new = new("forecast_reports.increments")
    deaths_new = new("forecast_deaths.increments")
    (isnothing(cases_new) || isnothing(deaths_new)) && throw(
        ArgumentError(
            "forecast_reported reads the joint's draws; for a single-stream " *
                "fit use `forecast_stream`."
        )
    )
    df = DataFrame(
        cases_cum = round(Int, obs_cases) .+ cases_new,
        deaths_cum = round(Int, obs_deaths) .+ deaths_new,
        cases_new = cases_new,
        deaths_new = deaths_new,
        infections_new = new("forecast_infections"),
        onsets_new = new("forecast_onsets"),
        deaths_latent_new = new("forecast_latent_deaths"),
        rt_forecast = at("forecast_rt")
    )
    conf = new("forecast_confirmed.increments")
    if obs_confirmed !== missing && !isnothing(conf)
        df.confirmed_cum = round(Int, obs_confirmed) .+ conf
        df.confirmed_new = conf
    end
    conf_deaths = new("forecast_confirmed_deaths.increments")
    if obs_confirmed_deaths !== missing && !isnothing(conf_deaths)
        df.confirmed_deaths_cum = round(Int, obs_confirmed_deaths) .+
            conf_deaths
        df.confirmed_deaths_new = conf_deaths
    end
    iso = at("forecast_isolation.obs")
    if !isnothing(iso)
        df.bed_demand = round.(Int, at("forecast_bed_demand"))
        df.isolation_level = round.(Int, iso)
        df.bed_shortfall = max.(
            df.bed_demand .- round.(Int, at("forecast_bed_capacity")), 0
        )
        df.admissions_fc = round.(Int, at("forecast_admissions.obs"))
        df.incare_deaths_fc = at("forecast_incare_deaths.increments")
        df.ruleouts_fc = at("forecast_ruleouts.increments")
    end
    rec = new("forecast_recovered.increments")
    if !isnothing(rec)
        base = obs_recovered === missing ?
            round.(Int, _draws(pp, :expected_recovered_T)) :
            round(Int, obs_recovered)
        df.recovered_cum = base .+ rec
        obs_recovered !== missing && (df.recovered_new = rec)
    end
    onsets = forecast_onsets(pp; horizon = h, strict = false)
    if !isnothing(onsets)
        for col in propertynames(onsets)
            col in propertynames(df) && continue
            df[!, col] = onsets[!, col]
        end
    end
    return df
end

"""
    forecast_onsets(pp; horizon = 7, obs_value = missing, strict = true)
        -> DataFrame

Nowcast and forecast of symptom onsets and of the reporting triangle that
observes them, one row per draw, read from the posterior-predictive draws
`pp` of a fit carrying the triangle ([`forecast_draws`](@ref)). It
separates cases that have had their onset but are not yet reported from
cases whose onset has not happened. The model runs the fitted reporting
hazard, ascertainment and scan noise forward to a future vintage
`horizon` days past the cut-off, which must be a whole number of weeks
(see [`onset_forecast_model`](@ref)). Columns, all per draw:

  - `onsets_to_date`: symptom onsets that have happened by the cut-off.
  - `onset_reports_to_date`: of those, the number the triangle should
    already have printed.
  - `onsets_unreported`: their difference, the reporting backlog and the
    cases ascertainment will never pick up.
  - `onsets_new`: onsets over the horizon.
  - `onset_reports_backfill`: reports arriving over the horizon of onsets
    on or before the cut-off.
  - `onset_reports_future`: reports arriving over the horizon of onsets
    after the cut-off.
  - `onset_reports_new`: the drawn new reported count, the scored
    forecast, rounded and floored at zero as the observed increment is.
  - `onset_reports_cum`: `obs_value + onset_reports_new` when `obs_value`
    (the triangle's own total at the cut-off) is supplied, so the forecast
    can be plotted on the observed scale.

With `strict = false` it returns `nothing` rather than erroring when the
draws carry no onset forecast or `horizon` is not a future vintage.
"""
function forecast_onsets(
        pp; horizon::Integer = 7, obs_value::Union{Real, Missing} = missing,
        strict::Bool = true
    )
    h = Int(horizon)
    H = _check_horizon(pp, h)
    draws = _forecast_vectors(pp, "forecast_onset_reports.increments")
    j = _vintage_index(h, H)
    if isnothing(draws) || isnothing(j)
        strict || return nothing
        throw(
            ArgumentError(
                isnothing(draws) ?
                    "these draws carry no onset forecast; the fit has no " *
                    "onset triangle." :
                    "horizon $h is not a forecast vintage; the onset " *
                    "forecast is drawn a whole number of weeks ahead."
            )
        )
    end
    cum = _draw_vectors(pp, :cumulative_onsets)
    to_date = [c[length(c) - H] for c in cum]
    reported = _draws(pp, :forecast_onset_reports_to_date)
    df = DataFrame(
        onsets_to_date = to_date,
        onset_reports_to_date = reported,
        onsets_unreported = to_date .- reported,
        onsets_new = _forecast_new(pp, "forecast_onsets", h),
        onset_reports_backfill = [
            x[j] for x in _draw_vectors(pp, :forecast_onset_reports_backfill)
        ],
        onset_reports_future = [
            x[j] for x in _draw_vectors(pp, :forecast_onset_reports_future)
        ],
        onset_reports_new = [max(round(Int, x[j]), 0) for x in draws]
    )
    obs_value === missing ||
        (df.onset_reports_cum = round(Int, obs_value) .+ df.onset_reports_new)
    return df
end

## The streams `forecast_stream` reads: the future-count key each composer
## draws and whether the stream is a level read on the last day rather than
## a count summed over the horizon.
const _STREAM_FORECAST = Dict{Symbol, NamedTuple}(
    :reported_cases => (key = "forecast_reports.increments", level = false),
    :suspected_deaths => (key = "forecast_deaths.increments", level = false),
    :confirmed_cases => (key = "forecast_confirmed.increments", level = false),
    :confirmed_deaths => (
        key = "forecast_confirmed_deaths.increments", level = false,
    ),
    :recovered => (key = "forecast_recovered.increments", level = false),
    :exports => (key = "forecast_exports.counts", level = false),
    :isolation_beds => (key = "forecast_isolation.obs", level = true),
    :onset_reports => (key = "forecast_onset_reports.increments", level = false),
)

"""
    forecast_stream(pp, stream; horizon = 7) -> Vector

The forecast of one observed stream `horizon` days past the cut-off, one
value per draw, from the posterior-predictive draws `pp` of either the joint
or the stream's own single-stream fit ([`forecast_draws`](@ref)). Every
composer draws its streams' future counts under the same names, so one
reader serves both.

`stream` is one of `:reported_cases`, `:suspected_deaths`,
`:confirmed_cases`, `:confirmed_deaths`, `:recovered`, `:exports`,
`:isolation_beds` and `:onset_reports`. The incident streams return the new
count over the horizon, the convention [`forecast_archive`](@ref) scores.
`:isolation_beds` returns the reported occupancy on the last day.
`:onset_reports` returns the new count the triangle prints over the horizon
([`forecast_onsets`](@ref)), which needs `horizon` to be a whole number of
weeks.
"""
function forecast_stream(pp, stream::Symbol; horizon::Integer = 7)
    spec = get(_STREAM_FORECAST, stream, nothing)
    isnothing(spec) && throw(
        ArgumentError(
            "forecast_stream: unknown stream `:$stream`; expected one of " *
                join(sort!([":$s" for s in keys(_STREAM_FORECAST)]), ", ")
        )
    )
    h = Int(horizon)
    _check_horizon(pp, h)
    stream === :onset_reports &&
        return forecast_onsets(pp; horizon = h).onset_reports_new
    _has_key(pp, Symbol(spec.key)) || throw(
        ArgumentError(
            "forecast_stream: these draws carry no `$(spec.key)`; the fit " *
                "does not include `:$stream`."
        )
    )
    spec.level && return round.(Int, _forecast_at(pp, spec.key, h))
    return _forecast_new(pp, spec.key, h)
end

"""
    PROVINCE_FORECAST_METHOD

The `method` [`province_forecast_archive`](@ref) records on each row, and the
only one the release scoring scores. It names the version of
[`forecast_provinces`](@ref), so archives of an earlier projection are never
scored alongside it.
"""
const PROVINCE_FORECAST_METHOD = "predict"

"""
    forecast_provinces(pp; horizon = 7, n_patches, patch_labels) -> DataFrame

The per-province forecast `horizon` days past the cut-off, read from the
posterior-predictive draws `pp` of the patch joint
([`forecast_draws`](@ref)). Each province's reproduction number continues
the fitted national walk and the province's correlated, mean-reverting
deviation, and each province renews with the fitted importation. Its
confirmed cases and confirmed deaths are the national forecast counts split
each week by the fitted province compositions, with the province's own
delays, relative ascertainment and severity, so the provinces add up to the
national forecast in every draw (see [`bvd_joint`](@ref)).

Returns one row per province and draw, with columns `patch`, `province`
(its label in `patch_labels`), `draw`, `infections_new`, `rt_forecast` (the
province's reproduction number on the last day), and `confirmed_new` and
`confirmed_deaths_new` when the fit carries the matching composition. The
counts need `horizon` to be a whole number of weeks.
"""
function forecast_provinces(
        pp;
        horizon::Integer = 7,
        n_patches::Integer = length(PROVINCE_NAMES),
        patch_labels::AbstractVector = PROVINCE_LABELS
    )
    h = Int(horizon)
    H = _check_horizon(pp, h)
    np = min(n_patches, length(patch_labels))
    inf = _forecast_vectors(pp, "forecast_infections_patch")
    (isnothing(inf) || length(first(inf)) != n_patches * H) && throw(
        ArgumentError(
            "these draws carry no $(n_patches)-patch forecast; draw them " *
                "from `bvd_joint` with $(n_patches) patches."
        )
    )
    rt = _draw_vectors(pp, :forecast_rt_patch)
    ## Weekly province counts, one column per future vintage, summed to
    ## the vintage `h` days on.
    j = _vintage_index(h, H)
    function weekly(key)
        v = _forecast_vectors(pp, key)
        isnothing(v) && return nothing
        isnothing(j) && throw(
            ArgumentError(
                "horizon $h is not a forecast vintage; the province " *
                    "counts are split a whole number of weeks ahead."
            )
        )
        return [vec(sum(reshape(x, n_patches, :)[:, 1:j]; dims = 2)) for x in v]
    end
    conf = weekly("forecast_province_confirmed")
    conf_deaths = weekly("forecast_province_deaths")
    nd = length(inf)
    out = DataFrame(
        patch = Int[], province = String[], draw = Int[],
        infections_new = Float64[], rt_forecast = Float64[]
    )
    for i in 1:nd
        I = reshape(inf[i], n_patches, :)
        R = reshape(rt[i], n_patches, :)
        for p in 1:np
            push!(
                out, (
                    p, String(patch_labels[p]), i,
                    sum(@view I[p, 1:h]), R[p, h],
                )
            )
        end
    end
    isnothing(conf) ||
        (out.confirmed_new = [conf[i][p] for i in 1:nd for p in 1:np])
    isnothing(conf_deaths) || (
        out.confirmed_deaths_new = [
            conf_deaths[i][p] for i in 1:nd for p in 1:np
        ]
    )
    return out
end

"""
    forecast_archive(fcs; made_date, thin = 1) -> DataFrame

Long-format archive of one or more [`forecast_reported`](@ref) results made
from a single cut-off, for the observed streams scored across releases. `fcs`
is an iterable of `(horizon, fc)` pairs and `made_date` is the cut-off `Date`
the forecasts were made from. Returns one row per `(stream, horizon, draw)`
with columns `made_date`, `horizon`, `target_date` (`made_date` plus the
horizon), `stream`, `draw` and `value`.

Only the incident and level quantities are archived: `confirmed cases` and
`confirmed deaths` new over the horizon, `recovered` new over the horizon,
`onset reports` new over the horizon, and the supply-limited `isolation beds`
occupancy. The cumulative totals are revised across data vintages, so they
are not archived. The reporting triangle's own total is revised by more than
the rest, since the ≈4% per-scan level error alone moves it by tens of cases
and the printed total falls between consecutive vintages (2531 to 2523, and
2018 to 1996), which late reporting cannot produce.

When the forecast carries the confirmed/suspect ward split
(`confirmed_occupancy` / `suspect_occupancy`, the occupancy partitioned by
the cut-off confirmed share), the two ward occupancy levels are archived too
as `treatment beds` (confirmed) and `isolation beds (suspected)`, each scored
against its own Tableau 6 occupancy sub-stock. The total `isolation beds`
occupancy stays as its own stream. Nothing produces those two columns yet:
[`forecast_reported`](@ref) projects the total occupancy alone, so the two
ward entries are reachable only from a forecast a caller has partitioned
itself, which today is the tests. Streams a forecast does not carry are
skipped, so this costs a real forecast nothing. `thin` keeps every
`thin`-th draw so the archive stays compact when it is saved as a release
asset.
"""
function forecast_archive(fcs; made_date::Date, thin::Integer = 1)
    streams = (
        (:confirmed_new, "confirmed cases"),
        (:confirmed_deaths_new, "confirmed deaths"),
        (:recovered_new, "recovered"),
        (:isolation_level, "isolation beds"),
        (:confirmed_occupancy, "treatment beds"),
        (:suspect_occupancy, "isolation beds (suspected)"),
        (:onset_reports_new, "onset reports"),
    )
    out = DataFrame(
        made_date = Date[], horizon = Int[], target_date = Date[],
        stream = String[], draw = Int[], value = Float64[]
    )
    for (horizon, fc) in fcs
        h = Int(horizon)
        target = made_date + Day(h)
        for (col, label) in streams
            col in propertynames(fc) || continue
            vals = fc[!, col]
            for (d, i) in enumerate(1:thin:length(vals))
                push!(
                    out,
                    (made_date, h, target, label, d, Float64(vals[i]))
                )
            end
        end
    end
    return out
end

"""
    province_forecast_archive(pp, fcs; made_date, thin = 1) -> DataFrame

Long-format archive of one or more per-province forecasts made from a single
cut-off, in the [`forecast_archive`](@ref) schema plus a `province` column.
`pp` is the patch joint's posterior-predictive draws
([`forecast_draws`](@ref)), `fcs` an iterable of `(horizon, fc)` pairs and
`made_date` the cut-off `Date`. Returns one row per
`(province, stream, horizon, draw)` with columns `made_date`, `horizon`,
`target_date`, `province`, `stream`, `draw`, `value` and `method`.
`method` is [`PROVINCE_FORECAST_METHOD`](@ref), so scoring can tell these
rows from earlier archives: the hand-written projection wrote
`"projection"`, and the share split before it wrote no `method` column.

The two incident streams the spatial tables report are archived, under the
same `confirmed cases` and `confirmed deaths` labels `forecast_archive`
gives the national streams. `province` is the patch
key from [`PROVINCE_NAMES`](@ref), which [`PROVINCE_MEMBERS`](@ref) maps to
the source provinces a patch pools, so a scorer can build each patch's truth
from the per-province histories in the same release's `observations.toml`.

Each province's values are the [`forecast_provinces`](@ref) forecast from
`pp` at that horizon. An `fc` that is already a province forecast is
archived as it is. A national [`forecast_reported`](@ref) result is replaced
by the province forecast from the same draws. `thin` keeps
every `thin`-th draw so the archive stays compact as a release asset.
"""
function province_forecast_archive(
        pp, fcs; made_date::Date,
        n_patches::Integer = length(PROVINCE_NAMES),
        patch_labels::AbstractVector = PROVINCE_NAMES,
        thin::Integer = 1
    )
    np = min(n_patches, length(patch_labels))
    out = DataFrame(
        made_date = Date[], horizon = Int[], target_date = Date[],
        province = String[], stream = String[], draw = Int[],
        value = Float64[], method = String[]
    )
    for (horizon, fc) in fcs
        h = Int(horizon)
        target = made_date + Day(h)
        for (label, province, vals) in _province_forecast_draws(
                pp, fc, np, patch_labels; horizon = h
            )
            for (d, i) in enumerate(1:thin:length(vals))
                push!(
                    out,
                    (
                        made_date, h, target, province, label, d,
                        Float64(vals[i]), PROVINCE_FORECAST_METHOD,
                    )
                )
            end
        end
    end
    return out
end

"""
Summarise a [`forecast_reported`](@ref) result into a `DataFrame` with
one row per confirmed stream (laboratory-confirmed cases and confirmed
deaths) and quantity (cumulative total by the cut-off plus the horizon, or
new this week), reporting the equal-tailed 30/60/90% credible interval
endpoints (`lower_90 … upper_90`) used by the other summary tables. The
suspected reported-case and suspected-death streams are no longer reported,
so they are not shown as forecast targets.
"""
function forecast_table(fc::DataFrame; digits::Integer = 0)
    _row(
        label,
        quantity,
        draws
    ) = begin
        s = posterior_summary(draws)
        (
            stream = label, quantity = quantity,
            lower_90 = round(s.lo90; digits), lower_60 = round(s.lo60; digits),
            lower_30 = round(s.lo30; digits), upper_30 = round(s.hi30; digits),
            upper_60 = round(s.hi60; digits), upper_90 = round(s.hi90; digits),
        )
    end
    streams = Tuple{String, Symbol, Symbol}[]
    :confirmed_cum in propertynames(fc) && push!(
        streams,
        ("DRC confirmed cases", :confirmed_cum, :confirmed_new)
    )
    :confirmed_deaths_cum in propertynames(fc) && push!(
        streams,
        (
            "DRC confirmed deaths", :confirmed_deaths_cum,
            :confirmed_deaths_new,
        )
    )
    rows = NamedTuple[]
    for (label, cum, new) in streams
        push!(rows, _row(label, "cumulative by T+7", fc[!, cum]))
        push!(rows, _row(label, "new this week", fc[!, new]))
    end
    ## Both are levels at the horizon. The gap between them is the shortfall.
    :bed_demand in propertynames(fc) && push!(
        rows,
        _row("DRC isolation beds", "demand at T+7", fc[!, :bed_demand])
    )
    :isolation_level in propertynames(fc) && push!(
        rows,
        _row("DRC isolation beds", "occupancy at T+7", fc[!, :isolation_level])
    )
    ## One-week-ahead daily isolation/treatment flows (a single-day rate at
    ## the horizon, not a cumulative total).
    :admissions_fc in propertynames(fc) && push!(
        rows,
        _row("DRC isolation admissions", "daily at T+7", fc[!, :admissions_fc])
    )
    :incare_deaths_fc in propertynames(fc) && push!(
        rows,
        _row("DRC in-care deaths", "daily at T+7", fc[!, :incare_deaths_fc])
    )
    :ruleouts_fc in propertynames(fc) && push!(
        rows,
        _row("DRC isolation rule-outs", "daily at T+7", fc[!, :ruleouts_fc])
    )
    if :recovered_cum in propertynames(fc)
        push!(
            rows,
            _row(
                "DRC recovered among confirmed", "cumulative by T+7",
                fc[!, :recovered_cum]
            )
        )
        :recovered_new in propertynames(fc) && push!(
            rows,
            _row(
                "DRC recovered among confirmed", "new this week",
                fc[!, :recovered_new]
            )
        )
    end
    return _prettify(DataFrame(rows))
end

"""
    onset_forecast_table(fc; digits = 0) -> DataFrame

Summarise a [`forecast_onsets`](@ref) result into the report's usual
30/60/90% credible-interval table, one row per quantity, ordered so the
nowcast is read before the forecast. The `Quantity` labels name what each
row is a statement about rather than which column it came from, since
"onsets" and "onsets reported" are different numbers:

  - `symptom onsets to date` and `of those, reported by T`: the nowcast
    pair. Their difference is the next row.
  - `onsets not yet reported at T`: already happened, not yet in the
    figure (reporting backlog and never-ascertained cases together, see
    [`forecast_onsets`](@ref)).
  - `reports this week of onsets before T`: the backlog the coming week
    should clear.
  - `reports this week of onsets after T`: reports of cases that have not
    yet had their symptom onset.
  - `new onset reports this week`: the replicated total of the two, the
    quantity scored against the triangle's own increment.
  - `new symptom onsets this week`: the latent forecast, unobserved.
"""
function onset_forecast_table(fc::DataFrame; digits::Integer = 0)
    _row(quantity, draws) = begin
        s = posterior_summary(draws)
        (
            quantity = quantity,
            lower_90 = round(s.lo90; digits), lower_60 = round(s.lo60; digits),
            lower_30 = round(s.lo30; digits), upper_30 = round(s.hi30; digits),
            upper_60 = round(s.hi60; digits), upper_90 = round(s.hi90; digits),
        )
    end
    rows = [
        _row("symptom onsets to date", fc.onsets_to_date),
        _row("of those, reported by T", fc.onset_reports_to_date),
        _row("onsets not yet reported at T", fc.onsets_unreported),
        _row(
            "reports this week of onsets before T",
            fc.onset_reports_backfill
        ),
        _row("reports this week of onsets after T", fc.onset_reports_future),
        _row("new onset reports this week", fc.onset_reports_new),
        _row("new symptom onsets this week", fc.onsets_new),
    ]
    return _prettify(DataFrame(rows))
end

"""
Validate a [`forecast_reported`](@ref) projection against the counts that
were later observed. `observed` is a `NamedTuple` mapping each stream's
cumulative column (`:confirmed_cum`, `:cases_cum`, …) to its observed
cumulative count at the forecast target date. `baseline` maps the same
columns to the cumulative count at the forecast origin (default `0`). A
stream is scored only when its cumulative column is in the forecast and a key
for it is in `observed`. The scored streams are the reported cases, suspected
deaths, laboratory-confirmed cases, confirmed deaths and recovered (the same
set [`plot_forecast_vs_truth`](@ref) draws), so the call site can pass the
same `observed`/`baseline` NamedTuples it builds for the plot.

Each scored stream gets two rows, mirroring the plot's two panels and the
`Quantity` split of [`forecast_table`](@ref): a `cumulative by T+7` row
scoring the projected cumulative against `observed`, and a `new this week`
row scoring the projected new count against
`max(observed − breaks − baseline, 0)`. `breaks` is keyed like `observed`
and carries each stream's retrospective harmonisation correction over the
forecast window (see [`confirmed_break_correction`](@ref)), defaulting to
zero. It comes out of both truths, because it lands in the reported
cumulative without having been notified in that week and a projection cannot
contain it.

When `isolation` (the observed bed occupancy at the target date) is supplied
and the forecast carries the beds, the projected supply-limited occupancy is
scored against it too as a single level row. Returns a `DataFrame` with the
observed count, the equal-tailed 30/60/90% predictive intervals (the same
endpoints as the other summary tables), and whether the observed count falls
inside the 90% interval.

Note that at a one-week-back freeze the bed capacity is weakly informed (the
reported occupancy rate starts only on 9 June), so the projected bed
occupancy rides the capacity random walk back to the freeze date and its
interval is wide.
"""
function forecast_vs_truth(
        fc::DataFrame;
        observed::NamedTuple, baseline::NamedTuple = NamedTuple(),
        breaks::NamedTuple = NamedTuple(),
        isolation::Union{Real, Missing} = missing,
        digits::Integer = 0
    )
    _row(label, quantity, draws, obs) = begin
        s = posterior_summary(draws)
        lo = round(s.lo90; digits)
        hi = round(s.hi90; digits)
        (
            stream = label, quantity = quantity, observed = round(obs; digits),
            lower_90 = lo, lower_60 = round(s.lo60; digits),
            lower_30 = round(s.lo30; digits), upper_30 = round(s.hi30; digits),
            upper_60 = round(s.hi60; digits), upper_90 = hi,
            within_90 = lo <= obs <= hi ? "yes" : "no",
        )
    end
    specs = (
        (:cases_cum, :cases_new, "DRC reported cases"),
        (:deaths_cum, :deaths_new, "DRC suspected deaths"),
        (:confirmed_cum, :confirmed_new, "DRC confirmed cases"),
        (:confirmed_deaths_cum, :confirmed_deaths_new, "DRC confirmed deaths"),
        (:recovered_cum, :recovered_new, "DRC recovered among confirmed"),
    )
    rows = NamedTuple[]
    for (cumcol, newcol, label) in specs
        (cumcol in propertynames(fc) && haskey(observed, cumcol)) || continue
        brk = float(get(breaks, cumcol, 0))
        obs_cum = float(observed[cumcol]) - brk
        push!(rows, _row(label, "cumulative by T+7", fc[!, cumcol], obs_cum))
        newcol in propertynames(fc) || continue
        obs_new = max(obs_cum - float(get(baseline, cumcol, 0)), 0.0)
        push!(rows, _row(label, "new this week", fc[!, newcol], obs_new))
    end
    isolation !== missing && :isolation_level in propertynames(fc) &&
        push!(
        rows, _row(
            "DRC isolation beds", "occupancy at T+7",
            fc[!, :isolation_level], isolation
        )
    )
    return _prettify(DataFrame(rows))
end

## Does the chain carry `key`? `FlexiChains` throws rather than returning
## `nothing` for an absent key, so the lookup is probed.
function _has_key(chn, key)
    try
        chn[key]
        return true
    catch
        return false
    end
end

## Draws of the first key in `candidates` the chain carries, or `nothing`
## when it carries none. Single-stream fits bind each stream's submodel under
## a composer-level name (`cases_state.expected_reports`), and the joint
## carries those nested names alongside its own un-prefixed aliases
## (`expected_reports_T`), so one ordered list serves both fit kinds.
function _resolve_draws(chn, candidates)
    for key in candidates
        _has_key(chn, key) && return _draws(chn, key)
    end
    return nothing
end
