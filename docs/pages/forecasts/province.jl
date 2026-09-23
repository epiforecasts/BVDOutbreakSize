# # Province forecasts
#
# This page splits the one-week-ahead national forecast on the [forecasts](@ref "Forecasts") page by province.
# Each province's count is the national draw times that province's modelled share at the most recent spatial vintage.
# The share is held at that value over the week.
# How well the model reproduces each province's share is on the [in-sample checks](@ref province-compositions) page.
# How the split has scored against what each province went on to report is in the [forecast by province](@ref "Forecast by province") evaluation.
#
# The spatial tables report confirmed cases and confirmed deaths, so the split covers those two streams.
# The symptom-onset curve is national only, so there is no province nowcast.
#
# This page is generated from
# [`docs/pages/forecasts/province.jl`](https://github.com/epiforecasts/BVDOutbreakSize/blob/main/docs/pages/forecasts/province.jl).
# The model code it calls is in
# [`src/`](https://github.com/epiforecasts/BVDOutbreakSize/tree/main/src).

#md # ```@raw html
#md # <details><summary>Load packages, data and fitted chains</summary>
#md # ```

## Shared setup: packages, observations, the fit registry and every model fit
## (loaded from the content-addressed cache). See `docs/pages/_setup.jl`.
using BVDOutbreakSize
include(joinpath(pkgdir(BVDOutbreakSize), "docs", "pages", "_setup.jl"))

#md # ```@raw html
#md # </details>
#md # ```

# ## One-week-ahead forecast by province
#
# The table and figure give the new confirmed cases and confirmed deaths expected in each province by $T + 7$.

#md # ```@raw html
#md # <details><summary>Generate the one-week-ahead forecast and its province split</summary>
#md # ```

## The same call as the national page, so the split is of the national
## forecast shown there.
forecast = forecast_reported(
    chn_joint;
    horizon = 7,
    obs_cases = obs.reported_cases,
    obs_deaths = obs.total_deaths,
    obs_confirmed = obs.confirmed_cases,
    obs_confirmed_deaths = obs.confirmed_deaths,
    obs_recovered = obs.recovered_cases
);
province_forecast = province_forecast_table(
    chn_joint, forecast;
    n_patches = N_PATCHES
);
province_forecast_fig = plot_province_forecast(
    chn_joint, forecast;
    n_patches = N_PATCHES
);

#md # ```@raw html
#md # </details>
#md # ```

#md # ```@raw html
#md # <details><summary>Province forecast summary table</summary>
#md # ```

MarkdownTable(province_forecast) #hide

#md # ```@raw html
#md # </details>
#md # ```

province_forecast_fig #hide

# ## Forecast for each province
#
# Each province below has its own summary table and forecast figure.
# The figure histograms the new count over the week, with the 90% predictive interval shaded.
# The dashed rule is the new count the province reported over its most recent week in the spatial tables, for comparison.

#md # ```@raw html
#md # <details><summary>Build the per-province tables and figures</summary>
#md # ```

## The spatial tables' most recent week in each province, clamped as the
## compositions are, so a downward revision reads as no new cases.
recent_cases = province_recent_counts(
    obs.province_confirmed_history, PROVINCE_NAMES, N_PATCHES
)
recent_deaths = province_recent_counts(
    obs.province_death_history, PROVINCE_NAMES, N_PATCHES
)
function forecast_province_observed(p)
    o = (;)
    recent_cases === nothing ||
        (o = merge(o, (; confirmed_new = recent_cases.counts[p])))
    recent_deaths === nothing ||
        (o = merge(o, (; confirmed_deaths_new = recent_deaths.counts[p])))
    return o
end
forecast_province_table(p) = MarkdownTable(
    province_forecast[
        province_forecast.Province .== PROVINCE_LABELS[p], :,
    ]
)
forecast_province_fig(p) = plot_province_forecast_detail(
    chn_joint, forecast;
    province = p, n_patches = N_PATCHES, observed = forecast_province_observed(p)
)
## The province blocks below are written out one per patch.
@assert N_PATCHES == 4 && PROVINCE_LABELS[1:4] ==
    ["Ituri", "Nord-Kivu", "Haut-Uele", "Other provinces"]
recent_note = recent_cases === nothing ?
    "The spatial tables carry no recent week to compare against." :
    string(
        "The observed week runs from ", grid_date(recent_cases.start_day),
        " to ", grid_date(recent_cases.last_day),
        ", the last spatial vintage, which can be earlier than the cut-off."
    );

#md # ```@raw html
#md # </details>
#md # ```

Markdown.parse(recent_note) #hide

#md # ```@raw html
#md # <details><summary>Ituri</summary>
#md # ```

forecast_province_table(1) #hide

#-

forecast_province_fig(1) #hide

#md # ```@raw html
#md # </details>
#md # ```

#md # ```@raw html
#md # <details><summary>Nord-Kivu</summary>
#md # ```

forecast_province_table(2) #hide

#-

forecast_province_fig(2) #hide

#md # ```@raw html
#md # </details>
#md # ```

#md # ```@raw html
#md # <details><summary>Haut-Uele</summary>
#md # ```

forecast_province_table(3) #hide

#-

forecast_province_fig(3) #hide

#md # ```@raw html
#md # </details>
#md # ```

# The other provinces are pooled into one patch, so their forecast is for the pool.

#md # ```@raw html
#md # <details><summary>Other provinces</summary>
#md # ```

forecast_province_table(4) #hide

#-

forecast_province_fig(4) #hide

#md # ```@raw html
#md # </details>
#md # ```

# ---
#
# The full analysis code, data and model definitions are in the
# [epiforecasts/BVDOutbreakSize](https://github.com/epiforecasts/BVDOutbreakSize)
# repository.
