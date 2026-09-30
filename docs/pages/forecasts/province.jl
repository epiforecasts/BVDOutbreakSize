# # Province forecasts
#
# This page gives the one-week-ahead forecast for each province from the joint model.
# The projection is defined in the [province forecast](@ref "Province forecast") Methods section.
# The national forecast is on the [forecasts](@ref "Forecasts") page.
# How these forecasts have scored against what each province went on to report is in the [forecast by province](@ref "Forecast by province") evaluation.
# How well the model reproduces each province's share is on the [in-sample checks](@ref province-compositions) page.

#md # ```@raw html
#md # <details><summary>Load packages, data and fitted chains</summary>
#md # ```

## Shared setup: packages, observations and the fit registry. See
## `docs/pages/_setup.jl`.
using BVDOutbreakSize
include(joinpath(pkgdir(BVDOutbreakSize), "docs", "pages", "_setup.jl"))
#-
## The fits this page reads, loaded from the cache here.
chn_joint = load_fit("joint");

#md # ```@raw html
#md # </details>
#md # ```

# ## Summary
#
# The table gives each province's forecast for the week to $T + 7$ as a 90% predictive interval, and the bullets under it compare the provinces draw by draw.
# The provinces add up to the national forecast.

#md # ```@raw html
#md # <details><summary>Project each province a week ahead</summary>
#md # ```

## Read from the same draws as the national forecast page, so the provinces
## add up to the national forecast shown there.
province_draws = fit_forecast("joint");
province_projection = forecast_provinces(
    province_draws; horizon = 7, n_patches = N_PATCHES
);
province_forecast = province_forecast_table(
    province_draws, province_projection;
    n_patches = N_PATCHES
);
province_forecast_fig = plot_province_forecast(
    province_draws, province_projection;
    n_patches = N_PATCHES
);
province_week_end = obs.cutoff + Day(7);
province_forecast_md = province_forecast_headline(
    province_projection; n_patches = N_PATCHES
);

#md # ```@raw html
#md # </details>
#md # ```

Markdown.parse(province_forecast_md) #hide

# ## One-week-ahead forecast by province
#
# Each panel is one forecast target, the provinces side by side, for the week to $(province_week_end).

province_forecast_fig #hide

# The maps shade each province by its patch's forecast, so the pooled patch shades all its provinces alike.
# The reproduction number is centred on one, and a province whose 90% predictive interval spans one is washed out.

#md # ```@raw html
#md # <details><summary>Build the forecast maps</summary>
#md # ```

province_forecast_draws(col) = [
    float.(province_projection[province_projection.patch .== p, col])
        for p in 1:N_PATCHES
]
province_forecast_map = plot_province_map(
    filter(
        !isnothing,
        [
            (;
                province_map_summary(province_forecast_draws(:rt_forecast))...,
                title = "Reproduction number at T+7", diverging_at = 1.0,
                colorbar_label = "R (median)",
            ),
            (;
                values = province_map_summary(
                    province_forecast_draws(:confirmed_new)
                ).values,
                title = "New confirmed cases by T+7",
                scale = CairoMakie.Makie.pseudolog10,
                colorbar_label = "Confirmed cases (median)",
            ),
            :isolation_level in propertynames(province_projection) ? (;
                    values = province_map_summary(
                        province_forecast_draws(:isolation_level)
                    ).values,
                    title = "Patients in isolation at T+7",
                    scale = CairoMakie.Makie.pseudolog10,
                    colorbar_label = "Patients (median)",
                ) : nothing,
        ]
    )
);

#md # ```@raw html
#md # </details>
#md # ```

province_forecast_map #hide

#md # ```@raw html
#md # <details><summary>Forecast summary table by province</summary>
#md # ```

MarkdownTable(province_forecast) #hide

#md # ```@raw html
#md # </details>
#md # ```

# ## Forecast for each province
#
# Each province has its own forecast figure, one histogram per target with its 90% predictive interval shaded.
# The dashed rule on the case and death panels is the count the province reported over its most recent week in the spatial tables.

#md # ```@raw html
#md # <details><summary>Build the per-province figures</summary>
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
forecast_province_fig(p) = plot_province_forecast_detail(
    province_draws, province_projection;
    province = p, n_patches = N_PATCHES,
    observed = forecast_province_observed(p)
)
## The province blocks below are written out one per patch.
@assert N_PATCHES == 4 && PROVINCE_LABELS[1:4] ==
    ["Ituri", "Nord-Kivu", "Haut-Uele", "Other provinces"]
## Cases and deaths come from separate tables, so each stream's observed
## week is stated on its own.
function forecast_province_recent_week(r, stream)
    r === nothing &&
        return "The spatial tables carry no recent week of $(stream)."
    return string(
        "The observed week of $(stream) runs from ",
        grid_date(r.start_day), " to ", grid_date(r.last_day), "."
    )
end
recent_note = join(
    [
        forecast_province_recent_week(recent_cases, "confirmed cases"),
        forecast_province_recent_week(recent_deaths, "confirmed deaths"),
        "Each ends at its last spatial vintage, which can be earlier than " *
            "the cut-off.",
    ], " "
);

#md # ```@raw html
#md # </details>
#md # ```

Markdown.parse(recent_note) #hide

# ### Ituri

forecast_province_fig(1) #hide

# ### Nord-Kivu

forecast_province_fig(2) #hide

# ### Haut-Uele

forecast_province_fig(3) #hide

# ### Other provinces
#
# The other provinces are pooled into one patch, so their forecast is for the pool.

forecast_province_fig(4) #hide

# ## Past forecasts against what was observed
#
# Only [province forecast](@ref "Province forecast") projections are shown, so the figure fills in as releases accumulate.
# Each panel is one province stream at one horizon.
# The x-axis is the cut-off each forecast was made from.
# Each forecast shows its median and 90% predictive interval, beside the persistence baseline and the count the province went on to report.
# A window holding a harmonisation-break day is left out, because that day's backfill is published for the country and not by province.

#md # ```@raw html
#md # <details><summary>Load the archived province forecasts and their outcomes</summary>
#md # ```

## Written by `scripts/score_releases.jl` from each release's
## `province_forecast.csv`, in the national overlay's schema. A missing
## file reads as an empty table, which the figure reports as nothing scored.
province_overlay_df = _release_data(
    joinpath("province", "forecast_overlay.csv"),
    (;
        release = String, made_date = Date, stream = String, horizon = Int,
        target_date = Date, fit = String, observed = Float64,
        median = Float64, lo30 = Float64, hi30 = Float64, lo60 = Float64,
        hi60 = Float64, lo90 = Float64, hi90 = Float64,
    )
)
province_overlay_fig = plot_forecast_overlay(
    scored_overlay(province_overlay_df);
    empty_message = "No projection forecast has been scored yet. " *
        "This fills in as releases accumulate."
);

#md # ```@raw html
#md # </details>
#md # ```

province_overlay_fig #hide

# ## Saving province map assets
#
# The [dashboard](@ref "Dashboard") map's province layer and its detail column read the per-province and national estimates and time series written here.

#md # ```@raw html
#md # <details><summary>Write the province map estimates and time series</summary>
#md # ```

## One row per province keyed as the health-zone geojson keys it: the
## patch's reproduction number at the cut-off, its confirmed cases over
## the coming week and its ascertainment, with the province's own confirmed
## cases and deaths. The national row reads the joint's national
## reproduction number and the provinces' summed forecast.
province_map_dir = joinpath(
    pkgdir(BVDOutbreakSize), "docs", "src", "summary_assets"
)
mkpath(province_map_dir)
CSV.write(
    joinpath(province_map_dir, "province_estimates.csv"),
    province_map_estimates(
        chn_joint, province_forecast_draws(:confirmed_new);
        n_patches = N_PATCHES,
        confirmed_history = obs.province_confirmed_history,
        death_history = obs.province_death_history, cutoff = obs.cutoff
    )
);
national_forecast_draws = let d = province_projection
    [float(sum(d.confirmed_new[d.draw .== i])) for i in sort(unique(d.draw))]
end
CSV.write(
    joinpath(province_map_dir, "national_estimates.csv"),
    province_map_estimates(
        [vec(Array(chn_joint[:R_T]))], [national_forecast_draws];
        confirmed_history = Dict("national" => obs.confirmed_history),
        death_history = Dict("national" => obs.confirmed_deaths_history),
        cutoff = obs.cutoff, patch_names = ["national"],
        patch_labels = ["National"],
        members = Dict("national" => ["national"])
    )
);
## The daily reproduction number and the weekly confirmed cases of each
## patch and of the country, keyed by the patch label.
province_ts_rt_args = (;
    n = obs.n, breakpoint = _BREAKPOINT, rt_start = _rt_start_plot,
    rt_walk_start = clamp(_BREAKPOINT - RT_WALK_LEAD, _rt_start_plot, obs.n),
    ramp = RT_INTERVENTION_RAMP,
)
province_ts_inc = province_increment_matrix(
    obs.province_confirmed_history, PROVINCE_NAMES, N_PATCHES
)
national_ts_inc = let c = obs.confirmed_history.counts
    isempty(c) ? zeros(Int, 1, 0) : reshape([c[1]; max.(diff(c), 0)], 1, :)
end
## Tag a long table with the series it holds.
_with_series(t, s) = (t[!, :series] .= s; t)
CSV.write(
    joinpath(province_map_dir, "province_timeseries.csv"),
    vcat(
        _with_series(
            rt_quantile_table(
                [
                    reconstruct_patch_rt(
                        chn_joint; n_patches = N_PATCHES,
                        province_ts_rt_args...
                    );
                    [reconstruct_rt(chn_joint; province_ts_rt_args...)]
                ],
                [PROVINCE_LABELS[1:N_PATCHES]; "National"];
                cutoff = obs.cutoff, n = obs.n, from = _rt_start_plot
            ), "rt"
        ),
        _with_series(
            vcat(
                weekly_count_table(
                    province_ts_inc.days, province_ts_inc.increments,
                    PROVINCE_LABELS[1:N_PATCHES];
                    cutoff = obs.cutoff, n = obs.n
                ),
                weekly_count_table(
                    obs.confirmed_history.days, national_ts_inc,
                    ["National"]; cutoff = obs.cutoff, n = obs.n
                )
            ), "cases"
        );
        cols = :union
    )
);

#md # ```@raw html
#md # </details>
#md # ```

# ---
#
# The full analysis code, data and model definitions are in the
# [epiforecasts/BVDOutbreakSize](https://github.com/epiforecasts/BVDOutbreakSize)
# repository.
