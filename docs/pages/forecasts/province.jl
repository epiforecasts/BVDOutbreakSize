# # Province forecasts
#
# This page projects each province a week ahead from the joint model's fit.
# Each province's renewal equation runs on past the cut-off, with the provinces still exchanging infections through importation.
# Its reproduction number follows the national trend's walk, and its deviation from that trend reverts toward zero at the fitted half-life.
# The confirmed counts start from the national daily rate at the cut-off times the province's modelled share at the most recent spatial vintage, and then grow with the province's own projected infections.
# The provinces are projected separately, so they need not add up to the national forecast on the [forecasts](@ref "Forecasts") page.
# How well the model reproduces each province's share is on the [in-sample checks](@ref province-compositions) page.
# The [forecast by province](@ref "Forecast by province") evaluation scores this projection.
#
# Releases before 24 September 2026 archived a different province forecast, the national forecast split by each province's share at the last spatial vintage and held over the week.
# Releases from that date archive this projection.
# Scores of past province forecasts, and the past forecasts shown below, therefore mix the two methods.
#
# The spatial tables report confirmed cases and confirmed deaths, so those are the observed streams projected.
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

# ## Province forecast summary
#
# The overall bullets compare the provinces.
# The detail under each province gives its own projection.

#md # ```@raw html
#md # <details><summary>Project each province a week ahead</summary>
#md # ```

province_projection = forecast_provinces(
    chn_joint;
    horizon = 7, n_patches = N_PATCHES
);
province_forecast = province_forecast_table(
    chn_joint, province_projection;
    n_patches = N_PATCHES
);
province_forecast_fig = plot_province_forecast(
    chn_joint, province_projection;
    n_patches = N_PATCHES
);
province_week_end = obs.cutoff + Day(7);
province_summary_markdown = let proj = province_projection
    draws(p, col) = float.(proj[proj.patch .== p, col])
    pct(x) = string(round(Int, 100 * x), "%")
    share_text(v) = string(
        "about ", pct(quantile(v, 0.5)), " (90% credible interval ",
        pct(quantile(v, 0.05)), " to ", pct(quantile(v, 0.95)), ")"
    )
    ## Overall: the provinces ranked by their median projection, each with
    ## its share of the provinces' combined projection and how often it
    ## projects the most.
    function overall(col, stream)
        shares = province_share_draws(proj, col; n_patches = N_PATCHES)
        top = [argmax([s[k] for s in shares]) for k in eachindex(shares[1])]
        order = sortperm(
            [quantile(draws(p, col), 0.5) for p in 1:N_PATCHES]; rev = true
        )
        lead = first(order)
        ranked = join(PROVINCE_LABELS[order], ", then ")
        split_text = join(
            [
                "$(PROVINCE_LABELS[p]) $(share_text(shares[p]))"
                    for p in order
            ], "; "
        )
        return [
            "- **Most $(stream):** $(ranked). " *
                "$(PROVINCE_LABELS[lead]) projects the most in " *
                "$(pct(count(==(lead), top) / length(top))) of draws.",
            "- **Share of $(stream):** $(split_text), of the provinces' " *
                "combined projection.",
        ]
    end
    function detail(p)
        return join(
            [
                "**$(PROVINCE_LABELS[p])**", "",
                "- New confirmed cases: " *
                    "$(median_interval_text(draws(p, :confirmed_new))).",
                "- New confirmed deaths: " *
                    "$(median_interval_text(draws(p, :confirmed_deaths_new))).",
                "- New infections, reported and unreported: " *
                    "$(median_interval_text(draws(p, :infections_new))).",
                "- Reproduction number on $(province_week_end): " *
                    "$(median_interval_text(draws(p, :rt_forecast); digits = 2)).",
            ], "\n"
        )
    end
    join(
        [
            "Projected counts are for the week to $(province_week_end).", "",
            overall(:confirmed_new, "new confirmed cases")...,
            overall(:confirmed_deaths_new, "new confirmed deaths")..., "",
            (detail(p) * "\n" for p in 1:N_PATCHES)...,
        ], "\n"
    )
end;

#md # ```@raw html
#md # </details>
#md # ```

Markdown.parse(province_summary_markdown) #hide

# ## One-week-ahead forecast by province
#
# The table and figure give the new confirmed cases and confirmed deaths expected in each province by $T + 7$.
# The table also gives each province's new infections and its reproduction number at $T + 7$.

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
    chn_joint, province_projection;
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

# ## Past forecasts against what was observed
#
# Each release since 15 September 2026 archives its province split, so this figure has few made dates so far.
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
    scored_overlay(province_overlay_df)
);

#md # ```@raw html
#md # </details>
#md # ```

province_overlay_fig #hide

# ---
#
# The full analysis code, data and model definitions are in the
# [epiforecasts/BVDOutbreakSize](https://github.com/epiforecasts/BVDOutbreakSize)
# repository.
