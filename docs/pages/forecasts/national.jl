# # Forecasts
#
# Every release forecasts each DRC stream a week ahead from the joint posterior, as the [one-week-ahead forecast](@ref "One-week-ahead forecast") Methods section describes.
# How these forecasts have scored against the data that arrived afterwards is on the [evaluation](@ref "Forecast evaluation") page.
# The split by province is on the [province forecasts](@ref "Province forecasts") page and the split by health zone on the [health-zone forecasts](@ref "Health-zone forecasts") page.

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
# These are the expected counts for the week after the cut-off.

#md # ```@raw html
#md # <details><summary>Generate the one-week-ahead forecast</summary>
#md # ```

forecast = forecast_reported(
    fit_forecast("joint");
    horizon = 7,
    obs_cases = obs.reported_cases,
    obs_deaths = obs.total_deaths,
    obs_confirmed = obs.confirmed_cases,
    obs_confirmed_deaths = obs.confirmed_deaths,
    obs_recovered = obs.recovered_cases
);
forecast_week_end = obs.cutoff + Day(7);
national_forecast_bullets = join(
    [
        "- **Confirmed cases:** $(median_interval_text(forecast.confirmed_new)) new laboratory-confirmed cases in the week to $(forecast_week_end).",
        "- **Confirmed deaths:** $(median_interval_text(forecast.confirmed_deaths_new)) new confirmed deaths over the same week.",
        "- **Infections:** $(median_interval_text(forecast.infections_new)) new infections, reported and unreported.",
        "- **Reproduction number:** $(median_interval_text(forecast.rt_forecast; digits = 2)) on $(forecast_week_end).",
    ], "\n"
);

#md # ```@raw html
#md # </details>
#md # ```

Markdown.parse(national_forecast_bullets) #hide

# ## One-week-ahead forecast results

#md # ```@raw html
#md # <details><summary>Summarise the one-week-ahead forecast</summary>
#md # ```

forecast_summary = forecast_table(forecast);

#md # ```@raw html
#md # </details>
#md # <figure class="bvd-figure">
#md # ```

forecast_summary #hide

#md # ```@raw html
#md # <figcaption>
#md # Cumulative and new expected counts a week after the cut-off for confirmed cases and deaths, the recovered total, and the isolation-bed levels and daily flows.
#md # </figcaption>
#md # </figure>
#md # ```

#md # ```@raw html
#md # <details><summary>One-week-ahead latent forecast plot</summary>
#md # ```

forecast_latent_fig = plot_forecast_latent(forecast);

#md # ```@raw html
#md # </details>
#md # <figure class="bvd-figure">
#md # ```

forecast_latent_fig #hide

#md # ```@raw html
#md # <figcaption>
#md # New infections, symptom onsets and deaths over the week after the cut-off, with the reproduction number left to evolve across it.
#md # </figcaption>
#md # </figure>
#md # ```

#md # ```@raw html
#md # <details><summary>One-week-ahead observed forecast plot</summary>
#md # ```

forecast_fig = plot_forecast(forecast);

#md # ```@raw html
#md # </details>
#md # <figure class="bvd-figure">
#md # ```

forecast_fig #hide

#md # ```@raw html
#md # <figcaption>
#md # The new count each observed stream adds over the week after the cut-off.
#md # The situation reports no longer update suspected cases and deaths, so their projections cannot be checked and the forecast evaluation leaves them out.
#md # </figcaption>
#md # </figure>
#md # ```

# The reported count of patients in isolation is the occupied-bed count, so it measures bed use limited by supply, and bed demand is the unobserved number who need a bed.
# The model carries one national bed capacity, so the national shortfall understates local unmet need (see [limitations](@ref "Limitations")).
# On 13 June Ituri was at 93.9% occupancy while Sud-Kivu was at 21.9%.

#md # ```@raw html
#md # <details><summary>One-week-ahead isolation-bed forecast plot</summary>
#md # ```

forecast_beds_fig = plot_forecast_beds(forecast);

#md # ```@raw html
#md # </details>
#md # <figure class="bvd-figure">
#md # ```

forecast_beds_fig #hide

#md # ```@raw html
#md # <figcaption>
#md # Projected isolation bed demand, the need under unconstrained supply, against the occupancy the beds can meet.
#md # The gap between them is the projected bed shortfall (right).
#md # </figcaption>
#md # </figure>
#md # ```

#md # ```@raw html
#md # <details><summary>One-week-ahead treatment-flow forecast plot</summary>
#md # ```

forecast_flows_fig = plot_forecast_flows(forecast);

#md # ```@raw html
#md # </details>
#md # <figure class="bvd-figure">
#md # ```

forecast_flows_fig #hide

#md # ```@raw html
#md # <figcaption>
#md # Daily isolation flows over the week after the cut-off: new admissions, in-care deaths and rule-outs.
#md # </figcaption>
#md # </figure>
#md # ```

# ## Symptom-onset nowcast and forecast results
#
# The [symptom-onset nowcast and forecast](@ref "Symptom-onset nowcast and forecast") Methods section describes how the onset projection is built.

#md # ```@raw html
#md # <details><summary>Generate the symptom-onset nowcast and forecast</summary>
#md # ```

onset_forecast = forecast_onsets(
    fit_forecast("joint");
    horizon = 7,
    obs_value = something(obs.onset_curve_history.last_total, 0)
);
onset_forecast_summary = onset_forecast_table(onset_forecast);

#md # ```@raw html
#md # </details>
#md # <figure class="bvd-figure">
#md # ```

onset_forecast_summary #hide

#md # ```@raw html
#md # <figcaption>
#md # The symptom-onset nowcast and forecast.
#md # The first three rows are the outbreak at the cut-off and the next three the coming week, so the two halves must not be added.
#md # </figcaption>
#md # </figure>
#md # ```

# "Onsets not yet reported at T" is not a backlog that will all arrive, because ascertainment does not reach one.
# It holds both the reporting backlog and the cases surveillance will never confirm.
# The smaller "reports this week of onsets before T" row is the part the coming week should clear.

#md # ```@raw html
#md # <details><summary>Symptom-onset nowcast and forecast plot</summary>
#md # ```

onset_forecast_fig = let
    fig = CairoMakie.Figure(; size = (960, 420))
    ## Two-line tick labels rather than rotated ones: the leftmost rotated
    ## label overhangs the axis and is clipped at the figure edge.
    ax1 = CairoMakie.Axis(
        fig[1, 1];
        title = "New onset reports over the coming week",
        ylabel = "cases", xticks = (
            1:4,
            [
                "already\nhappened", "not yet\nhappened", "sum of\nthe two",
                "as the next\nfigure reads it",
            ],
        )
    )
    ## The first three bars are latent, so the third is exactly the first
    ## two added. The fourth is that same sum replicated through the
    ## observation model, which is the scored quantity and the only one
    ## comparable to a digitised figure; it is wider by the read error,
    ## which is why the three latent bars are shown as well rather than a
    ## decomposition that appears not to add up.
    _latent_total = onset_forecast.onset_reports_backfill .+
        onset_forecast.onset_reports_future
    for (i, d, col) in (
            (
                1, onset_forecast.onset_reports_backfill,
                :mediumpurple,
            ),
            (2, onset_forecast.onset_reports_future, :mediumpurple),
            (3, _latent_total, :mediumpurple),
            (4, Float64.(onset_forecast.onset_reports_new), :slategray),
        )
        s = posterior_summary(d)
        CairoMakie.rangebars!(
            ax1, [Float64(i)], [s.lo90], [s.hi90];
            color = col, linewidth = 3
        )
        CairoMakie.rangebars!(
            ax1, [Float64(i)], [s.lo60], [s.hi60];
            color = col, linewidth = 8
        )
        CairoMakie.scatter!(
            ax1, [Float64(i)], [quantile(d, 0.5)];
            color = :black, markersize = 9
        )
    end
    ax2 = CairoMakie.Axis(
        fig[1, 2];
        title = "Symptom onsets by the cut-off",
        ylabel = "cases", xticks = (
            1:3,
            ["onsets\nto date", "reported\nby T", "not yet\nreported"],
        )
    )
    for (i, d) in enumerate(
            (
                onset_forecast.onsets_to_date,
                onset_forecast.onset_reports_to_date,
                onset_forecast.onsets_unreported,
            )
        )
        s = posterior_summary(d)
        CairoMakie.rangebars!(
            ax2, [Float64(i)], [s.lo90], [s.hi90];
            color = :seagreen, linewidth = 3
        )
        CairoMakie.rangebars!(
            ax2, [Float64(i)], [s.lo60], [s.hi60];
            color = :seagreen, linewidth = 8
        )
        CairoMakie.scatter!(
            ax2, [Float64(i)], [quantile(d, 0.5)];
            color = :black, markersize = 9
        )
    end
    ## The digitised total the "reported by T" bar is a model of, so the
    ## reader can see the fitted reported level against the figure itself.
    ismissing(obs.onset_curve_history.last_total) ||
        CairoMakie.hlines!(
        ax2,
        [Float64(obs.onset_curve_history.last_total)];
        color = :black, linestyle = :dash
    )
    fig
end;

#md # ```@raw html
#md # </details>
#md # <figure class="bvd-figure">
#md # ```

onset_forecast_fig #hide

#md # ```@raw html
#md # <figcaption>
#md # The coming week's new onset reports, split into onsets before and after the cut-off, with their sum and that sum after the observation model (left).
#md # Only the fourth bar is comparable to a digitised figure and only it is scored, and it is wider by the read error and the counting variation of new reports.
#md # The nowcast sets the onsets that have happened against the share the triangle has printed (right).
#md # </figcaption>
#md # </figure>
#md # ```
