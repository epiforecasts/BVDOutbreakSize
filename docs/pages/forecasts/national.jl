# # Forecasts
#
# Every release projects each DRC stream a week ahead from the joint
# posterior.
# The projection is a no-change forward run, defined in the
# [one-week-ahead forecast](@ref "One-week-ahead forecast") Methods section.
# How these forecasts have scored against the data that arrived afterwards is
# on the [evaluation](@ref "Forecast evaluation") page.
# The split by province is on the [province forecasts](@ref "Province forecasts") page.

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
# The expected counts for the week after the cut-off, from the forecast below.

#md # ```@raw html
#md # <details><summary>Generate the one-week-ahead forecast</summary>
#md # ```

forecast = forecast_reported(
    chn_joint;
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
#
# The table and figures below give the cumulative and new expected counts by $T + 7$ from the no-change projection defined in the [one-week-ahead forecast](@ref "One-week-ahead forecast") Methods section.
# The summary table reports the confirmed case and death streams, the recovered total and the isolation-bed levels and daily flows.
# The observed-forecast plot below additionally shows the suspected case and death streams, so every projected stream appears.
# The situation reports no longer update those two, so their projection cannot be checked against a later observation and the forecast validation leaves them out.

#md # ```@raw html
#md # <details><summary>Summarise the one-week-ahead forecast</summary>
#md # ```

forecast_summary = forecast_table(forecast);

#md # ```@raw html
#md # </details>
#md # ```

#md # ```@raw html
#md # <details><summary>One-week-ahead forecast summary table</summary>
#md # ```

forecast_summary #hide

#md # ```@raw html
#md # </details>
#md # ```

# The latent figure shows the new infections, symptom onsets and deaths over the horizon, with the reproduction number left to keep evolving across it.

#md # ```@raw html
#md # <details><summary>One-week-ahead latent forecast plot</summary>
#md # ```

forecast_latent_fig = plot_forecast_latent(forecast);

#md # ```@raw html
#md # </details>
#md # ```

forecast_latent_fig #hide

# The observed figure shows the new count each observed stream adds over the horizon: suspected cases, suspected deaths, laboratory-confirmed cases, confirmed deaths and recovered, one panel per stream the forecast carries.

#md # ```@raw html
#md # <details><summary>One-week-ahead observed forecast plot</summary>
#md # ```

forecast_fig = plot_forecast(forecast);

#md # ```@raw html
#md # </details>
#md # ```

forecast_fig #hide

# The bed figure shows the projected isolation/treatment-bed demand (the need a week ahead, under unconstrained supply) against the supply-limited occupancy the beds can actually meet.
# The gap between the two is the projected bed shortfall, shown in the right panel.
# The reported "Patients en isolement" count is the occupied-bed count (the report computes the "Taux d'occupation" as that count over the bed capacity), so isolation is bed usage, gated by supply.
# The demand is its unobserved counterpart, the number who need a bed.
# The model carries a single national bed capacity, so it cannot represent local saturation, and the national shortfall understates local unmet need.
# On 13 June Ituri was at 93.9% occupancy while Sud-Kivu was at 21.9%; beds free in one province cannot serve patients in another.

#md # ```@raw html
#md # <details><summary>One-week-ahead isolation-bed forecast plot</summary>
#md # ```

forecast_beds_fig = plot_forecast_beds(forecast);

#md # ```@raw html
#md # </details>
#md # ```

forecast_beds_fig #hide

# The flow figure projects the daily isolation/treatment flows a week ahead: new admissions, in-care deaths and rule-outs, each grown from its cut-off daily rate and replicated through the isolation dispersion.

#md # ```@raw html
#md # <details><summary>One-week-ahead treatment-flow forecast plot</summary>
#md # ```

forecast_flows_fig = plot_forecast_flows(forecast);

#md # ```@raw html
#md # </details>
#md # ```

forecast_flows_fig #hide

# ## Symptom-onset nowcast and forecast results
#
# The table below gives the onset stream's projection, built as described in the [symptom-onset nowcast and forecast](@ref "Symptom-onset nowcast and forecast") Methods section.
# The two halves must not be added together: the first three rows are the state of the outbreak at the cut-off, the next three the coming week.
#
# "Onsets not yet reported at T" is not a backlog that will all arrive, because ascertainment does not reach one.
# The row holds two things together: the reporting backlog, and the cases surveillance will never confirm.
# The "reports this week of onsets before T" row is the part of it the coming week should actually clear.
# It is the smaller number.

#md # ```@raw html
#md # <details><summary>Generate the symptom-onset nowcast and forecast</summary>
#md # ```

## `_onset_grid_start`/`_onset_grid_end` are the triangle's own grid, built
## from the observations in the shared setup.
onset_forecast = forecast_onsets(
    chn_joint;
    grid_start = _onset_grid_start, grid_end = _onset_grid_end,
    n = obs.n, horizon = 7,
    obs_value = something(obs.onset_curve_history.last_total, 0)
);
onset_forecast_summary = onset_forecast_table(onset_forecast);

#md # ```@raw html
#md # </details>
#md # ```

#md # ```@raw html
#md # <details><summary>Symptom-onset nowcast and forecast summary table</summary>
#md # ```

onset_forecast_summary #hide

#md # ```@raw html
#md # </details>
#md # ```

# The left panel splits the coming week's new onset reports into reports of onsets that had already happened by the cut-off and reports of onsets still to come, and shows their sum.
# The fourth bar is the same sum after it has been through the observation model, which is what the next vintage will actually print.
# It is wider than the sum it replicates by the two reads' error and the counting variation of the new reports.
# Only the fourth bar is comparable to a digitised figure, and only it is scored.
#
# The right panel puts the nowcast itself on the same axes, the onsets that have happened against the share of them the triangle has printed.

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
#md # ```

onset_forecast_fig #hide
