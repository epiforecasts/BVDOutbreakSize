```@meta
EditURL = "../examples/forecast.jl"
```

# Forecasts

Every release projects each DRC stream a week ahead from the joint
posterior.
The projection is a no-change forward run, defined in the
[one-week-ahead forecast](@ref "One-week-ahead forecast") Methods section.
How these forecasts have scored against the data that arrived afterwards is
on the [evaluation](@ref "Forecast evaluation") page.

```@raw html
<details><summary>Load packages, data and fitted chains</summary>
```

````julia
# Shared setup: packages, observations, the fit registry and every model fit
# (loaded from the content-addressed cache). See `docs/examples/_setup.jl`.
using BVDOutbreakSize
include(joinpath(pkgdir(BVDOutbreakSize), "docs", "examples", "_setup.jl"))
````

````
_release_data (generic function with 1 method)
````

```@raw html
</details>
```

## One-week-ahead forecast results

The table and figures below give the cumulative and new expected counts by $T + 7$ from the no-change projection defined in the [one-week-ahead forecast](@ref "One-week-ahead forecast") Methods section.
The summary table reports the confirmed case and death streams, the recovered total and the isolation-bed levels and daily flows.
The observed-forecast plot below additionally shows the suspected case and death streams, so every projected stream appears.
The situation reports no longer update those two, so their projection cannot be checked against a later observation and the forecast validation leaves them out.

```@raw html
<details><summary>Generate the one-week-ahead forecast</summary>
```

````julia
forecast = forecast_reported(
    chn_joint;
    horizon = 7,
    obs_cases = obs.reported_cases,
    obs_deaths = obs.total_deaths,
    obs_confirmed = obs.confirmed_cases,
    obs_confirmed_deaths = obs.confirmed_deaths,
    obs_recovered = obs.recovered_cases
);
forecast_summary = forecast_table(forecast);
````

```@raw html
</details>
```

```@raw html
<details><summary>One-week-ahead forecast summary table</summary>
```


```@raw html
<div><div style = "float: left;"><span>11×8 DataFrame</span></div><div style = "clear: both;"></div></div><div class = "data-frame" style = "overflow-x: scroll;"><table class = "data-frame" style = "margin-bottom: 6px;"><thead><tr class = "columnLabelRow"><th class = "stubheadLabel" style = "font-weight: bold; text-align: right;">Row</th><th style = "text-align: left;">Stream</th><th style = "text-align: left;">Quantity</th><th style = "text-align: left;">Lower 90%</th><th style = "text-align: left;">Lower 60%</th><th style = "text-align: left;">Lower 30%</th><th style = "text-align: left;">Upper 30%</th><th style = "text-align: left;">Upper 60%</th><th style = "text-align: left;">Upper 90%</th></tr><tr class = "columnLabelRow"><th class = "stubheadLabel" style = "font-weight: bold; text-align: right;"></th><th title = "String" style = "text-align: left;">String</th><th title = "String" style = "text-align: left;">String</th><th title = "Float64" style = "text-align: left;">Float64</th><th title = "Float64" style = "text-align: left;">Float64</th><th title = "Float64" style = "text-align: left;">Float64</th><th title = "Float64" style = "text-align: left;">Float64</th><th title = "Float64" style = "text-align: left;">Float64</th><th title = "Float64" style = "text-align: left;">Float64</th></tr></thead><tbody><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">1</td><td style = "text-align: left;">DRC confirmed cases</td><td style = "text-align: left;">cumulative by T+7</td><td style = "text-align: right;">7959.0</td><td style = "text-align: right;">8069.0</td><td style = "text-align: right;">8138.0</td><td style = "text-align: right;">8277.0</td><td style = "text-align: right;">8377.0</td><td style = "text-align: right;">8555.0</td></tr><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">2</td><td style = "text-align: left;">DRC confirmed cases</td><td style = "text-align: left;">new this week</td><td style = "text-align: right;">484.0</td><td style = "text-align: right;">594.0</td><td style = "text-align: right;">663.0</td><td style = "text-align: right;">802.0</td><td style = "text-align: right;">902.0</td><td style = "text-align: right;">1080.0</td></tr><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">3</td><td style = "text-align: left;">DRC confirmed deaths</td><td style = "text-align: left;">cumulative by T+7</td><td style = "text-align: right;">3831.0</td><td style = "text-align: right;">3866.0</td><td style = "text-align: right;">3890.0</td><td style = "text-align: right;">3933.0</td><td style = "text-align: right;">3963.0</td><td style = "text-align: right;">4028.0</td></tr><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">4</td><td style = "text-align: left;">DRC confirmed deaths</td><td style = "text-align: left;">new this week</td><td style = "text-align: right;">226.0</td><td style = "text-align: right;">261.0</td><td style = "text-align: right;">285.0</td><td style = "text-align: right;">328.0</td><td style = "text-align: right;">358.0</td><td style = "text-align: right;">423.0</td></tr><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">5</td><td style = "text-align: left;">DRC isolation beds</td><td style = "text-align: left;">demand at T+7</td><td style = "text-align: right;">839.0</td><td style = "text-align: right;">994.0</td><td style = "text-align: right;">1092.0</td><td style = "text-align: right;">1270.0</td><td style = "text-align: right;">1402.0</td><td style = "text-align: right;">1743.0</td></tr><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">6</td><td style = "text-align: left;">DRC isolation beds</td><td style = "text-align: left;">occupancy at T+7</td><td style = "text-align: right;">721.0</td><td style = "text-align: right;">867.0</td><td style = "text-align: right;">967.0</td><td style = "text-align: right;">1144.0</td><td style = "text-align: right;">1277.0</td><td style = "text-align: right;">1411.0</td></tr><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">7</td><td style = "text-align: left;">DRC isolation admissions</td><td style = "text-align: left;">daily at T+7</td><td style = "text-align: right;">125.0</td><td style = "text-align: right;">148.0</td><td style = "text-align: right;">166.0</td><td style = "text-align: right;">199.0</td><td style = "text-align: right;">224.0</td><td style = "text-align: right;">288.0</td></tr><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">8</td><td style = "text-align: left;">DRC in-care deaths</td><td style = "text-align: left;">daily at T+7</td><td style = "text-align: right;">16.0</td><td style = "text-align: right;">21.0</td><td style = "text-align: right;">25.0</td><td style = "text-align: right;">32.0</td><td style = "text-align: right;">37.0</td><td style = "text-align: right;">49.0</td></tr><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">9</td><td style = "text-align: left;">DRC isolation rule-outs</td><td style = "text-align: left;">daily at T+7</td><td style = "text-align: right;">73.0</td><td style = "text-align: right;">88.0</td><td style = "text-align: right;">97.0</td><td style = "text-align: right;">114.0</td><td style = "text-align: right;">126.0</td><td style = "text-align: right;">160.0</td></tr><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">10</td><td style = "text-align: left;">DRC recovered among confirmed</td><td style = "text-align: left;">cumulative by T+7</td><td style = "text-align: right;">1923.0</td><td style = "text-align: right;">1947.0</td><td style = "text-align: right;">1963.0</td><td style = "text-align: right;">1990.0</td><td style = "text-align: right;">2011.0</td><td style = "text-align: right;">2050.0</td></tr><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">11</td><td style = "text-align: left;">DRC recovered among confirmed</td><td style = "text-align: left;">new this week</td><td style = "text-align: right;">125.0</td><td style = "text-align: right;">149.0</td><td style = "text-align: right;">165.0</td><td style = "text-align: right;">192.0</td><td style = "text-align: right;">213.0</td><td style = "text-align: right;">252.0</td></tr></tbody></table></div>
```

```@raw html
</details>
```

The latent figure shows the new infections, symptom onsets and deaths over the horizon, with the reproduction number left to keep evolving across it.

```@raw html
<details><summary>One-week-ahead latent forecast plot</summary>
```

````julia
forecast_latent_fig = plot_forecast_latent(forecast);
````

```@raw html
</details>
```

![](forecast-17.png)

The observed figure shows the new count each observed stream adds over the horizon: suspected cases, suspected deaths, laboratory-confirmed cases, confirmed deaths and recovered, one panel per stream the forecast carries.

```@raw html
<details><summary>One-week-ahead observed forecast plot</summary>
```

````julia
forecast_fig = plot_forecast(forecast);
````

```@raw html
</details>
```

![](forecast-22.png)

The bed figure shows the projected isolation/treatment-bed demand (the need a week ahead, under unconstrained supply) against the supply-limited occupancy the beds can actually meet.
The gap between the two is the projected bed shortfall, shown in the right panel.
The reported "Patients en isolement" count is the occupied-bed count (the report computes the "Taux d'occupation" as that count over the bed capacity), so isolation is bed usage, gated by supply.
The demand is its unobserved counterpart, the number who need a bed.
The model carries a single national bed capacity, so it cannot represent local saturation, and the national shortfall understates local unmet need.
On 13 June Ituri was at 93.9% occupancy while Sud-Kivu was at 21.9%; beds free in one province cannot serve patients in another.

```@raw html
<details><summary>One-week-ahead isolation-bed forecast plot</summary>
```

````julia
forecast_beds_fig = plot_forecast_beds(forecast);
````

```@raw html
</details>
```

![](forecast-27.png)

The flow figure projects the daily isolation/treatment flows a week ahead: new admissions, in-care deaths and rule-outs, each grown from its cut-off daily rate and replicated through the isolation dispersion.

```@raw html
<details><summary>One-week-ahead treatment-flow forecast plot</summary>
```

````julia
forecast_flows_fig = plot_forecast_flows(forecast);
````

```@raw html
</details>
```

![](forecast-32.png)

The forecast split by province is below, for the two streams the spatial tables report.
Each province's count is the national draw times that province's modelled share at the most recent spatial vintage, so the split is held at its current value over the week.

```@raw html
<details><summary>Province forecast split</summary>
```

````julia
province_forecast_fig = plot_province_forecast(
    chn_joint, forecast;
    n_patches = N_PATCHES
);
province_forecast = province_forecast_table(
    chn_joint, forecast;
    n_patches = N_PATCHES
);
````

```@raw html
</details>
```

![](forecast-37.png)

## Symptom-onset nowcast and forecast results

The table below gives the onset stream's projection, built as described in the [symptom-onset nowcast and forecast](@ref "Symptom-onset nowcast and forecast") Methods section.
The two halves must not be added together: the first three rows are the state of the outbreak at the cut-off, the next three the coming week.

"Onsets not yet reported at T" is not a backlog that will all arrive, because ascertainment does not reach one.
The row holds two things together: the reporting backlog, and the cases surveillance will never confirm.
The "reports this week of onsets before T" row is the part of it the coming week should actually clear.
It is the smaller number.

```@raw html
<details><summary>Generate the symptom-onset nowcast and forecast</summary>
```

````julia
# `_onset_grid_start`/`_onset_grid_end` are the triangle's own grid, built
# from the observations in the shared setup.
onset_forecast = forecast_onsets(
    chn_joint;
    grid_start = _onset_grid_start, grid_end = _onset_grid_end,
    n = obs.n, horizon = 7,
    obs_value = something(obs.onset_curve_history.last_total, 0)
);
onset_forecast_summary = onset_forecast_table(onset_forecast);
````

```@raw html
</details>
```

```@raw html
<details><summary>Symptom-onset nowcast and forecast summary table</summary>
```


```@raw html
<div><div style = "float: left;"><span>7×7 DataFrame</span></div><div style = "clear: both;"></div></div><div class = "data-frame" style = "overflow-x: scroll;"><table class = "data-frame" style = "margin-bottom: 6px;"><thead><tr class = "columnLabelRow"><th class = "stubheadLabel" style = "font-weight: bold; text-align: right;">Row</th><th style = "text-align: left;">Quantity</th><th style = "text-align: left;">Lower 90%</th><th style = "text-align: left;">Lower 60%</th><th style = "text-align: left;">Lower 30%</th><th style = "text-align: left;">Upper 30%</th><th style = "text-align: left;">Upper 60%</th><th style = "text-align: left;">Upper 90%</th></tr><tr class = "columnLabelRow"><th class = "stubheadLabel" style = "font-weight: bold; text-align: right;"></th><th title = "String" style = "text-align: left;">String</th><th title = "Float64" style = "text-align: left;">Float64</th><th title = "Float64" style = "text-align: left;">Float64</th><th title = "Float64" style = "text-align: left;">Float64</th><th title = "Float64" style = "text-align: left;">Float64</th><th title = "Float64" style = "text-align: left;">Float64</th><th title = "Float64" style = "text-align: left;">Float64</th></tr></thead><tbody><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">1</td><td style = "text-align: left;">symptom onsets to date</td><td style = "text-align: right;">9427.0</td><td style = "text-align: right;">10636.0</td><td style = "text-align: right;">11597.0</td><td style = "text-align: right;">13635.0</td><td style = "text-align: right;">15139.0</td><td style = "text-align: right;">18534.0</td></tr><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">2</td><td style = "text-align: left;">of those, reported by T</td><td style = "text-align: right;">5341.0</td><td style = "text-align: right;">5581.0</td><td style = "text-align: right;">5735.0</td><td style = "text-align: right;">6019.0</td><td style = "text-align: right;">6197.0</td><td style = "text-align: right;">6517.0</td></tr><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">3</td><td style = "text-align: left;">onsets not yet reported at T</td><td style = "text-align: right;">3677.0</td><td style = "text-align: right;">4811.0</td><td style = "text-align: right;">5707.0</td><td style = "text-align: right;">7727.0</td><td style = "text-align: right;">9126.0</td><td style = "text-align: right;">12551.0</td></tr><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">4</td><td style = "text-align: left;">reports this week of onsets before T</td><td style = "text-align: right;">305.0</td><td style = "text-align: right;">357.0</td><td style = "text-align: right;">391.0</td><td style = "text-align: right;">455.0</td><td style = "text-align: right;">499.0</td><td style = "text-align: right;">581.0</td></tr><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">5</td><td style = "text-align: left;">reports this week of onsets after T</td><td style = "text-align: right;">100.0</td><td style = "text-align: right;">134.0</td><td style = "text-align: right;">158.0</td><td style = "text-align: right;">208.0</td><td style = "text-align: right;">244.0</td><td style = "text-align: right;">352.0</td></tr><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">6</td><td style = "text-align: left;">new onset reports this week</td><td style = "text-align: right;">270.0</td><td style = "text-align: right;">434.0</td><td style = "text-align: right;">535.0</td><td style = "text-align: right;">704.0</td><td style = "text-align: right;">806.0</td><td style = "text-align: right;">1021.0</td></tr><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">7</td><td style = "text-align: left;">new symptom onsets this week</td><td style = "text-align: right;">662.0</td><td style = "text-align: right;">856.0</td><td style = "text-align: right;">1002.0</td><td style = "text-align: right;">1314.0</td><td style = "text-align: right;">1561.0</td><td style = "text-align: right;">2221.0</td></tr></tbody></table></div>
```

```@raw html
</details>
```

The left panel splits the coming week's new onset reports into reports of onsets that had already happened by the cut-off and reports of onsets still to come, and shows their sum.
The fourth bar is the same sum after it has been through the observation model, which is what the next vintage will actually print.
It is much wider than the sum it replicates, and the gap is the rescan.
At a printed total of a couple of thousand cases the per-scan level error alone is worth tens of cases either way, well above the epidemic uncertainty on a week of new reports.
Only the fourth bar is comparable to a digitised figure, and only it is scored.

The right panel puts the nowcast itself on the same axes, the onsets that have happened against the share of them the triangle has printed.

```@raw html
<details><summary>Symptom-onset nowcast and forecast plot</summary>
```

````julia
onset_forecast_fig = let
    fig = CairoMakie.Figure(; size = (960, 420))
    # Two-line tick labels rather than rotated ones: the leftmost rotated
    # label overhangs the axis and is clipped at the figure edge.
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
    # The first three bars are latent, so the third is exactly the first
    # two added. The fourth is that same sum replicated through the
    # observation model, which is the scored quantity and the only one
    # comparable to a digitised figure; it is wider by the per-scan level
    # error, which is why the three latent bars are shown as well rather
    # than a decomposition that appears not to add up.
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
    # The digitised total the "reported by T" bar is a model of, so the
    # reader can see the fitted reported level against the figure itself.
    ismissing(obs.onset_curve_history.last_total) ||
        CairoMakie.hlines!(
        ax2,
        [Float64(obs.onset_curve_history.last_total)];
        color = :black, linestyle = :dash
    )
    fig
end;
````

```@raw html
</details>
```

![](forecast-49.png)

