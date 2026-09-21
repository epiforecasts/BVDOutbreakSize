```@meta
EditURL = "../examples/evaluation.jl"
```

# Forecast evaluation

How the forecasts on the [forecasts](@ref "Forecasts") page have scored
against the data that arrived afterwards.
Scoring is the continuous ranked probability score against a persistence
baseline, defined in the
[forecast scoring](@ref "Forecast scoring against a persistence baseline")
Methods section.

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

## Forecast validation

How last week's forecast held up against the data since observed, using the frozen re-fit and one-week projection defined in [forecast-versus-frozen evaluation](@ref "Forecast-versus-frozen evaluation").
Only the streams the situation reports are still updating are validated here.
A stream that has stopped being reported carries a cumulative total that repeats its last reported value, so there is no observation for the past week to score against.
The frozen fit also conditions on the isolation beds, so the projected bed occupancy is scored against the beds held a week later.
The bed validation is weak at a one-week-back freeze.
The reported occupancy rate starts only on 9 June, so the capacity has no implied-capacity anchor and rides its random walk back to the freeze date.
Like the scores further down, the confirmed new-count rows here take out any retrospective harmonisation step the week contained.
Such a step reattaches records notified earlier, so it is not something the forecast was predicting.
The cumulative rows are scored against the published total, harmonisation included.

```@raw html
<details><summary>Fit one week back and validate the one-week-ahead forecast</summary>
```

````julia
# frozen_lastweek and frozen_lastweek_streams are computed in the setup
# block above.
# `obs_recovered` is passed so the frozen fit's forecast carries a
# `recovered_new` column (materialised only when the recovered origin is
# given), letting the recovered stream be scored against the observed count
# below like the other streams.
# The onset grid is the one the FROZEN fit saw, not the live one, so the
# validation forecast carries an `onset reports` row scored on the triangle
# the frozen fit was actually fitted to.
_val_onset_days = frozen_lastweek.o.onset_curve_history.onset_days
_val_grid_start = isempty(_val_onset_days) ? nothing :
    minimum(_val_onset_days)
_val_grid_end = isnothing(_val_grid_start) ? nothing :
    max(
        maximum(frozen_lastweek.o.onset_curve_history.report_days),
        _val_grid_start
    )
validation_forecast = forecast_reported(
    frozen_lastweek.chn;
    horizon = 7,
    obs_cases = frozen_lastweek.o.reported_cases,
    obs_deaths = frozen_lastweek.o.total_deaths,
    obs_confirmed = frozen_lastweek.o.confirmed_cases,
    obs_confirmed_deaths = frozen_lastweek.o.confirmed_deaths,
    obs_recovered = frozen_lastweek.o.recovered_cases,
    grid_n = frozen_lastweek.o.n,
    onset_grid_start = _val_grid_start, onset_grid_end = _val_grid_end
);

# Each frozen individual (single-stream) fit's own one-week-ahead new-count
# forecast at the same cut-off as `frozen_lastweek`, from
# [`forecast_stream`](@ref) (the same per-stream forecaster
# `stream_forecasts.csv` uses), so the validation plots below can show the
# individual fit alongside the joint rather than the joint alone. Recovered
# has no individual fit and is absent here, as it is throughout this report.
# Only the still-reported streams are fitted at the validation cut-off, so
# a stream the situation reports have stopped updating is absent from
# `frozen_lastweek_streams` and carries no individual series here.
function _validation_individual_new(sid, stream::Symbol, obs_field)
    haskey(frozen_lastweek_streams, sid) || return nothing
    f = frozen_lastweek_streams[sid]
    bp = f.o.n - f.o.who_first_sitrep_days
    return Float64.(
        forecast_stream(
            f.chn, stream; horizon = 7,
            obs_value = getproperty(f.o, obs_field), n = f.o.n, breakpoint = bp,
            rt_start = 1, rt_walk_start = 1
        )
    )
end
validation_individual = NamedTuple(
    k => v
        for (k, v) in pairs(
            (;
                cases_new = _validation_individual_new(
                    "cases", :reported_cases, :reported_cases
                ),
                deaths_new = _validation_individual_new(
                    "deaths", :suspected_deaths, :total_deaths
                ),
                confirmed_new = _validation_individual_new(
                    "confirmed", :confirmed_cases, :confirmed_cases
                ),
                confirmed_deaths_new = _validation_individual_new(
                    "confirmed_deaths", :confirmed_deaths, :confirmed_deaths
                ),
            )
        )
        if !isnothing(v)
)
# The frozen individual (treatment-only) fit's own bed-occupancy forecast,
# anchored on the beds occupied at ITS OWN cut-off (the frozen fit's own
# `o`, not the current `obs`), matching how the joint frozen forecast is
# itself anchored.
# `nothing` when the beds have stopped being reported, so the treatment fit
# is absent; the bed panel then draws the joint alone.
# A `let` block, not a bare `if`: a top-level `if` shares the script's
# global scope, so its working names would leak into the rest of the page.
validation_individual_isolation = let
    if haskey(frozen_lastweek_streams, "treatment")
        tf = frozen_lastweek_streams["treatment"]
        beds = isempty(tf.o.isolation_history.counts) ? 0.0 :
            Float64(tf.o.isolation_history.counts[end])
        Float64.(
            forecast_stream(
                tf.chn, :isolation_beds; horizon = 7,
                obs_value = beds, n = tf.o.n,
                breakpoint = tf.o.n - tf.o.who_first_sitrep_days,
                rt_start = 1, rt_walk_start = 1
            )
        )
    else
        nothing
    end
end

# The observed beds at the current cut-off (the forecast target), so the
# frozen-fit bed forecast is scored against what the beds actually held.
# Held back once the beds stop being reported, since the last count would
# then be carried forward rather than observed at the target date.
_obs_beds = stream_reporting(obs, :isolation_beds) ?
    obs.isolation_history.counts[end] : missing
# Same observed/baseline keying as the plot below, so the table covers every
# fitted count stream (cumulative and new-count rows) plus the bed level.
# A harmonisation-break day between the frozen cut-off and the current one
# puts records into the confirmed cumulative that were never notified in that
# week, so the new-count truth carries a step the forecast was never
# predicting. Take it out, the same correction `score_releases.jl` applies.
# Grid days are relative to a seeding date fixed by the genetic tmrca, so the
# frozen fit's own `n` and the current `obs.n` index the same grid.
validation_breaks = (
    confirmed_cum = confirmed_break_correction(
        obs, frozen_lastweek.o.n, obs.n
    ),
    confirmed_deaths_cum = confirmed_break_correction(
        obs, frozen_lastweek.o.n, obs.n; deaths = true
    ),
)

# Observed cumulative at the target date per stream, keyed by the forecast's
# cumulative column; `baseline` is each stream's origin cumulative (the
# frozen cut-off), so the new count is scored against observed minus origin,
# less any harmonisation the window carries (see `validation_breaks`). Both
# the table and the plot below take the still-reported streams
# (`reporting_cum_cols`, from the setup block): a stream the situation
# reports have stopped updating has an origin and a target reading the same
# repeated total, so its cumulative truth is stale and its new-count truth is
# a guaranteed zero.
validation_observed = (
    cases_cum = obs.reported_cases,
    deaths_cum = obs.total_deaths,
    confirmed_cum = obs.confirmed_cases,
    confirmed_deaths_cum = obs.confirmed_deaths,
    recovered_cum = obs.recovered_cases,
)
validation_baseline = (
    cases_cum = frozen_lastweek.o.reported_cases,
    deaths_cum = frozen_lastweek.o.total_deaths,
    confirmed_cum = frozen_lastweek.o.confirmed_cases,
    confirmed_deaths_cum = frozen_lastweek.o.confirmed_deaths,
    recovered_cum = frozen_lastweek.o.recovered_cases,
)

validation_table = forecast_vs_truth(
    validation_forecast;
    observed = keep_streams(validation_observed, reporting_cum_cols),
    baseline = keep_streams(validation_baseline, reporting_cum_cols),
    breaks = validation_breaks,
    isolation = _obs_beds
);
````

```@raw html
</details>
```

```@raw html
<details><summary>Forecast-versus-observed validation table</summary>
```


| Stream | Quantity | Observed | Lower 90% | Lower 60% | Lower 30% | Upper 30% | Upper 60% | Upper 90% | Within 90% PI |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| DRC confirmed cases | cumulative by T+7 | 7475 | 7420 | 7509 | 7584 | 7728 | 7813 | 8020 | yes |
| DRC confirmed cases | new this week | 533 | 478 | 567 | 642 | 786 | 871 | 1078 | yes |
| DRC confirmed deaths | cumulative by T+7 | 3605 | 3571 | 3608 | 3634 | 3677 | 3711 | 3772 | yes |
| DRC confirmed deaths | new this week | 256 | 222 | 259 | 285 | 328 | 362 | 423 | yes |
| DRC recovered among confirmed | cumulative by T+7 | 1798 | 1812 | 1953 | 2023 | 2126 | 2191 | 2399 | no |
| DRC recovered among confirmed | new this week | 151 | 165 | 306 | 376 | 479 | 544 | 752 | no |
| DRC isolation beds | occupancy at T+7 | 905 | 650 | 803 | 885 | 1048 | 1183 | 1403 | yes |


```@raw html
</details>
```

The observation panels histogram the one-week-ahead forecast made from the frozen fit: a cumulative and a new-count panel for each still-reported count stream the forecast carries.
The 90% predictive interval is shaded, and the count observed by the current cut-off is a dashed black rule.
Where a stream has its own individual (single-stream) fit, that fit's forecast from the same frozen cut-off is overlaid as a dotted step outline on the joint's own histogram bins.

```@raw html
<details><summary>Forecast-versus-observed plot</summary>
```

````julia
validation_fig = plot_forecast_vs_truth(
    validation_forecast;
    observed = keep_streams(validation_observed, reporting_cum_cols),
    baseline = keep_streams(validation_baseline, reporting_cum_cols),
    breaks = validation_breaks,
    individual = keep_streams(validation_individual, reporting_cum_cols)
);
````

```@raw html
</details>
```

![](evaluation-17.png)

The bed panel scores last week's projected occupancy against the beds occupied now (the dashed rule), with the individual (treatment-only) fit's own projection overlaid as a dotted step outline on the joint's own histogram bins.

```@raw html
<details><summary>Bed forecast-versus-observed plot</summary>
```

````julia
validation_beds_fig = plot_forecast_beds_vs_truth(
    validation_forecast;
    isolation = _obs_beds, individual = validation_individual_isolation
);
````

```@raw html
</details>
```

![](evaluation-22.png)

The latent quantities are not observed, so they are scored distribution against distribution: what the frozen fit forecast for the past week's new infections, onsets and deaths against what the current fit now estimates for the same window.

```@raw html
<details><summary>Forecast-versus-now latent plot</summary>
```

````julia
# Current fit's draws of the new latent counts over the past week, the last
# seven days of each cumulative-trajectory deterministic.
function _now_new(chn, key)
    mat = chn[key]
    trajs = [collect(v) for v in vec(collect(mat))]
    return Float64[t[end] - t[max(1, length(t) - 7)] for t in trajs]
end
now_latent = (;
    infections_new = _now_new(chn_joint, :cumulative_infections),
    onsets_new = _now_new(chn_joint, :cumulative_onsets),
    deaths_latent_new = _now_new(chn_joint, :cumulative_expected_deaths),
)

validation_latent_fig = plot_forecast_vs_truth_latent(
    validation_forecast; now = now_latent
);
````

```@raw html
</details>
```

![](evaluation-27.png)

### Forecast by province

The one-week-ahead forecast split by province, scored against what each province went on to report.
Each province's forecast is the national draw times its modelled share at the frozen fit's most recent spatial vintage, multiplied draw by draw so the interval carries the correlation between the two rather than treating a province's share as independent of the national total.
The share is held over the horizon, which is the assumption the width does not express: a province whose share is moving is scored as though it were not.
Every release's archived split is scored against what has since been observed in [Forecast by province across releases](@ref "Forecast by province across releases").

```@raw html
<details><summary>Province forecast against observed</summary>
```

````julia
# Per-province cumulative confirmed cases and deaths at the frozen cut-off
# and at the current one, so the truth for the week is their difference.
# Read off the same increment matrices the compositions are scored on, so
# the clamped revision is treated identically on both sides.
province_truth = let
    cur_c = province_increment_matrix(
        obs.province_confirmed_history,
        PROVINCE_NAMES, N_PATCHES
    )
    cur_d = province_increment_matrix(
        obs.province_death_history,
        PROVINCE_NAMES, N_PATCHES
    )
    froz_c = province_increment_matrix(
        frozen_lastweek.o.province_confirmed_history,
        PROVINCE_NAMES, N_PATCHES
    )
    froz_d = province_increment_matrix(
        frozen_lastweek.o.province_death_history, PROVINCE_NAMES, N_PATCHES
    )
    (;
        observed = vec(sum(cur_c.increments; dims = 2)),
        baseline = vec(sum(froz_c.increments; dims = 2)),
        death_observed = vec(sum(cur_d.increments; dims = 2)),
        death_baseline = vec(sum(froz_d.increments; dims = 2)),
    )
end

province_validation_table = province_forecast_vs_truth(
    frozen_lastweek.chn, validation_forecast;
    observed = province_truth.observed,
    baseline = province_truth.baseline,
    death_observed = province_truth.death_observed,
    death_baseline = province_truth.death_baseline,
    n_patches = N_PATCHES
);
````

```@raw html
</details>
```


| Province | Stream | Lower 90% | Upper 90% | Observed | Within 90% PI |
| --- | --- | ---: | ---: | ---: | --- |
| Ituri | Confirmed cases | 265 | 607 | 284 | true |
| Ituri | Confirmed deaths | 120 | 235 | 141 | true |
| Nord-Kivu | Confirmed cases | 170 | 434 | 211 | true |
| Nord-Kivu | Confirmed deaths | 82 | 167 | 103 | true |
| Haut-Uele | Confirmed cases | 22 | 61 | 28 | true |
| Haut-Uele | Confirmed deaths | 12 | 28 | 8 | false |
| Other provinces | Confirmed cases | 3 | 10 | 10 | false |
| Other provinces | Confirmed deaths | 1 | 4 | 4 | false |


### Streams no longer reported

The situation reports have stopped updating some of the streams the model fits, listed with the date each was last reported below.
The panels show what the frozen fit projected for those streams over the same week, without an observed rule, since the count they would be scored against has not moved since the stream stopped.

```@raw html
<details><summary>Forecast for the streams no longer reported</summary>
```

````julia
# The last-reported date per stopped stream, and the frozen fit's own
# projection for them. `plot_forecast` draws a panel per new-count column
# the frame carries, so passing the stopped streams' columns alone gives the
# projection without the fabricated truth rule the validation figure would
# otherwise draw against a repeated total.
validation_stopped_streams = let s = stream_report_status(obs),
        ids = [stream_id(c) for c in stopped_cum_cols]

    keep = [r.stream in ids for r in eachrow(s)]
    DataFrame(
        "Stream" => s[keep, :label],
        "Last reported" => s[keep, :last_date]
    )
end
_stopped_new_cols = [
    c
        for c in new_cols(stopped_cum_cols)
        if c in propertynames(validation_forecast)
]
validation_stopped_fig = plot_forecast(
    validation_forecast[!, _stopped_new_cols]
);
````

```@raw html
</details>
```

![](evaluation-37.png)

## Forecast scoring across releases

Every release's saved one- to four-week-ahead forecast is scored against the data observed since, against a persistence baseline and, where one exists, the stream's own individual fit as well as the joint.
The tables in this section are the joint model's, one row per stream.
See [forecast scoring against a persistence baseline](@ref "Forecast scoring against a persistence baseline") for how the scores, the relative skill and the baseline are built.
Recovered has no individual fit of its own, so its comparison is the baseline against the joint only.
Reported cases and suspected deaths stopped being updated by the situation reports partway through the outbreak, and exports' confirmed-detection series is anchored to an earlier cut-off.
Exports therefore contributes no scored forecast, and reported cases and suspected deaths each rest on exactly one matched forecast, a single window rather than a settled sample.

Only a minority of the daily releases examined contribute a row to the table below, each a reconstruction of an earlier model version rather than the current fit.
Only the newest few releases carry the current model's own individual-stream forecasts, and the backfilled reconstructions carry none at all.
Every row also rests on one to a handful of matched forecasts, shown as its own count rather than rounded away.

Five things are excluded from the scores here and in the frozen section below, each for a stated reason rather than for scoring badly, and the second of them applies to the frozen section alone.

- One whole reconstruction (`results-v1.6.0`): its chain forecasts a near-zero median at every horizon and stream, with the upper predictive tail occasionally reaching five- and six-digit values, which is the signature of a chain that failed to sample rather than a forecast.
- Frozen section only: the confirmed-death rows of the fourteen frozen reconstructions cut between 16 July and 15 August 2026, whose forecaster could not project that stream from its own trajectory and floored it at zero, so each carries a one-week median of exactly zero against an observed 250 to 370. Reconstructions cut after that window project the stream normally.
- An onset window containing a vintage whose reread total falls, since its increment is not what the situation reports added.
- A stream that carries no persistence baseline, which is what makes a window scoreable at all.
- A province window holding a harmonisation-break day, since that day's backfill is published for the country and not by province.

Nothing is dropped from the archive itself; `data/forecast_scores*.csv` and `data/forecast_overlay*.csv` record everything that was scored.

The symptom-onset stream is scored on the new reported count each vintage adds rather than on its level, because every vintage rereads the whole figure.
Its printed total therefore moves with the scan error as well as with late reporting.
On fourteen vintages the reread total falls, which a cumulative onset curve cannot do, and on many others it repeats unchanged.
The fit absorbs that with a per-vintage scan level; the scored truth cannot, since it is the increment between the vintages at the two ends of a window.
A window containing a falling vintage is therefore left unscored, the rule province windows holding a harmonisation-break day already follow.
It bites hardest at the longer horizons, a four-week window being more likely to contain a reread than a one-week one: the frozen onset row keeps three of its twenty-nine windows, all at one week, and the cross-release row six of thirty-eight.
Read the onset row's skill against the baseline rather than its coverage, and read it as resting on a handful of windows.

```@raw html
<details><summary>Load and summarise the cross-release forecast scores</summary>
```

````julia
forecast_scores_df = _release_data(
    "forecast_scores.csv",
    (;
        release = String, made_date = Date, stream = String, horizon = Int,
        target_date = Date, fit = String, crps = Float64,
        log_crps = Float64, dispersion = Float64, overprediction = Float64,
        underprediction = Float64, coverage_50 = Float64,
        coverage_90 = Float64,
        bias = Float64, n_samples = Int,
        log_rel_to_baseline = Float64,
    )
)
forecast_overlay_df = _release_data(
    "forecast_overlay.csv",
    (;
        release = String, made_date = Date, stream = String, horizon = Int,
        target_date = Date, fit = String, observed = Float64,
        median = Float64, lo30 = Float64, hi30 = Float64, lo60 = Float64,
        hi60 = Float64, lo90 = Float64, hi90 = Float64,
    )
)
# The digitised onset triangle's own per-vintage total, as calendar dates,
# and the windows it makes unscoreable. A vintage that rereads the figure
# lower gives a window an increment the reports did not add, so the window
# is dropped from the scores and the overlay alike (see
# `drop_rescanned_onset_windows`). The current triangle is used to judge
# every release, since the scored truth is read off this one series.
_onset_vintage_dates = grid_date.(obs.onset_report_history.days)
_onset_vintage_totals = obs.onset_report_history.counts
_drop_rescanned(tbl) = drop_rescanned_onset_windows(
    tbl; vintage_dates = _onset_vintage_dates,
    vintage_totals = _onset_vintage_totals
)
forecast_scores_df = _drop_rescanned(forecast_scores_df)
forecast_overlay_df = _drop_rescanned(forecast_overlay_df)

# One row per (stream, fit) pooled over every horizon and release. The
# by-horizon and by-release detail tables carry the same columns at a finer
# grain (see src/scoring.jl). Every fit is kept here, since the
# relative-skill figure below compares the roles against each other. The
# tables rendered in this section select the joint role, and the individual
# fits are tabulated in their own section.
forecast_score_overview_table = forecast_score_overview(forecast_scores_df)
forecast_score_by_horizon_table = forecast_score_by_horizon(forecast_scores_df)
forecast_score_by_release_table = forecast_score_by_release(forecast_scores_df)

joint_score_overview_table = select_fit_role(
    forecast_score_overview_table, "joint"
)
joint_score_by_horizon_table = select_fit_role(
    forecast_score_by_horizon_table, "joint"
)
# The trailing `;` on this last assignment matters: without it, this whole
# setup chunk's last statement (the DataFrame it assigns) is Literate's
# implicitly displayed "result" for the chunk, on top of the deliberate
# display further down -- and a bare DataFrame is html-showable, so it
# goes out as a second, undisplayed-in-source `@raw html` block that (for
# a table this size) can itself hit the PCRE limit described above.
joint_score_by_release_table = select_fit_role(
    forecast_score_by_release_table, "joint"
);
````

```@raw html
</details>
```

The headline pools every horizon and release into one row per stream for the joint model: the mean CRPS and its decomposition, coverage, bias, and the relative skill against the persistence baseline, on both the natural and the log scale.
Each row also carries relative skill against the stream's own individual fit where one exists.


| stream | fit | n | crps | rel_to_baseline | log_crps | log_rel_to_baseline | rel_to_individual | log_rel_to_individual | dispersion | overprediction | underprediction | coverage_50 | coverage_90 | bias |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| confirmed cases | joint | 63 | 1258.97 | 7.28 | 0.242 | 1.77 | 0.8 | 0.74 | 1193.71 | 64.73 | 0.53 | 0.84 | 1 | 0.21 |
| confirmed deaths | joint | 61 | 1339.07 | 11.3 | 0.601 | 3.12 | 0.92 | 1.88 | 1212.65 | 6.16 | 120.25 | 0.36 | 0.98 | -0.48 |
| isolation beds | joint | 58 | 72.27 | 1.44 | 0.095 | 1.37 | 0.92 | 0.72 | 41.3 | 23.3 | 7.66 | 0.57 | 0.93 | 0.2 |
| onset reports | joint | 3 | 80.6 | 0.57 | 0.213 | 0.28 | 0.73 | 0.59 | 55.16 | 0.07 | 25.37 | 0.67 | 1 | -0.29 |
| recovered | joint | 11 | 139.4 | 1.68 | 0.376 | 1.09 | missing | missing | 120.02 | 6.31 | 13.07 | 1 | 1 | -0.18 |


The same relative skill against the baseline, by horizon: one panel per stream, one series per fit role, on a log-scaled skill axis with the reference line at one.

````julia
forecast_relative_skill_fig = plot_forecast_relative_skill(
    forecast_score_by_horizon_table
);

````
![](evaluation-45.png)

What that error is made of, by horizon: the mean CRPS split into its width, its overprediction and its underprediction, one stacked bar per horizon and fit role.

````julia
forecast_crps_by_horizon_fig = plot_forecast_crps_by_horizon(
    forecast_score_by_horizon_table
);

````
![](evaluation-47.png)

```@raw html
<details><summary>Scores by horizon</summary>
```


| stream | horizon | fit | n | crps | rel_to_baseline | log_crps | log_rel_to_baseline | rel_to_individual | log_rel_to_individual | dispersion | overprediction | underprediction | coverage_50 | coverage_90 | bias |
| --- | ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| confirmed cases | 7 | joint | 20 | 109.98 | 2.97 | 0.25 | 3.22 | 0.83 | 0.73 | 86.97 | 22.06 | 0.96 | 0.8 | 1 | 0.23 |
| confirmed cases | 14 | joint | 18 | 236.72 | 3.21 | 0.242 | 3 | 0.79 | 0.73 | 191.45 | 44.53 | 0.74 | 0.89 | 1 | 0.16 |
| confirmed cases | 21 | joint | 14 | 487.32 | 2.03 | 0.225 | 1.22 | 0.72 | 0.71 | 401.44 | 85.83 | 0.05 | 0.86 | 1 | 0.19 |
| confirmed cases | 28 | joint | 11 | 6002.93 | 12.07 | 0.252 | 0.91 | 0.86 | 0.82 | 5854.39 | 148.54 | 0 | 0.82 | 1 | 0.26 |
| confirmed deaths | 7 | joint | 20 | 90.41 | 3.69 | 0.647 | 5.95 | 1.35 | 1.8 | 35.28 | 0.4 | 54.73 | 0.4 | 0.95 | -0.49 |
| confirmed deaths | 14 | joint | 17 | 207.74 | 3.2 | 0.585 | 4.06 | 1.29 | 2.12 | 87.02 | 2.05 | 118.67 | 0.29 | 1 | -0.51 |
| confirmed deaths | 21 | joint | 13 | 408.31 | 2.44 | 0.596 | 2.42 | 1 | 2.05 | 218.11 | 8.33 | 181.86 | 0.15 | 1 | -0.51 |
| confirmed deaths | 28 | joint | 11 | 6457.76 | 20.55 | 0.548 | 1.53 | 0.64 | 1.53 | 6268.32 | 20.43 | 169.01 | 0.64 | 1 | -0.4 |
| isolation beds | 7 | joint | 18 | 65.92 | 1.59 | 0.086 | 1.47 | 0.74 | 0.56 | 30.66 | 27.48 | 7.78 | 0.44 | 0.94 | 0.29 |
| isolation beds | 14 | joint | 16 | 73.47 | 1.49 | 0.093 | 1.35 | 0.96 | 0.72 | 43.28 | 24.37 | 5.81 | 0.56 | 0.94 | 0.21 |
| isolation beds | 21 | joint | 13 | 75.62 | 1.46 | 0.099 | 1.4 | 1.1 | 0.89 | 46.79 | 19.48 | 9.35 | 0.69 | 0.92 | 0.12 |
| isolation beds | 28 | joint | 11 | 76.94 | 1.22 | 0.107 | 1.26 | 1.01 | 0.88 | 49.36 | 19.41 | 8.17 | 0.64 | 0.91 | 0.11 |
| onset reports | 7 | joint | 3 | 80.6 | 0.57 | 0.213 | 0.28 | 0.73 | 0.59 | 55.16 | 0.07 | 25.37 | 0.67 | 1 | -0.29 |
| recovered | 7 | joint | 3 | 26.96 | 2.03 | 0.319 | 1.67 | missing | missing | 20.93 | 1.28 | 4.76 | 1 | 1 | -0.18 |
| recovered | 14 | joint | 3 | 89.7 | 1.95 | 0.393 | 1.45 | missing | missing | 68.01 | 5.81 | 15.87 | 1 | 1 | -0.23 |
| recovered | 21 | joint | 3 | 297.16 | 3.24 | 0.422 | 1.15 | missing | missing | 269.9 | 16.05 | 11.2 | 1 | 1 | -0.09 |
| recovered | 28 | joint | 2 | 145.99 | 0.64 | 0.367 | 0.56 | missing | missing | 121.86 | 0 | 24.13 | 1 | 1 | -0.26 |


```@raw html
</details>
```

The same relative skill against the baseline, release by release, so a run of releases that lost to the baseline reads as a run rather than as an average.

````julia
forecast_skill_by_cutoff_fig = plot_forecast_skill_by_cutoff(
    forecast_score_by_release_table;
    title = "Relative skill against the baseline, by release"
);

````
![](evaluation-52.png)

```@raw html
<details><summary>Scores by release</summary>
```


| made_date | stream | fit | n | crps | rel_to_baseline | log_crps | log_rel_to_baseline | rel_to_individual | log_rel_to_individual | dispersion | overprediction | underprediction | coverage_50 | coverage_90 | bias |
| --- | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 2026-06-07 | confirmed cases | joint | 2 | 363.91 | 14.18 | 0.913 | 11.98 | missing | missing | 341.06 | 22.86 | 0 | 1 | 1 | 0.17 |
| 2026-06-10 | confirmed cases | joint | 3 | 534.87 | 6.69 | 0.791 | 4.97 | missing | missing | 435.67 | 99.2 | 0 | 1 | 1 | 0.27 |
| 2026-07-01 | confirmed cases | joint | 4 | 15862.78 | 74.74 | 0.583 | 2.84 | missing | missing | 15377.31 | 485.47 | 0 | 0 | 1 | 0.66 |
| 2026-07-06 | confirmed cases | joint | 4 | 421.28 | 1.31 | 0.269 | 0.95 | missing | missing | 311.85 | 109.43 | 0 | 0.5 | 1 | 0.49 |
| 2026-07-08 | confirmed cases | joint | 4 | 390.98 | 1.23 | 0.253 | 0.97 | missing | missing | 239.5 | 151.47 | 0 | 0.25 | 1 | 0.52 |
| 2026-07-23 | confirmed cases | joint | 4 | 335.57 | 0.86 | 0.169 | 0.66 | 1.02 | 0.82 | 299.27 | 35.68 | 0.62 | 1 | 1 | 0.13 |
| 2026-07-25 | confirmed cases | joint | 4 | 222.5 | 0.98 | 0.134 | 0.99 | 0.64 | 0.63 | 218.04 | 4.46 | 0 | 1 | 1 | 0.07 |
| 2026-07-26 | confirmed cases | joint | 4 | 235.69 | 1.07 | 0.141 | 1.06 | 0.79 | 0.67 | 229.65 | 5.68 | 0.36 | 1 | 1 | 0.06 |
| 2026-07-27 | confirmed cases | joint | 4 | 247.97 | 1.22 | 0.142 | 1.07 | 0.93 | 0.83 | 240.08 | 7.88 | 0 | 1 | 1 | 0.11 |
| 2026-07-31 | confirmed cases | joint | 4 | 259.64 | 1.85 | 0.149 | 1.69 | 0.64 | 0.68 | 235.73 | 23.76 | 0.15 | 1 | 1 | 0.13 |
| 2026-08-01 | confirmed cases | joint | 4 | 231.26 | 1.97 | 0.139 | 1.99 | 0.69 | 0.66 | 226.84 | 4.31 | 0.11 | 1 | 1 | 0.05 |
| 2026-08-02 | confirmed cases | joint | 4 | 278.45 | 2.12 | 0.145 | 1.75 | 0.94 | 0.85 | 263.77 | 14.68 | 0 | 1 | 1 | 0.14 |
| 2026-08-03 | confirmed cases | joint | 4 | 247.58 | 1.29 | 0.14 | 0.96 | 0.85 | 0.8 | 237.88 | 9.7 | 0 | 1 | 1 | 0.08 |
| 2026-08-04 | confirmed cases | joint | 3 | 153.72 | 3.53 | 0.123 | 3.34 | 0.63 | 0.68 | 152.58 | 1.15 | 0 | 1 | 1 | 0.04 |
| 2026-08-07 | confirmed cases | joint | 3 | 151.02 | 2.42 | 0.132 | 1.99 | 0.67 | 0.71 | 146.61 | 0 | 4.41 | 1 | 1 | -0.1 |
| 2026-08-11 | confirmed cases | joint | 2 | 98.54 | 2.69 | 0.122 | 2.85 | 0.62 | 0.69 | 91.08 | 0 | 7.46 | 1 | 1 | -0.18 |
| 2026-08-15 | confirmed cases | joint | 2 | 172.35 | 3.26 | 0.166 | 2.6 | 1.09 | 0.82 | 130.32 | 42.03 | 0 | 1 | 1 | 0.34 |
| 2026-08-17 | confirmed cases | joint | 2 | 202.19 | 2.31 | 0.187 | 1.84 | 1.4 | 0.89 | 125.12 | 77.07 | 0 | 1 | 1 | 0.4 |
| 2026-08-22 | confirmed cases | joint | 1 | 110.35 | 3.57 | 0.169 | 3.09 | 0.98 | 0.84 | 70.46 | 39.89 | 0 | 1 | 1 | 0.48 |
| 2026-08-24 | confirmed cases | joint | 1 | 77.26 | 3.26 | 0.127 | 2.92 | 0.87 | 0.82 | 33.95 | 43.31 | 0 | 0 | 1 | 0.63 |
| 2026-06-07 | confirmed deaths | joint | 1 | 41.49 | 1.87 | 0.903 | 2.55 | missing | missing | 36.94 | 0 | 4.55 | 1 | 1 | -0.2 |
| 2026-06-10 | confirmed deaths | joint | 2 | 64.93 | 2.27 | 0.735 | 3.05 | missing | missing | 62.7 | 0.99 | 1.24 | 1 | 1 | -0.03 |
| 2026-07-01 | confirmed deaths | joint | 4 | 17196 | 80.45 | 0.549 | 1.21 | missing | missing | 17103.11 | 92.89 | 0 | 1 | 1 | 0.28 |
| 2026-07-06 | confirmed deaths | joint | 4 | 197.2 | 0.89 | 0.871 | 2.29 | missing | missing | 110.35 | 0 | 86.85 | 0.25 | 1 | -0.49 |
| 2026-07-08 | confirmed deaths | joint | 4 | 177.59 | 0.78 | 0.713 | 1.84 | missing | missing | 84.29 | 0 | 93.3 | 0.25 | 1 | -0.61 |
| 2026-07-23 | confirmed deaths | joint | 4 | 271.49 | 1.67 | 0.484 | 2.45 | 0.39 | 1.08 | 172.61 | 0 | 98.88 | 1 | 1 | -0.45 |
| 2026-07-25 | confirmed deaths | joint | 4 | 241.83 | 1.98 | 0.504 | 3.41 | 1.02 | 2.15 | 106.46 | 0 | 135.37 | 0 | 1 | -0.57 |
| 2026-07-26 | confirmed deaths | joint | 4 | 268.46 | 2.31 | 0.554 | 4.16 | 1.05 | 2.17 | 101.59 | 0 | 166.87 | 0 | 1 | -0.64 |
| 2026-07-27 | confirmed deaths | joint | 4 | 249.47 | 2.29 | 0.534 | 3.99 | 0.8 | 2.12 | 109.3 | 0 | 140.17 | 0 | 1 | -0.59 |
| 2026-07-31 | confirmed deaths | joint | 4 | 274.55 | 2.31 | 0.56 | 3.96 | 0.88 | 1.94 | 108.79 | 0 | 165.76 | 0.25 | 1 | -0.58 |
| 2026-08-01 | confirmed deaths | joint | 4 | 313.48 | 2.69 | 0.683 | 4.66 | 1.12 | 2.7 | 94.48 | 0 | 219.01 | 0 | 1 | -0.68 |
| 2026-08-02 | confirmed deaths | joint | 4 | 258.17 | 2.85 | 0.507 | 4.38 | 1.01 | 1.83 | 128.19 | 0 | 129.98 | 0.25 | 1 | -0.53 |
| 2026-08-03 | confirmed deaths | joint | 4 | 268.89 | 1.96 | 0.58 | 2.91 | 1.2 | 2.37 | 123.56 | 0 | 145.33 | 0.25 | 1 | -0.55 |
| 2026-08-04 | confirmed deaths | joint | 3 | 253.16 | 5.21 | 0.682 | 9.44 | 1.51 | 2.5 | 66 | 0 | 187.16 | 0 | 1 | -0.7 |
| 2026-08-07 | confirmed deaths | joint | 3 | 320.06 | 4.83 | 0.923 | 8 | 1.94 | 2.84 | 53.21 | 0 | 266.85 | 0 | 0.67 | -0.84 |
| 2026-08-11 | confirmed deaths | joint | 2 | 208.92 | 10.85 | 0.849 | 19.33 | 1.81 | 2.66 | 48.56 | 0 | 160.36 | 0 | 1 | -0.8 |
| 2026-08-15 | confirmed deaths | joint | 2 | 107.78 | 2.58 | 0.323 | 3.74 | 0.99 | 1.02 | 79.26 | 0 | 28.52 | 1 | 1 | -0.34 |
| 2026-08-17 | confirmed deaths | joint | 2 | 105.06 | 2.11 | 0.307 | 2.72 | 0.97 | 1.13 | 83.42 | 0 | 21.64 | 1 | 1 | -0.3 |
| 2026-08-22 | confirmed deaths | joint | 1 | 62.71 | 2.4 | 0.306 | 3.48 | 0.78 | 0.63 | 43.58 | 0 | 19.13 | 1 | 1 | -0.35 |
| 2026-08-24 | confirmed deaths | joint | 1 | 17.41 | 1.57 | 0.055 | 1.49 | 0.18 | 0.1 | 15 | 2.41 | 0 | 1 | 1 | 0.22 |
| 2026-07-01 | isolation beds | joint | 4 | 90.19 | 1.48 | 0.138 | 1.52 | missing | missing | 13.48 | 0 | 76.71 | 0 | 0 | -0.99 |
| 2026-07-06 | isolation beds | joint | 4 | 42.87 | 0.62 | 0.066 | 0.64 | missing | missing | 21.07 | 0 | 21.8 | 0 | 1 | -0.64 |
| 2026-07-08 | isolation beds | joint | 4 | 54.75 | 1.54 | 0.073 | 1.54 | missing | missing | 11.67 | 43.08 | 0 | 0 | 1 | 0.8 |
| 2026-07-23 | isolation beds | joint | 4 | 80.06 | 2.72 | 0.101 | 2.63 | 1.44 | 1.14 | 29.39 | 50.67 | 0 | 0 | 1 | 0.58 |
| 2026-07-25 | isolation beds | joint | 4 | 80.84 | 2.19 | 0.105 | 2.16 | 1.56 | 1.28 | 32.65 | 48.19 | 0 | 0.5 | 1 | 0.53 |
| 2026-07-26 | isolation beds | joint | 4 | 73.83 | 1.97 | 0.097 | 1.9 | 1.2 | 0.96 | 36.96 | 36.87 | 0 | 0.5 | 1 | 0.46 |
| 2026-07-27 | isolation beds | joint | 4 | 71.58 | 2.47 | 0.093 | 2.35 | 1.2 | 0.95 | 39.52 | 32.06 | 0 | 0.5 | 1 | 0.44 |
| 2026-07-31 | isolation beds | joint | 4 | 71.05 | 1.68 | 0.088 | 1.63 | 0.86 | 0.67 | 51.69 | 19.36 | 0 | 0.75 | 1 | 0.34 |
| 2026-08-01 | isolation beds | joint | 4 | 72 | 1.12 | 0.094 | 1.09 | 0.81 | 0.61 | 64.51 | 7.48 | 0.02 | 1 | 1 | 0.18 |
| 2026-08-02 | isolation beds | joint | 4 | 83.69 | 1.63 | 0.102 | 1.47 | 0.99 | 0.72 | 52.67 | 31.02 | 0 | 0.75 | 1 | 0.39 |
| 2026-08-03 | isolation beds | joint | 4 | 78.83 | 1.85 | 0.099 | 1.72 | 0.99 | 0.74 | 60.42 | 18.41 | 0 | 1 | 1 | 0.3 |
| 2026-08-04 | isolation beds | joint | 3 | 102.79 | 1.98 | 0.141 | 1.8 | 1.19 | 0.99 | 55.77 | 47.03 | 0 | 0.67 | 1 | 0.46 |
| 2026-08-07 | isolation beds | joint | 3 | 74.19 | 1.29 | 0.103 | 1.39 | 0.69 | 0.59 | 57.75 | 0 | 16.44 | 1 | 1 | -0.27 |
| 2026-08-11 | isolation beds | joint | 2 | 49.12 | 0.37 | 0.068 | 0.34 | 0.54 | 0.44 | 48.55 | 0 | 0.57 | 1 | 1 | -0.07 |
| 2026-08-15 | isolation beds | joint | 2 | 60.39 | 0.78 | 0.068 | 0.7 | 0.5 | 0.36 | 57.4 | 2.99 | 0 | 1 | 1 | 0.13 |
| 2026-08-17 | isolation beds | joint | 2 | 73.47 | 1.56 | 0.083 | 1.33 | 0.75 | 0.52 | 57.83 | 15.64 | 0 | 1 | 1 | 0.32 |
| 2026-08-22 | isolation beds | joint | 1 | 37.09 | 0.81 | 0.04 | 0.75 | 0.24 | 0.17 | 35.32 | 1.77 | 0 | 1 | 1 | 0.13 |
| 2026-08-24 | isolation beds | joint | 1 | 58.73 | 2.43 | 0.065 | 2.18 | 0.44 | 0.3 | 36.03 | 22.7 | 0 | 1 | 1 | 0.44 |
| 2026-08-02 | onset reports | joint | 1 | 52.82 | 0.5 | 0.14 | 0.31 | 0.65 | 0.49 | 43.98 | 0 | 8.84 | 1 | 1 | -0.29 |
| 2026-08-11 | onset reports | joint | 1 | 122.26 | 0.51 | 0.373 | 0.23 | 0.86 | 0.61 | 55 | 0 | 67.26 | 0 | 1 | -0.6 |
| 2026-08-15 | onset reports | joint | 1 | 66.71 | 0.83 | 0.127 | 0.69 | 0.62 | 0.7 | 66.49 | 0.22 | 0 | 1 | 1 | 0.04 |
| 2026-07-01 | recovered | joint | 3 | 312.69 | 6.52 | 0.446 | 1.86 | missing | missing | 289.55 | 23.14 | 0 | 1 | 1 | 0.27 |
| 2026-07-06 | recovered | joint | 4 | 85.76 | 0.79 | 0.407 | 0.85 | missing | missing | 68.39 | 0 | 17.37 | 1 | 1 | -0.33 |
| 2026-07-08 | recovered | joint | 4 | 63.08 | 0.75 | 0.292 | 1 | missing | missing | 44.51 | 0 | 18.57 | 1 | 1 | -0.38 |


```@raw html
</details>
```

Forecasts made at each release against the value observed since, one panel per stream and horizon, the observed value in black.
The median and 90% interval are coloured by fit role: the persistence baseline, the stream's individual fit and the joint.
The x-axis is the date each forecast was made, so an incident stream's observed window pairs unambiguously with the forecast that made it.
Each panel's axis is cropped to a small multiple of what that stream actually reached, so one very wide interval cannot squash every other series flat.
An interval or median too wide for the panel is clamped at the top and marked with an open triangle rather than silently cut off.

```@raw html
<details><summary>Forecasts-versus-now overlay</summary>
```

````julia
forecast_overlay_fig = plot_forecast_overlay(
    scored_overlay(forecast_overlay_df)
);
````

```@raw html
</details>
```

![](evaluation-60.png)

### Forecast by province across releases

The archived provincial split of each release's forecast, scored against what each province went on to report, with a window holding a harmonisation-break day left unscored because that day's backfill is published for the country and not by province.

```@raw html
<details><summary>Load and summarise the province forecast scores</summary>
```

````julia
province_scores_df = _release_data(
    "province_forecast_scores.csv",
    (;
        release = String, made_date = Date, stream = String, horizon = Int,
        target_date = Date, fit = String, crps = Float64,
        log_crps = Float64, dispersion = Float64, overprediction = Float64,
        underprediction = Float64, coverage_50 = Float64,
        coverage_90 = Float64,
        bias = Float64, n_samples = Int,
        log_rel_to_baseline = Float64,
    )
)
# The joint patch model is the only model that forecasts the provinces, so
# there is no individual single-stream fit to compare against and `fit` is
# single-valued by construction. Both are dropped rather than rendered as
# columns that cannot vary.
#
# See the comment above `joint_score_by_release_table`'s assignment for why
# this setup chunk's last statement needs a trailing `;`.
province_score_overview_display = drop_degenerate_fit_column(
    drop_individual_fit_columns(forecast_score_overview(province_scores_df))
)
province_score_by_horizon_display = drop_degenerate_fit_column(
    drop_individual_fit_columns(
        forecast_score_by_horizon(province_scores_df)
    )
);
````

```@raw html
</details>
```


| stream | fit | n | crps | rel_to_baseline | log_crps | log_rel_to_baseline | dispersion | overprediction | underprediction | coverage_50 | coverage_90 | bias |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |


```@raw html
<details><summary>Province scores by horizon</summary>
```


| stream | horizon | fit | n | crps | rel_to_baseline | log_crps | log_rel_to_baseline | dispersion | overprediction | underprediction | coverage_50 | coverage_90 | bias |
| --- | ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |


```@raw html
</details>
```

## Frozen-fit forecast evaluation

The current model, frozen at earlier data cut-offs (see [Forecast-versus-frozen evaluation](@ref "Forecast-versus-frozen evaluation")), is scored the same way as the cross-release forecasts above, against the same persistence baseline.
The tables in this section are the frozen joint model's, one row per stream.
Each stream's own frozen fit is also scored where one exists, at the one-week-back cut-off and for still-reported streams only, and is carried by the skill figures rather than by the tables.
Every other cut-off carries the joint alone.
The May cut-offs predate the first reported bed occupancy and the first reported recoveries, so those windows are left unscored rather than scored against a series that had not started.
The baseline carries a weaker data-vintage guarantee than the cross-release one, since its snapshot was taken weeks after the frozen cut-off and can hold later revisions to earlier days (see [forecast scoring against a persistence baseline](@ref "Forecast scoring against a persistence baseline")).

```@raw html
<details><summary>Load and summarise the frozen-fit forecast scores</summary>
```

````julia
frozen_scores_df = _release_data(
    "forecast_scores_frozen.csv",
    (;
        release = String, made_date = Date, stream = String, horizon = Int,
        target_date = Date, fit = String, crps = Float64,
        log_crps = Float64, dispersion = Float64, overprediction = Float64,
        underprediction = Float64, coverage_50 = Float64,
        coverage_90 = Float64,
        bias = Float64, n_samples = Int,
        log_rel_to_baseline = Float64,
    )
)
frozen_overlay_df = _release_data(
    "forecast_overlay_frozen.csv",
    (;
        release = String, made_date = Date, stream = String, horizon = Int,
        target_date = Date, fit = String, observed = Float64,
        median = Float64, lo30 = Float64, hi30 = Float64, lo60 = Float64,
        hi60 = Float64, lo90 = Float64, hi90 = Float64,
    )
)
# Rows a superseded forecaster produced are dropped before anything is
# summarised or drawn, from the scores and the overlay alike, so the
# tables and the figures rest on one set of rows. See
# `drop_superseded_forecasts` for the one exclusion in force and why.
frozen_scores_df = _drop_rescanned(
    drop_superseded_forecasts(frozen_scores_df)
)
frozen_overlay_df = _drop_rescanned(
    drop_superseded_forecasts(frozen_overlay_df)
)

# The frozen joint carries `FROZEN_FIT`, so it is named as the joint role
# here and compared against each stream's own frozen fit. A release
# published before the archive had a `fit` column carries joint rows only.
frozen_score_overview_table = forecast_score_overview(
    frozen_scores_df; joint_fit = FROZEN_FIT
)
frozen_score_by_horizon_table = forecast_score_by_horizon(
    frozen_scores_df; joint_fit = FROZEN_FIT
)
frozen_score_by_release_table = forecast_score_by_release(
    frozen_scores_df; joint_fit = FROZEN_FIT
)

# The tables show the frozen joint alone, as the cross-release tables show
# the joint alone, so a row reads as one model at one cut-off rather than a
# stream interleaving two fits. The archive's single-stream frozen fits stay
# in the scored data and in the figures, which compare the roles against
# each other. Selecting the joint role leaves `fit` single-valued, so it is
# dropped and the model named in the prose instead.
_frozen_joint_only(tbl) = drop_degenerate_fit_column(
    select_fit_role(tbl, "joint")
)
frozen_score_overview_display = _frozen_joint_only(
    frozen_score_overview_table
)
frozen_score_by_horizon_display = _frozen_joint_only(
    frozen_score_by_horizon_table
)
frozen_score_by_release_display = _frozen_joint_only(
    frozen_score_by_release_table
)

# One row per release for the cut-offs more than one release forecast.
# See the comment above `joint_score_by_release_table`'s assignment for why
# this setup chunk's last statement needs a trailing `;`.
frozen_score_by_vintage_table = forecast_score_by_vintage(
    frozen_scores_df; joint_fit = FROZEN_FIT
)
frozen_score_by_vintage_display = _frozen_joint_only(
    frozen_score_by_vintage_table
);
````

```@raw html
</details>
```


| stream | n | crps | rel_to_baseline | log_crps | log_rel_to_baseline | rel_to_individual | log_rel_to_individual | dispersion | overprediction | underprediction | coverage_50 | coverage_90 | bias |
| --- | ---: | ---: | ---: | ---: | ---: | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| confirmed cases | 140 | 290.68 | 1.91 | 0.358 | 1.06 | missing | missing | 180.95 | 93.89 | 15.85 | 0.61 | 0.86 | 0.2 |
| confirmed deaths | 89 | 37.35 | 0.71 | 0.528 | 0.69 | missing | missing | 26.55 | 2.86 | 7.95 | 0.6 | 0.84 | -0.29 |
| isolation beds | 123 | 82.13 | 0.74 | 0.135 | 0.56 | missing | missing | 35.22 | 28.38 | 18.53 | 0.42 | 0.96 | 0.13 |
| onset reports | 3 | 98.55 | 0.75 | 0.25 | 0.4 | missing | missing | 56.63 | 9.34 | 32.57 | 0.67 | 1 | -0.08 |
| recovered | 55 | 199.18 | 2.45 | 0.453 | 1.89 | missing | missing | 178 | 17.52 | 3.66 | 0.93 | 0.95 | 0.07 |


The same relative skill against the baseline, by horizon, for the frozen cut-offs.

````julia
frozen_relative_skill_fig = plot_forecast_relative_skill(
    frozen_score_by_horizon_table
);

````
![](evaluation-75.png)

What that error is made of, by horizon, as in the cross-release section above.

````julia
frozen_crps_by_horizon_fig = plot_forecast_crps_by_horizon(
    frozen_score_by_horizon_table;
    title = "CRPS decomposition by horizon, frozen cut-offs"
);

````
![](evaluation-77.png)

```@raw html
<details><summary>Scores by horizon</summary>
```


| stream | horizon | n | crps | rel_to_baseline | log_crps | log_rel_to_baseline | rel_to_individual | log_rel_to_individual | dispersion | overprediction | underprediction | coverage_50 | coverage_90 | bias |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| confirmed cases | 7 | 68 | 91.5 | 1.06 | 0.396 | 0.77 | missing | missing | 42.44 | 16.43 | 32.62 | 0.59 | 0.76 | -0.07 |
| confirmed cases | 14 | 32 | 247.75 | 3.01 | 0.262 | 2.62 | missing | missing | 170.56 | 77.19 | 0 | 0.75 | 0.97 | 0.42 |
| confirmed cases | 21 | 30 | 615.99 | 3.03 | 0.415 | 2.14 | missing | missing | 340.44 | 275.55 | 0 | 0.43 | 0.9 | 0.54 |
| confirmed cases | 28 | 10 | 806.56 | 1.2 | 0.232 | 0.66 | missing | missing | 677.55 | 129.01 | 0 | 0.9 | 1 | 0.31 |
| confirmed deaths | 7 | 54 | 20.14 | 0.79 | 0.726 | 0.79 | missing | missing | 6.2 | 0.89 | 13.06 | 0.39 | 0.74 | -0.54 |
| confirmed deaths | 14 | 18 | 36.62 | 0.92 | 0.193 | 0.65 | missing | missing | 33.6 | 2.89 | 0.12 | 1 | 1 | 0.07 |
| confirmed deaths | 21 | 17 | 92.79 | 0.62 | 0.256 | 0.33 | missing | missing | 83.71 | 9.07 | 0 | 0.82 | 1 | 0.12 |
| isolation beds | 7 | 34 | 63.25 | 1.61 | 0.098 | 1.03 | missing | missing | 25.35 | 37.9 | 0 | 0.56 | 0.85 | 0.45 |
| isolation beds | 14 | 32 | 72.45 | 1.45 | 0.116 | 0.98 | missing | missing | 36.18 | 36.27 | 0 | 0.69 | 1 | 0.39 |
| isolation beds | 21 | 30 | 89.03 | 0.54 | 0.15 | 0.42 | missing | missing | 40.15 | 22.23 | 26.65 | 0.2 | 1 | -0.13 |
| isolation beds | 28 | 27 | 109.7 | 0.51 | 0.189 | 0.43 | missing | missing | 41.01 | 13.9 | 54.79 | 0.19 | 1 | -0.27 |
| onset reports | 7 | 3 | 98.55 | 0.75 | 0.25 | 0.4 | missing | missing | 56.63 | 9.34 | 32.57 | 0.67 | 1 | -0.08 |
| recovered | 7 | 17 | 79.54 | 2.82 | 0.471 | 2.35 | missing | missing | 46.7 | 29.92 | 2.93 | 0.82 | 0.88 | 0.13 |
| recovered | 14 | 15 | 142.54 | 3.67 | 0.441 | 2.98 | missing | missing | 111.42 | 26.81 | 4.31 | 0.93 | 0.93 | 0.08 |
| recovered | 21 | 13 | 204.13 | 1.81 | 0.44 | 1.54 | missing | missing | 196.69 | 1.24 | 6.2 | 1 | 1 | -0.02 |
| recovered | 28 | 10 | 481.11 | 2.47 | 0.459 | 1.18 | missing | missing | 476.8 | 3.67 | 0.65 | 1 | 1 | 0.05 |


```@raw html
</details>
```

The same relative skill against the baseline, cut-off by cut-off, pooled over the horizons each cut-off forecast.

````julia
frozen_skill_by_cutoff_fig = plot_forecast_skill_by_cutoff(
    frozen_score_by_release_table;
    xlabel = "Frozen cut-off",
    title = "Relative skill against the baseline, by frozen cut-off"
);

````
![](evaluation-82.png)

```@raw html
<details><summary>Scores by frozen cut-off</summary>
```


| made_date | stream | n | crps | rel_to_baseline | log_crps | log_rel_to_baseline | rel_to_individual | log_rel_to_individual | dispersion | overprediction | underprediction | coverage_50 | coverage_90 | bias |
| --- | --- | ---: | ---: | ---: | ---: | ---: | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 2026-05-23 | confirmed cases | 17 | 26.93 | 0.33 | 0.168 | 0.28 | missing | missing | 24.8 | 0 | 2.13 | 1 | 1 | -0.16 |
| 2026-05-27 | confirmed cases | 17 | 141.27 | 0.78 | 0.982 | 0.78 | missing | missing | 13.56 | 0 | 127.71 | 0 | 0.18 | -0.92 |
| 2026-06-08 | confirmed cases | 51 | 373.42 | 7.25 | 0.398 | 4 | missing | missing | 179.35 | 194.07 | 0 | 0.41 | 0.92 | 0.61 |
| 2026-07-16 | confirmed cases | 4 | 294.68 | 0.46 | 0.168 | 0.29 | missing | missing | 283.46 | 10.98 | 0.24 | 1 | 1 | 0.07 |
| 2026-07-18 | confirmed cases | 4 | 361.26 | 0.78 | 0.185 | 0.51 | missing | missing | 333.85 | 27.41 | 0 | 1 | 1 | 0.15 |
| 2026-07-19 | confirmed cases | 4 | 439.62 | 1.07 | 0.209 | 0.7 | missing | missing | 351.83 | 87.79 | 0 | 1 | 1 | 0.3 |
| 2026-07-20 | confirmed cases | 4 | 385.62 | 0.85 | 0.185 | 0.54 | missing | missing | 332.6 | 53.02 | 0 | 1 | 1 | 0.22 |
| 2026-07-24 | confirmed cases | 4 | 422.92 | 1.35 | 0.192 | 0.96 | missing | missing | 343.7 | 79.22 | 0 | 1 | 1 | 0.32 |
| 2026-07-25 | confirmed cases | 4 | 460.72 | 2.05 | 0.209 | 1.55 | missing | missing | 356.82 | 103.91 | 0 | 1 | 1 | 0.34 |
| 2026-07-26 | confirmed cases | 4 | 315.56 | 1.48 | 0.164 | 1.29 | missing | missing | 284.59 | 30.97 | 0 | 1 | 1 | 0.24 |
| 2026-07-27 | confirmed cases | 4 | 346.39 | 1.71 | 0.176 | 1.32 | missing | missing | 306.97 | 39.42 | 0 | 1 | 1 | 0.27 |
| 2026-07-28 | confirmed cases | 4 | 396.78 | 1.37 | 0.156 | 0.86 | missing | missing | 361.52 | 35.26 | 0 | 1 | 1 | 0.22 |
| 2026-07-31 | confirmed cases | 4 | 498.38 | 3.85 | 0.213 | 2.58 | missing | missing | 369.58 | 128.8 | 0 | 0.75 | 1 | 0.38 |
| 2026-08-04 | confirmed cases | 3 | 184.19 | 4.39 | 0.141 | 4.09 | missing | missing | 176.62 | 7.57 | 0 | 1 | 1 | 0.15 |
| 2026-08-08 | confirmed cases | 3 | 190.48 | 3.05 | 0.148 | 2.25 | missing | missing | 185.93 | 3.59 | 0.96 | 1 | 1 | 0.01 |
| 2026-08-10 | confirmed cases | 3 | 184.42 | 3.56 | 0.146 | 2.72 | missing | missing | 180.35 | 1.65 | 2.42 | 1 | 1 | -0.03 |
| 2026-08-15 | confirmed cases | 2 | 302.38 | 5.58 | 0.254 | 3.87 | missing | missing | 179 | 123.39 | 0 | 0 | 1 | 0.55 |
| 2026-08-17 | confirmed cases | 2 | 223.28 | 3.04 | 0.217 | 2.44 | missing | missing | 79.72 | 143.56 | 0 | 0 | 1 | 0.75 |
| 2026-08-18 | confirmed cases | 1 | 178.4 | 1.89 | 0.285 | 1.73 | missing | missing | 47.44 | 130.96 | 0 | 0 | 0 | 0.94 |
| 2026-08-23 | confirmed cases | 1 | 196.62 | 6.65 | 0.306 | 5.73 | missing | missing | 40.33 | 156.29 | 0 | 0 | 0 | 0.91 |
| 2026-05-23 | confirmed deaths | 17 | 10.23 | 0.41 | 0.479 | 0.34 | missing | missing | 3.3 | 0 | 6.93 | 0.18 | 1 | -0.61 |
| 2026-05-27 | confirmed deaths | 17 | 33.55 | 1.03 | 1.595 | 1.32 | missing | missing | 1.62 | 0 | 31.92 | 0 | 0.18 | -0.96 |
| 2026-06-08 | confirmed deaths | 51 | 47.72 | 0.7 | 0.223 | 0.49 | missing | missing | 43.12 | 3.68 | 0.92 | 0.94 | 1 | -0.02 |
| 2026-08-17 | confirmed deaths | 2 | 35.61 | 0.72 | 0.069 | 0.62 | missing | missing | 24.22 | 11.39 | 0 | 1 | 1 | 0.38 |
| 2026-08-18 | confirmed deaths | 1 | 37.68 | 0.64 | 0.124 | 0.65 | missing | missing | 14.8 | 22.88 | 0 | 0 | 1 | 0.72 |
| 2026-08-23 | confirmed deaths | 1 | 37.36 | 1.72 | 0.121 | 1.63 | missing | missing | 16.57 | 20.79 | 0 | 0 | 1 | 0.62 |
| 2026-06-08 | isolation beds | 68 | 74.09 | 0.44 | 0.15 | 0.38 | missing | missing | 38.21 | 2.4 | 33.48 | 0.5 | 1 | -0.21 |
| 2026-07-16 | isolation beds | 4 | 73.22 | 2.64 | 0.096 | 2.53 | missing | missing | 20.01 | 53.22 | 0 | 0 | 1 | 0.6 |
| 2026-07-18 | isolation beds | 4 | 77.78 | 2.76 | 0.103 | 2.61 | missing | missing | 18.71 | 59.07 | 0 | 0 | 1 | 0.61 |
| 2026-07-19 | isolation beds | 4 | 100.01 | 3.47 | 0.132 | 3.31 | missing | missing | 22.38 | 77.63 | 0 | 0 | 1 | 0.64 |
| 2026-07-20 | isolation beds | 4 | 82.23 | 3.2 | 0.106 | 3.05 | missing | missing | 19.89 | 62.33 | 0 | 0.25 | 1 | 0.61 |
| 2026-07-24 | isolation beds | 4 | 76.61 | 2.77 | 0.096 | 2.59 | missing | missing | 10.14 | 66.47 | 0 | 0 | 1 | 0.73 |
| 2026-07-25 | isolation beds | 4 | 106.43 | 2.71 | 0.136 | 2.62 | missing | missing | 13.84 | 92.59 | 0 | 0 | 0.75 | 0.73 |
| 2026-07-26 | isolation beds | 4 | 94.88 | 2.44 | 0.122 | 2.32 | missing | missing | 22.61 | 72.27 | 0 | 0.25 | 0.75 | 0.63 |
| 2026-07-27 | isolation beds | 4 | 85.56 | 2.9 | 0.113 | 2.81 | missing | missing | 38.92 | 46.65 | 0 | 0.5 | 1 | 0.48 |
| 2026-07-28 | isolation beds | 4 | 153.08 | 2.77 | 0.201 | 2.48 | missing | missing | 41.15 | 111.94 | 0 | 0.5 | 0.75 | 0.67 |
| 2026-07-31 | isolation beds | 4 | 119.06 | 3.04 | 0.139 | 2.78 | missing | missing | 41.36 | 77.7 | 0 | 0.5 | 0.75 | 0.57 |
| 2026-08-04 | isolation beds | 3 | 114.97 | 2.2 | 0.158 | 2.01 | missing | missing | 55.15 | 59.82 | 0 | 0.67 | 0.67 | 0.51 |
| 2026-08-08 | isolation beds | 3 | 54.61 | 0.87 | 0.07 | 0.86 | missing | missing | 44.63 | 9.34 | 0.65 | 0.67 | 1 | 0.2 |
| 2026-08-10 | isolation beds | 3 | 46.18 | 0.87 | 0.062 | 0.87 | missing | missing | 45.96 | 0.02 | 0.2 | 1 | 1 | 0 |
| 2026-08-15 | isolation beds | 2 | 83.21 | 1.1 | 0.088 | 0.93 | missing | missing | 59.78 | 23.43 | 0 | 1 | 1 | 0.36 |
| 2026-08-17 | isolation beds | 2 | 91.77 | 3.55 | 0.1 | 3.04 | missing | missing | 54.51 | 37.26 | 0 | 0.5 | 1 | 0.54 |
| 2026-08-18 | isolation beds | 1 | 80.15 | 2.76 | 0.094 | 2.38 | missing | missing | 36.09 | 44.06 | 0 | 0 | 1 | 0.67 |
| 2026-08-23 | isolation beds | 1 | 110.68 | 3.63 | 0.122 | 3.38 | missing | missing | 35.11 | 75.57 | 0 | 0 | 1 | 0.85 |
| 2026-07-31 | onset reports | 1 | 78.09 | 1.13 | 0.188 | 1.07 | missing | missing | 54.83 | 23.26 | 0 | 1 | 1 | 0.36 |
| 2026-08-10 | onset reports | 1 | 144.5 | 0.6 | 0.441 | 0.3 | missing | missing | 46.78 | 0 | 97.72 | 0 | 1 | -0.73 |
| 2026-08-15 | onset reports | 1 | 73.06 | 0.86 | 0.121 | 0.55 | missing | missing | 68.29 | 4.77 | 0 | 1 | 1 | 0.15 |
| 2026-07-16 | recovered | 4 | 139.34 | 1.1 | 0.469 | 1.09 | missing | missing | 135.5 | 0 | 3.84 | 1 | 1 | -0.1 |
| 2026-07-18 | recovered | 4 | 149.61 | 2.29 | 0.455 | 2.24 | missing | missing | 146.94 | 2.67 | 0 | 1 | 1 | 0.1 |
| 2026-07-19 | recovered | 4 | 181.47 | 2.24 | 0.482 | 2.11 | missing | missing | 179.35 | 2.12 | 0 | 1 | 1 | 0.07 |
| 2026-07-20 | recovered | 4 | 171.81 | 2.13 | 0.417 | 1.93 | missing | missing | 169.63 | 2.18 | 0 | 1 | 1 | 0.06 |
| 2026-07-24 | recovered | 4 | 179.37 | 2.26 | 0.409 | 1.93 | missing | missing | 175.36 | 4.01 | 0 | 1 | 1 | 0.13 |
| 2026-07-25 | recovered | 4 | 177.7 | 2.04 | 0.379 | 1.64 | missing | missing | 173.52 | 4.18 | 0 | 1 | 1 | 0.12 |
| 2026-07-26 | recovered | 4 | 139.99 | 1.51 | 0.381 | 1.57 | missing | missing | 138.99 | 0 | 1 | 1 | 1 | -0.05 |
| 2026-07-27 | recovered | 4 | 219.87 | 2.19 | 0.483 | 1.81 | missing | missing | 218.38 | 0 | 1.49 | 1 | 1 | -0.06 |
| 2026-07-28 | recovered | 4 | 544.07 | 3.75 | 0.369 | 0.83 | missing | missing | 543.25 | 0.5 | 0.33 | 1 | 1 | 0.02 |
| 2026-07-31 | recovered | 4 | 249.16 | 2.39 | 0.402 | 1.48 | missing | missing | 242.91 | 6.25 | 0 | 1 | 1 | 0.12 |
| 2026-08-04 | recovered | 3 | 118.93 | 2.56 | 0.444 | 3.02 | missing | missing | 110.06 | 0 | 8.87 | 1 | 1 | -0.19 |
| 2026-08-08 | recovered | 3 | 98.48 | 1.64 | 0.407 | 2.14 | missing | missing | 82.99 | 0 | 15.49 | 1 | 1 | -0.25 |
| 2026-08-10 | recovered | 3 | 119.75 | 1.97 | 0.524 | 2.58 | missing | missing | 85.84 | 0 | 33.9 | 1 | 1 | -0.32 |
| 2026-08-15 | recovered | 2 | 157.51 | 9.82 | 0.402 | 5.52 | missing | missing | 132.72 | 24.79 | 0 | 1 | 1 | 0.31 |
| 2026-08-17 | recovered | 2 | 309.21 | 12.68 | 0.773 | 6.65 | missing | missing | 57.82 | 251.39 | 0 | 0 | 0 | 0.92 |
| 2026-08-18 | recovered | 1 | 209.96 | 7 | 0.79 | 4.69 | missing | missing | 34.12 | 175.84 | 0 | 0 | 0 | 0.9 |
| 2026-08-23 | recovered | 1 | 190.8 | 16.84 | 0.679 | 9.4 | missing | missing | 43.12 | 147.68 | 0 | 0 | 1 | 0.87 |


```@raw html
</details>
```

### Frozen skill by release

Skill at each cut-off more than one release forecast, one point per release rather than pooled across releases.
Releases run in the order they were cut, evenly spaced rather than to calendar scale.

```@raw html
<details><summary>Frozen skill per release</summary>
```

````julia
frozen_skill_by_vintage_fig = plot_forecast_skill_by_vintage(
    frozen_score_by_vintage_table
);
````

```@raw html
</details>
```

![](evaluation-90.png)

```@raw html
<details><summary>Frozen scores by release</summary>
```


| stream | release | release_date | n | crps | rel_to_baseline | log_crps | log_rel_to_baseline | rel_to_individual | log_rel_to_individual | dispersion | overprediction | underprediction | coverage_50 | coverage_90 | bias |
| --- | --- | --- | ---: | ---: | ---: | ---: | ---: | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| confirmed cases | results-1359 | 2026-07-16 | 5 | 242.1 | 2.79 | 0.454 | 1.03 | missing | missing | 121.44 | 95.46 | 25.2 | 0.6 | 1 | 0.07 |
| confirmed cases | results-1391 | 2026-07-18 | 5 | 242.1 | 2.97 | 0.454 | 1.07 | missing | missing | 121.44 | 95.46 | 25.2 | 0.6 | 1 | 0.07 |
| confirmed cases | results-1394 | 2026-07-19 | 5 | 242.1 | 2.93 | 0.454 | 1.06 | missing | missing | 121.44 | 95.46 | 25.2 | 0.6 | 1 | 0.07 |
| confirmed cases | results-1441 | 2026-07-20 | 5 | 245.99 | 2.99 | 0.466 | 1.06 | missing | missing | 113.62 | 106.92 | 25.45 | 0.2 | 0.8 | 0.16 |
| confirmed cases | results-1446 | 2026-07-24 | 5 | 245.99 | 2.84 | 0.466 | 1.04 | missing | missing | 113.62 | 106.92 | 25.45 | 0.2 | 0.8 | 0.16 |
| confirmed cases | results-v1.12.0 | 2026-07-25 | 5 | 245.99 | 2.82 | 0.466 | 1.07 | missing | missing | 113.62 | 106.92 | 25.45 | 0.2 | 0.8 | 0.16 |
| confirmed cases | results-1462 | 2026-07-26 | 5 | 247.61 | 3.04 | 0.458 | 1.06 | missing | missing | 122.27 | 99.46 | 25.89 | 0.6 | 0.8 | 0.13 |
| confirmed cases | results-1469 | 2026-07-27 | 5 | 247.61 | 2.99 | 0.458 | 1.08 | missing | missing | 122.27 | 99.46 | 25.89 | 0.6 | 0.8 | 0.13 |
| confirmed cases | results-1479 | 2026-07-28 | 5 | 247.61 | 2.93 | 0.458 | 1.05 | missing | missing | 122.27 | 99.46 | 25.89 | 0.6 | 0.8 | 0.13 |
| confirmed cases | results-1489 | 2026-07-31 | 5 | 247.61 | 3.03 | 0.458 | 1.08 | missing | missing | 122.27 | 99.46 | 25.89 | 0.6 | 0.8 | 0.13 |
| confirmed cases | results-v1.13.2 | 2026-08-04 | 5 | 247.61 | 2.94 | 0.458 | 1.07 | missing | missing | 122.27 | 99.46 | 25.89 | 0.6 | 0.8 | 0.13 |
| confirmed cases | results-1505 | 2026-08-08 | 5 | 247.61 | 3.07 | 0.458 | 1.08 | missing | missing | 122.27 | 99.46 | 25.89 | 0.6 | 0.8 | 0.13 |
| confirmed cases | results-1527 | 2026-08-10 | 5 | 247.61 | 3.01 | 0.458 | 1.06 | missing | missing | 122.27 | 99.46 | 25.89 | 0.6 | 0.8 | 0.13 |
| confirmed cases | results-1573 | 2026-08-15 | 5 | 256.51 | 3 | 0.471 | 1.06 | missing | missing | 131.58 | 99.04 | 25.89 | 0.4 | 0.8 | 0.13 |
| confirmed cases | results-1591 | 2026-08-17 | 5 | 304.06 | 3.73 | 0.509 | 1.19 | missing | missing | 85.78 | 190.69 | 27.59 | 0.2 | 0.6 | 0.3 |
| confirmed cases | results-v1.15.0 | 2026-08-18 | 5 | 304.06 | 3.64 | 0.509 | 1.18 | missing | missing | 85.78 | 190.69 | 27.59 | 0.2 | 0.6 | 0.3 |
| confirmed cases | results-v1.16.0 | 2026-08-23 | 5 | 318.57 | 3.86 | 0.513 | 1.21 | missing | missing | 95.65 | 195.73 | 27.19 | 0.2 | 0.4 | 0.3 |
| confirmed deaths | results-1359 | 2026-07-16 | 5 | 36.61 | 0.68 | 0.553 | 0.69 | missing | missing | 28.04 | 0.07 | 8.49 | 0.6 | 0.8 | -0.39 |
| confirmed deaths | results-1391 | 2026-07-18 | 5 | 36.61 | 0.69 | 0.553 | 0.69 | missing | missing | 28.04 | 0.07 | 8.49 | 0.6 | 0.8 | -0.39 |
| confirmed deaths | results-1394 | 2026-07-19 | 5 | 36.61 | 0.7 | 0.553 | 0.71 | missing | missing | 28.04 | 0.07 | 8.49 | 0.6 | 0.8 | -0.39 |
| confirmed deaths | results-1441 | 2026-07-20 | 5 | 37.95 | 0.72 | 0.579 | 0.73 | missing | missing | 28.76 | 0.44 | 8.76 | 0.6 | 0.8 | -0.35 |
| confirmed deaths | results-1446 | 2026-07-24 | 5 | 37.95 | 0.71 | 0.579 | 0.71 | missing | missing | 28.76 | 0.44 | 8.76 | 0.6 | 0.8 | -0.35 |
| confirmed deaths | results-v1.12.0 | 2026-07-25 | 5 | 37.95 | 0.73 | 0.579 | 0.74 | missing | missing | 28.76 | 0.44 | 8.76 | 0.6 | 0.8 | -0.35 |
| confirmed deaths | results-1462 | 2026-07-26 | 5 | 37.32 | 0.71 | 0.635 | 0.8 | missing | missing | 27.48 | 0 | 9.84 | 0.6 | 0.8 | -0.43 |
| confirmed deaths | results-1469 | 2026-07-27 | 5 | 37.32 | 0.72 | 0.635 | 0.81 | missing | missing | 27.48 | 0 | 9.84 | 0.6 | 0.8 | -0.43 |
| confirmed deaths | results-1479 | 2026-07-28 | 5 | 37.32 | 0.72 | 0.635 | 0.8 | missing | missing | 27.48 | 0 | 9.84 | 0.6 | 0.8 | -0.43 |
| confirmed deaths | results-1489 | 2026-07-31 | 5 | 37.32 | 0.71 | 0.635 | 0.8 | missing | missing | 27.48 | 0 | 9.84 | 0.6 | 0.8 | -0.43 |
| confirmed deaths | results-v1.13.2 | 2026-08-04 | 5 | 37.32 | 0.71 | 0.635 | 0.8 | missing | missing | 27.48 | 0 | 9.84 | 0.6 | 0.8 | -0.43 |
| confirmed deaths | results-1505 | 2026-08-08 | 5 | 37.32 | 0.72 | 0.635 | 0.81 | missing | missing | 27.48 | 0 | 9.84 | 0.6 | 0.8 | -0.43 |
| confirmed deaths | results-1527 | 2026-08-10 | 5 | 37.32 | 0.71 | 0.635 | 0.8 | missing | missing | 27.48 | 0 | 9.84 | 0.6 | 0.8 | -0.43 |
| confirmed deaths | results-1573 | 2026-08-15 | 5 | 37.79 | 0.72 | 0.631 | 0.79 | missing | missing | 27.84 | 0.14 | 9.82 | 0.6 | 0.8 | -0.4 |
| confirmed deaths | results-1591 | 2026-08-17 | 5 | 36.94 | 0.71 | 0.288 | 0.36 | missing | missing | 20.9 | 12.14 | 3.9 | 0.6 | 1 | 0 |
| confirmed deaths | results-v1.15.0 | 2026-08-18 | 5 | 36.94 | 0.71 | 0.288 | 0.37 | missing | missing | 20.9 | 12.14 | 3.9 | 0.6 | 1 | 0 |
| confirmed deaths | results-v1.16.0 | 2026-08-23 | 5 | 38.99 | 0.73 | 0.27 | 0.33 | missing | missing | 24.17 | 11.59 | 3.22 | 0.6 | 1 | 0.04 |
| isolation beds | results-1359 | 2026-07-16 | 4 | 77.62 | 0.45 | 0.158 | 0.4 | missing | missing | 41.55 | 1.54 | 34.54 | 0.5 | 1 | -0.22 |
| isolation beds | results-1391 | 2026-07-18 | 4 | 77.62 | 0.46 | 0.158 | 0.41 | missing | missing | 41.55 | 1.54 | 34.54 | 0.5 | 1 | -0.22 |
| isolation beds | results-1394 | 2026-07-19 | 4 | 77.62 | 0.46 | 0.158 | 0.41 | missing | missing | 41.55 | 1.54 | 34.54 | 0.5 | 1 | -0.22 |
| isolation beds | results-1441 | 2026-07-20 | 4 | 77.32 | 0.46 | 0.158 | 0.41 | missing | missing | 38.49 | 1.02 | 37.81 | 0.5 | 1 | -0.22 |
| isolation beds | results-1446 | 2026-07-24 | 4 | 77.32 | 0.43 | 0.158 | 0.38 | missing | missing | 38.49 | 1.02 | 37.81 | 0.5 | 1 | -0.22 |
| isolation beds | results-v1.12.0 | 2026-07-25 | 4 | 77.32 | 0.46 | 0.158 | 0.4 | missing | missing | 38.49 | 1.02 | 37.81 | 0.5 | 1 | -0.22 |
| isolation beds | results-1462 | 2026-07-26 | 4 | 71.12 | 0.43 | 0.144 | 0.38 | missing | missing | 38.04 | 3.53 | 29.56 | 0.5 | 1 | -0.2 |
| isolation beds | results-1469 | 2026-07-27 | 4 | 71.12 | 0.43 | 0.144 | 0.38 | missing | missing | 38.04 | 3.53 | 29.56 | 0.5 | 1 | -0.2 |
| isolation beds | results-1479 | 2026-07-28 | 4 | 71.12 | 0.43 | 0.144 | 0.38 | missing | missing | 38.04 | 3.53 | 29.56 | 0.5 | 1 | -0.2 |
| isolation beds | results-1489 | 2026-07-31 | 4 | 71.12 | 0.44 | 0.144 | 0.39 | missing | missing | 38.04 | 3.53 | 29.56 | 0.5 | 1 | -0.2 |
| isolation beds | results-v1.13.2 | 2026-08-04 | 4 | 71.12 | 0.41 | 0.144 | 0.36 | missing | missing | 38.04 | 3.53 | 29.56 | 0.5 | 1 | -0.2 |
| isolation beds | results-1505 | 2026-08-08 | 4 | 71.12 | 0.41 | 0.144 | 0.36 | missing | missing | 38.04 | 3.53 | 29.56 | 0.5 | 1 | -0.2 |
| isolation beds | results-1527 | 2026-08-10 | 4 | 71.12 | 0.41 | 0.144 | 0.36 | missing | missing | 38.04 | 3.53 | 29.56 | 0.5 | 1 | -0.2 |
| isolation beds | results-1573 | 2026-08-15 | 4 | 76.5 | 0.44 | 0.154 | 0.38 | missing | missing | 35.78 | 2.15 | 38.57 | 0.5 | 1 | -0.22 |
| isolation beds | results-1591 | 2026-08-17 | 4 | 73.46 | 0.43 | 0.147 | 0.37 | missing | missing | 36.1 | 2.13 | 35.23 | 0.5 | 1 | -0.2 |
| isolation beds | results-v1.15.0 | 2026-08-18 | 4 | 73.46 | 0.43 | 0.147 | 0.37 | missing | missing | 36.1 | 2.13 | 35.23 | 0.5 | 1 | -0.2 |
| isolation beds | results-v1.16.0 | 2026-08-23 | 4 | 73.46 | 0.44 | 0.147 | 0.38 | missing | missing | 35.25 | 2.05 | 36.16 | 0.5 | 1 | -0.2 |


```@raw html
</details>
```

The frozen forecasts made at each cut-off against the value observed since, one panel per stream and horizon, the observed value in black.
Each panel carries the frozen forecast and the persistence baseline, coloured as in the cross-release overlay above, and the x-axis is the cut-off each forecast was made from.

```@raw html
<details><summary>Frozen-fit forecasts-versus-now overlay</summary>
```

````julia
# The frozen joint and the persistence baseline only. A single-stream
# frozen fit exists at the one-week-back cut-off alone, so its series lands
# on one made date of a panel spanning every cut-off, overplotting the
# joint point it sits beside rather than reading as a second series. The
# single-stream frozen fits are compared against the joint in the skill
# figures and in the validation plot at that cut-off.
frozen_overlay_fig = plot_forecast_overlay(
    scored_overlay(
        vcat(
            select_fit_role(frozen_overlay_df, "joint"),
            select_fit_role(frozen_overlay_df, "baseline")
        )
    )
);
````

```@raw html
</details>
```

![](evaluation-98.png)

The frozen re-fits below freeze the renewal data to an earlier cut-off and re-fit, so that a change driven by newer data can be distinguished from one driven by a change of method.
Each uses the full headline settings: 1000 draws across two chains.

```@raw html
<details><summary>Freeze the renewal data to a cut-off and re-fit</summary>
```

````julia
# Frozen re-fits and released_df are prepared in the setup block above.
````

```@raw html
</details>
```

## Individual fits against the baseline

This section carries the same cross-release forecast scoring as [Forecast scoring across releases](@ref "Forecast scoring across releases") above, for each stream's own individual fit rather than the joint, against the same persistence baseline.

```@raw html
<details><summary>Individual-fit rows of the cross-release scores</summary>
```

````julia
# The relative skill against a stream's individual fit is only ever
# computed on the joint model's row, so on these rows it is missing by
# construction and the column is dropped rather than shown empty.
individual_score_overview_table = drop_individual_fit_columns(
    select_fit_role(forecast_score_overview_table, "individual")
)
individual_score_by_horizon_table = drop_individual_fit_columns(
    select_fit_role(forecast_score_by_horizon_table, "individual")
)
# See the comment above `joint_score_by_release_table`'s assignment for why
# this setup chunk's last statement needs a trailing `;`.
individual_score_by_release_table = drop_individual_fit_columns(
    select_fit_role(forecast_score_by_release_table, "individual")
);
````

```@raw html
</details>
```


| stream | fit | n | crps | rel_to_baseline | log_crps | log_rel_to_baseline | dispersion | overprediction | underprediction | coverage_50 | coverage_90 | bias |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| confirmed cases | confirmed | 46 | 278.45 | 1.78 | 0.194 | 1.79 | 269.74 | 6.36 | 2.35 | 1 | 1 | 0 |
| confirmed deaths | confirmed_deaths | 46 | 264.4 | 2.7 | 0.297 | 2.27 | 220.56 | 40.31 | 3.53 | 0.93 | 1 | 0.06 |
| isolation beds | treatment | 46 | 81.57 | 1.67 | 0.133 | 2.01 | 69.65 | 0.13 | 11.79 | 0.93 | 1 | -0.21 |
| onset reports | onsets | 3 | 110.47 | 0.78 | 0.359 | 0.47 | 77.35 | 2.38 | 30.74 | 0.67 | 1 | -0.22 |


The same relative skill against the baseline, by horizon, one panel per stream (dataset), for each stream's own individual fit.

````julia
individual_relative_skill_fig = plot_forecast_relative_skill(
    individual_score_by_horizon_table;
    empty_message = "Empty: no release old enough for its targets to " *
        "have been observed carries an individual-stream " *
        "forecast. Not a missing forecast."
);

````
![](evaluation-109.png)

```@raw html
<details><summary>Scores by horizon</summary>
```


| stream | horizon | fit | n | crps | rel_to_baseline | log_crps | log_rel_to_baseline | dispersion | overprediction | underprediction | coverage_50 | coverage_90 | bias |
| --- | ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| confirmed cases | 7 | confirmed | 15 | 98.21 | 2.47 | 0.186 | 2.64 | 93.59 | 0.55 | 4.07 | 1 | 1 | -0.08 |
| confirmed cases | 14 | confirmed | 13 | 213.74 | 2.95 | 0.189 | 2.89 | 208.87 | 1.35 | 3.53 | 1 | 1 | -0.02 |
| confirmed cases | 21 | confirmed | 10 | 371.29 | 1.75 | 0.195 | 1.45 | 365.22 | 5.96 | 0.11 | 1 | 1 | 0.06 |
| confirmed cases | 28 | confirmed | 8 | 605.51 | 1.37 | 0.219 | 1 | 579.63 | 25.89 | 0 | 1 | 1 | 0.13 |
| confirmed deaths | 7 | confirmed_deaths | 15 | 73.26 | 2.67 | 0.294 | 3.17 | 62.47 | 2.18 | 8.61 | 1 | 1 | -0.11 |
| confirmed deaths | 14 | confirmed_deaths | 13 | 168.91 | 3 | 0.272 | 2.93 | 153.71 | 12.67 | 2.52 | 0.92 | 1 | 0.03 |
| confirmed deaths | 21 | confirmed_deaths | 10 | 343.65 | 2.6 | 0.296 | 1.9 | 298.26 | 45.33 | 0.06 | 0.9 | 1 | 0.16 |
| confirmed deaths | 28 | confirmed_deaths | 8 | 678.91 | 2.68 | 0.345 | 1.47 | 528.5 | 150.41 | 0 | 0.88 | 1 | 0.31 |
| isolation beds | 7 | treatment | 15 | 88.75 | 2.44 | 0.151 | 2.94 | 78.16 | 0.31 | 10.28 | 1 | 1 | -0.14 |
| isolation beds | 14 | treatment | 13 | 80.44 | 1.63 | 0.134 | 1.97 | 73.76 | 0.09 | 6.59 | 1 | 1 | -0.14 |
| isolation beds | 21 | treatment | 10 | 74.4 | 1.42 | 0.117 | 1.68 | 62.14 | 0 | 12.25 | 0.9 | 1 | -0.23 |
| isolation beds | 28 | treatment | 8 | 78.9 | 1.19 | 0.12 | 1.38 | 56.43 | 0 | 22.47 | 0.75 | 1 | -0.41 |
| onset reports | 7 | onsets | 3 | 110.47 | 0.78 | 0.359 | 0.47 | 77.35 | 2.38 | 30.74 | 0.67 | 1 | -0.22 |


```@raw html
</details>
```

```@raw html
<details><summary>Scores by release</summary>
```


| made_date | stream | fit | n | crps | rel_to_baseline | log_crps | log_rel_to_baseline | dispersion | overprediction | underprediction | coverage_50 | coverage_90 | bias |
| --- | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 2026-07-23 | confirmed cases | confirmed | 4 | 327.81 | 0.84 | 0.206 | 0.81 | 319.43 | 4.84 | 3.55 | 1 | 1 | -0.03 |
| 2026-07-25 | confirmed cases | confirmed | 4 | 346.24 | 1.53 | 0.213 | 1.57 | 343.88 | 2.3 | 0.06 | 1 | 1 | 0.02 |
| 2026-07-26 | confirmed cases | confirmed | 4 | 299.07 | 1.36 | 0.211 | 1.59 | 291.39 | 0.16 | 7.51 | 1 | 1 | -0.12 |
| 2026-07-27 | confirmed cases | confirmed | 4 | 267.62 | 1.32 | 0.171 | 1.3 | 262.03 | 5.38 | 0.21 | 1 | 1 | 0.04 |
| 2026-07-31 | confirmed cases | confirmed | 4 | 402.87 | 2.88 | 0.219 | 2.49 | 383.82 | 19.05 | 0 | 1 | 1 | 0.09 |
| 2026-08-01 | confirmed cases | confirmed | 4 | 336.31 | 2.87 | 0.209 | 3 | 322.98 | 13.33 | 0 | 1 | 1 | 0.08 |
| 2026-08-02 | confirmed cases | confirmed | 4 | 296.98 | 2.26 | 0.171 | 2.05 | 286.15 | 6.72 | 4.11 | 1 | 1 | -0.03 |
| 2026-08-03 | confirmed cases | confirmed | 4 | 291.8 | 1.52 | 0.175 | 1.2 | 281.42 | 10.08 | 0.3 | 1 | 1 | 0.05 |
| 2026-08-04 | confirmed cases | confirmed | 3 | 243.12 | 5.59 | 0.181 | 4.92 | 234.75 | 8.32 | 0.05 | 1 | 1 | 0.05 |
| 2026-08-07 | confirmed cases | confirmed | 3 | 225.88 | 3.62 | 0.187 | 2.8 | 218.25 | 2.14 | 5.5 | 1 | 1 | -0.06 |
| 2026-08-11 | confirmed cases | confirmed | 2 | 159.22 | 4.35 | 0.178 | 4.16 | 154.72 | 2.82 | 1.67 | 1 | 1 | 0 |
| 2026-08-15 | confirmed cases | confirmed | 2 | 158.81 | 3 | 0.203 | 3.17 | 154.14 | 0 | 4.67 | 1 | 1 | -0.1 |
| 2026-08-17 | confirmed cases | confirmed | 2 | 144.46 | 1.65 | 0.209 | 2.06 | 136.61 | 0 | 7.84 | 1 | 1 | -0.15 |
| 2026-08-22 | confirmed cases | confirmed | 1 | 113.08 | 3.65 | 0.201 | 3.65 | 111.45 | 1.63 | 0 | 1 | 1 | 0.05 |
| 2026-08-24 | confirmed cases | confirmed | 1 | 88.92 | 3.75 | 0.156 | 3.58 | 82.49 | 6.43 | 0 | 1 | 1 | 0.14 |
| 2026-07-23 | confirmed deaths | confirmed_deaths | 4 | 701.54 | 4.31 | 0.448 | 2.27 | 383.29 | 318.25 | 0 | 0.25 | 1 | 0.56 |
| 2026-07-25 | confirmed deaths | confirmed_deaths | 4 | 237.78 | 1.95 | 0.234 | 1.58 | 214.28 | 23.5 | 0 | 1 | 1 | 0.18 |
| 2026-07-26 | confirmed deaths | confirmed_deaths | 4 | 256.46 | 2.2 | 0.255 | 1.91 | 231.65 | 24.81 | 0 | 1 | 1 | 0.13 |
| 2026-07-27 | confirmed deaths | confirmed_deaths | 4 | 310.27 | 2.84 | 0.252 | 1.88 | 280.84 | 29.43 | 0 | 1 | 1 | 0.21 |
| 2026-07-31 | confirmed deaths | confirmed_deaths | 4 | 312.36 | 2.63 | 0.289 | 2.05 | 292.21 | 20.08 | 0.07 | 1 | 1 | 0.08 |
| 2026-08-01 | confirmed deaths | confirmed_deaths | 4 | 280.92 | 2.41 | 0.253 | 1.73 | 261.46 | 19.46 | 0 | 1 | 1 | 0.1 |
| 2026-08-02 | confirmed deaths | confirmed_deaths | 4 | 256.05 | 2.83 | 0.276 | 2.38 | 238.23 | 16.88 | 0.94 | 1 | 1 | 0.03 |
| 2026-08-03 | confirmed deaths | confirmed_deaths | 4 | 224.5 | 1.64 | 0.245 | 1.23 | 212.2 | 10.83 | 1.46 | 1 | 1 | 0.02 |
| 2026-08-04 | confirmed deaths | confirmed_deaths | 3 | 168.08 | 3.46 | 0.272 | 3.77 | 164.58 | 0.37 | 3.13 | 1 | 1 | -0.08 |
| 2026-08-07 | confirmed deaths | confirmed_deaths | 3 | 165.34 | 2.49 | 0.324 | 2.81 | 156.28 | 0 | 9.06 | 1 | 1 | -0.17 |
| 2026-08-11 | confirmed deaths | confirmed_deaths | 2 | 115.74 | 6.01 | 0.319 | 7.27 | 102.89 | 0 | 12.85 | 1 | 1 | -0.22 |
| 2026-08-15 | confirmed deaths | confirmed_deaths | 2 | 109.25 | 2.62 | 0.319 | 3.69 | 100.91 | 0 | 8.35 | 1 | 1 | -0.18 |
| 2026-08-17 | confirmed deaths | confirmed_deaths | 2 | 108.79 | 2.18 | 0.271 | 2.4 | 106.19 | 0 | 2.6 | 1 | 1 | -0.12 |
| 2026-08-22 | confirmed deaths | confirmed_deaths | 1 | 80.43 | 3.08 | 0.487 | 5.52 | 59.71 | 0 | 20.72 | 1 | 1 | -0.26 |
| 2026-08-24 | confirmed deaths | confirmed_deaths | 1 | 94.67 | 8.56 | 0.566 | 15.27 | 46.89 | 0 | 47.78 | 1 | 1 | -0.49 |
| 2026-07-23 | isolation beds | treatment | 4 | 55.73 | 1.89 | 0.089 | 2.3 | 47.79 | 0 | 7.94 | 1 | 1 | -0.21 |
| 2026-07-25 | isolation beds | treatment | 4 | 51.68 | 1.4 | 0.082 | 1.68 | 47.14 | 0.01 | 4.52 | 1 | 1 | -0.14 |
| 2026-07-26 | isolation beds | treatment | 4 | 61.67 | 1.64 | 0.1 | 1.98 | 54.82 | 0.23 | 6.63 | 1 | 1 | -0.12 |
| 2026-07-27 | isolation beds | treatment | 4 | 59.47 | 2.05 | 0.098 | 2.47 | 54.67 | 0 | 4.8 | 1 | 1 | -0.18 |
| 2026-07-31 | isolation beds | treatment | 4 | 82.76 | 1.95 | 0.132 | 2.43 | 68.17 | 0 | 14.6 | 0.75 | 1 | -0.28 |
| 2026-08-01 | isolation beds | treatment | 4 | 88.62 | 1.38 | 0.153 | 1.78 | 73.94 | 0 | 14.67 | 0.75 | 1 | -0.25 |
| 2026-08-02 | isolation beds | treatment | 4 | 84.37 | 1.65 | 0.141 | 2.03 | 69.7 | 0 | 14.67 | 0.75 | 1 | -0.28 |
| 2026-08-03 | isolation beds | treatment | 4 | 79.93 | 1.88 | 0.134 | 2.32 | 71.91 | 0 | 8.01 | 1 | 1 | -0.21 |
| 2026-08-04 | isolation beds | treatment | 3 | 86.12 | 1.66 | 0.142 | 1.81 | 84.49 | 1.63 | 0 | 1 | 1 | 0.05 |
| 2026-08-07 | isolation beds | treatment | 3 | 107.68 | 1.87 | 0.174 | 2.35 | 83.53 | 0 | 24.15 | 1 | 1 | -0.33 |
| 2026-08-11 | isolation beds | treatment | 2 | 91.8 | 0.69 | 0.156 | 0.78 | 86.73 | 0 | 5.07 | 1 | 1 | -0.13 |
| 2026-08-15 | isolation beds | treatment | 2 | 120.51 | 1.56 | 0.19 | 1.95 | 94.26 | 0 | 26.26 | 1 | 1 | -0.32 |
| 2026-08-17 | isolation beds | treatment | 2 | 98.59 | 2.09 | 0.158 | 2.54 | 90.9 | 0 | 7.69 | 1 | 1 | -0.18 |
| 2026-08-22 | isolation beds | treatment | 1 | 157.32 | 3.43 | 0.241 | 4.45 | 103.02 | 0 | 54.3 | 1 | 1 | -0.37 |
| 2026-08-24 | isolation beds | treatment | 1 | 134.63 | 5.57 | 0.216 | 7.18 | 100.67 | 0 | 33.96 | 1 | 1 | -0.36 |
| 2026-08-02 | onset reports | onsets | 1 | 81.88 | 0.77 | 0.283 | 0.63 | 63.32 | 0 | 18.56 | 1 | 1 | -0.28 |
| 2026-08-11 | onset reports | onsets | 1 | 141.47 | 0.59 | 0.61 | 0.37 | 67.8 | 0 | 73.67 | 0 | 1 | -0.55 |
| 2026-08-15 | onset reports | onsets | 1 | 108.06 | 1.35 | 0.183 | 1 | 100.93 | 7.13 | 0 | 1 | 1 | 0.17 |


```@raw html
</details>
```

## Saving forecast results

The one-week-back validation forecast, in the same archive format as the
release forecast, so the frozen "last week versus now" forecast is recorded
as a release asset alongside the forecast it is scored against.

```@raw html
<details><summary>Write forecast outputs</summary>
```

````julia
output_dir = get(
    ENV, "BVD_OUTPUT_DIR",
    joinpath(pkgdir(BVDOutbreakSize), "output")
)
mkpath(output_dir)
CSV.write(
    joinpath(output_dir, "forecast_validation.csv"),
    forecast_archive(
        [(7, validation_forecast)];
        made_date = frozen_lastweek.o.cutoff, thin = 5
    )
)
````

````
"/home/seabbs/code/seabbs/BVDOutbreakSize/worktrees/reporting/output/forecast_validation.csv"
````

```@raw html
</details>
```

