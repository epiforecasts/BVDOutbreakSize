# # Forecast evaluation
#
# How the forecasts on the [forecasts](@ref "Forecasts") page have scored against the data that arrived afterwards.
# How the outbreak-size and reproduction-number estimates moved from release to release is on the [in-sample checks](@ref "Estimate evolution across releases") page.
# Scoring is the continuous ranked probability score against a persistence baseline, defined in the [forecast scoring](@ref "Forecast scoring against a persistence baseline") Methods section.
# The same scoring by province is on the [province forecast evaluation](@ref "Province forecast evaluation") page and by health zone on the [health-zone forecast evaluation](@ref "Health-zone forecast evaluation") page.

#md # ```@raw html
#md # <details><summary>Load packages, data and fitted chains</summary>
#md # ```

## Shared setup: packages, observations and the fit registry. See
## `docs/pages/_setup.jl`.
using BVDOutbreakSize
include(joinpath(pkgdir(BVDOutbreakSize), "docs", "pages", "_setup.jl"))
#-
## The fits this page reads, loaded from the cache here.
chn_joint = load_fit("joint")
frozen_lastweek = load_fit("frozen_validation")
frozen_lastweek_streams = frozen_validation_stream_fits();

#md # ```@raw html
#md # </details>
#md # ```

# ## Summary
#
# The overall bullets come first, then a table of each stream's scores, from the sections below.
# Relative skill is the model's CRPS over the persistence baseline's, so a value below one beats the baseline.
# Coverage is the fraction of forecasts whose observed value falls inside the 90% predictive interval, nominally 0.9.

#md # ```@eval
#md # using Markdown, BVDOutbreakSize
#md # dir = joinpath(pkgdir(BVDOutbreakSize), "docs", "src", "summary_assets")
#md # Markdown.parse(read(joinpath(dir, "evaluation_forecast_national.md"), String))
#md # ```

# ## Forecast validation
#
# How last week's forecast held up against the data since observed, from the frozen refit in [forecast-versus-frozen evaluation](@ref "Forecast-versus-frozen evaluation").
# Only the streams the situation reports still update are validated, since a stopped stream repeats its last total and leaves nothing to score.
# The projected bed occupancy is scored against the beds held a week later, but this check is weak at a one-week-back freeze.
# The reported occupancy rate starts only on 9 June, so the capacity has no implied-capacity anchor and rides its random walk back to the freeze date.

#md # ```@raw html
#md # <details><summary>Fit one week back and validate the one-week-ahead forecast</summary>
#md # ```

## frozen_lastweek and frozen_lastweek_streams are computed in the setup
## block above, and `validation_forecast_from` is defined there.
validation_forecast = validation_forecast_from("frozen_validation");

## Each frozen individual (single-stream) fit's own one-week-ahead new-count
## forecast at the same cut-off as `frozen_lastweek`, from
## [`forecast_stream`](@ref) (the same per-stream forecaster
## `stream_forecasts.csv` uses), so the validation plots below can show the
## individual fit alongside the joint rather than the joint alone. Recovered
## has no individual fit and is absent here, as it is throughout this report.
## Only the still-reported streams are fitted at the validation cut-off, so
## a stream the situation reports have stopped updating is absent from
## `frozen_lastweek_streams` and carries no individual series here.
function _validation_individual_new(sid, stream::Symbol)
    haskey(frozen_lastweek_streams, sid) || return nothing
    return Float64.(
        forecast_stream(
            fit_forecast("frozen_validation_$sid"), stream; horizon = 7
        )
    )
end
validation_individual = NamedTuple(
    k => v
        for (k, v) in pairs(
            (;
                cases_new = _validation_individual_new(
                    "cases", :reported_cases
                ),
                deaths_new = _validation_individual_new(
                    "deaths", :suspected_deaths
                ),
                confirmed_new = _validation_individual_new(
                    "confirmed", :confirmed_cases
                ),
                confirmed_deaths_new = _validation_individual_new(
                    "confirmed_deaths", :confirmed_deaths
                ),
            )
        )
        if !isnothing(v)
)
## The frozen individual (treatment-only) fit's own bed-occupancy forecast.
## `nothing` when the beds have stopped being reported, so the treatment fit
## is absent; the bed panel then draws the joint alone.
validation_individual_isolation = haskey(frozen_lastweek_streams, "treatment") ?
    Float64.(
        forecast_stream(
            fit_forecast("frozen_validation_treatment"), :isolation_beds;
            horizon = 7
        )
    ) : nothing

## The observed beds at the current cut-off (the forecast target), so the
## frozen-fit bed forecast is scored against what the beds actually held.
## Held back once the beds stop being reported, since the last count would
## then be carried forward rather than observed at the target date.
_obs_beds = stream_reporting(obs, :isolation_beds) ?
    obs.isolation_history.counts[end] : missing
## Same observed/baseline keying as the plot below, so the table covers every
## fitted count stream (cumulative and new-count rows) plus the bed level.
## A harmonisation-break day between the frozen cut-off and the current one
## puts records into the confirmed cumulative that were never notified in that
## week, so the new-count truth carries a step the forecast was never
## predicting. Take it out, the same correction `score_releases.jl` applies.
## Grid days are relative to a seeding date fixed by the genetic tmrca, so the
## frozen fit's own `n` and the current `obs.n` index the same grid.
validation_breaks = (
    confirmed_cum = confirmed_break_correction(
        obs, frozen_lastweek.o.n, obs.n
    ),
    confirmed_deaths_cum = confirmed_break_correction(
        obs, frozen_lastweek.o.n, obs.n; deaths = true
    ),
)

## Observed cumulative at the target date per stream, keyed by the forecast's
## cumulative column; `baseline` is each stream's origin cumulative (the
## frozen cut-off), so the new count is scored against observed minus origin,
## less any harmonisation the window carries (see `validation_breaks`). Both
## the table and the plot below take the still-reported streams
## (`reporting_cum_cols`, from the setup block): a stream the situation
## reports have stopped updating has an origin and a target reading the same
## repeated total, so its cumulative truth is stale and its new-count truth is
## a guaranteed zero.
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

#md # ```@raw html
#md # </details>
#md # ```

#md # ```@raw html
#md # <details><summary>Forecast-versus-observed validation table</summary>
#md # ```

## `MarkdownTable` rather than a bare table expression: a DataFrame is #src
## `text/html`-showable, Literate prefers that mime, and the `@raw html` #src
## block it writes crosses Documenter's raw-block regex limit once the #src
## table grows. `MarkdownTable` is markdown-showable and not #src
## html-showable, so the table goes out as an ordinary markdown table #src
## rather than a fixed-width block of printed output. See its docstring #src
## for the mechanism. The same treatment is applied to every DataFrame #src
## display in the report pages. #src
#md # ```@raw html
#md # <figure class="bvd-figure">
#md # ```

MarkdownTable(validation_table) #hide

#md # ```@raw html
#md # <figcaption>
#md # Last week's frozen forecast against the data since observed, for the streams still reported.
#md # As in the scores below, the confirmed new-count rows take out any retrospective harmonisation step in the week, since it reattaches records notified earlier.
#md # The cumulative rows are scored against the published total, harmonisation included.
#md # </figcaption>
#md # </figure>
#md # </details>
#md # ```

#md # ```@raw html
#md # <details><summary>Forecast-versus-observed plot</summary>
#md # ```

validation_fig = plot_forecast_vs_truth(
    validation_forecast;
    observed = keep_streams(validation_observed, reporting_cum_cols),
    baseline = keep_streams(validation_baseline, reporting_cum_cols),
    breaks = validation_breaks,
    individual = keep_streams(validation_individual, reporting_cum_cols)
);

#md # ```@raw html
#md # </details>
#md # <figure class="bvd-figure">
#md # ```

validation_fig #hide

#md # ```@raw html
#md # <figcaption>
#md # The one-week-ahead forecast from the frozen fit for each still-reported count stream, cumulative and new, with the 90% predictive interval shaded and the count observed by the cut-off as a dashed black rule.
#md # A stream with its own single-stream fit carries that fit's forecast from the same cut-off as a dotted step outline.
#md # </figcaption>
#md # </figure>
#md # ```

#md # ```@raw html
#md # <details><summary>Bed forecast-versus-observed plot</summary>
#md # ```

validation_beds_fig = plot_forecast_beds_vs_truth(
    validation_forecast;
    isolation = _obs_beds, individual = validation_individual_isolation
);

#md # ```@raw html
#md # </details>
#md # <figure class="bvd-figure">
#md # ```

validation_beds_fig #hide

#md # ```@raw html
#md # <figcaption>
#md # Last week's projected bed occupancy against the beds occupied now (dashed), with the treatment-only fit's projection as a dotted step outline.
#md # </figcaption>
#md # </figure>
#md # ```

#md # ```@raw html
#md # <details><summary>Forecast-versus-now latent plot</summary>
#md # ```

## Current fit's draws of the new latent counts over the past week, the last
## seven days of each cumulative-trajectory deterministic.
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

#md # ```@raw html
#md # </details>
#md # <figure class="bvd-figure">
#md # ```

validation_latent_fig #hide

#md # ```@raw html
#md # <figcaption>
#md # The frozen fit's forecast of the past week's new infections, onsets and deaths against the current fit's estimate for the same window.
#md # These quantities are not observed, so they are scored distribution against distribution.
#md # </figcaption>
#md # </figure>
#md # ```

# ### Streams no longer reported

#md # ```@raw html
#md # <details><summary>Forecast for the streams no longer reported</summary>
#md # ```

## The last-reported date per stopped stream, and the frozen fit's own
## projection for them. `plot_forecast` draws a panel per new-count column
## the frame carries, so passing the stopped streams' columns alone gives the
## projection without the fabricated truth rule the validation figure would
## otherwise draw against a repeated total.
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

#md # ```@raw html
#md # </details>
#md # ```

## See the comment above `validation_table`'s display for why this wraps #src
## the table in `MarkdownTable` instead of showing it directly. #src
#md # ```@raw html
#md # <figure class="bvd-figure">
#md # ```

MarkdownTable(validation_stopped_streams) #hide

#md # ```@raw html
#md # <figcaption>
#md # The streams the situation reports no longer update, with the date each was last reported.
#md # </figcaption>
#md # </figure>
#md # ```

#-

#md # ```@raw html
#md # <figure class="bvd-figure">
#md # ```

validation_stopped_fig #hide

#md # ```@raw html
#md # <figcaption>
#md # The frozen fit's projection for these streams over the same week, with no observed rule, since their counts have not moved since they stopped.
#md # </figcaption>
#md # </figure>
#md # ```

# ## Forecast scoring across releases
#
# Every release's saved one- to four-week-ahead forecast is scored against the data observed since, against a persistence baseline and, where one exists, the stream's own individual fit.
# The tables in this section are the joint model's, one row per stream.
# Recovered has no individual fit of its own, so it is compared with the baseline only.
# Reported cases and suspected deaths stopped being updated partway through the outbreak, and the exports detection series is anchored to an earlier cut-off.
# Exports therefore has no scored forecast, and reported cases and suspected deaths each rest on a single matched forecast.
#
# Only a minority of the daily releases examined contribute a row, each a reconstruction of an earlier model version rather than the current fit.
# Only the newest few releases carry the current model's individual-stream forecasts.
# Each row shows the number of matched forecasts it rests on, from one to a handful.
#
# Four things are excluded from the scores, each for a stated reason rather than for scoring badly:
#
# - One whole reconstruction (`results-v1.6.0`), whose chain forecasts a near-zero median at every horizon and stream with an occasional five- or six-digit upper tail, the signature of a chain that failed to sample.
# - In the frozen section only, the confirmed-death rows of the fourteen frozen reconstructions cut between 16 July and 15 August 2026.
#   Their forecaster floored that stream at zero, giving a one-week median of zero against an observed 250 to 370.
# - An onset window containing a vintage whose reread total falls, since its increment is not what the situation reports added.
# - A stream that carries no persistence baseline, which is what makes a window scoreable at all.
#
# Nothing is dropped from the archive itself, and the score and overlay files in the data directory record everything that was scored.
#
# The symptom-onset stream is scored on the new reported count each vintage adds, because every vintage rereads the whole figure.
# Its printed total therefore moves with the scan error as well as with late reporting, and on fourteen vintages the reread total falls.
# A window containing a falling vintage is left unscored, as the [province scores](@ref "Forecast by province across releases") do for a window holding a harmonisation-break day.
# This bites hardest at the longer horizons: the frozen onset row keeps three of its twenty-nine windows, all at one week, and the cross-release row six of thirty-eight.
# Read the onset row's skill against the baseline rather than its coverage, and as resting on a handful of windows.

#md # ```@raw html
#md # <details><summary>Load and summarise the cross-release forecast scores</summary>
#md # ```

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
## The digitised onset triangle's own per-vintage total, as calendar dates,
## and the windows it makes unscoreable. A vintage that rereads the figure
## lower gives a window an increment the reports did not add, so the window
## is dropped from the scores and the overlay alike (see
## `drop_rescanned_onset_windows`). The current triangle is used to judge
## every release, since the scored truth is read off this one series.
_onset_vintage_dates = grid_date.(obs.onset_report_history.days)
_onset_vintage_totals = obs.onset_report_history.counts
_drop_rescanned(tbl) = drop_rescanned_onset_windows(
    tbl; vintage_dates = _onset_vintage_dates,
    vintage_totals = _onset_vintage_totals
)
forecast_scores_df = _drop_rescanned(forecast_scores_df)
forecast_overlay_df = _drop_rescanned(forecast_overlay_df)

## One row per (stream, fit) pooled over every horizon and release. The
## by-horizon and by-release detail tables carry the same columns at a finer
## grain (see src/scoring.jl). Every fit is kept here, since the
## relative-skill figure below compares the roles against each other. The
## tables rendered in this section select the joint role, and the individual
## fits are tabulated in their own section.
forecast_score_overview_table = forecast_score_overview(forecast_scores_df)
forecast_score_by_horizon_table = forecast_score_by_horizon(forecast_scores_df)
forecast_score_by_release_table = forecast_score_by_release(forecast_scores_df)

joint_score_overview_table = select_fit_role(
    forecast_score_overview_table, "joint"
)
joint_score_by_horizon_table = select_fit_role(
    forecast_score_by_horizon_table, "joint"
)
## The trailing `;` on this last assignment matters: without it, this whole
## setup chunk's last statement (the DataFrame it assigns) is Literate's
## implicitly displayed "result" for the chunk, on top of the deliberate
## display further down -- and a bare DataFrame is html-showable, so it
## goes out as a second, undisplayed-in-source `@raw html` block that (for
## a table this size) can itself hit the PCRE limit described above.
joint_score_by_release_table = select_fit_role(
    forecast_score_by_release_table, "joint"
);

#md # ```@raw html
#md # </details>
#md # <figure class="bvd-figure">
#md # ```

MarkdownTable(joint_score_overview_table) #hide

#md # ```@raw html
#md # <figcaption>
#md # Scores pooled over every horizon and release, one row per stream: the mean CRPS and its decomposition, coverage, bias, and relative skill against the persistence baseline on the natural and log scales.
#md # Each row also carries relative skill against the stream's own individual fit where one exists.
#md # </figcaption>
#md # </figure>
#md # ```

forecast_relative_skill_fig = plot_forecast_relative_skill(
    forecast_score_by_horizon_table
);

#md # ```@raw html
#md # <figure class="bvd-figure">
#md # ```

forecast_relative_skill_fig #hide

#md # ```@raw html
#md # <figcaption>
#md # Relative skill against the baseline by horizon, one panel per stream and one series per fit role, on a log scale with the reference line at one.
#md # </figcaption>
#md # </figure>
#md # ```

forecast_crps_by_horizon_fig = plot_forecast_crps_by_horizon(
    forecast_score_by_horizon_table
);

#md # ```@raw html
#md # <figure class="bvd-figure">
#md # ```

forecast_crps_by_horizon_fig #hide

#md # ```@raw html
#md # <figcaption>
#md # The mean CRPS by horizon split into its width, overprediction and underprediction, one stacked bar per horizon and fit role.
#md # </figcaption>
#md # </figure>
#md # <details><summary>Scores by horizon</summary>
#md # ```

MarkdownTable(joint_score_by_horizon_table) #hide

#md # ```@raw html
#md # </details>
#md # ```

forecast_skill_by_cutoff_fig = plot_forecast_skill_by_cutoff(
    forecast_score_by_release_table;
    title = "Relative skill against the baseline, by release"
);

#md # ```@raw html
#md # <figure class="bvd-figure">
#md # ```

forecast_skill_by_cutoff_fig #hide

#md # ```@raw html
#md # <figcaption>
#md # Relative skill against the baseline release by release, so a run of releases that lost to the baseline reads as a run rather than as an average.
#md # </figcaption>
#md # </figure>
#md # <details><summary>Scores by release</summary>
#md # ```

MarkdownTable(joint_score_by_release_table) #hide

#md # ```@raw html
#md # </details>
#md # ```

#md # ```@raw html
#md # <details><summary>Forecasts-versus-now overlay</summary>
#md # ```

forecast_overlay_fig = plot_forecast_overlay(
    scored_overlay(forecast_overlay_df)
);

#md # ```@raw html
#md # </details>
#md # <figure class="bvd-figure">
#md # ```

forecast_overlay_fig #hide

#md # ```@raw html
#md # <figcaption>
#md # Forecasts made at each release against the value observed since (black), one panel per stream and horizon by the date each forecast was made, with the median and 90% interval coloured by fit role.
#md # Each panel's axis is cropped to a small multiple of what the stream reached, and an interval or median past it is clamped and marked with an open triangle.
#md # </figcaption>
#md # </figure>
#md # ```

# ## Frozen-fit forecast evaluation
#
# The current model, frozen at earlier data cut-offs (see [forecast-versus-frozen evaluation](@ref "Forecast-versus-frozen evaluation")), is scored as the cross-release forecasts above, against the same baseline.
# The tables in this section are the frozen joint model's, one row per stream.
# Each still-reported stream's own frozen fit is scored at the one-week-back cut-off, and confirmed cases and confirmed deaths also at each May cut-off.
# These single-stream fits appear in the skill figures rather than the tables.
# The May cut-offs predate the first reported bed occupancy and recoveries, so those windows are left unscored.
# The baseline snapshot was taken weeks after the frozen cut-off, so it can hold later revisions to earlier days (see [forecast scoring against a persistence baseline](@ref "Forecast scoring against a persistence baseline")).

#md # ```@raw html
#md # <details><summary>Load and summarise the frozen-fit forecast scores</summary>
#md # ```

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
## Rows a superseded forecaster produced are dropped before anything is
## summarised or drawn, from the scores and the overlay alike, so the
## tables and the figures rest on one set of rows. See
## `drop_superseded_forecasts` for the one exclusion in force and why.
frozen_scores_df = _drop_rescanned(
    drop_superseded_forecasts(frozen_scores_df)
)
frozen_overlay_df = _drop_rescanned(
    drop_superseded_forecasts(frozen_overlay_df)
)

## The frozen joint carries `FROZEN_FIT`, so it is named as the joint role
## here and compared against each stream's own frozen fit. A release
## published before the archive had a `fit` column carries joint rows only.
frozen_score_overview_table = forecast_score_overview(
    frozen_scores_df; joint_fit = FROZEN_FIT
)
frozen_score_by_horizon_table = forecast_score_by_horizon(
    frozen_scores_df; joint_fit = FROZEN_FIT
)
frozen_score_by_release_table = forecast_score_by_release(
    frozen_scores_df; joint_fit = FROZEN_FIT
)

## The tables show the frozen joint alone, as the cross-release tables show
## the joint alone, so a row reads as one model at one cut-off rather than a
## stream interleaving two fits. The archive's single-stream frozen fits stay
## in the scored data and in the figures, which compare the roles against
## each other. Selecting the joint role leaves `fit` single-valued, so it is
## dropped and the model named in the prose instead.
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

## One row per release for the cut-offs more than one release forecast.
## See the comment above `joint_score_by_release_table`'s assignment for why
## this setup chunk's last statement needs a trailing `;`.
frozen_score_by_vintage_table = forecast_score_by_vintage(
    frozen_scores_df; joint_fit = FROZEN_FIT
)
frozen_score_by_vintage_display = _frozen_joint_only(
    frozen_score_by_vintage_table
);

#md # ```@raw html
#md # </details>
#md # ```

MarkdownTable(frozen_score_overview_display) #hide

#-

frozen_relative_skill_fig = plot_forecast_relative_skill(
    frozen_score_by_horizon_table
);

#md # ```@raw html
#md # <figure class="bvd-figure">
#md # ```

frozen_relative_skill_fig #hide

#md # ```@raw html
#md # <figcaption>
#md # Relative skill against the baseline by horizon for the frozen cut-offs.
#md # </figcaption>
#md # </figure>
#md # ```

frozen_crps_by_horizon_fig = plot_forecast_crps_by_horizon(
    frozen_score_by_horizon_table;
    title = "CRPS decomposition by horizon, frozen cut-offs"
);

#md # ```@raw html
#md # <figure class="bvd-figure">
#md # ```

frozen_crps_by_horizon_fig #hide

#md # ```@raw html
#md # <figcaption>
#md # The mean CRPS by horizon for the frozen cut-offs, split as in the cross-release section above.
#md # </figcaption>
#md # </figure>
#md # <details><summary>Scores by horizon</summary>
#md # ```

MarkdownTable(frozen_score_by_horizon_display) #hide

#md # ```@raw html
#md # </details>
#md # ```

frozen_skill_by_cutoff_fig = plot_forecast_skill_by_cutoff(
    frozen_score_by_release_table;
    xlabel = "Frozen cut-off",
    title = "Relative skill against the baseline, by frozen cut-off"
);

#md # ```@raw html
#md # <figure class="bvd-figure">
#md # ```

frozen_skill_by_cutoff_fig #hide

#md # ```@raw html
#md # <figcaption>
#md # Relative skill against the baseline cut-off by cut-off, pooled over the horizons each cut-off forecast.
#md # </figcaption>
#md # </figure>
#md # <details><summary>Scores by frozen cut-off</summary>
#md # ```

MarkdownTable(frozen_score_by_release_display) #hide

#md # ```@raw html
#md # </details>
#md # ```

# ### Frozen skill by release

#md # ```@raw html
#md # <details><summary>Frozen skill per release</summary>
#md # ```

frozen_skill_by_vintage_fig = plot_forecast_skill_by_vintage(
    frozen_score_by_vintage_table
);

#md # ```@raw html
#md # </details>
#md # <figure class="bvd-figure">
#md # ```

frozen_skill_by_vintage_fig #hide

#md # ```@raw html
#md # <figcaption>
#md # Skill at each cut-off that more than one release forecast, one point per release rather than pooled.
#md # Releases run in the order they were cut, evenly spaced rather than to calendar scale.
#md # </figcaption>
#md # </figure>
#md # <details><summary>Frozen scores by release</summary>
#md # ```

MarkdownTable(frozen_score_by_vintage_display) #hide

#md # ```@raw html
#md # </details>
#md # ```

#md # ```@raw html
#md # <details><summary>Frozen-fit forecasts-versus-now overlay</summary>
#md # ```

## The frozen joint and the persistence baseline only. A single-stream
## frozen fit sits at the same made dates as the joint, so drawing it here
## would overplot the joint point beside it. The single-stream frozen fits
## are compared against the joint in the skill figures and in the
## validation plot at the one-week-back cut-off.
frozen_overlay_fig = plot_forecast_overlay(
    scored_overlay(
        vcat(
            select_fit_role(frozen_overlay_df, "joint"),
            select_fit_role(frozen_overlay_df, "baseline")
        )
    )
);

#md # ```@raw html
#md # </details>
#md # <figure class="bvd-figure">
#md # ```

frozen_overlay_fig #hide

#md # ```@raw html
#md # <figcaption>
#md # The frozen forecasts made at each cut-off against the value observed since (black), one panel per stream and horizon by the cut-off each was made from.
#md # Each panel carries the frozen forecast and the persistence baseline, coloured as in the cross-release overlay above.
#md # </figcaption>
#md # </figure>
#md # ```

# The frozen re-fits below freeze the renewal data to an earlier cut-off and re-fit, so that a change driven by newer data can be distinguished from one driven by a change of method.
# The one-week-back re-fit uses the headline sampler settings, and the re-fits at earlier cut-offs use the single-population frozen settings given in the methods.

#md # ```@raw html
#md # <details><summary>Freeze the renewal data to a cut-off and re-fit</summary>
#md # ```

## Frozen re-fits and released_df are prepared in the setup block above.

#md # ```@raw html
#md # </details>
#md # ```

# ## Individual fits against the baseline
#
# This section scores each stream's own individual fit across releases, as [forecast scoring across releases](@ref "Forecast scoring across releases") does the joint, against the same persistence baseline.

#md # ```@raw html
#md # <details><summary>Individual-fit rows of the cross-release scores</summary>
#md # ```

## The relative skill against a stream's individual fit is only ever
## computed on the joint model's row, so on these rows it is missing by
## construction and the column is dropped rather than shown empty.
individual_score_overview_table = drop_individual_fit_columns(
    select_fit_role(forecast_score_overview_table, "individual")
)
individual_score_by_horizon_table = drop_individual_fit_columns(
    select_fit_role(forecast_score_by_horizon_table, "individual")
)
## See the comment above `joint_score_by_release_table`'s assignment for why
## this setup chunk's last statement needs a trailing `;`.
individual_score_by_release_table = drop_individual_fit_columns(
    select_fit_role(forecast_score_by_release_table, "individual")
);

#md # ```@raw html
#md # </details>
#md # ```

MarkdownTable(individual_score_overview_table) #hide

#-

individual_relative_skill_fig = plot_forecast_relative_skill(
    individual_score_by_horizon_table;
    empty_message = "Empty: no release old enough for its targets to " *
        "have been observed carries an individual-stream " *
        "forecast. Not a missing forecast."
);

#md # ```@raw html
#md # <figure class="bvd-figure">
#md # ```

individual_relative_skill_fig #hide

#md # ```@raw html
#md # <figcaption>
#md # Relative skill against the baseline by horizon for each stream's own individual fit, one panel per stream.
#md # </figcaption>
#md # </figure>
#md # <details><summary>Scores by horizon</summary>
#md # ```

MarkdownTable(individual_score_by_horizon_table) #hide

#md # ```@raw html
#md # </details>
#md # ```

#md # ```@raw html
#md # <details><summary>Scores by release</summary>
#md # ```

MarkdownTable(individual_score_by_release_table) #hide

#md # ```@raw html
#md # </details>
#md # ```

# ## Saving forecast outputs
#
# The one-week-back validation forecast is written in the release forecast's archive format, so it is recorded as a release asset alongside the forecast it is scored against.

#md # ```@raw html
#md # <details><summary>Write forecast outputs</summary>
#md # ```

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

#md # ```@raw html
#md # </details>
#md # ```

#md # ```@raw html
#md # <details><summary>Write the summary bullets</summary>
#md # ```

## The bullets under the summary heading at the top of the page. They read
## tables built further down, so they are written here and read back when
## the site is assembled.
evaluation_forecast_national_summary = let
    scored(tbl) = filter(r -> !ismissing(r.rel_to_baseline), tbl)
    function beat(label, tbl)
        rows = scored(tbl)
        size(rows, 1) == 0 &&
            return string("- **", label, ":** no scored forecasts yet.")
        return string(
            "- **", label, ":** beat the baseline on ",
            count(<(1), rows.rel_to_baseline), " of ", size(rows, 1),
            " streams."
        )
    end
    frozen_joint = select_fit_role(frozen_score_overview_table, "joint")
    overall = [
        beat("Joint model across releases", joint_score_overview_table),
        beat("Frozen joint model", frozen_joint),
    ]
    ## One row per stream and evaluation, the cross-release scores first.
    rows(lead, tbl) = let r = scored(tbl)
        DataFrame(
            "Stream" => r.stream,
            "Evaluation" => fill(lead, size(r, 1)),
            "Relative skill" => finite_round.(r.rel_to_baseline),
            "90% coverage" => finite_round.(r.coverage_90),
            "Bias" => finite_round.(r.bias),
            "Forecasts" => r.n,
        )
    end
    table = vcat(
        rows("Across releases", joint_score_overview_table),
        rows("Frozen fit", frozen_joint)
    )
    join(
        vcat(
            [join(overall, "\n")],
            size(table, 1) == 0 ? String[] :
                [markdown_table(table; missing_text = "n/a")]
        ), "\n\n"
    )
end
dashboard_dir = joinpath(
    pkgdir(BVDOutbreakSize), "docs", "src", "summary_assets"
)
mkpath(dashboard_dir)
write(
    joinpath(dashboard_dir, "evaluation_forecast_national.md"),
    evaluation_forecast_national_summary
);

#md # ```@raw html
#md # </details>
#md # ```
