# # Forecast evaluation
#
# How the forecasts on the [forecasts](@ref "Forecasts") page have scored against the data that arrived afterwards.
# It also tracks how the outbreak-size and reproduction-number estimates moved from release to release.
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
frozen_lastweek_streams = frozen_validation_stream_fits()
chn_exports = load_fit("exports")
chn_deaths = load_fit("deaths")
chn_cases = load_fit("cases")
chn_confirmed = load_fit("confirmed")
chn_confirmed_deaths = load_fit("confirmed_deaths")
chn_treatment = load_fit("treatment")
chn_onsets = load_fit("onsets")
frozen_by_cutoff = frozen_fits_by_cutoff()
frozen_C(c) = vec(Array(frozen_by_cutoff[c].chn[:C_T]))
## Every frozen fit is a full joint fit, so it carries the same walk base
## chn_joint does.
frozen_R0(c) = r0_walk_draws(frozen_by_cutoff[c].chn);

#md # ```@raw html
#md # </details>
#md # ```

# ## Summary
#
# The overall bullets come first, then each stream's own scores, from the sections further down this page.
# Relative skill is the model's CRPS over the persistence baseline's, so a value below one beats the baseline.
# Coverage is the fraction of forecasts whose observed value falls inside the 90% predictive interval, nominally 0.9.

#md # ```@eval
#md # using Markdown, BVDOutbreakSize
#md # dir = joinpath(pkgdir(BVDOutbreakSize), "docs", "src", "summary_assets")
#md # Markdown.parse(read(joinpath(dir, "evaluation_forecast_national.md"), String))
#md # ```

# ## Forecast validation
#
# How last week's forecast held up against the data since observed, using the frozen re-fit and one-week projection defined in [forecast-versus-frozen evaluation](@ref "Forecast-versus-frozen evaluation").
# Only the streams the situation reports are still updating are validated here.
# A stream that has stopped being reported carries a cumulative total that repeats its last reported value, so there is no observation for the past week to score against.
# The frozen fit also conditions on the isolation beds, so the projected bed occupancy is scored against the beds held a week later.
# The bed validation is weak at a one-week-back freeze.
# The reported occupancy rate starts only on 9 June, so the capacity has no implied-capacity anchor and rides its random walk back to the freeze date.
# Like the scores further down, the confirmed new-count rows here take out any retrospective harmonisation step the week contained.
# Such a step reattaches records notified earlier, so it is not something the forecast was predicting.
# The cumulative rows are scored against the published total, harmonisation included.

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
MarkdownTable(validation_table) #hide

#md # ```@raw html
#md # </details>
#md # ```

# The observation panels histogram the one-week-ahead forecast made from the frozen fit: a cumulative and a new-count panel for each still-reported count stream the forecast carries.
# The 90% predictive interval is shaded, and the count observed by the current cut-off is a dashed black rule.
# Where a stream has its own individual (single-stream) fit, that fit's forecast from the same frozen cut-off is overlaid as a dotted step outline on the joint's own histogram bins.

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
#md # ```

validation_fig #hide

# The bed panel scores last week's projected occupancy against the beds occupied now (the dashed rule), with the individual (treatment-only) fit's own projection overlaid as a dotted step outline on the joint's own histogram bins.

#md # ```@raw html
#md # <details><summary>Bed forecast-versus-observed plot</summary>
#md # ```

validation_beds_fig = plot_forecast_beds_vs_truth(
    validation_forecast;
    isolation = _obs_beds, individual = validation_individual_isolation
);

#md # ```@raw html
#md # </details>
#md # ```

validation_beds_fig #hide

# The latent quantities are not observed, so they are scored distribution against distribution: what the frozen fit forecast for the past week's new infections, onsets and deaths against what the current fit now estimates for the same window.

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
#md # ```

validation_latent_fig #hide

# ### Streams no longer reported
#
# The situation reports have stopped updating some of the streams the model fits, listed with the date each was last reported below.
# The panels show what the frozen fit projected for those streams over the same week, without an observed rule, since the count they would be scored against has not moved since the stream stopped.

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
MarkdownTable(validation_stopped_streams) #hide

#-

validation_stopped_fig #hide

# ## Forecast scoring across releases
#
# Every release's saved one- to four-week-ahead forecast is scored against the data observed since, against a persistence baseline and, where one exists, the stream's own individual fit as well as the joint.
# The tables in this section are the joint model's, one row per stream.
# See [forecast scoring against a persistence baseline](@ref "Forecast scoring against a persistence baseline") for how the scores, the relative skill and the baseline are built.
# Recovered has no individual fit of its own, so its comparison is the baseline against the joint only.
# Reported cases and suspected deaths stopped being updated by the situation reports partway through the outbreak, and exports' confirmed-detection series is anchored to an earlier cut-off.
# Exports therefore contributes no scored forecast, and reported cases and suspected deaths each rest on exactly one matched forecast, a single window rather than a settled sample.
#
# Only a minority of the daily releases examined contribute a row to the table below, each a reconstruction of an earlier model version rather than the current fit.
# Only the newest few releases carry the current model's own individual-stream forecasts, and the backfilled reconstructions carry none at all.
# Every row also rests on one to a handful of matched forecasts, shown as its own count rather than rounded away.
#
# Four things are excluded from the scores here and in the frozen section below, each for a stated reason rather than for scoring badly, and the second of them applies to the frozen section alone.
#
# - One whole reconstruction (`results-v1.6.0`): its chain forecasts a near-zero median at every horizon and stream, with the upper predictive tail occasionally reaching five- and six-digit values, which is the signature of a chain that failed to sample rather than a forecast.
# - Frozen section only: the confirmed-death rows of the fourteen frozen reconstructions cut between 16 July and 15 August 2026, whose forecaster could not project that stream from its own trajectory and floored it at zero, so each carries a one-week median of exactly zero against an observed 250 to 370. Reconstructions cut after that window project the stream normally.
# - An onset window containing a vintage whose reread total falls, since its increment is not what the situation reports added.
# - A stream that carries no persistence baseline, which is what makes a window scoreable at all.
#
# Nothing is dropped from the archive itself; `data/forecast_scores*.csv` and `data/forecast_overlay*.csv` record everything that was scored.
#
# The symptom-onset stream is scored on the new reported count each vintage adds rather than on its level, because every vintage rereads the whole figure.
# Its printed total therefore moves with the scan error as well as with late reporting.
# On fourteen vintages the reread total falls, which a cumulative onset curve cannot do, and on many others it repeats unchanged.
# The scored truth cannot absorb that, since it is the increment between the vintages at the two ends of a window.
# A window containing a falling vintage is therefore left unscored, the rule the [province scores](@ref "Forecast by province across releases") already apply to a window holding a harmonisation-break day.
# It bites hardest at the longer horizons, a four-week window being more likely to contain a reread than a one-week one: the frozen onset row keeps three of its twenty-nine windows, all at one week, and the cross-release row six of thirty-eight.
# Read the onset row's skill against the baseline rather than its coverage, and read it as resting on a handful of windows.

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
#md # ```

# The headline pools every horizon and release into one row per stream for the joint model: the mean CRPS and its decomposition, coverage, bias, and the relative skill against the persistence baseline, on both the natural and the log scale.
# Each row also carries relative skill against the stream's own individual fit where one exists.

MarkdownTable(joint_score_overview_table) #hide

# The same relative skill against the baseline, by horizon: one panel per stream, one series per fit role, on a log-scaled skill axis with the reference line at one.

forecast_relative_skill_fig = plot_forecast_relative_skill(
    forecast_score_by_horizon_table
);

forecast_relative_skill_fig #hide

# What that error is made of, by horizon: the mean CRPS split into its width, its overprediction and its underprediction, one stacked bar per horizon and fit role.

forecast_crps_by_horizon_fig = plot_forecast_crps_by_horizon(
    forecast_score_by_horizon_table
);

forecast_crps_by_horizon_fig #hide

#md # ```@raw html
#md # <details><summary>Scores by horizon</summary>
#md # ```

MarkdownTable(joint_score_by_horizon_table) #hide

#md # ```@raw html
#md # </details>
#md # ```

# The same relative skill against the baseline, release by release, so a run of releases that lost to the baseline reads as a run rather than as an average.

forecast_skill_by_cutoff_fig = plot_forecast_skill_by_cutoff(
    forecast_score_by_release_table;
    title = "Relative skill against the baseline, by release"
);

forecast_skill_by_cutoff_fig #hide

#md # ```@raw html
#md # <details><summary>Scores by release</summary>
#md # ```

MarkdownTable(joint_score_by_release_table) #hide

#md # ```@raw html
#md # </details>
#md # ```

# Forecasts made at each release against the value observed since, one panel per stream and horizon, the observed value in black.
# The median and 90% interval are coloured by fit role: the persistence baseline, the stream's individual fit and the joint.
# The x-axis is the date each forecast was made, so an incident stream's observed window pairs unambiguously with the forecast that made it.
# Each panel's axis is cropped to a small multiple of what that stream actually reached, so one very wide interval cannot squash every other series flat.
# An interval or median too wide for the panel is clamped at the top and marked with an open triangle rather than silently cut off.

#md # ```@raw html
#md # <details><summary>Forecasts-versus-now overlay</summary>
#md # ```

forecast_overlay_fig = plot_forecast_overlay(
    scored_overlay(forecast_overlay_df)
);

#md # ```@raw html
#md # </details>
#md # ```

forecast_overlay_fig #hide

# ## Frozen-fit forecast evaluation
#
# The current model, frozen at earlier data cut-offs (see [Forecast-versus-frozen evaluation](@ref "Forecast-versus-frozen evaluation")), is scored the same way as the cross-release forecasts above, against the same persistence baseline.
# The tables in this section are the frozen joint model's, one row per stream.
# Each stream's own frozen fit is also scored where one exists, at the one-week-back cut-off and for still-reported streams only, and is carried by the skill figures rather than by the tables.
# Every other cut-off carries the joint alone.
# The May cut-offs predate the first reported bed occupancy and the first reported recoveries, so those windows are left unscored rather than scored against a series that had not started.
# The baseline carries a weaker data-vintage guarantee than the cross-release one, since its snapshot was taken weeks after the frozen cut-off and can hold later revisions to earlier days (see [forecast scoring against a persistence baseline](@ref "Forecast scoring against a persistence baseline")).

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

# The same relative skill against the baseline, by horizon, for the frozen cut-offs.

frozen_relative_skill_fig = plot_forecast_relative_skill(
    frozen_score_by_horizon_table
);

frozen_relative_skill_fig #hide

# What that error is made of, by horizon, as in the cross-release section above.

frozen_crps_by_horizon_fig = plot_forecast_crps_by_horizon(
    frozen_score_by_horizon_table;
    title = "CRPS decomposition by horizon, frozen cut-offs"
);

frozen_crps_by_horizon_fig #hide

#md # ```@raw html
#md # <details><summary>Scores by horizon</summary>
#md # ```

MarkdownTable(frozen_score_by_horizon_display) #hide

#md # ```@raw html
#md # </details>
#md # ```

# The same relative skill against the baseline, cut-off by cut-off, pooled over the horizons each cut-off forecast.

frozen_skill_by_cutoff_fig = plot_forecast_skill_by_cutoff(
    frozen_score_by_release_table;
    xlabel = "Frozen cut-off",
    title = "Relative skill against the baseline, by frozen cut-off"
);

frozen_skill_by_cutoff_fig #hide

#md # ```@raw html
#md # <details><summary>Scores by frozen cut-off</summary>
#md # ```

MarkdownTable(frozen_score_by_release_display) #hide

#md # ```@raw html
#md # </details>
#md # ```

# ### Frozen skill by release
#
# Skill at each cut-off more than one release forecast, one point per release rather than pooled across releases.
# Releases run in the order they were cut, evenly spaced rather than to calendar scale.

#md # ```@raw html
#md # <details><summary>Frozen skill per release</summary>
#md # ```

frozen_skill_by_vintage_fig = plot_forecast_skill_by_vintage(
    frozen_score_by_vintage_table
);

#md # ```@raw html
#md # </details>
#md # ```

frozen_skill_by_vintage_fig #hide

#md # ```@raw html
#md # <details><summary>Frozen scores by release</summary>
#md # ```

MarkdownTable(frozen_score_by_vintage_display) #hide

#md # ```@raw html
#md # </details>
#md # ```

# The frozen forecasts made at each cut-off against the value observed since, one panel per stream and horizon, the observed value in black.
# Each panel carries the frozen forecast and the persistence baseline, coloured as in the cross-release overlay above, and the x-axis is the cut-off each forecast was made from.

#md # ```@raw html
#md # <details><summary>Frozen-fit forecasts-versus-now overlay</summary>
#md # ```

## The frozen joint and the persistence baseline only. A single-stream
## frozen fit exists at the one-week-back cut-off alone, so its series lands
## on one made date of a panel spanning every cut-off, overplotting the
## joint point it sits beside rather than reading as a second series. The
## single-stream frozen fits are compared against the joint in the skill
## figures and in the validation plot at that cut-off.
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
#md # ```

frozen_overlay_fig #hide

# The frozen re-fits below freeze the renewal data to an earlier cut-off and re-fit, so that a change driven by newer data can be distinguished from one driven by a change of method.
# Each uses the full headline settings: 1000 draws across two chains.

#md # ```@raw html
#md # <details><summary>Freeze the renewal data to a cut-off and re-fit</summary>
#md # ```

## Frozen re-fits and released_df are prepared in the setup block above.

#md # ```@raw html
#md # </details>
#md # ```

# ## Individual fits against the baseline
#
# This section carries the same cross-release forecast scoring as [Forecast scoring across releases](@ref "Forecast scoring across releases") above, for each stream's own individual fit rather than the joint, against the same persistence baseline.

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

# The same relative skill against the baseline, by horizon, one panel per stream (dataset), for each stream's own individual fit.

individual_relative_skill_fig = plot_forecast_relative_skill(
    individual_score_by_horizon_table;
    empty_message = "Empty: no release old enough for its targets to " *
        "have been observed carries an individual-stream " *
        "forecast. Not a missing forecast."
);

individual_relative_skill_fig #hide

#md # ```@raw html
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

# ## Estimate evolution across releases
#
# How the outbreak-size estimate has moved as situation reports accrued, three series on one calendar axis.
# The estimate published at each release is in blue, drawn as a median with nested 30/60/90% interval bars because each release is its own fit rather than one continuous model.
# The current model frozen at earlier cut-offs is in red.
# The current model on current data is the green band, drawn day by day so the latest estimate reads against the earlier points.
# Dotted vertical rules mark the release dates.
# The published series switches from a closed-form integral model to a renewal model on 7 June, so a step there can reflect the change of method rather than of data.

#md # ```@raw html
#md # <details><summary>Released estimates and the current-model frozen re-fits</summary>
#md # ```

## Released median and 30/60/90% intervals per release, from
## `data/released_estimates.csv`. Each tuple is
## `(date, median, lo30, hi30, lo60, hi60, lo90, hi90)`.
release_evolution = [
    (
        string(r.date), r.median, r.lo30, r.hi30, r.lo60, r.hi60,
        r.lo90, r.hi90,
    ) for r in eachrow(released_df)
]

## The current model frozen at earlier cut-offs, each its own discrete
## estimate: the McCabe-matched cut-offs (20, 23, 27 May), the 8 June Chamla
## confirmed-case anchor and the one-week-back validation fit
## (`frozen_lastweek`, at `validation_cutoff`). Each tuple carries the
## median and 30/60/90% credible bounds from the frozen draws; `round_fn`
## rounds to a whole count for outbreak size, and is passed through
## unrounded for a continuous quantity such as R0.
function _ci369(xs; round_fn = x -> round(Int, x))
    q(p) = round_fn(quantile(xs, p))
    return (q(0.5), q(0.35), q(0.65), q(0.2), q(0.8), q(0.05), q(0.95))
end
frozen_by_cutoff[validation_cutoff] = frozen_lastweek
## The cut-offs every frozen fit above was made at, shared by the
## outbreak-size and R0 by-release overlays below.
_frozen_matched_cutoffs = sort(
    union(
        frozen_cutoffs,
        [validation_cutoff, default_chamla_cutoff()]
    )
)
frozen_matched = [(c, _ci369(frozen_C(c))...) for c in _frozen_matched_cutoffs]

## The current-data, current-model estimate as the cumulative-infection
## trajectory over the day grid (one calendar date per grid day, day 1 is
## the seeding date), summarised by per-day 30/60/90% credible bounds. This
## is the same latent quantity the cumulative-trajectory figure shows, so
## the current estimate rises over time on the release-date axis instead of
## sitting flat. Drawn against calendar dates, it lines up with the
## release and frozen points.
infection_trajectory = let
    mat = chn_joint[:cumulative_infections]
    trajs = [collect(v) for v in vec(collect(mat))]
    ## Only over the comparison window — from the earliest release date to the
    ## cut-off — not back to the seeding date.
    start_day = obs.n - value(obs.cutoff - Date(release_evolution[1][1]))
    days = max(start_day, 1):obs.n
    dates = [obs.seeding + Day(d - 1) for d in days]
    q(d, p) = quantile(Float64[t[d] for t in trajs], p)
    (
        dates,
        [q(d, 0.35) for d in days], [q(d, 0.65) for d in days],
        [q(d, 0.2) for d in days], [q(d, 0.8) for d in days],
        [q(d, 0.05) for d in days], [q(d, 0.95) for d in days],
    )
end

evolution_fig = plot_estimate_evolution(
    release_evolution;
    renewal = frozen_matched,
    renewal_label = "Current model frozen at earlier cut-offs",
    trajectory = infection_trajectory,
    title = "Outbreak-size estimate as data accrued"
);

#md # ```@raw html
#md # </details>
#md # ```

evolution_fig #hide

# ### Reproduction number by release
#
# The reproduction number estimated at each release, the same kind of release-by-release picture as the outbreak-size evolution above.
# Each release's cut-off reproduction number $R_T$ is drawn as a discrete estimate, a median with nested 30/60/90% interval bars.
# The current fit's daily $R_t$ over its established window is drawn as the continuous band, and $R_t = 1$ is marked.
# The reproduction-number axis is fixed at three across this figure and the by-dataset one below, with an interval running past it clamped and marked with an open triangle.

#md # ```@raw html
#md # <details><summary>Reproduction number per release with the current-fit band</summary>
#md # ```

rt_release_df = CSV.read(
    joinpath(pkgdir(BVDOutbreakSize), "data", "rt_by_release.csv"), DataFrame
)
rt_release = [
    (
        string(r.date), r.median, r.lo30, r.hi30, r.lo60, r.hi60,
        r.lo90, r.hi90,
    ) for r in eachrow(rt_release_df)
]

## The current fit's daily Rt over its established window, summarised per day
## into a 30/60/90% band, reusing the same walk reconstruction the Rt figure
## uses so the band lines up with the per-release points on the calendar axis.
## The band is drawn only from the first release date onward, so it spans the
## same window as the per-release estimates rather than extending back to the
## renewal start. The first release day is the earliest date in
## `rt_by_release.csv` as a grid day; the walk is still reconstructed from the
## renewal start `_rt_start_plot` (the model knot grid) and the window is
## clamped into the reconstructed range so the quantiles never hit masked days.
rt_release_trajectory = let
    rt_walk_start = clamp(_BREAKPOINT - RT_WALK_LEAD, _rt_start_plot, obs.n)
    mat = reconstruct_rt(
        chn_joint; n = obs.n, breakpoint = _BREAKPOINT,
        rt_start = _rt_start_plot, rt_walk_start = rt_walk_start,
        ramp = RT_INTERVENTION_RAMP
    )
    first_release_day = clamp(
        value(minimum(rt_release_df.date) - obs.seeding) + 1,
        _rt_start_plot, obs.n
    )
    days = first_release_day:obs.n
    dates = [obs.seeding + Day(d - 1) for d in days]
    q(d, p) = quantile(collect(skipmissing(@view mat[:, d])), p)
    (
        dates,
        [q(d, 0.35) for d in days], [q(d, 0.65) for d in days],
        [q(d, 0.2) for d in days], [q(d, 0.8) for d in days],
        [q(d, 0.05) for d in days], [q(d, 0.95) for d in days],
    )
end

## One fixed reproduction-number axis across both figures below. The
## estimates sit around one and the widest stream's 90% upper pulls a free
## axis past four, which flattens every panel onto the lower quarter of its
## range. An interval past the crop is clamped and marked, so nothing is
## silently cut. The basic reproduction number keeps its own axis: it sits
## around two with tails near four, and the same crop would clip the
## estimates themselves.
const _RT_AXIS_MAX = 3.0

rt_evolution_fig = plot_estimate_evolution(
    rt_release;
    trajectory = rt_release_trajectory,
    ylabel = "Reproduction number",
    title = "Reproduction number as data accrued",
    released_label = "Released estimate (per project release)",
    trajectory_label = "Current model, current data",
    refline = 1.0,
    ymax = _RT_AXIS_MAX
);

#md # ```@raw html
#md # </details>
#md # ```

rt_evolution_fig #hide

# ### Reproduction number by release and dataset
#
# The same release-by-release reproduction number split into one panel per dataset, so each dataset's history reads against the others and against the joint.
# Panels share a calendar axis and the fixed reproduction-number range, and $R_t = 1$ is marked.
# Each release's cut-off value is a median with nested 30/60/90% interval bars.
# A dataset the report also fits on its own carries that fit's current-model band behind its points, built as in the overview above.
# Confirmed deaths carries no band, so its panel shows release points alone.
# Only the most recent releases published these per-dataset estimates, so every panel spans a much shorter window than the overview above rather than a different history.

#md # ```@raw html
#md # <details><summary>Reproduction number per release by fit</summary>
#md # ```

## Schema of the per-release, per-fit estimate tables written by
## scripts/score_releases.jl from each release's stream_estimates.csv.
_by_stream_schema = (;
    release = String, date = Date, fit = String,
    median = Float64, lo30 = Float64, hi30 = Float64, lo60 = Float64,
    hi60 = Float64, lo90 = Float64, hi90 = Float64,
)

## Fits in a fixed order, the joint first, so the panels do not reshuffle
## between builds. Labels match the per-stream table on the in-sample page.
## Recovered is absent because it has no individual fit.
_fit_order = [
    "joint", "cases", "deaths", "confirmed", "confirmed_deaths",
    "treatment", "onsets", "exports",
]
_fit_labels = Dict(
    "joint" => "joint", "cases" => "cases (DRC)",
    "deaths" => "deaths (DRC)", "confirmed" => "confirmed (DRC)",
    "confirmed_deaths" => "confirmed deaths (DRC)",
    "treatment" => "isolation (DRC)", "onsets" => "onsets (DRC)",
    "exports" => "exports"
)

## Group a per-fit estimate table into the label => tuples pairs the faceted
## plot takes, keyed on the date so the mixed release tag shapes
## (`results-v1.9.0` and `results-1243`) never reach the axis.
function _fit_groups(df)
    return [
        get(_fit_labels, f, f) =>
            [
            (
                string(r.date), r.median, r.lo30, r.hi30, r.lo60, r.hi60,
                r.lo90, r.hi90,
            ) for r in eachrow(df) if r.fit == f
        ]
            for f in _fit_order
    ]
end

## Per-fit reproduction-number trajectory, reconstructing the walk exactly as
## `plot_rt_streams` does per stream.
function _stream_rt_trajectory(chn, dates; rt_start, rt_walk_start)
    mat = reconstruct_rt(
        chn; n = obs.n, breakpoint = _BREAKPOINT,
        rt_start = rt_start, rt_walk_start = rt_walk_start,
        ramp = RT_INTERVENTION_RAMP
    )
    first_date = isempty(dates) ? obs.seeding : minimum(dates)
    first_day = clamp(value(first_date - obs.seeding) + 1, rt_start, obs.n)
    days = first_day:obs.n
    ds = [obs.seeding + Day(d - 1) for d in days]
    q(d, p) = quantile(collect(skipmissing(@view mat[:, d])), p)
    return (
        ds,
        [q(d, 0.35) for d in days], [q(d, 0.65) for d in days],
        [q(d, 0.2) for d in days], [q(d, 0.8) for d in days],
        [q(d, 0.05) for d in days], [q(d, 0.95) for d in days],
    )
end

## The single-stream chains and their renewal-walk starts, keyed on the fit
## id the per-release tables use. Both the joint walk start and the day-1
## per-stream starts are the ones the per-stream implied-Rt figure on the
## in-sample page uses, so the bands here match it. Confirmed deaths has no trajectory
## here: its panel still draws its release points alone.
_rt_walk_start_joint = clamp(
    _BREAKPOINT - RT_WALK_LEAD, _rt_start_plot, obs.n
)
_stream_chains = (
    "joint" => (;
        chn = chn_joint, rt_start = _rt_start_plot,
        rt_walk_start = _rt_walk_start_joint,
    ),
    "cases" => (; chn = chn_cases, rt_start = 1, rt_walk_start = 1),
    "deaths" => (; chn = chn_deaths, rt_start = 1, rt_walk_start = 1),
    "confirmed" => (; chn = chn_confirmed, rt_start = 1, rt_walk_start = 1),
    "confirmed_deaths" => (;
        chn = chn_confirmed_deaths, rt_start = 1,
        rt_walk_start = 1,
    ),
    "treatment" => (; chn = chn_treatment, rt_start = 1, rt_walk_start = 1),
    "onsets" => (; chn = chn_onsets, rt_start = 1, rt_walk_start = 1),
    "exports" => (; chn = chn_exports, rt_start = 1, rt_walk_start = 1),
)

## Build a fit label => trajectory dictionary from a per-release table,
## restricted to the fits `_stream_chains` names. A fit with no row in `df`
## gets no trajectory, so its panel still draws its release points alone.
function _rt_trajectories(df)
    trajs = Dict{String, Any}()
    for (fid, cfg) in _stream_chains
        fdates = df.date[df.fit .== fid]
        isempty(fdates) && continue
        trajs[get(_fit_labels, fid, fid)] = _stream_rt_trajectory(
            cfg.chn, fdates; rt_start = cfg.rt_start,
            rt_walk_start = cfg.rt_walk_start
        )
    end
    return trajs
end

rt_stream_df = _release_data(
    "rt_by_release_by_stream.csv",
    _by_stream_schema
)
rt_stream_fig = plot_evolution_by_group(
    _fit_groups(rt_stream_df);
    trajectories = _rt_trajectories(rt_stream_df),
    ylabel = "Reproduction number",
    title = "Reproduction number as data accrued, by dataset",
    released_label = "Released estimate (per release)",
    refline = 1.0,
    ymax = _RT_AXIS_MAX,
    empty_note = "No per-dataset reproduction numbers saved yet."
);

#md # ```@raw html
#md # </details>
#md # ```

rt_stream_fig #hide

# ### Basic reproduction number by release
#
# The basic reproduction number $R_0$ estimated at each release, the initial-transmission counterpart of the reproduction number above, before the time-varying decline.
# Released estimates are blue and the current model frozen at earlier cut-offs is red, each a median with nested 30/60/90% interval bars.
# The current fit sits behind both as a flat band, and $R_0 = 1$ is marked.
# Releases only began publishing this quantity recently, so the short blue history reflects that rather than any failed release.

#md # ```@raw html
#md # <details><summary>Basic reproduction number per release with frozen re-fits and the current-fit band</summary>
#md # ```

## Per-release R0 points from r0_by_release.csv, read through the typed
## fallback so a missing or header-only file (until a release carries
## `rt_state.log_R0` in its posterior draws) does not break the build. The
## schema mirrors rt_by_release.csv.
_r0_schema = (;
    release = String, date = Date, median = Float64,
    lo30 = Float64, hi30 = Float64, lo60 = Float64, hi60 = Float64,
    lo90 = Float64, hi90 = Float64,
)
r0_release_df = _release_data("r0_by_release.csv", _r0_schema)
r0_release = [
    (
        string(r.date), r.median, r.lo30, r.hi30, r.lo60, r.hi60,
        r.lo90, r.hi90,
    ) for r in eachrow(r0_release_df)
]

## The current model frozen at earlier cut-offs, one discrete estimate per
## cut-off, reusing the same frozen fits `frozen_matched` above already
## computed. No extra fits are run. Each tuple carries the median and
## 30/60/90% credible bounds of that frozen fit's own R0 draws, unrounded
## since R0 is continuous.
frozen_r0_matched = [
    (c, _ci369(frozen_R0(c); round_fn = identity)...)
        for c in _frozen_matched_cutoffs
]

## The current fit's R0 posterior is a single distribution rather than a
## daily series, so it summarises into a flat 30/60/90% reference band. The
## window runs from the earliest mark on the axis, the first frozen cut-off
## or release point, to the current cut-off, so the band reads behind both
## series rather than only their recent end.
r0_reference = let
    draws = r0_walk_draws(chn_joint)
    q(p) = quantile(draws, p)
    first_date = min(
        minimum(Date.(_frozen_matched_cutoffs)),
        isempty(r0_release_df.date) ? obs.cutoff :
            minimum(r0_release_df.date)
    )
    dates = [first_date, obs.cutoff]
    (
        dates, fill(q(0.35), 2), fill(q(0.65), 2), fill(q(0.2), 2),
        fill(q(0.8), 2), fill(q(0.05), 2), fill(q(0.95), 2),
    )
end

r0_evolution_fig = plot_estimate_evolution(
    r0_release;
    renewal = frozen_r0_matched,
    renewal_label = "Current model frozen at earlier cut-offs",
    trajectory = r0_reference,
    ylabel = "Basic reproduction number",
    title = "Basic reproduction number as data accrued",
    released_label = "Released estimate (per project release)",
    trajectory_label = "Current model, current data",
    refline = 1.0
);

#md # ```@raw html
#md # </details>
#md # ```

r0_evolution_fig #hide

# ### Basic reproduction number by release and dataset
#
# The basic reproduction number estimated at each release, one panel per fit, the by-dataset counterpart of the figure above.
# Panels share a calendar axis and a y range, and $R_0 = 1$ is marked.
# Each release is a median with nested 30/60/90% interval bars.
# Every fit the report runs on its own also carries a current-model reference band.

#md # ```@raw html
#md # <details><summary>Basic reproduction number per release by fit</summary>
#md # ```

## Per-fit R0 flat reference band, the by-dataset counterpart of
## `r0_reference` above, a single distribution rather than a daily walk, so
## each fit's band is flat across its own release window. `r0_walk_draws`
## probes for the walk base, so a single-stream model built without its own
## renewal walk drops its band instead of erroring.
function _r0_stream_trajectory(chn, dates)
    draws = r0_walk_draws(chn)
    isnothing(draws) && return nothing
    q(p) = quantile(draws, p)
    first_date = isempty(dates) ? obs.seeding : minimum(dates)
    ds = [first_date, obs.cutoff]
    return (
        ds, fill(q(0.35), 2), fill(q(0.65), 2), fill(q(0.2), 2),
        fill(q(0.8), 2), fill(q(0.05), 2), fill(q(0.95), 2),
    )
end

## Build a fit label => trajectory dictionary from a per-release R0 table,
## restricted to the fits `_stream_chains` names, the same restriction the
## reproduction-number-by-dataset trajectories use. A fit with no row in
## `df`, or whose chain carries no walk base, gets no trajectory, so its
## panel still draws its release points alone.
function _r0_trajectories(df)
    trajs = Dict{String, Any}()
    for (fid, cfg) in _stream_chains
        fdates = df.date[df.fit .== fid]
        isempty(fdates) && continue
        traj = _r0_stream_trajectory(cfg.chn, fdates)
        isnothing(traj) || (trajs[get(_fit_labels, fid, fid)] = traj)
    end
    return trajs
end

r0_stream_df = _release_data(
    "r0_by_release_by_stream.csv",
    _by_stream_schema
)
r0_stream_fig = plot_evolution_by_group(
    _fit_groups(r0_stream_df);
    trajectories = _r0_trajectories(r0_stream_df),
    ylabel = "Basic reproduction number",
    title = "Basic reproduction number as data accrued, by dataset",
    released_label = "Released estimate (per release)",
    refline = 1.0,
    empty_note = "No per-dataset basic reproduction numbers saved yet."
);

#md # ```@raw html
#md # </details>
#md # ```

r0_stream_fig #hide

# ## Saving forecast outputs
#
# The one-week-back validation forecast, in the same archive format as the release forecast, so the frozen "last week versus now" forecast is recorded as a release asset alongside the forecast it is scored against.

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
    fmt(x) = ismissing(x) || !isfinite(x) ? "n/a" :
        string(round(x; digits = 2))
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
    function block(lead, tbl)
        rows = scored(tbl)
        size(rows, 1) == 0 && return nothing
        return join(
            vcat(
                [string("**", lead, "**"), ""],
                [
                    string(
                        "- ", r.stream, ": relative skill ",
                        fmt(r.rel_to_baseline), ", 90% coverage ",
                        fmt(r.coverage_90), ", bias ", fmt(r.bias), " over ",
                        r.n, " forecasts."
                    )
                        for r in eachrow(rows)
                ]
            ), "\n"
        )
    end
    frozen_joint = select_fit_role(frozen_score_overview_table, "joint")
    overall = [
        beat("Joint model across releases", joint_score_overview_table),
        beat("Frozen joint model", frozen_joint),
    ]
    blocks = filter(
        !isnothing,
        [
            block("Across releases", joint_score_overview_table),
            block("Frozen fits", frozen_joint),
        ]
    )
    join(vcat([join(overall, "\n")], blocks), "\n\n")
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
