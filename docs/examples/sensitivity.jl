# # Sensitivity and comparison analyses
#
# This page continues the [main analysis](analysis.md) from the
# one-week-ahead forecast onward: how last week's forecast held up, how every
# release's forecasts scored against what has since been observed, the
# reproduction number by release, the outbreak size each data stream implies on
# its own, how the estimate has evolved across releases, comparisons with
# McCabe et al. and Chamla et al., and the delay and tree-prior sensitivity
# re-fits. The tree-prior
# sensitivity compares the Skygrid and Exponential growth TMRCA estimates
# [mbalaplacide2026](@cite) (~1.1E-3 subs/site/year, 139 BDBV genomes).
# It renders from the same fitted chains as the main analysis, loaded through
# the shared setup, so no model is re-fit here beyond the frozen and (gated)
# sensitivity re-fits.

#md # ```@raw html
#md # <details><summary>Load packages, data and fitted chains</summary>
#md # ```

## Shared setup: packages, observations, the fit registry and every model fit
## (loaded from the content-addressed cache). See docs/examples/_setup.jl.
using BVDOutbreakSize
include(joinpath(pkgdir(BVDOutbreakSize), "docs", "examples", "_setup.jl"))

#md # ```@raw html
#md # </details>
#md # ```

# ## Forecast validation (last week versus now)
#
# How last week's forecast held up against the data since observed, using the
# frozen re-fit and one-week projection defined in the methods
# [forecast-versus-frozen evaluation](@ref
# "Forecast-versus-frozen evaluation"). The frozen fit also conditions on
# the isolation beds, so the projected bed occupancy is scored against the
# beds held a week later. The bed validation is weak at a one-week-back freeze:
# the reported occupancy rate starts only on 9 June, so the capacity has no
# implied-capacity anchor and rides its random walk back to the freeze date,
# widening the projected bed interval.

#md # ```@raw html
#md # <details><summary>Fit one week back and validate the one-week-ahead forecast</summary>
#md # ```

## frozen_lastweek and frozen_lastweek_streams are computed in the setup
## block above.
## `obs_recovered` is passed so the frozen fit's forecast carries a
## `recovered_new` column (materialised only when the recovered origin is
## given), letting the recovered stream be scored against the observed count
## below like the other streams.
validation_forecast = forecast_reported(frozen_lastweek.chn;
    horizon = 7,
    obs_cases = frozen_lastweek.o.reported_cases,
    obs_deaths = frozen_lastweek.o.total_deaths,
    obs_confirmed = frozen_lastweek.o.confirmed_cases,
    obs_confirmed_deaths = frozen_lastweek.o.confirmed_deaths,
    obs_recovered = frozen_lastweek.o.recovered_cases);

## Each frozen individual (single-stream) fit's own one-week-ahead new-count
## forecast at the same cut-off as `frozen_lastweek`, from
## [`forecast_stream`](@ref) (the same per-stream forecaster
## `stream_forecasts.csv` uses), so the validation plots below can show the
## individual fit alongside the joint rather than the joint alone. Recovered
## has no individual fit and is absent here, as it is throughout this report.
function _validation_individual_new(sid, stream::Symbol, obs_value)
    f = frozen_lastweek_streams[sid]
    bp = f.o.n - f.o.who_first_sitrep_days
    return Float64.(forecast_stream(f.chn, stream; horizon = 7,
        obs_value = obs_value, n = f.o.n, breakpoint = bp,
        rt_start = 1, rt_walk_start = 1))
end
validation_individual = (;
    cases_new = _validation_individual_new(
        "cases", :reported_cases,
        frozen_lastweek_streams["cases"].o.reported_cases),
    deaths_new = _validation_individual_new(
        "deaths", :suspected_deaths,
        frozen_lastweek_streams["deaths"].o.total_deaths),
    confirmed_new = _validation_individual_new(
        "confirmed", :confirmed_cases,
        frozen_lastweek_streams["confirmed"].o.confirmed_cases),
    confirmed_deaths_new = _validation_individual_new(
        "confirmed_deaths", :confirmed_deaths,
        frozen_lastweek_streams["confirmed_deaths"].o.confirmed_deaths))
## The frozen individual (treatment-only) fit's own bed-occupancy forecast,
## anchored on the beds occupied at ITS OWN cut-off (the frozen fit's own
## `o`, not the current `obs`), matching how the joint frozen forecast is
## itself anchored.
_validation_iso_obs = let f = frozen_lastweek_streams["treatment"]
    isempty(f.o.isolation_history.counts) ? 0.0 :
    Float64(f.o.isolation_history.counts[end])
end
validation_individual_isolation = _validation_individual_new(
    "treatment", :isolation_beds, _validation_iso_obs)

## The observed beds at the current cut-off (the forecast target), so the
## frozen-fit bed forecast is scored against what the beds actually held.
_obs_beds = isempty(obs.isolation_history.counts) ? missing :
            obs.isolation_history.counts[end]
## Same observed/baseline keying as the plot below, so the table covers every
## fitted count stream (cumulative and new-count rows) plus the bed level.
validation_table = forecast_vs_truth(validation_forecast;
    observed = (cases_cum = obs.reported_cases,
        deaths_cum = obs.total_deaths,
        confirmed_cum = obs.confirmed_cases,
        confirmed_deaths_cum = obs.confirmed_deaths,
        recovered_cum = obs.recovered_cases),
    baseline = (cases_cum = frozen_lastweek.o.reported_cases,
        deaths_cum = frozen_lastweek.o.total_deaths,
        confirmed_cum = frozen_lastweek.o.confirmed_cases,
        confirmed_deaths_cum = frozen_lastweek.o.confirmed_deaths,
        recovered_cum = frozen_lastweek.o.recovered_cases),
    isolation = _obs_beds);

#md # ```@raw html
#md # </details>
#md # ```

#md # ```@raw html
#md # <details><summary>Forecast-versus-observed validation table</summary>
#md # ```

validation_table #hide

#md # ```@raw html
#md # </details>
#md # ```

# The observation panels histogram the one-week-ahead forecast made from the
# frozen fit, a cumulative and a new-count panel for each fitted count stream
# the forecast carries (reported cases, suspected deaths, confirmed cases,
# confirmed deaths and recovered), with the 90% predictive interval shaded and
# the count actually observed by the current cut-off drawn as a dashed black
# rule. Each stream draws only when the forecast carries its column, so a fit
# observing fewer streams shows fewer panels.
# Where a stream has its own individual (single-stream) fit, that fit's own
# forecast, made from the same frozen cut-off, is overlaid as a dotted density
# alongside the joint's histogram, so both are visible rather than the joint
# alone; recovered has no individual fit and draws the joint alone.

#md # ```@raw html
#md # <details><summary>Forecast-versus-observed plot</summary>
#md # ```

## Observed cumulative at the target date per stream, keyed by the forecast's
## cumulative column; `baseline` is each stream's origin cumulative (the frozen
## cut-off), so the new-count panel is scored against observed minus origin.
validation_fig = plot_forecast_vs_truth(validation_forecast;
    observed = (cases_cum = obs.reported_cases,
        deaths_cum = obs.total_deaths,
        confirmed_cum = obs.confirmed_cases,
        confirmed_deaths_cum = obs.confirmed_deaths,
        recovered_cum = obs.recovered_cases),
    baseline = (cases_cum = frozen_lastweek.o.reported_cases,
        deaths_cum = frozen_lastweek.o.total_deaths,
        confirmed_cum = frozen_lastweek.o.confirmed_cases,
        confirmed_deaths_cum = frozen_lastweek.o.confirmed_deaths,
        recovered_cum = frozen_lastweek.o.recovered_cases),
    individual = validation_individual);

#md # ```@raw html
#md # </details>
#md # ```

validation_fig #hide

# The bed panel scores last week's projected occupancy against the beds
# occupied now (the dashed rule), the individual (treatment-only) fit's own
# projection overlaid as a dotted density alongside the joint.

#md # ```@raw html
#md # <details><summary>Bed forecast-versus-observed plot</summary>
#md # ```

validation_beds_fig = plot_forecast_beds_vs_truth(validation_forecast;
    isolation = _obs_beds, individual = validation_individual_isolation);

#md # ```@raw html
#md # </details>
#md # ```

validation_beds_fig #hide

# The latent quantities are not observed, so they are scored distribution
# against distribution: what the frozen fit forecast for the past week's new
# infections, onsets and deaths against what the current fit now estimates
# for the same window.

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
    deaths_latent_new = _now_new(chn_joint, :cumulative_expected_deaths))

validation_latent_fig = plot_forecast_vs_truth_latent(
    validation_forecast; now = now_latent);

#md # ```@raw html
#md # </details>
#md # ```

validation_latent_fig #hide

# ## Forecast scoring across releases
#
# How every release's saved one- to four-week-ahead forecast scored against the data observed since, against a persistence baseline and, where one exists, the stream's own individual fit as well as the joint (see [forecast scoring against a persistence baseline](@ref "Forecast scoring against a persistence baseline") for how the scores, the relative skill and the baseline are built).
# Recovered has no individual fit of its own, so its comparison is the baseline against the joint only.
# Reported cases and suspected deaths stopped being updated by the situation reports partway through the outbreak, and exports' confirmed-detection series is anchored to an earlier cut-off, so exports contributes no scored forecast at all and reported cases and suspected deaths each rest on exactly one matched forecast, a single window rather than a settled sample.
#
# Only a minority of the daily releases examined across the outbreak contribute a row to the table below, and every one of them is a reconstruction of an earlier version of the model rather than the model as it stands today, so the table below is not a verdict on the current fit.
# One reconstruction is dropped from scoring entirely: its chain forecasts a near-zero median at every horizon and stream, with the upper predictive tail occasionally reaching five- and six-digit values, which is the signature of a chain that failed to sample properly rather than a genuine forecast (see `scripts/score_releases.jl`, `_FAILED_RECONSTRUCTIONS`).
# The releases that do carry the current model's own individual-stream forecasts are too recent for their targets to be observed yet, which is why the comparison against each stream's individual fit is still empty.
# It will populate once those targets resolve.
# Every row also rests on one to a handful of matched forecasts, shown as its own count rather than rounded away, so a ratio here should be read as an early signal rather than a settled result.

#md # ```@raw html
#md # <details><summary>Load and summarise the cross-release forecast scores</summary>
#md # ```

## scripts/score_releases.jl writes these after the fits. The committed files
## are header-only until a release carries the asset, so the common path
## reads a real file to a zero-row frame; the typed `schema` is the fallback
## for a file that is absent entirely, since CSV.read throws on a missing
## path and would take the whole docs build with it.
function _release_data(name, schema::NamedTuple)
    path = joinpath(pkgdir(BVDOutbreakSize), "data", name)
    isfile(path) && return CSV.read(path, DataFrame)
    return DataFrame([k => T[] for (k, T) in pairs(schema)])
end

forecast_scores_df = _release_data("forecast_scores.csv",
    (; release = String, made_date = Date, stream = String, horizon = Int,
        target_date = Date, fit = String, crps = Float64,
        log_crps = Float64, dispersion = Float64, overprediction = Float64,
        underprediction = Float64, coverage_50 = Float64,
        coverage_90 = Float64,
        bias = Float64, n_samples = Int,
        log_rel_to_baseline = Float64))
forecast_overlay_df = _release_data("forecast_overlay.csv",
    (; release = String, made_date = Date, stream = String, horizon = Int,
        target_date = Date, fit = String, observed = Float64,
        median = Float64, lo30 = Float64, hi30 = Float64, lo60 = Float64,
        hi60 = Float64, lo90 = Float64, hi90 = Float64))
## The headline table, one row per (stream, fit) pooled over every horizon
## and release; the by-horizon and by-release detail tables carry the same
## columns at a finer grain (see src/scoring.jl).
forecast_score_overview_table = forecast_score_overview(forecast_scores_df)
forecast_score_by_horizon_table = forecast_score_by_horizon(forecast_scores_df)
forecast_score_by_release_table = forecast_score_by_release(forecast_scores_df)

#md # ```@raw html
#md # </details>
#md # ```

# The headline pools every horizon and release into one row per stream and fit: the mean CRPS and its decomposition, coverage and bias, and the relative skill against the persistence baseline for every fit and, on the joint row where the stream has one, against its own individual fit, on both the natural and the log scale (column definitions in [forecast scoring against a persistence baseline](@ref "Forecast scoring against a persistence baseline")).

forecast_score_overview_table #hide

# The same relative skill against the baseline, by horizon: one panel per stream, one series per fit role, on a log-scaled skill axis with the reference line at one.
# A fit that beats the baseline on average but not at every cut-off is visible as a series that crosses the line rather than sitting under it throughout.

forecast_relative_skill_fig = plot_forecast_relative_skill(
    forecast_score_by_horizon_table);

forecast_relative_skill_fig #hide

# The same columns as a table, broken out by horizon, and again broken out by release and averaged across horizons, are behind the two dropdowns below.

#md # ```@raw html
#md # <details><summary>Scores by horizon</summary>
#md # ```

forecast_score_by_horizon_table #hide

#md # ```@raw html
#md # </details>
#md # ```

#md # ```@raw html
#md # <details><summary>Scores by release</summary>
#md # ```

forecast_score_by_release_table #hide

#md # ```@raw html
#md # </details>
#md # ```

# Forecasts made at each release against the value observed since, one panel per stream and horizon, the observed value in black.
# The median and 90% interval are coloured by fit role: the persistence baseline, the stream's individual fit and the joint.
# The x-axis is the date each forecast was made, so an incident stream's observed window pairs unambiguously with the forecast that made it.
# Each panel's axis is cropped to a small multiple of what that stream actually reached, so one very wide interval cannot squash every other series flat; an interval or median too wide for the panel is clamped at the top and marked with an open triangle rather than silently cut off.

#md # ```@raw html
#md # <details><summary>Forecasts-versus-now overlay</summary>
#md # ```

forecast_overlay_fig = plot_forecast_overlay(forecast_overlay_df);

#md # ```@raw html
#md # </details>
#md # ```

forecast_overlay_fig #hide

# ## Frozen-fit forecast evaluation
#
# The current model, frozen at earlier data cut-offs (see [Forecast-versus-frozen evaluation](@ref "Forecast-versus-frozen evaluation")), scored the same way as the cross-release forecasts above, against the same persistence baseline.
# The May cut-offs predate the first reported bed occupancy and the first reported recoveries, so those windows are dropped rather than scored against a series that had not started.
# The frozen baseline's data-vintage guarantee is weaker than the cross-release one above: see [forecast scoring against a persistence baseline](@ref "Forecast scoring against a persistence baseline") for the residual gap this leaves for a frozen made_date reused weeks ahead of its snapshot's own cut-off.

#md # ```@raw html
#md # <details><summary>Load and summarise the frozen-fit forecast scores</summary>
#md # ```

frozen_scores_df = _release_data("forecast_scores_frozen.csv",
    (; release = String, made_date = Date, stream = String, horizon = Int,
        target_date = Date, fit = String, crps = Float64,
        log_crps = Float64, dispersion = Float64, overprediction = Float64,
        underprediction = Float64, coverage_50 = Float64,
        coverage_90 = Float64,
        bias = Float64, n_samples = Int,
        log_rel_to_baseline = Float64))
frozen_overlay_df = _release_data("forecast_overlay_frozen.csv",
    (; release = String, made_date = Date, stream = String, horizon = Int,
        target_date = Date, fit = String, observed = Float64,
        median = Float64, lo30 = Float64, hi30 = Float64, lo60 = Float64,
        hi60 = Float64, lo90 = Float64, hi90 = Float64))
## The frozen evaluation never carries an individual single-stream fit
## (it scores only the joint model at past cut-offs), so the individual-fit
## comparison columns are dropped rather than shown as a column of missing.
frozen_score_overview_table = drop_individual_fit_columns(
    forecast_score_overview(frozen_scores_df))
frozen_score_by_horizon_table = drop_individual_fit_columns(
    forecast_score_by_horizon(frozen_scores_df))
frozen_score_by_release_table = drop_individual_fit_columns(
    forecast_score_by_release(frozen_scores_df))

## The frozen evaluation scores the current JOINT model alone at every
## cut-off, so `fit` is single-valued ("joint") by construction in every
## one of these tables, not just for the releases scored so far. Dropped
## from the DISPLAY tables below as a degenerate column (see
## `drop_degenerate_fit_column`); the `..._table` frames above keep `fit`
## and still feed the relative-skill plots, which read it to colour each
## series.
frozen_score_overview_display = drop_degenerate_fit_column(
    frozen_score_overview_table)
frozen_score_by_horizon_display = drop_degenerate_fit_column(
    frozen_score_by_horizon_table)
frozen_score_by_release_display = drop_degenerate_fit_column(
    frozen_score_by_release_table)

#md # ```@raw html
#md # </details>
#md # ```

# The headline frozen table, one row per stream pooled across cut-offs and horizons, against the persistence baseline on both scales, with the same CRPS decomposition, coverage and bias columns described above.
# Every row is the current JOINT model; the frozen evaluation carries no other model, so the `fit` column is omitted here.

frozen_score_overview_display #hide

# The same relative skill against the baseline, by horizon, for the frozen cut-offs.
# The frozen evaluation scores one fit at each cut-off, so each panel carries a single series, drawn in the joint colour.

frozen_relative_skill_fig = plot_forecast_relative_skill(
    frozen_score_by_horizon_table);

frozen_relative_skill_fig #hide

# The same columns as a table, broken out by horizon, and again broken out by frozen cut-off, are behind the two dropdowns below.
# Both are the current joint model only, the same as the headline table above.

#md # ```@raw html
#md # <details><summary>Scores by horizon</summary>
#md # ```

frozen_score_by_horizon_display #hide

#md # ```@raw html
#md # </details>
#md # ```

#md # ```@raw html
#md # <details><summary>Scores by frozen cut-off</summary>
#md # ```

frozen_score_by_release_display #hide

#md # ```@raw html
#md # </details>
#md # ```

# The frozen forecasts made at each cut-off against the value observed since, one panel per stream and horizon, the observed value in black.
# The frozen forecast draws in the joint colour against the persistence baseline.
# The x-axis is the frozen cut-off each forecast was made from.
# The same per-panel axis crop and overflow marker described above applies here too.

#md # ```@raw html
#md # <details><summary>Frozen-fit forecasts-versus-now overlay</summary>
#md # ```

frozen_overlay_fig = plot_forecast_overlay(frozen_overlay_df);

#md # ```@raw html
#md # </details>
#md # ```

frozen_overlay_fig #hide

# The frozen re-fits below freeze the renewal data to an earlier cut-off
# and re-fit, so that a change driven by newer data can be distinguished from
# one driven by a change of method.
# Each uses the full headline settings (1000 draws across two chains),
# reusing the frozen-fit helper defined above.

#md # ```@raw html
#md # <details><summary>Freeze the renewal data to a cut-off and re-fit</summary>
#md # ```

## Frozen re-fits and released_df are prepared in the setup block above.

#md # ```@raw html
#md # </details>
#md # ```

# ## Individual fits against the baseline
#
# The same cross-release forecast scoring as [Forecast scoring across
# releases](@ref "Forecast scoring across releases") above, filtered to each
# stream's own individual single-stream fit rather than the joint, against
# the same persistence baseline.
# Recovered has no individual fit, so it does not appear here.
# `forecast_score_overview_table` and its by-horizon and by-release
# counterparts already carry one row per (stream, fit) for every fit
# scored, joint and individual alike, so this section reuses those same
# tables filtered to the individual-fit rows rather than building its own.
# As noted in [Forecast scoring across releases](@ref "Forecast scoring
# across releases") above, the releases that carry the current model's own
# individual-stream forecasts are too recent for their targets to be
# observed yet, so every table and figure in this section is currently
# empty; it will populate once those targets resolve.

#md # ```@raw html
#md # <details><summary>Individual-fit rows of the cross-release scores</summary>
#md # ```

individual_score_overview_table = forecast_score_overview_table[
    forecast_score_overview_table.fit .!= "joint", :]
individual_score_by_horizon_table = forecast_score_by_horizon_table[
    forecast_score_by_horizon_table.fit .!= "joint", :]
individual_score_by_release_table = forecast_score_by_release_table[
    forecast_score_by_release_table.fit .!= "joint", :]

#md # ```@raw html
#md # </details>
#md # ```

individual_score_overview_table #hide

# The same relative skill against the baseline, by horizon, one panel per stream (dataset), for each stream's own individual fit.

individual_relative_skill_fig = plot_forecast_relative_skill(
    individual_score_by_horizon_table;
    empty_message = "Empty: the releases that carry the current model's " *
                    "own individual-stream forecasts are too recent for their " *
                    "targets to be observed yet. Not a missing forecast.");

individual_relative_skill_fig #hide

# The same columns as a table, broken out by horizon, and again broken out by release, are behind the two dropdowns below.

#md # ```@raw html
#md # <details><summary>Scores by horizon</summary>
#md # ```

individual_score_by_horizon_table #hide

#md # ```@raw html
#md # </details>
#md # ```

#md # ```@raw html
#md # <details><summary>Scores by release</summary>
#md # ```

individual_score_by_release_table #hide

#md # ```@raw html
#md # </details>
#md # ```

# ## Outbreak size estimated by each data stream
#
# Each data stream constrains the latent outbreak size differently.
# The table below puts the posteriors over the infection count side by side,
# the single-stream fits and the joint, to show what each stream implies on
# its own and what the joint adds.

#md # ```@raw html
#md # <details><summary>Per-stream infection-count table</summary>
#md # ```

streams_C_table = streams_table(
    "exports" => posterior_C_exports,
    "deaths (DRC)" => posterior_C_deaths,
    "cases (DRC)" => posterior_C_cases,
    "confirmed (DRC)" => posterior_C_confirmed,
    "isolation (DRC)" => posterior_C_treatment,
    "joint" => posterior_C_joint);

#md # ```@raw html
#md # </details>
#md # ```

streams_C_table #hide

# The first figure shows each single-stream fit's cumulative-infection
# trajectory projected to the cut-off, with a dotted rule in each stream's
# colour marking where its data stops and the ribbon beyond it becomes a
# forward projection.

#md # ```@raw html
#md # <details><summary>Per-stream projected-trajectory plot</summary>
#md # ```

## Per-draw cumulative-infection trajectory carried by each single-stream
## fit out to the cut-off on day `n`, so streams whose data ends earlier are
## still projected to today.
function _cuminf(chn)
    mat = chn[:cumulative_infections]
    return [collect(v) for v in vec(collect(mat))]
end
## Grid day a stream's data last reports, used for the dotted rule. The
## suspected case and death histories freeze at 26 May; exports and confirmed
## run to the cut-off.
_last_day(days) = isempty(days) ? nothing : maximum(days)

stream_traj_fig = plot_stream_trajectories(
    [
        (; label = "exports", trajs = _cuminf(chn_exports),
            last_day = _last_day(vcat(obs.export_case_days,
                obs.export_death_days)), colour = :seagreen),
        (; label = "deaths (DRC)", trajs = _cuminf(chn_deaths),
            last_day = _last_day(obs.deaths_history.days),
            colour = :firebrick),
        (; label = "cases (DRC)", trajs = _cuminf(chn_cases),
            last_day = _last_day(obs.reported_history.days),
            colour = :steelblue),
        (; label = "confirmed (DRC)", trajs = _cuminf(chn_confirmed),
            last_day = _last_day(obs.confirmed_history.days),
            colour = :goldenrod),
        (; label = "isolation (DRC)", trajs = _cuminf(chn_treatment),
            last_day = _last_day(obs.isolation_history.days),
            colour = :darkorange)];
    n = obs.n, seeding = obs.seeding);

#md # ```@raw html
#md # </details>
#md # ```

stream_traj_fig #hide

# The second figure is the posterior density of each fit's cumulative
# infection count at the cut-off; the x-axis is scaled to a multiple of the
# joint-fit 90% upper bound so the bulk of the streams stays visible rather
# than being flattened by the wide, ill-defined confirmed-only tail.

#md # ```@raw html
#md # <details><summary>Cut-off infection-count density plot</summary>
#md # ```

## Scale the x-axis to twice the joint-fit 90% upper bound, so the joint and
## the streams that track it read clearly while the confirmed-only tail runs
## off the axis rather than dominating it.
density_xmax = 2.0 * quantile(posterior_C_joint, 0.95)

cumulative_density_fig = plot_cumulative_cases(
    "exports" => posterior_C_exports,
    "deaths (DRC)" => posterior_C_deaths,
    "cases (DRC)" => posterior_C_cases,
    "confirmed (DRC)" => posterior_C_confirmed,
    "isolation (DRC)" => posterior_C_treatment,
    "joint" => posterior_C_joint;
    scenarios = [], xmax = density_xmax);

#md # ```@raw html
#md # </details>
#md # ```

cumulative_density_fig #hide

# ## Outbreak size by release and dataset
#
# The outbreak size $C_T$ estimated at each release, one panel per fit, the per-dataset analogue of the overall release evolution below.
# Each dataset constrains the latent size differently, so the panels show how far each fit's history sits from the joint's as data accrued.
# A dataset the report fits on its own also carries that fit's current-model outbreak-size band, the per-dataset analogue of the current-data trajectory below, so a panel spanning only one saved release still reads against its own current estimate.
# Confirmed deaths has no band of its own here, so its panel carries its release points alone.
# The per-dataset release history begins at the release that first published these per-fit estimates, so some panels open later than the overall outbreak-size picture below.
# Each panel keeps its own y range rather than a shared one, since the datasets span very different scales.

#md # ```@raw html
#md # <details><summary>Outbreak size per release, faceted by fit</summary>
#md # ```

## Schema of the per-release, per-fit estimate tables written by
## scripts/score_releases.jl from each release's stream_estimates.csv. Shared
## by this section and the reproduction-number-by-dataset section below.
_by_stream_schema = (; release = String, date = Date, fit = String,
    median = Float64, lo30 = Float64, hi30 = Float64, lo60 = Float64,
    hi60 = Float64, lo90 = Float64, hi90 = Float64)

## Fits in a fixed order, the joint first, so the panels do not reshuffle
## between builds. Labels match the per-stream table above (in "Outbreak
## size estimated by each data stream"). Recovered is absent because it
## has no individual fit.
_fit_order = ["joint", "cases", "deaths", "confirmed", "confirmed_deaths",
    "treatment", "exports"]
_fit_labels = Dict("joint" => "joint", "cases" => "cases (DRC)",
    "deaths" => "deaths (DRC)", "confirmed" => "confirmed (DRC)",
    "confirmed_deaths" => "confirmed deaths (DRC)",
    "treatment" => "isolation (DRC)", "exports" => "exports")

## Group a per-fit estimate table into the label => tuples pairs the faceted
## plot takes, keyed on the date so the mixed release tag shapes
## (`results-v1.9.0` and `results-1243`) never reach the axis.
function _fit_groups(df)
    return [get(_fit_labels, f, f) =>
                [(string(r.date), r.median, r.lo30, r.hi30, r.lo60, r.hi60,
                     r.lo90, r.hi90) for r in eachrow(df) if r.fit == f]
            for f in _fit_order]
end

## Current-data, current-model cumulative-infection trajectory for one
## single-stream fit, the per-dataset analogue of the joint trajectory built
## for the release-evolution section below. `dates` is that fit's own saved
## release dates, so the band starts at the same release its per-fit panel's
## points begin at.
function _stream_size_trajectory(chn, dates)
    mat = chn[:cumulative_infections]
    trajs = [collect(v) for v in vec(collect(mat))]
    first_date = isempty(dates) ? obs.seeding : minimum(dates)
    start_day = obs.n - value(obs.cutoff - first_date)
    days = max(start_day, 1):obs.n
    ds = [obs.seeding + Day(d - 1) for d in days]
    q(d, p) = quantile(Float64[t[d] for t in trajs], p)
    (ds,
        [q(d, 0.35) for d in days], [q(d, 0.65) for d in days],
        [q(d, 0.20) for d in days], [q(d, 0.80) for d in days],
        [q(d, 0.05) for d in days], [q(d, 0.95) for d in days])
end

## Same, for the reproduction number, reconstructing the walk exactly as
## `plot_rt_streams` does per stream.
function _stream_rt_trajectory(chn, dates; rt_start, rt_walk_start)
    mat = reconstruct_rt(chn; n = obs.n, breakpoint = _BREAKPOINT,
        rt_start = rt_start, rt_walk_start = rt_walk_start,
        ramp = RT_INTERVENTION_RAMP)
    first_date = isempty(dates) ? obs.seeding : minimum(dates)
    first_day = clamp(value(first_date - obs.seeding) + 1, rt_start, obs.n)
    days = first_day:obs.n
    ds = [obs.seeding + Day(d - 1) for d in days]
    q(d, p) = quantile(collect(skipmissing(@view mat[:, d])), p)
    (ds,
        [q(d, 0.35) for d in days], [q(d, 0.65) for d in days],
        [q(d, 0.20) for d in days], [q(d, 0.80) for d in days],
        [q(d, 0.05) for d in days], [q(d, 0.95) for d in days])
end

## The single-stream chains and their renewal-walk starts, keyed on the fit
## id the per-release tables use. The joint walks from `RT_WALK_LEAD` days
## before the first situation report; the per-stream fits walk from day 1,
## matching the per-stream implied-Rt figure further below. Confirmed deaths
## has no trajectory here: its panel still draws its release points alone.
_stream_rt_walk_start = clamp(
    _BREAKPOINT - RT_WALK_LEAD, _rt_start_plot, obs.n)
_stream_chains = (
    "joint" => (; chn = chn_joint, rt_start = _rt_start_plot,
        rt_walk_start = _stream_rt_walk_start),
    "cases" => (; chn = chn_cases, rt_start = 1, rt_walk_start = 1),
    "deaths" => (; chn = chn_deaths, rt_start = 1, rt_walk_start = 1),
    "confirmed" => (; chn = chn_confirmed, rt_start = 1, rt_walk_start = 1),
    "treatment" => (; chn = chn_treatment, rt_start = 1, rt_walk_start = 1),
    "exports" => (; chn = chn_exports, rt_start = 1, rt_walk_start = 1))

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
            rt_walk_start = cfg.rt_walk_start)
    end
    return trajs
end
function _size_trajectories(df)
    trajs = Dict{String, Any}()
    for (fid, cfg) in _stream_chains
        fdates = df.date[df.fit .== fid]
        isempty(fdates) && continue
        trajs[get(_fit_labels, fid, fid)] = _stream_size_trajectory(
            cfg.chn, fdates)
    end
    return trajs
end

size_stream_df = _release_data("size_by_release_by_stream.csv",
    _by_stream_schema)
## Outbreak size spans very different scales across datasets (exports, a
## handful of detections, against the joint's tens of thousands of
## infections), so each panel keeps its own y range rather than sharing one
## that a single wide-scale dataset would set for every other panel.
size_stream_fig = plot_evolution_by_group(_fit_groups(size_stream_df);
    trajectories = _size_trajectories(size_stream_df),
    ylabel = "Cumulative infections",
    title = "Outbreak-size estimate as data accrued, by dataset",
    released_label = "Released estimate (per release)",
    shared_yrange = false,
    empty_note = "No per-dataset outbreak sizes saved yet.");

#md # ```@raw html
#md # </details>
#md # ```

size_stream_fig #hide

# ## Estimate evolution across releases
#
# How the outbreak-size estimate moved as the situation reports accrued.
#
# The project publishes a tagged results release at each data cut-off
# (<https://github.com/epiforecasts/BVDOutbreakSize/releases>), bundling
# the posterior draws and input data.
# The released series, in blue, is the project's published estimate at each
# release: the closed-form integral model up to v1.3.0, then the renewal
# model from v1.4.0 on.
# Each release is its own fit, so it is drawn as a discrete estimate, a
# median with nested 30/60/90% interval bars, rather than a ribbon.
# The current model frozen at earlier cut-offs is drawn in red as discrete
# estimates: at the cut-offs matched to the McCabe et al. reports (20, 23
# and 27 May), Chamla et al.'s 8 June confirmed-case calibration date, and
# one week before the current cut-off. These are the frozen joint fits
# already computed for the matched-in-time comparison, the Chamla comparison
# and the forecast validation, reused here, so no extra fits are run.
# The current-data, current-model estimate is drawn in green as the
# cumulative-infection trajectory over time, a single fit shown across the
# period so the latest estimate reads against the earlier ones.
# Each release date is marked with a dotted vertical rule.

#md # ```@raw html
#md # <details><summary>Released estimates and the current-model frozen re-fits</summary>
#md # ```

## Released median and 30/60/90% intervals per release, from
## `data/released_estimates.csv`. Each tuple is
## `(date, median, lo30, hi30, lo60, hi60, lo90, hi90)`.
release_evolution = [(string(r.date), r.median, r.lo30, r.hi30, r.lo60, r.hi60,
                         r.lo90, r.hi90) for r in eachrow(released_df)]

## The current model frozen at earlier cut-offs, each its own discrete
## estimate: the matched-McCabe cut-offs (20, 23, 27 May) already computed
## for the matched-in-time comparison below, the 8 June Chamla
## confirmed-case anchor computed for the Chamla comparison, and the
## one-week-back validation fit (`frozen_lastweek`, at `validation_cutoff`)
## already computed for the forecast validation above. All are reused here so
## the current-model estimate at those earlier cut-offs reads against the
## released overlay, including a recent point one week before the cut-off.
## No extra fits are run. Each tuple carries the median and 30/60/90%
## credible bounds from the frozen `C_T` draws.
function _ci369(xs)
    q(p) = round(Int, quantile(xs, p))
    (q(0.5), q(0.35), q(0.65), q(0.20), q(0.80), q(0.05), q(0.95))
end
frozen_by_cutoff[validation_cutoff] = frozen_lastweek
frozen_matched = [(c, _ci369(frozen_C(c))...)
                  for c in sort(union(frozen_cutoffs,
    [validation_cutoff, default_chamla_cutoff()]))]

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
    (dates,
        [q(d, 0.35) for d in days], [q(d, 0.65) for d in days],
        [q(d, 0.20) for d in days], [q(d, 0.80) for d in days],
        [q(d, 0.05) for d in days], [q(d, 0.95) for d in days])
end

evolution_fig = plot_estimate_evolution(release_evolution;
    renewal = frozen_matched,
    renewal_label = "Current model frozen at earlier cut-offs",
    trajectory = infection_trajectory,
    title = "Outbreak-size estimate as data accrued");

#md # ```@raw html
#md # </details>
#md # ```

evolution_fig #hide

# ## Reproduction number estimated by each data stream
#
# The reproduction number each stream implies on its own, one panel per stream with the joint fit overlaid in grey as the reference.

#md # ```@raw html
#md # <details><summary>Per-stream implied-Rt plot</summary>
#md # ```

## The per-stream fits walk Rt from day 1 (the default `rt_start`), while the
## joint walks from `RT_WALK_LEAD` days before the first situation report; the
## shared `display_start` is the joint renewal start so every stream reads over
## the same established window. `ramp` matches the joint Rt figure.
_rt_walk_start_joint = clamp(_BREAKPOINT - RT_WALK_LEAD, _rt_start_plot, obs.n);
stream_rt_fig = plot_rt_streams(
    [
        (; label = "exports", chn = chn_exports, rt_start = 1,
            rt_walk_start = 1, colour = :seagreen),
        (; label = "deaths (DRC)", chn = chn_deaths, rt_start = 1,
            rt_walk_start = 1, colour = :firebrick),
        (; label = "cases (DRC)", chn = chn_cases, rt_start = 1,
            rt_walk_start = 1, colour = :steelblue),
        (; label = "confirmed (DRC)", chn = chn_confirmed, rt_start = 1,
            rt_walk_start = 1, colour = :goldenrod),
        (; label = "isolation (DRC)", chn = chn_treatment, rt_start = 1,
            rt_walk_start = 1, colour = :darkorange)];
    joint = (; label = "joint", chn = chn_joint, rt_start = _rt_start_plot,
        rt_walk_start = _rt_walk_start_joint),
    n = obs.n, breakpoint = _BREAKPOINT,
    as_of_date = string(obs.cutoff), seeding = obs.seeding,
    display_start = _rt_start_plot, ramp = RT_INTERVENTION_RAMP);

#md # ```@raw html
#md # </details>
#md # ```

stream_rt_fig #hide

# ## Reproduction number by release and dataset
#
# The same per-release reproduction number, one panel per fit, so each dataset's history reads against the others and against the joint.
# Panels share a calendar axis and a y range, $R_t = 1$ is marked, and each release's cut-off reproduction number is a median with nested 30/60/90% interval bars.
# A fit with no saved estimates is dropped rather than drawn empty.
# A dataset the report fits on its own also carries that fit's current-model reproduction number as a band, reconstructed the same way as the single overall fit below, so a panel spanning only one release still reads against its own current estimate rather than sitting alone.
# Confirmed deaths has no band of its own here, so its panel carries its release points alone.
# The per-dataset release history begins at the release that first published these per-fit estimates, so some panels open later than the overall reproduction-number picture below.

#md # ```@raw html
#md # <details><summary>Reproduction number per release by fit</summary>
#md # ```

## Schema, fit order/labels and grouping and trajectory helpers are shared
## with, and already defined by, the outbreak-size-by-dataset section above.
rt_stream_df = _release_data("rt_by_release_by_stream.csv",
    _by_stream_schema)
rt_stream_fig = plot_evolution_by_group(_fit_groups(rt_stream_df);
    trajectories = _rt_trajectories(rt_stream_df),
    ylabel = "Reproduction number",
    title = "Reproduction number as data accrued, by dataset",
    released_label = "Released estimate (per release)",
    refline = 1.0,
    empty_note = "No per-dataset reproduction numbers saved yet.");

#md # ```@raw html
#md # </details>
#md # ```

rt_stream_fig #hide

# ## Reproduction number by release
#
# The reproduction number estimated at each release, the same kind of release-by-release picture as the outbreak-size evolution above.
# Each release's cut-off reproduction number $R_T$ is drawn as a discrete estimate, a median with nested 30/60/90% interval bars.
# The current fit's daily $R_t$ over its established window is drawn as the continuous band, and $R_t = 1$ is marked.

#md # ```@raw html
#md # <details><summary>Reproduction number per release with the current-fit band</summary>
#md # ```

rt_release_df = CSV.read(
    joinpath(pkgdir(BVDOutbreakSize), "data", "rt_by_release.csv"), DataFrame)
rt_release = [(string(r.date), r.median, r.lo30, r.hi30, r.lo60, r.hi60,
                  r.lo90, r.hi90) for r in eachrow(rt_release_df)]

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
    mat = reconstruct_rt(chn_joint; n = obs.n, breakpoint = _BREAKPOINT,
        rt_start = _rt_start_plot, rt_walk_start = rt_walk_start,
        ramp = RT_INTERVENTION_RAMP)
    first_release_day = clamp(
        value(minimum(rt_release_df.date) - obs.seeding) + 1,
        _rt_start_plot, obs.n)
    days = first_release_day:obs.n
    dates = [obs.seeding + Day(d - 1) for d in days]
    q(d, p) = quantile(collect(skipmissing(@view mat[:, d])), p)
    (dates,
        [q(d, 0.35) for d in days], [q(d, 0.65) for d in days],
        [q(d, 0.20) for d in days], [q(d, 0.80) for d in days],
        [q(d, 0.05) for d in days], [q(d, 0.95) for d in days])
end

rt_evolution_fig = plot_estimate_evolution(rt_release;
    trajectory = rt_release_trajectory,
    ylabel = "Reproduction number",
    title = "Reproduction number as data accrued",
    released_label = "Released estimate (per project release)",
    trajectory_label = "Current model, current data",
    refline = 1.0);

#md # ```@raw html
#md # </details>
#md # ```

rt_evolution_fig #hide

# ## Basic reproduction number by release
#
# The basic reproduction number $R_0$, the renewal walk's starting value before the time-varying decline, estimated at each release, the initial-transmission analogue of the cut-off reproduction number above.
# Each release's $R_0$ posterior is drawn as a discrete estimate, a median with nested 30/60/90% interval bars.
# The current fit's $R_0$ posterior is the flat reference band across the release window, and $R_0 = 1$ is marked.
# This quantity is only saved by a release once the release started publishing it, so the discrete points begin at the release that introduced it rather than at the start of the release history: the underlying posterior column landed only a few releases before this page was last regenerated, so only the most recent releases carry a point here.
# This is a genuinely sparse start, not a dropped or failed release: every release scored above produced a forecast and an $R_T$ estimate without error, and the gap closes release by release from here as more releases carry the column.

#md # ```@raw html
#md # <details><summary>Basic reproduction number per release with the current-fit band</summary>
#md # ```

## Per-release R0 points from r0_by_release.csv, read through the typed
## fallback so a missing or header-only file (until a release carries
## `rt_state.log_R0` in its posterior draws) does not break the build. The
## schema mirrors rt_by_release.csv.
_r0_schema = (; release = String, date = Date, median = Float64,
    lo30 = Float64, hi30 = Float64, lo60 = Float64, hi60 = Float64,
    lo90 = Float64, hi90 = Float64)
r0_release_df = _release_data("r0_by_release.csv", _r0_schema)
r0_release = [(string(r.date), r.median, r.lo30, r.hi30, r.lo60, r.hi60,
                  r.lo90, r.hi90) for r in eachrow(r0_release_df)]

## The current fit's R0 posterior is the exponential of the renewal walk's log
## base `rt_state.log_R0`, a single distribution rather than a daily series.
## Summarise it into a flat 30/60/90% reference band spanning the release
## window (first release date to the cut-off), so it reads across the
## per-release points like the time-varying band in the Rt figure. When no
## release has been scored yet the window falls back to the seeding date.
r0_reference = let
    draws = exp.(vec(Array(chn_joint[Symbol("rt_state.log_R0")])))
    q(p) = quantile(draws, p)
    first_date = isempty(r0_release_df.date) ? obs.seeding :
                 minimum(r0_release_df.date)
    dates = [first_date, obs.cutoff]
    (dates, fill(q(0.35), 2), fill(q(0.65), 2), fill(q(0.20), 2),
        fill(q(0.80), 2), fill(q(0.05), 2), fill(q(0.95), 2))
end

r0_evolution_fig = plot_estimate_evolution(r0_release;
    trajectory = r0_reference,
    ylabel = "Basic reproduction number",
    title = "Basic reproduction number as data accrued",
    released_label = "Released estimate (per project release)",
    trajectory_label = "Current model, current data",
    refline = 1.0);

#md # ```@raw html
#md # </details>
#md # ```

r0_evolution_fig #hide

# ## Comparison with McCabe et al.
#
# Our model is a discrete-time renewal model with a time-varying
# reproduction number and every data stream fitted jointly.
# McCabe et al. published their estimates as scenarios at fixed
# situation-report cut-offs, each scenario carrying a 95% confidence
# interval.
# We show all three, the 18 May report, the 20 May update and the 27 May
# Lancet publication, as one panel each, with their intervals kept.
# Within a panel each method and scenario family is a single line, carrying its
# sweep over the nuisance assumptions: the case-fatality ratio, the geographic
# window and the doubling time.
# The geographic-spread scenarios come from exported cases and travel volume.
# Their back-calculation-from-deaths scenarios differ between the reports,
# since the 18 May report used 88 reported deaths and the 20 May update 131,
# with a corrected set of case-fatality ratios.
# McCabe's scenarios estimate cumulative cases at their report dates, though
# their report is not fully explicit about whether this is symptomatic cases
# or all infections.
# We take the like-for-like quantity to be our cumulative symptom onsets, the
# symptomatic cases, on the same dates, rather than the latent infections
# (which include the not-yet-symptomatic) or our current cut-off total.
# We read our value off the joint fit's cumulative-onset trajectory at the
# grid day for each report date and show it with its credible interval, the
# 18 May report against our 18 May value, the 20 May update against our
# 20 May value, and the 27 May Lancet publication against our 27 May value,
# so each scenario sits beside our estimate for the date it was made.

#md # ```@raw html
#md # <details><summary>McCabe scenarios with uncertainty against our estimates</summary>
#md # ```

function _ci90row(xs)
    (round(Int, quantile(xs, 0.5)),
        round(Int, quantile(xs, 0.05)),
        round(Int, quantile(xs, 0.95)))
end

## Our modelled cumulative symptom onsets on a McCabe report date, read off
## the joint fit's per-draw `cumulative_onsets` trajectory. The grid runs to
## the cut-off on day `n`, so the day-index for a date is `n` minus the days
## from that date back to the cut-off (`grid_day("2026-06-07") = n`,
## `"2026-05-20") = n - 18`, `"2026-05-18") = n - 20`).
_onset_trajs = let mat = chn_joint[:cumulative_onsets]
    [collect(v) for v in vec(collect(mat))]
end
## Inverse of `grid_date(day) = obs.cutoff - Day(obs.n - day)`: the day-index
## whose calendar date is `date`, using `value` (imported above) for the
## day count rather than the non-exported `Dates.date2epochdays`.
_grid_day(date) = obs.n - value(obs.cutoff - Date(date))
function _ours_on(date)
    d = _grid_day(date)
    _ci90row(Float64[t[d] for t in _onset_trajs])
end

## Our matched cumulative-onset estimate for each report date, keyed by date so
## it lands beside that vintage's scenarios in its own panel.
mccabe_ours = Dict(
    "2026-05-18" => _ours_on("2026-05-18"),
    "2026-05-20" => _ours_on("2026-05-20"),
    "2026-05-27" => _ours_on("2026-05-27"))

## One panel per report date; within a panel each method-and-family is one row,
## with the case-fatality / window / doubling-time sweep dodged onto that single
## line, so the ~40 scenarios keep their intervals without becoming ~40 rows.
matched_comparison_fig = plot_scenario_comparison(REPORT_SCENARIOS_CI;
    ours = mccabe_ours,
    date_titles = ["2026-05-18" => "18 May report",
        "2026-05-20" => "20 May update",
        "2026-05-27" => "27 May (Lancet)"],
    xlabel = "Cumulative cases");

#md # ```@raw html
#md # </details>
#md # ```

matched_comparison_fig #hide

# The McCabe scenarios are outbreak-size estimates, the same quantity our
# renewal model and the released integral model report.
# Their 95% confidence intervals come from exact negative-binomial counts
# for the geographic-spread method and a Poisson likelihood profile for the
# back-calculation from deaths.

#md # ```@raw html
#md # <details><summary>Frozen-fit C_T intervals (kept for the CSV export, not shown)</summary>
#md # ```

## The estimate-evolution figure above already shows how the size estimate
## shifts as data accrues, so the side-by-side frozen-fit table is no longer
## rendered in the report; it is kept only to populate the published
## `frozen_matched_cutoffs.csv` export.
frozen_streams_table = streams_table(
    "frozen 20 May" => frozen_C("2026-05-20"),
    "frozen 23 May" => frozen_C("2026-05-23"),
    "frozen 27 May" => frozen_C("2026-05-27"),
    "frozen 8 June" => frozen_C(default_chamla_cutoff()),
    "current data" => posterior_C_joint);

#md # ```@raw html
#md # </details>
#md # ```

# ## Comparison with Chamla et al.
#
# A second group, Chamla et al. [chamla2026](@cite) at the World Health
# Organization Regional Office for Africa, published a stochastic compartmental
# model of the same outbreak on 25 June 2026.
# Their model is a discrete-time susceptible-exposed-infectious-recovered-dead
# ensemble, recalibrated by simulation filtering to the laboratory-confirmed
# case series, anchored on the 598 confirmed cases reported by 8 June, then run
# forward to project the confirmed-case trajectory under a low, central and
# high transmissibility scenario.
#
# Their published quantity is the cumulative confirmed-case count, with the
# reporting fraction held at one, so it does not adjust for the cases that are
# infected but never laboratory-confirmed.
# This is a different quantity from the cumulative cases this analysis and
# McCabe et al. estimate, which include the unconfirmed and unascertained, and
# it sits below them: it is a floor on the true size rather than an estimate of
# it.
# The like-for-like comparison is therefore against our own confirmed-case
# projection, not against our cumulative infection count.
#
# We compare forward projections rather than refitting to their assumptions.
# We take our fit frozen at 8 June, the exact date of their confirmed-case
# calibration anchor, and roll its confirmed-case stream forward to the dates
# Chamla report with the same machinery as the one-week-ahead forecast.
# Setting our projection, their projection and the confirmed cases observed
# since on one timeline shows how each projection has held up against the data.

#md # ```@raw html
#md # <details><summary>Project the 8 June fit forward and assemble the Chamla comparison</summary>
#md # ```

## The 8 June frozen joint fit matches Chamla's confirmed-case calibration
## anchor exactly and carries the confirmed-case testing history through then,
## so we roll its confirmed-case stream forward with the one-week-ahead forecast
## machinery to the dates Chamla report.
chamla_anchor = frozen_by_cutoff["2026-06-08"]

## Our projected cumulative confirmed cases at a horizon of `h` days past the
## 8 June cut-off: a forward `forecast_reported` run (its reproduction number
## left to keep evolving), summarised as (median, 5%, 95%).
function _our_confirmed_h(h)
    fc = forecast_reported(chamla_anchor.chn;
        horizon = h,
        obs_cases = chamla_anchor.o.reported_cases,
        obs_deaths = chamla_anchor.o.total_deaths,
        obs_confirmed = chamla_anchor.o.confirmed_cases,
        obs_confirmed_deaths = chamla_anchor.o.confirmed_deaths)
    return _ci90row(float.(fc.confirmed_cum))
end

## Our projection at Chamla's forward report dates (10 and 24 June, week 12):
## the anchor day is the fitted confirmed total at 8 June, each later date a
## forward forecast. Reused for the matched-date table and the week-12 figure.
chamla_fan = map(["2026-06-08", "2026-06-10", "2026-06-24"]) do d
    h = value(Date(d) - chamla_anchor.cutoff)
    row = h == 0 ?
          (chamla_anchor.o.confirmed_cases, chamla_anchor.o.confirmed_cases,
        chamla_anchor.o.confirmed_cases) : _our_confirmed_h(h)
    (d, row...)
end
_fan_at(date) =
    let r = first(x for x in chamla_fan if x[1] == date)
        (r[2], r[3], r[4])
    end
ours_10jun = _fan_at("2026-06-10")
ours_24jun = _fan_at("2026-06-24")

## Observed confirmed cases over the comparison window: the daily cumulative
## series read off the chain's grid from 18 May (Chamla's first projected point)
## to the cut-off.
chamla_obs_series = let
    ds = [grid_date(d) for d in obs.confirmed_history.days]
    cs = obs.confirmed_history.counts
    [(string(ds[i]), cs[i]) for i in eachindex(ds) if ds[i] >= Date("2026-05-18")]
end

## Chamla's central confirmed-case projection over the comparison window; their
## later, far-larger horizons are noted in the text rather than plotted so the
## window stays legible.
chamla_central_window = CHAMLA_CONFIRMED_CENTRAL[1:4]

chamla_projection_fig = plot_projection_comparison(;
    external = chamla_central_window,
    ours = chamla_fan,
    observed = chamla_obs_series,
    external_label = "Chamla et al. central (R₀=1.71)",
    ours_label = "Our projection (from 8 June)",
    observed_label = "Observed confirmed",
    title = "Confirmed-case projections versus observed, from mid-May");

#md # ```@raw html
#md # </details>
#md # ```

chamla_projection_fig #hide

# By 24 June their central scenario projected just under a thousand confirmed
# cases, and their low and high scenarios ranged from roughly 870 to 1360.
# The figure below sets that week-12 scenario spread beside our 8 June
# projection for the same date and the confirmed count observed by the cut-off,
# so each reads against their three scenarios at a glance.

#md # ```@raw html
#md # <details><summary>Week-12 (24 June) scenario spread against ours and observed</summary>
#md # ```

chamla_w12_rows = vcat(
    [(label, m, lo, hi) for (label, m, lo, hi) in CHAMLA_CONFIRMED_W12],
    [("Our projection (from 8 June)", ours_24jun...)],
    [("Observed by 23 June cut-off", obs.confirmed_cases,
        obs.confirmed_cases, obs.confirmed_cases)])
chamla_w12_groups = vcat(fill("Chamla et al. scenarios", 3),
    ["Our projection"], ["Observed"])

chamla_w12_fig = plot_estimate_comparison(chamla_w12_rows;
    xlabel = "Cumulative confirmed cases by 24 June",
    groups = chamla_w12_groups,
    group_colours = ["Chamla et al. scenarios" => :steelblue,
        "Our projection" => :firebrick,
        "Observed" => :black]);

#md # ```@raw html
#md # </details>
#md # ```

chamla_w12_fig #hide

# The matched-date numbers behind these figures are in the dropdown below, with
# the observed column taken to the 23 June cut-off.

#md # ```@raw html
#md # <details><summary>Matched-date projection numbers (10 and 24 June)</summary>
#md # ```

chamla_comparison_table = let
    fmt(t) = string(t[1], " (", t[2], "–", t[3], ")")
    central(date) =
        let r = first(x for x in CHAMLA_CONFIRMED_CENTRAL
            if x[1] == date)
            fmt((r[2], r[3], r[4]))
        end
    DataFrame(
        "Date" => ["10 June", "24 June"],
        "Chamla central (90% PI)" => [central("2026-06-10"),
            central("2026-06-24")],
        "Our projection (90% CrI)" => [fmt(ours_10jun), fmt(ours_24jun)],
        "Observed confirmed" => [
            string(freeze_observations("2026-06-10").confirmed_cases),
            string(obs.confirmed_cases) * " (23 June)"])
end;

chamla_comparison_table #hide

#md # ```@raw html
#md # </details>
#md # ```

# Beyond the comparison window their central scenario continues to roughly 8200
# confirmed cases by mid-September, with the high scenario far higher; those
# longer projections are not set against data here.
# Their confirmed-case projections are a floor on the outbreak size, so they sit
# below our cumulative infection count, which adds the unconfirmed and
# unascertained cases on top.

# ## Reproduction number behind the projection
#
# The forward projection above is carried by the reproduction-number trajectory
# our 8 June fit estimated, a quantity we report in its own right rather than as
# a comparison.
# The figure shows that trajectory, the time-varying reproduction number from
# the renewal walk with its credible intervals, as the fit saw it at 8 June.
# It declines over the weeks leading to the cut-off, and that decline is what
# bends the projected trajectory away from sustained early growth.

#md # ```@raw html
#md # <details><summary>Reproduction number as estimated by the 8 June fit</summary>
#md # ```

## Reconstruct the reproduction-number trajectory the 8 June fit estimated,
## mirroring the current-data R_t figure but with the frozen vintage's own grid,
## breakpoint and renewal start.
chamla_rt_obs = chamla_anchor.o
chamla_rt_breakpoint = chamla_rt_obs.n - chamla_rt_obs.who_first_sitrep_days
chamla_rt_start = clamp(
    chamla_rt_obs.n - round(Int, chamla_rt_obs.tmrca_days) + RENEWAL_START_LEAD,
    1, chamla_rt_obs.n)
chamla_rt_fig = plot_rt(chamla_anchor.chn;
    n = chamla_rt_obs.n, breakpoint = chamla_rt_breakpoint,
    rt_start = chamla_rt_start,
    rt_walk_start = clamp(chamla_rt_breakpoint - RT_WALK_LEAD,
        chamla_rt_start, chamla_rt_obs.n),
    as_of_date = string(chamla_rt_obs.cutoff),
    seeding = chamla_rt_obs.seeding, ramp = RT_INTERVENTION_RAMP);

#md # ```@raw html
#md # </details>
#md # ```

chamla_rt_fig #hide

# ## Delay sensitivity
#
# The death stream dates the outbreak from how far deaths lag symptom onset,
# so the assumed onset-to-death delay sets the implied infection count.
# The baseline uses the hospital-pathway delay from the Isiro 2012 line-list
# reanalysis (onset to admission then admission to death, implied mean about
# 12 d).
# We re-fit the joint model under the community-pathway delay from the same
# reanalysis, the delay for deaths that occur in the community without a
# recorded admission, which is shorter (implied mean about 8 d).
# Both pathways come from the line list, so this varies the actual delay
# assumption rather than an arbitrary scenario.
# The re-fit uses the full headline settings (1000 draws across two chains).
#
# The infection count to date shifts with the assumed delay, and the
# table and overlaid densities below show how far.

#md # ```@raw html
#md # <details><summary>Re-fit the joint under the community-pathway onset-to-death delay</summary>
#md # ```

## The sensitivity re-fits (community-delay variant) are
## defined in the fit registry (`docs/fits/registry.jl`) and loaded through the cache
## (when enabled) in the setup block above.
posterior_C_community_delay = RUN_SENSITIVITY ?
                              vec(Array(chn_joint_community_delay[:C_T])) : nothing

#md # ```@raw html
#md # </details>
#md # ```

#md # ```@raw html
#md # <details><summary>Delay-sensitivity infection-count table</summary>
#md # ```

delay_sensitivity_table = RUN_SENSITIVITY ?
                          streams_table("baseline (hospital pathway)" => posterior_C_joint,
    "community pathway" => posterior_C_community_delay) :
                          Markdown.md"_Delay sensitivity analysis not shown in this build._"

#md # ```@raw html
#md # </details>
#md # ```

delay_sensitivity_table #hide

#md # ```@raw html
#md # <details><summary>Delay-sensitivity infection-count density plot</summary>
#md # ```

delay_sensitivity_fig = RUN_SENSITIVITY ?
                        plot_cumulative_cases(
    "baseline (hospital pathway)" => posterior_C_joint,
    "community pathway" => posterior_C_community_delay; scenarios = []) :
                        Markdown.md"_Delay sensitivity analysis not shown in this build._"

#md # ```@raw html
#md # </details>
#md # ```

delay_sensitivity_fig #hide

# ## Tree-prior sensitivity
#
# The outbreak-age estimate depends on the coalescent tree prior assumed
# in the BEAST X analysis. The baseline uses the more flexible Skygrid
# non-parametric model, which dates the common ancestor to 15 March 2026
# ($95\%$ HPD 09 Feb -- 12 Apr). The report also fits an Exponential
# growth tree prior, which dates the common ancestor about a week earlier
# to 08 March 2026 ($95\%$ HPD 01 Feb -- 05 Apr)
# [mbalaplacide2026](@cite). Both priors give very similar evolutionary
# rates ($\sim 1.1\times10^{-3}$ subs/site/year). We re-fit the joint
# model under the Exponential growth TMRCA and compare the infection
# count to date and the outbreak age.

#md # ```@raw html
#md # <details><summary>Re-fit the joint under the Exponential growth tree prior</summary>
#md # ```

## The Exponential-growth re-fit (and its `tmrca_days` offset) is defined in the fit
## registry (`docs/fits/registry.jl`) and loaded through the cache (when enabled) in the
## setup block above.
posterior_C_exp_growth = RUN_SENSITIVITY ?
                         vec(Array(chn_joint_exp_growth_clock[:C_T])) : nothing
T_skygrid = vec(Array(chn_joint[:T]))
T_exp_growth = RUN_SENSITIVITY ? vec(Array(chn_joint_exp_growth_clock[:T])) : nothing

#md # ```@raw html
#md # </details>
#md # ```

# The infection count to date under the two tree priors, side by side.
# A slightly earlier common ancestor (Exponential growth) permits a
# marginally older outbreak, though the difference is small because the
# evolutionary rates are nearly identical.

#md # ```@raw html
#md # <details><summary>Tree-prior infection-count table</summary>
#md # ```

clock_sensitivity_C_table = RUN_SENSITIVITY ?
                            streams_table("Skygrid (baseline)" => posterior_C_joint,
    "Exponential growth" => posterior_C_exp_growth) :
                            Markdown.md"_Tree-prior sensitivity analysis not shown in this build._"

#md # ```@raw html
#md # </details>
#md # ```

clock_sensitivity_C_table #hide

#md # ```@raw html
#md # <details><summary>Tree-prior infection-count density plot</summary>
#md # ```

clock_sensitivity_C_fig = RUN_SENSITIVITY ?
                          plot_cumulative_cases("Skygrid (baseline)" => posterior_C_joint,
    "Exponential growth" => posterior_C_exp_growth; scenarios = []) :
                          Markdown.md"_Tree-prior sensitivity analysis not shown in this build._"

#md # ```@raw html
#md # </details>
#md # ```

clock_sensitivity_C_fig #hide

# The outbreak age, the number of days from seeding to the cut-off, under
# the two tree priors.

#md # ```@raw html
#md # <details><summary>Tree-prior outbreak-age table</summary>
#md # ```

clock_sensitivity_T_table = RUN_SENSITIVITY ?
                            streams_table("Skygrid (baseline)" => T_skygrid,
    "Exponential growth" => T_exp_growth; digits = 0) :
                            Markdown.md"_Tree-prior sensitivity analysis not shown in this build._"

#md # ```@raw html
#md # </details>
#md # ```

clock_sensitivity_T_table #hide

#md # ```@raw html
#md # <details><summary>Tree-prior outbreak-age density plot</summary>
#md # ```

clock_sensitivity_T_fig = RUN_SENSITIVITY ?
                          plot_density_overlay("Skygrid (baseline)" => T_skygrid,
    "Exponential growth" => T_exp_growth;
    xlabel = "Outbreak age (days before cut-off)",
    title = "Posterior outbreak age by tree prior") :
                          Markdown.md"_Tree-prior sensitivity analysis not shown in this build._"

#md # ```@raw html
#md # </details>
#md # ```

clock_sensitivity_T_fig #hide

# ## Saving sensitivity results
#
# The stream-comparison and frozen-fit tables and the per-stream reproduction
# number figure are written to the shared `output/` directory (the main
# analysis writes the rest), so the combined release and summary dashboard
# pick up both pages' outputs.

#md # ```@raw html
#md # <details><summary>Write sensitivity outputs</summary>
#md # ```

output_dir = get(ENV, "BVD_OUTPUT_DIR",
    joinpath(pkgdir(BVDOutbreakSize), "output"))
mkpath(output_dir)
CSV.write(joinpath(output_dir, "cumulative_cases_by_stream.csv"),
    streams_C_table)
CSV.write(joinpath(output_dir, "frozen_matched_cutoffs.csv"),
    frozen_streams_table)

## The one-week-back validation forecast, in the same archive format as the
## release forecast, so the frozen "last week versus now" forecast is recorded
## as a release asset alongside the forecast it is scored against.
CSV.write(joinpath(output_dir, "forecast_validation.csv"),
    forecast_archive([(7, validation_forecast)];
        made_date = frozen_lastweek.o.cutoff, thin = 5))

## The per-stream reproduction-number figure for the summary dashboard; the
## main analysis writes the other three dashboard figures.
dashboard_dir = joinpath(
    pkgdir(BVDOutbreakSize), "docs", "src", "summary_assets")
mkpath(dashboard_dir)
CairoMakie.save(joinpath(dashboard_dir, "rt_streams.png"), stream_rt_fig)

#md # ```@raw html
#md # </details>
#md # ```
