## Tests for the shared stream registry in src/data.jl: which streams the
## observation loader carries, when each last reported, and whether it was
## still reporting at the cut-off. This is the one rule the report pages
## and the release scorer partition streams on.

@testitem "stream_id resolves every stream vocabulary" begin
    using BVDOutbreakSize: stream_id
    ## The canonical id, the observation-set field, the scoring label and
    ## the forecast column names all name the same stream.
    @test stream_id(:confirmed_cases) == :confirmed_cases
    @test stream_id(:confirmed_history) == :confirmed_cases
    @test stream_id("confirmed cases") == :confirmed_cases
    @test stream_id(:confirmed_cum) == :confirmed_cases
    @test stream_id(:confirmed_new) == :confirmed_cases
    @test stream_id(:cases_cum) == :suspected_cases
    @test stream_id("isolation beds (suspected)") == :suspect_beds
    @test stream_id(:export_case_days) == :exports
    err = try
        stream_id(:not_a_stream)
    catch e
        e
    end
    @test err isa ErrorException
    @test occursin("confirmed_cases", err.msg)
end

@testitem "stream_forecast_columns names the forecast columns" begin
    using BVDOutbreakSize: stream_forecast_columns
    ## Derived from the registry's forecast_prefix, so the column names are
    ## not repeated per consumer, and resolvable from any vocabulary.
    @test stream_forecast_columns(:confirmed_cases) ==
          (cum = :confirmed_cum, new = :confirmed_new)
    @test stream_forecast_columns(:confirmed_cum) ==
          (cum = :confirmed_cum, new = :confirmed_new)
    @test stream_forecast_columns("suspected deaths") ==
          (cum = :deaths_cum, new = :deaths_new)
    ## A stream the forecast does not carry has no columns.
    @test isnothing(stream_forecast_columns(:tests_analysed))
end

@testitem "stream_last_date is the stream's last vintage" begin
    using Dates: Date, Day
    using BVDOutbreakSize: stream_last_date
    n = 40
    cutoff = Date(2026, 7, 15)
    obs = (; cutoff = cutoff, n = n,
        reported_history = (; days = [5, 12], counts = [50.0, 90.0]),
        deaths_history = (; days = Int[], counts = Float64[]))
    @test stream_last_date(obs, :suspected_cases) == cutoff - Day(n - 12)
    ## A history with no vintages, and a stream this observation set does
    ## not carry at all, both have no last reported date.
    @test ismissing(stream_last_date(obs, :suspected_deaths))
    @test ismissing(stream_last_date(obs, :confirmed_cases))
end

@testitem "stream_last_date reads exports from the dated detections" begin
    using Dates: Date, Day
    using BVDOutbreakSize: stream_last_date
    n = 40
    cutoff = Date(2026, 7, 15)
    ## The exports are a dated list of detections rather than a series of
    ## vintages, and the later of the import and import-death detections
    ## is the last thing that stream reported.
    obs = (; cutoff = cutoff, n = n,
        export_case_days = [8, 15, 20], export_death_days = [9, 22])
    @test stream_last_date(obs, :exports) == cutoff - Day(n - 22)
end

@testitem "stream_reporting turns over at the grace boundary" begin
    using Dates: Date, Day
    using BVDOutbreakSize: stream_reporting, STREAM_REPORTING_GRACE_DAYS
    n = 40
    cutoff = Date(2026, 7, 15)
    grace = STREAM_REPORTING_GRACE_DAYS
    mk(day) = (; cutoff = cutoff, n = n,
        reported_history = (; days = [day], counts = [90.0]))
    @test stream_reporting(mk(n), :suspected_cases)
    @test stream_reporting(mk(n - grace), :suspected_cases)
    @test !stream_reporting(mk(n - grace - 1), :suspected_cases)
    ## A stream with no vintages is not reporting.
    empty_obs = (; cutoff = cutoff, n = n,
        reported_history = (; days = Int[], counts = Float64[]))
    @test !stream_reporting(empty_obs, :suspected_cases)
end

@testitem "stream_report_status covers the loaded manifest" begin
    using DataFrames: nrow
    using BVDOutbreakSize: load_observations, stream_report_status,
                           OBSERVATION_STREAMS
    obs = load_observations()
    status = stream_report_status(obs)
    ## Every registry entry names a field the loader really returns, so a
    ## rename in either place fails here rather than silently dropping a
    ## stream.
    @test all(e -> hasproperty(obs, e.field), OBSERVATION_STREAMS)
    @test nrow(status) == length(OBSERVATION_STREAMS)
    @test names(status) ==
          ["stream", "label", "last_date", "reporting", "days_since"]
    @test all(d -> d <= obs.cutoff, skipmissing(status.last_date))
    @test maximum(skipmissing(status.last_date)) == obs.cutoff
    ## The split the report's posterior predictive checks partition on.
    reporting = Set(status.stream[status.reporting])
    for s in (:confirmed_cases, :confirmed_deaths, :recovered,
        :isolation_beds, :tests_analysed_daily, :onset_reports)
        @test s in reporting
    end
    for s in (:suspected_cases, :suspected_deaths, :suspected_daily,
        :suspected_daily_deaths, :tests_analysed, :treatment_admissions,
        :treatment_deaths, :treatment_ruleouts, :treatment_absconded,
        :treatment_beds, :suspect_beds)
        @test !(s in reporting)
    end
end

@testitem "grid_date maps the grid ends to the seeding and cut-off" begin
    using BVDOutbreakSize: load_observations
    using BVDOutbreakSize: grid_date
    obs = load_observations()
    @test grid_date(obs, obs.n) == obs.cutoff
    @test grid_date(obs, 1) == obs.seeding
end
