## Tests for `confirmed_break_correction` and its use in the in-report
## forecast validation. On a listed break day the reported cumulative jumps
## by more than that day's notifications, so a raw cumulative difference is
## not the count that was notified across the week.
##
## `test_break_day_correction.jl` covers the same quantity as the release
## scoring reads it, through `scripts/score_releases.jl`'s stream dispatch
## onto `confirmed_break_steps`.

@testitem "confirmed_break_correction sums net minus gross in the window" begin
    using BVDOutbreakSize: confirmed_break_correction

    ## A two-vintage history with one break day at day 5. The step is 100
    ## against a printed 24h count of 30, so 70 is retrospective.
    obs = (confirmed_break_days = [5],
        confirmed_break_gross_cases = [30],
        confirmed_break_gross_deaths = [4],
        confirmed_history = (days = [3, 5], counts = [200, 300]),
        confirmed_deaths_history = (days = [3, 5], counts = [50, 70]))

    @test confirmed_break_correction(obs, 4, 6) == 70.0
    @test confirmed_break_correction(obs, 4, 6; deaths = true) == 16.0

    ## The window is half open on the left, so a break day on the origin
    ## belongs to the previous window, not this one.
    @test confirmed_break_correction(obs, 5, 9) == 0.0
    @test confirmed_break_correction(obs, 1, 4) == 0.0
    @test confirmed_break_correction(obs, 1, 5) == 70.0

    ## A break day that fell before the validation window contributes
    ## nothing to it, for cases and for deaths alike, so a week holding no
    ## harmonisation is scored on its raw counts.
    @test confirmed_break_correction(obs, 6, 9) == 0.0
    @test confirmed_break_correction(obs, 6, 9; deaths = true) == 0.0
end

@testitem "confirmed_break_correction is zero without a listed break" begin
    using BVDOutbreakSize: confirmed_break_correction

    none = (confirmed_break_days = Int[],
        confirmed_break_gross_cases = Int[],
        confirmed_break_gross_deaths = Int[],
        confirmed_history = (days = [3, 5], counts = [200, 300]),
        confirmed_deaths_history = (days = [3, 5], counts = [50, 70]))
    @test confirmed_break_correction(none, 0, 99) == 0.0

    ## A break day the history carries no vintage for has not arrived by this
    ## cut-off, so it contributes nothing rather than erroring.
    unarrived = (confirmed_break_days = [7],
        confirmed_break_gross_cases = [30],
        confirmed_break_gross_deaths = [4],
        confirmed_history = (days = [3, 5], counts = [200, 300]),
        confirmed_deaths_history = (days = [3, 5], counts = [50, 70]))
    @test confirmed_break_correction(unarrived, 0, 99) == 0.0

    ## An empty history has no step to correct.
    empty_hist = (confirmed_break_days = [5],
        confirmed_break_gross_cases = [30],
        confirmed_break_gross_deaths = [4],
        confirmed_history = (days = Int[], counts = Int[]),
        confirmed_deaths_history = (days = Int[], counts = Int[]))
    @test confirmed_break_correction(empty_hist, 0, 99) == 0.0
end

@testitem "forecast_vs_truth subtracts the break from both truths" begin
    using DataFrames: DataFrame
    using BVDOutbreakSize: forecast_vs_truth

    fc = DataFrame(confirmed_cum = collect(400.0:1.0:499.0),
        confirmed_new = collect(100.0:1.0:199.0))
    observed = (confirmed_cum = 500.0,)
    baseline = (confirmed_cum = 300.0,)

    plain = forecast_vs_truth(fc; observed = observed, baseline = baseline)
    corrected = forecast_vs_truth(fc; observed = observed,
        baseline = baseline, breaks = (confirmed_cum = 70.0,))

    newrow(df) = df[df[!, "Quantity"] .== "new this week", "Observed"][1]
    cumrow(df) = df[df[!, "Quantity"] .== "cumulative by T+7", "Observed"][1]

    @test newrow(plain) == 200.0
    @test newrow(corrected) == 130.0
    ## The projected cumulative is the origin plus the projected new count,
    ## so it cannot contain a reattachment either. Both truths drop it, and
    ## the new count stays the difference between them and the origin.
    @test cumrow(plain) == 500.0
    @test cumrow(corrected) == 430.0
    @test cumrow(corrected) - 300.0 == newrow(corrected)
end

@testitem "plot_forecast_vs_truth takes the break out of both truths" begin
    using DataFrames: DataFrame
    using CairoMakie: Axis
    using BVDOutbreakSize: plot_forecast_vs_truth

    fc = DataFrame(confirmed_cum = collect(400.0:1.0:499.0),
        confirmed_new = collect(100.0:1.0:199.0))
    indiv = (confirmed_new = collect(90.0:1.0:189.0),)
    origin = (confirmed_cum = 300.0,)

    plain = plot_forecast_vs_truth(fc; observed = (confirmed_cum = 500.0,),
        baseline = origin, individual = indiv)
    corrected = plot_forecast_vs_truth(fc;
        observed = (confirmed_cum = 500.0,), baseline = origin,
        individual = indiv, breaks = (confirmed_cum = 70.0,))
    ## A break is exactly a reduction of the reported cumulative, so passing
    ## it matches passing a total already net of it. The overlay's origin is
    ## the frozen baseline either way and does not move with it.
    presubtracted = plot_forecast_vs_truth(fc;
        observed = (confirmed_cum = 430.0,), baseline = origin,
        individual = indiv)

    cum_limits(fig) = [c for c in fig.content if c isa Axis][1].limits[]
    @test cum_limits(corrected) == cum_limits(presubtracted)
    ## The cumulative panel does move, which is the change: it used to score
    ## the reported total against a projection that could not contain the
    ## harmonisation.
    @test cum_limits(plain) != cum_limits(corrected)
end
