## Smoke tests for the one-week-ahead forecast. Builds a tiny
## synthetic chain carrying the parameters `forecast_reported` reads,
## then checks the returned DataFrame contract.

@testsnippet ForecastFixtures begin
    using Turing: @model, sample, Prior
    using Distributions: Beta, Normal, truncated
    import FlexiChains
    using BVDOutbreakSize: deaths_only_model, bvd_joint

    ## Synthetic prior carrying every parameter name that
    ## `forecast_reported` reads: :r, :expected_reports_T,
    ## :expected_deaths_T, :expected_infections_T, :R_T,
    ## :expected_confirmed_T, :expected_confirmed_deaths_T, :k.
    @model function _forecast_test()
        r ~ truncated(Normal(0.05, 0.01); lower = 1e-3)
        inv_sqrt_k ~ truncated(Normal(0.5, 0.2); lower = 1e-3)
        k := 1.0 / (inv_sqrt_k^2 + eps(typeof(inv_sqrt_k)))
        expected_reports_T ~ truncated(Normal(300.0, 50.0); lower = 1.0)
        expected_deaths_T ~ truncated(Normal(15.0, 3.0); lower = 1.0)
        expected_infections_T ~ truncated(Normal(800.0, 100.0); lower = 1.0)
        R_T ~ truncated(Normal(1.5, 0.3); lower = 1e-3)
        expected_confirmed_T ~ truncated(Normal(120.0, 20.0); lower = 1.0)
        expected_confirmed_deaths_T ~ truncated(Normal(8.0, 2.0); lower = 0.5)
        return nothing
    end

    _forecast_chain(n) = sample(
        _forecast_test(), Prior(), n;
        chain_type = FlexiChains.VNChain, progress = false
    )
end

@testitem "forecast_reported returns the documented columns" tags=[:slow] setup=[ForecastFixtures] begin
    using DataFrames: DataFrame, nrow
    using BVDOutbreakSize: forecast_reported

    chn=_forecast_chain(200)
    fc=forecast_reported(chn;
        horizon = 7,
        obs_cases = 905,
        obs_deaths = 18,
        obs_confirmed = 210,
        obs_confirmed_deaths = 17)

    @test fc isa DataFrame
    @test nrow(fc) == 200
    cols=[:cases_cum, :deaths_cum, :confirmed_cum, :confirmed_deaths_cum,
        :cases_new, :deaths_new, :confirmed_new, :confirmed_deaths_new,
        :infections_new, :rt_forecast]
    @test all(c -> c in propertynames(fc), cols)
    @test all(fc.infections_new .>= 0)
    @test all(fc.cases_cum .>= 0)
    @test all(fc.deaths_cum .>= 0)
    @test all(fc.confirmed_cum .>= 0)
    @test all(fc.confirmed_deaths_cum .>= 0)
    @test all(fc.cases_new .>= 0)
    @test all(fc.deaths_new .>= 0)
    ## New-this-week cannot exceed the cumulative forecast.
    @test all(fc.cases_new .<= fc.cases_cum)
    @test all(fc.deaths_new .<= fc.deaths_cum)
    @test all(fc.confirmed_new .<= fc.confirmed_cum)
    ## Confirmed deaths are a thinning of suspected deaths, so cannot exceed
    ## the forecast cumulative suspected deaths.
    @test all(fc.confirmed_deaths_cum .<= fc.deaths_cum)
end

@testitem "forecast_table has expected rows and columns" tags=[:slow] setup=[ForecastFixtures] begin
    using DataFrames: DataFrame, nrow
    using BVDOutbreakSize: forecast_reported, forecast_table

    chn=_forecast_chain(200)
    fc=forecast_reported(chn;
        horizon = 7,
        obs_cases = 905,
        obs_deaths = 18,
        obs_confirmed = 210,
        obs_confirmed_deaths = 17)

    tbl=forecast_table(fc)
    @test tbl isa DataFrame
    ## Four streams x two quantities (cumulative, new this week).
    @test nrow(tbl) == 8
    @test names(tbl) ==
          ["Stream", "Quantity", "Lower 90%", "Lower 60%", "Lower 30%",
        "Upper 30%", "Upper 60%", "Upper 90%"]
    @test Set(tbl[!, "Quantity"]) ==
          Set(["cumulative by T+7", "new this week"])
end

@testitem "forecast_vs_truth compares forecast to observed counts" tags=[:slow] setup=[ForecastFixtures] begin
    using DataFrames: DataFrame, nrow
    using BVDOutbreakSize: forecast_reported, forecast_vs_truth

    chn=_forecast_chain(200)
    fc=forecast_reported(chn;
        horizon = 7,
        obs_cases = 905,
        obs_deaths = 18,
        obs_confirmed = 210,
        obs_confirmed_deaths = 17)

    tbl=forecast_vs_truth(fc;
        cases = 1000, deaths = 25, confirmed = 260, confirmed_deaths = 20)

    @test tbl isa DataFrame
    @test nrow(tbl) == 4
    @test names(tbl) ==
          ["Stream", "Observed", "Lower 90%", "Lower 60%", "Lower 30%",
        "Upper 30%", "Upper 60%", "Upper 90%", "Within 90% PI"]
    @test Set(tbl[!, "Stream"]) ==
          Set(["DRC reported cases", "DRC deaths",
        "DRC confirmed cases", "DRC confirmed deaths"])

    for row in eachrow(tbl)
        covered=row["Lower 90%"]<=row.Observed<=row["Upper 90%"]
        @test row["Within 90% PI"] == (covered ? "yes" : "no")
    end
end
