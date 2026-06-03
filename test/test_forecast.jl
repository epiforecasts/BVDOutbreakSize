## Smoke tests for the one-week-ahead forecast. Builds a tiny
## synthetic chain carrying the parameters `forecast_reported` reads,
## then checks the returned DataFrame contract.

@testsnippet ForecastFixtures begin
    using Turing: Turing, @model, sample, Prior
    using Distributions: Beta, Normal, truncated
    import FlexiChains
    using BVDOutbreakSize: bvd_joint

    ## Synthetic prior carrying every parameter name that
    ## `forecast_reported` reads. `include_lab = true` adds the
    ## lab-turnaround delay and PCR sensitivity draws so the
    ## confirmed-cases columns are populated.
    @model function _forecast_test(; include_lab::Bool = false)
        r ~ truncated(Normal(0.05, 0.01); lower = 1e-3)
        T ~ truncated(Normal(100.0, 10.0); lower = 1.0)
        CFR ~ Beta(6.0, 14.0)
        α ~ truncated(Normal(4.3, 0.5); lower = 0.5)
        θ ~ truncated(Normal(2.6, 0.3); lower = 0.2)
        w ~ truncated(Normal(15.0, 2.0); lower = 1.0)
        p_drc ~ Beta(2.0, 6.0)
        p_uganda ~ Beta(2.0, 6.0)
        inv_sqrt_k ~ truncated(Normal(0.0, 1.0); lower = 1e-3)
        k := 1.0 / (inv_sqrt_k^2 + eps(typeof(inv_sqrt_k)))
        α_rep ~ truncated(Normal(4.0, 0.5); lower = 0.5)
        θ_rep ~ truncated(Normal(3.0, 0.3); lower = 0.2)
        λ_bg ~ truncated(Normal(0.0, 10.0); lower = 0)
        if include_lab
            ## Severe-first confirmed-stream draws the forecast reads:
            ## sensitivity, specificity, the q-curve shape (q0, qinf,
            ## decay) and the forwarded fraction.
            s_test ~ Beta(15.0, 2.0)
            spec_test ~ Beta(50.0, 1.5)
            q0 ~ Beta(20.0, 1.5)
            qinf ~ Beta(6.0, 6.0)
            decay_scale ~ truncated(Normal(0.0, 10.0); lower = 0)
            τ_forward ~ Beta(5.0, 2.0)
            ## Receipt-delay Gamma and the cut-off analysis capacity the
            ## capacity-limited confirmed forecast reads.
            α_recv ~ truncated(Normal(2.0, 1.0); lower = 0.1)
            θ_recv ~ truncated(Normal(1.5, 0.75); lower = 0.1)
            capacity_cutoff ~ truncated(Normal(150.0, 30.0); lower = 1.0)
        end
        return nothing
    end

    _forecast_chain(n;
        include_lab::Bool = false) = sample(_forecast_test(; include_lab), Prior(), n;
        chain_type = FlexiChains.VNChain, progress = false)
end

@testitem "forecast_reported returns the documented columns" tags=[:slow] setup=[ForecastFixtures] begin
    using DataFrames: DataFrame, nrow
    using BVDOutbreakSize: forecast_reported

    chn=_forecast_chain(200)

    fc=forecast_reported(chn;
        horizon = 7,
        daily_travellers = 1871,
        source_population = 4_392_200,
        obs_cases = 514,
        obs_deaths = 136,
        obs_exports = 2)

    @test fc isa DataFrame
    @test nrow(fc) == 200
    cols=[:cases_cum, :deaths_cum, :exports_cum,
        :cases_new, :deaths_new, :exports_new]
    @test all(c -> c in propertynames(fc), cols)
    ## Counts are non-negative integers.
    @test all(fc.cases_cum .>= 0)
    @test all(fc.deaths_cum .>= 0)
    @test all(fc.exports_cum .>= 0)
    @test all(fc.cases_new .>= 0)
    ## New-this-week cannot exceed the cumulative forecast.
    @test all(fc.cases_new .<= fc.cases_cum)
    @test all(fc.deaths_new .<= fc.deaths_cum)
    ## No lab-turnaround draws in this fixture → no confirmed columns.
    @test !(:confirmed_cum in propertynames(fc))
end

@testitem "forecast_reported adds confirmed columns when lab delay sampled" tags=[:slow] setup=[ForecastFixtures] begin
    using DataFrames: DataFrame
    using BVDOutbreakSize: forecast_reported

    chn=_forecast_chain(200; include_lab = true)
    fc=forecast_reported(chn;
        horizon = 7, daily_travellers = 1871,
        source_population = 4_392_200,
        obs_cases = 514, obs_deaths = 136, obs_exports = 2,
        obs_confirmed = 33, obs_analysed = 211)

    @test :confirmed_cum in propertynames(fc)
    @test :confirmed_new in propertynames(fc)
    @test all(fc.confirmed_cum .>= 0)
    @test all(fc.confirmed_new .>= 0)
    @test all(fc.confirmed_new .<= fc.confirmed_cum)
end

@testitem "forecast_table and plot_forecast" tags=[:slow] setup=[
    ForecastFixtures, HeadlessMakie] begin
    using DataFrames: DataFrame, nrow
    using BVDOutbreakSize: forecast_reported, forecast_table, plot_forecast

    chn=_forecast_chain(200)
    fc=forecast_reported(chn;
        horizon = 7, daily_travellers = 1871,
        source_population = 4_392_200,
        obs_cases = 514, obs_deaths = 136, obs_exports = 2)

    tbl=forecast_table(fc)
    @test tbl isa DataFrame
    ## Three streams × two quantities (cumulative, new this week).
    @test nrow(tbl) == 6
    @test names(tbl) ==
          ["Stream", "Quantity", "Lower 90%", "Lower 60%", "Lower 30%",
        "Upper 30%", "Upper 60%", "Upper 90%"]
    @test Set(tbl[!, "Quantity"]) ==
          Set(["cumulative by T+7", "new this week"])

    fig=plot_forecast(fc)
    @test fig !== nothing
end

@testitem "forecast_vs_truth compares cumulative forecast to observed" tags=[:slow] setup=[
    ForecastFixtures, HeadlessMakie] begin
    using DataFrames: DataFrame, nrow
    using BVDOutbreakSize: forecast_reported, forecast_vs_truth, plot_forecast_vs_truth

    chn=_forecast_chain(200)
    fc=forecast_reported(chn;
        horizon = 7, daily_travellers = 1871,
        source_population = 4_392_200,
        obs_cases = 514, obs_deaths = 136, obs_exports = 2)

    tbl=forecast_vs_truth(fc;
        cases = 600, deaths = 150, exports = 3)

    @test tbl isa DataFrame
    ## One row per stream (cases, deaths, exports).
    @test nrow(tbl) == 3
    @test names(tbl) ==
          ["Stream", "Observed", "Lower 90%", "Lower 60%", "Lower 30%",
        "Upper 30%", "Upper 60%", "Upper 90%", "Within 90% PI"]
    @test Set(tbl[!, "Stream"]) ==
          Set(["DRC reported cases", "DRC deaths", "Uganda exports"])

    ## Coverage flag agrees with the reported interval endpoints.
    for row in eachrow(tbl)
        covered=row["Lower 90%"]<=row.Observed<=row["Upper 90%"]
        @test row["Within 90% PI"] == (covered ? "yes" : "no")
    end

    fig=plot_forecast_vs_truth(fc; cases = 600, deaths = 150, exports = 3)
    @test fig !== nothing
end

@testitem "forecast BVD means apply the incubation onset_fraction" begin
    using Distributions: Gamma
    using BVDOutbreakSize: onset_rescale
    import BVDOutbreakSize as B

    r = 0.05
    Th = 107.0
    α, θ, CFR = 4.3, 2.6, 0.3
    α_rep, θ_rep, p_drc, λ_bg = 4.0, 3.0, 0.2, 1.5
    s_test, spec_test = 0.9, 0.97
    q0, qinf, decay_scale, τ_forward = 0.95, 0.4, 5.0, 0.6
    t_report = Th - 8.0
    os = onset_rescale(Gamma(3.0, 2.1), r)

    ## Deaths are purely BVD-driven, so the whole mean scales by
    ## onset_fraction.
    d1 = B._forecast_deaths_mean(r, Th, α, θ, CFR; onset_fraction = os)
    d0 = B._forecast_deaths_mean(r, Th, α, θ, CFR)
    @test d1 ≈ os * d0

    ## Received backlog (receipt-delayed), positivity and the
    ## capacity-limited analysed increment are positive and finite; the
    ## lab analyses at most the available backlog.
    f_receipt = Gamma(2.0, 1.5)
    recv = B._forecast_received(r, Th, α_rep, θ_rep, p_drc, λ_bg, τ_forward,
        f_receipt; onset_fraction = os)
    p_pos = B._forecast_positivity(Th, t_report, s_test, spec_test, q0, qinf,
        decay_scale)
    Δa = B._forecast_analysed_increment(recv, 40.0, 50.0, 7.0)
    @test recv > 0 && isfinite(recv)
    @test 0 < p_pos < 1
    @test 0 < Δa <= recv

    ## Reported cases scale only the BVD part; the constant background
    ## λ_bg·Th is unscaled, so the scaled mean is strictly above os·mean.
    rc1 = B._forecast_cases_mean(r, Th, α_rep, θ_rep, p_drc, λ_bg;
        onset_fraction = os)
    rc0 = B._forecast_cases_mean(r, Th, α_rep, θ_rep, p_drc, λ_bg)
    bg_cum = λ_bg * Th
    bvd0 = rc0 - bg_cum
    @test rc1 ≈ os * bvd0 + bg_cum
    @test rc1 > os * rc0

    ## Default onset_fraction = 1 recovers the pre-infection-layer mean.
    @test B._forecast_deaths_mean(r, Th, α, θ, CFR) ≈ d0
end

@testitem "committed-deaths counterfactual applies onset_fraction" begin
    using Distributions: Gamma
    using BVDOutbreakSize: onset_rescale
    import BVDOutbreakSize as B

    r, T, α, θ, CFR = 0.05, 100.0, 4.3, 2.6, 0.3
    os = onset_rescale(Gamma(3.0, 2.1), r)
    base = B._committed_deaths_one(r, T, α, θ, CFR)
    scaled = B._committed_deaths_one(r, T, α, θ, CFR; onset_fraction = os)
    @test scaled ≈ os * base
end

@testitem "forecast_vs_truth_trajectory scores the daily trajectory" tags=[:slow] setup=[
    ForecastFixtures] begin
    using DataFrames: DataFrame, nrow
    using BVDOutbreakSize: forecast_vs_truth_trajectory

    chn=_forecast_chain(200)
    ## 16 May is the fit cut-off; 18-22 May are the post-cut-off
    ## vintages that should be scored.
    dates=["2026-05-16", "2026-05-18", "2026-05-20", "2026-05-22"]
    cases=[336, 516, 672, 872]
    deaths=[88, 131, 160, 204]

    tbl=forecast_vs_truth_trajectory(chn;
        dates = dates, cases = cases, deaths = deaths,
        snapshot_date = "2026-05-16",
        daily_travellers = 1871, source_population = 4_392_200,
        baseline_cases = 336, baseline_deaths = 88)

    @test tbl isa DataFrame
    ## Two streams x three post-snapshot vintages; the snapshot date
    ## itself (horizon 0) is excluded.
    @test nrow(tbl) == 6
    @test !("2026-05-16" in tbl[!, "Date"])
    @test Set(tbl[!, "Date"]) ==
          Set(["2026-05-18", "2026-05-20", "2026-05-22"])
    @test Set(tbl[!, "Stream"]) ==
          Set(["DRC reported cases", "DRC deaths"])
    @test sort(unique(tbl[!, "Horizon (days)"])) == [2, 4, 6]

    ## Coverage flag agrees with the reported 90% endpoints.
    for row in eachrow(tbl)
        covered=row["Lower 90%"]<=row.Observed<=row["Upper 90%"]
        @test row["Within 90% PI"] == (covered ? "yes" : "no")
    end

    ## Vintages at or before the cut-off yield no rows.
    empty=forecast_vs_truth_trajectory(chn;
        dates = ["2026-05-14", "2026-05-16"], cases = [200, 336],
        deaths = [50, 88], snapshot_date = "2026-05-16",
        daily_travellers = 1871, source_population = 4_392_200)
    @test nrow(empty) == 0

    ## Mismatched input lengths error.
    @test_throws ErrorException forecast_vs_truth_trajectory(chn;
        dates = ["2026-05-18"], cases = [1, 2], deaths = [1],
        snapshot_date = "2026-05-16",
        daily_travellers = 1871, source_population = 4_392_200)
end
