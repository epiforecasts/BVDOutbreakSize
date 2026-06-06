## Smoke tests for the one-week-ahead forecast. Builds a tiny
## synthetic chain carrying the parameters `forecast_reported` reads,
## then checks the returned DataFrame contract.

@testsnippet ForecastFixtures begin
    using Turing: Turing, @model, sample, Prior
    using Distributions: Beta, Normal, LogNormal, truncated
    import FlexiChains
    using BVDOutbreakSize: bvd_joint

    ## Synthetic prior carrying every parameter name that
    ## `forecast_reported` reads. `include_lab = true` adds the
    ## lab-turnaround delay and PCR sensitivity draws so the
    ## confirmed-cases columns are populated.
    @model function _forecast_test(; include_lab::Bool = false,
            delay::Bool = false)
        r ~ truncated(Normal(0.05, 0.01); lower = 1e-3)
        T ~ truncated(Normal(100.0, 10.0); lower = 1.0)
        ## Latent cumulative infections C(T) = 2^m, the chain symbol the
        ## forecast and counterfactual read for new/committed infections.
        m ~ truncated(Normal(13.0, 0.5); lower = 1.0)
        cumulative_infections := 2.0^m
        CFR ~ Beta(6.0, 14.0)
        α ~ truncated(Normal(4.3, 0.5); lower = 0.5)
        θ ~ truncated(Normal(2.6, 0.3); lower = 0.2)
        ## Window-mechanism chains carry `w`; delay-mechanism chains do
        ## not. Delay-mechanism export forecasts build the
        ## infection→detection delay from the incubation (α_inc / θ_inc)
        ## and DRC report (α_rep / θ_rep) draws.
        if delay
            α_inc ~ truncated(Normal(1.1, 0.3); lower = 0.3)
            θ_inc ~ truncated(Normal(5.7, 1.0); lower = 0.5)
        else
            w ~ truncated(Normal(15.0, 2.0); lower = 1.0)
        end
        p_drc ~ Beta(2.0, 6.0)
        p_uganda ~ Beta(2.0, 6.0)
        inv_sqrt_k ~ truncated(Normal(0.0, 1.0); lower = 1e-3)
        k := 1.0 / (inv_sqrt_k^2 + eps(typeof(inv_sqrt_k)))
        α_rep ~ truncated(Normal(4.0, 0.5); lower = 0.5)
        θ_rep ~ truncated(Normal(3.0, 0.3); lower = 0.2)
        λ_bg ~ truncated(Normal(0.0, 10.0); lower = 0)
        if include_lab
            ## Confirmed-stream draws the forecast reads: sensitivity,
            ## specificity and the q-curve shape (q0, qinf, decay).
            ## Forwarding is handled by the lab queue, so there is no
            ## forwarded fraction.
            s_test ~ Beta(15.0, 2.0)
            spec_test ~ Beta(50.0, 1.5)
            q0 ~ Beta(20.0, 1.5)
            qinf ~ Beta(6.0, 6.0)
            decay_scale ~ truncated(Normal(0.0, 10.0); lower = 0)
            ## Receipt-delay Gamma and the cut-off analysis capacity the
            ## capacity-limited confirmed forecast reads.
            α_recv ~ truncated(Normal(2.0, 1.0); lower = 0.1)
            θ_recv ~ truncated(Normal(1.5, 0.75); lower = 0.1)
            capacity_cutoff ~ truncated(Normal(150.0, 30.0); lower = 1.0)
            ## Confirmed-death stream draws: the case→death testing factor
            ## (scaling the shared capacity), deaths drift and suspect-death
            ## background; share the case-lab sensitivity / specificity /
            ## capacity / receipt delay and the CFR / onset-to-death delay.
            death_factor ~ LogNormal(0.0, 0.45)
            p_deaths ~ Beta(6.0, 2.0)
            λ_bg_death ~ truncated(Normal(0.0, 5.0); lower = 0)
        end
        return nothing
    end

    _forecast_chain(n;
        include_lab::Bool = false, delay::Bool = false) = sample(
        _forecast_test(; include_lab, delay), Prior(), n;
        chain_type = FlexiChains.VNChain, progress = false)
end

@testitem "forecast_reported returns the documented columns" tags=[:slow] setup=[ForecastFixtures] begin
    using DataFrames: DataFrame, nrow
    using BVDOutbreakSize: forecast_reported

    chn=_forecast_chain(200)

    fc=forecast_reported(chn;
        horizon = 7,
        daily_travellers = 1871,
        source_population = 4_392_200)

    @test fc isa DataFrame
    @test nrow(fc) == 200
    cols=[:infections_new, :bvd_deaths_new]
    @test all(c -> c in propertynames(fc), cols)
    ## Counts are non-negative, finite integers.
    @test all(fc.infections_new .>= 0)
    @test all(fc.bvd_deaths_new .>= 0)
    @test all(isfinite, fc.infections_new)
    @test all(isfinite, fc.bvd_deaths_new)
    ## Suspected/reported streams are dropped.
    @test !(:cases_new in propertynames(fc))
    @test !(:deaths_new in propertynames(fc))
    @test !(:tests_new in propertynames(fc))
    ## Exports default off.
    @test !(:exports_new in propertynames(fc))
    ## No lab-turnaround draws in this fixture → no confirmed columns.
    @test !(:confirmed_new in propertynames(fc))
end

@testitem "forecast_reported adds exports when forecast_exports=true" tags=[:slow] setup=[ForecastFixtures] begin
    using DataFrames: DataFrame
    using BVDOutbreakSize: forecast_reported, forecast_table, plot_forecast

    chn=_forecast_chain(200)
    fc=forecast_reported(chn;
        horizon = 7, daily_travellers = 1871,
        source_population = 4_392_200,
        forecast_exports = true)

    @test :exports_new in propertynames(fc)
    @test all(fc.exports_new .>= 0)
    ## Default (false) drops the export column.
    fc0=forecast_reported(chn;
        horizon = 7, daily_travellers = 1871,
        source_population = 4_392_200)
    @test !(:exports_new in propertynames(fc0))
    ## The table and plot must not assume the export column.
    tbl=forecast_table(fc0)
    @test !any(occursin.("Uganda", string.(tbl[!, 1])))
    @test plot_forecast(fc0) !== nothing
end

@testitem "forecast_reported works on the onset-to-detection delay chain" tags=[:slow] setup=[ForecastFixtures] begin
    using DataFrames: DataFrame
    using BVDOutbreakSize: forecast_reported

    ## Delay-mechanism chain carries no `w` (exports reuse the DRC report
    ## delay), so this exercises the non-window branch of
    ## `forecast_reported` that the window fixture never reaches.
    chn=_forecast_chain(50; delay = true)

    fc=forecast_reported(chn;
        horizon = 7,
        daily_travellers = 1871,
        source_population = 4_392_200,
        forecast_exports = true)

    @test fc isa DataFrame
    @test :exports_new in propertynames(fc)
    @test all(fc.exports_new .>= 0)
end

@testitem "forecast_reported adds confirmed columns when lab delay sampled" tags=[:slow] setup=[ForecastFixtures] begin
    using DataFrames: DataFrame
    using BVDOutbreakSize: forecast_reported

    chn=_forecast_chain(200; include_lab = true)
    fc=forecast_reported(chn;
        horizon = 7, daily_travellers = 1871,
        source_population = 4_392_200,
        obs_confirmed = 33, obs_analysed = 211)

    @test :confirmed_new in propertynames(fc)
    @test all(fc.confirmed_new .>= 0)
    @test all(isfinite, fc.confirmed_new)
end

@testitem "forecast_vs_truth scores confirmed cases and deaths" tags=[:slow] setup=[
    ForecastFixtures, HeadlessMakie] begin
    using DataFrames: DataFrame, nrow
    using BVDOutbreakSize: forecast_reported, forecast_vs_truth,
                           plot_forecast_vs_truth

    chn=_forecast_chain(200; include_lab = true)
    fc=forecast_reported(chn;
        horizon = 6, daily_travellers = 1871,
        source_population = 4_392_200,
        obs_confirmed = 210, obs_confirmed_deaths = 17,
        obs_analysed = 755, forecast_exports = false)
    @test :confirmed_new in propertynames(fc)
    @test :confirmed_deaths_new in propertynames(fc)

    tbl=forecast_vs_truth(fc;
        confirmed = 381, confirmed_deaths = 64,
        baseline_confirmed = 210, baseline_confirmed_deaths = 17)
    @test tbl isa DataFrame
    ## Cumulative and new view for each of confirmed cases and deaths.
    @test nrow(tbl) == 4
    @test names(tbl) ==
          ["Quantity", "Observed", "Lower 90%", "Lower 60%",
        "Lower 30%", "Upper 30%", "Upper 60%", "Upper 90%",
        "Within 90% PI"]
    @test Set(tbl[!, "Quantity"]) == Set([
        "Confirmed cases (DRC), cumulative", "Confirmed cases (DRC), new",
        "Confirmed deaths (DRC), cumulative", "Confirmed deaths (DRC), new"])

    ## A truth inside the predicted 90% interval is flagged covered; one
    ## outside is not. Read the cumulative-cases row's interval and probe
    ## both sides of it with a re-scored table.
    row=findfirst(==("Confirmed cases (DRC), cumulative"),
        tbl[!, "Quantity"])
    lo=tbl[row, "Lower 90%"]
    hi=tbl[row, "Upper 90%"]
    inside=(lo+hi)/2
    covered=forecast_vs_truth(fc;
        confirmed = inside, confirmed_deaths = 17,
        baseline_confirmed = 210, baseline_confirmed_deaths = 17)
    @test covered[row, "Within 90% PI"] == "yes"
    missed=forecast_vs_truth(fc;
        confirmed = hi+1_000, confirmed_deaths = 17,
        baseline_confirmed = 210, baseline_confirmed_deaths = 17)
    @test missed[row, "Within 90% PI"] == "no"

    fig=plot_forecast_vs_truth(fc;
        confirmed = 381, confirmed_deaths = 64,
        baseline_confirmed = 210, baseline_confirmed_deaths = 17)
    @test fig !== nothing
end

@testitem "forecast_vs_truth_trajectory scores each post-cutoff vintage" tags=[:slow] setup=[ForecastFixtures] begin
    using DataFrames: DataFrame
    using BVDOutbreakSize: forecast_vs_truth_trajectory

    chn=_forecast_chain(100; include_lab = true)
    traj=forecast_vs_truth_trajectory(chn;
        dates = ["2026-05-29", "2026-05-31", "2026-06-03"],
        confirmed = [263, 321, 381],
        confirmed_deaths = [42, 48, 64],
        snapshot_date = "2026-05-28",
        baseline_confirmed = 210, baseline_confirmed_deaths = 17,
        baseline_analysed = 755)
    @test traj isa DataFrame
    ## Two confirmed quantities at each of three post-cut-off vintages.
    @test "Horizon (days)" in names(traj)
    @test "Within 90% PI" in names(traj)
    @test Set(traj[!, "Quantity"]) ==
          Set(["Confirmed cases (DRC)", "Confirmed deaths (DRC)"])
    @test all(traj[!, "Horizon (days)"] .> 0)
end

@testitem "forecast_table and plot_forecast" tags=[:slow] setup=[
    ForecastFixtures, HeadlessMakie] begin
    using DataFrames: DataFrame, nrow
    using BVDOutbreakSize: forecast_reported, forecast_table, plot_forecast

    chn=_forecast_chain(200)
    fc=forecast_reported(chn;
        horizon = 7, daily_travellers = 1871,
        source_population = 4_392_200)

    tbl=forecast_table(fc)
    @test tbl isa DataFrame
    ## Two trusted quantities (infections, BVD deaths) without lab draws.
    @test nrow(tbl) == 2
    @test names(tbl) ==
          ["Quantity", "Lower 90%", "Lower 60%", "Lower 30%",
        "Upper 30%", "Upper 60%", "Upper 90%"]
    @test Set(tbl[!, "Quantity"]) ==
          Set(["New infections", "New BVD deaths (DRC)"])

    fig=plot_forecast(fc)
    @test fig !== nothing
end

@testitem "predict_committed returns the four committed totals" tags=[:slow] setup=[
    ForecastFixtures, HeadlessMakie] begin
    using DataFrames: DataFrame, nrow
    using BVDOutbreakSize: predict_committed, plot_committed

    chn=_forecast_chain(200; include_lab = true)
    obs_confirmed=33
    obs_confirmed_deaths=11
    df=predict_committed(chn;
        obs_confirmed = obs_confirmed,
        obs_confirmed_deaths = obs_confirmed_deaths,
        obs_analysed = 211,
        horizon_days = 365)

    @test df isa DataFrame
    @test nrow(df) == 200
    cols=[:infections, :bvd_deaths, :confirmed_cases, :confirmed_deaths]
    @test all(c -> c in propertynames(df), cols)
    @test all(isfinite, df.infections)
    @test all(isfinite, df.bvd_deaths)
    @test all(isfinite, df.confirmed_cases)
    @test all(isfinite, df.confirmed_deaths)
    ## Committed infections equal the cumulative infection draws C(T) = 2^m.
    ci=vec(Array(chn[:cumulative_infections]))
    @test df.infections ≈ float.(ci)
    ## Committed BVD deaths = CFR · C(T) and never exceed the infected pool.
    @test all(df.bvd_deaths .>= 0)
    @test all(df.bvd_deaths .<= df.infections)
    ## Committed confirmed totals can only add to what was already confirmed.
    @test all(df.confirmed_cases .>= obs_confirmed)
    @test all(df.confirmed_deaths .>= obs_confirmed_deaths)

    fig=plot_committed(df;
        obs_confirmed = obs_confirmed,
        obs_confirmed_deaths = obs_confirmed_deaths)
    @test fig !== nothing

    ## Without lab args the lab columns are dropped.
    df0=predict_committed(chn)
    @test :infections in propertynames(df0)
    @test :bvd_deaths in propertynames(df0)
    @test !(:confirmed_cases in propertynames(df0))
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
    q0, qinf, decay_scale = 0.95, 0.4, 5.0
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
    recv = B._forecast_received(r, Th, α_rep, θ_rep, p_drc, λ_bg,
        f_receipt; onset_fraction = os)
    p_pos = B._forecast_positivity(r, Th, α_rep, θ_rep, p_drc, λ_bg, s_test,
        spec_test; onset_fraction = os)
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

@testitem "committed confirmed-deaths helper adds to the observed count" begin
    using Distributions: Gamma
    using BVDOutbreakSize: onset_rescale
    import BVDOutbreakSize as B

    r, T, α, θ, CFR = 0.05, 100.0, 4.3, 2.6, 0.3
    p_deaths, λ_bg_death = 1.0, 0.5
    α_recv, θ_recv = 2.0, 1.5
    s, spec = 0.9, 0.97
    os = onset_rescale(Gamma(3.0, 2.1), r)
    obs_cd = 11.0
    ## Over a long horizon the committed total adds the drained suspect-death
    ## backlog (shared receipt delay) to the already-observed confirmed
    ## deaths; the death factor (capacity scaling) does not bind.
    committed = B._committed_confirmed_deaths_one(r, T, T + 365, α, θ, CFR,
        p_deaths, λ_bg_death, α_recv, θ_recv, s, spec, obs_cd, 0;
        onset_fraction = os)
    @test committed >= obs_cd
    @test isfinite(committed)
end
