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
    ## :expected_confirmed_T, :expected_confirmed_deaths_T, :k. The
    ## cumulative-onset and cumulative-death trajectories are absent, so the
    ## latent onset and death forecasts fall back to the scalar proxies.
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

@testitem "forecast_reported returns the documented columns" tags=[
    :slow
] setup=[ForecastFixtures] begin
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
        :infections_new, :onsets_new, :deaths_latent_new, :rt_forecast]
    @test all(c -> c in propertynames(fc), cols)
    @test all(fc.infections_new .>= 0)
    @test all(fc.onsets_new .>= 0)
    @test all(fc.deaths_latent_new .>= 0)
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

@testitem "forecast_reported projects isolation beds and recovered" tags=[
    :slow
] begin
    using Turing: @model, sample, Prior
    using Distributions: Normal, truncated
    import FlexiChains
    using DataFrames: DataFrame
    using BVDOutbreakSize: forecast_reported, forecast_table

    ## Synthetic chain that also carries the isolation and recovered
    ## deterministics, so the forecast projects the bed occupancy and the
    ## cumulative recovered using their own dispersions.
    @model function _forecast_iso_test()
        r ~ truncated(Normal(0.05, 0.01); lower = 1e-3)
        inv_sqrt_k ~ truncated(Normal(0.5, 0.2); lower = 1e-3)
        k := 1.0 / (inv_sqrt_k^2 + eps(typeof(inv_sqrt_k)))
        expected_reports_T ~ truncated(Normal(300.0, 50.0); lower = 1.0)
        expected_deaths_T ~ truncated(Normal(15.0, 3.0); lower = 1.0)
        expected_infections_T ~ truncated(Normal(800.0, 100.0); lower = 1.0)
        R_T ~ truncated(Normal(1.5, 0.3); lower = 1e-3)
        expected_confirmed_T ~ truncated(Normal(120.0, 20.0); lower = 1.0)
        expected_confirmed_deaths_T ~ truncated(Normal(8.0, 2.0); lower = 0.5)
        expected_bed_demand_T ~ truncated(Normal(600.0, 80.0); lower = 1.0)
        bed_capacity ~ truncated(Normal(430.0, 40.0); lower = 1.0)
        isolation_dispersion ~ truncated(Normal(10.0, 3.0); lower = 1.0)
        expected_recovered_T ~ truncated(Normal(32.0, 8.0); lower = 1.0)
        recovered_dispersion ~ truncated(Normal(10.0, 3.0); lower = 1.0)
        return nothing
    end
    chn = sample(_forecast_iso_test(), Prior(), 200;
        chain_type = FlexiChains.VNChain, progress = false)
    fc = forecast_reported(chn; horizon = 7,
        obs_cases = 905, obs_deaths = 18,
        obs_confirmed = 210, obs_confirmed_deaths = 17,
        obs_recovered = 32)
    @test :bed_demand in propertynames(fc)
    @test :isolation_level in propertynames(fc)
    @test :recovered_cum in propertynames(fc)
    @test :recovered_new in propertynames(fc)
    @test all(fc.isolation_level .>= 0)
    @test all(fc.bed_demand .>= 0)
    ## Supply-limited occupancy never exceeds the projected demand.
    @test all(fc.isolation_level .<= fc.bed_demand)
    ## Occupancy is min(demand, capacity): with demand (~600) above capacity
    ## (~430) it reaches the capacity, where the old exponential soft cap would
    ## have plateaued well below it (~0.75 of capacity).
    @test maximum(fc.isolation_level) > 400
    @test all(fc.recovered_cum .>= 0)
    @test all(fc.recovered_new .<= fc.recovered_cum)
    tbl = forecast_table(fc)
    @test "DRC isolation beds" in tbl[!, "Stream"]
    @test "DRC recovered" in tbl[!, "Stream"]
    @test "demand at T+7" in tbl[!, "Quantity"]
    @test "occupancy at T+7" in tbl[!, "Quantity"]
    ## The bed forecast is validated against an observed occupancy when one is
    ## supplied; without it the beds are not scored.
    using BVDOutbreakSize: forecast_vs_truth
    vt = forecast_vs_truth(fc;
        observed = (confirmed_cum = 210, confirmed_deaths_cum = 17),
        isolation = 359)
    @test "DRC isolation beds" in vt[!, "Stream"]
    vt0 = forecast_vs_truth(fc;
        observed = (confirmed_cum = 210, confirmed_deaths_cum = 17))
    @test "DRC isolation beds" ∉ vt0[!, "Stream"]
end

@testitem "forecast_table has expected rows and columns" tags=[:slow] setup=[
    ForecastFixtures
] begin
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
    ## Two confirmed streams x two quantities (cumulative, new this week).
    @test nrow(tbl) == 4
    @test names(tbl) ==
          ["Stream", "Quantity", "Lower 90%", "Lower 60%", "Lower 30%",
        "Upper 30%", "Upper 60%", "Upper 90%"]
    @test Set(tbl[!, "Stream"]) ==
          Set(["DRC confirmed cases", "DRC confirmed deaths"])
    @test Set(tbl[!, "Quantity"]) ==
          Set(["cumulative by T+7", "new this week"])
end

@testitem "forecast_vs_truth compares forecast to observed counts" tags=[
    :slow
] setup=[ForecastFixtures] begin
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
        observed = (confirmed_cum = 260, confirmed_deaths_cum = 20))

    @test tbl isa DataFrame
    ## Two confirmed streams x two quantities (cumulative, new this week).
    @test nrow(tbl) == 4
    @test names(tbl) ==
          ["Stream", "Quantity", "Observed", "Lower 90%", "Lower 60%",
        "Lower 30%", "Upper 30%", "Upper 60%", "Upper 90%", "Within 90% PI"]
    @test Set(tbl[!, "Stream"]) ==
          Set(["DRC confirmed cases", "DRC confirmed deaths"])
    @test Set(tbl[!, "Quantity"]) ==
          Set(["cumulative by T+7", "new this week"])

    for row in eachrow(tbl)
        covered=row["Lower 90%"]<=row.Observed<=row["Upper 90%"]
        @test row["Within 90% PI"] == (covered ? "yes" : "no")
    end
end

@testitem "forecast_vs_truth covers all streams and guards on observed" begin
    using Random: MersenneTwister
    using DataFrames: DataFrame, nrow
    using BVDOutbreakSize: forecast_vs_truth
    rng = MersenneTwister(51)
    n = 300
    fc = DataFrame(
        cases_cum = rand(rng, 50:150, n), cases_new = rand(rng, 0:30, n),
        deaths_cum = rand(rng, 40:100, n), deaths_new = rand(rng, 0:20, n),
        confirmed_cum = rand(rng, 20:80, n),
        confirmed_new = rand(rng, 0:15, n),
        confirmed_deaths_cum = rand(rng, 1:20, n),
        confirmed_deaths_new = rand(rng, 0:5, n),
        recovered_cum = rand(rng, 10:60, n),
        recovered_new = rand(rng, 0:10, n),
        isolation_level = rand(rng, 250:400, n)
    )
    ## Every count stream supplied an observed cumulative gives two rows each
    ## (cumulative + new); the beds add one level row when an occupancy is
    ## supplied. Ten count rows + one beds row.
    tbl = forecast_vs_truth(fc;
        observed = (cases_cum = 140, deaths_cum = 90, confirmed_cum = 70,
            confirmed_deaths_cum = 18, recovered_cum = 55),
        baseline = (confirmed_cum = 40,),
        isolation = 359)
    @test nrow(tbl) == 11
    @test "Quantity" in names(tbl)
    @test Set(tbl[!, "Stream"]) == Set([
        "DRC reported cases", "DRC suspected deaths", "DRC confirmed cases",
        "DRC confirmed deaths", "DRC recovered", "DRC isolation beds"])
    ## A stream present in the frame but absent from `observed` is skipped, and
    ## without an observed occupancy the beds row is dropped.
    tbl2 = forecast_vs_truth(fc;
        observed = (confirmed_cum = 70, confirmed_deaths_cum = 18))
    @test nrow(tbl2) == 4
    @test "DRC reported cases" ∉ tbl2[!, "Stream"]
    @test "DRC isolation beds" ∉ tbl2[!, "Stream"]
    ## baseline shifts the new-count observed: the confirmed "new this week"
    ## row is scored against observed - baseline = 70 - 40 = 30.
    conf_new = tbl[(tbl.Stream .== "DRC confirmed cases") .&
                   (tbl.Quantity .== "new this week"), :]
    @test only(conf_new.Observed) == 30
    ## Absent baseline defaults to zero, so the new-count observed is the full
    ## observed cumulative.
    tbl3 = forecast_vs_truth(fc; observed = (confirmed_cum = 70,))
    c3 = tbl3[(tbl3.Stream .== "DRC confirmed cases") .&
              (tbl3.Quantity .== "new this week"), :]
    @test only(c3.Observed) == 70
end

@testitem "forecast_archive returns tidy long scored streams" tags=[
    :slow
] setup=[ForecastFixtures] begin
    using DataFrames: DataFrame, nrow
    using Dates: Date, Day
    using BVDOutbreakSize: forecast_reported, forecast_archive

    chn=_forecast_chain(200)
    made=Date("2026-06-07")
    fcs=[(h,
             forecast_reported(chn;
                 horizon = h,
                 obs_cases = 905, obs_deaths = 18,
                 obs_confirmed = 210, obs_confirmed_deaths = 17))
         for h in (7, 14)]
    arch=forecast_archive(fcs; made_date = made, thin = 2)

    @test arch isa DataFrame
    @test names(arch) ==
          ["made_date", "horizon", "target_date", "stream", "draw", "value"]
    ## Only the incident confirmed streams are carried by this chain (recovered
    ## and beds are absent, so skipped), across the two horizons.
    @test Set(arch.stream) == Set(["confirmed cases", "confirmed deaths"])
    @test Set(arch.horizon) == Set([7, 14])
    @test all(arch.made_date .== made)
    ## target_date is made_date plus the horizon.
    @test all(arch.target_date .== arch.made_date .+ Day.(arch.horizon))
    ## Thinning keeps every second draw: 200 / 2 = 100 per (stream, horizon).
    sub=arch[(arch.stream .== "confirmed cases") .& (arch.horizon .== 7), :]
    @test nrow(sub) == 100
    @test all(arch.value .>= 0)
end

@testitem "forecast_archive carries the ward-bed occupancy split" begin
    using DataFrames: DataFrame
    using Dates: Date
    using BVDOutbreakSize: forecast_archive

    ## A forecast frame carrying the confirmed/suspect ward split (as
    ## `forecast_reported` emits once the chain carries the confirmed in-care
    ## sub-stock) archives the two ward occupancy levels alongside the total,
    ## each under its own scored label. A frame without the split carries the
    ## total alone, so the hook is dormant until the columns appear.
    n = 40
    occ = collect(1:n)
    conf_occ = fld.(occ, 2)
    fc = DataFrame(isolation_level = occ, confirmed_occupancy = conf_occ,
        suspect_occupancy = occ .- conf_occ)
    arch = forecast_archive([(7, fc)]; made_date = Date("2026-06-13"))
    @test "isolation beds" in arch.stream
    @test "treatment beds" in arch.stream
    @test "isolation beds (suspected)" in arch.stream
    ## The ward occupancy values round-trip and partition the total.
    _vals(s) = sort(arch[arch.stream .== s, :], :draw).value
    tot = _vals("isolation beds")
    tre = _vals("treatment beds")
    iso = _vals("isolation beds (suspected)")
    @test tre .+ iso == tot

    ## A frame without the split archives only the total.
    arch0 = forecast_archive([(7, DataFrame(isolation_level = occ))];
        made_date = Date("2026-06-13"))
    @test Set(arch0.stream) == Set(["isolation beds"])
end

@testitem "forecast cumulative streams never fall below the cut-off" tags=[
    :slow
] begin
    using Turing: @model, sample, Prior
    using Distributions: Normal, truncated
    import FlexiChains
    using Statistics: median
    using BVDOutbreakSize: forecast_reported

    ## A declining chain (growth rate r < 0, R_T < 1): the regime where a
    ## stock-scaling projection (`cumulative_T * exp(r * horizon)`) can
    ## shrink the cumulative below the observed cut-off, an impossible
    ## decreasing cumulative. The outbreak age `:T` is carried so the daily
    ## incidence at the cut-off is inferred from the cumulative total.
    @model function _forecast_decline_test()
        r ~ truncated(Normal(-0.05, 0.01); upper = -1e-3)
        inv_sqrt_k ~ truncated(Normal(0.5, 0.2); lower = 1e-3)
        k := 1.0 / (inv_sqrt_k^2 + eps(typeof(inv_sqrt_k)))
        T := 100.0
        expected_reports_T ~ truncated(Normal(900.0, 80.0); lower = 1.0)
        expected_deaths_T ~ truncated(Normal(40.0, 6.0); lower = 1.0)
        expected_infections_T ~ truncated(Normal(2000.0, 200.0); lower = 1.0)
        R_T ~ truncated(Normal(0.7, 0.1); lower = 1e-3, upper = 1.0)
        expected_confirmed_T ~ truncated(Normal(800.0, 60.0); lower = 1.0)
        expected_confirmed_deaths_T ~ truncated(Normal(30.0, 5.0); lower = 0.5)
        return nothing
    end
    chn = sample(_forecast_decline_test(), Prior(), 400;
        chain_type = FlexiChains.VNChain, progress = false)

    obs_cases, obs_deaths, obs_confirmed = 1205, 60, 1203
    fc(h) = forecast_reported(chn; horizon = h,
        obs_cases = obs_cases, obs_deaths = obs_deaths,
        obs_confirmed = obs_confirmed, obs_confirmed_deaths = 35)
    f7, f21 = fc(7), fc(21)

    ## Every draw's projected cumulative stays at or above the observed
    ## cut-off, even though the outbreak is shrinking.
    @test all(f7.cases_cum .>= obs_cases)
    @test all(f7.deaths_cum .>= obs_deaths)
    @test all(f7.confirmed_cum .>= obs_confirmed)
    @test all(f21.confirmed_cum .>= obs_confirmed)

    ## The cumulative grows with the horizon (more new counts accrue) rather
    ## than shrinking, which a stock-scaling projection would do here.
    @test median(f21.cases_cum) >= median(f7.cases_cum)
    @test median(f21.deaths_cum) >= median(f7.deaths_cum)
    @test median(f21.confirmed_cum) >= median(f7.confirmed_cum)
end

## Per-stream forecasts from either fit kind. Two synthetic chains stand in
## for the two shapes `forecast_stream` must serve: a top-level chain naming
## the joint's un-prefixed aliases, and a nested chain naming the
## single-stream composers' submodel-bound deterministics. The nested chain
## carries no `r`, `R_T` or `T`, so it exercises the reconstruction path.
@testsnippet StreamFixtures begin
    using Turing: @model, sample, Prior
    using Distributions: Normal, truncated, product_distribution
    import FlexiChains
    using BVDOutbreakSize: knot_days

    ## Grid length and intervention day shared by the nested chain and the
    ## `forecast_stream` calls that reconstruct its walk.
    const STREAM_N = 60
    const STREAM_BREAK = 30.0
    ## One innovation per walk step, matching `rt_walk_model` at these
    ## settings, so `reconstruct_rt` accepts the chain.
    const STREAM_NZ = length(knot_days(STREAM_N; week = 7, start = 1)) - 1

    ## Joint-shaped: every quantity un-prefixed, as `bvd_joint`'s `:=`
    ## aliases expose them.
    @model function _stream_toplevel()
        r ~ truncated(Normal(0.05, 0.01); lower = 1e-3)
        R_T ~ truncated(Normal(1.5, 0.3); lower = 1e-3)
        T ~ truncated(Normal(100.0, 5.0); lower = 10.0)
        k_cases ~ truncated(Normal(10.0, 3.0); lower = 1.0)
        k_deaths ~ truncated(Normal(10.0, 3.0); lower = 1.0)
        k_confirmed ~ truncated(Normal(10.0, 3.0); lower = 1.0)
        k_confirmed_deaths ~ truncated(Normal(10.0, 3.0); lower = 1.0)
        isolation_dispersion ~ truncated(Normal(10.0, 3.0); lower = 1.0)
        expected_reports_T ~ truncated(Normal(900.0, 50.0); lower = 1.0)
        expected_deaths_T ~ truncated(Normal(40.0, 5.0); lower = 1.0)
        expected_confirmed_T ~ truncated(Normal(210.0, 20.0); lower = 1.0)
        expected_confirmed_deaths_T ~ truncated(Normal(17.0, 3.0); lower = 0.5)
        expected_exports_T ~ truncated(Normal(12.0, 3.0); lower = 0.5)
        expected_bed_demand_T ~ truncated(Normal(600.0, 80.0); lower = 1.0)
        bed_capacity ~ truncated(Normal(430.0, 40.0); lower = 1.0)
        return nothing
    end

    ## Single-stream-shaped: each stream's expected count under its composer
    ## binding, one scalar `dispersion_state.k`, and the walk and generation
    ## interval the reconstruction reads. `growth_state.T` is the cryptic
    ## duration only, and there is no `r`, `R_T` or `T`.
    ##
    ## Confirmed cases are the exception to the nested naming: their
    ## submodel keeps `expected_confirmed` on a plain `=` (required for
    ## Enzyme's AD), so `confirmed_only_model` aliases the cut-off count
    ## un-prefixed as `expected_confirmed_T`, exactly as the joint does.
    @model function _stream_nested()
        var"growth_state.T" ~ truncated(Normal(40.0, 5.0); lower = 1.0)
        var"growth_state.r" ~ truncated(Normal(0.2, 0.02); lower = 1e-3)
        var"rt_state.log_R0" ~ Normal(log(1.8), 0.1)
        var"rt_state.sigma_rw" ~ truncated(Normal(0.1, 0.02); lower = 1e-3)
        var"rt_state.intervention_effect" ~ Normal(-0.2, 0.05)
        var"rt_state.z" ~ product_distribution(fill(Normal(0, 1), STREAM_NZ))
        var"gi_state.α" ~ truncated(Normal(2.71, 0.1); lower = 0.1)
        var"gi_state.θ" ~ truncated(Normal(5.65, 0.2); lower = 0.1)
        var"dispersion_state.k" ~ truncated(Normal(10.0, 3.0); lower = 1.0)
        var"cases_state.expected_reports" ~
        truncated(Normal(900.0, 50.0); lower = 1.0)
        var"deaths_state.expected_deaths_T" ~
        truncated(Normal(40.0, 5.0); lower = 1.0)
        expected_confirmed_T ~ truncated(Normal(210.0, 20.0); lower = 1.0)
        var"confirmed_deaths_state.expected_confirmed_deaths" ~
        truncated(Normal(17.0, 3.0); lower = 0.5)
        var"exports_state.expected_exports_T" ~
        truncated(Normal(12.0, 3.0); lower = 0.5)
        var"treatment_state.expected_bed_demand" ~
        truncated(Normal(600.0, 80.0); lower = 1.0)
        var"treatment_state.expected_isolation" ~
        truncated(Normal(400.0, 20.0); lower = 1.0)
        ## `bed_utilisation := occ_T / C_T`, so a capacity near 430 given the
        ## occupancy above.
        var"treatment_state.bed_utilisation" ~
        truncated(Normal(0.93, 0.02); lower = 0.1, upper = 1.0)
        var"treatment_state.disp_state.k" ~
        truncated(Normal(10.0, 3.0); lower = 1.0)
        return nothing
    end

    _toplevel_chain(n) = sample(_stream_toplevel(), Prior(), n;
        chain_type = FlexiChains.VNChain, progress = false)
    _nested_chain(n) = sample(_stream_nested(), Prior(), n;
        chain_type = FlexiChains.VNChain, progress = false)

    const STREAM_OBS = Dict(
        :reported_cases => 905, :suspected_deaths => 40,
        :confirmed_cases => 210, :confirmed_deaths => 17,
        :exports => 12, :isolation_beds => 359)
    const STREAM_ALL = collect(keys(STREAM_OBS))
end

@testitem "forecast_stream covers every stream from a joint-shaped chain" tags=[
    :slow
] setup=[StreamFixtures] begin
    using BVDOutbreakSize: forecast_stream

    chn=_toplevel_chain(200)
    for s in STREAM_ALL
        fc=forecast_stream(chn, s; horizon = 7, obs_value = STREAM_OBS[s])
        @test fc isa Vector
        @test length(fc) == 200
        @test all(fc .>= 0)
    end
    ## Beds are a supply-limited level, so the occupancy cannot exceed the
    ## capacity (~430) however far the demand (~600) is projected.
    beds=forecast_stream(chn, :isolation_beds; horizon = 7,
        obs_value = 359)
    @test maximum(beds) <= 600
    @test maximum(beds) > 300
end

@testitem "forecast_stream projects every stream from a nested chain" tags=[
    :slow
] setup=[StreamFixtures] begin
    using BVDOutbreakSize: forecast_stream

    chn=_nested_chain(200)
    for s in STREAM_ALL
        fc=forecast_stream(chn, s; horizon = 7, obs_value = STREAM_OBS[s],
            n = STREAM_N, breakpoint = STREAM_BREAK)
        @test fc isa Vector
        @test length(fc) == 200
        @test all(fc .>= 0)
    end
    beds=forecast_stream(chn, :isolation_beds; horizon = 7,
        obs_value = 359, n = STREAM_N, breakpoint = STREAM_BREAK)
    @test maximum(beds) > 300
end

@testitem "forecast_stream incident streams grow with the horizon" tags=[
    :slow
] setup=[StreamFixtures] begin
    using Statistics: median
    using BVDOutbreakSize: forecast_stream

    ## More new counts accrue over a longer horizon, on both chain shapes.
    chn=_toplevel_chain(400)
    for s in (:reported_cases, :confirmed_cases, :exports)
        f7=forecast_stream(chn, s; horizon = 7, obs_value = STREAM_OBS[s])
        f21=forecast_stream(chn, s; horizon = 21, obs_value = STREAM_OBS[s])
        @test median(f21) >= median(f7)
    end
    nchn=_nested_chain(400)
    f7=forecast_stream(nchn, :confirmed_cases; horizon = 7,
        obs_value = 210, n = STREAM_N, breakpoint = STREAM_BREAK)
    f21=forecast_stream(nchn, :confirmed_cases; horizon = 21,
        obs_value = 210, n = STREAM_N, breakpoint = STREAM_BREAK)
    @test median(f21) >= median(f7)
end

@testitem "forecast_stream rejects unknown and unfitted streams" tags=[
    :slow
] setup=[StreamFixtures] begin
    using BVDOutbreakSize: forecast_stream

    chn=_toplevel_chain(50)
    @test_throws ArgumentError forecast_stream(chn, :not_a_stream;
        horizon = 7, obs_value = 1)
    ## A nested chain carries no `r` or `R_T`, so without the grid length and
    ## breakpoint needed to rebuild them the projection is an error rather
    ## than a silent fallback to the wrong growth rate.
    nchn=_nested_chain(50)
    @test_throws ArgumentError forecast_stream(nchn, :reported_cases;
        horizon = 7, obs_value = 905)
end

@testitem "every composer carries its own stream's forecast keys" tags=[
    :slow
] begin
    using Turing: sample, Prior
    import FlexiChains
    using BVDOutbreakSize: cases_only_model, deaths_only_model,
                           confirmed_only_model, confirmed_deaths_only_model,
                           treatment_only_model, exports_only_model,
                           _STREAM_SPEC, _resolve_draws, _bed_capacity

    ## The fixtures above pin `forecast_stream`'s key resolution against
    ## hand-written chains, which cannot catch `_STREAM_SPEC` naming a key
    ## the real model never exposes. Sample each single-stream composer
    ## from the prior and check its own stream resolves off a real chain.
    m_cases = cases_only_model(40, missing)
    m_deaths = deaths_only_model(33, missing;
        deaths_history = (; days = [13, 18, 23], counts = [131, 204, 246]))
    m_confirmed = confirmed_only_model(40, 8;
        confirmed_history = (; days = [20, 40], counts = [3, 8]),
        lab_history = (; days = [20, 40], counts = [5, 9]))
    m_confirmed_deaths = confirmed_deaths_only_model(40, 17, 246;
        deaths_history = (; days = [20, 40], counts = [120, 246]))
    fits = [(:reported_cases, m_cases), (:suspected_deaths, m_deaths),
        (:confirmed_cases, m_confirmed),
        (:confirmed_deaths, m_confirmed_deaths),
        (:isolation_beds, treatment_only_model(33)),
        (:exports, exports_only_model(40, 2))]
    for (stream, fit) in fits
        chn = sample(fit, Prior(), 10; chain_type = FlexiChains.VNChain,
            progress = false)
        spec = _STREAM_SPEC[stream]
        @test !isnothing(_resolve_draws(chn, spec.expected))
        if spec.noise === :nb
            @test !isnothing(_resolve_draws(chn, spec.dispersion))
        end
        ## The treatment composer exposes no `bed_capacity` key, so the
        ## supply limit must still be recoverable from the ratio.
        if stream === :isolation_beds
            @test !isnothing(_bed_capacity(chn))
        end
    end
end

@testitem "forecast_stream recovers bed capacity on a standalone fit" tags=[
    :slow
] begin
    using Turing: @model, sample, Prior
    using Distributions: Normal, truncated, product_distribution
    using Statistics: mean
    import FlexiChains
    using BVDOutbreakSize: forecast_stream, _bed_capacity, knot_days

    ## A standalone treatment fit exposes no capacity key, so the cut-off
    ## capacity is recovered from `expected_isolation / bed_utilisation`.
    ## Fix a known capacity and check the recovery returns it: the occupancy
    ## cancels, so the answer must not depend on the occupancy level.
    nz = length(knot_days(60; week = 7, start = 1)) - 1
    @model function _cap_test(occ, capacity)
        var"growth_state.T" ~ truncated(Normal(40.0, 5.0); lower = 1.0)
        var"rt_state.log_R0" ~ Normal(log(1.8), 0.1)
        var"rt_state.sigma_rw" ~ truncated(Normal(0.1, 0.02); lower = 1e-3)
        var"rt_state.intervention_effect" ~ Normal(-0.2, 0.05)
        var"rt_state.z" ~ product_distribution(fill(Normal(0, 1), nz))
        var"gi_state.α" ~ truncated(Normal(2.71, 0.1); lower = 0.1)
        var"gi_state.θ" ~ truncated(Normal(5.65, 0.2); lower = 0.1)
        var"treatment_state.disp_state.k" ~
        truncated(Normal(10.0, 3.0); lower = 1.0)
        var"treatment_state.expected_bed_demand" ~
        truncated(Normal(600.0, 80.0); lower = 1.0)
        var"treatment_state.expected_isolation" := occ
        var"treatment_state.bed_utilisation" := occ / capacity
        return nothing
    end

    ## A mid-outbreak occupancy and a near-zero early one both recover the
    ## same fixed capacity.
    for occ in (400.0, 1e-8)
        chn = sample(_cap_test(occ, 430.0), Prior(), 50;
            chain_type = FlexiChains.VNChain, progress = false)
        cap = _bed_capacity(chn)
        @test all(c -> isapprox(c, 430.0; rtol = 1e-6), cap)
        ## The supply limit binds: demand (~600) is capped at the capacity.
        beds = forecast_stream(chn, :isolation_beds; horizon = 7,
            obs_value = 359, n = 60, breakpoint = 30.0)
        @test maximum(beds) <= 430
    end
end

@testitem "forecast_stream beds reproduce forecast_reported's occupancy" tags=[
    :slow
] begin
    using Turing: @model, sample, Prior
    using Distributions: Normal, truncated
    using Statistics: median, mean
    import FlexiChains
    using BVDOutbreakSize: forecast_reported, forecast_stream

    ## The joint's isolation forecast must come out the same whether it is
    ## taken from `forecast_reported` (which projects every stream at once)
    ## or from `forecast_stream` (one stream), since both read the same
    ## demand, capacity and isolation dispersion. The replicates draw from
    ## different points of the RNG stream, so the distributions are compared
    ## rather than the draws.
    @model function _iso_parity()
        r ~ truncated(Normal(0.05, 0.01); lower = 1e-3)
        k ~ truncated(Normal(10.0, 3.0); lower = 1.0)
        expected_reports_T ~ truncated(Normal(900.0, 50.0); lower = 1.0)
        expected_deaths_T ~ truncated(Normal(40.0, 5.0); lower = 1.0)
        expected_infections_T ~ truncated(Normal(800.0, 100.0); lower = 1.0)
        R_T ~ truncated(Normal(1.5, 0.3); lower = 1e-3)
        expected_confirmed_T ~ truncated(Normal(210.0, 20.0); lower = 1.0)
        expected_confirmed_deaths_T ~ truncated(Normal(17.0, 3.0); lower = 0.5)
        expected_bed_demand_T ~ truncated(Normal(600.0, 80.0); lower = 1.0)
        bed_capacity ~ truncated(Normal(430.0, 40.0); lower = 1.0)
        isolation_dispersion ~ truncated(Normal(10.0, 3.0); lower = 1.0)
        return nothing
    end
    chn = sample(_iso_parity(), Prior(), 2000;
        chain_type = FlexiChains.VNChain, progress = false)

    for h in (7, 21)
        fr = forecast_reported(chn; horizon = h, obs_cases = 905,
            obs_deaths = 40, obs_confirmed = 210,
            obs_confirmed_deaths = 17).isolation_level
        fs = forecast_stream(chn, :isolation_beds; horizon = h,
            obs_value = 359)
        @test isapprox(median(fr), median(fs); rtol = 0.05)
        @test isapprox(mean(fr), mean(fs); rtol = 0.05)
    end
end

@testitem "forecast_stream reconstructs cut-off R_T at the model's ramp" tags=[
    :slow
] setup=[StreamFixtures] begin
    using BVDOutbreakSize: _cutoff_rt, reconstruct_rt

    ## A single-stream chain carries no top-level `R_T`, so `_cutoff_rt`
    ## rebuilds it from the walk. The reconstruction is ramp-sensitive, and
    ## the model's ramp is 21 (`sigmoid_ramp` / `rt_walk_model`), not
    ## `reconstruct_rt`'s lighter 14 default. Pin that `_cutoff_rt` uses 21:
    ## it must equal the ramp = 21 cut-off column and differ from ramp = 14,
    ## so a regression back to the default would fail here.
    chn=_nested_chain(200)
    got=_cutoff_rt(chn; n = STREAM_N, breakpoint = STREAM_BREAK,
        rt_start = 1, rt_walk_start = 1)
    rt21=reconstruct_rt(chn; n = STREAM_N, breakpoint = STREAM_BREAK,
        ramp = 21.0)
    rt14=reconstruct_rt(chn; n = STREAM_N, breakpoint = STREAM_BREAK,
        ramp = 14.0)
    cut21=[rt21[i, STREAM_N] for i in axes(rt21, 1)]
    cut14=[rt14[i, STREAM_N] for i in axes(rt14, 1)]

    ## Uses the model ramp exactly.
    @test got == cut21
    ## The ramp genuinely moves the cut-off R_T, so the wrong default is a
    ## real, detectable error rather than a harmless relabelling.
    @test maximum(abs.(cut21 .- cut14)) > 1e-6
    @test got != cut14
end
