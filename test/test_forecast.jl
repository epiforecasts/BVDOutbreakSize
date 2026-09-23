## The forecast readers. Every forecast is the fitted model run past its
## cut-off (`forecast_draws`), so these items draw from the fixture joint of
## test_forecast_horizon.jl and check that the national, per-stream, onset
## and province readers return the documented shapes off one set of draws.

@testsnippet ForecastDraws begin
    using BVDOutbreakSize: forecast_draws
    using Turing: sample, Prior
    import FlexiChains
    using Random: Xoshiro

    const ND = 12
    const HMAX = 14
    const OBS = (;
        obs_cases = 905, obs_deaths = 18, obs_confirmed = 40,
        obs_confirmed_deaths = 5, obs_recovered = 6,
    )
    fitted = patch_joint()
    chn = sample(
        Xoshiro(21), fitted, Prior(), ND;
        chain_type = FlexiChains.VNChain, progress = false
    )
    pp = forecast_draws(fitted, chn; horizon = HMAX)
    draws(key) = [collect(v) for v in vec(collect(pp[Symbol(key)]))]
end

@testitem "forecast_reported reads the documented columns off the draws" setup = [
    HorizonFixtures, ForecastDraws,
] begin
    using DataFrames: DataFrame, nrow
    using BVDOutbreakSize: forecast_reported

    fc = forecast_reported(pp; horizon = 7, OBS...)
    @test fc isa DataFrame
    @test nrow(fc) == ND
    cols = [
        :cases_cum, :deaths_cum, :confirmed_cum, :confirmed_deaths_cum,
        :cases_new, :deaths_new, :confirmed_new, :confirmed_deaths_new,
        :infections_new, :onsets_new, :deaths_latent_new, :rt_forecast,
        :bed_demand, :isolation_level, :bed_shortfall, :admissions_fc,
        :incare_deaths_fc, :ruleouts_fc, :recovered_cum, :recovered_new,
        :onsets_to_date, :onset_reports_new,
    ]
    @test all(c -> c in propertynames(fc), cols)
    ## Each new count is the model's own future daily counts summed over the
    ## week, and each cumulative the observed cut-off plus it.
    @test fc.confirmed_new == sum.(x -> x[1:7], draws("forecast_confirmed.increments"))
    @test fc.cases_new == sum.(x -> x[1:7], draws("forecast_reports.increments"))
    @test fc.confirmed_cum == OBS.obs_confirmed .+ fc.confirmed_new
    @test fc.recovered_cum == OBS.obs_recovered .+ fc.recovered_new
    @test fc.isolation_level ==
        round.(Int, getindex.(draws("forecast_isolation.obs"), 7))
    @test fc.rt_forecast == getindex.(draws("forecast_rt"), 7)
    for c in (:cases_new, :deaths_new, :confirmed_new, :confirmed_deaths_new)
        @test all(>=(0), fc[!, c])
    end
    ## A longer horizon adds the second week's counts to the same draws.
    fc14 = forecast_reported(pp; horizon = 14, OBS...)
    @test all(fc14.confirmed_new .>= fc.confirmed_new)
    @test_throws ArgumentError forecast_reported(pp; horizon = 15, OBS...)
    ## Without a confirmed origin the confirmed columns are left out.
    bare = forecast_reported(pp; obs_cases = 905, obs_deaths = 18)
    @test !(:confirmed_cum in propertynames(bare))
end

@testitem "forecast_draws is reproducible and keeps the fitted draws" setup = [
    HorizonFixtures, ForecastDraws,
] begin
    again = forecast_draws(fitted, chn; horizon = HMAX)
    @test draws("forecast_confirmed.increments") ==
        [collect(v) for v in vec(collect(again[Symbol("forecast_confirmed.increments")]))]
    @test vec(Array(pp[:C_T])) == vec(Array(chn[:C_T]))
end

@testitem "forecast_stream reads each stream off the joint's draws" setup = [
    HorizonFixtures, ForecastDraws,
] begin
    using BVDOutbreakSize: forecast_stream, forecast_reported, forecast_onsets

    fc = forecast_reported(pp; horizon = 7, OBS...)
    @test forecast_stream(pp, :confirmed_cases) == fc.confirmed_new
    @test forecast_stream(pp, :reported_cases) == fc.cases_new
    @test forecast_stream(pp, :suspected_deaths) == fc.deaths_new
    @test forecast_stream(pp, :confirmed_deaths) == fc.confirmed_deaths_new
    @test forecast_stream(pp, :recovered) == fc.recovered_new
    @test forecast_stream(pp, :isolation_beds) == fc.isolation_level
    @test forecast_stream(pp, :onset_reports) == fc.onset_reports_new
    ex = forecast_stream(pp, :exports; horizon = 14)
    @test length(ex) == ND && all(>=(0), ex)
    @test_throws ArgumentError forecast_stream(pp, :not_a_stream)
end

@testitem "forecast_stream reads a single-stream fit's own draws" setup = [
    HorizonFixtures,
] begin
    using BVDOutbreakSize: forecast_draws, forecast_stream
    using Turing: sample, Prior
    import FlexiChains
    using Random: Xoshiro

    for (name, stream) in (
            "cases" => :reported_cases, "deaths" => :suspected_deaths,
            "confirmed" => :confirmed_cases,
            "confirmed deaths" => :confirmed_deaths,
            "treatment" => :isolation_beds, "onsets" => :onset_reports,
            "exports" => :exports,
        )
        m = Dict(composers())[name]
        c = sample(
            Xoshiro(3), m, Prior(), 5;
            chain_type = FlexiChains.VNChain, progress = false
        )
        v = forecast_stream(forecast_draws(m, c; horizon = 7), stream)
        @test length(v) == 5
        @test all(>=(0), v)
        ## A single-stream fit forecasts only its own stream.
        stream === :reported_cases && @test_throws ArgumentError forecast_stream(
            forecast_draws(m, c; horizon = 7), :confirmed_cases
        )
    end
end

@testitem "forecast_onsets splits the reported total at the cut-off" setup = [
    HorizonFixtures, ForecastDraws,
] begin
    using BVDOutbreakSize: forecast_onsets

    fc = forecast_onsets(pp; horizon = 14, obs_value = 100)
    @test fc.onsets_unreported ≈ fc.onsets_to_date .- fc.onset_reports_to_date
    @test all(>=(0), fc.onset_reports_new)
    @test fc.onset_reports_cum == 100 .+ fc.onset_reports_new
    ## The backfill and the future reports are the two parts of the mean
    ## increment, the reports of onsets before and after the cut-off.
    @test all(>=(-1.0e-8), fc.onset_reports_backfill)
    @test all(>=(-1.0e-8), fc.onset_reports_future)
    ## The triangle is forecast a whole number of weeks ahead only.
    @test_throws ArgumentError forecast_onsets(pp; horizon = 10)
end

@testitem "forecast_table has expected rows and columns" setup = [
    HorizonFixtures, ForecastDraws,
] begin
    using DataFrames: nrow
    using BVDOutbreakSize: forecast_reported, forecast_table

    tbl = forecast_table(forecast_reported(pp; horizon = 7, OBS...))
    @test "Stream" in names(tbl)
    @test Set(tbl[!, "Quantity"]) ⊇ Set(["cumulative by T+7", "new this week"])
end

@testitem "forecast_archive returns tidy long scored streams" setup = [
    HorizonFixtures, ForecastDraws,
] begin
    using DataFrames: nrow
    using Dates: Date, Day
    using BVDOutbreakSize: forecast_reported, forecast_archive

    made = Date("2026-06-07")
    fcs = [(h, forecast_reported(pp; horizon = h, OBS...)) for h in (7, 14)]
    arch = forecast_archive(fcs; made_date = made, thin = 2)
    @test names(arch) ==
        ["made_date", "horizon", "target_date", "stream", "draw", "value"]
    @test Set(arch.stream) == Set(
        [
            "confirmed cases", "confirmed deaths", "recovered",
            "isolation beds", "onset reports",
        ]
    )
    @test all(arch.target_date .== arch.made_date .+ Day.(arch.horizon))
    sub = arch[(arch.stream .== "confirmed cases") .& (arch.horizon .== 7), :]
    @test nrow(sub) == cld(ND, 2)
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
    tbl = forecast_vs_truth(
        fc;
        observed = (
            cases_cum = 140, deaths_cum = 90, confirmed_cum = 70,
            confirmed_deaths_cum = 18, recovered_cum = 55,
        ),
        baseline = (confirmed_cum = 40,),
        isolation = 359
    )
    @test nrow(tbl) == 11
    @test "Quantity" in names(tbl)
    @test Set(tbl[!, "Stream"]) == Set(
        [
            "DRC reported cases", "DRC suspected deaths", "DRC confirmed cases",
            "DRC confirmed deaths", "DRC recovered among confirmed",
            "DRC isolation beds",
        ]
    )
    ## A stream present in the frame but absent from `observed` is skipped, and
    ## without an observed occupancy the beds row is dropped.
    tbl2 = forecast_vs_truth(
        fc;
        observed = (confirmed_cum = 70, confirmed_deaths_cum = 18)
    )
    @test nrow(tbl2) == 4
    @test "DRC reported cases" ∉ tbl2[!, "Stream"]
    @test "DRC isolation beds" ∉ tbl2[!, "Stream"]
    ## baseline shifts the new-count observed: the confirmed "new this week"
    ## row is scored against observed - baseline = 70 - 40 = 30.
    conf_new = tbl[(tbl.Stream .== "DRC confirmed cases") .& (tbl.Quantity .== "new this week"), :]
    @test only(conf_new.Observed) == 30
    ## Absent baseline defaults to zero, so the new-count observed is the full
    ## observed cumulative.
    tbl3 = forecast_vs_truth(fc; observed = (confirmed_cum = 70,))
    c3 = tbl3[(tbl3.Stream .== "DRC confirmed cases") .& (tbl3.Quantity .== "new this week"), :]
    @test only(c3.Observed) == 70
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
    fc = DataFrame(
        isolation_level = occ, confirmed_occupancy = conf_occ,
        suspect_occupancy = occ .- conf_occ
    )
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
    arch0 = forecast_archive(
        [(7, DataFrame(isolation_level = occ))];
        made_date = Date("2026-06-13")
    )
    @test Set(arch0.stream) == Set(["isolation beds"])
end
