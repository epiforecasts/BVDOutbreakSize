## Tests for the per-province forecast, read off the posterior-predictive
## draws of the fixture patch joint (see test_forecast_horizon.jl).

@testitem "forecast_provinces adds up to the national forecast" setup = [
    HorizonFixtures, ForecastDraws,
] begin
    using BVDOutbreakSize: forecast_provinces, forecast_reported

    for h in (7, 14)
        fc = forecast_provinces(pp; horizon = h, n_patches = NP)
        nat = forecast_reported(pp; horizon = h, OBS...)
        @test size(fc, 1) == NP * ND
        @test names(fc) == [
            "patch", "province", "draw", "infections_new", "rt_forecast",
            "confirmed_new", "confirmed_deaths_new",
        ]
        for d in 1:ND
            rows = fc.draw .== d
            @test sum(fc.confirmed_new[rows]) == nat.confirmed_new[d]
            @test sum(fc.confirmed_deaths_new[rows]) ==
                nat.confirmed_deaths_new[d]
        end
        ## The provinces' infections are the national ones split by patch.
        @test all(isfinite, fc.rt_forecast)
        @test [sum(fc.infections_new[fc.draw .== d]) for d in 1:ND] ≈
            nat.infections_new
    end
    ## The counts are split a week at a time.
    @test_throws ArgumentError forecast_provinces(pp; horizon = 10, n_patches = NP)
end

@testitem "province summaries read the draws or a projection alike" setup = [
    HorizonFixtures, ForecastDraws,
] begin
    using BVDOutbreakSize: forecast_provinces, forecast_reported,
        province_forecast_table

    nat = forecast_reported(pp; horizon = 14, OBS...)
    proj = forecast_provinces(pp; horizon = 14, n_patches = NP)
    tbl = province_forecast_table(pp, nat; horizon = 14, n_patches = NP)
    @test tbl == province_forecast_table(pp, proj; horizon = 14, n_patches = NP)
    @test "New confirmed cases by T+14" in tbl.Quantity
    ## Draws with no province forecast cannot be read by province.
    single = forecast_draws(
        joint(), sample(
            Xoshiro(2), joint(), Prior(), 3;
            chain_type = FlexiChains.VNChain, progress = false
        ); horizon = 7
    )
    @test_throws ArgumentError province_forecast_table(single, nat; n_patches = NP)
end

@testitem "province_forecast_archive archives the draws" setup = [
    HorizonFixtures, ForecastDraws,
] begin
    using Dates: Date, Day
    using DataFrames: nrow
    using BVDOutbreakSize: forecast_provinces, forecast_reported,
        province_forecast_archive, PROVINCE_NAMES

    made = Date("2026-09-06")
    fcs = [(h, forecast_reported(pp; horizon = h, OBS...)) for h in (7, 14)]
    arch = province_forecast_archive(pp, fcs; made_date = made, n_patches = NP)
    @test names(arch) == [
        "made_date", "horizon", "target_date", "province", "stream",
        "draw", "value", "method",
    ]
    ## The method is recorded so scoring can leave out older archives.
    @test all(==("predict"), arch.method)
    @test nrow(arch) == NP * 2 * 2 * ND
    @test Set(arch.province) == Set(PROVINCE_NAMES[1:NP])
    @test all(arch.target_date .== made .+ Day.(arch.horizon))
    for h in (7, 14)
        proj = forecast_provinces(
            pp; horizon = h, n_patches = NP, patch_labels = PROVINCE_NAMES
        )
        got = arch[
            (arch.horizon .== h) .& (arch.province .== PROVINCE_NAMES[1]) .&
                (arch.stream .== "confirmed cases"), :value,
        ]
        @test got == float.(proj[proj.patch .== 1, :confirmed_new])
    end
end

@testitem "province_forecast_vs_truth scores the forecast" setup = [
    HorizonFixtures, ForecastDraws,
] begin
    using Statistics: quantile
    using BVDOutbreakSize: forecast_provinces, forecast_reported,
        province_forecast_vs_truth

    kw = (;
        observed = [900, 100, 30, 5], baseline = [800, 60, 20, 5],
        death_observed = [40, 20, 5, 1], death_baseline = [20, 10, 5, 1],
        n_patches = NP,
    )
    nat = forecast_reported(pp; horizon = 7, OBS...)
    df = province_forecast_vs_truth(pp, nat; kw...)
    @test size(df, 1) == 2 * NP
    cases = df[df[!, "Stream"] .== "Confirmed cases", :]
    @test cases[!, "Observed"] == [100, 40, 10, 0]
    proj = forecast_provinces(pp; n_patches = NP)
    v = float.(proj[proj.patch .== 1, :confirmed_new])
    @test cases[1, "Upper 90%"] == round(quantile(v, 0.95))
    @test province_forecast_vs_truth(pp, proj; kw...) == df
end

@testitem "province forecast summaries read a projection frame" begin
    using DataFrames: DataFrame
    using BVDOutbreakSize: province_forecast_table, PROVINCE_LABELS

    ## A projection frame carries each province's own draws, so the table
    ## reads them directly and needs no shares on the chain.
    fc = DataFrame(
        patch = repeat(1:2; inner = 100),
        province = repeat(PROVINCE_LABELS[1:2]; inner = 100),
        draw = repeat(1:100; outer = 2),
        confirmed_new = vcat(fill(10.0, 100), fill(50.0, 100))
    )
    tbl = province_forecast_table((;), fc; n_patches = 2)
    @test tbl.Province == PROVINCE_LABELS[1:2]
    @test tbl[!, "Lower 90%"] == [10.0, 50.0]
    @test tbl[!, "Upper 90%"] == [10.0, 50.0]

    ## The latent quantities follow each province's observed streams, the
    ## reproduction number to two decimals.
    fc.infections_new = vcat(fill(30.0, 100), fill(90.0, 100))
    fc.rt_forecast = vcat(fill(1.234, 100), fill(0.876, 100))
    tbl = province_forecast_table((;), fc; n_patches = 2)
    @test tbl.Province == repeat(PROVINCE_LABELS[1:2]; inner = 3)
    @test tbl.Quantity == repeat(
        [
            "New confirmed cases by T+7", "New infections by T+7",
            "Reproduction number at T+7",
        ], 2
    )
    @test tbl[!, "Upper 90%"] == [10.0, 30.0, 1.23, 50.0, 90.0, 0.88]
end


@testitem "plot_province_forecast captions the province forecast" setup = [
    HeadlessMakie,
] begin
    using DataFrames: DataFrame
    using CairoMakie: Makie as Mk
    using BVDOutbreakSize: plot_province_forecast, PROVINCE_LABELS

    fc = DataFrame(
        patch = repeat(1:2; inner = 50),
        province = repeat(PROVINCE_LABELS[1:2]; inner = 50),
        draw = repeat(1:50; outer = 2),
        confirmed_new = collect(1.0:100.0)
    )
    fig = plot_province_forecast((;), fc; n_patches = 2)
    captions = [x.text[] for x in fig.content if x isa Mk.Label]
    @test any(c -> occursin("add up to the national forecast", c), captions)
    @test !any(c -> occursin("modelled share", c), captions)
end


@testitem "province_share_draws splits each draw's total" begin
    using DataFrames: DataFrame
    using BVDOutbreakSize: province_share_draws

    fc = DataFrame(
        patch = [1, 2, 1, 2], draw = [1, 1, 2, 2],
        confirmed_new = [30.0, 10.0, 0.0, 0.0]
    )
    sh = province_share_draws(fc, :confirmed_new; n_patches = 2)
    ## Shares are per draw, and a draw with nothing projected anywhere has
    ## no share to give, so it is left out rather than divided by zero.
    @test sh == [[0.75], [0.25]]
end

@testitem "forecast_provinces splits the isolation and bed forecasts" setup = [
    HorizonFixtures,
] begin
    using BVDOutbreakSize: forecast_draws, forecast_provinces,
        forecast_reported, _draw_vectors
    using Turing: sample, Prior
    import FlexiChains
    using Random: Xoshiro

    ## The fixture patch joint with province occupancy and bed rows, so the
    ## forecast carries their splits.
    care = (;
        days = [34, 34, 34, 40, 40, 40], patches = [1, 2, 3, 1, 2, 3],
        counts = [15, 7, 3, 18, 7, 3],
    )
    beds = (; days = [30, 30, 30], patches = [1, 2, 3], counts = [40, 15, 5])
    m = patch_joint(; province_isolation = care, province_capacity = beds)
    nd = 8
    chn = sample(
        Xoshiro(3), m, Prior(), nd;
        chain_type = FlexiChains.VNChain, progress = false
    )
    pp = forecast_draws(m, chn; horizon = 14)
    fc = forecast_provinces(pp; horizon = 7, n_patches = NP)
    nat = forecast_reported(
        pp; horizon = 7, obs_cases = 905, obs_deaths = 18, obs_confirmed = 40
    )
    capacity = _draw_vectors(pp, :forecast_bed_capacity)
    shortfall = _draw_vectors(pp, :forecast_bed_shortfall)
    @test "isolation_level" in names(fc) && "bed_capacity" in names(fc)
    @test "admissions_new" in names(fc)
    national = _draw_vectors(pp, Symbol("forecast_admissions.obs"))
    province = _draw_vectors(pp, Symbol("forecast_admissions.province"))
    for d in 1:nd
        rows = fc.draw .== d
        ## The provinces' patients in isolation add up to the national
        ## occupancy forecast on the same day, and their beds to the
        ## national beds.
        @test round(Int, sum(fc.isolation_level[rows])) ==
            nat.isolation_level[d]
        @test sum(fc.bed_capacity[rows]) ≈ capacity[d][7]
        ## The provinces' admissions add up to the national admissions on
        ## every day.
        @test vec(sum(reshape(province[d], NP, :); dims = 1)) ≈ national[d]
        @test sum(fc.admissions_new[rows]) ≈ sum(national[d][1:7])
        ## The national shortfall column is the forecast shortfall that day.
        @test nat.bed_shortfall[d] == round(shortfall[d][7])
        @test all(>=(0), shortfall[d])
    end
end

@testitem "the province forecast runs the fitted stock forward" begin
    using BVDOutbreakSize: treatment_forecast_model, cutoff_occupancy
    using Turing: returned
    using Random: Xoshiro

    ## Three patches over twelve days, forecast from day 11. Patch 1 holds
    ## 60 of the 100 occupied against 30 beds, so it is 30 above them and
    ## nothing caps it.
    n = 12
    fd = 11:12
    admit_patch = [fill(40.0, 1, n); fill(5.0, 1, n); fill(1.0, 1, n)]
    demand_patch = [fill(60.0, 1, n); fill(20.0, 1, n); fill(20.0, 1, n)]
    capacity_patch = [fill(30.0, 1, n); fill(40.0, 1, n); fill(30.0, 1, n)]
    cut = cutoff_occupancy(100.0, demand_patch[:, 10], capacity_patch[:, 10])
    @test cut.beds == capacity_patch[:, 10]
    @test cut.occupancy ≈ [60.0, 20.0, 20.0]
    @test cut.shortfall ≈ [30.0, 0.0, 0.0]
    flow = fill(1.0, n)
    state = (;
        occupancy_mean = fill(100.0, n), admit_patch, demand_patch,
        capacity_patch, deaths_daily = flow, recover_daily = flow,
        ruleout_daily = flow, abscond_daily = flow, demand = fill(100.0, n),
        admit_daily = vec(sum(admit_patch; dims = 1)),
        occupancy_split_rho = 0.05, admissions_split_rho = 0.05, k_flow = 10.0,
    )
    cap = (; days = [5, 10], counts = [100, 100])
    m = treatment_forecast_model(state, fd, cap, 10.0)
    for seed in 1:20
        r = returned(m, rand(Xoshiro(seed), m))
        ## The means are the fitted model's: the occupancy shared by
        ## demand, the beds the patch capacity and the admissions the
        ## patch admissions, uncapped.
        @test r.held ≈ [fill(60.0, 1, 2); fill(20.0, 1, 2); fill(20.0, 1, 2)]
        @test r.beds == capacity_patch[:, fd]
        @test r.admissions ≈ fill(46.0, 2)
        @test r.isolation ≈ fill(100.0, 2)
        ## Drawn occupancy and admissions sum over the provinces to the
        ## national counts. The national occupancy is censored at the
        ## national beds, as the fitted occupancy is; the admissions are not
        ## censored.
        for d in (r.occupancy_draws, r.admission_draws)
            prov = reshape(d.province, 3, :)
            @test vec(sum(prov; dims = 1)) == round.(Int, d.obs)
        end
        @test all(r.occupancy_draws.obs .<= vec(sum(r.beds; dims = 1)))
        ## Patch 1's 60 stand 30 above its 30 beds, while the national
        ## demand of 100 fills the national 100 beds with none short.
        @test r.province_shortfall ≈ [fill(30.0, 1, 2); zeros(2, 2)]
        @test r.bed_shortfall ≈ [0.0, 0.0]
    end
    ## Admissions are not held below any free beds: a patch over its beds
    ## still draws around its modelled admissions.
    draws = [
        reshape(returned(m, rand(Xoshiro(s), m)).admission_draws.province, 3, :)[1, 1]
            for s in 1:200
    ]
    @test sum(draws) / length(draws) > 30
    ## Half the occupancy fits every patch's beds.
    half = treatment_forecast_model(
        merge(state, (; occupancy_mean = fill(50.0, n))), fd, cap, 10.0
    )
    @test returned(half, rand(Xoshiro(1), half)).province_shortfall ==
        zeros(3, 2)
    ## National demand above the national beds is the national shortfall.
    over_demand = treatment_forecast_model(
        merge(state, (; demand = fill(130.0, n))), fd, cap, 10.0
    )
    @test returned(over_demand, rand(Xoshiro(2), over_demand)).bed_shortfall ≈
        [30.0, 30.0]
    ## A one-day horizon keeps one column per patch.
    one_day = treatment_forecast_model(state, 11:11, cap, 10.0)
    r1 = returned(one_day, rand(Xoshiro(6), one_day))
    @test size(r1.held) == (3, 1)
    @test r1.province_shortfall ≈ reshape([30.0, 0.0, 0.0], 3, 1)
    ## With no recorded capacity nothing is censored.
    nocap = treatment_forecast_model(
        state, fd, (; days = Int[], counts = Int[]), 10.0
    )
    over = [
        any(
            reshape(
                returned(nocap, rand(Xoshiro(s), nocap)).occupancy_draws.province,
                3, :
            )[1, :] .> 30
        ) for s in 1:20
    ]
    @test any(over)
end

@testitem "the province forecast keeps the fitted national spread" begin
    using BVDOutbreakSize: treatment_forecast_model
    using Turing: returned
    using Random: Xoshiro
    using Statistics: mean, var

    ## The fit scores the national occupancy and admissions as one negative
    ## binomial each and the provinces as a split of them, so the forecast's
    ## national counts carry the national NB variance μ + μ²/k, not the
    ## smaller variance of a sum of independent province draws. The
    ## occupancy takes the stock dispersion `k` and the admissions the flow
    ## dispersion `k_flow`, as on the fitted days.
    n = 3
    fd = 3:3
    k = 20.0
    k_flow = 8.0
    μ = 700.0
    level = [0.6, 0.2, 0.2]
    patch = reshape(μ .* level, 3, 1) .* ones(1, n)
    flow = fill(1.0, n)
    state = (;
        occupancy_mean = fill(μ, n), admit_patch = patch,
        demand_patch = patch, capacity_patch = fill(1.0e5, 3, n),
        deaths_daily = flow, recover_daily = flow, ruleout_daily = flow,
        abscond_daily = flow, demand = fill(μ, n), admit_daily = fill(μ, n),
        occupancy_split_rho = 0.05, admissions_split_rho = 0.05, k_flow,
    )
    m = treatment_forecast_model(state, fd, (; days = Int[], counts = Int[]), k)
    rng = Xoshiro(20261022)
    draws = [returned(m, rand(rng, m)) for _ in 1:4000]
    for (field, kd) in ((:occupancy_draws, k), (:admission_draws, k_flow))
        nb_var = μ + μ^2 / kd
        indep_var = μ + sum((μ .* level) .^ 2) / kd
        nat = [only(getfield(d, field).obs) for d in draws]
        @test mean(nat) ≈ μ rtol = 0.02
        @test var(nat) ≈ nb_var rtol = 0.1
        @test var(nat) > 1.5 * indep_var
        ## The split follows the fitted shares and adds up to the national
        ## count.
        prov = reduce(hcat, [getfield(d, field).province for d in draws])
        @test vec(sum(prov; dims = 1)) == round.(Int, nat)
        @test vec(mean(prov; dims = 2)) ./ mean(nat) ≈ level rtol = 0.03
    end
end
