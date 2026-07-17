## Smoke tests for the plotting functions. We do not compare
## pixels — we only check that each call returns a renderable object
## without throwing. CairoMakie is activated headless via the
## HeadlessMakie testsnippet.

@testsnippet HeadlessMakie begin
    using CairoMakie
    CairoMakie.activate!(type = "png")
end

@testitem "plot_cumulative_cases returns a figure-grid" setup=[HeadlessMakie] begin
    using Random: MersenneTwister
    using BVDOutbreakSize: plot_cumulative_cases
    rng = MersenneTwister(4)
    a = randn(rng, 300) .* 50 .+ 400
    b = randn(rng, 300) .* 80 .+ 600
    fg = plot_cumulative_cases("fit A" => a, "fit B" => b; xmax = 1_500)
    @test fg !== nothing
    # AlgebraOfGraphics.draw returns a FigureGrid wrapping a Makie Figure.
    @test fg.figure isa CairoMakie.Makie.Figure
end

@testitem "plot_density_overlay returns a figure-grid" setup=[HeadlessMakie] begin
    using Random: MersenneTwister
    using BVDOutbreakSize: plot_density_overlay
    rng = MersenneTwister(14)
    a = randn(rng, 300) .* 5 .+ 50
    b = randn(rng, 300) .* 5 .+ 40
    fg = plot_density_overlay("fit A" => a, "fit B" => b;
        xlabel = "Seeding time", title = "by clock rate")
    @test fg !== nothing
    @test fg.figure isa CairoMakie.Makie.Figure
end

@testitem "plot_posterior_predictive returns a Makie figure" setup=[HeadlessMakie] begin
    using Random: MersenneTwister
    using BVDOutbreakSize: plot_posterior_predictive
    rng = MersenneTwister(5)
    pp_exports = rand(rng, 0:10, 500)
    pp_deaths = rand(rng, 0:5, 500)
    fig = plot_posterior_predictive(pp_exports, pp_deaths, 3, 1)
    @test fig isa CairoMakie.Makie.Figure
end

@testitem "plot_posterior_predictive lays out four streams" setup=[HeadlessMakie] begin
    using Random: MersenneTwister
    using BVDOutbreakSize: plot_posterior_predictive
    rng = MersenneTwister(15)
    fig = plot_posterior_predictive(
        rand(rng, 0:10, 400), rand(rng, 0:60, 400), 3, 40;
        pp_cases = rand(rng, 0:30, 400), obs_cases = 20,
        pp_exports_deaths = rand(rng, 0:3, 400), obs_exports_deaths = 1,
        pp_confirmed_deaths = rand(rng, 0:30, 400),
        obs_confirmed_deaths = 17)
    @test fig isa CairoMakie.Makie.Figure
end

@testitem "plot_prior_predictive returns a Makie figure" setup=[HeadlessMakie] begin
    using Random: MersenneTwister
    using BVDOutbreakSize: plot_prior_predictive
    rng = MersenneTwister(6)
    pp_exports = rand(rng, 0:10, 500)
    pp_deaths = rand(rng, 0:5, 500)
    fig = plot_prior_predictive(pp_exports, pp_deaths, 3, 1)
    @test fig isa CairoMakie.Makie.Figure
end

@testitem "plot_posterior_predictive_grid lays out four columns" setup=[HeadlessMakie] begin
    using Random: MersenneTwister
    using BVDOutbreakSize
    rng = MersenneTwister(17)
    streams = (; exports = rand(rng, 0:10, 300),
        exports_deaths = rand(rng, 0:3, 300),
        deaths = rand(rng, 0:60, 300),
        cases = rand(rng, 0:30, 300))
    observed = (; exports = 2, exports_deaths = 1,
        deaths = 40, cases = 20)
    fig = BVDOutbreakSize.plot_posterior_predictive_grid(;
        individual = streams, joint = streams, observed = observed)
    @test fig isa CairoMakie.Makie.Figure
end

@testitem "plot_posterior_predictive_grid floats vector draws" setup=[HeadlessMakie] begin
    ## predict returns vector-valued observations (e.g. per-vintage
    ## total_deaths) as a Vector{Vector{Int}}; rendering must not throw
    ## isfinite(::Vector{Int}) under Makie 0.24 / AoG 0.12.
    using Random: MersenneTwister
    using BVDOutbreakSize
    rng = MersenneTwister(23)
    vecdraws(n, r) = [rand(rng, r, 1) for _ in 1:n]
    streams = (; exports = rand(rng, 0:10, 300),
        exports_deaths = rand(rng, 0:3, 300),
        deaths = vecdraws(300, 0:60),
        cases = vecdraws(300, 0:30))
    observed = (; exports = 2, exports_deaths = 1,
        deaths = 40, cases = 20)
    fig = BVDOutbreakSize.plot_posterior_predictive_grid(;
        individual = streams, joint = streams, observed = observed)
    @test fig isa CairoMakie.Makie.Figure
    ## Saving forces Makie to compute data limits, where the
    ## isfinite regression manifested.
    path = tempname() * ".png"
    CairoMakie.save(path, fig)
    @test isfile(path)
end

@testitem "plot_pair returns a renderable object" setup=[HeadlessMakie] begin
    using Distributions: Normal
    using Turing: @model, sample, Prior
    import FlexiChains
    using BVDOutbreakSize: plot_pair

    ## kept: plot_pair only needs two named parameters; using a real
    ## BVD model would bloat the test without changing what it checks.
    @model function _plot_model()
        a ~ Normal(0.0, 1.0)
        b ~ Normal(2.0, 0.5)
    end

    chn = sample(_plot_model(), Prior(), 200;
        chain_type = FlexiChains.VNChain, progress = false)
    obj = plot_pair(chn, [:a, :b]; thin = 4)
    @test obj !== nothing
end

@testitem "plot_pair overlays a prior series" setup=[HeadlessMakie] begin
    using Distributions: Normal
    using Turing: @model, sample, Prior
    import FlexiChains
    using BVDOutbreakSize: plot_pair

    ## kept: plot_pair only needs two named parameters; using a real
    ## BVD model would bloat the test without changing what it checks.
    @model function _plot_model()
        a ~ Normal(0.0, 1.0)
        b ~ Normal(2.0, 0.5)
    end

    chn = sample(_plot_model(), Prior(), 200;
        chain_type = FlexiChains.VNChain, progress = false)
    obj = plot_pair(chn, [:a, :b]; thin = 4, prior = chn)
    @test obj !== nothing
end

@testitem "plot_correlation_heatmap returns a Makie figure" setup=[HeadlessMakie] begin
    using Distributions: Normal
    using Turing: @model, sample, Prior
    import FlexiChains
    using BVDOutbreakSize: plot_correlation_heatmap

    @model function _corr_model()
        a ~ Normal(0.0, 1.0)
        b ~ Normal(2.0, 0.5)
        c ~ Normal(-1.0, 2.0)
    end

    chn = sample(_corr_model(), Prior(), 200;
        chain_type = FlexiChains.VNChain, progress = false)
    fig = plot_correlation_heatmap(chn, [:a, :b, :c];
        labels = Dict(:a => "A", :b => "B"))
    @test fig isa CairoMakie.Makie.Figure
end

@testitem "plot_stream_pairs returns a renderable object" setup=[HeadlessMakie] begin
    using Random: MersenneTwister
    using BVDOutbreakSize: plot_stream_pairs
    rng = MersenneTwister(7)
    modelled = (; cases = randn(rng, 300) .* 50 .+ 1000,
        deaths = randn(rng, 300) .* 20 .+ 250)
    observed = (; cases = 1077.0, deaths = 246.0)
    obj = plot_stream_pairs(modelled, observed)
    @test obj !== nothing
end

@testitem "plot_estimate_comparison returns a Makie figure" setup=[HeadlessMakie] begin
    using BVDOutbreakSize: plot_estimate_comparison
    rows = [
        ("Source A", 313, 39, 870),
        ("Source B", 501, 402, 612),
        ("Our model", 240, 150, 400)
    ]
    fig = plot_estimate_comparison(rows)
    @test fig isa CairoMakie.Makie.Figure
end

@testitem "plot_estimate_evolution returns a Makie figure" setup=[HeadlessMakie] begin
    using BVDOutbreakSize: plot_estimate_evolution
    ## Each tuple is (date, median, lo30, hi30, lo60, hi60, lo90, hi90).
    rows = [
        ("2026-05-18", 925, 765, 1095, 628, 1378, 438, 2234),
        ("2026-05-23", 1364, 1142, 1688, 915, 2128, 656, 3385),
        ("2026-05-28", 3510, 3135, 3969, 2750, 4602, 2231, 6103)
    ]
    ## Released series only, drawn as discrete per-date marks.
    @test plot_estimate_evolution(rows) isa CairoMakie.Makie.Figure
    ## With the discrete frozen-fit renewal marks and the current-data,
    ## current-model trajectory ribbon over the date grid.
    renewal = [
        ("2026-05-20", 1666, 1400, 2000, 1100, 2400, 900, 2900),
        ("2026-05-23", 1900, 1600, 2300, 1300, 2700, 1000, 3300),
        ("2026-05-28", 4000, 3500, 4600, 3000, 5200, 2500, 6000)
    ]
    ## `trajectory` is `(dates, lo30, hi30, lo60, hi60, lo90, hi90)`.
    trajectory = (
        ["2026-05-20", "2026-05-23", "2026-05-28"],
        [1500, 1800, 3800], [2100, 2400, 4400],
        [1200, 1400, 3200], [2500, 2800, 5000],
        [900, 1000, 2600], [3000, 3400, 5800])
    fig = plot_estimate_evolution(rows; renewal = renewal,
        trajectory = trajectory)
    @test fig isa CairoMakie.Makie.Figure

    ## Marks sharing a date (two released estimates at one cut-off, and a
    ## release plus its frozen re-fit) must dodge without error.
    same_date = [
        ("2026-05-18", 972, 813, 1170, 668, 1510, 478, 2437),
        ("2026-05-18", 925, 765, 1095, 628, 1378, 438, 2234),
        ("2026-06-07", 4509, 4250, 5176, 3845, 5931, 3213, 7665),
        ("2026-06-07", 4161, 3498, 5018, 2958, 6822, 2345, 12421)
    ]
    same_date_renewal = [
        ("2026-05-18", 1100, 900, 1300, 750, 1500, 600, 1900)
    ]
    @test plot_estimate_evolution(same_date; renewal = same_date_renewal) isa
          CairoMakie.Makie.Figure
end

@testitem "plot_evolution_by_group empty and filled" setup=[HeadlessMakie] begin
    using BVDOutbreakSize: plot_evolution_by_group
    ## Every group empty is the state before any per-dataset estimate is
    ## saved: the guard must return the note figure, not throw.
    empty_groups = ["joint" => NamedTuple[], "cases" => NamedTuple[]]
    @test plot_evolution_by_group(empty_groups) isa CairoMakie.Makie.Figure
    ## Each tuple is (date, median, lo30, hi30, lo60, hi60, lo90, hi90).
    joint = [
        ("2026-06-21", 1.2, 1.1, 1.3, 1.0, 1.4, 0.9, 1.6),
        ("2026-07-08", 1.5, 1.4, 1.6, 1.3, 1.8, 1.1, 2.1)
    ]
    cases = [
        ("2026-06-21", 1.6, 1.4, 1.8, 1.2, 2.0, 1.0, 2.4),
        ("2026-07-08", 1.9, 1.7, 2.1, 1.5, 2.4, 1.2, 2.9)
    ]
    ## A group with no estimates is dropped rather than drawn as an empty
    ## panel, the recovered case where no individual fit exists.
    groups = ["joint" => joint, "cases" => cases, "recovered" => NamedTuple[]]
    @test plot_evolution_by_group(groups; refline = 1.0) isa
          CairoMakie.Makie.Figure
end

@testitem "plot_forecast_overlay empty and filled" setup=[HeadlessMakie] begin
    using BVDOutbreakSize: plot_forecast_overlay
    using DataFrames: DataFrame
    using Dates: Date, Day
    ## A zero-row overlay is the state before any release stores a forecast:
    ## the guard must return the note figure, not throw.
    empty = DataFrame(stream = String[], made_date = Date[], horizon = Int[],
        target_date = Date[], fit = String[], observed = Float64[],
        median = Float64[], lo90 = Float64[], hi90 = Float64[])
    @test plot_forecast_overlay(empty) isa CairoMakie.Makie.Figure
    ## Filled with the fit roles: confirmed_cases carries baseline, its
    ## individual fit and the joint; recovered carries only baseline and the
    ## joint, so its individual role is simply absent.
    rows = NamedTuple[]
    spec = ["confirmed_cases" => ["baseline", "confirmed", "joint"],
        "recovered" => ["baseline", "joint"]]
    for (stream, fits) in spec, md in [Date(2026, 6, 21), Date(2026, 6, 28)],
        h in [7, 14], fit in fits
        med = 20.0 + h
        push!(rows,
            (; stream = stream, made_date = md, horizon = h,
                target_date = md + Day(h), fit = fit, observed = 18.0 + h,
                median = med, lo90 = med * 0.7, hi90 = med * 1.4))
    end
    @test plot_forecast_overlay(DataFrame(rows)) isa CairoMakie.Makie.Figure
end

@testitem "plot_cumulative_trajectories returns a Makie figure" setup=[HeadlessMakie] begin
    using Random: MersenneTwister
    using Dates: Date
    import FlexiChains
    using BVDOutbreakSize: plot_cumulative_trajectories
    rng = MersenneTwister(21)
    ndraws = 60
    n = 40
    ## Each cumulative trajectory deterministic is a draws×chains matrix of
    ## per-draw monotone vectors.
    _traj() = reshape(
        [cumsum(abs.(randn(rng, n))) for _ in 1:ndraws], ndraws, 1)
    chn = FlexiChains.FlexiChain{Symbol}(ndraws, 1,
        Dict(
            FlexiChains.Parameter(:cumulative_infections) => _traj(),
            FlexiChains.Parameter(:cumulative_onsets) => _traj(),
            FlexiChains.Parameter(:cumulative_expected_deaths) => _traj()))
    fig = plot_cumulative_trajectories(chn; n = n,
        seeding = Date("2026-02-23"))
    @test fig isa CairoMakie.Makie.Figure
end

@testitem "plot_start_date_pair returns a Makie figure" setup=[HeadlessMakie] begin
    using Random: MersenneTwister
    import FlexiChains
    using BVDOutbreakSize: plot_start_date_pair
    rng = MersenneTwister(16)
    n = 200
    vals = hcat(abs.(randn(rng, n)) .+ 7, abs.(randn(rng, n)) .* 30)
    ## :doubling_time replaces the removed :τ parameter
    chn = FlexiChains.FlexiChain{Symbol}(n,
        1,
        Dict(
            FlexiChains.Parameter(:doubling_time) => reshape(vals[:, 1], n, 1),
            FlexiChains.Parameter(:T) => reshape(vals[:, 2], n, 1)))
    fig = plot_start_date_pair(chn; as_of_date = "2026-05-20")
    @test fig isa CairoMakie.Makie.Figure
end

@testitem "plot_rt reconstructs and returns a Makie figure" setup=[HeadlessMakie] begin
    using Random: MersenneTwister
    using Dates: Date
    import FlexiChains
    using BVDOutbreakSize: plot_rt, knot_days, RT_INTERVENTION_RAMP
    rng = MersenneTwister(17)
    ndraws = 120
    n = 95
    nz = length(knot_days(n; week = 7)) - 1
    ## Vector-valued `rt_state.z`: one innovation vector per draw, stored as
    ## a draws×chains matrix of vectors (as the predictive chain returns it).
    zcol = reshape([randn(rng, nz) for _ in 1:ndraws], ndraws, 1)
    chn = FlexiChains.FlexiChain{Symbol}(ndraws, 1,
        Dict(
            FlexiChains.Parameter(Symbol("rt_state.log_R0")) => reshape(
                log.(1.0 .+ abs.(randn(rng, ndraws))), ndraws, 1),
            FlexiChains.Parameter(Symbol("rt_state.sigma_rw")) => reshape(
                abs.(randn(rng, ndraws)) .* 0.02, ndraws, 1),
            FlexiChains.Parameter(Symbol("rt_state.intervention_effect")) => reshape(
                -abs.(randn(rng, ndraws)) .* 0.3, ndraws, 1),
            FlexiChains.Parameter(Symbol("rt_state.z")) => zcol,
            FlexiChains.Parameter(:T) =>
                reshape(abs.(randn(rng, ndraws)) .* 10 .+ 40, ndraws, 1)))
    fig = plot_rt(chn; n = n, breakpoint = n - 11,
        as_of_date = "2026-05-28", seeding = Date("2026-02-23"),
        ramp = RT_INTERVENTION_RAMP)
    @test fig isa CairoMakie.Makie.Figure
end

@testitem "plot_rt_streams overlays streams and joint" setup=[HeadlessMakie] begin
    using Random: MersenneTwister
    using Dates: Date
    import FlexiChains
    using BVDOutbreakSize: plot_rt_streams, knot_days, RT_INTERVENTION_RAMP
    rng = MersenneTwister(31)
    ndraws = 80
    n = 95
    ## Build one chain reconstructable by `reconstruct_rt`, parameterised by
    ## the walk start so the joint (walk from the breakpoint lead) and the
    ## per-stream fits (walk from day 1) each get an innovation vector of the
    ## right length.
    function make_chain(walk_start)
        nz = length(knot_days(n; week = 7, start = walk_start)) - 1
        zcol = reshape([randn(rng, nz) for _ in 1:ndraws], ndraws, 1)
        FlexiChains.FlexiChain{Symbol}(ndraws, 1,
            Dict(
                FlexiChains.Parameter(Symbol("rt_state.log_R0")) => reshape(
                    log.(1.0 .+ abs.(randn(rng, ndraws))), ndraws, 1),
                FlexiChains.Parameter(Symbol("rt_state.sigma_rw")) => reshape(
                    abs.(randn(rng, ndraws)) .* 0.02, ndraws, 1),
                FlexiChains.Parameter(
                    Symbol("rt_state.intervention_effect")) => reshape(
                    -abs.(randn(rng, ndraws)) .* 0.3, ndraws, 1),
                FlexiChains.Parameter(Symbol("rt_state.z")) => zcol))
    end
    breakpoint = n - 11
    joint_walk = breakpoint - 14
    streams = [
        (; label = "cases", chn = make_chain(1), rt_start = 1,
            rt_walk_start = 1, colour = :steelblue),
        (; label = "deaths", chn = make_chain(1), rt_start = 1,
            rt_walk_start = 1, colour = :firebrick)]
    joint = (; label = "joint", chn = make_chain(joint_walk),
        rt_start = joint_walk, rt_walk_start = joint_walk)
    fig = plot_rt_streams(streams; joint = joint, n = n,
        breakpoint = breakpoint, as_of_date = "2026-05-28",
        seeding = Date("2026-02-23"), display_start = joint_walk,
        ramp = RT_INTERVENTION_RAMP)
    @test fig isa CairoMakie.Makie.Figure
end

@testitem "plot_vintage_conditional_ppc returns a Makie figure" setup=[HeadlessMakie] begin
    using Random: MersenneTwister
    using BVDOutbreakSize: plot_vintage_conditional_ppc
    rng = MersenneTwister(21)
    dates = ["2026-05-18", "2026-05-19", "2026-05-20",
        "2026-05-21", "2026-05-22", "2026-05-23"]
    ## Per-draw per-bin increment vectors, as the predictive chain
    ## returns them (here a plain vector of draws).
    reps = [rand(rng, 1:30, length(dates)) for _ in 1:150]
    observed = cumsum([18, 9, 12, 7, 6, 5])
    fig = plot_vintage_conditional_ppc([
        (; title = "Suspected", dates = dates,
            replicates = reps, observed = observed, colour = :steelblue),
        (; title = "Confirmed", dates = dates,
            replicates = reps, observed = observed)])
    @test fig isa CairoMakie.Makie.Figure
end

@testitem "plot_vintage_conditional_ppc draws a daily (non-cumulative) panel" setup=[HeadlessMakie] begin
    using Random: MersenneTwister
    using BVDOutbreakSize: plot_vintage_conditional_ppc
    rng = MersenneTwister(22)
    ## The daily new-suspect inflow: per-day counts, NOT cumulated. With
    ## `cumulative = false` the observed are the raw daily counts and each
    ## replicate is its own daily count (no running baseline).
    dates = ["2026-06-04", "2026-06-05", "2026-06-06", "2026-06-07"]
    reps = [rand(rng, 80:180, length(dates)) for _ in 1:150]
    observed = [153, 119, 117, 94]
    fig = plot_vintage_conditional_ppc([
        (; title = "New suspects/day", dates = dates,
        replicates = reps, observed = observed,
        colour = :slateblue, cumulative = false)])
    @test fig isa CairoMakie.Makie.Figure
end

@testitem "plot_vintage_incidence_ppc returns a Makie figure" setup=[HeadlessMakie] begin
    using Random: MersenneTwister
    using BVDOutbreakSize: plot_vintage_incidence_ppc
    rng = MersenneTwister(23)
    dates = ["2026-05-18", "2026-05-19", "2026-05-20",
        "2026-05-21", "2026-05-22", "2026-05-23"]
    reps = [rand(rng, 1:30, length(dates)) for _ in 1:150]
    ## A cumulative panel (observed is the running total; the incidence view
    ## differences it) and a non-cumulative daily panel (observed already a
    ## per-vintage count).
    cum_observed = cumsum([18, 9, 12, 7, 6, 5])
    daily_observed = [18, 9, 12, 7, 6, 5]
    fig = plot_vintage_incidence_ppc([
        (; title = "Suspected", dates = dates,
            replicates = reps, observed = cum_observed, colour = :steelblue),
        (; title = "New suspects/day", dates = dates,
            replicates = reps, observed = daily_observed,
            cumulative = false)])
    @test fig isa CairoMakie.Makie.Figure
end

@testitem "plot_cfr_prior returns a Makie figure" setup=[HeadlessMakie] begin
    using Distributions: Beta
    using BVDOutbreakSize: plot_cfr_prior
    prior = Beta(6.6, 13.4)
    fig = plot_cfr_prior(prior)
    @test fig isa CairoMakie.Makie.Figure
end

@testitem "plot_no_onward_deaths returns a Makie figure" setup=[HeadlessMakie] begin
    using Random: MersenneTwister
    using DataFrames: DataFrame
    using BVDOutbreakSize: plot_no_onward_deaths
    rng = MersenneTwister(31)
    df = DataFrame(
        delta_deaths = abs.(randn(rng, 300)) .* 5,
        total_projected = abs.(randn(rng, 300)) .* 5 .+ 55
    )
    fig = plot_no_onward_deaths(df; obs_deaths = 55)
    @test fig isa CairoMakie.Makie.Figure
end

@testitem "plot_forecast returns a Makie figure" setup=[HeadlessMakie] begin
    using Random: MersenneTwister
    using DataFrames: DataFrame
    using BVDOutbreakSize: plot_forecast
    rng = MersenneTwister(32)
    n = 300
    fc = DataFrame(
        cases_cum = rand(rng, 50:150, n),
        deaths_cum = rand(rng, 40:100, n),
        confirmed_cum = rand(rng, 20:80, n),
        confirmed_deaths_cum = rand(rng, 1:20, n),
        cases_new = rand(rng, 0:30, n),
        deaths_new = rand(rng, 0:20, n),
        confirmed_new = rand(rng, 0:15, n),
        confirmed_deaths_new = rand(rng, 0:5, n),
        infections_new = abs.(randn(rng, n)) .* 500,
        rt_forecast = 1.0 .+ abs.(randn(rng, n)) .* 0.5
    )
    fig = plot_forecast(fc)
    @test fig isa CairoMakie.Makie.Figure
end

@testitem "plot_forecast_beds returns a Makie figure" setup=[HeadlessMakie] begin
    using Random: MersenneTwister
    using DataFrames: DataFrame
    using BVDOutbreakSize: plot_forecast_beds
    rng = MersenneTwister(33)
    n = 300
    ## Demand exceeds the supply-limited occupancy, so the shortfall panel is
    ## non-degenerate.
    demand = rand(rng, 400:900, n)
    occ = min.(demand, rand(rng, 300:450, n))
    fc = DataFrame(bed_demand = demand, isolation_level = occ)
    fig = plot_forecast_beds(fc)
    @test fig isa CairoMakie.Makie.Figure
end

@testitem "plot_forecast_beds_vs_truth scores beds against observed" setup=[HeadlessMakie] begin
    using Random: MersenneTwister
    using DataFrames: DataFrame
    using BVDOutbreakSize: plot_forecast_beds_vs_truth
    rng = MersenneTwister(34)
    fc = DataFrame(isolation_level = rand(rng, 250:400, 300))
    @test plot_forecast_beds_vs_truth(fc; isolation = 359) isa
          CairoMakie.Makie.Figure
    ## A missing observed value returns an empty figure rather than erroring.
    @test plot_forecast_beds_vs_truth(fc; isolation = missing) isa
          CairoMakie.Makie.Figure
end

@testitem "plot_forecast_latent returns a Makie figure" setup=[HeadlessMakie] begin
    using Random: MersenneTwister
    using DataFrames: DataFrame
    using BVDOutbreakSize: plot_forecast_latent
    rng = MersenneTwister(34)
    n = 300
    fc = DataFrame(
        infections_new = abs.(randn(rng, n)) .* 500,
        onsets_new = abs.(randn(rng, n)) .* 300,
        deaths_latent_new = abs.(randn(rng, n)) .* 30,
        rt_forecast = 1.0 .+ abs.(randn(rng, n)) .* 0.5
    )
    fig = plot_forecast_latent(fc)
    @test fig isa CairoMakie.Makie.Figure
end

@testitem "plot_forecast_vs_truth_latent returns a Makie figure" setup=[HeadlessMakie] begin
    using Random: MersenneTwister
    using DataFrames: DataFrame
    using BVDOutbreakSize: plot_forecast_vs_truth_latent
    rng = MersenneTwister(35)
    n = 300
    fc = DataFrame(
        infections_new = abs.(randn(rng, n)) .* 500,
        onsets_new = abs.(randn(rng, n)) .* 300,
        deaths_latent_new = abs.(randn(rng, n)) .* 30
    )
    now = (;
        infections_new = abs.(randn(rng, n)) .* 600,
        onsets_new = abs.(randn(rng, n)) .* 350,
        deaths_latent_new = abs.(randn(rng, n)) .* 35)
    fig = plot_forecast_vs_truth_latent(fc; now = now)
    @test fig isa CairoMakie.Makie.Figure
end

@testitem "plot_forecast_vs_truth returns a Makie figure" setup=[HeadlessMakie] begin
    using Random: MersenneTwister
    using DataFrames: DataFrame
    using BVDOutbreakSize: plot_forecast_vs_truth
    rng = MersenneTwister(33)
    n = 300
    fc = DataFrame(
        confirmed_cum = rand(rng, 20:80, n),
        confirmed_deaths_cum = rand(rng, 1:20, n),
        confirmed_new = rand(rng, 0:15, n),
        confirmed_deaths_new = rand(rng, 0:5, n)
    )
    fig = plot_forecast_vs_truth(fc;
        confirmed = 70, confirmed_deaths = 18)
    @test fig isa CairoMakie.Makie.Figure
    ## Confirmed-deaths column absent: the confirmed-deaths panel is dropped
    ## without error, leaving the confirmed-cases panel alone.
    fc_cases = DataFrame(
        confirmed_cum = rand(rng, 20:80, n),
        confirmed_new = rand(rng, 0:15, n))
    fig2 = plot_forecast_vs_truth(fc_cases;
        confirmed = 70, confirmed_deaths = 18)
    @test fig2 isa CairoMakie.Makie.Figure
end

@testitem "plot_projection_comparison returns a Makie figure" setup=[HeadlessMakie] begin
    using BVDOutbreakSize: plot_projection_comparison, CHAMLA_CONFIRMED_CENTRAL
    ## External projection from the packaged Chamla central trajectory; our
    ## projection as a dated fan (ribbon) including a zero-width anchor; observed
    ## a dated value series, all with dates out of order to exercise sorting.
    ours = [("2026-06-24", 1200, 800, 1700), ("2026-05-27", 250, 250, 250),
        ("2026-06-10", 700, 500, 950), ("2026-06-03", 430, 330, 560),
        ("2026-06-17", 930, 660, 1300)]
    observed = [("2026-06-08", 598), ("2026-05-27", 250),
        ("2026-06-23", 1118), ("2026-06-15", 850)]
    fig = plot_projection_comparison(;
        external = CHAMLA_CONFIRMED_CENTRAL[1:4],
        ours = ours, observed = observed)
    @test fig isa CairoMakie.Makie.Figure
end

@testitem "plot_scenario_comparison facets the published scenarios" setup=[HeadlessMakie] begin
    using BVDOutbreakSize: plot_scenario_comparison, REPORT_SCENARIOS_CI
    ## The real scenario set exercises the parser, the dodge of the swept level,
    ## and the geographic/back-calc block layout (18 May has no geographic row).
    ours = Dict("2026-05-18" => (520, 320, 860),
        "2026-05-20" => (760, 470, 1180),
        "2026-05-27" => (1250, 720, 2050))
    fig = plot_scenario_comparison(REPORT_SCENARIOS_CI; ours = ours)
    @test fig isa CairoMakie.Makie.Figure
    ## Renders without an `ours` overlay too (every panel still draws).
    @test plot_scenario_comparison(REPORT_SCENARIOS_CI) isa
          CairoMakie.Makie.Figure
end
