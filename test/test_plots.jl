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
    rows = [
        ("2026-05-18", 925, 419, 2075),
        ("2026-05-23", 1364, 680, 3137),
        ("2026-05-28", 3510, 2196, 6325)
    ]
    ## Released series only.
    @test plot_estimate_evolution(rows) isa CairoMakie.Makie.Figure
    ## With a renewal frozen-fit ribbon and a scenario range band.
    renewal = [
        ("2026-05-20", 1666, 900, 2900),
        ("2026-05-23", 1900, 1000, 3300)
    ]
    fig = plot_estimate_evolution(rows; renewal = renewal,
        scenario_range = (235, 1386))
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
    using BVDOutbreakSize: plot_rt, knot_days
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
            FlexiChains.Parameter(Symbol("rt_state.log_R0_seed")) => reshape(
                log.(1.0 .+ abs.(randn(rng, ndraws)) .* 0.2), ndraws, 1),
            FlexiChains.Parameter(Symbol("rt_state.sigma_rw")) => reshape(
                abs.(randn(rng, ndraws)) .* 0.02, ndraws, 1),
            FlexiChains.Parameter(Symbol("rt_state.intervention_effect")) => reshape(
                -abs.(randn(rng, ndraws)) .* 0.3, ndraws, 1),
            FlexiChains.Parameter(Symbol("rt_state.z")) => zcol,
            FlexiChains.Parameter(:T) => reshape(abs.(randn(rng, ndraws)) .* 10 .+ 40, ndraws, 1)))
    fig = plot_rt(chn; n = n, breakpoint = n - 11,
        as_of_date = "2026-05-28", seeding = Date("2026-02-23"))
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

@testitem "plot_forecast_vs_truth returns a Makie figure" setup=[HeadlessMakie] begin
    using Random: MersenneTwister
    using DataFrames: DataFrame
    using BVDOutbreakSize: plot_forecast_vs_truth
    rng = MersenneTwister(33)
    n = 300
    fc = DataFrame(
        cases_cum = rand(rng, 50:150, n),
        deaths_cum = rand(rng, 40:100, n),
        exports_cum = rand(rng, 2:10, n),
        cases_new = rand(rng, 0:30, n),
        deaths_new = rand(rng, 0:20, n),
        exports_new = rand(rng, 0:5, n)
    )
    fig = plot_forecast_vs_truth(fc;
        cases = 130, deaths = 80, exports = 7)
    @test fig isa CairoMakie.Makie.Figure
end
