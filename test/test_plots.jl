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
        pp_exports_deaths = rand(rng, 0:3, 400), obs_exports_deaths = 1)
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

@testitem "plot_start_date_pair returns a Makie figure" setup=[HeadlessMakie] begin
    using Random: MersenneTwister
    import FlexiChains
    using BVDOutbreakSize: plot_start_date_pair
    rng = MersenneTwister(16)
    n = 200
    vals = hcat(abs.(randn(rng, n)) .+ 7, abs.(randn(rng, n)) .* 30)
    chn = FlexiChains.FlexiChain{Symbol}(n,
        1,
        Dict(
            FlexiChains.Parameter(:τ) => reshape(vals[:, 1], n, 1),
            FlexiChains.Parameter(:T) => reshape(vals[:, 2], n, 1)))
    fig = plot_start_date_pair(chn; as_of_date = "2026-05-20")
    @test fig isa CairoMakie.Makie.Figure
end

@testitem "plot_vintage_conditional_ppc figure" setup=[HeadlessMakie] begin
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

@testitem "conditional one-step-ahead predictive construction" begin
    ## The conditional one-step-ahead cumulative at vintage v is the
    ## observed previous cumulative plus the drawn increment: ŷ_v =
    ## y_{v-1} + Δ_v with y_0 = 0. We reproduce the construction used by
    ## `plot_vintage_conditional_ppc` on a tiny synthetic example and
    ## check it against the closed form, including the bin-1 baseline 0.
    observed = [18, 27, 39, 46]          # cumulative y_1..y_n
    n = length(observed)
    obs_cum = float.(observed)
    obs_prev = [v == 1 ? 0.0 : obs_cum[v - 1] for v in 1:n]
    @test obs_prev == [0.0, 18.0, 27.0, 39.0]   # y_{v-1}, y_0 = 0

    increments = [3.0, 5.0, 8.0, 2.0]    # Δ_1..Δ_n for one draw
    cond = obs_prev .+ increments
    ## Bin 1 conditions on baseline 0, the rest on the observed previous
    ## cumulative; nothing here is a running sum of modelled increments.
    @test cond[1] == 0.0 + increments[1]
    for v in 2:n
        @test cond[v] == observed[v - 1] + increments[v]
    end
    ## Sanity: the conditional reconstruction is not the unconditional
    ## cumulative sum of the increments unless the data matched exactly.
    @test cond != cumsum(increments)
end
