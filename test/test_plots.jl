## Smoke tests for the plotting functions. We do not compare
## pixels — we only check that each call returns a renderable object
## without throwing. CairoMakie is activated headless via the
## HeadlessMakie testsnippet.

@testsnippet HeadlessMakie begin
    using CairoMakie
    CairoMakie.activate!(type = "png")
end

@testitem "draw-vector plots return a renderable figure" setup = [
    HeadlessMakie,
] begin
    using Random: MersenneTwister
    using BVDOutbreakSize: plot_cumulative_cases, plot_density_overlay,
        plot_posterior_predictive, plot_prior_predictive
    rng = MersenneTwister(4)
    a = randn(rng, 300) .* 50 .+ 400
    b = randn(rng, 300) .* 80 .+ 600

    ## AlgebraOfGraphics.draw returns a FigureGrid wrapping a Makie Figure.
    for fg in (
            plot_cumulative_cases("fit A" => a, "fit B" => b; xmax = 1_500),
            plot_density_overlay(
                "fit A" => a, "fit B" => b;
                xlabel = "Seeding time", title = "by clock rate"
            ),
        )
        @test fg.figure isa CairoMakie.Makie.Figure
    end

    ## The prior and posterior predictive panels take the same arguments.
    pp_exports = rand(rng, 0:10, 500)
    pp_deaths = rand(rng, 0:5, 500)
    for f in (plot_posterior_predictive, plot_prior_predictive)
        @test f(pp_exports, pp_deaths, 3, 1) isa CairoMakie.Makie.Figure
    end

    ## The optional streams add panels to the posterior predictive.
    fig = plot_posterior_predictive(
        rand(rng, 0:10, 400), rand(rng, 0:60, 400), 3, 40;
        pp_cases = rand(rng, 0:30, 400), obs_cases = 20,
        pp_exports_deaths = rand(rng, 0:3, 400), obs_exports_deaths = 1,
        pp_confirmed_deaths = rand(rng, 0:30, 400),
        obs_confirmed_deaths = 17
    )
    @test fig isa CairoMakie.Makie.Figure
end

@testitem "plot_posterior_predictive_grid lays out four columns" setup = [
    HeadlessMakie,
] begin
    using Random: MersenneTwister
    using BVDOutbreakSize
    rng = MersenneTwister(17)
    streams = (;
        exports = rand(rng, 0:10, 300),
        exports_deaths = rand(rng, 0:3, 300),
        deaths = rand(rng, 0:60, 300),
        cases = rand(rng, 0:30, 300),
    )
    observed = (;
        exports = 2, exports_deaths = 1,
        deaths = 40, cases = 20,
    )
    fig = BVDOutbreakSize.plot_posterior_predictive_grid(;
        individual = streams, joint = streams, observed = observed
    )
    @test fig isa CairoMakie.Makie.Figure
end

@testitem "plot_posterior_predictive_grid floats vector draws" setup = [
    HeadlessMakie,
] begin
    ## predict returns vector-valued observations (e.g. per-vintage
    ## total_deaths) as a Vector{Vector{Int}}; rendering must not throw
    ## isfinite(::Vector{Int}) under Makie 0.24 / AoG 0.12.
    using Random: MersenneTwister
    using BVDOutbreakSize
    rng = MersenneTwister(23)
    vecdraws(n, r) = [rand(rng, r, 1) for _ in 1:n]
    streams = (;
        exports = rand(rng, 0:10, 300),
        exports_deaths = rand(rng, 0:3, 300),
        deaths = vecdraws(300, 0:60),
        cases = vecdraws(300, 0:30),
    )
    observed = (;
        exports = 2, exports_deaths = 1,
        deaths = 40, cases = 20,
    )
    fig = BVDOutbreakSize.plot_posterior_predictive_grid(;
        individual = streams, joint = streams, observed = observed
    )
    @test fig isa CairoMakie.Makie.Figure
    ## Saving forces Makie to compute data limits, where the
    ## isfinite regression manifested.
    path = tempname() * ".png"
    CairoMakie.save(path, fig)
    @test isfile(path)
end

@testitem "plot_pair returns a renderable object" setup = [HeadlessMakie] begin
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

    chn = sample(
        _plot_model(), Prior(), 200;
        chain_type = FlexiChains.VNChain, progress = false
    )
    obj = plot_pair(chn, [:a, :b]; thin = 4)
    @test obj !== nothing
end

@testitem "plot_pair overlays a prior series" setup = [HeadlessMakie] begin
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

    chn = sample(
        _plot_model(), Prior(), 200;
        chain_type = FlexiChains.VNChain, progress = false
    )
    obj = plot_pair(chn, [:a, :b]; thin = 4, prior = chn)
    @test obj !== nothing
end

@testitem "plot_correlation_heatmap returns a Makie figure" setup = [
    HeadlessMakie,
] begin
    using Distributions: Normal
    using Turing: @model, sample, Prior
    import FlexiChains
    using BVDOutbreakSize: plot_correlation_heatmap

    @model function _corr_model()
        a ~ Normal(0.0, 1.0)
        b ~ Normal(2.0, 0.5)
        c ~ Normal(-1.0, 2.0)
    end

    chn = sample(
        _corr_model(), Prior(), 200;
        chain_type = FlexiChains.VNChain, progress = false
    )
    fig = plot_correlation_heatmap(
        chn, [:a, :b, :c];
        labels = Dict(:a => "A", :b => "B")
    )
    @test fig isa CairoMakie.Makie.Figure
end

@testitem "plot_correlation_heatmap takes named draw vectors" setup = [
    HeadlessMakie,
] begin
    using Random: MersenneTwister
    using BVDOutbreakSize: plot_correlation_heatmap

    ## Per-province quantities are vector deterministics, so the page hands
    ## the heatmap one draw vector per named quantity rather than a chain.
    rng = MersenneTwister(3)
    a = randn(rng, 300)
    draws = (; a = a, b = a .+ 0.1 .* randn(rng, 300), c = randn(rng, 300))
    fig = plot_correlation_heatmap(draws; labels = Dict(:a => "A"))
    @test fig isa CairoMakie.Makie.Figure
end

@testitem "plot_pair takes named draw vectors" setup = [HeadlessMakie] begin
    using Random: MersenneTwister
    using BVDOutbreakSize: plot_pair

    rng = MersenneTwister(4)
    draws = (; a = randn(rng, 200), b = randn(rng, 200))
    @test plot_pair(draws; thin = 2) !== nothing
end

@testitem "plot_stream_pairs returns a renderable object" setup = [
    HeadlessMakie,
] begin
    using Random: MersenneTwister
    using BVDOutbreakSize: plot_stream_pairs
    rng = MersenneTwister(7)
    modelled = (;
        cases = randn(rng, 300) .* 50 .+ 1000,
        deaths = randn(rng, 300) .* 20 .+ 250,
    )
    observed = (; cases = 1077.0, deaths = 246.0)
    obj = plot_stream_pairs(modelled, observed)
    @test obj !== nothing
end

@testitem "plot_estimate_comparison returns a Makie figure" setup = [
    HeadlessMakie,
] begin
    using BVDOutbreakSize: plot_estimate_comparison
    rows = [
        ("Source A", 313, 39, 870),
        ("Source B", 501, 402, 612),
        ("Our model", 240, 150, 400),
    ]
    fig = plot_estimate_comparison(rows)
    @test fig isa CairoMakie.Makie.Figure
end

@testitem "plot_estimate_evolution returns a Makie figure" setup = [
    HeadlessMakie,
] begin
    using BVDOutbreakSize: plot_estimate_evolution
    ## Each tuple is (date, median, lo30, hi30, lo60, hi60, lo90, hi90).
    rows = [
        ("2026-05-18", 925, 765, 1095, 628, 1378, 438, 2234),
        ("2026-05-23", 1364, 1142, 1688, 915, 2128, 656, 3385),
        ("2026-05-28", 3510, 3135, 3969, 2750, 4602, 2231, 6103),
    ]
    ## Released series only, drawn as discrete per-date marks.
    @test plot_estimate_evolution(rows) isa CairoMakie.Makie.Figure
    ## With the discrete frozen-fit renewal marks and the current-data,
    ## current-model trajectory ribbon over the date grid.
    renewal = [
        ("2026-05-20", 1666, 1400, 2000, 1100, 2400, 900, 2900),
        ("2026-05-23", 1900, 1600, 2300, 1300, 2700, 1000, 3300),
        ("2026-05-28", 4000, 3500, 4600, 3000, 5200, 2500, 6000),
    ]
    ## `trajectory` is `(dates, lo30, hi30, lo60, hi60, lo90, hi90)`.
    trajectory = (
        ["2026-05-20", "2026-05-23", "2026-05-28"],
        [1500, 1800, 3800], [2100, 2400, 4400],
        [1200, 1400, 3200], [2500, 2800, 5000],
        [900, 1000, 2600], [3000, 3400, 5800],
    )
    fig = plot_estimate_evolution(
        rows; renewal = renewal,
        trajectory = trajectory
    )
    @test fig isa CairoMakie.Makie.Figure

    ## Marks sharing a date (two released estimates at one cut-off, and a
    ## release plus its frozen re-fit) must dodge without error.
    same_date = [
        ("2026-05-18", 972, 813, 1170, 668, 1510, 478, 2437),
        ("2026-05-18", 925, 765, 1095, 628, 1378, 438, 2234),
        ("2026-06-07", 4509, 4250, 5176, 3845, 5931, 3213, 7665),
        ("2026-06-07", 4161, 3498, 5018, 2958, 6822, 2345, 12421),
    ]
    same_date_renewal = [
        ("2026-05-18", 1100, 900, 1300, 750, 1500, 600, 1900),
    ]
    @test plot_estimate_evolution(same_date; renewal = same_date_renewal) isa
        CairoMakie.Makie.Figure
end

@testitem "plot_estimate_evolution widens a degenerate trajectory" setup = [
    HeadlessMakie,
] begin
    using BVDOutbreakSize: plot_estimate_evolution
    rows = [
        ("2026-06-21", 1.2, 1.1, 1.3, 1.0, 1.4, 0.9, 1.6),
        ("2026-07-08", 1.5, 1.4, 1.6, 1.3, 1.8, 1.1, 2.1),
    ]
    ## Both trajectory dates collapse to the same day, the basic
    ## reproduction number's case where its release history begins at the
    ## current release: a naive band would have zero width and vanish.
    traj = (
        ["2026-07-08", "2026-07-08"], [1.4, 1.4], [1.6, 1.6],
        [1.3, 1.3], [1.7, 1.7], [1.1, 1.1], [1.9, 1.9],
    )
    fig = plot_estimate_evolution(rows; trajectory = traj)
    @test fig isa CairoMakie.Makie.Figure
    ax = only(x for x in fig.content if x isa CairoMakie.Makie.Axis)
    bands = [p for p in ax.scene.plots if p isa CairoMakie.Makie.Band]
    @test length(bands) == 3
    xs = first.(bands[1][1][])
    @test maximum(xs) > minimum(xs)
end

@testitem "plot_evolution_by_group empty and filled" setup = [HeadlessMakie] begin
    using BVDOutbreakSize: plot_evolution_by_group
    ## Every group empty is the state before any per-dataset estimate is
    ## saved: the guard must return the note figure, not throw.
    empty_groups = ["joint" => NamedTuple[], "cases" => NamedTuple[]]
    @test plot_evolution_by_group(empty_groups) isa CairoMakie.Makie.Figure
    ## Each tuple is (date, median, lo30, hi30, lo60, hi60, lo90, hi90).
    joint = [
        ("2026-06-21", 1.2, 1.1, 1.3, 1.0, 1.4, 0.9, 1.6),
        ("2026-07-08", 1.5, 1.4, 1.6, 1.3, 1.8, 1.1, 2.1),
    ]
    cases = [
        ("2026-06-21", 1.6, 1.4, 1.8, 1.2, 2.0, 1.0, 2.4),
        ("2026-07-08", 1.9, 1.7, 2.1, 1.5, 2.4, 1.2, 2.9),
    ]
    ## A group with no estimates is dropped rather than drawn as an empty
    ## panel, the recovered case where no individual fit exists.
    groups = ["joint" => joint, "cases" => cases, "recovered" => NamedTuple[]]
    @test plot_evolution_by_group(groups; refline = 1.0) isa
        CairoMakie.Makie.Figure
end

@testitem "plot_evolution_by_group draws a per-group trajectory" setup = [
    HeadlessMakie,
] begin
    using BVDOutbreakSize: plot_evolution_by_group
    joint = [
        ("2026-06-21", 1.2, 1.1, 1.3, 1.0, 1.4, 0.9, 1.6),
        ("2026-07-08", 1.5, 1.4, 1.6, 1.3, 1.8, 1.1, 2.1),
    ]
    ## A single release, the case a per-dataset history usually starts from
    ## (R22): the group's own trajectory is what makes a one-point panel
    ## readable.
    cases = [("2026-07-08", 1.9, 1.7, 2.1, 1.5, 2.4, 1.2, 2.9)]
    groups = ["joint" => joint, "cases" => cases]

    ## Both trajectory dates collapse to the single release date, the
    ## degenerate case a single-release history produces.
    traj = (
        ["2026-07-08", "2026-07-08"], [1.8, 1.8], [2.0, 2.0],
        [1.6, 1.6], [2.2, 2.2], [1.3, 1.3], [2.6, 2.6],
    )
    fig = plot_evolution_by_group(groups; trajectories = Dict("cases" => traj))
    @test fig isa CairoMakie.Makie.Figure

    axes = [x for x in fig.content if x isa CairoMakie.Makie.Axis]
    cases_ax = only(a for a in axes if a.title[] == "cases")
    joint_ax = only(a for a in axes if a.title[] == "joint")
    cases_bands = [
        p for p in cases_ax.scene.plots
            if p isa CairoMakie.Makie.Band
    ]
    ## The trajectory is drawn in the "cases" panel as a real, non-zero-width
    ## band even though both its dates are the same calendar day.
    @test length(cases_bands) == 3
    xs = first.(cases_bands[1][1][])
    @test maximum(xs) > minimum(xs)
    ## The "joint" panel has no trajectory of its own, so it draws no band.
    @test !any(p -> p isa CairoMakie.Makie.Band, joint_ax.scene.plots)
end

@testitem "plot_evolution_by_group honours shared_yrange=false" setup = [
    HeadlessMakie,
] begin
    using BVDOutbreakSize: plot_evolution_by_group

    ## "joint" reaches far higher than "exports", the outbreak-size case where
    ## a shared axis would squash the small-scale panel to a hairline.
    joint = [
        (
            "2026-07-08", 50000.0, 48000.0, 52000.0,
            46000.0, 54000.0, 44000.0, 56000.0,
        ),
    ]
    exports = [("2026-07-08", 4.0, 3.0, 5.0, 2.0, 6.0, 1.0, 7.0)]
    groups = ["joint" => joint, "exports" => exports]

    fig = plot_evolution_by_group(groups; shared_yrange = false)
    axes = [x for x in fig.content if x isa CairoMakie.Makie.Axis]
    joint_ax = only(a for a in axes if a.title[] == "joint")
    exports_ax = only(a for a in axes if a.title[] == "exports")
    joint_ylim = joint_ax.limits[][2]
    exports_ylim = exports_ax.limits[][2]
    ## Each panel's own upper limit tracks its own data, so the small-scale
    ## panel is not squashed by the large-scale one.
    @test joint_ylim[2] > 10000
    @test exports_ylim[2] < 100

    ## The default (shared_yrange = true) instead gives every panel the
    ## same, joint-dominated limit.
    shared_fig = plot_evolution_by_group(groups)
    shared_axes = [
        x for x in shared_fig.content
            if x isa CairoMakie.Makie.Axis
    ]
    shared_joint = only(a for a in shared_axes if a.title[] == "joint")
    shared_exports = only(a for a in shared_axes if a.title[] == "exports")
    @test shared_joint.limits[][2] == shared_exports.limits[][2]
end

@testitem "plot_forecast_overlay empty and filled" setup = [HeadlessMakie] begin
    using BVDOutbreakSize: plot_forecast_overlay
    using DataFrames: DataFrame
    using Dates: Date, Day
    ## A zero-row overlay is the state before any release stores a forecast:
    ## the guard must return the note figure, not throw.
    empty = DataFrame(
        stream = String[], made_date = Date[], horizon = Int[],
        target_date = Date[], fit = String[], observed = Float64[],
        median = Float64[], lo90 = Float64[], hi90 = Float64[]
    )
    @test plot_forecast_overlay(empty) isa CairoMakie.Makie.Figure
    ## Filled with the fit roles: confirmed_cases carries baseline, its
    ## individual fit and the joint; recovered carries only baseline and the
    ## joint, so its individual role is simply absent.
    rows = NamedTuple[]
    spec = [
        "confirmed_cases" => ["baseline", "confirmed", "joint"],
        "recovered" => ["baseline", "joint"],
    ]
    for (stream, fits) in spec, md in [Date(2026, 6, 21), Date(2026, 6, 28)],
            h in [7, 14], fit in fits
        med = 20.0 + h
        push!(
            rows,
            (;
                stream = stream, made_date = md, horizon = h,
                target_date = md + Day(h), fit = fit, observed = 18.0 + h,
                median = med, lo90 = med * 0.7, hi90 = med * 1.4,
            )
        )
    end
    @test plot_forecast_overlay(DataFrame(rows)) isa CairoMakie.Makie.Figure
end

@testitem "plot_forecast_overlay crops the panel and places the overflow
    marker" setup = [HeadlessMakie] begin
    using BVDOutbreakSize: plot_forecast_overlay
    using DataFrames: DataFrame
    using Dates: Date, Day
    ## One stream and horizon, with a joint interval far wider than the
    ## observed value or either fit's median: the case that can squash
    ## every series in the panel to a line near the bottom.
    md = Date(2026, 6, 21)
    rows = [
        (;
            stream = "confirmed cases", made_date = md, horizon = 7,
            target_date = md + Day(7), fit = "baseline", observed = 21.0,
            median = 20.0, lo90 = 16.0, hi90 = 24.0,
        ),
        (;
            stream = "confirmed cases", made_date = md, horizon = 7,
            target_date = md + Day(7), fit = "joint", observed = 21.0,
            median = 22.0, lo90 = 5.0, hi90 = 6000.0,
        ),
    ]
    fig = plot_forecast_overlay(DataFrame(rows))
    ax = only(x for x in fig.content if x isa CairoMakie.Makie.Axis)
    _, ylims = ax.limits[]
    cap = ylims[2]

    ## Cropped to three times the larger of the observed value (21) and the
    ## largest median (22), not to the joint's much wider 90% interval.
    @test cap ≈ 3.0 * 22.0
    @test cap < 100.0

    ## The marker must sit strictly below the axis's own upper limit, not
    ## coincident with it, or CairoMakie's plot-area clipping cuts the
    ## triangle in half. Makie converts a marker symbol into a path before
    ## it reaches the plot, so the triangle is matched against the
    ## converted path rather than against the `:utriangle` symbol.
    tri = CairoMakie.Makie.convert_attribute(
        :utriangle,
        CairoMakie.Makie.key"marker"(), CairoMakie.Makie.key"scatter"()
    )
    scatters = [x for x in ax.scene.plots if x isa CairoMakie.Makie.Scatter]
    marker = only(s for s in scatters if s.marker[] == tri)
    marker_y = only(unique(last.(marker[1][])))
    @test marker_y < cap
    @test marker_y > 0.9 * cap
end

@testitem "plot_forecast_overlay tick density scales with made dates" setup = [
    HeadlessMakie,
] begin
    using BVDOutbreakSize: plot_forecast_overlay
    using DataFrames: DataFrame
    using Dates: Date, Day
    dates = [Date(2026, 5, 1) + Day(3 * i) for i in 0:14]

    row(md) = (;
        stream = "confirmed cases", made_date = md, horizon = 7,
        target_date = md + Day(7), fit = "joint", observed = 20.0,
        median = 20.0, lo90 = 15.0, hi90 = 25.0,
    )

    ## Few made dates: every one gets its own tick.
    few = DataFrame(row.(dates[1:5]))
    fig_few = plot_forecast_overlay(few)
    ax_few = only(x for x in fig_few.content if x isa CairoMakie.Makie.Axis)
    @test length(ax_few.xticks[][1]) == 5

    ## A moderate release history still gets a tick per made date: the made
    ## dates are evenly spaced rather than placed to calendar scale, so
    ## fifteen labels sit apart from each other rather than overprinting.
    many = DataFrame(row.(dates))
    fig_many = plot_forecast_overlay(many)
    ax_many = only(x for x in fig_many.content if x isa CairoMakie.Makie.Axis)
    @test length(ax_many.xticks[][1]) == length(dates)

    ## A long enough history is thinned, well above the four ticks a fixed
    ## quarter-split would leave but below one per date.
    longer = [Date(2026, 5, 1) + Day(3 * i) for i in 0:39]
    fig_long = plot_forecast_overlay(DataFrame(row.(longer)))
    ax_long = only(x for x in fig_long.content if x isa CairoMakie.Makie.Axis)
    nticks = length(ax_long.xticks[][1])
    @test nticks > 4
    @test nticks < length(longer)
end

@testitem "plot_forecast_overlay draws a frozen fit as joint" setup = [
    HeadlessMakie,
] begin
    using BVDOutbreakSize: plot_forecast_overlay, FROZEN_FIT, BASELINE_FIT
    using DataFrames: DataFrame
    using Dates: Date, Day
    ## A frozen row is the joint model re-fit at a past cut-off, so it
    ## takes the joint colour and legend entry rather than the individual
    ## role's, no individual fit having run at these cut-offs.
    md = Date(2026, 6, 21)
    rows = [
        (;
            stream = "confirmed cases", made_date = md, horizon = 7,
            target_date = md + Day(7), fit = BASELINE_FIT, observed = 21.0,
            median = 20.0, lo90 = 16.0, hi90 = 24.0,
        ),
        (;
            stream = "confirmed cases", made_date = md, horizon = 7,
            target_date = md + Day(7), fit = FROZEN_FIT, observed = 21.0,
            median = 22.0, lo90 = 18.0, hi90 = 26.0,
        ),
    ]
    fig = plot_forecast_overlay(DataFrame(rows))
    leg = only(x for x in fig.content if x isa CairoMakie.Legend)
    labels = [e.label[] for (_, es) in leg.entrygroups[] for e in es]
    @test "joint" in labels
    @test !("individual" in labels)

    ax = only(x for x in fig.content if x isa CairoMakie.Makie.Axis)
    scatters = [x for x in ax.scene.plots if x isa CairoMakie.Makie.Scatter]
    joint_colour = CairoMakie.Makie.to_color(:firebrick)
    individual_colour = CairoMakie.Makie.to_color(:steelblue)
    @test any(s -> s.color[] == joint_colour, scatters)
    @test !any(s -> s.color[] == individual_colour, scatters)
end

@testitem "plot_forecast_relative_skill empty and filled" setup = [
    HeadlessMakie,
] begin
    using BVDOutbreakSize: plot_forecast_relative_skill
    using DataFrames: DataFrame

    empty = DataFrame(
        stream = String[], horizon = Int[], fit = String[],
        rel_to_baseline = Union{Missing, Float64}[]
    )
    @test plot_forecast_relative_skill(empty) isa CairoMakie.Makie.Figure

    rows = NamedTuple[]
    for h in [7, 14, 21], (fit, val) in [("confirmed", 0.8), ("joint", 1.3)]

        push!(
            rows,
            (;
                stream = "confirmed cases", horizon = h, fit = fit,
                rel_to_baseline = val + 0.05 * h,
            )
        )
    end
    ## Recovered has no individual fit, and one guarded cell is missing
    ## (R20), neither of which should break the render.
    push!(
        rows,
        (;
            stream = "recovered", horizon = 7, fit = "joint",
            rel_to_baseline = 1.1,
        )
    )
    push!(
        rows,
        (;
            stream = "recovered", horizon = 14, fit = "joint",
            rel_to_baseline = missing,
        )
    )
    df = DataFrame(rows)
    fig = plot_forecast_relative_skill(df)
    @test fig isa CairoMakie.Makie.Figure
    naxes = count(x -> x isa CairoMakie.Makie.Axis, fig.content)
    @test naxes == 2

    ## The log-scale column plots the same way through `value_col`.
    df2 = copy(df)
    df2.log_rel_to_baseline = df2.rel_to_baseline
    @test plot_forecast_relative_skill(
        df2;
        value_col = :log_rel_to_baseline
    ) isa CairoMakie.Makie.Figure
end

@testitem "plot_forecast_relative_skill draws a frozen fit as joint" setup = [
    HeadlessMakie,
] begin
    using BVDOutbreakSize: plot_forecast_relative_skill, FROZEN_FIT
    using DataFrames: DataFrame

    ## The frozen evaluation scores the joint model alone, so its single
    ## series belongs to the joint role rather than reading as an
    ## individual fit that was never run.
    rows = [
        (;
            stream = "confirmed cases", horizon = h, fit = FROZEN_FIT,
            rel_to_baseline = 0.8 + 0.02 * h,
        ) for h in [7, 14, 21]
    ]
    fig = plot_forecast_relative_skill(DataFrame(rows))
    leg = only(x for x in fig.content if x isa CairoMakie.Legend)
    labels = [e.label[] for (_, es) in leg.entrygroups[] for e in es]
    @test labels == ["joint"]
end

@testitem "plot_cumulative_trajectories returns a Makie figure" setup = [
    HeadlessMakie,
] begin
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
        [cumsum(abs.(randn(rng, n))) for _ in 1:ndraws], ndraws, 1
    )
    chn = FlexiChains.FlexiChain{Symbol}(
        ndraws, 1,
        Dict(
            FlexiChains.Parameter(:cumulative_infections) => _traj(),
            FlexiChains.Parameter(:cumulative_onsets) => _traj(),
            FlexiChains.Parameter(:cumulative_expected_deaths) => _traj()
        )
    )
    fig = plot_cumulative_trajectories(
        chn; n = n,
        seeding = Date("2026-02-23")
    )
    @test fig isa CairoMakie.Makie.Figure
end

@testitem "plot_start_date_pair returns a Makie figure" setup = [
    HeadlessMakie,
] begin
    using Random: MersenneTwister
    import FlexiChains
    using BVDOutbreakSize: plot_start_date_pair
    rng = MersenneTwister(16)
    n = 200
    vals = hcat(abs.(randn(rng, n)) .+ 7, abs.(randn(rng, n)) .* 30)
    ## :doubling_time replaces the removed :τ parameter
    chn = FlexiChains.FlexiChain{Symbol}(
        n,
        1,
        Dict(
            FlexiChains.Parameter(:doubling_time) => reshape(vals[:, 1], n, 1),
            FlexiChains.Parameter(:T) => reshape(vals[:, 2], n, 1)
        )
    )
    fig = plot_start_date_pair(chn; as_of_date = "2026-05-20")
    @test fig isa CairoMakie.Makie.Figure
end

@testitem "plot_rt reconstructs and returns a Makie figure" setup = [
    HeadlessMakie,
] begin
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
    chn = FlexiChains.FlexiChain{Symbol}(
        ndraws, 1,
        Dict(
            FlexiChains.Parameter(Symbol("rt_state.log_R0")) => reshape(
                log.(1.0 .+ abs.(randn(rng, ndraws))), ndraws, 1
            ),
            FlexiChains.Parameter(Symbol("rt_state.sigma_rw")) => reshape(
                abs.(randn(rng, ndraws)) .* 0.02, ndraws, 1
            ),
            FlexiChains.Parameter(Symbol("rt_state.intervention_effect")) =>
                reshape(-abs.(randn(rng, ndraws)) .* 0.3, ndraws, 1),
            FlexiChains.Parameter(Symbol("rt_state.z")) => zcol,
            FlexiChains.Parameter(:T) =>
                reshape(abs.(randn(rng, ndraws)) .* 10 .+ 40, ndraws, 1)
        )
    )
    fig = plot_rt(
        chn; n = n, breakpoint = n - 11,
        as_of_date = "2026-05-28", seeding = Date("2026-02-23"),
        ramp = RT_INTERVENTION_RAMP
    )
    @test fig isa CairoMakie.Makie.Figure
end

@testitem "plot_rt_streams overlays streams and joint" setup = [
    HeadlessMakie,
] begin
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
        FlexiChains.FlexiChain{Symbol}(
            ndraws, 1,
            Dict(
                FlexiChains.Parameter(Symbol("rt_state.log_R0")) => reshape(
                    log.(1.0 .+ abs.(randn(rng, ndraws))), ndraws, 1
                ),
                FlexiChains.Parameter(Symbol("rt_state.sigma_rw")) => reshape(
                    abs.(randn(rng, ndraws)) .* 0.02, ndraws, 1
                ),
                FlexiChains.Parameter(
                    Symbol("rt_state.intervention_effect")
                ) => reshape(
                    -abs.(randn(rng, ndraws)) .* 0.3, ndraws, 1
                ),
                FlexiChains.Parameter(Symbol("rt_state.z")) => zcol
            )
        )
    end
    breakpoint = n - 11
    joint_walk = breakpoint - 14
    streams = [
        (;
            label = "cases", chn = make_chain(1), rt_start = 1,
            rt_walk_start = 1, colour = :steelblue,
        ),
        (;
            label = "deaths", chn = make_chain(1), rt_start = 1,
            rt_walk_start = 1, colour = :firebrick,
        ),
    ]
    joint = (;
        label = "joint", chn = make_chain(joint_walk),
        rt_start = joint_walk, rt_walk_start = joint_walk,
    )
    fig = plot_rt_streams(
        streams; joint = joint, n = n,
        breakpoint = breakpoint, as_of_date = "2026-05-28",
        seeding = Date("2026-02-23"), display_start = joint_walk,
        ramp = RT_INTERVENTION_RAMP
    )
    @test fig isa CairoMakie.Makie.Figure
end

@testitem "vintage PPC panels render cumulative and daily views" setup = [
    HeadlessMakie,
] begin
    using Random: MersenneTwister
    using BVDOutbreakSize: plot_vintage_conditional_ppc,
        plot_vintage_incidence_ppc
    rng = MersenneTwister(21)
    dates = [
        "2026-05-18", "2026-05-19", "2026-05-20",
        "2026-05-21", "2026-05-22", "2026-05-23",
    ]
    ## Per-draw per-bin increment vectors, as the predictive chain returns
    ## them. A cumulative panel takes a running total as `observed`; a
    ## `cumulative = false` panel takes the raw per-day counts.
    reps = [rand(rng, 1:30, length(dates)) for _ in 1:150]
    daily = [18, 9, 12, 7, 6, 5]
    cumulative = cumsum(daily)
    panels = [
        (;
            title = "Suspected", dates = dates,
            replicates = reps, observed = cumulative, colour = :steelblue,
        ),
        (;
            title = "Confirmed", dates = dates,
            replicates = reps, observed = cumulative,
        ),
        (;
            title = "New suspects/day", dates = dates,
            replicates = reps, observed = daily,
            colour = :slateblue, cumulative = false,
        ),
    ]
    for f in (plot_vintage_conditional_ppc, plot_vintage_incidence_ppc)
        @test f(panels) isa CairoMakie.Makie.Figure
    end
end

@testitem "_vintage_ticks labels about one vintage a week" begin
    using Dates: Date, Day
    using BVDOutbreakSize: _vintage_ticks
    ## A long, irregular vintage grid: daily from 14 May with a handful of
    ## days missing, as the situation reports skip days.
    all_dates = [Date("2026-05-14") + Day(i) for i in 0:100]
    skipped = Set(
        Date.(
            [
                "2026-05-24", "2026-05-25", "2026-06-30",
                "2026-07-01", "2026-08-11",
            ]
        )
    )
    dates = string.(filter(!in(skipped), all_dates))
    pos, labels = _vintage_ticks(dates)
    @test length(pos) == length(labels)
    @test issorted(pos)
    @test allunique(pos)
    @test all(1 .<= pos .<= length(dates))
    @test labels == [dates[i] for i in pos]
    ## The last vintage is always labelled, and consecutive labels are at
    ## least a week apart.
    @test last(pos) == length(dates)
    gaps = diff(Date.(labels))
    @test all(g -> g >= Day(7), gaps)
    ## A weekly cadence over a 100-day span is around 15 labels, far below
    ## the ~93 vintages the axis would otherwise carry.
    @test 10 <= length(pos) <= 18
end

@testitem "_vintage_ticks widens the step on a long series" begin
    using Dates: Date, Day
    using BVDOutbreakSize: _vintage_ticks
    ## Three years of weekly vintages. A strict weekly rule would put ~157
    ## labels on the axis, so the step widens in whole weeks instead.
    dates = string.([Date("2026-05-14") + Day(7i) for i in 0:156])
    pos, labels = _vintage_ticks(dates)
    @test length(pos) <= 18
    @test last(pos) == length(dates)
    gaps = diff(Date.(labels))
    @test all(g -> g >= Day(7), gaps)
    @test all(g -> g.value % 7 == 0, gaps)
end

@testitem "_vintage_ticks keeps every label on a short series" begin
    using Dates: Date, Day
    using BVDOutbreakSize: _vintage_ticks
    ## Six consecutive daily vintages span under two weeks, so a weekly
    ## rule would leave two labels. Every vintage is kept instead.
    dates = string.([Date("2026-05-18") + Day(i) for i in 0:5])
    pos, labels = _vintage_ticks(dates)
    @test pos == collect(1:6)
    @test labels == dates
    @test _vintage_ticks(String[]) == (Int[], String[])
end

@testitem "vintage PPC plots take a long series and an empty group" setup = [
    HeadlessMakie,
] begin
    using Dates: Date, Day
    using Random: MersenneTwister
    using BVDOutbreakSize: plot_vintage_conditional_ppc,
        plot_vintage_incidence_ppc
    rng = MersenneTwister(24)
    dates = string.([Date("2026-05-14") + Day(i) for i in 0:100])
    reps = [rand(rng, 1:30, length(dates)) for _ in 1:60]
    observed = cumsum(rand(rng, 1:30, length(dates)))
    panel = (;
        title = "Confirmed", dates = dates,
        replicates = reps, observed = observed,
    )
    @test plot_vintage_conditional_ppc([panel]) isa CairoMakie.Makie.Figure
    @test plot_vintage_incidence_ppc([panel]) isa CairoMakie.Makie.Figure
    ## A stream group with no members (an early cut-off where nothing has
    ## stopped reporting yet) renders as a blank figure.
    @test plot_vintage_conditional_ppc(NamedTuple[]) isa
        CairoMakie.Makie.Figure
    @test plot_vintage_incidence_ppc(NamedTuple[]) isa CairoMakie.Makie.Figure
end

@testitem "plot_cfr_prior returns a Makie figure" setup = [HeadlessMakie] begin
    using Distributions: Beta
    using BVDOutbreakSize: plot_cfr_prior
    prior = Beta(6.6, 13.4)
    fig = plot_cfr_prior(prior)
    @test fig isa CairoMakie.Makie.Figure
end

@testitem "plot_no_onward_deaths returns a Makie figure" setup = [
    HeadlessMakie,
] begin
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

@testitem "plot_forecast returns a Makie figure" setup = [HeadlessMakie] begin
    using Random: MersenneTwister
    using DataFrames: DataFrame
    using BVDOutbreakSize: plot_forecast
    rng = MersenneTwister(32)
    n = 300
    naxes(fig) = count(x -> x isa CairoMakie.Makie.Axis, fig.content)
    ## A frame carrying every observed count stream draws one panel per stream:
    ## reported cases, suspected deaths, confirmed cases, confirmed deaths and
    ## recovered.
    fc = DataFrame(
        cases_cum = rand(rng, 50:150, n),
        deaths_cum = rand(rng, 40:100, n),
        confirmed_cum = rand(rng, 20:80, n),
        confirmed_deaths_cum = rand(rng, 1:20, n),
        cases_new = rand(rng, 0:30, n),
        deaths_new = rand(rng, 0:20, n),
        confirmed_new = rand(rng, 0:15, n),
        confirmed_deaths_new = rand(rng, 0:5, n),
        recovered_new = rand(rng, 0:10, n),
        infections_new = abs.(randn(rng, n)) .* 500,
        rt_forecast = 1.0 .+ abs.(randn(rng, n)) .* 0.5
    )
    fig = plot_forecast(fc)
    @test fig isa CairoMakie.Makie.Figure
    @test naxes(fig) == 5
    ## A single-stream fit that carries only the confirmed columns still
    ## draws just those two panels (backward-compatible with the
    ## confirmed-only frame).
    fc_conf = DataFrame(
        confirmed_new = rand(rng, 0:15, n),
        confirmed_deaths_new = rand(rng, 0:5, n)
    )
    fig_conf = plot_forecast(fc_conf)
    @test fig_conf isa CairoMakie.Makie.Figure
    @test naxes(fig_conf) == 2
    ## A frame carrying no observed count stream returns an empty figure rather
    ## than erroring on the panel-grid layout.
    fig_empty = plot_forecast(DataFrame(rt_forecast = rand(rng, n)))
    @test fig_empty isa CairoMakie.Makie.Figure
    @test naxes(fig_empty) == 0
end

@testitem "plot_forecast_beds returns a Makie figure" setup = [
    HeadlessMakie,
] begin
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

@testitem "plot_forecast_beds_vs_truth scores beds against observed" setup = [
    HeadlessMakie,
] begin
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
    ## An individual-fit forecast overlays a second (dashed) density without
    ## erroring, and without changing the figure type.
    indiv = rand(rng, 200:380, 250)
    @test plot_forecast_beds_vs_truth(
        fc; isolation = 359,
        individual = indiv
    ) isa CairoMakie.Makie.Figure
    ## A degenerate (single-valued) individual sample is skipped rather than
    ## erroring inside `density!`, which needs more than one distinct value.
    @test plot_forecast_beds_vs_truth(
        fc; isolation = 359,
        individual = fill(300.0, 10)
    ) isa CairoMakie.Makie.Figure
    ## An empty individual sample is likewise skipped without erroring.
    @test plot_forecast_beds_vs_truth(
        fc; isolation = 359,
        individual = Float64[]
    ) isa CairoMakie.Makie.Figure
end

@testitem "plot_forecast_latent returns a Makie figure" setup = [
    HeadlessMakie,
] begin
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

@testitem "plot_forecast_latent clips the reproduction number at zero" setup = [
    HeadlessMakie,
] begin
    using Random: MersenneTwister
    using DataFrames: DataFrame
    using BVDOutbreakSize: plot_forecast_latent
    rng = MersenneTwister(36)
    n = 400
    ## A right-skewed forecast reproduction number sitting close to zero, so
    ## the Gaussian kernel reaches past the smallest draw and the estimator
    ## itself spans negative values.
    fc = DataFrame(
        infections_new = abs.(randn(rng, n)) .* 500,
        onsets_new = abs.(randn(rng, n)) .* 300,
        deaths_latent_new = abs.(randn(rng, n)) .* 30,
        rt_forecast = 0.05 .+ abs.(randn(rng, n)) .* 0.3
    )
    fig = plot_forecast_latent(fc)
    axs = [x for x in fig.content if x isa CairoMakie.Makie.Axis]
    rt_label = "Forecast reproduction number (DRC)"
    ax = only(a for a in axs if a.xlabel[] == rt_label)
    dens = only(p for p in ax.scene.plots if p isa CairoMakie.Makie.Density)
    ## The estimator keeps its full support, mass below zero included, and
    ## the axis is what crops it, so the test fails if either half is lost.
    @test CairoMakie.Makie.data_limits(dens).origin[1] < 0
    xlims, _ = ax.limits[]
    @test xlims == (0.0, nothing)
    CairoMakie.Makie.update_state_before_display!(fig)
    @test ax.finallimits[].origin[1] == 0.0
end

@testitem "plot_forecast_vs_truth_latent returns a Makie figure" setup = [
    HeadlessMakie,
] begin
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
        deaths_latent_new = abs.(randn(rng, n)) .* 35,
    )
    fig = plot_forecast_vs_truth_latent(fc; now = now)
    @test fig isa CairoMakie.Makie.Figure
end

@testitem "plot_forecast_vs_truth returns a Makie figure" setup = [
    HeadlessMakie,
] begin
    using Random: MersenneTwister
    using DataFrames: DataFrame
    using BVDOutbreakSize: plot_forecast_vs_truth
    rng = MersenneTwister(33)
    n = 300
    naxes(fig) = count(x -> x isa CairoMakie.Makie.Axis, fig.content)
    fc = DataFrame(
        confirmed_cum = rand(rng, 20:80, n),
        confirmed_deaths_cum = rand(rng, 1:20, n),
        confirmed_new = rand(rng, 0:15, n),
        confirmed_deaths_new = rand(rng, 0:5, n)
    )
    ## Two confirmed streams supplied observed cumulatives: two columns of a
    ## cumulative-and-new panel each, four axes.
    fig = plot_forecast_vs_truth(
        fc;
        observed = (confirmed_cum = 70, confirmed_deaths_cum = 18)
    )
    @test fig isa CairoMakie.Makie.Figure
    @test naxes(fig) == 4
    ## A stream whose cumulative column is absent is dropped without error even
    ## when an observed value is supplied, leaving the confirmed-cases panel
    ## alone (two axes).
    fc_cases = DataFrame(
        confirmed_cum = rand(rng, 20:80, n),
        confirmed_new = rand(rng, 0:15, n)
    )
    fig2 = plot_forecast_vs_truth(
        fc_cases;
        observed = (confirmed_cum = 70, confirmed_deaths_cum = 18)
    )
    @test fig2 isa CairoMakie.Makie.Figure
    @test naxes(fig2) == 2
    ## All five scored streams present with observed cumulatives, one with a
    ## baseline: five columns of two panels give ten axes.
    fc_all = DataFrame(
        cases_cum = rand(rng, 50:150, n),
        cases_new = rand(rng, 0:30, n),
        deaths_cum = rand(rng, 40:100, n),
        deaths_new = rand(rng, 0:20, n),
        confirmed_cum = rand(rng, 20:80, n),
        confirmed_new = rand(rng, 0:15, n),
        confirmed_deaths_cum = rand(rng, 1:20, n),
        confirmed_deaths_new = rand(rng, 0:5, n),
        recovered_cum = rand(rng, 10:60, n),
        recovered_new = rand(rng, 0:10, n)
    )
    fig3 = plot_forecast_vs_truth(
        fc_all;
        observed = (
            cases_cum = 140, deaths_cum = 90, confirmed_cum = 70,
            confirmed_deaths_cum = 18, recovered_cum = 55,
        ),
        baseline = (confirmed_cum = 40,)
    )
    @test fig3 isa CairoMakie.Makie.Figure
    @test naxes(fig3) == 10
    ## A stream present in the frame but absent from `observed` stays absent, so
    ## here recovered adds no panel and only the confirmed streams are drawn.
    fig4 = plot_forecast_vs_truth(
        fc_all;
        observed = (confirmed_cum = 70, confirmed_deaths_cum = 18)
    )
    @test fig4 isa CairoMakie.Makie.Figure
    @test naxes(fig4) == 4
    ## No observed values at all yields an empty figure rather than erroring.
    fig5 = plot_forecast_vs_truth(fc_all; observed = NamedTuple())
    @test fig5 isa CairoMakie.Makie.Figure
    @test naxes(fig5) == 0

    ## An individual-fit forecast overlays a dashed density on both panels
    ## of the streams it covers, without adding axes (the panel count is
    ## driven by `observed`, not by which streams carry an individual
    ## overlay) and without erroring on a stream `individual` has no entry
    ## for (recovered, which has no individual fit).
    fig6 = plot_forecast_vs_truth(
        fc_all;
        observed = (
            cases_cum = 140, deaths_cum = 90, confirmed_cum = 70,
            confirmed_deaths_cum = 18, recovered_cum = 55,
        ),
        baseline = (confirmed_cum = 40,),
        individual = (
            cases_new = rand(rng, 0:30, n),
            confirmed_new = rand(rng, 0:15, n),
        )
    )
    @test fig6 isa CairoMakie.Makie.Figure
    @test naxes(fig6) == 10

    ## A degenerate (single-valued) individual sample for one stream is
    ## skipped rather than erroring inside `density!`.
    fig7 = plot_forecast_vs_truth(
        fc_cases;
        observed = (confirmed_cum = 70, confirmed_deaths_cum = 18),
        individual = (confirmed_new = fill(5.0, 10),)
    )
    @test fig7 isa CairoMakie.Makie.Figure
    @test naxes(fig7) == 2
end

@testitem "forecast validation splits off the streams that stopped" setup = [
    HeadlessMakie,
] begin
    using Random: MersenneTwister
    using Dates: Date, Day
    using DataFrames: DataFrame
    using BVDOutbreakSize: plot_forecast, plot_forecast_vs_truth,
        stream_reporting, stream_forecast_columns
    rng = MersenneTwister(35)
    n = 300
    naxes(fig) = count(x -> x isa CairoMakie.Makie.Axis, fig.content)
    fc = DataFrame(
        cases_cum = rand(rng, 50:150, n), cases_new = rand(rng, 0:30, n),
        deaths_cum = rand(rng, 40:100, n), deaths_new = rand(rng, 0:20, n),
        confirmed_cum = rand(rng, 20:80, n),
        confirmed_new = rand(rng, 0:15, n),
        confirmed_deaths_cum = rand(rng, 1:20, n),
        confirmed_deaths_new = rand(rng, 0:5, n)
    )
    ## The suspected streams stopped being reported 60 days before the
    ## cut-off; the confirmed streams run to it.
    grid = 90
    cutoff = Date(2026, 8, 22)
    obs = (;
        cutoff = cutoff, n = grid,
        reported_history = (; days = [10, grid - 60], counts = [50.0, 90.0]),
        deaths_history = (; days = [10, grid - 60], counts = [5.0, 9.0]),
        confirmed_history = (; days = [10, grid], counts = [20.0, 80.0]),
        confirmed_deaths_history = (;
            days = [10, grid],
            counts = [1.0, 18.0],
        ),
    )
    observed = (
        cases_cum = 1077, deaths_cum = 246, confirmed_cum = 70,
        confirmed_deaths_cum = 18,
    )
    cum_cols = (
        :cases_cum, :deaths_cum, :confirmed_cum,
        :confirmed_deaths_cum,
    )
    reporting = Tuple(c for c in cum_cols if stream_reporting(obs, c))
    stopped = Tuple(c for c in cum_cols if !stream_reporting(obs, c))
    @test reporting == (:confirmed_cum, :confirmed_deaths_cum)
    @test stopped == (:cases_cum, :deaths_cum)
    ## The validation figure draws only the streams still reported, so the
    ## two stale streams contribute no dashed truth rule: two columns of a
    ## cumulative and a new panel each.
    kept = NamedTuple(k => v for (k, v) in pairs(observed) if k in reporting)
    fig = plot_forecast_vs_truth(fc; observed = kept)
    @test naxes(fig) == 4
    ## Their projection is kept as its own figure instead, one panel per
    ## stopped stream and no observation drawn.
    stopped_new = [stream_forecast_columns(c).new for c in stopped]
    fig_stopped = plot_forecast(fc[!, stopped_new])
    @test fig_stopped isa CairoMakie.Makie.Figure
    @test naxes(fig_stopped) == 2
end

@testitem "plot_projection_comparison returns a Makie figure" setup = [
    HeadlessMakie,
] begin
    using BVDOutbreakSize: plot_projection_comparison, CHAMLA_CONFIRMED_CENTRAL
    ## External projection from the packaged Chamla central trajectory; our
    ## projection as a dated fan (ribbon) including a zero-width anchor;
    ## observed a dated value series, all with dates out of order to
    ## exercise sorting.
    ours = [
        ("2026-06-24", 1200, 800, 1700), ("2026-05-27", 250, 250, 250),
        ("2026-06-10", 700, 500, 950), ("2026-06-03", 430, 330, 560),
        ("2026-06-17", 930, 660, 1300),
    ]
    observed = [
        ("2026-06-08", 598), ("2026-05-27", 250),
        ("2026-06-23", 1118), ("2026-06-15", 850),
    ]
    fig = plot_projection_comparison(;
        external = CHAMLA_CONFIRMED_CENTRAL[1:4],
        ours = ours, observed = observed
    )
    @test fig isa CairoMakie.Makie.Figure
end

@testitem "plot_scenario_comparison facets the published scenarios" setup = [
    HeadlessMakie,
] begin
    using BVDOutbreakSize: plot_scenario_comparison, REPORT_SCENARIOS_CI
    ## The real scenario set exercises the parser, the dodge of the swept level,
    ## and the geographic/back-calc block layout (18 May has no geographic row).
    ours = Dict(
        "2026-05-18" => (520, 320, 860),
        "2026-05-20" => (760, 470, 1180),
        "2026-05-27" => (1250, 720, 2050)
    )
    fig = plot_scenario_comparison(REPORT_SCENARIOS_CI; ours = ours)
    @test fig isa CairoMakie.Makie.Figure
    ## Renders without an `ours` overlay too (every panel still draws).
    @test plot_scenario_comparison(REPORT_SCENARIOS_CI) isa
        CairoMakie.Makie.Figure
end

@testitem "vintage PPC plots label an occupancy census as a level" setup = [
    HeadlessMakie,
] begin
    using Random: MersenneTwister
    using BVDOutbreakSize: plot_vintage_conditional_ppc,
        plot_vintage_incidence_ppc
    rng = MersenneTwister(24)
    ## Bed occupancy is a census stock: a level at the end of each report
    ## day, not a count of new events. It shares `cumulative = false` with
    ## the genuine per-day flows, so without its own `ylabel` it would be
    ## labelled "Daily count" here and "New per vintage" in the incidence
    ## view, both of which read as an accumulating total.
    dates = ["2026-06-04", "2026-06-05", "2026-06-06", "2026-06-07"]
    reps = [rand(rng, 200:300, length(dates)) for _ in 1:80]
    occupancy = [258, 267, 283, 260]
    flow = [153, 119, 117, 94]
    stock = (;
        title = "Patients in isolation", dates = dates,
        replicates = reps, observed = occupancy, cumulative = false,
        ylabel = "Beds occupied",
    )
    daily = (;
        title = "New suspects/day", dates = dates,
        replicates = reps, observed = flow, cumulative = false,
    )
    ylabels(fig) = [
        ax.ylabel[] for ax in fig.content
            if ax isa CairoMakie.Makie.Axis
    ]

    ## The census panel keeps its own label in both views; the per-day flow
    ## beside it keeps each view's default.
    cond = plot_vintage_conditional_ppc([stock, daily])
    @test ylabels(cond) == ["Beds occupied", "Daily count"]

    inc = plot_vintage_incidence_ppc([stock, daily])
    @test ylabels(inc) == ["Beds occupied", ""]

    ## A cumulative panel is untouched by the override: it still names the
    ## running total in the first column only.
    cumulative = (;
        title = "Confirmed", dates = dates,
        replicates = reps, observed = cumsum(flow),
    )
    @test ylabels(plot_vintage_conditional_ppc([cumulative])) ==
        ["Cumulative count"]
    @test ylabels(plot_vintage_incidence_ppc([cumulative])) ==
        ["New per vintage"]
end

@testitem "onset_nowcast_draws narrows as the reporting delay grows" begin
    using Random: MersenneTwister
    using Statistics: quantile
    using BVDOutbreakSize: onset_nowcast_draws
    rng = MersenneTwister(4242)
    D, gs, ge = 21, 1, 60
    ndraws = 200
    ## A fitted hazard with posterior spread in every component, so the
    ## nowcast interval has something to be wide about at a short delay.
    hazard = (;
        logit_h0 = [fill(-1.4 + 0.2 * randn(rng), D) for _ in 1:ndraws],
        γ = [zeros(ge) for _ in 1:ndraws],
        alpha = [
            fill(0.4 + 0.05 * randn(rng), ge - gs + 1)
                for _ in 1:ndraws
        ],
    )
    onsets = [fill(150.0 + 30 * randn(rng), ge) for _ in 1:ndraws]
    ## One onset day per delay, all carrying the same observed count so the
    ## only thing separating them is how much reporting has happened.
    delays = [0, 5, 10, D - 1]
    days = [40, 39, 38, 37]
    observed = fill(50.0, length(days))
    draws = onset_nowcast_draws(
        days, observed, delays, onsets, hazard;
        grid_start = gs
    )
    @test length(draws) == length(days)
    @test all(length(d) == ndraws for d in draws)
    ## Never below what is already reported.
    @test all(all(d .>= 50.0 - 1.0e-9) for d in draws)
    ## The interval collapses onto the observed count once the delay has
    ## run out, and widens monotonically as the delay shortens.
    width(d) = quantile(d, 0.95) - quantile(d, 0.05)
    ws = width.(draws)
    @test ws[end] < 1.0e-6
    @test all(diff(ws) .< 0)
    @test all(isapprox.(draws[end], 50.0; atol = 1.0e-6))
    @test_throws ErrorException onset_nowcast_draws(
        days, observed[1:2],
        delays, onsets, hazard; grid_start = gs
    )
    ## Onsets and hazard must be the same fit's draws, paired one to one.
    @test_throws ErrorException onset_nowcast_draws(
        days, observed, delays,
        onsets[1:(ndraws - 1)], hazard; grid_start = gs
    )
    ## A day off the end of the onset series is named rather than left to a
    ## `BoundsError` from inside the draw loop.
    @test_throws ErrorException onset_nowcast_draws(
        [ge + 1], [1.0], [0],
        onsets, hazard; grid_start = gs
    )
    ## `target_delays` stops the prediction at a given delay: nothing
    ## outstanding when it is the delay already reached, and no more than
    ## the eventual total when it is the end of the delay axis.
    same = onset_nowcast_draws(
        days, observed, delays, onsets, hazard;
        grid_start = gs, target_delays = delays
    )
    @test all(all(isapprox.(d, 50.0; atol = 1.0e-9)) for d in same)
    full = onset_nowcast_draws(
        days, observed, delays, onsets, hazard;
        grid_start = gs, target_delays = fill(D - 1, length(days))
    )
    @test all(all(full[k] .<= draws[k] .+ 1.0e-9) for k in eachindex(days))
    @test_throws ErrorException onset_nowcast_draws(
        days, observed, delays,
        onsets, hazard; grid_start = gs, target_delays = delays[1:2]
    )
end

@testitem "plot_onset_nowcast_grid returns a Makie figure" setup = [HeadlessMakie] begin
    using Random: MersenneTwister
    using Dates: Date, Day
    using BVDOutbreakSize: plot_onset_nowcast_grid
    rng = MersenneTwister(4243)
    function panel(title, n)
        dates = [Date("2026-07-01") + Day(i) for i in 0:(n - 1)]
        observed = [40.0 + 10 * randn(rng) for _ in dates]
        nowcast = [
            observed[k] .+ abs.(randn(rng, 120)) .* k
                for k in eachindex(dates)
        ]
        return (; title, dates, observed, nowcast, latest = observed .+ 5)
    end
    fig = plot_onset_nowcast_grid(
        [
            panel("2026-08-01", 30),
            panel("2026-08-08", 34),
        ]
    )
    @test fig isa CairoMakie.Makie.Figure
    ## No digitised snapshots: a blank figure rather than an empty grid.
    @test plot_onset_nowcast_grid([]) isa CairoMakie.Makie.Figure
    p = panel("2026-08-15", 12)
    @test_throws ErrorException plot_onset_nowcast_grid(
        [
            (;
                p.title, p.dates,
                observed = p.observed[1:5], p.nowcast, p.latest,
            ),
        ]
    )
end

@testitem "onset_report_delay_pmf sums to one and stays nonnegative" begin
    using BVDOutbreakSize: onset_report_delay_pmf
    using Random: MersenneTwister

    logit_h0 = randn(MersenneTwister(1), 28) .- 1.0
    γ = 0.3 .* randn(MersenneTwister(2), 10)
    for t in (1, 5, 10)
        pmf = onset_report_delay_pmf(logit_h0, γ, t, 1)
        @test length(pmf) == 28
        @test all(>=(0), pmf)
        @test sum(pmf) ≈ 1.0 atol = 1.0e-8
    end
    ## `t` outside `[grid_start, grid_start + length(γ) - 1]` holds the
    ## calendar effect flat at the walk's nearest edge rather than erroring.
    pmf_before = onset_report_delay_pmf(logit_h0, γ, -5, 1)
    pmf_edge = onset_report_delay_pmf(logit_h0, γ, 1, 1)
    @test pmf_before ≈ pmf_edge
end

@testitem "onset_report_delay_moments matches a hand-computed pmf" begin
    ## D = 2, a flat hazard (logit_h0 both zero) and no calendar effect: the
    ## un-normalised cdf is the truncated-geometric survival product, hand
    ## computable, and the pmf is that normalised by its own last entry
    ## (`onset_report_delay_pmf`'s guard against underflow).
    using BVDOutbreakSize: onset_report_delay_moments, onset_report_delay_pmf
    using StatsFuns: logistic

    logit_h0 = [0.0, 0.0]
    γ = [0.0]
    h = logistic(0.0)
    cdf0 = 1 - (1 - h)
    cdf1 = 1 - (1 - h)^2
    p0 = cdf0 / cdf1
    p1 = (cdf1 - cdf0) / cdf1
    pmf = onset_report_delay_pmf(logit_h0, γ, 1, 1)
    @test pmf ≈ [p0, p1] atol = 1.0e-8

    m = onset_report_delay_moments(logit_h0, γ, 1, 1)
    @test m.mean ≈ 0 * pmf[1] + 1 * pmf[2] atol = 1.0e-8
    @test m.sd ≈
        sqrt(max(0^2 * pmf[1] + 1^2 * pmf[2] - m.mean^2, 0)) atol = 1.0e-8
end

@testitem "onset_level_predictive_draws samples through the missing branch" setup = [
    HeadlessMakie,
] begin
    using BVDOutbreakSize: onset_level_predictive_draws, onsets_only_model,
        reconstruct_onset_hazard, onset_report_F
    using Turing: Prior, sample
    using Random: MersenneTwister

    oc = (;
        onset_days = [10, 11, 12, 13, 10, 11, 12, 13, 14],
        report_days = [15, 15, 15, 15, 20, 20, 20, 20, 20],
        prev_report_days = [0, 0, 0, 0, 15, 15, 15, 15, 0],
        increments = [2, 3, 1, 0, 1, 2, 3, 4, 5],
    )
    n = 40
    chn = sample(
        onsets_only_model(n; onset_curve_history = oc), Prior(), 20;
        progress = false
    )
    grid_start = minimum(oc.onset_days)
    grid_end = maximum(oc.report_days)
    hz = reconstruct_onset_hazard(chn; grid_start, grid_end)
    daily = [
        (v = collect(t); vcat(v[1], diff(v)))
            for t in vec(collect(chn[:cumulative_onsets]))
    ]
    scan_level = [collect(v) for v in vec(collect(chn[:onset_scan_level]))]
    noise_scale = [collect(v) for v in vec(collect(chn[:onset_noise_scale]))]
    u = 12

    draws = onset_level_predictive_draws(
        u, daily, hz, scan_level, noise_scale, 1;
        grid_start, alpha_grid_start = grid_start, n_rep = 5,
        rng = MersenneTwister(11)
    )
    @test length(draws) == 20 * 5
    @test all(isfinite, draws)

    ## Same seed, same replicate: the model call is deterministic given an
    ## explicit `rng`.
    draws_again = onset_level_predictive_draws(
        u, daily, hz, scan_level, noise_scale, 1;
        grid_start, alpha_grid_start = grid_start, n_rep = 5,
        rng = MersenneTwister(11)
    )
    @test draws == draws_again

    ## `target_delay` reduces the truncation the level is read at: for every
    ## draw, the target mean at a short delay is no larger than at the
    ## walk's asymptote (the eventual total, `target_delay = nothing`).
    for i in eachindex(hz.logit_h0)
        a = hz.alpha[i]
        α = a[clamp(u - grid_start + 1, 1, length(a))]
        f_short = onset_report_F(3, hz.logit_h0[i], hz.γ[i], u, grid_start, α)
        f_full = onset_report_F(27, hz.logit_h0[i], hz.γ[i], u, grid_start, α)
        @test f_short <= f_full + 1.0e-10
    end

    ## A vintage index outside the fitted range falls back to a scan
    ## multiplier of one and a noise scale of zero rather than indexing out
    ## of bounds.
    draws_oob = onset_level_predictive_draws(
        u, daily, hz, scan_level, noise_scale, 999;
        grid_start, alpha_grid_start = grid_start,
        rng = MersenneTwister(12)
    )
    @test length(draws_oob) == 20 * 4
    @test all(isfinite, draws_oob)

    @test_throws ErrorException onset_level_predictive_draws(
        999, daily, hz, scan_level, noise_scale, 1;
        grid_start, alpha_grid_start = grid_start
    )
end

@testitem "plot_onset_delay_profile returns a Makie figure" setup = [
    HeadlessMakie,
] begin
    using BVDOutbreakSize: onsets_only_model, reconstruct_onset_hazard,
        plot_onset_delay_profile
    using Turing: Prior, sample
    using Dates: Date

    oc = (;
        onset_days = [10, 11, 12, 13, 10, 11, 12, 13, 14],
        report_days = [15, 15, 15, 15, 20, 20, 20, 20, 20],
        prev_report_days = [0, 0, 0, 0, 15, 15, 15, 15, 0],
        increments = [2, 3, 1, 0, 1, 2, 3, 4, 5],
    )
    chn = sample(
        onsets_only_model(40; onset_curve_history = oc), Prior(), 10;
        progress = false
    )
    grid_start = minimum(oc.onset_days)
    grid_end = maximum(oc.report_days)
    hz = reconstruct_onset_hazard(chn; grid_start, grid_end)
    fig = plot_onset_delay_profile(
        hz; grid_start, grid_end, seeding = Date("2026-01-01")
    )
    @test fig isa CairoMakie.Makie.Figure
end

@testitem "plot_onset_level_band returns a Makie figure" setup = [
    HeadlessMakie,
] begin
    using BVDOutbreakSize: plot_onset_level_band
    using Dates: Date, Day
    using Random: MersenneTwister

    rng = MersenneTwister(77)
    dates = [Date("2026-06-01") + Day(i) for i in 0:9]
    observed = [10.0 + i for i in 0:9]
    draws = [observed[k] .+ randn(rng, 200) for k in eachindex(dates)]

    fig = plot_onset_level_band(dates, observed, draws; title = "test")
    @test fig isa CairoMakie.Makie.Figure
    ## An empty series is a blank figure, not an error.
    @test plot_onset_level_band(Date[], Float64[], []; title = "empty") isa
        CairoMakie.Makie.Figure
    @test_throws ErrorException plot_onset_level_band(
        dates, observed, draws[1:3]; title = "mismatch"
    )
end

@testitem "_composition_predictive: allocates the observed total, wider with rho" begin
    using Statistics: mean, std
    using BVDOutbreakSize: _composition_predictive

    ## The predictive path must reproduce the composition's own generative
    ## step: the patches partition each vintage's observed total exactly, so
    ## their predicted shares sum to one draw by draw.
    nd = 400
    shares = [0.6 0.5; 0.3 0.35; 0.1 0.15]
    ms = [shares for _ in 1:nd]
    totals = [200, 300]
    preds = _composition_predictive(ms, fill(0.05, nd), totals, 2)
    @test length(preds) == 3
    for d in 1:nd, i in 1:2

        @test sum(preds[p][d][i] for p in 1:3) ≈ 1 atol = 1.0e-12
    end
    ## Centred on the expected share: the allocation is unbiased, so only the
    ## scatter around it is new.
    for p in 1:3, i in 1:2

        @test mean(preds[p][d][i] for d in 1:nd) ≈ shares[p, i] atol = 0.03
    end
    ## The overdispersion is what the band is for. A larger rho must scatter
    ## the predicted shares further, or the band says nothing the expected
    ## ribbon did not.
    tight = _composition_predictive(ms, fill(0.001, nd), totals, 2)
    loose = _composition_predictive(ms, fill(0.3, nd), totals, 2)
    for p in 1:2
        @test std(t[1] for t in loose[p]) > 2 * std(t[1] for t in tight[p])
    end
    ## A vintage with no observed cases has no split to predict.
    empty_v = _composition_predictive(ms, fill(0.05, nd), [0, 300], 2)
    @test all(isnan(t[1]) for t in empty_v[1])
    @test all(!isnan(t[2]) for t in empty_v[1])
end

@testitem "plot_province_composition_ppc: predictive band behind the expected" setup = [
    HeadlessMakie,
] begin
    using Dates: Date
    using BVDOutbreakSize: plot_province_composition_ppc

    nd = 120
    shares = [0.7 0.7 0.6 0.6; 0.2 0.2 0.25 0.25; 0.1 0.1 0.15 0.15]
    obs = [70 70 60 60; 20 20 25 25; 10 10 15 15]
    days = [7, 14, 21, 28]
    seeding = Date(2025, 8, 1)
    chn = (;
        province_shares = [shares for _ in 1:nd],
        province_composition_rho = fill(0.2, nd),
    )
    fig = plot_province_composition_ppc(
        chn; share_key = :province_shares,
        obs_increments = obs, days, seeding, n_patches = 3
    )
    @test fig isa CairoMakie.Makie.Figure
    axes = [x for x in fig.content if x isa CairoMakie.Makie.Axis]
    @test length(axes) == 3
    ## Three predictive bands behind three expected-share bands per panel.
    for ax in axes
        @test count(p -> p isa CairoMakie.Makie.Band, ax.scene.plots) == 6
    end
    ## The death composition reads its own overdispersion, so a chain
    ## carrying only the case one draws no band on the death panels.
    death = (;
        province_death_shares = [shares for _ in 1:nd],
        province_death_composition_rho = fill(0.2, nd),
    )
    dfig = plot_province_composition_ppc(
        death;
        share_key = :province_death_shares, obs_increments = obs, days,
        seeding, n_patches = 3
    )
    dax = first(x for x in dfig.content if x isa CairoMakie.Makie.Axis)
    @test count(p -> p isa CairoMakie.Makie.Band, dax.scene.plots) == 6
    ## A chain predating the deterministic still plots the expected share.
    plain = (; province_shares = [shares for _ in 1:nd])
    pfig = plot_province_composition_ppc(
        plain; share_key = :province_shares,
        obs_increments = obs, days, seeding, n_patches = 3
    )
    pax = first(x for x in pfig.content if x isa CairoMakie.Makie.Axis)
    @test count(p -> p isa CairoMakie.Makie.Band, pax.scene.plots) == 3
end

@testitem "plot_patch_summary: one panel per quantity, one interval per province" setup = [
    HeadlessMakie,
] begin
    using Random: MersenneTwister
    using BVDOutbreakSize: plot_patch_summary, PROVINCE_LABELS

    rng = MersenneTwister(11)
    nd = 200
    np = 3
    draws(centre) = [centre .+ 0.1 .* randn(rng, np) for _ in 1:nd]
    base = (;
        C_T_patch = draws([900.0, 300.0, 80.0]),
        R_T_patch = draws([1.2, 0.9, 0.7]),
        infections_T_patch = draws([40.0, 12.0, 3.0]),
        delta_patch = draws([0.2, -0.05, -0.15]),
    )
    fig = plot_patch_summary(base, np)
    @test fig isa CairoMakie.Makie.Figure
    axes = [x for x in fig.content if x isa CairoMakie.Makie.Axis]
    @test length(axes) == 4
    for ax in axes
        ## Three nested bars and a median dot per province.
        @test count(p -> p isa CairoMakie.Makie.Lines, ax.scene.plots) ==
            3 * np
        @test count(p -> p isa CairoMakie.Makie.Scatter, ax.scene.plots) == np
        ## The provinces share the axis, one named tick each, rather than
        ## being stacked on one position.
        positions, labels = ax.xticks[]
        @test positions == Float64.(1:np)
        @test labels == String.(PROVINCE_LABELS[1:np])
    end
    ## A reference rule only where the quantity has one: the reproduction
    ## number against one, the log-Rt deviation against zero.
    @test [
        count(p -> p isa CairoMakie.Makie.HLines, ax.scene.plots)
            for ax in axes
    ] == [0, 1, 0, 1]
    ## The optional quantities gain a panel each when the chain carries them,
    ## matching the seven `patch_summary_table` reports.
    full = (;
        base..., log_rt_contrast = draws([0.0, -0.3, -0.5]),
        region_drift_sd = draws([0.05, 0.04, 0.06]),
        province_ascertainment = draws([1.4, 0.8, 0.6]),
    )
    ffig = plot_patch_summary(full, np)
    @test count(x -> x isa CairoMakie.Makie.Axis, ffig.content) == 7
    ## A chain that is not from `bvd_joint` says so rather than failing deep
    ## inside the draw lookup.
    @test_throws ErrorException plot_patch_summary((; base.C_T_patch), np)
end

@testitem "reconstruct_patch_rt: the deviation knots unflatten per province" begin
    using Random: MersenneTwister
    import FlexiChains
    using BVDOutbreakSize: reconstruct_patch_rt, reconstruct_rt, knot_days,
        interpolate_knots, RT_INTERVENTION_RAMP

    rng = MersenneTwister(5)
    nd, n, np = 40, 60, 3
    walk_start = 12
    days = knot_days(n; week = 7, start = walk_start)
    nb = length(days)
    P = FlexiChains.Parameter
    ## Deviations that differ along both axes, so a row/column-major mix-up
    ## in the reshape moves values rather than permuting equal ones.
    dev(p, k) = 0.1 * p + 0.03 * k
    rt_keys = Dict(
        P(Symbol("rt_state.log_R0")) => reshape(fill(log(1.5), nd), nd, 1),
        P(Symbol("rt_state.sigma_rw")) => reshape(fill(0.05, nd), nd, 1),
        P(Symbol("rt_state.intervention_effect")) => reshape(
            fill(-0.3, nd), nd, 1
        ),
        P(Symbol("rt_state.z")) => reshape(
            [randn(rng, nb - 1) for _ in 1:nd], nd, 1
        )
    )
    chain(knots) = FlexiChains.FlexiChain{Symbol}(
        nd, 1,
        merge(rt_keys, Dict(P(:delta_knots) => reshape(knots, nd, 1)))
    )
    args = (;
        n, breakpoint = n - 11, n_patches = np, rt_start = walk_start,
        rt_walk_start = walk_start, week = 7, ramp = RT_INTERVENTION_RAMP,
    )

    national = reconstruct_rt(
        chain([zeros(np * nb) for _ in 1:nd]);
        n, breakpoint = n - 11, rt_start = walk_start,
        rt_walk_start = walk_start, week = 7, ramp = RT_INTERVENTION_RAMP
    )
    ## Zero deviations leave the national trajectory untouched in every
    ## panel, mask included: nothing rescales a province to the trend.
    flat = reconstruct_patch_rt(chain([zeros(np * nb) for _ in 1:nd]); args...)
    @test length(flat) == np
    for p in 1:np
        @test isequal(flat[p], national)
    end
    ## Days before the renewal start are masked in every province exactly as
    ## they are nationally, rather than filled with the walk base.
    @test all(ismissing, flat[1][:, 1:(walk_start - 1)])
    @test all(!ismissing, flat[1][:, walk_start:n])

    ## The matrix reaches the chain flattened column-major, so patch `p` at
    ## knot `k` sits at `(k - 1) * np + p`.
    knots = [vec([dev(p, k) for p in 1:np, k in 1:nb]) for _ in 1:nd]
    rt = reconstruct_patch_rt(chain(knots); args...)
    for p in 1:np
        daily = interpolate_knots([dev(p, k) for k in 1:nb], days, n)
        @test all(
            rt[p][i, d] ≈ national[i, d] * exp(daily[d])
                for i in 1:nd, d in walk_start:n
        )
    end
    ## Province 3 sits above province 1 throughout, since its deviation is
    ## larger at every knot. A transposed reshape would not preserve that.
    @test all(rt[3][i, d] > rt[1][i, d] for i in 1:nd, d in walk_start:n)

    ## A chain sampled with no spatial structure says so rather than
    ## silently repeating the national trajectory in every panel.
    nodev = FlexiChains.FlexiChain{Symbol}(nd, 1, rt_keys)
    @test_throws ErrorException reconstruct_patch_rt(nodev; args...)
    ## A knot vector that does not match `n_patches` by the knot grid is a
    ## mismatched `rt_walk_start` or patch count, not a usable deviation.
    short = [zeros(np * nb - 1) for _ in 1:nd]
    @test_throws ErrorException reconstruct_patch_rt(chain(short); args...)
end

@testitem "plot_rt_patches: one panel per province on a shared axis" setup = [
    HeadlessMakie,
] begin
    using Random: MersenneTwister
    using Dates: Date
    import FlexiChains
    using BVDOutbreakSize: plot_rt_patches, knot_days, PROVINCE_LABELS,
        RT_INTERVENTION_RAMP

    rng = MersenneTwister(23)
    nd, n, np = 40, 60, 3
    walk_start = 12
    nb = length(knot_days(n; week = 7, start = walk_start))
    P = FlexiChains.Parameter
    ## Province 3 runs well above province 1, so a per-panel autoscale would
    ## give the panels different y-limits.
    knots = [vec([0.4 * p for p in 1:np, _ in 1:nb]) for _ in 1:nd]
    chn = FlexiChains.FlexiChain{Symbol}(
        nd, 1,
        Dict(
            P(Symbol("rt_state.log_R0")) => reshape(fill(log(1.5), nd), nd, 1),
            P(Symbol("rt_state.sigma_rw")) => reshape(fill(0.05, nd), nd, 1),
            P(Symbol("rt_state.intervention_effect")) => reshape(
                fill(-0.3, nd), nd, 1
            ),
            P(Symbol("rt_state.z")) => reshape(
                [randn(rng, nb - 1) for _ in 1:nd], nd, 1
            ),
            P(:delta_knots) => reshape(knots, nd, 1)
        )
    )

    fig = plot_rt_patches(
        chn; n, breakpoint = n - 11,
        as_of_date = "2026-05-28", seeding = Date("2026-02-23"),
        n_patches = np, rt_start = walk_start, rt_walk_start = walk_start,
        ramp = RT_INTERVENTION_RAMP
    )
    @test fig isa CairoMakie.Makie.Figure
    axes = [x for x in fig.content if x isa CairoMakie.Makie.Axis]
    @test length(axes) == np
    @test [ax.title[] for ax in axes] == String.(PROVINCE_LABELS[1:np])
    for ax in axes
        ## The national reference behind the province: three ribbons each.
        @test count(p -> p isa CairoMakie.Makie.Band, ax.scene.plots) == 6
        ## No-growth line, plus the breakpoint, the end of the scale-up and
        ## the cut-off as vertical rules.
        @test count(p -> p isa CairoMakie.Makie.HLines, ax.scene.plots) == 1
        @test count(p -> p isa CairoMakie.Makie.VLines, ax.scene.plots) == 3
    end
    ## The provinces are compared rather than each rescaled to its own range,
    ## so every panel carries the same y-limits, from zero.
    ylims = [ax.limits[][2] for ax in axes]
    @test all(==(ylims[1]), ylims)
    @test ylims[1][1] == 0
    @test ylims[1][2] >= 2.5
    ## The panels do differ in what they draw, so a per-panel autoscale
    ## would have given them different limits.
    function panel_top(ax)
        lim = CairoMakie.Makie.data_limits(ax.scene)
        return lim.origin[2] + lim.widths[2]
    end
    @test panel_top(axes[3]) > panel_top(axes[1])

    ## The panels wrap onto a second row when the column count does not
    ## divide the provinces.
    wide = plot_rt_patches(
        chn; n, breakpoint = n - 11,
        as_of_date = "2026-05-28", seeding = Date("2026-02-23"),
        n_patches = np, rt_start = walk_start, rt_walk_start = walk_start,
        ramp = RT_INTERVENTION_RAMP, ncols = 2
    )
    @test count(x -> x isa CairoMakie.Makie.Axis, wide.content) == np
end

@testitem "_patch_daily: reads the column-major per-province matrix" begin
    import FlexiChains
    using BVDOutbreakSize: _patch_daily

    np, n, nd = 3, 5, 4
    ## Day `t` of patch `p` sits at `(t - 1) * np + p`, so a value encoding
    ## both indices catches a transposed read.
    flat = [Float64[100 * p + t for t in 1:n for p in 1:np] for _ in 1:nd]
    P = FlexiChains.Parameter
    chn = FlexiChains.FlexiChain{Symbol}(
        nd, 1,
        Dict(P(:infections_patch) => reshape(flat, nd, 1))
    )
    out = _patch_daily(chn, :infections_patch, np, n)
    @test length(out) == np
    for p in 1:np
        @test length(out[p]) == nd
        @test all(v == Float64[100 * p + t for t in 1:n] for v in out[p])
    end
    ## A chain without the deterministic, and one whose matrix does not match
    ## the patch count by the day count, both name what is wrong.
    @test_throws ErrorException _patch_daily(chn, :importation_patch, np, n)
    @test_throws ErrorException _patch_daily(chn, :infections_patch, np, n + 1)
end

@testitem "plot_infections_patches: the cumulative row sums the daily one" setup = [
    HeadlessMakie,
] begin
    using Dates: Date
    import FlexiChains
    using BVDOutbreakSize: plot_infections_patches, PROVINCE_LABELS
    Mk = CairoMakie.Makie

    np, n, nd = 3, 20, 30
    ## A constant daily rate per province, so the cumulative panel must end
    ## at `n` times it and the daily panel must stay flat.
    rate(p) = p == 1 ? 1.0 : 0.5
    flat = [Float64[rate(p) for t in 1:n for p in 1:np] for _ in 1:nd]
    P = FlexiChains.Parameter
    chn = FlexiChains.FlexiChain{Symbol}(
        nd, 1,
        Dict(P(:infections_patch) => reshape(flat, nd, 1))
    )
    fig = plot_infections_patches(
        chn; n, seeding = Date("2026-02-23"),
        n_patches = np
    )
    @test fig isa Mk.Figure
    axes = [x for x in fig.content if x isa Mk.Axis]
    ## A daily and a cumulative panel per province, the province named once.
    @test length(axes) == 2 * np
    @test [ax.ylabel[] for ax in axes] ==
        repeat(["Daily infections", "Cumulative infections"], np)
    @test [ax.title[] for ax in axes] ==
        vcat([[String(PROVINCE_LABELS[p]), ""] for p in 1:np]...)
    for (i, ax) in enumerate(axes)
        @test count(p -> p isa Mk.Band, ax.scene.plots) == 3
        p = cld(i, 2)
        ylo = Mk.data_limits(ax.scene).origin[2]
        yhi = ylo + Mk.data_limits(ax.scene).widths[2]
        if isodd(i)
            ## Daily: flat at the province's own rate.
            @test ylo ≈ rate(p)
            @test yhi ≈ rate(p)
        else
            ## Cumulative: the running sum of the panel above it.
            @test ylo ≈ rate(p)
            @test yhi ≈ n * rate(p)
        end
    end
end

@testitem "plot_imports_patches: one panel per province" setup = [
    HeadlessMakie,
] begin
    using Dates: Date
    import FlexiChains
    using BVDOutbreakSize: plot_imports_patches, PROVINCE_LABELS
    Mk = CairoMakie.Makie

    np, n, nd = 3, 20, 30
    ## Arrivals that differ by province and rise over the window, so a panel
    ## reading the wrong slice of the flattened matrix shows it.
    flat = [Float64[p * t for t in 1:n for p in 1:np] for _ in 1:nd]
    P = FlexiChains.Parameter
    chn = FlexiChains.FlexiChain{Symbol}(
        nd, 1,
        Dict(P(:importation_patch) => reshape(flat, nd, 1))
    )
    fig = plot_imports_patches(
        chn; n, seeding = Date("2026-02-23"),
        n_patches = np
    )
    @test fig isa Mk.Figure
    axes = [x for x in fig.content if x isa Mk.Axis]
    @test length(axes) == np
    @test [ax.title[] for ax in axes] == String.(PROVINCE_LABELS[1:np])
    ## Each panel carries its own y-axis, reaching that province's own peak.
    for (p, ax) in enumerate(axes)
        @test count(q -> q isa Mk.Band, ax.scene.plots) == 3
        lim = Mk.data_limits(ax.scene)
        @test lim.origin[2] ≈ p * 1.0
        @test lim.origin[2] + lim.widths[2] ≈ p * n
    end
end

@testitem "plot_forecast_flows: a panel per flow stream carried" setup = [
    HeadlessMakie,
] begin
    using DataFrames: DataFrame
    using Random: MersenneTwister
    using BVDOutbreakSize: plot_forecast_flows
    Mk = CairoMakie.Makie

    rng = MersenneTwister(3)
    draws = 400
    full = DataFrame(
        admissions_fc = rand(rng, 0:20, draws),
        incare_deaths_fc = rand(rng, 0:5, draws),
        ruleouts_fc = rand(rng, 0:40, draws)
    )
    fig = plot_forecast_flows(full)
    axes = [x for x in fig.content if x isa Mk.Axis]
    @test length(axes) == 3
    @test [ax.xlabel[] for ax in axes] == [
        "New isolation admissions (DRC)",
        "New in-care deaths (DRC)", "New rule-outs (DRC)",
    ]
    for ax in axes
        ## A histogram of the predictive with its 90% interval shaded.
        @test count(p -> p isa Mk.Hist, ax.scene.plots) == 1
        @test count(p -> p isa Mk.VSpan, ax.scene.plots) == 1
    end
    ## Only the streams the forecast carries get a panel.
    part = plot_forecast_flows(full[!, [:admissions_fc, :ruleouts_fc]])
    @test [ax.xlabel[] for ax in part.content if ax isa Mk.Axis] ==
        ["New isolation admissions (DRC)", "New rule-outs (DRC)"]
    ## A forecast without the flow streams draws nothing rather than an
    ## empty grid of axes.
    none = plot_forecast_flows(DataFrame(cases_fc = [1, 2, 3]))
    @test none isa Mk.Figure
    @test count(x -> x isa Mk.Axis, none.content) == 0
end

@testitem "plot_stream_calibration: streams against their nominal level" setup = [
    HeadlessMakie,
] begin
    using Random: MersenneTwister
    using BVDOutbreakSize: plot_stream_calibration, stream_calibration
    Mk = CairoMakie.Makie

    rng = MersenneTwister(3)
    observed = [10, 20, 30, 40]
    ## One calibrated stream and one that predicts far too high, so the bias
    ## panel has a sign to place.
    panels = [
        (;
            title = "cases", observed = observed,
            replicates = [10 .+ rand(rng, -5:5, 4) for _ in 1:500],
        ),
        (;
            title = "deaths", observed = observed,
            replicates = [fill(40.0, 4) for _ in 1:200],
        ),
    ]
    tbl = stream_calibration(panels)
    fig = plot_stream_calibration(tbl)
    @test fig isa Mk.Figure
    axes = [x for x in fig.content if x isa Mk.Axis]
    @test length(axes) == 2
    ## Streams share the row labels, the first stream at the top, and only
    ## the left panel names them.
    pos, labs = axes[1].yticks[]
    @test pos == [2, 1]
    @test labs == ["cases", "deaths"]
    @test axes[2].yticks[][2] == ["", ""]
    ## Coverage is a fraction, so the left panel is pinned to [0, 1] and
    ## carries the two nominal reference lines.
    @test axes[1].limits[][1] == (0, 1)
    @test count(p -> p isa Mk.VLines, axes[1].scene.plots) == 2
    @test count(p -> p isa Mk.Scatter, axes[1].scene.plots) == 2
    ## The bias panel is symmetric about zero and clears the largest bias.
    lo, hi = axes[2].limits[][1]
    @test lo ≈ -hi
    @test hi >= maximum(abs, tbl[!, "Bias"])
    @test count(p -> p isa Mk.VLines, axes[2].scene.plots) == 1
end

@testitem "plot_province_forecast draws one panel per stream" setup = [
    HeadlessMakie,
] begin
    using DataFrames: DataFrame
    using CairoMakie: Makie as Mk
    using BVDOutbreakSize: plot_province_forecast, PROVINCE_LABELS

    ## Deterministic shares against a spread national forecast, so each
    ## province's interval is a known fraction of the national one.
    nd = 200
    shares = [0.8 0.75; 0.15 0.2; 0.05 0.05]
    chn = (;
        province_shares = [shares for _ in 1:nd],
        province_death_shares = [shares for _ in 1:nd],
    )
    v = collect(range(50.0, 150.0; length = nd))
    fc = DataFrame(confirmed_new = v, confirmed_deaths_new = v ./ 5)
    fig = plot_province_forecast(chn, fc; n_patches = 3)
    @test fig isa Mk.Figure

    axes = [x for x in fig.content if x isa Mk.Axis]
    @test length(axes) == 2
    @test axes[1].title[] == "New confirmed cases by T+7"
    @test axes[2].title[] == "New confirmed deaths by T+7"
    ## Provinces are the shared category axis, in patch order.
    pos, labs = axes[1].xticks[]
    @test collect(pos) == [1, 2, 3]
    @test labs == PROVINCE_LABELS[1:3]
    ## Three nested interval bars and a median dot per province.
    for ax in axes
        @test count(p -> p isa Mk.Lines, ax.scene.plots) == 9
        @test count(p -> p isa Mk.Scatter, ax.scene.plots) == 3
        ## Counts are non-negative, so the panel is read against zero.
        @test ax.limits[][2][1] == 0
    end

    ## A forecast carrying one stream draws that panel alone, and one
    ## carrying neither returns an empty figure rather than erroring.
    one = plot_province_forecast(
        chn,
        DataFrame(confirmed_new = v); n_patches = 3
    )
    @test length([x for x in one.content if x isa Mk.Axis]) == 1
    none = plot_province_forecast(
        chn, DataFrame(cases_new = v);
        n_patches = 3
    )
    @test isempty([x for x in none.content if x isa Mk.Axis])
end

@testitem "plot_province_forecast_detail draws one province's streams" setup = [
    HeadlessMakie,
] begin
    using DataFrames: DataFrame
    using CairoMakie: Makie as Mk
    using Statistics: quantile
    using BVDOutbreakSize: plot_province_forecast_detail, plot_forecast,
        PROVINCE_LABELS

    nd = 200
    shares = [0.8 0.75; 0.15 0.2; 0.05 0.05]
    chn = (;
        province_shares = [shares for _ in 1:nd],
        province_death_shares = [shares for _ in 1:nd],
    )
    v = collect(range(50.0, 150.0; length = nd))
    fc = DataFrame(confirmed_new = v, confirmed_deaths_new = v ./ 5)

    fig = plot_province_forecast_detail(chn, fc; province = 2, n_patches = 3)
    @test fig isa Mk.Figure
    axes = [x for x in fig.content if x isa Mk.Axis]
    @test length(axes) == 2
    ## Each panel names the stream and the province it is a split for.
    @test axes[1].xlabel[] == "New confirmed cases ($(PROVINCE_LABELS[2]))"
    @test axes[2].xlabel[] == "New confirmed deaths ($(PROVINCE_LABELS[2]))"
    ## The histogram is that province's share of the national draws, so its
    ## 90% band is 0.2 times the national band for the deaths.
    band = only(p for p in axes[2].scene.plots if p isa Mk.VSpan)
    @test band[1][][1] ≈ 0.2 * quantile(v ./ 5, 0.05)
    ## Each panel takes the colour of the matching national panel.
    national = [
        x for x in plot_forecast(fc).content if x isa Mk.Axis
    ]
    hist_colour(ax) = only(p for p in ax.scene.plots if p isa Mk.Hist).color[]
    @test hist_colour.(axes) == hist_colour.(national)
    ## No observed week, no reference rule.
    @test all(ax -> !any(p -> p isa Mk.VLines, ax.scene.plots), axes)

    ## An observed week is drawn as a rule, and the axis widens to hold an
    ## observation beyond the forecast's upper tail.
    obs_fig = plot_province_forecast_detail(
        chn, fc; province = 2, n_patches = 3,
        observed = (; confirmed_new = 500, confirmed_deaths_new = 3)
    )
    obs_axes = [x for x in obs_fig.content if x isa Mk.Axis]
    for ax in obs_axes
        @test count(p -> p isa Mk.VLines, ax.scene.plots) == 1
    end
    @test obs_axes[1].limits[][1][2] >= 500

    ## A forecast without the deaths column draws the cases panel alone.
    one = plot_province_forecast_detail(
        chn, DataFrame(confirmed_new = v); province = 1, n_patches = 3
    )
    @test length([x for x in one.content if x isa Mk.Axis]) == 1

    ## A province outside the patches is an error, not an empty figure.
    @test_throws ArgumentError plot_province_forecast_detail(
        chn, fc; province = 4, n_patches = 3
    )
end

@testitem "plot_evolution_by_group clamps and marks past ymax" setup = [
    HeadlessMakie,
] begin
    using BVDOutbreakSize: plot_evolution_by_group
    ## One release whose 90% upper runs far past the crop and one that does
    ## not, the case a fixed reproduction-number axis exists for.
    groups = [
        "cases" => [
            ("2026-08-01", 1.2, 1.0, 1.4, 0.9, 1.6, 0.6, 9.0),
            ("2026-08-08", 1.1, 0.9, 1.3, 0.8, 1.5, 0.5, 2.0),
        ],
    ]
    fig = plot_evolution_by_group(groups; refline = 1.0, ymax = 3.0)
    ax = only(x for x in fig.content if x isa CairoMakie.Makie.Axis)
    @test ax.limits[][2] == (0, 3.0)
    ## The overflow marker sits just inside the crop, so the clipping does
    ## not cut it in half.
    ## Makie converts a marker symbol into a path before it reaches the
    ## plot, so the triangle is matched against the converted path.
    tri = CairoMakie.Makie.convert_attribute(
        :utriangle,
        CairoMakie.Makie.key"marker"(), CairoMakie.Makie.key"scatter"()
    )
    markers = [
        p for p in ax.scene.plots
            if p isa CairoMakie.Makie.Scatter && p.marker[] == tri
    ]
    @test length(markers) == 1
    @test all(pt -> pt[2] < 3.0, only(markers)[1][])
    ## Without a crop the axis still sizes itself to the widest interval.
    free = plot_evolution_by_group(groups; refline = 1.0)
    free_ax = only(x for x in free.content if x isa CairoMakie.Makie.Axis)
    @test free_ax.limits[][2][2] > 9.0
end

@testitem "plot_stream_trajectories crops to ymax" setup = [HeadlessMakie] begin
    using BVDOutbreakSize: plot_stream_trajectories
    using Dates: Date
    n = 20
    ## One stream on the scale of the source population and one on the scale
    ## of the outbreak, the pair that flattens a free axis.
    streams = [
        (;
            label = "wide", trajs = [cumsum(fill(1.0e6, n)) for _ in 1:20],
            last_day = 10, colour = :seagreen,
        ),
        (;
            label = "narrow", trajs = [cumsum(fill(10.0, n)) for _ in 1:20],
            last_day = 15, colour = :steelblue,
        ),
    ]
    fig = plot_stream_trajectories(
        streams; n = n, seeding = Date(2026, 3, 1), ymax = 1.0e4
    )
    ax = only(x for x in fig.content if x isa CairoMakie.Makie.Axis)
    @test ax.limits[][2] == (0, 1.0e4)
    ## Makie converts a marker symbol into a path before it reaches the
    ## plot, so the triangle is matched against the converted path.
    tri = CairoMakie.Makie.convert_attribute(
        :utriangle,
        CairoMakie.Makie.key"marker"(), CairoMakie.Makie.key"scatter"()
    )
    markers = [
        p for p in ax.scene.plots
            if p isa CairoMakie.Makie.Scatter && p.marker[] == tri
    ]
    ## The wide stream leaves the axis and is marked; the narrow one never
    ## reaches it.
    @test length(markers) == 1
    free = plot_stream_trajectories(streams; n = n, seeding = Date(2026, 3, 1))
    free_ax = only(x for x in free.content if x isa CairoMakie.Makie.Axis)
    @test free_ax.limits[][2][2] > 1.0e6
end

@testitem "plot_forecast_crps_by_horizon empty and filled" setup = [
    HeadlessMakie,
] begin
    using BVDOutbreakSize: plot_forecast_crps_by_horizon, FROZEN_FIT
    using DataFrames: DataFrame
    schema = (;
        stream = String[], horizon = Int[], fit = String[],
        dispersion = Float64[], overprediction = Float64[],
        underprediction = Float64[],
    )
    @test plot_forecast_crps_by_horizon(DataFrame(schema)) isa
        CairoMakie.Makie.Figure
    rows = NamedTuple[]
    for s in ("confirmed cases", "confirmed deaths"), h in (7, 14),
            f in (FROZEN_FIT, "confirmed_deaths")
        push!(
            rows,
            (;
                stream = s, horizon = h, fit = f, dispersion = 40.0,
                overprediction = 30.0, underprediction = 30.0 + h,
            )
        )
    end
    fig = plot_forecast_crps_by_horizon(DataFrame(rows))
    axes = [x for x in fig.content if x isa CairoMakie.Makie.Axis]
    @test length(axes) == 2
end

@testitem "plot_forecast_skill_by_cutoff empty and filled" setup = [
    HeadlessMakie,
] begin
    using BVDOutbreakSize: plot_forecast_skill_by_cutoff, FROZEN_FIT
    using DataFrames: DataFrame
    using Dates: Date
    schema = (;
        stream = String[], made_date = Date[], fit = String[],
        rel_to_baseline = Float64[],
    )
    @test plot_forecast_skill_by_cutoff(DataFrame(schema)) isa
        CairoMakie.Makie.Figure
    rows = NamedTuple[]
    for d in (Date(2026, 7, 16), Date(2026, 8, 4)),
            f in (FROZEN_FIT, "confirmed")
        push!(
            rows,
            (;
                stream = "confirmed cases", made_date = d, fit = f,
                rel_to_baseline = 1.5,
            )
        )
    end
    fig = plot_forecast_skill_by_cutoff(DataFrame(rows))
    ax = only(x for x in fig.content if x isa CairoMakie.Makie.Axis)
    ## One slot per made date, whichever fits carry it.
    @test ax.limits[][1] == (0.5, 2.5)
end
