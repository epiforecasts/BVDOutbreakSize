## The province summary scatter of weekly confirmed cases against the
## reproduction number, built from synthetic trajectories and counts.

@testsnippet ProvinceCasesRt begin
    using Dates: Date, Day
    using DataFrames: DataFrame, Not
    ## Two patches on a 28-day grid with weekly vintages. Draw `i` of
    ## patch 1 on day `d` is `(i - 1) / 100 + d / 10`, so its median is
    ## `0.5 + d / 10` and its 90% interval is 0.45 either side. Patch 2 is
    ## patch 1 plus one, with no draws on day 14.
    n = 28
    cutoff = Date(2026, 9, 28)
    days = [7, 14, 21, 28]
    increments = [5 3 4 6; 0 1 0 2]
    a = Matrix{Union{Missing, Float64}}(
        [(i - 1) / 100 + d / 10 for i in 1:101, d in 1:n]
    )
    b = Matrix{Union{Missing, Float64}}(a .+ 1)
    b[:, 14] .= missing
    rt = [a, b]
    ## The forecast's confirmed cases run 0 to 100 and its reproduction
    ## number 0 to 1, plus the patch index.
    forecast = DataFrame(
        patch = repeat(1:2; inner = 101), draw = repeat(1:101; outer = 2),
        confirmed_new = repeat(collect(0.0:100.0); outer = 2),
        rt_forecast = vcat(collect(0.0:0.01:1.0), collect(1.0:0.01:2.0))
    )
end

@testitem "province_cases_rt_table pairs weekly cases with Rt" setup = [
    ProvinceCasesRt,
] begin
    using BVDOutbreakSize: province_cases_rt_table
    t = province_cases_rt_table(
        rt, days, increments; cutoff, n, weeks = 3,
        patch_labels = ["A", "B"], forecast
    )
    obs = t[t.kind .== "observed", :]
    ## Weeks end on days 14, 21 and 28; the week to day 7 starts before the
    ## first vintage and patch 2 has no draws on day 14.
    @test obs.province == ["A", "A", "A", "B", "B"]
    @test obs.date == cutoff .- Day.([14, 7, 0, 7, 0])
    @test obs.cases == [3, 4, 6, 0, 2]
    @test obs.cases_lower == obs.cases == obs.cases_upper
    @test obs.rt_median ≈ [1.9, 2.6, 3.3, 3.6, 4.3]
    @test obs.rt_lower ≈ obs.rt_median .- 0.45
    @test obs.rt_upper ≈ obs.rt_median .+ 0.45
    fc = t[t.kind .== "forecast", :]
    @test fc.patch == [1, 2]
    @test fc.date == [cutoff + Day(7), cutoff + Day(7)]
    @test fc.cases ≈ [50, 50]
    @test fc.cases_lower ≈ [5, 5]
    @test fc.cases_upper ≈ [95, 95]
    @test fc.rt_median ≈ [0.5, 1.5]
    @test fc.rt_lower ≈ [0.05, 1.05]
    ## Each patch's rows run in date order, the forecast last.
    @test t.kind[t.patch .== 1] == [fill("observed", 3); "forecast"]
    ## The weeks end at the last vintage, not the cut-off.
    late = province_cases_rt_table(
        [a[:, 1:27], b[:, 1:27]], days[1:3], increments[:, 1:3];
        cutoff = cutoff - Day(1), n = 27, weeks = 1,
        patch_labels = ["A", "B"]
    )
    @test late.date == fill(cutoff - Day(7), 2)
    @test late.cases == [4, 0]
    @test_throws ArgumentError province_cases_rt_table(
        rt, days, increments; cutoff, n, patch_labels = ["A", "B"],
        forecast = forecast[:, Not(:rt_forecast)]
    )
    @test_throws ArgumentError province_cases_rt_table(
        rt, days, increments[1:1, :]; cutoff, n, patch_labels = ["A", "B"]
    )
end

@testitem "plot_province_cases_rt draws trajectories and forecasts" setup = [
    ProvinceCasesRt,
] begin
    using CairoMakie
    using CairoMakie: Makie as Mk
    using BVDOutbreakSize: province_cases_rt_table, plot_province_cases_rt
    CairoMakie.activate!(type = "png")
    ## Makie converts a marker or line style symbol before it reaches the
    ## plot, so the plots are matched against the converted values.
    conv(v, k, p) = Mk.convert_attribute(v, Mk.Key{k}(), Mk.Key{p}())
    diamond = conv(:diamond, :marker, :scatter)
    dot = conv(:dot, :linestyle, :lines)
    dash = conv(:dash, :linestyle, :lines)
    labels = ["A", "B"]
    t = province_cases_rt_table(
        rt, days, increments; cutoff, n, weeks = 3,
        patch_labels = labels, forecast
    )
    fig = plot_province_cases_rt(t; patch_labels = labels)
    @test fig isa Mk.Figure
    ax = only(x for x in fig.content if x isa Mk.Axis)
    @test ax.xscale[] === Mk.pseudolog10
    plots = ax.scene.plots
    ## One line at R = 1 and one at the median of the latest weekly counts.
    @test only(p for p in plots if p isa Mk.HLines)[1][] == [1.0]
    @test only(p for p in plots if p isa Mk.VLines)[1][] == [4.0]
    ## A dotted trajectory and a dashed forecast link per province.
    lines = [p for p in plots if p isa Mk.Lines]
    @test count(p -> p.linestyle[] == dot, lines) == 2
    @test count(p -> p.linestyle[] == dash, lines) == 2
    ## Rt bars on each latest week, Rt and case bars on each forecast.
    @test count(p -> p isa Mk.Rangebars, plots) == 6
    @test count(p -> p isa Mk.Scatter && p.marker[] == diamond, plots) == 2
    @test sort([only(p.text[]) for p in plots if p isa Mk.Text]) == labels
    leg = only(x for x in fig.content if x isa Mk.Legend)
    entries = [e.label[] for e in leg.entrygroups[][1][2]]
    @test "Forecast week" in entries
    ## Without a forecast the diamonds and their legend entry are gone, and
    ## a stated reference moves the vertical line.
    obs_only = t[t.kind .== "observed", :]
    fig2 = plot_province_cases_rt(
        obs_only; patch_labels = labels, case_reference = 10
    )
    ax2 = only(x for x in fig2.content if x isa Mk.Axis)
    @test only(p for p in ax2.scene.plots if p isa Mk.VLines)[1][] == [10.0]
    @test !any(
        p -> p isa Mk.Scatter && p.marker[] == diamond, ax2.scene.plots
    )
    leg2 = only(x for x in fig2.content if x isa Mk.Legend)
    @test !("Forecast week" in [e.label[] for e in leg2.entrygroups[][1][2]])
    @test any(
        x -> x isa Mk.Label && occursin("10 confirmed cases", x.text[]),
        fig2.content
    )
    ## An empty table says so rather than drawing empty axes.
    empty_fig = plot_province_cases_rt(obs_only[1:0, :]; patch_labels = labels)
    @test !any(x -> x isa Mk.Axis, empty_fig.content)
end
