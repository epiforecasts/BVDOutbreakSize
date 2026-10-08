## The summary scatter of weekly confirmed cases against the reproduction
## number, built from synthetic trajectories and counts.

@testsnippet CasesRt begin
    using Dates: Date, Day
    using DataFrames: DataFrame, Not, nrow
    ## Two areas on a 28-day grid with weekly vintages. Draw `i` of area 1
    ## on day `d` is `(i - 1) / 100 + d / 10`, so its median is
    ## `0.5 + d / 10` and its 90% interval is 0.45 either side. Area 2 is
    ## area 1 plus one, with no draws on day 14.
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

@testitem "cases_rt_table pairs weekly cases with Rt" setup = [CasesRt] begin
    using BVDOutbreakSize: cases_rt_table
    t = cases_rt_table(
        rt, days, increments; cutoff, n, weeks = 3, areas = ["A", "B"],
        forecast
    )
    obs = t[t.kind .== "observed", :]
    ## Weeks end on days 14, 21 and 28; the week to day 7 starts before the
    ## first vintage and area 2 has no draws on day 14.
    @test obs.area == ["A", "A", "A", "B", "B"]
    @test obs.patch == [1, 1, 1, 2, 2]
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
    ## Each area's rows run in date order, the forecast last.
    @test t.kind[t.area .== "A"] == [fill("observed", 3); "forecast"]
    ## Areas keep their given order, not alphabetical, and take their
    ## patch from `patch`.
    z = cases_rt_table(
        rt, days, increments; cutoff, n, areas = ["Z", "Y"],
        patch = [3, 3]
    )
    @test unique(z.area) == ["Z", "Y"]
    @test all(==(3), z.patch)
    @test nrow(z) == 4
    ## The weeks end at the last vintage, not the cut-off.
    late = cases_rt_table(
        [a[:, 1:27], b[:, 1:27]], days[1:3], increments[:, 1:3];
        cutoff = cutoff - Day(1), n = 27, weeks = 1, areas = ["A", "B"]
    )
    @test late.date == fill(cutoff - Day(7), 2)
    @test late.cases == [4, 0]
    @test_throws ArgumentError cases_rt_table(
        rt, days, increments; cutoff, n, areas = ["A", "B"],
        forecast = forecast[:, Not(:rt_forecast)]
    )
    @test_throws ArgumentError cases_rt_table(
        rt, days, increments[1:1, :]; cutoff, n, areas = ["A", "B"]
    )
end

@testitem "plot_cases_rt draws intervals, the week before and forecasts" setup = [
    CasesRt,
] begin
    using CairoMakie
    using CairoMakie: Makie as Mk
    using BVDOutbreakSize: cases_rt_table, plot_cases_rt
    CairoMakie.activate!(type = "png")
    ## Makie converts a marker or line style symbol before it reaches the
    ## plot, so the plots are matched against the converted values.
    conv(v, k, p) = Mk.convert_attribute(v, Mk.Key{k}(), Mk.Key{p}())
    triangle = conv(:utriangle, :marker, :scatter)
    dash = conv(:dash, :linestyle, :lines)
    labels = ["A", "B"]
    t = cases_rt_table(
        rt, days, increments; cutoff, n, weeks = 2, areas = labels, forecast
    )
    fig = plot_cases_rt(t; patch_labels = labels)
    @test fig isa Mk.Figure
    ax = only(x for x in fig.content if x isa Mk.Axis)
    @test ax.xscale[] === Mk.pseudolog10
    plots = ax.scene.plots
    ## One line at R = 1 and one at the median of the latest weekly counts.
    @test only(p for p in plots if p isa Mk.HLines)[1][] == [1.0]
    @test only(p for p in plots if p isa Mk.VLines)[1][] == [4.0]
    ## A faded triangle for the week before, joined to the point, and a
    ## dashed link to the forecast, for each area.
    @test count(p -> p isa Mk.Scatter && p.marker[] == triangle, plots) == 2
    lines = [p for p in plots if p isa Mk.Lines]
    @test count(p -> p.linestyle[] == dash, lines) == 2
    @test length(lines) == 4
    ## Each point carries crossed bars spanning its 90% intervals: the
    ## reproduction number for every row drawn, and the cases only for the
    ## forecast, since an observed count has no interval.
    bars = [p for p in plots if p isa Mk.Rangebars]
    vertical = [p for p in bars if p.direction[] == :y]
    horizontal = [p for p in bars if p.direction[] == :x]
    @test length(vertical) == 6
    @test length(horizontal) == 2
    latest = t[(t.kind .== "observed") .& (t.date .== cutoff), :]
    fc = t[t.kind .== "forecast", :]
    ## Makie stores each bar as a (position, low, high) triple.
    span(p) = (only(p[1][])[2], only(p[1][])[3])
    spans = sort(span.(vertical))
    for r in eachrow(vcat(latest, fc))
        @test any(s -> s[1] ≈ r.rt_lower && s[2] ≈ r.rt_upper, spans)
    end
    @test sort(span.(horizontal)) ==
        sort(collect(zip(fc.cases_lower, fc.cases_upper)))
    ## Point size is fixed, so it encodes nothing.
    points = [
        p for p in plots
            if p isa Mk.Scatter && p.marker[] != triangle &&
            p.strokecolor[] == Mk.to_color(:white)
    ]
    @test length(points) == 2
    @test allequal(p.markersize[] for p in points)
    @test sort([only(p.text[]) for p in plots if p isa Mk.Text]) == labels
    leg = only(x for x in fig.content if x isa Mk.Legend)
    entries = [e.label[] for e in leg.entrygroups[][1][2]]
    @test "Week ahead (forecast)" in entries
    @test "Week before" in entries
    @test "90% credible interval" in entries
    ## Without a forecast its entry is gone; `top` keeps the areas with the
    ## most cases over the weeks shown.
    obs_only = t[t.kind .== "observed", :]
    fig2 = plot_cases_rt(
        obs_only; patch_labels = labels, top = 1, unit = "health zone"
    )
    ax2 = only(x for x in fig2.content if x isa Mk.Axis)
    @test [only(p.text[]) for p in ax2.scene.plots if p isa Mk.Text] == ["A"]
    leg2 = only(x for x in fig2.content if x isa Mk.Legend)
    @test !(
        "Week ahead (forecast)" in [e.label[] for e in leg2.entrygroups[][1][2]]
    )
    ## An empty table says so rather than drawing empty axes.
    empty_fig = plot_cases_rt(obs_only[1:0, :]; patch_labels = labels)
    @test !any(x -> x isa Mk.Axis, empty_fig.content)
end
