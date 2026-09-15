## Smoke tests for the health-zone figures. Each builds its figure from
## small synthetic inputs and checks the layout it promises, without
## comparing pixels. The map tests read a three-zone geojson written to a
## temporary file, so they do not depend on the packaged polygons.

@testsnippet ZoneGeojson begin
    using CairoMakie
    CairoMakie.activate!(type = "png")

    ## Two provinces: P holds a unit square, a square to its right whose
    ## shared edge carries an extra vertex (the T-junction the dissolve
    ## must drop), and a two-part zone with a hole; Q holds one square
    ## below them. The feature without a `zone` property keys by its name.
    function write_zone_geojson()
        sq(x0, y0, x1, y1) = [[x0, y0], [x1, y0], [x1, y1], [x0, y1], [x0, y0]]
        ft(props, geom) = """{"type": "Feature", "properties": $props,
            "geometry": $geom}"""
        poly(rings) = """{"type": "Polygon", "coordinates": $rings}"""
        a = ft("""{"zone": "alpha", "label": "Alpha", "province": "P"}""",
            poly([sq(0, 0, 1, 1)]))
        b = ft("""{"zone": "", "label": "Beta Two", "province": "P"}""",
            poly([[[1, 0], [2, 0], [2, 1], [1, 1], [1, 0.5], [1, 0]]]))
        c = ft("""{"zone": "gamma", "label": "Gamma", "province": "P"}""",
            """{"type": "MultiPolygon", "coordinates": [
                [$(sq(0, 1, 1, 2)), $(sq(0.3, 1.3, 0.6, 1.6))],
                [$(sq(3, 3, 3.5, 3.5))]]}""")
        d = ft("""{"zone": "delta", "label": "Delta", "province": "Q"}""",
            poly([sq(0, -1, 1, 0)]))
        e = ft("""{"zone": "omega", "label": "Omega", "province": "Z"}""",
            poly([sq(5, 5, 6, 6)]))
        path = tempname() * ".geojson"
        fts = join((a, b, c, d, e), ", ")
        write(path, """{"type": "FeatureCollection", "features": [$fts]}""")
        return path
    end
end

@testitem "zone_key folds a zone name to its snake_case key" begin
    using BVDOutbreakSize: zone_key
    @test zone_key("Bili (Bas-Uele)") == "bili_bas_uele"
    @test zone_key("Miti Murhesa") == "miti_murhesa"
    @test zone_key("Bas-Uélé") == "bas_uele"
    @test zone_key("  Aru ") == "aru"
end

@testitem "load_health_zones_geojson keys, filters and dissolves" setup=[
    ZoneGeojson
] begin
    using BVDOutbreakSize: load_health_zones_geojson, _province_outlines,
                           _zone_table
    path = write_zone_geojson()
    zones = load_health_zones_geojson(path; provinces = ["P", "Q"])
    ## One record per zone. The province outside the filter is dropped;
    ## the empty `zone` property falls back to the folded label.
    @test [z.zone for z in zones] == ["alpha", "beta_two", "gamma", "delta"]
    @test zones[2].label == "Beta Two"
    @test [z.province for z in zones] == ["P", "P", "P", "Q"]
    @test [length(z.polygons) for z in zones] == [1, 1, 2, 1]
    ## The hole is kept, and the centroid sits on the largest part.
    @test length(zones[3].polygons[1].interiors) == 1
    @test zones[3].centroid ≈ CairoMakie.Point2f(0.5, 1.5)
    @test zones[1].centroid ≈ CairoMakie.Point2f(0.5, 0.5)
    geo = _zone_table(zones)
    @test geo.key == ["alpha", "beta_two", "gamma", "delta"]
    ## Dissolving province P: the shared edge between alpha and beta drops
    ## although beta splits it at an extra vertex, the alpha-gamma edge
    ## drops, and the hole ring survives as part of the outline. Province
    ## Q keeps all four sides, since delta borders alpha across a province
    ## line. Each segment is two points.
    seg = _province_outlines(geo)
    @test iseven(length(seg))
    ## P: alpha 2 outer sides + beta 3 sides + gamma main 3 sides + hole 4
    ## + gamma island 4; Q: 4.
    @test length(seg) ÷ 2 == 2 + 3 + 3 + 4 + 4 + 4
    @test_throws ErrorException load_health_zones_geojson(path;
        provinces = ["none"])
end

@testitem "plot_zone_map colours the matched zones and greys the rest" setup=[
    ZoneGeojson
] begin
    using CairoMakie: Makie as Mk
    using BVDOutbreakSize: plot_zone_map
    path = write_zone_geojson()
    fig = plot_zone_map([1.6, 0.7], ["alpha", "gamma"]; geojson = path,
        provinces = ["P", "Q"], diverging_at = 1.0, lower = [1.2, 0.4],
        upper = [2.1, 1.1],
        title = "R", colorbar_label = "R", label_top = 1)
    @test fig isa Mk.Figure
    axes = [x for x in fig.content if x isa Mk.Axis]
    @test length(axes) == 1
    @test count(x -> x isa Mk.Colorbar, fig.content) == 1
    polys = [p for p in axes[1].scene.plots if p isa Mk.Poly]
    ## Grey layer for the unmatched zones, the coloured layer, and the
    ## wash over the one zone whose interval straddles the centre.
    @test length(polys) == 3
    @test length(polys[1][1][]) == 2
    @test length(polys[2][1][]) == 3
    @test length(polys[3][1][]) == 2
    @test count(p -> p isa Mk.LineSegments, axes[1].scene.plots) == 1
    ## One labelled zone, drawn as halo plus text.
    @test count(p -> p isa Mk.Text, axes[1].scene.plots) == 2
    ## A sequential map on a log scale with a zone the geojson lacks,
    ## which is warned about rather than fatal.
    seq = @test_logs (:warn, r"not in the geojson") plot_zone_map(
        [3, 30, 300], ["alpha", "beta_two", "nowhere"]; geojson = path,
        provinces = ["P", "Q"], scale = log10, colorbar_label = "Cases")
    @test seq isa Mk.Figure
    ## Values and zones must pair off.
    @test_throws ErrorException plot_zone_map([1.0], ["alpha", "gamma"];
        geojson = path, provinces = ["P", "Q"])
end

@testitem "plot_zone_map_panels shares the polygons across panels" setup=[
    ZoneGeojson
] begin
    using CairoMakie: Makie as Mk
    using BVDOutbreakSize: plot_zone_map_panels
    path = write_zone_geojson()
    zones = ["alpha", "beta_two", "delta"]
    fig = plot_zone_map_panels(
        [
            (values = [1.4, 0.8, 1.0], zones, title = "R", diverging_at = 1.0),
            (values = [0, 5, 20], zones, title = "Forecast",
                scale = Mk.pseudolog10),
            (values = [10, 50, 200], zones, title = "Cumulative",
                scale = log10)];
        geojson = path, provinces = ["P", "Q"], title = "Zones", ncols = 2)
    @test fig isa Mk.Figure
    axes = [x for x in fig.content if x isa Mk.Axis]
    @test length(axes) == 3
    @test [ax.title[] for ax in axes] == ["R", "Forecast", "Cumulative"]
    @test count(x -> x isa Mk.Colorbar, fig.content) == 3
end

@testitem "plot_zone_map draws the packaged health zones" setup=[
    ZoneGeojson
] begin
    using CairoMakie: Makie as Mk
    using BVDOutbreakSize: plot_zone_map, load_health_zones_geojson,
                           zone_geojson_path
    ## The packaged polygons ship with the zone data; a checkout without
    ## them skips this rather than failing on the synthetic-geojson tests.
    if isfile(zone_geojson_path())
        zones = load_health_zones_geojson()
        keys = [z.zone for z in zones]
        @test length(keys) > 100
        @test allunique(keys)
        fig = plot_zone_map([1.2, 0.9], keys[1:2]; label_top = 2)
        @test fig isa Mk.Figure
    else
        @test_skip isfile(zone_geojson_path())
    end
end

@testitem "plot_rt_zones facets the top zones grouped by patch" setup=[
    ZoneGeojson
] begin
    using CairoMakie: Makie as Mk
    using Random: MersenneTwister
    using BVDOutbreakSize: plot_rt_zones
    rng = MersenneTwister(3)
    nd, n, nz = 30, 40, 6
    function traj(level)
        m = Matrix{Union{Missing, Float64}}(missing, nd, n)
        m[:, 8:n] .= level .* exp.(0.1 .* randn(rng, nd, n - 7))
        return m
    end
    rt = [traj(l) for l in (1.5, 0.8, 1.1, 2.0, 0.6, 1.0)]
    labels = ["Z$i" for i in 1:nz]
    patch = [1, 2, 1, 2, 4, 1]
    cum = [50, 500, 5, 200, 80, 20]
    fig = plot_rt_zones(rt, labels, patch; as_of_date = "2026-09-10",
        cumulative = cum, top = 4, patch_rt = [traj(1.0), traj(1.2),
            traj(1.0), traj(0.9)], ncols = 2)
    @test fig isa Mk.Figure
    axes = [x for x in fig.content if x isa Mk.Axis]
    ## The four largest zones, grouped by patch and by size within it.
    @test [ax.title[] for ax in axes] == ["Z1", "Z2", "Z4", "Z5"]
    for ax in axes
        ## Patch reference behind the zone: three ribbons each, plus the
        ## patch median line, the no-growth line and the cut-off.
        @test count(p -> p isa Mk.Band, ax.scene.plots) == 6
        @test count(p -> p isa Mk.Lines, ax.scene.plots) == 1
        @test count(p -> p isa Mk.HLines, ax.scene.plots) == 1
        @test count(p -> p isa Mk.VLines, ax.scene.plots) == 1
    end
    ylims = [ax.limits[][2] for ax in axes]
    @test all(==(ylims[1]), ylims)
    @test ylims[1][1] == 0
    @test count(x -> x isa Mk.Legend, fig.content) == 1
    ## Without a patch reference or a cut-off, and with explicit dates,
    ## every zone in the order given.
    using Dates: Date, Day
    dates = [Date(2026, 6, 1) + Day(d - 1) for d in 1:n]
    plain = plot_rt_zones(rt, labels, patch; dates, top = 10)
    paxes = [x for x in plain.content if x isa Mk.Axis]
    @test length(paxes) == nz
    @test count(p -> p isa Mk.Band, paxes[1].scene.plots) == 3
    @test count(p -> p isa Mk.VLines, paxes[1].scene.plots) == 0
    ## A reference fit per zone adds its 90% band and dashed median.
    ref = plot_rt_zones(rt, labels, patch; dates, top = 2,
        reference_rt = [traj(1.0) for _ in 1:nz], reference_label = "Frozen")
    rax = first(x for x in ref.content if x isa Mk.Axis)
    @test count(p -> p isa Mk.Band, rax.scene.plots) == 4
    @test count(p -> p isa Mk.Lines, rax.scene.plots) == 1
    leg = first(x for x in ref.content if x isa Mk.Legend)
    @test any(e -> e.label[] == "Frozen", leg.entrygroups[][1][2])
    @test_throws ErrorException plot_rt_zones(rt, labels, patch)
    @test_throws ErrorException plot_rt_zones(rt, labels[1:2], patch;
        as_of_date = "2026-09-10")
end

@testitem "plot_zone_shares draws the observed over the modelled" setup=[
    ZoneGeojson
] begin
    using CairoMakie: Makie as Mk
    using Random: MersenneTwister
    using Dates: Date, Day
    using BVDOutbreakSize: plot_zone_shares
    rng = MersenneTwister(5)
    nd, nv, nz = 40, 6, 5
    dates = [Date(2026, 7, 1) + Day(7 * (v - 1)) for v in 1:nv]
    shares = [clamp.(s .+ 0.05 .* randn(rng, nd, nv), 0, 1)
              for s in (0.5, 0.3, 0.1, 0.6, 0.4)]
    obs = [0.5 0.4 0.6 0.5 0.5 0.5; 0.3 0.3 0.2 0.3 NaN 0.3;
           0.1 0.1 0.1 0.1 0.1 0.1; 0.6 0.7 0.5 0.6 0.6 0.6;
           0.4 0.3 0.5 0.4 0.4 0.4]
    fig = plot_zone_shares(shares, obs, dates, ["A", "B", "C", "D", "E"];
        zone_patch = [1, 1, 1, 2, 2], top = 3, pred_draws = shares)
    @test fig isa Mk.Figure
    axes = [x for x in fig.content if x isa Mk.Axis]
    ## The three largest observed shares, grouped by patch.
    @test [ax.title[] for ax in axes] == ["A", "D", "E"]
    for ax in axes
        @test count(p -> p isa Mk.Band, ax.scene.plots) == 6
        @test count(p -> p isa Mk.Scatter, ax.scene.plots) == 1
        @test ax.limits[][2][1] == 0
    end
    plain = plot_zone_shares(shares, obs, dates, ["A", "B", "C", "D", "E"])
    pax = first(x for x in plain.content if x isa Mk.Axis)
    @test count(p -> p isa Mk.Band, pax.scene.plots) == 3
    @test_throws ErrorException plot_zone_shares(shares, obs[1:2, :], dates,
        ["A", "B", "C", "D", "E"])
end

@testitem "plot_zone_forecast ranks the zones and marks the observed" setup=[
    ZoneGeojson
] begin
    using CairoMakie: Makie as Mk
    using DataFrames: DataFrame
    using Random: MersenneTwister
    using BVDOutbreakSize: plot_zone_forecast
    rng = MersenneTwister(7)
    draws = [round.(Int, m .* exp.(0.3 .* randn(rng, 100)))
             for m in (30, 5, 12, 40, 8)]
    labels = ["A", "B", "C", "D", "E"]
    patch = [1, 2, 1, 2, 4]
    fig = plot_zone_forecast(draws, labels, patch;
        observed = [35, missing, 10, 50, 9], top = 4)
    @test fig isa Mk.Figure
    ax = first(x for x in fig.content if x isa Mk.Axis)
    ## The four largest medians, grouped by patch: A, C then D, E.
    pos, labs = ax.yticks[]
    @test labs == ["A", "C", "D", "E"]
    @test collect(pos) == [4, 3, 2, 1]
    ## Two bars per zone, a median dot per zone, one observed layer.
    @test count(p -> p isa Mk.Lines, ax.scene.plots) == 8
    @test count(p -> p isa Mk.Scatter, ax.scene.plots) == 5
    @test ax.limits[][1][1] < 0 < ax.limits[][1][2]
    ## The table form, without the 50% interval or the observed.
    tbl = DataFrame(label = labels, patch = patch,
        median = [30.0, 5, 12, 40, 8], lo90 = [20.0, 3, 8, 30, 5],
        hi90 = [45.0, 8, 18, 55, 12])
    tfig = plot_zone_forecast(tbl; top = 10)
    tax = first(x for x in tfig.content if x isa Mk.Axis)
    @test count(p -> p isa Mk.Lines, tax.scene.plots) == 5
    @test count(p -> p isa Mk.Scatter, tax.scene.plots) == 5
    @test_throws ErrorException plot_zone_forecast(
        DataFrame(label = labels, median = [1.0, 2, 3, 4, 5]))
    ## The report's zone forecast layout: labels under `zone`, the patch
    ## as its label, the 60% interval as the thick bar, and the appended
    ## "Patch total" rows left out.
    rep = DataFrame(zone = ["A", "B", "Patch total", "C"],
        patch = ["Ituri", "Nord-Kivu", "Ituri", "Ituri"],
        central_estimate = [30.0, 5, 42, 12], lower_90 = [20.0, 3, 30, 8],
        lower_60 = [24.0, 4, 35, 10], upper_60 = [36.0, 6, 50, 14],
        upper_90 = [45.0, 8, 60, 18])
    rfig = plot_zone_forecast(rep)
    rax = first(x for x in rfig.content if x isa Mk.Axis)
    @test rax.yticks[][2] == ["A", "C", "B"]
    @test count(p -> p isa Mk.Lines, rax.scene.plots) == 6
    @test_throws ErrorException plot_zone_forecast(
        DataFrame(zone = ["A"], patch = ["Mars"], median = [1.0],
        lo90 = [0.5], hi90 = [2.0]))
end

@testitem "zone_summary_table summarises per-zone draws" begin
    using BVDOutbreakSize: zone_summary_table
    draws = [collect(1.0:100.0), collect(10.0:10.0:1000.0)]
    tbl = zone_summary_table(draws, ["A", "B"], [1, 2];
        observed = [50, missing])
    @test names(tbl) == ["label", "patch", "median", "lo90", "hi90", "lo50",
        "hi50", "observed"]
    @test tbl.median == [50.5, 505.0]
    @test tbl.lo90[1] < tbl.lo50[1] < tbl.median[1] < tbl.hi50[1] < tbl.hi90[1]
    @test ismissing(tbl.observed[2])
    @test !("observed" in names(zone_summary_table(draws, ["A", "B"], [1, 2])))
end

@testitem "plot_zone_comparison offsets each fit within a zone's row" setup=[
    ZoneGeojson
] begin
    using CairoMakie: Makie as Mk
    using DataFrames: DataFrame
    using BVDOutbreakSize: plot_zone_comparison
    labels = ["A", "B", "C", "D"]
    patch = [1, 2, 1, 4]
    tbl(m) = DataFrame(label = labels, patch = patch, median = m,
        lo90 = m .* 0.7, hi90 = m .* 1.4, lo50 = m .* 0.9, hi50 = m .* 1.1)
    live = tbl([1.5, 0.8, 1.2, 1.0])
    ## The frozen fit lacks the 50% columns and one zone.
    frozen = DataFrame(label = labels[1:3], patch = patch[1:3],
        median = [1.3, 0.9, 1.1], lo90 = [1.0, 0.6, 0.8],
        hi90 = [1.7, 1.3, 1.5])
    fig = plot_zone_comparison(["Live" => live, "Frozen" => frozen];
        top = 3, reference_line = 1.0, xlabel = "R")
    @test fig isa Mk.Figure
    ax = first(x for x in fig.content if x isa Mk.Axis)
    ## The three largest live medians, grouped by patch: A, C then D.
    @test ax.yticks[][2] == ["A", "C", "D"]
    ## Live: two bars and a marker per zone; frozen: one bar and a marker
    ## for the two of those zones it carries; one reference rule.
    @test count(p -> p isa Mk.Lines, ax.scene.plots) == 3 * 2 + 2
    @test count(p -> p isa Mk.Scatter, ax.scene.plots) == 3 + 2
    @test count(p -> p isa Mk.VLines, ax.scene.plots) == 1
    ## Ranking by the second series instead, and the guards.
    alt = plot_zone_comparison(["Live" => live, "Frozen" => frozen];
        top = 2, rank_by = 2)
    aax = first(x for x in alt.content if x isa Mk.Axis)
    @test aax.yticks[][2] == ["A", "C"]
    @test_throws ErrorException plot_zone_comparison(Pair{String, DataFrame}[])
    @test_throws ErrorException plot_zone_comparison(
        ["x" => DataFrame(label = labels, median = [1.0, 2, 3, 4])])
end

@testitem "plot_zone_ranking sorts by P(R > 1), greys level-only" setup=[
    ZoneGeojson
] begin
    using CairoMakie: Makie as Mk
    using DataFrames: DataFrame, Not
    using BVDOutbreakSize: plot_zone_ranking
    ov = DataFrame(label = ["A", "B", "C", "D"], patch = [1, 2, 1, 4],
        rt_median = [1.4, 0.9, 1.1, 1.0], rt_lo90 = [1.1, 0.6, 0.8, 0.7],
        rt_hi90 = [1.8, 1.2, 1.5, 1.4], p_rt_above_one = [0.95, 0.3, 0.7, 0.5],
        walking = [true, true, false, true])
    fig = plot_zone_ranking(ov)
    @test fig isa Mk.Figure
    axes = [x for x in fig.content if x isa Mk.Axis]
    @test length(axes) == 2
    pos, labs = axes[1].yticks[]
    @test labs == ["A", "C", "D", "B"]
    @test count(p -> p isa Mk.Scatter, axes[1].scene.plots) == 4
    @test count(p -> p isa Mk.Lines, axes[2].scene.plots) == 4
    @test count(p -> p isa Mk.VLines, axes[2].scene.plots) == 1
    ## The level-only zone is drawn hollow; the legend says so.
    leg = first(x for x in fig.content if x isa Mk.Legend)
    @test any(e -> occursin("Level only", e.label[]), leg.entrygroups[][1][2])
    ## Without a `walking` column every zone is its own walk, and
    ## `max_zones` truncates the ranking.
    short = plot_zone_ranking(ov[:, Not(:walking)]; max_zones = 2)
    sax = first(x for x in short.content if x isa Mk.Axis)
    @test sax.yticks[][2] == ["A", "C"]
    @test_throws ErrorException plot_zone_ranking(ov[:, Not(:rt_median)])
    ## The report's overview layout: `zone`, the patch label, and
    ## `p_R_above_1` with NaN for a zone below the reporting floor.
    rep = DataFrame(zone = ["A", "B", "C"],
        patch = ["Ituri", "Nord-Kivu", "Other provinces"],
        rt_median = [1.4, 0.9, NaN], rt_lo90 = [1.1, 0.6, NaN],
        rt_hi90 = [1.8, 1.2, NaN], p_R_above_1 = [0.95, 0.3, NaN],
        walking = [true, true, false])
    rfig = plot_zone_ranking(rep)
    rax = first(x for x in rfig.content if x isa Mk.Axis)
    @test rax.yticks[][2] == ["A", "B"]
end
