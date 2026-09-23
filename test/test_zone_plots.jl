## Smoke tests for the health-zone figures on small synthetic inputs,
## checking the layout each promises. The map tests read a synthetic
## geojson written to a temporary file.

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
        a = ft(
            """{"zone": "alpha", "label": "Alpha", "province": "P"}""",
            poly([sq(0, 0, 1, 1)])
        )
        b = ft(
            """{"zone": "", "label": "Beta Two", "province": "P"}""",
            poly([[[1, 0], [2, 0], [2, 1], [1, 1], [1, 0.5], [1, 0]]])
        )
        c = ft(
            """{"zone": "gamma", "label": "Gamma", "province": "P"}""",
            """{"type": "MultiPolygon", "coordinates": [
            [$(sq(0, 1, 1, 2)), $(sq(0.3, 1.3, 0.6, 1.6))],
            [$(sq(3, 3, 3.5, 3.5))]]}"""
        )
        d = ft(
            """{"zone": "delta", "label": "Delta", "province": "Q"}""",
            poly([sq(0, -1, 1, 0)])
        )
        e = ft(
            """{"zone": "omega", "label": "Omega", "province": "Z"}""",
            poly([sq(5, 5, 6, 6)])
        )
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

@testitem "load_health_zones_geojson keys, filters and dissolves" setup = [
    ZoneGeojson,
] begin
    using BVDOutbreakSize: load_health_zones_geojson, _province_outlines,
        _zone_table, _feature_polygons
    path = write_zone_geojson()
    zones = load_health_zones_geojson(path; provinces = ["P", "Q"])
    ## The province outside the filter is dropped; the empty `zone`
    ## property falls back to the folded label.
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
    @test_throws ErrorException load_health_zones_geojson(
        path;
        provinces = ["none"]
    )
    ## Only polygons and multipolygons are drawn.
    @test_throws ErrorException _feature_polygons(
        Dict("type" => "Point", "coordinates" => [0.0, 0.0])
    )
end

@testitem "zone colourbar ticks thin to the round values in range" begin
    using CairoMakie: Makie as Mk
    using BVDOutbreakSize: _zone_cbar_ticks

    ## An identity scale keeps Makie's own ticks.
    @test _zone_cbar_ticks((0.0, 10.0), identity) === Mk.automatic
    ## A count scale from zero: round values, none between zero and one.
    ticks, labels = _zone_cbar_ticks((0.0, 30.0), Mk.pseudolog10)
    @test 0.0 in ticks
    @test !any(t -> 0 < t < 1, ticks)
    @test length(ticks) <= 7
    @test labels == string.(round.(Int, ticks))
    ## A range too wide for any mantissa set falls back to Makie's ticks.
    @test _zone_cbar_ticks((1.0e-4, 1.0e8), log10) === Mk.automatic
end

@testitem "plot_zone_map colours the matched zones and greys the rest" setup = [
    ZoneGeojson,
] begin
    using CairoMakie: Makie as Mk
    using BVDOutbreakSize: plot_zone_map
    path = write_zone_geojson()
    fig = plot_zone_map(
        [1.6, 0.7], ["alpha", "gamma"]; geojson = path,
        provinces = ["P", "Q"], diverging_at = 1.0, lower = [1.2, 0.4],
        upper = [2.1, 1.1],
        title = "R", colorbar_label = "R", label_top = 1
    )
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
        provinces = ["P", "Q"], scale = log10, colorbar_label = "Cases"
    )
    @test seq isa Mk.Figure
    ## Values and zones must pair off.
    @test_throws ErrorException plot_zone_map(
        [1.0], ["alpha", "gamma"];
        geojson = path, provinces = ["P", "Q"]
    )
end

@testitem "plot_zone_map_panels shares the polygons across panels" setup = [
    ZoneGeojson,
] begin
    using CairoMakie: Makie as Mk
    using BVDOutbreakSize: plot_zone_map_panels
    path = write_zone_geojson()
    zones = ["alpha", "beta_two", "delta"]
    fig = plot_zone_map_panels(
        [
            (values = [1.4, 0.8, 1.0], zones, title = "R", diverging_at = 1.0),
            (
                values = [0, 5, 20], zones, title = "Forecast",
                scale = Mk.pseudolog10,
            ),
            (
                values = [10, 50, 200], zones, title = "Cumulative",
                scale = log10,
            ),
        ];
        geojson = path, provinces = ["P", "Q"], title = "Zones", ncols = 2
    )
    @test fig isa Mk.Figure
    axes = [x for x in fig.content if x isa Mk.Axis]
    @test length(axes) == 3
    @test [ax.title[] for ax in axes] == ["R", "Forecast", "Cumulative"]
    @test count(x -> x isa Mk.Colorbar, fig.content) == 3
end

@testitem "plot_zone_map draws the packaged health zones" setup = [
    ZoneGeojson,
] begin
    using CairoMakie: Makie as Mk
    using BVDOutbreakSize: plot_zone_map, load_health_zones_geojson,
        zone_geojson_path
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
