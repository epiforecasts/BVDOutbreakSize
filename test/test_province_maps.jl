## Tests for the province maps, which colour the health-zone polygons by
## the patch each zone's province belongs to. They read a synthetic geojson
## whose provinces carry real patch keys.

@testsnippet ProvinceGeojson begin
    using CairoMakie
    CairoMakie.activate!(type = "png")

    ## Ituri holds two unit squares side by side, Nord-Kivu one below them,
    ## Tshopo (pooled into `other`) one to the left and Haut-Uele one above.
    function write_province_geojson()
        sq(x0, y0, x1, y1) = [[x0, y0], [x1, y0], [x1, y1], [x0, y1], [x0, y0]]
        ft(zone, prov, x0, y0) = """{"type": "Feature",
        "properties": {"zone": "$zone", "label": "$zone",
        "province": "$prov"}, "geometry": {"type": "Polygon",
        "coordinates": [$(sq(x0, y0, x0 + 1, y0 + 1))]}}"""
        fts = join(
            (
                ft("a", "ituri", 0, 0), ft("b", "ituri", 1, 0),
                ft("c", "nord_kivu", 0, -1), ft("d", "tshopo", -1, 0),
                ft("e", "haut_uele", 0, 1),
            ), ", "
        )
        path = tempname() * ".geojson"
        write(path, """{"type": "FeatureCollection", "features": [$fts]}""")
        return path
    end
    provinces = ["Ituri", "Nord-Kivu", "Haut-Uele", "Tshopo"]
end

@testitem "province_map_summary gives the median and 90% interval" begin
    using BVDOutbreakSize: province_map_summary
    draws = [collect(1.0:101.0), collect(201.0:301.0)]
    s = province_map_summary(draws)
    @test s.values == [51.0, 251.0]
    @test s.lower ≈ [6.0, 206.0]
    @test s.upper ≈ [96.0, 296.0]
    wide = province_map_summary(draws; level = 0.5)
    @test wide.lower ≈ [26.0, 226.0]
    @test wide.upper ≈ [76.0, 276.0]
    ## From a chain's per-patch deterministic, one vector per draw.
    chn = (; R_T_patch = [[d, 10 + d, 20 + d] for d in 1.0:101.0])
    c = province_map_summary(chn, :R_T_patch, 2)
    @test c.values == [51.0, 61.0]
    @test c.lower ≈ [6.0, 16.0]
end

@testitem "province_zone_values spreads each patch over its zones" setup = [
    ProvinceGeojson,
] begin
    using BVDOutbreakSize: province_zone_values
    path = write_province_geojson()
    ## Patches in PROVINCE_NAMES order: ituri, nord_kivu, haut_uele, other.
    z = province_zone_values(
        [1.0, 2.0, 3.0, 4.0]; lower = [0.5, 1.5, 2.5, 3.5],
        upper = [1.5, 2.5, 3.5, 4.5], geojson = path, provinces
    )
    @test z.zones == ["a", "b", "c", "d", "e"]
    @test z.values == [1.0, 1.0, 2.0, 4.0, 3.0]
    @test z.lower == [0.5, 0.5, 1.5, 3.5, 2.5]
    @test z.upper == [1.5, 1.5, 2.5, 4.5, 3.5]
    ## Without bounds the panel carries no bounds, so a sequential map
    ## is not washed out.
    nob = province_zone_values(
        [1.0, 2.0, 3.0, 4.0]; geojson = path, provinces
    )
    @test !haskey(nob, :lower) && !haskey(nob, :upper)
    ## Fewer values than patches leaves the later patches' zones out.
    three = province_zone_values([1.0, 2.0, 3.0]; geojson = path, provinces)
    @test three.zones == ["a", "b", "c", "e"]
    @test_throws ErrorException province_zone_values(
        [1.0, 2.0, 3.0, 4.0]; lower = [1.0], upper = [2.0, 3.0, 4.0, 5.0],
        geojson = path, provinces
    )
end

@testitem "plot_province_map labels provinces rather than zones" setup = [
    ProvinceGeojson,
] begin
    using CairoMakie: Makie as Mk
    using BVDOutbreakSize: plot_province_map
    path = write_province_geojson()
    fig = plot_province_map(
        [
            (values = [10.0, 100.0, 1.0, 5.0], title = "Infections", scale = log10),
            (
                values = [1.3, 0.8, 1.1, 0.9], lower = [1.1, 0.6, 0.7, 0.5],
                upper = [1.5, 0.95, 1.4, 1.2], title = "R",
                diverging_at = 1.0,
            ),
        ];
        geojson = path, provinces, ncols = 2
    )
    @test fig isa Mk.Figure
    axes = [x for x in fig.content if x isa Mk.Axis]
    @test [ax.title[] for ax in axes] == ["Infections", "R"]
    @test count(x -> x isa Mk.Colorbar, fig.content) == 2
    for ax in axes
        texts = [p for p in ax.scene.plots if p isa Mk.Text]
        ## One halo and one label layer, each naming the four provinces.
        @test length(texts) == 2
        @test sort(texts[2].text[]) == sort(provinces)
    end
    ## Haut-Uele and Other provinces straddle one, so their zones (e, d)
    ## are washed out on the R panel.
    polys = [p for p in axes[2].scene.plots if p isa Mk.Poly]
    @test length(polys) == 2
    @test length(polys[2][1][]) == 2
    ## A single panel is one axis.
    one = plot_province_map(
        [1.0, 2.0, 3.0, 4.0]; geojson = path, provinces, title = "One"
    )
    @test count(x -> x isa Mk.Axis, one.content) == 1
end

@testitem "the packaged zones each fall in one patch" begin
    using BVDOutbreakSize: load_health_zones_geojson, province_zone_values,
        PROVINCE_NAMES
    zones = load_health_zones_geojson()
    @test allunique([z.zone for z in zones])
    z = province_zone_values(Float64.(eachindex(PROVINCE_NAMES)))
    @test sort(z.zones) == sort([z.zone for z in zones])
    @test sort(unique(z.values)) == Float64.(eachindex(PROVINCE_NAMES))
end
