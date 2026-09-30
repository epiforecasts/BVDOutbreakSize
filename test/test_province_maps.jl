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

@testitem "province_map_estimates gives one row per source province" begin
    using BVDOutbreakSize: province_map_estimates
    using Dates: Date
    rt = [collect(0.5:0.01:1.5), fill(2.0, 101)]
    fc = [collect(0.0:100.0), collect(100.0:200.0)]
    hist(c) = (; days = [1, 5], counts = c)
    confirmed = Dict("a" => hist([3, 7]), "c" => hist([1, 2]))
    deaths = Dict("a" => hist([0, 1]))
    est = province_map_estimates(
        rt, fc; confirmed_history = confirmed, death_history = deaths,
        cutoff = Date(2026, 9, 1), patch_names = ["x", "y"],
        patch_labels = ["X", "Pool"],
        members = Dict("x" => ["a"], "y" => ["b", "c"])
    )
    @test est.province == ["a", "b", "c"]
    @test est.patch == ["X", "Pool", "Pool"]
    @test est.pooled == [0, 1, 1]
    @test isequal(est.cases, [7, missing, 2])
    @test isequal(est.deaths, [1, missing, missing])
    @test est.R_T_median ≈ [1.0, 2.0, 2.0]
    @test est.R_T_lower ≈ [0.55, 2.0, 2.0]
    @test est.R_T_upper ≈ [1.45, 2.0, 2.0]
    @test est.p_rt_above_one ≈ [50 / 101, 1.0, 1.0]
    @test est.forecast_median ≈ [50.0, 150.0, 150.0]
    @test est.forecast_lower ≈ [5.0, 105.0, 105.0]
    @test est.forecast_upper ≈ [95.0, 195.0, 195.0]
    @test all(==("2026-09-01"), est.as_of)
    ## From a chain's per-patch reproduction number, one vector per draw.
    chn = (; R_T_patch = [[rt[1][i], rt[2][i], 9.0] for i in 1:101])
    from_chn = province_map_estimates(
        chn, fc; n_patches = 2, confirmed_history = confirmed,
        death_history = deaths, cutoff = Date(2026, 9, 1),
        patch_names = ["x", "y"], patch_labels = ["X", "Pool"],
        members = Dict("x" => ["a"], "y" => ["b", "c"])
    )
    @test isequal(from_chn, est)
end

@testitem "province_map_estimates adds ascertainment when given" begin
    using BVDOutbreakSize: province_map_estimates
    using Dates: Date
    d = [collect(0.5:0.01:1.5), fill(2.0, 101)]
    kw = (;
        confirmed_history = Dict(), death_history = Dict(),
        cutoff = Date(2026, 9, 1), patch_names = ["x", "y"],
        patch_labels = ["X", "Pool"],
        members = Dict("x" => ["a"], "y" => ["b", "c"]),
    )
    est = province_map_estimates(d, d; ascertainment = d, kw...)
    @test est.ascertainment_median ≈ [1.0, 2.0, 2.0]
    @test est.ascertainment_lower ≈ [0.55, 2.0, 2.0]
    @test est.ascertainment_upper ≈ [1.45, 2.0, 2.0]
    @test !(
        :ascertainment_median in propertynames(
            province_map_estimates(d, d; kw...)
        )
    )
end

@testitem "rt_quantile_table gives daily quantiles in long format" begin
    using BVDOutbreakSize: rt_quantile_table
    using Dates: Date
    ## Two areas over four days, 101 draws each. The second is unreported
    ## (NaN) before day 3, and neither is written before `from`.
    a = repeat(collect(0.0:0.01:1.0), 1, 4) .+ [0 1 2 3]
    b = copy(a)
    b[:, 1:2] .= NaN
    t = rt_quantile_table(
        [a, b], ["A", "B"]; cutoff = Date(2026, 9, 4), n = 4, from = 2
    )
    @test t.area == ["A", "A", "A", "B", "B"]
    @test t.date == Date.(
        [
            "2026-09-02", "2026-09-03", "2026-09-04",
            "2026-09-03", "2026-09-04",
        ]
    )
    @test t.median ≈ [1.5, 2.5, 3.5, 2.5, 3.5]
    @test t.lower_90 ≈ [1.05, 2.05, 3.05, 2.05, 3.05]
    @test t.upper_90 ≈ [1.95, 2.95, 3.95, 2.95, 3.95]
    @test t.lower_50 ≈ [1.25, 2.25, 3.25, 2.25, 3.25]
    @test t.upper_50 ≈ [1.75, 2.75, 3.75, 2.75, 3.75]
end

@testitem "weekly_count_table sums increments into weeks ending at the cut-off" begin
    using BVDOutbreakSize: weekly_count_table
    using Dates: Date
    ## Vintages on days 2, 9, 10, 16 and 20 of a 20-day grid. Weeks end on
    ## days 20, 13 and 6; the week (−1, 6] starts before the first vintage,
    ## whose increment lumps in everything before it, so it is left out.
    days = [2, 9, 10, 16, 20]
    inc = [5 1 2 3 4; 0 0 1 0 0]
    t = weekly_count_table(
        days, inc, ["A", "B"]; cutoff = Date(2026, 9, 20), n = 20,
        weeks = 3
    )
    @test t.area == ["A", "A", "B", "B"]
    @test t.date == Date.(
        [
            "2026-09-13", "2026-09-20",
            "2026-09-13", "2026-09-20",
        ]
    )
    @test t.count == [3, 7, 1, 0]
end
