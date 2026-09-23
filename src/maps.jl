## Health-zone and province maps. The zone map code (from ZONE_MAP_PROVINCES
## to plot_zone_map_panels) is shared with the health-zone model.

export ZONE_MAP_PROVINCES, zone_key, zone_geojson_path,
    load_health_zones_geojson, plot_zone_map, plot_zone_map_panels

## Health-zone figures. Each takes plain vectors, matrices or a DataFrame
## rather than a chain.

"""
    ZONE_MAP_PROVINCES

The provinces whose health zones the maps draw. Zones outside them are
dropped when the geojson is read.
"""
const ZONE_MAP_PROVINCES = [
    "Ituri", "Nord-Kivu", "Sud-Kivu", "Haut-Uele",
    "Tshopo", "Bas-Uele", "Sud-Ubangi",
]

## Fill colours by patch index: the patch figures' three, then the pooled
## fourth patch.
const _ZONE_PATCH_COLOURS = [:firebrick, :steelblue, :seagreen, :darkorange]

"""
$(TYPEDSIGNATURES)

Path of the packaged health-zone geojson, `src/assets/health_zones.geojson`.
"""
zone_geojson_path() = joinpath(
    pkgdir(@__MODULE__), "src", "assets",
    "health_zones.geojson"
)

"""
$(TYPEDSIGNATURES)

Fold a zone name to the snake_case key the zone data and the geojson share:
accents stripped, lower-cased, and every run of non-alphanumerics collapsed
to one underscore, so `"Bili (Bas-Uele)"` becomes `"bili_bas_uele"`.
"""
function zone_key(name::AbstractString)
    s = Base.Unicode.normalize(String(name); stripmark = true, casefold = true)
    s = replace(s, r"[^a-z0-9]+" => "_")
    return String(strip(s, '_'))
end

## One geojson ring as Makie points, without the closing repeat of the
## first point.
function _ring_points(ring)
    pts = CairoMakie.Point2f[(Float32(p[1]), Float32(p[2])) for p in ring]
    length(pts) > 1 && pts[end] == pts[1] && pop!(pts)
    return pts
end

## A geojson polygon (exterior ring, then holes) as one Makie polygon.
function _zone_polygon(coords)
    return CairoMakie.Makie.Polygon(
        _ring_points(coords[1]),
        Vector{CairoMakie.Point2f}[_ring_points(r) for r in coords[2:end]]
    )
end

## Every part of a feature's geometry as Makie polygons.
function _feature_polygons(geom)
    t = geom["type"]
    t == "Polygon" && return [_zone_polygon(geom["coordinates"])]
    t == "MultiPolygon" &&
        return [_zone_polygon(c) for c in geom["coordinates"]]
    return error("zone geojson: unsupported geometry type `$t`.")
end

## Signed shoelace area and area-weighted centroid of one ring.
function _ring_centroid(pts)
    n = length(pts)
    a = cx = cy = 0.0
    for i in 1:n
        x1, y1 = pts[i]
        x2, y2 = pts[mod1(i + 1, n)]
        w = x1 * y2 - x2 * y1
        a += w
        cx += (x1 + x2) * w
        cy += (y1 + y2) * w
    end
    a == 0 && return (mean(p[1] for p in pts), mean(p[2] for p in pts), 0.0)
    return (cx / (3a), cy / (3a), abs(a) / 2)
end

## Centroid of a zone's largest part, where its label goes.
function _zone_centroid(polys)
    best = (0.0, 0.0, -1.0)
    for p in polys
        c = _ring_centroid(
            CairoMakie.Makie.GeometryBasics.coordinates(p.exterior)
        )
        c[3] > best[3] && (best = c)
    end
    return CairoMakie.Point2f(best[1], best[2])
end

"""
$(TYPEDSIGNATURES)

Read the health zones of `provinces` (matched through [`zone_key`](@ref))
from a geojson, one NamedTuple per zone of `zone` (the key), `label`,
`province`, `polygons` (one Makie polygon per part, holes included) and
`centroid` (of the largest part). A feature is named by its `label` (or
`nom`) property and keyed by its `zone` property when that is non-empty,
otherwise by [`zone_key`](@ref) of the name.
"""
function load_health_zones_geojson(
        path::AbstractString = zone_geojson_path();
        provinces::AbstractVector = ZONE_MAP_PROVINCES
    )
    gj = JSON.parsefile(path)
    keep = Set(zone_key(String(p)) for p in provinces)
    zones = @NamedTuple{
        zone::String, label::String, province::String,
        polygons::Vector{CairoMakie.Makie.Polygon},
        centroid::CairoMakie.Point2f,
    }[]
    for f in gj["features"]
        props = f["properties"]
        prov = something(get(props, "province", nothing), "")
        zone_key(prov) in keep || continue
        name = something(
            get(props, "label", nothing),
            get(props, "nom", nothing), ""
        )
        k = something(get(props, "zone", nothing), "")
        isempty(k) && (k = zone_key(name))
        polys = _feature_polygons(f["geometry"])
        push!(
            zones,
            (;
                zone = String(k), label = String(name),
                province = String(prov), polygons = polys,
                centroid = _zone_centroid(polys),
            )
        )
    end
    isempty(zones) && error("zone geojson `$path` holds no zone in $provinces.")
    return zones
end

## The zone records as one column per field, the form the map code indexes.
function _zone_table(zones::AbstractVector{<:NamedTuple})
    return (;
        key = [z.zone for z in zones], label = [z.label for z in zones],
        province = [z.province for z in zones],
        polygons = [z.polygons for z in zones],
        centroid = [z.centroid for z in zones],
    )
end

## Every edge of every ring of the zones `idx`, as point pairs.
function _zone_edges(geo, idx)
    edges = Tuple{CairoMakie.Point2f, CairoMakie.Point2f}[]
    for z in idx, poly in geo.polygons[z]

        for ring in (poly.exterior, poly.interiors...)
            pts = CairoMakie.Makie.GeometryBasics.coordinates(ring)
            m = length(pts)
            for i in 1:m
                push!(edges, (pts[i], pts[mod1(i + 1, m)]))
            end
        end
    end
    return edges
end

## Whether `p` lies inside the ring `pts`, by the even-odd rule.
function _in_ring(p, pts)
    inside = false
    j = length(pts)
    for i in eachindex(pts)
        xi, yi = pts[i]
        xj, yj = pts[j]
        if (yi > p[2]) != (yj > p[2]) &&
                p[1] < (xj - xi) * (p[2] - yi) / (yj - yi) + xi
            inside = !inside
        end
        j = i
    end
    return inside
end

## Whether `p` lies inside the polygon `poly`, holes excluded.
function _in_polygon(p, poly)
    coords = CairoMakie.Makie.GeometryBasics.coordinates
    return _in_ring(p, coords(poly.exterior)) &&
        !any(r -> _in_ring(p, coords(r)), poly.interiors)
end

## Province boundaries by dissolving the zone polygons: an edge lies on the
## outline when the points `δ` degrees either side of its midpoint fall one
## inside a zone of the province and one in none. Returns the segments as
## consecutive point pairs for `linesegments!`.
function _province_outlines(geo; δ::Real = 2.0e-3)
    out = CairoMakie.Point2f[]
    for prov in unique(geo.province)
        idx = findall(==(prov), geo.province)
        polys = [poly for z in idx for poly in geo.polygons[z]]
        boxes = [
            let r = CairoMakie.Makie.GeometryBasics.coordinates(
                    poly.exterior
                )
                (extrema(p[1] for p in r), extrema(p[2] for p in r))
            end
                for poly in polys
        ]
        inside(q) = any(eachindex(polys)) do k
            (bx, by) = boxes[k]
            bx[1] <= q[1] <= bx[2] && by[1] <= q[2] <= by[2] &&
                _in_polygon(q, polys[k])
        end
        for (a, b) in _zone_edges(geo, idx)
            len = hypot(b[1] - a[1], b[2] - a[2])
            len == 0 && continue
            nx, ny = -(b[2] - a[2]) / len * δ, (b[1] - a[1]) / len * δ
            mx, my = (a[1] + b[1]) / 2, (a[2] + b[2]) / 2
            inside((mx + nx, my + ny)) == inside((mx - nx, my - ny)) ||
                push!(out, a, b)
        end
    end
    return out
end

## Every part of the listed zones flattened to one polygon list, with the
## zone each part belongs to.
function _zone_parts(geo, idx)
    parts = CairoMakie.Makie.Polygon[]
    owner = Int[]
    for z in idx, p in geo.polygons[z]

        push!(parts, p)
        push!(owner, z)
    end
    return parts, owner
end

## Colour range for one map: symmetric about `diverging_at` on the colour
## scale, or the span of the values drawn (widened by one when constant).
function _zone_colorrange(v, diverging_at, scale)
    if diverging_at === nothing
        lo, hi = extrema(v)
        return lo == hi ? (lo, hi + 1) : (lo, hi)
    end
    h = maximum(abs.(scale.(v) .- scale(diverging_at)); init = 0.0)
    h == 0 && (h = 1.0)
    inv = CairoMakie.Makie.inverse_transform(scale)
    return (inv(scale(diverging_at) - h), inv(scale(diverging_at) + h))
end

## Colourbar ticks for a transformed colour scale: the round values inside
## the range, thinned to at most seven. A range from zero is a count scale.
## Ticks between zero and one are left out. An identity scale keeps
## Makie's automatic ticks.
function _zone_cbar_ticks(crange, scale)
    scale === identity && return CairoMakie.Makie.automatic
    lo, hi = crange
    for mantissas in ((1, 1.5, 2, 3, 5, 7), (1, 2, 5), (1, 3), (1,))
        ticks = [
            0.0; vec(
                [
                    round(m * 10.0^e; sigdigits = 3)
                        for m in mantissas, e in -3:7
                ]
            )
        ]
        ticks = [
            t
                for t in ticks
                if lo <= t <= hi && isfinite(scale(t)) &&
                !(lo <= 0 && 0 < t < 1)
        ]
        if length(ticks) <= 7
            labels = [
                t == round(t) ? string(round(Int, t)) : string(t)
                    for t in ticks
            ]
            return (ticks, labels)
        end
    end
    return CairoMakie.Makie.automatic
end

## The first `top` zones of `order` whose centroid is not within `gap`
## degrees (east-west, north-south) of a label already placed.
function _zone_label_pick(
        centroids, order, top::Integer;
        gap = (1.0, 0.35)
    )
    picked = Int[]
    for z in order
        length(picked) >= top && break
        c = centroids[z]
        clash = any(picked) do q
            abs(centroids[q][1] - c[1]) < gap[1] &&
                abs(centroids[q][2] - c[2]) < gap[2]
        end
        clash || push!(picked, z)
    end
    return picked
end

## One choropleth into `gl[1, 1]` with its colourbar in `gl[1, 2]`. `values`
## and `zones` pair a number with a zone key. The geo is the loaded polygon
## table and `outline` its dissolved province boundaries.
function _zone_map_panel!(
        gl, geo, outline, values, zones;
        title::AbstractString = "", colormap = nothing, colorrange = nothing,
        diverging_at = nothing, scale = identity, lower = nothing,
        upper = nothing, label_top::Integer = 10,
        colorbar_label::AbstractString = "", missing_colour = :grey85,
        fontsize = 11
    )
    length(values) == length(zones) || error(
        "plot_zone_map: $(length(values)) values for $(length(zones)) zones."
    )
    idx = Dict(k => i for (i, k) in enumerate(geo.key))
    nz = length(geo.key)
    vals = fill(NaN, nz)
    muted = falses(nz)
    unmatched = String[]
    for (j, z) in enumerate(zones)
        i = get(idx, String(z), nothing)
        if i === nothing
            push!(unmatched, String(z))
            continue
        end
        vals[i] = Float64(values[j])
        if diverging_at !== nothing && lower !== nothing && upper !== nothing
            muted[i] = lower[j] <= diverging_at <= upper[j]
        end
    end
    isempty(unmatched) || @warn "plot_zone_map: zones not in the geojson " *
        "are not drawn" unmatched
    present = findall(!isnan, vals)
    absent = findall(isnan, vals)
    cmap = colormap === nothing ?
        (diverging_at === nothing ? :Blues : CairoMakie.Reverse(:RdBu)) :
        colormap
    crange = colorrange === nothing ?
        (
            isempty(present) ? (0.0, 1.0) :
            _zone_colorrange(vals[present], diverging_at, scale)
        ) :
        colorrange
    ax = Axis(gl[1, 1]; title, aspect = CairoMakie.DataAspect())
    CairoMakie.hidedecorations!(ax)
    CairoMakie.hidespines!(ax)
    if !isempty(absent)
        parts, _ = _zone_parts(geo, absent)
        CairoMakie.poly!(
            ax, parts; color = missing_colour,
            strokecolor = :white, strokewidth = 0.4
        )
    end
    if !isempty(present)
        parts, owner = _zone_parts(geo, present)
        CairoMakie.poly!(
            ax, parts; color = vals[owner], colormap = cmap,
            colorrange = crange, colorscale = scale, strokecolor = :white,
            strokewidth = 0.6
        )
        ## Wash out the zones whose interval straddles the centre.
        mi = findall(muted)
        if !isempty(mi)
            mparts, _ = _zone_parts(geo, mi)
            CairoMakie.poly!(
                ax, mparts; color = (:white, 0.55),
                strokecolor = :white, strokewidth = 0.6
            )
        end
    end
    CairoMakie.linesegments!(ax, outline; color = :black, linewidth = 1.0)
    ## Label the largest values, or the furthest from the centre on a
    ## diverging map.
    rank = diverging_at === nothing ? vals[present] :
        abs.(scale.(vals[present]) .- scale(diverging_at))
    top = _zone_label_pick(
        geo.centroid, present[sortperm(rank; rev = true)],
        label_top
    )
    if !isempty(top)
        CairoMakie.text!(
            ax, geo.centroid[top]; text = geo.label[top],
            fontsize, align = (:center, :center), color = :white,
            strokecolor = :white, strokewidth = 3.0
        )
        CairoMakie.text!(
            ax, geo.centroid[top]; text = geo.label[top],
            fontsize, align = (:center, :center), color = :black
        )
    end
    CairoMakie.Colorbar(
        gl[1, 2]; colormap = cmap, colorrange = crange,
        scale, label = colorbar_label, width = 12,
        ticks = _zone_cbar_ticks(crange, scale)
    )
    return ax
end

"""
$(TYPEDSIGNATURES)

Choropleth of the health zones, coloured by `values`, one per key in
`zones`. Zones absent from `zones` are grey, province boundaries are black
outlines dissolved from the zone polygons, and the `label_top` zones with
the largest values are labelled at their centroids.

With `diverging_at` the colour range is symmetric about that value on a
red-blue colormap with a neutral centre. Without it the map is sequential
in blue over the range drawn. `scale` (`log10`, or
`CairoMakie.Makie.pseudolog10` when zeros are present) is applied before
colouring and the colourbar shows the untransformed values. With
`diverging_at`, per-zone `lower` and `upper` bounds wash out the zones whose
interval straddles the centre. `colormap` and `colorrange` override the
defaults, and `geojson` is read by [`load_health_zones_geojson`](@ref) for
`provinces`.
"""
function plot_zone_map(
        values::AbstractVector, zones::AbstractVector;
        geojson::AbstractString = zone_geojson_path(),
        provinces::AbstractVector = ZONE_MAP_PROVINCES,
        title::AbstractString = "", colormap = nothing, colorrange = nothing,
        lower = nothing, upper = nothing, label_top::Integer = 10,
        diverging_at = nothing, scale = identity,
        colorbar_label::AbstractString = "", size = (760, 540)
    )
    geo = _zone_table(load_health_zones_geojson(geojson; provinces))
    outline = _province_outlines(geo)
    fig = Figure(; size)
    gl = fig[1, 1] = CairoMakie.GridLayout()
    _zone_map_panel!(
        gl, geo, outline, values, zones; title, colormap,
        colorrange, diverging_at, scale, lower, upper, label_top,
        colorbar_label
    )
    return fig
end

"""
$(TYPEDSIGNATURES)

Several choropleths of the same zones side by side, one per entry of
`panels`, wrapping after `ncols`. Each panel is a NamedTuple with `values`
and `zones`, plus any of the [`plot_zone_map`](@ref) keywords (`title`,
`colormap`, `colorrange`, `diverging_at`, `scale`, `lower`, `upper`,
`colorbar_label`, `label_top`). The polygons are read once and shared.
"""
function plot_zone_map_panels(
        panels::AbstractVector{<:NamedTuple};
        geojson::AbstractString = zone_geojson_path(),
        provinces::AbstractVector = ZONE_MAP_PROVINCES, ncols::Integer = 3,
        label_top::Integer = 10, title::AbstractString = "",
        panel_size = (520, 400)
    )
    geo = _zone_table(load_health_zones_geojson(geojson; provinces))
    outline = _province_outlines(geo)
    np = length(panels)
    nc = min(ncols, np)
    nr = cld(np, nc)
    fig = Figure(;
        size = (
            panel_size[1] * nc, panel_size[2] * nr +
                (isempty(title) ? 0 : 40),
        )
    )
    for (k, panel) in enumerate(panels)
        r, c = cld(k, nc), mod1(k, nc)
        gl = fig[r, c] = CairoMakie.GridLayout()
        opts = Base.structdiff(panel, NamedTuple{(:values, :zones)})
        _zone_map_panel!(
            gl, geo, outline, panel.values, panel.zones;
            label_top, opts...
        )
    end
    isempty(title) || CairoMakie.Label(
        fig[0, 1:nc], title; fontsize = 16,
        font = :bold
    )
    return fig
end
