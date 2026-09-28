#!/usr/bin/env python3
"""Print the population-weighted centre of each source province.

Reads the INRB-UMIE build of the DRC health-zone map
(``build/drc_health_zones.geojson`` in
https://github.com/INRB-UMIE/BDBV2026-Data, the MoH ``DRC_Health_zones``
shapefile from the Humanitarian Data Exchange with WorldPop population
counts attached per zone) and prints, for each province in
``PROVINCE_SOURCE_NAMES`` order, the WorldPop-weighted mean of its zones'
centroids as the ``PROVINCE_SOURCE_CENTRES`` literal in src/constants.jl.
A zone's centroid is the area-weighted centroid of its largest ring.

Needs only the Python standard library.

Usage:

    python3 scripts/build_province_centres.py path/to/drc_health_zones.geojson
"""

import json
import sys

# PROVINCE_SOURCE_NAMES order, keyed by the shapefile's province name.
PROVINCES = [
    ("ituri", "Ituri"),
    ("nord_kivu", "Nord-Kivu"),
    ("sud_kivu", "Sud-Kivu"),
    ("haut_uele", "Haut-Uele"),
    ("tshopo", "Tshopo"),
    ("bas_uele", "Bas-Uele"),
    ("sud_ubangi", "Sud-Ubangi"),
]


def ring_area_centroid(ring):
    """Signed area and centroid (lat, lon) of a closed lon/lat ring."""
    a = cx = cy = 0.0
    for (x0, y0), (x1, y1) in zip(ring, ring[1:]):
        c = x0 * y1 - x1 * y0
        a += c
        cx += (x0 + x1) * c
        cy += (y0 + y1) * c
    a /= 2
    return abs(a), cy / (6 * a), cx / (6 * a)


def zone_centroid(geometry):
    """Centroid (lat, lon) of a zone's largest outer ring."""
    if geometry["type"] == "Polygon":
        polygons = [geometry["coordinates"]]
    else:
        polygons = geometry["coordinates"]
    rings = [ring_area_centroid(p[0]) for p in polygons]
    _, lat, lon = max(rings)
    return lat, lon


def main(path):
    features = json.load(open(path))["features"]
    rows = []
    for key, name in PROVINCES:
        zones = [
            f for f in features if f["properties"]["province"] == name
        ]
        if not zones:
            sys.exit(f"no zones for {name}")
        total = lat = lon = 0.0
        for f in zones:
            pop = f["properties"]["worldpop"]["pop_count"]["pop_count"]
            zlat, zlon = zone_centroid(f["geometry"])
            total += pop
            lat += pop * zlat
            lon += pop * zlon
        rows.append((key, len(zones), total, lat / total, lon / total))
    print("const PROVINCE_SOURCE_CENTRES = [")
    for key, n, total, lat, lon in rows:
        print(
            f"    ({lat:.4f}, {lon:.4f}),  "
            f"# {key}, {n} zones, WorldPop {total / 1e6:.2f}M"
        )
    print("]")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    main(sys.argv[1])
