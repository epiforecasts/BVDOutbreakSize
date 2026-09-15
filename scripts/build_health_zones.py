#!/usr/bin/env python3
"""Build data/health_zones.csv and data/health_zones.geojson.

Reads the INRB-UMIE build of the DRC health-zone map
(``build/drc_health_zones.geojson`` in
https://github.com/INRB-UMIE/BDBV2026-Data, the MoH ``DRC_Health_zones``
shapefile from the Humanitarian Data Exchange with WorldPop population
counts attached per zone) and the zone keys in the
``[zone_confirmed_history]`` block of data/observations.toml, and writes:

- ``data/health_zones.csv``: one row per zone key that appears in the
  block (the unallocated rows excluded) with its display label, province
  key, WorldPop population, polygon centroid and the shapefile's zscode.
- ``data/health_zones.geojson``: the map features of the seven affected
  provinces, geometry rounded to four decimals and simplified with
  Douglas-Peucker at 0.005 degrees, properties reduced to ``zone`` (our
  key, empty for a zone with no confirmed case), ``label``, ``province``
  and ``zscode``.

Every zone key must match a feature; the script exits non-zero listing
any that do not. Spelling differences between the SitRep names and the
shapefile's ``nom`` are resolved in ``ALIASES``. Rings are simplified
independently, so a shared boundary can open a hairline gap or overlap;
the file is for display, not for topology.

Needs only the Python standard library.

Usage:

    python3 scripts/build_health_zones.py path/to/drc_health_zones.geojson
"""

import csv
import json
import math
import re
import sys
import unicodedata
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
MANIFEST = ROOT / "data" / "observations.toml"
OUT_CSV = ROOT / "data" / "health_zones.csv"
OUT_GEOJSON = ROOT / "data" / "health_zones.geojson"

## Shapefile PROVINCE attribute -> manifest province key, in patch order.
PROVINCES = {
    "Ituri": "ituri",
    "Nord-Kivu": "nord_kivu",
    "Sud-Kivu": "sud_kivu",
    "Haut-Uele": "haut_uele",
    "Tshopo": "tshopo",
    "Bas-Uele": "bas_uele",
    "Sud-Ubangi": "sud_ubangi",
}

## Zone key -> folded shapefile name, where the two differ.
ALIASES = {
    "lubunga": "lubunga tshopo",
    "wanie_rukula": "wanierukula",
}

## Zone key -> display label, where the SitRep tables hyphenate a name the
## shapefile spaces. Every other label is the shapefile's own name.
LABELS = {
    "lubunga": "Lubunga",
    "makiso_kisangani": "Makiso-Kisangani",
    "miti_murhesa": "Miti-Murhesa",
    "nia_nia": "Nia-Nia",
    "wanie_rukula": "Wanie-Rukula",
}

DECIMALS = 4
TOLERANCE = 0.005


def fold(name):
    """Accents stripped, lower case, non-letters collapsed to one space."""
    s = unicodedata.normalize("NFKD", name)
    s = "".join(c for c in s if not unicodedata.combining(c)).lower()
    return re.sub(r"[^a-z]+", " ", s).strip()


def zone_keys(manifest):
    """(province key, zone key) pairs from the zone_confirmed_history block."""
    text = manifest.read_text(encoding="utf-8")
    block = re.search(r"^\[zone_confirmed_history\]\n(.*?)(?=^\[|\Z)", text,
        re.S | re.M)
    if block is None:
        sys.exit("no [zone_confirmed_history] block in " + str(manifest))
    keys = re.findall(r"^([a-z_]+)\.([a-z_]+) = ", block.group(1), re.M)
    return [(p, z) for p, z in keys if z != "unallocated"]


def ring_area_centroid(ring):
    """Signed area and centroid of a closed ring by the shoelace formula."""
    a = cx = cy = 0.0
    for (x0, y0), (x1, y1) in zip(ring, ring[1:]):
        w = x0 * y1 - x1 * y0
        a += w
        cx += (x0 + x1) * w
        cy += (y0 + y1) * w
    if a == 0:
        xs = [p[0] for p in ring]
        ys = [p[1] for p in ring]
        return 0.0, sum(xs) / len(xs), sum(ys) / len(ys)
    a *= 0.5
    return a, cx / (6 * a), cy / (6 * a)


def polygons(geometry):
    if geometry["type"] == "Polygon":
        return [geometry["coordinates"]]
    if geometry["type"] == "MultiPolygon":
        return geometry["coordinates"]
    raise ValueError(geometry["type"])


def centroid(geometry):
    """Centroid of the largest outer ring, as (lat, lon)."""
    best = None
    for poly in polygons(geometry):
        a, cx, cy = ring_area_centroid(poly[0])
        if best is None or abs(a) > best[0]:
            best = (abs(a), cx, cy)
    return best[2], best[1]


def perpendicular_distance(p, a, b):
    (px, py), (ax, ay), (bx, by) = p, a, b
    dx, dy = bx - ax, by - ay
    if dx == 0 and dy == 0:
        return math.hypot(px - ax, py - ay)
    t = ((px - ax) * dx + (py - ay) * dy) / (dx * dx + dy * dy)
    t = max(0.0, min(1.0, t))
    return math.hypot(px - (ax + t * dx), py - (ay + t * dy))


def douglas_peucker(points, tol):
    """Iterative Douglas-Peucker on an open polyline."""
    keep = [False] * len(points)
    keep[0] = keep[-1] = True
    stack = [(0, len(points) - 1)]
    while stack:
        i, j = stack.pop()
        if j <= i + 1:
            continue
        dmax, imax = 0.0, i
        for k in range(i + 1, j):
            d = perpendicular_distance(points[k], points[i], points[j])
            if d > dmax:
                dmax, imax = d, k
        if dmax > tol:
            keep[imax] = True
            stack.append((i, imax))
            stack.append((imax, j))
    return [p for p, k in zip(points, keep) if k]


def simplify_ring(ring, tol):
    """Simplify a closed ring, keeping it closed and at least a triangle."""
    pts = ring[:-1] if ring[0] == ring[-1] else ring
    if len(pts) < 4:
        return ring
    ## Split at the point farthest from the first so the closing edge is
    ## not a fixed endpoint of the whole ring.
    far = max(range(1, len(pts)),
        key=lambda k: math.hypot(pts[k][0] - pts[0][0], pts[k][1] - pts[0][1]))
    first = douglas_peucker(pts[: far + 1], tol)
    second = douglas_peucker(pts[far:] + [pts[0]], tol)
    out = first[:-1] + second[:-1]
    if len(out) < 3:
        return ring
    return out + [out[0]]


def round_ring(ring):
    out = []
    for x, y in ring:
        p = [round(x, DECIMALS), round(y, DECIMALS)]
        if not out or out[-1] != p:
            out.append(p)
    if out[0] != out[-1]:
        out.append(out[0])
    return out


def simplify(geometry):
    polys = []
    for poly in polygons(geometry):
        rings = []
        for ring in poly:
            r = round_ring(simplify_ring(ring, TOLERANCE))
            if len(r) >= 4:
                rings.append(r)
        if rings:
            polys.append(rings)
    if len(polys) == 1:
        return {"type": "Polygon", "coordinates": polys[0]}
    return {"type": "MultiPolygon", "coordinates": polys}


def main(path):
    with open(path, encoding="utf-8") as f:
        source = json.load(f)
    features = [f for f in source["features"]
                if f["properties"].get("province") in PROVINCES]
    by_name = {}
    for f in features:
        p = f["properties"]
        by_name[(PROVINCES[p["province"]], fold(p["nom"]))] = f

    keys = zone_keys(MANIFEST)
    matched = {}
    missing = []
    for prov, zone in keys:
        name = ALIASES.get(zone, zone.replace("_", " "))
        f = by_name.get((prov, name))
        if f is None:
            missing.append(f"{prov}.{zone}")
        else:
            matched[(prov, zone)] = f
    if missing:
        sys.exit("zone keys with no feature in the shapefile: "
                 + ", ".join(missing))
    order = {p: i for i, p in enumerate(PROVINCES.values())}
    keys.sort(key=lambda k: (order[k[0]], k[1]))

    with open(OUT_CSV, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(["zone", "label", "province", "population", "lat", "lon",
            "zscode"])
        for prov, zone in keys:
            p = matched[(prov, zone)]["properties"]
            lat, lon = centroid(matched[(prov, zone)]["geometry"])
            pop = p["worldpop"]["pop_count"]["pop_count"]
            w.writerow([zone, LABELS.get(zone, p["nom"]), prov,
                round(pop), f"{lat:.4f}", f"{lon:.4f}", p["zscode"]])

    key_of = {id(f): f"{prov}.{zone}" for (prov, zone), f in matched.items()}
    out = []
    for f in sorted(features, key=lambda f: (order[PROVINCES[f["properties"]["province"]]], f["properties"]["nom"])):
        p = f["properties"]
        zone = key_of.get(id(f), "")
        zone = zone.split(".", 1)[1] if zone else ""
        out.append({
            "type": "Feature",
            "properties": {
                "zone": zone,
                "label": LABELS.get(zone, p["nom"]) if zone else p["nom"],
                "province": PROVINCES[p["province"]],
                "zscode": p["zscode"],
            },
            "geometry": simplify(f["geometry"]),
        })
    with open(OUT_GEOJSON, "w", encoding="utf-8") as f:
        json.dump({"type": "FeatureCollection", "features": out}, f,
            ensure_ascii=False, separators=(",", ":"))
        f.write("\n")
    size = OUT_GEOJSON.stat().st_size
    print(f"{len(keys)} zones matched, {len(out)} features written, "
          f"{size / 1e6:.2f} MB")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    main(sys.argv[1])
