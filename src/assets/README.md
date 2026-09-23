# Map assets

Files the report's maps read.
They sit outside `data/` because every file under `data/` is hashed into the fit cache key, and nothing here is a model input.

## `health_zones.geojson`

Simplified boundaries of every health zone in the seven affected provinces (Ituri, Nord-Kivu, Sud-Kivu, Haut-Uele, Tshopo, Bas-Uele and Sud-Ubangi), 167 zones.
Read by `load_health_zones_geojson()` and drawn by `plot_zone_map`, `plot_zone_map_panels` and `plot_province_map`.
Its properties are `zone` (the manifest zone key, empty for a zone with no confirmed case), `label`, `province` and `zscode`.

It is built by `scripts/build_health_zones.py`, which is on the `zone-mainstream` branch (#779), from the INRB-UMIE build of the DRC health-zone map (`build/drc_health_zones.geojson` in <https://github.com/INRB-UMIE/BDBV2026-Data>, build of 21 September 2026, commit `f489a92`).
That build uses the Ministry of Health `DRC_Health_zones` shapefile (Zones de santé, DSNIS) from the Humanitarian Data Exchange, <https://data.humdata.org/dataset/drc-health-data>.
Coordinates are rounded to four decimals and each ring is simplified on its own by Douglas-Peucker at 0.005 degrees, so a shared boundary can show a hairline gap.
The file is for display, not for topology.

The HDX dataset is published under the Open Database License (ODbL 1.0), <https://opendatacommons.org/licenses/odbl/1-0/>.
This file is a derived database of it and is made available under the same licence, with the attribution above.
The file carries none of the INRB-UMIE build's epidemiological or population values.
