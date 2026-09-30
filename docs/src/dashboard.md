```@raw html
---
aside: false
pageClass: bvd-dashboard
---
```

# Dashboard

The map shades each health zone or province by its reproduction number at the cut-off, its confirmed cases forecast over the coming week, its confirmed cases to date or the chance its reproduction number exceeds one.
The health-zone layer adds each zone's share of its patch and the chance of at least $K$ cases, with filters for zones with or without recent cases and for zones whose reproduction number is modelled separately.
The province layer takes the reproduction number and forecast from the joint model's patches, so the pooled provinces share one estimate and are hatched.
Hover over an area for its estimate, or click it for a side column with its key numbers, its reproduction number over time and its weekly confirmed cases against next week's forecast.
With nothing selected the column shows the national figures, and the table sorts the areas by any column.
The detail is on the [province estimates](estimates/province.md) and [health-zone estimates](estimates/zone.md) pages.

```@raw html
<p><a href="zone_map/index.html" target="_blank" rel="noopener">Open the map in a new tab</a>.</p>
<iframe class="bvd-dashboard-map" src="zone_map/index.html" title="Health-zone and province map" loading="lazy"></iframe>
```
