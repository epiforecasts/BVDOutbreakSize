```@raw html
---
aside: false
pageClass: bvd-dashboard
---
```

# Dashboard

The map shades each health zone or province by its reproduction number at the cut-off, its forecast over the coming week, its count to date, the chance its reproduction number exceeds one or its case fatality ratio.
The forecast, the count to date and the weekly chart show confirmed cases or confirmed deaths, chosen from a menu beside the views.
The health-zone model forecasts confirmed cases only, so health zones carry no forecast of confirmed deaths.
The case fatality ratio is confirmed deaths over confirmed cases to date, and provinces and the country also give the model's case fatality ratio.
The health-zone layer adds each zone's share of its patch and the chance of at least $K$ confirmed cases, with filters for zones with or without recent cases and for zones whose reproduction number is modelled separately.
The province layer takes the reproduction number, forecasts and modelled case fatality ratio from the joint model's patches, so the pooled provinces share one estimate and are hatched.
The treatment-centre view is by province only, since health zones report no treatment figures, and choosing it opens the province layer.
It shades each province by its patients in isolation, its beds or the share of its beds in use.
Each is shown as last reported by the province, as modelled at the cut-off or as forecast a week ahead, with the modelled figures for the pooled provinces again shared.
Hover over an area for its estimate, or click it for a side column with its key numbers, its reproduction number over time and its weekly counts against next week's forecast.
With nothing selected the column shows the national figures, and the table sorts the areas by any column.
The detail is on the [province estimates](estimates/province.md) and [health-zone estimates](estimates/zone.md) pages.

```@raw html
<p><a href="zone_map/index.html" target="_blank" rel="noopener">Open the map in a new tab</a>.</p>
<iframe class="bvd-dashboard-map" src="zone_map/index.html" title="Health-zone and province map" loading="lazy"></iframe>
```
