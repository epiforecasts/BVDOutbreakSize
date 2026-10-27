```@raw html
---
aside: false
pageClass: bvd-dashboard
---
```

# Dashboard

The map shades each health zone or province by the estimate chosen above it.
Hover over an area for its value, or click it for its detail.
A date slider under the title steps the map back from the cut-off, and its play button moves it a week at a time.
The reproduction number, the counts and the reported treatment-centre figures follow it from their daily or weekly series.
The estimates without a series, such as the forecasts, show the latest saved cut-off on or before the date, or the current cut-off where no earlier one is saved.
The colour classes are the cut-off's at every date, so colours compare across dates.
The detail charts show the last 4 or 8 weeks before the cut-off or the whole series, and mark the map's date.
The compare menu sets an earlier data cut-off's estimates beside the current ones in the tooltips and the detail.
The [methods](methods.md) and the [province](estimates/province.md) and [health-zone](estimates/zone.md) estimates explain each figure.

```@raw html
<p><a href="zone_map/index.html" target="_blank" rel="noopener">Open the map in a new tab</a>.</p>
<iframe class="bvd-dashboard-map" src="zone_map/index.html" title="Health-zone and province map" loading="lazy"></iframe>
```
