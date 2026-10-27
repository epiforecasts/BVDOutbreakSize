```@raw html
---
aside: false
pageClass: bvd-dashboard
---
```

# Dashboard

The map shades each health zone or province by the estimate chosen above it.
Hover over an area for its value, or click it for its detail.
On the reproduction-number metrics a slider under the controls steps the map back through earlier dates, and its play button moves it a week at a time.
The chance that the reproduction number exceeds one is estimated at the cut-off only, so earlier dates show the median reproduction number.
The detail charts show the last 4 or 8 weeks before the cut-off or the whole series.
The compare menu sets an earlier data cut-off's estimates beside the current ones in the tooltips and the detail.
The [methods](methods.md) and the [province](estimates/province.md) and [health-zone](estimates/zone.md) estimates explain each figure.

```@raw html
<p><a href="zone_map/index.html" target="_blank" rel="noopener">Open the map in a new tab</a>.</p>
<iframe class="bvd-dashboard-map" src="zone_map/index.html" title="Health-zone and province map" loading="lazy"></iframe>
```
