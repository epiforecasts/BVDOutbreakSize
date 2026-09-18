# Spatial results

Where the outbreak sits, health zone by health zone.
These estimates come from the same fits as the [Summary](summary.md) and the full [Analysis](analysis.md), and refresh whenever the data updates.
See the [Analysis](analysis.md) page for the methods behind the health-zone model.

## Health zones

Each patch's infections split across its health zones, conditional on the province model.
The maps colour each affected zone by its reproduction number at the cut-off, the two bounds of the 90% interval on its forecast confirmed cases over the coming week, and its confirmed cases to date.
A zone whose reproduction-number interval straddles one is washed towards white, and a zone with too few infections for an estimate is grey.

![Health zones at the cut-off](summary_assets/zone_rt_map.png)

The one-week-ahead confirmed-case forecast for the fifteen zones with the largest forecasts, each the national forecast times the zone's projected share.

![One-week-ahead forecast by health zone](summary_assets/zone_forecast.png)

The zone ranking and the composition check are in the [health-zone estimates](analysis.md#Health-zone-estimates) of the analysis, and the validation against last week's observed split is in the [zone forecast validation](sensitivity.md#Forecast-by-health-zone) of the sensitivity page.
The same estimates are on the interactive map below.

## Health-zone map

Each affected health zone coloured by its current reproduction number, with the seven-day confirmed-case forecast, the confirmed cases to date and the zone's share of its patch's infections available from the switcher.
Hover over a zone for its estimate and 90% credible interval, click it for every number, or open the table view for a sortable list.
The map needs a browser, so it appears only on the documentation site.

```@raw html
<p><a href="zone_map/index.html" target="_blank" rel="noopener">Open the map in a new tab</a>.</p>
<iframe src="zone_map/index.html" title="Health-zone map" loading="lazy" style="width:100%;height:640px;border:1px solid var(--vp-c-divider);border-radius:8px;background:var(--vp-c-bg)"></iframe>
```

---

For the full results, methods and code see the [Analysis](analysis.md) page and the [epiforecasts/BVDOutbreakSize](https://github.com/epiforecasts/BVDOutbreakSize) repository.
