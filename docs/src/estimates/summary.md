# Summary dashboard

```@eval
using Markdown, BVDOutbreakSize, Dates
include(joinpath(pkgdir(BVDOutbreakSize), "docs", "front_matter.jl"))
dir = joinpath(pkgdir(BVDOutbreakSize), "docs", "src", "summary_assets")
cutoff = Date(strip(read(joinpath(dir, "cutoff.md"), String)))
Markdown.parse(report_dates(cutoff))
```

See [Methods](../methods.md) for the model and the definition of each quantity.

## Headline estimates

```@eval
using Markdown, BVDOutbreakSize
dir = joinpath(pkgdir(BVDOutbreakSize), "docs", "src", "summary_assets")
Markdown.parse(read(joinpath(dir, "headline.md"), String))
```

```@raw html
<figure class="bvd-figure">
<figcaption>Outbreak size and timing, with 30%, 60% and 90% credible intervals.</figcaption>
```

```@eval
using Markdown, BVDOutbreakSize
dir = joinpath(pkgdir(BVDOutbreakSize), "docs", "src", "summary_assets")
Markdown.parse(read(joinpath(dir, "headline_counts.md"), String))
```

```@raw html
</figure>
<figure class="bvd-figure">
<figcaption>Growth and severity, with 30%, 60% and 90% credible intervals.</figcaption>
```

```@eval
using Markdown, BVDOutbreakSize
dir = joinpath(pkgdir(BVDOutbreakSize), "docs", "src", "summary_assets")
Markdown.parse(read(joinpath(dir, "headline_rates.md"), String))
```

```@raw html
</figure>
```

### By province

```@raw html
<figure class="bvd-figure">
<figcaption>Province estimates with 90% credible intervals.</figcaption>
```

```@eval
using Markdown, BVDOutbreakSize
dir = joinpath(pkgdir(BVDOutbreakSize), "docs", "src", "summary_assets")
Markdown.parse(read(joinpath(dir, "provinces_summary.md"), String))
```

```@raw html
</figure>
<figure class="bvd-figure">
```

![Map of modelled infections to date by province](../summary_assets/province_map.png)

```@raw html
<figcaption>Modelled infections to date, shaded by the posterior median on a log scale.</figcaption>
</figure>
```

![Confirmed cases and reproduction number by province](../summary_assets/cases_rt_provinces.png)

### By health zone

```@raw html
<figure class="bvd-figure">
```

![Health zones at the cut-off](../summary_assets/zone_rt_map.png)

```@raw html
<figcaption>Each zone's reproduction number at the cut-off, its confirmed cases over the coming week (90% predictive interval) and its confirmed cases to date.
A reproduction number whose 90% interval spans one is washed towards white.</figcaption>
</figure>
```

```@raw html
<figure class="bvd-figure">
<figcaption>Zones with the most confirmed cases over the past two weeks, with 90% credible and predictive intervals.</figcaption>
```

```@eval
using Markdown, BVDOutbreakSize
dir = joinpath(pkgdir(BVDOutbreakSize), "docs", "src", "summary_assets")
Markdown.parse(read(joinpath(dir, "zone_headline.md"), String))
```

```@raw html
</figure>
<figure class="bvd-figure">
```

![Confirmed cases and reproduction number by health zone](../summary_assets/cases_rt_zones.png)

```@raw html
<figcaption>The twelve zones with the most confirmed cases over the past two weeks.</figcaption>
</figure>
```

### Health-zone forecast

```@eval
using Markdown, BVDOutbreakSize
dir = joinpath(pkgdir(BVDOutbreakSize), "docs", "src", "summary_assets")
Markdown.parse(read(joinpath(dir, "zone_forecast.md"), String))
```

### Fit diagnostics

```@eval
using Markdown, BVDOutbreakSize
dir = joinpath(pkgdir(BVDOutbreakSize), "docs", "src", "summary_assets")
Markdown.parse(read(joinpath(dir, "diagnostics_summary.md"), String))
```

```@raw html
<details><summary>Expand: how the fits behind these numbers sampled</summary>
```

```@eval
using Markdown, BVDOutbreakSize
dir = joinpath(pkgdir(BVDOutbreakSize), "docs", "src", "summary_assets")
Markdown.parse(read(joinpath(dir, "diagnostics.md"), String))
```

```@raw html
</details>
```

## Reproduction number trajectory

```@raw html
<figure class="bvd-figure">
```

![Estimated reproduction number over time](../summary_assets/rt.png)

```@raw html
<figcaption>National daily reproduction number with 30%, 60% and 90% credible ribbons.</figcaption>
</figure>
```

![Estimated reproduction number over time by province](../summary_assets/rt_provinces.png)

## Infections over time

```@raw html
<figure class="bvd-figure">
```

![Estimated cumulative infections, onsets and deaths over time](../summary_assets/infections.png)

```@raw html
<figcaption>Modelled cumulative infections, onsets and deaths with 30%, 60% and 90% credible ribbons.</figcaption>
</figure>
```
