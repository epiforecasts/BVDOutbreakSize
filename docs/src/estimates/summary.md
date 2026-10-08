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

### Outbreak size and timing

```@eval
using Markdown, BVDOutbreakSize
dir = joinpath(pkgdir(BVDOutbreakSize), "docs", "src", "summary_assets")
Markdown.parse(read(joinpath(dir, "headline_counts.md"), String))
```

### Growth and severity

```@eval
using Markdown, BVDOutbreakSize
dir = joinpath(pkgdir(BVDOutbreakSize), "docs", "src", "summary_assets")
Markdown.parse(read(joinpath(dir, "headline_rates.md"), String))
```

Equal-tailed 30%, 60% and 90% credible intervals from the joint posterior.

### By province

```@eval
using Markdown, BVDOutbreakSize
dir = joinpath(pkgdir(BVDOutbreakSize), "docs", "src", "summary_assets")
Markdown.parse(read(joinpath(dir, "provinces_summary.md"), String))
```

Equal-tailed 90% credible intervals. Confirmed cases are the reported count over the latest week of situation reports.

![Modelled infections to date by province](../summary_assets/province_map.png)

Modelled infections to date by province, shaded by the posterior median on a log scale.

![Confirmed cases and reproduction number by province](../summary_assets/cases_rt_provinces.png)

### By health zone

```@eval
using Markdown, BVDOutbreakSize
dir = joinpath(pkgdir(BVDOutbreakSize), "docs", "src", "summary_assets")
Markdown.parse(read(joinpath(dir, "zone_headline.md"), String))
```

The zones with the most confirmed cases over the past two weeks, with equal-tailed 90% credible intervals.

![Confirmed cases and reproduction number by health zone](../summary_assets/cases_rt_zones.png)

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

## Reproduction number

![Estimated reproduction number over time](../summary_assets/rt.png)

The daily reproduction number nationally and by province, with 30%, 60% and 90% credible ribbons and the national trajectory in grey behind each province.

![Estimated reproduction number over time by province](../summary_assets/rt_provinces.png)

## Infections over time

![Estimated cumulative infections, onsets and deaths over time](../summary_assets/infections.png)

Modelled cumulative infections, symptom onsets and deaths, with 30%, 60% and 90% credible ribbons.
