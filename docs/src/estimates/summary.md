# Summary dashboard

```@eval
using Markdown, BVDOutbreakSize, Dates
include(joinpath(pkgdir(BVDOutbreakSize), "docs", "front_matter.jl"))
dir = joinpath(pkgdir(BVDOutbreakSize), "docs", "src", "summary_assets")
cutoff = Date(strip(read(joinpath(dir, "cutoff.md"), String)))
Markdown.parse(report_dates(cutoff) * "\n\n" * readme_abstract())
```

This page summarises the headline results.
See the [National](national.md) and [Provinces](province.md) pages for the estimates at each level.
See [Forecasts](../forecasts/national.md) for the week ahead.
See [In-sample](../evaluation/insample.md) for how the model fits the data and [Forecast](../evaluation/forecast.md) evaluation for how past forecasts scored.
See [Methods](../methods.md) for the model, [Limitations](../limitations.md) for its caveats and [Sensitivity](../sensitivity.md) for the sensitivity analyses.

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

All intervals are equal-tailed 30%, 60% and 90% credible intervals from the joint posterior.

### By province

The model runs one renewal equation per province and fits the national streams against the summed provinces, so the national count above is the sum of the provinces.
Each range is an equal-tailed 90% credible interval.
The reproduction number and the relative ascertainment are read together, because the per-province case data identify only their product.
Each province's own estimates are in the [detail by province](province.md#Detail-by-province).
Modelled infections by province, the per-province parameter detail and the composition checks are on the same page.

```@eval
using Markdown, BVDOutbreakSize
dir = joinpath(pkgdir(BVDOutbreakSize), "docs", "src", "summary_assets")
Markdown.parse(read(joinpath(dir, "provinces.md"), String))
```

### Fit diagnostics

```@raw html
<details><summary>Expand: how the fits behind these numbers sampled</summary>
```

R-hat sets the spread within each chain against the spread across chains, and a value near one says the chains agree.
The bulk effective sample size is the number of independent draws the chains are worth, counted for the parameter where that count is lowest.
A divergent transition is a step the sampler could not take accurately.
The [Sensitivity](../sensitivity.md) page breaks these numbers down by parameter.

```@eval
using Markdown, BVDOutbreakSize
dir = joinpath(pkgdir(BVDOutbreakSize), "docs", "src", "summary_assets")
Markdown.parse(read(joinpath(dir, "diagnostics.md"), String))
```

```@raw html
</details>
```

## Estimated reproduction number

The time-varying reproduction number R(t), the average number of further infections caused by each infection.
A value above one means the outbreak is growing.

![Estimated reproduction number over time](../summary_assets/rt.png)

The same trajectory by province, with the national one in grey behind each panel.
A panel tracking grey says that province moves with the national trend.

![Estimated reproduction number over time by province](../summary_assets/rt_provinces.png)

## Infections over time

Modelled cumulative infections, symptom onsets and deaths.
These are the underlying outbreak, upstream of the testing and reporting that produce the observed counts, so they are larger than the reported cases.

![Estimated cumulative infections, onsets and deaths over time](../summary_assets/infections.png)

## Estimate variation by data stream

The [outbreak size](../sensitivity.md#Outbreak-size-estimated-by-each-data-stream) and the [reproduction number](../sensitivity.md#Reproduction-number-estimated-by-each-data-stream) each data stream implies on its own are on the sensitivity page.

---

For the full results, methods and code see the [National](national.md) page and the [epiforecasts/BVDOutbreakSize](https://github.com/epiforecasts/BVDOutbreakSize) repository.
