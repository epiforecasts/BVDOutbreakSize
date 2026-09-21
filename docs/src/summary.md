# Summary dashboard

A one-page overview of the headline results for readers with limited time.
Every number, table and figure on this page is produced by the same model fit as the full [Analysis](analysis.md) and refreshes whenever the data updates.
See the [Analysis](analysis.md) page for the methods, assumptions and supporting detail behind each result, and the [Sensitivity](sensitivity.md) page for the forecast validation, the outbreak size implied by each data stream, the comparisons with McCabe et al. and Chamla et al., and the delay and tree-prior sensitivity analyses.
Which streams each vintage carries, and which are frozen, is recorded in the inclusion rules in `data/README.md`.

```@eval
using Markdown, BVDOutbreakSize
dir = joinpath(pkgdir(BVDOutbreakSize), "docs", "src", "summary_assets")
Markdown.parse("**Data as of:** " * read(joinpath(dir, "cutoff.md"), String))
```

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

The model runs one renewal equation per province and fits the national streams against the summed provinces, so the national count above is the sum of the three below.
Each cell is a median with a 90% credible interval.
The reproduction number and the relative ascertainment are read together, because the per-province case data identify only their product.
Modelled infections by province and the per-province parameter detail are in the [joint model estimates](analysis.md#Joint-model-estimates) and the [reproduction number over time](analysis.md#Reproduction-number-over-time) on the analysis page.

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
The [Sensitivity](sensitivity.md) page breaks these numbers down by parameter.

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

![Estimated reproduction number over time](summary_assets/rt.png)

The same trajectory by province, with the national one in grey behind each panel.
A panel tracking grey says that province moves with the national trend.

![Estimated reproduction number over time by province](summary_assets/rt_provinces.png)

## Infections over time

Modelled cumulative infections, symptom onsets and deaths.
These are the underlying outbreak, upstream of the testing and reporting that produce the observed counts, so they are larger than the reported cases.

![Estimated cumulative infections, onsets and deaths over time](summary_assets/infections.png)

## Reported cases: model versus observed

Modelled reported cases against the observed reported cases over time, a check that the fit reproduces what was seen on the ground.

![Modelled versus observed reported cases over time](summary_assets/reported_cases.png)

## Where each data stream points

The outbreak size and the reproduction number each data stream implies on its own, fitted to that stream alone, are on the sensitivity page: [outbreak size](sensitivity.md#Outbreak-size-estimated-by-each-data-stream) and [reproduction number](sensitivity.md#Reproduction-number-estimated-by-each-data-stream).
Agreement between the streams supports the joint estimate.
Disagreement shows where they pull in different directions.

---

For the full results, methods and code see the [Analysis](analysis.md) page and the [epiforecasts/BVDOutbreakSize](https://github.com/epiforecasts/BVDOutbreakSize) repository.
