# Summary dashboard

A one-page overview of the headline results for readers with limited time.
Every number, table and figure on this page is produced by the same model
fit as the full [Analysis](analysis.md) and refreshes whenever the data updates.
See the [Analysis](analysis.md) page for the methods, assumptions and supporting
detail behind each result, and the [Sensitivity](sensitivity.md) page for the
forecast validation, the outbreak size implied by each data stream, the
comparisons with McCabe et al. and Chamla et al., and the delay and tree-prior sensitivity analyses.

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

All intervals are equal-tailed 30%, 60% and 90% credible intervals from the
joint posterior.

## Estimated reproduction number

The time-varying reproduction number R(t), the average number of further
infections caused by each infection.
A value above one means the outbreak is growing.

![Estimated reproduction number over time](summary_assets/rt.png)

## Infections over time

Modelled cumulative infections, symptom onsets and deaths.
These are the underlying outbreak, upstream of the testing and reporting that
produce the observed counts, so they are larger than the reported cases.

![Estimated cumulative infections, onsets and deaths over time](summary_assets/infections.png)

## Reported cases: model versus observed

Modelled reported cases against the observed reported cases over time, a check
that the fit reproduces what was seen on the ground.

![Modelled versus observed reported cases over time](summary_assets/reported_cases.png)

## Reproduction number by data stream

The reproduction number each data stream implies on its own, fitted to that
stream alone. Agreement between the streams supports the joint estimate;
disagreement shows where they pull in different directions.

![Reproduction number implied by each data stream](summary_assets/rt_streams.png)

---

For the full results, methods and code see the [Analysis](analysis.md) page and the
[epiforecasts/BVDOutbreakSize](https://github.com/epiforecasts/BVDOutbreakSize)
repository.
