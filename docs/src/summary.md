# Summary

A one-page overview of the headline results for readers with limited time.
Every number, table and figure on this page is produced by the same model
fit as the full [Analysis](@ref) and refreshes whenever the data updates.
See the [Analysis](@ref) page for the methods, assumptions and supporting
detail behind each result.

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

---

For the full results, methods and code see the [Analysis](@ref) page and the
[epiforecasts/BVDOutbreakSize](https://github.com/epiforecasts/BVDOutbreakSize)
repository.
