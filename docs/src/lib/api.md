# [API overview](@id api-overview)

The code is a Julia package, `BVDOutbreakSize`, and the report pages are
literate scripts that call it.
This section documents what the package offers, grouped by the job each part
does.

## Public and internal

Everything the package exports is public and documented here.
The exported surface is what the report pages, the fit registry and the
scripts use, and it is what a reader should build on.
Anything not exported is an implementation detail, kept on the
[internals](internals.md) page, and may change without notice.

## The groups

The groups follow the order in which a fit runs.

| Group | What it covers |
|---|---|
| [Data and constants](data.md) | Loading and freezing the situation-report data, and the fixed quantities |
| [Renewal and delays](renewal.md) | The renewal recursion, the delay convolutions and the growth-rate conversions |
| [Priors and latent submodels](priors.md) | Turing submodels for the epidemiological and surveillance processes |
| [Observation models](observations.md) | One submodel per reported data stream |
| [Joint and single-stream models](joint.md) | The composers that assemble a fittable model |
| [Fitting](fitting.md) | Running NUTS, initialising it, and the progress callbacks |
| [Summaries and diagnostics](summaries.md) | Posterior tables and sampler diagnostics |
| [Forecasts, scoring and counterfactuals](forecasts.md) | Projecting forward, scoring against what arrived, and the counterfactuals |
| [Plotting](plotting.md) | Every figure in the report |

## Full index

```@index
Modules = [BVDOutbreakSize]
```
