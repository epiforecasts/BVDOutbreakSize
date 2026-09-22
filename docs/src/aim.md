# Aim and origins

The aim is a transparent estimate of how far the outbreak has already grown, and of how much each published data stream contributes to that estimate.
Most infections are not yet reported, so the current size has to be inferred from the surveillance data that are available.
What the estimate can and cannot support is set out in the [limitations](limitations.md).

## Origins of this work

This work began as a replication of the [mccabe2026](@citet) report.
It has since evolved into a real-time joint Bayesian estimate of the current outbreak size.
The model is a discrete-time renewal process with a time-varying reproduction number, fitted to more of the available data streams than the original.
What is shared with that report and what has changed is set out component by component in the [comparison with published estimates](@ref "Comparison with published estimates").
The [comparison with McCabe et al.](@ref "Comparison with McCabe et al.") sets the current estimates against theirs.
