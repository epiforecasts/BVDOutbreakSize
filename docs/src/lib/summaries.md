# Summaries and diagnostics

Turning a fitted chain into the tables the report prints, and measuring how the sampler behaved.
The diagnostics cover R-hat, effective sample size and divergences, both overall and per parameter.
Point prevalence at the cut-off weights each draw's daily infections by how likely an infection of each age is still incubating or undetected.

## Index

```@index
Pages = ["summaries.md"]
```

## Reference

```@autodocs
Modules = [BVDOutbreakSize]
Pages = ["summaries.jl", "diagnostics.jl", "prevalence.jl"]
Private = false
```
