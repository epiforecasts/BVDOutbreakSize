# Forecasts, scoring and counterfactuals

Projecting each stream forward from the posterior, scoring those forecasts against what was later observed, and the counterfactual and delay-corrected case-fatality calculations.
Scoring is continuous ranked probability score and its decomposition, with a persistence baseline for relative skill.

## Index

```@index
Pages = ["forecasts.md"]
```

## Reference

```@autodocs
Modules = [BVDOutbreakSize]
Pages = ["forecast.jl", "scoring.jl", "counterfactual.jl", "confirmed_cfr.jl"]
Private = false
```
