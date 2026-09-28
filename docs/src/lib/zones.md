# Health zones

The health-zone stage melded onto the joint fit: the model, the fixed inputs it reads from a parent chain, and the post-processing that turns its draws into shares, reproduction numbers, forecasts and scores.
The model and its inputs are in `models/zone.jl` and the post-processing of a zone chain is in `zone.jl`.

## Index

```@index
Pages = ["zones.md"]
```

## Reference

```@autodocs
Modules = [BVDOutbreakSize]
Pages = ["models/zone.jl", "zone.jl"]
Private = false
```
