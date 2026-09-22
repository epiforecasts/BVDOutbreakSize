# Data and constants

Loading the situation-report data, freezing it to a past cut-off, and the fixed quantities the model treats as known.
`load_observations` reads `data/observations.toml` and is the entry point every fit starts from.
The constants carry the province geography, the travel volumes and the sampler settings.

## Index

```@index
Pages = ["data.md"]
```

## Reference

```@autodocs
Modules = [BVDOutbreakSize]
Pages = ["constants.jl", "data.jl", "onset_curve.jl"]
Private = false
```
