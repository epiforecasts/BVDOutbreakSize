# Observation models

One submodel per reported data stream, mapping the latent infections onto what surveillance recorded.
Each takes the latent trajectory and the stream's observations and adds that stream's likelihood.

## Index

```@index
Pages = ["observations.md"]
```

## Reference

```@autodocs
Modules = [BVDOutbreakSize]
Pages = ["models/observations.jl", "models/observation_distributions.jl"]
Private = false
```
