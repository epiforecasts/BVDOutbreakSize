# Counterfactual: outbreak size with no onward transmission from the
# infections already seeded by the cut-off. With the renewal this is the
# expected eventual deaths from the current cohort: every infection
# present by the cut-off still dies with probability CFR, so the eventual
# deaths are `CFR · C_T`. The deaths already expected by the cut-off are
# `expected_deaths_T`, so the committed future deaths in the
# onset-to-death tail are `ΔD = CFR · C_T − expected_deaths_T`, and the
# total projected cumulative deaths under the counterfactual are
# `obs_deaths + ΔD`.

"""
Per-draw projection of cumulative deaths under the counterfactual that
every onward transmission stops at the cut-off. Reads `:CFR`, `:C_T` and
`:expected_deaths_T` from the posterior `chn` and forms the committed
future deaths

```math
\\Delta D = \\mathrm{CFR} \\cdot C_T - \\mathbb{E}[D_T],
```

with `C_T` the cumulative infections and `E[D_T]` the deaths already
expected by the cut-off, returning a `DataFrame` with one row per draw:

- `:delta_deaths`     additional future expected deaths beyond `obs_deaths`
- `:total_projected`  `obs_deaths + delta_deaths`

`obs_deaths` is the number of deaths already observed at the cut-off
(e.g. `obs.total_deaths` from the bundled observations).
"""
function predict_no_onward_deaths(chn; obs_deaths::Real)
    CFR = _draws(chn, :CFR)
    C_T = _draws(chn, :C_T)
    expected_deaths_T = _draws(chn, :expected_deaths_T)

    delta = max.(CFR .* C_T .- expected_deaths_T, 0.0)
    total = float(obs_deaths) .+ delta
    return DataFrame(delta_deaths = delta, total_projected = total)
end
