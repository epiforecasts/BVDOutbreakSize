# Future-expected-deaths counterfactual: a lower bound on future deaths
# under the assumption that every onward transmission stops at time `T`.
# The cohort already infected by `T` still contributes future expected
# deaths in the onset-to-death tail; per draw,
#
#     ΔD = CFR · ∫_0^T r · exp(r · s) · (1 - F_d(T - s)) ds,
#
# with `F_d` the Gamma(α, θ) CDF. Total projected cumulative deaths
# under the no-onward-transmission counterfactual is `obs_deaths + ΔD`.

function _committed_deaths_one(r, T, α, θ, CFR;
        onset_fraction = 1.0, alg = DEATH_INTEGRAL_ALG)
    delay_dist = Gamma(α, θ)
    scale = _delay_scale(delay_dist)
    g = let r = r, T = T, delay_dist = delay_dist
        s -> r * exp(r * s) * ccdf(delay_dist, T - s)
    end
    ## `onset_fraction` (the incubation mgf) maps the latent infection
    ## trajectory onto onsets before the onset-to-death survival tail,
    ## matching `deaths_model`.
    return onset_fraction * CFR * integrate(g, zero(T), T, scale; alg)
end

"""
Per-draw projection of cumulative deaths under the counterfactual
that every onward transmission stops at time `T`. Reads `:r, :T, :α,
:θ, :CFR` from the posterior `chn` and integrates

```math
\\Delta D = \\mathrm{CFR} \\cdot \\int_0^T r\\,\\exp(r\\,s)
            \\,\\bigl(1 - F_d(T - s)\\bigr)\\,ds,
```

with `F_d` the Gamma(α, θ) onset-to-death CDF, returning a
`DataFrame` with one row per draw:

- `:delta_deaths`     additional future expected deaths beyond `obs_deaths`
- `:total_projected`  `obs_deaths + delta_deaths`

`obs_deaths` is the number of deaths already observed at time `T`
(e.g. `obs.total_deaths` from the bundled observations). `alg` is the
quadrature scheme used for ΔD; defaults to `GaussLegendre(n = 64)`,
matching the rest of the package.
"""
function predict_no_onward_deaths(chn;
        obs_deaths::Real,
        alg = DEATH_INTEGRAL_ALG)
    r_draws = _draws(chn, :r)
    T_draws = _draws(chn, :T)
    α_draws = _draws(chn, :α)
    θ_draws = _draws(chn, :θ)
    CFR_draws = _draws(chn, :CFR)
    ## Incubation draws map the latent infection trajectory onto onsets;
    ## absent on chains predating the infection layer, where
    ## `onset_fraction = 1` recovers the previous behaviour.
    has_incubation = all(haskey_chain(chn, n) for n in (:α_inc, :θ_inc))
    α_inc = has_incubation ? _draws(chn, :α_inc) : nothing
    θ_inc = has_incubation ? _draws(chn, :θ_inc) : nothing

    n = length(r_draws)
    delta = Vector{Float64}(undef, n)
    @inbounds for i in 1:n
        os = has_incubation ?
             onset_rescale(Gamma(α_inc[i], θ_inc[i]), r_draws[i]) : 1.0
        delta[i] = _committed_deaths_one(r_draws[i], T_draws[i],
            α_draws[i], θ_draws[i],
            CFR_draws[i]; onset_fraction = os, alg)
    end
    total = float(obs_deaths) .+ delta
    return DataFrame(delta_deaths = delta, total_projected = total)
end
