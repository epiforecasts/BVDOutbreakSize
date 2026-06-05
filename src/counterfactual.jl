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

# Flat-infection cumulative trajectory for the counterfactual: exponential
# growth `exp(r·s)` up to the cut-off `T`, then held constant at `C(T)`
# afterwards (no onward transmission). The lab keeps draining the backlog
# accrued from the already-infected pool plus the continuing non-BVD
# background.
_flat_cumulative(r, T) =
    let r = r, T = T
        s -> exp(r * min(s, T))
    end

# Committed suspect-case backlog at horizon time `Th` under flat infections:
# the BVD onset-to-report pool (the flat infection trajectory convolved with
# the report delay, mapped to onsets by `onset_fraction`) plus the
# continuing constant-rate non-BVD background `λ_bg·Th`. Matches the suspect
# pool `N_susp` of `confirmed_cases_model` but with infections plateaued
# after `T`.
function _committed_nsusp_one(r, T, Th, α_rep, θ_rep, p_drc, λ_bg;
        onset_fraction = 1.0)
    f_rep = Gamma(α_rep, θ_rep)
    C = _flat_cumulative(r, T)
    bvd = onset_fraction * p_drc * delay_convolution(C, Th, f_rep)
    return max(bvd, zero(bvd)) + max(λ_bg * Th, zero(λ_bg)), max(bvd, zero(bvd))
end

# Committed laboratory-confirmed cases under the no-onward-transmission
# counterfactual: the whole suspect-case backlog that will ever be forwarded
# and received drains over the long horizon, so committed confirmed cases
# are the already-observed confirmed plus the forwarded-and-received backlog
# beyond what has already been analysed, weighted by the horizon positivity.
# The capacity limit is not applied: over a year the lab clears the finite
# committed backlog, so the committed total is capacity-independent. Returns
# `obs_confirmed + p_pos · max(received(Th) − obs_analysed, 0)`.
function _committed_confirmed_one(r, T, Th, α_rep, θ_rep, p_drc, λ_bg,
        τ_forward, α_recv, θ_recv, s_test, spec_test, obs_confirmed,
        obs_analysed; onset_fraction = 1.0)
    f_recv = Gamma(α_recv, θ_recv)
    C = _flat_cumulative(r, T)
    f_rep = Gamma(α_rep, θ_rep)
    ## Suspect-case pool as a function of time, BVD onsets (flat infections)
    ## plus continuing background, then receipt-delay convolved and forwarded.
    nsusp = let r = r, T = T, f_rep = f_rep, p_drc = p_drc, onset_fraction = onset_fraction,
        λ_bg = λ_bg, C = C

        u -> u <= zero(u) ? zero(u) :
             onset_fraction * p_drc * delay_convolution(C, u, f_rep) +
             max(λ_bg * u, zero(λ_bg))
    end
    received = τ_forward * max(delay_convolution(nsusp, Th, f_recv),
        zero(Th))
    ## Composition-linked positivity at the horizon: tested BVD share is the
    ## suspect-pool composition φ = μ_BVD / (μ_BVD + μ_bg).
    ntot, μ_bvd = _committed_nsusp_one(r, T, Th, α_rep, θ_rep, p_drc, λ_bg;
        onset_fraction)
    q = ntot > zero(ntot) ? clamp(μ_bvd / ntot, zero(ntot), one(ntot)) :
        zero(ntot)
    p_pos = clamp(s_test * q + (one(q) - spec_test) * (one(q) - q),
        zero(q), one(q))
    committed_backlog = max(received - obs_analysed, zero(received))
    return float(obs_confirmed) + p_pos * committed_backlog
end

# Committed laboratory-confirmed deaths under the counterfactual: the
# already-observed confirmed deaths plus the forwarded-positive increment of
# the suspect-death backlog drained over the long horizon, with infections
# held FLAT at `C(T)` after `T` (no onward transmission). Mirrors the
# suspect-death forward of `_forecast_confirmed_deaths_increment` but on the
# flat infection trajectory, so the BVD-death pool plateaus rather than
# growing exponentially over the year. The death pool is
# `os · p_deaths · CFR · ∫₀^u C_flat(s) f_death(u−s) ds` plus the continuing
# constant-rate background `λ_bg_death·u`; positivity uses the BVD share of
# the suspect-death pool at the horizon.
function _committed_confirmed_deaths_one(r, T, Th, α, θ, CFR, p_deaths,
        λ_bg_death, τ_death, s, spec, obs_confirmed_deaths;
        onset_fraction = 1.0, alg = DEATH_INTEGRAL_ALG)
    f_death = Gamma(α, θ)
    C = _flat_cumulative(r, T)
    bvd_death(u) = u <= zero(u) ? zero(u) :
                   onset_fraction * p_deaths * CFR *
                   delay_convolution(C, u, f_death)
    nsusp(u) = max(bvd_death(u), zero(u)) +
               max(λ_bg_death * u, zero(λ_bg_death))
    ΔN = max(nsusp(Th) - nsusp(T), zero(Th))
    denom = nsusp(Th)
    q_d = denom > zero(denom) ?
          clamp(max(bvd_death(Th), zero(Th)) / denom, zero(Th), one(Th)) :
          zero(Th)
    p_pos = clamp(s * q_d + (one(q_d) - spec) * (one(q_d) - q_d),
        zero(q_d), one(q_d))
    inc = max(τ_death * p_pos * ΔN, zero(ΔN))
    return float(obs_confirmed_deaths) + inc
end

"""
Per-draw committed totals under the counterfactual that every onward
transmission stops at the cut-off `T`: infections stay flat at the current
outbreak size `C(T)` and the already-infected pool continues to die and be
confirmed over the `horizon_days` horizon (default one year). Reads the
posterior `chn` and returns a `DataFrame` with one row per draw and
columns:

- `:infections` — committed infections, the current outbreak size
  `C(T) = 2^m` (`:cumulative_infections`); flat, since no new infections
  occur under the counterfactual.
- `:bvd_deaths` — committed true BVD deaths, `CFR · C(T)`. Every infection
  eventually contributes its death, so over a one-year horizon essentially
  all delayed onset-to-death events occur and the committed toll is the
  closed form `CFR · C(T)` (equivalently the realised deaths plus the
  committed-deaths tail of [`predict_no_onward_deaths`](@ref)).
- `:confirmed_cases` — committed laboratory-confirmed cases when the chain
  carries the lab parameters (`:s_test`, `:spec_test`, `:τ_forward`,
  `:α_recv`, `:θ_recv`) and `obs_confirmed` / `obs_analysed` are supplied.
  Computed by draining the committed suspect-case backlog (flat infections,
  continuing background) over the horizon: `obs_confirmed + p_pos ·
  max(received(T + horizon) − obs_analysed, 0)`. The capacity limit is not
  applied because over a year the lab clears the finite committed backlog.
- `:confirmed_deaths` — committed laboratory-confirmed deaths when the chain
  additionally carries `:τ_death`, `:p_deaths`, `:λ_bg_death` and
  `obs_confirmed_deaths` is supplied: the observed confirmed deaths plus the
  forwarded-positive suspect-death backlog drained over the horizon.

`horizon_days` sets how far the lab drains; the default `365` is long
enough that the committed lab outcomes have effectively saturated. `alg` is
the quadrature scheme, defaulting to `GaussLegendre(n = 64)`.
"""
function predict_committed(chn;
        obs_confirmed::Union{Real, Missing} = missing,
        obs_confirmed_deaths::Union{Real, Missing} = missing,
        obs_analysed::Union{Real, Missing} = missing,
        horizon_days::Real = 365,
        alg = DEATH_INTEGRAL_ALG)
    r = _draws(chn, :r)
    T = _draws(chn, :T)
    CFR = _draws(chn, :CFR)
    cumulative_infections = _draws(chn, :cumulative_infections)
    has_incubation = all(haskey_chain(chn, n) for n in (:α_inc, :θ_inc))
    α_inc = has_incubation ? _draws(chn, :α_inc) : nothing
    θ_inc = has_incubation ? _draws(chn, :θ_inc) : nothing

    ## Lab parameters for the committed confirmed-case backlog. The
    ## forwarding fraction is `τ_forward_out` on the production queue chain
    ## and `τ_forward` on prior/test chains; resolve whichever is present.
    fwd_key = _forward_key(chn)
    has_lab = fwd_key !== nothing &&
              all(haskey_chain(chn, n)
              for n in (:s_test, :spec_test, :α_recv, :θ_recv,
                  :α_rep, :θ_rep, :p_drc, :λ_bg)) &&
              obs_confirmed !== missing && obs_analysed !== missing
    s_test = has_lab ? _draws(chn, :s_test) : nothing
    spec_test = has_lab ? _draws(chn, :spec_test) : nothing
    τ_forward = has_lab ? _draws(chn, fwd_key) : nothing
    α_recv = has_lab ? _draws(chn, :α_recv) : nothing
    θ_recv = has_lab ? _draws(chn, :θ_recv) : nothing
    α_rep = has_lab ? _draws(chn, :α_rep) : nothing
    θ_rep = has_lab ? _draws(chn, :θ_rep) : nothing
    p_drc = has_lab ? _draws(chn, :p_drc) : nothing
    λ_bg = has_lab ? _draws(chn, :λ_bg) : nothing

    has_lab_deaths = has_lab &&
                     all(haskey_chain(chn, n)
                     for n in (:τ_death, :p_deaths, :λ_bg_death, :α, :θ)) &&
                     obs_confirmed_deaths !== missing
    α = has_lab_deaths ? _draws(chn, :α) : nothing
    θ = has_lab_deaths ? _draws(chn, :θ) : nothing
    τ_death = has_lab_deaths ? _draws(chn, :τ_death) : nothing
    p_deaths = has_lab_deaths ? _draws(chn, :p_deaths) : nothing
    λ_bg_death = has_lab_deaths ? _draws(chn, :λ_bg_death) : nothing

    n = length(r)
    infections = float.(cumulative_infections)
    ## Closed form: over a one-year horizon essentially all delayed deaths in
    ## the already-infected pool occur, so committed BVD deaths = CFR · C(T).
    bvd_deaths = CFR .* infections
    confirmed_cases = has_lab ? Vector{Float64}(undef, n) : nothing
    confirmed_deaths = has_lab_deaths ? Vector{Float64}(undef, n) : nothing

    @inbounds for i in 1:n
        Th = T[i] + horizon_days
        os = has_incubation ?
             onset_rescale(Gamma(α_inc[i], θ_inc[i]), r[i]) : 1.0
        if has_lab
            confirmed_cases[i] = _committed_confirmed_one(r[i], T[i], Th,
                α_rep[i], θ_rep[i], p_drc[i], λ_bg[i], τ_forward[i],
                α_recv[i], θ_recv[i], s_test[i], spec_test[i],
                obs_confirmed, obs_analysed; onset_fraction = os)
        end
        if has_lab_deaths
            confirmed_deaths[i] = _committed_confirmed_deaths_one(r[i], T[i],
                Th, α[i], θ[i], CFR[i], p_deaths[i], λ_bg_death[i],
                τ_death[i], s_test[i], spec_test[i], obs_confirmed_deaths;
                onset_fraction = os, alg)
        end
    end

    df = DataFrame(infections = infections, bvd_deaths = bvd_deaths)
    has_lab && (df.confirmed_cases = confirmed_cases)
    has_lab_deaths && (df.confirmed_deaths = confirmed_deaths)
    return df
end
