# Health-zone model, the second stage of a two-stage Markov melding. Stage
# one is the four-patch joint fit. Stage two runs the patch renewal over the
# 62 health zones on the parent's grid, with the parent's reproduction
# numbers, seed and growth rate melded through a Gaussian summary of its
# posterior, the parent's weekly patch infections as a melding term on the
# zone sums, and the per-zone confirmed-case tables as a within-patch
# composition. `zone_fit_inputs` builds everything from a parent chain,
# `zone_forward` is the pure forward pass the model and the render share,
# `bvd_zone` is the Turing model and `fit_zone` fits it.

## Truncation lags of the parent's delay PMFs, the defaults of
## `patch_infection_model` (generation interval, incubation) and
## `lab_delay_model` (receipt).
const ZONE_GI_NMAX = cdf_nmax(Gamma(2.71, 5.65))
const ZONE_INCUBATION_NMAX = cdf_nmax(lognormal_meansd(6.3, 3.5))
const ZONE_RECEIPT_NMAX = cdf_nmax(lognormal_meansd(4.5, 4.0))

## Chain keys of the parent quantities the zone stage reads. The joint
## attaches its latent submodel unprefixed, so the growth, generation
## interval and incubation blocks sit under their submodel names.
const _ZONE_PARENT_KEYS = (
    gi_alpha = Symbol("gi_state.α"),
    gi_theta = Symbol("gi_state.θ"),
    inc_mean = Symbol("inc_state.delay_mean"),
    inc_sd = Symbol("inc_state.delay_sd"),
    receipt_mean = Symbol("confirmed_state.receipt_state.d.delay_mean"),
    receipt_sd = Symbol("confirmed_state.receipt_state.d.delay_sd"),
    infections = :infections_patch,
    death_confirmation = :onset_to_death_confirmation_pmf,
    seed = Symbol("growth_state.C_T"),
    r = Symbol("growth_state.r"))

## --- Pure building blocks ------------------------------------------------

"""
$(TYPEDSIGNATURES)

Weekly deviation knots `(n_zones × n_knots)`. The first knot is the level
`σ_level (z_level − mean_p z_level)`, centred within each patch. Later knots
follow an AR(1) toward zero with retention `φ`: a walking zone adds
`σ_δ (z − mean_{W_p} z)` centred over the walking zones of its patch, a
level-only zone decays along `φ^{k−1} δ(1)`. `walk_index[z]` is the zone's
position among the `n_walking` walking zones (zero otherwise) and
`z_drift` holds the innovations knot by knot at
`(k − 2) n_walking + walk_index[z]`.
"""
function zone_deviation_knots(z_level::AbstractVector, z_drift::AbstractVector,
        σ_level::Real, σ_δ::Real, φ::Real,
        patch_ranges::AbstractVector{<:UnitRange},
        walking::AbstractVector{Bool}, walk_index::AbstractVector{<:Integer},
        n_walking::Integer, n_knots::Integer)
    Tp = promote_type(eltype(z_level), eltype(z_drift),
        typeof(float(σ_level)), typeof(float(σ_δ)), typeof(float(φ)))
    nz = length(z_level)
    δ = zeros(Tp, nz, n_knots)
    @inbounds for zs in patch_ranges
        isempty(zs) && continue
        m = zero(Tp)
        for z in zs
            m += z_level[z]
        end
        m /= length(zs)
        for z in zs
            δ[z, 1] = σ_level * (z_level[z] - m)
        end
    end
    @inbounds for k in 2:n_knots
        off = (k - 2) * n_walking
        for zs in patch_ranges
            isempty(zs) && continue
            nw = 0
            s = zero(Tp)
            for z in zs
                walking[z] || continue
                nw += 1
                s += z_drift[off + walk_index[z]]
            end
            bar = nw > 0 ? s / nw : zero(Tp)
            for z in zs
                innov = walking[z] ?
                        σ_δ * (z_drift[off + walk_index[z]] - bar) : zero(Tp)
                δ[z, k] = φ * δ[z, k - 1] + innov
            end
        end
    end
    return δ
end

"""
$(TYPEDSIGNATURES)

The linear map from knot values to the days `t0 … n`, as an `(n_days ×
n_knots)` weight matrix `W` with `δ_daily = W δ_knotsᵀ`. Each column is
[`interpolate_knots`](@ref) applied to a unit vector.
"""
function zone_interpolation_weights(knots::AbstractVector{<:Integer},
        t0::Integer, n::Integer)
    nb = length(knots)
    nd = n - t0 + 1
    W = zeros(Float64, nd, nb)
    unit = zeros(Float64, nb)
    for k in 1:nb
        fill!(unit, 0.0)
        unit[k] = 1.0
        daily = interpolate_knots(unit, knots, n)
        for j in 1:nd
            W[j, k] = daily[t0 + j - 1]
        end
    end
    return W
end

"""
$(TYPEDSIGNATURES)

The delay convolution over `nd` days as an `(nd × nd)` lower-triangular
Toeplitz matrix `F` with `F[j, j − s] = f_s` (lag 0 on the diagonal), so
the convolution of every zone's infections is one product `F Iᵀ`.
"""
function zone_delay_operator(f::AbstractVector, nd::Integer)
    F = zeros(Float64, nd, nd)
    for j in 1:nd, s in 0:min(length(f) - 1, j - 1)

        F[j, j - s] = f[s + 1]
    end
    return F
end

"""
$(TYPEDSIGNATURES)

Window indicator `(n × n_windows)` for the windows `(d_{v−1}, d_v]` of the
grid days `days`, the first opening on day one, so a daily series binned
into the windows is one product with this matrix.
"""
function zone_window_operator(days::AbstractVector{<:Integer}, n::Integer)
    nv = length(days)
    B = zeros(Float64, n, nv)
    for v in 1:nv
        lo = v == 1 ? 1 : Int(days[v - 1]) + 1
        for t in lo:min(Int(days[v]), n)
            B[t, v] = 1.0
        end
    end
    return B
end

"""
$(TYPEDSIGNATURES)

Patch membership indicator `(n_patches × n_zones)`.
"""
function zone_membership(patch_ranges::AbstractVector{<:UnitRange},
        nz::Integer)
    P = zeros(Float64, length(patch_ranges), nz)
    for (p, zs) in enumerate(patch_ranges), z in zs

        P[p, z] = 1.0
    end
    return P
end

"""
$(TYPEDSIGNATURES)

Log probability mass of counts `y` under a Dirichlet-multinomial with
concentration vector `α`, `Σ y` trials:

```math
\\log \\frac{N!\\,Γ(A)}{Γ(N + A)}
+ \\sum_i \\log \\frac{Γ(y_i + α_i)}{Γ(α_i)\\, y_i!},
\\qquad A = \\sum_i α_i.
```
"""
function dirichlet_multinomial_logpdf(y::AbstractVector{<:Integer},
        α::AbstractVector)
    Tp = float(eltype(α))
    N = 0
    A = zero(Tp)
    @inbounds for i in eachindex(y, α)
        N += y[i]
        A += α[i]
    end
    lp = loggamma(Tp(N + 1)) + loggamma(A) - loggamma(N + A)
    @inbounds for i in eachindex(y, α)
        lp += loggamma(y[i] + α[i]) - loggamma(α[i]) - loggamma(Tp(y[i] + 1))
    end
    return lp
end

"""
$(TYPEDSIGNATURES)

The composition log likelihood summed over the scored cells: for cell `c`,
patch `cell_patch[c]` at vintage `cell_vintage[c]` with allocated total
`cell_total[c]`, the observed zone counts in `counts` `(n_zones ×
n_vintages)` follow a Dirichlet-multinomial with concentration `κ π`,
where `π` is the modelled expected reports `C` of that patch's zones
normalised within the patch. `cell_const[c]` is the parameter-free part
of the mass, `log N! − Σ log y!`.
"""
function zone_composition_logpdf(counts::AbstractMatrix{<:Integer},
        C::AbstractMatrix, cell_patch::AbstractVector{<:Integer},
        cell_vintage::AbstractVector{<:Integer},
        cell_total::AbstractVector{<:Integer},
        cell_const::AbstractVector{<:Real},
        patch_ranges::AbstractVector{<:UnitRange}, κ::Real)
    Tp = promote_type(eltype(C), typeof(float(κ)))
    lp = zero(Tp)
    lgκ = loggamma(κ)
    @inbounds for c in eachindex(cell_patch)
        v = cell_vintage[c]
        zs = patch_ranges[cell_patch[c]]
        tot = zero(Tp)
        for z in zs
            tot += safe_rate(C[z, v])
        end
        acc = lgκ - loggamma(cell_total[c] + κ) + cell_const[c]
        for z in zs
            y = counts[z, v]
            y == 0 && continue
            α = κ * safe_rate(C[z, v]) / tot
            acc += loggamma(y + α) - loggamma(α)
        end
        lp += acc
    end
    return lp
end

## `log N! − Σ_z log y_z!` for the counts in one column of `counts` over the
## zones `zs`: the parameter-free part of the Dirichlet-multinomial mass.
function _zone_cell_const(counts::AbstractMatrix{<:Integer}, zs, v::Integer)
    N = 0
    acc = 0.0
    @inbounds for z in zs
        y = counts[z, v]
        N += y
        acc -= loggamma(y + 1.0)
    end
    return acc + loggamma(N + 1.0)
end

"""
$(TYPEDSIGNATURES)

The melding term: `log S ~ Normal(m, s)` summed over the kept cells, with
`S = weekly[cell_patch[c], cell_week[c]]` the zone stage's weekly patch
sum of infections and `m`, `s` the parent's posterior mean and sd of the
same log sum. The sum is floored at machine epsilon before the log.
"""
function zone_melding_logpdf(weekly::AbstractMatrix,
        cell_patch::AbstractVector{<:Integer},
        cell_week::AbstractVector{<:Integer},
        m::AbstractVector{<:Real}, s::AbstractVector{<:Real})
    Tp = float(eltype(weekly))
    lp = zero(Tp)
    @inbounds for c in eachindex(cell_patch)
        x = log(safe_rate(weekly[cell_patch[c], cell_week[c]]))
        z = (x - m[c]) / s[c]
        lp -= z * z / 2 + log(s[c]) + log(2π) / 2
    end
    return lp
end

"""
$(TYPEDSIGNATURES)

Daily log reproduction number of every patch `(n_patches × n)` from the
melded parent vector `φ`, whose first `n_patches × n_knots` entries are
the patch log reproduction numbers at the `knots`, flattened
column-major, interpolated linearly and held flat outside the knot span
([`interpolate_knots`](@ref)).
"""
function zone_parent_log_rt(φ::AbstractVector, knots::AbstractVector{<:Integer},
        np::Integer, n::Integer)
    K = length(knots)
    Tp = eltype(φ)
    out = zeros(Tp, np, n)
    @inbounds for p in 1:np
        daily = interpolate_knots(view(φ, p:np:(np * K)), knots, n)
        for t in 1:n
            out[p, t] = daily[t]
        end
    end
    return out
end

"""
$(TYPEDSIGNATURES)

Zone reproduction numbers `(n_zones × n)`, `R_z(t) = R_{p(z)}(t)
exp(δ_z(t))`, from the patch log reproduction numbers `(n_patches × n)`
and the daily deviations `(n × n_zones)`.
"""
function zone_rt_matrix(log_Rp::AbstractMatrix, δ_daily::AbstractMatrix,
        patch_of_zone::AbstractVector{<:Integer})
    n, nz = size(δ_daily)
    Tp = promote_type(eltype(log_Rp), eltype(δ_daily))
    R = zeros(Tp, nz, n)
    @inbounds for z in 1:nz
        p = patch_of_zone[z]
        for t in 1:n
            R[z, t] = exp(log_Rp[p, t] + δ_daily[t, z])
        end
    end
    return R
end

"""
$(TYPEDSIGNATURES)

Seed share of every zone: a within-patch softmax of `2 (z − mean)` over the
zones of each seeded patch (`seed_zones` lists them and `seed_ranges`
their positions in `z_seed`, patch by patch), times that patch's share
`patch_share[p]` of the parent's seed. Zones outside the seeded patches
are zero.
"""
function zone_seed_shares(z_seed::AbstractVector, patch_share::AbstractVector,
        seed_zones::AbstractVector{<:Integer},
        seed_ranges::AbstractVector{<:UnitRange},
        seed_patches::AbstractVector{<:Integer}, nz::Integer;
        scale::Real = 2.0)
    Tp = promote_type(eltype(z_seed), eltype(patch_share), typeof(float(scale)))
    out = zeros(Tp, nz)
    @inbounds for (i, rs) in enumerate(seed_ranges)
        isempty(rs) && continue
        m = zero(Tp)
        for j in rs
            m += z_seed[j]
        end
        m /= length(rs)
        tot = zero(Tp)
        for j in rs
            e = exp(scale * (z_seed[j] - m))
            out[seed_zones[j]] = e
            tot += e
        end
        share = patch_share[seed_patches[i]]
        for j in rs
            out[seed_zones[j]] = share * out[seed_zones[j]] / tot
        end
    end
    return out
end

"""
$(TYPEDSIGNATURES)

Seed matrix `(n_zones × L)`: each zone's share of the parent's cryptic
curve [`seed_infections`](@ref)`(seed, r, L)`.
"""
function zone_seed_matrix(seed_shares::AbstractVector, seed::Real, r::Real,
        L::Integer)
    curve = seed_infections(seed, r, L)
    nz = length(seed_shares)
    Tp = promote_type(eltype(seed_shares), eltype(curve))
    S = zeros(Tp, nz, L)
    @inbounds for z in 1:nz
        s = seed_shares[z]
        for j in 1:L
            S[z, j] = s * curve[j]
        end
    end
    return S
end

"""
$(TYPEDSIGNATURES)

Importation intensity per origin zone and day `(n_zones × n)`,
`min(level_z exp(effect ramp_t), 1)`, the parent's construction.
"""
function zone_mixing_intensity(level::AbstractVector, effect::Real,
        ramp::AbstractVector)
    nz = length(level)
    n = length(ramp)
    Tp = promote_type(eltype(level), typeof(float(effect)), eltype(ramp))
    ε = zeros(Tp, nz, n)
    @inbounds for t in 1:n
        e = exp(effect * ramp[t])
        for z in 1:nz
            ε[z, t] = min(level[z] * e, one(Tp))
        end
    end
    return ε
end

## Per-origin intensity levels `ε_bar exp(σ_ε (z_ε − mean z_ε))`.
function _zone_mixing_levels(ε_bar::Real, σ_ε::Real, z_ε::AbstractVector)
    Tp = promote_type(typeof(float(ε_bar)), typeof(float(σ_ε)), eltype(z_ε))
    nz = length(z_ε)
    m = zero(Tp)
    @inbounds for z in 1:nz
        m += z_ε[z]
    end
    m /= nz
    out = zeros(Tp, nz)
    @inbounds for z in 1:nz
        out[z] = ε_bar * exp(σ_ε * (z_ε[z] - m))
    end
    return out
end

## Patch seed shares under the parent's uncoupled rule: the primary takes
## `1 / (1 + Σ f)` and secondary patch `p` takes `f_{p−1} / (1 + Σ f)`.
function _zone_patch_seed_shares(seed_fraction::AbstractVector, np::Integer)
    Tp = float(eltype(seed_fraction))
    out = zeros(Tp, np)
    if isempty(seed_fraction)
        out[1] = one(Tp)
        return out
    end
    denom = one(Tp) + sum(seed_fraction)
    out[1] = one(Tp) / denom
    @inbounds for p in 2:np
        out[p] = seed_fraction[p - 1] / denom
    end
    return out
end

"""
$(TYPEDSIGNATURES)

Share of each zone in its patch's infections on day `t`, the patch sum
floored at `floatmin`.
"""
function zone_patch_shares(infections::AbstractMatrix,
        patch_ranges::AbstractVector{<:UnitRange}, t::Integer)
    Tp = eltype(infections)
    nz = size(infections, 1)
    out = zeros(Tp, nz)
    @inbounds for zs in patch_ranges
        tot = zero(Tp)
        for z in zs
            tot += infections[z, t]
        end
        tot = max(tot, floatmin(Tp))
        for z in zs
            out[z] = infections[z, t] / tot
        end
    end
    return out
end

"""
$(TYPEDSIGNATURES)

One forward pass of the zone model from its fixed data `zd` (the
`model_data` of [`zone_fit_inputs`](@ref)): the melded parent vector `φ`,
the deviation knots `(n_zones × n_knots)`, the seed shares, and the
importation levels per origin zone with the ramp effect (`nothing` and
zero when mixing is off). Builds the zone reproduction numbers, runs
[`patch_infections`](@ref) over the zones, convolves the reporting delay
and bins the vintage windows, and sums the patches by week. Called by the
model and, per draw, by the render. Returns `(; Rt, infections,
importation, reports, increments, patch_sums, weekly)`.
"""
function zone_forward(zd, φ::AbstractVector, δ_knots::AbstractMatrix,
        seed_shares::AbstractVector,
        ε_level::Union{Nothing, AbstractVector}, ε_effect::Real)
    np = length(zd.patch_ranges)
    K = length(zd.knots)
    log_Rp = zone_parent_log_rt(φ, zd.knots, np, zd.n)
    δ_daily = zd.interp * transpose(δ_knots)
    Rt = zone_rt_matrix(log_Rp, δ_daily, zd.patch_of_zone)
    seeds = zone_seed_matrix(seed_shares, exp(φ[np * K + 1]), φ[np * K + 2],
        zd.L)
    st = if ε_level === nothing
        patch_infections(Rt, zd.g, seeds, zd.kernel, 0.0)
    else
        patch_infections(Rt, zd.g, seeds, zd.kernel,
            zone_mixing_intensity(ε_level, ε_effect, zd.ramp))
    end
    I = st.infections
    reports = zd.report_matrix * transpose(I)
    increments = transpose(reports) * zd.bin_matrix
    patch_sums = zd.membership * I
    weekly = patch_sums * zd.week_matrix
    return (; Rt, infections = I, importation = st.importation, reports,
        increments, patch_sums, weekly)
end

## Dirichlet-multinomial concentration `(1 − ρ) / ρ`, with `ρ` and `1 − ρ`
## floored at machine epsilon so a proposal at either end of `[0, 1]` gives
## a finite, positive `κ`.
_zone_kappa(ρ::Real) = safe_rate(one(ρ) - ρ) / safe_rate(ρ)

## --- The model -----------------------------------------------------------

"""
Health-zone model, stage two of the melding. `zd` is the `model_data` of
[`zone_fit_inputs`](@ref): the Gaussian summary of the parent's melded
vector, the parent's weekly patch sums, the delay PMFs, the zone kernel,
the count matrix and its scored cells, the grid and the walking mask, all
plain arrays.

### Parameters

```math
\\begin{aligned}
η &\\sim N(0, I), & φ &= \\hat μ + \\hat L η \\\\
σ_L &\\sim N^+(0, 0.3),\\; z^L_z \\sim N(0, 1), &
δ_z(1) &= σ_L (z^L_z − \\bar z^L_p) \\\\
h &\\sim \\mathrm{LogNormal}(\\log 42, 0.6), & φ_δ &= 2^{−w/h} \\\\
σ_δ &\\sim N^+(0, 0.1),\\; z^δ_{z,k} \\sim N(0, 1), &
δ_z(k) &= φ_δ\\, δ_z(k−1) + σ_δ (z^δ_{z,k} − \\bar z^δ_{W_p,k}) \\\\
z^s_z &\\sim N(0, 1), & s_z &= \\mathrm{softmax}_p(2 (z^s_z − \\bar z^s_p)) \\\\
ε̄ &\\sim \\mathrm{Beta}(1, 100),\\; σ_ε \\sim N^+(0, 0.5),\\;
z^ε_z \\sim N(0, 1),\\; β_ε \\sim N(0, 0.5) \\\\
ρ &\\sim N^+(0, 0.1) \\text{ on } [0, 1], & κ &= (1 − ρ)/ρ
\\end{aligned}
```

`φ` holds the patch log reproduction numbers at the parent's knots, the
log seed and the growth rate; `melding = :cut` fixes it at `\\hat μ`. The
innovations exist for the walking zones `W_p` only. With `mixing` the
seed is split over the primary patch's zones and every other zone starts
empty; without it each secondary patch takes a sampled fraction of the
seed, as the parent does against an all-zero kernel, and `ε ≡ 0`.

### Likelihood

[`zone_forward`](@ref) gives the zone infections. The observed zone
increments of each patch and vintage follow a Dirichlet-multinomial on the
allocated total with concentration `κ π` ([`zone_composition_logpdf`](@ref))
and the weekly patch sums of zone infections follow the parent's posterior
of the same sums ([`zone_melding_logpdf`](@ref)), one `@addlogprob!` each.
`deaths = true` adds one Dirichlet-multinomial on the cumulative allocated
zone deaths at the last vintage through the death-confirmation PMF.

### Deterministics

`delta_knots_zone` (flattened `(n_zones × n_knots)`), `delta_T_zone`,
`R_T_zone`, `share_T_zone`, `seed_share_zone`, `mixing_epsilon_zone`
(at the cut-off), `mixing_level_zone`, `mixing_effect_zone`,
`parent_logR_knots` (flattened `(n_patches × n_knots)`), `parent_C_T`,
`parent_r`, `patch_sum_T`, `region_sd_zone`, `region_drift_sd_zone`,
`region_halflife_zone` and `composition_rho_zone`. Daily trajectories are
rebuilt from these by [`zone_forward`](@ref).
"""
@model function bvd_zone(zd;
        mixing::Bool = true,
        deaths::Bool = false,
        melding::Symbol = :sampled,
        region_sd_prior = truncated(Normal(0, 0.3); lower = 0),
        region_drift_sd_prior = truncated(Normal(0, 0.1); lower = 0),
        region_halflife_prior = LogNormal(log(42), 0.6),
        rho_prior = truncated(Normal(0, 0.1); lower = 0, upper = 1),
        importation_epsilon_prior = Beta(1, 100),
        importation_sd_prior = truncated(Normal(0, 0.5); lower = 0),
        importation_effect_prior = Normal(0, 0.5),
        seed_fraction_prior = LogNormal(log(0.05), 1.0),
        offset_prior = Normal(0, 1))
    nz = size(zd.counts, 1)
    np = length(zd.patch_ranges)
    K = length(zd.knots)
    nφ = np * K
    φ = if melding === :sampled
        η ~ product_distribution(fill(offset_prior, length(zd.phi_mean)))
        zd.phi_mean .+ zd.phi_chol * η
    else
        zd.phi_mean
    end
    σ_level ~ region_sd_prior
    z_level ~ product_distribution(fill(offset_prior, nz))
    δ_halflife ~ region_halflife_prior
    σ_δ ~ region_drift_sd_prior
    ρ ~ rho_prior
    ## Sampled only when used, or they would be prior-only dimensions.
    n_drift = zd.n_walking * (K - 1)
    z_drift = if n_drift > 0
        z_drift_ ~ product_distribution(fill(offset_prior, n_drift))
        z_drift_
    else
        Float64[]
    end
    ## Seed layout. Coupled zones are seeded in the primary patch only;
    ## uncoupled ones take the parent's per-patch fractions.
    seed_zones = mixing ? zd.seed_zones_coupled : zd.seed_zones_uncoupled
    seed_ranges = mixing ? zd.seed_ranges_coupled : zd.seed_ranges_uncoupled
    seed_patches = mixing ? zd.seed_patches_coupled :
                   zd.seed_patches_uncoupled
    z_seed ~ product_distribution(fill(offset_prior, length(seed_zones)))
    patch_share = if mixing || np == 1
        _zone_patch_seed_shares(Float64[], np)
    else
        seed_fraction ~ product_distribution(
            fill(seed_fraction_prior, np - 1))
        _zone_patch_seed_shares(seed_fraction, np)
    end
    seed_shares = zone_seed_shares(z_seed, patch_share, seed_zones,
        seed_ranges, seed_patches, nz)
    ε_level = if mixing
        ε_bar ~ importation_epsilon_prior
        σ_ε ~ importation_sd_prior
        z_ε ~ product_distribution(fill(offset_prior, nz))
        _zone_mixing_levels(ε_bar, σ_ε, z_ε)
    else
        nothing
    end
    ε_effect = if mixing
        β_ε ~ importation_effect_prior
        β_ε
    else
        0.0
    end
    φ_δ = exp2(-zd.week / δ_halflife)
    δ_knots = zone_deviation_knots(z_level, z_drift, σ_level, σ_δ, φ_δ,
        zd.patch_ranges, zd.walking, zd.walk_index, zd.n_walking, K)
    fw = zone_forward(zd, φ, δ_knots, seed_shares, ε_level, ε_effect)
    κ = _zone_kappa(ρ)
    @addlogprob! zone_composition_logpdf(zd.counts, fw.increments,
        zd.cell_patch, zd.cell_vintage, zd.cell_total, zd.cell_const,
        zd.patch_ranges, κ)
    @addlogprob! zone_melding_logpdf(fw.weekly, zd.meld_patch, zd.meld_week,
        zd.meld_mean, zd.meld_sd)
    if deaths
        D = fw.infections * zd.death_weights
        @addlogprob! zone_composition_logpdf(zd.death_counts,
            reshape(D, nz, 1), zd.death_cell_patch, zd.death_cell_vintage,
            zd.death_cell_total, zd.death_cell_const, zd.patch_ranges, κ)
    end
    n = zd.n
    R_T_zone := fw.Rt[:, n]
    share_T_zone := zone_patch_shares(fw.infections, zd.patch_ranges, n)
    patch_sum_T := fw.patch_sums[:, n]
    delta_knots_zone := vec(δ_knots)
    delta_T_zone := δ_knots[:, K]
    seed_share_zone := seed_shares
    parent_logR_knots := φ[1:nφ]
    parent_C_T := exp(φ[nφ + 1])
    parent_r := φ[nφ + 2]
    region_sd_zone := σ_level
    region_drift_sd_zone := σ_δ
    region_halflife_zone := δ_halflife
    composition_rho_zone := ρ
    if mixing
        mixing_level_zone := ε_level
        mixing_effect_zone := ε_effect
        mixing_epsilon_zone := zone_mixing_intensity(ε_level, ε_effect,
            zd.ramp)[:, n]
    end
    return (; fw.Rt, fw.infections, fw.increments, fw.weekly, δ_knots,
        seed_shares, φ)
end

## --- Parent summaries -----------------------------------------------------

## Draws of `key` from the parent chain (per-draw vectors with
## `vectors = true`), with a clear error naming the key when the chain does
## not carry it.
function _zone_parent_draws(chn, key::Symbol; vectors::Bool = false)
    _has_key(chn, key) || error(
        "zone_fit_inputs: the parent chain carries no `$(key)`; it must be " *
        "a `bvd_joint` chain sampled with the patch structure on.")
    return vectors ? _draw_vectors(chn, key) : _draws(chn, key)
end

## Generation-interval PMF (lag 1) of one parent draw, discretised as the
## parent did.
function _zone_gi_pmf(α::Real, θ::Real)
    pmf = discretise_censored(Gamma(α, θ), ZONE_GI_NMAX)
    return pmf[2:end] ./ sum(pmf[2:end])
end

## Mean over draws of a per-draw PMF.
function _zone_mean_pmf(build, ndraws::Integer)
    acc = build(1)
    for i in 2:ndraws
        acc = acc .+ build(i)
    end
    return acc ./ ndraws
end

"""
$(TYPEDSIGNATURES)

The parent's posterior-mean delay PMFs: the generation interval `g` (lag
1), the infection-to-confirmed-report PMF `f = incubation ⊛ receipt` (lag
0) and the infection-to-confirmed-death PMF `death_pmf` (empty when the
chain carries no death-confirmation delay). Returns `(; g, f, death_pmf)`.
"""
function zone_parent_inputs(chn)
    keys_ = _ZONE_PARENT_KEYS
    α = _zone_parent_draws(chn, keys_.gi_alpha)
    θ = _zone_parent_draws(chn, keys_.gi_theta)
    inc_m = _zone_parent_draws(chn, keys_.inc_mean)
    inc_s = _zone_parent_draws(chn, keys_.inc_sd)
    rec_m = _zone_parent_draws(chn, keys_.receipt_mean)
    rec_s = _zone_parent_draws(chn, keys_.receipt_sd)
    ndraws = length(α)
    g = _zone_mean_pmf(i -> _zone_gi_pmf(α[i], θ[i]), ndraws)
    inc_pmf(i) = discretise_censored(lognormal_meansd(inc_m[i], inc_s[i]),
        ZONE_INCUBATION_NMAX)
    rec_pmf(i) = discretise_censored(lognormal_meansd(rec_m[i], rec_s[i]),
        ZONE_RECEIPT_NMAX)
    f = _zone_mean_pmf(i -> convolve_pmf(inc_pmf(i), rec_pmf(i)), ndraws)
    death_pmf = if _has_key(chn, keys_.death_confirmation)
        dc = _draw_vectors(chn, keys_.death_confirmation)
        _zone_mean_pmf(i -> convolve_pmf(inc_pmf(i), Float64.(dc[i])), ndraws)
    else
        Float64[]
    end
    return (; g, f, death_pmf)
end

"""
$(TYPEDSIGNATURES)

Ledoit–Wolf shrinkage of the sample covariance of `X` (one draw per row)
toward the scaled identity `μ I`, `μ = tr(S) / p`. Returns the shrunk
covariance and the shrinkage intensity.
"""
function ledoit_wolf(X::AbstractMatrix)
    ndraw, p = size(X)
    Xc = X .- (sum(X; dims = 1) ./ ndraw)
    S = transpose(Xc) * Xc ./ ndraw
    μ = tr(S) / p
    d2 = sum(abs2, S - μ * I) / p
    b2 = 0.0
    for i in 1:ndraw
        x = view(Xc, i, :)
        b2 += sum(abs2, x * transpose(x) - S)
    end
    b2 /= ndraw^2 * p
    shrink = d2 > 0 ? clamp(b2 / d2, 0.0, 1.0) : 1.0
    Σ = shrink * μ * Matrix{Float64}(I, p, p) + (1 - shrink) * S
    return Σ, shrink
end

"""
$(TYPEDSIGNATURES)

Gaussian summary of the parent's melded vector `φ = (log R_p at the
knots, log seed, r)` and its weekly patch sums. The patch log reproduction
numbers come from [`reconstruct_patch_rt`](@ref) read at the knot days,
the seed and growth rate from the growth submodel's `C_T` and `r`. The
covariance is shrunk by [`ledoit_wolf`](@ref) and factored. The weekly
windows are `[1, knots[1]]` then `(knots[k−1], knots[k]]`; a window whose
mean patch sum is below `floor` infections is left out of the melding
cells. Returns `(; phi_mean, phi_cov, phi_chol, phi_sd, phi_labels,
meld_patch, meld_week, meld_mean, meld_sd, shrinkage)`.
"""
function zone_parent_summary(chn; n::Integer, np::Integer,
        knots::AbstractVector{<:Integer}, breakpoint, rt_start::Integer,
        rt_walk_start::Integer, week::Integer, ramp::Real,
        patch_names::AbstractVector, floor::Real = 1.0)
    keys_ = _ZONE_PARENT_KEYS
    rt = reconstruct_patch_rt(chn; n, breakpoint, n_patches = np, rt_start,
        rt_walk_start, week, ramp)
    seed = _zone_parent_draws(chn, keys_.seed)
    r = _zone_parent_draws(chn, keys_.r)
    ndraw = length(seed)
    K = length(knots)
    X = zeros(Float64, ndraw, np * K + 2)
    for i in 1:ndraw, p in 1:np, k in 1:K
        v = rt[p][i, knots[k]]
        ismissing(v) && error(
            "zone_parent_summary: the parent's reproduction number is " *
            "undefined on knot day $(knots[k]); the knots must sit inside " *
            "its established window.")
        X[i, p + np * (k - 1)] = log(safe_rate(v))
    end
    X[:, np * K + 1] .= log.(safe_rate.(seed))
    X[:, np * K + 2] .= r
    phi_mean = vec(sum(X; dims = 1)) ./ ndraw
    phi_cov, shrinkage = ledoit_wolf(X)
    phi_chol = Matrix(cholesky(Symmetric(phi_cov)).L)
    phi_sd = sqrt.(diag(phi_cov))
    phi_labels = vcat(
        ["$(patch_names[p]) knot $k" for k in 1:K for p in 1:np],
        ["log C_T", "r"])
    ## Weekly patch sums.
    infs = _zone_parent_draws(chn, keys_.infections; vectors = true)
    length(infs[1]) == np * n || error(
        "zone_fit_inputs: the parent's `infections_patch` holds " *
        "$(length(infs[1])) entries but $np patches by $n days is " *
        "$(np * n); the parent must be fitted to the same data cut-off.")
    W = zone_window_operator(knots, n)
    sums = zeros(Float64, ndraw, np, K)
    for i in 1:ndraw
        sums[i, :, :] = reshape(Float64.(infs[i]), np, n) * W
    end
    meld_patch = Int[]
    meld_week = Int[]
    meld_mean = Float64[]
    meld_sd = Float64[]
    for k in 1:K, p in 1:np

        s = view(sums, :, p, k)
        mean(s) >= floor || continue
        ls = log.(safe_rate.(s))
        push!(meld_patch, p)
        push!(meld_week, k)
        push!(meld_mean, mean(ls))
        push!(meld_sd, max(std(ls), 1e-3))
    end
    return (; phi_mean, phi_cov, phi_chol, phi_sd, phi_labels, meld_patch,
        meld_week, meld_mean, meld_sd, shrinkage)
end

## --- Fixed inputs from the parent chain ---------------------------------

## Cumulative count of `zone` in `prov` at the last vintage of `history`, or
## zero when the zone is absent from the block.
function _zone_last_cumulative(history, prov::AbstractString,
        zone::AbstractString)
    haskey(history, prov) || return 0
    haskey(history[prov], zone) || return 0
    c = history[prov][zone].counts
    return isempty(c) ? 0 : Int(c[end])
end

"""
$(TYPEDSIGNATURES)

Gravity kernel over the zones, [`province_importation_kernel`](@ref) on
the zone populations and haversine centroid distances, so every column's
off-diagonal sum is `1 − N_z / N` over all zones across patches.
"""
function zone_kernel(pops::AbstractVector, coords::AbstractVector)
    return province_importation_kernel(pops;
        distances = province_distance_matrix(coords))
end

## The package's health-zone metadata, or `nothing` when the file is absent.
function _default_health_zones()
    path = joinpath(@__DIR__, "..", "..", "data", "health_zones.csv")
    isfile(path) || return nothing
    return load_health_zones(path)
end

## Metadata row per zone; an error names the first zone without one.
function _zone_metadata_rows(zones, zone_province, zone_names)
    zones === nothing && error(
        "zone_fit_inputs: the health-zone metadata (data/health_zones.csv) " *
        "is needed for the zone kernel; pass `zones`.")
    lookup = Dict((r.province, r.zone) => r for r in zones)
    rows = [get(lookup, (zone_province[z], zone_names[z]), nothing)
            for z in eachindex(zone_names)]
    z = findfirst(isnothing, rows)
    z === nothing || error(
        "zone_fit_inputs: no population and centroid for zone " *
        "`$(zone_province[z]).$(zone_names[z])` in the health-zone metadata.")
    return rows
end

## The zone indices, ranges into the seed vector and patch ids of the
## seeded patches: the primary patch alone, or every non-empty patch.
function _zone_seed_layout(patch_ranges, patches)
    zones = Int[]
    ranges = UnitRange{Int}[]
    ids = Int[]
    for p in patches
        zs = patch_ranges[p]
        isempty(zs) && continue
        start = length(zones) + 1
        append!(zones, zs)
        push!(ranges, start:length(zones))
        push!(ids, p)
    end
    return (; zones, ranges, patches = ids)
end

"""
$(TYPEDSIGNATURES)

Everything the zone model and its render need, built once from the
parent chain and the observations.

Zones are the health-zone rows of `obs.zone_confirmed_history`
([`zone_increment_matrix`](@ref)), nested in the parent's patches of
`patch_names` and ordered patch by patch; `unallocated` pseudo-rows are
never units. Per vintage the zone increments are the consecutive-vintage
differences clamped at zero and the allocated patch total their sum; a
vintage on which a patch's unallocated count falls in either zone block
is zeroed for that patch and listed in `excluded`. Cells with a zero
allocated total are not scored. The walking set is the zones whose
cumulative confirmed count at the cut-off is at least `walk_threshold`,
in patches with at least two such zones.

The grid is the parent's: `n` days, the renewal start `rt_start` (the
seed length), the walk start `rt_walk_start` and the weekly knots from it
([`knot_days`](@ref)), the intervention `breakpoint` and its ramp. The
defaults derive them from `obs` as the joint does. The parent summary
([`zone_parent_summary`](@ref)) and delay PMFs
([`zone_parent_inputs`](@ref)) are read from `parent_chain`, and the zone
kernel ([`zone_kernel`](@ref)) from the `zones` metadata.

Returns a named tuple whose `model_data` field is the plain-array input
of [`bvd_zone`](@ref), alongside the zone keys and labels, the patch index
and labels, the vintage days and dates, the counts, the walking mask, the
grid, the parent summary and the kernel. `t0` is the walk start, the
first day the zone reproduction numbers are plotted from.
"""
function zone_fit_inputs(parent_chain, obs;
        walk_threshold::Integer = 30,
        week::Integer = 7,
        breakpoint = hasproperty(obs, :who_first_sitrep_days) ?
                     obs.n - obs.who_first_sitrep_days : missing,
        rt_start::Integer = hasproperty(obs, :tmrca_days) ?
                            clamp(
            obs.n - round(Int, obs.tmrca_days) +
            RENEWAL_START_LEAD, 1, obs.n) : 1,
        rt_walk_start::Integer = ismissing(breakpoint) ? rt_start :
                                 clamp(round(Int, breakpoint) - RT_WALK_LEAD,
            rt_start, obs.n),
        ramp::Real = RT_INTERVENTION_RAMP,
        melding_floor::Real = 1.0,
        zones = _default_health_zones(),
        patch_names::AbstractVector = PROVINCE_NAMES,
        patch_labels::AbstractVector = PROVINCE_LABELS)
    hasproperty(obs, :zone_confirmed_history) || error(
        "zone_fit_inputs: `obs` carries no `zone_confirmed_history`; load " *
        "the observations with a manifest that has the zone tables.")
    n = obs.n
    np = length(patch_names)
    ## Zone units and counts, patch by patch.
    death_history = hasproperty(obs, :zone_death_history) ?
                    obs.zone_death_history : Dict{String, Any}()
    reattribution = mergewith(vcat,
        zone_reattribution_days(obs.zone_confirmed_history),
        zone_reattribution_days(death_history))
    per_patch = zone_increment_matrix(obs.zone_confirmed_history, patch_names;
        reattribution)
    excluded = [(; patch = String(entry.patch),
                    date = obs.seeding + Day(d - 1))
                for entry in per_patch for d in entry.excluded]
    isempty(per_patch) && error(
        "zone_fit_inputs: the observations carry no zone histories.")
    zone_province = String[]
    zone_names = String[]
    patch_of_zone = Int[]
    patch_ranges = UnitRange{Int}[]
    days = Int[]
    blocks = Matrix{Int}[]
    for (p, entry) in enumerate(per_patch)
        start = length(zone_names) + 1
        for (prov, zone) in entry.zones
            push!(zone_province, prov)
            push!(zone_names, zone)
            push!(patch_of_zone, p)
        end
        push!(patch_ranges, start:length(zone_names))
        isempty(entry.zones) && continue
        if isempty(days)
            days = copy(entry.days)
        elseif entry.days != days
            error("zone_fit_inputs: patch `$(entry.patch)` is reported on " *
                  "different vintage days to the first patch with zones; " *
                  "the zone tables must share one `dates` array.")
        end
        push!(blocks, entry.increments)
    end
    nz = length(zone_names)
    nz > 0 || error("zone_fit_inputs: no health zones in the observations.")
    nv = length(days)
    counts = reduce(vcat, blocks)
    size(counts) == (nz, nv) || error(
        "zone_fit_inputs: the zone count blocks do not stack to " *
        "$(nz) zones by $(nv) vintages.")
    zone_keys = [zone_province[z] * "." * zone_names[z] for z in 1:nz]
    rows = _zone_metadata_rows(zones, zone_province, zone_names)
    labels = [r.label for r in rows]
    kernel = zone_kernel(Float64[r.population for r in rows],
        [(r.lat, r.lon) for r in rows])
    ## Scored cells: every (patch, vintage) with a positive allocated total.
    cell_patch = Int[]
    cell_vintage = Int[]
    cell_total = Int[]
    cell_const = Float64[]
    for v in 1:nv, (p, zs) in enumerate(patch_ranges)

        isempty(zs) && continue
        N = sum(@view counts[zs, v])
        N > 0 || continue
        push!(cell_patch, p)
        push!(cell_vintage, v)
        push!(cell_total, N)
        push!(cell_const, _zone_cell_const(counts, zs, v))
    end
    ## Cumulative confirmed at the cut-off, from the manifest rather than the
    ## clamped increments, and the walking set.
    cumulative = [_zone_last_cumulative(obs.zone_confirmed_history,
                      zone_province[z], zone_names[z]) for z in 1:nz]
    eligible = cumulative .>= walk_threshold
    walking = falses(nz)
    for zs in patch_ranges
        count(eligible[zs]) >= 2 || continue
        walking[zs] .= eligible[zs]
    end
    walk_index = zeros(Int, nz)
    n_walking = 0
    for z in 1:nz
        walking[z] || continue
        n_walking += 1
        walk_index[z] = n_walking
    end
    ## The parent's grid and summaries.
    knots = knot_days(n; week, start = rt_walk_start)
    parent = zone_parent_inputs(parent_chain)
    summary = zone_parent_summary(parent_chain; n, np, knots, breakpoint,
        rt_start, rt_walk_start, week, ramp, patch_names,
        floor = melding_floor)
    coupled = _zone_seed_layout(patch_ranges, 1:1)
    uncoupled = _zone_seed_layout(patch_ranges, 1:np)
    ## Optional death composition: cumulative allocated deaths at the last
    ## vintage, zones absent from the death table counting zero.
    death_counts = reshape(
        [_zone_last_cumulative(death_history, zone_province[z],
             zone_names[z]) for z in 1:nz], nz, 1)
    death_cell_patch = Int[]
    death_cell_total = Int[]
    death_cell_const = Float64[]
    for (p, zs) in enumerate(patch_ranges)
        isempty(zs) && continue
        N = sum(@view death_counts[zs, 1])
        N > 0 || continue
        push!(death_cell_patch, p)
        push!(death_cell_total, N)
        push!(death_cell_const, _zone_cell_const(death_counts, zs, 1))
    end
    death_pmf = isempty(parent.death_pmf) ? [1.0] : parent.death_pmf
    model_data = (; counts, cell_patch, cell_vintage, cell_total, cell_const,
        days, n, L = rt_start, knots, week, g = parent.g, f = parent.f,
        patch_ranges, patch_of_zone,
        walking = collect(walking), walk_index, n_walking,
        phi_mean = summary.phi_mean, phi_chol = summary.phi_chol,
        summary.meld_patch, summary.meld_week, summary.meld_mean,
        summary.meld_sd,
        ramp = sigmoid_ramp(n, breakpoint; ramp), kernel,
        interp = zone_interpolation_weights(knots, 1, n),
        report_matrix = zone_delay_operator(parent.f, n),
        bin_matrix = zone_window_operator(days, n),
        membership = zone_membership(patch_ranges, nz),
        week_matrix = zone_window_operator(knots, n),
        seed_zones_coupled = coupled.zones,
        seed_ranges_coupled = coupled.ranges,
        seed_patches_coupled = coupled.patches,
        seed_zones_uncoupled = uncoupled.zones,
        seed_ranges_uncoupled = uncoupled.ranges,
        seed_patches_uncoupled = uncoupled.patches,
        death_counts, death_cell_patch,
        death_cell_vintage = ones(Int, length(death_cell_patch)),
        death_cell_total, death_cell_const,
        death_weights = vec(sum(zone_delay_operator(death_pmf, n); dims = 1)))
    dates = [obs.seeding + Day(d - 1) for d in days]
    return (; model_data, zone_keys, zone_labels = labels, zone_province,
        zone_names, patch_of_zone, patch_ranges,
        patch_names = collect(String, patch_names),
        patch_labels = collect(String, patch_labels),
        days, dates, counts, cumulative, excluded, walking = collect(walking),
        walk_threshold, knots, t0 = rt_walk_start, n, week, breakpoint,
        rt_start, rt_walk_start, ramp, seeding = obs.seeding,
        cutoff = obs.cutoff, g = parent.g, f = parent.f, death_pmf, kernel,
        summary.phi_mean, summary.phi_sd, summary.phi_cov, summary.phi_labels,
        summary.shrinkage, summary.meld_patch, summary.meld_week,
        summary.meld_mean, summary.meld_sd)
end

## --- Fitting -------------------------------------------------------------

## Iterations-by-chains matrix of a sampler statistic, or `nothing`.
function _zone_stat(chn, name::Symbol)
    for e in FlexiChains.extras(chn)
        e.name === name || continue
        return Float64.(coalesce.(chn[e], NaN))
    end
    return nothing
end

"""
$(TYPEDSIGNATURES)

Fit the health-zone model [`bvd_zone`](@ref) melded from `parent_chain`, a
`bvd_joint` chain with the patch structure on, to the zone tables in
`obs`. Builds the inputs with [`zone_fit_inputs`](@ref) and runs
[`nuts_sample`](@ref) with its default initialisation at the given
settings. Logs the adapted step size, the fraction of iterations at the
tree-depth cap and the divergences per chain after. `mixing`, `deaths`
and `melding` select the variants; `walk_threshold`, `week`, `zones`,
`patch_names` and `patch_labels` pass through to
[`zone_fit_inputs`](@ref) and any other keyword to [`nuts_sample`](@ref).
Returns the chain.
"""
function fit_zone(parent_chain, obs;
        samples::Integer = 500,
        chains::Integer = 2,
        n_adapts::Integer = min(200, samples ÷ 2),
        target_accept::Real = 0.85,
        max_depth::Integer = 10,
        seed::Integer = 20260518,
        callback = nothing,
        mixing::Bool = true,
        deaths::Bool = false,
        melding::Symbol = :sampled,
        walk_threshold::Integer = 30,
        week::Integer = 7,
        zones = _default_health_zones(),
        patch_names::AbstractVector = PROVINCE_NAMES,
        patch_labels::AbstractVector = PROVINCE_LABELS,
        kwargs...)
    melding in (:sampled, :cut) || error(
        "fit_zone: `melding` must be :sampled or :cut, got $(repr(melding)).")
    inputs = zone_fit_inputs(parent_chain, obs; walk_threshold, week, zones,
        patch_names, patch_labels)
    zd = inputs.model_data
    deaths && isempty(zd.death_cell_patch) &&
        error(
            "fit_zone: deaths = true but the observations carry no allocated " *
            "zone deaths.")
    model = bvd_zone(zd; mixing, deaths, melding)
    n_zones = length(inputs.zone_keys)
    n_walking = count(inputs.walking)
    n_knots = length(inputs.knots)
    n_cells = length(zd.cell_patch)
    n_melding = length(zd.meld_patch)
    @info("fit_zone: starting", n_zones, n_walking, n_knots, n_cells,
        n_melding, mixing, deaths, melding)
    t = time()
    chn = nuts_sample(model; samples, chains, n_adapts, target_accept,
        max_depth, seed, callback, kwargs...)
    elapsed = time() - t
    steps = _zone_stat(chn, :step_size)
    depth = _zone_stat(chn, :tree_depth)
    div = _zone_stat(chn, :numerical_error)
    per_chain(m, f) = m === nothing ? missing :
                      [f(view(m, :, c)) for c in 1:size(m, 2)]
    minutes = round(elapsed / 60; digits = 1)
    step_size = per_chain(steps, x -> x[end])
    depth_cap_fraction = per_chain(depth, x -> mean(x .>= max_depth))
    divergences = per_chain(div, x -> Int(sum(x)))
    @info("fit_zone: finished", minutes, step_size, depth_cap_fraction,
        divergences)
    return chn
end
