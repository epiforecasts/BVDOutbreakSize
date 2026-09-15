# Health-zone composition model, the second stage of a two-stage Markov
# melding. Stage one is the four-patch joint fit; stage two takes its
# posterior as fixed inputs (a cut, nothing feeds back) and fits only the
# per-zone confirmed-case tables, as a within-patch composition. The
# functions here build the fixed inputs from a parent chain
# (`zone_fit_inputs`), run the share renewal and its observation model as
# plain functions the model and the render both call, define the Turing
# model (`bvd_zone`) and fit it (`fit_zone`).

## Truncation lags of the stage-1 delay PMFs the zone model reuses. They
## mirror the defaults of `patch_infection_model` (generation interval,
## incubation) and `lab_delay_model` (receipt), so the posterior-mean PMFs
## built here are discretised exactly as the parent fit discretised them.
const ZONE_GI_NMAX = cdf_nmax(Gamma(2.71, 5.65))
const ZONE_INCUBATION_NMAX = cdf_nmax(lognormal_meansd(6.3, 3.5))
const ZONE_RECEIPT_NMAX = cdf_nmax(lognormal_meansd(4.5, 4.0))

## Chain keys of the parent quantities the cut reads. The joint attaches its
## latent submodel unprefixed, so the generation interval and incubation
## sit under their submodel names, and the receipt delay under the
## confirmed stream's.
const _ZONE_PARENT_KEYS = (
    gi_alpha = Symbol("gi_state.α"),
    gi_theta = Symbol("gi_state.θ"),
    inc_mean = Symbol("inc_state.delay_mean"),
    inc_sd = Symbol("inc_state.delay_sd"),
    receipt_mean = Symbol("confirmed_state.receipt_state.d.delay_mean"),
    receipt_sd = Symbol("confirmed_state.receipt_state.d.delay_sd"),
    infections = :infections_patch,
    death_confirmation = :onset_to_death_confirmation_pmf,
    C_T = :C_T)

## --- Pure building blocks ------------------------------------------------

"""
$(TYPEDSIGNATURES)

Initial zone shares at the grid start: a within-patch softmax of the
centred standard-normal draws `z_w` at a fixed `scale`,
`w_z = exp(scale (z_z − mean_p z)) / Σ_p`. The centring removes the flat
direction of the softmax, so the prior is proper on the shares. Returns a
vector over the zones, summing to one within every patch range.
"""
function zone_initial_shares(z_w::AbstractVector,
        patch_ranges::AbstractVector{<:UnitRange}, scale::Real)
    Tp = promote_type(eltype(z_w), typeof(float(scale)))
    w = zeros(Tp, length(z_w))
    @inbounds for zs in patch_ranges
        isempty(zs) && continue
        m = zero(Tp)
        for z in zs
            m += z_w[z]
        end
        m /= length(zs)
        tot = zero(Tp)
        for z in zs
            w[z] = exp(scale * (z_w[z] - m))
            tot += w[z]
        end
        for z in zs
            w[z] /= tot
        end
    end
    return w
end

"""
$(TYPEDSIGNATURES)

Weekly deviation knots `(n_zones × n_knots)` of the zone log-transmission
deviations. The first knot is the level `σ_level (z_level − mean_p z_level)`,
centred within each patch. Later knots follow an AR(1) toward zero with
retention `φ`: a walking zone adds an innovation
`σ_δ (z − mean_{W_p} z)` centred over the walking zones of its patch, and a
level-only zone decays along the AR mean path `φ^{k−1} δ(1)`. The patch
sum is therefore zero at every knot. `walk_index[z]` is the zone's position
among the `n_walking` walking zones (zero for a level-only zone) and
`z_drift` holds the innovations knot by knot, `(k − 2) n_walking +
walk_index[z]`.
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
                innov = walking[z] ? σ_δ * (z_drift[off + walk_index[z]] - bar) :
                        zero(Tp)
                δ[z, k] = φ * δ[z, k - 1] + innov
            end
        end
    end
    return δ
end

"""
$(TYPEDSIGNATURES)

Interpolate the deviation knots `(n_zones × n_knots)` placed on grid days
`knots` onto the days `t0 … n`, returning an `(n_days × n_zones)` matrix
with `n_days = n − t0 + 1`. Linear between knots and flat outside them, as
[`interpolate_knots`](@ref).
"""
function zone_interpolate_knots(δ_knots::AbstractMatrix,
        knots::AbstractVector{<:Integer}, t0::Integer, n::Integer)
    nz = size(δ_knots, 1)
    nd = n - t0 + 1
    out = zeros(eltype(δ_knots), nd, nz)
    @inbounds for z in 1:nz
        daily = interpolate_knots(view(δ_knots, z, :), knots, n)
        for j in 1:nd
            out[j, z] = daily[t0 + j - 1]
        end
    end
    return out
end

"""
$(TYPEDSIGNATURES)

The parts of the zone renewal and its delays that involve only patch
infections before the grid start `t0`, which are constant under the cut.
Every zone holds its initial share over those days, so each term enters a
zone's series multiplied by `w_z(t0)` only.

Returns, with `I_bar` `(n_patches × n)` and PMFs `g` (lag 1) and `f`
(lag 0):

- `force_pre[p, t] = Σ_{s ≥ 1, t − s < t0} g_s Ī_p(t − s)`, the pre-`t0`
  force of infection on every day;
- `report_pre[p, t] = Σ_{s ≥ 0, t − s < t0} f_s Ī_p(t − s)`, the
  pre-`t0` contribution to the daily expected reports;
- `report_pre_cum[p, t] = Σ_{t' ≤ min(t, t0 − 1)} report_pre[p, t']`, the
  expected reports accrued over the days before `t0`, so a vintage window
  opening before the grid start bins them by difference;
- `infections_pre[p] = Σ_{t < t0} Ī_p(t)`, the patch infections before
  the grid start.
"""
function zone_fixed_terms(I_bar::AbstractMatrix, g::AbstractVector,
        f::AbstractVector, t0::Integer)
    np, n = size(I_bar)
    Tp = promote_type(eltype(I_bar), eltype(g), eltype(f))
    force_pre = zeros(Tp, np, n)
    report_pre = zeros(Tp, np, n)
    report_pre_cum = zeros(Tp, np, n)
    infections_pre = zeros(Tp, np)
    @inbounds for p in 1:np
        for t in 1:n
            acc = zero(Tp)
            for s in 1:min(length(g), t - 1)
                t - s < t0 || continue
                acc += g[s] * I_bar[p, t - s]
            end
            force_pre[p, t] = acc
            acc = zero(Tp)
            for s in 0:min(length(f) - 1, t - 1)
                t - s < t0 || continue
                acc += f[s + 1] * I_bar[p, t - s]
            end
            report_pre[p, t] = acc
        end
        run = zero(Tp)
        for t in 1:n
            t < t0 && (run += report_pre[p, t])
            report_pre_cum[p, t] = run
        end
        for t in 1:min(t0 - 1, n)
            infections_pre[p] += I_bar[p, t]
        end
    end
    return (; force_pre, report_pre, report_pre_cum, infections_pre)
end

"""
$(TYPEDSIGNATURES)

The share renewal over the days `t0 … n`. Each zone's force of infection is
its own past infections through the generation interval `g` (lag 1), with
the days before `t0` entering as the initial share times the patch's
pre-`t0` force (`force_pre`, from [`zone_fixed_terms`](@ref)):

```math
Λ_z(t) = w_z(t_0)\\, Λ^{pre}_p(t) + \\sum_{s ≥ 1,\\; t − s ≥ t_0} g_s I_z(t − s),
\\qquad u_z(t) = e^{δ_z(t)} Λ_z(t),
\\qquad w_z(t) = u_z(t) / \\sum_{z' ∈ p} u_{z'}(t),
\\qquad I_z(t) = Ī_p(t)\\, w_z(t).
```

`δ_daily` is `(n_days × n_zones)` from [`zone_interpolate_knots`](@ref),
`w0` the initial shares and `I_bar` the fixed patch infections. Only the
share denominator is floored, so a patch with almost no infections keeps
its initial split rather than collapsing to uniform.

With `kernel` (an `(n_zones × n_zones)` matrix, column-stochastic within
each patch and zero across patches) and per-patch `ε`, a fraction of each
zone's force is redistributed within its patch before normalising,
`v = (1 − ε_p) u + ε_p K u`. Pass `kernel = nothing` to skip it.

Returns `(; shares, forces, infections)`, each `(n_days × n_zones)`.
"""
function zone_share_renewal(I_bar::AbstractMatrix, g::AbstractVector,
        δ_daily::AbstractMatrix, w0::AbstractVector,
        patch_ranges::AbstractVector{<:UnitRange}, t0::Integer,
        force_pre::AbstractMatrix;
        kernel::Union{Nothing, AbstractMatrix} = nothing,
        ε::Union{Nothing, AbstractVector} = nothing)
    nd, nz = size(δ_daily)
    Tp = promote_type(eltype(I_bar), eltype(g), eltype(δ_daily),
        eltype(w0), eltype(force_pre),
        ε === nothing ? Float64 : eltype(ε),
        kernel === nothing ? Float64 : eltype(kernel))
    I = zeros(Tp, nd, nz)
    Λ = zeros(Tp, nd, nz)
    W = zeros(Tp, nd, nz)
    u = zeros(Tp, nz)
    v = zeros(Tp, nz)
    L = length(g)
    floor_ = floatmin(Tp)
    @inbounds for j in 1:nd
        t = t0 + j - 1
        smax = min(L, j - 1)
        for p in eachindex(patch_ranges)
            zs = patch_ranges[p]
            isempty(zs) && continue
            fp = force_pre[p, t]
            for z in zs
                acc = w0[z] * fp
                for s in 1:smax
                    acc += g[s] * I[j - s, z]
                end
                Λ[j, z] = acc
                u[z] = exp(δ_daily[j, z]) * acc
            end
            if kernel === nothing || ε === nothing
                for z in zs
                    v[z] = u[z]
                end
            else
                e = ε[p]
                for z in zs
                    mixed = zero(Tp)
                    for q in zs
                        mixed += kernel[z, q] * u[q]
                    end
                    v[z] = (one(Tp) - e) * u[z] + e * mixed
                end
            end
            tot = zero(Tp)
            for z in zs
                tot += v[z]
            end
            tot = max(tot, floor_)
            Ip = I_bar[p, t]
            for z in zs
                W[j, z] = v[z] / tot
                I[j, z] = Ip * W[j, z]
            end
        end
    end
    return (; shares = W, forces = Λ, infections = I)
end

"""
$(TYPEDSIGNATURES)

Daily expected confirmed reports per zone over the days `t0 … n`, the zone
infections `(n_days × n_zones)` carried through the infection-to-report
PMF `f` (lag 0), plus the initial share times the patch's pre-`t0`
contribution (`report_pre`, from [`zone_fixed_terms`](@ref)). Returns an
`(n_days × n_zones)` matrix.
"""
function zone_reports(infections::AbstractMatrix, f::AbstractVector,
        w0::AbstractVector, patch_ranges::AbstractVector{<:UnitRange},
        t0::Integer, report_pre::AbstractMatrix)
    nd, nz = size(infections)
    Tp = promote_type(eltype(infections), eltype(f), eltype(w0),
        eltype(report_pre))
    r = zeros(Tp, nd, nz)
    F = length(f)
    @inbounds for j in 1:nd
        t = t0 + j - 1
        smax = min(F - 1, j - 1)
        for p in eachindex(patch_ranges)
            zs = patch_ranges[p]
            rp = report_pre[p, t]
            for z in zs
                acc = w0[z] * rp
                for s in 0:smax
                    acc += f[s + 1] * infections[j - s, z]
                end
                r[j, z] = acc
            end
        end
    end
    return r
end

"""
$(TYPEDSIGNATURES)

Bin the daily expected reports `(n_days × n_zones)` from
[`zone_reports`](@ref) into the vintage windows `(d_{v−1}, d_v]` given by
the grid days `days`, the first window opening on day one. Days before the
grid start `t0` contribute the initial share times the patch's pre-`t0`
accrual (`report_pre_cum`, from [`zone_fixed_terms`](@ref)). Returns the
`(n_zones × n_vintages)` expected increments.
"""
function zone_report_increments(reports::AbstractMatrix, w0::AbstractVector,
        patch_ranges::AbstractVector{<:UnitRange},
        days::AbstractVector{<:Integer}, t0::Integer,
        report_pre_cum::AbstractMatrix)
    nd, nz = size(reports)
    nv = length(days)
    Tp = promote_type(eltype(reports), eltype(w0), eltype(report_pre_cum))
    C = zeros(Tp, nz, nv)
    @inbounds for v in 1:nv
        lo = v == 1 ? 1 : Int(days[v - 1]) + 1
        hi = Int(days[v])
        for p in eachindex(patch_ranges)
            zs = patch_ranges[p]
            ## Accrued before the grid start, by difference of the cumulative.
            pre = zero(Tp)
            if lo < t0
                top = min(hi, t0 - 1)
                pre = report_pre_cum[p, top] -
                      (lo > 1 ? report_pre_cum[p, lo - 1] : zero(Tp))
            end
            jlo = max(lo, t0) - t0 + 1
            jhi = hi - t0 + 1
            for z in zs
                acc = w0[z] * pre
                for j in jlo:min(jhi, nd)
                    acc += reports[j, z]
                end
                C[z, v] = acc
            end
        end
    end
    return C
end

"""
$(TYPEDSIGNATURES)

Log probability mass of counts `y` under a Dirichlet-multinomial with
concentration vector `α`, `Σ y` trials:

```math
\\log \\frac{N!\\,Γ(A)}{Γ(N + A)} + \\sum_i \\log \\frac{Γ(y_i + α_i)}{Γ(α_i)\\, y_i!},
\\qquad A = \\sum_i α_i.
```

Plain arithmetic and `loggamma`, so it differentiates under Mooncake.
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
normalised within the patch. `cell_const[c]` is the count-only part of the
mass, `log N! − Σ log y!`, precomputed since it carries no parameter. A
zero count contributes nothing to the sum over zones, so those terms are
skipped.
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

One forward pass of the zone model from its fixed data `zd` (the
`model_data` of [`zone_fit_inputs`](@ref)) and a draw's deviation knots
`(n_zones × n_knots)`, initial shares `w0` and per-patch mixing fractions
`ε` (`nothing` when mixing is off): the daily deviations, the share
renewal and the binned expected reports. This is the function the model
evaluates and the render re-runs per draw, so the two never diverge.
Returns `(; shares, forces, infections, reports, increments)`.
"""
function zone_forward(zd, δ_knots::AbstractMatrix, w0::AbstractVector,
        ε::Union{Nothing, AbstractVector})
    δ_daily = zone_interpolate_knots(δ_knots, zd.knots, zd.t0, zd.n)
    return zone_forward_daily(zd, δ_daily, w0, ε)
end

"""
$(TYPEDSIGNATURES)

[`zone_forward`](@ref) from daily deviations `(n_days × n_zones)` rather
than knots, for the forecast projection, which continues the deviations
on their AR mean path day by day.
"""
function zone_forward_daily(zd, δ_daily::AbstractMatrix, w0::AbstractVector,
        ε::Union{Nothing, AbstractVector})
    kernel = ε === nothing ? nothing : zd.mixing_kernel
    st = zone_share_renewal(zd.I_bar, zd.g, δ_daily, w0, zd.patch_ranges,
        zd.t0, zd.force_pre; kernel, ε)
    reports = zone_reports(st.infections, zd.f, w0, zd.patch_ranges,
        zd.t0, zd.report_pre)
    increments = zone_report_increments(reports, w0, zd.patch_ranges,
        zd.days, zd.t0, zd.report_pre_cum)
    return (; st.shares, st.forces, st.infections, reports, increments)
end

## Cumulative zone infections to the last grid day: the pre-`t0` patch
## infections at the initial share plus the renewal's own days.
function _zone_cumulative_infections(infections::AbstractMatrix,
        w0::AbstractVector, patch_ranges::AbstractVector{<:UnitRange},
        infections_pre::AbstractVector)
    nd, nz = size(infections)
    Tp = promote_type(eltype(infections), eltype(w0), eltype(infections_pre))
    cum = zeros(Tp, nz)
    @inbounds for p in eachindex(patch_ranges)
        zs = patch_ranges[p]
        for z in zs
            acc = w0[z] * infections_pre[p]
            for j in 1:nd
                acc += infections[j, z]
            end
            cum[z] = acc
        end
    end
    return cum
end

## Implied zone reproduction number on grid row `j`, `I_z / Λ_z`, `NaN` where
## the zone's cumulative infections are below `floor` (the ratio of two
## near-zero numbers says nothing there).
function _zone_rt_at(infections::AbstractMatrix, forces::AbstractMatrix,
        cum::AbstractVector, j::Integer, floor::Real)
    Tp = promote_type(eltype(infections), eltype(forces))
    nz = size(infections, 2)
    out = zeros(Tp, nz)
    nan = convert(Tp, NaN)
    @inbounds for z in 1:nz
        r = infections[j, z] / max(forces[j, z], floatmin(Tp))
        out[z] = cum[z] >= floor ? r : nan
    end
    return out
end

## Shares at the knot days as an `(n_zones × n_knots)` matrix.
function _zone_shares_at_knots(shares::AbstractMatrix,
        knots::AbstractVector{<:Integer}, t0::Integer)
    nz = size(shares, 2)
    out = zeros(eltype(shares), nz, length(knots))
    @inbounds for (k, d) in enumerate(knots), z in 1:nz

        out[z, k] = shares[d - t0 + 1, z]
    end
    return out
end

## --- The model -----------------------------------------------------------

"""
Health-zone composition model, stage two of the melding. `zd` is the
`model_data` of [`zone_fit_inputs`](@ref): the fixed stage-1 inputs (patch
infections `Ī_p`, the generation-interval and infection-to-report PMFs),
the zone count matrix and the scored cells, the knot days and the walking
mask, all plain arrays. Nothing in it is sampled here.

### Parameters

```math
\\begin{aligned}
z^w_z &\\sim N(0, 1), & w_z(t_0) &= \\mathrm{softmax}_p\\bigl(2\\,(z^w_z − \\bar z^w_p)\\bigr) \\\\
σ_L &\\sim N^+(0, 0.3),\\; z^L_z \\sim N(0, 1), & δ_z(1) &= σ_L (z^L_z − \\bar z^L_p) \\\\
h &\\sim \\mathrm{LogNormal}(\\log 42, 0.6), & φ &= 2^{−7/h} \\\\
σ_δ &\\sim N^+(0, 0.1),\\; z^δ_{z,k} \\sim N(0, 1), & δ_z(k) &= φ\\, δ_z(k−1) + σ_δ (z^δ_{z,k} − \\bar z^δ_{W_p,k}) \\\\
ρ &\\sim N^+(0, 0.1) \\text{ on } [0, 1], & κ &= (1 − ρ)/ρ
\\end{aligned}
```

The innovations exist for the walking zones `W_p` only (cumulative
confirmed cases at the cut-off at or above the threshold, and at least two
such zones in the patch); a level-only zone decays along the AR mean path.
Each block is one `~` over a `product_distribution`.

### Likelihood

The share renewal [`zone_share_renewal`](@ref) gives each zone's
infections as its share of the fixed patch infections, the delay
[`zone_reports`](@ref) and binning [`zone_report_increments`](@ref) give
the expected confirmed reports per vintage window, and the observed zone
increments of each patch and vintage follow a Dirichlet-multinomial on the
allocated total with concentration `κ π`
([`zone_composition_logpdf`](@ref)), as one `@addlogprob!` over every
cell.

`mixing = true` adds one shared within-patch mixing fraction per patch,
`ε_p ~ Beta(1, 100)`, redistributing force through the gravity kernel
`zd.mixing_kernel`. `deaths = true` adds one Dirichlet-multinomial on the
cumulative allocated zone deaths at the last vintage, through the
infection-to-death-confirmation PMF, with no extra multiplier. Both are off
by default.

### Deterministics

Flattened column-major where a matrix: `delta_knots_zone` `(n_zones ×
n_knots)`, `share_knots_zone` `(n_zones × n_knots)`, `delta_T_zone`,
`share_T_zone`, `share_start_zone` (the initial shares `w_z(t_0)`),
`R_T_zone` (the implied `I_z / Λ_z` at the cut-off, `NaN` below ten
cumulative zone infections), `region_sd_zone` (`σ_L`),
`region_drift_sd_zone` (`σ_δ`), `region_halflife_zone` (`h`),
`composition_rho_zone` (`ρ`) and, with mixing, `mixing_epsilon_zone`. Daily
trajectories are rebuilt from these by [`zone_forward`](@ref).
"""
@model function bvd_zone(zd;
        mixing::Bool = false,
        deaths::Bool = false,
        share_scale::Real = 2.0,
        rt_floor::Real = 10.0,
        week::Integer = 7,
        region_sd_prior = truncated(Normal(0, 0.3); lower = 0),
        region_drift_sd_prior = truncated(Normal(0, 0.1); lower = 0),
        region_halflife_prior = LogNormal(log(42), 0.6),
        rho_prior = truncated(Normal(0, 0.1); lower = 0, upper = 1),
        mixing_prior = Beta(1, 100),
        offset_prior = Normal(0, 1))
    nz = size(zd.counts, 1)
    np = length(zd.patch_ranges)
    K = length(zd.knots)
    z_w ~ product_distribution(fill(offset_prior, nz))
    σ_level ~ region_sd_prior
    z_level ~ product_distribution(fill(offset_prior, nz))
    δ_halflife ~ region_halflife_prior
    σ_δ ~ region_drift_sd_prior
    z_drift ~ product_distribution(
        fill(offset_prior, max(zd.n_walking * (K - 1), 1)))
    ρ ~ rho_prior
    ## The mixing fractions are sampled only when the kernel is used; against
    ## `mixing = false` they would be prior-only dimensions.
    if mixing
        ε_mix ~ product_distribution(fill(mixing_prior, np))
    else
        ε_mix = nothing
    end
    w0 = zone_initial_shares(z_w, zd.patch_ranges, share_scale)
    φ = exp2(-week / δ_halflife)
    δ_knots = zone_deviation_knots(z_level, z_drift, σ_level, σ_δ, φ,
        zd.patch_ranges, zd.walking, zd.walk_index, zd.n_walking, K)
    fw = zone_forward(zd, δ_knots, w0, ε_mix)
    κ = (1 - ρ) / ρ
    @addlogprob! zone_composition_logpdf(zd.counts, fw.increments,
        zd.cell_patch, zd.cell_vintage, zd.cell_total, zd.cell_const,
        zd.patch_ranges, κ)
    if deaths
        ## One composition of the cumulative allocated deaths at the last
        ## vintage, through the death-confirmation delay.
        death_daily = zone_reports(fw.infections, zd.death_pmf, w0,
            zd.patch_ranges, zd.t0, zd.death_pre)
        D = zone_report_increments(death_daily, w0, zd.patch_ranges,
            zd.death_days, zd.t0, zd.death_pre_cum)
        @addlogprob! zone_composition_logpdf(zd.death_counts, D,
            zd.death_cell_patch, zd.death_cell_vintage,
            zd.death_cell_total, zd.death_cell_const, zd.patch_ranges, κ)
    end
    nd = zd.n - zd.t0 + 1
    cum = _zone_cumulative_infections(fw.infections, w0, zd.patch_ranges,
        zd.infections_pre)
    R_T_zone := _zone_rt_at(fw.infections, fw.forces, cum, nd, rt_floor)
    delta_knots_zone := vec(δ_knots)
    delta_T_zone := δ_knots[:, K]
    share_T_zone := fw.shares[nd, :]
    share_start_zone := w0
    share_knots_zone := vec(_zone_shares_at_knots(fw.shares, zd.knots, zd.t0))
    region_sd_zone := σ_level
    region_drift_sd_zone := σ_δ
    region_halflife_zone := δ_halflife
    composition_rho_zone := ρ
    if mixing
        mixing_epsilon_zone := ε_mix
    end
    return (; shares = fw.shares, forces = fw.forces,
        infections = fw.infections, increments = fw.increments,
        δ_knots, w0, cum)
end

## --- Fixed inputs from the parent chain ---------------------------------

## Draws of `key` from the parent chain, with a clear error naming the key
## when the chain does not carry it.
function _zone_parent_draws(chn, key::Symbol)
    _has_key(chn, key) || error(
        "zone_fit_inputs: the parent chain carries no `$(key)`; it must be " *
        "a `bvd_joint` chain sampled with the patch structure on.")
    return _draws(chn, key)
end

## Index of the parent draw whose `C_T` sits at quantile `q`, or `nothing`
## for the posterior mean.
function _zone_parent_index(chn, parent_summary::Symbol)
    parent_summary === :mean && return nothing
    q = parent_summary === :draw_low ? 0.05 :
        parent_summary === :draw_high ? 0.95 :
        error("zone_fit_inputs: `parent_summary` must be :mean, :draw_low " *
              "or :draw_high, got $(repr(parent_summary)).")
    C_T = _zone_parent_draws(chn, _ZONE_PARENT_KEYS.C_T)
    target = quantile(C_T, q)
    return argmin(abs.(C_T .- target))
end

## Generation-interval PMF (lag 1) of one parent draw, discretised as the
## parent did.
function _zone_gi_pmf(α::Real, θ::Real)
    pmf = discretise_censored(Gamma(α, θ), ZONE_GI_NMAX)
    return pmf[2:end] ./ sum(pmf[2:end])
end

## Mean over draws of a per-draw PMF, or that of a single draw.
function _zone_mean_pmf(build, ndraws::Integer, idx)
    idx === nothing || return build(idx)
    acc = build(1)
    for i in 2:ndraws
        acc = acc .+ build(i)
    end
    return acc ./ ndraws
end

"""
$(TYPEDSIGNATURES)

The stage-1 quantities the zone model conditions on, read from a
`bvd_joint` parent chain: the patch infections `Ī_p(t)` `(n_patches × n)`,
the generation-interval PMF `g` (lag 1), the infection-to-confirmed-report
PMF `f = incubation ⊛ receipt` (lag 0, the delay the parent's per-province
composition applies to infections) and the infection-to-confirmed-death PMF
`death_pmf = incubation ⊛ onset-to-death ⊛ receipt`.

With `parent_summary = :mean` each is a posterior mean: the infections are
the exponential of the mean log infections per day, the PMFs the mean of
the per-draw PMFs. With `:draw_low` or `:draw_high` every quantity comes
from the single draw whose national `C_T` sits nearest its 5th or 95th
percentile, for the feedback sensitivity. Returns
`(; I_bar, g, f, death_pmf, draw)`.
"""
function zone_parent_inputs(chn; parent_summary::Symbol = :mean)
    keys_ = _ZONE_PARENT_KEYS
    idx = _zone_parent_index(chn, parent_summary)
    infs = _draw_vectors(chn, keys_.infections)
    ndraws = length(infs)
    logI = if idx === nothing
        acc = log.(safe_rate.(Float64.(infs[1])))
        for i in 2:ndraws
            acc .+= log.(safe_rate.(Float64.(infs[i])))
        end
        acc ./ ndraws
    else
        log.(safe_rate.(Float64.(infs[idx])))
    end
    α = _zone_parent_draws(chn, keys_.gi_alpha)
    θ = _zone_parent_draws(chn, keys_.gi_theta)
    inc_m = _zone_parent_draws(chn, keys_.inc_mean)
    inc_s = _zone_parent_draws(chn, keys_.inc_sd)
    rec_m = _zone_parent_draws(chn, keys_.receipt_mean)
    rec_s = _zone_parent_draws(chn, keys_.receipt_sd)
    g = _zone_mean_pmf(i -> _zone_gi_pmf(α[i], θ[i]), ndraws, idx)
    inc_pmf(i) = discretise_censored(lognormal_meansd(inc_m[i], inc_s[i]),
        ZONE_INCUBATION_NMAX)
    rec_pmf(i) = discretise_censored(lognormal_meansd(rec_m[i], rec_s[i]),
        ZONE_RECEIPT_NMAX)
    f = _zone_mean_pmf(i -> convolve_pmf(inc_pmf(i), rec_pmf(i)), ndraws, idx)
    death_pmf = if _has_key(chn, keys_.death_confirmation)
        dc = _draw_vectors(chn, keys_.death_confirmation)
        _zone_mean_pmf(i -> convolve_pmf(inc_pmf(i), Float64.(dc[i])),
            ndraws, idx)
    else
        Float64[]
    end
    return (; log_infections = logI, g, f, death_pmf, draw = idx)
end

## Cumulative count of `zone` in `prov` at the last vintage of `history`, or
## zero when the zone is absent from the block.
function _zone_last_cumulative(history, prov::AbstractString,
        zone::AbstractString)
    haskey(history, prov) || return 0
    haskey(history[prov], zone) || return 0
    c = history[prov][zone].counts
    return isempty(c) ? 0 : Int(c[end])
end

## Within-patch gravity kernel over the zones, column-stochastic within each
## patch and zero across patches, from the zone populations and centroids.
function _zone_mixing_kernel(pops::AbstractVector, coords::AbstractVector,
        patch_ranges::AbstractVector{<:UnitRange})
    nz = length(pops)
    K = zeros(Float64, nz, nz)
    for zs in patch_ranges
        length(zs) >= 2 || continue
        sub = province_importation_kernel(pops[zs];
            distances = province_distance_matrix(coords[zs]))
        for (j, q) in enumerate(zs)
            s = sum(@view sub[:, j])
            s > 0 || continue
            for (i, z) in enumerate(zs)
                K[z, q] = sub[i, j] / s
            end
        end
    end
    return K
end

"""
$(TYPEDSIGNATURES)

Everything the zone model and its render need, built once from the
parent chain and the observations: the fixed stage-1 inputs
([`zone_parent_inputs`](@ref)), the zone units and their data, the grid
and the sampler's starting point.

Zones are the health-zone rows of `obs.zone_confirmed_history`
([`zone_increment_matrix`](@ref)), nested in the stage-1 patches of
[`PROVINCE_NAMES`](@ref) and ordered patch by patch; `unallocated`
pseudo-rows are never units. Per vintage the zone increments are the
consecutive-vintage differences clamped at zero, the first vintage's
increment its cumulative, and the allocated patch total their sum. Cells
with a zero allocated total are dropped here, so the model scores only
positive totals. The walking set is the zones whose cumulative confirmed
count at the cut-off is at least `walk_threshold`, in patches with at least
two such zones (a single walking zone's centred innovations would vanish).
The grid starts `lead_days` before the first vintage and carries weekly
knots ([`knot_days`](@ref)) to the cut-off.

`zones` is the metadata table of [`load_health_zones`](@ref), read from the
package data by default, for labels and, with mixing, populations and
centroids. Returns a named tuple whose `model_data` field is the plain-array
input of [`bvd_zone`](@ref), alongside the zone keys and labels, the patch
index and labels, the vintage days and dates, the cumulative counts, the
walking mask, the knots, the grid start and the initial-share start values.
"""
function zone_fit_inputs(parent_chain, obs;
        parent_summary::Symbol = :mean,
        walk_threshold::Integer = 30,
        lead_days::Integer = 42,
        week::Integer = 7,
        zones = _default_health_zones(),
        patch_names::AbstractVector = PROVINCE_NAMES,
        patch_labels::AbstractVector = PROVINCE_LABELS)
    hasproperty(obs, :zone_confirmed_history) || error(
        "zone_fit_inputs: `obs` carries no `zone_confirmed_history`; load " *
        "the observations with a manifest that has the zone tables.")
    n = obs.n
    parent = zone_parent_inputs(parent_chain; parent_summary)
    np = length(patch_names)
    length(parent.log_infections) == np * n || error(
        "zone_fit_inputs: the parent's `infections_patch` holds " *
        "$(length(parent.log_infections)) entries but $np patches by $n " *
        "days is $(np * n); the parent must be fitted to the same data " *
        "cut-off.")
    I_bar = exp.(reshape(parent.log_infections, np, n))
    ## Zone units and counts, patch by patch.
    per_patch = zone_increment_matrix(obs.zone_confirmed_history, patch_names)
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
    labels = _zone_labels(zones, zone_province, zone_names)
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
    ## Grid and knots.
    t0 = clamp(days[1] - lead_days, 1, n)
    knots = knot_days(n; week, start = t0)
    fixed = zone_fixed_terms(I_bar, parent.g, parent.f, t0)
    ## Optional death composition: cumulative allocated deaths at the last
    ## vintage, zones absent from the death table counting zero.
    death_history = hasproperty(obs, :zone_death_history) ?
                    obs.zone_death_history : Dict{String, Any}()
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
    death_fixed = zone_fixed_terms(I_bar, parent.g, death_pmf, t0)
    ## Within-patch mixing kernel, when the metadata covers every zone.
    mixing_kernel = _zone_mixing_kernel_or_zeros(zones, zone_province,
        zone_names, patch_ranges)
    ## Data-informed starting point for the initial shares: the log observed
    ## cumulative share at the first vintage with a pseudo-count, centred
    ## within patch and divided by the softmax scale.
    z_w_start = zeros(Float64, nz)
    for zs in patch_ranges
        isempty(zs) && continue
        lg = [log(counts[z, 1] + 0.5) for z in zs]
        z_w_start[zs] .= (lg .- (sum(lg) / length(lg))) ./ 2
    end
    model_data = (; counts, cell_patch, cell_vintage, cell_total, cell_const,
        days, I_bar, g = parent.g, f = parent.f, patch_ranges, knots, t0, n,
        walking = collect(walking), walk_index, n_walking,
        fixed.force_pre, fixed.report_pre, fixed.report_pre_cum,
        fixed.infections_pre, mixing_kernel,
        death_counts, death_cell_patch,
        death_cell_vintage = ones(Int, length(death_cell_patch)),
        death_cell_total, death_cell_const, death_days = [n],
        death_pmf, death_pre = death_fixed.report_pre,
        death_pre_cum = death_fixed.report_pre_cum)
    dates = [obs.seeding + Day(d - 1) for d in days]
    return (; model_data, zone_keys, zone_labels = labels, zone_province,
        zone_names, patch_of_zone, patch_ranges,
        patch_names = collect(String, patch_names),
        patch_labels = collect(String, patch_labels),
        days, dates, counts, cumulative, walking = collect(walking),
        walk_threshold, knots, t0, n, week,
        seeding = obs.seeding, cutoff = obs.cutoff, z_w_start,
        I_bar, g = parent.g, f = parent.f, death_pmf,
        parent_summary, parent_draw = parent.draw)
end

## The package's health-zone metadata, or `nothing` when the file is absent.
function _default_health_zones()
    path = joinpath(@__DIR__, "..", "..", "data", "health_zones.csv")
    isfile(path) || return nothing
    return load_health_zones(path)
end

## Display label per zone from the metadata, falling back to the key.
function _zone_labels(zones, zone_province::AbstractVector,
        zone_names::AbstractVector)
    fallback(z) = titlecase(replace(zone_names[z], "_" => " "))
    zones === nothing && return [fallback(z) for z in eachindex(zone_names)]
    lookup = Dict((r.province, r.zone) => r.label for r in zones)
    return [get(lookup, (zone_province[z], zone_names[z]), fallback(z))
            for z in eachindex(zone_names)]
end

## The mixing kernel when every zone has metadata, else an all-zero matrix
## (mixing is then refused by `fit_zone`).
function _zone_mixing_kernel_or_zeros(zones, zone_province, zone_names,
        patch_ranges)
    nz = length(zone_names)
    zones === nothing && return zeros(Float64, nz, nz)
    lookup = Dict((r.province, r.zone) => r for r in zones)
    rows = [get(lookup, (zone_province[z], zone_names[z]), nothing)
            for z in 1:nz]
    any(isnothing, rows) && return zeros(Float64, nz, nz)
    pops = Float64[r.population for r in rows]
    coords = [(r.lat, r.lon) for r in rows]
    return _zone_mixing_kernel(pops, coords, patch_ranges)
end

## --- Fitting -------------------------------------------------------------

"""
$(TYPEDSIGNATURES)

Per-chain starting points for [`bvd_zone`](@ref), in unconstrained space.
The constrained start is the data-informed point of the specification
(`δ = 0`, the initial shares from the observed first-vintage cumulative,
`σ_L = 0.2`, `σ_δ = 0.05`, `h = 42`, `ρ = 0.05`, and `ε = 0.01` with
mixing), built through a `VarInfo`, linked, then jittered per chain with
`N(0, jitter²)` noise so the chains do not share one start. Returns
`(; inits, x0, logp)`, the `InitFromVector` strategies, the unjittered
unconstrained vector and each chain's initial log joint density.
"""
function zone_initial_params(model, inputs; chains::Integer = 2,
        seed::Integer = 20260518, jitter::Real = 0.1, mixing::Bool = false)
    zd = inputs.model_data
    nz = length(inputs.z_w_start)
    K = length(zd.knots)
    params = (; z_w = inputs.z_w_start, σ_level = 0.2, z_level = zeros(nz),
        δ_halflife = 42.0, σ_δ = 0.05,
        z_drift = zeros(max(zd.n_walking * (K - 1), 1)), ρ = 0.05)
    if mixing
        params = merge(params,
            (; ε_mix = fill(0.01, length(zd.patch_ranges))))
    end
    rng = MersenneTwister(seed)
    vi = VarInfo(rng, model, InitFromParams(params))
    vi_linked = link(vi, model)
    x0 = collect(vi_linked[:])
    ldf = LogDensityFunction(model, getlogjoint, vi_linked)
    xs = [x0 .+ jitter .* randn(rng, length(x0)) for _ in 1:chains]
    logp = [LogDensityProblems.logdensity(ldf, x) for x in xs]
    inits = [InitFromVector(x, ldf) for x in xs]
    return (; inits, x0, logp)
end

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
`obs`. Builds the fixed inputs with [`zone_fit_inputs`](@ref), starts every
chain at the data-informed point of [`zone_initial_params`](@ref) and runs
[`nuts_sample`](@ref) with a diagonal metric at the given settings. Logs
each chain's initial log density before sampling and its adapted step size,
divergences and the fraction of iterations at the tree-depth cap after.
`mixing`, `deaths` and `parent_summary` select the sensitivity variants;
`zones`, `patch_names` and `patch_labels` pass through to
[`zone_fit_inputs`](@ref) and any other keyword to [`nuts_sample`](@ref).
Returns the chain.
"""
function fit_zone(parent_chain, obs;
        samples::Integer = 600,
        chains::Integer = 2,
        n_adapts::Integer = 400,
        target_accept::Real = 0.8,
        max_depth::Integer = 8,
        seed::Integer = 20260518,
        callback = nothing,
        mixing::Bool = false,
        deaths::Bool = false,
        parent_summary::Symbol = :mean,
        walk_threshold::Integer = 30,
        lead_days::Integer = 42,
        jitter::Real = 0.1,
        zones = _default_health_zones(),
        patch_names::AbstractVector = PROVINCE_NAMES,
        patch_labels::AbstractVector = PROVINCE_LABELS,
        kwargs...)
    inputs = zone_fit_inputs(parent_chain, obs; parent_summary,
        walk_threshold, lead_days, zones, patch_names, patch_labels)
    zd = inputs.model_data
    mixing && !any(!iszero, zd.mixing_kernel) &&
        error(
            "fit_zone: mixing needs a population and centroid for every zone " *
            "in the health-zone metadata.")
    deaths && isempty(zd.death_cell_patch) &&
        error(
            "fit_zone: deaths = true but the observations carry no allocated " *
            "zone deaths.")
    model = bvd_zone(zd; mixing, deaths)
    start = zone_initial_params(model, inputs; chains, seed, jitter, mixing)
    @info "fit_zone: starting" zones=length(inputs.zone_keys) walking=count(
        inputs.walking) knots=length(inputs.knots) cells=length(
        zd.cell_patch) dimension=length(start.x0) initial_logp=start.logp
    t = time()
    chn = nuts_sample(model; samples, chains, n_adapts, target_accept,
        max_depth, seed, callback, init = start.inits, kwargs...)
    elapsed = time() - t
    steps = _zone_stat(chn, :step_size)
    depth = _zone_stat(chn, :tree_depth)
    div = _zone_stat(chn, :numerical_error)
    per_chain(m, f) = m === nothing ? missing :
                      [f(view(m, :, c)) for c in 1:size(m, 2)]
    @info "fit_zone: finished" minutes=round(elapsed/60; digits = 1) step_size=per_chain(
        steps, x->x[end]) depth_cap_fraction=per_chain(
        depth, x->mean(x .>= max_depth)) divergences=per_chain(
        div, x->Int(sum(x)))
    return chn
end
