# Shared Gauss-Legendre quadrature: a single reusable integrator backs
# every forward integral in the package (the gamma onset-to-death
# convolution for deaths, the at-risk person-time integral for exports,
# the deaths-among-exports convolution, and the counterfactual and
# forecast integrals). Also defines `ExportDeathDelay`, which carries
# a precomputed onset-to-death CDF on a fixed grid.

# Reduce an integrand `g` over the reference domain `[-1, 1]` against the
# `alg` Gauss-Legendre rule. The package only ever integrates on `[-1, 1]`
# (each method folds the change of variables into `g`), so the SciML
# `scale = (ub - lb) / 2 = 1` and `shift = 0`, leaving the bare weighted
# sum `Σ wᵢ g(xᵢ)`. Calling `g` directly here — rather than through a
# SciML `IntegralProblem`/`solve` parameter boundary — lets Julia
# specialise on `g`'s concrete return type, so the result type is inferred
# (AD tangents and `Dual`s propagate) and the loop allocates nothing. The
# accumulation order matches the reference `gauss_legendre`
# (`w₁ g(x₁) + …`), so values are reproduced bit-for-bit. `alg.nodes` and
# `alg.weights` are the very nodes/weights `gausslegendre(n)` built into
# the algorithm object, so the rule is identical.
@inline function _gl_reduce(g::G, alg) where {G}
    nodes = alg.nodes
    weights = alg.weights
    @inbounds acc = weights[1] * g(nodes[1])
    @inbounds for i in 2:length(nodes)
        acc += weights[i] * g(nodes[i])
    end
    return acc
end

"""
Integrate a scalar function `f` over `[lo, hi]` by Gauss-Legendre
quadrature. The reference domain `[-1, 1]` is mapped onto `[lo, hi]`,
so any forward integral in the package can share one integrator. Returns
zero when `hi <= lo`. The default scheme is [`CUMULATIVE_INTEGRAL_ALG`](@ref).
"""
function integrate(f::F, lo, hi; alg = CUMULATIVE_INTEGRAL_ALG) where {F}
    hi <= lo && return zero(hi - lo)
    halfwidth = (hi - lo) / 2
    g = let f = f, halfwidth = halfwidth, lo = lo
        u -> f(halfwidth * (u + one(u)) + lo)
    end
    return halfwidth * _gl_reduce(g, alg)
end

# Change of variables clustering the reference nodes towards `hi`. With
# `d = hi - s` measured back from the upper limit and `v = (u + 1) / 2`,
# the map `d = span · vᵖ` sends the dense end of the reference grid to
# `s = hi` and stretches the sparse end out to `s = lo`. The Jacobian
# `dd/du = span · p · vᵖ⁻¹ / 2` is folded into the integrand so a single
# fixed-node reduction covers the whole domain.
function _clustered_kernel(u, f, hi, span, expo)
    v = (u + one(u)) / 2
    d = span * v^expo
    return f(hi - d) * (span * expo * v^(expo - one(expo)) / 2)
end

"""
Integrate a scalar function `f` over `[lo, hi]` by Gauss-Legendre
quadrature with the nodes clustered towards `hi`, resolving features of
size `scale` near that limit. Use for onset-to-death convolutions, whose
integrand mass piles up against the cut-off `hi` over a window set by the
sampled delay scale: the uniform [`integrate`](@ref) spread across a wide
`[lo, hi]` would resolve that peak too coarsely. The whole domain is
still covered (no truncation), so a delay wider than `[lo, hi]` is
integrated exactly; when `scale ≥ hi - lo` the map reduces to the uniform
method. Returns zero when `hi <= lo`. The default scheme is
[`DEATH_INTEGRAL_ALG`](@ref).
"""
function integrate(f::F, lo, hi, scale; alg = DEATH_INTEGRAL_ALG) where {F}
    hi <= lo && return zero(hi - lo)
    span = hi - lo
    # `expo` places ~half the nodes within `scale` of `hi`; `expo = 1`
    # (uniform) once the feature is as wide as the domain. Fall back to
    # the uniform rule for a degenerate scale so a non-finite delay
    # moment cannot produce a `NaN` exponent on an AD proposal.
    (isfinite(scale) && scale > zero(scale)) ||
        return integrate(f, lo, hi; alg)
    expo = max(one(span), log(span / scale) / log(oftype(span, 2)))
    # The clustered Jacobian is folded into the kernel, so the reduction
    # carries no outer `halfwidth` factor (unlike the uniform method).
    g = let f = f, hi = hi, span = span, expo = expo
        u -> _clustered_kernel(u, f, hi, span, expo)
    end
    return _gl_reduce(g, alg)
end

# Clustering scale for a delay distribution, `mean + K·std`: the width of
# the window near the cut-off where the convolution integrand has mass.
# Scales with the sampled `std`, so gradients flow through it on the
# AD-supported parameter path used by the clustered `integrate`.
_delay_scale(dist) = mean(dist) + DELAY_SUPPORT_K * std(dist)

"""
Integrate a cumulative-incidence trajectory `cumulative` (a callable
`C(s)`) over `[lo, hi]`. Backs the at-risk person-time export integral
``\\int_{T-w}^{T} C(s)\\, ds``. Uses [`CUMULATIVE_INTEGRAL_ALG`](@ref).
"""
function integrate_cumulative(cumulative, lo, hi; alg = CUMULATIVE_INTEGRAL_ALG)
    return integrate(cumulative, lo, hi; alg)
end

"""
Deaths-among-exports convolution
``\\int_{lo}^{hi} C(s)\\, F_d(T - s)\\, ds`` with `F_d` the `delay_dist`
onset-to-death CDF. The CDF is itself written as the inner integral of
the density, ``F_d(x) = \\int_0^x f_d(u)\\,du``, so the whole expression
differentiates through the density alone.
Previously, the reverse-mode AD backend did not support the gamma CDF
shape-parameter derivative). Outer and inner integrals both use
[`CUMULATIVE_INTEGRAL_ALG`](@ref).
"""
function integrate_exports_deaths(cumulative, delay_dist, lo, hi, T;
        alg = CUMULATIVE_INTEGRAL_ALG)
    cdf_to(x) = integrate(u -> pdf(delay_dist, u), zero(x), x; alg)
    g = let cumulative = cumulative, T = T
        s -> cumulative(s) * cdf_to(T - s)
    end
    return integrate(g, lo, hi; alg)
end

"""
Deaths-among-exports convolution specialised for `Gamma` delay
distributions. Same expression as the generic method —
``\\int_{lo}^{hi} C(s)\\, F_d(T - s)\\, ds`` — but the onset-to-death
CDF is evaluated in closed form via [`_gamma_cdf`](@ref) rather than
re-integrated from the density at each outer quadrature node. Drops one
entire quadrature level (just the outer integral remains) and
sidesteps the singularity that fixed-node Gauss-Legendre hits at the
head-form ``\\int_0^y t^{\\alpha-1} e^{-t} \\, du`` for `α < 1`. The
shape-parameter derivative the AD backend needed for `cdf(::Gamma, ·)`
is supplied by `_gamma_cdf`'s rrule (see `src/gamma_cdf.jl`).
"""
function integrate_exports_deaths(cumulative, delay_dist::Gamma, lo, hi, T;
        alg = CUMULATIVE_INTEGRAL_ALG)
    α, θ = delay_dist.α, delay_dist.θ
    g = let cumulative = cumulative, T = T, α = α, θ = θ
        s -> cumulative(s) * _gamma_cdf(α, θ, T - s)
    end
    return integrate(g, lo, hi; alg)
end

"""
Onset-to-death delay carrying its CDF `F_d` precomputed on an evenly
spaced grid over `[0, gmax]`, built once and reused across every outer
node and bin edge of the deaths-among-exports convolution rather than
re-integrating the density at each (see [`integrate_exports_deaths`](@ref)).
`gmax` must cover the largest delay argument reached, the detection window
`w`. Pass one in place of the raw delay distribution to select this path by
dispatch; the distribution methods are unchanged and remain the reference.

`F_d` is a cumulative trapezoid of the density, so the convolution still
differentiates through the density alone (the AD backend lacks the Gamma
CDF shape derivative). The density is never evaluated at `0`, whose Gamma
shape derivative is `0·log 0 = NaN` under AD; for `f_d(0) = 0` (Gamma
shape > 1) treating it as zero is exact.
"""
struct ExportDeathDelay{D, T}
    dist::D
    gmax::T
    dx::T
    F::Vector{T}
end

function ExportDeathDelay(dist, gmax::Real;
        npts::Integer = EXPORT_DELAY_GRID_POINTS)
    g = float(gmax)
    dx = g / (npts - 1)
    Tt = typeof(pdf(dist, dx) * dx)
    F = Vector{Tt}(undef, npts)
    F[1] = zero(Tt)
    prev = zero(Tt)
    @inbounds for i in 2:npts
        fx = pdf(dist, (i - 1) * dx)
        F[i] = F[i - 1] + (prev + fx) * dx / 2
        prev = fx
    end
    return ExportDeathDelay(dist, g, dx, F)
end

# Linear interpolation of the precomputed CDF: F_d(0) = 0 and flat past
# `gmax` (all delay mass within the window has accumulated by then).
@inline function _cdf_to(ed::ExportDeathDelay, y)
    y <= zero(y) && return zero(eltype(ed.F))
    y >= ed.gmax && return @inbounds ed.F[end]
    pos = y / ed.dx
    i = floor(Int, pos) + 1
    frac = pos - (i - 1)
    return @inbounds ed.F[i] + frac * (ed.F[i + 1] - ed.F[i])
end

"""
Deaths-among-exports convolution using a precomputed [`ExportDeathDelay`](@ref):
``\\int_{lo}^{hi} C(s)\\, F_d(T - s)\\, ds`` with `F_d` interpolated off the
grid rather than re-integrated at each node. Mathematically identical to
the distribution method up to the grid resolution.
"""
function integrate_exports_deaths(cumulative, ed::ExportDeathDelay, lo, hi, T;
        alg = CUMULATIVE_INTEGRAL_ALG)
    g = let cumulative = cumulative, T = T, ed = ed
        s -> cumulative(s) * _cdf_to(ed, T - s)
    end
    return integrate(g, lo, hi; alg)
end

"""
Precomputed BVD-suspected trajectory carrying the unit-ascertainment
trajectory ``\\mu_{BVD,0}(s_j) = \\int_0^{s_j} e^{r u} f_{rep}(s_j - u)\\,du``
evaluated on a fixed Gauss-Legendre node set over `[0, T_max]`. The same
nodes back the outer quadrature of the confirmed-cases convolution at every
daily bin edge, so the expensive per-incidence-time piece is built once
per draw and reused across all `T_k`. Mirrors [`ExportDeathDelay`](@ref):
a struct that precomputes the expensive per-incidence-time piece of the
integrand, dispatched to by a specialised [`delay_convolution`](@ref)
method.

The DRC ascertainment fraction is *not* baked into the precomputation:
the per-bin random-effect ascertainment of the daily likelihood
(see [`daily_ascertainment_model`](@ref)) multiplies the convolution
output per bin, after `delay_convolution(d, t_edges, f_lab)` returns the
unit-ascertainment cumulative ``I_{lab,0}(T_k)``. Pass the
unit-ascertainment trajectory through and apply `p_drc_t` at the bin
level.

The clustered Gauss-Legendre map of [`integrate`](@ref) is replaced by
a uniform-grid map so the same `bvd_at_nodes` vector serves every bin
edge. Resolution scales with `npts`, which defaults to the
[`DEATH_INTEGRAL_ALG`](@ref) node count.
"""
struct DailyBVDTrajectory{T}
    nodes::Vector{T}           # incidence-time samples s_j ∈ [0, T_max]
    weights::Vector{T}         # quadrature weights (scaled to T_max / 2)
    bvd_at_nodes::Vector{T}    # μ_BVD,0(s_j) — propagates the AD type
    w_bvd::Vector{T}           # weights[j] · bvd_at_nodes[j], edge-independent
    T_max::T
end

"""
Build a [`DailyBVDTrajectory`](@ref) precomputation. `T_max` is the
upper edge of the bin grid (typically the latent seeding-to-cut-off
time `T`), `r` is the growth rate, and `f_rep` is the onset-to-report
kernel. Allocates one quadrature node / weight pair scaled to
`[0, T_max]`, then evaluates `μ_BVD,0(s_j)` (the unit-ascertainment
trajectory) at each node.

DRC ascertainment is applied separately by the daily likelihood, after
the convolution returns, so the precomputation is independent of any
per-bin or pooled `p_drc`.

The Gamma-specialised [`delay_convolution`](@ref) is used for the
per-node BVD evaluations, so each is closed-form rather than a nested
quadrature.
"""
function DailyBVDTrajectory(T_max::Real, r::Real, f_rep;
        npts::Integer = 64)
    Tt = promote_type(typeof(float(T_max)), typeof(r))
    raw_nodes, raw_weights = FastGaussQuadrature.gausslegendre(npts)
    half = T_max / 2
    nodes = Tt[half * (u + one(u)) for u in raw_nodes]
    weights = Tt[half * w for w in raw_weights]
    bvd = Tt[delay_convolution(one(Tt), r, s, f_rep) for s in nodes]
    w_bvd = weights .* bvd
    return DailyBVDTrajectory{Tt}(nodes, weights, bvd, w_bvd,
        convert(Tt, T_max))
end

"""
Cumulative confirmed-cases convolution evaluated at a vector of bin
edges `t_edges`, reusing the BVD-suspected node evaluations carried by
`d`:

```math
I_{lab}(T_k) = \\int_0^{T_k} \\mu_{BVD}(s)\\, f_{lab}(T_k - s)\\, ds.
```

The outer Gauss-Legendre nodes are fixed across `[0, T_max]` so each
``\\mu_{BVD}(s_j)`` is precomputed once. For every bin edge, the
integrand at each node is multiplied by `f_lab(T_k - s_j)` (zero past
the upper limit) and reduced against the shared weights. Returns one
cumulative value per edge.
"""
function delay_convolution(d::DailyBVDTrajectory, t_edges::AbstractVector,
        f_lab)
    Tt = eltype(d.bvd_at_nodes)
    n = length(t_edges)
    out = Vector{Tt}(undef, n)
    for k in 1:n
        T_k = t_edges[k]
        acc = zero(Tt)
        @inbounds for j in eachindex(d.nodes)
            δ = T_k - d.nodes[j]
            δ <= zero(δ) && continue
            ## `w_bvd[j]` folds in the edge-independent weight × μ_BVD,0
            ## product, leaving only the per-edge lab-delay density.
            acc += d.w_bvd[j] * pdf(f_lab, δ)
        end
        out[k] = acc
    end
    return out
end

"""
[`delay_convolution`](@ref) over a [`DailyBVDTrajectory`](@ref)
specialised to a `Gamma` lab-delay distribution. The Gamma density
``f(δ) = δ^{α-1} e^{-δ/θ} / (Γ(α)\\, θ^{α})`` splits into a
node-dependent shape ``δ^{α-1} e^{-δ/θ}`` and the constant
``1 / (Γ(α)\\, θ^{α})``. The constant is evaluated once per call — so
`loggamma(α)` is computed once rather than at every (edge, node) pair —
and factored back in per edge. Numerically identical to the generic
`pdf`-based method.
"""
function delay_convolution(d::DailyBVDTrajectory, t_edges::AbstractVector,
        f_lab::Gamma)
    α, θ = f_lab.α, f_lab.θ
    αm1 = α - one(α)
    invθ = inv(θ)
    Z = exp(-loggamma(α) - α * log(θ))   # 1 / (Γ(α) θ^α)
    Tt = promote_type(eltype(d.bvd_at_nodes), typeof(Z), typeof(invθ))
    n = length(t_edges)
    out = Vector{Tt}(undef, n)
    for k in 1:n
        T_k = t_edges[k]
        acc = zero(Tt)
        @inbounds for j in eachindex(d.nodes)
            δ = T_k - d.nodes[j]
            δ <= zero(δ) && continue
            acc += d.w_bvd[j] * exp(αm1 * log(δ) - δ * invθ)
        end
        out[k] = acc * Z
    end
    return out
end

"""
Precomputed at-risk export trajectory: the growth curve ``e^{r u}`` on a
fixed Gauss-Legendre node set over ``[0, t_{last}]``, built once per draw
and reused across every export bin edge. The export streams are the
travel-gated analogue of the DRC convolutions: the expected detected
exports by `t` is a convolution of the epidemic with the detection
survival,

```math
\\int_0^t e^{r u}\\, S_{det}(t - u)\\, du,
    \\qquad S_{det}(x) = 1 - F_{det}(x),
```

and the deaths-among-exports weight that integrand further by the
infection→death CDF ``F_{death}``. Mirrors [`DailyBVDTrajectory`](@ref):
the per-incidence-time piece (here just the weighted growth) is evaluated
once at the shared nodes, and each bin edge reuses it, multiplying by the
edge-specific survival kernel and skipping nodes past the edge. The travel
rate `q`, ascertainment `p` and `CFR` are applied by the caller, so the
trajectory carries person-time, not counts.
"""
struct ExportRiskTrajectory{T}
    nodes::Vector{T}      # incidence-time samples u_j ∈ [0, t_last]
    w_growth::Vector{T}   # weights[j] · e^{r·u_j}, edge-independent
    r::T                  # growth rate, for the exact growth integral
    t_last::T
end

"""
Build an [`ExportRiskTrajectory`](@ref) over `[0, t_last]` for growth rate
`r`. Lays down `npts` Gauss-Legendre nodes scaled to the interval and
precomputes the weighted growth ``w_j\\, e^{r u_j}`` at each, the
edge-independent piece reused across every export edge. `npts = 32`
resolves the at-risk and deaths-among-exports integrals to a relative
error below `1e-4` across the supported growth range (verified against a
1024-node reference), well under the Poisson/NegBinomial noise on the
O(1) observed export counts, while halving the per-gradient gamma-CDF
evaluations in [`export_at_risk`](@ref) / [`export_death_at_risk`](@ref)
relative to the previous 64.
"""
function ExportRiskTrajectory(t_last::Real, r::Real; npts::Integer = 32)
    Tt = promote_type(typeof(float(t_last)), typeof(r))
    raw_nodes, raw_weights = FastGaussQuadrature.gausslegendre(npts)
    half = t_last / 2
    nodes = Tt[half * (u + one(u)) for u in raw_nodes]
    w_growth = Tt[half * w * exp(r * s)
                  for (w, s) in zip(raw_weights, nodes)]
    return ExportRiskTrajectory{Tt}(nodes, w_growth, convert(Tt, r),
        convert(Tt, t_last))
end

# Cumulative at-risk export person-time ∫₀^{t_k} e^{r u} S_det(t_k − u) du
# at each bin edge. Written as the exact growth integral minus the
# growth⊛F_det convolution: S_det = 1 − F_det does *not* vanish at lag 0
# (survival is 1 there), so convolving it directly and truncating the
# fixed grid at each edge would put a discontinuity at an interior node
# and the uniform rule would resolve it poorly. F_det *does* vanish at
# lag 0, so the convolved integrand is continuous and the shared grid is
# accurate; the e^{r u} part is integrated in closed form.
"""
At-risk export person-time ``\\int_0^{t_k} e^{r u}\\, S_{det}(t_k - u)\\,
du`` at each edge in `t_edges`, over an [`ExportRiskTrajectory`](@ref)
`traj` with a `Gamma` infection→detection delay `f_det`. Evaluated as the
closed-form growth integral minus the ``e^{r u} \\otimes F_{det}``
convolution. The delay-convolution form of
[`expected_exports_delay`](@ref); multiply by `p · q` for the expected
detected exports.
"""
function export_at_risk(traj::ExportRiskTrajectory,
        t_edges::AbstractVector, f_det::Gamma)
    α, θ = f_det.α, f_det.θ
    r = traj.r
    Tt = eltype(traj.w_growth)
    out = Vector{Tt}(undef, length(t_edges))
    @inbounds for k in eachindex(t_edges)
        tk = t_edges[k]
        tk <= zero(tk) && (out[k] = zero(Tt); continue)
        removed = zero(Tt)
        for j in eachindex(traj.nodes)
            δ = tk - traj.nodes[j]
            δ <= zero(δ) && continue
            removed += traj.w_growth[j] * _gamma_cdf(α, θ, δ)
        end
        ## Clamp to finite non-negative: at an extreme `r` proposal both
        ## the growth integral and `removed` overflow to `Inf`, so their
        ## difference is `NaN`; left unclamped it would reach `Poisson` and
        ## throw under AD (the prior path clamped the same way).
        val = _exp_cumulative_integral(r, zero(tk), tk) - removed
        out[k] = isfinite(val) ? max(val, zero(Tt)) : zero(Tt)
    end
    return out
end

"""
At-risk deaths-among-exports weight ``\\int_0^{t_k} e^{r u}\\,
S_{det}(t_k - u)\\, F_{death}(t_k - u)\\, du`` at each edge in `t_edges`,
over an [`ExportRiskTrajectory`](@ref) `traj` with `Gamma`
infection→detection (`f_det`) and infection→death (`f_death`) delays. The
delay-convolution form of [`expected_exports_deaths_delay`](@ref);
multiply by `CFR · p · q` for the expected export deaths.
"""
function export_death_at_risk(traj::ExportRiskTrajectory,
        t_edges::AbstractVector, f_det::Gamma, f_death::Gamma)
    αd, θd = f_det.α, f_det.θ
    αx, θx = f_death.α, f_death.θ
    Tt = eltype(traj.w_growth)
    out = Vector{Tt}(undef, length(t_edges))
    @inbounds for k in eachindex(t_edges)
        tk = t_edges[k]
        acc = zero(Tt)
        for j in eachindex(traj.nodes)
            δ = tk - traj.nodes[j]
            δ <= zero(δ) && continue
            acc += traj.w_growth[j] *
                   (one(δ) - _gamma_cdf(αd, θd, δ)) *
                   _gamma_cdf(αx, θx, δ)
        end
        ## Clamp to finite non-negative against `Inf`/`NaN` from extreme
        ## `r` proposals (see `export_at_risk`).
        out[k] = isfinite(acc) ? max(acc, zero(Tt)) : zero(Tt)
    end
    return out
end
