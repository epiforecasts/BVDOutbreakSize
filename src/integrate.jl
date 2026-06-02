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
