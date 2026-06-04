"""
Wrapper around `cdf(Gamma(α, θ), x)` as a 3-argument scalar
scalar function to attach reverse-mode rule.

Instead we attach our own analytic rrule:

```math
\\begin{aligned}
\\partial_x F      &= f(x; \\alpha), \\\\
\\partial_\\theta F &= -\\frac{x}{\\theta}\\, f(x; \\alpha), \\\\
\\partial_\\alpha F &= -\\psi(\\alpha)\\, P(\\alpha, y)
    + \\frac{1}{\\Gamma(\\alpha)} \\int_0^y t^{\\alpha-1} e^{-t} \\log t \\, dt,
\\end{aligned}
```

with `y = x/θ`, `P` the regularized lower incomplete gamma and `ψ` is the digamma function.
"""
_gamma_cdf(α, θ, x) = cdf(Gamma(α, θ), x)

"""
Integral of the Gamma(α, θ) CDF over `[0, x]`,

```math
\\int_0^x F(v;\\, α, θ)\\, dv = x\\, F(x;\\, α, θ) - α θ\\, F(x;\\, α + 1, θ),
```

from integration by parts (the second term is
``\\int_0^x v\\, f(v)\\,dv = α θ\\, F(x;\\, α + 1, θ)``). Backs the
background-tested volume of the confirmed-cases lab stream: a non-BVD
background arriving at constant rate, convolved against the lab-delay
CDF and right-truncated at the testing cut-off, reduces to this
integral. Two [`_gamma_cdf`](@ref) evaluations replace a full
quadrature, and the α / θ derivatives flow through the existing
`_gamma_cdf` rrule with no extra rule. Returns zero for `x ≤ 0`.
"""
function _gamma_cdf_integral(α, θ, x)
    x <= zero(x) &&
        return zero(float(promote_type(typeof(α), typeof(θ), typeof(x))))
    return x * _gamma_cdf(α, θ, x) - α * θ * _gamma_cdf(α + one(α), θ, x)
end

"""
Severe-first BVD share of the tested pool: a high→baseline exponential
relaxation. The laboratory tests the most-likely-BVD (severest / most
obvious) suspects first, so the BVD fraction `q` of the analysed batch
starts at the early severe-cluster share `q0` and decays toward the
broad-pool baseline `q∞` as testing widens:

```math
q(c) = q_\\infty + (q_0 - q_\\infty)\\, e^{-c / \\text{scale}},
```

with `c ≥ 0` the elapsed time since testing (reporting) onset and `scale`
the relaxation timescale (days). At `c = 0` (the first vintage) `q = q0`
(near 1 — the obvious-BVD cluster, so positivity `≈ s` identifies the
sensitivity); as `c → ∞` `q → q∞` and STAYS there (it does not overshoot
to zero, unlike a saturating coverage curve). The baseline `q∞` is tied to
the count-implied BVD composition of the suspect pool in
[`confirmed_cases_model`](@ref), pinning the plateau positivity to the
outbreak size; `q0` and `scale` are the [`test_selection_model`](@ref)
shape parameters. AD-safe: smooth in `q0`, `q∞`, `c`, `scale`.
"""
function severe_first_share(q0, q∞, c, scale)
    R = float(promote_type(typeof(q0), typeof(q∞), typeof(c), typeof(scale)))
    cc = max(convert(R, c), zero(R))
    sc = max(convert(R, scale), eps(R))
    w = exp(-cc / sc)
    q = convert(R, q∞) + (convert(R, q0) - convert(R, q∞)) * w
    return clamp(q, zero(R), one(R))
end

"""
Series sum of term derivatives for `∂_α P(α, z)`, using the
absolutely-convergent Kummer expansion

```math
\\begin{aligned}
P(\\alpha, z)                   &= z^{\\alpha}\\, e^{-z}
    \\sum_{n=0}^{\\infty} \\frac{z^n}{\\Gamma(\\alpha + n + 1)}, \\\\
\\partial_\\alpha P(\\alpha, z) &= \\log(z)\\, P(\\alpha, z)
    - z^{\\alpha}\\, e^{-z}
    \\sum_{n=0}^{\\infty} \\frac{\\psi(\\alpha + n + 1)\\, z^n}
                               {\\Gamma(\\alpha + n + 1)}.
\\end{aligned}
```

The digamma factor advances by the recurrence
`ψ(α + n + 1) = ψ(α + n) + 1/(α + n)`, so each iteration costs only a
division, an add and two multiplies — no per-iteration gamma or
digamma calls. Stan's `grad_reg_inc_gamma` uses this same series as
its small-`x` branch; the Julia formulation here is the direct port
from EpiAware/CensoredDistributions PR #250.
"""
function _grad_p_a_series(a, z; rtol = 1e-14, maxiter = 10_000)
    z <= zero(z) && return zero(a) * zero(z) #type promotion
    # avoid recalculating the same digamma values across iterations
    log_term0 = a * log(z) - z - loggamma(a + 1)
    term = exp(log_term0)
    ψ = digamma(a + 1)
    P = term
    S = term * ψ
    for n in 1:maxiter
        term *= z / (a + n)
        ψ += 1 / (a + n)
        P += term
        S += term * ψ
        # convergence check: both the P and S series must have converged to
        # ensure the final result is accurate to rtol.
        abs(term * ψ) <= rtol * abs(S) &&
            abs(term) <= rtol * abs(P) && break
    end
    return log(z) * P - S
end

"""
Compute the partial derivatives of the gamma CDF with respect to α, θ, and x.
"""
function _gamma_cdf_partials(α, θ, x)
    R = float(promote_type(typeof(α), typeof(θ), typeof(x)))
    y = x / θ
    y <= zero(y) && return zero(R), zero(R), zero(R)
    f = pdf(Gamma(α, θ), x)
    df_dx = f
    df_dθ = -y * f
    df_dα = _grad_p_a_series(α, y)
    return df_dα, df_dθ, df_dx
end

"""
ChainRulesCore rrule for the gamma CDF, using the above partials.

NB: the `NoTangent()` is for the function argument itself, which is not a callable/functor.
Using standard pullback convention where seed combines with the transposed Jacobian, although in this case
all the partials are real scalars.
"""
function ChainRulesCore.rrule(::typeof(_gamma_cdf),
        α::Real, θ::Real, x::Real)
    y = _gamma_cdf(α, θ, x)
    dα, dθ, dx = _gamma_cdf_partials(α, θ, x)
    function _gamma_cdf_pullback(ȳ)
        return (NoTangent(), dα' * ȳ, dθ' * ȳ, dx' * ȳ)
    end
    return y, _gamma_cdf_pullback
end

# Generate reverse-mode rules for Mooncake AD
Mooncake.@from_rrule Mooncake.DefaultCtx Tuple{
    typeof(_gamma_cdf), Float64, Float64, Float64}
