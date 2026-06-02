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
Saturating ramp non-BVD background, anchored to surveillance/reporting
onset. Carries the baseline rate `λ0`, the ramp amplitude `Δλ`, the
scale-up timescale `scale` and the reporting-onset elapsed time
`t_report` (the latent time at which case-finding begins; before it
there is no surveillance so no background suspects accrue). The per-day
non-BVD suspected-case arrival rate, with the ramp clock starting at
`t_report`,

```math
\\lambda_{bg}(t) = \\begin{cases}
0, & t \\le t_{report}, \\\\
\\lambda_0 + \\Delta\\lambda\\,
    \\bigl(1 - e^{-(t - t_{report})/\\text{scale}}\\bigr),
    & t > t_{report},
\\end{cases}
```

rises from `λ0` towards `λ0 + Δλ` over the surveillance scale-up window
`scale` once reporting has begun. Its cumulative is the integral from
`t_report`,

```math
\\mu_{bg}(t) = \\int_{t_{report}}^{t} \\lambda_{bg}(u)\\, du
  = \\lambda_0\\,\\Delta t + \\Delta\\lambda\\,
    \\bigl(\\Delta t - \\text{scale}
    + \\text{scale}\\, e^{-\\Delta t/\\text{scale}}\\bigr),
\\quad \\Delta t = t - t_{report},
```

`μ_bg(t) = 0` for `t ≤ t_report`. Anchoring to reporting onset (rather
than seeding) lets the ramp still be climbing across the late lab
vintages, so the broadening suspected-case definition pulls the observed
per-test positivity *down* there; a seeding-anchored ramp has long since
saturated by the cut-off and cannot. `Δλ = 0` *and* `t_report = 0`
recover the constant-rate background `λ_bg(t) = λ0`, `μ_bg(t) = λ0·t`,
exactly. Construct with [`background_ramp`](@ref); the rate is
[`bg_rate`](@ref) and the cumulative [`bg_cumulative`](@ref).
"""
struct BackgroundRamp{T}
    λ0::T
    Δλ::T
    scale::T
    t_report::T
end

"""
Build a [`BackgroundRamp`](@ref) from the baseline rate `λ0`, ramp
amplitude `Δλ`, fixed scale-up timescale `scale` and reporting-onset
elapsed time `t_report` (default `0`, recovering the seeding-anchored
clock), promoting to a common element type so AD tangents flow through
`λ0`, `Δλ` and `t_report`.
"""
function background_ramp(λ0, Δλ, scale, t_report = zero(scale))
    T = float(promote_type(typeof(λ0), typeof(Δλ), typeof(scale),
        typeof(t_report)))
    return BackgroundRamp{T}(convert(T, λ0), convert(T, Δλ),
        convert(T, scale), convert(T, t_report))
end

"""
Instantaneous non-BVD background rate `λ_bg(t)` of a
[`BackgroundRamp`](@ref). Zero before reporting onset `t_report`; for
`t > t_report` returns `λ0` exactly when `Δλ = 0`.
"""
function bg_rate(b::BackgroundRamp, t)
    Δt = t - b.t_report
    Δt <= zero(Δt) && return zero(b.λ0 * Δt)
    return b.λ0 + b.Δλ * (one(Δt) - exp(-Δt / b.scale))
end

"""
Cumulative non-BVD background `μ_bg(t)` of a [`BackgroundRamp`](@ref),
the integral of [`bg_rate`](@ref) from the reporting onset `t_report` to
`t`. `μ_bg(t) = 0` for `t ≤ t_report`; with `Δλ = 0` it reduces to
`λ0·(t − t_report)`.
"""
function bg_cumulative(b::BackgroundRamp, t)
    Δt = t - b.t_report
    Δt <= zero(Δt) && return zero(b.λ0 * Δt)
    return b.λ0 * Δt +
           b.Δλ * (Δt - b.scale + b.scale * exp(-Δt / b.scale))
end

"""
Cumulative background *tested* volume of a [`BackgroundRamp`](@ref):
the time-varying background arrivals (from reporting onset `t_report`)
convolved against the lab-delay CDF `F_lab` and right-truncated at the
testing cut-off `t`,

```math
\\int_{t_{report}}^{t} \\lambda_{bg}(u)\\, F_{lab}(t - u)\\, du
  = (\\lambda_0 + \\Delta\\lambda)\\, G - \\Delta\\lambda\\, H,
```

where, writing the substitution `v = u - t_report` and
`Δt = t - t_report`,
``G = \\int_{0}^{Δt} F_{lab}(Δt - v)\\,dv`` is the closed-form
[`_gamma_cdf_integral`](@ref) at `Δt`, and
``H = \\int_{0}^{Δt} e^{-v/\\text{scale}}\\, F_{lab}(Δt - v)\\,dv`` is
the ramp's exponential weighting, evaluated numerically with the shared
[`integrate`](@ref) quadrature. With `Δλ = 0` only the closed-form term
`λ0·G` remains, reproducing the constant-rate background tested volume
exactly when `t_report = 0`. Returns zero for `t ≤ t_report`.
"""
function bg_tested_integral(b::BackgroundRamp{T}, α, θ, t;
        alg = DEATH_INTEGRAL_ALG) where {T}
    R = float(promote_type(T, typeof(α), typeof(θ), typeof(t)))
    Δt = t - b.t_report
    Δt <= zero(Δt) && return zero(R)
    G = _gamma_cdf_integral(α, θ, Δt)
    iszero(b.Δλ) && return convert(R, b.λ0 * G)
    ## H = ∫₀^{Δt} e^{-v/scale} F_lab(Δt-v) dv, the ramp's exponential
    ## weighting of the lab-delay CDF (clock measured from reporting
    ## onset). No closed form, so reuse the shared quadrature; the
    ## integrand peaks near v = 0 (small delay argument Δt - v ≈ Δt) and
    ## decays, so the uniform rule on [0, Δt] is adequate.
    H = integrate(zero(Δt), Δt; alg) do v
        exp(-v / b.scale) * _gamma_cdf(α, θ, Δt - v)
    end
    return convert(R, (b.λ0 + b.Δλ) * G - b.Δλ * H)
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
