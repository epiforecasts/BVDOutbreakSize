# Expected-count integrals built on top of the shared `integrate`
# helpers: expected cumulative deaths, expected detected exports, and
# expected deaths-among-exports.

"""
Infection-to-onset rescale under exponential growth. The latent
cumulative *infections* trajectory is ``C(s) = e^{r s}``; convolving it
with the incubation-period density `f` gives the cumulative symptom
onsets

```math
\\int_0^{s} e^{r u}\\, f(s - u)\\, du = e^{r s}\\, M_f(-r),
```

a constant rescale of the same exponential by the incubation
moment-generating function ``M_f(-r)``. Onsets lag infections, so the
factor lies in ``(0, 1]`` and `onset_rescale(f, 0) = 1`. Applied as the
`onset_fraction` of the onset-driven observation submodels so every
downstream delay acts on onsets rather than infections.

This is a timing factor, not a biological proportion: at the cut-off it
is the share of infections that have *already* progressed to symptom
onset, the rest still being within their incubation period. It is not an
asymptomatic fraction — every infection becomes a symptomatic case
eventually. It is fully derived from the sampled incubation period and
the growth rate, so it carries no prior of its own.
"""
onset_rescale(delay_dist, r) = mgf(delay_dist, -r)

"""
Delay-convolved cumulative-incidence count by time `T` under
exponential growth at rate `r`:

```math
\\text{scale} \\cdot
    \\int_0^T e^{r s}\\, f(T - s)\\, ds,
```

with `f` the `delay_dist` density (per-event delay from incidence to
the observed event class — e.g. onset-to-death for deaths,
infection-to-report for suspected cases). The integrand returns zero
past `T` so the convolution support is respected. Integrated with the
clustered [`integrate`](@ref) so the quadrature nodes pack near `T`,
where `f(T − s)` has mass, over a window set by the delay scale (see
[`DELAY_SUPPORT_K`](@ref)). Uses [`DEATH_INTEGRAL_ALG`](@ref).
"""
function delay_convolution(scale, r, T, delay_dist; alg = DEATH_INTEGRAL_ALG)
    s_lo = zero(T)
    s_scale = _delay_scale(delay_dist)
    g = let r = r, T = T, delay_dist = delay_dist
        s -> begin
            d = T - s
            d <= 0 ? zero(r) : exp(r * s) * pdf(delay_dist, d)
        end
    end
    return scale * integrate(g, s_lo, T, s_scale; alg)
end

"""
[`delay_convolution`](@ref) specialised to the Gamma family, using the
analytic closed form

```math
\\text{scale} \\cdot
    \\int_0^T e^{r s}\\, f(T - s)\\, ds
  = \\text{scale} \\cdot e^{r T} \\cdot M(-r) \\cdot F(T (1 + \\theta r)),
```

where ``f`` is the Gamma(α, θ) density, ``M`` its moment-generating
function and ``F`` its CDF.
"""
function delay_convolution(scale, r, T, delay_dist::Gamma)
    α, θ = delay_dist.α, delay_dist.θ
    return scale * exp(r * T) * mgf(delay_dist, -r) *
           _gamma_cdf(α, θ, T * (1 + θ * r))
end

"""
Variant accepting an arbitrary cumulative-incidence trajectory
`cumulative(s)` in place of the ``e^{r s}`` proxy:

```math
\\int_0^T \\text{cumulative}(s)\\, f(T - s)\\, ds.
```

Used by the confirmed-cases likelihood, where `cumulative` is the
reported-cases submodel's BVD convolution as a function of time.
"""
function delay_convolution(cumulative::Function, T, delay_dist;
        alg = DEATH_INTEGRAL_ALG)
    scale = _delay_scale(delay_dist)
    g = let cumulative = cumulative, T = T, delay_dist = delay_dist
        s -> begin
            d = T - s
            d <= 0 ? zero(T) : cumulative(s) * pdf(delay_dist, d)
        end
    end
    return integrate(g, zero(T), T, scale; alg)
end

"""
Expected cumulative detected exports by elapsed time `t`, clamped to be
strictly positive and finite:

```math
\\mathbb{E}[\\text{exports}(t)] = p \\cdot q \\cdot
    \\int_{t-w}^{t} C(s)\\, ds,
```

with `cumulative` the trajectory ``C(s)``, `p` the detection
probability, `q` the per-capita travel rate, and `window` the detection
window ``w``. Backs both the exports count likelihood (evaluated at
`t = T`) and the first-export-detection survival term (evaluated at an
earlier `t`).

Pass `r` (the exponential growth rate) to use the exact closed form
``\\int_a^b e^{r s}\\,ds = (e^{r b} - e^{r a})/r`` instead of numerical
quadrature — faster and with a cleaner gradient under AD, used by the
model where the growth trajectory is the exponential default. With
`r = nothing` (the default) the window integral is evaluated
numerically via [`CUMULATIVE_INTEGRAL_ALG`](@ref) for an arbitrary
`cumulative`.
"""
function expected_exports(cumulative, p, q, t, window;
        r = nothing, alg = CUMULATIVE_INTEGRAL_ALG)
    window_start = max(t - window, zero(t))
    integral = r === nothing ?
               integrate_cumulative(cumulative, window_start, t; alg) :
               _exp_cumulative_integral(r, window_start, t)
    raw = p * q * integral
    return isfinite(raw) ? max(raw, eps(typeof(raw))) : eps(typeof(raw))
end

## Exact ∫_a^b exp(r·s) ds for an exponential cumulative-incidence
## trajectory. Written via `expm1` so it stays accurate as r → 0;
## algebraically equal to (exp(r·b) − exp(r·a)) / r.
_exp_cumulative_integral(r, a, b) = exp(r * a) * expm1(r * (b - a)) / r

"""
Expected cumulative deaths among detected exports by elapsed time `t`,
clamped to be strictly positive and finite:

```math
\\mathbb{E}[D_{\\text{uganda}}(t)] = \\mathrm{CFR} \\cdot
    p \\cdot q \\cdot
    \\int_{t-w}^{t} C(s)\\, F_d(t - s)\\, ds,
```

with `cumulative` the trajectory ``C(s)``, `delay_dist` the onset-to-death
distribution (CDF ``F_d``), `p` the detection probability, `q` the
per-capita travel rate, and `window` the detection window ``w``. Backs
the binned-Poisson export-death likelihood, evaluated at the daily bin
edges. Uses [`CUMULATIVE_INTEGRAL_ALG`](@ref).
"""
function expected_exports_deaths(cumulative, delay_dist, CFR, p, q,
        t, window; alg = CUMULATIVE_INTEGRAL_ALG)
    window_start = max(t - window, zero(t))
    integral = integrate_exports_deaths(
        cumulative, delay_dist, window_start, t, t; alg)
    raw = CFR * p * q * integral
    return isfinite(raw) ? max(raw, eps(typeof(raw))) : eps(typeof(raw))
end
