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
Moment-matched Gamma approximation to the sum of two independent Gamma
delays `a` and `b`. The convolution `a ⊕ b` is not Gamma in general, so
its mean ``\\mu = \\mu_a + \\mu_b`` and variance ``v = v_a + v_b`` are
matched to a single `Gamma(μ²/v, v/μ)`. Used to build the Uganda-export
infection→detection delay as incubation ⊕ onset-to-report, so the
at-risk export window runs from infection (capturing pre-symptomatic
travel) through to detection abroad. Plain arithmetic on the gamma
moments keeps the result differentiable under AD.
"""
function combined_delay(a::Gamma, b::Gamma)
    m = mean(a) + mean(b)
    v = std(a)^2 + std(b)^2
    return Gamma(m^2 / v, v / m)
end

"""
Expected cumulative detected exports by elapsed time `t` under an
explicit infection→detection delay, clamped to be strictly positive and
finite. The delay-survival generalisation of [`expected_exports`](@ref):
instead of a rectangular detection window of fixed width `w`, a case
stays at risk of being exported and detected abroad only until the
infection→detection delay `f_det` has elapsed. The exports stream is
travel-gated, so the at-risk clock starts at infection (the traveller
moves during incubation, pre-symptomatic); `f_det` is therefore the full
infection→detection delay (incubation ⊕ onset-to-report, see
[`combined_delay`](@ref)), not an onset-to-report delay rescaled by the
incubation moment. Writing the cumulative infections as ``C(s) = e^{r s}``
and the infections that have already completed infection→detection by `s`
as the convolution
``\\text{detected}(s) = \\int_0^s e^{r u}\\, f_{det}(s-u)\\,du``, the
at-risk prevalence is the difference

```math
\\text{prevalence}(s) = e^{r s} - \\text{detected}(s),
```

and the expected detected exports accumulate the per-day per-capita
travel rate `q` over the at-risk person-time

```math
\\mathbb{E}[\\text{exports}(t)] = p\\, q
    \\int_0^t \\text{prevalence}(s)\\, ds
  = p\\, q
    \\left(\\int_0^t e^{r s}\\, ds - \\int_0^t \\text{detected}(s)\\, ds\\right).
```

`r` is the exponential growth rate, `p` the detection probability and
`q` the per-day per-capita travel rate (so the result is in cases). The
``\\int_0^t e^{r s}\\, ds`` term is the exact closed form
[`_exp_cumulative_integral`](@ref); the inner convolution uses the Gamma
closed form of [`delay_convolution`](@ref), so only the outer
``\\int_0^t \\text{detected}(s)\\, ds`` is quadrature and the path stays
AD-friendly. As the infection→detection delay collapses to a point mass
at `w` the survival becomes the top-hat ``1\\{t-s<w\\}`` and this reduces
exactly to the McCabe window form [`expected_exports`](@ref).
"""
function expected_exports_delay(r, p, q, t, f_det;
        alg = CUMULATIVE_INTEGRAL_ALG)
    cum_integral = _exp_cumulative_integral(r, zero(t), t)
    removed = integrate(s -> delay_convolution(one(t), r, s, f_det),
        zero(t), t; alg)
    raw = p * q * (cum_integral - removed)
    return isfinite(raw) ? max(raw, eps(typeof(raw))) : eps(typeof(raw))
end

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

"""
Expected cumulative deaths among detected exports by elapsed time `t`
under an explicit infection→detection delay, clamped to be strictly
positive and finite. The delay-survival generalisation of
[`expected_exports_deaths`](@ref): McCabe et al. integrate the
onset-to-death CDF over the top-hat detection window `[t-w, t]`, i.e. a
top-hat detection survival ``S(t-s) = 1\\{t-s<w\\}``. Here the indicator
is replaced by the infection→detection survival
``\\overline{F}_{det}(t-s)``, so

```math
\\mathbb{E}[D_{\\text{uganda}}(t)] = \\mathrm{CFR}\\, p\\, q\\,
    \\int_0^t C(s)\\, \\overline{F}_{det}(t-s)\\, F_{death}(t-s)\\, ds,
```

with `cumulative` the infection trajectory ``C(s) = e^{r s}``, `f_det`
the infection→detection Gamma delay (incubation ⊕ onset-to-report, see
[`combined_delay`](@ref); survival ``\\overline{F}_{det}``), `f_death`
the onset-to-death Gamma delay (CDF ``F_{death}``), `CFR` the
case-fatality ratio, `p` the detection probability and `q` the per-day
per-capita travel rate (so the result is in cases). Because the at-risk
clock starts at infection, the detection survival is 1 at age 0 (a
just-infected traveller is certainly not yet detected). The survival and
the death CDF are both evaluated through the Gamma closed form
[`_gamma_cdf`](@ref), which carries the reverse-mode rule for the shape
derivative, so the integrand is closed-form and only the outer integral
is quadrature (clustered near `t` where the integrand has mass). As the
infection→detection delay collapses to a point mass at `w` the survival
becomes the top-hat and this reduces exactly to
[`expected_exports_deaths`](@ref). Uses [`DEATH_INTEGRAL_ALG`](@ref).
"""
function expected_exports_deaths_delay(cumulative, f_det::Gamma,
        f_death::Gamma, CFR, p, q, t; alg = DEATH_INTEGRAL_ALG)
    scale = _delay_scale(f_det) + _delay_scale(f_death)
    g = let cumulative = cumulative, t = t, f_det = f_det, f_death = f_death
        s -> begin
            a = t - s
            a <= zero(a) ? zero(a) :
            cumulative(s) *
            (one(a) - _gamma_cdf(f_det.α, f_det.θ, a)) *
            _gamma_cdf(f_death.α, f_death.θ, a)
        end
    end
    integral = integrate(g, zero(t), t, scale; alg)
    raw = CFR * p * q * integral
    return isfinite(raw) ? max(raw, eps(typeof(raw))) : eps(typeof(raw))
end
