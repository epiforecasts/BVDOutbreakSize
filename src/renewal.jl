# Discrete-time renewal primitives. Pure, allocation-light, and
# AD-transparent (output element types are promoted from the inputs) so
# they differentiate cleanly under Mooncake inside a Turing model. Delays
# are applied by daily convolution.

"""
NaN / Inf-safe positive rate. The renewal recursion can transiently
overflow on extreme NUTS warmup proposals, giving a non-finite expected
count. A plain `max(x, eps)` would propagate the NaN
(`max(NaN, eps) = NaN`) and trip the Poisson / NegativeBinomial domain
check.
"""
@inline safe_rate(x) = _safe_rate_on(x) ? x : eps(typeof(x))

## Whether `safe_rate` passes `x` through rather than flooring it, which is
## where its slope is one.
@inline _safe_rate_on(x) = isfinite(x) && x > eps(typeof(x))

"""
LogNormal with the given `mean` and standard deviation `sd`, by moment
matching `var = mean^2 (exp(σ^2) − 1)`. The inputs are passed through
[`safe_rate`](@ref) first so a NaN-prone warmup proposal cannot push
`σ = sqrt(log1p(·))` into NaN territory and trip the LogNormal domain
check. Used by every delay submodel so a delay is parameterised by its
mean and SD rather than the log-scale parameters.
"""
function lognormal_meansd(mean, sd)
    m = safe_rate(mean)
    s = safe_rate(sd)
    σ2 = log1p((s / m)^2)
    μ = log(m) - σ2 / 2
    return LogNormal(μ, sqrt(σ2))
end

"""
Daily probability mass function for the continuous delay `dist` over lags
`0, 1, …, nmax`, discretised by double interval censoring (uniform primary
event over a one-day window, then unit-interval censoring of the secondary
event, truncated at `nmax`), as `CensoredDistributions.double_interval_censored`
defines it. The truncation holds the CDF at one from `nmax` on, so the lag
`nmax` entry is zero. For a LogNormal or Gamma delay the CDF differentiates
cleanly under Mooncake, so this is AD-safe. Extreme warmup proposals that
drive the total to a non-finite or zero value fall back to a uniform PMF, so
the downstream convolution stays finite (the proposal is still rejected
through its low log-likelihood). Returns a vector whose element type follows
the delay parameters.
"""
function discretise_censored(dist, nmax::Integer)
    return _pmf_from_cdfs(_primary_censored_cdfs(dist, nmax), dist, nmax)
end

## Delays with a closed-form primary-censored CDF under a uniform primary.
const _AnalyticalDelay = Union{Gamma, LogNormal, Distributions.Weibull}

## Primary-censored CDF `F₊(b)` at the boundaries `b = 1, …, nmax`, with a
## `Uniform(0, 1)` primary event.
##
## For the analytical delays `F₊(b) = H(b) − H(b − 1)`, where
## `H(t) = t·F(t) − M(t)` and `M` is the delay's partial first moment. The
## analytical CDF under a `Uniform(0, t)` primary, evaluated at `t`, is
## `H(t) / t`, so each `H` costs one delay-CDF endpoint and each boundary
## reuses its neighbour's.
function _primary_censored_cdfs(dist::_AnalyticalDelay, nmax::Integer)
    solver = AnalyticalSolver()
    pc = Vector{float(Distributions.partype(dist))}(undef, nmax)
    H_prev = zero(eltype(pc))
    for i in 1:nmax
        t = float(i)
        H = t * primarycensored_cdf(dist, Uniform(zero(t), t), t, solver)
        pc[i] = H - H_prev
        H_prev = H
    end
    return pc
end

function _primary_censored_cdfs(dist, nmax::Integer)
    pc = primary_censored(dist, Uniform(0.0, 1.0))
    return [cdf(pc, t) for t in 1.0:nmax]
end

## Lag masses from the boundary CDFs `F₊(1), …, F₊(nmax)`: lag `d < nmax`
## has `F₊(d + 1) − F₊(d)` (with `F₊(0) = 0`) and lag `nmax` has zero, then
## the masses are normalised by their sum, which is the truncation at `nmax`.
function _pmf_from_cdfs(pc, dist, nmax::Integer)
    z0 = zero(eltype(pc))
    raw = Vector{eltype(pc)}(undef, nmax + 1)
    pc_prev = z0
    for i in 1:nmax
        raw[i] = max(pc[i] - pc_prev, z0)
        pc_prev = pc[i]
    end
    raw[nmax + 1] = z0
    s = sum(raw)
    if !isfinite(s) || s <= zero(s)
        z = zero(pdf(dist, oneunit(float(nmax))))
        return fill(one(z) / (nmax + 1), nmax + 1)
    end
    return raw ./ s
end

"""
    cdf_nmax(dist; q = 0.98, cap = 120, minlag = 5)

Maximum lag for discretising a delay `dist`: the smallest integer covering
`q` of its CDF (the `q`th quantile rounded up), clamped to `[minlag, cap]`.
Sizing the truncation by the distribution keeps a consistent tail mass across
every delay. It is evaluated once outside the Turing model when each delay
submodel is constructed, so the PMF length is fixed and AD-safe.
"""
function cdf_nmax(dist; q::Real = 0.98, cap::Integer = 120, minlag::Integer = 5)
    return clamp(ceil(Int, quantile(dist, q)), minlag, cap)
end

"""
Exponential growth rate `r` implied by a reproduction number `R` and a
generation-interval PMF `g` (indexed from lag 1): the root of the
Euler–Lotka identity `R · Σ_s g_s e^{−r s} = 1`. Starts from the
small-`r` approximation `r ≈ (R − 1) / (R · ḡ)` with `ḡ` the mean
generation time, or below `R = 1` from `log(R) / ḡ`, then takes Newton
steps on `log(R Σ_s g_s e^{−r s}) = 0` until a step is below `1e-12`, at
most 20. That function is convex and falling in `r`, and both starts sit
below the root, so Newton climbs to it without overshooting. The second
start is closer as `R` falls, and the log form stays close to linear there,
so a near-zero `R` from a depleted pool gives a finite rate. `R = 0` has no
finite root and gives `-Inf`. Its Mooncake rule differentiates the root
through the implicit function theorem. Mirrors the `R_to_r` seeding helper
in EpiAware.jl and the implied-growth initialisation in the EpiNow2 Stan
model.
"""
function euler_lotka_r(R, g::AbstractVector)
    Tp = promote_type(typeof(float(R)), eltype(g))
    ḡ = zero(Tp)
    @inbounds for i in eachindex(g)
        ḡ += g[i] * i
    end
    R > zero(R) || return Tp(-Inf)
    r = R < one(R) ? log(R) / ḡ : (R - one(R)) / (R * ḡ)
    for _ in 1:20
        G, dG = euler_lotka_sums(r, g)
        ## f(r) = log(R·G), f'(r) = −dG / G.
        Δ = log(R * G) * G / dG
        r += Δ
        abs(Δ) < 1.0e-12 && break
    end
    return r
end

## `G = Σ_s g_s e^{−r s}` and `dG = Σ_s s g_s e^{−r s}`, with `e^{−r s}` a
## running product of one `exp`.
function euler_lotka_sums(r, g::AbstractVector)
    Tp = promote_type(typeof(float(r)), eltype(g))
    q = exp(-r)
    e = one(Tp)
    G = zero(Tp)
    dG = zero(Tp)
    @inbounds for i in eachindex(g)
        e *= q
        G += g[i] * e
        dG += g[i] * i * e
    end
    return G, dG
end

"""
Reproduction number `R` implied by an exponential growth rate `r` and a
generation-interval PMF `g` (indexed from lag 1), the forward Euler–Lotka
relation `R = 1 / Σ_s g_s e^{−r s}`. The inverse of [`euler_lotka_r`](@ref),
so a prior can be placed on the growth rate and the reproduction number
derived from it under the model's generation interval. `e^{−r s}` is a
running product of one `exp`, so the gradient tapes one `exp`, not one per
lag.
"""
function r_to_R0(r, g::AbstractVector)
    Tp = promote_type(typeof(float(r)), eltype(g))
    q = exp(-r)
    e = one(Tp)
    G = zero(Tp)
    @inbounds for i in eachindex(g)
        e *= q
        G += g[i] * e
    end
    return one(Tp) / G
end

"""
Doubling time `log(2) / r` implied by an exponential growth rate `r`.
Returns a non-finite value as `r` crosses zero, matching the limit of an
unbounded doubling time at zero growth.
"""
doubling_time(r) = log(oftype(float(r), 2)) / r

"""
Seed the first `len` days of the infection trajectory as exponential
growth `I_t = I0 · e^{r (t − len)}` at the implied growth rate `r` (see
[`euler_lotka_r`](@ref)), so the seeding window is pinned at `I0` on
its last day and tails off backwards. This is the initialisation used by
EpiNow2 and EpiAware.jl. Placing the whole seed on a single day instead
injects a transient the renewal recursion has to relax away from. Returns a
length-`len` vector whose element type follows `I0` and `r`.
"""
function seed_infections(I0, r, len::Integer)
    Tp = promote_type(typeof(float(I0)), typeof(float(r)))
    seed = Vector{Tp}(undef, len)
    @inbounds for j in 1:len
        seed[j] = I0 * exp(r * (j - len))
    end
    return seed
end

"""
Renewal-start seed magnitude for the two-phase renewal: the daily infection
incidence the analytic cryptic phase reaches on the renewal-start day.

It is an incidence, not a cumulative count. The cryptic total over the
`renewal_start` grid days is larger by more than an order of magnitude, so
reading the seed as cumulative would shift `m` by several generations.

The cryptic phase runs from the origin to the renewal start (≈ the genetic
TMRCA day, off the renewal grid), spanning `T = m · G` days for `m` generations
at mean generation interval `G`. A single daily infection at the origin grows
over it at the cryptic rate `r` to

```math
\\text{seed\\_at\\_renewal\\_start} = C_T = e^{r T},
```

which the renewal grows forward under `R_t`, so the realised cut-off size stays
data-driven while the prior fixes only the renewal-start scale.

The magnitude is referenced to the origin, not the cut-off. A cut-off-referenced
size would put `r` into the seed and the renewal growth in opposing directions,
cancelling for a fixed realised size and leaving a flat ridge along which `R0`
could slide. Referenced to the origin they compound instead, so `r` does enter
the seed.

`C_T_prior` is returned unchanged. The helper names the quantity at the
seeding call site.
"""
@inline function seed_at_renewal_start(C_T_prior)
    return C_T_prior
end

"""
Daily latent infections from the renewal equation with weak susceptible
depletion. Each day's renewal force `R_t Σ_{s ≥ 1} I_{t−s} g_s`, with
generation-interval PMF `g` (indexed from lag 1) and per-day reproduction
numbers `Rt` (length `n`), is taken as a rate on a pool of `N`:

```math
x_t = \\frac{R_t}{N} \\sum_{s \\ge 1} I_{t-s} g_s, \\qquad
I_t = S_{t-1}\\left(1 - e^{-x_t}\\right), \\qquad
S_t = S_{t-1} e^{-x_t}.
```

While the pool is large against the outbreak, `I_t ≈ (S_{t−1} / N) R_t f_t`,
the renewal scaled by the susceptible fraction. `I_t` never exceeds the pool
left, so an overflowing force takes the pool rather than giving `Inf`.
`N` is a population size, not an estimate of who can be reached. The models
pass census counts ([`PROVINCE_POPULATIONS`](@ref)), and at that scale it is
a light bound.
`Rt` is therefore the reproduction number in a fully susceptible population.

A pre-computed `seed` of length `L < n` fills the first `L` days (see
[`seed_infections`](@ref)) and is drawn from the pool before the recursion
runs for days `L+1 … n`. Returns the length-`n` infection trajectory. The
output element type is promoted from `Rt`, `g`, `seed` and `N`.

!!! note "Multi-patch analogue"
    See [`patch_infections`](@ref) for the meta-population extension
    with between-patch importation.
"""
function renewal_infections(
        Rt::AbstractVector, g::AbstractVector,
        seed::AbstractVector, N::Real
    )
    return renewal_infections_with_state(Rt, g, seed, N).infections
end

"""
The renewal trajectory with the state it was built from, as
`(; infections, force, susceptible)`: the per-day force
`Σ_{s ≥ 1} I_{t−s} g_s` and the pool left at the end of each day, which is
the pool after the seed on the seeded days. [`renewal_infections`](@ref)
returns the first; the derivative rule needs the others, which it would
otherwise have to rebuild from a copy of this loop.

Each day's force is one `dot` of the most recent infections with the
generation interval reversed, a single BLAS call on float arrays.
"""
function renewal_infections_with_state(
        Rt::AbstractVector, g::AbstractVector,
        seed::AbstractVector, N::Real
    )
    n = length(Rt)
    L = length(seed)
    G = length(g)
    Tp = promote_type(eltype(Rt), eltype(g), eltype(seed), typeof(float(N)))
    I = zeros(Tp, n)
    force = zeros(Tp, n)
    S = zeros(Tp, n)
    @inbounds for j in 1:min(L, n)
        I[j] = seed[j]
    end
    pool = convert(Tp, _pool_after_seed(N, view(I, 1:min(L, n))))
    @inbounds for j in 1:min(L, n)
        S[j] = pool
    end
    rg = reverse(g)
    for t in (L + 1):n
        k = min(t - 1, G)
        f = dot(view(rg, (G - k + 1):G), view(I, (t - k):(t - 1)))
        force[t] = f
        x = Rt[t] * f / N
        I[t] = -pool * expm1(-x)
        pool *= exp(-x)
        S[t] = pool
    end
    return (; infections = I, force, susceptible = S)
end

"""
    pool_fraction(cumulative, N)

Share of the pool `N` left susceptible at the end of each day, `1 − C_t / N`
for the cumulative infections `C_t`, seed included, floored at zero. Exact
for the depletion in [`renewal_infections`](@ref), where each day's
infections are what leaves the pool and a seed larger than `N` leaves none.
"""
pool_fraction(cumulative::AbstractVector, N::Real) =
    max.(1 .- cumulative ./ N, 0)

"""
    adjusted_rt(Rt, fraction, days)

The reproduction number net of depletion on each of `days`, `R_t` times the
susceptible `fraction` at the end of the day before
([`pool_fraction`](@ref)). This is the number of infections each
infection causes given the pool left, where `Rt` is the number in a fully
susceptible population. Day one takes the full pool.
"""
function adjusted_rt(Rt::AbstractVector, fraction::AbstractVector, days)
    return [
        @inbounds(Rt[t] * (t > 1 ? fraction[t - 1] : one(eltype(fraction))))
            for t in days
    ]
end

## The pool left once the seed is drawn from `N`. The seed can itself
## overflow on extreme warmup proposals, so the pool is floored at zero.
@inline function _pool_after_seed(N, seed)
    left = N - sum(seed)
    return left > zero(left) ? left : zero(left)
end

## --- Multi-patch (meta-population) renewal primitives --------------------

## Importation intensity of origin `q` on day `t`. A scalar applies to every
## origin and every day; a matrix carries one level per origin over time.
@inline _eps(e::Real, q::Integer, t::Integer) = e
@inline _eps(e::AbstractMatrix, q::Integer, t::Integer) = @inbounds e[q, t]

## Renewal force of patch `p` on day `t`, `Σ_{s ≥ 1} I[p, t − s] g[s]` over
## the lags inside the grid. Shared by `patch_infections` and its rule.
## `@simd` lets the sum reassociate, so it can differ from a sequential sum
## in the last bits.
@inline function _patch_force(I::AbstractMatrix, g::AbstractVector, p, t)
    f = zero(eltype(I))
    @inbounds @simd for s in 1:min(t - 1, length(g))
        f += I[p, t - s] * g[s]
    end
    return f
end

"""
    patch_infections(Rt_matrix, g, seeds_matrix, importation_kernel, epsilon, N)

Multi-patch (meta-population) renewal with between-patch importation and
weak susceptible depletion. Each patch `p` generates infections from its
own renewal force on a shared daily grid, and importation relocates a share
of each day's generated infections through a kernel `K`:

```math
G_{p,t} = R_{p,t} \\sum_{s \\ge 1} I_{p,t-s}\\, g_s, \\qquad
y_{p,t} = \\Bigl(1 - \\varepsilon_{p,t} \\sum_{q \\ne p} K_{q,p}\\Bigr) G_{p,t}
    + \\sum_{q \\ne p} \\varepsilon_{q,t} K_{p,q}\\, G_{q,t}.
```

`y_{p,t}` is taken as a rate on the patch's own pool of `N[p]` as in
[`renewal_infections`](@ref): `I_{p,t} = S_{p,t−1}(1 − e^{−y_{p,t} / N_p})`
and `S_{p,t} = S_{p,t−1} e^{−y_{p,t} / N_p}`.

# Arguments

- `Rt_matrix`: `n_patches x n_days` matrix whose `[p, t]` entry is the
  reproduction number in patch `p` on day `t`. Each row is one patch's
  daily `R_t` trajectory.
- `g`: shared generation-interval PMF (indexed from lag 1, so `g[1]` is
  the probability of a one-day generation interval). Same for all patches.
- `seeds_matrix`: `n_patches x L` matrix whose `[p, :]` row is the
  pre-computed seed infection trajectory for patch `p` (see
  [`seed_infections`](@ref)). The seed fills days `1 ... L` and the renewal
  recursion begins on day `L+1`.
- `importation_kernel`: `n_patches x n_patches` matrix `K` where
  `K[p, q]` is the share of patch `q`'s transmission that lands in patch
  `p` rather than at home. The diagonal should be zero, and each column's
  off-diagonal sum times `epsilon` must be at most one, so a patch cannot
  export more transmission than it generates. Both hold for
  [`province_importation_kernel`](@ref) at any `epsilon` in `[0, 1]`.
- `epsilon`: importation intensity of each origin, a scalar for every
  origin and day or an `n_patches x n_days` matrix. Importation
  is a transfer on the day it happens, so the origin patch is debited exactly
  what the destination patches are credited and coupling is never a source of
  infections. It is not conserved across days. The destination grows at its
  own reproduction number, so relocating infections from a fast patch to a
  slow one lowers the national total and the reverse raises it. With one
  shared reproduction number the transfer cancels exactly. Depletion acts
  on the destination after the transfer, so where a pool binds the
  destination takes in less than the origin sent.
- `N`: population of each patch, the pool it depletes. The seed is drawn
  from it first.

# Returns

`(; infections, importation)`. `infections` is a matrix of shape
`(n_patches, n_days)` where row `p` is the daily infection trajectory for
patch `p`. The first `L` days are copied from `seeds_matrix` and the
remaining days are the renewal recursion with importation. `importation` is
the matching matrix of infections each patch received from the others, the
arrivals term alone before depletion rather than the net of arrivals and
departures. The element type is promoted from all input types.
AD-transparent under Mooncake.
"""
function patch_infections(
        Rt_matrix::AbstractMatrix, g::AbstractVector,
        seeds_matrix::AbstractMatrix, importation_kernel::AbstractMatrix,
        epsilon::Union{Real, AbstractMatrix}, N::AbstractVector
    )
    st = patch_infections_with_state(
        Rt_matrix, g, seeds_matrix, importation_kernel, epsilon, N
    )
    return (; st.infections, st.importation)
end

"""
The [`patch_infections`](@ref) trajectories with the state they were built
from: `susceptible`, the pool each patch has left at the end of each day
(the pool after its seed on the seeded days), and `rate`, each day's
depletion rate `y_{p,t} / N_p`. The derivative rule reads both.
"""
function patch_infections_with_state(
        Rt_matrix::AbstractMatrix, g::AbstractVector,
        seeds_matrix::AbstractMatrix, importation_kernel::AbstractMatrix,
        epsilon::Union{Real, AbstractMatrix}, N::AbstractVector
    )
    np, n = size(Rt_matrix)
    L = size(seeds_matrix, 2)
    Tp = promote_type(
        eltype(Rt_matrix), eltype(g), eltype(seeds_matrix),
        eltype(importation_kernel),
        epsilon isa Real ? typeof(float(epsilon)) : eltype(epsilon),
        float(eltype(N))
    )
    I = zeros(Tp, np, n)
    imports = zeros(Tp, np, n)
    S = zeros(Tp, np, n)
    rate = zeros(Tp, np, n)
    outflow = _patch_outflow(Tp, importation_kernel, np)
    pool = zeros(Tp, np)
    @inbounds for p in 1:np
        for j in 1:min(L, n)
            I[p, j] = seeds_matrix[p, j]
        end
        pool[p] = _pool_after_seed(N[p], view(I, p, 1:min(L, n)))
        for j in 1:min(L, n)
            S[p, j] = pool[p]
        end
    end
    gen = zeros(Tp, np)
    @inbounds for t in (L + 1):n
        ## What each patch generates today from its own renewal force.
        for p in 1:np
            gen[p] = Rt_matrix[p, t] * _patch_force(I, g, p, t)
        end
        ## Importation redistributes transmission rather than adding to it. A
        ## fraction `epsilon * K[p, q]` of what `q` generates is realised in
        ## `p` instead of at home, so `q` is debited exactly what the
        ## destinations are credited. The intensity belongs to the origin, so
        ## `p` is debited at its own rate and credited at each sender's.
        for p in 1:np
            arrivals = zero(Tp)
            for q in 1:np
                q == p && continue
                arrivals += _eps(epsilon, q, t) *
                    importation_kernel[p, q] * gen[q]
            end
            imports[p, t] = arrivals
            y = (one(Tp) - _eps(epsilon, p, t) * outflow[p]) * gen[p] +
                arrivals
            x = y / N[p]
            rate[p, t] = x
            I[p, t] = -pool[p] * expm1(-x)
            pool[p] *= exp(-x)
            S[p, t] = pool[p]
        end
    end
    return (; infections = I, importation = imports, susceptible = S, rate)
end

## What each of the first `np` origins sends away per unit of its own
## generated infections: the importation kernel's off-diagonal column sums
## over those patches, constant in time. The kernel may cover more patches
## than the model runs, so `np` comes from the caller.
function _patch_outflow(::Type{T}, K::AbstractMatrix, np::Integer) where {T}
    outflow = zeros(T, np)
    @inbounds for q in 1:np, r in 1:np
        r == q && continue
        outflow[q] += K[r, q]
    end
    return outflow
end

"""
Convolve a daily trajectory `x` (infections or onsets) with a delay PMF
`delay` (indexed from lag 0), returning the expected daily counts of the
delayed event on the same daily grid: entry `t` sums `x[t−d] · delay[d+1]`
over lags `d` that stay in range. Maps infections to onsets, onsets to
deaths, onsets to reports and onsets to detected exports. Type-stable and
AD-transparent.

Each lag adds one scaled, shifted copy of `x` with `axpy!`, a single BLAS
call on float arrays. BLAS skips a lag whose weight is exactly zero, so on
float arrays an `Inf` or `NaN` in `x` does not reach the days that lag
feeds; other element types add `0 · x` there and propagate it.
"""
function convolve_delay(x::AbstractVector, delay::AbstractVector)
    n = length(x)
    y = zeros(promote_type(eltype(x), eltype(delay)), n)
    for d in 1:min(length(delay), n)
        axpy!(delay[d], view(x, 1:(n - d + 1)), view(y, d:n))
    end
    return y
end

"""
Discrete convolution of two delay PMFs `a` and `b` (both indexed from
delay 0 at element 1), giving the PMF of the summed delay `a ⊕ b`. The
result has length `length(a) + length(b) - 1`. Its mass equals
`sum(a) * sum(b)`, so normalised inputs give a normalised output. Used to
build the infection→detection delay (incubation ⊕ onset-to-detection) and
the infection→death delay (incubation ⊕ onset-to-death) for the exports
streams from their component PMFs. Type-stable and AD-transparent.
"""
function convolve_pmf(a::AbstractVector, b::AbstractVector)
    (isempty(a) || isempty(b)) &&
        return zeros(promote_type(eltype(a), eltype(b)), 0)
    na = length(a)
    nb = length(b)
    Tp = promote_type(eltype(a), eltype(b))
    y = zeros(Tp, na + nb - 1)
    @inbounds for i in 1:na, j in 1:nb
        y[i + j - 1] += a[i] * b[j]
    end
    return y
end

"""
Day indices of the weekly reproduction-number knots over an `n`-day grid.
The first knot sits on day `start` (default 1) and the last knot on day
`n`, with regular knots every `week` days, so a knot is pinned to `start`
and to the end of the grid. With `start > 1` the reproduction number is
held flat (at the first knot's value) for all days before `start`
([`interpolate_knots`](@ref) clamps below the first knot), so the random
walk only varies `R_t` from `start` onward. This fixes `R_t` over the
pre-establishment seeding window before the genetic TMRCA bound. Returns a
sorted vector of unique day indices.
"""
function knot_days(n::Integer; week::Integer = 7, start::Integer = 1)
    n <= 1 && return [1]
    s = clamp(Int(start), 1, n)
    days = collect(s:week:n)
    days[end] == n || push!(days, n)
    return days
end

"""
    ForecastHorizon(days)

The number of days a model runs past its cut-off to forecast. Passed as a
composer's `forecast` keyword (see [`with_horizon`](@ref)). The latent grid
then runs to day `n + days`, the walks draw fresh innovations for the
future knots, and each stream draws its future counts as `missing`
observations for `predict` to generate. The default `forecast = nothing`
is the fitted model. It is a type rather than an integer so the fitted
model compiles with no forecast code in it.
"""
struct ForecastHorizon
    days::Int
    function ForecastHorizon(days::Integer)
        days >= 1 || throw(
            ArgumentError("a forecast horizon needs at least one day, got $days")
        )
        return new(Int(days))
    end
end

"Days a model runs past its cut-off: `0` for the fitted model."
horizon_days(::Nothing) = 0
horizon_days(f::ForecastHorizon) = f.days

"""
Knot days after the cut-off `n` for a forecast of `horizon` days: one every
`week` days and one on the last day. [`knot_days`](@ref) pins a knot to the
cut-off, so these continue the fitted knots without moving any of them.
They are also the future vintages the weekly forecast quantities are
binned to. Empty when `horizon` is zero.
"""
function future_knot_days(n::Integer, horizon::Integer; week::Integer = 7)
    horizon <= 0 && return Int[]
    days = collect((n + week):week:(n + horizon))
    (isempty(days) || days[end] != n + horizon) && push!(days, n + horizon)
    return days
end

"The first `n` entries of `x`, or `x` itself when it has no more."
upto(x::AbstractVector, n::Integer) = length(x) <= n ? x : x[1:n]

"""
Smooth intervention ramp over an `n`-day grid: the logistic curve
`1 / (1 + e^{−(t − day) / ramp})` for each day `t`, rising from ≈0 well
before `day` to ≈1 well after, with `ramp` setting the transition width
in days. Multiplied by a sampled effect size and added to log-`R_t`, this
gives an intervention (e.g. the first WHO situation report) a gradual
ramped effect on transmission rather than an instantaneous step. Returns
a length-`n` `Float64` vector. `day = missing` gives an all-zero ramp (no
intervention). Type-stable and AD-transparent in the effect size it
multiplies.

Split into two dispatches on `day`'s concrete type rather than one method
with a runtime `ismissing(day)` branch. Mooncake cannot build a reverse rule
for the branching form (`TypeError: non-boolean (Missing) used in boolean
context`). Dispatch resolves `day`'s type at the call site, so each method
body is differentiated on its own concretely-typed slot.
"""
function sigmoid_ramp(n::Integer, day::Missing; ramp::Real = RT_INTERVENTION_RAMP)
    return zeros(Float64, n)
end

function sigmoid_ramp(n::Integer, day::Real; ramp::Real = RT_INTERVENTION_RAMP)
    return Float64[logistic((t - day) / ramp) for t in 1:n]
end

"""
Gaussian random walk in centred form, shifted by a known offset: entry `i`
is the walk's level `i` plus `offset[i]`. The walk starts at `start`, and
each level is `Normal` about the one before it with SD `σ`. The shifted
values themselves are the sampled coordinates. `σ` is floored at `eps` so a
zero draw still gives a proper density. [`rt_walk_model`](@ref) draws its
weekly log-`R_t` knots from it, with the intervention ramp as the offset.
"""
struct RandomWalkVector{
        S <: Real, T <: Real, V <: AbstractVector{<:Real},
    } <: Distributions.ContinuousMultivariateDistribution
    "Level the first step is taken from."
    start::S
    "Per-step SD."
    σ::T
    "Shift added to each level."
    offset::V
end

Base.length(d::RandomWalkVector) = length(d.offset)
Base.eltype(::Type{<:RandomWalkVector}) = Float64

function Distributions._logpdf(d::RandomWalkVector, x::AbstractVector)
    T = float(
        promote_type(
            typeof(d.start), typeof(d.σ), eltype(d.offset), eltype(x)
        )
    )
    σ = d.σ + eps(typeof(float(d.σ)))
    c = -log(σ) - T(log(2π) / 2)
    s = zero(T)
    prev = d.start
    @inbounds for i in eachindex(x, d.offset)
        level = x[i] - d.offset[i]
        s += c - ((level - prev) / σ)^2 / 2
        prev = level
    end
    return s
end

function Distributions._rand!(
        rng::AbstractRNG, d::RandomWalkVector, x::AbstractVector{<:Real}
    )
    σ = d.σ + eps(typeof(float(d.σ)))
    prev = d.start
    @inbounds for i in eachindex(x, d.offset)
        prev += σ * randn(rng)
        x[i] = prev + d.offset[i]
    end
    return x
end

## The values are unconstrained, so a sampled walk links through the
## identity.
VectorBijectors.from_linked_vec(::RandomWalkVector) =
    VectorBijectors.TypedIdentity()
VectorBijectors.to_linked_vec(::RandomWalkVector) =
    VectorBijectors.TypedIdentity()
VectorBijectors.linked_vec_length(d::RandomWalkVector) = length(d)
VectorBijectors.linked_optic_vec(d::RandomWalkVector) =
    VectorBijectors.optic_vec(d)

"""
Outbreak age in days: the elapsed time from the model-implied seeding
day to the cut-off (day `n`), where the seeding day is the smooth
crossing at which cumulative infections first reach one. The crossing is
linearly interpolated between the two days that bracket a cumulative of
one, so it is a continuous function of the trajectory. Before the
trajectory reaches one it returns `n` (the full grid). Used for the
seeding-date plots and the genetic-TMRCA bound.
"""
function seeding_age(cumulative::AbstractVector, n::Integer)
    Tp = float(eltype(cumulative))
    one_ = one(Tp)
    cumulative[end] < one_ && return Tp(n)
    j = 1
    @inbounds while j < length(cumulative) && cumulative[j] < one_
        j += 1
    end
    ## j is the first day at or above one. Interpolate within [j-1, j].
    if j == 1
        cross = one(Tp)
    else
        lo = cumulative[j - 1]
        hi = cumulative[j]
        frac = hi == lo ? zero(Tp) : (one_ - lo) / (hi - lo)
        cross = (j - 1) + frac
    end
    return Tp(n) - cross
end

"""
Linearly interpolate the knot values `knot_vals`, placed on the day
indices `days`, onto the full daily grid `1:n`, returning the length-`n`
series. Piecewise-linear between bracketing knots, so the series bends
only at the knots and is otherwise straight. Applied on the log-`R_t`
scale, this gives weekly random-walk knots with within-week linear
interpolation. Type-stable and AD-transparent (the output element type
follows `knot_vals`).

Outside the knot span the series is held flat at the nearest knot value
rather than extrapolated (the interpolation fraction is clamped to
`[0, 1]`). This lets the reproduction-number walk start at a day `> 1` and
hold `R_t` flat at the established `R0` over every earlier day, rather than
running the first segment's slope backwards off the start of the grid.
"""
function interpolate_knots(
        knot_vals::AbstractVector,
        days::AbstractVector{<:Integer}, n::Integer
    )
    Tp = eltype(knot_vals)
    out = Vector{Tp}(undef, n)
    nb = length(days)
    ## A single knot spans no segment to interpolate over (and `days[b + 1]`
    ## would read out of bounds), so the window holds flat at that knot's value.
    if nb == 1
        fill!(out, knot_vals[1])
        return out
    end
    @inbounds for t in 1:n
        b, frac = _knot_bracket(days, t, Tp)
        out[t] = knot_vals[b] + frac * (knot_vals[b + 1] - knot_vals[b])
    end
    return out
end

## The knot segment `b` that day `t` falls in, between `days[b]` and
## `days[b + 1]`, and the fraction of the way along it. The fraction is
## clamped to `[0, 1]` so days outside the knot span hold flat at the nearest
## knot instead of extrapolating the end segment. Needs at least two knots.
@inline function _knot_bracket(
        days::AbstractVector{<:Integer}, t::Integer, ::Type{T}
    ) where {T}
    nb = length(days)
    b = 1
    @inbounds while b < nb - 1 && t > days[b + 1]
        b += 1
    end
    @inbounds d0, d1 = days[b], days[b + 1]
    frac = d1 == d0 ? zero(T) :
        clamp(T(t - d0) / T(d1 - d0), zero(T), one(T))
    return b, frac
end

"""
Partially pooled relative multiplier over the units of each group, with
log contrasts that sum to zero within the group,

```math
m_{u \\in g} = \\exp\\bigl(\\sigma\\, (Q_g z_g)_u\\bigr),
\\qquad z_g \\sim N(0, I_{n_g - 1}).
```

`Q_g` is the sum-to-zero basis of group `g` ([`sum_to_zero_basis`](@ref)),
so the log multipliers have the distribution of `n_g` independent
`N(0, σ²)` draws centred within the group, from the `n_g - 1` directions a
composition over the group can see. The level stays with the group and as
`σ` shrinks the composition weights its units by the modelled quantity
alone. This is the construction of [`province_composition_model`](@ref),
applied with one group per patch over its health zones in
[`bvd_zone`](@ref).

`z` holds each group's draws in turn, [`relative_multiplier_dims`](@ref) in
all. A group of one unit takes no draw and its multiplier is one.

With `Q` the block-diagonal basis of [`relative_multiplier_basis`](@ref)
the log multipliers are the one product `σ Q z`. Pass `Q` in place of
`groups` to build it once rather than on every call.
"""
function relative_multiplier(
        z::AbstractVector, σ::Real,
        groups::AbstractVector{<:UnitRange}
    )
    return relative_multiplier(z, σ, relative_multiplier_basis(groups))
end

function relative_multiplier(z::AbstractVector, σ::Real, Q::AbstractMatrix)
    length(z) == size(Q, 2) || throw(
        DimensionMismatch(
            "relative_multiplier: $(length(z)) draws for " *
                "$(size(Q, 2)) contrasts"
        )
    )
    return exp.(σ .* (Q * z))
end

"""
Block-diagonal sum-to-zero basis of [`relative_multiplier`](@ref),
`(n_units × relative_multiplier_dims(groups))`.

Group `g` holds its basis `Q_g` ([`sum_to_zero_basis`](@ref)) on its own
rows and its own [`relative_multiplier_dims`](@ref) columns, in the order of
`groups`. Every other entry is zero, so a unit in a group of one, or in no
group, has a zero row. `n_units` is the last unit of any group.
"""
function relative_multiplier_basis(groups::AbstractVector{<:UnitRange})
    nu = maximum((last(us) for us in groups if !isempty(us)); init = 0)
    Q = zeros(nu, relative_multiplier_dims(groups))
    off = 0
    for us in groups
        k = length(us) - 1
        k >= 1 || continue
        Q[us, (off + 1):(off + k)] = sum_to_zero_basis(k + 1)
        off += k
    end
    return Q
end

"Number of draws [`relative_multiplier`](@ref) takes over `groups`."
relative_multiplier_dims(groups::AbstractVector{<:UnitRange}) =
    sum((max(length(us) - 1, 0) for us in groups); init = 0)

"""
Mean-reverting AR(1) deviation knots over the health zones of
[`bvd_zone`](@ref), the construction of the province deviations in
[`patch_rt_model`](@ref) applied within each group. `groups` are the index
ranges the deviations sum to zero within, one per patch.

Group `g` of `n_g` units draws its level on the `n_g - 1` directions of its
sum-to-zero basis `Q_g` ([`sum_to_zero_basis`](@ref)), and its innovations
on the `n^W_g - 1` directions of the basis `Q^W_g` over its walking units
`W_g`, zero on the other units,

```math
δ_g(1) = σ_L Q_g A_g z^L_g,
\\qquad
δ_g(k) = φ\\, δ_g(k-1) + σ_{δ,g} Q^W_g A^W_g z^δ_{g,k},
```

through [`sum_to_zero_knots`](@ref), so every group sums to zero at every
knot. `level_factors[g]` is `A_g`, the lower
Cholesky factor of `Q_gᵀ C_g Q_g` for the units' correlation `C_g`, and
`drift_factors[g]` is `A^W_g`, the same over the walking units
([`zone_correlation_factors`](@ref)). An empty list, or a group past its
end, takes `A = I`. As `Q_g Q_gᵀ = I - J / n_g`, the level's covariance is

```math
σ_L^2\\, Q_g A_g A_g^\\top Q_g^\\top = σ_L^2\\, P_g C_g P_g,
\\qquad P_g = I - J / n_g,
```

the covariance of `σ_L L z` centred within the group for `L Lᵀ = C_g` and
`z ~ N(0, I_{n_g})`, from `n_g - 1` draws instead of `n_g`. The innovations
are the same over the walking units, and a level-only unit decays along
the mean path `φ^{k-1} δ_u(1)`.

`σ_δ` holds one drift scale per group. `z_level` holds each group's level
draws in turn and `z_drift` the innovations knot by knot, each knot's
block holding each group's in turn, [`deviation_knot_dims`](@ref) of each.
Returns an `(n_units × n_knots)` matrix.
"""
function deviation_knots(
        z_level::AbstractVector, z_drift::AbstractVector,
        σ_level::Real, σ_δ::AbstractVector, φ::Real,
        groups::AbstractVector{<:UnitRange},
        level_factors::AbstractVector, drift_factors::AbstractVector,
        walking::AbstractVector{Bool}, n_knots::Integer
    )
    dims = deviation_knot_dims(groups, walking)
    length(z_level) == dims.level || throw(
        DimensionMismatch(
            "deviation_knots: $(length(z_level)) level draws for " *
                "$(dims.level) directions"
        )
    )
    length(z_drift) == dims.drift * (n_knots - 1) || throw(
        DimensionMismatch(
            "deviation_knots: $(length(z_drift)) innovation draws for " *
                "$(dims.drift) directions over $(n_knots - 1) knots"
        )
    )
    Z = reshape(z_drift, dims.drift, n_knots - 1)
    Tp = promote_type(
        eltype(z_level), eltype(z_drift), typeof(float(σ_level)),
        eltype(σ_δ), typeof(float(φ)),
        _factor_eltype(level_factors), _factor_eltype(drift_factors)
    )
    δ = zeros(Tp, length(walking), n_knots)
    off_level = 0
    off_drift = 0
    for (g, us) in enumerate(groups)
        n = length(us)
        n >= 2 || continue
        pos = [a for (a, u) in enumerate(us) if walking[u]]
        k = max(length(pos) - 1, 0)
        F_level = sum_to_zero_factor(
            sum_to_zero_basis(n), σ_level, _group_factor(level_factors, g)
        )
        F_drift = sum_to_zero_factor(
            _embedded_basis(n, pos), σ_δ[g], _group_factor(drift_factors, g)
        )
        knots = sum_to_zero_knots(
            F_level, F_drift, z_level[(off_level + 1):(off_level + n - 1)],
            Z[(off_drift + 1):(off_drift + k), :], φ
        )
        @inbounds for j in 1:n_knots, (a, u) in enumerate(us)
            δ[u, j] = knots[a, j]
        end
        off_level += n - 1
        off_drift += k
    end
    return δ
end

"""
Draws [`deviation_knots`](@ref) takes over `groups`: `level`, the
`n_g - 1` level directions summed over the groups, and `drift`, the
`n^W_g - 1` innovation directions per knot summed over the groups with a
walking unit.
"""
function deviation_knot_dims(
        groups::AbstractVector{<:UnitRange}, walking::AbstractVector{Bool}
    )
    drift = sum(
        (max(count(u -> walking[u], us) - 1, 0) for us in groups); init = 0
    )
    return (; level = relative_multiplier_dims(groups), drift)
end

## The sum-to-zero basis over the units at positions `pos` of a group of
## `n`, zero on its other rows.
function _embedded_basis(n::Integer, pos::AbstractVector{<:Integer})
    Q = zeros(n, max(length(pos) - 1, 0))
    length(pos) >= 2 || return Q
    Q[pos, :] = sum_to_zero_basis(length(pos))
    return Q
end

_group_factor(factors::AbstractVector, g::Integer) =
    g <= length(factors) ? factors[g] : nothing

## Element type of a list of correlation factors, ignoring the `nothing`
## entries that stand for independent draws.
function _factor_eltype(factors::AbstractVector)
    T = Float64
    for F in factors
        F === nothing && continue
        T = promote_type(T, eltype(F))
    end
    return T
end

"""
Derive the implied national reproduction number from a summed infection
trajectory by inverting the renewal equation:

    Rt_national(t) = I_total(t) / sum_s I_total(t-s) * g_s

This reconstructs what a single-patch model would estimate as the national
Rt from the aggregated infection count. The first day is set to zero (no
prior infections to divide by). Days where the force of infection is zero
(no prior infections) also return zero. AD-transparent under Mooncake
(only arithmetic and `@inbounds` loops).

A model that needs the reproduction number on one day only should call
[`implied_national_Rt_at`](@ref), which this is the trajectory form of.
"""
function implied_national_Rt(
        infections_total::AbstractVector,
        g::AbstractVector
    )
    n = length(infections_total)
    Tp = promote_type(eltype(infections_total), eltype(g))
    Rt = zeros(Tp, n)
    for t in 2:n
        Rt[t] = implied_national_Rt_at(infections_total, g, t)
    end
    return Rt
end

"""
    implied_national_Rt_at(infections_total, g, t)

The implied national reproduction number on a single day `t`, the one entry
[`implied_national_Rt`](@ref) would put at index `t`. Day one and any day whose
force of infection is zero give zero, as they do there.

The model reports the aggregate reproduction number at the cut-off alone, and
building the whole trajectory to read its last entry would put `n` divisions
and `n` force sums on the gradient tape for one number.
"""
function implied_national_Rt_at(
        infections_total::AbstractVector,
        g::AbstractVector, t::Integer
    )
    Tp = promote_type(eltype(infections_total), eltype(g))
    t <= 1 && return zero(Tp)
    force = zero(Tp)
    kmax = min(t - 1, length(g))
    @inbounds for s in 1:kmax
        force += infections_total[t - s] * g[s]
    end
    force > zero(Tp) || return zero(Tp)
    return @inbounds(safe_rate(infections_total[t])) / force
end
