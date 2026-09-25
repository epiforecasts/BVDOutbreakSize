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
generation-interval PMF `g` (indexed from lag 1), solving the
Euler–Lotka identity `R · Σ_s g_s e^{−r s} = 1`. Starts from the
small-`r` approximation `r ≈ (R − 1) / (R · ḡ)` with `ḡ` the mean
generation time, then refines with `steps` Newton iterations. The loop
is unrolled over a fixed step count and uses only arithmetic and `exp`,
so it is AD-transparent under Mooncake. Mirrors the `R_to_r` seeding
helper in EpiAware.jl and the implied-growth initialisation in the
EpiNow2 Stan model.
"""
function euler_lotka_r(R, g::AbstractVector; steps::Integer = 2)
    Tp = promote_type(typeof(float(R)), eltype(g))
    ḡ = zero(Tp)
    @inbounds for i in eachindex(g)
        ḡ += g[i] * i
    end
    r = (R - one(R)) / (R * ḡ)
    @inbounds for _ in 1:steps
        G = zero(Tp)
        dG = zero(Tp)
        for i in eachindex(g)
            e = exp(-r * i)
            G += g[i] * e
            dG += g[i] * i * e
        end
        ## f(r) = R·G − 1, f'(r) = −R·dG.
        r = r - (R * G - one(R)) / (-R * dG)
    end
    return r
end

"""
Reproduction number `R` implied by an exponential growth rate `r` and a
generation-interval PMF `g` (indexed from lag 1), the forward Euler–Lotka
relation `R = 1 / Σ_s g_s e^{−r s}`. The inverse of [`euler_lotka_r`](@ref),
so a prior can be placed on the growth rate and the reproduction number
derived from it under the model's generation interval. Uses only arithmetic
and `exp`, so it is AD-transparent under Mooncake.
"""
function r_to_R0(r, g::AbstractVector)
    Tp = promote_type(typeof(float(r)), eltype(g))
    G = zero(Tp)
    @inbounds for i in eachindex(g)
        G += g[i] * exp(-r * i)
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
depletion,

```math
I_t = S_{t-1} \\left(1 - e^{-R_t F_t / N}\\right), \\qquad
F_t = \\sum_{s \\ge 1} I_{t-s} g_s, \\qquad S_t = S_{t-1} - I_t,
```

with generation-interval PMF `g` (indexed from lag 1), per-day reproduction
numbers `Rt` (length `n`), a pre-computed `seed` of length `L < n` filling the
first `L` days (see [`seed_infections`](@ref)) and a susceptible
`population` `N`. The susceptibles start at `N` less the seed, floored at
zero. While `S_{t-1} ≈ N` and `R_t F_t ≪ N` this is the plain renewal
`I_t = R_t F_t` scaled by the share still susceptible, and it keeps the
cumulative infections below `N` whatever `Rt` does. The recursion runs for
days `L+1 … n`, so `Rt[1]` (used to imply the seeding growth) and the seed
are mutually consistent. Returns the length-`n` infection trajectory. The
output element type is promoted from `Rt`, `g`, `seed` and `population`.

!!! note "Multi-patch analogue"
    See [`patch_infections`](@ref) for the meta-population extension
    with between-patch importation.
"""
function renewal_infections(
        Rt::AbstractVector, g::AbstractVector,
        seed::AbstractVector, population::Real
    )
    return first(renewal_infections_with_force(Rt, g, seed, population))
end

"""
The renewal trajectory, the per-day force of infection it was built from
and the susceptibles left after each day, as `(infections, force,
susceptible)`. [`renewal_infections`](@ref) returns the first; the
derivative rule needs the other two, which it would otherwise have to
rebuild from a copy of this loop. `susceptible` is zero on the seed days
before the last.

Each day's force is one `dot` of the most recent infections with the
generation interval reversed, a single BLAS call on float arrays.
"""
function renewal_infections_with_force(
        Rt::AbstractVector, g::AbstractVector,
        seed::AbstractVector, population::Real
    )
    n = length(Rt)
    L = length(seed)
    G = length(g)
    Tp = promote_type(
        eltype(Rt), eltype(g), eltype(seed), typeof(float(population))
    )
    I = zeros(Tp, n)
    force = zeros(Tp, n)
    susceptible = zeros(Tp, n)
    @inbounds for j in 1:min(L, n)
        I[j] = seed[j]
    end
    S = _initial_susceptible(Tp, population, seed)
    L >= 1 && L <= n && (susceptible[L] = S)
    rg = reverse(g)
    for t in (L + 1):n
        k = min(t - 1, G)
        f = dot(view(rg, (G - k + 1):G), view(I, (t - k):(t - 1)))
        force[t] = f
        I[t] = S * -expm1(-Rt[t] * f / population)
        S -= I[t]
        susceptible[t] = S
    end
    return I, force, susceptible
end

## Susceptibles on the last seed day: the population less the seed, floored
## at zero so a seed larger than the population leaves nothing to infect.
## Shared by the renewal kernels and their rules.
@inline function _initial_susceptible(::Type{T}, population, seed) where {T}
    return max(T(population) - T(sum(seed)), zero(T))
end

## --- Multi-patch (meta-population) renewal primitives --------------------

"""
    importation_from_kernel(K, I_prev, epsilon)

Per-patch importation into each of `n_patches` patches on a single day,
given the `n_patches x n_patches` importation kernel `K`, the previous
day's infections per patch `I_prev` (length `n_patches`), and the
importation intensity `epsilon`.

```math
\\text{importation}_p = \\varepsilon \\sum_{q} K_{p,q} I_{q,t-1}
```

`K[p, q]` is the per-capita daily travel rate from patch `q` to patch `p`
(the first index is the destination). Diagonal entries should be zero (no
self-importation). Each entry is unitless (a rate per day per traveller in
the source patch).

Returns a length-`n_patches` vector of imported infections expected on the
current day. AD-transparent under Mooncake.
"""
function importation_from_kernel(
        K::AbstractMatrix, I_prev::AbstractVector,
        epsilon::Real
    )
    np = size(K, 1)
    Tp = promote_type(eltype(K), eltype(I_prev), typeof(float(epsilon)))
    imp = zeros(Tp, np)
    @inbounds for p in 1:np
        acc = zero(Tp)
        for q in 1:np
            acc += K[p, q] * I_prev[q]
        end
        imp[p] = epsilon * acc
    end
    return imp
end

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
    patch_infections(
        Rt_matrix, g, seeds_matrix, importation_kernel, epsilon, populations
    )

Multi-patch (meta-population) renewal with between-patch importation and
weak susceptible depletion. Each patch `p` follows a modified renewal
equation on a shared daily grid:

```math
\\lambda_{p,t} = R_{p,t}\\, \\sum_{s \\ge 1} I_{p,t-s}\\, g_s\\;+\\;\\text{importation}_{p,t},
\\qquad
I_{p,t} = S_{p,t-1} \\left(1 - e^{-\\lambda_{p,t} / N_p}\\right)
```

as in [`renewal_infections`](@ref). The susceptibles `S_{p,t}` start at the
population `N_p` less the patch's seed and fall by each day's infections, and
the importation term couples patches through a kernel `K`:

```math
\\text{importation}_{p,t} =
    \\varepsilon \\sum_{q} K_{p,q}\\, I_{q,t-1}.
```

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
- `epsilon`: importation intensity, scaling the whole kernel. Importation
  is a transfer on the day it happens, so the origin patch is debited exactly
  what the destination patches are credited and coupling is never a source of
  infections. It is not conserved across days. The destination grows at its
  own reproduction number, so relocating infections from a fast patch to a
  slow one lowers the national total and the reverse raises it. With one
  shared reproduction number the transfer cancels exactly.
- `populations`: length-`n_patches` vector of susceptible populations `N_p`.

# Returns

`(; infections, importation)`. `infections` is a matrix of shape
`(n_patches, n_days)` where row `p` is the daily infection trajectory for
patch `p`. The first `L` days are copied from `seeds_matrix` and the
remaining days are the renewal recursion with importation. `importation` is
the matching matrix of infections each patch received from the others, the
arrivals term alone rather than the net of arrivals and departures, before
depletion at the destination. The
element type is promoted from all input types. AD-transparent under Mooncake.
"""
function patch_infections(
        Rt_matrix::AbstractMatrix, g::AbstractVector,
        seeds_matrix::AbstractMatrix, importation_kernel::AbstractMatrix,
        epsilon::Union{Real, AbstractMatrix}, populations::AbstractVector
    )
    st = patch_infections_with_state(
        Rt_matrix, g, seeds_matrix, importation_kernel, epsilon, populations
    )
    return (; st.infections, st.importation)
end

"""
[`patch_infections`](@ref) with the per-patch pressure `λ_{p,t}` and the
susceptibles left after each day as well, as `(; infections, importation,
pressure, susceptible)`. The derivative rule reads the last two rather than
rebuilding them. `susceptible` is zero on the seed days before the last.
"""
function patch_infections_with_state(
        Rt_matrix::AbstractMatrix, g::AbstractVector,
        seeds_matrix::AbstractMatrix, importation_kernel::AbstractMatrix,
        epsilon::Union{Real, AbstractMatrix}, populations::AbstractVector
    )
    np, n = size(Rt_matrix)
    L = size(seeds_matrix, 2)
    Tp = promote_type(
        eltype(Rt_matrix), eltype(g), eltype(seeds_matrix),
        eltype(importation_kernel),
        epsilon isa Real ? typeof(float(epsilon)) : eltype(epsilon),
        typeof(float(first(populations)))
    )
    I = zeros(Tp, np, n)
    imports = zeros(Tp, np, n)
    pressure = zeros(Tp, np, n)
    susceptible = zeros(Tp, np, n)
    outflow = _patch_outflow(Tp, importation_kernel, np)
    S = zeros(Tp, np)
    @inbounds for p in 1:np
        for j in 1:min(L, n)
            I[p, j] = seeds_matrix[p, j]
        end
        S[p] = _initial_susceptible(Tp, populations[p], view(seeds_matrix, p, :))
        L >= 1 && L <= n && (susceptible[p, L] = S[p])
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
            λ = (one(Tp) - _eps(epsilon, p, t) * outflow[p]) * gen[p] +
                arrivals
            pressure[p, t] = λ
            I[p, t] = S[p] * -expm1(-λ / populations[p])
            S[p] -= I[p, t]
            susceptible[p, t] = S[p]
        end
    end
    return (;
        infections = I, importation = imports, pressure, susceptible,
    )
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
