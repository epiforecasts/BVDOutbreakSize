# Discrete-time renewal primitives. Pure, allocation-light, and
# AD-transparent (output element types are promoted from the inputs) so
# they differentiate cleanly under Mooncake inside a Turing model. These
# back the renewal architecture: a generating infection process whose
# expected infections every downstream stream consumes, with delays
# applied by daily convolution rather than continuous-time integrals.

"""
NaN / Inf-safe positive rate. The renewal recursion can transiently
overflow on extreme NUTS warmup proposals (large `R_t` compounding),
giving a non-finite expected count. A plain `max(x, eps)` would
propagate the NaN (`max(NaN, eps) = NaN`) and trip the Poisson /
NegativeBinomial domain check.
"""
@inline function safe_rate(x)
    return isfinite(x) ? max(x, eps(typeof(x))) : eps(typeof(x))
end

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
Daily probability mass function for the continuous delay `dist` over
lags `0, 1, …, nmax`, discretised by double interval censoring (uniform
primary event over a one-day window, then unit-interval censoring of the
secondary event) via
`CensoredDistributions.double_interval_censored`. This is the
discrete analogue of the continuous onset-to-event densities used by the
integral model, and is the discretisation route the renewal convolutions
rely on. For a LogNormal primary the CDF differentiates cleanly under
Mooncake, so this is AD-safe. Extreme warmup proposals that drive the
quadrature to a non-finite or zero total fall back to a uniform PMF, so
the downstream convolution stays finite (the proposal is still rejected
through its low log-likelihood). Returns a vector whose element type
follows the delay parameters.
"""
function discretise_censored(dist, nmax::Integer)
    dic = double_interval_censored(dist; interval = 1.0, upper = float(nmax))
    return _pmf_from_dic(dic, dist, nmax)
end

## Function barrier: `double_interval_censored` returns a `Union` of solver
## types, so it is inferred abstractly at the call site above. Isolating
## the PMF loop in its own method lets it specialise on the concrete `dic`
## type, making the ~`nmax` censored-CDF evaluations type-stable under AD.
@inline function _pmf_from_dic(dic, dist, nmax::Integer)
    ## Differencing one CDF path, not `nmax + 1` overlapping `pdf` calls. The
    ## interval-censored lag-`d` mass is `cdf(dic, d+1) − cdf(dic, d)`, and
    ## `pdf(dic, d)` computes exactly that pair, so the old
    ## `[pdf(dic, d) for d in 0:nmax]` evaluated every interior integer
    ## boundary CDF twice. Evaluating the boundary CDFs once over `0:nmax+1`
    ## and taking adjacent differences halves the censored-CDF evaluations the
    ## Mooncake reverse pass walks (each CDF is the expensive part: a
    ## primary-censored, truncation-normalised incomplete-gamma / Normal-CDF
    ## call), and is numerically identical to the `pdf` differences it
    ## replaces. `IntervalCensored`'s `cdf` floors to the interval, so at an
    ## integer boundary it returns the same inner CDF the `pdf` pair reads,
    ## with the same below-minimum→0 / at-maximum→1 edge handling.
    ##
    ## CensoredDistributions' batched `pdf(dic, 0:nmax)` also evaluates each
    ## boundary CDF once (via a `Dict` cache) and is value-identical, but its
    ## cache path is not Mooncake-differentiable (it hits a bitcast-to-a-
    ## differentiable-type error in the reverse pass), so the gradient hot
    ## path stays on this plain-array CDF difference, which Mooncake handles
    ## cleanly.
    c = [cdf(dic, float(b)) for b in 0:(nmax + 1)]
    z0 = zero(eltype(c))
    raw = [max(c[i + 1] - c[i], z0) for i in 1:(nmax + 1)]
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
Sizing the truncation by the distribution rather than a hand-set constant
keeps a consistent tail mass (98% by default) across every delay. This is a
deterministic function of the prior-centre distribution, evaluated once
outside the Turing model when each delay submodel is constructed, so the PMF
length is fixed and AD-safe.
"""
function cdf_nmax(dist; q::Real = 0.98, cap::Integer = 120, minlag::Integer = 5)
    clamp(ceil(Int, quantile(dist, q)), minlag, cap)
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
EpiNow2 Stan model, replacing the doubling-time parameterisation of the
integral model.
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
relation `R = 1 / Σ_s g_s e^{−r s}`. The inverse of [`euler_lotka_r`](@ref):
where that solves `R · Σ_s g_s e^{−r s} = 1` for `r` given `R`, this returns
`R` directly for a given `r`, so the prior can be placed on the growth rate
and the reproduction number derived from it under the model's generation
interval. Uses only arithmetic and `exp`, so it is AD-transparent under
Mooncake.
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
Doubling time `log(2) / r` implied by an exponential growth rate `r`,
the renewal model's reported analogue of the integral model's sampled
doubling time. Returns a non-finite value as `r` crosses zero, matching
the limit of an unbounded doubling time at zero growth.
"""
doubling_time(r) = log(oftype(float(r), 2)) / r

"""
Seed the first `len` days of the infection trajectory as exponential
growth `I_t = I0 · e^{r (t − len)}` at the implied growth rate `r` (see
[`euler_lotka_r`](@ref)), so the seeding window is pinned at `I0` on
its last day and tails off backwards. This is the model-based
initialisation used by EpiNow2 and EpiAware.jl in place of placing the
whole seed on a single day, which would otherwise inject a transient the
renewal recursion has to relax away from. Returns a length-`len` vector
whose element type follows `I0` and `r`.
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
Renewal-start seed magnitude for the two-phase renewal: the cumulative
infection count reached by the analytic cryptic phase at the renewal-start
day. The renewal is two-phase: an analytic exponential cryptic phase from
the origin to the renewal start (≈ the genetic TMRCA day, off the renewal
grid), then the renewal recursion on `[renewal_start, cut-off]`. The
doubling count `m` counts the doublings during the cryptic phase
(origin → renewal start), so the cryptic phase grows a single import to

```math
\\text{seed\\_at\\_renewal\\_start} = 2^m,
```

at the renewal start, independent of the growth rate `r`. The rate `r`
shapes the cryptic exponential history feeding the recursion just before
the renewal start (see [`seed_infections`](@ref)), but the magnitude at the
renewal start is fixed by `m` alone. This deliberately keeps `r` (hence
`R0`) out of the seed magnitude. A cut-off-referenced size
`2^m e^{-rτ_obs}` would put `r` into both the seed and the renewal growth,
so the two would cancel for a fixed realised size and leave a flat ridge
along which `R0` could slide freely. With the magnitude `r`-free, the
renewal grows `2^m` forward over the observation window under the
time-varying `R_t`, so the realised cut-off size is data-driven through
`R_t` while the prior fixes only the renewal-start scale. The argument is
`C_T_prior = 2^m`, returned unchanged and kept as a named helper for the
seeding call site.
"""
@inline function seed_at_renewal_start(C_T_prior)
    return C_T_prior
end

"""
Daily latent infections from the renewal equation
`I_t = R_t Σ_{s ≥ 1} I_{t−s} g_s`, with generation-interval PMF `g`
(indexed from lag 1), per-day reproduction numbers `Rt` (length `n`) and
a pre-computed `seed` of length `L < n` filling the first `L` days (see
[`seed_infections`](@ref)). The recursion runs for days `L+1 … n`, so
`Rt[1]` (used to imply the seeding growth) and the seed are mutually
consistent. Returns the length-`n` infection trajectory. The output
element type is promoted from `Rt`, `g` and `seed`.

!!! note "Multi-patch analogue"
    See [`patch_infections`](@ref) for the meta-population extension
    with between-patch importation.
"""
function renewal_infections(Rt::AbstractVector, g::AbstractVector,
        seed::AbstractVector)
    n = length(Rt)
    L = length(seed)
    Tp = promote_type(eltype(Rt), eltype(g), eltype(seed))
    I = zeros(Tp, n)
    @inbounds for j in 1:min(L, n)
        I[j] = seed[j]
    end
    @inbounds for t in (L + 1):n
        force = zero(Tp)
        kmax = min(t - 1, length(g))
        for s in 1:kmax
            force += I[t - s] * g[s]
        end
        I[t] = Rt[t] * force
    end
    return I
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
current day. Pure and AD-transparent: only arithmetic and `@inbounds` loops,
no allocations of tracked containers.
"""
function importation_from_kernel(K::AbstractMatrix, I_prev::AbstractVector,
        epsilon::Real)
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

"""
    patch_infections(Rt_matrix, g, seeds_matrix, importation_kernel, epsilon)

Multi-patch (meta-population) renewal with between-patch importation.
Each patch `p` follows a modified renewal equation on a shared daily grid:

```math
I_{p,t} = R_{p,t}\\, \\sum_{s \\ge 1} I_{p,t-s}\\, g_s\\;+\\;\\text{importation}_{p,t}
```

where the importation term couples patches through a kernel `K`:

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
- `epsilon`: importation intensity, scaling the whole kernel. Coupling
  conserves infections: the origin patch is debited exactly what the
  destination patches are credited, so the national total is the same as
  it would be with `epsilon = 0`. Importation changes where infections
  occur, never how many.

# Returns

Matrix `I` of shape `(n_patches, n_days)` where row `p` is the daily
infection trajectory for patch `p`. The first `L` days are copied from
`seeds_matrix`; the remaining days are the renewal recursion with
importation. The element type is promoted from all input types.

# AD transparency

Uses only basic arithmetic and `@inbounds` loops. No `push!`, `append!`,
closures that capture mutated variables, or other constructs that would
obscure Mooncake's AD reverse pass. The importation is computed inline
in each day's patch loop (no closure allocation).
"""
function patch_infections(Rt_matrix::AbstractMatrix, g::AbstractVector,
        seeds_matrix::AbstractMatrix, importation_kernel::AbstractMatrix,
        epsilon::Real)
    np, n = size(Rt_matrix)
    L = size(seeds_matrix, 2)
    Tp = promote_type(eltype(Rt_matrix), eltype(g), eltype(seeds_matrix),
        eltype(importation_kernel), typeof(float(epsilon)))
    I = zeros(Tp, np, n)
    @inbounds for p in 1:np
        for j in 1:min(L, n)
            I[p, j] = seeds_matrix[p, j]
        end
    end
    gen = zeros(Tp, np)
    @inbounds for t in (L + 1):n
        ## What each patch generates today from its own renewal force.
        for p in 1:np
            force = zero(Tp)
            kmax = min(t - 1, length(g))
            for s in 1:kmax
                force += I[p, t - s] * g[s]
            end
            gen[p] = Rt_matrix[p, t] * force
        end
        ## Importation redistributes that transmission rather than adding to
        ## it: a fraction `epsilon * K[p, q]` of what `q` generates is realised
        ## in `p` instead of at home, so `q` is debited exactly what the
        ## destinations are credited and the national total is untouched by
        ## coupling. Crediting the destination without debiting the origin, as
        ## an earlier version did, makes coupling a net source of infections --
        ## every patch's total rises and none falls -- and the surplus
        ## compounds through the renewal into the national cumulative total.
        for p in 1:np
            outflow = zero(Tp)
            arrivals = zero(Tp)
            for q in 1:np
                q == p && continue
                outflow += importation_kernel[q, p]
                arrivals += importation_kernel[p, q] * gen[q]
            end
            I[p, t] = (one(Tp) - epsilon * outflow) * gen[p] +
                      epsilon * arrivals
        end
    end
    return I
end

"""
As [`patch_infections`](@ref), but anchored so that the reproduction number
implied by the summed patches is exactly `national_rt` on every day.

`Rt_matrix` is read for its per-patch ratios only. On each day the patch
reproduction numbers are scaled by one common factor chosen so that the
force-weighted mean matches `national_rt[t]`, which is the quantity
[`implied_national_Rt`](@ref) recovers from the national trajectory.

This is what makes `mu(t)` in [`patch_rt_model`](@ref) mean what the model
says it means. The deviations are centred as `sum_p delta_p = 0`, so `mu` is
the geometric mean of the patch reproduction numbers, while the national
epidemic grows at the force-weighted arithmetic mean. The arithmetic mean is
the larger of the two, and the gap widens as the dominant patch pulls away, so
without anchoring the country grows faster than the trend the molecular-clock
prior constrains, and the surplus compounds over the whole renewal window into
the national cumulative total. Anchoring makes the deviations pure contrasts
between provinces and leaves the national level to `mu` alone.

Returns `(; infections, Rt_matrix)`, the second being the realised per-patch
reproduction numbers after scaling, which is what should be reported.
"""
function patch_infections_anchored(Rt_matrix::AbstractMatrix,
        g::AbstractVector, seeds_matrix::AbstractMatrix,
        importation_kernel::AbstractMatrix, epsilon::Real,
        national_rt::AbstractVector)
    np, n = size(Rt_matrix)
    L = size(seeds_matrix, 2)
    Tp = promote_type(eltype(Rt_matrix), eltype(g), eltype(seeds_matrix),
        eltype(importation_kernel), typeof(float(epsilon)),
        eltype(national_rt))
    I = zeros(Tp, np, n)
    Rt_realised = zeros(Tp, np, n)
    @inbounds for p in 1:np
        for j in 1:min(L, n)
            I[p, j] = seeds_matrix[p, j]
        end
        for t in 1:min(L, n)
            Rt_realised[p, t] = Rt_matrix[p, t]
        end
    end
    force = zeros(Tp, np)
    gen = zeros(Tp, np)
    @inbounds for t in (L + 1):n
        total_force = zero(Tp)
        weighted = zero(Tp)
        for p in 1:np
            f = zero(Tp)
            kmax = min(t - 1, length(g))
            for s in 1:kmax
                f += I[p, t - s] * g[s]
            end
            force[p] = f
            total_force += f
            weighted += Rt_matrix[p, t] * f
        end
        ## The common scale that puts the force-weighted mean on `national_rt`.
        ## With no force yet, or a degenerate weighted sum, leave the ratios
        ## alone rather than dividing by zero; those days carry no infections
        ## so the choice does not affect the trajectory.
        scale = (total_force > 0 && weighted > 0) ?
                national_rt[t] * total_force / weighted : one(Tp)
        for p in 1:np
            Rt_realised[p, t] = scale * Rt_matrix[p, t]
            gen[p] = Rt_realised[p, t] * force[p]
        end
        for p in 1:np
            outflow = zero(Tp)
            arrivals = zero(Tp)
            for q in 1:np
                q == p && continue
                outflow += importation_kernel[q, p]
                arrivals += importation_kernel[p, q] * gen[q]
            end
            I[p, t] = (one(Tp) - epsilon * outflow) * gen[p] +
                      epsilon * arrivals
        end
    end
    return (; infections = I, Rt_matrix = Rt_realised)
end

"""
Convolve a daily trajectory `x` (infections or onsets) with a delay PMF
`delay` (indexed from lag 0), returning the expected daily counts of the
delayed event on the same daily grid: entry `t` sums `x[t−d] · delay[d+1]`
over lags `d` that stay in range. This is the discrete convolution that
replaces the continuous onset-to-event integrals of the integral model,
and maps infections to onsets, onsets to deaths, onsets to reports and
onsets to detected exports. Type-stable and AD-transparent.

The scalar double loop is deliberate: a vectorised lag-AXPY rewrite
(`scripts/bench_convolve.jl`) was no faster under Mooncake (≈0.95x, since
the reverse pass over the scalar loop is already efficient), so the simpler
form stays. The per-gradient cost the joint actually pays sits in the delay
discretisation (the censored-distribution CDF evaluations), which the delay
`nmax` trims target instead.
"""
function convolve_delay(x::AbstractVector, delay::AbstractVector)
    n = length(x)
    Tp = promote_type(eltype(x), eltype(delay))
    y = zeros(Tp, n)
    @inbounds for t in 1:n
        acc = zero(Tp)
        dmax = min(t - 1, length(delay) - 1)
        for d in 0:dmax
            acc += x[t - d] * delay[d + 1]
        end
        y[t] = acc
    end
    return y
end

"""
Survival-weighted convolution for an occupancy (prevalence) stream. Given a
daily admission series `x` and a length-of-stay PMF `los` (indexed from lag
0, so `los[1] = P(LOS = 0)`), return the daily occupancy

```math
\\text{occupancy}(t) = \\sum_{\\tau \\ge 0} x_{t-\\tau}\\, S(\\tau),
\\qquad S(\\tau) = P(\\text{LOS} \\ge \\tau),
```

so an admission on day `s` occupies a bed on days `s, s+1, …` until it is
discharged: the admission day is always counted (`S(0) = 1`) and a stay of
`LOS` days contributes to `LOS + 1` daily occupancies. This is the
prevalence analogue of [`convolve_delay`](@ref)'s incidence convolution:
the length-of-stay survival replaces the onset-to-event delay PMF, turning
an inflow into a stock. The survival weights are the reverse-cumulative
tail sums of the PMF (`S(τ) = Σ_{j ≥ τ} los[j]`), so for a normalised PMF
`S(0) = 1`, and the occupancy is `convolve_delay(x, S)`. Type-stable and
AD-transparent. The element type follows the inputs.
"""
function convolve_survival(x::AbstractVector, los::AbstractVector)
    L = length(los)
    surv = similar(los)
    acc = zero(eltype(los))
    @inbounds for i in L:-1:1
        acc += los[i]
        surv[i] = acc
    end
    return convolve_delay(x, surv)
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
Smooth intervention ramp over an `n`-day grid: the logistic curve
`1 / (1 + e^{−(t − day) / ramp})` for each day `t`, rising from ≈0 well
before `day` to ≈1 well after, with `ramp` setting the transition width
in days. Multiplied by a sampled effect size and added to log-`R_t`, this
gives an intervention (e.g. the first WHO situation report) a gradual
ramped effect on transmission rather than an instantaneous step. Returns
a length-`n` `Float64` vector. `day = missing` gives an all-zero ramp (no
intervention). Type-stable and AD-transparent in the effect size it
multiplies.
"""
function sigmoid_ramp(n::Integer, day::Union{Missing, Real};
        ramp::Real = RT_INTERVENTION_RAMP)
    ismissing(day) && return zeros(Float64, n)
    return Float64[logistic((t - day) / ramp) for t in 1:n]
end

"""
Outbreak age in days: the elapsed time from the model-implied seeding
day to the cut-off (day `n`), where the seeding day is the smooth
crossing at which cumulative infections first reach one. The crossing is
linearly interpolated between the two days that bracket a cumulative of
one, so it is a continuous function of the trajectory. Before the
trajectory reaches one it returns `n` (the full grid). The renewal
analogue of the integral model's sampled outbreak age `T`, used for the
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
(the interpolation fraction is clamped to `[0, 1]`), not extrapolated: days
before the first knot take `knot_vals[1]` and days after the last take
`knot_vals[end]`. This is what lets the reproduction-number walk start at a
day `> 1` and hold `R_t` flat at the established `R0` (the first knot value)
over every earlier day, rather than running the first segment's slope
backwards off the start of the grid.
"""
function interpolate_knots(knot_vals::AbstractVector,
        days::AbstractVector{<:Integer}, n::Integer)
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
        b = 1
        while b < nb - 1 && t > days[b + 1]
            b += 1
        end
        d0 = days[b]
        d1 = days[b + 1]
        ## Clamp the fraction to `[0, 1]` so days outside the knot span hold
        ## flat at the nearest knot instead of extrapolating the end segment.
        frac = d1 == d0 ? zero(Tp) :
               clamp(Tp(t - d0) / Tp(d1 - d0), zero(Tp), one(Tp))
        out[t] = knot_vals[b] + frac * (knot_vals[b + 1] - knot_vals[b])
    end
    return out
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
"""
function implied_national_Rt(infections_total::AbstractVector,
        g::AbstractVector)
    n = length(infections_total)
    Tp = promote_type(eltype(infections_total), eltype(g))
    Rt = zeros(Tp, n)
    @inbounds for t in 2:n
        force = zero(Tp)
        kmax = min(t - 1, length(g))
        for s in 1:kmax
            force += infections_total[t - s] * g[s]
        end
        if force > zero(Tp)
            Rt[t] = safe_rate(infections_total[t]) / force
        end
    end
    return Rt
end
