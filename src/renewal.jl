# Discrete-time renewal primitives. Pure, allocation-light, and
# AD-transparent (output element types are promoted from the inputs) so
# they differentiate cleanly under Mooncake inside a Turing model. These
# back the renewal architecture: a generating infection process whose
# expected infections every downstream stream consumes, with delays
# applied by daily convolution rather than continuous-time integrals.

"""
NaN / Inf-safe positive rate. The renewal recursion can transiently
overflow on extreme NUTS warmup proposals (large `R_t` compounding),
giving a non-finite expected count; a plain `max(x, eps)` would
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
Mooncake, so this is AD-safe; extreme warmup proposals that drive the
quadrature to a non-finite or zero total fall back to a uniform PMF so
the downstream convolution stays finite (the proposal is still rejected
through its low log-likelihood). Returns a vector whose element type
follows the delay parameters.
"""
function discretise_censored(dist, nmax::Integer)
    dic = double_interval_censored(dist; interval = 1.0, upper = float(nmax))
    raw = [pdf(dic, float(d)) for d in 0:nmax]
    s = sum(raw)
    if !isfinite(s) || s <= zero(s)
        z = zero(pdf(dist, oneunit(float(nmax))))
        return fill(one(z) / (nmax + 1), nmax + 1)
    end
    return raw ./ s
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
Doubling time `log(2) / r` implied by an exponential growth rate `r`,
the renewal model's reported analogue of the integral model's sampled
doubling time. Returns a non-finite value as `r` crosses zero, matching
the limit of an unbounded doubling time at zero growth.
"""
doubling_time(r) = log(oftype(float(r), 2)) / r

"""
Seed the first `len` days of the infection trajectory as exponential
growth `I_t = I0 · e^{r (t − len)}` at the implied growth rate `r` (see
[`euler_lotka_r`](@ref)), so the seeding window is anchored at `I0` on
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
Daily latent infections from the renewal equation
`I_t = R_t Σ_{s ≥ 1} I_{t−s} g_s`, with generation-interval PMF `g`
(indexed from lag 1), per-day reproduction numbers `Rt` (length `n`) and
a pre-computed `seed` of length `L < n` filling the first `L` days (see
[`seed_infections`](@ref)). The recursion runs for days `L+1 … n`, so
`Rt[1]` (used to imply the seeding growth) and the seed are mutually
consistent. Returns the length-`n` infection trajectory; the output
element type is promoted from `Rt`, `g` and `seed`.
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

"""
Convolve a daily trajectory `x` (infections or onsets) with a delay PMF
`delay` (indexed from lag 0), returning the expected daily counts of the
delayed event on the same daily grid: entry `t` sums `x[t−d] · delay[d+1]`
over lags `d` that stay in range. This is the discrete convolution that
replaces the continuous onset-to-event integrals of the integral model,
and maps infections to onsets, onsets to deaths, onsets to reports and
onsets to detected exports. Type-stable and AD-transparent.

The scalar double loop is deliberate: a vectorised lag-AXPY rewrite
(`scripts/bench_convolve.jl`) was no faster under Mooncake (≈0.95x — the
reverse pass over the scalar loop is already efficient), so the simpler
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
Day indices of the weekly reproduction-number knots over an `n`-day grid.
Knot 1 sits on day 1 and the last knot on day `n`, with regular knots
every `week` days, so a knot is always pinned to each end of the grid.
Returns a sorted vector of unique day indices.
"""
function knot_days(n::Integer; week::Integer = 7)
    n <= 1 && return [1]
    days = collect(1:week:n)
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
a length-`n` `Float64` vector; `day = missing` gives an all-zero ramp (no
intervention). Type-stable and AD-transparent in the effect size it
multiplies.
"""
function sigmoid_ramp(n::Integer, day::Union{Missing, Real};
        ramp::Real = 14.0)
    ismissing(day) && return zeros(Float64, n)
    return Float64[logistic((t - day) / ramp) for t in 1:n]
end

"""
Outbreak age in days: the elapsed time from the model-implied seeding
day to the cut-off (day `n`), where the seeding day is the smooth
crossing at which cumulative infections first reach one. The crossing is
linearly interpolated between the two days that bracket a cumulative of
one, so it is a continuous function of the trajectory; before the
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
    ## j is the first day at or above one; interpolate within [j-1, j].
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
only at the knots and is otherwise straight; applied on the log-`R_t`
scale this gives weekly random-walk knots with within-week linear
interpolation. Type-stable and AD-transparent (the output element type
follows `knot_vals`).
"""
function interpolate_knots(knot_vals::AbstractVector,
        days::AbstractVector{<:Integer}, n::Integer)
    Tp = eltype(knot_vals)
    out = Vector{Tp}(undef, n)
    nb = length(days)
    @inbounds for t in 1:n
        b = 1
        while b < nb - 1 && t > days[b + 1]
            b += 1
        end
        d0 = days[b]
        d1 = days[b + 1]
        frac = d1 == d0 ? zero(Tp) : Tp(t - d0) / Tp(d1 - d0)
        out[t] = knot_vals[b] + frac * (knot_vals[b + 1] - knot_vals[b])
    end
    return out
end
