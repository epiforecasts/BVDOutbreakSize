# Observation distributions, one family at a time: the scalar `safe_*`
# constructor, which guards a NegativeBinomial, BetaBinomial or Student-t
# against extreme NUTS proposals, its guards, and its vector distributions.
# A vector distribution scores a whole observed vector with a single `~`
# and draws a `missing` one with the same statement. Its `logpdf` sums the
# guarded scalar terms in one loop, and `src/mooncake_rules.jl` gives that
# `logpdf` a closed-form Mooncake rule. Its `rand` draws each entry from
# the scalar distribution. The vector is one variable, so a sampled vector
# is stored under one whole-vector key (`<prefix>.increments`,
# `<prefix>.obs`).

## `x` where it is finite and positive, otherwise `fallback`, and whether
## `x` was kept. The domain guard of `safe_nbinomial`'s dispersion and
## `safe_studentt`'s scale and degrees of freedom.
@inline function _positive_or(x, fallback)
    on = isfinite(x) && x > zero(x)
    return (on ? x : fallback), on
end

## Count draws in floating point, floored to an `Int` that saturates at
## `typemax(Int)`. The Distributions samplers convert to `Int` and throw
## `InexactError` once a mean passes it, which a forecast from a loosely
## constrained fit can reach. Following ComposableTuringIDModels'
## `SafePoisson` and `SafeNegativeBinomial`, the Gamma–Poisson mixture is
## drawn in floating point and floored safely, saturating rather than
## widening to `BigInt` because every draw lands in an `Int` container.

## A Poisson count of mean `λ` as a float. Below `2^62` this is the
## Distributions draw. Above it the sampler's own `Int` conversion can
## overflow, so the draw takes the normal limit, whose relative error there is
## below `1e-9`.
function _poisson_count(rng::AbstractRNG, λ::Real)
    λ < 2.0^62 && return float(rand(rng, Poisson(λ)))
    return max(round(λ + sqrt(λ) * randn(rng)), zero(float(λ)))
end

## A negative binomial count of shape `r` and success probability `p` as a
## float, by the Gamma–Poisson mixture `rand(::NegativeBinomial)` uses, so a
## mean within range draws the same count from the same generator state.
function _nbinomial_count(rng::AbstractRNG, d::NegativeBinomial)
    isone(d.p) && return 0.0
    return _poisson_count(rng, rand(rng, Gamma(d.r, (1 - d.p) / d.p)))
end

## `x` floored to an `Int`, saturating at `typemax(Int)`.
_saturated_count(x::Real) = x < 2.0^63 ? floor(Int, x) : typemax(Int)

"""
    SafePoisson(λ)

`Poisson` with mean `λ` whose draws saturate at `typemax(Int)` rather than
throwing `InexactError` past it. `logpdf` is the `Poisson` one. The draw
follows ComposableTuringIDModels' `SafePoisson`.
"""
struct SafePoisson{T <: Real} <: Distributions.DiscreteUnivariateDistribution
    λ::T
end

Base.minimum(::SafePoisson) = 0
Base.maximum(::SafePoisson) = Inf
Distributions.insupport(::SafePoisson, x::Real) = isinteger(x) && x >= 0
function Distributions.logpdf(d::SafePoisson, x::Real)
    return logpdf(Poisson(d.λ; check_args = false), x)
end
function Base.rand(rng::AbstractRNG, d::SafePoisson)
    return _saturated_count(_poisson_count(rng, d.λ))
end

"""
NaN / Inf-safe `NegativeBinomial` constructor parameterised by mean `μ`
and dispersion `k`, with clamping on the success probability so extreme
NUTS proposals during warmup do not trip the distribution domain check.
Shared by the count-stream observation submodels.
"""
function safe_nbinomial(k, μ)
    g = _nbinomial_params(k, μ)
    return NegativeBinomial(g.r, g.p)
end

## The guarded parameters `safe_nbinomial(k, μ)` builds from: the dispersion
## `r`, the floored mean `m` and the success probability `p`, with whether
## `p` fell inside its clamp. `r` is guarded as well as `p`, because
## `NegativeBinomial(0, p)` throws `DomainError: r > 0`, which aborts a
## gradient rather than rejecting the step.
@inline function _nbinomial_params(k, μ)
    r, _ = _positive_or(k, eps(typeof(k)))
    m = max(μ, eps(typeof(μ)))
    p_raw = r / (r + m)
    lo = eps(typeof(r))
    hi = one(r) - lo
    p_on = isfinite(p_raw) && !(p_raw > hi) && !(p_raw < lo)
    p = isfinite(p_raw) ? clamp(p_raw, lo, hi) : lo
    return (; r, m, p, p_on)
end

"""
Independent [`safe_nbinomial`](@ref) counts, entry `i` about the
[`safe_rate`](@ref) of `μ[i]` with the shared dispersion `k`. `logpdf` sums
the entries' terms in one loop.
"""
struct NegBinomialVector{T <: Real, V <: AbstractVector{<:Real}} <:
    Distributions.DiscreteMultivariateDistribution
    "Shared dispersion."
    k::T
    "Per-entry means."
    μ::V
end

Base.length(d::NegBinomialVector) = length(d.μ)
Base.eltype(::Type{<:NegBinomialVector}) = Int

function Distributions._logpdf(d::NegBinomialVector, x::AbstractVector)
    s = zero(float(promote_type(typeof(d.k), eltype(d.μ))))
    @inbounds for i in eachindex(d.μ, x)
        s += logpdf(safe_nbinomial(d.k, safe_rate(d.μ[i])), x[i])
    end
    return s
end

function Distributions._rand!(
        rng::AbstractRNG, d::NegBinomialVector, x::AbstractVector{<:Real}
    )
    @inbounds for i in eachindex(x, d.μ)
        x[i] = _saturated_count(
            _nbinomial_count(rng, safe_nbinomial(d.k, safe_rate(d.μ[i])))
        )
    end
    return x
end

"""
Log-probability that a draw from the `NegativeBinomial` `d` is at least `x`,
the censored tail of a count at its ceiling. Equal to `logccdf(d, x - 1)`,
computed as `1 − I_p(r, x)` from `SpecialFunctions.beta_inc`, whose
plain-Julia body Mooncake differentiates. `logccdf` calls Rmath's
`pnbinom` through a `ccall`, which Mooncake cannot. A tail too small for a
normal float is summed term by term in log space from `x` up, until a term
adds less than `exp(-40)` of the total. A fractional `x` takes the tail from
the next count up.
"""
function nbinomial_logtail(d::NegativeBinomial, x::Real)
    x > 0 || return zero(float(typeof(d.p)))
    u = ceil(x)
    I, J = beta_inc(d.r, u, d.p)
    I < J && return log1p(-I)
    J >= floatmin(J) && return log(J)
    s = logpdf(d, u)
    t = logpdf(d, u + 1)
    while t > s - 40
        s = logaddexp(s, t)
        u += 1
        t = logpdf(d, u + 1)
    end
    return s
end

"""
Independent right-censored [`safe_nbinomial`](@ref) counts, entry `i` about
`μ[i]` and censored at `upper[i]`, both through [`safe_rate`](@ref), with
the shared dispersion `k`. Built by
`censored(NegBinomialVector(k, μ); upper)`. A count at its ceiling scores
the censored tail [`nbinomial_logtail`](@ref).
"""
struct CensoredNegBinomialVector{
        T <: Real, V <: AbstractVector{<:Real}, U <: AbstractVector{<:Real},
    } <: Distributions.DiscreteMultivariateDistribution
    "Shared dispersion."
    k::T
    "Per-entry means."
    μ::V
    "Per-entry censoring ceilings."
    upper::U
end

"""
    censored(d::NegBinomialVector; upper::AbstractVector)

Right-censor each entry of `d` at the matching `upper`, as `censored` does
for one `NegativeBinomial`. Returns a
[`CensoredNegBinomialVector`](@ref).
"""
function Distributions.censored(d::NegBinomialVector; upper::AbstractVector)
    return CensoredNegBinomialVector(d.k, d.μ, upper)
end

Base.length(d::CensoredNegBinomialVector) = length(d.μ)
## A draw above a ceiling returns the ceiling, which need not be a whole
## number (`admission_headroom`), so censored draws are floats.
Base.eltype(::Type{<:CensoredNegBinomialVector}) = Float64

## A count below its ceiling scores the uncensored `logpdf`, so those go
## through one `NegBinomialVector` together. A count at its ceiling scores
## the censored tail `nbinomial_logtail`, and one above it `-Inf`.
function Distributions._logpdf(
        d::CensoredNegBinomialVector, x::AbstractVector
    )
    upper = safe_rate.(d.upper)
    below = x .< upper
    s = logpdf(NegBinomialVector(d.k, d.μ[below]), x[below])
    @inbounds for i in findall(!, below)
        s += x[i] > upper[i] ? oftype(s, -Inf) :
            nbinomial_logtail(safe_nbinomial(d.k, safe_rate(d.μ[i])), x[i])
    end
    return s
end

function Distributions._rand!(
        rng::AbstractRNG, d::CensoredNegBinomialVector,
        x::AbstractVector{<:Real}
    )
    @inbounds for i in eachindex(x, d.μ, d.upper)
        count = _saturated_count(
            _nbinomial_count(rng, safe_nbinomial(d.k, safe_rate(d.μ[i])))
        )
        x[i] = min(count, safe_rate(d.upper[i]))
    end
    return x
end

"""
NaN / Inf-safe overdispersed `Binomial` (`BetaBinomial`) constructor
parameterised by the trial count `n`, the mean positive probability `p`
and an intra-window overdispersion `ρ ∈ (0, 1)`. With concentration
`s = (1 − ρ)/ρ`, `α = s·p` and `β = s·(1 − p)`, the `BetaBinomial(n, α, β)`
has mean `n·p` and variance `n·p·(1 − p)·(1 + (n − 1)·ρ)`, so `ρ → 0`
recovers the plain `Binomial(n, p)` and larger `ρ` inflates the variance
above it. This carries the day-to-day laboratory batching and within-window
positivity heterogeneity a pooled per-window `p` does not. `ρ` is floored
away from `0` and `p` clamped into `(0, 1)` so the distribution stays
defined under extreme NUTS proposals. Shared by the confirmed-positives
windows.
"""
function safe_betabinomial(n::Integer, p, ρ)
    T = float(promote_type(typeof(p), typeof(ρ)))
    g = _betabinomial_shapes(p, _betabinomial_concentration(T, ρ).s)
    return BetaBinomial(n, g.α, g.β)
end

## The concentration `s = (1 − ρc) / ρc` of `safe_betabinomial`, with the
## clamped `ρc` and whether `ρ` fell inside its clamp. `ρ` is floored at
## 1e-6 (capping `s` at ≈1e6) so a near-zero draw stays a well-conditioned
## near-Binomial, and capped below 1 so `s` stays positive.
@inline function _betabinomial_concentration(::Type{T}, ρ) where {T}
    ρ_lo = T(1.0e-6)
    ρ_hi = one(T) - ρ_lo
    ρ_on = isfinite(ρ) && !(ρ < ρ_lo) && !(ρ > ρ_hi)
    ρc = isfinite(ρ) ? clamp(T(ρ), ρ_lo, ρ_hi) : ρ_lo
    return (; s = (one(T) - ρc) / ρc, ρc, ρ_on)
end

## The shapes `α = s·pc` and `β = s·(1 − pc)` of `safe_betabinomial` at
## concentration `s`, each floored at `eps`, with the mean `pc` clamped into
## `(0, 1)`. Also returns whether `p` fell inside its clamp and whether each
## shape sits above its floor.
@inline function _betabinomial_shapes(p, s::T) where {T}
    lo = eps(T)
    p_on = isfinite(p) && !(p < lo) && !(p > one(T) - lo)
    pc = isfinite(p) ? clamp(T(p), lo, one(T) - lo) : one(T) / 2
    sα = s * pc
    sβ = s * (one(T) - pc)
    return (;
        α = max(sα, lo), β = max(sβ, lo), pc, p_on,
        α_on = sα > lo, β_on = sβ > lo,
    )
end

"""
Independent [`safe_betabinomial`](@ref) counts, entry `i` out of
`trials[i]` with mean probability `p[i]` and the shared overdispersion
`ρ`. `logpdf` sums the entries' terms in one loop.
"""
struct BetaBinomialVector{
        I <: AbstractVector{<:Integer}, P <: AbstractVector{<:Real}, R <: Real,
    } <: Distributions.DiscreteMultivariateDistribution
    "Per-entry trial counts."
    trials::I
    "Per-entry mean probabilities."
    p::P
    "Shared overdispersion."
    ρ::R
end

Base.length(d::BetaBinomialVector) = length(d.p)
Base.eltype(::Type{<:BetaBinomialVector}) = Int

function Distributions._logpdf(d::BetaBinomialVector, x::AbstractVector)
    s = zero(float(promote_type(eltype(d.p), typeof(d.ρ))))
    @inbounds for i in eachindex(d.trials, d.p, x)
        s += logpdf(safe_betabinomial(d.trials[i], d.p[i], d.ρ), x[i])
    end
    return s
end

function Distributions._rand!(
        rng::AbstractRNG, d::BetaBinomialVector, x::AbstractVector{<:Real}
    )
    @inbounds for i in eachindex(x, d.trials, d.p)
        x[i] = rand(rng, safe_betabinomial(d.trials[i], d.p[i], d.ρ))
    end
    return x
end

"""
    safe_studentt(μ, σ, ν)

NaN / Inf-safe location-scale Student-t distribution `μ + σ · Tν`, built
from `Distributions.TDist(ν)` via the affine-combination operators. `σ` is
floored away from zero and non-finite values, mirroring
[`safe_nbinomial`](@ref)'s domain guard. A non-positive or non-finite `ν`
falls back to `4`, the caller's own default, rather than to the smallest
value `TDist` accepts: `TDist(1)` is Cauchy, so a floor at the domain edge
would turn a bad degrees-of-freedom argument into a likelihood with no mean
or variance.

Used by [`onset_reporting_model`](@ref) to score the reporting-triangle
increments, which are frequently negative (a later scan reads fewer cases
at some onset date than an earlier one, from digitisation noise rather than
a real reporting reversal) and so cannot take a count distribution.
"""
function safe_studentt(μ::Real, σ::Real, ν::Real)
    σc, _ = _studentt_scale(σ)
    νc, _ = _studentt_dof(ν)
    return μ + σc * TDist(νc)
end

## `safe_studentt`'s guards on the scale and the degrees of freedom, each
## returning the guarded value and whether the argument was kept.
_studentt_scale(σ) = _positive_or(σ, eps(typeof(float(σ))))
_studentt_dof(ν) = _positive_or(ν, oftype(float(ν), 4))

## One cell's term of a `StudentTVector`'s `logpdf`, given the normalising
## constant `c`, `(ν + 1) / 2` and the guarded `ν`, with the pieces its
## gradient reads: the standardised residual `z`, `log1p(z² / ν)` and the
## guarded scale with whether it was kept.
@inline function _studentt_cell(c, νp12, νc, μ, σ, x)
    σc, σ_on = _studentt_scale(σ)
    z = (x - μ) / σc
    l1 = log1p(z^2 / νc)
    return (; ℓ = (c - νp12 * l1) - log(σc), z, l1, σc, σ_on)
end

"""
Independent [`safe_studentt`](@ref) values, entry `i` about `μ[i]` with
scale `σ[i]` and the shared degrees of freedom `ν`. `logpdf` sums the
cells' terms in one loop, with the normalising constant evaluated once.
"""
struct StudentTVector{
        V <: AbstractVector{<:Real}, S <: AbstractVector{<:Real}, T <: Real,
    } <: Distributions.ContinuousMultivariateDistribution
    "Per-entry locations."
    μ::V
    "Per-entry scales."
    σ::S
    "Shared degrees of freedom."
    ν::T
end

Base.length(d::StudentTVector) = length(d.μ)
Base.eltype(::Type{<:StudentTVector}) = Float64

function Distributions._logpdf(d::StudentTVector, x::AbstractVector)
    T = float(promote_type(eltype(d.μ), eltype(d.σ), typeof(d.ν)))
    νc, _ = _studentt_dof(d.ν)
    νp12 = (νc + 1) / 2
    ## `TDist`'s log-density at zero is its normalising constant, so each
    ## cell's term below is the location-scale `logpdf` evaluated in full.
    c = logpdf(TDist(νc), zero(T))
    s = zero(T)
    @inbounds for i in eachindex(d.μ, d.σ, x)
        s += _studentt_cell(c, νp12, νc, d.μ[i], d.σ[i], x[i]).ℓ
    end
    return s
end

function Distributions._rand!(
        rng::AbstractRNG, d::StudentTVector, x::AbstractVector{<:Real}
    )
    @inbounds for i in eachindex(x, d.μ, d.σ)
        x[i] = rand(rng, safe_studentt(d.μ[i], d.σ[i], d.ν))
    end
    return x
end

## A `StudentTVector` is unconstrained, like each `safe_studentt` entry, so
## a sampled one links through the identity.
VectorBijectors.from_linked_vec(::StudentTVector) =
    VectorBijectors.TypedIdentity()
VectorBijectors.to_linked_vec(::StudentTVector) =
    VectorBijectors.TypedIdentity()
VectorBijectors.linked_vec_length(d::StudentTVector) = length(d)
VectorBijectors.linked_optic_vec(d::StudentTVector) =
    VectorBijectors.optic_vec(d)

"""
Independent counts split on their trial counts: entry `i` with
`trials[i] > 0` is a [`safe_betabinomial`](@ref) of `trials[i]` with mean
probability `p[i]` and overdispersion `ρ`, and every other entry a
[`safe_nbinomial`](@ref) about `μ[i]` with dispersion `k`. `logpdf` is a
[`BetaBinomialVector`](@ref) over the first group plus a
[`NegBinomialVector`](@ref) over the second.
"""
struct SplitCountVector{
        B <: BetaBinomialVector, N <: NegBinomialVector,
    } <: Distributions.DiscreteMultivariateDistribution
    "Entries with a trial count."
    anchored::Vector{Int}
    "Entries without one."
    unanchored::Vector{Int}
    "Distribution of the anchored entries."
    binomial::B
    "Distribution of the unanchored entries."
    negbinomial::N
end

function SplitCountVector(
        trials::AbstractVector{<:Integer}, p::AbstractVector, ρ::Real,
        k::Real, μ::AbstractVector
    )
    a = findall(>(0), trials)
    u = findall(<=(0), trials)
    return SplitCountVector(
        a, u, BetaBinomialVector(trials[a], p[a], ρ),
        NegBinomialVector(k, μ[u])
    )
end

Base.length(d::SplitCountVector) = length(d.anchored) + length(d.unanchored)
Base.eltype(::Type{<:SplitCountVector}) = Int

function Distributions._logpdf(d::SplitCountVector, x::AbstractVector)
    return logpdf(d.binomial, x[d.anchored]) +
        logpdf(d.negbinomial, x[d.unanchored])
end

function Distributions._rand!(
        rng::AbstractRNG, d::SplitCountVector, x::AbstractVector{<:Real}
    )
    x[d.anchored] = rand(rng, d.binomial)
    x[d.unanchored] = rand(rng, d.negbinomial)
    return x
end
