# Vector-valued observation distributions: one per stream kernel, so a
# submodel scores a whole observed vector with a single `~` and draws a
# `missing` one with the same statement. `logpdf` is the summed helper in
# `models/observations.jl`, which carries a closed-form Mooncake rule
# (`src/ad_rules.jl`), and `rand` draws each entry from the matching scalar
# `safe_*` distribution. The vector is one variable, so a sampled vector is
# stored under one whole-vector key (`<prefix>.increments`, `<prefix>.obs`).

"""
Independent [`safe_nbinomial`](@ref) counts, entry `i` about the
[`safe_rate`](@ref) of `μ[i]` with the shared dispersion `k`. `logpdf` is
[`nbinomial_loglik`](@ref).
"""
struct NegBinomialVector{T <: Real, V <: AbstractVector{<:Real}} <:
    Distributions.DiscreteMultivariateDistribution
    "Shared dispersion."
    k::T
    "Per-entry means."
    μ::V
end

"""
Independent right-censored [`safe_nbinomial`](@ref) counts, entry `i` about
`μ[i]` and censored at `upper[i]`, both through [`safe_rate`](@ref), with
the shared dispersion `k`. `logpdf` is
[`censored_nbinomial_loglik`](@ref).
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
Independent [`safe_studentt`](@ref) values, entry `i` about `μ[i]` with
scale `σ[i]` and the shared degrees of freedom `ν`. `logpdf` is
[`studentt_loglik`](@ref).
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

Base.length(
    d::Union{NegBinomialVector, CensoredNegBinomialVector, StudentTVector}
) = length(d.μ)
Base.eltype(::Type{<:NegBinomialVector}) = Int
## A draw above a ceiling returns the ceiling, which need not be a whole
## number (`admission_headroom`), so censored draws are floats.
Base.eltype(::Type{<:CensoredNegBinomialVector}) = Float64
Base.eltype(::Type{<:StudentTVector}) = Float64

function Distributions._logpdf(d::NegBinomialVector, x::AbstractVector)
    return nbinomial_loglik(d.k, d.μ, x)
end
function Distributions._logpdf(
        d::CensoredNegBinomialVector, x::AbstractVector
    )
    return censored_nbinomial_loglik(d.k, d.μ, d.upper, x)
end
function Distributions._logpdf(d::StudentTVector, x::AbstractVector)
    return studentt_loglik(d.μ, d.σ, x, d.ν)
end

function Distributions._rand!(
        rng::AbstractRNG, d::NegBinomialVector, x::AbstractVector{<:Real}
    )
    @inbounds for i in eachindex(x, d.μ)
        x[i] = rand(rng, safe_nbinomial(d.k, safe_rate(d.μ[i])))
    end
    return x
end
function Distributions._rand!(
        rng::AbstractRNG, d::CensoredNegBinomialVector,
        x::AbstractVector{<:Real}
    )
    @inbounds for i in eachindex(x, d.μ, d.upper)
        x[i] = rand(
            rng, censored(
                safe_nbinomial(d.k, safe_rate(d.μ[i]));
                upper = safe_rate(d.upper[i])
            )
        )
    end
    return x
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
