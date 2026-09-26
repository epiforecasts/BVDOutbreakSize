# Sum-to-zero deviations across patches, parameterised on the `n - 1`
# free directions of a sum-to-zero vector. Plain functions of their
# inputs, so they differentiate under Mooncake inside a Turing model.

"""
Orthonormal basis of the sum-to-zero subspace of `R^n`, as an
`n × (n - 1)` matrix `Q` with `Qᵀ Q = I` and `Qᵀ 1 = 0`.

Column `j` is the Helmert contrast of the first `j` entries against entry
`j + 1`,

```math
Q_{ij} = \\frac{1}{\\sqrt{j (j + 1)}} \\;(i \\le j), \\qquad
Q_{j+1, j} = \\frac{-j}{\\sqrt{j (j + 1)}},
```

zero below. Any vector `Q y` sums to zero, and `Q Qᵀ = I - J / n` is the
centring projector, so `σ Q z` with `z ~ N(0, I_{n-1})` has exactly the
distribution of a standard-normal `n`-vector scaled by `σ` and centred.
It does so with `n - 1` draws, one per direction a sum-to-zero vector can
move in. With `n = 1` the basis is empty (`1 × 0`).
"""
function sum_to_zero_basis(n::Integer)
    n >= 1 || throw(ArgumentError("sum_to_zero_basis: n = $n < 1"))
    Q = zeros(n, n - 1)
    for j in 1:(n - 1)
        c = 1 / sqrt(j * (j + 1))
        for i in 1:j
            Q[i, j] = c
        end
        Q[j + 1, j] = -j * c
    end
    return Q
end

"""
Loading matrix `F = Q diag(s) L` `(n × (n - 1))` of a sum-to-zero vector
`δ = F z` with `z ~ N(0, I_{n-1})`.

`Q` is [`sum_to_zero_basis`](@ref), `s` the scales of the `n - 1` basis
directions (a vector, or one scalar for an exchangeable vector) and `L` a
lower-triangular Cholesky factor of their `(n - 1) × (n - 1)` correlation,
or `nothing` for none. The covariance of `δ` is then

```math
\\Sigma = Q \\, \\mathrm{diag}(s) \\, L L^\\top \\mathrm{diag}(s) \\, Q^\\top,
```

a full covariance of the sum-to-zero vector with its `n (n - 1) / 2` free
parameters and no more. With a scalar `s` and no `L`, `Σ = s² (I - J / n)`,
the covariance of a centred vector of independent `N(0, s²)` draws.
"""
function sum_to_zero_factor(Q::AbstractMatrix, s, L = nothing)
    n, k = size(Q)
    T = promote_type(
        eltype(Q), eltype(s), isnothing(L) ? Bool : eltype(L)
    )
    F = zeros(T, n, k)
    @inbounds for i in 1:n, j in 1:k
        acc = zero(T)
        if isnothing(L)
            acc = Q[i, j] * _stz_scale(s, j)
        else
            for m in j:k
                acc += Q[i, m] * _stz_scale(s, m) * L[m, j]
            end
        end
        F[i, j] = acc
    end
    return F
end

@inline _stz_scale(s::Real, ::Integer) = s
@inline _stz_scale(s::AbstractVector, j::Integer) = @inbounds s[j]

"""
Lower-triangular `k × k` Bartlett factor with diagonal `d` and strictly
lower entries `o`, filled row by row.

With `d[j] ~ Chi(ν - j + 1)` and `o ~ N(0, I)` the product `A Aᵀ` is a
`Wishart(ν, I_k)` draw ([Bartlett decomposition](https://en.wikipedia.org/wiki/Wishart_distribution#Bartlett_decomposition)).
That prior is invariant under any rotation of the `k` directions, so with
[`sum_to_zero_factor`](@ref)`(Q, c, A)` the implied prior on the
sum-to-zero vector is the same whatever order its entries come in. `A`
is the unique Cholesky factor of `A Aᵀ`, so its `k (k + 1) / 2` entries
are exactly the free parameters of the covariance.
"""
function bartlett_factor(d::AbstractVector, o::AbstractVector)
    k = length(d)
    length(o) == k * (k - 1) ÷ 2 || throw(
        DimensionMismatch(
            "bartlett_factor: $(length(o)) lower entries for a $k × $k factor"
        )
    )
    A = zeros(promote_type(eltype(d), eltype(o)), k, k)
    m = 0
    @inbounds for i in 1:k
        A[i, i] = d[i]
        for j in 1:(i - 1)
            m += 1
            A[i, j] = o[m]
        end
    end
    return A
end

"""
Sum-to-zero vector `F z` for a loading matrix `F` from
[`sum_to_zero_factor`](@ref) and `n - 1` standard-normal draws `z`.
"""
function sum_to_zero(F::AbstractMatrix, z::AbstractVector)
    n, k = size(F)
    length(z) == k || throw(
        DimensionMismatch("sum_to_zero: $(length(z)) draws for $k directions")
    )
    T = promote_type(eltype(F), eltype(z))
    δ = zeros(T, n)
    @inbounds for i in 1:n
        acc = zero(T)
        for j in 1:k
            acc += F[i, j] * z[j]
        end
        δ[i] = acc
    end
    return δ
end

"""
Per-entry standard deviations and correlation matrix of the sum-to-zero
vector `F z`, from its covariance `Σ = F Fᵀ`. Returns `(; sd, cor)`.

The correlations of a sum-to-zero vector are constrained: its entries
cannot all be positively correlated, since `Σ 1 = 0` makes each row of `Σ`
sum to zero. With equal standard deviations each entry's correlations with
the others average `-1 / (n - 1)`, and they all equal it only when the
vector is exchangeable. For `n = 3` the three standard deviations determine the
three correlations exactly, `cor_{12} = (sd_3² - sd_1² - sd_2²) /
(2 sd_1 sd_2)`, and for larger `n` they constrain them. With `n = 1` the
vector is identically zero and the correlation is reported as one.
"""
function sum_to_zero_moments(F::AbstractMatrix)
    n, k = size(F)
    T = eltype(F)
    Σ = zeros(T, n, n)
    @inbounds for i in 1:n, j in 1:i
        acc = zero(T)
        for m in 1:k
            acc += F[i, m] * F[j, m]
        end
        Σ[i, j] = acc
        Σ[j, i] = acc
    end
    sd = [sqrt(Σ[i, i]) for i in 1:n]
    cor = ones(T, n, n)
    @inbounds for i in 1:n, j in 1:n
        i == j && continue
        d = sd[i] * sd[j]
        cor[i, j] = d > 0 ? Σ[i, j] / d : zero(T)
    end
    return (; sd, cor)
end
