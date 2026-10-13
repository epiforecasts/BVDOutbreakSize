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

"""
    flow_pair_basis(n)

Orthonormal bases, to rounding, of the double-centred flows between `n`
patches, as `(; symmetric, antisymmetric)` matrices of size `n² × d` whose
columns reshape to `n × n` flow matrices.

A double-centred flow matrix has a zero diagonal and zero row and column
sums. The symmetric part moves `q → p` and `p → q` together and has
`n (n - 3) / 2` directions. The antisymmetric part moves them in opposite
directions (circulations) and has `(n - 1)(n - 2) / 2`.
"""
function flow_pair_basis(n::Integer)
    n >= 1 || throw(ArgumentError("flow_pair_basis: n = $n < 1"))
    idx(p, q) = (q - 1) * n + p
    unit(ks...) = (v = zeros(n^2); foreach(k -> v[k] += 1, ks); v)
    centred = Vector{Float64}[]
    for p in 1:n
        push!(centred, unit(idx(p, p)))
        push!(centred, unit((idx(p, q) for q in 1:n)...))
        push!(centred, unit((idx(q, p) for q in 1:n)...))
    end
    ## `x_pq = x_qp` for the symmetric part, `x_pq = -x_qp` for the other.
    swapped(sgn) = [
        unit(idx(p, q)) .+ sgn .* unit(idx(q, p)) for q in 1:n for p in 1:(q - 1)
    ]
    symmetric = _orthonormal_complement(vcat(centred, swapped(-1)), n^2)
    antisymmetric = _orthonormal_complement(vcat(centred, swapped(1)), n^2)
    ## The diagonal is zero up to rounding; set it exactly.
    for p in 1:n
        symmetric[idx(p, p), :] .= 0
        antisymmetric[idx(p, p), :] .= 0
    end
    return (; symmetric, antisymmetric)
end

## Orthonormal basis (`m × d`) of the complement of `constraints` in `R^m`,
## by Gram-Schmidt in plain loops, since Mooncake has no rule for an SVD.
function _orthonormal_complement(constraints, m::Integer; tol = 1.0e-9)
    span = Vector{Float64}[]
    function project_out!(v)
        for b in span
            v .-= dot(b, v) .* b
        end
        nv = sqrt(dot(v, v))
        nv > tol || return false
        push!(span, v ./ nv)
        return true
    end
    foreach(c -> project_out!(copy(c)), constraints)
    n_constraints = length(span)
    for k in 1:m
        v = zeros(m)
        v[k] = 1
        project_out!(v)
    end
    d = length(span) - n_constraints
    basis = zeros(m, d)
    for j in 1:d
        basis[:, j] = span[n_constraints + j]
    end
    return basis
end

"""
    importation_flow_factor(Q, B, σ_dest, σ_flow, ρ_flow)

Non-centred factor `L` (`n² × (n - 1 + d)`) of the multivariate normal
prior on the log importation flows `ℓ = log K + W`. Gravity is the prior
mean of the flows, and the deviation `W = L ζ`, with `ζ ~ N(0, I)`, is their
data-driven departure from it:

```math
\\operatorname{vec}(\\ell) \\sim
    \\mathrm{MvNormal}(\\operatorname{vec}(\\log K),\\ \\Sigma),
\\qquad
\\Sigma = L L^\\top = \\sigma_{\\text{dest}}^2 C_{\\text{dest}}
    + \\sigma_{\\text{flow}}^2 R(\\rho_{\\text{flow}}),
```

before each origin's column is rescaled to its gravity total
([`destination_weighted_kernel`](@ref)). Row `(q - 1) n + p` is the flow
`q → p`, and the diagonal rows are zero.

The first `n - 1` columns are the destination effect, `σ_dest` times the
sum-to-zero basis `Q` ([`sum_to_zero_basis`](@ref)) repeated over the
origins, so `C_dest` has entries `δ_{pp'} - 1/n` between flows into `p`
and `p'`. The other `d` columns span the double-centred flows of `B`
([`flow_pair_basis`](@ref)), the symmetric directions first, scaled so
`R(ρ)` has unit diagonal on the flows and `ρ` between each flow and its
reverse: a correlation matrix. `Σ` is singular, since it leaves out the
shift shared within an origin's column that the rescale cancels, so it has
no Cholesky factor and `L` is a rectangular square root.
"""
function importation_flow_factor(
        Q::AbstractMatrix, B::NamedTuple, σ_dest, σ_flow, ρ_flow
    )
    S, A = B.symmetric, B.antisymmetric
    n, k = size(Q)
    ds, da = size(S, 2), size(A, 2)
    size(S, 1) == size(A, 1) == n^2 || throw(
        DimensionMismatch(
            "importation_flow_factor: a $n-patch basis and flow bases " *
                "with $(size(S, 1)) and $(size(A, 1)) rows"
        )
    )
    pairs = n * (n - 1) ÷ 2
    T = promote_type(
        eltype(Q), eltype(S), typeof(σ_dest), typeof(σ_flow), typeof(ρ_flow)
    )
    ## `M / d` of each part, so the two scales give each flow variance
    ## `σ_flow²` and correlation `ρ_flow` with its reverse.
    s_sym = ds > 0 ? σ_flow * sqrt((1 + ρ_flow) * pairs / ds) : zero(T)
    s_anti = da > 0 ? σ_flow * sqrt((1 - ρ_flow) * pairs / da) : zero(T)
    L = zeros(T, n^2, k + ds + da)
    @inbounds for q in 1:n, p in 1:n
        p == q && continue
        i = (q - 1) * n + p
        for j in 1:k
            L[i, j] = σ_dest * Q[p, j]
        end
        for j in 1:ds
            L[i, k + j] = s_sym * S[i, j]
        end
        for j in 1:da
            L[i, k + ds + j] = s_anti * A[i, j]
        end
    end
    return L
end

"""
Mean-reverting AR(1) knots of a sum-to-zero vector, one column per knot,
`size(Z, 2) + 1` in all,

```math
δ(1) = F_L z_L, \\qquad δ(k) = φ\\, δ(k - 1) + F_δ Z_{k-1},
```

with `F_L` and `F_δ` loading matrices from [`sum_to_zero_factor`](@ref),
`z_L` the level's draws and `Z` the innovation draws, one column per later
knot. A drift loading with no columns leaves the knots decaying along
`φ^{k-1} δ(1)`.
"""
function sum_to_zero_knots(
        F_level::AbstractMatrix, F_drift::AbstractMatrix,
        z_level::AbstractVector, Z::AbstractMatrix, φ::Real
    )
    T = promote_type(
        eltype(F_level), eltype(F_drift), eltype(z_level), eltype(Z),
        typeof(φ)
    )
    knots = zeros(T, size(F_level, 1), size(Z, 2) + 1)
    lvl = sum_to_zero(F_level, z_level)
    @inbounds for i in eachindex(lvl)
        knots[i, 1] = lvl[i]
    end
    return sum_to_zero_ar1!(knots, F_drift, Z, φ, 1)
end

"""
Fill knots `k0 + 1, …, k0 + m` of `knots` by the AR(1) of
[`sum_to_zero_knots`](@ref) from knot `k0`, with the `m` columns of `Z` as
the innovation draws through the loading `F`.
"""
function sum_to_zero_ar1!(
        knots::AbstractMatrix, F::AbstractMatrix, Z::AbstractMatrix,
        φ::Real, k0::Integer
    )
    innovations = F * Z
    @inbounds for j in axes(Z, 2), i in axes(knots, 1)
        knots[i, k0 + j] = φ * knots[i, k0 + j - 1] + innovations[i, j]
    end
    return knots
end
