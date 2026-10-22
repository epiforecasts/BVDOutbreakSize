# Native `Mooncake.rrule!!` methods for the `renewal.jl` kernels and the
# observation kernels in `models/observations.jl`. Each is a loop over a
# daily grid or an observation vector, so left to the backend every
# iteration's intermediates go on the tape; these replace that with a
# closed-form adjoint of the same shape as the forward loop. The kernels
# that are ComposableRecurrences operators bring their own rules.
#
# Differentiable arguments are restricted to `Array{<:IEEEFloat}` and
# `IEEEFloat` scalars, plus a view of one into such an array where a call
# site passes a row of a per-patch matrix. Anything else falls through to
# Mooncake's derived rule.
#
# These earn their place only while the backend's own derivation is worse.
# `task benchmark-rules` times both arms; if the gap has closed after a
# backend upgrade, delete this file rather than maintain it.

using Mooncake: CoDual, NoFData, NoRData, primal, tangent
using SpecialFunctions: digamma

## An array argument a rule accepts: a float array, or a view into one (a
## row of a per-patch matrix). A view's tangent is its parent's tangent, so
## `_tangent_array` views that the same way to accumulate into it.
const _FloatArray = Array{<:Mooncake.IEEEFloat}
const _FloatView = SubArray{<:Mooncake.IEEEFloat, 1, <:_FloatArray}
const _FloatVec = Union{_FloatArray, _FloatView}

_tangent_array(x::CoDual{<:_FloatArray}) = tangent(x)
function _tangent_array(x::CoDual{<:_FloatView})
    return view(tangent(x).data.parent, primal(x).indices...)
end

Mooncake.@is_primitive(
    Mooncake.MinimalCtx,
    Tuple{
        typeof(interpolate_knots), _FloatVec, Array{<:Integer}, Integer,
    },
)
## Data-only helpers: lookups over the observation grid and the recorded
## histories, with no sampled quantity among their arguments. They pass no
## derivative, so the backend need not tape them.
Mooncake.@zero_derivative(
    Mooncake.MinimalCtx, Tuple{typeof(admission_headroom), Vararg}
)
Mooncake.@zero_derivative(
    Mooncake.MinimalCtx, Tuple{typeof(censoring_cap), Vararg}
)
## Work that reaches only reported quantities (see `_detached`).
Mooncake.@zero_derivative(Mooncake.MinimalCtx, Tuple{typeof(_detached), Vararg})

function Mooncake.rrule!!(
        ::CoDual{typeof(interpolate_knots)},
        knot_vals::CoDual{<:_FloatVec},
        days::CoDual{<:Array{<:Integer}}, n::CoDual{<:Integer}
    )
    kp = primal(knot_vals)
    dayp = primal(days)
    np = primal(n)
    k̄ = _tangent_array(knot_vals)
    out = interpolate_knots(kp, dayp, np)
    ō = zero(out)
    function interpolate_knots_pullback!!(::NoRData)
        nb = length(dayp)
        if nb == 1
            @inbounds for t in 1:np
                k̄[1] += ō[t]
            end
            return NoRData(), NoRData(), NoRData(), NoRData()
        end
        Tf = eltype(ō)
        @inbounds for t in 1:np
            b, frac = _knot_bracket(dayp, t, Tf)
            g = ō[t]
            k̄[b] += (one(Tf) - frac) * g
            k̄[b + 1] += frac * g
        end
        return NoRData(), NoRData(), NoRData(), NoRData()
    end
    return CoDual(out, ō), interpolate_knots_pullback!!
end

Mooncake.@is_primitive(
    Mooncake.MinimalCtx,
    Tuple{
        typeof(euler_lotka_r), Mooncake.IEEEFloat,
        Array{<:Mooncake.IEEEFloat},
    },
)

## At the root `r` of `log R + log G(r) = 0`, with `G = Σ_s g_s e^{−r s}`
## and `dG = Σ_s s g_s e^{−r s}`, the implicit function theorem gives
##
##     ∂r/∂R = G / (R dG),    ∂r/∂g_s = e^{−r s} / dG.
##
## `R ≤ 0` returns `-Inf` and passes no derivative.
function Mooncake.rrule!!(
        ::CoDual{typeof(euler_lotka_r)}, R::CoDual{<:Mooncake.IEEEFloat},
        g::CoDual{<:Array{<:Mooncake.IEEEFloat}}
    )
    Rp = primal(R)
    gp = primal(g)
    ḡ = tangent(g)
    r = euler_lotka_r(Rp, gp)
    function euler_lotka_r_pullback!!(r̄)
        Rp > zero(Rp) || return NoRData(), zero(Rp), NoRData()
        G, dG = euler_lotka_sums(r, gp)
        a = r̄ / dG
        q = exp(-r)
        e = one(r)
        @inbounds for i in eachindex(gp)
            e *= q
            ḡ[i] += a * e
        end
        return NoRData(), oftype(Rp, a * G / Rp), NoRData()
    end
    return CoDual(r, NoFData()), euler_lotka_r_pullback!!
end

Mooncake.@is_primitive(
    Mooncake.MinimalCtx,
    Tuple{
        typeof(abscond_thinned), Array{<:Mooncake.IEEEFloat},
        Mooncake.IEEEFloat,
    },
)

## `y[i] = pmf[i] · s^(i−1)` with `s = 1 − κ`. The scalar `κ` has no fdata,
## so its cotangent goes back as rdata.
function Mooncake.rrule!!(
        ::CoDual{typeof(abscond_thinned)},
        pmf::CoDual{<:Array{<:Mooncake.IEEEFloat}},
        κ::CoDual{<:Mooncake.IEEEFloat}
    )
    pp = primal(pmf)
    κp = primal(κ)
    p̄ = tangent(pmf)
    y = abscond_thinned(pp, κp)
    ȳ = zero(y)
    function abscond_thinned_pullback!!(::NoRData)
        s = one(κp) - κp
        ## `pw` is `s^(i−1)` and `dpw` its derivative in `s`, both carried
        ## as running products so no power is taken.
        pw = one(s)
        dpw = zero(s)
        s̄ = zero(s)
        @inbounds for i in eachindex(pp)
            p̄[i] += ȳ[i] * pw
            s̄ += ȳ[i] * pp[i] * dpw
            dpw = dpw * s + pw
            pw *= s
        end
        return NoRData(), NoRData(), oftype(κp, -s̄)
    end
    return CoDual(y, ȳ), abscond_thinned_pullback!!
end

Mooncake.@is_primitive(
    Mooncake.MinimalCtx,
    Tuple{
        typeof(clinical_stay_survival), Array{<:Mooncake.IEEEFloat},
        Array{<:Mooncake.IEEEFloat}, Mooncake.IEEEFloat,
    },
)

function Mooncake.rrule!!(
        ::CoDual{typeof(clinical_stay_survival)},
        death_pmf::CoDual{<:Array{<:Mooncake.IEEEFloat}},
        recover_pmf::CoDual{<:Array{<:Mooncake.IEEEFloat}},
        cfr::CoDual{<:Mooncake.IEEEFloat}
    )
    dp = primal(death_pmf)
    rp = primal(recover_pmf)
    c = primal(cfr)
    d̄ = tangent(death_pmf)
    r̄ = tangent(recover_pmf)
    S = clinical_stay_survival(dp, rp, c)
    S̄ = zero(S)
    function clinical_stay_survival_pullback!!(::NoRData)
        ## `S[d] = 1 − Σ_{i ≤ d} (c · dp[i] + (1 − c) · rp[i])`, so entry
        ## `i` of either PMF feeds every survival weight at or above `i`:
        ## its adjoint is the reverse cumulative sum of `S̄`, scaled by its
        ## mixture weight.
        run = zero(eltype(S̄))
        c̄ = zero(eltype(S̄))
        @inbounds for i in length(S̄):-1:1
            run += S̄[i]
            pd = zero(eltype(dp))
            pr = zero(eltype(rp))
            if i <= length(dp)
                pd = dp[i]
                d̄[i] -= c * run
            end
            if i <= length(rp)
                pr = rp[i]
                r̄[i] -= (one(c) - c) * run
            end
            c̄ -= run * (pd - pr)
        end
        return NoRData(), NoRData(), NoRData(), convert(typeof(c), c̄)
    end
    return CoDual(S, S̄), clinical_stay_survival_pullback!!
end


Mooncake.@is_primitive(
    Mooncake.MinimalCtx,
    Tuple{
        typeof(incare_census), Array{<:Mooncake.IEEEFloat},
        Array{<:Mooncake.IEEEFloat}, Array{<:Mooncake.IEEEFloat},
        Mooncake.IEEEFloat, Array{<:Mooncake.IEEEFloat},
    },
)

## Adjoint of the in-care census. Per day, with `c = clamp(x, 0, O_bvd)`,
## `u = max(D − c, 0)` and `s = max(D + Δ − c, 0)`,
##
##     c̄ = ȳ_conf − s̄ [s side] − ū [u side],   ū = κ · ā[t + 1],
##     D̄ += s̄ [s side] + ū [u side],   Δ̄ += s̄ [s side] + t̄ot,
##     D̄ += t̄ot,   κ̄ += ā[t + 1] · u,
##
## and `c̄` goes to `x` or to `O_bvd` by the side the clamp took.
function Mooncake.rrule!!(
        ::CoDual{typeof(incare_census)},
        demand::CoDual{<:Array{<:Mooncake.IEEEFloat}},
        O_bvd::CoDual{<:Array{<:Mooncake.IEEEFloat}},
        O_conf_raw::CoDual{<:Array{<:Mooncake.IEEEFloat}},
        κ::CoDual{<:Mooncake.IEEEFloat},
        offset::CoDual{<:Array{<:Mooncake.IEEEFloat}}
    )
    κp = primal(κ)
    D̄ = tangent(demand)
    B̄ = tangent(O_bvd)
    X̄ = tangent(O_conf_raw)
    Δ̄ = tangent(offset)
    y, unconf, flags = _incare_census(
        Val(true), primal(demand), primal(O_bvd), primal(O_conf_raw), κp,
        primal(offset)
    )
    ȳ = map(zero, y)
    function incare_census_pullback!!(::NoRData)
        Tf = eltype(y.total)
        κ̄ = zero(Tf)
        n = length(flags)
        @inbounds for t in 1:n
            f = flags[t]
            gs = f & _CEN_SUSP != 0 ? ȳ.suspect[t] : zero(Tf)
            ā = t < n ? ȳ.abscond[t + 1] : zero(Tf)
            κ̄ += ā * unconf[t]
            gu = f & _CEN_UNCONF != 0 ? κp * ā : zero(Tf)
            gtot = ȳ.total[t] + gs
            D̄[t] += gtot + gu
            Δ̄[t] += gtot
            gc = ȳ.confirmed[t] - gs - gu
            if f & _CEN_CONF_X != 0
                X̄[t] += gc
            elseif f & _CEN_CONF_HI != 0
                B̄[t] += gc
            end
        end
        return NoRData(), NoRData(), NoRData(), NoRData(),
            convert(typeof(κp), κ̄), NoRData()
    end
    return CoDual(y, ȳ), incare_census_pullback!!
end

# Onset-reporting kernels from `models/observations.jl`. Each is a loop over
# the (delay, onset date) grid or the scored cells, so the derived rule tapes
# every logistic, product and division in it.

Mooncake.@is_primitive(
    Mooncake.MinimalCtx,
    Tuple{
        typeof(onset_report_cdf_table), Array{<:Mooncake.IEEEFloat},
        Array{<:Mooncake.IEEEFloat}, Integer, Integer, Integer,
    },
)
Mooncake.@is_primitive(
    Mooncake.MinimalCtx,
    Tuple{
        typeof(onset_report_anchor_series), Matrix{<:Mooncake.IEEEFloat},
        Integer, Array{<:Mooncake.IEEEFloat},
    },
)
Mooncake.@is_primitive(
    Mooncake.MinimalCtx,
    Tuple{
        typeof(onset_report_moments), Matrix{<:Mooncake.IEEEFloat},
        Integer, Array{<:Mooncake.IEEEFloat}, Integer,
        Array{<:Mooncake.IEEEFloat}, Array{<:Integer}, Array{<:Integer},
        Array{<:Integer},
    },
)
## Derivative of `safe_rate`: the identity above its floor, flat on it.
_safe_rate_slope(x) = _safe_rate_on(x) ? one(x) : zero(x)

## Adjoint of one `_onset_columns!` column, with `c̄` the cotangent of
## `1 - surv`. `1 - surv[j]` has derivative `surv[j] · h[i]` in the logit
## `x_i` for every `i ≤ j`, so a column's adjoint is
## `x̄_i = h_i · Σ_{j ≥ i} c̄_j surv_j`, one backward running sum. Each `x̄_i`
## lands on `logit_h0[i]` and on the (clamped) calendar day it read from `γ`.
function _onset_column_adjoint!(
        l̄::AbstractVector, γ̄::AbstractVector, c̄::AbstractVector,
        h::AbstractVector, surv::AbstractVector, u::Integer,
        grid_start::Integer
    )
    D = length(h)
    ng = length(γ̄)
    acc = zero(eltype(l̄))
    @inbounds for j in (D - 1):-1:0
        acc += c̄[j + 1] * surv[j + 1]
        x̄ = h[j + 1] * acc
        l̄[j + 1] += x̄
        γ̄[clamp(u + j - grid_start + 1, 1, ng)] += x̄
    end
    return nothing
end

function Mooncake.rrule!!(
        ::CoDual{typeof(onset_report_cdf_table)},
        logit_h0::CoDual{<:Array{<:Mooncake.IEEEFloat}},
        γ::CoDual{<:Array{<:Mooncake.IEEEFloat}},
        grid_start::CoDual{<:Integer}, u_lo::CoDual{<:Integer},
        u_hi::CoDual{<:Integer}
    )
    lp = primal(logit_h0)
    γp = primal(γ)
    gs = primal(grid_start)
    lo = Int(primal(u_lo))
    l̄ = tangent(logit_h0)
    γ̄ = tangent(γ)
    ## The forward pass keeps each column's hazards and survival products,
    ## so the pullback reuses them rather than re-evaluating the logistic.
    D = length(lp)
    nu = max(Int(primal(u_hi)) - lo + 1, 0)
    T = promote_type(eltype(lp), eltype(γp))
    H = Matrix{T}(undef, D, nu)
    S = Matrix{T}(undef, D, nu)
    _onset_columns!(H, S, lp, γp, gs, lo)
    table = one(T) .- S
    out = Mooncake.zero_fcodual(table)
    t̄ = tangent(out)
    function onset_report_cdf_table_pullback!!(::NoRData)
        ## `table = 1 - surv`, so the column adjoint is `t̄` itself.
        @inbounds for k in 1:nu
            _onset_column_adjoint!(
                l̄, γ̄, view(t̄, :, k), view(H, :, k), view(S, :, k),
                lo + k - 1, gs
            )
        end
        return ntuple(_ -> NoRData(), 6)
    end
    return out, onset_report_cdf_table_pullback!!
end

function Mooncake.rrule!!(
        ::CoDual{typeof(onset_report_anchor_series)},
        cdf_table::CoDual{<:Matrix{<:Mooncake.IEEEFloat}},
        u_lo::CoDual{<:Integer}, a::CoDual{<:Array{<:Mooncake.IEEEFloat}}
    )
    cp = primal(cdf_table)
    ap = primal(a)
    lo = Int(primal(u_lo))
    c̄ = tangent(cdf_table)
    ā = tangent(a)
    y = onset_report_anchor_series(cp, lo, ap)
    out = Mooncake.zero_fcodual(y)
    ȳ = tangent(out)
    ## `out[k] = S_k / den_k` with `S_k = Σ_d (c_d - c_{d-1}) a_{u+d}` and
    ## `den_k = safe_rate(c_{D-1})`, so `c_d` enters `S_k` through two
    ## neighbouring anchor days and `c_{D-1}` also through the denominator.
    function onset_report_anchor_series_pullback!!(::NoRData)
        D = size(cp, 1)
        na = length(ap)
        D == 0 && return NoRData(), NoRData(), NoRData(), NoRData()
        @inbounds for k in axes(cp, 2)
            g = ȳ[k]
            u = lo + k - 1
            cD = cp[D, k]
            invden = inv(safe_rate(cD))
            gi = g * invden
            c_prev = zero(eltype(cp))
            S = zero(eltype(cp))
            for d in 0:(D - 1)
                cd = cp[d + 1, k]
                ad = ap[clamp(u + d, 1, na)]
                a_next = d < D - 1 ? ap[clamp(u + d + 1, 1, na)] : zero(ad)
                c̄[d + 1, k] += gi * (ad - a_next)
                ā[clamp(u + d, 1, na)] += gi * (cd - c_prev)
                S += (cd - c_prev) * ad
                c_prev = cd
            end
            c̄[D, k] -= gi * S * invden * _safe_rate_slope(cD)
        end
        return NoRData(), NoRData(), NoRData(), NoRData()
    end
    return out, onset_report_anchor_series_pullback!!
end

function Mooncake.rrule!!(
        ::CoDual{typeof(onset_report_moments)},
        cdf_table::CoDual{<:Matrix{<:Mooncake.IEEEFloat}},
        u_lo::CoDual{<:Integer},
        onsets::CoDual{<:Array{<:Mooncake.IEEEFloat}},
        grid_start::CoDual{<:Integer},
        alpha::CoDual{<:Array{<:Mooncake.IEEEFloat}},
        onset_idx::CoDual{<:Array{<:Integer}},
        cur_report_idx::CoDual{<:Array{<:Integer}},
        prev_report_idx::CoDual{<:Array{<:Integer}}
    )
    cp = primal(cdf_table)
    lo = Int(primal(u_lo))
    op = primal(onsets)
    gs = Int(primal(grid_start))
    alp = primal(alpha)
    oi = primal(onset_idx)
    ci = primal(cur_report_idx)
    pri = primal(prev_report_idx)
    c̄ = tangent(cdf_table)
    ō = tangent(onsets)
    ᾱ = tangent(alpha)
    y = onset_report_moments(cp, lo, op, gs, alp, oi, ci, pri)
    out = Mooncake.zero_fcodual(y)
    ȳ = tangent(out)
    ## Per cell, `level = onset · α · num / den` for the current and the
    ## previous report date, and `means = level_cur - level_prev`, so the
    ## `means` cotangent folds into the two level cotangents first.
    function onset_report_moments_pullback!!(::NoRData)
        D = size(cp, 1)
        n = length(op)
        na = length(alp)
        T = eltype(cp)
        @inbounds for i in eachindex(oi)
            gm = ȳ.means[i]
            gc = ȳ.level_cur[i] + gm
            gp = ȳ.level_prev[i] - gm
            u = oi[i]
            k = u - lo + 1
            in_range = u >= 1 && u <= n
            onset_rate = in_range ? op[u] : zero(T)
            ia = clamp(u - gs + 1, 1, na)
            α = alp[ia]
            δc = ci[i] - u
            δp = pri[i] - u
            jc = _onset_delay_row(δc, D)
            jp = _onset_delay_row(δp, D)
            num_c = jc == 0 ? zero(T) : cp[jc, k]
            num_p = jp == 0 ? zero(T) : cp[jp, k]
            cD = D > 0 ? cp[D, k] : zero(T)
            sden = safe_rate(cD)
            r = (num_c * gc + num_p * gp) / sden
            in_range && (ō[u] += α * r)
            ᾱ[ia] += onset_rate * r
            w = onset_rate * α / sden
            jc > 0 && (c̄[jc, k] += w * gc)
            jp > 0 && (c̄[jp, k] += w * gp)
            D > 0 && (
                c̄[D, k] -= onset_rate * α * r / sden *
                    _safe_rate_slope(cD)
            )
        end
        return ntuple(_ -> NoRData(), 9)
    end
    return out, onset_report_moments_pullback!!
end

# Log densities of the vector observation distributions in
# `models/observation_distributions.jl`. Each rule takes the distribution
# itself: the tangent of an array field is read from the distribution's
# fdata and accumulated in place, and the cotangent of a scalar field goes
# back in its rdata.

## `logpdf(NegBinomialVector(k, μ), x)` in one pass. With `r = k` and
## `p = r / (r + μ)` each count `x` contributes
##
##     ∂ℓ/∂r = log p + ψ(r + x) − ψ(r),    ∂ℓ/∂p = r / p − x / (1 − p),
##
## chained through `p` to `μ` and `k`. The guards are the ones `safe_rate`
## and `safe_nbinomial` apply: a floored `k` or `μ`, or a clamped `p`,
## passes no derivative. A term that is not finite adds none either.
function _negbinomial_vector_grad(
        k::T, modelled::AbstractVector,
        obs::AbstractVector
    ) where {T}
    _, k_on = _positive_or(k, eps(T))
    s = zero(T)
    dk = zero(T)
    dμ = zeros(T, length(modelled))
    @inbounds for i in eachindex(modelled, obs)
        μ = modelled[i]
        x = obs[i]
        (; r, m, p, p_on) = _nbinomial_params(k, safe_rate(μ))
        ℓ = logpdf(NegativeBinomial(r, p), x)
        s += ℓ
        isfinite(ℓ) || continue
        ∂r = iszero(x) ? log(p) : log(p) + digamma(r + x) - digamma(r)
        if p_on
            ∂p = r / p - x / (one(T) - p)
            den = (r + m)^2
            ∂r += ∂p * m / den
            _safe_rate_on(μ) && (dμ[i] = -∂p * r / den)
        end
        dk += ∂r
    end
    return s, (k_on ? dk : zero(T)), dμ
end

Mooncake.@is_primitive(
    Mooncake.MinimalCtx,
    Tuple{
        typeof(Distributions._logpdf),
        NegBinomialVector{<:Mooncake.IEEEFloat, <:_FloatArray},
        Array{<:Integer},
    },
)

function Mooncake.rrule!!(
        ::CoDual{typeof(Distributions._logpdf)},
        d::CoDual{<:NegBinomialVector{<:Mooncake.IEEEFloat, <:_FloatArray}},
        obs::CoDual{<:Array{<:Integer}}
    )
    dist = primal(d)
    s, dk, dμ = _negbinomial_vector_grad(dist.k, dist.μ, primal(obs))
    μ̄ = tangent(d).data.μ
    function negbinomial_vector_pullback!!(s̄)
        μ̄ .+= s̄ .* dμ
        return NoRData(), Mooncake.RData((k = s̄ * dk, μ = NoRData())),
            NoRData()
    end
    return CoDual(s, NoFData()), negbinomial_vector_pullback!!
end

## The increments a `StudentTVector` scores: counts in the package data,
## floats in simulated triangles. Float increments take the adjoint
## `−∂ℓ/∂μ`.
const _StudentTObs = Union{Array{<:Integer}, _FloatArray}
const _StudentTVec = StudentTVector{
    <:_FloatArray, <:_FloatArray, <:Mooncake.IEEEFloat,
}

Mooncake.@is_primitive(
    Mooncake.MinimalCtx,
    Tuple{typeof(Distributions._logpdf), _StudentTVec, _StudentTObs},
)

## `logpdf(StudentTVector(μ, σ, ν), x)` in one pass. With
## `z = (x − μ) / σ` and `g = (ν + 1) z / (ν + z²)` each cell contributes
##
##     ∂ℓ/∂μ = g / σ,    ∂ℓ/∂σ = (g z − 1) / σ,
##     ∂ℓ/∂ν = (ψ((ν + 1)/2) − ψ(ν/2) − 1/ν
##              − log1p(z²/ν) + g z / ν) / 2,
##
## and `∂ℓ/∂x = −∂ℓ/∂μ`. The guards are `safe_studentt`'s: a floored `σ`
## or a defaulted `ν` passes no derivative. A cell whose term is not finite
## adds none either.
function _studentt_vector_grad(
        means::AbstractVector, sds::AbstractVector,
        obs::AbstractVector, ν::Real
    )
    T = float(promote_type(eltype(means), eltype(sds), typeof(ν)))
    νc, ν_on = _studentt_dof(ν)
    νp12 = (νc + 1) / 2
    c = logpdf(TDist(νc), zero(T))
    ∂ν_c = (digamma(νp12) - digamma(νc / 2) - 1 / νc) / 2
    s = zero(T)
    dν = zero(T)
    dμ = zeros(T, length(means))
    dσ = zeros(T, length(means))
    @inbounds for i in eachindex(means, sds, obs)
        (; ℓ, z, l1, σc, σ_on) = _studentt_cell(
            c, νp12, νc, means[i], sds[i], obs[i]
        )
        s += ℓ
        isfinite(ℓ) || continue
        g = (νc + 1) * z / (νc + z^2)
        dμ[i] = g / σc
        σ_on && (dσ[i] = (g * z - 1) / σc)
        dν += ∂ν_c - l1 / 2 + g * z / (2 * νc)
    end
    return s, (ν_on ? dν : zero(T)), dμ, dσ
end

function Mooncake.rrule!!(
        ::CoDual{typeof(Distributions._logpdf)},
        d::CoDual{<:_StudentTVec},
        obs::CoDual{<:_StudentTObs}
    )
    dist = primal(d)
    s, dν, dμ, dσ = _studentt_vector_grad(
        dist.μ, dist.σ, primal(obs), dist.ν
    )
    μ̄ = tangent(d).data.μ
    σ̄ = tangent(d).data.σ
    x̄ = tangent(obs)
    ## Integer increments carry no tangent.
    function studentt_vector_pullback!!(s̄)
        μ̄ .+= s̄ .* dμ
        σ̄ .+= s̄ .* dσ
        x̄ isa _FloatArray && (x̄ .-= s̄ .* dμ)
        d̄ = Mooncake.RData((μ = NoRData(), σ = NoRData(), ν = s̄ * dν))
        return NoRData(), d̄, NoRData()
    end
    return CoDual(s, NoFData()), studentt_vector_pullback!!
end

const _BetaBinomialVec = BetaBinomialVector{
    <:AbstractVector{<:Integer}, <:_FloatArray, <:Mooncake.IEEEFloat,
}

Mooncake.@is_primitive(
    Mooncake.MinimalCtx,
    Tuple{
        typeof(Distributions._logpdf), _BetaBinomialVec,
        AbstractVector{<:Integer},
    },
)

## `logpdf(BetaBinomialVector(n, p, ρ), x)` in one pass. With
## concentration `c = (1 − ρ) / ρ`, `α = c·p` and `β = c·(1 − p)`, each
## count `x` of `n` trials contributes
##
##     ∂ℓ/∂α = ψ(x + α) − ψ(α) + ψ(α + β) − ψ(n + α + β),
##     ∂ℓ/∂β = ψ(n − x + β) − ψ(β) + ψ(α + β) − ψ(n + α + β),
##
## chained through `α` and `β` to `p` and `c`, and through `c` to `ρ`. The
## guards are `safe_betabinomial`'s: a clamped `p` or `ρ`, or an `α` or `β`
## held at its floor, passes no derivative. A term that is not finite adds
## none either.
function _betabinomial_vector_grad(
        trials::AbstractVector, p::AbstractVector, ρ, obs::AbstractVector
    )
    T = float(promote_type(eltype(p), typeof(ρ)))
    conc = _betabinomial_concentration(T, ρ)
    c = conc.s
    s = zero(T)
    dc = zero(T)
    dp = zeros(T, length(p))
    @inbounds for i in eachindex(trials, p, obs)
        n = trials[i]
        x = obs[i]
        (; α, β, pc, p_on, α_on, β_on) = _betabinomial_shapes(p[i], c)
        ℓ = logpdf(BetaBinomial(n, α, β), x)
        s += ℓ
        isfinite(ℓ) || continue
        ψ_ab = digamma(α + β) - digamma(n + α + β)
        ∂α = α_on ? digamma(x + α) - digamma(α) + ψ_ab : zero(T)
        ∂β = β_on ? digamma(n - x + β) - digamma(β) + ψ_ab : zero(T)
        p_on && (dp[i] = c * (∂α - ∂β))
        dc += ∂α * pc + ∂β * (one(T) - pc)
    end
    return s, (conc.ρ_on ? -dc / conc.ρc^2 : zero(T)), dp
end

function Mooncake.rrule!!(
        ::CoDual{typeof(Distributions._logpdf)},
        d::CoDual{<:_BetaBinomialVec},
        obs::CoDual{<:AbstractVector{<:Integer}}
    )
    dist = primal(d)
    ## The gradient is taken in the forward pass, since a caller may
    ## overwrite `p` in place before the pullback runs.
    s, dρ, dp = _betabinomial_vector_grad(
        dist.trials, dist.p, dist.ρ, primal(obs)
    )
    p̄ = tangent(d).data.p
    function betabinomial_vector_pullback!!(s̄)
        p̄ .+= s̄ .* dp
        d̄ = Mooncake.RData((trials = NoRData(), p = NoRData(), ρ = s̄ * dρ))
        return NoRData(), d̄, NoRData()
    end
    return CoDual(s, NoFData()), betabinomial_vector_pullback!!
end

## The health-zone window weights the model rebuilds on every draw's
## delays, from the cumulative PMF.
Mooncake.@is_primitive(
    Mooncake.MinimalCtx,
    Tuple{
        typeof(zone_window_weights), Vector{<:Mooncake.IEEEFloat},
        Vector{<:Integer}, Vector{<:Integer}, Integer,
    },
)

function Mooncake.rrule!!(
        ::CoDual{typeof(zone_window_weights)},
        f::CoDual{<:Vector{<:Mooncake.IEEEFloat}},
        starts::CoDual{<:Vector{<:Integer}},
        stops::CoDual{<:Vector{<:Integer}}, nd::CoDual{<:Integer}
    )
    fp, a_, b_ = primal(f), primal(starts), primal(stops)
    f̄ = tangent(f)
    y = zone_window_weights(fp, a_, b_, primal(nd))
    ȳ = zero(y)
    function zone_window_weights_pullback!!(::NoRData)
        ## Each weight is a difference of two entries of the cumulative PMF
        ## `cf`, so the adjoint gathers onto `cf` and then reverses the
        ## cumulative sum onto `f`.
        L = length(fp)
        c̄ = zeros(eltype(fp), L)
        @inbounds for v in eachindex(a_)
            a, b = a_[v], b_[v]
            a <= b || continue
            for k in max(1, a - L + 1):b
                c̄[min(b - k, L - 1) + 1] += ȳ[k, v]
                a - 1 - k < 0 || (c̄[min(a - 1 - k, L - 1) + 1] -= ȳ[k, v])
            end
        end
        acc = zero(eltype(fp))
        @inbounds for i in L:-1:1
            acc += c̄[i]
            f̄[i] += acc
        end
        return NoRData(), NoRData(), NoRData(), NoRData(), NoRData()
    end
    return CoDual(y, ȳ), zone_window_weights_pullback!!
end
