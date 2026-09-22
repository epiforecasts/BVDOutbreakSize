# Native `Mooncake.rrule!!` methods for the `renewal.jl` kernels. Each is a
# loop over a daily grid, so left to the backend every iteration's
# intermediates go on the tape; these replace that with a closed-form
# adjoint of the same shape as the forward loop.
#
# Signatures are restricted to `Array{<:IEEEFloat}`, what every call site
# passes. Anything else falls through to Mooncake's derived rule.
#
# These earn their place only while the backend's own derivation is worse.
# `task benchmark-rules` times both arms; if the gap has closed after a
# backend upgrade, delete this file rather than maintain it.

using Mooncake: CoDual, NoRData, primal, tangent, zero_fcodual

Mooncake.@is_primitive(
    Mooncake.MinimalCtx,
    Tuple{
        typeof(convolve_delay), Array{<:Mooncake.IEEEFloat},
        Array{<:Mooncake.IEEEFloat},
    },
)
Mooncake.@is_primitive(
    Mooncake.MinimalCtx,
    Tuple{
        typeof(convolve_survival), Array{<:Mooncake.IEEEFloat},
        Array{<:Mooncake.IEEEFloat},
    },
)
Mooncake.@is_primitive(
    Mooncake.MinimalCtx,
    Tuple{
        typeof(convolve_pmf), Array{<:Mooncake.IEEEFloat},
        Array{<:Mooncake.IEEEFloat},
    },
)
Mooncake.@is_primitive(
    Mooncake.MinimalCtx,
    Tuple{
        typeof(interpolate_knots), Array{<:Mooncake.IEEEFloat},
        Array{<:Integer}, Integer,
    },
)
Mooncake.@is_primitive(
    Mooncake.MinimalCtx,
    Tuple{
        typeof(renewal_infections), Array{<:Mooncake.IEEEFloat},
        Array{<:Mooncake.IEEEFloat}, Array{<:Mooncake.IEEEFloat},
    },
)

## Adjoint of the daily delay convolution, shared by `convolve_delay` and
## `convolve_survival`. A convolution's pullback is the matching
## correlation, one pass over the same `(t, d)` pairs:
##
##     y[t]    = Σ_d x[t−d] · w[d+1]
##     x̄[s]   += Σ_d ȳ[s+d] · w[d+1]
##     w̄[d+1] += Σ_t ȳ[t] · x[t−d]
function _convolve_delay_adjoint(
        ȳ::AbstractVector, x::AbstractVector,
        delay::AbstractVector
    )
    n = length(x)
    D = length(delay)
    x̄ = zeros(promote_type(eltype(ȳ), eltype(delay)), n)
    d̄ = zeros(promote_type(eltype(ȳ), eltype(x)), D)
    @inbounds for t in 1:min(n, length(ȳ))
        g = ȳ[t]
        dmax = min(t - 1, D - 1)
        for d in 0:dmax
            x̄[t - d] += g * delay[d + 1]
            d̄[d + 1] += g * x[t - d]
        end
    end
    return x̄, d̄
end

function Mooncake.rrule!!(
        ::CoDual{typeof(convolve_delay)},
        x::CoDual{<:Array{<:Mooncake.IEEEFloat}},
        delay::CoDual{<:Array{<:Mooncake.IEEEFloat}}
    )
    xp = primal(x)
    dp = primal(delay)
    x̄ = tangent(x)
    d̄ = tangent(delay)
    y = convolve_delay(xp, dp)
    ȳ = zero(y)
    function convolve_delay_pullback!!(::NoRData)
        Δx, Δd = _convolve_delay_adjoint(ȳ, xp, dp)
        x̄ .+= Δx
        d̄ .+= Δd
        return NoRData(), NoRData(), NoRData()
    end
    return CoDual(y, ȳ), convolve_delay_pullback!!
end

function Mooncake.rrule!!(
        ::CoDual{typeof(convolve_survival)},
        x::CoDual{<:Array{<:Mooncake.IEEEFloat}},
        los::CoDual{<:Array{<:Mooncake.IEEEFloat}}
    )
    xp = primal(x)
    lp = primal(los)
    x̄ = tangent(x)
    l̄ = tangent(los)
    surv = survival_weights(lp)
    y = convolve_delay(xp, surv)
    ȳ = zero(y)
    function convolve_survival_pullback!!(::NoRData)
        Δx, s̄ = _convolve_delay_adjoint(ȳ, xp, surv)
        x̄ .+= Δx
        ## `surv[i] = Σ_{j ≥ i} los[j]`, so `los[j]` feeds every survival
        ## weight at or below `j`: the adjoint is the forward cumulative
        ## sum of the survival adjoint.
        run = zero(eltype(s̄))
        @inbounds for j in eachindex(s̄)
            run += s̄[j]
            l̄[j] += run
        end
        return NoRData(), NoRData(), NoRData()
    end
    return CoDual(y, ȳ), convolve_survival_pullback!!
end

function Mooncake.rrule!!(
        ::CoDual{typeof(convolve_pmf)},
        a::CoDual{<:Array{<:Mooncake.IEEEFloat}},
        b::CoDual{<:Array{<:Mooncake.IEEEFloat}}
    )
    ap = primal(a)
    bp = primal(b)
    ā = tangent(a)
    b̄ = tangent(b)
    y = convolve_pmf(ap, bp)
    ȳ = zero(y)
    function convolve_pmf_pullback!!(::NoRData)
        na = length(ap)
        nb = length(bp)
        @inbounds for i in 1:na, j in 1:nb

            g = ȳ[i + j - 1]
            ā[i] += g * bp[j]
            b̄[j] += g * ap[i]
        end
        return NoRData(), NoRData(), NoRData()
    end
    return CoDual(y, ȳ), convolve_pmf_pullback!!
end

function Mooncake.rrule!!(
        ::CoDual{typeof(interpolate_knots)},
        knot_vals::CoDual{<:Array{<:Mooncake.IEEEFloat}},
        days::CoDual{<:Array{<:Integer}}, n::CoDual{<:Integer}
    )
    kp = primal(knot_vals)
    dayp = primal(days)
    np = primal(n)
    k̄ = tangent(knot_vals)
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
            b = 1
            while b < nb - 1 && t > dayp[b + 1]
                b += 1
            end
            d0 = dayp[b]
            d1 = dayp[b + 1]
            frac = d1 == d0 ? zero(Tf) :
                clamp(Tf(t - d0) / Tf(d1 - d0), zero(Tf), one(Tf))
            g = ō[t]
            k̄[b] += (one(Tf) - frac) * g
            k̄[b + 1] += frac * g
        end
        return NoRData(), NoRData(), NoRData(), NoRData()
    end
    return CoDual(out, ō), interpolate_knots_pullback!!
end

function Mooncake.rrule!!(
        ::CoDual{typeof(renewal_infections)},
        Rt::CoDual{<:Array{<:Mooncake.IEEEFloat}},
        g::CoDual{<:Array{<:Mooncake.IEEEFloat}},
        seed::CoDual{<:Array{<:Mooncake.IEEEFloat}}
    )
    Rtp = primal(Rt)
    gp = primal(g)
    seedp = primal(seed)
    R̄ = tangent(Rt)
    ḡ = tangent(g)
    s̄ = tangent(seed)
    n = length(Rtp)
    L = length(seedp)
    ## The forward pass is the model's own, which also hands back the
    ## per-day force the pullback needs, so the recursion is defined once.
    I, force = renewal_infections_with_force(Rtp, gp, seedp)
    Ī = zero(I)
    function renewal_infections_pullback!!(::NoRData)
        ## Sequential recursion, so the walk runs backwards: each day's
        ## adjoint must land on the lagged infections before those days
        ## are read. `Ī` is Mooncake's buffer, so accumulate into a copy.
        acc = copy(Ī)
        @inbounds for t in n:-1:(L + 1)
            it = acc[t]
            R̄[t] += it * force[t]
            f̄ = it * Rtp[t]
            kmax = min(t - 1, length(gp))
            for s in 1:kmax
                acc[t - s] += f̄ * gp[s]
                ḡ[s] += f̄ * I[t - s]
            end
        end
        @inbounds for j in 1:min(L, n)
            s̄[j] += acc[j]
        end
        return NoRData(), NoRData(), NoRData(), NoRData()
    end
    return CoDual(I, Ī), renewal_infections_pullback!!
end
