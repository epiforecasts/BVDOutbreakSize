# Hand-written reverse-mode derivative rules for the numeric kernels in
# `renewal.jl`.
#
# The rules are plain `ChainRulesCore.rrule` methods with explicit
# pullbacks, so they carry nothing backend-specific. Mooncake picks them up
# through the `Mooncake.@from_rrule` registrations at the bottom of this
# file; the Enzyme extension imports the same methods with
# `Enzyme.@import_rrule`, so both backends differentiate one derivation.
#
# Each kernel is a loop over a daily grid. Left to a backend, every
# iteration's intermediates go on the tape. The rules below replace that
# with a closed-form adjoint loop of the same shape as the forward one.

## Adjoint of the daily delay convolution, shared by `convolve_delay` and
## `convolve_survival`. For
##
##     y[t] = Σ_d x[t−d] · delay[d+1]    (d ≥ 0, t−d ≥ 1)
##
## the two input adjoints are
##
##     x̄[s]          += Σ_d ȳ[s+d] · delay[d+1]
##     delaybar[d+1] += Σ_t ȳ[t] · x[t−d]
##
## which is one pass over the same `(t, d)` pairs the forward loop walks.
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

function ChainRulesCore.rrule(
        ::typeof(convolve_delay), x::AbstractVector,
        delay::AbstractVector
    )
    y = convolve_delay(x, delay)
    function convolve_delay_pullback(Δy)
        ȳ = ChainRulesCore.unthunk(Δy)
        x̄, d̄ = _convolve_delay_adjoint(ȳ, x, delay)
        return ChainRulesCore.NoTangent(), x̄, d̄
    end
    return y, convolve_delay_pullback
end

function ChainRulesCore.rrule(
        ::typeof(convolve_survival), x::AbstractVector,
        los::AbstractVector
    )
    L = length(los)
    surv = similar(los)
    acc = zero(eltype(los))
    @inbounds for i in L:-1:1
        acc += los[i]
        surv[i] = acc
    end
    y = convolve_delay(x, surv)
    function convolve_survival_pullback(Δy)
        ȳ = ChainRulesCore.unthunk(Δy)
        x̄, s̄ = _convolve_delay_adjoint(ȳ, x, surv)
        ## `surv[i] = Σ_{j ≥ i} los[j]`, so `los[j]` feeds every survival
        ## weight at or below `j`: the adjoint is the forward cumulative
        ## sum of the survival adjoint.
        l̄ = zeros(eltype(s̄), L)
        run = zero(eltype(s̄))
        @inbounds for j in 1:L
            run += s̄[j]
            l̄[j] = run
        end
        return ChainRulesCore.NoTangent(), x̄, l̄
    end
    return y, convolve_survival_pullback
end

function ChainRulesCore.rrule(
        ::typeof(convolve_pmf), a::AbstractVector,
        b::AbstractVector
    )
    y = convolve_pmf(a, b)
    function convolve_pmf_pullback(Δy)
        ȳ = ChainRulesCore.unthunk(Δy)
        na = length(a)
        nb = length(b)
        ā = zeros(promote_type(eltype(ȳ), eltype(b)), na)
        b̄ = zeros(promote_type(eltype(ȳ), eltype(a)), nb)
        @inbounds for i in 1:na, j in 1:nb

            g = ȳ[i + j - 1]
            ā[i] += g * b[j]
            b̄[j] += g * a[i]
        end
        return ChainRulesCore.NoTangent(), ā, b̄
    end
    return y, convolve_pmf_pullback
end

function ChainRulesCore.rrule(
        ::typeof(interpolate_knots), knot_vals::AbstractVector,
        days::AbstractVector{<:Integer}, n::Integer
    )
    out = interpolate_knots(knot_vals, days, n)
    function interpolate_knots_pullback(Δout)
        ō = ChainRulesCore.unthunk(Δout)
        nb = length(days)
        k̄ = zeros(eltype(ō), nb)
        if nb == 1
            @inbounds for t in 1:n
                k̄[1] += ō[t]
            end
            return (
                ChainRulesCore.NoTangent(), k̄,
                ChainRulesCore.NoTangent(), ChainRulesCore.NoTangent(),
            )
        end
        Tf = eltype(ō)
        @inbounds for t in 1:n
            b = 1
            while b < nb - 1 && t > days[b + 1]
                b += 1
            end
            d0 = days[b]
            d1 = days[b + 1]
            frac = d1 == d0 ? zero(Tf) :
                clamp(Tf(t - d0) / Tf(d1 - d0), zero(Tf), one(Tf))
            g = ō[t]
            k̄[b] += (one(Tf) - frac) * g
            k̄[b + 1] += frac * g
        end
        return (
            ChainRulesCore.NoTangent(), k̄,
            ChainRulesCore.NoTangent(), ChainRulesCore.NoTangent(),
        )
    end
    return out, interpolate_knots_pullback
end

function ChainRulesCore.rrule(
        ::typeof(renewal_infections), Rt::AbstractVector,
        g::AbstractVector, seed::AbstractVector
    )
    n = length(Rt)
    L = length(seed)
    Tp = promote_type(eltype(Rt), eltype(g), eltype(seed))
    I = zeros(Tp, n)
    force = zeros(Tp, n)
    @inbounds for j in 1:min(L, n)
        I[j] = seed[j]
    end
    @inbounds for t in (L + 1):n
        f = zero(Tp)
        kmax = min(t - 1, length(g))
        for s in 1:kmax
            f += I[t - s] * g[s]
        end
        force[t] = f
        I[t] = Rt[t] * f
    end
    function renewal_infections_pullback(ΔI)
        ## The recursion is sequential, so the reverse pass walks the days
        ## backwards, pushing each day's infection adjoint onto the lagged
        ## infections it was built from before those days are themselves read.
        Ī = collect(float.(ChainRulesCore.unthunk(ΔI)))
        Tb = eltype(Ī)
        R̄ = zeros(Tb, n)
        ḡ = zeros(Tb, length(g))
        s̄ = zeros(Tb, L)
        @inbounds for t in n:-1:(L + 1)
            it = Ī[t]
            R̄[t] += it * force[t]
            f̄ = it * Rt[t]
            kmax = min(t - 1, length(g))
            for s in 1:kmax
                Ī[t - s] += f̄ * g[s]
                ḡ[s] += f̄ * I[t - s]
            end
        end
        @inbounds for j in 1:min(L, n)
            s̄[j] = Ī[j]
        end
        return ChainRulesCore.NoTangent(), R̄, ḡ, s̄
    end
    return I, renewal_infections_pullback
end

## --- Mooncake registration ----------------------------------------------
##
## Mooncake is the package default, so the rules above are wired into it
## here rather than in an extension. The signatures are restricted to
## `Array{<:IEEEFloat}` arguments: that is what every model call site
## passes, and it keeps the tangent types Mooncake builds for the rule
## concrete. Anything else falls through to Mooncake's own derived rule,
## which still differentiates the plain Julia body.
Mooncake.@from_rrule(
    Mooncake.DefaultCtx,
    Tuple{
        typeof(convolve_delay), Array{<:Mooncake.IEEEFloat},
        Array{<:Mooncake.IEEEFloat},
    },
)
Mooncake.@from_rrule(
    Mooncake.DefaultCtx,
    Tuple{
        typeof(convolve_survival), Array{<:Mooncake.IEEEFloat},
        Array{<:Mooncake.IEEEFloat},
    },
)
Mooncake.@from_rrule(
    Mooncake.DefaultCtx,
    Tuple{
        typeof(convolve_pmf), Array{<:Mooncake.IEEEFloat},
        Array{<:Mooncake.IEEEFloat},
    },
)
Mooncake.@from_rrule(
    Mooncake.DefaultCtx,
    Tuple{
        typeof(interpolate_knots), Array{<:Mooncake.IEEEFloat},
        Array{<:Integer}, Integer,
    },
)
Mooncake.@from_rrule(
    Mooncake.DefaultCtx,
    Tuple{
        typeof(renewal_infections), Array{<:Mooncake.IEEEFloat},
        Array{<:Mooncake.IEEEFloat}, Array{<:Mooncake.IEEEFloat},
    },
)
