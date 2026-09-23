# Native `Mooncake.rrule!!` methods for the `renewal.jl` kernels and the
# observation kernels in `models/observations.jl`. Each is a loop over a
# daily grid or an observation vector, so left to the backend every
# iteration's intermediates go on the tape; these replace that with a
# closed-form adjoint of the same shape as the forward loop.
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
    Tuple{typeof(convolve_delay), _FloatVec, _FloatVec},
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
        typeof(interpolate_knots), _FloatVec, Array{<:Integer}, Integer,
    },
)
Mooncake.@is_primitive(
    Mooncake.MinimalCtx,
    Tuple{
        typeof(renewal_infections), Array{<:Mooncake.IEEEFloat},
        Array{<:Mooncake.IEEEFloat}, Array{<:Mooncake.IEEEFloat},
    },
)
Mooncake.@is_primitive(
    Mooncake.MinimalCtx,
    Tuple{
        typeof(patch_infections), Matrix{<:Mooncake.IEEEFloat},
        Vector{<:Mooncake.IEEEFloat}, Matrix{<:Mooncake.IEEEFloat},
        Matrix{<:Mooncake.IEEEFloat}, Matrix{<:Mooncake.IEEEFloat},
    },
)
Mooncake.@is_primitive(
    Mooncake.MinimalCtx,
    Tuple{
        typeof(nbinomial_loglik), Mooncake.IEEEFloat,
        Array{<:Mooncake.IEEEFloat}, Array{<:Integer},
    },
)

## Adjoint of the daily delay convolution `convolve_delay`. A convolution's
## pullback is the matching correlation, one pass over the same `(t, d)`
## pairs:
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
        x::CoDual{<:_FloatVec}, delay::CoDual{<:_FloatVec}
    )
    xp = primal(x)
    dp = primal(delay)
    x̄ = _tangent_array(x)
    d̄ = _tangent_array(delay)
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

Mooncake.@is_primitive(
    Mooncake.MinimalCtx,
    Tuple{
        typeof(abscond_thinned), Array{<:Mooncake.IEEEFloat},
        Mooncake.IEEEFloat,
    },
)
Mooncake.@is_primitive(
    Mooncake.MinimalCtx,
    Tuple{
        typeof(abscond_thinned_flows), Array{<:Mooncake.IEEEFloat},
        Array{<:Mooncake.IEEEFloat}, Array{<:Mooncake.IEEEFloat},
        Array{<:Mooncake.IEEEFloat}, Mooncake.IEEEFloat,
        Array{<:Mooncake.IEEEFloat},
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

## Adjoint of the confirmation-stopped abscond thinning. Per cohort `t` the
## forward walk is
##
##     u[d] = u[d−1] · (1 − h[t+d−1]),    s[d] = s[d−1] · (1 − κ u[d])
##     outk[t+d] += admk[t] · pmfk[d+1] · s[d]
##
## The pullback recomputes `u` and `s` for one cohort into scratch buffers,
## then walks the cohort backwards carrying `s̄` and `ū`. It walks the
## cohorts the forward pass skips (both admissions zero) too: their output
## is zero but its derivative in the admissions is not.
function Mooncake.rrule!!(
        ::CoDual{typeof(abscond_thinned_flows)},
        adm1::CoDual{<:Array{<:Mooncake.IEEEFloat}},
        pmf1::CoDual{<:Array{<:Mooncake.IEEEFloat}},
        adm2::CoDual{<:Array{<:Mooncake.IEEEFloat}},
        pmf2::CoDual{<:Array{<:Mooncake.IEEEFloat}},
        κ::CoDual{<:Mooncake.IEEEFloat},
        conf_hazard::CoDual{<:Array{<:Mooncake.IEEEFloat}}
    )
    a1p, p1p = primal(adm1), primal(pmf1)
    a2p, p2p = primal(adm2), primal(pmf2)
    κp = primal(κ)
    hp = primal(conf_hazard)
    ## `abscond_thinned_flow` passes the same arrays twice, so these tangents
    ## can alias: the pullback only ever adds into them.
    ā1, p̄1 = tangent(adm1), tangent(pmf1)
    ā2, p̄2 = tangent(adm2), tangent(pmf2)
    h̄ = tangent(conf_hazard)
    y1, y2 = abscond_thinned_flows(a1p, p1p, a2p, p2p, κp, hp)
    ȳ1 = zero(y1)
    ȳ2 = zero(y2)
    function abscond_thinned_flows_pullback!!(::NoRData)
        n = length(a1p)
        nmax1 = length(p1p)
        nmax2 = length(p2p)
        T = eltype(y1)
        one_T = one(T)
        u = Vector{T}(undef, max(nmax1, nmax2))
        s = similar(u)
        κ̄ = zero(T)
        @inbounds for t in 1:n
            a1 = a1p[t]
            a2 = a2p[t]
            dmax1 = min(nmax1 - 1, n - t)
            dmax2 = min(nmax2 - 1, n - t)
            dm = max(dmax1, dmax2)
            ## `u[d+1]`, `s[d+1]` hold the survivals at cohort age `d`.
            ā1[t] += ȳ1[t] * p1p[1]
            ā2[t] += ȳ2[t] * p2p[1]
            p̄1[1] += ȳ1[t] * a1
            p̄2[1] += ȳ2[t] * a2
            u[1] = one_T
            s[1] = one_T
            for d in 1:dm
                u[d + 1] = u[d] * (one_T - hp[t + d - 1])
                s[d + 1] = s[d] * (one_T - κp * u[d + 1])
            end
            ūc = zero(T)
            s̄c = zero(T)
            for d in dm:-1:1
                sd = s[d + 1]
                if d <= dmax1
                    g = ȳ1[t + d]
                    ā1[t] += g * p1p[d + 1] * sd
                    p̄1[d + 1] += g * a1 * sd
                    s̄c += g * a1 * p1p[d + 1]
                end
                if d <= dmax2
                    g = ȳ2[t + d]
                    ā2[t] += g * p2p[d + 1] * sd
                    p̄2[d + 1] += g * a2 * sd
                    s̄c += g * a2 * p2p[d + 1]
                end
                ## s[d] = s[d−1] · (1 − κ u[d])
                ud = u[d + 1]
                sprev = s[d]
                κ̄ -= s̄c * sprev * ud
                ūc -= s̄c * sprev * κp
                s̄c *= one_T - κp * ud
                ## u[d] = u[d−1] · (1 − h[t+d−1])
                hd = hp[t + d - 1]
                h̄[t + d - 1] -= ūc * u[d]
                ūc *= one_T - hd
            end
        end
        return NoRData(), NoRData(), NoRData(), NoRData(), NoRData(),
            oftype(κp, κ̄), NoRData()
    end
    return CoDual((y1, y2), (ȳ1, ȳ2)), abscond_thinned_flows_pullback!!
end

Mooncake.@is_primitive(
    Mooncake.MinimalCtx,
    Tuple{
        typeof(two_clock_confirmed), Array{<:Mooncake.IEEEFloat},
        Array{<:Mooncake.IEEEFloat}, Array{<:Mooncake.IEEEFloat},
    },
)
Mooncake.@is_primitive(
    Mooncake.MinimalCtx,
    Tuple{
        typeof(clinical_stay_survival), Array{<:Mooncake.IEEEFloat},
        Array{<:Mooncake.IEEEFloat}, Mooncake.IEEEFloat,
    },
)

function Mooncake.rrule!!(
        ::CoDual{typeof(two_clock_confirmed)},
        A_bvd::CoDual{<:Array{<:Mooncake.IEEEFloat}},
        conf_hazard::CoDual{<:Array{<:Mooncake.IEEEFloat}},
        S_clin::CoDual{<:Array{<:Mooncake.IEEEFloat}}
    )
    Ap = primal(A_bvd)
    hp = primal(conf_hazard)
    Sp = primal(S_clin)
    Ā = tangent(A_bvd)
    h̄ = tangent(conf_hazard)
    S̄ = tangent(S_clin)
    O = two_clock_confirmed(Ap, hp, Sp)
    Ō = zero(O)
    function two_clock_confirmed_pullback!!(::NoRData)
        n = length(Ap)
        L = length(Sp)
        one_T = one(eltype(O))
        ## Day `t` walks cohorts `u = t, t-1, …` with the unconfirmed
        ## product `p` extended by `(1 − h[u])` after each is scored:
        ##
        ##     O[t] += A[u] · (1 − p) · S[t−u+1],   p ← p · (1 − h[u])
        ##
        ## The pullback replays each day's products into `p_at`, then walks
        ## the cohorts back in the opposite order carrying the adjoint of
        ## `p`, so no factor is divided out.
        p_at = Vector{eltype(O)}(undef, L)
        @inbounds for t in 1:n
            g = Ō[t]
            iszero(g) && continue
            umin = max(1, t - L + 1)
            p = one_T
            for u in t:-1:umin
                p_at[t - u + 1] = p
                p *= (one_T - hp[u])
            end
            p̄ = zero(one_T)
            for u in umin:t
                d = t - u
                pk = p_at[d + 1]
                s = Sp[d + 1]
                h̄[u] -= p̄ * pk
                Ā[u] += g * (one_T - pk) * s
                S̄[d + 1] += g * Ap[u] * (one_T - pk)
                p̄ = p̄ * (one_T - hp[u]) - g * Ap[u] * s
            end
        end
        return NoRData(), NoRData(), NoRData(), NoRData()
    end
    return CoDual(O, Ō), two_clock_confirmed_pullback!!
end

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
        typeof(accumulate_occupancy), Array{<:Mooncake.IEEEFloat},
        Array{<:Mooncake.IEEEFloat}, Array{<:Mooncake.IEEEFloat},
        Array{<:Mooncake.IEEEFloat}, Array{<:Mooncake.IEEEFloat},
        Mooncake.IEEEFloat, Array{<:Mooncake.IEEEFloat},
    },
)

## Branch flags of the occupancy balance, one bit per `max`/`clamp` side.
## Ties go to the second argument of `max` and to the bound of `clamp`,
## lower bound first, the sides Mooncake's own rules for them take.
const _OCC_UNCONF = 0x01
const _OCC_DENOM = 0x02
const _OCC_BVD = 0x04
const _OCC_BG = 0x08
const _OCC_CONF_X = 0x10
const _OCC_CONF_HI = 0x20
const _OCC_SUSP = 0x40

## The forward balance of `accumulate_occupancy`, also recording the
## non-case stock and which side of each `max` and `clamp` was taken.
function _accumulate_occupancy_taped(
        A_bvd, A_bg, deaths, recover, ruleout, κ, conf_hazard
    )
    n = length(A_bvd)
    T = promote_type(
        eltype(A_bvd), eltype(A_bg), eltype(deaths),
        eltype(recover), eltype(ruleout), typeof(κ), eltype(conf_hazard)
    )
    demand = Vector{T}(undef, n)
    O_bvd = Vector{T}(undef, n)
    O_conf = Vector{T}(undef, n)
    O_susp = Vector{T}(undef, n)
    abscond = Vector{T}(undef, n)
    O_bg = Vector{T}(undef, n)
    flags = Vector{UInt8}(undef, n)
    z = zero(T)
    Obvd_prev = z
    Obg_prev = z
    Oconf_prev = z
    Osusp_prev = z
    ε = eps(T)
    @inbounds for t in 1:n
        bvd_out = deaths[t] + recover[t]
        ab = κ * Osusp_prev
        x_u = Obvd_prev - Oconf_prev
        unconf = max(x_u, z)
        denom = max(Osusp_prev, ε)
        ab_bvd = ab * (unconf / denom)
        ab_bg = ab * (Obg_prev / denom)
        x_bvd = Obvd_prev + A_bvd[t] - bvd_out - ab_bvd
        Obvd_t = max(x_bvd, z)
        x_bg = Obg_prev + A_bg[t] - ruleout[t] - ab_bg
        Obg_t = max(x_bg, z)
        Dt = Obvd_t + Obg_t
        conf_in = conf_hazard[t] * unconf
        share = Obvd_prev > z ? Oconf_prev / Obvd_prev : z
        x_conf = Oconf_prev + conf_in - bvd_out * share
        Oconf_t = clamp(x_conf, z, Obvd_t)
        x_susp = Dt - Oconf_t
        Osusp_t = max(x_susp, z)
        f = 0x00
        x_u > z && (f |= _OCC_UNCONF)
        Osusp_prev > ε && (f |= _OCC_DENOM)
        x_bvd > z && (f |= _OCC_BVD)
        x_bg > z && (f |= _OCC_BG)
        if x_conf > z
            f |= x_conf < Obvd_t ? _OCC_CONF_X : _OCC_CONF_HI
        end
        x_susp > z && (f |= _OCC_SUSP)
        demand[t] = Dt
        O_bvd[t] = Obvd_t
        O_conf[t] = Oconf_t
        O_susp[t] = Osusp_t
        abscond[t] = ab
        O_bg[t] = Obg_t
        flags[t] = f
        Obvd_prev = Obvd_t
        Obg_prev = Obg_t
        Oconf_prev = Oconf_t
        Osusp_prev = Osusp_t
    end
    return (; demand, O_bvd, O_conf, O_susp, abscond), O_bg, flags
end

function Mooncake.rrule!!(
        ::CoDual{typeof(accumulate_occupancy)},
        A_bvd::CoDual{<:Array{<:Mooncake.IEEEFloat}},
        A_bg::CoDual{<:Array{<:Mooncake.IEEEFloat}},
        deaths::CoDual{<:Array{<:Mooncake.IEEEFloat}},
        recover::CoDual{<:Array{<:Mooncake.IEEEFloat}},
        ruleout::CoDual{<:Array{<:Mooncake.IEEEFloat}},
        κ::CoDual{<:Mooncake.IEEEFloat},
        conf_hazard::CoDual{<:Array{<:Mooncake.IEEEFloat}}
    )
    κp = primal(κ)
    hp = primal(conf_hazard)
    dp = primal(deaths)
    rp = primal(recover)
    Āb = tangent(A_bvd)
    Āg = tangent(A_bg)
    d̄ = tangent(deaths)
    r̄ = tangent(recover)
    ō = tangent(ruleout)
    h̄ = tangent(conf_hazard)
    y, O_bg, flags = _accumulate_occupancy_taped(
        primal(A_bvd), primal(A_bg), dp, rp, primal(ruleout), κp, hp
    )
    ȳ = map(zero, y)
    function accumulate_occupancy_pullback!!(::NoRData)
        Tf = eltype(y.demand)
        z = zero(Tf)
        ε = eps(Tf)
        κ̄ = zero(Tf)
        ## Adjoints of the day-`t` stocks, carried back from day `t + 1`.
        cb = z
        cg = z
        cc = z
        cs = z
        @inbounds for t in length(flags):-1:1
            f = flags[t]
            Pb = t > 1 ? y.O_bvd[t - 1] : z
            Pg = t > 1 ? O_bg[t - 1] : z
            Pc = t > 1 ? y.O_conf[t - 1] : z
            Ps = t > 1 ? y.O_susp[t - 1] : z
            bvd_out = dp[t] + rp[t]
            ab = κp * Ps
            unconf = f & _OCC_UNCONF != 0 ? Pb - Pc : z
            denom = f & _OCC_DENOM != 0 ? Ps : ε
            share = Pb > z ? Pc / Pb : z
            gOb = ȳ.O_bvd[t] + cb
            gOg = cg
            gOc = ȳ.O_conf[t] + cc
            gOs = ȳ.O_susp[t] + cs
            gD = ȳ.demand[t]
            gab = ȳ.abscond[t]
            ## `O_susp = max(D − O_conf, 0)`
            if f & _OCC_SUSP != 0
                gD += gOs
                gOc -= gOs
            end
            ## `O_conf = clamp(x_conf, 0, O_bvd)`
            gx_c = z
            if f & _OCC_CONF_X != 0
                gx_c = gOc
            elseif f & _OCC_CONF_HI != 0
                gOb += gOc
            end
            nb = z
            nc = gx_c
            gbo = -gx_c * share
            if Pb > z
                gshare = -gx_c * bvd_out
                nc += gshare / Pb
                nb -= gshare * share / Pb
            end
            h̄[t] += gx_c * unconf
            gu = gx_c * hp[t]
            ## `D = O_bvd + O_bg`, each `max(x, 0)` of its balance
            gOb += gD
            gOg += gD
            gx_b = f & _OCC_BVD != 0 ? gOb : z
            gx_g = f & _OCC_BG != 0 ? gOg : z
            nb += gx_b
            Āb[t] += gx_b
            gbo -= gx_b
            ng = gx_g
            Āg[t] += gx_g
            ō[t] -= gx_g
            ## Abscond split `ab · (unconf / denom)` and `ab · (O_bg / denom)`
            qu = unconf / denom
            qg = Pg / denom
            gab -= gx_b * qu + gx_g * qg
            gqu = -gx_b * ab
            gqg = -gx_g * ab
            gu += gqu / denom
            ng += gqg / denom
            gden = -(gqu * qu + gqg * qg) / denom
            ns = f & _OCC_DENOM != 0 ? gden : z
            if f & _OCC_UNCONF != 0
                nb += gu
                nc -= gu
            end
            κ̄ += gab * Ps
            ns += gab * κp
            d̄[t] += gbo
            r̄[t] += gbo
            cb = nb
            cg = ng
            cc = nc
            cs = ns
        end
        return NoRData(), NoRData(), NoRData(), NoRData(), NoRData(),
            NoRData(), convert(typeof(κp), κ̄), NoRData()
    end
    return CoDual(y, ȳ), accumulate_occupancy_pullback!!
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
Mooncake.@is_primitive(
    Mooncake.MinimalCtx,
    Tuple{
        typeof(onset_report_expected_total), Array{<:Mooncake.IEEEFloat},
        Array{<:Mooncake.IEEEFloat}, Array{<:Mooncake.IEEEFloat}, Integer,
        Array{<:Mooncake.IEEEFloat}, Integer,
    },
)

## Derivative of `safe_rate`: the identity above its floor, flat on it.
_safe_rate_slope(x) = (isfinite(x) && x > eps(typeof(x))) ? one(x) : zero(x)

## Hazards and survival products of one onset date's delay column, as the
## forward loop builds them. `surv[j + 1]` is the product up to and
## including delay `j`.
function _onset_column!(
        h::AbstractVector, surv::AbstractVector,
        logit_h0::AbstractVector, γ::AbstractVector, u::Integer,
        grid_start::Integer
    )
    D = length(logit_h0)
    ng = length(γ)
    T = eltype(surv)
    s = one(T)
    @inbounds for j in 0:(D - 1)
        gi = clamp(u + j - grid_start + 1, 1, ng)
        hj = logistic(logit_h0[j + 1] + γ[gi])
        s *= (one(T) - hj)
        h[j + 1] = hj
        surv[j + 1] = s
    end
    return nothing
end

## `1 - surv[j]` has derivative `surv[j] · h[i]` in the logit `x_i` for every
## `i ≤ j`, so a column's adjoint is `x̄_i = h_i · Σ_{j ≥ i} c̄_j surv_j`, one
## backward running sum. Each `x̄_i` lands on `logit_h0[i]` and on the
## (clamped) calendar day it read from `γ`.
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
    @inbounds for k in 1:nu
        _onset_column!(view(H, :, k), view(S, :, k), lp, γp, lo + k - 1, gs)
    end
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
            jc = (δc < 0 || D == 0) ? 0 : min(Int(δc), D - 1) + 1
            jp = (δp < 0 || D == 0) ? 0 : min(Int(δp), D - 1) + 1
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

function Mooncake.rrule!!(
        ::CoDual{typeof(onset_report_expected_total)},
        onsets::CoDual{<:Array{<:Mooncake.IEEEFloat}},
        logit_h0::CoDual{<:Array{<:Mooncake.IEEEFloat}},
        γ::CoDual{<:Array{<:Mooncake.IEEEFloat}},
        grid_start::CoDual{<:Integer},
        alpha::CoDual{<:Array{<:Mooncake.IEEEFloat}},
        as_of::CoDual{<:Integer}
    )
    op = primal(onsets)
    lp = primal(logit_h0)
    γp = primal(γ)
    gs = Int(primal(grid_start))
    alp = primal(alpha)
    t = Int(primal(as_of))
    ō = tangent(onsets)
    l̄ = tangent(logit_h0)
    γ̄ = tangent(γ)
    ᾱ = tangent(alpha)
    ## Each term is `onsets[u] · α_u · G_u` with `G = (1 - surv_jn) / den`
    ## and `den = safe_rate(1 - surv_{D-1})`: the numerator's column adjoint
    ## is the cotangent at `jn`, the denominator's at `D - 1`. The forward
    ## pass keeps each onset date's column for the pullback.
    D = length(lp)
    n = length(op)
    na = length(alp)
    T = promote_type(eltype(op), eltype(lp), eltype(γp), eltype(alp))
    ge = min(t, n)
    H = Matrix{T}(undef, D, max(ge, 0))
    S = Matrix{T}(undef, D, max(ge, 0))
    total = zero(T)
    @inbounds for u in 1:ge
        _onset_column!(view(H, :, u), view(S, :, u), lp, γp, u, gs)
        δ = t - u
        α = alp[clamp(u - gs + 1, 1, na)]
        jn = min(δ, D - 1)
        num = (δ < 0 || D == 0) ? zero(T) : one(T) - S[jn + 1, u]
        den = one(T) - (D == 0 ? one(T) : S[D, u])
        total += op[u] * (α * (num / safe_rate(den)))
    end
    function onset_report_expected_total_pullback!!(ȳ::Mooncake.IEEEFloat)
        D == 0 && return ntuple(_ -> NoRData(), 7)
        c̄ = zeros(T, D)
        @inbounds for u in 1:ge
            δ = t - u
            ia = clamp(u - gs + 1, 1, na)
            α = alp[ia]
            jn = min(δ, D - 1)
            num = δ < 0 ? zero(T) : one(T) - S[jn + 1, u]
            cD = one(T) - S[D, u]
            sden = safe_rate(cD)
            G = num / sden
            ō[u] += ȳ * α * G
            ᾱ[ia] += ȳ * op[u] * G
            w = ȳ * op[u] * α / sden
            fill!(c̄, zero(T))
            δ >= 0 && (c̄[jn + 1] += w)
            c̄[D] -= w * G * _safe_rate_slope(cD)
            _onset_column_adjoint!(
                l̄, γ̄, c̄, view(H, :, u), view(S, :, u), u, gs
            )
        end
        return ntuple(_ -> NoRData(), 7)
    end
    return CoDual(total, NoFData()), onset_report_expected_total_pullback!!
end

function Mooncake.rrule!!(
        ::CoDual{typeof(patch_infections)},
        Rt::CoDual{<:Matrix{<:Mooncake.IEEEFloat}},
        g::CoDual{<:Vector{<:Mooncake.IEEEFloat}},
        seeds::CoDual{<:Matrix{<:Mooncake.IEEEFloat}},
        K::CoDual{<:Matrix{<:Mooncake.IEEEFloat}},
        ε::CoDual{<:Matrix{<:Mooncake.IEEEFloat}}
    )
    Rp = primal(Rt)
    gp = primal(g)
    Kp = primal(K)
    εp = primal(ε)
    R̄ = tangent(Rt)
    ḡ = tangent(g)
    s̄ = tangent(seeds)
    K̄ = tangent(K)
    ε̄ = tangent(ε)
    np, n = size(Rp)
    L = size(primal(seeds), 2)
    y = patch_infections(Rp, gp, primal(seeds), Kp, εp)
    out = Mooncake.zero_fcodual(y)
    Ī = tangent(out).infections
    Ā = tangent(out).importation
    I = y.infections
    function patch_infections_pullback!!(::NoRData)
        Tf = eltype(I)
        outflow = zeros(Tf, np)
        @inbounds for q in 1:np, r in 1:np
            r == q && continue
            outflow[q] += Kp[r, q]
        end
        ## Each day's force and generated infections are rebuilt from the
        ## forward infections rather than stored. The walk runs backwards so
        ## a day's adjoint is complete before it is pushed onto earlier days,
        ## and it accumulates into a copy since `Ī` is Mooncake's buffer.
        acc = copy(Ī)
        force = zeros(Tf, np)
        gen = zeros(Tf, np)
        ḡen = zeros(Tf, np)
        ōut = zeros(Tf, np)
        @inbounds for t in n:-1:(L + 1)
            kmax = min(t - 1, length(gp))
            for p in 1:np
                f = zero(Tf)
                for s in 1:kmax
                    f += I[p, t - s] * gp[s]
                end
                force[p] = f
                gen[p] = Rp[p, t] * f
                ḡen[p] = zero(Tf)
            end
            ## I[p, t] = (1 - ε[p, t] outflow[p]) gen[p] + arrivals[p], and
            ## the importation output is arrivals[p] alone.
            for p in 1:np
                a = acc[p, t]
                ā = a + Ā[p, t]
                ḡen[p] += a * (one(Tf) - εp[p, t] * outflow[p])
                ε̄[p, t] -= a * outflow[p] * gen[p]
                ōut[p] -= a * εp[p, t] * gen[p]
                for q in 1:np
                    q == p && continue
                    ε̄[q, t] += ā * Kp[p, q] * gen[q]
                    K̄[p, q] += ā * εp[q, t] * gen[q]
                    ḡen[q] += ā * εp[q, t] * Kp[p, q]
                end
            end
            for p in 1:np
                R̄[p, t] += ḡen[p] * force[p]
                f̄ = ḡen[p] * Rp[p, t]
                for s in 1:kmax
                    acc[p, t - s] += f̄ * gp[s]
                    ḡ[s] += f̄ * I[p, t - s]
                end
            end
        end
        @inbounds for q in 1:np, r in 1:np
            r == q && continue
            K̄[r, q] += ōut[q]
        end
        @inbounds for p in 1:np, j in 1:min(L, n)
            s̄[p, j] += acc[p, j]
        end
        return ntuple(_ -> NoRData(), 6)
    end
    return out, patch_infections_pullback!!
end

## Value and gradient of `nbinomial_loglik` in one pass. With `r = k` and
## `p = r / (r + μ)` each count `x` contributes
##
##     ∂ℓ/∂r = log p + ψ(r + x) − ψ(r),    ∂ℓ/∂p = r / p − x / (1 − p),
##
## chained through `p` to `μ` and `k`. The guards in `safe_rate` and
## `safe_nbinomial` are mirrored: a floored `k` or `μ`, or a clamped `p`,
## passes no derivative. A term that is not finite adds none either.
function _nbinomial_loglik_grad(
        k::T, modelled::AbstractVector,
        obs::AbstractVector
    ) where {T}
    lo = eps(T)
    hi = one(T) - lo
    k_on = isfinite(k) && k > zero(k)
    r = k_on ? k : lo
    s = zero(T)
    dk = zero(T)
    dμ = zeros(T, length(modelled))
    @inbounds for i in eachindex(modelled, obs)
        μ = modelled[i]
        x = obs[i]
        μ_on = isfinite(μ) && μ > lo
        m = μ_on ? μ : lo
        p_raw = r / (r + m)
        p = isfinite(p_raw) ? clamp(p_raw, lo, hi) : lo
        ℓ = logpdf(NegativeBinomial(r, p), x)
        s += ℓ
        isfinite(ℓ) || continue
        ∂r = iszero(x) ? log(p) : log(p) + digamma(r + x) - digamma(r)
        if isfinite(p_raw) && !(p_raw > hi) && !(p_raw < lo)
            ∂p = r / p - x / (one(T) - p)
            den = (r + m)^2
            ∂r += ∂p * m / den
            μ_on && (dμ[i] = -∂p * r / den)
        end
        dk += ∂r
    end
    return s, (k_on ? dk : zero(T)), dμ
end

function Mooncake.rrule!!(
        ::CoDual{typeof(nbinomial_loglik)},
        k::CoDual{<:Mooncake.IEEEFloat},
        modelled::CoDual{<:Array{<:Mooncake.IEEEFloat}},
        obs::CoDual{<:Array{<:Integer}}
    )
    s, dk, dμ = _nbinomial_loglik_grad(
        primal(k), primal(modelled), primal(obs)
    )
    μ̄ = tangent(modelled)
    ## `k` is a scalar, so its adjoint goes back as rdata rather than into
    ## a tangent buffer.
    function nbinomial_loglik_pullback!!(s̄)
        μ̄ .+= s̄ .* dμ
        return NoRData(), s̄ * dk, NoRData(), NoRData()
    end
    return CoDual(s, NoFData()), nbinomial_loglik_pullback!!
end

## The increments `studentt_loglik` scores: counts in the package data, floats
## in simulated triangles. Float increments take the adjoint `−∂ℓ/∂μ`.
const _StudentTObs = Union{Array{<:Integer}, _FloatArray}

Mooncake.@is_primitive(
    Mooncake.MinimalCtx,
    Tuple{
        typeof(studentt_loglik), Array{<:Mooncake.IEEEFloat},
        Array{<:Mooncake.IEEEFloat}, _StudentTObs, Mooncake.IEEEFloat,
    },
)

## Value and gradient of `studentt_loglik` in one pass. With
## `z = (x − μ) / σ` and `g = (ν + 1) z / (ν + z²)` each cell contributes
##
##     ∂ℓ/∂μ = g / σ,    ∂ℓ/∂σ = (g z − 1) / σ,
##     ∂ℓ/∂ν = (ψ((ν + 1)/2) − ψ(ν/2) − 1/ν − log1p(z²/ν) + g z / ν) / 2,
##
## and `∂ℓ/∂x = −∂ℓ/∂μ`. The guards in `safe_studentt` are mirrored: a
## floored `σ` or a defaulted `ν` passes no derivative. A cell whose term is
## not finite adds none either.
function _studentt_loglik_grad(
        means::AbstractVector, sds::AbstractVector,
        obs::AbstractVector, ν::Real
    )
    T = float(promote_type(eltype(means), eltype(sds), typeof(ν)))
    ν_on = isfinite(ν) && ν > zero(ν)
    νc = ν_on ? ν : oftype(float(ν), 4)
    νp12 = (νc + 1) / 2
    c = logpdf(TDist(νc), zero(T))
    ∂ν_c = (digamma(νp12) - digamma(νc / 2) - 1 / νc) / 2
    s = zero(T)
    dν = zero(T)
    dμ = zeros(T, length(means))
    dσ = zeros(T, length(means))
    @inbounds for i in eachindex(means, sds, obs)
        σ = sds[i]
        σ_on = isfinite(σ) && σ > zero(σ)
        σc = σ_on ? σ : eps(typeof(float(σ)))
        z = (obs[i] - means[i]) / σc
        l1 = log1p(z^2 / νc)
        ℓ = (c - νp12 * l1) - log(σc)
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
        ::CoDual{typeof(studentt_loglik)},
        means::CoDual{<:Array{<:Mooncake.IEEEFloat}},
        sds::CoDual{<:Array{<:Mooncake.IEEEFloat}},
        obs::CoDual{<:_StudentTObs},
        ν::CoDual{<:Mooncake.IEEEFloat}
    )
    s, dν, dμ, dσ = _studentt_loglik_grad(
        primal(means), primal(sds), primal(obs), primal(ν)
    )
    μ̄ = tangent(means)
    σ̄ = tangent(sds)
    x̄ = tangent(obs)
    ## `ν` is a scalar, so its adjoint goes back as rdata. Integer
    ## increments carry no tangent.
    function studentt_loglik_pullback!!(s̄)
        μ̄ .+= s̄ .* dμ
        σ̄ .+= s̄ .* dσ
        x̄ isa _FloatArray && (x̄ .-= s̄ .* dμ)
        return NoRData(), NoRData(), NoRData(), NoRData(), s̄ * dν
    end
    return CoDual(s, NoFData()), studentt_loglik_pullback!!
end

Mooncake.@is_primitive(
    Mooncake.MinimalCtx,
    Tuple{
        typeof(betabinomial_loglik), AbstractVector{<:Integer},
        Array{<:Mooncake.IEEEFloat}, Mooncake.IEEEFloat,
        AbstractVector{<:Integer},
    },
)

## Value and gradient of `betabinomial_loglik` in one pass. With
## concentration `c = (1 − ρ) / ρ`, `α = c·p` and `β = c·(1 − p)`, each
## count `x` of `n` trials contributes
##
##     ∂ℓ/∂α = ψ(x + α) − ψ(α) + ψ(α + β) − ψ(n + α + β),
##     ∂ℓ/∂β = ψ(n − x + β) − ψ(β) + ψ(α + β) − ψ(n + α + β),
##
## chained through `α` and `β` to `p` and `c`, and through `c` to `ρ`. The
## guards in `safe_betabinomial` are mirrored: a clamped `p` or `ρ`, or an
## `α` or `β` held at its floor, passes no derivative. A term that is not
## finite adds none either.
function _betabinomial_loglik_grad(
        trials::AbstractVector, p::AbstractVector, ρ, obs::AbstractVector
    )
    T = float(promote_type(eltype(p), typeof(ρ)))
    lo = eps(T)
    ρ_lo = T(1.0e-6)
    ρ_on = isfinite(ρ) && !(ρ < ρ_lo) && !(ρ > one(T) - ρ_lo)
    ρc = isfinite(ρ) ? clamp(T(ρ), ρ_lo, one(T) - ρ_lo) : ρ_lo
    c = (one(T) - ρc) / ρc
    s = zero(T)
    dc = zero(T)
    dp = zeros(T, length(p))
    @inbounds for i in eachindex(trials, p, obs)
        n = trials[i]
        x = obs[i]
        q = p[i]
        d = safe_betabinomial(n, q, ρ)
        ℓ = logpdf(d, x)
        s += ℓ
        isfinite(ℓ) || continue
        α, β = d.α, d.β
        pc = isfinite(q) ? clamp(T(q), lo, one(T) - lo) : one(T) / 2
        ψ_ab = digamma(α + β) - digamma(n + α + β)
        ∂α = c * pc > lo ? digamma(x + α) - digamma(α) + ψ_ab : zero(T)
        ∂β = c * (one(T) - pc) > lo ?
            digamma(n - x + β) - digamma(β) + ψ_ab : zero(T)
        p_on = isfinite(q) && !(q < lo) && !(q > one(T) - lo)
        p_on && (dp[i] = c * (∂α - ∂β))
        dc += ∂α * pc + ∂β * (one(T) - pc)
    end
    return s, (ρ_on ? -dc / ρc^2 : zero(T)), dp
end

function Mooncake.rrule!!(
        ::CoDual{typeof(betabinomial_loglik)},
        trials::CoDual{<:AbstractVector{<:Integer}},
        p::CoDual{<:Array{<:Mooncake.IEEEFloat}},
        ρ::CoDual{<:Mooncake.IEEEFloat},
        obs::CoDual{<:AbstractVector{<:Integer}}
    )
    ## The gradient is taken in the forward pass, since a caller may
    ## overwrite `p` in place before the pullback runs.
    s, dρ, dp = _betabinomial_loglik_grad(
        primal(trials), primal(p), primal(ρ), primal(obs)
    )
    p̄ = tangent(p)
    ## `ρ` is a scalar, so its adjoint goes back as rdata rather than into
    ## a tangent buffer.
    function betabinomial_loglik_pullback!!(s̄)
        p̄ .+= s̄ .* dp
        return NoRData(), NoRData(), NoRData(), s̄ * dρ, NoRData()
    end
    return CoDual(s, NoFData()), betabinomial_loglik_pullback!!
end
