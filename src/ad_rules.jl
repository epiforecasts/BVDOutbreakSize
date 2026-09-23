# Native `Mooncake.rrule!!` methods for the `renewal.jl` kernels and the
# observation kernels in `models/observations.jl`. Each is a loop over a
# daily grid, so left to the backend every iteration's intermediates go on
# the tape; these replace that with a closed-form adjoint of the same shape
# as the forward loop.
#
# Signatures are restricted to `Array{<:IEEEFloat}` and `IEEEFloat`
# scalars, what every call site passes. Anything else falls through to
# Mooncake's derived rule.
#
# These earn their place only while the backend's own derivation is worse.
# `task benchmark-rules` times both arms; if the gap has closed after a
# backend upgrade, delete this file rather than maintain it.

using Mooncake: CoDual, NoRData, primal, tangent

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
