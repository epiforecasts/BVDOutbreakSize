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
        Mooncake.IEEEFloat,
    },
)
Mooncake.@is_primitive(
    Mooncake.MinimalCtx,
    Tuple{
        typeof(patch_infections), Matrix{<:Mooncake.IEEEFloat},
        Vector{<:Mooncake.IEEEFloat}, Matrix{<:Mooncake.IEEEFloat},
        Matrix{<:Mooncake.IEEEFloat}, Matrix{<:Mooncake.IEEEFloat},
        Vector{<:Mooncake.IEEEFloat},
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
        ## The forward adds `w[d] · x[1:n−d+1]` to `y[d:n]` for each lag, so
        ## the pullback is the matching correlation, one lag at a time:
        ##
        ##     x̄[1:n−d+1] += w[d] · ȳ[d:n]
        ##     w̄[d]       += ȳ[d:n] · x[1:n−d+1]
        ##
        ## `x̄` gathers into a contiguous buffer first, since `x` may be a
        ## strided matrix row.
        n = length(xp)
        Δx = zeros(eltype(x̄), n)
        for d in 1:min(length(dp), n)
            ȳd = view(ȳ, d:n)
            axpy!(dp[d], ȳd, view(Δx, 1:(n - d + 1)))
            d̄[d] += dot(ȳd, view(xp, 1:(n - d + 1)))
        end
        x̄ .+= Δx
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
        ## `y[j:j+na−1]` holds `b[j] · a` for each `j`, so
        ##
        ##     ā    += b[j] · ȳ[j:j+na−1]
        ##     b̄[j] += ȳ[j:j+na−1] · a
        na = length(ap)
        for j in eachindex(bp)
            ȳj = view(ȳ, j:(j + na - 1))
            axpy!(bp[j], ȳj, ā)
            b̄[j] += dot(ȳj, ap)
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
            b, frac = _knot_bracket(dayp, t, Tf)
            g = ō[t]
            k̄[b] += (one(Tf) - frac) * g
            k̄[b + 1] += frac * g
        end
        return NoRData(), NoRData(), NoRData(), NoRData()
    end
    return CoDual(out, ō), interpolate_knots_pullback!!
end

## Per day, with `S⁻` the pool before it and `x = R f / N`,
##
##     I = S⁻ (1 − e^{−x}),    S = S⁻ e^{−x},
##
## so walking back `x̄ = e^{−x} S⁻ (Ī − S̄)` and the pool before the day
## collects `S̄⁻ = Ī (1 − e^{−x}) + S̄ e^{−x}`. The pool after the seed is
## `N − Σ seed`, which hands its adjoint to `N` and, negated, to each seed.
## The scalar `N` has no fdata, so its cotangent goes back as rdata.
function Mooncake.rrule!!(
        ::CoDual{typeof(renewal_infections)},
        Rt::CoDual{<:Array{<:Mooncake.IEEEFloat}},
        g::CoDual{<:Array{<:Mooncake.IEEEFloat}},
        seed::CoDual{<:Array{<:Mooncake.IEEEFloat}},
        N::CoDual{<:Mooncake.IEEEFloat}
    )
    Rtp = primal(Rt)
    gp = primal(g)
    seedp = primal(seed)
    Np = primal(N)
    R̄ = tangent(Rt)
    ḡ = tangent(g)
    s̄ = tangent(seed)
    n = length(Rtp)
    L = length(seedp)
    ## The forward pass is the model's own, which also hands back the
    ## per-day force and pool the pullback needs, so the recursion is
    ## defined once.
    st = renewal_infections_with_state(Rtp, gp, seedp, Np)
    I, force, S = st.infections, st.force, st.susceptible
    Ī = zero(I)
    function renewal_infections_pullback!!(::NoRData)
        ## Sequential recursion, so the walk runs backwards: each day's
        ## adjoint must land on the lagged infections before those days
        ## are read. `Ī` is Mooncake's buffer, so accumulate into a copy.
        acc = copy(Ī)
        S̄ = zero(Np)
        N̄ = zero(Np)
        @inbounds for t in n:-1:(L + 1)
            x = Rtp[t] * force[t] / Np
            e = exp(-x)
            x̄ = e * S[t - 1] * (acc[t] - S̄)
            S̄ = -acc[t] * expm1(-x) + S̄ * e
            R̄[t] += x̄ * force[t] / Np
            N̄ -= x̄ * x / Np
            f̄ = x̄ * Rtp[t] / Np
            kmax = min(t - 1, length(gp))
            for s in 1:kmax
                acc[t - s] += f̄ * gp[s]
                ḡ[s] += f̄ * I[t - s]
            end
        end
        m = min(L, n)
        live = Np - sum(view(seedp, 1:m)) > zero(Np)
        live && (N̄ += S̄)
        @inbounds for j in 1:m
            s̄[j] += acc[j] - (live ? S̄ : zero(Np))
        end
        return NoRData(), NoRData(), NoRData(), NoRData(), N̄
    end
    return CoDual(I, Ī), renewal_infections_pullback!!
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
    y, O_bg, flags = _accumulate_occupancy(
        Val(true), primal(A_bvd), primal(A_bg), dp, rp, primal(ruleout), κp,
        hp
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

function Mooncake.rrule!!(
        ::CoDual{typeof(patch_infections)},
        Rt::CoDual{<:Matrix{<:Mooncake.IEEEFloat}},
        g::CoDual{<:Vector{<:Mooncake.IEEEFloat}},
        seeds::CoDual{<:Matrix{<:Mooncake.IEEEFloat}},
        K::CoDual{<:Matrix{<:Mooncake.IEEEFloat}},
        ε::CoDual{<:Matrix{<:Mooncake.IEEEFloat}},
        N::CoDual{<:Vector{<:Mooncake.IEEEFloat}}
    )
    Rp = primal(Rt)
    gp = primal(g)
    Kp = primal(K)
    εp = primal(ε)
    Np = primal(N)
    seedsp = primal(seeds)
    R̄ = tangent(Rt)
    ḡ = tangent(g)
    s̄ = tangent(seeds)
    K̄ = tangent(K)
    ε̄ = tangent(ε)
    N̄ = tangent(N)
    np, n = size(Rp)
    L = size(seedsp, 2)
    st = patch_infections_with_state(Rp, gp, seedsp, Kp, εp, Np)
    y = (; st.infections, st.importation)
    out = Mooncake.zero_fcodual(y)
    Ī = tangent(out).infections
    Ā = tangent(out).importation
    I, S, rate = st.infections, st.susceptible, st.rate
    function patch_infections_pullback!!(::NoRData)
        Tf = eltype(I)
        outflow = _patch_outflow(Tf, Kp, np)
        ## Each day's force and generated infections are rebuilt from the
        ## forward infections rather than stored. The walk runs backwards so
        ## a day's adjoint is complete before it is pushed onto earlier days,
        ## and it accumulates into a copy since `Ī` is Mooncake's buffer.
        acc = copy(Ī)
        force = zeros(Tf, np)
        gen = zeros(Tf, np)
        ḡen = zeros(Tf, np)
        ōut = zeros(Tf, np)
        ȳ = zeros(Tf, np)
        S̄ = zeros(Tf, np)
        @inbounds for t in n:-1:(L + 1)
            kmax = min(t - 1, length(gp))
            for p in 1:np
                f = _patch_force(I, gp, p, t)
                force[p] = f
                gen[p] = Rp[p, t] * f
                ḡen[p] = zero(Tf)
            end
            ## Depletion, as in the single-patch rule, on the rate
            ## `x = y / N[p]` with `y` the force after the transfer.
            for p in 1:np
                x = rate[p, t]
                e = exp(-x)
                x̄ = e * S[p, t - 1] * (acc[p, t] - S̄[p])
                S̄[p] = -acc[p, t] * expm1(-x) + S̄[p] * e
                ȳ[p] = x̄ / Np[p]
                N̄[p] -= x̄ * x / Np[p]
            end
            ## y[p] = (1 - ε[p, t] outflow[p]) gen[p] + arrivals[p], and
            ## the importation output is arrivals[p] alone.
            for p in 1:np
                a = ȳ[p]
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
        ## Each pool after the seed is `N[p] − Σ seeds[p, :]`.
        m = min(L, n)
        @inbounds for p in 1:np
            live = Np[p] - sum(view(seedsp, p, 1:m)) > zero(Tf)
            live && (N̄[p] += S̄[p])
            for j in 1:m
                s̄[p, j] += acc[p, j] - (live ? S̄[p] : zero(Tf))
            end
        end
        return ntuple(_ -> NoRData(), 7)
    end
    return out, patch_infections_pullback!!
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
