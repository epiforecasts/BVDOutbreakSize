## Correctness gate for the hand-written reverse rules in `src/ad_rules.jl`.
##
## A wrong derivative is worse than a slow one, and a hand-written rule
## replaces the one place a backend would otherwise derive it for itself,
## so nothing else in the suite would notice it drifting. Each kernel is
## checked twice:
##
##   1. the pullback against central finite differences, on both (or all
##      three) of its differentiable arguments;
##   2. the Mooncake gradient with the rule registered against the Mooncake
##      gradient of an unregistered clone of the same function body, which
##      is the gradient the backend derives without the rule.
##
## The second check is what proves the rule is both active and faithful:
## the clone shares the body but not the `Tuple{typeof(f), ...}` signature
## the rule is registered against, so Mooncake differentiates it the old
## way.
##
## Tagged `:ad` with the rest of the gradient items.

@testitem "AD rules: kernel pullbacks match finite differences" tags = [:ad] begin
    using Random: seed!
    using FiniteDifferences: central_fdm, grad
    using BVDOutbreakSize: BVDOutbreakSize, convolve_delay, convolve_survival,
        convolve_pmf, interpolate_knots, renewal_infections
    CRC = BVDOutbreakSize.ChainRulesCore

    fdm = central_fdm(5, 1)
    seed!(20260518)

    @testset "convolve_delay" begin
        x = abs.(randn(40)) .+ 0.5
        w = abs.(randn(15)) .+ 0.1
        w ./= sum(w)
        ȳ = randn(40)
        y, pb = CRC.rrule(convolve_delay, x, w)
        @test y == convolve_delay(x, w)
        _, x̄, w̄ = pb(ȳ)
        fx, fw = grad(fdm, (a, b) -> sum(ȳ .* convolve_delay(a, b)), x, w)
        @test x̄ ≈ fx rtol = 1.0e-7
        @test w̄ ≈ fw rtol = 1.0e-7
    end

    @testset "convolve_survival" begin
        x = abs.(randn(40)) .+ 0.5
        l = abs.(randn(12)) .+ 0.1
        l ./= sum(l)
        ȳ = randn(40)
        y, pb = CRC.rrule(convolve_survival, x, l)
        @test y == convolve_survival(x, l)
        _, x̄, l̄ = pb(ȳ)
        fx, fl = grad(fdm, (a, b) -> sum(ȳ .* convolve_survival(a, b)), x, l)
        @test x̄ ≈ fx rtol = 1.0e-7
        @test l̄ ≈ fl rtol = 1.0e-7
    end

    @testset "convolve_pmf" begin
        a = abs.(randn(14)) .+ 0.1
        b = abs.(randn(9)) .+ 0.1
        ȳ = randn(22)
        y, pb = CRC.rrule(convolve_pmf, a, b)
        @test y == convolve_pmf(a, b)
        _, ā, b̄ = pb(ȳ)
        fa, fb = grad(fdm, (u, v) -> sum(ȳ .* convolve_pmf(u, v)), a, b)
        @test ā ≈ fa rtol = 1.0e-7
        @test b̄ ≈ fb rtol = 1.0e-7
    end

    @testset "interpolate_knots" begin
        ## Every knot layout the model builds: a regular weekly grid, a
        ## two-knot span, an irregular grid, and the single-knot window
        ## that takes the `nb == 1` branch.
        for days in (collect(1:7:40), [1, 40], [5, 12, 19, 40], [1])
            n = 40
            kv = randn(length(days))
            ō = randn(n)
            out, pb = CRC.rrule(interpolate_knots, kv, days, n)
            @test out == interpolate_knots(kv, days, n)
            k̄ = pb(ō)[2]
            fk = only(
                grad(fdm, v -> sum(ō .* interpolate_knots(v, days, n)), kv)
            )
            @test k̄ ≈ fk rtol = 1.0e-7
            ## The knot days and the grid length carry no derivative.
            @test pb(ō)[3] isa CRC.NoTangent
            @test pb(ō)[4] isa CRC.NoTangent
        end
    end

    @testset "renewal_infections" begin
        Rt = abs.(randn(40)) .* 0.3 .+ 1.2
        g = abs.(randn(12)) .+ 0.1
        g ./= sum(g)
        seed_vec = abs.(randn(7)) .+ 1.0
        Ī = randn(40)
        Ī_in = copy(Ī)
        I, pb = CRC.rrule(renewal_infections, Rt, g, seed_vec)
        @test I == renewal_infections(Rt, g, seed_vec)
        _, R̄, ḡ, s̄ = pb(Ī)
        fR, fg, fs = grad(
            fdm,
            (a, b, c) -> sum(Ī .* renewal_infections(a, b, c)),
            Rt, g, seed_vec
        )
        @test R̄ ≈ fR rtol = 1.0e-6
        @test ḡ ≈ fg rtol = 1.0e-6
        @test s̄ ≈ fs rtol = 1.0e-6
        ## The pullback must not consume the incoming cotangent in place:
        ## it accumulates the recursion into a working copy of it.
        @test Ī == Ī_in
    end
end

@testitem "AD rules: Mooncake with the rule matches Mooncake without it" tags = [
    :ad,
] begin
    using Random: seed!
    using Mooncake: Mooncake
    using BVDOutbreakSize: convolve_delay, convolve_survival, convolve_pmf,
        interpolate_knots, renewal_infections

    ## Unregistered clones. Same bodies as `src/renewal.jl`, different
    ## signatures, so no rule fires and Mooncake derives the gradient the
    ## way it did before `src/ad_rules.jl` existed.
    function cd_ref(x::AbstractVector, delay::AbstractVector)
        n = length(x)
        Tp = promote_type(eltype(x), eltype(delay))
        y = zeros(Tp, n)
        @inbounds for t in 1:n
            acc = zero(Tp)
            dmax = min(t - 1, length(delay) - 1)
            for d in 0:dmax
                acc += x[t - d] * delay[d + 1]
            end
            y[t] = acc
        end
        return y
    end
    function cs_ref(x::AbstractVector, los::AbstractVector)
        L = length(los)
        surv = similar(los)
        acc = zero(eltype(los))
        @inbounds for i in L:-1:1
            acc += los[i]
            surv[i] = acc
        end
        return cd_ref(x, surv)
    end
    function cp_ref(a::AbstractVector, b::AbstractVector)
        na, nb = length(a), length(b)
        Tp = promote_type(eltype(a), eltype(b))
        y = zeros(Tp, na + nb - 1)
        @inbounds for i in 1:na, j in 1:nb

            y[i + j - 1] += a[i] * b[j]
        end
        return y
    end
    function ik_ref(
            kv::AbstractVector, days::AbstractVector{<:Integer},
            n::Integer
        )
        Tp = eltype(kv)
        out = Vector{Tp}(undef, n)
        nb = length(days)
        nb == 1 && return fill!(out, kv[1])
        @inbounds for t in 1:n
            b = 1
            while b < nb - 1 && t > days[b + 1]
                b += 1
            end
            d0 = days[b]
            d1 = days[b + 1]
            frac = d1 == d0 ? zero(Tp) :
                clamp(Tp(t - d0) / Tp(d1 - d0), zero(Tp), one(Tp))
            out[t] = kv[b] + frac * (kv[b + 1] - kv[b])
        end
        return out
    end
    function ri_ref(Rt::AbstractVector, g::AbstractVector, seed::AbstractVector)
        n = length(Rt)
        L = length(seed)
        Tp = promote_type(eltype(Rt), eltype(g), eltype(seed))
        I = zeros(Tp, n)
        @inbounds for j in 1:min(L, n)
            I[j] = seed[j]
        end
        @inbounds for t in (L + 1):n
            force = zero(Tp)
            kmax = min(t - 1, length(g))
            for s in 1:kmax
                force += I[t - s] * g[s]
            end
            I[t] = Rt[t] * force
        end
        return I
    end

    function mgrad(f, args...)
        rule = Mooncake.build_rrule(f, args...)
        _, g = Mooncake.value_and_gradient!!(rule, f, args...)
        return g[2:end]
    end

    seed!(20260518)
    x = abs.(randn(40)) .+ 0.5
    w = abs.(randn(15)) .+ 0.1
    w ./= sum(w)
    ȳ = randn(40)
    @test all(
        isapprox.(
            mgrad((a, b) -> sum(ȳ .* convolve_delay(a, b)), x, w),
            mgrad((a, b) -> sum(ȳ .* cd_ref(a, b)), x, w);
            rtol = 1.0e-10
        )
    )

    l = abs.(randn(12)) .+ 0.1
    l ./= sum(l)
    @test all(
        isapprox.(
            mgrad((a, b) -> sum(ȳ .* convolve_survival(a, b)), x, l),
            mgrad((a, b) -> sum(ȳ .* cs_ref(a, b)), x, l);
            rtol = 1.0e-10
        )
    )

    a = abs.(randn(14)) .+ 0.1
    b = abs.(randn(9)) .+ 0.1
    z̄ = randn(22)
    @test all(
        isapprox.(
            mgrad((u, v) -> sum(z̄ .* convolve_pmf(u, v)), a, b),
            mgrad((u, v) -> sum(z̄ .* cp_ref(u, v)), a, b);
            rtol = 1.0e-10
        )
    )

    days = collect(1:7:40)
    kv = randn(length(days))
    ō = randn(40)
    @test all(
        isapprox.(
            mgrad(v -> sum(ō .* interpolate_knots(v, days, 40)), kv),
            mgrad(v -> sum(ō .* ik_ref(v, days, 40)), kv);
            rtol = 1.0e-10
        )
    )

    Rt = abs.(randn(40)) .* 0.3 .+ 1.2
    g = abs.(randn(12)) .+ 0.1
    g ./= sum(g)
    seed_vec = abs.(randn(7)) .+ 1.0
    Ī = randn(40)
    @test all(
        isapprox.(
            mgrad(
                (p, q, r) -> sum(Ī .* renewal_infections(p, q, r)),
                Rt, g, seed_vec
            ),
            mgrad((p, q, r) -> sum(Ī .* ri_ref(p, q, r)), Rt, g, seed_vec);
            rtol = 1.0e-10
        )
    )
end
