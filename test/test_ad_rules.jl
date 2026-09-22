## Correctness gate for the hand-written reverse rules in `src/ad_rules.jl`.
##
## A wrong derivative is worse than a slow one, and a hand-written rule
## replaces the one place a backend would otherwise derive it for itself,
## so nothing else in the suite would notice it drifting. Each kernel is
## checked twice:
##
##   1. the pullback against central finite differences, on both (or all
##      three) of its differentiable arguments;
##   2. the Mooncake gradient, which the rule serves, against ForwardDiff
##      on the same objective.
##
## The second check proves the rule is both active and faithful.
## ForwardDiff does not consult Mooncake's rule table, so it differentiates
## the kernel body itself, and it is exact, so the two must agree to
## round-off rather than to a finite-difference tolerance.
##
## Tagged `:ad` with the rest of the gradient items.

@testitem "AD rules: kernel pullbacks match finite differences" tags = [:ad] begin
    using Random: seed!
    using FiniteDifferences: central_fdm, grad
    using Mooncake: Mooncake, NoRData, primal, tangent, zero_fcodual
    using BVDOutbreakSize: convolve_delay, convolve_survival,
        convolve_pmf, interpolate_knots, renewal_infections

    ## Drive one native rule: build zero-tangent coduals for the arguments,
    ## seed the output tangent with the cotangent, run the pullback, and
    ## hand back the argument tangents it accumulated into along with the
    ## rdata tuple it returned.
    function run_rule(f, args...; cotangent)
        codual_args = map(zero_fcodual, args)
        out, pb = Mooncake.rrule!!(zero_fcodual(f), codual_args...)
        tangent(out) .= cotangent
        rdata = pb(NoRData())
        return primal(out), map(tangent, codual_args), rdata, tangent(out)
    end

    fdm = central_fdm(5, 1)
    seed!(20260518)

    @testset "convolve_delay" begin
        x = abs.(randn(40)) .+ 0.5
        w = abs.(randn(15)) .+ 0.1
        w ./= sum(w)
        ȳ = randn(40)
        y, (x̄, w̄), _ = run_rule(convolve_delay, x, w; cotangent = ȳ)
        @test y == convolve_delay(x, w)
        fx, fw = grad(fdm, (a, b) -> sum(ȳ .* convolve_delay(a, b)), x, w)
        @test x̄ ≈ fx rtol = 1.0e-7
        @test w̄ ≈ fw rtol = 1.0e-7
    end

    @testset "convolve_survival" begin
        x = abs.(randn(40)) .+ 0.5
        l = abs.(randn(12)) .+ 0.1
        l ./= sum(l)
        ȳ = randn(40)
        y, (x̄, l̄), _ = run_rule(convolve_survival, x, l; cotangent = ȳ)
        @test y == convolve_survival(x, l)
        fx, fl = grad(fdm, (a, b) -> sum(ȳ .* convolve_survival(a, b)), x, l)
        @test x̄ ≈ fx rtol = 1.0e-7
        @test l̄ ≈ fl rtol = 1.0e-7
    end

    @testset "convolve_pmf" begin
        a = abs.(randn(14)) .+ 0.1
        b = abs.(randn(9)) .+ 0.1
        ȳ = randn(22)
        y, (ā, b̄), _ = run_rule(convolve_pmf, a, b; cotangent = ȳ)
        @test y == convolve_pmf(a, b)
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
            out, tangents, rdata = run_rule(
                interpolate_knots, kv, days, n; cotangent = ō
            )
            @test out == interpolate_knots(kv, days, n)
            k̄ = tangents[1]
            fk = only(
                grad(fdm, v -> sum(ō .* interpolate_knots(v, days, n)), kv)
            )
            @test k̄ ≈ fk rtol = 1.0e-7
            ## The knot days and the grid length carry no derivative, and
            ## neither does the function itself.
            @test length(rdata) == 4
            @test all(r -> r isa NoRData, rdata)
        end
    end

    @testset "renewal_infections" begin
        Rt = abs.(randn(40)) .* 0.3 .+ 1.2
        g = abs.(randn(12)) .+ 0.1
        g ./= sum(g)
        seed_vec = abs.(randn(7)) .+ 1.0
        Ī = randn(40)
        I, (R̄, ḡ, s̄), _, out_tangent = run_rule(
            renewal_infections, Rt, g, seed_vec; cotangent = Ī
        )
        @test I == renewal_infections(Rt, g, seed_vec)
        fR, fg, fs = grad(
            fdm,
            (a, b, c) -> sum(Ī .* renewal_infections(a, b, c)),
            Rt, g, seed_vec
        )
        @test R̄ ≈ fR rtol = 1.0e-6
        @test ḡ ≈ fg rtol = 1.0e-6
        @test s̄ ≈ fs rtol = 1.0e-6
        ## The pullback must not consume the output tangent in place: it
        ## accumulates the recursion into a working copy of it.
        @test out_tangent == Ī
    end
end

@testitem "AD rules: Mooncake with the rule matches ForwardDiff" tags = [
    :ad,
] begin
    using Random: seed!
    using ForwardDiff: ForwardDiff
    using Mooncake: Mooncake
    using BVDOutbreakSize: convolve_delay, convolve_survival, convolve_pmf,
        interpolate_knots, renewal_infections

    ## Mooncake's gradient of `f` with respect to each of its arguments.
    ## The registered rule fires here, so this is the gradient the model
    ## actually sees.
    function mgrad(f, args...)
        rule = Mooncake.build_rrule(f, args...)
        _, g = Mooncake.value_and_gradient!!(rule, f, args...)
        return g[2:end]
    end

    ## ForwardDiff's gradient of the same objective, one argument at a time.
    ## ForwardDiff never consults Mooncake's rule table, so it differentiates
    ## the kernel body in `src/renewal.jl` itself.
    function fgrad(f, args...)
        return ntuple(length(args)) do i
            ForwardDiff.gradient(args[i]) do v
                f(ntuple(j -> j == i ? v : args[j], length(args))...)
            end
        end
    end

    ## Both are exact methods on the same body, so they agree to round-off.
    function check_grads(f, args...)
        m = mgrad(f, args...)
        d = fgrad(f, args...)
        for (mi, di) in zip(m, d)
            @test mi ≈ di rtol = 1.0e-10
        end
        return nothing
    end

    seed!(20260518)

    @testset "convolve_delay" begin
        x = abs.(randn(40)) .+ 0.5
        w = abs.(randn(15)) .+ 0.1
        w ./= sum(w)
        ȳ = randn(40)
        check_grads((a, b) -> sum(ȳ .* convolve_delay(a, b)), x, w)
    end

    @testset "convolve_survival" begin
        x = abs.(randn(40)) .+ 0.5
        l = abs.(randn(12)) .+ 0.1
        l ./= sum(l)
        ȳ = randn(40)
        check_grads((a, b) -> sum(ȳ .* convolve_survival(a, b)), x, l)
    end

    @testset "convolve_pmf" begin
        a = abs.(randn(14)) .+ 0.1
        b = abs.(randn(9)) .+ 0.1
        z̄ = randn(22)
        check_grads((u, v) -> sum(z̄ .* convolve_pmf(u, v)), a, b)
    end

    @testset "interpolate_knots" begin
        ## The knot days and the grid length carry no derivative, so they
        ## stay closed over rather than differentiated.
        for days in (collect(1:7:40), [1, 40], [5, 12, 19, 40], [1])
            n = 40
            kv = randn(length(days))
            ō = randn(n)
            check_grads(v -> sum(ō .* interpolate_knots(v, days, n)), kv)
        end
    end

    @testset "renewal_infections" begin
        Rt = abs.(randn(40)) .* 0.3 .+ 1.2
        g = abs.(randn(12)) .+ 0.1
        g ./= sum(g)
        seed_vec = abs.(randn(7)) .+ 1.0
        Ī = randn(40)
        check_grads(
            (p, q, r) -> sum(Ī .* renewal_infections(p, q, r)),
            Rt, g, seed_vec
        )
    end
end
