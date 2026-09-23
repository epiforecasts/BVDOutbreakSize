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
        convolve_pmf, interpolate_knots, renewal_infections,
        abscond_thinned, abscond_thinned_flows, patch_infections,
        nbinomial_loglik, studentt_loglik

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

    @testset "abscond_thinned" begin
        pmf = abs.(randn(15)) .+ 0.1
        pmf ./= sum(pmf)
        κ = 0.07
        ȳ = randn(15)
        y, (p̄, _), rdata = run_rule(abscond_thinned, pmf, κ; cotangent = ȳ)
        @test y == abscond_thinned(pmf, κ)
        fp, fκ = grad(fdm, (p, k) -> sum(ȳ .* abscond_thinned(p, k)), pmf, κ)
        @test p̄ ≈ fp rtol = 1.0e-7
        ## The scalar rate has no fdata, so its cotangent is the rdata.
        @test rdata[3] ≈ fκ rtol = 1.0e-7
        @test rdata[1] isa NoRData && rdata[2] isa NoRData
    end

    @testset "abscond_thinned_flows" begin
        ## Unequal schedule lengths so both branches of the cohort walk run,
        ## and leading zero admissions so the forward pass skips cohorts.
        n = 40
        adm1 = [zeros(3); abs.(randn(n - 3)) .+ 0.5]
        adm2 = [zeros(3); abs.(randn(n - 3)) .+ 0.5]
        pmf1 = abs.(randn(15)) .+ 0.1
        pmf1 ./= sum(pmf1)
        pmf2 = abs.(randn(9)) .+ 0.1
        pmf2 ./= sum(pmf2)
        κ = 0.07
        h = 0.05 .+ 0.3 .* rand(n)
        ȳ1 = randn(n)
        ȳ2 = randn(n)
        args = (adm1, pmf1, adm2, pmf2, κ, h)
        codual_args = map(zero_fcodual, args)
        out, pb = Mooncake.rrule!!(
            zero_fcodual(abscond_thinned_flows), codual_args...
        )
        tangent(out)[1] .= ȳ1
        tangent(out)[2] .= ȳ2
        rdata = pb(NoRData())
        @test primal(out) == abscond_thinned_flows(args...)
        fgs = grad(
            fdm,
            (a, p, b, q, k, c) -> begin
                y1, y2 = abscond_thinned_flows(a, p, b, q, k, c)
                sum(ȳ1 .* y1) + sum(ȳ2 .* y2)
            end,
            args...
        )
        for i in (1, 2, 3, 4, 6)
            @test tangent(codual_args[i]) ≈ fgs[i] rtol = 1.0e-7
        end
        @test rdata[6] ≈ fgs[5] rtol = 1.0e-7
        @test all(i -> rdata[i] isa NoRData, (1, 2, 3, 4, 5, 7))
    end

    @testset "patch_infections" begin
        ## Three patches, a daily importation intensity per origin, and both
        ## outputs carrying a cotangent.
        np, n, L = 3, 40, 7
        Rt = abs.(randn(np, n)) .* 0.3 .+ 1.0
        g = abs.(randn(12)) .+ 0.1
        g ./= sum(g)
        seeds = abs.(randn(np, L)) .+ 1.0
        K = rand(np, np) .* 0.2
        K[[1, 5, 9]] .= 0.0
        ε = rand(np, n) .* 0.5
        Ī = randn(np, n)
        Ā = randn(np, n)
        args = (Rt, g, seeds, K, ε)
        codual_args = map(zero_fcodual, args)
        out, pb = Mooncake.rrule!!(
            zero_fcodual(patch_infections), codual_args...
        )
        tangent(out).infections .= Ī
        tangent(out).importation .= Ā
        rdata = pb(NoRData())
        @test primal(out) == patch_infections(args...)
        @test all(r -> r isa NoRData, rdata)
        objective(a...) = let r = patch_infections(a...)
            sum(Ī .* r.infections) + sum(Ā .* r.importation)
        end
        fd = grad(fdm, objective, args...)
        for (t, f) in zip(map(tangent, codual_args), fd)
            @test t ≈ f rtol = 1.0e-6
        end
        @test tangent(out).infections == Ī
    end

    @testset "nbinomial_loglik" begin
        ## Scalar output and a scalar `k`, so the cotangent goes in as the
        ## pullback's argument and `k`'s adjoint comes back as rdata.
        μ = exp.(4 .+ 0.5 .* randn(60))
        x = rand(0:120, 60)
        x[3] = 0
        for k in (8.3, 0.7, 150.0)
            s̄ = randn()
            codual_μ = zero_fcodual(μ)
            out, pb = Mooncake.rrule!!(
                zero_fcodual(nbinomial_loglik), zero_fcodual(k),
                codual_μ, zero_fcodual(x)
            )
            @test primal(out) == nbinomial_loglik(k, μ, x)
            rdata = pb(s̄)
            fk, fμ = grad(fdm, (a, b) -> s̄ * nbinomial_loglik(a, b, x), k, μ)
            @test rdata[2] ≈ fk rtol = 1.0e-7
            @test tangent(codual_μ) ≈ fμ rtol = 1.0e-7
            ## The function and the counts carry no derivative.
            @test rdata[1] isa NoRData && rdata[3] isa NoRData
            @test rdata[4] isa NoRData
        end
    end

    @testset "studentt_loglik" begin
        ## Scalar output and a scalar `ν`, whose adjoint comes back as rdata.
        ## Increments are integers of either sign, as scanned, or floats, as
        ## simulated, which take a tangent of their own.
        μ = 20 .* randn(60)
        σ = exp.(1 .+ 0.5 .* randn(60))
        x = round.(Int, μ .+ 3 .* σ .* randn(60))
        xf = μ .+ 3 .* σ .* randn(60)
        for (ν, c) in ((4.0, x), (1.5, x), (40.0, x), (4.0, xf))
            s̄ = randn()
            codual_μ = zero_fcodual(μ)
            codual_σ = zero_fcodual(σ)
            codual_c = zero_fcodual(c)
            out, pb = Mooncake.rrule!!(
                zero_fcodual(studentt_loglik), codual_μ, codual_σ,
                codual_c, zero_fcodual(ν)
            )
            @test primal(out) == studentt_loglik(μ, σ, c, ν)
            rdata = pb(s̄)
            fμ, fσ, fν = grad(
                fdm, (a, b, d) -> s̄ * studentt_loglik(a, b, c, d), μ, σ, ν
            )
            @test tangent(codual_μ) ≈ fμ rtol = 1.0e-7
            @test tangent(codual_σ) ≈ fσ rtol = 1.0e-7
            @test rdata[5] ≈ fν rtol = 1.0e-7
            @test all(r -> r isa NoRData, rdata[1:4])
            if c isa Vector{Float64}
                @test tangent(codual_c) ≈ -tangent(codual_μ)
            end
        end
    end
end

@testitem "AD rules: Mooncake with the rule matches ForwardDiff" tags = [
    :ad,
] begin
    using Random: seed!
    using ForwardDiff: ForwardDiff
    using Mooncake: Mooncake
    using BVDOutbreakSize: convolve_delay, convolve_survival, convolve_pmf,
        interpolate_knots, renewal_infections, abscond_thinned,
        abscond_thinned_flow, abscond_thinned_flows, patch_infections,
        nbinomial_loglik, studentt_loglik, censored_nbinomial_loglik,
        NegBinomialVector, StudentTVector
    using Distributions: logpdf, censored

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
    ## the kernel body in `src/renewal.jl` or `src/models/` itself.
    function fgrad(f, args...)
        return ntuple(length(args)) do i
            fi(v) = f(ntuple(j -> j == i ? v : args[j], length(args))...)
            args[i] isa Real ? ForwardDiff.derivative(fi, args[i]) :
                ForwardDiff.gradient(fi, args[i])
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

    @testset "abscond_thinned" begin
        pmf = abs.(randn(15)) .+ 0.1
        pmf ./= sum(pmf)
        ȳ = randn(15)
        check_grads((p, k) -> sum(ȳ .* abscond_thinned(p, k)), pmf, 0.07)
    end

    @testset "abscond_thinned_flows" begin
        ## Leading zero admissions: the forward pass skips those cohorts, and
        ## ForwardDiff still sees their derivative in the admissions.
        n = 40
        adm1 = [zeros(3); abs.(randn(n - 3)) .+ 0.5]
        adm2 = [zeros(3); abs.(randn(n - 3)) .+ 0.5]
        h = 0.05 .+ 0.3 .* rand(n)
        ȳ1 = randn(n)
        ȳ2 = randn(n)
        obj = (a, p, b, q, k, c) -> begin
            y1, y2 = abscond_thinned_flows(a, p, b, q, k, c)
            sum(ȳ1 .* y1) + sum(ȳ2 .* y2)
        end
        ## Each schedule in turn the longer, and one longer than the series
        ## so every cohort is truncated.
        for (l1, l2) in ((15, 9), (9, 15), (50, 12))
            pmf1 = abs.(randn(l1)) .+ 0.1
            pmf1 ./= sum(pmf1)
            pmf2 = abs.(randn(l2)) .+ 0.1
            pmf2 ./= sum(pmf2)
            check_grads(obj, adm1, pmf1, adm2, pmf2, 0.07, h)
        end
        ## The single-flow wrapper passes the same arrays as both flows, so
        ## the rule's tangents alias.
        pmf = abs.(randn(15)) .+ 0.1
        pmf ./= sum(pmf)
        check_grads(
            (a, p, k, c) -> sum(ȳ1 .* abscond_thinned_flow(a, p, k, c)),
            adm1, pmf, 0.07, h
        )
    end

    ## The per-patch models pass a row of a matrix, a view rather than an
    ## `Array`. The rule must fire on it, not only agree with it.
    fires(sig) = Mooncake.is_primitive(
        Mooncake.MinimalCtx, Mooncake.ReverseMode, sig,
        Base.get_world_counter()
    )

    @testset "convolve_delay on a matrix row" begin
        M = abs.(randn(3, 40)) .+ 0.5
        w = abs.(randn(15)) .+ 0.1
        w ./= sum(w)
        Ȳ = randn(3, 40)
        @test fires(
            Tuple{typeof(convolve_delay), typeof(view(M, 1, :)), typeof(w)}
        )
        check_grads(
            (A, b) -> sum(
                p -> sum(Ȳ[p, :] .* convolve_delay(view(A, p, :), b)), 1:3
            ),
            M, w
        )
    end

    @testset "interpolate_knots on a matrix row" begin
        days = collect(1:7:40)
        n = 40
        δ = randn(3, length(days))
        Ō = randn(3, n)
        @test fires(
            Tuple{
                typeof(interpolate_knots), typeof(view(δ, 1, :)),
                typeof(days), typeof(n),
            }
        )
        check_grads(
            A -> sum(
                p -> sum(Ō[p, :] .* interpolate_knots(view(A, p, :), days, n)),
                1:3
            ),
            δ
        )
    end

    @testset "patch_infections" begin
        np, n, L = 3, 40, 7
        Rt = abs.(randn(np, n)) .* 0.3 .+ 1.0
        g = abs.(randn(12)) .+ 0.1
        g ./= sum(g)
        seeds = abs.(randn(np, L)) .+ 1.0
        K = rand(np, np) .* 0.2
        K[[1, 5, 9]] .= 0.0
        ε = rand(np, n) .* 0.5
        Ī = randn(np, n)
        Ā = randn(np, n)
        check_grads(
            (a...) -> let r = patch_infections(a...)
                sum(Ī .* r.infections) + sum(Ā .* r.importation)
            end,
            Rt, g, seeds, K, ε
        )
    end

    @testset "nbinomial_loglik" begin
        μ = exp.(4 .+ 0.5 .* randn(60))
        x = rand(0:120, 60)
        x[3] = 0
        ## A mean at the `safe_rate` floor passes no derivative. ForwardDiff's
        ## own quotient rule cancels there when `k` is small, so the floored
        ## means are checked at a moderate `k` only.
        μ_floor = vcat(μ, [0.0, 1.0e-18])
        x_floor = vcat(x, [3, 40])
        for (k, m, c) in ((8.3, μ_floor, x_floor), (0.4, μ, x), (250.0, μ, x))
            m_grad = mgrad((a, b) -> nbinomial_loglik(a, b, c), k, m)
            @test m_grad[1] ≈ ForwardDiff.derivative(
                a -> nbinomial_loglik(a, m, c), k
            ) rtol = 1.0e-10
            @test m_grad[2] ≈ ForwardDiff.gradient(
                b -> nbinomial_loglik(k, b, c), m
            ) rtol = 1.0e-10
        end
        ## `k` large enough that `p` clamps at `1 − eps` for every count:
        ## the clamp passes no derivative to the means.
        m_grad = mgrad((b) -> nbinomial_loglik(1.0e20, b, x), μ)
        @test all(iszero, m_grad[1])
    end

    @testset "studentt_loglik" begin
        μ = 20 .* randn(60)
        σ = exp.(1 .+ 0.5 .* randn(60))
        x = round.(Int, μ .+ 3 .* σ .* randn(60))
        ## The call site's argument types, so the rule serves the model.
        @test fires(
            Tuple{
                typeof(studentt_loglik), typeof(μ), typeof(σ), typeof(x),
                Float64,
            }
        )
        ## Scales at the `safe_studentt` floor (negative, non-finite) pass no
        ## derivative to the scale but still pass one to the mean.
        μ_floor = vcat(μ, [-2.0, 1.0, 0.5])
        σ_floor = vcat(σ, [-1.0, NaN, Inf])
        x_floor = vcat(x, [1, -4, 2])
        for ν in (4.0, 1.5, 40.0)
            check_grads(
                (a, b, d) -> studentt_loglik(a, b, x_floor, d),
                μ_floor, σ_floor, ν
            )
        end
        ## ForwardDiff's `Dual` comparison reads a scale of exactly zero as
        ## above the floor, so a zero scale is checked in the mean alone.
        μ0, σ0, x0 = vcat(μ_floor, 3.0), vcat(σ_floor, 0.0), vcat(x_floor, 5)
        g = mgrad((a, b) -> studentt_loglik(a, b, x0, 4.0), μ0, σ0)
        @test all(iszero, g[2][61:64])
        @test all(!iszero, g[1][61:64])
        @test g[1] ≈ ForwardDiff.gradient(
            a -> studentt_loglik(a, σ0, x0, 4.0), μ0
        ) rtol = 1.0e-10
        ## Float increments take a tangent too.
        xf = μ .+ 3 .* σ .* randn(60)
        check_grads(studentt_loglik, μ, σ, xf, 4.0)
        ## A defaulted `ν` passes no derivative.
        @test iszero(mgrad(d -> studentt_loglik(μ, σ, x, d), -1.0)[1])
    end

    @testset "vector distributions" begin
        ## `logpdf` of each vector distribution is its summed helper, so the
        ## helper's rule fires through it: the gradient is the rule's own,
        ## bit for bit, and matches ForwardDiff.
        μ = exp.(4 .+ 0.5 .* randn(60))
        x = rand(0:120, 60)
        nb(a, b) = logpdf(NegBinomialVector(a, b), x)
        @test mgrad(nb, 8.3, μ) == mgrad(
            (a, b) -> nbinomial_loglik(a, b, x), 8.3, μ
        )
        check_grads(nb, 8.3, μ)
        ## Counts below their ceilings only: Mooncake has no rule for the
        ## censored tail's Rmath call.
        up = fill(1.0e6, 60)
        @test mgrad(
            (a, b) -> logpdf(censored(NegBinomialVector(a, b); upper = up), x),
            8.3, μ
        ) == mgrad(
            (a, b) -> censored_nbinomial_loglik(a, b, up, x), 8.3, μ
        )
        m = 20 .* randn(60)
        σ = exp.(1 .+ 0.5 .* randn(60))
        y = round.(Int, m .+ 3 .* σ .* randn(60))
        st(a, b, d) = logpdf(StudentTVector(a, b, d), y)
        @test mgrad(st, m, σ, 4.0) == mgrad(
            (a, b, d) -> studentt_loglik(a, b, y, d), m, σ, 4.0
        )
        check_grads(st, m, σ, 4.0)
    end
end

@testitem "AD rules: treatment cohort kernels match ForwardDiff" tags = [
    :ad,
] begin
    using Random: seed!
    using ForwardDiff: ForwardDiff
    using Mooncake: Mooncake
    using BVDOutbreakSize: clinical_stay_survival, two_clock_confirmed

    ## Mooncake's gradient, with the registered rule, against ForwardDiff
    ## on the kernel body. `cfr` is a scalar, so its ForwardDiff reference
    ## is a derivative rather than a gradient.
    function mgrad(f, args...)
        rule = Mooncake.build_rrule(f, args...)
        _, g = Mooncake.value_and_gradient!!(rule, f, args...)
        return g[2:end]
    end
    function fgrad(f, args...)
        return ntuple(length(args)) do i
            fi = v -> f(ntuple(j -> j == i ? v : args[j], length(args))...)
            args[i] isa Real ? ForwardDiff.derivative(fi, args[i]) :
                ForwardDiff.gradient(fi, args[i])
        end
    end
    function check_grads(f, args...)
        for (mi, di) in zip(mgrad(f, args...), fgrad(f, args...))
            @test mi ≈ di rtol = 1.0e-10
        end
        return nothing
    end

    seed!(20260923)

    @testset "clinical_stay_survival" begin
        ## Unequal PMF lengths, either way round, exercise the zero tail.
        for (Ld, Lr) in ((18, 12), (10, 16), (14, 14))
            dp = abs.(randn(Ld)) .+ 0.1
            dp ./= 1.05 * sum(dp)
            rp = abs.(randn(Lr)) .+ 0.1
            rp ./= 1.02 * sum(rp)
            S̄ = randn(max(Ld, Lr))
            check_grads(
                (a, b, c) -> sum(S̄ .* clinical_stay_survival(a, b, c)),
                dp, rp, 0.37
            )
        end
    end

    @testset "two_clock_confirmed" begin
        ## A survival longer than the grid stops the cohort walk at day 1
        ## rather than at the support.
        for (n, L) in ((40, 15), (10, 15))
            A = abs.(randn(n)) .* 10 .+ 1
            h = rand(n) .* 0.3
            S = sort(rand(L); rev = true)
            Ō = randn(n)
            check_grads(
                (a, b, c) -> sum(Ō .* two_clock_confirmed(a, b, c)),
                A, h, S
            )
        end
    end

    @testset "composed as the treatment model calls them" begin
        n = 40
        A = abs.(randn(n)) .* 10 .+ 1
        h = rand(n) .* 0.3
        dp = abs.(randn(18)) .+ 0.1
        dp ./= 1.05 * sum(dp)
        rp = abs.(randn(12)) .+ 0.1
        rp ./= 1.02 * sum(rp)
        Ō = randn(n)
        check_grads(
            (a, b, d, r, c) -> sum(
                Ō .* two_clock_confirmed(
                    a, b, clinical_stay_survival(d, r, c)
                )
            ),
            A, h, dp, rp, 0.37
        )
    end
end

@testitem "AD rules: accumulate_occupancy matches ForwardDiff" tags = [
    :ad,
] begin
    using Random: seed!
    using ForwardDiff: ForwardDiff
    using Mooncake: Mooncake, NoRData, primal, tangent, zero_fcodual
    using BVDOutbreakSize: accumulate_occupancy, convolve_delay,
        _accumulate_occupancy_taped, _OCC_CONF_HI

    ## Admissions, discharge schedules and a confirmation hazard shaped like
    ## the treatment-flow model's. `dmult` above one discharges more than
    ## was admitted, so the stocks hit their zero floors. No BVD admissions
    ## on the last day under a high flat hazard `hflat` confirms more than
    ## the stock holds, so the confirmed clamp takes its upper side there.
    ## Only the last day, as the next day's `max(O_bvd − O_conf, 0)` would
    ## sit on its tie, where ForwardDiff and Mooncake take different sides.
    function occupancy_inputs(
            n; dmult = 1.0, κ = 0.02, stop = n, hflat = nothing
        )
        days = 1:n
        A_bvd = @. 20 * exp(-((days - n / 2) / 20)^2) + 0.5
        A_bvd[(stop + 1):end] .= 0
        A_bg = @. 30 * exp(-((days - n / 2) / 30)^2) + 2.0
        w(L) = (p = abs.(randn(L)) .+ 0.1; p ./ sum(p))
        deaths = dmult .* convolve_delay(0.4 .* A_bvd, w(15))
        recover = dmult .* convolve_delay(0.6 .* A_bvd, w(20))
        ruleout = dmult .* convolve_delay(A_bg, w(8))
        h = isnothing(hflat) ? 0.1 .+ 0.3 .* abs.(sin.(days ./ 7)) :
            fill(hflat, n)
        return A_bvd, A_bg, deaths, recover, ruleout, κ, h
    end

    ## Mooncake with the rule against ForwardDiff on the kernel body, all
    ## seven arguments at once with every output weighted.
    function check_occupancy(args)
        n = length(args[1])
        W = [randn(n) for _ in 1:5]
        loss = (a, b, c, d, e, k, h) -> begin
            o = accumulate_occupancy(a, b, c, d, e, k, h)
            sum(W[1] .* o.demand) + sum(W[2] .* o.O_bvd) +
                sum(W[3] .* o.O_conf) + sum(W[4] .* o.O_susp) +
                sum(W[5] .* o.abscond)
        end
        rule = Mooncake.build_rrule(loss, args...)
        _, mg = Mooncake.value_and_gradient!!(rule, loss, args...)
        lens = map(length, args)
        x = reduce(vcat, map(a -> a isa Number ? [a] : a, args))
        unpack(v) = let i = Ref(0)
            map(args, lens) do a, L
                seg = v[(i[] + 1):(i[] + L)]
                i[] += L
                a isa Number ? only(seg) : seg
            end
        end
        fg = unpack(ForwardDiff.gradient(v -> loss(unpack(v)...), x))
        for (mi, di) in zip(mg[2:end], fg)
            @test mi ≈ di rtol = 1.0e-10
        end
        return nothing
    end

    seed!(20260923)
    scenarios = (
        (; dmult = 1.0, κ = 0.02), (; dmult = 1.6, κ = 0.08),
        (; κ = 0.02, stop = 59, hflat = 0.9),
    )
    for kw in scenarios
        args = occupancy_inputs(60; kw...)
        ## The rule's forward pass is its own copy of the balance, so it
        ## must reproduce the model's outputs bit for bit.
        y = accumulate_occupancy(args...)
        ỹ, _, flags = _accumulate_occupancy_taped(args...)
        @test all(k -> getfield(y, k) == getfield(ỹ, k), keys(y))
        haskey(kw, :hflat) && @test any(f -> f & _OCC_CONF_HI != 0, flags)
        ## The rule fires: the output tangent is the rule's `NamedTuple`.
        out, _ = Mooncake.rrule!!(
            zero_fcodual(accumulate_occupancy), map(zero_fcodual, args)...
        )
        @test tangent(out) isa NamedTuple
        check_occupancy(args)
    end
end

@testitem "AD rules: onset-reporting kernels match ForwardDiff" tags = [
    :ad,
] begin
    using Random: seed!
    using ForwardDiff: ForwardDiff
    using Mooncake: Mooncake
    using BVDOutbreakSize: onset_report_cdf_table,
        onset_report_anchor_series, onset_report_moments,
        onset_report_expected_total

    ## Same pairing as above: Mooncake runs the registered rule, ForwardDiff
    ## the kernel body in `src/models/observations.jl`.
    function check_grads(f, args...)
        rule = Mooncake.build_rrule(f, args...)
        v, g = Mooncake.value_and_gradient!!(rule, f, args...)
        ## The table and total rules rebuild the forward loop themselves.
        @test v == f(args...)
        for i in eachindex(args)
            d = ForwardDiff.gradient(args[i]) do v
                f(ntuple(j -> j == i ? v : args[j], length(args))...)
            end
            ## `atol` covers the constant anchor, whose table gradient is
            ## zero up to round-off.
            @test g[i + 1] ≈ d rtol = 1.0e-10 atol = 1.0e-12
        end
        return nothing
    end

    seed!(20260923)
    D = 12
    gs = 5
    ge = 40
    n = 45
    nu = ge - gs + 1
    lh = randn(D) .- 2
    ## Shorter than the table's reach, so the calendar index clamps at the
    ## top end as it does for onset dates near `grid_end`.
    γ = 0.3 .* randn(nu)
    onsets = abs.(randn(n)) .* 20
    alpha = abs.(randn(nu)) .* 0.3
    tab = onset_report_cdf_table(lh, γ, gs, gs, ge)

    @testset "onset_report_cdf_table" begin
        W = randn(D, nu)
        check_grads(
            (l, g) -> sum(W .* onset_report_cdf_table(l, g, gs, gs, ge)),
            lh, γ
        )
    end

    @testset "onset_report_anchor_series" begin
        w = randn(nu)
        ## A daily anchor series and the length-1 constant default.
        for a in (abs.(randn(n)) .* 0.3, [0.15])
            check_grads(
                (t, v) -> sum(w .* onset_report_anchor_series(t, gs, v)),
                tab, a
            )
        end
    end

    @testset "onset_report_moments" begin
        ## Cells from a first vintage (previous report date before every
        ## onset, so a negative delay), from later vintages, and past the
        ## delay support on both report dates.
        oi = Int[]
        ri = Int[]
        pri = Int[]
        for (R, Rp) in ((gs + 4, 0), (22, gs + 4), (31, 22), (ge, 31))
            for u in max(R - D + 1, gs):R
                push!(oi, u)
                push!(ri, R)
                push!(pri, Rp)
            end
        end
        append!(oi, [gs, gs + 1])
        append!(ri, [ge, ge])
        append!(pri, [ge - 1, gs + D + 3])
        w1, w2, w3 = randn(length(oi)), randn(length(oi)), randn(length(oi))
        check_grads(
            (t, o, al) -> begin
                r = onset_report_moments(t, gs, o, gs, al, oi, ri, pri)
                sum(w1 .* r.means) + sum(w2 .* r.level_cur) +
                    sum(w3 .* r.level_prev)
            end,
            tab, onsets, alpha
        )
    end

    @testset "onset_report_expected_total" begin
        ## Onset dates before `grid_start` clamp both `γ` and `alpha`.
        check_grads(
            (o, l, g, al) -> onset_report_expected_total(o, l, g, gs, al, n),
            onsets, lh, γ, alpha
        )
    end
end

@testitem "AD rules: betabinomial_loglik passes Mooncake's test_rule" tags = [
    :ad,
] begin
    using Random: Xoshiro
    using Mooncake: Mooncake
    using Mooncake.TestUtils: test_rule
    using BVDOutbreakSize: betabinomial_loglik

    ## The composition hands its observed row as a `Vector{Int}` or, from a
    ## matrix, a view. The rule must fire on both.
    fires(sig) = Mooncake.is_primitive(
        Mooncake.MinimalCtx, Mooncake.ReverseMode, sig,
        Base.get_world_counter()
    )
    row = view(zeros(Int, 2, 2), 1, :)
    for obs_type in (Vector{Int}, typeof(row))
        @test fires(
            Tuple{
                typeof(betabinomial_loglik), Vector{Int}, Vector{Float64},
                Float64, obs_type,
            }
        )
    end

    rng = Xoshiro(20260923)
    n = rand(rng, 0:300, 30)
    p = 0.05 .+ 0.9 .* rand(rng, 30)
    x = [rand(rng, 0:n[i]) for i in 1:30]
    ## A window with no trials.
    n[3] = 0
    x[3] = 0
    ## Probabilities outside `[0, 1]` sit in the clamps' flat region, so they
    ## pass no derivative and finite differences agree.
    p_clamp = copy(p)
    p_clamp[[1, 2]] .= (-0.5, 1.5)
    x_clamp = copy(x)
    x_clamp[1] = 0
    x_clamp[2] = n[2]
    M = zeros(Int, 2, 30)
    M[1, :] .= x
    ## `ρ` below its floor or above its cap is clamped and passes no
    ## derivative.
    for (q, ρ, obs) in (
            (p, 0.05, x), (p, 0.3, x), (p, 0.8, view(M, 1, :)),
            (p_clamp, 0.05, x_clamp), (p, 1.0e-7, x), (p, 1.5, x),
        )
        test_rule(
            rng, betabinomial_loglik, n, q, ρ, obs;
            is_primitive = true, mode = Mooncake.ReverseMode
        )
    end
end
