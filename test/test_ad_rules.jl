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
        abscond_thinned, abscond_thinned_flows

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
end

@testitem "AD rules: Mooncake with the rule matches ForwardDiff" tags = [
    :ad,
] begin
    using Random: seed!
    using ForwardDiff: ForwardDiff
    using Mooncake: Mooncake
    using BVDOutbreakSize: convolve_delay, convolve_survival, convolve_pmf,
        interpolate_knots, renewal_infections, abscond_thinned,
        abscond_thinned_flow, abscond_thinned_flows

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
