## Correctness and speed gates for the hand-written reverse rules in
## `src/ad_rules.jl`.
##
## A wrong derivative is worse than a slow one, and a hand-written rule
## replaces the one place a backend would otherwise derive it for itself,
## so nothing else in the suite would notice it drifting. `ADRuleCases`
## lists every rule with the argument types its call sites pass and the
## edge cases its guards take. Mooncake's `test_rule` checks each case
## against finite differences with `is_primitive = true`, which also proves
## the rule fires for that signature. The speed item times one
## production-sized case per rule against Mooncake's own derivation.
##
## Tagged `:ad` with the rest of the gradient items, and the speed item
## also `:ad_perf`.

@testsnippet ADRuleCases begin
    using Random: Xoshiro
    using BVDOutbreakSize: convolve_delay, convolve_pmf, interpolate_knots,
        knot_days, renewal_infections, patch_infections, nbinomial_loglik,
        abscond_thinned, abscond_thinned_flows, two_clock_confirmed,
        clinical_stay_survival, accumulate_occupancy,
        onset_report_cdf_table, onset_report_anchor_series,
        onset_report_moments, onset_report_expected_total, studentt_loglik,
        betabinomial_loglik, onset_vintage_indices, censoring_cap,
        admission_headroom

    ## A positive PMF of length `L` with total mass `mass`.
    pmf(rng, L; mass = 1.0) = (p = rand(rng, L) .+ 0.1; p .* (mass / sum(p)))

    ## Admissions, discharges and a confirmation hazard shaped like the
    ## treatment-flow model's. `dmult` above one discharges more than was
    ## admitted, so the stocks sit on their zero floors. No BVD admissions on
    ## the last day under a high flat hazard `hflat` confirms more than the
    ## stock holds, so the confirmed clamp takes its upper side there. It is
    ## the last day only, because the next day's `max(O_bvd − O_conf, 0)`
    ## would sit on its tie, where finite differences straddle the kink.
    function occupancy_args(
            rng, n; dmult = 1.0, κ = 0.02, stop = n, hflat = nothing
        )
        days = 1:n
        A_bvd = @. 20 * exp(-((days - n / 2) / 20)^2) + 0.5
        A_bvd[(stop + 1):end] .= 0
        A_bg = @. 30 * exp(-((days - n / 2) / 30)^2) + 2.0
        deaths = dmult .* convolve_delay(0.4 .* A_bvd, pmf(rng, 15))
        recover = dmult .* convolve_delay(0.6 .* A_bvd, pmf(rng, 20))
        ruleout = dmult .* convolve_delay(A_bg, pmf(rng, 8))
        h = isnothing(hflat) ? 0.1 .+ 0.3 .* abs.(sin.(days ./ 7)) :
            fill(hflat, n)
        return (A_bvd, A_bg, deaths, recover, ruleout, κ, h)
    end

    ## Onset-reporting inputs over a `D`-day delay support and the onset
    ## grid `gs:ge`. `γ` is shorter than the table's reach, so the calendar
    ## index clamps at the top end as it does near `grid_end`. The scored
    ## cells come from a first vintage (previous report date before every
    ## onset, so a negative delay), from later vintages, and past the delay
    ## support on both report dates.
    function onset_args(rng; D, gs, ge, n)
        nu = ge - gs + 1
        lh = randn(rng, D) .- 2
        γ = 0.3 .* randn(rng, nu)
        tab = onset_report_cdf_table(lh, γ, gs, gs, ge)
        oi, ri, pri = Int[], Int[], Int[]
        vintages = round.(Int, range(gs + 4, ge; length = 4))
        for (R, Rp) in zip(vintages, [0; vintages[1:(end - 1)]])
            for u in max(R - D + 1, gs):R
                push!(oi, u)
                push!(ri, R)
                push!(pri, Rp)
            end
        end
        append!(oi, [gs, gs + 1])
        append!(ri, [ge, ge])
        append!(pri, [ge - 1, gs + D + 3])
        return (;
            lh, γ, tab, oi, ri, pri,
            onsets = abs.(randn(rng, n)) .* 20,
            alpha = abs.(randn(rng, nu)) .* 0.3,
        )
    end

    ## Every native rule with the argument types its call sites pass: float
    ## arrays, matrix-row views, integer observations and scalar κ, cfr, k,
    ## ρ and ν. `perf` marks the one production-sized case per rule the
    ## speed item times.
    function rule_cases(rng)
        cases = NamedTuple[]
        add!(note, f, args...; perf = false) = push!(
            cases, (; name = "$(nameof(f)): $note", f, args, perf)
        )

        M = rand(rng, 3, 40) .+ 0.5
        add!("D = 15", convolve_delay, rand(rng, 40) .+ 0.5, pmf(rng, 15))
        add!("n < D", convolve_delay, rand(rng, 5), pmf(rng, 12))
        add!("D = 1", convolve_delay, rand(rng, 20), [1.0])
        add!("empty series", convolve_delay, Float64[], pmf(rng, 5))
        w0 = pmf(rng, 10)
        w0[[1, 4]] .= 0
        add!("zero-weight lags", convolve_delay, rand(rng, 30), w0)
        add!("matrix row", convolve_delay, view(M, 2, :), pmf(rng, 15))
        add!(
            "n = 220, D = 35", convolve_delay, rand(rng, 220) .+ 0.5,
            pmf(rng, 35); perf = true
        )

        add!("14 ⊕ 9", convolve_pmf, pmf(rng, 14), pmf(rng, 9))
        add!("1 ⊕ 1", convolve_pmf, [0.4], [0.7])
        add!("empty ⊕ 5", convolve_pmf, Float64[], pmf(rng, 5))
        add!("45 ⊕ 45", convolve_pmf, pmf(rng, 45), pmf(rng, 45); perf = true)

        ## Every knot layout the model builds: a regular weekly grid, a
        ## two-knot span, an irregular grid, and the single-knot window.
        for days in (knot_days(40), [1, 40], [5, 12, 19, 40], [1])
            add!(
                "knots $days", interpolate_knots, randn(rng, length(days)),
                days, 40
            )
        end
        δ = randn(rng, 3, length(knot_days(40)))
        add!("matrix row", interpolate_knots, view(δ, 2, :), knot_days(40), 40)
        add!(
            "n = 220", interpolate_knots, randn(rng, length(knot_days(220))),
            knot_days(220), 220; perf = true
        )

        Rt(n) = abs.(randn(rng, n)) .* 0.3 .+ 1.1
        add!(
            "G = 12", renewal_infections, Rt(40), pmf(rng, 12),
            rand(rng, 7) .+ 1
        )
        add!(
            "G > n", renewal_infections, Rt(10), pmf(rng, 15),
            rand(rng, 3) .+ 1
        )
        add!(
            "seed covers the grid", renewal_infections, Rt(6),
            pmf(rng, 4), rand(rng, 6) .+ 1
        )
        add!(
            "n = 220, G = 35", renewal_infections, Rt(220), pmf(rng, 35),
            rand(rng, 14) .+ 1; perf = true
        )

        ## No self-importation, and a daily importation intensity per patch.
        function patch_args(np, n, L, G)
            K = rand(rng, np, np) .* 0.2
            foreach(p -> K[p, p] = 0, 1:np)
            return (
                abs.(randn(rng, np, n)) .* 0.3 .+ 1.0, pmf(rng, G),
                rand(rng, np, L) .+ 1.0, K, rand(rng, np, n) .* 0.5,
            )
        end
        add!("3 patches", patch_infections, patch_args(3, 40, 7, 12)...)
        add!("1 patch", patch_infections, patch_args(1, 30, 5, 10)...)
        add!("G > n", patch_infections, patch_args(3, 10, 3, 15)...)
        add!(
            "seed covers the grid", patch_infections,
            patch_args(2, 5, 7, 4)...
        )
        add!(
            "3 patches, n = 220", patch_infections,
            patch_args(3, 220, 14, 35)...; perf = true
        )

        ## Floored inputs sit inside the flat region, away from the kink, so
        ## the finite-difference steps stay on one side: means below the
        ## `safe_rate` floor and a negative `k` on its fallback. A `k` large
        ## enough to clamp `p` is too large for a finite-difference step to
        ## move, so the guarded-inputs item checks that clamp instead.
        μ = exp.(4 .+ 0.5 .* randn(rng, 60))
        x = rand(rng, 0:120, 60)
        x[3] = 0
        for k in (8.3, 0.7, 150.0)
            add!("k = $k", nbinomial_loglik, k, μ, x)
        end
        add!(
            "floored means", nbinomial_loglik, 8.3, [μ; -1.0; -5.0],
            [x; 3; 0]
        )
        add!("fallback k", nbinomial_loglik, -1.0, μ, x)
        add!("empty", nbinomial_loglik, 8.3, Float64[], Int[])
        add!(
            "n = 220", nbinomial_loglik, 8.3,
            exp.(4 .+ 0.5 .* randn(rng, 220)), rand(rng, 0:120, 220);
            perf = true
        )

        add!("κ = 0.07", abscond_thinned, pmf(rng, 15), 0.07)
        add!("κ = 0", abscond_thinned, pmf(rng, 15), 0.0)
        add!("L = 45", abscond_thinned, pmf(rng, 45), 0.05; perf = true)

        ## Leading zero admissions skip those cohorts. The first two cases
        ## swap which schedule is longer. In the third the first schedule is
        ## longer than the series, so every cohort it admits is truncated.
        adm(n) = [zeros(3); abs.(randn(rng, n - 3)) .+ 0.5]
        for (l1, l2) in ((15, 9), (9, 15), (50, 12))
            add!(
                "schedules $l1 and $l2", abscond_thinned_flows, adm(40),
                pmf(rng, l1), adm(40), pmf(rng, l2), 0.07,
                0.05 .+ 0.3 .* rand(rng, 40)
            )
        end
        ## The forward pass skips a cohort with no admissions, so Mooncake's
        ## own derivation passes those admissions no derivative where the
        ## rule does. The timed case admits on every day, so the two agree.
        add!(
            "n = 220", abscond_thinned_flows, rand(rng, 220) .+ 0.5,
            pmf(rng, 35), rand(rng, 220) .+ 0.5, pmf(rng, 45), 0.05,
            0.05 .+ 0.3 .* rand(rng, 220);
            perf = true
        )

        ## Zero admissions, hazards of exactly 0 and 1, and a survival longer
        ## than the grid, which stops the cohort walk at day 1.
        survival(L) = sort(rand(rng, L); rev = true)
        A = abs.(randn(rng, 40)) .* 10 .+ 1
        A[[2, 9]] .= 0
        h = rand(rng, 40) .* 0.3
        h[[5, 6]] .= (0.0, 1.0)
        add!("n = 40, L = 15", two_clock_confirmed, A, h, survival(15))
        add!(
            "n < L", two_clock_confirmed, A[1:10], h[1:10], survival(15)
        )
        add!(
            "n = 220, L = 45", two_clock_confirmed,
            abs.(randn(rng, 220)) .* 10 .+ 1, rand(rng, 220) .* 0.3,
            survival(45); perf = true
        )

        ## Unequal PMF lengths, either way round, exercise the zero tail.
        for (Ld, Lr, cfr) in ((18, 12, 0.37), (10, 16, 0.0), (14, 14, 1.0))
            add!(
                "lengths $Ld and $Lr, cfr = $cfr", clinical_stay_survival,
                pmf(rng, Ld; mass = 0.95), pmf(rng, Lr; mass = 0.98), cfr
            )
        end
        add!(
            "lengths 45 and 45", clinical_stay_survival,
            pmf(rng, 45; mass = 0.95), pmf(rng, 45; mass = 0.98), 0.37;
            perf = true
        )

        for (note, kw) in (
                ("balanced", (;)),
                ("floored stocks", (; dmult = 1.6, κ = 0.08)),
                ("confirmed clamp", (; stop = 59, hflat = 0.9)),
            )
            add!(note, accumulate_occupancy, occupancy_args(rng, 60; kw...)...)
        end
        add!(
            "n = 220", accumulate_occupancy, occupancy_args(rng, 220)...;
            perf = true
        )

        gs, ge = 5, 40
        o = onset_args(rng; D = 12, gs, ge, n = 45)
        big = onset_args(rng; D = 28, gs = 120, ge = 220, n = 220)
        add!("D = 12", onset_report_cdf_table, o.lh, o.γ, gs, gs, ge)
        add!("D = 1", onset_report_cdf_table, o.lh[1:1], o.γ, gs, gs, ge)
        add!(
            "D = 28", onset_report_cdf_table, big.lh, big.γ, 120, 120, 220;
            perf = true
        )
        ## A daily anchor series and the length-1 constant default.
        add!("daily anchor", onset_report_anchor_series, o.tab, gs, o.alpha)
        add!("constant anchor", onset_report_anchor_series, o.tab, gs, [0.15])
        add!("D = 1", onset_report_anchor_series, o.tab[1:1, :], gs, o.alpha)
        add!(
            "D = 28", onset_report_anchor_series, big.tab, 120,
            abs.(randn(rng, 220)) .* 0.3; perf = true
        )
        add!(
            "D = 12", onset_report_moments, o.tab, gs, o.onsets, gs,
            o.alpha, o.oi, o.ri, o.pri
        )
        add!(
            "no cells", onset_report_moments, o.tab, gs, o.onsets, gs,
            o.alpha, Int[], Int[], Int[]
        )
        add!(
            "D = 28", onset_report_moments, big.tab, 120, big.onsets, 120,
            big.alpha, big.oi, big.ri, big.pri; perf = true
        )
        ## Onset dates before `grid_start` clamp both `γ` and `alpha`, and a
        ## cut-off before `grid_start` reads only clamped days.
        add!(
            "D = 12", onset_report_expected_total, o.onsets, o.lh, o.γ, gs,
            o.alpha, 45
        )
        add!(
            "D = 1", onset_report_expected_total, o.onsets, o.lh[1:1], o.γ,
            gs, o.alpha, 45
        )
        add!(
            "D = 0", onset_report_expected_total, o.onsets, Float64[], o.γ,
            gs, o.alpha, 45
        )
        add!(
            "cut-off before grid_start", onset_report_expected_total,
            o.onsets, o.lh, o.γ, gs, o.alpha, gs - 2
        )
        add!(
            "D = 28", onset_report_expected_total, big.onsets, big.lh,
            big.γ, 120, big.alpha, 220; perf = true
        )

        ## Increments are integers of either sign, as scanned, or floats, as
        ## simulated. Scales below the floor pass no derivative to the scale,
        ## and a negative `ν` takes the default. A non-finite scale is guarded
        ## too, but finite differences cannot step from it.
        m = 20 .* randn(rng, 60)
        σ = exp.(1 .+ 0.5 .* randn(rng, 60))
        y = round.(Int, m .+ 3 .* σ .* randn(rng, 60))
        for ν in (4.0, 1.5, 40.0)
            add!("ν = $ν", studentt_loglik, m, σ, y, ν)
        end
        add!(
            "float increments", studentt_loglik, m, σ,
            m .+ 3 .* σ .* randn(rng, 60), 4.0
        )
        add!(
            "floored scales", studentt_loglik, [m; -2.0; 1.0], [σ; -1.0; -3.0],
            [y; 1; -4], 4.0
        )
        add!("default ν", studentt_loglik, m, σ, y, -1.0)
        add!(
            "n = 220", studentt_loglik, 20 .* randn(rng, 220),
            exp.(1 .+ 0.5 .* randn(rng, 220)), rand(rng, -40:40, 220), 4.0;
            perf = true
        )

        ## A window with no trials, probabilities outside `[0, 1]` in the
        ## clamps' flat region, `ρ` below its floor and above its cap, and
        ## observations as a matrix row.
        n = rand(rng, 0:300, 30)
        p = 0.05 .+ 0.9 .* rand(rng, 30)
        k = [rand(rng, 0:n[i]) for i in 1:30]
        n[3] = 0
        k[3] = 0
        p_clamp = copy(p)
        p_clamp[[1, 2]] .= (-0.5, 1.5)
        k_clamp = copy(k)
        k_clamp[[1, 2]] .= (0, n[2])
        K = zeros(Int, 2, 30)
        K[1, :] .= k
        for (note, q, ρ, obs) in (
                ("ρ = 0.05", p, 0.05, k), ("ρ = 0.3", p, 0.3, k),
                ("matrix row", p, 0.8, view(K, 1, :)),
                ("clamped p", p_clamp, 0.05, k_clamp),
                ("ρ below its floor", p, -0.5, k),
                ("ρ above its cap", p, 1.5, k),
            )
            add!(note, betabinomial_loglik, n, q, ρ, obs)
        end
        nb = rand(rng, 0:300, 20)
        add!(
            "20 vintages", betabinomial_loglik, nb,
            0.05 .+ 0.9 .* rand(rng, 20), 0.05,
            [rand(rng, 0:t) for t in nb]; perf = true
        )

        ## The data-only helpers pass no derivative. Their inputs are the
        ## integer day indices and counts the histories carry, including a
        ## `missing` observation vector and a capacity history with no
        ## counts, as the predictive generator passes.
        add!("scored cells", onset_vintage_indices, o.ri, o.pri)
        add!(
            "production-sized", onset_vintage_indices, big.ri, big.pri;
            perf = true
        )
        days = collect(100:219)
        counts = rand(rng, 150:400, 120)
        cap = (; days = collect(100:3:219), counts = rand(rng, 300:500, 40))
        occ = (; days, counts)
        add!("120 days", censoring_cap, days, counts, cap; perf = true)
        add!("missing counts", censoring_cap, days, missing, cap)
        add!(
            "no recorded capacity", censoring_cap, days, counts,
            (; days = Int[], counts = Int[])
        )
        admitted = rand(rng, 0:40, 120)
        add!(
            "120 days", admission_headroom, days, admitted, cap, occ;
            perf = true
        )
        add!("missing counts", admission_headroom, days, missing, cap, occ)
        return cases
    end
end

@testitem "AD rules: every rule passes Mooncake's test_rule" tags = [
    :ad,
] setup = [ADRuleCases] begin
    using Mooncake: Mooncake, ReverseMode
    using Mooncake.TestUtils: test_rule
    using BVDOutbreakSize: abscond_thinned_flow, clinical_stay_survival,
        two_clock_confirmed

    rng = Xoshiro(20260923)
    for c in rule_cases(rng)
        @testset "$(c.name)" begin
            test_rule(
                rng, c.f, c.args...; is_primitive = true, mode = ReverseMode,
                rtol = 1.0e-6, atol = 1.0e-8
            )
        end
    end

    ## Two compositions the model calls through the rules: the single-flow
    ## wrapper, which passes the same arrays as both flows so the rule's
    ## tangents alias, and the treatment model's survival into the
    ## confirmed stock.
    adm = [zeros(3); abs.(randn(rng, 37)) .+ 0.5]
    pmf15 = pmf(rng, 15)
    h = 0.05 .+ 0.3 .* rand(rng, 40)
    test_rule(
        rng, abscond_thinned_flow, adm, pmf15, 0.07, h;
        is_primitive = false, mode = ReverseMode, rtol = 1.0e-6, atol = 1.0e-8
    )
    test_rule(
        rng,
        (a, b, d, r, c) -> two_clock_confirmed(
            a, b, clinical_stay_survival(d, r, c)
        ),
        abs.(randn(rng, 40)) .* 10 .+ 1, h, pmf(rng, 18; mass = 0.95),
        pmf(rng, 12; mass = 0.98), 0.37;
        is_primitive = false, mode = ReverseMode, rtol = 1.0e-6, atol = 1.0e-8
    )
end

@testitem "AD rules: each rule beats Mooncake's own derivation" tags = [
    :ad, :ad_perf,
] setup = [ADRuleCases] begin
    using Mooncake: Mooncake, DefaultCtx, MinimalCtx, Mode, ReverseMode,
        MooncakeInterpreter, build_rrule, get_interpreter, randn_tangent,
        zero_codual
    using Printf: @printf
    using BVDOutbreakSize: BVDOutbreakSize

    ## A context that sees every primitive the default one does except the
    ## rules `src/ad_rules.jl` registers, so Mooncake derives those kernels
    ## from their source as it would with the `ad_rules` preference off.
    struct NoPackageRulesCtx end
    function registered_here(M, sig)
        m = try
            which(
                Mooncake._is_primitive,
                Tuple{Type{MinimalCtx}, Type{M}, Type{sig}}
            )
        catch
            return false
        end
        return parentmodule(m) === BVDOutbreakSize
    end
    function Mooncake.is_primitive(
            ::Type{NoPackageRulesCtx}, M::Type{<:Mode}, sig, world::UInt
        )
        @nospecialize sig
        registered_here(M, sig) && return false
        return Mooncake.is_primitive(DefaultCtx, M, sig, world)
    end

    ## Value and argument tangents of one pullback. `__value_and_pullback!!`
    ## is the call a prepared `Mooncake.Cache` makes, here on either rule.
    pullback!(rule, ȳ, cx) = Mooncake.__value_and_pullback!!(rule, ȳ, cx...)

    ## Fastest time per call over batches of about 20 µs.
    function fastest(run; seconds = 0.3)
        run()
        evals = max(1, round(Int, 2.0e-5 / max(@elapsed(run()), 1.0e-9)))
        best = Inf
        stop = time() + seconds
        while time() < stop
            best = min(best, @elapsed(foreach(_ -> run(), 1:evals)) / evals)
        end
        return best
    end

    agrees(a::Union{Tuple, NamedTuple}, b) = all(map(agrees, a, b))
    agrees(a::Union{Real, AbstractArray{<:Real}}, b) = isapprox(
        a, b; rtol = 1.0e-8
    )
    agrees(a, b) = a == b

    ## Under coverage instrumentation the timings are not those of a fit, so
    ## the ratios are reported but not asserted there.
    instrumented = Base.JLOptions().code_coverage != 0
    rng = Xoshiro(20260923)
    for c in filter(c -> c.perf, rule_cases(rng))
        fx = (c.f, c.args...)
        sig = Tuple{map(Core.Typeof, fx)...}
        with = build_rrule(get_interpreter(ReverseMode), sig)
        without = build_rrule(
            MooncakeInterpreter(NoPackageRulesCtx, ReverseMode), sig
        )
        ȳ = randn_tangent(rng, c.f(c.args...))
        @testset "$(c.name)" begin
            ## The derived arm bypasses the rule and gives the same answer.
            @test with === Mooncake.rrule!!
            @test !(without isa typeof(Mooncake.rrule!!))
            @test agrees(
                pullback!(with, ȳ, map(zero_codual, fx)),
                pullback!(without, ȳ, map(zero_codual, fx))
            )
            cx = map(zero_codual, fx)
            t_rule = fastest(() -> pullback!(with, ȳ, cx))
            t_derived = fastest(() -> pullback!(without, ȳ, cx))
            ratio = t_rule / t_derived
            @printf(
                "%-40s rule %8.2f µs  derived %8.2f µs  ratio %.3f\n",
                c.name, 1.0e6 * t_rule, 1.0e6 * t_derived, ratio
            )
            instrumented || @test ratio <= 0.8
        end
    end
end

@testitem "AD rules: pullbacks leave the output tangent as given" tags = [
    :ad,
] begin
    using Random: Xoshiro
    using Mooncake: Mooncake, NoRData, tangent, zero_fcodual
    using BVDOutbreakSize: renewal_infections, patch_infections

    ## The renewal pullbacks walk the recursion backwards through a working
    ## copy of the output tangent rather than consuming it in place.
    rng = Xoshiro(20260923)
    out, pb = Mooncake.rrule!!(
        zero_fcodual(renewal_infections),
        map(zero_fcodual, (rand(rng, 40) .+ 1, rand(rng, 12), rand(rng, 7)))...
    )
    Ī = randn(rng, 40)
    tangent(out) .= Ī
    pb(NoRData())
    @test tangent(out) == Ī

    np, n = 3, 40
    K = rand(rng, np, np) .* 0.2
    K[[1, 5, 9]] .= 0
    args = (
        rand(rng, np, n) .+ 1, rand(rng, 12), rand(rng, np, 7) .+ 1, K,
        rand(rng, np, n) .* 0.5,
    )
    out, pb = Mooncake.rrule!!(
        zero_fcodual(patch_infections), map(zero_fcodual, args)...
    )
    Ī = randn(rng, np, n)
    tangent(out).infections .= Ī
    pb(NoRData())
    @test tangent(out).infections == Ī
end

@testitem "AD rules: guarded inputs pass no derivative" tags = [:ad] begin
    using Random: Xoshiro
    using Mooncake: Mooncake
    using BVDOutbreakSize: nbinomial_loglik, studentt_loglik

    function mgrad(f, args...)
        rule = Mooncake.build_rrule(f, args...)
        return Mooncake.value_and_gradient!!(rule, f, args...)[2][2:end]
    end

    rng = Xoshiro(20260923)
    μ = exp.(4 .+ 0.5 .* randn(rng, 60))
    x = rand(rng, 0:120, 60)
    ## `k` large enough that `p` clamps at `1 − eps` for every count.
    @test all(iszero, only(mgrad(b -> nbinomial_loglik(1.0e20, b, x), μ)))

    m = 20 .* randn(rng, 60)
    σ = exp.(1 .+ 0.5 .* randn(rng, 60))
    y = round.(Int, m .+ 3 .* σ .* randn(rng, 60))
    ## A scale at or below its floor, or infinite, passes no derivative to
    ## the scale but still passes one to the mean.
    m0, σ0, y0 = [m; 3.0; -2.0; 0.5], [σ; 0.0; -1.0; Inf], [y; 5; 1; 2]
    gm, gσ = mgrad((a, b) -> studentt_loglik(a, b, y0, 4.0), m0, σ0)
    @test all(iszero, gσ[61:63])
    @test all(!iszero, gm[61:63])
    ## A defaulted `ν` passes no derivative.
    @test iszero(only(mgrad(d -> studentt_loglik(m, σ, y, d), -1.0)))
end

@testitem "AD rules: the vector distributions score through the rules" tags = [
    :ad,
] begin
    using Random: Xoshiro
    using Mooncake: Mooncake
    using Distributions: logpdf, censored
    using BVDOutbreakSize: nbinomial_loglik, studentt_loglik,
        betabinomial_loglik, censored_nbinomial_loglik, NegBinomialVector,
        StudentTVector, BetaBinomialVector

    function mgrad(f, args...)
        rule = Mooncake.build_rrule(f, args...)
        return Mooncake.value_and_gradient!!(rule, f, args...)[2][2:end]
    end

    ## `logpdf` of each vector distribution is its summed helper, so the
    ## helper's rule fires through it and the gradient is the rule's own,
    ## bit for bit.
    rng = Xoshiro(20260923)
    μ = exp.(4 .+ 0.5 .* randn(rng, 60))
    x = rand(rng, 0:120, 60)
    @test mgrad((a, b) -> logpdf(NegBinomialVector(a, b), x), 8.3, μ) ==
        mgrad((a, b) -> nbinomial_loglik(a, b, x), 8.3, μ)
    ## Every seventh count sits at its ceiling and scores the tail.
    up = [i % 7 == 0 ? float(x[i]) : 1.0e6 for i in 1:60]
    @test mgrad(
        (a, b) -> logpdf(censored(NegBinomialVector(a, b); upper = up), x),
        8.3, μ
    ) == mgrad((a, b) -> censored_nbinomial_loglik(a, b, up, x), 8.3, μ)
    m = 20 .* randn(rng, 60)
    σ = exp.(1 .+ 0.5 .* randn(rng, 60))
    y = round.(Int, m .+ 3 .* σ .* randn(rng, 60))
    @test mgrad((a, b, d) -> logpdf(StudentTVector(a, b, d), y), m, σ, 4.0) ==
        mgrad((a, b, d) -> studentt_loglik(a, b, y, d), m, σ, 4.0)
    n = rand(rng, 0:300, 30)
    k = [rand(rng, 0:t) for t in n]
    p = 0.05 .+ 0.9 .* rand(rng, 30)
    @test mgrad((q, r) -> logpdf(BetaBinomialVector(n, q, r), k), p, 0.05) ==
        mgrad((q, r) -> betabinomial_loglik(n, q, r, k), p, 0.05)
end

@testitem "AD rules: the occupancy rule records the model's balance" tags = [
    :ad,
] setup = [ADRuleCases] begin
    using BVDOutbreakSize: _accumulate_occupancy, _OCC_CONF_HI

    ## The recording forward pass the rule runs gives the model's outputs
    ## bit for bit, and the clamp scenario reaches the confirmed cap.
    rng = Xoshiro(20260923)
    for kw in ((;), (; dmult = 1.6, κ = 0.08), (; stop = 59, hflat = 0.9))
        args = occupancy_args(rng, 60; kw...)
        y = accumulate_occupancy(args...)
        ỹ, _, flags = _accumulate_occupancy(Val(true), args...)
        @test all(k -> getfield(y, k) == getfield(ỹ, k), keys(y))
        haskey(kw, :hflat) && @test any(f -> f & _OCC_CONF_HI != 0, flags)
    end
end

@testitem "AD: the censored NegativeBinomial tail passes Mooncake's test_rule" tags = [
    :ad,
] begin
    using Random: Xoshiro
    using Mooncake: Mooncake
    using Mooncake.TestUtils: test_rule
    using Distributions: NegativeBinomial, Normal
    using Turing: @model, filldist, to_submodel, DynamicPPL
    using LogDensityProblems: logdensity_and_gradient
    using ADTypes: AutoMooncake
    using BVDOutbreakSize: nbinomial_logtail, censored_nbinomial_loglik,
        censored_occupancy_model

    rng = Xoshiro(20260923)
    ## Tails near one (demand far above the ceiling), moderate and far, one
    ## too small for a normal float, a dispersion below one, and a zero
    ## count, whose tail is one. A small step keeps `p` inside `(0, 1)`.
    for (r, μ, u) in (
            (8.3, 400.0, 40), (8.3, 50.0, 40), (8.3, 20.0, 60),
            (0.7, 100.0, 40), (30.0, 200.0, 120), (2.0, 5.0, 60),
            (200.0, 5.0, 300), (8.3, 50.0, 0),
        )
        test_rule(
            rng, nbinomial_logtail, NegativeBinomial(r, r / (r + μ)), u;
            is_primitive = false, mode = Mooncake.ReverseMode,
            max_fd_step = 1.0e-5
        )
    end

    ## The summed likelihood with counts at their ceilings. The ceilings are
    ## captured as integers, so finite differences never move a count across
    ## its ceiling.
    μ = exp.(4 .+ 0.5 .* randn(rng, 20))
    x = rand(rng, 0:120, 20)
    up = [i % 4 == 0 ? x[i] : 10^6 for i in 1:20]
    test_rule(
        rng, (a, b) -> censored_nbinomial_loglik(a, b, float.(up), x), 8.3, μ;
        is_primitive = false, mode = Mooncake.ReverseMode
    )

    ## Through the model: one count at its ceiling gives a finite gradient.
    @model function occupancy(obs, ceilings)
        lk ~ Normal(2, 0.5)
        lμ ~ filldist(Normal(3.5, 0.5), length(obs))
        x ~ to_submodel(
            censored_occupancy_model(exp.(lμ), ceilings, obs, exp(lk))
        )
    end
    model = occupancy([20, 40, 31], [60.0, 40.0, 50.0])
    ldf = DynamicPPL.LogDensityFunction(
        model, DynamicPPL.getlogjoint, DynamicPPL.VarInfo(model);
        adtype = AutoMooncake(; config = nothing)
    )
    lp, g = logdensity_and_gradient(ldf, [2.0, 3.4, 3.7, 3.5])
    @test isfinite(lp)
    @test all(isfinite, g) && !iszero(g[3])
end
