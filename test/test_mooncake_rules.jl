## Correctness and speed gates for the hand-written reverse rules in
## `src/mooncake_rules.jl`.
##
## A wrong derivative is worse than a slow one, and a hand-written rule
## replaces the one place a backend would otherwise derive it for itself,
## so nothing else in the suite would notice it drifting. `ADRuleCases`
## lists every rule with the argument types its call sites pass and the
## edge cases its guards take. Mooncake's `test_rule` checks each case
## against finite differences with `is_primitive = true`, which also proves
## the rule fires for that signature. The speed item times one
## production-sized case per rule against a process that loads the package
## with the `mooncake_rules` preference off, and the joint item compares
## the production joint's gradient with that process.
##
## Tagged `:ad` with the rest of the gradient items, and the speed item
## also `:ad_perf`.

@testsnippet ADRuleCases begin
    using Random: Xoshiro
    using Distributions: _logpdf
    using BVDOutbreakSize: convolve_delay, convolve_pmf, interpolate_knots,
        knot_days, renewal_infections, patch_infections, NegBinomialVector,
        abscond_thinned, abscond_thinned_flows, two_clock_confirmed,
        clinical_stay_survival, accumulate_occupancy, incare_census,
        onset_report_cdf_table, onset_report_anchor_series,
        onset_report_moments, StudentTVector,
        BetaBinomialVector, censoring_cap, admission_headroom

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

    ## In-care census inputs cycling through five kinds of day, each clear of
    ## its tie: the two-clock stock inside `[0, O_bvd]`, above `O_bvd`, below
    ## zero, above the demand, and above a total a negative offset takes to
    ## zero.
    function census_args(rng, n; κ = 0.02)
        D = 50 .+ 100 .* rand(rng, n)
        O_bvd = 0.6 .* D
        x = O_bvd .* (0.2 .+ 0.6 .* rand(rng, n))
        Δ = zeros(n)
        for t in 1:n
            r = t % 5
            r == 1 && (x[t] = 1.3 * O_bvd[t])
            r == 2 && (x[t] = -5.0)
            r == 3 && (D[t] = 0.5 * x[t])
            r == 4 && (Δ[t] = -D[t])
        end
        return (D, O_bvd, x, κ, Δ)
    end

    ## Stick-breaking rows for groups of the given sizes: each group's shares
    ## sum to one, and its counts are a random split of a random total.
    function stick_args(rng, sizes; ρ = 0.05)
        groups = reduce(
            vcat, [fill(g, k) for (g, k) in enumerate(sizes)]; init = Int[]
        )
        shares = reduce(
            vcat, [(w = rand(rng, k) .+ 0.2; w ./ sum(w)) for k in sizes];
            init = Float64[]
        )
        counts = rand(rng, 0:60, length(groups))
        return (groups, counts, shares, ρ)
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
        ## A distribution's rule is named after the distribution.
        label(f, args) = f === _logpdf ? nameof(typeof(first(args))) :
            nameof(f)
        add!(note, f, args...; perf = false) = push!(
            cases, (; name = "$(label(f, args)): $note", f, args, perf)
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
            add!("k = $k", _logpdf, NegBinomialVector(k, μ), x)
        end
        add!(
            "floored means", _logpdf,
            NegBinomialVector(8.3, [μ; -1.0; -5.0]), [x; 3; 0]
        )
        add!("fallback k", _logpdf, NegBinomialVector(-1.0, μ), x)
        add!("empty", _logpdf, NegBinomialVector(8.3, Float64[]), Int[])
        add!(
            "n = 220", _logpdf,
            NegBinomialVector(8.3, exp.(4 .+ 0.5 .* randn(rng, 220))),
            rand(rng, 0:120, 220);
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
        add!("every branch", incare_census, census_args(rng, 60)...)
        add!("one day", incare_census, census_args(rng, 1)...)
        add!("no days", incare_census, census_args(rng, 0)...)
        add!(
            "n = 220", incare_census, census_args(rng, 220)...; perf = true
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

        ## Increments are integers of either sign, as scanned, or floats, as
        ## simulated. Scales below the floor pass no derivative to the scale,
        ## and a negative `ν` takes the default. A non-finite scale is guarded
        ## too, but finite differences cannot step from it.
        m = 20 .* randn(rng, 60)
        σ = exp.(1 .+ 0.5 .* randn(rng, 60))
        y = round.(Int, m .+ 3 .* σ .* randn(rng, 60))
        for ν in (4.0, 1.5, 40.0)
            add!("ν = $ν", _logpdf, StudentTVector(m, σ, ν), y)
        end
        add!(
            "float increments", _logpdf, StudentTVector(m, σ, 4.0),
            m .+ 3 .* σ .* randn(rng, 60)
        )
        add!(
            "floored scales", _logpdf,
            StudentTVector([m; -2.0; 1.0], [σ; -1.0; -3.0], 4.0), [y; 1; -4]
        )
        add!("default ν", _logpdf, StudentTVector(m, σ, -1.0), y)
        add!(
            "n = 220", _logpdf,
            StudentTVector(
                20 .* randn(rng, 220), exp.(1 .+ 0.5 .* randn(rng, 220)), 4.0
            ),
            rand(rng, -40:40, 220);
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
            add!(note, _logpdf, BetaBinomialVector(n, q, ρ), obs)
        end
        nb = rand(rng, 0:300, 20)
        add!(
            "20 vintages", _logpdf,
            BetaBinomialVector(nb, 0.05 .+ 0.9 .* rand(rng, 20), 0.05),
            [rand(rng, 0:t) for t in nb]; perf = true
        )

        ## The data-only helpers pass no derivative. Their inputs are the
        ## integer day indices and counts the histories carry, including a
        ## `missing` observation vector and a capacity history with no
        ## counts, as the predictive generator passes.
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

@testitem "Mooncake rules: every rule passes Mooncake's test_rule" tags = [
    :ad,
] setup = [ADRuleCases] begin
    using Mooncake: Mooncake, ReverseMode
    using Mooncake.TestUtils: test_rule
    using BVDOutbreakSize: abscond_thinned_flow, clinical_stay_survival,
        two_clock_confirmed, stick_breaking_loglik

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

    ## The stick-breaking composition, derived around the BetaBinomial rule:
    ## groups of one to four rows, a group with no counts, and each share's
    ## clamp and tail floor. In the last group the second share exceeds the
    ## tail the first leaves, so its conditional share clamps at one, the
    ## tail after it sits on its floor, and the next share clamps against
    ## that floor.
    g, y, sh, _ = stick_args(rng, [4, 4, 4])
    y[1:4] .= 0
    sh[9:12] .= (0.6, 0.7, 0.2, 0.1)
    for args in (
            stick_args(rng, [1, 2, 3, 4, 3]), (g, y, sh, 0.3),
            stick_args(rng, [4, 4]; ρ = -0.5), (Int[], Int[], Float64[], 0.05),
        )
        test_rule(
            rng, stick_breaking_loglik, args...;
            is_primitive = false, mode = ReverseMode, rtol = 1.0e-6,
            atol = 1.0e-8
        )
    end
end

@testsnippet RulesOff begin
    using Serialization: serialize, deserialize
    using BVDOutbreakSize: BVDOutbreakSize
    include(joinpath(@__DIR__, "rules_comparison.jl"))

    ## Runs `code` in a Julia process that loads the package with the
    ## `mooncake_rules` preference off, with `input` bound, and returns the
    ## value it leaves in `result`. The preference comes from a temporary
    ## environment at the end of the load path, so the checkout is untouched.
    function rules_off(code, input)
        dir = mktempdir()
        uuid = Base.PkgId(BVDOutbreakSize).uuid
        write(
            joinpath(dir, "Project.toml"),
            "[deps]\nBVDOutbreakSize = \"$uuid\"\n"
        )
        write(
            joinpath(dir, "LocalPreferences.toml"),
            "[BVDOutbreakSize]\nmooncake_rules = false\n"
        )
        inp, out = joinpath(dir, "in.jls"), joinpath(dir, "out.jls")
        serialize(inp, input)
        arm = joinpath(@__DIR__, "rules_comparison.jl")
        script = """
        using Serialization: serialize, deserialize
        include($(repr(arm)))
        input = deserialize($(repr(inp)))
        $code
        serialize($(repr(out)), (; rules_loaded = rules_loaded(), result))
        """
        sep = Sys.iswindows() ? ";" : ":"
        cmd = addenv(
            `$(Base.julia_cmd()) --threads=1 --project=$(Base.active_project()) -e $script`,
            "JULIA_LOAD_PATH" => join(["@", dir, "@stdlib"], sep)
        )
        run(cmd)
        return deserialize(out)
    end

    agrees(a::Union{Tuple, NamedTuple}, b; rtol) = all(
        map((x, y) -> agrees(x, y; rtol), a, b)
    )
    agrees(a::Mooncake.Tangent, b; rtol) = agrees(a.fields, b.fields; rtol)
    agrees(a::Real, b; rtol) = isapprox(a, b; rtol)
    agrees(a::AbstractArray{<:Real}, b; rtol) = isapprox(a, b; rtol)
    agrees(a::AbstractArray, b; rtol) = size(a) == size(b) &&
        all(agrees(x, y; rtol) for (x, y) in zip(a, b))
    ## Any other struct field by field, and anything else exactly.
    function agrees(a, b; rtol)
        isstructtype(typeof(a)) && fieldcount(typeof(a)) > 0 || return a == b
        return all(
            agrees(getfield(a, i), getfield(b, i); rtol)
                for i in 1:fieldcount(typeof(a))
        )
    end
end

@testitem "Mooncake rules: the joint matches the rules-off derivation" tags = [
    :ad,
] setup = [RulesOff] begin
    ## The production joint's log density and gradient with the rules and
    ## without them, at seeded prior points.
    model, ldf = production_ldf()
    xs = prior_points(model)
    on = joint_values(ldf, xs)
    off = rules_off(
        "result = joint_values(last(production_ldf()), input)", xs
    )
    @test rules_loaded()
    @test !off.rules_loaded
    for (a, b) in zip(on, off.result)
        @test isfinite(a.lp)
        @test isapprox(a.lp, b.lp; rtol = 1.0e-8)
        ## Each component, with a floor scaled to the largest for entries
        ## near zero.
        floor = 1.0e-10 * maximum(abs, b.g)
        @test all(@. abs(a.g - b.g) <= 1.0e-8 * abs(b.g) + floor)
    end
end

@testitem "Mooncake rules: each rule beats the rules-off derivation" tags = [
    :ad, :ad_perf,
] setup = [ADRuleCases, RulesOff] begin
    using Printf: @printf

    ## One production-sized case per rule, timed as a value and pullback
    ## with the rules loaded and in a process without them.
    cases = [
        (; c.name, c.f, c.args)
            for c in filter(c -> c.perf, rule_cases(Xoshiro(20260923)))
    ]
    on = pullback_times(cases)
    off = rules_off("result = pullback_times(input)", cases)
    @test !off.rules_loaded

    ## Under coverage instrumentation the timings are not those of a fit, so
    ## the ratios are reported but not asserted there.
    instrumented = Base.JLOptions().code_coverage != 0
    for (a, b) in zip(on, off.result)
        @testset "$(a.name)" begin
            @test agrees(a.result, b.result; rtol = 1.0e-8)
            ratio = a.time / b.time
            @printf(
                "%-40s rule %8.2f µs  derived %8.2f µs  ratio %.3f\n",
                a.name, 1.0e6 * a.time, 1.0e6 * b.time, ratio
            )
            instrumented || @test ratio <= 0.8
        end
    end
end

@testitem "Mooncake rules: pullbacks leave the output tangent as given" tags = [
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

@testitem "Mooncake rules: guarded inputs pass no derivative" tags = [:ad] begin
    using Random: Xoshiro
    using Mooncake: Mooncake
    using Distributions: logpdf
    using BVDOutbreakSize: NegBinomialVector, StudentTVector

    function mgrad(f, args...)
        rule = Mooncake.build_rrule(f, args...)
        return Mooncake.value_and_gradient!!(rule, f, args...)[2][2:end]
    end

    rng = Xoshiro(20260923)
    μ = exp.(4 .+ 0.5 .* randn(rng, 60))
    x = rand(rng, 0:120, 60)
    ## `k` large enough that `p` clamps at `1 − eps` for every count.
    @test all(
        iszero, only(mgrad(b -> logpdf(NegBinomialVector(1.0e20, b), x), μ))
    )

    m = 20 .* randn(rng, 60)
    σ = exp.(1 .+ 0.5 .* randn(rng, 60))
    y = round.(Int, m .+ 3 .* σ .* randn(rng, 60))
    ## A scale at or below its floor, or infinite, passes no derivative to
    ## the scale but still passes one to the mean.
    m0, σ0, y0 = [m; 3.0; -2.0; 0.5], [σ; 0.0; -1.0; Inf], [y; 5; 1; 2]
    gm, gσ = mgrad((a, b) -> logpdf(StudentTVector(a, b, 4.0), y0), m0, σ0)
    @test all(iszero, gσ[61:63])
    @test all(!iszero, gm[61:63])
    ## A defaulted `ν` passes no derivative.
    @test iszero(only(mgrad(d -> logpdf(StudentTVector(m, σ, d), y), -1.0)))
end

@testitem "Mooncake rules: the joint's likelihood rules fire" tags = [:ad] begin
    using Mooncake: Mooncake, MinimalCtx, ReverseMode
    using Distributions: _logpdf
    using BVDOutbreakSize: load_observations, joint_fit_args,
        default_breakpoint, incare_census, NegBinomialVector,
        StudentTVector, BetaBinomialVector

    ## The argument types the streams pass on the production data.
    obs = load_observations()
    h = joint_fit_args(
        obs; breakpoint = default_breakpoint(obs)
    ).onset_curve_history
    F = Vector{Float64}
    fires(sig) = Mooncake.is_primitive(
        MinimalCtx, ReverseMode, sig, Base.get_world_counter()
    )
    ## The onset cells score the history's increments.
    @test fires(
        Tuple{
            typeof(_logpdf), StudentTVector{F, F, Float64},
            typeof(h.increments),
        }
    )
    ## The treatment model's census takes its stocks and offset as float
    ## vectors. The count streams score integer increments about float
    ## means, and the composition's rows are the vectors
    ## `stick_breaking_loglik` builds.
    @test fires(Tuple{typeof(incare_census), F, F, F, Float64, F})
    @test fires(
        Tuple{
            typeof(_logpdf), NegBinomialVector{Float64, F}, Vector{Int},
        }
    )
    @test fires(
        Tuple{
            typeof(_logpdf), BetaBinomialVector{Vector{Int}, F, Float64},
            Vector{Int},
        }
    )
end

@testitem "Mooncake rules: a vector distribution scores as its entries" tags = [
    :ad,
] begin
    using Random: Xoshiro
    using Mooncake: Mooncake
    using Distributions: logpdf, censored
    using BVDOutbreakSize: safe_nbinomial, safe_studentt, safe_betabinomial,
        safe_rate, nbinomial_logtail, NegBinomialVector, StudentTVector,
        BetaBinomialVector, SplitCountVector

    function mgrad(f, args...)
        rule = Mooncake.build_rrule(f, args...)
        return Mooncake.value_and_gradient!!(rule, f, args...)
    end
    ## Value and gradient agree with the sum of the scalar terms, which
    ## Mooncake derives for itself, entry by entry.
    function agrees(a, b)
        return isapprox(a[1], b[1]; rtol = 1.0e-12) && all(
            map((x, y) -> isapprox(x, y; rtol = 1.0e-10), a[2][2:end], b[2][2:end])
        )
    end

    rng = Xoshiro(20260923)
    μ = exp.(4 .+ 0.5 .* randn(rng, 60))
    x = rand(rng, 0:120, 60)
    nb(k, m, i) = safe_nbinomial(k, safe_rate(m[i]))
    @test agrees(
        mgrad((a, b) -> logpdf(NegBinomialVector(a, b), x), 8.3, μ),
        mgrad((a, b) -> sum(logpdf(nb(a, b, i), x[i]) for i in 1:60), 8.3, μ)
    )
    ## Every seventh count sits at its ceiling and scores the tail.
    up = [i % 7 == 0 ? float(x[i]) : 1.0e6 for i in 1:60]
    function censored_terms(a, b)
        return sum(
            x[i] < up[i] ? logpdf(nb(a, b, i), x[i]) :
                nbinomial_logtail(nb(a, b, i), x[i]) for i in 1:60
        )
    end
    @test agrees(
        mgrad(
            (a, b) -> logpdf(censored(NegBinomialVector(a, b); upper = up), x),
            8.3, μ
        ),
        mgrad(censored_terms, 8.3, μ)
    )
    m = 20 .* randn(rng, 60)
    σ = exp.(1 .+ 0.5 .* randn(rng, 60))
    y = round.(Int, m .+ 3 .* σ .* randn(rng, 60))
    @test agrees(
        mgrad((a, b, d) -> logpdf(StudentTVector(a, b, d), y), m, σ, 4.0),
        mgrad(
            (a, b, d) -> sum(
                logpdf(safe_studentt(a[i], b[i], d), y[i]) for i in 1:60
            ), m, σ, 4.0
        )
    )
    n = rand(rng, 0:300, 30)
    k = [rand(rng, 0:t) for t in n]
    p = 0.05 .+ 0.9 .* rand(rng, 30)
    bb(q, r, i) = logpdf(safe_betabinomial(n[i], q[i], r), k[i])
    @test agrees(
        mgrad((q, r) -> logpdf(BetaBinomialVector(n, q, r), k), p, 0.05),
        mgrad((q, r) -> sum(bb(q, r, i) for i in 1:30), p, 0.05)
    )
    ## The split vector scores each entry by whether it has trials.
    n[[4, 9, 17]] .= 0
    μs = exp.(3 .+ 0.5 .* randn(rng, 30))
    function split_terms(q, r, s, m)
        return sum(
            n[i] > 0 ? bb(q, r, i) : logpdf(nb(s, m, i), k[i]) for i in 1:30
        )
    end
    @test agrees(
        mgrad(
            (q, r, s, m) -> logpdf(SplitCountVector(n, q, r, s, m), k),
            p, 0.05, 8.3, μs
        ),
        mgrad(split_terms, p, 0.05, 8.3, μs)
    )
end

@testitem "Mooncake rules: the occupancy rule records the model's balance" tags = [
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

@testitem "Mooncake rules: the census rule records the model's census" tags = [
    :ad,
] setup = [ADRuleCases] begin
    using BVDOutbreakSize: _incare_census, _CEN_CONF_X, _CEN_CONF_HI,
        _CEN_UNCONF, _CEN_SUSP

    ## The recording forward pass gives the census bit for bit, and the test
    ## inputs take every side of each clamp and floor.
    args = census_args(Xoshiro(20260923), 60)
    y = incare_census(args...)
    ỹ, _, flags = _incare_census(Val(true), args...)
    @test all(k -> getfield(y, k) == getfield(ỹ, k), keys(y))
    for bit in (_CEN_CONF_X, _CEN_CONF_HI, _CEN_UNCONF, _CEN_SUSP)
        @test any(f -> f & bit != 0, flags)
        @test any(f -> f & bit == 0, flags)
    end
end

@testitem "AD: the censored NegativeBinomial tail passes Mooncake's test_rule" tags = [
    :ad,
] begin
    using Random: Xoshiro
    using Mooncake: Mooncake
    using Mooncake.TestUtils: test_rule
    using Distributions: NegativeBinomial, Normal, logpdf, censored
    using Turing: @model, filldist, to_submodel, DynamicPPL
    using LogDensityProblems: logdensity_and_gradient
    using ADTypes: AutoMooncake
    using BVDOutbreakSize: nbinomial_logtail, NegBinomialVector,
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
        rng,
        (a, b) -> logpdf(
            censored(NegBinomialVector(a, b); upper = float.(up)), x
        ), 8.3, μ;
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

@testitem "reporting guard: true only where `:=` values are recorded" begin
    using BVDOutbreakSize: _reporting
    using Distributions: Normal
    using Turing: @model, DynamicPPL
    using .DynamicPPL: OnlyAccsVarInfo, ParamsWithStats, InitFromPrior,
        ldf_accs, getlogjoint_internal, @varname

    @model function probe()
        x ~ Normal()
        if _reporting(__varinfo__)
            y := x + 1
        end
        return _reporting(__varinfo__)
    end
    m = probe()

    ## A chain row records `y`; the gradient's accumulators and `m()` do not
    ## run the guarded branch.
    pws = ParamsWithStats(InitFromPrior(), m)
    @test haskey(pws.params, @varname(y))
    @test pws.params[@varname(y)] == pws.params[@varname(x)] + 1
    @test !_reporting(OnlyAccsVarInfo(ldf_accs(getlogjoint_internal)))
    @test m() === false
end
