## The package kernels written for speed against the textbook loops in
## `test/reference_kernels.jl`, on values. Gradients are `test_rule`'s job in
## `test/test_ad_rules.jl`.

@testitem "reference kernels: fast kernels match their textbook loops" setup = [
    ADRuleCases,
] begin
    using BVDOutbreakSize: renewal_infections_with_force,
        onset_report_expected_total, stick_breaking_loglik
    include(joinpath(@__DIR__, "reference_kernels.jl"))

    agrees(a::Union{Tuple, NamedTuple}, b) = all(map(agrees, a, b))
    agrees(a, b) = isapprox(a, b; rtol = 1.0e-12)

    rng = Xoshiro(20260923)
    cases = Any[]
    add!(note, f, ref, args...) = push!(cases, (; note, f, ref, args))

    for (n, D) in ((0, 5), (1, 1), (5, 1), (10, 3), (7, 12), (93, 40))
        add!(
            "convolve_delay n = $n, D = $D", convolve_delay,
            ref_convolve_delay, rand(rng, n), pmf(rng, D)
        )
    end
    w0 = pmf(rng, 10)
    w0[[1, 4]] .= 0
    add!(
        "convolve_delay with zero-weight lags", convolve_delay,
        ref_convolve_delay, rand(rng, 30), w0
    )
    add!(
        "convolve_delay of a zero series", convolve_delay,
        ref_convolve_delay, zeros(30), pmf(rng, 10)
    )
    add!(
        "convolve_delay on a matrix row", convolve_delay, ref_convolve_delay,
        view(rand(rng, 3, 40), 2, :), pmf(rng, 15)
    )
    for (na, nb) in ((14, 9), (1, 1), (45, 45), (0, 4))
        add!(
            "convolve_pmf $na ⊕ $nb", convolve_pmf, ref_convolve_pmf,
            pmf(rng, na), pmf(rng, nb)
        )
    end
    for (n, G, L) in ((40, 12, 7), (10, 15, 3), (6, 4, 6), (220, 35, 14))
        add!(
            "renewal n = $n, G = $G, L = $L", renewal_infections_with_force,
            ref_renewal_infections_with_force, rand(rng, n) .+ 0.8,
            pmf(rng, G), rand(rng, L) .+ 1
        )
    end
    ## Uncoupled, one patch, the seed covering the grid, and the production
    ## shape with a generation interval longer than the seed.
    for (np, n, L, G, scale) in (
            (3, 40, 7, 12, 0.0), (1, 30, 5, 10, 0.5), (2, 6, 6, 4, 0.5),
            (3, 220, 14, 35, 0.5),
        )
        K = rand(rng, np, np) .* 0.2
        foreach(p -> K[p, p] = 0, 1:np)
        add!(
            "patch_infections np = $np, n = $n, L = $L, G = $G, ε × $scale",
            patch_infections, ref_patch_infections,
            rand(rng, np, n) .+ 0.8, pmf(rng, G), rand(rng, np, L) .+ 1, K,
            scale .* rand(rng, np, n)
        )
    end
    for (note, kw) in (
            ("balanced", (;)), ("floored stocks", (; dmult = 1.6, κ = 0.08)),
            ("confirmed clamp", (; stop = 59, hflat = 0.9)),
            ("no admissions", (; stop = 0)),
        )
        add!(
            "accumulate_occupancy $note", accumulate_occupancy,
            ref_accumulate_occupancy, occupancy_args(rng, 60; kw...)...
        )
    end
    for n in (0, 1, 60)
        add!(
            "incare_census n = $n", incare_census, ref_incare_census,
            census_args(rng, n)...
        )
    end
    ## Hazards of 0 and 1 at the ends of the delay support, a one-day
    ## support, an empty onset range, and `γ` shorter than the table's reach.
    o = onset_args(rng; D = 12, gs = 5, ge = 40, n = 45)
    lh01 = [-800.0; o.lh[2:(end - 1)]; 800.0]
    for (note, lh, lo, hi) in (
            ("D = 12", o.lh, 5, 40), ("hazards 0 and 1", lh01, 5, 40),
            ("D = 1", o.lh[1:1], 5, 40), ("empty range", o.lh, 9, 8),
        )
        add!(
            "onset_report_cdf_table $note", onset_report_cdf_table,
            ref_onset_report_cdf_table, lh, o.γ, 5, lo, hi
        )
    end
    for (note, lh, as_of) in (
            ("D = 12", o.lh, 45), ("hazards 0 and 1", lh01, 45),
            ("cut-off before the grid", o.lh, 3), ("D = 1", o.lh[1:1], 45),
            ("cut-off inside the support", o.lh, 8),
        )
        add!(
            "onset_report_expected_total $note", onset_report_expected_total,
            ref_onset_report_expected_total, o.onsets, lh, o.γ, 5, o.alpha,
            as_of
        )
    end
    ## Positive and negative means, the sentinel and out-of-range previous
    ## scans, and no cells.
    for (note, lratio) in (("positive means", 0.6), ("negative means", 1.6))
        add!(
            "onset_scanned_cells $note", onset_scanned_cells,
            ref_onset_scanned_cells, scanned_args(rng, 40, 6; lratio)...
        )
    end
    add!(
        "onset_scanned_cells no cells", onset_scanned_cells,
        ref_onset_scanned_cells, Float64[], Float64[], [1.0], Int[], Int[],
        Int[], 2.1
    )
    μ = exp.(4 .+ 0.5 .* randn(rng, 60))
    x = rand(rng, 0:120, 60)
    for k in (0.7, 8.3, 150.0)
        add!(
            "nbinomial_loglik k = $k", nbinomial_loglik, ref_nbinomial_loglik,
            k, μ, x
        )
    end
    m = 20 .* randn(rng, 60)
    σ = exp.(1 .+ 0.5 .* randn(rng, 60))
    y = round.(Int, m .+ 3 .* σ .* randn(rng, 60))
    for ν in (1.5, 4.0, 40.0)
        add!(
            "studentt_loglik ν = $ν", studentt_loglik, ref_studentt_loglik,
            m, σ, y, ν
        )
    end
    n = rand(rng, 0:300, 30)
    n[3] = 0
    for ρ in (0.05, 0.3, 0.8)
        add!(
            "betabinomial_loglik ρ = $ρ", betabinomial_loglik,
            ref_betabinomial_loglik, n, 0.05 .+ 0.9 .* rand(rng, 30), ρ,
            [rand(rng, 0:t) for t in n]
        )
    end

    for (note, sizes) in (
            ("4 × 20", fill(4, 20)), ("ragged", [1, 2, 3, 4, 3]),
            ("one group of one", [1]), ("no rows", Int[]),
        )
        add!(
            "stick_breaking_loglik $note", stick_breaking_loglik,
            ref_stick_breaking_loglik, stick_args(rng, sizes; ρ = 0.3)...
        )
    end
    for c in cases
        @testset "$(c.note)" begin
            @test agrees(c.f(c.args...), c.ref(c.args...))
        end
    end

    ## The one intended difference. BLAS skips a lag whose weight is exactly
    ## zero, so the Inf on day 2 reaches only day 2 (lag 0) and day 4
    ## (lag 2), where the textbook loop gives `0 · Inf = NaN` on day 3.
    x_inf = [1.0, Inf, 2.0, 3.0]
    w_inf = [0.5, 0.0, 0.5]
    @test isnan(ref_convolve_delay(x_inf, w_inf)[3])
    @test isequal(convolve_delay(x_inf, w_inf), [0.5, Inf, 1.5, Inf])

    ## An empty delay support gives a zero total. With every hazard zero,
    ## a date whose reports are all in counts in full, where the share off
    ## the floored denominator would be zero, and the unsettled dates add
    ## nothing.
    @test onset_report_expected_total(o.onsets, Float64[], o.γ, 5, o.alpha, 45) ==
        0
    lh0 = fill(-800.0, 12)
    α(u) = o.alpha[clamp(u - 5 + 1, 1, length(o.alpha))]
    @test onset_report_expected_total(o.onsets, lh0, o.γ, 5, o.alpha, 45) ≈
        sum(o.onsets[u] * α(u) for u in 1:34) rtol = 1.0e-12
end
