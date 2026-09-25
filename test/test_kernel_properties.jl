## Properties the kernels with hand-written rules must satisfy. Their
## derivatives are `test_rule`'s job in `test/test_mooncake_rules.jl`.

@testsnippet KernelProperties begin
    using Random: Xoshiro
    ## A positive PMF of length `L` with total mass `mass`.
    pmf(rng, L; mass = 1.0) = (p = rand(rng, L) .+ 0.1; p .* (mass / sum(p)))
end

@testitem "convolve_delay: shift, linearity and mass" setup = [
    KernelProperties,
] begin
    using BVDOutbreakSize: convolve_delay
    rng = Xoshiro(1)
    x = rand(rng, 40)
    @test convolve_delay(x, [0.0, 1.0]) == [0.0; x[1:(end - 1)]]
    w, v = pmf(rng, 8), pmf(rng, 8)
    y = rand(rng, 40)
    @test convolve_delay(2x .+ 3y, w) ≈
        2convolve_delay(x, w) .+ 3convolve_delay(y, w)
    @test convolve_delay(x, 2w .+ 3v) ≈
        2convolve_delay(x, w) .+ 3convolve_delay(x, v)
    ## A series that is zero for the last `D - 1` days keeps all its mass.
    xz = [x[1:30]; zeros(10)]
    @test sum(convolve_delay(xz, w)) ≈ sum(xz)
    ## Day `t < D` holds only the lags that reach back to day 1.
    z = convolve_delay(x, w)
    @test all(t -> z[t] ≈ sum(w[d + 1] * x[t - d] for d in 0:(t - 1)), 1:7)
    ## A matrix row gives the same as its copy.
    M = rand(rng, 3, 40)
    @test convolve_delay(view(M, 2, :), w) ≈ convolve_delay(M[2, :], w)
    ## On floats a lag whose weight is exactly zero carries nothing, so an
    ## `Inf` on day 2 reaches only day 2 (lag 0) and day 4 (lag 2).
    @test isequal(
        convolve_delay([1.0, Inf, 2.0, 3.0], [0.5, 0.0, 0.5]),
        [0.5, Inf, 1.5, Inf]
    )
end

@testitem "convolve_pmf: mass, identity and symmetry" setup = [
    KernelProperties,
] begin
    using BVDOutbreakSize: convolve_pmf
    rng = Xoshiro(2)
    a, b = pmf(rng, 14; mass = 0.9), pmf(rng, 9; mass = 0.8)
    c = convolve_pmf(a, b)
    @test length(c) == 22
    @test sum(c) ≈ 0.9 * 0.8
    @test convolve_pmf(a, [1.0]) == a
    @test c ≈ convolve_pmf(b, a)
    @test isempty(convolve_pmf(Float64[], b))
end

@testitem "renewal_infections: seeds kept, geometric growth" setup = [
    KernelProperties,
] begin
    using BVDOutbreakSize: renewal_infections, renewal_infections_with_state
    seed = [1.0, 2.0, 4.0]
    ## A one-day generation interval multiplies each day by `R`.
    ## A pool so large the renewal is undepleted, exactly so as a power of two.
    I = renewal_infections(fill(1.5, 12), [1.0], seed, 2.0^900)
    @test I[1:3] == seed
    @test I[4:end] ≈ 4.0 .* 1.5 .^ (1:9)
    ## Each day after the seed is `R_t` times its force.
    rng = Xoshiro(3)
    R = rand(rng, 40) .+ 0.5
    st = renewal_infections_with_state(
        R, pmf(rng, 12), rand(rng, 7), 2.0^900
    )
    @test st.infections[8:end] == R[8:end] .* st.force[8:end]
end

@testitem "renewal_infections: depletion bounds new infections by the pool" setup = [
    KernelProperties,
] begin
    using BVDOutbreakSize: renewal_infections, renewal_infections_with_state
    rng = Xoshiro(5)
    g, seed = pmf(rng, 12), rand(rng, 7) .+ 1
    N = 500.0
    ## Each day takes `S (1 - e^{-x})` of the pool `S`, with `x = R f / N`,
    ## and leaves `S e^{-x}` for the next.
    R = rand(rng, 40) .+ 0.8
    st = renewal_infections_with_state(R, g, seed, N)
    x = R[8:40] .* st.force[8:40] ./ N
    pools = (N - sum(seed)) .* exp.(-cumsum(x))
    before = [N - sum(seed); pools[1:(end - 1)]]
    @test st.infections[8:40] ≈ before .* (1 .- exp.(-x))
    @test st.susceptible[8:40] ≈ pools
    ## New infections never exceed the pool left after the seed.
    I = renewal_infections(fill(3.0, 60), g, seed, N)
    @test all(>=(0), I)
    @test sum(I[8:end]) <= N - sum(seed) + 1.0e-8
    ## A proposal whose force overflows takes the whole pool rather than `Inf`.
    I = renewal_infections(fill(1.0e308, 30), g, seed, N)
    @test all(isfinite, I)
    @test sum(I[8:end]) ≈ N - sum(seed)
end

@testitem "pool_fraction: the pool left is the population less the cumulative" setup = [
    KernelProperties,
] begin
    using BVDOutbreakSize: renewal_infections_with_state, pool_fraction,
        adjusted_rt
    rng = Xoshiro(6)
    N = 800.0
    R = rand(rng, 50) .+ 1.0
    st = renewal_infections_with_state(R, pmf(rng, 10), rand(rng, 5) .+ 1, N)
    frac = pool_fraction(cumsum(st.infections), N)
    @test N .* frac[5:end] ≈ st.susceptible[5:end]
    ## Net of depletion, each day's reproduction number uses the pool the
    ## day before, and the first day the full pool.
    @test adjusted_rt(R, frac, 2:50) ≈ R[2:50] .* frac[1:49]
    @test adjusted_rt(R, frac, 1:1) == R[1:1]
end

@testitem "renewal_infections: depletion is light at census scale" setup = [
    KernelProperties,
] begin
    using BVDOutbreakSize: renewal_infections, PROVINCE_POPULATIONS
    using Distributions: Gamma
    using BVDOutbreakSize: discretise_censored
    ## Cumulative infections on the day the undepleted renewal first reaches
    ## `target`, with and without the pool `N`.
    function at_size(N, target)
        gi = discretise_censored(Gamma(2.71, 5.65), 40)
        g = gi[2:end] ./ sum(gi[2:end])
        R, seed = fill(2.0, 400), fill(10.0, length(g))
        free = cumsum(renewal_infections(R, g, seed, 2.0^900))
        t = findfirst(>=(target), free)
        return cumsum(renewal_infections(R, g, seed, N))[t] / free[t]
    end
    ## The outbreak's current size against the affected provinces barely
    ## moves, while a million-infection tail in Ituri is visibly trimmed.
    @test at_size(float(sum(PROVINCE_POPULATIONS)), 15_000) > 0.99
    @test at_size(float(PROVINCE_POPULATIONS[1]), 1.0e6) < 0.9
end

@testitem "abscond thinning: nests the plain convolution" setup = [
    KernelProperties,
] begin
    using BVDOutbreakSize: abscond_thinned, abscond_thinned_flows,
        convolve_delay
    rng = Xoshiro(4)
    adm, w, h = rand(rng, 40) .* 10, pmf(rng, 15), 0.3 .* rand(rng, 40)
    adm2, w2 = rand(rng, 40) .* 5, pmf(rng, 20)
    ## No absconding is the plain convolution; with no confirmation every
    ## cohort thins at `κ` per day of age.
    o1, o2 = abscond_thinned_flows(adm, w, adm2, w2, 0.0, h)
    @test o1 ≈ convolve_delay(adm, w)
    @test o2 ≈ convolve_delay(adm2, w2)
    @test first(abscond_thinned_flows(adm, w, adm2, w2, 0.1, zeros(40))) ≈
        convolve_delay(adm, abscond_thinned(w, 0.1))
    ## Absconding only removes, and a confirmation hazard slows it.
    t1, _ = abscond_thinned_flows(adm, w, adm2, w2, 0.1, h)
    @test all(>=(0), t1)
    @test all(t1 .<= o1 .+ 1.0e-12)
    @test all(
        t1 .>= first(abscond_thinned_flows(adm, w, adm2, w2, 0.1, zeros(40))) .-
            1.0e-12
    )
end

@testitem "two_clock_confirmed: the confirmation clock's limits" setup = [
    KernelProperties,
] begin
    using BVDOutbreakSize: two_clock_confirmed, convolve_delay
    rng = Xoshiro(5)
    A = rand(rng, 30) .* 10
    S = reverse(cumsum(pmf(rng, 12)))
    @test all(iszero, two_clock_confirmed(A, zeros(30), S))
    ## A hazard of one confirms every cohort the day after admission, so
    ## the confirmed stock is the whole stock less that day's admissions.
    @test two_clock_confirmed(A, ones(30), S) ≈
        convolve_delay(A, S) .- A .* S[1]
end

@testitem "incare_census: a four-day case" begin
    using BVDOutbreakSize: incare_census
    ## Day 3's confirmed stock exceeds the BVD stock and clamps to it, and
    ## day 4's is negative and clamps to zero.
    D = [10.0, 12.0, 8.0, 5.0]
    O_bvd = [6.0, 7.0, 4.0, 3.0]
    x = [2.0, 3.0, 5.0, -1.0]
    Δ = [0.0, 1.0, -9.0, 0.0]
    c = incare_census(D, O_bvd, x, 0.1, Δ)
    @test c.confirmed == [2.0, 3.0, 4.0, 0.0]
    @test c.abscond ≈ [0.0, 0.8, 0.9, 0.4]
    @test c.total == [10.0, 13.0, -1.0, 5.0]
    @test c.suspect == [8.0, 10.0, 0.0, 5.0]
end

@testitem "onset_report_cdf_table: hazard limits and monotone in delay" setup = [
    KernelProperties,
] begin
    using BVDOutbreakSize: onset_report_cdf_table
    rng = Xoshiro(6)
    γ = 0.2 .* randn(rng, 40)
    lh = randn(rng, 12) .- 1
    tab = onset_report_cdf_table(lh, γ, 5, 5, 40)
    @test size(tab) == (12, 36)
    @test all(>=(0), diff(tab; dims = 1))
    @test all(x -> 0 <= x <= 1, tab)
    @test all(iszero, onset_report_cdf_table(fill(-800.0, 12), γ, 5, 5, 40))
    @test all(isone, onset_report_cdf_table(fill(800.0, 12), γ, 5, 5, 40))
    @test size(onset_report_cdf_table(lh, γ, 5, 9, 8)) == (12, 0)
end

@testitem "stick_breaking_loglik: one conditional BetaBinomial per row" begin
    using BVDOutbreakSize: stick_breaking_loglik, safe_betabinomial
    using Distributions: logpdf
    ## Two groups: counts 5, 3, 2 of 10 at shares 0.5, 0.3, 0.2, and 4, 6
    ## of 10 at 0.4, 0.6. The last row of each group is its remainder.
    ρ = 0.1
    expected = logpdf(safe_betabinomial(10, 0.5, ρ), 5) +
        logpdf(safe_betabinomial(5, 0.3 / 0.5, ρ), 3) +
        logpdf(safe_betabinomial(10, 0.4, ρ), 4)
    @test stick_breaking_loglik(
        [1, 1, 1, 2, 2], [5, 3, 2, 4, 6], [0.5, 0.3, 0.2, 0.4, 0.6], ρ
    ) ≈ expected
    @test stick_breaking_loglik([1], [7], [1.0], ρ) == 0
end

@testitem "onset_report_expected_total: the sum of each date's reported share" setup = [
    KernelProperties,
] begin
    using BVDOutbreakSize: onset_report_expected_total, onset_report_F
    rng = Xoshiro(8)
    onsets = abs.(randn(rng, 45)) .* 20
    γ = 0.3 .* randn(rng, 36)
    alpha = abs.(randn(rng, 36)) .* 0.3
    α(u) = alpha[clamp(u - 5 + 1, 1, 36)]
    ## An empty delay support reports nothing.
    @test onset_report_expected_total(onsets, Float64[], γ, 5, alpha, 45) == 0
    ## Each date contributes its onsets times its reported share, settled
    ## dates included, also when the report probability sits on the
    ## `safe_rate` floor.
    for lh in (randn(rng, 12) .- 2, fill(-40.0, 12), fill(-800.0, 12))
        by_date = sum(
            onsets[u] * onset_report_F(45 - u, lh, γ, u, 5, α(u)) for u in 1:45
        )
        @test onset_report_expected_total(onsets, lh, γ, 5, alpha, 45) ≈
            by_date rtol = 1.0e-12
    end
end
