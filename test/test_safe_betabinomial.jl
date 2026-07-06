## Tests for `safe_betabinomial`, the guarded overdispersed Binomial used by
## the confirmed-positives windows (`src/models/observations.jl`). The
## confirmed positives carry extra-Binomial variation, so scoring them with a
## plain Binomial gives predictive intervals that are far too tight. The
## BetaBinomial inflates the variance while `ρ → 0` recovers the Binomial and
## the mean is unchanged.

@testitem "safe_betabinomial: mean matches the Binomial mean" begin
    using BVDOutbreakSize: safe_betabinomial
    using Distributions: BetaBinomial, mean

    ## Mean is `n·p` regardless of the overdispersion.
    for ρ in (1e-6, 0.01, 0.1, 0.3)
        d = safe_betabinomial(200, 0.3, ρ)
        @test d isa BetaBinomial
        @test mean(d) ≈ 200 * 0.3 rtol = 1e-6
    end
end

@testitem "safe_betabinomial: overdispersion inflates the variance" begin
    using BVDOutbreakSize: safe_betabinomial
    using Distributions: Binomial, var

    n, p = 200, 0.3
    binom_var = var(Binomial(n, p))
    ## `ρ → 0` collapses onto the Binomial variance.
    @test var(safe_betabinomial(n, p, 1e-6)) ≈ binom_var rtol = 1e-3
    ## Larger `ρ` strictly inflates the variance above the Binomial, and it
    ## grows monotonically in `ρ`.
    v_lo = var(safe_betabinomial(n, p, 0.02))
    v_hi = var(safe_betabinomial(n, p, 0.1))
    @test v_lo > binom_var
    @test v_hi > v_lo
    ## The inflation matches `1 + (n − 1)·ρ` for the intra-class correlation.
    @test var(safe_betabinomial(n, p, 0.05)) ≈
          binom_var * (1 + (n - 1) * 0.05) rtol = 1e-6
end

@testitem "safe_betabinomial: extreme inputs stay valid" begin
    using BVDOutbreakSize: safe_betabinomial
    using Distributions: BetaBinomial, params, mean

    ## Non-finite / out-of-range positivity and overdispersion are clamped so
    ## the distribution stays defined under extreme NUTS proposals.
    for (p, ρ) in ((NaN, 0.05), (1.5, 0.05), (-0.2, 0.05),
        (0.3, 0.0), (0.3, NaN), (0.3, 1.0), (0.3, -0.1))
        d = safe_betabinomial(100, p, ρ)
        @test d isa BetaBinomial
        α, β = params(d)[2], params(d)[3]
        @test α > 0 && β > 0
        @test isfinite(mean(d))
    end
end

@testitem "confirmed_overdispersion_model returns a unit-interval ρ" begin
    using BVDOutbreakSize: confirmed_overdispersion_model
    using Turing: returned
    using Random: MersenneTwister
    using Statistics: mean

    draws = [returned(confirmed_overdispersion_model(),
                 rand(MersenneTwister(i), confirmed_overdispersion_model())).ρ
             for i in 1:500]
    @test all(d -> 0 < d < 1, draws)
    ## Default `Beta(1, 24)` favours a small overdispersion (mean ≈ 0.04).
    @test 0.01 < mean(draws) < 0.1
end
