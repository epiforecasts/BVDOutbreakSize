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
    for ρ in (1.0e-6, 0.01, 0.1, 0.3)
        d = safe_betabinomial(200, 0.3, ρ)
        @test d isa BetaBinomial
        @test mean(d) ≈ 200 * 0.3 rtol = 1.0e-6
    end
end

@testitem "safe_betabinomial: overdispersion inflates the variance" begin
    using BVDOutbreakSize: safe_betabinomial
    using Distributions: Binomial, var

    n, p = 200, 0.3
    binom_var = var(Binomial(n, p))
    ## `ρ → 0` collapses onto the Binomial variance.
    @test var(safe_betabinomial(n, p, 1.0e-6)) ≈ binom_var rtol = 1.0e-3
    ## Larger `ρ` strictly inflates the variance above the Binomial, and it
    ## grows monotonically in `ρ`.
    v_lo = var(safe_betabinomial(n, p, 0.02))
    v_hi = var(safe_betabinomial(n, p, 0.1))
    @test v_lo > binom_var
    @test v_hi > v_lo
    ## The inflation matches `1 + (n − 1)·ρ` for the intra-class
    ## correlation.
    @test var(safe_betabinomial(n, p, 0.05)) ≈
        binom_var * (1 + (n - 1) * 0.05) rtol = 1.0e-6
end

@testitem "safe_betabinomial: extreme inputs stay valid" begin
    using BVDOutbreakSize: safe_betabinomial
    using Distributions: BetaBinomial, params, mean

    ## Non-finite / out-of-range positivity and overdispersion are clamped so
    ## the distribution stays defined under extreme NUTS proposals.
    for (p, ρ) in (
            (NaN, 0.05), (1.5, 0.05), (-0.2, 0.05),
            (0.3, 0.0), (0.3, NaN), (0.3, 1.0), (0.3, -0.1),
        )
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

    draws = [
        returned(
            confirmed_overdispersion_model(),
            rand(MersenneTwister(i), confirmed_overdispersion_model())
        ).ρ
            for i in 1:500
    ]
    @test all(d -> 0 < d < 1, draws)
    ## Default `Beta(1, 24)` favours a small overdispersion (mean ≈ 0.04).
    @test 0.01 < mean(draws) < 0.1
end

@testitem "BetaBinomialVector matches one BetaBinomial per count" begin
    using BVDOutbreakSize: BetaBinomialVector, safe_betabinomial,
        province_composition_model
    using Distributions: logpdf, product_distribution
    using Random: Xoshiro
    using Turing: DynamicPPL, returned

    trials = [120, 0, 300, 0, 45]
    p = [0.2, 0.3, 0.05, 0.1, 0.6]
    x = [30, 0, 12, 0, 20]
    ρ = 0.03
    d = BetaBinomialVector(trials, p, ρ)
    product = product_distribution(
        [safe_betabinomial(trials[i], p[i], ρ) for i in 1:5]
    )
    @test length(d) == 5
    @test logpdf(d, x) ≈ logpdf(product, x)
    ## The same draws as the product it stands in for, from the same seed.
    @test rand(Xoshiro(7), d) == rand(Xoshiro(7), product)
    @test rand(Xoshiro(7), d, 3) == rand(Xoshiro(7), product, 3)

    ## Composition: each observed patch row is one summed term, equal to the
    ## row's `product_distribution` of BetaBinomials.
    rng = Xoshiro(3)
    modelled = exp.([3.0, 1.5, 0.5] .+ 0.3 .* randn(rng, 3, 8))
    obs = round.(Int, modelled)
    model = province_composition_model(obs, modelled)
    draw = rand(Xoshiro(2), model)
    state = returned(model, draw)
    function stick_breaking(shares, rho)
        remaining = vec(sum(obs; dims = 1))
        tail = ones(8)
        total = 0.0
        for q in 1:2
            p_cond = clamp.(shares[q, :] ./ tail, 0.0, 1.0)
            total += logpdf(
                product_distribution(
                    [
                        safe_betabinomial(max(remaining[i], 0), p_cond[i], rho)
                            for i in 1:8
                    ]
                ), obs[q, :]
            )
            remaining .-= obs[q, :]
            tail .= max.(tail .- shares[q, :], 1.0e-10)
        end
        return total
    end
    @test DynamicPPL.loglikelihood(model, draw) ≈
        stick_breaking(state.shares, state.rho)
    ## The predictive path samples each free row under one key and fills the
    ## last row with the remainder.
    predictive = province_composition_model(missing, modelled)
    rows = filter(
        v -> DynamicPPL.getsym(v) == :obs_increments,
        collect(keys(DynamicPPL.VarInfo(Xoshiro(1), predictive)))
    )
    @test length(rows) == 2
    drawn = returned(predictive, rand(Xoshiro(1), predictive)).obs_increments
    totals = [round(Int, max(sum(modelled[:, i]), 0.0)) for i in 1:8]
    @test vec(sum(drawn; dims = 1)) == totals
end
