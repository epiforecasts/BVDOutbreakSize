## Tests for the weighted generation-interval prior. The weight raises the
## prior density to the power `w`, so a log-density difference between two
## points scales by `w`.

@testitem "generation_interval_model: the prior weight scales the density" begin
    using BVDOutbreakSize: generation_interval_model, cdf_nmax
    using Turing: DynamicPPL, condition, loglikelihood, @varname
    using Distributions: Gamma

    nmax = cdf_nmax(Gamma(2.71, 5.65))
    ## Conditioning puts the prior terms in the likelihood, the only terms
    ## the submodel has.
    lp(m, α, θ) = loglikelihood(
        condition(m, Dict(@varname(α) => α, @varname(θ) => θ)),
        DynamicPPL.VarInfo()
    )
    Δ(m) = lp(m, 1.4, 2.3) - lp(m, 2.71, 5.65)

    base = generation_interval_model(nmax)
    @test Δ(generation_interval_model(nmax; prior_weight = 1)) ≈ Δ(base)
    for w in (4, 24, 178)
        @test Δ(generation_interval_model(nmax; prior_weight = w)) ≈
            w * Δ(base)
    end
    weighted = generation_interval_model(nmax; prior_weight = 178)
    @test Set(keys(DynamicPPL.VarInfo(weighted))) ==
        Set([@varname(α), @varname(θ)])
end

@testitem "patch_infection_model: GI prior weighted by the renewal span" begin
    using BVDOutbreakSize: generation_interval_model, patch_infection_model,
        cdf_nmax
    using Turing: DynamicPPL, condition, loglikelihood, @varname
    using Distributions: Gamma
    using Random: Xoshiro

    n, n_patches, rt_start = 60, 4, 11
    nmax = cdf_nmax(Gamma(2.71, 5.65))
    m = patch_infection_model(n, n_patches; rt_start)
    vi = DynamicPPL.VarInfo(Xoshiro(3), m)
    lp(α, θ) = loglikelihood(
        condition(
            m, Dict(
                @varname(gi_state.α) => α, @varname(gi_state.θ) => θ
            )
        ), vi
    )
    gi = generation_interval_model(nmax)
    lp_gi(α, θ) = loglikelihood(
        condition(gi, Dict(@varname(α) => α, @varname(θ) => θ)),
        DynamicPPL.VarInfo()
    )
    @test lp(1.4, 2.3) - lp(2.71, 5.65) ≈
        (n - rt_start) * (lp_gi(1.4, 2.3) - lp_gi(2.71, 5.65))
end

@testitem "infection_model: GI prior weighted by the renewal span" begin
    using BVDOutbreakSize: generation_interval_model, infection_model,
        cdf_nmax
    using Turing: DynamicPPL, condition, loglikelihood, @varname
    using Distributions: Gamma
    using Random: Xoshiro

    n, rt_start = 60, 11
    nmax = cdf_nmax(Gamma(2.71, 5.65))
    m = infection_model(n; rt_start)
    vi = DynamicPPL.VarInfo(Xoshiro(3), m)
    lp(α, θ) = loglikelihood(
        condition(
            m, Dict(
                @varname(gi_state.α) => α, @varname(gi_state.θ) => θ
            )
        ), vi
    )
    gi = generation_interval_model(nmax)
    lp_gi(α, θ) = loglikelihood(
        condition(gi, Dict(@varname(α) => α, @varname(θ) => θ)),
        DynamicPPL.VarInfo()
    )
    @test lp(1.4, 2.3) - lp(2.71, 5.65) ≈
        (n - rt_start) * (lp_gi(1.4, 2.3) - lp_gi(2.71, 5.65))
end
