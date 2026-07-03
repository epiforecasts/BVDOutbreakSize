## Tests for the latent onset-to-sample delay constraint from the NEJM DRC
## 2026 BVD cohort (Akilimali et al. 2026, doi:10.1056/NEJMc2608070). The
## submodel is a double-censored Gamma whose mean/SD priors carry the SI
## Table S3 confirmed-positive fit and its uncertainty; the constraint ties
## it to the joint's onset-to-confirmation convolution.

@testitem "onset_to_sample_model: double-censored Gamma PMF is valid" begin
    using Distributions: truncated, Normal
    using BVDOutbreakSize: onset_to_sample_model
    using Turing: returned
    using Random: MersenneTwister

    nmax = 60
    m = onset_to_sample_model(nmax;
        mean_prior = truncated(Normal(7.39, 2.09); lower = 1),
        sd_prior = truncated(Normal(7.95, 2.0); lower = 1))
    out = returned(m, rand(MersenneTwister(1), m))

    @test length(out.pmf) == nmax + 1
    @test all(out.pmf .>= 0)
    @test isapprox(sum(out.pmf), 1.0; atol = 1e-8)
    ## The sampled delay mean/SD are surfaced and positive.
    @test out.mean > 0
    @test out.sd > 0
end

@testitem "delay_match_logweight rewards a matching convolution" begin
    using BVDOutbreakSize: delay_match_logweight

    target = [0.1, 0.3, 0.4, 0.2]
    close = [0.12, 0.28, 0.4, 0.2]
    far = [0.4, 0.3, 0.2, 0.1]

    ## Closer modelled PMF scores higher, and the exact match is the best.
    @test delay_match_logweight(target, close, 10) >
          delay_match_logweight(target, far, 10)
    @test delay_match_logweight(target, target, 1) >
          delay_match_logweight(target, close, 1)
    ## Weight scales the term linearly (it is N pseudo-observations).
    @test isapprox(delay_match_logweight(target, close, 20),
        2 * delay_match_logweight(target, close, 10); rtol = 1e-12)
    ## Tolerates differing lengths by matching the shared support.
    @test isfinite(delay_match_logweight(target, [0.1, 0.3, 0.6], 5))
end

@testitem "nejm_onset_to_sample carries the SI Table S3 confirmed fit" begin
    using Distributions: mean
    using BVDOutbreakSize: nejm_onset_to_sample

    cfg = nejm_onset_to_sample()
    @test cfg.n_obs == 129
    ## Priors centred on the recovered confirmed-positive Gamma moments.
    @test isapprox(mean(cfg.mean_prior), 7.39; atol = 0.2)
    @test isapprox(mean(cfg.sd_prior), 7.95; atol = 0.3)
end

@testitem "onset_to_sample composes into bvd_joint via kwarg" tags=[:slow] begin
    using Turing: sample, Prior
    import FlexiChains
    using BVDOutbreakSize: bvd_joint, nejm_onset_to_sample

    chn = sample(
        bvd_joint(40, 2, 18; onset_to_sample = nejm_onset_to_sample()),
        Prior(), 50;
        chain_type = FlexiChains.VNChain, progress = false)

    mean_draws = vec(Array(chn[:onset_to_sample_mean]))
    @test length(mean_draws) == 50
    @test all(isfinite, mean_draws)
    @test all(mean_draws .> 0)
end
