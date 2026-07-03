## Tests for the onset-to-sample delay constraint from the NEJM DRC 2026 BVD
## cohort (Akilimali et al. 2026, doi:10.1056/NEJMc2608070). The cohort delay
## is a Gamma whose mean is fixed to the reported value and whose SD is the
## single free parameter, inferred by fitting the Gamma median to the reported
## median; the delay is tied to the joint's onset-to-confirmation convolution
## by an N-weighted cross-entropy.

@testitem "nejm_onset_to_sample returns the delay configuration" begin
    using Distributions: Distribution, minimum
    using BVDOutbreakSize: nejm_onset_to_sample

    cfg = nejm_onset_to_sample()
    @test cfg.n_obs == 129
    @test cfg.mean_obs == 7.4
    @test cfg.median_obs == 4.8
    ## median_sd is the reported median 95% CrI half-width over 1.96.
    @test isapprox(cfg.median_sd, (7.84 - 3.46) / 2 / 1.96; atol = 1e-9)
    ## The SD is the free parameter, drawn from a positive, weakly-informative
    ## prior.
    @test cfg.sd_prior isa Distribution
    @test minimum(cfg.sd_prior) >= 0
end

@testitem "gamma_median_wh reproduces the median at the reported SD" begin
    using BVDOutbreakSize: gamma_median_wh

    ## At mean 7.4 and SD ~7.95 (the Gamma implied by the reported summaries)
    ## the Wilson-Hilferty median is within a few percent of the reported 4.8 d.
    @test isapprox(gamma_median_wh(7.4, 7.95), 4.8; atol = 0.2)
    ## Monotone: a larger SD (more right-skew) lowers the median below the mean.
    @test gamma_median_wh(7.4, 9.0) < gamma_median_wh(7.4, 6.0) < 7.4
end

@testitem "onset_to_sample_model returns a valid discretised PMF" begin
    using BVDOutbreakSize: onset_to_sample_model, nejm_onset_to_sample

    cfg = nejm_onset_to_sample()
    nmax = 60
    ## Evaluating the submodel draws the SD from its prior and builds the PMF.
    out = onset_to_sample_model(nmax; mean_obs = cfg.mean_obs,
        median_obs = cfg.median_obs, median_sd = cfg.median_sd,
        sd_prior = cfg.sd_prior)()

    @test length(out.pmf) == nmax + 1
    @test all(out.pmf .>= 0)
    @test isapprox(sum(out.pmf), 1.0; atol = 1e-8)
    @test out.mean == cfg.mean_obs
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

@testitem "onset_to_sample_model infers the SD from the median" tags=[:slow] begin
    using Turing: sample, NUTS
    using Statistics: mean
    using BVDOutbreakSize: onset_to_sample_model, nejm_onset_to_sample,
                           gamma_median_wh

    cfg = nejm_onset_to_sample()
    chn = sample(
        onset_to_sample_model(60; mean_obs = cfg.mean_obs,
            median_obs = cfg.median_obs, median_sd = cfg.median_sd,
            sd_prior = cfg.sd_prior),
        NUTS(200, 0.8), 400; progress = false)
    sd = vec(chn[:sd])
    ## The fitted delay reproduces the reported median (within the WH accuracy
    ## and the reported uncertainty), so the SD is inferred, not assigned.
    @test isapprox(mean(gamma_median_wh.(cfg.mean_obs, sd)), 4.8; atol = 0.4)
    @test mean(sd) > 5.0
end

@testitem "onset_to_sample carried by bvd_joint by default" tags=[:slow] begin
    using Turing: sample, Prior
    import FlexiChains
    using BVDOutbreakSize: bvd_joint

    ## Default-on: the joint constructs and samples with the latent onset-to-
    ## sample delay (its SD sampled and its median fitted) and the N-weighted
    ## cross-entropy tie.
    chn = sample(
        bvd_joint(40, 2, 18), Prior(), 50;
        chain_type = FlexiChains.VNChain, progress = false)
    r0 = vec(Array(chn[:R0]))
    @test length(r0) == 50
    @test all(isfinite, r0)

    ## Explicitly dropping the term with `nothing` also constructs and samples.
    chn0 = sample(
        bvd_joint(40, 2, 18; onset_to_sample = nothing), Prior(), 10;
        chain_type = FlexiChains.VNChain, progress = false)
    @test length(vec(Array(chn0[:R0]))) == 10
end
