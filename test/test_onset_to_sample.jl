## Tests for the fixed onset-to-sample delay constraint from the NEJM DRC
## 2026 BVD cohort (Akilimali et al. 2026, doi:10.1056/NEJMc2608070). The
## target is a fixed double-censored Gamma pinned by the cohort's reported
## mean and median; it is tied to the joint's onset-to-confirmation
## convolution by an N-weighted cross-entropy, with the cohort's uncertainty
## carried by the sample size N rather than any prior spread.

@testitem "nejm_onset_to_sample pins the Gamma from mean and median" begin
    using Distributions: Gamma, mean, median, std, quantile
    using BVDOutbreakSize: nejm_onset_to_sample

    cfg = nejm_onset_to_sample()
    @test cfg.n_obs == 129
    dist = Gamma(cfg.shape, cfg.scale)
    ## The target reproduces the two directly-reported pinning summaries.
    @test isapprox(mean(dist), 7.4; atol = 0.05)
    @test isapprox(median(dist), 4.8; atol = 0.1)
    ## The SD is a derived property near 8 d, not an input; it is consistent
    ## with the reported quartiles 1.81 / 10.23 d (sanity check).
    @test isapprox(std(dist), 7.95; atol = 0.6)
    @test isapprox(quantile(dist, 0.25), 1.81; atol = 0.7)
    @test isapprox(quantile(dist, 0.75), 10.23; atol = 1.0)
end

@testitem "onset_to_sample_model builds a valid fixed-target PMF" begin
    using BVDOutbreakSize: onset_to_sample_model, nejm_onset_to_sample

    cfg = nejm_onset_to_sample()
    nmax = 60
    out = onset_to_sample_model(nmax; shape = cfg.shape, scale = cfg.scale)

    @test length(out.pmf) == nmax + 1
    @test all(out.pmf .>= 0)
    @test isapprox(sum(out.pmf), 1.0; atol = 1e-8)
    ## The surfaced mean/SD are the fixed Gamma's moments, not sampled draws.
    @test isapprox(out.mean, 7.4; atol = 0.05)
    @test isapprox(out.sd, 7.95; atol = 0.6)
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

@testitem "onset_to_sample carried by bvd_joint by default" tags=[:slow] begin
    using Turing: sample, Prior
    import FlexiChains
    using BVDOutbreakSize: bvd_joint

    ## Default-on: the joint constructs and samples with the fixed-target
    ## cross-entropy term, which adds to the log-density without introducing
    ## any latent onset-to-sample parameter.
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
