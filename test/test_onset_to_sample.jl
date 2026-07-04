## Tests for the onset-to-sample delay constraint from the NEJM DRC 2026 BVD
## cohort (Akilimali et al. 2026, doi:10.1056/NEJMc2608070). The reported
## mean and median are continuous (epidist censoring-corrected), so the
## confirmed onset→report→receipt convolution is grounded on its own
## continuous mean (sum of the leg means) and median (Wilson-Hilferty of the
## summed variances) by soft Normal observations.

@testitem "nejm_onset_to_sample returns the cohort summaries and SEs" begin
    using BVDOutbreakSize: nejm_onset_to_sample

    cfg = nejm_onset_to_sample()
    @test cfg.mean_obs == 7.4
    @test cfg.median_obs == 4.8
    ## SEs are the reported 95% CrI half-widths over 1.96.
    @test isapprox(cfg.mean_se, (13.5 - 5.3) / 2 / 1.96; atol = 1e-9)
    @test isapprox(cfg.median_se, (7.84 - 3.46) / 2 / 1.96; atol = 1e-9)
end

@testitem "gamma_median_wh: continuous median under the mean" begin
    using BVDOutbreakSize: gamma_median_wh

    ## At mean 7.4 and SD ~7.95 (the Gamma matching the reported summaries) the
    ## Wilson-Hilferty continuous median is within a few percent of 4.8.
    @test isapprox(gamma_median_wh(7.4, 7.95), 4.8; atol = 0.2)
    ## Monotone: a larger SD (more right-skew) lowers the median below the mean.
    @test gamma_median_wh(7.4, 9.0) < gamma_median_wh(7.4, 6.0) < 7.4
end

@testitem "onset_to_sample_logweight peaks at the reported summaries" begin
    using BVDOutbreakSize: onset_to_sample_logweight, nejm_onset_to_sample,
                           gamma_median_wh

    cfg = nejm_onset_to_sample()
    ## Split the reported mean/median across two legs whose convolution
    ## reproduces them exactly: mean = μ_rep + μ_rec, and an SD whose WH median
    ## equals the reported median. Solve the SD from the WH relation.
    r = (cfg.median_obs / cfg.mean_obs)^(1 / 3)
    sd_star = 3 * cfg.mean_obs * sqrt(1 - r)          # WH-median SD at 4.8
    @test isapprox(gamma_median_wh(cfg.mean_obs, sd_star), cfg.median_obs;
        atol = 1e-6)
    ## Two legs summing to mean_obs and var sd_star^2.
    rep_m, rec_m = 3.0, cfg.mean_obs - 3.0
    rep_v = 0.4 * sd_star^2
    rec_v = sd_star^2 - rep_v
    at_target = onset_to_sample_logweight(rep_m, sqrt(rep_v), rec_m,
        sqrt(rec_v), cfg)
    ## Perturbing either leg away from the matched mean lowers the log-weight.
    off_mean = onset_to_sample_logweight(rep_m + 2.0, sqrt(rep_v), rec_m,
        sqrt(rec_v), cfg)
    off_spread = onset_to_sample_logweight(rep_m, sqrt(rep_v * 3), rec_m,
        sqrt(rec_v * 3), cfg)
    @test at_target > off_mean
    @test at_target > off_spread
    @test at_target ≈ 0 atol = 1e-6      # both residuals zero at the target
end

@testitem "onset_to_sample carried by bvd_joint by default" tags=[:slow] begin
    using Turing: sample, Prior
    import FlexiChains
    using BVDOutbreakSize: bvd_joint

    ## Default-on: the joint constructs and samples with the onset-to-sample
    ## mean/median soft constraint, which adds to the log-density without any
    ## new latent parameter. The modelled convolution mean is exposed.
    chn = sample(
        bvd_joint(40, 2, 18), Prior(), 50;
        chain_type = FlexiChains.VNChain, progress = false)
    r0 = vec(Array(chn[:R0]))
    @test length(r0) == 50
    @test all(isfinite, r0)
    @test all(isfinite, vec(Array(chn[:onset_to_sample_mean])))

    ## Explicitly dropping the term with `nothing` also constructs and samples.
    chn0 = sample(
        bvd_joint(40, 2, 18; onset_to_sample = nothing), Prior(), 10;
        chain_type = FlexiChains.VNChain, progress = false)
    @test length(vec(Array(chn0[:R0]))) == 10
end
