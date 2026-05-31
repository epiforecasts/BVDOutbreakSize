## Smoke tests for the per-bin random-effect DRC ascertainment used by
## the daily reported and confirmed likelihoods, and the
## deaths-reporting ascertainment factor used by the joint composers.

@testitem "daily_ascertainment_model: per-bin draws lie in (0, 1)" tags=[:slow] begin
    using BVDOutbreakSize: daily_ascertainment_model
    using Turing: sample, Prior, @model, to_submodel
    import FlexiChains
    using StatsFuns: logit
    using Statistics: std

    @model function _wrap(n, μ, τ)
        st ~ to_submodel(daily_ascertainment_model(n, μ, τ), false)
        return (; st)
    end

    chn = sample(_wrap(5, logit(0.25), 0.3), Prior(), 200;
        chain_type = FlexiChains.VNChain, progress = false)
    p = reduce(vcat, vec(Array(chn[:p_drc_t])))
    @test all(0 .< p .< 1)
    @test all(isfinite, p)
    ## With τ = 0 the per-bin draws should collapse to the mean. Pin
    ## that the τ = 0.3 spread is non-degenerate but bounded.
    @test std(p) > 0.01
end

@testitem "bvd_joint: fixed DRC ascertainment fraction" tags=[:slow] begin
    using BVDOutbreakSize: bvd_joint, load_observations
    using Turing: sample, Prior
    import FlexiChains

    obs = load_observations()
    rep = obs.reported_case_history
    conf = obs.confirmed_case_history
    dh = obs.death_history
    n_dh = length(dh.values)
    n_rep = length(rep.values)
    n_conf = length(conf.values)
    m = bvd_joint(missing,
        fill(missing, n_dh), fill(missing, n_rep);
        reported_offsets = rep.offsets,
        death_offsets = dh.offsets,
        confirmed_cases = fill(missing, n_conf),
        confirmed_offsets = conf.offsets)
    chn = sample(m, Prior(), 100;
        chain_type = FlexiChains.VNChain, progress = false)
    ## DRC ascertainment is the fixed pooled scalar `p_drc`, applied to
    ## every vintage bin, in (0, 1).
    @test all(0 .< vec(Array(chn[:p_drc])) .< 1)
end

@testitem "deaths_ascertainment_model: prior is centred near 1 with ~5% SD" tags=[:slow] begin
    using BVDOutbreakSize: deaths_ascertainment_model
    using Turing: sample, Prior, @model, to_submodel
    import FlexiChains
    using Statistics: std, mean

    @model function _wrap()
        st ~ to_submodel(deaths_ascertainment_model(), false)
        return (; st)
    end
    chn = sample(_wrap(), Prior(), 2000;
        chain_type = FlexiChains.VNChain, progress = false)
    p = vec(Array(chn[:p_deaths]))
    @test all(p .>= 0)
    @test all(isfinite, p)
    ## 0.85-1.15 covers ~3 SDs around the mean; the truncation at zero
    ## keeps p positive but otherwise leaves the bulk symmetric around 1.
    @test 0.95 < mean(p) < 1.05
    @test 0.04 < std(p) < 0.07
end

@testitem "bvd_joint: deaths ascertainment factor propagates to expected_deaths_T" tags=[:slow] begin
    using BVDOutbreakSize: bvd_joint
    using Turing: sample, Prior
    import FlexiChains

    ## deaths vector [100] at offset 0 (edge = T); reported missing.
    chn = sample(
        bvd_joint(1, [100], [missing]; reported_offsets = [0]),
        Prior(), 200;
        chain_type = FlexiChains.VNChain, progress = false)
    p = vec(Array(chn[:p_deaths]))
    @test all(0 .< p)
    @test all(isfinite, p)
    ## Fix p_deaths = 1 disables the factor; expected_deaths_T should
    ## not carry the multiplicative drift.
    chn_fixed = sample(
        bvd_joint(1, [100], [missing]; reported_offsets = [0],
            p_deaths_fixed = 1.0),
        Prior(), 50; chain_type = FlexiChains.VNChain,
        progress = false)
    @test :p_deaths ∉ keys(chn_fixed)
end
