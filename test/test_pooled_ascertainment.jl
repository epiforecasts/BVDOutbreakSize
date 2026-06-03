## Smoke tests for the per-stream DRC / Uganda ascertainment submodel.
## Exercises the real `pooled_ascertainment_model` from
## `src/models/priors.jl`, which samples each fraction independently from
## its own Beta prior (the name is retained for API compatibility).

@testsnippet PooledFixtures begin
    using Turing: @model, to_submodel
    using BVDOutbreakSize: pooled_ascertainment_model

    @model function _pooled_test_compose()
        asc ~ to_submodel(pooled_ascertainment_model(), false)
        p_drc_outer := asc.p_drc
        p_uganda_outer := asc.p_uganda
        return asc
    end
end

@testitem "ascertainment prior draws produce p ∈ (0, 1)" tags=[:slow] setup=[PooledFixtures] begin
    using Turing: sample, Prior
    import FlexiChains
    chn=sample(pooled_ascertainment_model(), Prior(), 200;
        chain_type = FlexiChains.VNChain, progress = false)
    p_drc=vec(Array(chn[:p_drc]))
    p_uganda=vec(Array(chn[:p_uganda]))
    @test length(p_drc) == 200
    @test length(p_uganda) == 200
    @test all(0 .< p_drc .< 1)
    @test all(0 .< p_uganda .< 1)
    @test all(isfinite, p_drc)
    @test all(isfinite, p_uganda)
end

@testitem "ascertainment composes via to_submodel" tags=[:slow] setup=[PooledFixtures] begin
    using Turing: sample, Prior
    import FlexiChains
    chn=sample(_pooled_test_compose(), Prior(), 100;
        chain_type = FlexiChains.VNChain, progress = false)
    p_drc=vec(Array(chn[:p_drc_outer]))
    p_uganda=vec(Array(chn[:p_uganda_outer]))
    @test length(p_drc) == 100
    @test all(0 .< p_drc .< 1)
    @test all(0 .< p_uganda .< 1)
end

@testitem "bvd_joint defaults to pooled ascertainment" tags=[:slow] begin
    using BVDOutbreakSize: bvd_joint, load_observations
    using Turing: sample, Prior, @varname
    import FlexiChains

    obs = load_observations()
    m = bvd_joint(obs.n, obs.exported_cases, obs.total_deaths,
        obs.reported_cases, obs.exports_deaths, obs.confirmed_cases,
        obs.tests_analysed;
        confirmed_deaths = obs.confirmed_deaths,
        deaths_history = obs.deaths_history,
        reported_history = obs.reported_history,
        confirmed_history = obs.confirmed_history,
        confirmed_deaths_history = obs.confirmed_deaths_history,
        lab_history = obs.lab_history,
        tests_received_history = obs.tests_received_history,
        breakpoint = obs.n - obs.who_first_sitrep_days,
        tmrca_days = obs.tmrca_days)
    chn = sample(m, Prior(), 50;
        chain_type = FlexiChains.VNChain, progress = false)
    ## The composer default is the non-centred two-group pooled
    ## hierarchy, so the shared hyperparameters `μ_logit` and the
    ## pooling SD `τ_logit` are sampled parameters of the joint model.
    params = FlexiChains.parameters(chn)
    @test @varname(τ_logit) in params
    @test @varname(μ_logit) in params
    @test all(vec(Array(chn[:τ_logit])) .>= 0)
    @test all(0 .< vec(Array(chn[:p_drc])) .< 1)
    @test all(0 .< vec(Array(chn[:p_uganda])) .< 1)
end
