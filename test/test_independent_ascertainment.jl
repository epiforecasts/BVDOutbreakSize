## Smoke tests for the independent (unpooled) DRC / Uganda ascertainment
## submodel. Exercises the real `independent_ascertainment_model` from
## `src/models/priors.jl`, the composer default.

@testsnippet IndependentFixtures begin
    using Turing: @model, to_submodel
    using BVDOutbreakSize: independent_ascertainment_model

    @model function _independent_test_compose()
        asc ~ to_submodel(independent_ascertainment_model(), false)
        p_drc_outer := asc.p_drc
        p_uganda_outer := asc.p_uganda
        return asc
    end
end

@testitem "independent_ascertainment prior draws produce p ∈ (0, 1)" tags=[:slow] setup=[IndependentFixtures] begin
    using Turing: sample, Prior
    import FlexiChains
    chn=sample(independent_ascertainment_model(), Prior(), 200;
        chain_type = FlexiChains.VNChain, progress = false)
    p_drc=vec(Array(chn[:p_drc]))
    p_uganda=vec(Array(chn[:p_uganda]))
    @test length(p_drc) == 200
    @test length(p_uganda) == 200
    @test all(0 .< p_drc .< 1)
    @test all(0 .< p_uganda .< 1)
    @test all(isfinite, p_drc)
    @test all(isfinite, p_uganda)
    ## No shared hyperparameter: the two fractions carry their own
    ## logit-scale draws and are not tied through a pooling SD.
    @test !(:τ_logit in FlexiChains.parameters(chn))
    @test !(:μ_logit in FlexiChains.parameters(chn))
end

@testitem "independent_ascertainment composes via to_submodel" tags=[:slow] setup=[IndependentFixtures] begin
    using Turing: sample, Prior
    import FlexiChains
    chn=sample(_independent_test_compose(), Prior(), 100;
        chain_type = FlexiChains.VNChain, progress = false)
    p_drc=vec(Array(chn[:p_drc_outer]))
    p_uganda=vec(Array(chn[:p_uganda_outer]))
    @test length(p_drc) == 100
    @test all(0 .< p_drc .< 1)
    @test all(0 .< p_uganda .< 1)
end

@testitem "independent_ascertainment priors are separable" tags=[:slow] setup=[IndependentFixtures] begin
    using Turing: sample, Prior
    using Distributions: Normal
    using StatsFuns: logit
    import FlexiChains
    ## Distinct DRC / Uganda centres come through on their own fractions.
    chn=sample(
        independent_ascertainment_model(
            drc_prior = Normal(logit(0.05), 0.1),
            uganda_prior = Normal(logit(0.8), 0.1)),
        Prior(), 400; chain_type = FlexiChains.VNChain, progress = false)
    p_drc=vec(Array(chn[:p_drc]))
    p_uganda=vec(Array(chn[:p_uganda]))
    @test 0.02 < sum(p_drc) / length(p_drc) < 0.12
    @test 0.7 < sum(p_uganda) / length(p_uganda) < 0.9
end
