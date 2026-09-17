## Smoke tests for the independent DRC / Uganda ascertainment submodel.
## Exercises the real `independent_ascertainment_model` from
## `src/models/priors.jl`, retained as the sensitivity alternative to the
## composer-default pooled hierarchy.

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

@testitem "independent_ascertainment prior draws produce p ∈ (0, 1)" tags=[
    :slow
] setup=[IndependentFixtures] begin
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
end

@testitem "independent_ascertainment composes via to_submodel" tags=[
    :slow
] setup=[IndependentFixtures] begin
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
