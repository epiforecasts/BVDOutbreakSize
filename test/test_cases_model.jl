## Smoke tests for the reported-cases stream composer.
## Exercises the real `cases_only_model` from `src/models/joint.jl`.

@testitem "cases_only_model prior draws are finite and non-negative" tags=[:slow] begin
    using Turing: sample, Prior
    import FlexiChains
    using BVDOutbreakSize: cases_only_model

    chn = sample(
        cases_only_model(40, missing),
        Prior(), 100;
        chain_type = FlexiChains.VNChain, progress = false
    )

    C_T = vec(Array(chn[:C_T]))
    @test length(C_T) == 100
    @test all(isfinite, C_T)
    @test all(C_T .> 0)

    ## The DRC ascertainment fraction is internal to the cases composer
    ## and surfaces under the ascertainment submodel prefix.
    p_drc = vec(Array(chn[Symbol("asc_state.p_drc")]))
    @test all(0 .< p_drc .< 1)
end

@testitem "cases_only_model conditioned on an observation stays positive" tags=[:slow] begin
    using Turing: sample, Prior
    import FlexiChains
    using BVDOutbreakSize: cases_only_model

    chn = sample(
        cases_only_model(40, 50),
        Prior(), 100;
        chain_type = FlexiChains.VNChain, progress = false
    )
    C_T = vec(Array(chn[:C_T]))
    @test length(C_T) == 100
    @test all(isfinite, C_T)
    @test all(C_T .> 0)
end
