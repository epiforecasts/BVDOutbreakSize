## Smoke tests for the per-vintage reported-cases likelihood exercised
## through `cases_only_model` with a history `(; days, counts)`.

@testitem "reported_cases: prior draws are finite and non-negative" tags=[:slow] begin
    using Turing: sample, Prior
    import FlexiChains
    using BVDOutbreakSize: cases_only_model

    history = (; days = [13, 18, 23], counts = [340, 516, 905])
    chn = sample(
        cases_only_model(23, missing; reported_history = history),
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

@testitem "reported_cases: tiny fit with history and total stays positive" tags=[:slow] begin
    using Turing: sample, Prior
    import FlexiChains
    using BVDOutbreakSize: cases_only_model

    history = (; days = [13, 18, 23], counts = [340, 516, 905])
    chn = sample(
        cases_only_model(23, 905; reported_history = history),
        Prior(), 100;
        chain_type = FlexiChains.VNChain, progress = false
    )
    C_T = vec(Array(chn[:C_T]))
    @test length(C_T) == 100
    @test all(isfinite, C_T)
    @test all(C_T .> 0)
end
