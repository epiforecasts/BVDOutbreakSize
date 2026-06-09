## Smoke tests for predict_no_onward_deaths. A prior-predictive draw from
## the joint model provides the bare `:CFR`, `:C_T`, and
## `:expected_deaths_T` deterministics the function reads (the joint
## composer re-exposes these; the single-stream composers expose only
## `C_T`).

@testitem "predict_no_onward_deaths returns the documented columns" tags=[:slow] begin
    using DataFrames: DataFrame, nrow
    using Turing: sample, Prior
    import FlexiChains
    using BVDOutbreakSize: bvd_joint, predict_no_onward_deaths

    chn = sample(
        bvd_joint(40, missing, missing),
        Prior(), 100;
        chain_type = FlexiChains.VNChain, progress = false
    )
    obs_deaths = 18
    df = predict_no_onward_deaths(chn; obs_deaths = obs_deaths)

    @test df isa DataFrame
    @test sort(names(df)) == sort(["delta_deaths", "total_projected"])
    @test nrow(df) == 100

    @test all(isfinite, df.delta_deaths)
    @test all(isfinite, df.total_projected)
    @test all(df.delta_deaths .>= 0)
    @test all(df.total_projected .>= obs_deaths)
    @test maximum(
        abs.(df.total_projected .- (obs_deaths .+ df.delta_deaths))
    ) < 1e-8
end
