## End-to-end smoke test for nuts_sample on a trivial one-parameter
## model. 50 draws × 2 chains is enough to exercise the wiring without
## blowing the CI budget. We accept whichever sample container Turing
## returns by default (a FlexiChains.VNChain here) and check shape +
## finite draws.

@testitem "nuts_sample returns a sample container with finite draws" tags=[
    :slow
] begin
    using Distributions: Normal
    using Turing: @model
    using BVDOutbreakSize: nuts_sample

    ## kept: a trivial one-parameter Gaussian is the cheapest target
    ## that still exercises every NUTS wiring path.
    @model function _nuts_model()
        x ~ Normal(0.0, 1.0)
    end

    chn = nuts_sample(_nuts_model(); samples = 50, chains = 2)
    @test chn !== nothing

    # Extract the draws of `x` from the FlexiChains.VNChain.
    xs = vec(Array(chn[:x]))
    @test length(xs) == 50 * 2
    @test all(isfinite, xs)
end
