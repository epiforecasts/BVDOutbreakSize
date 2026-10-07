## End-to-end smoke test for nuts_sample on a trivial one-parameter
## model. 50 draws × 2 chains is enough to exercise the wiring without
## blowing the CI budget. We accept whichever sample container Turing
## returns by default (a FlexiChains.VNChain here) and check shape +
## finite draws.

@testitem "nuts_sample returns a sample container with finite draws" tags = [
    :slow,
] begin
    using Distributions: Normal
    using Turing: @model
    using Turing.DynamicPPL: InitFromPrior, InitFromUniform
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

    ## One initialisation strategy per chain, and the refusal when the
    ## vector and the chain count disagree.
    per_chain = nuts_sample(
        _nuts_model(); samples = 50, chains = 2,
        init = [InitFromPrior(), InitFromUniform()]
    )
    @test length(vec(Array(per_chain[:x]))) == 50 * 2
    @test_throws ArgumentError nuts_sample(
        _nuts_model(); samples = 50, chains = 2,
        init = [InitFromPrior()]
    )
end

@testitem "nuts_sample adapts a dense mass matrix when asked" tags = [
    :slow,
] begin
    using Distributions: MvNormal
    using Turing: @model
    using BVDOutbreakSize: nuts_sample

    ## Two strongly correlated parameters: the case a dense metric is for.
    @model function _correlated_model()
        x ~ MvNormal([0.0, 0.0], [1.0 0.95; 0.95 1.0])
    end

    chn = nuts_sample(
        _correlated_model(); samples = 100, chains = 2, metric = :dense
    )
    xs = Array(chn[:x])
    @test all(isfinite, reduce(vcat, vec(xs)))
    @test_throws ArgumentError nuts_sample(
        _correlated_model(); samples = 10, chains = 1, metric = :unit
    )
end
