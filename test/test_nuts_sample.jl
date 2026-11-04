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
    using BVDOutbreakSize: BVDOutbreakSize, nuts_sample

    ## kept: a trivial one-parameter Gaussian is the cheapest target
    ## that exercises every NUTS wiring path.
    @model function _nuts_model()
        x ~ Normal(0.0, 1.0)
    end

    ## BLAS runs on one thread while the chains sample, and the caller's
    ## thread count comes back afterwards.
    BLAS = BVDOutbreakSize.BLAS
    before = BLAS.get_num_threads()
    during = Int[]
    record(args...; kwargs...) = push!(during, BLAS.get_num_threads())
    chn = nuts_sample(
        _nuts_model(); samples = 50, chains = 2, callback = record
    )
    @test chn !== nothing
    @test !isempty(during)
    @test all(==(1), during)
    @test BLAS.get_num_threads() == before

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
