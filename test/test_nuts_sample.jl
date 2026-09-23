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

@testitem "nuts_sample builds plain NUTS when no term buffer is set" begin
    using Turing: NUTS
    using BVDOutbreakSize: default_adtype
    using BVDOutbreakSize: BVDOutbreakSize as B

    adtype = default_adtype()
    alg = B.nuts_algorithm(;
        n_adapts = 500, target_accept = 0.8, max_depth = 12, adtype,
        term_buffer = nothing
    )
    @test alg isa NUTS
    @test alg == NUTS(500, 0.8; max_depth = 12, adtype)
end

@testitem "a term buffer sets the terminal step-size window" begin
    using Turing: NUTS, Turing
    using BVDOutbreakSize: default_adtype
    using BVDOutbreakSize: BVDOutbreakSize as B
    TI = Turing.Inference
    AHMC = TI.AHMC

    adtype = default_adtype()
    metric = AHMC.DiagEuclideanMetric(3)
    alg(tb; n_adapts = 500) = B.nuts_algorithm(;
        n_adapts, target_accept = 0.7, max_depth = 12, adtype,
        term_buffer = tb
    )
    ref_alg = NUTS(500, 0.7; max_depth = 12, adtype)
    ref = TI.AHMCAdaptor(ref_alg, metric, 500; ϵ = 0.1)

    ## The default window reproduces the adaptor plain NUTS builds.
    a50 = TI.AHMCAdaptor(alg(50), metric, 500; ϵ = 0.1)
    @test a50 isa AHMC.StanHMCAdaptor
    @test ref.term_buffer == 50
    @test a50.term_buffer == 50
    @test a50.init_buffer == ref.init_buffer
    @test a50.window_size == ref.window_size
    @test a50.state.window_start == ref.state.window_start
    @test a50.state.window_end == ref.state.window_end
    @test a50.state.window_splits == ref.state.window_splits
    @test a50.ssa.state.ϵ == ref.ssa.state.ϵ
    @test a50.ssa.δ == ref.ssa.δ

    a200 = TI.AHMCAdaptor(alg(200), metric, 500; ϵ = 0.1)
    @test a200.term_buffer == 200
    @test a200.state.window_end == 300
    @test maximum(a200.state.window_splits) <= 300

    ## The trajectory and metric match plain NUTS.
    @test repr(TI.make_ahmc_kernel(alg(200), 0.1)) ==
        repr(TI.make_ahmc_kernel(ref_alg, 0.1))
    @test TI.getmetricT(alg(200)) == TI.getmetricT(ref_alg)

    ## No warmup means no adaptation, as for plain NUTS.
    @test TI.AHMCAdaptor(alg(200; n_adapts = 0), metric, 0; ϵ = 0.1) isa
        AHMC.Adaptation.NoAdaptation

    @test_throws ArgumentError alg(-1)
end

@testitem "nuts_sample runs with a term buffer" tags = [:slow] begin
    using Distributions: Normal
    using Turing: @model
    using BVDOutbreakSize: nuts_sample

    @model function _nuts_tb_model()
        x ~ Normal(0.0, 1.0)
    end

    chn = nuts_sample(
        _nuts_tb_model(); samples = 50, chains = 2, n_adapts = 40,
        term_buffer = 20
    )
    xs = vec(Array(chn[:x]))
    @test length(xs) == 50 * 2
    @test all(isfinite, xs)
end
