## Tests for the opt-in Pathfinder NUTS initialisation (the weak-dependency
## extension `ext/BVDOutbreakSizePathfinderExt.jl`). Pathfinder is a test-only
## dependency here; `using Pathfinder` activates the extension.

@testitem "pathfinder_init returns one InitFromParams per chain" begin
    using BVDOutbreakSize: pathfinder_init, default_adtype
    using Pathfinder
    using Turing
    using Turing.DynamicPPL: InitFromParams
    using Random: MersenneTwister

    ## A vector parameter (`z`) mirrors the joint's non-centred blocks, so this
    ## checks `InitFromParams` is built for vector-valued variable names too.
    @model function _m()
        a ~ Normal(0, 1)
        z ~ filldist(Normal(0, 1), 3)
        s ~ truncated(Normal(0, 1); lower = 0)
        y ~ Normal(a + sum(z), s + 0.1)
    end

    inits = pathfinder_init(_m(), 3; nruns = 4, ndraws = 100,
        adtype = default_adtype(), rng = MersenneTwister(1))
    @test length(inits) == 3
    @test all(x -> x isa InitFromParams, inits)
end

@testitem "nuts_sample pathfinder-init runs and returns a chain" begin
    using BVDOutbreakSize: nuts_sample
    using Pathfinder
    using Turing
    import FlexiChains

    @model function _m()
        a ~ Normal(0, 1)
        z ~ filldist(Normal(0, 1), 3)
        s ~ truncated(Normal(0, 1); lower = 0)
        y ~ Normal(a + sum(z), s + 0.1)
    end

    ## Pathfinder-initialised fit returns the same chain type as a prior-init
    ## fit of the same model and settings.
    chn = nuts_sample(_m(); samples = 40, chains = 2, pathfinder = true,
        progress = false)
    @test chn isa FlexiChains.VNChain
    @test size(chn, 1) == 40

    chn_prior = nuts_sample(_m(); samples = 40, chains = 2, progress = false)
    @test chn_prior isa FlexiChains.VNChain
end
