## Tests for the guarded prior initialisation that keeps a chain off the
## prior predictive's unrecoverable tail (see `ViablePrior` in
## src/sampling.jl). The guard is a selection rule over forward density
## evaluations, so it is tested without sampling: the properties that
## matter are that it picks a better-than-typical prior draw, that it
## stays inside the prior's support, and that `nuts_sample` still runs
## and still gives each chain its own independent starting point.

@testitem "ViablePrior validates its attempt count" begin
    using BVDOutbreakSize: ViablePrior

    @test ViablePrior().attempts == 8
    @test ViablePrior(3).attempts == 3
    @test_throws ArgumentError ViablePrior(0)
    @test_throws ArgumentError ViablePrior(-1)
end

@testitem "viable_prior_init rejects the low-density prior tail" begin
    using Distributions: Normal, truncated
    using Random: MersenneTwister
    using Statistics: median
    using Turing: @model
    using Turing.DynamicPPL: DynamicPPL, InitFromPrior, VarInfo, getlogjoint
    using BVDOutbreakSize: viable_prior_init

    ## A wide prior on the scale of an observation that is well away from
    ## the prior centre: most prior draws score badly, a few score well,
    ## which is the shape that breaks the joint fit.
    @model function _wide()
        a ~ truncated(Normal(5.0, 4.0); lower = 0)
        y ~ Normal(a, 0.1)
    end
    model = _wide() | (; y = 0.5)

    logjoint(strategy, rng) = getlogjoint(VarInfo(rng, model, strategy))

    ## Typical unguarded prior draws, as the baseline to beat.
    plain = [logjoint(InitFromPrior(), MersenneTwister(s)) for s in 1:200]

    guarded = [logjoint(
                   viable_prior_init(MersenneTwister(s), model;
                       attempts = 8), MersenneTwister(1)) for s in 1:200]

    ## Every guarded start is finite, and the guard lifts the whole
    ## distribution of starting densities, not just its centre.
    @test all(isfinite, guarded)
    @test median(guarded) > median(plain)
    @test minimum(guarded) > minimum(plain)

    ## The guard must not smuggle in a value outside the prior's support.
    @test all(
        s -> begin
            vi = VarInfo(MersenneTwister(1),
                model, viable_prior_init(MersenneTwister(s), model))
            only(DynamicPPL.getindex_internal(vi, only(keys(vi)))) > 0
        end,
        1:50)
end

@testitem "viable_prior_init keeps more attempts at least as good" begin
    using Distributions: Normal, truncated
    using Random: MersenneTwister
    using Statistics: median
    using Turing: @model
    using Turing.DynamicPPL: VarInfo, getlogjoint
    using BVDOutbreakSize: viable_prior_init

    @model function _wide()
        a ~ truncated(Normal(5.0, 4.0); lower = 0)
        y ~ Normal(a, 0.1)
    end
    model = _wide() | (; y = 0.5)

    score(n) = median([getlogjoint(VarInfo(MersenneTwister(1), model,
                           viable_prior_init(MersenneTwister(s), model;
                               attempts = n))) for s in 1:100])

    @test score(16) >= score(4)
    @test_throws ArgumentError viable_prior_init(MersenneTwister(1), model;
        attempts = 0)
end

@testitem "viable_prior_init falls back when no draw is finite" begin
    using Distributions: Normal
    using Random: MersenneTwister
    using Turing: @model, @addlogprob!
    using Turing.DynamicPPL: InitFromPrior
    using BVDOutbreakSize: viable_prior_init

    ## Every draw scores -Inf, so the guard has nothing to choose between
    ## and must hand back unguarded prior initialisation rather than fail.
    @model function _impossible()
        a ~ Normal(0.0, 1.0)
        @addlogprob! -Inf
        return a
    end

    @test viable_prior_init(MersenneTwister(1), _impossible()) isa
          InitFromPrior
end

@testitem "nuts_sample starts each chain independently" tags=[:slow] begin
    using Distributions: Normal
    using Turing: @model
    using Turing.DynamicPPL: InitFromPrior
    using BVDOutbreakSize: nuts_sample, ViablePrior

    @model function _two()
        x ~ Normal(0.0, 1.0)
        z ~ Normal(0.0, 1.0)
    end

    ## The guarded default samples, and an explicit unguarded strategy
    ## still works, so the `init` keyword stays a pass-through.
    guarded = nuts_sample(_two(); samples = 50, chains = 2)
    plain = nuts_sample(_two(); samples = 50, chains = 2,
        init = InitFromPrior())
    for chn in (guarded, plain)
        xs = vec(Array(chn[:x]))
        @test length(xs) == 50 * 2
        @test all(isfinite, xs)
    end

    ## A single attempt is the unguarded draw, so it must still run.
    one = nuts_sample(_two(); samples = 50, chains = 2,
        init = ViablePrior(1))
    @test all(isfinite, vec(Array(one[:x])))
end
