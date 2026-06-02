## Tests for the streaming-callback layer on `nuts_sample`
## (`src/sampling.jl`): the dependency-free `progress_callback`, the
## optional TensorBoard stub `tensorboard_callback`, and that
## `nuts_sample` forwards an arbitrary AbstractMCMC callback. These run
## without TensorBoardLogger installed, so the stub-error path is what CI
## exercises. Kept fast: trivial one-parameter Gaussian, few draws.

@testitem "nuts_sample forwards a callback to sample" tags=[:slow] begin
    using Distributions: Normal
    using Turing: @model
    using BVDOutbreakSize: nuts_sample

    @model function _cb_model()
        x ~ Normal(0.0, 1.0)
    end

    calls = Threads.Atomic{Int}(0)
    cb = (rng, model, sampler, transition, state, iteration;
        kwargs...) -> Threads.atomic_add!(calls, 1)

    nuts_sample(_cb_model(); samples = 20, chains = 1, callback = cb)
    ## One callback per post-warmup draw (warmup steps also fire), so the
    ## count is at least the number of retained draws.
    @test calls[] >= 20
end

@testitem "progress_callback writes the expected lines" tags=[:slow] begin
    using Distributions: Normal
    using Turing: @model
    using BVDOutbreakSize: nuts_sample, progress_callback

    @model function _pc_model()
        x ~ Normal(0.0, 1.0)
    end

    path = tempname()
    cb = progress_callback(; path = path, every = 5)
    nuts_sample(_pc_model(); samples = 50, chains = 1, callback = cb)

    @test isfile(path)
    lines = readlines(path)
    ## At least the retained draws (50) at one line per `every`=5 steps;
    ## warmup steps also fire, so allow a lower bound rather than equality.
    @test length(lines) >= 50 ÷ 5
    @test all(
        l -> occursin("iteration=", l) &&
             occursin("lp=", l) &&
             occursin("divergences=", l), lines)
    rm(path; force = true)
end

@testitem "tensorboard_callback errors without TensorBoardLogger" begin
    using BVDOutbreakSize: tensorboard_callback
    ## Mirrors the Enzyme stub pattern: the stub in src/sampling.jl
    ## throws until the extension is loaded. TensorBoardLogger is not a
    ## test dependency, so this is the active method here.
    @test_throws ErrorException tensorboard_callback("logs/run")
    try
        tensorboard_callback("logs/run")
    catch err
        @test occursin("TensorBoardLogger", err.msg)
    end
end
