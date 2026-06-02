## Tests for the streaming-callback layer on `nuts_sample`
## (`src/sampling.jl`): the dependency-free `progress_callback`, the
## optional TensorBoard `tensorboard_callback`, and that `nuts_sample`
## forwards an arbitrary AbstractMCMC callback. TensorBoardLogger is a
## test dependency (mirroring how Enzyme covers its extension), so the
## extension method in `ext/BVDOutbreakSizeTensorBoardLoggerExt.jl` is the
## active method here and is exercised directly. Kept fast: trivial
## one-parameter Gaussian, few draws.

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

@testitem "tensorboard_callback closure logs divergence diagnostics" begin
    using TensorBoardLogger: TBLogger
    using BVDOutbreakSize: tensorboard_callback
    ## Loading TensorBoardLogger activates the extension method (mirroring
    ## how `using Enzyme` activates `enzyme_adtype` in test_enzyme.jl), so
    ## this is the real method, not the src/sampling.jl error stub. Drive
    ## the closure directly with synthetic transitions to cover the
    ## divergence tally, the `every` gate and the scalar logging without a
    ## fit.
    ## A fresh, not-yet-existing path: `TBLogger` writes event files
    ## straight into `logdir` only when it does not already exist (an
    ## existing dir is bumped to a `_1` sibling), so build the path under
    ## a temp dir rather than passing `mktempdir()` itself.
    logdir = joinpath(mktempdir(), "run")
    cb = tensorboard_callback(logdir; every = 2)
    @test cb isa Function

    ## A NUTS-shaped transition: `lp` plus a `stat` exposing
    ## `numerical_error`, `step_size` and `tree_depth`.
    ok = (lp = -1.5,
        stat = (numerical_error = false, step_size = 0.3, tree_depth = 4))
    bad = (lp = -2.5,
        stat = (numerical_error = true, step_size = 0.1, tree_depth = 5))

    ## `every = 2`: iteration 1 is skipped for logging but still tallies a
    ## divergence; iterations 2 and 4 log. Never throws.
    @test cb(nothing, nothing, nothing, bad, nothing, 1) === nothing
    @test cb(nothing, nothing, nothing, ok, nothing, 2) === nothing
    @test cb(nothing, nothing, nothing, bad, nothing, 4) === nothing
    ## A malformed transition must be swallowed, not abort the fit.
    @test cb(nothing, nothing, nothing, (;), nothing, 6) === nothing

    files = readdir(logdir)
    @test any(startswith("events.out.tfevents"), files)
    @test any(f -> filesize(joinpath(logdir, f)) > 0, files)
end

@testitem "nuts_sample streams to TensorBoard via tensorboard_callback" tags=[:slow] begin
    using Distributions: Normal
    using Turing: @model
    using TensorBoardLogger: TBLogger
    using BVDOutbreakSize: nuts_sample, tensorboard_callback
    ## End-to-end: a tiny real fit driving the extension closure with
    ## genuine NUTS transitions, so the writer flushes event files.

    @model function _tb_model()
        x ~ Normal(0.0, 1.0)
    end

    logdir = joinpath(mktempdir(), "run")
    cb = tensorboard_callback(logdir; every = 1)
    nuts_sample(_tb_model(); samples = 20, chains = 1, callback = cb)

    files = readdir(logdir)
    @test any(startswith("events.out.tfevents"), files)
    @test any(f -> filesize(joinpath(logdir, f)) > 0, files)
end
