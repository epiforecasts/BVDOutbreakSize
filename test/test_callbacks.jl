## Tests for the streaming-callback layer on `nuts_sample`
## (`src/sampling.jl`): the dependency-free `progress_callback`, the
## optional TensorBoard `tensorboard_callback` (a thin wrapper over
## `AbstractMCMC.mcmc_callback`), and that `nuts_sample` forwards an
## arbitrary AbstractMCMC callback. Both streaming callbacks read step
## statistics through the `AbstractMCMC.ParamsWithStats` interface rather
## than poking transition fields, so they are exercised against genuine
## Turing transitions from a real (tiny) fit. TensorBoardLogger is a test
## dependency (mirroring how Enzyme covers its extension). Kept fast:
## trivial one-parameter Gaussian, few draws.

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
    ## One callback per post-warmup draw, so the count is at least the
    ## number of retained draws.
    @test calls[] >= 20
end

@testitem "progress_callback records real log-density" tags=[:slow] begin
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
    ## At least the retained draws (50) at one line per `every`=5 steps.
    @test length(lines) >= 50 ÷ 5
    @test all(
        l -> occursin("iteration=", l) &&
             occursin("lp=", l) &&
             occursin("divergences=", l), lines)
    ## The log-density must be a real number pulled from the transition
    ## statistics, not `missing` — this is what regresses if the stats
    ## interface stops exposing `logjoint` (the bug the old hand-rolled
    ## `transition.lp` access hit silently).
    lps = [match(r"lp=(\S+)", l).captures[1] for l in lines]
    @test any(s -> tryparse(Float64, s) !== nothing, lps)
    @test !all(==("missing"), lps)
    rm(path; force = true)
end

@testitem "tensorboard_callback streams grouped scalars and histograms" tags=[:slow] begin
    using Distributions: Normal
    using Turing: @model
    using TensorBoardLogger: TBLogger, map_summaries
    using BVDOutbreakSize: nuts_sample, tensorboard_callback
    ## Loading TensorBoardLogger activates the BVDOutbreakSize stub's real
    ## method (mirroring how `using Enzyme` activates `enzyme_adtype`).
    ## End-to-end: a tiny real fit drives the callback with genuine NUTS
    ## transitions, then the event files are read back to confirm the
    ## grouped tag layout and that histograms are logged alongside scalars.

    @model function _tb_model()
        x ~ Normal(0.0, 1.0)
    end

    cb = tensorboard_callback(joinpath(mktempdir(), "construct"))
    @test cb !== nothing

    logdir = joinpath(mktempdir(), "run")
    fit_cb = tensorboard_callback(logdir; every = 2)
    nuts_sample(_tb_model(); samples = 20, chains = 1, callback = fit_cb)

    files = readdir(logdir)
    @test any(startswith("events.out.tfevents"), files)
    @test any(f -> filesize(joinpath(logdir, f)) > 0, files)

    tags = Set{String}()
    map_summaries((name, step, val) -> push!(tags, name), logdir)
    ## Parameters and diagnostics land under their own grouped prefixes,
    ## the divergence flag is logged, and both scalar `value` traces and
    ## `distribution` histograms are written.
    @test any(t -> startswith(t, "params/"), tags)
    @test any(t -> startswith(t, "diagnostics/"), tags)
    @test "diagnostics/numerical_error/value" in tags
    @test any(t -> endswith(t, "/value"), tags)
    @test any(t -> endswith(t, "/distribution"), tags)
end
