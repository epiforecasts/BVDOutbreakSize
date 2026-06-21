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

@testitem "nuts_sample warmup=true streams adaptation steps" tags=[:slow] begin
    using Distributions: Normal
    using Turing: @model
    import Turing
    using BVDOutbreakSize: nuts_sample

    @model function _wu_model()
        x ~ Normal(0.0, 1.0)
    end

    ## Collect the step size the callback sees on each kept draw, pulled
    ## through the same `ParamsWithStats` interface the callbacks use.
    function step_sizes(; warmup)
        sizes = Float64[]
        lk = ReentrantLock()
        cb = function (rng, model, sampler, transition, state, iteration;
                kwargs...)
            pws = Turing.AbstractMCMC.ParamsWithStats(
                model, sampler, transition, state; stats = true)
            ss = get(pws.stats, :step_size, missing)
            ss isa Real && Base.@lock lk push!(sizes, ss)
            return nothing
        end
        nuts_sample(_wu_model(); samples = 100, chains = 1,
            callback = cb, warmup = warmup)
        return sizes
    end

    ## Adaptation retunes the step size every warmup step, so streaming
    ## warmup exposes many distinct step sizes; the post-warmup phase alone
    ## keeps it frozen (bar the single boundary draw as adaptation ends).
    @test length(unique(step_sizes(; warmup = true))) > 10
    @test length(unique(step_sizes(; warmup = false))) <= 3
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

@testitem "combined_callback composes survivors and prunes nothings" begin
    using BVDOutbreakSize: combined_callback
    ## All-`nothing` collapses so `nuts_sample` sees no callback.
    @test combined_callback(nothing, nothing) === nothing
    ## A single survivor is returned unwrapped.
    f = (args...; kwargs...) -> nothing
    @test combined_callback(f, nothing) === f
    ## Several survivors fire in order on each step.
    calls = Int[]
    a = (args...; kwargs...) -> push!(calls, 1)
    b = (args...; kwargs...) -> push!(calls, 2)
    cb = combined_callback(a, nothing, b)
    cb(nothing, nothing, nothing, nothing, nothing, 1)
    @test calls == [1, 2]
end

@testitem "fit_callback: BVD_FIT_LOG=none disables logging" begin
    using BVDOutbreakSize: fit_callback
    @test fit_callback("x"; spec = "none") === nothing
    ## Case-insensitive and whitespace-tolerant.
    @test fit_callback("x"; spec = " NONE ") === nothing
end

@testitem "fit_callback: progress streams a named log" tags=[:slow] begin
    using Distributions: Normal
    using Turing: @model
    using BVDOutbreakSize: fit_callback, nuts_sample

    @model function _fc_model()
        x ~ Normal(0.0, 1.0)
    end

    dir = mktempdir()
    cb = fit_callback("unit"; logdir = dir, spec = "progress")
    @test cb !== nothing
    nuts_sample(_fc_model(); samples = 20, chains = 1, callback = cb)
    ## The progress file is named after the fit; no TensorBoard run dir.
    @test isfile(joinpath(dir, "unit.log"))
    @test !isdir(joinpath(dir, "tensorboard", "unit"))
end

@testitem "fit_callback: all streams progress and TensorBoard" tags=[:slow] begin
    using Distributions: Normal
    using Turing: @model
    using TensorBoardLogger
    using BVDOutbreakSize: fit_callback, nuts_sample

    @model function _fc_all_model()
        x ~ Normal(0.0, 1.0)
    end

    dir = mktempdir()
    cb = fit_callback("unit"; logdir = dir, spec = "all")
    @test cb !== nothing
    nuts_sample(_fc_all_model(); samples = 20, chains = 1, callback = cb)
    @test isfile(joinpath(dir, "unit.log"))
    tbdir = joinpath(dir, "tensorboard", "unit")
    @test isdir(tbdir)
    @test any(startswith("events.out.tfevents"), readdir(tbdir))
end
