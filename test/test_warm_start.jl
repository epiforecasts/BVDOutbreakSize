## The warm-start sampler path: `FixedMetricNUTS` holds its metric through
## warm-up, `nuts_sample` with no seed is plain `NUTS`, and the fit layer's
## metric mapping (`docs/fits/warmstart.jl`).

@testitem "FixedMetricNUTS adapts the step size and keeps the metric" tags = [
    :slow,
] begin
    using Distributions: Normal, MvNormal
    using LinearAlgebra: I
    using Random: Xoshiro
    using Turing: @model, sample, NUTS
    using BVDOutbreakSize: FlexiChains
    import BVDOutbreakSize
    using BVDOutbreakSize: FixedMetricNUTS

    @model function _correlated()
        x ~ MvNormal([0.0, 0.0], [1.0 0.8; 0.8 1.0])
    end
    Minv = [1.0 0.8; 0.8 1.0]
    alg = FixedMetricNUTS(40, 0.8, Minv; ϵ = 0.5)
    chn = sample(
        Xoshiro(1), _correlated(), alg, 30; save_state = true,
        progress = false
    )
    state = only(FlexiChains.last_sampler_state(chn))
    ## The metric the sampler ends warm-up with is the one it was given.
    @test state.hamiltonian.metric.M⁻¹ == Minv
    @test state.adaptor isa BVDOutbreakSize._TracedStepSize
    ## One trace for the one chain, one entry per warm-up iteration, and the
    ## final step size is the averaged one.
    trace = only(alg.step_size_traces)
    @test length(trace) == 40
    @test exp(last(trace)[2]) ≈ state.kernel.τ.integrator.ϵ
    ## The chain carries the statistics `NUTS` gives it.
    ref = sample(Xoshiro(1), _correlated(), NUTS(40, 0.8), 30; progress = false)
    @test Set(FlexiChains.extras(chn)) == Set(FlexiChains.extras(ref))
    ## A metric of the wrong size is refused.
    @test_throws DimensionMismatch sample(
        Xoshiro(1), _correlated(), FixedMetricNUTS(10, 0.8, Matrix(1.0I, 3, 3)),
        5; progress = false
    )
end

@testitem "nuts_sample with no seed samples plain NUTS" tags = [:slow] begin
    using Distributions: Normal, LogNormal
    using Turing: @model
    import BVDOutbreakSize
    using BVDOutbreakSize: nuts_sample, NUTS

    @model function _two()
        x ~ Normal(0.0, 1.0)
        s ~ LogNormal(0.0, 0.5)
    end
    ## A warm start that declines to seed, and records what it was handed.
    struct _Declines
        seen::Base.RefValue{Any}
    end
    BVDOutbreakSize.warm_start_settings(::_Declines, model) = nothing
    function BVDOutbreakSize.warm_start_record(w::_Declines, model, chn, alg)
        w.seen[] = alg
        return nothing
    end
    w = _Declines(Ref{Any}(nothing))
    plain = nuts_sample(_two(); samples = 40, chains = 2)
    declined = nuts_sample(_two(); samples = 40, chains = 2, warm_start = w)
    @test Array(plain[:x]) == Array(declined[:x])
    @test Array(plain[:s]) == Array(declined[:s])
    @test w.seen[] isa NUTS
end

@testitem "linked_draws gives each draw in the linked space" tags = [
    :slow,
] begin
    using Distributions: Normal, LogNormal
    using Turing: @model
    using BVDOutbreakSize: nuts_sample, linked_draws, linked_layout

    @model function _two()
        x ~ Normal(0.0, 1.0)
        s ~ LogNormal(0.0, 0.5)
    end
    chn = nuts_sample(_two(); samples = 20, chains = 2)
    ld = linked_draws(_two(), chn)
    @test ld.layout == linked_layout(_two())
    @test first.(ld.layout) == ["x", "s"]
    @test size(ld.draws) == (40, 2)
    @test ld.chain_index == repeat(1:2; inner = 20)
    @test ld.draws[:, 1] ≈ vec(Array(chn[:x]))
    @test ld.draws[:, 2] ≈ log.(vec(Array(chn[:s])))
end

@testitem "warm-start metric maps the parent's draws by variable" begin
    using LinearAlgebra: diag, isposdef, Symmetric
    using Random: Xoshiro
    include(joinpath(@__DIR__, "..", "docs", "fits", "registry.jl"))

    rng = Xoshiro(3)
    ## A parent with `a` (2 entries) and `b` (3 entries) over two chains;
    ## the child keeps `a`, grows `b` to 4 and adds `c`.
    X = randn(rng, 400, 5) .* [1.0 2.0 0.5 0.5 0.5]
    X[:, 2] .+= 0.9 .* 2.0 .* X[:, 1]
    parent = (;
        layout = ["a" => 1:2, "b" => 3:5], draws = X,
        chain_index = repeat(1:2; inner = 200),
    )
    child = ["a" => 1:2, "b" => 3:6, "c" => 7:7]
    no_shrink = seeded_inverse_metric(parent, child; shrink = 0.0)
    M = no_shrink.inverse_metric
    @test size(M) == (7, 7)
    @test no_shrink.coverage ≈ 5 / 7
    ## Mapped entries carry the parent's within-chain covariance.
    Xc = copy(X)
    for c in 1:2
        rows = parent.chain_index .== c
        Xc[rows, :] .-= sum(Xc[rows, :]; dims = 1) ./ count(rows)
    end
    S = Xc' * Xc ./ (size(Xc, 1) - 1)
    @test M[1:5, 1:5] ≈ S
    ## The new entry of `b` takes `b`'s median variance, `c` the median of
    ## all, neither correlated with anything.
    @test M[6, 6] ≈ sort(diag(S)[3:5])[2]
    @test M[7, 7] ≈ sort(diag(S))[3]
    @test all(iszero, M[6, [1:5; 7]])
    @test all(iszero, M[7, 1:6])
    ## Full shrinkage keeps the variances and drops every correlation.
    D = seeded_inverse_metric(parent, child; shrink = 1.0).inverse_metric
    @test diag(D) ≈ diag(M)
    @test D[1, 2] == 0
    shrunk = seeded_inverse_metric(parent, child).inverse_metric
    @test isposdef(Symmetric(shrunk))
end

@testitem "warm-start parent checks" begin
    using Serialization: serialize
    include(joinpath(@__DIR__, "..", "docs", "fits", "registry.jl"))

    dir = mktempdir()
    function state(;
            key = "joint__aaaa", gate = :pass,
            models = warm_start_models_hash(), schema = WARM_START_SCHEMA
        )
        return (;
            schema, id = "joint", key, models_hash = models, gate,
            metric_kind = :adapted_diag, step_size = 0.01,
            layout = ["a" => 1:1], draws = randn(10, 1),
            chain_index = ones(Int, 10),
        )
    end
    put(s) = serialize(warm_start_file(dir, s.key), s)
    seed(key = "joint__bbbb") = fit_warm_start("joint", key; parent_dir = dir)
    withenv("BVD_WARMSTART" => nothing, "BVD_WARMSTART_RULE" => nothing) do
        @test seed().parent === nothing
        put(state())
        @test seed().parent.key == "joint__aaaa"
        ## The same key, a failed gate, a model change and an old schema
        ## each fall back to the full warm-up.
        @test seed("joint__aaaa").parent === nothing
        put(state(; gate = :fail))
        @test occursin("gate", seed().reason)
        put(state(; models = "x"))
        @test occursin("src/models", seed().reason)
        ## The coverage rule does not fall back on a model change.
        withenv("BVD_WARMSTART_RULE" => "coverage") do
            @test seed().parent !== nothing
        end
        put(state(; schema = "v0"))
        @test occursin("schema", seed().reason)
        ## Another fit id's state is not a parent.
        rm(warm_start_file(dir, "joint__aaaa"))
        put(state(; key = "joint_x__cccc"))
        @test seed().parent === nothing
    end
    withenv("BVD_WARMSTART" => "false") do
        @test warm_start_key() == "warm:off"
        @test seed().parent === nothing
    end
end
