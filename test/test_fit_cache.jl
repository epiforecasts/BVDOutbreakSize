## Tests for the content-addressed fit cache (`docs/fits/cache.jl`): a cache
## hit reuses the serialised result, a miss (or forced refit) runs the thunk,
## and the content hash changes when the inputs change.

@testitem "fit_or_load caches, reuses and refits" tags=[:quality] begin
    include(joinpath(@__DIR__, "..", "docs", "fits", "cache.jl"))

    dir = mktempdir()
    calls = Ref(0)
    thunk = () -> (calls[] += 1; (payload = "chain", n = 42))
    key = "demo__" * content_hash([@__FILE__]; extra = "settings")

    r1 = fit_or_load(key, thunk; cache_dir = dir)                # miss → fit
    r2 = fit_or_load(key, thunk; cache_dir = dir)                # hit → load
    r3 = fit_or_load(key, thunk; cache_dir = dir, refit = true)  # forced refit

    @test calls[] == 2                     # fit once + forced refit once
    @test r1 == r2 == r3
    @test isfile(joinpath(dir, key * ".jls"))
end

@testitem "fit_or_load strict mode errors on a miss instead of fitting" tags=[
    :quality
] begin
    include(joinpath(@__DIR__, "..", "docs", "fits", "cache.jl"))

    dir = mktempdir()
    calls = Ref(0)
    thunk = () -> (calls[] += 1; (payload = "chain", n = 42))
    key = "demo__" * content_hash([@__FILE__]; extra = "strict")

    ## A strict miss must throw (naming the key and dir) and never run the
    ## thunk, so a render can fail fast rather than silently refit the whole
    ## report.
    @test_throws Exception fit_or_load(
        key, thunk; cache_dir = dir, strict = true)
    @test calls[] == 0
    @test !isfile(joinpath(dir, key * ".jls"))

    ## Once the fit exists, strict mode loads it like a normal hit.
    fit_or_load(key, thunk; cache_dir = dir)      # populate (non-strict miss)
    r = fit_or_load(key, thunk; cache_dir = dir, strict = true)
    @test r == (payload = "chain", n = 42)
    @test calls[] == 1                            # only the populating fit ran
end

@testitem "every score_releases overlay is excluded from the fit hash" tags=[
    :quality
] begin
    include(joinpath(@__DIR__, "..", "docs", "fits", "registry.jl"))

    ## score_releases.jl runs in the render job (before rendering) and writes
    ## overlay CSVs into data/. Every such file must be in FIT_DATA_EXCLUDE, or
    ## the render's data-dir hash diverges from the fit matrix's, every fit key
    ## changes, and the render misses the whole cache (a 2h refit / strict-mode
    ## failure). Auto-derive the written files from the script so a new overlay
    ## that forgets the exclusion fails here instead of in CI.
    ## Matched on the script's own write path (`@__DIR__/../data/...`), not
    ## on any mention of the data directory: the script also reads a fit
    ## input from there (the digitised onset triangle), which must stay in
    ## the hash rather than be excluded from it.
    src = read(
        joinpath(@__DIR__, "..", "scripts", "score_releases.jl"), String)
    written = Set(m.captures[1]
    for m in eachmatch(
        r"@__DIR__,\s*\"\.\.\",\s*\"data\",\s*\"([\w.]+\.csv)\"", src))
    @test length(written) >= 4  # guards against a silent regex miss
    for f in written
        @test f in FIT_DATA_EXCLUDE
    end
end

@testitem "content hash reflects inputs" tags=[:quality] begin
    include(joinpath(@__DIR__, "..", "docs", "fits", "cache.jl"))

    h = content_hash([@__FILE__]; extra = "a")
    @test length(h) == 16
    ## The extra settings string is part of the hash.
    @test content_hash([@__FILE__]; extra = "a") == h
    @test content_hash([@__FILE__]; extra = "b") != h
    ## A missing source file hashes deterministically (does not throw).
    @test content_hash(["/no/such/file"]; extra = "a") isa String
    @test file_sha256("/no/such/file") == "absent"

    ## A directory tree hashes order-independently and picks up changes.
    d = mktempdir()
    write(joinpath(d, "b.csv"), "2")
    write(joinpath(d, "a.csv"), "1")
    t1 = tree_sha256(d)
    write(joinpath(d, "a.csv"), "1x")
    @test tree_sha256(d) != t1
end

@testitem "content hash can exclude non-input data files" tags=[:quality] begin
    include(joinpath(@__DIR__, "..", "docs", "fits", "cache.jl"))

    ## Excluding a file removes it from the tree digest; a fit-input CSV still
    ## contributes. This mirrors `fit_content_hash` ignoring the published
    ## `released_estimates.csv` overlay while genuine data changes still refit.
    d = mktempdir()
    write(joinpath(d, "observations.csv"), "1")
    write(joinpath(d, "released_estimates.csv"), "overlay")
    excl = ("released_estimates.csv",)

    ## The overlay is not part of the excluded digest.
    t_full = tree_sha256(d)
    t_excl = tree_sha256(d; exclude = excl)
    @test t_full != t_excl

    ## Changing the excluded file leaves the excluded digest unchanged.
    write(joinpath(d, "released_estimates.csv"), "overlay-v2")
    @test tree_sha256(d; exclude = excl) == t_excl

    ## Changing a fit-input file still changes the excluded digest.
    write(joinpath(d, "observations.csv"), "2")
    @test tree_sha256(d; exclude = excl) != t_excl

    ## The same guarantees hold through `content_hash`'s `data_exclude`.
    src = [@__FILE__]
    h = content_hash(src; data_dir = d, data_exclude = excl)
    write(joinpath(d, "released_estimates.csv"), "overlay-v3")
    @test content_hash(src; data_dir = d, data_exclude = excl) == h
    write(joinpath(d, "observations.csv"), "3")
    @test content_hash(src; data_dir = d, data_exclude = excl) != h
end

@testitem "validation fits follow the reporting status" tags=[:quality] begin
    using Dates
    using Dates: Date, Day

    include(joinpath(@__DIR__, "..", "docs", "fits", "registry.jl"))

    ## The validation panels take the still-reported streams, so a stream the
    ## situation reports have stopped updating must not be fitted at all: its
    ## fit feeds an overlay that is filtered out, and each one is a NUTS fit
    ## a docs build runs serially.
    cutoff = Date(2026, 7, 15)
    n = 60
    day(d) = n - Dates.value(cutoff - d)
    live = (; days = [day(cutoff - Day(1))], counts = [10.0])
    stale = (; days = [day(cutoff - Day(60))], counts = [10.0])
    obs = (; cutoff = cutoff, n = n,
        reported_history = stale, deaths_history = stale,
        confirmed_history = live, confirmed_deaths_history = live,
        isolation_history = live)

    @test validation_stream_ids(obs) ==
          ("confirmed", "confirmed_deaths", "treatment")

    ## A stream that starts being reported again comes back on its own.
    revived = merge(obs, (; reported_history = live))
    @test "cases" in validation_stream_ids(revived)

    ## Every id names a single-stream fit the registry can build.
    ids = fit_ids(load_observations(); run_sensitivity = false)
    for sid in validation_stream_ids(load_observations())
        @test "frozen_validation_$sid" in ids
        @test sid in ids
    end
end

@testitem "dependent fits follow their parents in the registry" tags=[
    :quality
] begin
    include(joinpath(@__DIR__, "..", "docs", "fits", "registry.jl"))

    obs = load_observations()
    specs = build_fit_specs(obs; run_sensitivity = true)

    ## Every parent a spec names is listed before it.
    seen = String[]
    for s in specs
        @test s.needs isa Vector{String}
        for parent in s.needs
            @test parent in seen
        end
        push!(seen, s.id)
    end
    @test validate_fit_specs(specs) === specs
    @test_throws Exception validate_fit_specs([
        (; id = "child", needs = ["parent"]),
        (; id = "parent", needs = String[])])

    ## The health-zone fits are the dependent stage and nothing else is.
    base = base_fit_ids(obs; run_sensitivity = true)
    dependent = dependent_fit_ids(obs; run_sensitivity = true)
    @test dependent == ["local", "local_frozen_validation", "local_mixing",
        "local_deaths", "local_parent_low", "local_parent_high"]
    @test dependent_fit_ids(obs; run_sensitivity = false) ==
          ["local", "local_frozen_validation"]
    @test "joint" in base
    @test "frozen_validation" in base
    @test isempty(intersect(base, dependent))
    @test sort(vcat(base, dependent)) ==
          sort(fit_ids(obs; run_sensitivity = true))
    @test fit_ids(obs; run_sensitivity = true, stage = :base) == base
    @test fit_ids(obs; run_sensitivity = true, stage = :dependent) ==
          dependent
    @test_throws Exception fit_ids(obs; stage = :nonsense)

    ## `BVD_FIT_STAGE` picks the stage for list.jl and all.jl.
    withenv("BVD_FIT_STAGE" => nothing) do
        @test fit_stage_env(:base) === :base
        @test fit_stage_env(:all) === :all
    end
    withenv("BVD_FIT_STAGE" => "dependent") do
        @test fit_stage_env(:base) === :dependent
    end
    withenv("BVD_FIT_STAGE" => "later") do
        @test_throws Exception fit_stage_env(:base)
    end
end

@testitem "dependent thunks load their parent from the cache" tags=[
    :quality
] begin
    using Dates: Date

    include(joinpath(@__DIR__, "..", "docs", "fits", "registry.jl"))

    obs = load_observations()
    dir = mktempdir()
    calls = Any[]
    ## A stand-in for `BVDOutbreakSize.fit_zone`, recording what it was
    ## handed and returning a fake chain.
    fake_zone = (parent, o; kwargs...) -> begin
        push!(calls, (; parent, o, kwargs...))
        (; zone = "chain", from = parent)
    end
    specs = build_fit_specs(obs; run_sensitivity = true,
        zone_fitter = fake_zone, cache_dir = dir)
    spec(id) = specs[findfirst(s -> s.id == id, specs)]

    withenv("BVD_FIT_LOG" => "none", "BVD_ZONE_SAMPLES" => nothing,
        "BVD_ZONE_WARMUP" => nothing) do
        ## A missing parent is an error, not a refit: the zone fitter must
        ## not run.
        @test_throws Exception spec("local").thunk()
        @test_throws Exception spec("local_frozen_validation").thunk()
        @test isempty(calls)

        ## A fake headline under the joint's own key reaches the zone fitter
        ## with the current observations and the zone sampler settings.
        fake_joint = (; payload = "joint chain")
        fit_or_load(fit_key("joint"), () -> fake_joint; cache_dir = dir)
        r = spec("local").thunk()
        @test r == (; zone = "chain", from = fake_joint)
        @test length(calls) == 1
        c = calls[1]
        @test c.parent == fake_joint
        @test c.o === obs
        @test c.samples == 600
        @test c.n_adapts == 400
        @test c.target_accept == 0.8
        @test c.max_depth == 8
        @test c.chains == 2
        @test c.callback === nothing

        ## The draw counts follow the zone overrides.
        withenv("BVD_ZONE_SAMPLES" => "50", "BVD_ZONE_WARMUP" => "25") do
            spec("local").thunk()
        end
        @test calls[2].samples == 50
        @test calls[2].n_adapts == 25

        ## The frozen dependent melds from the frozen parent's chain on that
        ## parent's observations and returns the `(; cutoff, o, chn)` shape.
        frozen_o = (; n = 42, cutoff = Date(2026, 9, 1))
        fake_frozen = (; cutoff = frozen_o.cutoff, o = frozen_o,
            chn = (; payload = "frozen chain"))
        fit_or_load(fit_key("frozen_validation"), () -> fake_frozen;
            cache_dir = dir)
        f = spec("local_frozen_validation").thunk()
        @test keys(f) == (:cutoff, :o, :chn)
        @test f.cutoff == frozen_o.cutoff
        @test f.o == frozen_o
        @test f.chn == (; zone = "chain", from = fake_frozen.chn)
        @test calls[3].parent == fake_frozen.chn
        @test calls[3].o == frozen_o

        ## Each zone sensitivity variant melds from the same joint and adds
        ## exactly its own switch to the zone fitter's keywords.
        variants = (("local_mixing", :mixing, true),
            ("local_deaths", :deaths, true),
            ("local_parent_low", :parent_summary, :draw_low),
            ("local_parent_high", :parent_summary, :draw_high))
        for (id, key, value) in variants
            spec(id).thunk()
            c = calls[end]
            @test c.parent == fake_joint
            @test c.o === obs
            @test getproperty(c, key) == value
            @test c.samples == 600
        end
        @test !haskey(calls[1], :mixing)
        @test !haskey(calls[1], :parent_summary)
    end

    ## Without an injected fitter the package function is looked up when the
    ## thunk runs.
    lazy = build_fit_specs(obs; run_sensitivity = false, cache_dir = dir)
    @test any(s -> s.id == "local", lazy)
    if !isdefined(BVDOutbreakSize, :fit_zone)
        @test_throws Exception default_zone_fitter()
    else
        @test default_zone_fitter() === BVDOutbreakSize.fit_zone
    end
end
