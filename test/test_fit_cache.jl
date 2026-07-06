## Tests for the content-addressed fit cache (`docs/fits/cache.jl`): a cache
## hit reuses the serialised result, a miss (or forced refit) runs the thunk,
## and the content hash changes when the inputs change.

@testitem "fit_or_load caches, reuses and refits" tags=[:quality] begin
    include(joinpath(@__DIR__, "..", "docs", "fits", "cache.jl"))

    dir = mktempdir()
    calls = Ref(0)
    thunk = () -> (calls[] += 1; (payload = "chain", n = 42))
    key = "demo__" * content_hash([@__FILE__]; extra = "settings")

    r1 = fit_or_load(key, thunk; cache_dir = dir)                 # miss → fit
    r2 = fit_or_load(key, thunk; cache_dir = dir)                 # hit  → load
    r3 = fit_or_load(key, thunk; cache_dir = dir, refit = true)   # forced refit

    @test calls[] == 2                     # fit once + forced refit once
    @test r1 == r2 == r3
    @test isfile(joinpath(dir, key * ".jls"))
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
