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
