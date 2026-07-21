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

@testitem "fit_or_load strict mode errors on a miss instead of fitting" tags=[:quality] begin
    include(joinpath(@__DIR__, "..", "docs", "fits", "cache.jl"))

    dir = mktempdir()
    calls = Ref(0)
    thunk = () -> (calls[] += 1; (payload = "chain", n = 42))
    key = "demo__" * content_hash([@__FILE__]; extra = "strict")

    ## A strict miss must throw (naming the key and dir) and never run the thunk,
    ## so a render can fail fast rather than silently refit the whole report.
    @test_throws Exception fit_or_load(key, thunk; cache_dir = dir, strict = true)
    @test calls[] == 0
    @test !isfile(joinpath(dir, key * ".jls"))

    ## Once the fit exists, strict mode loads it like a normal hit.
    fit_or_load(key, thunk; cache_dir = dir)          # populate (non-strict miss)
    r = fit_or_load(key, thunk; cache_dir = dir, strict = true)
    @test r == (payload = "chain", n = 42)
    @test calls[] == 1                                # only the populating fit ran
end

@testitem "every score_releases overlay is excluded from the fit hash" tags=[:quality] begin
    include(joinpath(@__DIR__, "..", "docs", "fits", "registry.jl"))

    ## score_releases.jl runs in the render job (before rendering) and writes
    ## overlay CSVs into data/. Every such file MUST be in FIT_DATA_EXCLUDE, or
    ## the render's data-dir hash diverges from the fit matrix's, every fit key
    ## changes, and the render misses the whole cache (a 2h refit / strict-mode
    ## failure). Auto-derive the written files from the script so a new overlay
    ## that forgets the exclusion fails here instead of in CI.
    src = read(joinpath(@__DIR__, "..", "scripts", "score_releases.jl"), String)
    written = Set(m.captures[1]
    for m in eachmatch(r"\"data\",\s*\"([\w.]+\.csv)\"", src))
    @test length(written) >= 4          # guards against the regex silently missing all
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
