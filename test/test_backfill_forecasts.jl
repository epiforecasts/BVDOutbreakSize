## Tests for the flag parsing in scripts/backfill_forecasts.jl. The script's
## `main` is guarded by a PROGRAM_FILE check, so including it here defines the
## helpers without running any fit.

@testitem "backfill parse_args reads only, keep and concurrency" begin
    using BVDOutbreakSize: BVDOutbreakSize
    include(joinpath(pkgdir(BVDOutbreakSize), "scripts",
        "backfill_forecasts.jl"))

    ## Defaults: all releases, worktrees removed, conservative concurrency.
    ## The concurrency default is asserted against its intended literal so a
    ## wrong DEFAULT_CONCURRENCY constant is caught rather than passing
    ## vacuously.
    d = parse_args(String[])
    @test d.only === nothing
    @test d.keep == false
    @test d.concurrency == 2
    @test DEFAULT_CONCURRENCY == 2

    ## Flags are read in any order.
    a = parse_args(["--only", "v1.8.0", "--keep", "--concurrency", "3"])
    @test a.only == "v1.8.0"
    @test a.keep == true
    @test a.concurrency == 3

    b = parse_args(["--concurrency", "4", "--only", "v1.7.0"])
    @test b.only == "v1.7.0"
    @test b.keep == false
    @test b.concurrency == 4
end

@testitem "backfill parse_args rejects a non-positive concurrency" begin
    using BVDOutbreakSize: BVDOutbreakSize
    include(joinpath(pkgdir(BVDOutbreakSize), "scripts",
        "backfill_forecasts.jl"))

    @test_throws ErrorException parse_args(["--concurrency", "0"])
end

@testitem "backfill parse_args rejects a non-integer concurrency" begin
    using BVDOutbreakSize: BVDOutbreakSize
    include(joinpath(pkgdir(BVDOutbreakSize), "scripts",
        "backfill_forecasts.jl"))

    @test_throws ErrorException parse_args(["--concurrency", "abc"])
end

@testitem "backfill parse_args rejects a dangling flag value" begin
    using BVDOutbreakSize: BVDOutbreakSize
    include(joinpath(pkgdir(BVDOutbreakSize), "scripts",
        "backfill_forecasts.jl"))

    ## A trailing flag with no value is a user error, not a silent no-op.
    @test_throws ErrorException parse_args(["--only"])
    @test_throws ErrorException parse_args(["--concurrency"])
end
