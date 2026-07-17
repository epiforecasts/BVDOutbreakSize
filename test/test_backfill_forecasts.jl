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

@testitem "backfill routes pre-registry renewal tags to their own driver" begin
    using BVDOutbreakSize: BVDOutbreakSize
    include(joinpath(pkgdir(BVDOutbreakSize), "scripts",
        "backfill_forecasts.jl"))

    ## The pre-registry renewal tags are exactly the ones `preregistry_driver`
    ## records a joint call for; a registry-era tag is not among them, so
    ## `backfill_one` keeps it on the registry path.
    @test PREREGISTRY_TAGS == ("v1.4.0", "v1.5.0", "v1.6.0")
    @test all(t -> t in PREREGISTRY_TAGS, ("v1.4.0", "v1.5.0", "v1.6.0"))
    @test !("v1.7.0" in PREREGISTRY_TAGS)
end

@testitem "backfill pre-registry driver source valid and schema-matched" begin
    using BVDOutbreakSize: BVDOutbreakSize
    include(joinpath(pkgdir(BVDOutbreakSize), "scripts",
        "backfill_forecasts.jl"))

    ## Source generation needs no fit: the driver is emitted as text and only
    ## run inside a worktree. Every pre-registry tag emits parseable Julia.
    for tag in PREREGISTRY_TAGS
        src = preregistry_driver("/tmp/forecast_$(tag).csv", tag)
        @test Meta.parseall(src) isa Expr
        ## Schema-identical to `backfill_driver`: same archive columns and the
        ## live build's every-fifth-draw thinning.
        @test occursin("made_date = Date[], horizon = Int[]", src)
        @test occursin("stream = String[], draw = Int[], value = Float64[]",
            src)
        @test occursin("1:$(THIN):length(vals)", src)
        ## The tag's own headline joint call, not a registry lookup.
        @test occursin("bvd_joint(", src)
        @test !occursin("build_fit_specs", src)
    end

    ## v1.6.0 adds the isolation/bed and recovered histories; v1.4.0/v1.5.0
    ## predate those columns, so their call carries neither.
    v16 = preregistry_driver("/tmp/forecast_v1.6.0.csv", "v1.6.0")
    @test occursin("isolation_history = obs.isolation_history", v16)
    @test occursin("recovered_history = obs.recovered_history", v16)
    for tag in ("v1.4.0", "v1.5.0")
        src = preregistry_driver("/tmp/forecast_$(tag).csv", tag)
        @test !occursin("isolation_history", src)
        @test !occursin("recovered_history", src)
    end

    ## `samples`/`chains` default to the production setting and can be lowered
    ## for a cheap smoke run; `target_accept` is never restated, so each tag's
    ## own `nuts_sample` default applies.
    prod = preregistry_driver("/tmp/f.csv", "v1.6.0")
    @test occursin("samples = 1000, chains = 2", prod)
    smoke = preregistry_driver("/tmp/f.csv", "v1.6.0"; samples = 30, chains = 1)
    @test occursin("samples = 30, chains = 1", smoke)
    @test !occursin("target_accept", prod)

    ## An unrecorded tag is a hard error, not a silently empty driver.
    @test_throws ErrorException preregistry_joint_call("v1.99.0", 10, 1)
end
