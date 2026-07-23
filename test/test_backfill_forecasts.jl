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

@testitem "backfill routes integral-era tags to the integral driver" begin
    using BVDOutbreakSize: BVDOutbreakSize
    include(joinpath(pkgdir(BVDOutbreakSize), "scripts",
        "backfill_forecasts.jl"))

    ## v1.2.0 and v1.3.0 are the integral-era tags `integral_driver`
    ## reconstructs inline; neither is a pre-registry renewal tag and neither
    ## carries a fit registry, so `backfill_one` reaches them only through the
    ## `INTEGRAL_TAGS` branch. `main` appends these tags to the run so
    ## `--only v1.2.0`/`--only v1.3.0` resolve at all. v1.0.0/v1.1.0 are not
    ## inlined, so they are absent here.
    @test INTEGRAL_TAGS == ("v1.2.0", "v1.3.0")
    for tag in INTEGRAL_TAGS
        @test !(tag in PREREGISTRY_TAGS)
        @test !has_registry(tag)
    end
    @test !("v1.0.0" in INTEGRAL_TAGS)
    @test !("v1.1.0" in INTEGRAL_TAGS)
end

@testitem "backfill integral driver source valid and schema-matched" begin
    using BVDOutbreakSize: BVDOutbreakSize
    include(joinpath(pkgdir(BVDOutbreakSize), "scripts",
        "backfill_forecasts.jl"))

    ## Source generation needs no fit: the driver is emitted as text and only
    ## run inside a worktree. Both inlined integral tags emit parseable Julia
    ## with the shared archive schema, the every-fifth-draw thinning, a chain
    ## serialised before the forecast, and the guarded stream superset (labels
    ## match the scorer's `STREAM_HISTORY`/`STREAM_ASSEMBLED`).
    for tag in INTEGRAL_TAGS
        src = integral_driver("/tmp/forecast_$(tag).csv", tag)
        @test Meta.parseall(src) isa Expr
        @test occursin("made_date = Date[], horizon = Int[]", src)
        @test occursin("stream = String[], draw = Int[]", src)
        @test occursin("value = Float64[]", src)
        @test occursin("1:$(THIN):length(vals)", src)
        ## The tag's own headline joint call, not a registry lookup.
        @test occursin("bvd_joint(obs.exported_cases", src)
        @test !occursin("build_fit_specs", src)
        @test occursin("serialize(chain_path, chn)", src)
        for label in ("confirmed cases", "confirmed deaths", "reported cases",
            "suspected deaths", "exports")
            @test occursin("\"$(label)\"", src)
        end
        @test occursin("col in propertynames(fc) || continue", src)
        ## `samples`/`chains` default to the production setting and can be
        ## lowered for a cheap smoke run.
        @test occursin("samples = 1000, chains = 2", src)
        smoke = integral_driver("/tmp/f.csv", tag; samples = 30, chains = 1)
        @test occursin("samples = 30, chains = 1", smoke)
    end

    ## v1.3.0 reproduces the analysis.jl-local `joint_obs` helper, forecasts
    ## the confirmed streams and restates its own `target_accept`. v1.2.0 has
    ## no local model block, forecasts the reported/suspected/export streams
    ## and takes its `nuts_sample` default (no `target_accept`).
    v13 = integral_driver("/tmp/f.csv", "v1.3.0")
    @test occursin("function joint_obs", v13)
    @test occursin("obs_confirmed = obs.confirmed_cases", v13)
    @test occursin("target_accept = 0.9", v13)
    v12 = integral_driver("/tmp/f.csv", "v1.2.0")
    @test !occursin("function joint_obs", v12)
    @test occursin("obs_cases = obs.reported_cases", v12)
    @test occursin("obs_exports = obs.exported_cases", v12)
    @test !occursin("target_accept", v12)

    ## An unrecorded integral tag is a hard error, not a silently empty driver.
    @test_throws ErrorException integral_fit_forecast("/tmp/f.csv",
        "v1.99.0", 10, 1)
end
