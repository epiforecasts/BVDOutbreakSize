@testitem "Aqua: Ambiguities" tags=[:quality] begin
    using Aqua, BVDOutbreakSize
    Aqua.test_ambiguities(BVDOutbreakSize)
end

@testitem "Aqua: unbound_args" tags=[:quality] begin
    using Aqua, BVDOutbreakSize
    Aqua.test_unbound_args(BVDOutbreakSize)
end

@testitem "Aqua: undefined_exports" tags=[:quality] begin
    using Aqua, BVDOutbreakSize
    Aqua.test_undefined_exports(BVDOutbreakSize)
end

@testitem "Aqua: project_extras" tags=[:quality] begin
    using Aqua, BVDOutbreakSize
    Aqua.test_project_extras(BVDOutbreakSize)
end

@testitem "Aqua: stale_deps" tags=[:quality] begin
    using Aqua, BVDOutbreakSize
    # DataFramesMeta and Chain aren't used in src/ — the data-wrangling
    # macros live in the scripts/ data tooling (confirm_insp_data.jl),
    # which runs under --project=. and needs the packages present.
    Aqua.test_stale_deps(BVDOutbreakSize; ignore = [:DataFramesMeta, :Chain])
end

@testitem "Aqua: deps_compat" tags=[:quality] begin
    using Aqua, BVDOutbreakSize
    Aqua.test_deps_compat(BVDOutbreakSize)
end

@testitem "Aqua: piracies" tags=[:quality] begin
    using Aqua, BVDOutbreakSize
    Aqua.test_piracies(BVDOutbreakSize)
end

@testitem "Aqua: persistent_tasks" tags=[:quality] begin
    using Aqua, BVDOutbreakSize
    using Test: collect_test_logs
    # Aqua's probe builds a throwaway package depending on this one,
    # precompiles it in a subprocess, and has that subprocess write a sentinel
    # file once `using BVDOutbreakSize` returns. It then allows `tmax` for the
    # subprocess to exit; one still alive after loading is the persistent task
    # the check looks for.
    #
    # The probe conflates two outcomes. When the subprocess exits WITHOUT
    # writing the sentinel — a failed or killed precompilation, nothing to do
    # with tasks — Aqua logs `<file> was not created, but precompilation
    # exited` and reports a persistent task regardless. `tmax` does not cover
    # that phase: the wait for the sentinel has no timeout at all, so raising
    # `tmax` cannot help. The durations seen, 1m39 and 8m48 against
    # `tmax = 600`, were the precompilation dying rather than a timeout.
    #
    # Rebuilding this package's dependency stack in a fresh project while this
    # process already holds it loaded is heavy enough that the subprocess dies
    # intermittently on the runners: two of three concluded
    # `Julia 1 - ubuntu-latest` runs on 27 July 2026, one of them on main.
    #
    # So the two outcomes are separated here. A genuine detection still fails.
    # A probe that could not run is skipped with its reason visible, rather
    # than reported as a package defect. Issue #495 tracks the root cause,
    # which this does not address.
    #
    # macOS is excluded outright: that runner leaves a background task alive
    # after loading the graphics stack (CairoMakie/Makie), so the check cannot
    # pass there.
    if Sys.isapple()
        @test_skip "persistent_tasks: macOS CI runner leaves a graphics task"
    else
        ## `Base.CoreLogging.Error` rather than `Logging.Error`: the
        ## constants are the same object, and `Logging` is not a
        ## declared dependency of the test project, so importing it
        ## resolves locally but not in a clean CI environment.
        logs, detected = collect_test_logs(
            min_level = Base.CoreLogging.Error) do
            Aqua.has_persistent_tasks(Base.PkgId(BVDOutbreakSize); tmax = 600)
        end
        ## collect_test_logs swallows the records, so re-emit them: the probe's
        ## own diagnostics are the only clue to why it failed.
        for rec in logs
            @warn "Aqua persistent_tasks probe" message = string(rec.message)
        end
        incomplete = any(logs) do rec
            occursin("was not created", string(rec.message))
        end
        if incomplete
            @test_skip "persistent_tasks: probe did not complete, so the " *
                       "presence of a task is unknown (#495)"
        else
            @test !detected
        end
    end
end
