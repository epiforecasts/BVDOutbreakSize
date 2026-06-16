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
    # Skipped on macOS: the macOS CI runner leaves a background task alive
    # after loading the heavy graphics stack (CairoMakie/Makie), so this
    # check fails there even at tmax = 60 s. It passes on Linux and Windows
    # and locally, so this is a runner-environment issue, not package code.
    # Keep the real check on the other platforms; raise tmax for the heavy
    # AD/plotting load. The Enzyme precompile-and-load step is slow on the
    # Julia 1.11 ubuntu runner (Enzyme alone takes ~3.5 min to build there),
    # so the load can exceed 120 s and the check then false-positives a
    # persistent task; it passes on Windows and Julia LTS. Give it a generous
    # timeout so the slow load is not mistaken for a lingering task.
    if Sys.isapple()
        @test_skip "persistent_tasks: macOS CI runner leaves a graphics task"
    else
        Aqua.test_persistent_tasks(BVDOutbreakSize; tmax = 600)
    end
end
