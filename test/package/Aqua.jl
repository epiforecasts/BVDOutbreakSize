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
    # DataFramesMeta isn't used in src/ — the data-wrangling macros
    # live in scripts/refresh_insp_data.jl, which runs under
    # --project=. and needs the package present.
    Aqua.test_stale_deps(BVDOutbreakSize; ignore = [:DataFramesMeta])
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
    ## Give the check a longer settle window (`tmax`): loading the package
    ## pulls heavy precompiled dependencies (Turing, Makie,
    ## CensoredDistributions), and on a loaded CI runner the default
    ## threshold can false-positive on still-finishing precompile work
    ## rather than a genuine lingering task.
    Aqua.test_persistent_tasks(BVDOutbreakSize; tmax = 30.0)
end
