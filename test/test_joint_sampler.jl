@testitem "joint sampler budget defaults and overrides" tags = [:quality] begin
    include(joinpath(@__DIR__, "..", "docs", "fits", "registry.jl"))

    unset = (
        "BVD_JOINT_SAMPLES" => nothing, "BVD_JOINT_WARMUP" => nothing,
        "BVD_JOINT_TARGET_ACCEPT" => nothing, "BVD_JOINT_MAX_DEPTH" => nothing,
    )
    withenv(unset...) do
        s = joint_sampler_args()
        @test s.samples == 1000
        @test s.n_adapts == 1000
        @test s.target_accept == 0.8
        @test s.max_depth == 10
    end
    withenv(
        "BVD_JOINT_SAMPLES" => "1200", "BVD_JOINT_WARMUP" => "400",
        "BVD_JOINT_TARGET_ACCEPT" => "0.85", "BVD_JOINT_MAX_DEPTH" => "11"
    ) do
        s = joint_sampler_args()
        @test s.samples == 1200
        @test s.n_adapts == 400
        @test s.target_accept == 0.85
        @test s.max_depth == 11
    end
end

@testitem "health-zone fits take a lower target acceptance" tags = [
    :quality,
] begin
    include(joinpath(@__DIR__, "..", "docs", "fits", "registry.jl"))

    unset = (
        "BVD_JOINT_SAMPLES" => nothing, "BVD_JOINT_WARMUP" => nothing,
        "BVD_JOINT_TARGET_ACCEPT" => nothing, "BVD_JOINT_MAX_DEPTH" => nothing,
        "BVD_ZONE_TARGET_ACCEPT" => nothing,
    )
    withenv(unset...) do
        j, z = joint_sampler_args(), zone_sampler_args()
        @test z.target_accept == 0.75
        @test j.target_accept == 0.8
        @test Base.structdiff(z, (; target_accept = 0)) ==
            Base.structdiff(j, (; target_accept = 0))
        specs = build_fit_specs(load_observations(); run_sensitivity = true)
        for s in specs
            zone = s.id in ("local", "local_frozen_validation")
            @test fit_sampler_args(s) == (zone ? z : j)
        end
    end
    ## The zone override moves only the zone target, and the joint's
    ## overrides still reach the zone fits.
    withenv(
        unset..., "BVD_ZONE_TARGET_ACCEPT" => "0.8",
        "BVD_JOINT_SAMPLES" => "300"
    ) do
        z = zone_sampler_args()
        @test z == joint_sampler_args()
        @test z.samples == 300
    end
end

@testitem "every fit shares the joint's sampler settings" tags = [
    :quality,
] begin
    ## Every fit draws from the one budget; a fit that spells its own
    ## `samples` out would drift from the others.
    root = joinpath(@__DIR__, "..")
    src = read(joinpath(root, "docs", "fits", "registry.jl"), String)
    ## The settings each fit samples at are checked by behaviour in
    ## `test_fit_cache.jl` and the zone target test above.
    @test !occursin("samples = samples", src)

    ## The scripts that fit outside the registry take the same helpers.
    for (f, helper) in (
            "fit_zone.jl" => "zone_sampler_args()",
            "zone_fit_report.jl" => "zone_sampler_args()",
            "recovery.jl" => "joint_sampler_args()",
        )
        @test occursin(helper, read(joinpath(root, "scripts", f), String))
    end
    rec = read(joinpath(root, "scripts", "recovery.jl"), String)
    @test !occursin("BVD_RECOVERY_WARMUP", rec)
    @test !occursin("BVD_RECOVERY_SAMPLES", rec)

    ## Every fit is diagnosed at the joint's tree-depth cap.
    include(joinpath(root, "docs", "fits", "registry.jl"))
    withenv("BVD_JOINT_MAX_DEPTH" => "12") do
        for id in ("joint", "deaths", "frozen_2026-05-20", "local")
            @test fit_max_depth(id) == 12
        end
    end
end

@testitem "the joint sampler settings are part of the fit cache key" tags = [
    :quality,
] begin
    include(joinpath(@__DIR__, "..", "docs", "fits", "registry.jl"))

    vars = (
        "BVD_JOINT_SAMPLES", "BVD_JOINT_WARMUP",
        "BVD_JOINT_TARGET_ACCEPT", "BVD_JOINT_MAX_DEPTH", "BVD_ZONE_TARGET_ACCEPT",
    )
    unset = Tuple(v => nothing for v in vars)
    specs = build_fit_specs(load_observations(); run_sensitivity = true)
    spec(id) = specs[findfirst(s -> s.id == id, specs)]
    key(id) = fit_key(spec(id))
    zone_ids = ("local", "local_frozen_validation")
    ids = Tuple(s.id for s in specs)
    base = withenv(unset...) do
        Dict(id => key(id) for id in ids)
    end

    ## The same settings give the same key, and so does a default set
    ## explicitly.
    withenv(unset...) do
        @test key("joint") == base["joint"]
    end
    withenv(
        "BVD_JOINT_SAMPLES" => "1000", "BVD_JOINT_WARMUP" => "1000",
        "BVD_JOINT_TARGET_ACCEPT" => "0.80", "BVD_JOINT_MAX_DEPTH" => "10"
    ) do
        @test key("joint") == base["joint"]
        @test key("sens_no_patches") == base["sens_no_patches"]
    end

    overrides = (
        "BVD_JOINT_SAMPLES" => "1200", "BVD_JOINT_WARMUP" => "400",
        "BVD_JOINT_TARGET_ACCEPT" => "0.70", "BVD_JOINT_MAX_DEPTH" => "12",
    )
    ## The health-zone fits take their own target acceptance, so the
    ## joint's leaves their keys alone.
    for ov in overrides
        withenv(unset..., ov) do
            for id in ids
                if first(ov) == "BVD_JOINT_TARGET_ACCEPT" && id in zone_ids
                    @test key(id) == base[id]
                else
                    @test key(id) != base[id]
                end
            end
        end
    end

    ## The zone target acceptance moves only the health-zone fits' keys.
    withenv(unset..., "BVD_ZONE_TARGET_ACCEPT" => "0.70") do
        for id in ids
            if id in zone_ids
                @test key(id) != base[id]
            else
                @test key(id) == base[id]
            end
        end
    end
end
