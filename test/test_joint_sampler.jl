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

@testitem "every fit shares the joint's sampler settings" tags = [
    :quality,
] begin
    ## Every fit draws from the one budget. Each sampler call in the
    ## registry splats the one helper; a fit that spells its own `samples`,
    ## `n_adapts` or `target_accept` out would drift from the others.
    root = joinpath(@__DIR__, "..")
    src = read(joinpath(root, "docs", "fits", "registry.jl"), String)
    ## Every `nuts_sample` call, plus the one health-zone fitter call.
    n_calls = count(r"nuts_sample\(", src) + count("return fitter(", src)
    @test n_calls >= 13
    @test count("joint_sampler_args()...", src) == n_calls
    @test !occursin("samples = samples", src)
    @test count("n_adapts = joint_warmup(", src) == 1
    @test count("target_accept = ", src) == 1
    @test !occursin("zone_sampler_args", src)

    ## The scripts that fit outside the registry take the same helper.
    for f in ("fit_zone.jl", "zone_fit_report.jl", "recovery.jl")
        @test occursin(
            "joint_sampler_args()", read(joinpath(root, "scripts", f), String)
        )
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
        "BVD_JOINT_TARGET_ACCEPT", "BVD_JOINT_MAX_DEPTH",
    )
    unset = Tuple(v => nothing for v in vars)
    ids = (
        "joint", "sens_no_patches", "frozen_validation",
        "sens_community_delay", "sens_exp_growth_clock", "deaths",
        "frozen_validation_cases", "frozen_2026-05-20", "local",
        "local_frozen_validation",
    )
    base = withenv(unset...) do
        Dict(id => fit_key(id) for id in ids)
    end

    ## The same settings give the same key, and so does a default set
    ## explicitly.
    withenv(unset...) do
        @test fit_key("joint") == base["joint"]
    end
    withenv(
        "BVD_JOINT_SAMPLES" => "1000", "BVD_JOINT_WARMUP" => "1000",
        "BVD_JOINT_TARGET_ACCEPT" => "0.80", "BVD_JOINT_MAX_DEPTH" => "10"
    ) do
        @test fit_key("joint") == base["joint"]
        @test fit_key("sens_no_patches") == base["sens_no_patches"]
    end

    ## Each override moves every fit's key, since every fit samples at the
    ## joint's settings.
    overrides = (
        "BVD_JOINT_SAMPLES" => "1200", "BVD_JOINT_WARMUP" => "400",
        "BVD_JOINT_TARGET_ACCEPT" => "0.70", "BVD_JOINT_MAX_DEPTH" => "12",
    )
    for ov in overrides
        withenv(unset..., ov) do
            for id in ids
                @test fit_key(id) != base[id]
            end
        end
    end
end
