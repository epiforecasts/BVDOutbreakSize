@testitem "joint sampler budget defaults and overrides" tags = [:quality] begin
    include(joinpath(@__DIR__, "..", "docs", "fits", "registry.jl"))

    unset = (
        "BVD_JOINT_SAMPLES" => nothing, "BVD_JOINT_WARMUP" => nothing,
        "BVD_JOINT_TARGET_ACCEPT" => nothing, "BVD_JOINT_MAX_DEPTH" => nothing,
    )
    withenv(unset...) do
        s = joint_sampler_args()
        @test s.samples == 1000
        @test s.n_adapts == 500
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

@testitem "the joint-model fits share one sampler budget" tags = [
    :quality,
] begin
    ## The headline and its control are the halves of the spatial
    ## sensitivity, and the validation and sensitivity joints are read
    ## against them, so all of them draw from the same budget. Each splats
    ## the one helper; a fit that spells its own `samples` or `n_adapts` out
    ## would drift from the others.
    src = read(
        joinpath(@__DIR__, "..", "docs", "fits", "registry.jl"), String
    )
    @test count("joint_sampler_args()...", src) == 3
    @test count("budget = patches ? joint_sampler_args()", src) == 1
    @test count("n_adapts = joint_warmup(", src) == 1
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
    base = withenv(unset...) do
        Dict(
            id => fit_key(id) for id in (
                    "joint", "sens_no_patches", "frozen_validation",
                    "sens_community_delay", "sens_exp_growth_clock", "deaths",
                )
        )
    end

    @test issubset(
        (
            "joint", "sens_no_patches", "frozen_validation",
            "sens_community_delay", "sens_exp_growth_clock",
        ),
        JOINT_SAMPLER_FITS
    )

    ## The same settings give the same key, and so does a default set
    ## explicitly.
    withenv(unset...) do
        @test fit_key("joint") == base["joint"]
    end
    withenv(
        "BVD_JOINT_SAMPLES" => "1000", "BVD_JOINT_WARMUP" => "500",
        "BVD_JOINT_TARGET_ACCEPT" => "0.80", "BVD_JOINT_MAX_DEPTH" => "10"
    ) do
        @test fit_key("joint") == base["joint"]
        @test fit_key("sens_no_patches") == base["sens_no_patches"]
    end

    ## Each override moves every joint-budget key and leaves the other fits
    ## alone.
    overrides = (
        "BVD_JOINT_SAMPLES" => "1200", "BVD_JOINT_WARMUP" => "400",
        "BVD_JOINT_TARGET_ACCEPT" => "0.70", "BVD_JOINT_MAX_DEPTH" => "12",
    )
    for ov in overrides
        withenv(unset..., ov) do
            for id in JOINT_SAMPLER_FITS
                @test fit_key(id) != base[id]
            end
            @test fit_key("deaths") == base["deaths"]
        end
    end
end
