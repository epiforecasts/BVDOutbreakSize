@testitem "joint sampler budget defaults and overrides" tags = [:quality] begin
    include(joinpath(@__DIR__, "..", "docs", "fits", "registry.jl"))

    unset = (
        "BVD_JOINT_SAMPLES" => nothing, "BVD_JOINT_WARMUP" => nothing,
        "BVD_JOINT_TARGET_ACCEPT" => nothing,
    )
    withenv(unset...) do
        s = joint_sampler_args()
        @test s.samples == 1000
        @test s.n_adapts == 500
        @test s.target_accept == 0.80
    end
    withenv(
        "BVD_JOINT_SAMPLES" => "1200", "BVD_JOINT_WARMUP" => "400",
        "BVD_JOINT_TARGET_ACCEPT" => "0.85"
    ) do
        s = joint_sampler_args()
        @test s.samples == 1200
        @test s.n_adapts == 400
        @test s.target_accept == 0.85
    end
end

@testitem "the headline joint and its control share one sampler budget" tags = [
    :quality,
] begin
    ## The two fits are the halves of the spatial sensitivity, so they must
    ## draw from the same budget. Both splat the one helper; a fit that spells
    ## its own `samples` or `n_adapts` out would drift from the other.
    src = read(
        joinpath(@__DIR__, "..", "docs", "fits", "registry.jl"), String
    )
    @test count("joint_sampler_args()...", src) == 2
    @test !occursin("n_adapts = joint_warmup(", src)
end
