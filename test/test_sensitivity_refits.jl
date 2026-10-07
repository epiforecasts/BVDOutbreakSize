@testitem "sensitivity re-fits differ from the headline only in overrides" tags = [
    :quality,
] begin
    include(joinpath(@__DIR__, "..", "docs", "fits", "registry.jl"))

    obs = load_observations()
    headline = headline_joint_args(obs)
    @test headline == merge(
        joint_fit_args(obs; breakpoint = default_breakpoint(obs)),
        patch_fit_args(obs)
    )
    @test headline.n_patches == length(PROVINCE_NAMES)

    overrides = sensitivity_overrides(obs)
    @test keys(overrides) == (:sens_community_delay, :sens_exp_growth_clock)
    for ov in overrides
        variant = headline_joint_args(obs; ov...)
        @test issetequal(keys(variant), union(keys(headline), keys(ov)))
        for k in keys(variant)
            @test isequal(variant[k], k in keys(ov) ? ov[k] : headline[k])
        end
    end
    @test overrides.sens_exp_growth_clock.tmrca_days != headline.tmrca_days

    withenv("BVD_RUN_SENSITIVITY" => "true") do
        ids = fit_ids(obs)
        for id in string.(keys(overrides))
            @test id in ids
            @test id in JOINT_SAMPLER_FITS
        end
    end
end
