## Per-province suspected-case detection in the joint: a zero spread is the
## model without it, the contrast reaches the likelihood, and the chain
## records each province's probability.

@testsnippet ProvinceDetectionJoint begin
    using BVDOutbreakSize
    using Turing: DynamicPPL, fix, @varname
    using Random: Xoshiro

    obs = load_observations()
    ## The production inputs, with or without the province care splits and
    ## the detection contrast.
    function detection_joint(; care = true, detection = true)
        spatial = patch_fit_args(obs)
        if !care
            spatial = merge(
                spatial, (; province_isolation = nothing, province_capacity = nothing)
            )
        end
        return bvd_joint(
            obs.n, obs.exported_cases, obs.total_deaths,
            obs.reported_cases, obs.exports_deaths,
            obs.confirmed_cases, obs.tests_analysed;
            joint_fit_args(obs; breakpoint = default_breakpoint(obs))...,
            spatial...,
            (detection ? (;) : (; suspected_ascertainment = nothing))...
        )
    end
    no_spread(z) = Dict(
        @varname(detection_state.τ_susp) => 0.0,
        @varname(detection_state.z_susp) => z
    )
end

@testitem "bvd_joint: a zero detection spread is the model without it" setup = [
    ProvinceDetectionJoint,
] begin
    z = [0.8, -1.3, 0.4]
    ## Without the province care splits the detection enters only through
    ## the national suspected cases and the case and laboratory
    ## compositions, where a common probability is `p_drc`, so the density
    ## matches at any confirmation contrast.
    off = detection_joint(; care = false, detection = false)
    vi = DynamicPPL.VarInfo(Xoshiro(11), off)
    @test !any(contains("detection_state"), string.(keys(vi)))
    on = fix(detection_joint(; care = false), no_spread(z))
    lp = DynamicPPL.logjoint(off, vi)
    @test isfinite(lp)
    @test DynamicPPL.logjoint(on, vi) ≈ lp rtol = 1.0e-10

    ## With the care splits the model without detection also splits the
    ## occupancy by the confirmation contrast. With that contrast at zero
    ## the two models agree.
    flat = Dict(@varname(composition_state.τ_asc) => 0.0)
    off_care = fix(detection_joint(; detection = false), flat)
    vi_care = DynamicPPL.VarInfo(Xoshiro(12), off_care)
    on_care = fix(detection_joint(), merge(no_spread(z), flat))
    lp_care = DynamicPPL.logjoint(off_care, vi_care)
    @test DynamicPPL.logjoint(on_care, vi_care) ≈ lp_care rtol = 1.0e-10
    ## A confirmation contrast moves the occupancy split only without the
    ## detection contrast.
    tilt = Dict(@varname(composition_state.τ_asc) => 0.5)
    @test !isapprox(
        DynamicPPL.logjoint(
            fix(detection_joint(; detection = false), tilt), vi_care
        ),
        DynamicPPL.logjoint(
            fix(detection_joint(), merge(no_spread(z), tilt)), vi_care
        );
        rtol = 1.0e-8
    )
    ## A detection spread reaches the likelihood.
    spread = Dict(
        @varname(detection_state.τ_susp) => 0.5,
        @varname(detection_state.z_susp) => z
    )
    @test !isapprox(
        DynamicPPL.logjoint(fix(detection_joint(; care = false), spread), vi),
        lp; rtol = 1.0e-8
    )
end

@testitem "bvd_joint: records each province's detection and ascertainment" setup = [
    ProvinceDetectionJoint,
] begin
    using BVDOutbreakSize: _draw_vectors, _draws
    using Turing: sample, Prior
    using StatsFuns: logit
    import FlexiChains
    chn = sample(
        Xoshiro(13), detection_joint(), Prior(), 5;
        chain_type = FlexiChains.VNChain, progress = false
    )
    p_drc = vec(_draws(chn, :p_drc))
    p = _draw_vectors(chn, :province_suspected_ascertainment)
    conf = _draw_vectors(chn, :province_confirmation_relative)
    asc = _draw_vectors(chn, :province_ascertainment)
    np = length(PROVINCE_NAMES)
    for i in eachindex(p)
        @test length(p[i]) == np
        @test all(a -> 0 < a < 1, p[i])
        ## A logit contrast about the national centre, summing to zero.
        @test sum(logit.(p[i]) .- logit(p_drc[i])) ≈ 0 atol = 1.0e-8
        ## The relative case ascertainment is detection times confirmation,
        ## centred on the log scale.
        @test sum(log, asc[i]) ≈ 0 atol = 1.0e-10
        @test asc[i] ./ asc[i][1] ≈ (p[i] .* conf[i]) ./ (p[i][1] * conf[i][1])
    end
end
