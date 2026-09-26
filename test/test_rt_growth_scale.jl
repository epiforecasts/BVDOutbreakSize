@testitem "rt_walk_model: gi_scale scales the walk on the log-R scale" begin
    using BVDOutbreakSize: rt_walk_model, knot_days
    using Turing: fix
    using Random: Xoshiro

    ## To first order `log R ≈ r G`, so the same growth-rate path needs
    ## log-R steps in proportion to the generation interval.
    n = 70
    nb = length(knot_days(n; start = 1))
    z = randn(Xoshiro(1), nb - 1)
    vals = (; sigma_rw = 0.1, z, intervention_effect = 0.0)
    base = log(1.5)
    logR(s) = log.(fix(rt_walk_model(n, base; gi_scale = s), vals)().Rt)
    @test logR(1.0) .- base ≈ 0.5 .* (logR(2.0) .- base)
    @test logR(0.5) .- base ≈ 0.5 .* (logR(1.0) .- base)
end

@testitem "patch_rt_model: gi_scale scales the provincial deviations" begin
    using BVDOutbreakSize: patch_rt_model
    using Turing: DynamicPPL, sample, Prior
    using Turing.DynamicPPL: returned
    using Random: Xoshiro

    ## The innovation covariance has trace `(gi_scale σ_drift)² (np - 1)`.
    n, np = 120, 4
    for s in (0.5, 2.0)
        m = patch_rt_model(
            n, np, log(1.5); rt_start = 20, breakpoint = 60.0, gi_scale = s
        )
        chn = sample(Xoshiro(3), m, Prior(), 50; progress = false)
        σs = vec(chn[DynamicPPL.VarName{:σ_drift}()])
        for (r, σ) in zip(vec(returned(m, chn)), σs)
            @test sum(abs2, r.drift_factor) ≈ (s * σ)^2 * (np - 1) rtol = 1.0e-10
        end
    end
end
