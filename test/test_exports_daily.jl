## Tests for the dated export-case likelihood. The real
## `exports_daily_delay_model` places the cumulative export expectation
## `expected_exports_delay` onto a daily detection grid: a continuous
## survival term over the pre-detection stretch followed by per-day
## Poisson counts, the same construction as `exports_deaths_delay_model`.
## The pre-detection survival term subsumes the first-export-detection
## timing bound, and the total expected (pre + daily increments) recovers
## the single-total cumulative expectation at the cut-off.

@testsnippet ExportsDailyFixtures begin
    using Distributions: Gamma, Normal, truncated
    using StatsFuns: logit, logistic
    using Turing: @model, to_submodel
    using BVDOutbreakSize: exports_daily_delay_model, combined_delay,
                           incubation_model, report_delay_model

    @model function _xc_growth()
        log_τ ~ Normal(log(14), 0.4)
        m ~ truncated(Normal(7.0, 2.5); lower = 0, upper = 13.0)
        τ := exp(log_τ)
        r := log(2) / τ
        T := m * τ
        C_T := 2.0^m
        cumulative = s -> exp(r * s)
        return (; log_τ, τ, r, m, T, C_T, cumulative)
    end

    @model function _xc_pooled(;
            mu_prior = Normal(logit(0.25), 1.0),
            tau_prior = truncated(Normal(0.0, 0.5); lower = 1e-4))
        μ_logit ~ mu_prior
        τ_logit ~ tau_prior
        z_uganda ~ Normal(0, 1)
        logit_p_uganda = μ_logit + τ_logit * z_uganda
        p_uganda := logistic(logit_p_uganda)
        return (; p_uganda)
    end

    @model function _xc_only(exported_cases_daily::AbstractVector;
            last_offset = 0)
        growth_state ~ to_submodel(_xc_growth(), false)
        asc_state ~ to_submodel(_xc_pooled(), false)
        incubation_state ~ to_submodel(incubation_model(), false)
        report_state ~ to_submodel(report_delay_model(), false)
        f_det = combined_delay(incubation_state.dist, report_state.dist)

        exports_state ~ to_submodel(
            exports_daily_delay_model(exported_cases_daily, growth_state,
                asc_state.p_uganda, f_det; last_offset = last_offset), false)

        cumulative_cases := growth_state.C_T
    end
end

@testitem "export_at_risk matches expected_exports_delay" begin
    using Distributions: Gamma
    using BVDOutbreakSize: ExportRiskTrajectory, export_at_risk,
                           expected_exports_delay

    r = 0.05
    f_det = Gamma(4.3, 2.6)
    edges = [20.0, 45.0, 70.0, 90.0]
    traj = ExportRiskTrajectory(maximum(edges), r)
    got = export_at_risk(traj, edges, f_det)
    ## p = q = 1 so the wrapper returns the bare at-risk person-time.
    for (i, t) in enumerate(edges)
        want = expected_exports_delay(r, 1.0, 1.0, t, f_det)
        @test got[i] ≈ want rtol = 1e-3
    end
end

@testitem "export_death_at_risk matches expected_exports_deaths_delay" begin
    using Distributions: Gamma
    using BVDOutbreakSize: ExportRiskTrajectory, export_death_at_risk,
                           expected_exports_deaths_delay

    r = 0.05
    cumulative = s -> exp(r * s)
    f_det = Gamma(4.3, 2.6)
    f_death = Gamma(4.3, 2.6)
    edges = [30.0, 60.0, 90.0]
    traj = ExportRiskTrajectory(maximum(edges), r)
    got = export_death_at_risk(traj, edges, f_det, f_death)
    for (i, t) in enumerate(edges)
        want = expected_exports_deaths_delay(cumulative, f_det, f_death,
            1.0, 1.0, 1.0, t)
        @test got[i] ≈ want rtol = 1e-2
    end
end

@testitem "export at-risk engines stay finite at extreme growth" begin
    using Distributions: Gamma
    using BVDOutbreakSize: ExportRiskTrajectory, export_at_risk,
                           export_death_at_risk

    f_det = Gamma(4.3, 2.6)
    f_death = Gamma(4.3, 2.6)
    edges = [30.0, 60.0, 90.0]
    ## A NUTS proposal can drive `r` high enough that exp(r·t) overflows;
    ## the growth integral and the removed convolution both go to Inf, so
    ## their difference is NaN. The engines must clamp, or the NaN reaches
    ## Poisson and throws under AD.
    for r in (0.5, 1.0, 2.0)
        traj = ExportRiskTrajectory(90.0, r)
        a = export_at_risk(traj, edges, f_det)
        d = export_death_at_risk(traj, edges, f_det, f_death)
        @test all(isfinite, a) && all(>=(0.0), a)
        @test all(isfinite, d) && all(>=(0.0), d)
    end
end

@testitem "export_at_risk is reverse-mode differentiable" tags=[:ad] begin
    using Distributions: Gamma
    using Mooncake: Mooncake
    using BVDOutbreakSize: ExportRiskTrajectory, export_at_risk

    edges = [20.0, 45.0, 70.0, 90.0]
    f(r) = begin
        traj = ExportRiskTrajectory(90.0, r)
        sum(export_at_risk(traj, edges, Gamma(4.3, 2.6)))
    end
    cache = Mooncake.prepare_gradient_cache(f, 0.05)
    g = Mooncake.value_and_gradient!!(cache, f, 0.05)[2][2]
    @test isfinite(g)
    @test g > 0   # more growth ⇒ more at-risk person-time
end

@testitem "exports daily means partition the cumulative expectation" begin
    using Distributions: Gamma
    using BVDOutbreakSize: expected_exports_delay, daily_increment_kernel

    r = 0.05
    p_uganda = 0.25
    q = 1871 / 4_392_200
    f_det = Gamma(4.3, 2.6)
    T = 90.0
    n = 16   # daily series from earliest detection to the cut-off

    Λ(t) = expected_exports_delay(r, p_uganda, q, t, f_det)
    pre = Λ(T - n)
    Λ_at_edges = [Λ(T - n + i) for i in 1:n]
    μ_day = daily_increment_kernel(Λ_at_edges, pre)

    ## Total person-time is conserved: the pre-detection survival weight
    ## plus the per-day increments recover the cumulative expectation at
    ## the cut-off, the quantity the scalar single-total likelihood uses.
    @test pre + sum(μ_day) ≈ Λ(T) rtol = 1e-10

    ## With a last_offset the series ends at t_last = T - last_offset, so
    ## the conserved total is the cumulative expectation at t_last, not T.
    last_offset = 3.0
    t_last = T - last_offset
    pre_l = Λ(t_last - n)
    edges_l = [Λ(t_last - n + i) for i in 1:n]
    μ_l = daily_increment_kernel(edges_l, pre_l)
    @test pre_l + sum(μ_l) ≈ Λ(t_last) rtol = 1e-10
    @test Λ(t_last) < Λ(T)
    @test all(μ_day .> 0)
    @test pre > 0
end

@testitem "exports_daily_delay_model prior draws produce valid counts" tags=[:slow] setup=[ExportsDailyFixtures] begin
    using Turing: sample, Prior
    import FlexiChains

    ## Earliest detection 15 days before the cut-off, then two more
    ## detections; zeros elsewhere (the Uganda import pattern).
    daily=[1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0]
    chn=sample(_xc_only(daily), Prior(), 200;
        chain_type = FlexiChains.VNChain, progress = false)

    expected=vec(Array(chn[:expected_exports_T]))
    @test length(expected) == 200
    @test all(isfinite, expected)
    @test all(expected .> 0)
end

@testitem "exports_daily_delay_model fits an all-zero series" tags=[:slow] setup=[ExportsDailyFixtures] begin
    using Turing: sample, Prior
    import FlexiChains

    chn=sample(_xc_only(zeros(Int, 16)), Prior(), 200;
        chain_type = FlexiChains.VNChain, progress = false)
    C=vec(Array(chn[:cumulative_cases]))
    @test length(C) == 200
    @test all(isfinite, C)
    @test all(C .> 0)
end

@testitem "exports_daily_delay_model last_offset stops the model early" tags=[:slow] setup=[ExportsDailyFixtures] begin
    using Turing: sample, Prior
    using Random: MersenneTwister
    import FlexiChains

    ## The same detection series fit to the cut-off vs stopped three days
    ## early (last reported import). Stopping early evaluates the expected
    ## exports at the earlier date, so the per-draw expectation is lower.
    ## Same seed pairs the parameter draws, so the only difference is the
    ## anchor date and the comparison is exact per draw (the expectation is
    ## heavy-tailed, so unpaired means are too noisy to compare).
    daily=[1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1]
    chn0=sample(MersenneTwister(1), _xc_only(daily; last_offset = 0),
        Prior(), 200; chain_type = FlexiChains.VNChain, progress = false)
    chn3=sample(MersenneTwister(1), _xc_only(daily; last_offset = 3),
        Prior(), 200; chain_type = FlexiChains.VNChain, progress = false)
    e0=vec(Array(chn0[:expected_exports_T]))
    e3=vec(Array(chn3[:expected_exports_T]))
    @test all(e3 .< e0)
end

@testitem "bvd_joint uses the dated export series when provided" tags=[:slow] begin
    using Turing: sample, Prior
    using BVDOutbreakSize: bvd_joint
    import FlexiChains

    ## With a dated detection series the joint switches to the
    ## time-resolved export likelihood. Pass a non-empty daily series (the
    ## scalar count is then ignored) and a first-detection delta: the
    ## dated likelihood's pre-detection term carries the timing bound, so
    ## the separate timing term must be a no-op and not double-count it.
    model=bvd_joint(missing, [120], [500], [0, 0, 1];
        exported_cases_daily = [1, 0, 0, 1],
        reported_offsets = [0],
        death_offsets = [0],
        first_export_detection_delta = 3)

    chn=sample(model, Prior(), 50;
        chain_type = FlexiChains.VNChain, progress = false)

    expected=vec(Array(chn[:expected_exports_T]))
    @test length(expected) == 50
    @test all(isfinite, expected)
    @test all(expected .> 0)
    cases=vec(Array(chn[:cumulative_cases]))
    @test all(isfinite, cases) && all(>(0), cases)
end
