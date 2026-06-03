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

    @model function _xc_only(exported_cases_daily::AbstractVector)
        growth_state ~ to_submodel(_xc_growth(), false)
        asc_state ~ to_submodel(_xc_pooled(), false)
        incubation_state ~ to_submodel(incubation_model(), false)
        report_state ~ to_submodel(report_delay_model(), false)
        f_det = combined_delay(incubation_state.dist, report_state.dist)

        exports_state ~ to_submodel(
            exports_daily_delay_model(exported_cases_daily, growth_state,
                asc_state.p_uganda, f_det), false)

        cumulative_cases := growth_state.C_T
    end
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
