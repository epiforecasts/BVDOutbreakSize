## Tests for the confirmation-process sensitivity integration and its
## single-assay downside variant. The headline `test_sensitivity_model` now
## credits the repeat-control confirmation process (Beta(38, 2), mean 0.95)
## rather than one analytical assay draw (issue #374). The downside variant
## re-fits the confirmed stream with the single-assay prior (Beta(10, 1.76),
## mean 0.85) through the existing `sensitivity` hook of `confirmed_cases_model`
## and the `confirmed` override of `bvd_joint`.

@testitem "headline process prior is higher and tighter than single-assay" tags=[:slow] begin
    using Turing: sample, Prior
    using Random: MersenneTwister
    using Statistics: mean, std
    using Distributions: Beta
    using BVDOutbreakSize: test_sensitivity_model

    headline = sample(MersenneTwister(1), test_sensitivity_model(),
        Prior(), 4_000; progress = false)
    single = sample(MersenneTwister(1),
        test_sensitivity_model(sensitivity_prior = Beta(10.0, 1.76)),
        Prior(), 4_000; progress = false)
    s_headline = vec(Array(headline[:s_test]))
    s_single = vec(Array(single[:s_test]))
    ## The headline default credits the process: higher mean, tighter spread,
    ## and still a valid probability.
    @test mean(s_headline) > mean(s_single)
    @test std(s_headline) < std(s_single)
    @test mean(s_headline) > 0.9
    @test all(0 .< s_headline .< 1)
end

@testitem "confirmed override with single-assay sensitivity runs in the joint" tags=[:slow] begin
    using Turing: sample, Prior
    import FlexiChains
    using Distributions: Beta
    using BVDOutbreakSize: load_observations, bvd_joint, confirmed_cases_model,
                           test_sensitivity_model, genetic_seeding_model

    obs = load_observations()
    breakpoint = obs.n - obs.who_first_sitrep_days

    ## Same override the registry's sens_ruleout_sensitivity fit uses: inject the
    ## single-assay downside prior into the confirmed stream.
    confirmed_single = (history, cc, onsets, k, p_drc, bg, tau, bvd;
        kwargs...) -> confirmed_cases_model(
        history, cc, onsets, k, p_drc, bg, tau, bvd;
        sensitivity = test_sensitivity_model(
            sensitivity_prior = Beta(10.0, 1.76)),
        kwargs...)

    m = bvd_joint(obs.n, obs.exported_cases, obs.total_deaths,
        obs.reported_cases, obs.exports_deaths, obs.confirmed_cases,
        obs.tests_analysed;
        confirmed_deaths = obs.confirmed_deaths,
        deaths_history = obs.deaths_history,
        reported_history = obs.reported_history,
        confirmed_history = obs.confirmed_history,
        confirmed_deaths_history = obs.confirmed_deaths_history,
        lab_history = obs.lab_history,
        lab_daily_history = obs.lab_daily_history,
        confirmed = confirmed_single,
        export_case_days = obs.export_case_days,
        export_death_days = obs.export_death_days,
        breakpoint = breakpoint,
        genetic = genetic_seeding_model,
        tmrca_days = obs.tmrca_days)
    chn = sample(m, Prior(), 30;
        chain_type = FlexiChains.VNChain, progress = false)
    C_T = vec(Array(chn[:C_T]))
    @test length(C_T) == 30
    @test all(isfinite, C_T)
    @test all(C_T .> 0)
end
