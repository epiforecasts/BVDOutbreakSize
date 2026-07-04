## Tests for the rule-out / confirmation-process sensitivity variant: the
## confirmed stream re-fit with a higher, tighter effective sensitivity prior
## (Beta(38, 2), mean 0.95) reflecting the repeat-control confirmation process
## rather than a single analytical assay draw (issue #374). The variant uses
## only the existing `sensitivity` hook of `confirmed_cases_model` and the
## `confirmed` override of `bvd_joint`, so it needs no headline-model change.

@testitem "process-sensitivity prior is higher and tighter than single-assay" tags=[:slow] begin
    using Turing: sample, Prior
    using Random: MersenneTwister
    using Statistics: mean, std
    using Distributions: Beta
    using BVDOutbreakSize: test_sensitivity_model

    base = sample(MersenneTwister(1), test_sensitivity_model(),
        Prior(), 4_000; progress = false)
    proc = sample(MersenneTwister(1),
        test_sensitivity_model(sensitivity_prior = Beta(38.0, 2.0)),
        Prior(), 4_000; progress = false)
    s_base = vec(Array(base[:s_test]))
    s_proc = vec(Array(proc[:s_test]))
    ## Higher mean, tighter spread, and still a valid probability.
    @test mean(s_proc) > mean(s_base)
    @test std(s_proc) < std(s_base)
    @test all(0 .< s_proc .< 1)
end

@testitem "confirmed override with process sensitivity runs in the joint" tags=[:slow] begin
    using Turing: sample, Prior
    import FlexiChains
    using Distributions: Beta
    using BVDOutbreakSize: load_observations, bvd_joint, confirmed_cases_model,
                           test_sensitivity_model, genetic_seeding_model

    obs = load_observations()
    breakpoint = obs.n - obs.who_first_sitrep_days

    ## Same override the registry's sens_ruleout_sensitivity fit uses: inject a
    ## higher process-sensitivity prior into the confirmed stream.
    confirmed_process = (history, cc, onsets, k, p_drc, bg, tau, bvd;
        kwargs...) -> confirmed_cases_model(
        history, cc, onsets, k, p_drc, bg, tau, bvd;
        sensitivity = test_sensitivity_model(
            sensitivity_prior = Beta(38.0, 2.0)),
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
        confirmed = confirmed_process,
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
