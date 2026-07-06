## Tests for the opt-in suspected-case reporting-effort multiplier
## ([`reporting_effort_walk_model`](@ref)) and its wiring into the suspected-case
## likelihood ([`reported_cases_model`](@ref)) and the joint composer.

@testitem "reporting effort: multiplier is one before onset and positive" begin
    using BVDOutbreakSize: reporting_effort_walk_model

    n = 40
    onset = 15
    ## Calling the model runs it under the prior and returns the effort series.
    ret = reporting_effort_walk_model(n; onset)()
    @test length(ret.effort) == n
    @test all(ret.effort .> 0)
    ## Before the window onset the multiplier is exactly one (a no-op there).
    @test all(ret.effort[1:(onset - 1)] .== 1)
end

@testitem "reporting effort: sigma to zero recovers constant effort one" begin
    using BVDOutbreakSize: reporting_effort_walk_model
    using Distributions: Dirac

    n = 30
    onset = 8
    ## Pin the innovation SD to zero: the effort must be exactly one over the
    ## whole window, whatever the innovations z are.
    ret = reporting_effort_walk_model(n; onset, sigma_prior = Dirac(0.0))()
    @test length(ret.effort) == n
    @test all(isapprox.(ret.effort, 1.0; atol = 1e-10))
end

@testitem "reporting effort: joint flag off leaves the likelihood unchanged" tags=[:slow] begin
    using Turing: sample, Prior
    using Random: MersenneTwister
    import FlexiChains
    using BVDOutbreakSize: load_observations, bvd_joint, genetic_seeding_model

    obs = load_observations()
    breakpoint = obs.n - obs.who_first_sitrep_days
    mk(flag) = bvd_joint(obs.n, obs.exported_cases, obs.total_deaths,
        obs.reported_cases, obs.exports_deaths, obs.confirmed_cases,
        obs.tests_analysed;
        confirmed_deaths = obs.confirmed_deaths,
        deaths_history = obs.deaths_history,
        reported_history = obs.reported_history,
        confirmed_history = obs.confirmed_history,
        confirmed_deaths_history = obs.confirmed_deaths_history,
        lab_history = obs.lab_history,
        lab_daily_history = obs.lab_daily_history,
        suspected_daily_history = obs.suspected_daily_history,
        isolation_history = obs.isolation_history,
        bed_capacity_history = obs.bed_capacity_history,
        export_case_days = obs.export_case_days,
        export_death_days = obs.export_death_days,
        breakpoint = breakpoint, background_re = true,
        suspected_reporting_effort = flag,
        genetic = genetic_seeding_model, tmrca_days = obs.tmrca_days)

    ks_off = string.(collect(keys(sample(MersenneTwister(1), mk(false),
        Prior(), 5; chain_type = FlexiChains.VNChain, progress = false))))
    ks_on = string.(collect(keys(sample(MersenneTwister(1), mk(true),
        Prior(), 5; chain_type = FlexiChains.VNChain, progress = false))))
    ## Off: no effort parameters in the chain. On: the effort walk adds them.
    @test !any(k -> occursin("effort", k), ks_off)
    @test any(k -> occursin("effort", k) || occursin("σ_eff", k), ks_on)
end

@testitem "reporting effort: joint flag on stays finite" tags=[:slow] begin
    using Turing: sample, Prior
    import FlexiChains
    using BVDOutbreakSize: load_observations, bvd_joint, genetic_seeding_model

    obs = load_observations()
    breakpoint = obs.n - obs.who_first_sitrep_days
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
        suspected_daily_history = obs.suspected_daily_history,
        export_case_days = obs.export_case_days,
        export_death_days = obs.export_death_days,
        breakpoint = breakpoint, background_re = true,
        suspected_reporting_effort = true,
        genetic = genetic_seeding_model, tmrca_days = obs.tmrca_days)
    chn = sample(m, Prior(), 20;
        chain_type = FlexiChains.VNChain, progress = false)
    C_T = vec(Array(chn[:C_T]))
    R_T = vec(Array(chn[:R_T]))
    @test all(isfinite, C_T) && all(C_T .> 0)
    @test all(isfinite, R_T)
end
