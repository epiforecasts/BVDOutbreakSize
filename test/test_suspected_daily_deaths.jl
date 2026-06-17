## Tests for the post-26 May daily new suspected deaths ("cas suspects du
## jour N (M deces)"), scored against the modelled daily suspected-death
## series at each report day, the deaths analogue of the daily new-suspect
## inflow. Exercised through `deaths_only_model` and `bvd_joint`.

@testitem "suspected daily deaths: conditioned fit stays positive" tags=[:slow] begin
    using Turing: sample, Prior
    import FlexiChains
    using BVDOutbreakSize: deaths_only_model

    ## Cumulative suspected-death history (frozen) plus the disjoint daily
    ## inflow on later days, supplied as observed counts.
    deaths_history = (; days = [13, 18, 23], counts = [131, 204, 246])
    suspected_daily_deaths_history = (; days = [30, 31, 32, 33],
        counts = [35, 41, 30, 50])
    chn = sample(
        deaths_only_model(33, missing; deaths_history,
            suspected_daily_deaths_history),
        Prior(), 100;
        chain_type = FlexiChains.VNChain, progress = false
    )
    C_T = vec(Array(chn[:C_T]))
    @test length(C_T) == 100
    @test all(isfinite, C_T)
    @test all(C_T .> 0)
end

@testitem "suspected daily deaths: predictive path samples the counts" tags=[:slow] begin
    using Turing: sample, Prior
    import FlexiChains
    using BVDOutbreakSize: deaths_only_model

    ## Days but no counts: the daily inflow is a predictive generator, so its
    ## per-day counts are sampled under the `deaths_state.suspected_daily_deaths`
    ## submodel rather than conditioned.
    deaths_history = (; days = [13, 18, 23], counts = [131, 204, 246])
    suspected_daily_deaths_history = (; days = [30, 31, 32, 33], counts = Int[])
    chn = sample(
        deaths_only_model(33, missing; deaths_history,
            suspected_daily_deaths_history),
        Prior(), 50;
        chain_type = FlexiChains.VNChain, progress = false
    )
    ks = string.(collect(keys(chn)))
    @test any(k -> occursin("suspected_daily_deaths", k), ks)
end

@testitem "suspected daily deaths: empty history is a no-op" tags=[:slow] begin
    using Turing: sample, Prior
    import FlexiChains
    using BVDOutbreakSize: deaths_only_model

    ## With no daily history the model matches the plain deaths composer: the
    ## inflow submodel scores nothing and adds no sampled keys.
    deaths_history = (; days = [13, 18, 23], counts = [131, 204, 246])
    chn = sample(
        deaths_only_model(33, missing; deaths_history),
        Prior(), 50;
        chain_type = FlexiChains.VNChain, progress = false
    )
    ks = string.(collect(keys(chn)))
    @test !any(k -> occursin("suspected_daily_deaths.increments", k), ks)
    C_T = vec(Array(chn[:C_T]))
    @test all(isfinite, C_T)
    @test all(C_T .> 0)
end

@testitem "suspected daily deaths: joint prior runs with the live data" tags=[:slow] begin
    using Turing: sample, Prior
    import FlexiChains
    using BVDOutbreakSize: load_observations, bvd_joint, genetic_seeding_model

    obs = load_observations()
    @test !isempty(obs.suspected_daily_deaths_history.counts)
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
        suspected_daily_deaths_history = obs.suspected_daily_deaths_history,
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
