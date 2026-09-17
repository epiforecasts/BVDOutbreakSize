## Tests for the recovered-among-confirmed stream ("cumul guéris"), an
## incidence (scaled-convolution) stream: survivors among the modelled daily
## confirmed cases, scaled by the recovery probability and lagged by a
## confirmation-to-recovery delay. The submodel is sampled standalone over a
## fixed daily confirmed series, and the stream is exercised through
## `bvd_joint`.

@testitem "recovered: conditioned fit tracks a finite recovered total" tags=[
    :slow
] begin
    using Turing: sample, Prior
    import FlexiChains
    using BVDOutbreakSize: recovered_model

    ## A fixed daily confirmed incidence and a cumulative recovered history
    ## on later days, supplied as observed counts.
    confirmed_daily = fill(8.0, 33)
    recovered_history = (; days = [30, 31, 32, 33], counts = [12, 19, 22, 30])
    chn = sample(
        recovered_model(recovered_history, 30, confirmed_daily, 0.3),
        Prior(), 100;
        chain_type = FlexiChains.VNChain, progress = false
    )
    er = vec(Array(chn[:expected_recovered]))
    ## `p_recover` is sampled inside the injected recovery submodel, so it is
    ## keyed under that submodel's prefix.
    pr = vec(Array(chn[Symbol("rec_state.p_recover")]))
    @test length(er) == 100
    @test all(isfinite, er)
    @test all(er .>= 0)
    ## The recovery probability is a fraction.
    @test all(0 .<= pr .<= 1)
end

@testitem "recovered: predictive path samples the increments" tags=[:slow] begin
    using Turing: sample, Prior
    import FlexiChains
    using BVDOutbreakSize: recovered_model

    ## Days but no counts: the increments are sampled under the
    ## `recovered_increments` submodel rather than conditioned.
    confirmed_daily = fill(8.0, 33)
    recovered_history = (; days = [30, 31, 32, 33], counts = Int[])
    chn = sample(
        recovered_model(recovered_history, missing, confirmed_daily, 0.3),
        Prior(), 50;
        chain_type = FlexiChains.VNChain, progress = false
    )
    ks = string.(collect(keys(chn)))
    @test any(k -> occursin("recovered_increments", k), ks)
end

@testitem "recovered: empty history is a no-op" tags=[:slow] begin
    using Turing: sample, Prior
    import FlexiChains
    using BVDOutbreakSize: recovered_model

    confirmed_daily = fill(8.0, 33)
    chn = sample(
        recovered_model((; days = Int[], counts = Int[]), missing,
            confirmed_daily, 0.3),
        Prior(), 50;
        chain_type = FlexiChains.VNChain, progress = false
    )
    ks = string.(collect(keys(chn)))
    @test !any(k -> occursin("recovered_increments.increments", k), ks)
end

@testitem "recovered: joint prior runs with the live data" tags=[:slow] begin
    using Turing: sample, Prior
    import FlexiChains
    using BVDOutbreakSize: load_observations, bvd_joint, genetic_seeding_model

    obs = load_observations()
    @test !isempty(obs.recovered_history.counts)
    breakpoint = obs.n - obs.who_first_sitrep_days
    m = bvd_joint(obs.n, obs.exported_cases, obs.total_deaths,
        obs.reported_cases, obs.exports_deaths, obs.confirmed_cases,
        obs.tests_analysed;
        confirmed_deaths = obs.confirmed_deaths,
        recovered_cases = obs.recovered_cases,
        deaths_history = obs.deaths_history,
        reported_history = obs.reported_history,
        confirmed_history = obs.confirmed_history,
        confirmed_deaths_history = obs.confirmed_deaths_history,
        lab_history = obs.lab_history,
        lab_daily_history = obs.lab_daily_history,
        suspected_daily_history = obs.suspected_daily_history,
        isolation_history = obs.isolation_history,
        recovered_history = obs.recovered_history,
        export_case_days = obs.export_case_days,
        export_death_days = obs.export_death_days,
        breakpoint = breakpoint,
        genetic = genetic_seeding_model,
        tmrca_days = obs.tmrca_days)
    chn = sample(m, Prior(), 30;
        chain_type = FlexiChains.VNChain, progress = false)
    rec = vec(Array(chn[:expected_recovered_T]))
    @test length(rec) == 30
    @test all(isfinite, rec)
    @test all(rec .>= 0)
end
