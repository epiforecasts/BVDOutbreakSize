@testitem "test_sensitivity_model samples a probability" begin
    using BVDOutbreakSize: test_sensitivity_model
    using Random: seed!
    ## A submodel's `(; ...)` return tuple (with the derived fields) comes
    ## from calling it; `rand(rng, model)` returns only sampled variables.
    seed!(1)
    s = test_sensitivity_model()().s_test
    @test 0 <= s <= 1
end

@testitem "lab_delay_model returns a normalised daily PMF" begin
    using BVDOutbreakSize: lab_delay_model
    using Random: seed!
    seed!(1)
    d = lab_delay_model(20)()
    @test all(>=(0), d.pmf)
    @test isapprox(sum(d.pmf), 1; atol = 1e-8)
    @test d.mean > 0
end

@testitem "test_positivity_model samples background and testing fraction" begin
    using BVDOutbreakSize: test_positivity_model
    using Random: seed!
    seed!(1)
    s = test_positivity_model()()
    @test s.λ_bg >= 0
    @test 0 <= s.τ_test <= 1
end

@testitem "confirmed_only_model conditions on the lab pipeline" begin
    using BVDOutbreakSize: confirmed_only_model
    using Turing.DynamicPPL: logjoint
    using Random: MersenneTwister

    m = confirmed_only_model(40, 27;
        confirmed_history = (; days = [18, 40], counts = [17, 27]),
        lab_history = (; days = [18, 40], counts = [12, 28]))
    draw = rand(MersenneTwister(1), m)
    @test isfinite(logjoint(m, draw))
end

@testitem "bvd_joint exposes lab-pipeline deterministics" tags=[:slow] begin
    using BVDOutbreakSize: bvd_joint, nuts_sample, load_observations
    using Statistics: mean

    obs = load_observations()
    m = bvd_joint(
        obs.n, obs.exported_cases, obs.total_deaths,
        obs.reported_cases, obs.exports_deaths, obs.confirmed_cases,
        obs.tests_analysed;
        confirmed_deaths = obs.confirmed_deaths,
        deaths_history = obs.deaths_history,
        reported_history = obs.reported_history,
        confirmed_history = obs.confirmed_history,
        confirmed_deaths_history = obs.confirmed_deaths_history,
        lab_history = obs.lab_history,
        lab_daily_history = obs.lab_daily_history,
        breakpoint = obs.n - obs.who_first_sitrep_days,
        tmrca_days = obs.tmrca_days)
    chn = nuts_sample(m; samples = 25, chains = 1, progress = false)
    for key in (:expected_confirmed_T, :expected_analysed_T,
        :tau_test, :lambda_bg, :suspected_positivity, :test_positivity,
        :death_ascertainment, :background_cfr, :tau_death,
        :death_composition, :death_confirmation,
        :expected_confirmed_deaths_T)
        v = vec(Array(chn[key]))
        @test all(isfinite, v)
    end
    @test all(0 .<= vec(Array(chn[:test_positivity])) .<= 1)
    @test all(0 .<= vec(Array(chn[:death_confirmation])) .<= 1)
    @test all(0 .<= vec(Array(chn[:death_composition])) .<= 1)
    @test all(0 .<= vec(Array(chn[:death_ascertainment])) .<= 1)
    ## `tau_death` here is the realised cut-off death-testing intensity
    ## (scaling * analysed / suspected), surfaced as a diagnostic. The death
    ## volume is capped at the suspected-death pool each day, so the realised
    ## intensity stays in [0, 1] (it can sit at 1 on a backlog day, where the
    ## uncapped intensity would have exceeded the suspected pool).
    @test all(0 .<= vec(Array(chn[:tau_death])) .<= 1)
end
