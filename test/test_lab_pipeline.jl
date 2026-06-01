@testitem "test_sensitivity_model samples a probability" begin
    using BVDOutbreakSize: test_sensitivity_model
    using Random: MersenneTwister
    s = rand(MersenneTwister(1), test_sensitivity_model()).s_test
    @test 0 <= s <= 1
end

@testitem "lab_delay_model returns a normalised daily PMF" begin
    using BVDOutbreakSize: lab_delay_model
    using Random: MersenneTwister
    d = rand(MersenneTwister(1), lab_delay_model(20))
    @test all(>=(0), d.pmf)
    @test isapprox(sum(d.pmf), 1; atol = 1e-8)
    @test d.mean > 0
end

@testitem "test_positivity_model samples background and testing fraction" begin
    using BVDOutbreakSize: test_positivity_model
    using Random: MersenneTwister
    s = rand(MersenneTwister(1), test_positivity_model())
    @test s.λ_bg >= 0
    @test 0 <= s.τ_test <= 1
end

@testitem "confirmed_only_model conditions on the lab pipeline" begin
    using BVDOutbreakSize: confirmed_only_model
    using Turing.DynamicPPL: logjoint
    using Random: MersenneTwister

    m = confirmed_only_model(40, 27;
        confirmed_history = (; days = [18, 40], counts = [17, 27]),
        lab_history = (; days = [18, 40], counts = [12, 27]),
        tests_analysed = 27)
    draw = rand(MersenneTwister(1), m)
    @test isfinite(logjoint(m, draw))
end

@testitem "bvd_joint exposes lab-pipeline deterministics" tags=[:slow] begin
    using BVDOutbreakSize: bvd_joint, nuts_sample, load_observations
    using Statistics: mean

    obs = load_observations()
    tests = obs.lab_history.counts[end]
    m = bvd_joint(
        obs.n, obs.exported_cases, obs.total_deaths,
        obs.reported_cases, obs.exports_deaths, obs.confirmed_cases, tests;
        deaths_history = obs.deaths_history,
        reported_history = obs.reported_history,
        confirmed_history = obs.confirmed_history,
        lab_history = obs.lab_history,
        breakpoint = obs.n - obs.who_first_sitrep_days,
        tmrca_days = obs.tmrca_days)
    chn = nuts_sample(m; samples = 25, chains = 1, progress = false)
    for key in (:expected_confirmed_T, :expected_tested_T, :s_test,
        :tau_test, :lambda_bg, :suspected_positivity, :test_positivity)
        v = vec(Array(chn[key]))
        @test all(isfinite, v)
    end
    @test all(0 .<= vec(Array(chn[:s_test])) .<= 1)
    @test all(0 .<= vec(Array(chn[:test_positivity])) .<= 1)
end
