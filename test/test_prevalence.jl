## Point prevalence: the state weights against hand counts under fixed
## delays, the weighted sum, the overview table, and a prior-chain smoke test.

@testsnippet FixedDelays begin
    ## A PMF putting all its mass on lag `d`.
    at(d) = [zeros(d); 1.0]
end

@testitem "infection_state_weights counts states under fixed delays" setup = [
    FixedDelays,
] begin
    using BVDOutbreakSize: infection_state_weights, point_prevalence
    ## One infection a day for 30 days, onset five days after infection,
    ## detection two days after onset, death four and recovery seven.
    flat = ones(30)
    w(; asc, cfr) = infection_state_weights(
        at(5), at(2), at(4), at(7);
        ascertainment = asc, cfr, horizon = 30
    )
    all_detected = w(; asc = 1.0, cfr = 0.5)
    ## Infected 0 to 4 days ago: incubating.
    @test point_prevalence(flat, all_detected.incubating) ≈ 5
    ## Onset 0 or 1 days ago and not yet detected.
    @test point_prevalence(flat, all_detected.symptomatic) ≈ 2
    ## Undetected cases stay until death (4 days) or recovery (7 days).
    @test point_prevalence(flat, w(; asc = 0.0, cfr = 1.0).symptomatic) ≈ 4
    @test point_prevalence(flat, w(; asc = 0.0, cfr = 0.0).symptomatic) ≈ 7
    @test point_prevalence(flat, w(; asc = 0.5, cfr = 1.0).symptomatic) ≈
        0.5 * 2 + 0.5 * 4
end

@testitem "infection_state_weights rejects invalid arguments" setup = [
    FixedDelays,
] begin
    using BVDOutbreakSize: infection_state_weights
    f(; kw...) = infection_state_weights(
        at(1), at(1), at(1), at(1);
        ascertainment = 0.5, cfr = 0.5, horizon = 10, kw...
    )
    @test_throws ArgumentError f(; horizon = 0)
    @test_throws ArgumentError f(; ascertainment = 1.5)
    @test_throws ArgumentError f(; cfr = -0.1)
end

@testitem "point_prevalence weights each infection by its age" begin
    using BVDOutbreakSize: point_prevalence
    ## Three infections two days before the last day, weights by age 0, 1, 2.
    @test point_prevalence([3.0, 0.0, 0.0], [0.1, 0.2, 0.5]) ≈ 1.5
    ## Weights beyond the series and a series beyond the weights are ignored.
    @test point_prevalence([2.0], [0.5, 1.0, 1.0]) ≈ 1.0
    @test point_prevalence([5.0, 1.0], [1.0]) ≈ 1.0
end

@testitem "prevalence_overview rates, sorting and table" begin
    using BVDOutbreakSize: prevalence_overview, prevalence_table
    using DataFrames: nrow
    ## Two areas, two draws each; the second area has the higher rate.
    states = (:incubating, :symptomatic, :community)
    prev = NamedTuple{states}(
        ntuple(_ -> [[1.0, 3.0], [10.0, 30.0]], 3)
    )
    ov = prevalence_overview(
        prev, [100_000, 200_000]; labels = ["a", "b"], patch = [1, 2]
    )
    @test ov.label == ["b", "a"]
    @test ov.community_median == [20.0, 2.0]
    @test ov.community_rate_median == [10.0, 2.0]
    @test_throws DimensionMismatch prevalence_overview(
        prev, [1]; labels = ["a", "b"], patch = [1, 2]
    )
    tab = prevalence_table(ov; max_rows = 1, area = "Zone")
    @test nrow(tab) == 1
    @test tab[1, "Zone"] == "b"
    @test tab[1, "People in the community"] == "20 (11–29)"
end

@testitem "patch_prevalence on a four-patch prior chain" tags = [:slow] begin
    using BVDOutbreakSize
    using Turing: sample, Prior
    using Random: Xoshiro
    import FlexiChains

    obs = load_observations()
    np = length(PROVINCE_NAMES)
    m = bvd_joint(
        obs.n,
        obs.exported_cases, obs.total_deaths, obs.reported_cases,
        obs.exports_deaths, obs.confirmed_cases, obs.tests_analysed;
        reported_history = obs.reported_history,
        isolation_history = obs.isolation_history,
        bed_capacity_history = obs.bed_capacity_history,
        n_patches = np,
        breakpoint = obs.who_first_sitrep_days,
        tmrca_days = obs.tmrca_days
    )
    chn = sample(
        Xoshiro(20261008), m, Prior(), 5;
        chain_type = FlexiChains.VNChain, progress = false
    )
    prev = patch_prevalence(chn, np)
    @test length(prev.incubating) == np
    @test all(length(v) == 5 for v in prev.community)
    for p in 1:np
        @test all(>=(0), prev.symptomatic[p])
        @test prev.community[p] ≈ prev.incubating[p] .+ prev.symptomatic[p]
    end
end
