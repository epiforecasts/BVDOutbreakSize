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

@testitem "zone prevalence parameters shift the province CFR by the zone's" begin
    using BVDOutbreakSize: _zone_prevalence_parameters, _ZONE_DELAY_KEYS
    using StatsFuns: logit, logistic
    joint = [
        (;
            inc_mean = 9.0, recovery_mean = 20.0,
            ascertainment = [0.4, 0.6], cfr = [0.3, 0.5],
        ),
        (;
            inc_mean = 9.5, recovery_mean = 21.0,
            ascertainment = [0.5, 0.7], cfr = [0.4, 0.6],
        ),
    ]
    ## Three zone draws, so the third cycles back to the first joint draw.
    θ = [Float64.(eachindex(_ZONE_DELAY_KEYS)) .+ j for j in 1:3]
    ## Logit contrasts; the second draw's third zone is far in the tail.
    sev = [[-0.5, 0.5, 0.0], [0.0, 0.0, 30.0], [1.0, -1.0, 0.0]]
    asc = [[1.0, 2.0, 0.5], [2.5, 1.0, 1.0], [1.0, 1.0, 1.0]]
    w = [[1.0, 1.0, 5.0], [3.0, 1.0, 1.0], [1.0, 2.0, 1.0]]
    params = _zone_prevalence_parameters(joint, θ, asc, sev, w, [1, 1, 2])
    ## Within each province the zones, weighted by their infections, keep
    ## the province's ascertainment, in proportion to their contrasts; a
    ## province's only zone takes its rate.
    @test params[1].ascertainment ≈ [0.4 / 1.5, 0.8 / 1.5, 0.6]
    @test params[2].ascertainment ≈ [0.5 * 2.5 * 4 / 8.5, 0.5 * 4 / 8.5, 0.7]
    @test params[3].ascertainment ≈ [0.4, 0.4, 0.6]
    shift(c, s) = logistic(logit(c) + s)
    @test params[1].cfr ≈ shift.([0.3, 0.3, 0.5], sev[1])
    @test params[2].cfr[1:2] ≈ [0.4, 0.4]
    @test params[3].cfr ≈ shift.([0.3, 0.3, 0.5], sev[3])
    ## A zero contrast is the province's ratio, and an extreme one stays a
    ## probability.
    @test params[1].cfr[3] ≈ 0.5
    @test all(c -> 0 < c < 1, params[2].cfr)
    ## The zone's own delays replace the joint's; the stay is the joint's.
    inc = findfirst(==(:inc_mean), _ZONE_DELAY_KEYS)
    @test [p.inc_mean for p in params] == [θ[j][inc] for j in 1:3]
    @test [p.recovery_mean for p in params] == [20.0, 21.0, 20.0]
end

@testitem "patch ascertainment keeps the national rate" begin
    using BVDOutbreakSize: _patch_ascertainment
    rel = [0.5, 1.0, 2.0]
    w = [100.0, 50.0, 10.0]
    asc = _patch_ascertainment(0.3, rel, w)
    ## Weighted by the infections, the patches detect at the national rate,
    ## in proportion to their relative ascertainment.
    @test sum(asc .* w) / sum(w) ≈ 0.3
    @test asc ./ asc[2] ≈ rel
    ## The cap at one, and no infections to weight by.
    @test _patch_ascertainment(0.9, [1.0, 3.0], [10.0, 10.0]) ≈ [0.45, 1.0]
    @test _patch_ascertainment(0.3, rel, zeros(3)) ≈ 0.3 .* rel
end

@testitem "province case-fatality ratios stay in (0, 1) at extreme severity" begin
    using BVDOutbreakSize: province_composition_model, _logit_shift,
        infection_state_weights
    using Turing: fix
    using Distributions: truncated, Normal

    ## The joint's province ratio is the national one shifted by the death
    ## composition's severity contrast on the logit scale, and the
    ## prevalence weights reject a ratio outside [0, 1]. A national ratio
    ## near one and a contrast far into the tail must still give ratios
    ## inside it, each the national ratio times its multiplier.
    modelled = [10.0 12.0; 4.0 5.0; 1.0 1.0]
    obs = [8 9; 3 4; 1 1]
    m(cfr) = province_composition_model(
        obs, modelled;
        severity_sd_prior = truncated(Normal(0, 0.15); lower = 0),
        severity_baseline = cfr
    )
    pmf = [0.2, 0.3, 0.3, 0.2]
    for cfr in (0.05, 0.5, 0.9), z in ([4.0, -4.0], [-4.0, 4.0])
        r = fix(m(cfr), (; τ_sev = 5.0, z_sev = z))()
        cfr_patch = _logit_shift(cfr, r.severity_log_odds)
        @test all(c -> 0 < c < 1, cfr_patch)
        @test cfr_patch ≈ cfr .* r.province_severity
        @test maximum(abs, r.severity_log_odds) > 10
        for c in cfr_patch
            w = infection_state_weights(
                pmf, pmf, pmf, pmf; ascertainment = 0.5, cfr = c, horizon = 10
            )
            @test all(>=(0), w.symptomatic)
        end
    end
    ## No contrast is the common-ratio model.
    r0 = fix(m(0.4), (; τ_sev = 0.0, z_sev = [1.5, -0.7]))()
    @test _logit_shift(0.4, r0.severity_log_odds) ≈ fill(0.4, 3)
    @test r0.province_severity ≈ ones(3)
end

@testitem "ascertainment weights lag the infections by the confirmation delay" begin
    using BVDOutbreakSize: _cutoff_confirmations
    ## Confirmation two days after infection: an area whose infections
    ## stopped two days before the cut-off still weighs, one infected only
    ## on the cut-off day does not yet.
    f = [0.0, 0.0, 1.0]
    w = _cutoff_confirmations([[5.0, 0.0, 0.0], [0.0, 0.0, 5.0]], f)
    @test w == [5.0, 0.0]
end

@testitem "prevalence_overview rates, sorting and table" begin
    using BVDOutbreakSize: prevalence_overview, prevalence_table,
        plot_prevalence_ranking
    import CairoMakie
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
    @test tab[1, "Number in the community"] == "20 (11–29)"
    @test plot_prevalence_ranking(
        ov; patch_labels = ["A", "B"], max_rows = 1
    ) isa CairoMakie.Makie.Figure
end

@testitem "patch and zone prevalence on prior chains" tags = [:slow] setup = [
    ZoneSynthetic,
] begin
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

    ## The zone stage of the synthetic two-patch setup, with this chain as
    ## the joint its zone draws take their stay, ascertainment and CFR from.
    syn = zone_synthetic()
    inputs = zone_inputs(syn)
    zchn = sample(
        Xoshiro(20261013), bvd_zone(inputs.model_data), Prior(), 8;
        chain_type = FlexiChains.VNChain, progress = false
    )
    zprev = zone_prevalence(zchn, inputs, chn)
    nz = length(inputs.zone_keys)
    @test length(zprev.community) == nz
    @test all(length(v) == 8 for v in zprev.community)
    for a in 1:nz
        @test all(>=(0), zprev.symptomatic[a])
        @test zprev.community[a] ≈ zprev.incubating[a] .+ zprev.symptomatic[a]
    end
end
