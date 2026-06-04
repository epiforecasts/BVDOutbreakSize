## Smoke tests for the laboratory-confirmed-deaths stream. Confirmed
## deaths (`Cumul décès parmi les confirmés`) are true BVD deaths whose
## post-mortem swab was tested positive. The increment thins the MODELLED
## BVD-death trajectory by `coverage_death · s`, where `s` is the shared
## confirmed-case sensitivity and `coverage_death` is the free death-
## specimen submission rate (no `s · q` composition: the specimen is BVD).
## Fast `Prior()` + small `predict` runs, no NUTS where avoidable.

@testitem "death_coverage_model: Beta(2,18) prior, mean ≈ 0.10" tags=[:slow] begin
    ## Post-mortem death-specimen coverage carries a weakly-informative
    ## Beta(2, 18) prior favouring low coverage (mean 0.10, support (0, 1)).
    ## A fixed RNG keeps the prior-mean check deterministic and non-flaky.
    using Turing: sample, Prior
    using Random: MersenneTwister
    import FlexiChains
    using Statistics: mean
    using Distributions: Beta
    using BVDOutbreakSize: death_coverage_model
    chn = sample(MersenneTwister(20260603), death_coverage_model(),
        Prior(), 20_000; chain_type = FlexiChains.VNChain, progress = false)
    c = vec(Array(chn[:coverage_death]))
    @test length(c) == 20_000
    @test all(0 .< c .< 1)
    ## Beta(2, 18) mean = 2 / 20 = 0.10; loose tolerance avoids flakiness.
    @test isapprox(mean(c), 0.10; atol = 0.01)

    ## The prior is overridable via the `coverage_prior` keyword.
    chn2 = sample(MersenneTwister(20260603),
        death_coverage_model(; coverage_prior = Beta(8.0, 2.0)),
        Prior(), 8_000; chain_type = FlexiChains.VNChain, progress = false)
    c2 = vec(Array(chn2[:coverage_death]))
    @test mean(c2) > mean(c)                          # override took effect
end

@testitem "bvd_death_trajectory: monotone non-negative" tags=[:slow] begin
    using Distributions: Gamma
    using BVDOutbreakSize: bvd_death_trajectory
    f_death = Gamma(4.0, 3.0)
    μ = bvd_death_trajectory(0.05, 0.4, f_death, [5.0, 20.0, 40.0];
        onset_fraction = 0.8)
    @test length(μ) == 3
    @test all(isfinite, μ)
    @test all(μ .>= 0)
    ## Cumulative expectation grows with elapsed time.
    @test issorted(μ)
    ## A zero-time edge has no backlog, so the trajectory is zero there.
    μ0 = bvd_death_trajectory(0.05, 0.4, f_death, [0.0])
    @test μ0[1] == 0
end

@testitem "confirmed_deaths_model: thins modelled BVD-death trajectory" tags=[:slow] begin
    ## The confirmed-death submodel thins the per-edge increment of the
    ## modelled cumulative BVD-death trajectory by `coverage_death · s`,
    ## observing each edge with a shared-`k` NegBinomial. Drive it through a
    ## `to_submodel` harness with a representative monotone trajectory and
    ## the shared sensitivity / dispersion passed in, matching how the joint
    ## wires the stream. A fixed RNG keeps the predictive draws deterministic.
    using Turing: sample, Prior, predict, @model, to_submodel
    using Random: MersenneTwister
    import FlexiChains
    using BVDOutbreakSize: confirmed_deaths_model

    ## Cumulative modelled BVD-death expectation at three confirmed-death
    ## edges; increments are 40, 60, 50 deaths.
    death_mu = [40.0, 100.0, 150.0]
    s = 0.75                                   # shared PCR sensitivity
    k = 5.0                                    # shared count dispersion
    cd_obs = Union{Missing, Int}[10, 12, 8]    # observed confirmed deaths

    @model function _cdeath_harness(confirmed_deaths, death_mu, s, k)
        cd_state ~ to_submodel(
            confirmed_deaths_model(confirmed_deaths, death_mu, s, k), false)
    end

    ## Conditioned fit: coverage stays a genuine probability and the tracked
    ## expected confirmed-death total is positive and finite.
    chn = sample(MersenneTwister(20260603),
        _cdeath_harness(cd_obs, death_mu, s, k), Prior(), 300;
        chain_type = FlexiChains.VNChain, progress = false)
    cov = vec(Array(chn[:coverage_death_out]))
    @test all(0 .< cov .< 1)
    et = vec(Array(chn[:expected_confirmed_deaths_total]))
    @test all(isfinite, et)
    @test all(et .> 0)

    ## Posterior-predictive confirmed deaths: one count per edge, all
    ## non-negative integers.
    pp = predict(MersenneTwister(20260603),
        _cdeath_harness(fill(missing, 3), death_mu, s, k), chn)
    cc = reduce(hcat, vec(Array(pp[:confirmed_deaths])))   # 3 edges × draws
    @test size(cc, 1) == 3
    @test all(cc .>= 0)
    @test all(cc .== round.(Int, cc))                     # integer counts
end

@testitem "confirmed_deaths_model: near-zero trajectory ⇒ ~zero deaths" tags=[:slow] begin
    ## With essentially no modelled BVD deaths the thinned NegBinomial mean
    ## collapses, so predictive confirmed deaths must sit at (or beside)
    ## zero regardless of coverage·s. A non-trivial trajectory by contrast
    ## produces a strictly larger expected confirmed-death total, confirming
    ## the count scales with the thinned increment.
    using Turing: sample, Prior, predict, @model, to_submodel
    using Random: MersenneTwister
    import FlexiChains
    using Statistics: mean, maximum
    using BVDOutbreakSize: confirmed_deaths_model

    s = 0.75
    k = 5.0

    @model function _cdeath_harness(confirmed_deaths, death_mu, s, k)
        cd_state ~ to_submodel(
            confirmed_deaths_model(confirmed_deaths, death_mu, s, k), false)
    end

    ## Near-zero modelled BVD-death trajectory: predictive deaths hug zero.
    tiny = [1e-9, 2e-9, 3e-9]
    chn0 = sample(MersenneTwister(20260603),
        _cdeath_harness(fill(missing, 3), tiny, s, k), Prior(), 400;
        chain_type = FlexiChains.VNChain, progress = false)
    et0 = vec(Array(chn0[:expected_confirmed_deaths_total]))
    @test all(et0 .< 1e-3)                       # thinned mean ≈ 0
    pp0 = predict(MersenneTwister(20260603),
        _cdeath_harness(fill(missing, 3), tiny, s, k), chn0)
    cc0 = reduce(hcat, vec(Array(pp0[:confirmed_deaths])))
    @test all(cc0 .>= 0)
    @test maximum(cc0) <= 1                       # essentially no deaths

    ## A substantial trajectory yields a strictly larger expected total.
    big = [40.0, 100.0, 150.0]
    chn1 = sample(MersenneTwister(20260603),
        _cdeath_harness(fill(missing, 3), big, s, k), Prior(), 400;
        chain_type = FlexiChains.VNChain, progress = false)
    et1 = vec(Array(chn1[:expected_confirmed_deaths_total]))
    @test mean(et1) > mean(et0)
end

@testitem "confirmed_deaths_only_model: prior draws finite, bounded" tags=[:slow] begin
    using Turing: sample, Prior
    import FlexiChains
    using BVDOutbreakSize: confirmed_deaths_only_model
    chn = sample(confirmed_deaths_only_model(246), Prior(), 200;
        chain_type = FlexiChains.VNChain, progress = false)
    cd = reduce(vcat, vec(Array(chn[:confirmed_deaths])))
    @test length(cd) == 200
    @test all(isfinite, cd)
    @test all(cd .>= 0)
    ## Coverage is a genuine probability.
    c = vec(Array(chn[:coverage_death]))
    @test all(0 .< c .< 1)
end

@testitem "confirmed_deaths_only_model: conditioning gives finite C(T)" tags=[:slow] begin
    using Turing: sample, Prior
    import FlexiChains
    using BVDOutbreakSize: confirmed_deaths_only_model
    chn = sample(confirmed_deaths_only_model(246, 17), Prior(), 200;
        chain_type = FlexiChains.VNChain, progress = false)
    C = vec(Array(chn[:cumulative_cases]))
    @test length(C) == 200
    @test all(isfinite, C)
    @test all(C .> 0)
end

@testitem "confirmed_deaths_only_model NUTS-fits the 28-May data" tags=[:slow] begin
    ## The single informative point (17 confirmed deaths) must fit: finite,
    ## positive C(T) and a bounded coverage. The free parameter is the
    ## death-specimen coverage, identified given the shared sensitivity.
    using BVDOutbreakSize: confirmed_deaths_only_model, nuts_sample
    chn = nuts_sample(confirmed_deaths_only_model(246, 17);
        samples = 50, chains = 2, seed = 1, progress = false)
    c = vec(Array(chn[:coverage_death]))
    @test all(isfinite, c)
    @test all(0 .< c .< 1)
end

@testitem "confirmed_deaths_only: PP recentres on 17 (not 38)" tags=[:slow] begin
    ## The coverage·s thinning must recover the observed 17 confirmed
    ## deaths. The old composition link (s·q) over-predicted (~38); fitting
    ## the 28-May point and drawing the posterior predictive must centre the
    ## replicate confirmed-death total near 17, well below 38.
    using Turing: predict
    using Statistics: median
    using BVDOutbreakSize: confirmed_deaths_only_model, nuts_sample
    chn = nuts_sample(confirmed_deaths_only_model(246, 17);
        samples = 200, chains = 2, seed = 1, progress = false)
    pp = predict(confirmed_deaths_only_model(246, missing), chn)
    cd = reduce(vcat, vec(Array(pp[:confirmed_deaths])))
    m = median(cd)
    @test all(cd .>= 0)
    ## Centred on the observed point, far from the composition-link 38.
    @test 8 <= m <= 30
    @test m < 38
end

@testitem "bvd_joint: confirmed-deaths stream on vs off" tags=[:slow] begin
    ## Enabling the stream conditioned on the real 28-May confirmed-death
    ## point must keep finite generated quantities, and leave the model
    ## unchanged when the stream is off. Predictive draws are non-negative.
    using Turing: sample, Prior, predict
    import FlexiChains
    using BVDOutbreakSize: bvd_joint, load_observations
    obs = load_observations()
    rh = obs.reported_case_history
    dh = obs.death_history
    ch = obs.confirmed_case_history
    function _inc(values)
        out = similar(values, Int)
        prev = 0
        for i in eachindex(values)
            out[i] = values[i] - prev
            prev = values[i]
        end
        return out
    end
    death_incr = _inc(dh.values)

    ## Stream off (default): finite cumulative cases.
    chn_off = sample(
        bvd_joint(obs.exported_cases, death_incr, _inc(rh.values),
            obs.export_deaths_daily;
            reported_offsets = rh.offsets, death_offsets = dh.offsets,
            confirmed_cases = _inc(ch.values),
            confirmed_offsets = ch.offsets,
            tests_analysed = obs.cumulative_tests_analysed),
        Prior(), 100;
        chain_type = FlexiChains.VNChain, progress = false)
    C_off = vec(Array(chn_off[:cumulative_cases]))
    @test all(isfinite, C_off)

    ## Stream on: a single cumulative confirmed-death total thinned from the
    ## modelled BVD-death trajectory at the cut-off.
    cdeath_cum = obs.confirmed_death_history.values[end]
    model_on = bvd_joint(obs.exported_cases, death_incr, _inc(rh.values),
        obs.export_deaths_daily;
        reported_offsets = rh.offsets, death_offsets = dh.offsets,
        confirmed_cases = _inc(ch.values),
        confirmed_offsets = ch.offsets,
        tests_analysed = obs.cumulative_tests_analysed,
        confirmed_deaths = Union{Missing, Int}[cdeath_cum],
        confirmed_death_offsets = [0])
    chn_on = sample(model_on, Prior(), 100;
        chain_type = FlexiChains.VNChain, progress = false)
    C_on = vec(Array(chn_on[:cumulative_cases]))
    @test all(isfinite, C_on)
    c = vec(Array(chn_on[:coverage_death]))
    @test all(0 .< c .< 1)

    ## Predictive confirmed deaths are finite and non-negative.
    model_pp = bvd_joint(obs.exported_cases, death_incr, _inc(rh.values),
        obs.export_deaths_daily;
        reported_offsets = rh.offsets, death_offsets = dh.offsets,
        confirmed_cases = _inc(ch.values),
        confirmed_offsets = ch.offsets,
        tests_analysed = obs.cumulative_tests_analysed,
        confirmed_deaths = Union{Missing, Int}[missing],
        confirmed_death_offsets = [0])
    pp = predict(model_pp, chn_on)
    cd = reduce(vcat, vec(Array(pp[:confirmed_deaths])))
    @test all(isfinite, cd)
    @test all(cd .>= 0)
end

@testitem "bvd_joint: confirmed-deaths stream NUTS-fits (AD init)" tags=[:slow] begin
    ## With the confirmed-deaths stream on, the thinned NegBinomial mean
    ## must stay positive and finite under extreme prior draws so a
    ## gradient-based fit can initialise and give finite cumulative cases.
    using BVDOutbreakSize: bvd_joint, load_observations, nuts_sample
    obs = load_observations()
    rh = obs.reported_case_history
    dh = obs.death_history
    ch = obs.confirmed_case_history
    sa = obs.tests_analysed_history
    sr = obs.tests_received_history
    function _inc(values)
        out = similar(values, Int)
        prev = 0
        for i in eachindex(values)
            out[i] = values[i] - prev
            prev = values[i]
        end
        return out
    end
    ## Per-vintage confirmed / analysed / received on the lab vintages, so
    ## the received NegBinomial is observed (not a discrete latent NUTS
    ## rejects) and the confirmed-deaths stream is conditioned alongside.
    idx = [findfirst(==(off), ch.offsets) for off in sa.offsets]
    keep = [i == 1 || sa.values[i] > sa.values[i - 1]
            for i in eachindex(sa.values)]
    aoff = collect(sa.offsets)[keep]
    analysed = Union{Missing, Int}[_inc(sa.values[keep])...]
    conf = Union{Missing, Int}[_inc([ch.values[i] for i in idx][keep])...]
    ridx = [findfirst(==(off), sr.offsets) for off in aoff]
    received = Union{Missing, Int}[_inc([sr.values[i] for i in ridx])...]
    model = bvd_joint(obs.exported_cases, _inc(dh.values),
        _inc(rh.values), obs.export_deaths_daily;
        reported_offsets = rh.offsets, death_offsets = dh.offsets,
        confirmed_cases = conf, confirmed_offsets = aoff,
        samples_analysed = analysed, samples_received = received,
        tests_analysed = obs.cumulative_tests_analysed,
        confirmed_deaths = Union{Missing, Int}[
        obs.confirmed_death_history.values[end]],
        confirmed_death_offsets = [0],
        first_export_detection_delta = obs.first_export_detection_delta)
    chn = nuts_sample(model; samples = 5, chains = 1, seed = 1,
        progress = false)
    C = vec(Array(chn[:cumulative_cases]))
    @test all(isfinite, C)
    @test all(C .> 0)
    c = vec(Array(chn[:coverage_death]))
    @test all(0 .< c .< 1)
end
