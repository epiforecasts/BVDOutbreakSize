## Smoke tests for the laboratory-confirmed-deaths stream. Confirmed
## deaths (`Cumul décès parmi les confirmés`) are true BVD deaths whose
## post-mortem swab was tested positive. The increment thins the MODELLED
## BVD-death trajectory by `coverage_death · s`, where `s` is the shared
## confirmed-case sensitivity and `coverage_death` is the free death-
## specimen submission rate (no `s · q` composition: the specimen is BVD).
## Fast `Prior()` + small `predict` runs, no NUTS where avoidable.

@testitem "death_coverage_model: prior bounded in (0, 1)" tags=[:slow] begin
    using Turing: sample, Prior
    import FlexiChains
    using BVDOutbreakSize: death_coverage_model
    chn = sample(death_coverage_model(), Prior(), 300;
        chain_type = FlexiChains.VNChain, progress = false)
    c = vec(Array(chn[:coverage_death]))
    @test length(c) == 300
    @test all(0 .< c .< 1)
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
