## Smoke tests for the laboratory-confirmed-deaths stream. Confirmed
## deaths (`Cumul décès parmi les confirmés`) are deaths that got
## laboratory-confirmed; the increment is a Binomial thinning of the
## suspected-death increment, with the thinning probability an INPUT from
## the confirmed-case BVD composition `q = μ_BVD/N_susp` times a tight
## enrichment scalar `m_death` (no death-specific denominator exists, so
## the probability is not a free curve). Fast `Prior()` + small `predict`
## runs, no NUTS where avoidable.

@testitem "bvd_count_composition: finite, bounded in (0, 1)" tags=[:slow] begin
    using Distributions: Gamma
    using BVDOutbreakSize: bvd_count_composition
    f_rep = Gamma(4.0, 3.0)
    q = bvd_count_composition(0.05, 0.3, 0.5, f_rep, [5.0, 20.0, 40.0];
        onset_fraction = 0.8)
    @test length(q) == 3
    @test all(isfinite, q)
    @test all(0 .<= q .<= 1)
    ## A zero-time edge has no backlog, so composition is zero there.
    q0 = bvd_count_composition(0.05, 0.3, 0.5, f_rep, [0.0])
    @test q0[1] == 0
    ## Extreme growth overflows the BVD backlog to Inf; the guarded ratio
    ## must still return a finite share (Inf/Inf would otherwise be NaN and
    ## propagate into a NaN Binomial gradient downstream).
    qbig = bvd_count_composition(5.0, 0.3, 0.5, f_rep, [200.0];
        onset_fraction = 1.0)
    @test all(isfinite, qbig)
    @test all(0 .<= qbig .<= 1)
end

@testitem "confirmed_deaths_only_model: prior draws finite, bounded" tags=[:slow] begin
    using Turing: sample, Prior
    import FlexiChains
    using BVDOutbreakSize: confirmed_deaths_only_model
    susp = 246
    chn = sample(confirmed_deaths_only_model(susp), Prior(), 200;
        chain_type = FlexiChains.VNChain, progress = false)
    cd = reduce(vcat, vec(Array(chn[:confirmed_deaths])))
    @test length(cd) == 200
    @test all(isfinite, cd)
    @test all(cd .>= 0)
    ## Bounded by the suspected-death denominator (single cumulative).
    @test all(cd .<= susp)
    ## The thinning probability at the cut-off is a genuine probability.
    p = vec(Array(chn[:p_death_conf_cutoff]))
    @test all(0 .<= p .<= 1)
    ## The enrichment scalar is centred on 1 (no enrichment).
    m = vec(Array(chn[:m_death]))
    @test all(m .> 0)
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
    ## The single informative point (17 confirmed deaths among 246
    ## suspected) must fit: finite, positive C(T) and a bounded thinning
    ## probability. Only the tight enrichment scalar is free.
    using BVDOutbreakSize: confirmed_deaths_only_model, nuts_sample
    chn = nuts_sample(confirmed_deaths_only_model(246, 17);
        samples = 50, chains = 2, seed = 1, progress = false)
    p = vec(Array(chn[:p_death_conf_cutoff]))
    @test all(0 .<= p .<= 1)
    m = vec(Array(chn[:m_death]))
    @test all(isfinite, m)
    @test all(m .> 0)
end

@testitem "bvd_joint: confirmed-deaths stream on vs off" tags=[:slow] begin
    ## Enabling the stream conditioned on the real 28-May confirmed-death
    ## point must keep finite generated quantities, and leave the model
    ## unchanged when the stream is off. Predictive draws are bounded by
    ## the suspected-death denominator.
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

    ## Stream on: a single cumulative confirmed-death total thinned from
    ## the cumulative suspected-death total at the cut-off.
    susp_cum = dh.values[end]
    cdeath_cum = obs.confirmed_death_history.values[end]
    model_on = bvd_joint(obs.exported_cases, death_incr, _inc(rh.values),
        obs.export_deaths_daily;
        reported_offsets = rh.offsets, death_offsets = dh.offsets,
        confirmed_cases = _inc(ch.values),
        confirmed_offsets = ch.offsets,
        tests_analysed = obs.cumulative_tests_analysed,
        confirmed_deaths = Union{Missing, Int}[cdeath_cum],
        confirmed_death_susp_increments = Union{Missing, Int}[susp_cum],
        confirmed_death_offsets = [0])
    chn_on = sample(model_on, Prior(), 100;
        chain_type = FlexiChains.VNChain, progress = false)
    C_on = vec(Array(chn_on[:cumulative_cases]))
    @test all(isfinite, C_on)
    p = vec(Array(chn_on[:p_death_conf_cutoff]))
    @test all(0 .<= p .<= 1)

    ## Predictive confirmed deaths are bounded by the suspected total.
    model_pp = bvd_joint(obs.exported_cases, death_incr, _inc(rh.values),
        obs.export_deaths_daily;
        reported_offsets = rh.offsets, death_offsets = dh.offsets,
        confirmed_cases = _inc(ch.values),
        confirmed_offsets = ch.offsets,
        tests_analysed = obs.cumulative_tests_analysed,
        confirmed_deaths = Union{Missing, Int}[missing],
        confirmed_death_susp_increments = Union{Missing, Int}[susp_cum],
        confirmed_death_offsets = [0])
    pp = predict(model_pp, chn_on)
    cd = reduce(vcat, vec(Array(pp[:confirmed_deaths])))
    @test all(0 .<= cd .<= susp_cum)
end

@testitem "bvd_joint: confirmed-deaths stream NUTS-fits (AD init)" tags=[:slow] begin
    ## Regression for the joint AD-init NaN: with the confirmed-deaths
    ## stream on, an extreme prior growth draw can overflow the
    ## count-composition `q`, giving a NaN `p_death_conf` whose Binomial
    ## reverse-mode gradient is NaN and kills NUTS initialisation. The
    ## guarded composition and the eps-clamped probability must let a
    ## gradient-based fit run and give finite cumulative cases.
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
        confirmed_death_susp_increments = Union{Missing, Int}[dh.values[end]],
        confirmed_death_offsets = [0],
        first_export_detection_delta = obs.first_export_detection_delta)
    chn = nuts_sample(model; samples = 5, chains = 1, seed = 1,
        progress = false)
    C = vec(Array(chn[:cumulative_cases]))
    @test all(isfinite, C)
    @test all(C .> 0)
    p = vec(Array(chn[:p_death_conf_cutoff]))
    @test all(0 .< p .< 1)
end
