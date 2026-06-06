## Tests for the laboratory-confirmed-deaths stream. Confirmed deaths
## (`Cumul décès parmi les confirmés`) are a genuine lab/positivity process
## on the post-mortem death specimens forwarded to the laboratory (issue
## #193): a fraction `τ_death` of the suspect-death backlog (BVD plus the
## non-BVD background `λ_bg_death`) is analysed, and its BVD share sets the
## positivity `p_pos_death = s·q_death + (1−spec)(1−q_death)` with `s`, `spec`
## imported shared from the confirmed-case lab pipeline. Fast `Prior()` +
## small `predict` runs, no NUTS where avoidable.

@testitem "death_background_model: half-normal mean ≈ 0.2" tags=[:slow] begin
    ## The non-BVD suspected-death background carries a half-normal
    ## `truncated(Normal(0, 0.25); lower = 0)` prior; its mean is the
    ## half-normal mean 0.25·√(2/π) ≈ 0.2. A fixed RNG keeps the check
    ## deterministic and non-flaky.
    using Turing: sample, Prior
    using Random: MersenneTwister
    import FlexiChains
    using Statistics: mean
    using Distributions: truncated, Normal
    using BVDOutbreakSize: death_background_model
    chn = sample(MersenneTwister(20260603), death_background_model(),
        Prior(), 20_000; chain_type = FlexiChains.VNChain, progress = false)
    λ = vec(Array(chn[:λ_bg_death]))
    @test length(λ) == 20_000
    @test all(λ .>= 0)
    @test isapprox(mean(λ), 0.2; atol = 0.02)

    ## The prior is overridable via the `lambda_prior` keyword.
    chn2 = sample(MersenneTwister(20260603),
        death_background_model(;
            lambda_prior = truncated(Normal(0.0, 1.0); lower = 0)),
        Prior(), 8_000; chain_type = FlexiChains.VNChain, progress = false)
    λ2 = vec(Array(chn2[:λ_bg_death]))
    @test mean(λ2) > mean(λ)                          # override took effect
end

@testitem "death_forward_model: Beta(2,8) prior, mean ≈ 0.2" tags=[:slow] begin
    ## The death-specimen forwarding fraction carries a weakly-informative
    ## Beta(2, 8) prior favouring low forwarding (mean 0.20, support (0, 1)).
    using Turing: sample, Prior
    using Random: MersenneTwister
    import FlexiChains
    using Statistics: mean
    using Distributions: Beta
    using BVDOutbreakSize: death_forward_model
    chn = sample(MersenneTwister(20260603), death_forward_model(),
        Prior(), 20_000; chain_type = FlexiChains.VNChain, progress = false)
    τ = vec(Array(chn[:τ_death]))
    @test length(τ) == 20_000
    @test all(0 .< τ .< 1)
    ## Beta(2, 8) mean = 2 / 10 = 0.20; loose tolerance avoids flakiness.
    @test isapprox(mean(τ), 0.20; atol = 0.01)

    ## The prior is overridable via the `fraction_prior` keyword.
    chn2 = sample(MersenneTwister(20260603),
        death_forward_model(; fraction_prior = Beta(8.0, 2.0)),
        Prior(), 8_000; chain_type = FlexiChains.VNChain, progress = false)
    τ2 = vec(Array(chn2[:τ_death]))
    @test mean(τ2) > mean(τ)                          # override took effect
end

@testitem "death_testing_model: LogNormal(log2,0.4), median ≈ 2" tags=[:slow] begin
    ## The death-testing enrichment factor carries a LogNormal(log 2, 0.4)
    ## prior centred on a 2× enrichment (deaths tested more than cases),
    ## strictly positive.
    using Turing: sample, Prior
    using Random: MersenneTwister
    import FlexiChains
    using Statistics: median
    using Distributions: LogNormal
    using BVDOutbreakSize: death_testing_model
    chn = sample(MersenneTwister(20260606), death_testing_model(),
        Prior(), 20_000; chain_type = FlexiChains.VNChain, progress = false)
    f = vec(Array(chn[:death_enrichment]))
    @test length(f) == 20_000
    @test all(f .> 0)
    ## LogNormal(log 2, ·) has median exp(log 2) = 2.
    @test isapprox(median(f), 2.0; atol = 0.05)

    ## The prior is overridable via the `factor_prior` keyword.
    chn2 = sample(MersenneTwister(20260606),
        death_testing_model(; factor_prior = LogNormal(log(4.0), 0.4)),
        Prior(), 8_000; chain_type = FlexiChains.VNChain, progress = false)
    f2 = vec(Array(chn2[:death_enrichment]))
    @test median(f2) > median(f)                       # override took effect
end

@testitem "confirmed_deaths_model: enrichment scales case rate" tags=[:slow] begin
    ## When the case forwarding rate `τ_forward` is supplied, the death rate
    ## is `clamp(τ_forward · death_enrichment, 0, 1)`: with the same RNG draws
    ## a fixed `τ_forward` and a 2×-centred enrichment must exceed the
    ## case rate, and the enrichment factor is exposed.
    using Turing: sample, Prior, @model, to_submodel
    using Random: MersenneTwister
    import FlexiChains
    using Statistics: median
    using BVDOutbreakSize: confirmed_deaths_model

    bvd_death = [30.0, 70.0, 100.0]
    nsusp_death = [40.0, 110.0, 170.0]
    s = 0.75; spec = 0.97; k = 5.0
    cd_obs = Union{Missing, Int}[6, 5, 4]
    τf = 0.3

    @model function _enrich_harness(cd, bvd, nsusp, s, spec, k, τf)
        cd_state ~ to_submodel(
            confirmed_deaths_model(cd, bvd, nsusp, s, spec, k;
                case_test_rate = τf), false)
    end

    chn = sample(MersenneTwister(20260606),
        _enrich_harness(cd_obs, bvd_death, nsusp_death, s, spec, k, τf),
        Prior(), 2_000; chain_type = FlexiChains.VNChain, progress = false)
    e = vec(Array(chn[:death_enrichment]))
    @test all(e .> 0)
    τd = vec(Array(chn[:τ_death_out]))
    @test all(0 .<= τd .<= 1)
    ## Death rate is the case rate enriched: with a 2×-centred factor the
    ## median death rate exceeds the case rate (before any clamp bites at
    ## these low rates).
    @test median(τd) > τf
    ## Where the product is below 1 it equals τ_forward · enrichment.
    unclamped = (τf .* e) .< 1
    @test all(isapprox.(τd[unclamped], τf .* e[unclamped]; atol = 1e-8))
end

@testitem "deaths_model: λ_bg_death lifts the suspected-death expectation" tags=[:slow] begin
    ## Adding a positive constant-rate background must raise the expected
    ## suspected-death total above the pure-BVD (λ_bg_death = 0) baseline by
    ## the linear background contribution, and the BVD-only trajectory must be
    ## unchanged. Drive the deaths submodel through a `to_submodel` harness
    ## with the growth / delay / CFR priors pinned via fixed RNG draws.
    using Turing: sample, Prior, @model, to_submodel
    using Random: MersenneTwister
    import FlexiChains
    using Statistics: mean
    using BVDOutbreakSize: deaths_model, exponential_growth_model

    @model function _deaths_harness(td, t_edges, λ)
        g ~ to_submodel(exponential_growth_model(), false)
        d ~ to_submodel(
            deaths_model(td, g, 5.0, t_edges; λ_bg_death = λ), false)
    end

    t_edges = [40.0, 80.0, 120.0]
    td = Union{Missing, Int}[20, 40, 60]
    chn0 = sample(MersenneTwister(7),
        _deaths_harness(td, t_edges, 0.0), Prior(), 1; progress = false,
        chain_type = FlexiChains.VNChain)
    chnλ = sample(MersenneTwister(7),
        _deaths_harness(td, t_edges, 0.5), Prior(), 1; progress = false,
        chain_type = FlexiChains.VNChain)
    e0 = vec(Array(chn0[:expected_deaths_T]))[1]
    eλ = vec(Array(chnλ[:expected_deaths_T]))[1]
    ## Same RNG ⇒ identical growth/CFR/delay draws, so the only difference is
    ## the background. The total expected suspected deaths at T rises by
    ## λ_bg_death · T = 0.5 · T (≈ tens of deaths), strictly positive.
    @test eλ > e0
end

@testitem "confirmed_deaths_model: lab/positivity on the suspect-death pool" tags=[:slow] begin
    ## The confirmed-death submodel forwards a fraction `τ_death` of the
    ## suspect-death backlog increment and weights it by the death-specimen
    ## positivity `s·q_death + (1−spec)(1−q_death)`, observing each edge with
    ## a shared-`k` NegBinomial. Drive it through a `to_submodel` harness with
    ## representative monotone BVD-death and total backlog trajectories and
    ## the shared sensitivity / specificity / dispersion passed in, matching
    ## how the joint wires the stream.
    using Turing: sample, Prior, predict, @model, to_submodel
    using Random: MersenneTwister
    import FlexiChains
    using BVDOutbreakSize: confirmed_deaths_model

    ## Cumulative modelled BVD deaths and total suspect deaths (BVD + bg) at
    ## three confirmed-death edges; the BVD share q_death = bvd / total falls
    ## as the background accumulates.
    bvd_death = [30.0, 70.0, 100.0]
    nsusp_death = [40.0, 110.0, 170.0]
    s = 0.75                                   # shared PCR sensitivity
    spec = 0.97                                # shared PCR specificity
    k = 5.0                                    # shared count dispersion
    cd_obs = Union{Missing, Int}[6, 5, 4]      # observed confirmed deaths

    @model function _cdeath_harness(cd, bvd, nsusp, s, spec, k)
        cd_state ~ to_submodel(
            confirmed_deaths_model(cd, bvd, nsusp, s, spec, k), false)
    end

    ## Conditioned fit: τ_death stays a genuine probability, the BVD share
    ## and positivity stay in (0, 1), and the tracked expected confirmed-death
    ## total is positive and finite.
    chn = sample(MersenneTwister(20260603),
        _cdeath_harness(cd_obs, bvd_death, nsusp_death, s, spec, k),
        Prior(), 300; chain_type = FlexiChains.VNChain, progress = false)
    τ = vec(Array(chn[:τ_death_out]))
    @test all(0 .< τ .< 1)
    qd = vec(Array(chn[:q_death_cutoff]))
    @test all(0 .<= qd .<= 1)
    pp_pos = vec(Array(chn[:p_pos_death_cutoff]))
    @test all(0 .< pp_pos .< 1)
    et = vec(Array(chn[:expected_confirmed_deaths_total]))
    @test all(isfinite, et)
    @test all(et .> 0)

    ## Posterior-predictive confirmed deaths: one count per edge, all
    ## non-negative integers.
    pp = predict(MersenneTwister(20260603),
        _cdeath_harness(fill(missing, 3), bvd_death, nsusp_death, s, spec, k),
        chn)
    cc = reduce(hcat, vec(Array(pp[:confirmed_deaths])))   # 3 edges × draws
    @test size(cc, 1) == 3
    @test all(cc .>= 0)
    @test all(cc .== round.(Int, cc))                     # integer counts
end

@testitem "confirmed_deaths_model: near-zero backlog ⇒ ~zero deaths" tags=[:slow] begin
    ## With essentially no suspect-death backlog the forwarded positive mean
    ## collapses, so predictive confirmed deaths must sit at (or beside) zero
    ## regardless of τ_death. A non-trivial backlog by contrast produces a
    ## strictly larger expected confirmed-death total, confirming the count
    ## scales with the forwarded positivity-weighted increment.
    using Turing: sample, Prior, predict, @model, to_submodel
    using Random: MersenneTwister
    import FlexiChains
    using Statistics: mean, maximum
    using BVDOutbreakSize: confirmed_deaths_model

    s = 0.75
    spec = 0.97
    k = 5.0

    @model function _cdeath_harness(cd, bvd, nsusp, s, spec, k)
        cd_state ~ to_submodel(
            confirmed_deaths_model(cd, bvd, nsusp, s, spec, k), false)
    end

    ## Near-zero suspect-death backlog: predictive deaths hug zero.
    tiny_bvd = [1e-9, 2e-9, 3e-9]
    tiny_n = [2e-9, 4e-9, 6e-9]
    chn0 = sample(MersenneTwister(20260603),
        _cdeath_harness(fill(missing, 3), tiny_bvd, tiny_n, s, spec, k),
        Prior(), 400; chain_type = FlexiChains.VNChain, progress = false)
    et0 = vec(Array(chn0[:expected_confirmed_deaths_total]))
    @test all(et0 .< 1e-3)                       # forwarded mean ≈ 0
    pp0 = predict(MersenneTwister(20260603),
        _cdeath_harness(fill(missing, 3), tiny_bvd, tiny_n, s, spec, k),
        chn0)
    cc0 = reduce(hcat, vec(Array(pp0[:confirmed_deaths])))
    @test all(cc0 .>= 0)
    @test maximum(cc0) <= 1                       # essentially no deaths

    ## A substantial backlog yields a strictly larger expected total.
    big_bvd = [30.0, 70.0, 100.0]
    big_n = [40.0, 110.0, 170.0]
    chn1 = sample(MersenneTwister(20260603),
        _cdeath_harness(fill(missing, 3), big_bvd, big_n, s, spec, k),
        Prior(), 400; chain_type = FlexiChains.VNChain, progress = false)
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
    ## Forwarding fraction is a genuine probability.
    τ = vec(Array(chn[:τ_death]))
    @test all(0 .< τ .< 1)
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
    ## positive C(T) and a bounded forwarding fraction. The free parameter is
    ## the death-specimen forwarding fraction, identified given the shared
    ## sensitivity / specificity through the suspect-death composition.
    using BVDOutbreakSize: confirmed_deaths_only_model, nuts_sample
    chn = nuts_sample(confirmed_deaths_only_model(246, 17);
        samples = 50, chains = 2, seed = 1, progress = false)
    τ = vec(Array(chn[:τ_death]))
    @test all(isfinite, τ)
    @test all(0 .< τ .< 1)
end

@testitem "confirmed_deaths_only: PP recentres on 17 (not over-predict)" tags=[:slow] begin
    ## The lab/positivity model must recover the observed 17 confirmed
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
    ## Background is sampled even when the confirmed-deaths stream is off
    ## (it feeds the suspected-death likelihood).
    λ_off = vec(Array(chn_off[:λ_bg_death]))
    @test all(λ_off .>= 0)

    ## Stream on: a single cumulative confirmed-death total via the
    ## lab/positivity process at the cut-off.
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
    ## Joint path derives the death rate as the case rate × enrichment, so the
    ## forwarding rate is exposed as `τ_death_out` (and the enrichment factor
    ## as `death_enrichment`), both genuine and bounded.
    τ = vec(Array(chn_on[:τ_death_out]))
    @test all(0 .<= τ .<= 1)
    e = vec(Array(chn_on[:death_enrichment]))
    @test all(e .> 0)

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
    ## With the confirmed-deaths stream on, the forwarded positivity-weighted
    ## mean must stay positive and finite under extreme prior draws so a
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
    ## Joint path: the derived death forwarding rate is exposed as
    ## `τ_death_out`, the case rate scaled by the death-testing enrichment.
    τ = vec(Array(chn[:τ_death_out]))
    @test all(0 .<= τ .<= 1)
    e = vec(Array(chn[:death_enrichment]))
    @test all(isfinite, e)
    @test all(e .> 0)
end
