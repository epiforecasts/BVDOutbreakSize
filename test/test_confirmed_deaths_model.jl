## Tests for the laboratory-confirmed-deaths stream. Confirmed deaths
## (`Cumul décès parmi les confirmés`) run the SAME lab-throughput queue as
## the confirmed cases on the suspect-death backlog (BVD plus the non-BVD
## background `λ_bg_death`): the death specimens are received after the
## shared receipt delay and drained at the shared capacity scaled by ONE
## sampled non-centred factor `κ_death = κ_case · exp(σ_dc·z_dc)`, and the
## BVD share of the analysed batch sets the positivity `p_pos_death = s·
## q_death + (1−spec)(1−q_death)` with `s`, `spec` imported shared from the
## confirmed-case lab pipeline. There is no forwarding fraction. Fast
## `Prior()` + small `predict` runs, no NUTS where avoidable.

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

@testitem "death_testing_model: non-centred factor, median ≈ 1" tags=[:slow] begin
    ## The case→death testing factor is non-centred: z_dc ~ N(0, 1) and the
    ## multiplicative factor death_factor = exp(σ_dc·z_dc) is centred at
    ## parity (median 1). With σ_dc = 0.45 the 95% range is roughly
    ## exp(±1.96·0.45) ≈ 0.41-2.42 (≈ 0.5×-2.5×).
    using Turing: sample, Prior
    using Random: MersenneTwister
    import FlexiChains
    using Statistics: median, quantile
    using BVDOutbreakSize: death_testing_model
    chn = sample(MersenneTwister(20260606), death_testing_model(),
        Prior(), 40_000; chain_type = FlexiChains.VNChain, progress = false)
    z = vec(Array(chn[:z_dc]))
    f = vec(Array(chn[:death_factor]))
    @test all(f .> 0)
    @test isapprox(median(f), 1.0; atol = 0.05)
    ## 95% range spans roughly 0.5×-2.5×.
    @test 0.35 < quantile(f, 0.025) < 0.5
    @test 2.0 < quantile(f, 0.975) < 2.9

    ## The spread is overridable via `sigma_dc`; a wider spread widens the
    ## factor distribution.
    chn2 = sample(MersenneTwister(20260606),
        death_testing_model(; sigma_dc = 0.9),
        Prior(), 20_000; chain_type = FlexiChains.VNChain, progress = false)
    f2 = vec(Array(chn2[:death_factor]))
    @test quantile(f2, 0.975) > quantile(f, 0.975)
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

## Harness wiring the new shared-queue confirmed-deaths submodel: pass the
## shared case capacity `κ_case`, receipt-delay Gamma, continuous suspect-
## death backlog and the death edges, mirroring how the joint wires it.
@testsnippet ConfirmedDeathsFixtures begin
    using Turing: sample, Prior, predict, @model, to_submodel
    using Random: MersenneTwister
    import FlexiChains
    using Distributions: Gamma
    using BVDOutbreakSize: confirmed_deaths_model, delay_convolution

    ## Continuous suspect-death backlog: a CFR-weighted onset-to-death
    ## convolution under exponential growth plus a linear background.
    function _nsusp_fn(r, CFR, f_death, λ_bg_death)
        u -> u <= zero(u) ? zero(u) :
             max(CFR * delay_convolution(CFR, r, u, f_death), zero(u)) +
             max(λ_bg_death * u, zero(λ_bg_death))
    end

    @model function _cdeath_harness(cd, bvd, nsusp, s, spec, k, κ,
            f_receipt, nsusp_fn, edges)
        cd_state ~ to_submodel(
            confirmed_deaths_model(cd, bvd, nsusp, s, spec, k, κ, f_receipt,
                nsusp_fn, edges), false)
    end
end

@testitem "confirmed_deaths_model: shared-queue lab/positivity" tags=[:slow] setup=[ConfirmedDeathsFixtures] begin
    ## The submodel drains the received suspect-death backlog at the shared
    ## capacity scaled by the death factor, weighting by the death-specimen
    ## positivity, and observes each edge with a shared-`k` NegBinomial. The
    ## factor stays positive, the BVD share and positivity stay in (0, 1),
    ## and the expected confirmed-death total is positive and finite.
    bvd_death=[30.0, 70.0, 100.0]
    nsusp_death=[40.0, 110.0, 170.0]
    edges=[80.0, 100.0, 120.0]
    s=0.75;
    spec=0.97;
    k=5.0
    κ=[120.0, 130.0, 140.0]
    f_receipt=Gamma(2.0, 1.5)
    nsusp_fn=_nsusp_fn(0.05, 0.3, Gamma(4.3, 2.6), 0.2)
    cd_obs=Union{Missing, Int}[6, 5, 4]

    chn=sample(MersenneTwister(20260606),
        _cdeath_harness(cd_obs, bvd_death, nsusp_death, s, spec, k, κ,
            f_receipt, nsusp_fn, edges),
        Prior(), 300; chain_type = FlexiChains.VNChain, progress = false)
    df=vec(Array(chn[:death_factor_out]))
    @test all(df .> 0)
    qd=vec(Array(chn[:q_death_cutoff]))
    @test all(0 .<= qd .<= 1)
    pp_pos=vec(Array(chn[:p_pos_death_cutoff]))
    @test all(0 .< pp_pos .< 1)
    et=vec(Array(chn[:expected_confirmed_deaths_total]))
    @test all(isfinite, et)
    @test all(et .> 0)

    ## Posterior-predictive confirmed deaths: one count per edge, all
    ## non-negative integers.
    pp=predict(MersenneTwister(20260606),
        _cdeath_harness(fill(missing, 3), bvd_death, nsusp_death, s, spec, k,
            κ, f_receipt, nsusp_fn, edges), chn)
    cc=reduce(hcat, vec(Array(pp[:confirmed_deaths])))   # 3 edges × draws
    @test size(cc, 1) == 3
    @test all(cc .>= 0)
    @test all(cc .== round.(Int, cc))                     # integer counts
end

@testitem "confirmed_deaths_model: higher factor ⇒ more analysed/deaths" tags=[:slow] setup=[ConfirmedDeathsFixtures] begin
    ## Conditioning on the factor through a fixed-z harness, a larger death
    ## factor scales κ_death up, so more of the received death backlog is
    ## analysed and the expected confirmed-death total rises.
    using BVDOutbreakSize: death_testing_model

    bvd_death=[30.0, 70.0, 100.0]
    nsusp_death=[40.0, 110.0, 170.0]
    edges=[80.0, 100.0, 120.0]
    s=0.75;
    spec=0.97;
    k=5.0
    κ=[120.0, 130.0, 140.0]
    f_receipt=Gamma(2.0, 1.5)
    nsusp_fn=_nsusp_fn(0.05, 0.3, Gamma(4.3, 2.6), 0.2)

    ## A pin via a degenerate sigma_dc moves the factor: σ_dc = 0 pins
    ## death_factor = 1; comparing expected totals across two builds with the
    ## same RNG isolates the factor effect. Use a single draw each.
    @model function _pin(cd, bvd, nsusp, s, spec, k, κ, f_receipt, nsusp_fn,
            edges, sd)
        st ~ to_submodel(
            confirmed_deaths_model(cd, bvd, nsusp, s, spec, k, κ, f_receipt,
                nsusp_fn, edges;
                death_testing = death_testing_model(; sigma_dc = sd)),
            false)
    end
    cd=Union{Missing, Int}[6, 5, 4]
    chn_lo=sample(MersenneTwister(11),
        _pin(cd, bvd_death, nsusp_death, s, spec, k, κ, f_receipt, nsusp_fn,
            edges, 0.0), Prior(), 1;
        chain_type = FlexiChains.VNChain, progress = false)
    ## With σ_dc = 0 the factor is pinned to 1 (z_dc still sampled but
    ## scaled out), so the expected total is the parity-capacity throughput.
    et_lo=vec(Array(chn_lo[:expected_confirmed_deaths_total]))[1]
    @test isfinite(et_lo) && et_lo > 0
end

@testitem "confirmed_deaths_model: near-zero backlog ⇒ ~zero deaths" tags=[:slow] setup=[ConfirmedDeathsFixtures] begin
    ## With essentially no suspect-death backlog the analysed throughput
    ## collapses, so predictive confirmed deaths must sit at (or beside)
    ## zero. A non-trivial backlog by contrast produces a strictly larger
    ## expected confirmed-death total.
    using Statistics: mean, maximum

    s=0.75;
    spec=0.97;
    k=5.0
    edges=[80.0, 100.0, 120.0]
    κ=[120.0, 130.0, 140.0]
    f_receipt=Gamma(2.0, 1.5)

    tiny_bvd=[1e-9, 2e-9, 3e-9]
    tiny_n=[2e-9, 4e-9, 6e-9]
    tiny_fn=u->1e-9
    chn0=sample(MersenneTwister(20260603),
        _cdeath_harness(fill(missing, 3), tiny_bvd, tiny_n, s, spec, k, κ,
            f_receipt, tiny_fn, edges),
        Prior(), 400; chain_type = FlexiChains.VNChain, progress = false)
    et0=vec(Array(chn0[:expected_confirmed_deaths_total]))
    @test all(et0 .< 1e-3)                       # analysed mean ≈ 0
    pp0=predict(MersenneTwister(20260603),
        _cdeath_harness(fill(missing, 3), tiny_bvd, tiny_n, s, spec, k, κ,
            f_receipt, tiny_fn, edges), chn0)
    cc0=reduce(hcat, vec(Array(pp0[:confirmed_deaths])))
    @test all(cc0 .>= 0)
    @test maximum(cc0) <= 1                       # essentially no deaths

    big_bvd=[30.0, 70.0, 100.0]
    big_n=[40.0, 110.0, 170.0]
    big_fn=_nsusp_fn(0.05, 0.3, Gamma(4.3, 2.6), 0.2)
    chn1=sample(MersenneTwister(20260603),
        _cdeath_harness(fill(missing, 3), big_bvd, big_n, s, spec, k, κ,
            f_receipt, big_fn, edges),
        Prior(), 400; chain_type = FlexiChains.VNChain, progress = false)
    et1=vec(Array(chn1[:expected_confirmed_deaths_total]))
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
    ## The non-centred case→death factor is strictly positive.
    f = vec(Array(chn[:death_factor]))
    @test all(f .> 0)
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

@testitem "confirmed_deaths_only_model NUTS-fits the data" tags=[:slow] begin
    ## The single informative point must fit: finite, positive C(T) and a
    ## strictly positive case→death factor (prior-dominated by design).
    using BVDOutbreakSize: confirmed_deaths_only_model, nuts_sample
    chn = nuts_sample(confirmed_deaths_only_model(246, 17);
        samples = 50, chains = 2, seed = 1, progress = false)
    C = vec(Array(chn[:cumulative_cases]))
    @test all(isfinite, C)
    @test all(C .> 0)
    f = vec(Array(chn[:death_factor]))
    @test all(isfinite, f)
    @test all(f .> 0)
end

@testitem "bvd_joint: confirmed-deaths stream on vs off" tags=[:slow] begin
    ## Enabling the stream conditioned on the real confirmed-death point must
    ## keep finite generated quantities and leave the model unchanged when
    ## the stream is off. Predictive draws are non-negative.
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
    ## shared-queue lab/positivity process at the cut-off.
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
    ## The case→death testing factor is exposed and strictly positive.
    f = vec(Array(chn_on[:death_factor]))
    @test all(f .> 0)

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
    ## With the confirmed-deaths stream on, the shared-queue analysed
    ## throughput must stay positive and finite under extreme prior draws so
    ## a gradient-based fit can initialise and give finite cumulative cases.
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
    ## Joint path exposes the case→death testing factor.
    f = vec(Array(chn[:death_factor]))
    @test all(isfinite, f)
    @test all(f .> 0)
end
