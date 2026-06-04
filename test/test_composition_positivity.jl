## Tests for the composition-linked confirmed positivity
## (`positivity_link = :composition`). The tested BVD share is the
## suspect-pool composition φ = μ_BVD/(μ_BVD+μ_bg) upsampled by a decaying
## severity enrichment δ0 (see `severity_enrichment_model`), so the
## confirmed/positivity data identify the non-BVD background λ_bg rather
## than a free severe-first curve.

@testitem "severity_enrichment_model: δ0≥0, decay≥0, mean ≈ 1.5" tags=[:slow] begin
    using Turing: sample, Prior
    using Random: MersenneTwister
    import FlexiChains
    using Statistics: mean
    using Distributions: truncated, Normal
    using BVDOutbreakSize: severity_enrichment_model
    chn = sample(MersenneTwister(20260603), severity_enrichment_model(),
        Prior(), 20_000; chain_type = FlexiChains.VNChain, progress = false)
    δ = vec(Array(chn[:δ0]))
    dec = vec(Array(chn[:decay_scale]))
    @test all(δ .>= 0)
    @test all(dec .>= 0)
    ## Half-normal-ish centre near 1.5 (targets ~0.75 early tested-BVD for a
    ## pool composition φ ≈ 0.4).
    @test isapprox(mean(δ), 1.5; atol = 0.2)

    ## Overridable prior.
    chn2 = sample(MersenneTwister(20260603),
        severity_enrichment_model(;
            logodds_prior = truncated(Normal(3.0, 0.5); lower = 0)),
        Prior(), 8_000; chain_type = FlexiChains.VNChain, progress = false)
    @test mean(vec(Array(chn2[:δ0]))) > mean(δ)
end

@testitem "composition positivity: enrichment upsamples BVD (q ≥ φ)" tags=[:slow] begin
    ## With δ0 ≥ 0 the tested BVD share q must sit at or above the pool
    ## composition φ (severity triage over-tests BVD), and the per-test
    ## positivity stays a valid probability. Driven through the confirmed
    ## submodel on a single representative vintage with a known composition.
    using Turing: sample, Prior, @model, to_submodel
    using Random: MersenneTwister
    import FlexiChains
    using BVDOutbreakSize: confirmed_cases_model, exponential_growth_model,
                           severity_enrichment_model
    using Distributions: Gamma

    ## A single cumulative vintage at the cut-off; the composition φ is set by
    ## λ_bg against the BVD onset trajectory the growth/r imply. `f_rep` is a
    ## fixed onset-to-report Gamma (the report-delay prior centre).
    @model function _harness(link)
        g ~ to_submodel(exponential_growth_model(), false)
        c ~ to_submodel(
            confirmed_cases_model(
                Union{Missing, Int}[120], Union{Missing, Int}[400],
                Union{Missing, Int}[missing], 400, g, 5.0,
                [0.3], 0.6, 0.7, Gamma(2.5, 4.5),
                [g.T], g.T;
                positivity_link = link,
                severity_enrichment = severity_enrichment_model(),
                q_random_effect = nothing), false)
    end

    chn = sample(MersenneTwister(7), _harness(:composition), Prior(), 400;
        chain_type = FlexiChains.VNChain, progress = false)
    q = vec(Array(chn[:q_cutoff]))
    φ = vec(Array(chn[:q_baseline_count]))
    pp = vec(Array(chn[:p_positive]))
    @test all(0 .<= q .<= 1)
    @test all(0 .< pp .< 1)
    ## Enrichment upsamples: q ≥ φ for (nearly) every draw (allow a tiny
    ## numerical slack from the logit/logistic round-trip).
    @test mean(q .>= φ .- 1e-6) > 0.99
end

@testitem "bvd_joint :composition builds, samples, identifies λ_bg path" tags=[:slow] begin
    ## The joint with composition-linked positivity must sample finite
    ## generated quantities under the prior and expose δ0 / λ_bg, and fit a
    ## few NUTS steps (composition mode is AD-differentiable).
    using Turing: sample, Prior
    import FlexiChains
    using BVDOutbreakSize: bvd_joint, load_observations, nuts_sample,
                           severity_enrichment_model
    obs = load_observations()
    rh = obs.reported_case_history
    dh = obs.death_history
    ch = obs.confirmed_case_history
    sa = obs.tests_analysed_history
    sr = obs.tests_received_history
    function _inc(v)
        out = similar(v, Int)
        prev = 0
        for i in eachindex(v)
            out[i] = v[i] - prev
            prev = v[i]
        end
        return out
    end
    idx = [findfirst(==(off), ch.offsets) for off in sa.offsets]
    keep = [i == 1 || sa.values[i] > sa.values[i - 1]
            for i in eachindex(sa.values)]
    aoff = collect(sa.offsets)[keep]
    analysed = Union{Missing, Int}[_inc(sa.values[keep])...]
    conf = Union{Missing, Int}[_inc([ch.values[i] for i in idx][keep])...]
    ridx = [findfirst(==(off), sr.offsets) for off in aoff]
    received = Union{Missing, Int}[_inc([sr.values[i] for i in ridx])...]
    model = bvd_joint(obs.exported_cases, _inc(dh.values), _inc(rh.values),
        obs.export_deaths_daily;
        reported_offsets = rh.offsets, death_offsets = dh.offsets,
        confirmed_cases = conf, confirmed_offsets = aoff,
        samples_analysed = analysed, samples_received = received,
        tests_analysed = obs.cumulative_tests_analysed,
        confirmed_positivity_link = :composition,
        first_export_detection_delta = obs.first_export_detection_delta)

    chn = sample(model, Prior(), 100;
        chain_type = FlexiChains.VNChain, progress = false)
    C = vec(Array(chn[:cumulative_cases]))
    @test all(isfinite, C)
    @test all(0 .<= vec(Array(chn[:p_positive])) .<= 1)
    @test all(vec(Array(chn[:δ0])) .>= 0)
    @test all(vec(Array(chn[:λ_bg])) .>= 0)

    ## A few NUTS steps: composition mode must initialise and stay finite.
    chn2 = nuts_sample(model; samples = 5, chains = 1, seed = 1,
        progress = false)
    @test all(isfinite, vec(Array(chn2[:cumulative_cases])))
end
