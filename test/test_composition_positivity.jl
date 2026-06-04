## Tests for the composition-linked confirmed positivity
## (`positivity_link = :composition`). The per-window tested BVD share is
## the suspect-pool composition φ_v = (p_drc·BVD)_v / ((p_drc·BVD)_v +
## λ_bg_v) upsampled by a decaying severity enrichment δ0 (see
## `severity_enrichment_model`), so the confirmed/positivity data identify
## the non-BVD background λ_bg rather than a free per-window random effect.

@testitem "severity_enrichment_model: δ0≥0, decay≥0, mean ≈ 1.5" begin
    using Turing: sample, Prior
    using Random: MersenneTwister
    import FlexiChains
    using Statistics: mean
    using Distributions: truncated, Normal
    using BVDOutbreakSize: severity_enrichment_model
    chn = sample(MersenneTwister(20260604), severity_enrichment_model(),
        Prior(), 20_000; chain_type = FlexiChains.VNChain, progress = false)
    δ = vec(Array(chn[:δ0]))
    dec = vec(Array(chn[:decay_scale]))
    @test all(δ .>= 0)
    @test all(dec .>= 0)
    ## Half-normal-ish centre near 1.5 (targets ~0.75 early tested-BVD for a
    ## pool composition φ ≈ 0.4).
    @test isapprox(mean(δ), 1.5; atol = 0.2)

    ## Overridable prior.
    chn2 = sample(MersenneTwister(20260604),
        severity_enrichment_model(;
            logodds_prior = truncated(Normal(3.0, 0.5); lower = 0)),
        Prior(), 8_000; chain_type = FlexiChains.VNChain, progress = false)
    @test mean(vec(Array(chn2[:δ0]))) > mean(δ)
end

@testitem "composition positivity: enrichment upsamples BVD (q ≥ φ)" begin
    ## With δ0 ≥ 0 the per-window tested BVD share must sit at or above the
    ## suspect-pool composition φ (severity triage over-tests BVD), and the
    ## per-window positivity stays a valid probability.
    using BVDOutbreakSize: confirmed_cases_model, reported_cases_model,
                           infection_model, onset_incidence_model
    using Turing: @model, to_submodel, returned
    using Random: MersenneTwister
    using Statistics: mean

    ## Local binning helper mirroring the model's `bin_increments`.
    function _bin(daily, days)
        out = zeros(eltype(daily), length(days))
        prev = 0
        for (i, d) in enumerate(days)
            hi = clamp(Int(d), 0, length(daily))
            lo = clamp(prev, 0, length(daily))
            out[i] = hi > lo ? sum(@view daily[(lo + 1):hi]) :
                     zero(eltype(daily))
            prev = hi
        end
        return out
    end

    @model function _latent(n)
        inf ~ to_submodel(infection_model(n), false)
        ons ~ to_submodel(onset_incidence_model(inf.infections), false)
        return ons.onsets
    end
    n = 40
    onsets = returned(_latent(n), rand(MersenneTwister(1), _latent(n)))

    @model function _rep(onsets)
        st ~ to_submodel(
            reported_cases_model(
                (; days = Int[], counts = Int[]), missing, onsets, 5.0, 0.3),
            false)
        return st
    end
    rep = returned(_rep(onsets), rand(MersenneTwister(3), _rep(onsets)))

    @model function _conf(onsets, rep)
        st ~ to_submodel(
            confirmed_cases_model(
                (; days = [20, 40], counts = [3, 8]), 8, onsets, 5.0, 0.3,
                rep.bg_daily, rep.τ_test, rep.bvd_reports_daily;
                lab_history = (; days = [20, 40], counts = [5, 9]),
                tests_received_history =
                (; days = [20, 40], counts = [6, 11]),
                positivity_link = :composition),
            false)
        return st
    end
    ## Reproduce the per-window composition φ from the same series the model
    ## consumes, so q ≥ φ can be checked directly per draw.
    qs = Float64[]
    φs = Float64[]
    for s in 1:200
        m = _conf(onsets, rep)
        st = returned(m, rand(MersenneTwister(100 + s), m))
        @test all(0 .<= st.p_pos .<= 1)
        @test 0 <= st.p_positive <= 1
        ## window day-indices are early then observed
        wd = vcat(st.windows.early_days, st.windows.obs_days)
        bvd_w = _bin(0.3 .* rep.bvd_reports_daily, wd)
        bg_w = _bin(rep.bg_daily, wd)
        for i in eachindex(wd)
            φ = bvd_w[i] / (bvd_w[i] + bg_w[i] + eps())
            push!(φs, φ)
            push!(qs, st.p_pos[i])
        end
    end
    ## Enrichment upsamples: q ≥ φ for (nearly) every window/draw.
    @test mean(qs .>= φs .- 1e-6) > 0.99
end

@testitem "bvd_joint :composition builds, samples, exposes δ0 / λ_bg" tags=[:slow] begin
    using Turing: sample, Prior
    import FlexiChains
    using BVDOutbreakSize: bvd_joint, load_observations, nuts_sample,
                           genetic_seeding_model
    obs = load_observations()
    BP = obs.n - obs.who_first_sitrep_days
    model = bvd_joint(
        obs.n, obs.exported_cases, obs.total_deaths,
        obs.reported_cases, obs.exports_deaths, obs.confirmed_cases,
        obs.tests_analysed;
        confirmed_deaths = obs.confirmed_deaths,
        deaths_history = obs.deaths_history,
        reported_history = obs.reported_history,
        confirmed_history = obs.confirmed_history,
        confirmed_deaths_history = obs.confirmed_deaths_history,
        lab_history = obs.lab_history,
        tests_received_history = obs.tests_received_history,
        breakpoint = BP, background_re = true,
        genetic = genetic_seeding_model, tmrca_days = obs.tmrca_days,
        confirmed_positivity_link = :composition)

    chn = sample(model, Prior(), 100;
        chain_type = FlexiChains.VNChain, progress = false)
    @test all(isfinite, vec(Array(chn[:C_T])))
    @test all(0 .<= vec(Array(chn[:test_positivity])) .<= 1)
    @test all(vec(Array(chn[:δ0])) .>= 0)
    @test all(vec(Array(chn[:lambda_bg])) .>= 0)

    ## A few NUTS steps: composition mode must initialise and stay finite.
    chn2 = nuts_sample(model; samples = 5, chains = 1, seed = 1,
        progress = false)
    @test all(isfinite, vec(Array(chn2[:C_T])))
end
