## Tests for the direct confirmed-case fit. The laboratory-confirmed case
## counts are scored as a single NegBinomial likelihood on the between-vintage
## increments of the modelled daily confirmed series `confirmed_daily`, binned
## over the confirmed-history vintage day grid (sharing the dispersion `k`), so
## the confirmed level and shape inform the renewal trajectory. Each confirmed
## datum is used exactly once; the separate analysed-volume likelihoods score
## the ANALYSED counts, not the confirmed counts.

@testitem "confirmed_cases_model fits confirmed increments directly" begin
    using BVDOutbreakSize: confirmed_only_model
    using Turing: logjoint
    using Random: MersenneTwister

    m = confirmed_only_model(40, 8;
        confirmed_history = (; days = [18, 30, 40], counts = [3, 6, 8]),
        lab_history = (; days = [20, 40], counts = [5, 9]))
    lp = logjoint(m, rand(MersenneTwister(1), m))
    @test isfinite(lp)
end

@testitem "confirmed level informs the joint via the confirmed fit" begin
    using BVDOutbreakSize: confirmed_cases_model, reported_cases_model,
                           infection_model, onset_incidence_model
    using Turing: @model, to_submodel, returned
    using Turing.DynamicPPL: logjoint
    using Random: MersenneTwister

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
    rep_state = returned(_rep(onsets), rand(MersenneTwister(3), _rep(onsets)))

    ## Two confirmed histories on the same vintage day grid differing only in
    ## the confirmed LEVEL. The confirmed increments enter the likelihood, so
    ## the log-joint must move at a fixed draw.
    function build(conf)
        @model function _conf(onsets, rep)
            st ~ to_submodel(
                confirmed_cases_model(
                    conf, conf.counts[end], onsets, 5.0, 0.3,
                    rep.bg_daily, rep.τ_test, rep.bvd_reports_daily;
                    lab_history = (; days = [20, 40], counts = [5, 9])),
                false)
            return st
        end
        return _conf(onsets, rep_state)
    end
    m_a = build((; days = [18, 30, 40], counts = [3, 6, 8]))
    m_b = build((; days = [18, 30, 40], counts = [3, 40, 90]))
    draw = rand(MersenneTwister(4), m_a)
    @test logjoint(m_a, draw) != logjoint(m_b, draw)
end

@testitem "confirmed counts are scored exactly once" begin
    using BVDOutbreakSize: confirmed_only_model
    using Turing
    using Turing.DynamicPPL: predict
    using Random: MersenneTwister

    hist = (; confirmed_history = (; days = [18, 30, 40], counts = [3, 6, 8]),
        lab_history = (; days = [20, 40], counts = [5, 9]),
        lab_daily_history = (; days = [35], counts = [4]))

    ## In generator mode `predict` emits exactly the sampled-data keys. The
    ## confirmed counts are scored once, via a single confirmed-increment
    ## NegBinomial likelihood, so that key is present while the retired
    ## per-window confirmed likelihoods (the observed-window Binomial and the
    ## early/late confirmed-increment terms) are absent — no double count.
    m_fit = confirmed_only_model(40, 8; hist...)
    chn = sample(MersenneTwister(1), m_fit, Prior(), 8; progress = false)
    m_gen = confirmed_only_model(40, missing; hist...)
    pp = predict(MersenneTwister(2), m_gen, chn)
    ks = [string(k) for k in keys(pp)]
    @test any(occursin("confirmed_increments.increments", k) for k in ks)
    @test !any(occursin("confirmed_positives.positives", k) for k in ks)
    @test !any(occursin("early_increments.increments", k) for k in ks)
    @test !any(occursin("late_increments.increments", k) for k in ks)
    ## The analysed-volume likelihoods (a different stream) are still scored.
    @test any(occursin("analysed_increments.increments", k) for k in ks)
    @test any(occursin("analysed_daily_increments.increments", k) for k in ks)
end

@testitem "confirmed fit respects the testing-onset baseline" begin
    using BVDOutbreakSize: confirmed_cases_model, reported_cases_model,
                           infection_model, onset_incidence_model
    using Turing: @model, to_submodel, returned
    using Turing.DynamicPPL: logjoint
    using Random: MersenneTwister

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
    rep_state = returned(_rep(onsets), rand(MersenneTwister(3), _rep(onsets)))

    ## The first confirmed vintage is the testing-onset baseline and is not
    ## scored, so two histories sharing the same later increments (+3, +2) but
    ## different baselines have an equal log-joint at a fixed draw.
    function build(conf)
        @model function _conf(onsets, rep)
            st ~ to_submodel(
                confirmed_cases_model(
                    conf, conf.counts[end], onsets, 5.0, 0.3,
                    rep.bg_daily, rep.τ_test, rep.bvd_reports_daily;
                    lab_history = (; days = [20, 40], counts = [5, 9])),
                false)
            return st
        end
        return _conf(onsets, rep_state)
    end
    m_a = build((; days = [18, 30, 40], counts = [3, 6, 8]))
    m_b = build((; days = [18, 30, 40], counts = [5, 8, 10]))
    draw = rand(MersenneTwister(4), m_a)
    @test logjoint(m_a, draw) == logjoint(m_b, draw)
end
