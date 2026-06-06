@testitem "confirmed_cases_model exposes lab-pipeline positivity" begin
    using BVDOutbreakSize: confirmed_cases_model, reported_cases_model,
                           infection_model, onset_incidence_model
    using Turing: @model, to_submodel, returned
    using Turing.DynamicPPL: VarInfo
    using Random: MersenneTwister

    ## Build daily onsets and the shared report kernel / background /
    ## testing fraction by running the latent + reported submodels and
    ## reading their `returned` named tuples.
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

    @model function _conf(onsets, rep)
        st ~ to_submodel(
            confirmed_cases_model(
                (; days = [20, 40], counts = [3, 8]), 8, onsets, 5.0, 0.3,
                rep.bg_daily, rep.τ_test, rep.bvd_reports_daily;
                lab_history = (; days = [20, 40], counts = [5, 9]),
                tests_received_history =
                (; days = [20, 40], counts = [6, 11])),
            false)
        return st
    end
    m = _conf(onsets, rep_state)
    st = returned(m, rand(MersenneTwister(4), m))
    @test 0 <= st.p_positive <= 1
    @test st.expected_confirmed >= 0
    @test st.expected_received >= 0
    @test all(0 .<= st.p_pos .<= 1)
end

@testitem "confirmed_cases_model composition link ties positivity to λ_bg" begin
    using BVDOutbreakSize: confirmed_cases_model, reported_cases_model,
                           infection_model, onset_incidence_model,
                           severity_enrichment_model
    using Turing: @model, to_submodel, returned
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

    ## `:composition` link: the tested share is the suspect-pool composition
    ## upsampled by the severity enrichment, not a free per-window effect.
    @model function _conf_comp(onsets, rep)
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
    m = _conf_comp(onsets, rep_state)
    st = returned(m, rand(MersenneTwister(4), m))
    @test 0 <= st.p_positive <= 1
    @test st.expected_confirmed >= 0
    @test all(isfinite, st.p_pos)
    @test all(0 .<= st.p_pos .<= 1)

    ## The severity-enrichment submodel constructs and stays positive.
    sv = returned(severity_enrichment_model(),
        rand(MersenneTwister(5), severity_enrichment_model()))
    @test sv.δ0 >= 0
    @test sv.decay_scale >= 0
end

@testitem "test_specificity_model returns a high-but-imperfect spec" begin
    using BVDOutbreakSize: test_specificity_model
    using Turing: returned
    using Random: MersenneTwister
    using Statistics: mean

    spec = returned(test_specificity_model(),
        rand(MersenneTwister(7), test_specificity_model())).spec
    @test 0 < spec < 1
    ## Beta(60, 2) mean ≈ 0.97: a small false-positive rate, never zero.
    draws = [returned(test_specificity_model(),
                 rand(MersenneTwister(i), test_specificity_model())).spec
             for i in 1:500]
    @test 0.9 < mean(draws) < 1.0
    @test all(d -> 0 < d < 1, draws)
end

@testitem "composition positivity carries a false-positive floor" begin
    using BVDOutbreakSize: confirmed_cases_model, reported_cases_model,
                           infection_model, onset_incidence_model,
                           test_specificity_model
    using Turing: @model, to_submodel, returned
    using Distributions: Beta
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

    ## With a near-degenerate specificity prior (spec ≈ 0.7), every window's
    ## positivity must sit above the false-positive floor `1 − spec` because
    ## the transform `p = s·q + (1 − spec)(1 − q)` adds a non-BVD share.
    @model function _conf_fp(onsets, rep)
        st ~ to_submodel(
            confirmed_cases_model(
                (; days = [20, 40], counts = [3, 8]), 8, onsets, 5.0, 0.3,
                rep.bg_daily, rep.τ_test, rep.bvd_reports_daily;
                lab_history = (; days = [20, 40], counts = [5, 9]),
                tests_received_history =
                (; days = [20, 40], counts = [6, 11]),
                positivity_link = :composition,
                specificity = test_specificity_model(;
                    specificity_prior = Beta(700.0, 300.0))),
            false)
        return st
    end
    m = _conf_fp(onsets, rep_state)
    st = returned(m, rand(MersenneTwister(11), m))
    @test all(0 .< st.p_pos .< 1)
    ## False-positive floor ≈ 1 − 0.7 = 0.3 minus headroom: positivity is
    ## bounded away from zero by the non-BVD share, identifying λ_bg.
    @test all(st.p_pos .> 0.05)
end
