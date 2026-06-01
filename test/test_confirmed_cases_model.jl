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
                rep.λ_bg, rep.τ_test, rep.bvd_reports_daily;
                lab_history = (; days = [20, 40], counts = [5, 8]),
                tests_analysed = 8),
            false)
        return st
    end
    m = _conf(onsets, rep_state)
    st = returned(m, rand(MersenneTwister(4), m))
    @test 0 <= st.p_positive <= 1
    @test st.expected_confirmed >= 0
    @test st.expected_tested >= 0
end
