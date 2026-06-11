## Tests for the collapsed laboratory volume: a single suspected->analysed
## volume fitted to the specimens-analysed series (`lab_history`) and reused
## as the denominator in the early and unanchored late windows. The received-specimen
## stream is not modelled.

@testitem "confirmed_cases_model fits the analysed volume" begin
    using BVDOutbreakSize: confirmed_cases_model, reported_cases_model,
                           infection_model, onset_incidence_model
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

    ## No received series is supplied: the single modelled volume is fitted
    ## to the analysed series (`lab_history`) alone.
    @model function _conf(onsets, rep)
        st ~ to_submodel(
            confirmed_cases_model(
                (; days = [20, 40], counts = [3, 8]), 8, onsets, 5.0, 0.3,
                rep.bg_daily, rep.τ_test, rep.bvd_reports_daily;
                lab_history = (; days = [20, 40], counts = [5, 9]),
                tests_analysed = 9),
            false)
        return st
    end
    m = _conf(onsets, rep_state)
    st = returned(m, rand(MersenneTwister(4), m))
    @test 0 <= st.p_positive <= 1
    @test st.expected_confirmed >= 0
    @test st.expected_analysed >= 0
    @test all(0 .<= st.p_pos .<= 1)
end

@testitem "the analysed series drives the laboratory likelihood" begin
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

    ## Two analysed histories that differ on the same vintage day grid; the
    ## fitted volume and the window denominators both come from this series,
    ## so changing it moves the log-joint at a fixed draw.
    function build(lab)
        @model function _conf(onsets, rep)
            st ~ to_submodel(
                confirmed_cases_model(
                    (; days = [20, 40], counts = [3, 8]), 8, onsets, 5.0, 0.3,
                    rep.bg_daily, rep.τ_test, rep.bvd_reports_daily;
                    lab_history = lab, tests_analysed = lab.counts[end]),
                false)
            return st
        end
        return _conf(onsets, rep_state)
    end

    m_a = build((; days = [20, 40], counts = [5, 9]))
    m_b = build((; days = [20, 40], counts = [50, 90]))
    draw = rand(MersenneTwister(4), m_a)
    @test logjoint(m_a, draw) != logjoint(m_b, draw)
end

@testitem "the 24h analysed volume is a fitted per-day stream" begin
    using BVDOutbreakSize: confirmed_only_model
    using Turing
    using Turing.DynamicPPL: predict
    using Random: MersenneTwister

    hist = (; confirmed_history = (; days = [20, 30, 35], counts = [5, 9, 14]),
        lab_history = (; days = [10, 20], counts = [12, 28]),
        lab_daily_history = (; days = [35], counts = [30]))
    m_fit = confirmed_only_model(40, 14; hist...)
    chn = sample(MersenneTwister(1), m_fit, Prior(), 8; progress = false)
    ## In generator mode the per-day volume is resampled, so `predict` emits
    ## the `analysed_daily_increments.increments` key — i.e. it is a modelled,
    ## fitted stream rather than only a conditioned denominator.
    m_gen = confirmed_only_model(40, missing; hist...)
    pp = predict(MersenneTwister(2), m_gen, chn)
    @test any(occursin("analysed_daily_increments.increments", string(k))
    for k in keys(pp))
end

@testitem "bvd_joint fits the analysed volume without a received stream" tags=[:slow] begin
    using BVDOutbreakSize: bvd_joint, load_observations
    using Turing.DynamicPPL: logjoint
    using Random: MersenneTwister

    obs = load_observations()
    m = bvd_joint(
        obs.n, obs.exported_cases, obs.total_deaths,
        obs.reported_cases, obs.exports_deaths, obs.confirmed_cases,
        obs.tests_analysed;
        confirmed_deaths = obs.confirmed_deaths,
        deaths_history = obs.deaths_history,
        reported_history = obs.reported_history,
        confirmed_history = obs.confirmed_history,
        confirmed_deaths_history = obs.confirmed_deaths_history,
        lab_history = obs.lab_history,
        lab_daily_history = obs.lab_daily_history,
        breakpoint = obs.n - obs.who_first_sitrep_days,
        tmrca_days = obs.tmrca_days)
    draw = rand(MersenneTwister(1), m)
    @test isfinite(logjoint(m, draw))
end
