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
                (; days = Int[], counts = Int[]), missing, onsets, 5.0, 0.3
            ),
            false
        )
        return st
    end
    rep_state = returned(_rep(onsets), rand(MersenneTwister(3), _rep(onsets)))

    @model function _conf(onsets, rep)
        st ~ to_submodel(
            confirmed_cases_model(
                (; days = [20, 40], counts = [3, 8]), 8, onsets, 5.0, 0.3,
                rep.bg_daily, rep.τ_test, rep.bvd_reports_daily;
                lab_history = (; days = [20, 40], counts = [5, 9])
            ),
            false
        )
        return st
    end
    m = _conf(onsets, rep_state)
    st = returned(m, rand(MersenneTwister(4), m))
    @test 0 <= st.p_positive <= 1
    @test st.expected_confirmed >= 0
    @test st.expected_analysed >= 0
    @test all(0 .<= st.p_pos .<= 1)
end

@testitem "confirmed_cases_model composition ties positivity to λ_bg" begin
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
                (; days = Int[], counts = Int[]), missing, onsets, 5.0, 0.3
            ),
            false
        )
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
                lab_history = (; days = [20, 40], counts = [5, 9])
            ),
            false
        )
        return st
    end
    m = _conf_comp(onsets, rep_state)
    st = returned(m, rand(MersenneTwister(4), m))
    @test 0 <= st.p_positive <= 1
    @test st.expected_confirmed >= 0
    @test all(isfinite, st.p_pos)
    @test all(0 .<= st.p_pos .<= 1)

    ## The severity-enrichment submodel constructs and stays positive.
    sv = returned(
        severity_enrichment_model(),
        rand(MersenneTwister(5), severity_enrichment_model())
    )
    @test sv.δ0 >= 0
    @test sv.decay_scale >= 0
end

@testitem "test_specificity_model returns a high-but-imperfect spec" begin
    using BVDOutbreakSize: test_specificity_model
    using Turing: returned
    using Random: MersenneTwister
    using Statistics: mean

    spec = returned(
        test_specificity_model(),
        rand(MersenneTwister(7), test_specificity_model())
    ).spec
    @test 0 < spec < 1
    ## Beta(60, 2) mean ≈ 0.97: a small false-positive rate, never zero.
    draws = [
        returned(
            test_specificity_model(),
            rand(MersenneTwister(i), test_specificity_model())
        ).spec
            for i in 1:500
    ]
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
                (; days = Int[], counts = Int[]), missing, onsets, 5.0, 0.3
            ),
            false
        )
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
                specificity = test_specificity_model(;
                    specificity_prior = Beta(700.0, 300.0)
                )
            ),
            false
        )
        return st
    end
    m = _conf_fp(onsets, rep_state)
    st = returned(m, rand(MersenneTwister(11), m))
    @test all(0 .< st.p_pos .< 1)
    ## False-positive floor ≈ 1 − 0.7 = 0.3 minus headroom: positivity is
    ## bounded away from zero by the non-BVD share, identifying λ_bg.
    @test all(st.p_pos .> 0.05)
end

@testitem "composition_positivity reads the BVD share of the pool window" begin
    using BVDOutbreakSize: composition_positivity

    ## With no enrichment and a perfect assay the positivity is the pool
    ## composition φ = BVD / pool, clamped to [lo, hi].
    bvd = [2.0, 0.0, 5.0, 1.0]
    pool = [8.0, 3.0, 5.0, 0.0]
    lo, hi = 1.0e-8, 1 - 1.0e-8
    p = composition_positivity(
        1:4, bvd, pool, zeros(4), 0.0, 1.0, 1.0, 1.0, lo, hi
    )
    @test p ≈ clamp.(bvd ./ (pool .+ lo), lo, hi)
    @test p[1] ≈ 0.25
    @test p[4] == hi
end

@testitem "composition_case_confirmation is true positives per BVD suspect" begin
    using BVDOutbreakSize: composition_positivity,
        composition_case_confirmation

    bvd = [2.0, 1.0, 5.0]
    pool = [8.0, 50.0, 6.0]
    c = [0.0, 10.0, 40.0]
    lo, hi = 1.0e-8, 1 - 1.0e-8
    s_test, spec = 0.9, 0.98
    ## With no enrichment a BVD suspect is tested as often as the pool
    ## average, so the per-case rate is the sensitivity alone.
    @test composition_case_confirmation(
        bvd, pool, c, 0.0, 1.0, s_test, lo, hi
    ) ≈ fill(s_test, 3)
    ## With enrichment it rises above the sensitivity, and times the BVD
    ## share it gives the positivity less its false positives.
    δ0, dscale = 1.5, 20.0
    r = composition_case_confirmation(bvd, pool, c, δ0, dscale, s_test, lo, hi)
    @test all(r .> s_test)
    p = composition_positivity(
        1:3, bvd, pool, c, δ0, dscale, s_test, spec, lo, hi
    )
    φ = bvd ./ (pool .+ lo)
    q = r .* φ ./ s_test
    @test p ≈ s_test .* q .+ (1 - spec) .* (1 .- q)
    ## A tiny BVD share does not blow the rate up through false positives.
    tiny = composition_case_confirmation(
        [1.0e-6], [100.0], [0.0], δ0, dscale, s_test, lo, hi
    )
    @test only(tiny) < s_test * exp(δ0) + 1.0e-6
end
