## Smoke tests for the truth-anchored reported and confirmed cases
## likelihoods. C(T) is the latent eventually-confirmed pool; reported
## counts add an independent non-BVD background to the BVD-suspected
## convolution. Confirmed counts are the same BVD-suspected convolution
## evaluated one more step through f_lab, multiplied by the PCR
## sensitivity. The `@model` blocks live in the literate walkthrough;
## we recreate the minimal set here.

@testsnippet ConfirmedCasesFixtures begin
    using Distributions: Beta, Gamma, NegativeBinomial, truncated, Normal,
                         mean, std, pdf
    using Turing: Turing, @model, sample, Prior, to_submodel
    import FlexiChains
    using BVDOutbreakSize: delay_convolution, integrate

    function _safe_nbinomial(k, μ)
        p_raw = k / (k + max(μ, eps(typeof(μ))))
        p = isfinite(p_raw) ?
            clamp(p_raw, eps(typeof(k)), one(k) - eps(typeof(k))) :
            eps(typeof(k))
        return NegativeBinomial(k, p)
    end

    @model function _confirmed_test_growth()
        log_τ ~ Normal(log(14), 0.4)
        m ~ truncated(Normal(7.0, 2.5); lower = 0, upper = 13.0)
        τ := exp(log_τ)
        r := log(2) / τ
        T := m * τ
        C_T := 2.0^m
        return (; log_τ, τ, r, m, T, C_T)
    end

    @model function _confirmed_test_background(;
            lambda_prior = truncated(Normal(0.0, 10.0); lower = 0))
        λ_bg ~ lambda_prior
        return (; λ_bg)
    end

    @model function _confirmed_test_dispersion()
        inv_sqrt_k ~ truncated(Normal(0.5, 0.2); lower = 1e-3, upper = 5.0)
        k := 1.0 / (inv_sqrt_k^2 + eps(typeof(inv_sqrt_k)))
        return (; k, inv_sqrt_k)
    end

    @model function _confirmed_test_ascertainment(; p_prior = Beta(2.0, 6.0))
        p_drc ~ p_prior
        return (; p_drc)
    end

    @model function _confirmed_test_report_delay()
        α_rep ~ truncated(Normal(4.0, 1.5); lower = 0)
        θ_rep ~ truncated(Normal(3.0, 1.0); lower = 0)
        return (; α = α_rep, θ = θ_rep, dist = Gamma(α_rep, θ_rep))
    end

    @model function _confirmed_test_lab_delay()
        α_lab ~ truncated(Normal(2.0, 1.0); lower = 0)
        θ_lab ~ truncated(Normal(1.5, 0.75); lower = 0)
        return (; α = α_lab, θ = θ_lab, dist = Gamma(α_lab, θ_lab))
    end

    @model function _confirmed_test_sensitivity()
        s_test ~ Beta(15.0, 2.0)
        return (; s_test)
    end

    @model function _confirmed_test_reported(
            reported_cases::Union{Missing, Integer},
            growth_state, k::Real, p_drc::Real,
            λ_bg::Real, f_rep::Gamma)
        r = growth_state.r
        T = growth_state.T
        bvd_reported_at = let r = r, p_drc = p_drc, f_rep = f_rep
            τ -> p_drc * delay_convolution(one(p_drc), r, τ, f_rep)
        end
        μ_BVD_raw = bvd_reported_at(T)
        μ_bg_raw = λ_bg * T
        μ_BVD := isfinite(μ_BVD_raw) ?
                 max(μ_BVD_raw, eps(typeof(μ_BVD_raw))) :
                 eps(typeof(μ_BVD_raw))
        μ_bg := isfinite(μ_bg_raw) ?
                max(μ_bg_raw, eps(typeof(μ_bg_raw))) :
                eps(typeof(μ_bg_raw))
        expected_reports := μ_BVD + μ_bg
        positivity := μ_BVD / expected_reports
        reported_cases ~ _safe_nbinomial(k, expected_reports)
        return (; expected_reports, reported_cases, positivity,
            bvd_reported_at)
    end

    @model function _confirmed_test_confirmed(
            confirmed_cases::Union{Missing, Integer},
            bvd_reported_at, growth_state, k::Real,
            s_test::Real, f_lab::Gamma)
        T = growth_state.T
        integrand = let bvd_reported_at = bvd_reported_at, f_lab = f_lab, T = T
            u -> begin
                d = T - u
                d <= 0 ? zero(T) :
                pdf(f_lab, u) * bvd_reported_at(d)
            end
        end
        scale_lab = mean(f_lab) + 10 * std(f_lab)
        raw_confirmed = s_test * integrate(integrand, zero(T), T, scale_lab)
        expected_confirmed := isfinite(raw_confirmed) ?
                              max(raw_confirmed, eps(typeof(raw_confirmed))) :
                              eps(typeof(raw_confirmed))
        confirmed_cases ~ _safe_nbinomial(k, expected_confirmed)
        return (; expected_confirmed)
    end

    @model function _confirmed_test_only(
            reported_cases::Union{Missing, Integer},
            confirmed_cases::Union{Missing, Integer};
            growth = _confirmed_test_growth(),
            background = _confirmed_test_background(),
            report_delay = _confirmed_test_report_delay(),
            lab_delay = _confirmed_test_lab_delay(),
            test_sensitivity = _confirmed_test_sensitivity(),
            dispersion = _confirmed_test_dispersion(),
            ascertainment = _confirmed_test_ascertainment())
        growth_state ~ to_submodel(growth, false)
        background_state ~ to_submodel(background, false)
        report_state ~ to_submodel(report_delay, false)
        lab_state ~ to_submodel(lab_delay, false)
        sensitivity_state ~ to_submodel(test_sensitivity, false)
        dispersion_state ~ to_submodel(dispersion, false)
        asc_state ~ to_submodel(ascertainment, false)
        k = dispersion_state.k
        p_drc = asc_state.p_drc
        λ_bg = background_state.λ_bg
        s_test = sensitivity_state.s_test
        f_rep = report_state.dist
        f_lab = lab_state.dist
        reported_state ~ to_submodel(
            _confirmed_test_reported(
                reported_cases, growth_state, k, p_drc, λ_bg, f_rep), false)
        confirmed_state ~ to_submodel(
            _confirmed_test_confirmed(
                confirmed_cases, reported_state.bvd_reported_at,
                growth_state, k, s_test, f_lab), false)
        cumulative_cases := growth_state.C_T
    end
end

@testitem "truth-anchored streams: prior draws are finite and non-negative" tags=[:slow] setup=[ConfirmedCasesFixtures] begin
    m=_confirmed_test_only(missing, missing)
    chn=sample(m, Prior(), 200;
        chain_type = FlexiChains.VNChain, progress = false)
    rc=vec(Array(chn[:reported_cases]))
    cc=vec(Array(chn[:confirmed_cases]))
    @test length(cc) == 200
    @test all(isfinite, rc)
    @test all(isfinite, cc)
    @test all(rc .>= 0)
    @test all(cc .>= 0)

    π=vec(Array(chn[:positivity]))
    @test all(0 .<= π .<= 1)
end

@testitem "truth-anchored streams: fit with both observations" tags=[:slow] setup=[ConfirmedCasesFixtures] begin
    chn=sample(_confirmed_test_only(516, 33), Prior(), 200;
        chain_type = FlexiChains.VNChain, progress = false)
    C=vec(Array(chn[:cumulative_cases]))
    @test length(C) == 200
    @test all(isfinite, C)
    @test all(C .> 0)
end

@testitem "confirmed_only_model fits and predicts lab streams" tags=[:slow] begin
    ## Exercises the real production submodels through the
    ## confirmed-and-tested-only composer: a prior fit conditioned on the
    ## confirmed observation must give finite C(T), and `predict` on the
    ## same composer with missing observations must populate the
    ## per-vintage `confirmed_cases` and `samples_received` predictive
    ## draws. Confirmed positives are a Binomial on the analysed
    ## denominator; received samples are a NegBinomial on the suspect
    ## backlog (the τ_forward stream).
    using Turing: sample, Prior, predict
    import FlexiChains
    using BVDOutbreakSize: confirmed_only_model
    chn = sample(confirmed_only_model(101, 211), Prior(), 200;
        chain_type = FlexiChains.VNChain, progress = false)
    C = vec(Array(chn[:cumulative_cases]))
    @test length(C) == 200
    @test all(isfinite, C)
    @test all(C .> 0)
    ## Per-test positivity exposed as a derived quantity in (0, 1).
    pos = vec(Array(chn[:p_positive]))
    @test all(0 .< pos .< 1)

    pp = predict(confirmed_only_model(missing, missing), chn)
    ## Single confirmed vintage stored as a length-1 vector per draw. The
    ## received counts are a conditioned queue input, not a generated
    ## quantity, so only the confirmed positives are predicted.
    cc = reduce(vcat, vec(Array(pp[:confirmed_cases])))
    @test all(cc .>= 0)
end

@testitem "confirmed_only_model NUTS-fits with received observed" tags=[:slow] begin
    ## Regression test: supplying `tests_received` conditions the
    ## received-count NegBinomial as an observation. Without it the
    ## received count is a sampled discrete latent, which fails the NUTS
    ## model check ("sampled from a discrete distribution"). The
    ## gradient-based fit must run and give finite, positive C(T).
    using BVDOutbreakSize: confirmed_only_model, nuts_sample
    chn = nuts_sample(confirmed_only_model(121, 403, 662);
        samples = 5, chains = 1, seed = 1, progress = false)
    C = vec(Array(chn[:cumulative_cases]))
    @test all(isfinite, C)
    @test all(C .> 0)
end

@testitem "bvd_joint initialises with tested + confirmed observations" tags=[:slow] begin
    ## Regression test for the Binomial init failure: the full
    ## per-vintage `bvd_joint` path (deaths + reported + confirmed
    ## histories, tests_analysed) must produce finite log-densities
    ## under prior draws so NUTS can find valid initial parameters.
    ## Exercises the full per-test positivity Binomial branch with the
    ## real submodels and the real observation history.
    using Turing: sample, Prior
    import FlexiChains
    using BVDOutbreakSize: bvd_joint, load_observations
    obs = load_observations()
    rh = obs.reported_case_history
    ch = obs.confirmed_case_history
    dh = obs.death_history
    ## Models observe between-vintage increments, not cumulative totals.
    function _inc(values)
        out = similar(values, Int)
        prev = 0
        for i in eachindex(values)
            out[i] = values[i] - prev
            prev = values[i]
        end
        return out
    end
    chn = sample(
        bvd_joint(obs.exported_cases, _inc(dh.values), _inc(rh.values),
            obs.export_deaths_daily;
            reported_offsets = rh.offsets,
            death_offsets = dh.offsets,
            confirmed_cases = _inc(ch.values),
            confirmed_offsets = ch.offsets,
            tests_analysed = obs.cumulative_tests_analysed,
            tests_offset = 0,
            first_export_detection_delta = obs.first_export_detection_delta),
        Prior(), 200;
        chain_type = FlexiChains.VNChain, progress = false)
    C = vec(Array(chn[:cumulative_cases]))
    @test length(C) == 200
    @test all(isfinite, C)
end
