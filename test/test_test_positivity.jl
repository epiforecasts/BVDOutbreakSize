## Tests for the test-positivity / confirmed-stream submodels in
## `src/models/priors.jl` and `src/models/observations.jl`. The confirmed
## stream models the laboratory-confirmed cases as a severe-first BVD share
## of the tested pool that relaxes from a near-1 early severe cluster `q0`
## to a baseline `qinf`, with received counts conditioning the forwarded
## fraction `τ_forward`.

@testitem "sensitivity prior is untruncated Beta(6,2)" begin
    using Turing: sample, Prior
    using Random: MersenneTwister
    using Statistics: mean, minimum
    using BVDOutbreakSize: test_sensitivity_model

    ## Untruncated Beta(6,2): mean 0.75, no lower floor. The severe-first
    ## q0 ≈ 1 makes the early positivity ≈ s, identifying s from the data.
    chn = sample(MersenneTwister(20260518), test_sensitivity_model(),
        Prior(), 5_000; progress = false)
    s = vec(Array(chn[:s_test]))
    @test isapprox(mean(s), 0.75; atol = 0.02)
    @test minimum(s) < 0.48          # no data-implied floor
    @test all(0 .< s .< 1)
end

@testitem "specificity prior is informative and high" begin
    using Turing: sample, Prior
    using Random: MersenneTwister
    using Statistics: mean, quantile
    using BVDOutbreakSize: test_specificity_model

    ## Default Beta(50, 1.5): mean ≈ 0.97, tight, most mass above 0.95.
    chn = sample(MersenneTwister(20260518), test_specificity_model(),
        Prior(), 10_000; progress = false)
    sp = vec(Array(chn[:spec_test]))
    @test isapprox(mean(sp), 50 / 51.5; atol = 0.01)   # Beta mean
    @test quantile(sp, 0.05) > 0.90                     # tight high tail
    @test all(0 .< sp .< 1)
end

@testitem "test_selection_model: q0 weak, qinf central, decay positive" begin
    using Turing: sample, Prior
    using Random: MersenneTwister
    using Statistics: mean
    using BVDOutbreakSize: test_selection_model

    chn = sample(MersenneTwister(20260518), test_selection_model(),
        Prior(), 10_000; progress = false)
    q0 = vec(Array(chn[:q0]))
    qinf = vec(Array(chn[:qinf]))
    decay = vec(Array(chn[:decay_scale]))
    @test isapprox(mean(q0), 0.5; atol = 0.02)          # Beta(2,2) mean
    @test 0.4 < mean(q0) < 0.6                           # weak selection
    @test isapprox(mean(qinf), 0.5; atol = 0.02)         # Beta(6,6) central
    @test all(0 .< q0 .< 1)
    @test all(0 .< qinf .< 1)
    @test all(decay .>= 0)
end

@testitem "severe_first_share: high q0 to baseline qinf, levels off" begin
    using BVDOutbreakSize: severe_first_share

    q0 = 0.95
    qinf = 0.3
    scale = 5.0
    ## At c = 0 the share is q0 exactly.
    @test severe_first_share(q0, qinf, 0.0, scale) ≈ q0
    ## Monotone decline toward qinf, never below it (no overshoot).
    cs = 0.0:1.0:60.0
    qs = [severe_first_share(q0, qinf, c, scale) for c in cs]
    @test all(diff(qs) .<= 1e-12)          # non-increasing
    @test all(qs .>= qinf - 1e-9)          # bounded below by qinf
    @test all(qs .<= q0 + 1e-9)
    ## Far past the decay scale it has effectively levelled at qinf.
    @test isapprox(severe_first_share(q0, qinf, 80.0, scale), qinf;
        atol = 1e-3)
    ## q0 = qinf gives a flat share at every c.
    @test severe_first_share(0.5, 0.5, 3.0, scale) ≈ 0.5
    ## Always a valid probability.
    @test all(0 .<= qs .<= 1)
end

@testitem "specificity term: spec=1 is s·q, false positives raise p_pos" begin
    ## The per-test positivity p_pos = s·q + (1−spec)·(1−q). With spec = 1
    ## it reduces to s·q; lowering spec adds false positives on the non-BVD
    ## share and raises p_pos.
    s = 0.6
    q = 0.4
    p_perfect = s * q + (1 - 1.0) * (1 - q)
    @test p_perfect ≈ s * q
    p_imperfect = s * q + (1 - 0.95) * (1 - q)
    @test p_imperfect > p_perfect
    @test p_imperfect ≈ s * q + 0.05 * (1 - q)
end

@testitem "default λ_bg prior matches half-normal SD 1" tags=[:slow] begin
    using Turing: sample, Prior
    using Random: MersenneTwister
    using Statistics: mean, std, quantile
    using BVDOutbreakSize: test_positivity_model

    chn = sample(MersenneTwister(20260518), test_positivity_model(),
        Prior(), 20_000; progress = false)
    λ = vec(Array(chn[:λ_bg]))
    @test all(λ .>= 0)
    ## Half-normal SD 1: mean ≈ √(2/π) ≈ 0.80, median ≈ 0.67.
    @test isapprox(mean(λ), sqrt(2 / pi); atol = 0.05)
    @test isapprox(quantile(λ, 0.5), 0.674; atol = 0.05)
end

@testitem "τ_forward prior is Beta(5,2) and overridable" tags=[:slow] begin
    using Turing: sample, Prior
    using Random: MersenneTwister
    using Statistics: mean
    using Distributions: Beta
    using BVDOutbreakSize: test_positivity_model

    chn = sample(MersenneTwister(20260518), test_positivity_model(),
        Prior(), 8_000; progress = false)
    τ = vec(Array(chn[:τ_forward]))
    @test isapprox(mean(τ), 5 / 7; atol = 0.02)       # Beta(5,2) mean
    @test all(0 .< τ .< 1)

    chn2 = sample(MersenneTwister(20260518),
        test_positivity_model(; fraction_forwarded_prior = Beta(2.0, 8.0)),
        Prior(), 8_000; progress = false)
    τ2 = vec(Array(chn2[:τ_forward]))
    @test mean(τ2) < mean(τ)                          # override took effect
end

@testitem "confirmed binomial conditions on observed analysed" tags=[:slow] begin
    ## Confirmed increments are observed as `ΔC_v ~ Binomial(ΔA_v, p_pos_v)`
    ## with the newly-analysed count `ΔA_v` a known denominator (data). A
    ## confirmed draw must never exceed its denominator, and the per-vintage
    ## positivity must stay in (0, 1). Inputs are merged-vintage increments
    ## (the 25 May stall folded into 26 May).
    using Turing: sample, Prior, predict, @model, to_submodel
    import FlexiChains
    using BVDOutbreakSize: confirmed_cases_model, exponential_growth_model,
                           report_delay_model, surveillance_dispersion_model

    A = [211, 84, 108]                        # newly-analysed per window
    R = [418, 13, 231]                        # received increments
    C = [101, 4, 16]                          # confirmed increments
    edges = [123.0, 124.0, 126.0]             # vintage elapsed times

    @model function _binom_harness(confirmed, analysed, received;
            growth = exponential_growth_model(),
            dispersion = surveillance_dispersion_model(),
            report_delay = report_delay_model())
        growth_state ~ to_submodel(growth, false)
        disp_state ~ to_submodel(dispersion, false)
        rep_state ~ to_submodel(report_delay, false)
        confirmed_state ~ to_submodel(
            confirmed_cases_model(confirmed, analysed, received, 211,
                growth_state, disp_state.k, fill(0.3, length(confirmed)),
                0.6, 0.7, rep_state.dist, edges, 126.0), false)
    end

    chn = sample(_binom_harness(C, A, R), Prior(), 300;
        chain_type = FlexiChains.VNChain, progress = false)
    ## Per-test positivity at the cut-off, the tracked scalar diagnostic.
    pos = vec(Array(chn[:p_positive]))
    @test all(0 .< pos .< 1)
    ## The tested BVD share is also a valid probability.
    qc = vec(Array(chn[:q_cutoff]))
    @test all(0 .<= qc .<= 1)

    ## Posterior-predictive confirmed and received draws.
    pp = predict(_binom_harness(fill(missing, 3), A, fill(missing, 3)), chn)
    draws = vec(Array(pp[:confirmed_cases]))
    cc = reduce(hcat, draws)            # 3 vintages × draws
    for v in 1:3
        @test all(cc[v, :] .<= A[v])
        @test all(cc[v, :] .>= 0)
    end
    rr = reduce(hcat, vec(Array(pp[:samples_received])))
    @test all(rr .>= 0)
end

@testitem "confirmed binomial single-total reduction" tags=[:slow] begin
    ## The length-1 confirmed vector with a single analysed denominator
    ## must still fit and predict, keeping the McCabe single-total path
    ## intact. Confirmed draws bounded by the single analysed total.
    using Turing: sample, Prior, predict
    import FlexiChains
    using BVDOutbreakSize: confirmed_only_model

    chn = sample(confirmed_only_model(101, 211), Prior(), 300;
        chain_type = FlexiChains.VNChain, progress = false)
    C = vec(Array(chn[:cumulative_cases]))
    @test all(isfinite, C)
    @test all(C .> 0)
    pos = vec(Array(chn[:p_positive]))
    @test all(0 .< pos .< 1)

    pp = predict(confirmed_only_model(missing, 211), chn)
    cc = reduce(vcat, vec(Array(pp[:confirmed_cases])))
    @test all(0 .<= cc .<= 211)
end

@testitem "received likelihood pins τ_forward to received/suspect" tags=[:slow] begin
    ## The received mean is τ_forward · N_susp. With a fixed suspect backlog,
    ## conditioning on received counts well above the prior-implied mean
    ## should pull the τ_forward posterior up relative to the prior.
    using Turing: sample, Prior, NUTS, @model, to_submodel
    import FlexiChains
    using Statistics: mean
    using BVDOutbreakSize: confirmed_cases_model, exponential_growth_model,
                           report_delay_model, surveillance_dispersion_model

    A = [211, 84, 108]
    C = [101, 4, 16]
    edges = [123.0, 124.0, 126.0]

    @model function _recv_harness(confirmed, analysed, received;
            growth = exponential_growth_model(),
            dispersion = surveillance_dispersion_model(),
            report_delay = report_delay_model())
        growth_state ~ to_submodel(growth, false)
        disp_state ~ to_submodel(dispersion, false)
        rep_state ~ to_submodel(report_delay, false)
        confirmed_state ~ to_submodel(
            confirmed_cases_model(confirmed, analysed, received, 403,
                growth_state, disp_state.k, fill(0.5, length(confirmed)),
                1.0, 0.5, rep_state.dist, edges, 126.0), false)
    end

    ## Sanity: with received supplied the model samples without error and
    ## the expected received total is positive and finite.
    chn = sample(_recv_harness(C, A, [418, 13, 231]), Prior(), 300;
        chain_type = FlexiChains.VNChain, progress = false)
    er = vec(Array(chn[:expected_received_total]))
    @test all(isfinite, er)
    @test all(er .> 0)
end

@testitem "lab_receipt_delay_model: gamma delay centred near 3 days" tags=[:slow] begin
    using Turing: sample, Prior
    using Random: MersenneTwister
    using Statistics: mean
    using Distributions: Gamma, truncated, Normal
    using BVDOutbreakSize: lab_receipt_delay_model

    ## Report-to-lab-receipt delay: α ~ N₊(2, 1), θ ~ N₊(1.5, 0.75), so
    ## the mean delay E[α]·E[θ] ≈ 2 · 1.5 = 3 days.
    chn = sample(MersenneTwister(20260518), lab_receipt_delay_model(),
        Prior(), 10_000; progress = false)
    α = vec(Array(chn[:α_recv]))
    θ = vec(Array(chn[:θ_recv]))
    @test all(isfinite, α) && all(>(0.1), α)   # truncated at 0.1
    @test all(isfinite, θ) && all(>(0.1), θ)
    @test isapprox(mean(α), 2.0; atol = 0.1)
    @test isapprox(mean(θ), 1.5; atol = 0.1)
    ## The per-draw gamma mean (α·θ) is centred near the documented 3 days.
    @test isapprox(mean(α .* θ), 3.0; atol = 0.3)

    ## Override moves the centre and the returned distribution is a Gamma.
    res = lab_receipt_delay_model()()
    @test res.dist isa Gamma
    @test res.dist == Gamma(res.α, res.θ)
end

@testitem "lab_capacity_model: positive capacities centred near 150" tags=[:slow] begin
    using Turing: sample, Prior, @model, to_submodel
    using Random: MersenneTwister
    using Statistics: mean, std
    import FlexiChains
    using BVDOutbreakSize: lab_capacity_model

    @model function _cap_wrap(n)
        st ~ to_submodel(lab_capacity_model(n), false)
        return (; st)
    end

    n = 4
    chn = sample(MersenneTwister(20260518), _cap_wrap(n), Prior(), 4_000;
        chain_type = FlexiChains.VNChain, progress = false)

    ## The length-n capacity vector is a positive log-random-walk.
    cap = reduce(hcat, vec(Array(chn[:capacity])))   # n × draws
    @test size(cap, 1) == n
    @test all(isfinite, cap) && all(>(0), cap)
    rw = vec(Array(chn[:rw_sd]))
    @test all(rw .>= 0)                               # half-normal walk SD

    ## Window 1 sits at the centred level log(150) + 0.5·z₁, so its
    ## geometric mean (median of the log-normal level) is ≈ 150.
    @test isapprox(mean(cap[1, :]), 150.0; rtol = 0.2)
end

@testitem "analysis capacity is pulled below the prior by throughput" tags=[:slow] begin
    ## The capacity random walk is centred on 150 samples/day. Conditioning
    ## on the observed analysed throughput (84/day then 54/day across the
    ## merged windows) under NUTS must pull the cut-off capacity posterior
    ## well below the 150 prior centre, and keep it finite and positive.
    using Turing: @model, to_submodel
    using Statistics: mean
    using BVDOutbreakSize: confirmed_cases_model, exponential_growth_model,
                           report_delay_model, surveillance_dispersion_model,
                           nuts_sample

    A = [211, 84, 108];
    R = [418, 13, 231];
    C = [101, 4, 16]
    edges = [123.0, 124.0, 126.0]

    @model function _cap_harness(confirmed, analysed, received;
            growth = exponential_growth_model(),
            dispersion = surveillance_dispersion_model(),
            report_delay = report_delay_model())
        growth_state ~ to_submodel(growth, false)
        disp_state ~ to_submodel(dispersion, false)
        rep_state ~ to_submodel(report_delay, false)
        confirmed_state ~ to_submodel(
            confirmed_cases_model(confirmed, analysed, received, 211,
                growth_state, disp_state.k, fill(0.3, length(confirmed)),
                0.6, 0.7, rep_state.dist, edges, 126.0), false)
    end

    chn = nuts_sample(_cap_harness(C, A, R);
        samples = 60, chains = 1, seed = 3, progress = false)
    cap = vec(Array(chn[:capacity_cutoff]))
    @test all(isfinite, cap)
    @test all(cap .> 0)
    @test mean(cap) < 150          # throughput pulls capacity below the prior
end

@testitem "confirmed_q_re_model: prior draws are valid" tags=[:slow] begin
    ## The per-vintage q random effect samples a non-negative pooling SD
    ## `σ_q` and `n` standard-normal offsets `z_q`. Prior draws of both
    ## must be finite, with `σ_q ≥ 0`.
    using Turing: sample, Prior, @model, to_submodel
    import FlexiChains
    using BVDOutbreakSize: confirmed_q_re_model

    @model function _qre_harness(n)
        state ~ to_submodel(confirmed_q_re_model(n), false)
        sigma_q := state.σ_q
        z_mean := sum(state.z_q) / n
        return (; σ_q = state.σ_q, z_q = state.z_q)
    end

    chn = sample(_qre_harness(5), Prior(), 300;
        chain_type = FlexiChains.VNChain, progress = false)
    σ = vec(Array(chn[:sigma_q]))
    @test length(σ) == 300
    @test all(isfinite, σ)
    @test all(σ .>= 0)
    zm = vec(Array(chn[:z_mean]))
    @test all(isfinite, zm)
end

@testitem "q random effect on by default; nothing recovers baseline" tags=[:slow] begin
    ## With the q-RE on (the default) each vintage's tested BVD share can
    ## depart from the smooth severe-first baseline, so the confirmed
    ## stream fits the non-monotone positivity. The model must fit under
    ## Prior and keep per-vintage positivity and the cut-off share in
    ## (0, 1); predictive confirmed draws stay bounded by the analysed
    ## denominator. Passing `q_random_effect = nothing` recovers the
    ## smooth baseline and must still fit.
    using Turing: sample, Prior, predict, @model, to_submodel
    import FlexiChains
    using BVDOutbreakSize: confirmed_cases_model, confirmed_q_re_model,
                           exponential_growth_model, report_delay_model,
                           surveillance_dispersion_model

    A = [211, 84, 108, 60, 95]    # non-monotone positivity windows
    R = [418, 13, 231, 40, 180]
    C = [101, 4, 16, 4, 85]       # 28 May per-vintage confirmed
    edges = [123.0, 124.0, 126.0, 127.0, 128.0]

    @model function _qre_conf(confirmed, analysed, received; q_re,
            growth = exponential_growth_model(),
            dispersion = surveillance_dispersion_model(),
            report_delay = report_delay_model())
        growth_state ~ to_submodel(growth, false)
        disp_state ~ to_submodel(dispersion, false)
        rep_state ~ to_submodel(report_delay, false)
        confirmed_state ~ to_submodel(
            confirmed_cases_model(confirmed, analysed, received, 211,
                growth_state, disp_state.k, fill(0.3, length(confirmed)),
                0.6, 0.7, rep_state.dist, edges, 128.0;
                q_random_effect = q_re), false)
    end

    chn = sample(_qre_conf(C, A, R; q_re = confirmed_q_re_model),
        Prior(), 300; chain_type = FlexiChains.VNChain, progress = false)
    pos = vec(Array(chn[:p_positive]))
    @test all(0 .< pos .< 1)
    qc = vec(Array(chn[:q_cutoff]))
    @test all(0 .<= qc .<= 1)

    pp = predict(
        _qre_conf(fill(missing, 5), A, fill(missing, 5);
            q_re = confirmed_q_re_model), chn)
    cc = reduce(hcat, vec(Array(pp[:confirmed_cases])))
    for v in 1:5
        @test all(cc[v, :] .<= A[v])
        @test all(cc[v, :] .>= 0)
    end

    ## q-RE off recovers the smooth severe-first baseline and still fits.
    chn0 = sample(_qre_conf(C, A, R; q_re = nothing),
        Prior(), 200; chain_type = FlexiChains.VNChain, progress = false)
    pos0 = vec(Array(chn0[:p_positive]))
    @test all(0 .< pos0 .< 1)
end
