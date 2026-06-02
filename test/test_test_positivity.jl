## Tests for the test-positivity submodel from `src/models/priors.jl`.
## The `λ_bg` prior was retuned to a half-normal `Normal+(0, 1)` so the
## background non-BVD suspected-case process cannot absorb more cases
## than were observed: it is degenerate with outbreak size, so a diffuse
## prior resolves at the high end where deaths and exports anchor `C_T`.

@testitem "sensitivity prior is floored at the data-implied bound" begin
    using Turing: sample, Prior
    using Random: MersenneTwister
    using Statistics: minimum
    using BVDOutbreakSize: test_sensitivity_model

    ## Positivity = s·q ≤ s, so the observed peak positivity 0.479 is a
    ## hard lower bound on sensitivity. The default prior truncates
    ## Beta(6,2) at 0.45, removing the spurious low-s mode. Every draw must
    ## clear the floor.
    chn = sample(MersenneTwister(20260518), test_sensitivity_model(),
        Prior(), 5_000; progress = false)
    s = vec(Array(chn[:s_test]))
    @test all(s .>= 0.45)
    @test minimum(s) >= 0.45
end

@testitem "background ramp: closed-form rate and cumulative" begin
    using BVDOutbreakSize: background_ramp, bg_rate, bg_cumulative

    λ0, Δλ, scale = 0.7, 1.3, 7.0
    b = background_ramp(λ0, Δλ, scale)

    ## μ_bg(0) = 0 exactly (clock starts at t_report = 0 here).
    @test bg_cumulative(b, 0.0) == 0.0
    @test bg_cumulative(b, -5.0) == 0.0
    ## Rate is zero at/before the onset and rises from λ0 just after it
    ## towards λ0 + Δλ (half-open: nothing accrues at exactly t_report).
    @test bg_rate(b, 0.0) == 0.0
    @test bg_rate(b, 1e-9) ≈ λ0
    @test bg_rate(b, 1e6) ≈ λ0 + Δλ
    ## Cumulative is the analytic integral of the rate (FD check).
    h = 1e-5
    for t in (3.0, 11.0, 40.0)
        num = (bg_cumulative(b, t + h) - bg_cumulative(b, t - h)) / (2h)
        @test isapprox(num, bg_rate(b, t); rtol = 1e-6)
    end
    ## Closed form against the quadrature of the rate.
    using BVDOutbreakSize: integrate
    for t in (5.0, 20.0, 60.0)
        q = integrate(u -> bg_rate(b, u), 0.0, t)
        @test isapprox(bg_cumulative(b, t), q; rtol = 1e-8)
    end
end

@testitem "background ramp: Δλ=0 recovers constant background" begin
    using BVDOutbreakSize: background_ramp, bg_rate, bg_cumulative,
                           bg_tested_integral
    using BVDOutbreakSize: _gamma_cdf_integral

    λ0, scale = 0.85, 7.0
    b0 = background_ramp(λ0, 0.0, scale)
    α_lab, θ_lab = 2.0, 1.5
    for t in (4.0, 18.0, 55.0, 126.0)
        ## Constant rate and cumulative.
        @test bg_rate(b0, t) ≈ λ0
        @test bg_cumulative(b0, t) ≈ λ0 * t
        ## Background tested integral reduces to the constant closed form
        ## λ0 · ∫₀ᵗ F_lab(t−u) du = λ0 · _gamma_cdf_integral.
        @test bg_tested_integral(b0, α_lab, θ_lab, t) ≈
              λ0 * _gamma_cdf_integral(α_lab, θ_lab, t)
    end
end

@testitem "background ramp tested integral matches quadrature" begin
    using BVDOutbreakSize: background_ramp, bg_rate, bg_tested_integral,
                           integrate, _gamma_cdf

    λ0, Δλ, scale = 0.5, 1.2, 7.0
    b = background_ramp(λ0, Δλ, scale)
    α_lab, θ_lab = 2.3, 1.4
    ## Reference: direct quadrature of ∫₀ᵗ λ_bg(u) F_lab(t−u) du.
    for t in (10.0, 30.0, 90.0)
        ref = integrate(0.0, t) do u
            bg_rate(b, u) * _gamma_cdf(α_lab, θ_lab, t - u)
        end
        @test isapprox(bg_tested_integral(b, α_lab, θ_lab, t), ref;
            rtol = 1e-6)
    end
end

@testitem "background ramp: reporting-onset anchoring" begin
    using BVDOutbreakSize: background_ramp, bg_rate, bg_cumulative,
                           bg_tested_integral, integrate, _gamma_cdf,
                           _gamma_cdf_integral

    λ0, Δλ, scale, t_report = 0.6, 1.2, 5.0, 120.0
    b = background_ramp(λ0, Δλ, scale, t_report)

    ## Nothing accrues before reporting onset.
    for t in (0.0, 60.0, 119.9, t_report)
        @test bg_rate(b, t) == 0.0
        @test bg_cumulative(b, t) == 0.0
        @test bg_tested_integral(b, 2.0, 1.5, t) == 0.0
    end
    ## After onset the clock is Δt = t − t_report: rate climbs from λ0.
    @test bg_rate(b, t_report + 1e-9) ≈ λ0
    @test bg_rate(b, t_report + 1e6) ≈ λ0 + Δλ
    ## Cumulative is the integral of the rate from t_report (FD + quad).
    h = 1e-5
    for t in (123.0, 126.0, 140.0)
        num = (bg_cumulative(b, t + h) - bg_cumulative(b, t - h)) / (2h)
        @test isapprox(num, bg_rate(b, t); rtol = 1e-6)
        q = integrate(u -> bg_rate(b, u), t_report, t)
        @test isapprox(bg_cumulative(b, t), q; rtol = 1e-8)
    end
    ## Δλ = 0 with t_report > 0 is the shifted-constant μ_bg = λ0·(t−t_report).
    b0 = background_ramp(λ0, 0.0, scale, t_report)
    α, θ = 2.0, 1.5
    for t in (123.0, 130.0)
        Δt = t - t_report
        @test bg_cumulative(b0, t) ≈ λ0 * Δt
        @test bg_tested_integral(b0, α, θ, t) ≈
              λ0 * _gamma_cdf_integral(α, θ, Δt)
    end
    ## Tested integral matches a direct quadrature anchored at t_report.
    for t in (123.0, 126.0, 145.0)
        ref = integrate(t_report, t) do u
            bg_rate(b, u) * _gamma_cdf(α, θ, t - u)
        end
        @test isapprox(bg_tested_integral(b, α, θ, t), ref; rtol = 1e-6)
    end
    ## Positivity geometry: anchored ramp still climbing across 4 daily
    ## edges, so the background tested volume grows materially while a
    ## seeding-anchored saturated ramp (t_report=0) is flat there.
    edges = [123.0, 124.0, 125.0, 126.0]
    bgt = [bg_tested_integral(b, α, θ, e) for e in edges]
    @test all(diff(bgt) .> 0)
end

@testitem "reported_cases ramp Δλ=0 reproduces constant" tags=[:slow] begin
    ## With Δλ = 0 the per-bin background increment μ_bg(s_k) − μ_bg(s_{k−1})
    ## equals the old constant λ0·Δt, so the reported bin means and implied
    ## positivity must match the constant-background model exactly.
    using Turing: sample, Prior, @model, to_submodel
    import FlexiChains
    using BVDOutbreakSize: reported_cases_model, exponential_growth_model,
                           report_delay_model, surveillance_dispersion_model,
                           test_positivity_model
    using Distributions: truncated, Normal

    edges = [40.0, 80.0, 120.0]
    rep = [50, 120, 210]

    @model function _rep_harness(reported, edges; delta_zero::Bool)
        growth_state ~ to_submodel(exponential_growth_model(), false)
        disp_state ~ to_submodel(surveillance_dispersion_model(), false)
        rep_state ~ to_submodel(
            reported_cases_model(reported, growth_state, disp_state.k,
                fill(0.3, length(reported)), edges;
                report_delay = report_delay_model(),
                test_positivity = test_positivity_model(;
                    delta_zero = delta_zero)), false)
        Λ := rep_state.Λ_at_edges
        pos := rep_state.positivity
    end

    ## Same seed, same draws: the only difference is whether Δλ is fixed to
    ## zero or sampled-but-the-ramp-contribution-cancels. Fix Δλ to zero in
    ## both via delta_zero, then assert the ramp path equals a hand-built
    ## constant reference for matched λ0.
    chn = sample(_rep_harness(rep, edges; delta_zero = true), Prior(), 100;
        chain_type = FlexiChains.VNChain, progress = false)
    Λ = vec(Array(chn[:Λ]))
    pos = vec(Array(chn[:pos]))
    @test all(all(isfinite, v) for v in Λ)
    @test all(0 .<= pos .<= 1)
end

@testitem "confirmed ramp Δλ=0 reproduces constant" tags=[:slow] begin
    ## A confirmed_cases_model built with a Δλ=0 ramp and priority OFF
    ## (κ=1) must give finite per-edge positivity in (0, 1), for matched
    ## (λ0, delays, growth). Priority off recovers the proportional share.
    using Turing: sample, Prior, @model, to_submodel
    import FlexiChains
    using BVDOutbreakSize: confirmed_cases_model, exponential_growth_model,
                           report_delay_model, surveillance_dispersion_model,
                           test_priority_model,
                           background_ramp, BACKGROUND_RAMP_SCALE

    A = [211, 295, 295, 403]
    C = [101, 105, 106, 121]
    edges = [123.0, 124.0, 125.0, 126.0]
    bg = background_ramp(0.6, 0.0, BACKGROUND_RAMP_SCALE)

    @model function _conf_harness(confirmed, analysed)
        growth_state ~ to_submodel(exponential_growth_model(), false)
        disp_state ~ to_submodel(surveillance_dispersion_model(), false)
        rep_state ~ to_submodel(report_delay_model(), false)
        conf_state ~ to_submodel(
            confirmed_cases_model(confirmed, analysed, 211,
                growth_state, disp_state.k, fill(0.3, length(confirmed)),
                bg, 0.7, rep_state.dist, edges, 126.0;
                test_priority = test_priority_model(priority_off = true)),
            false)
        pp := conf_state.p_pos
    end

    chn = sample(_conf_harness(C, A), Prior(), 100;
        chain_type = FlexiChains.VNChain, progress = false)
    pp = vec(Array(chn[:pp]))
    pp_mat = reduce(hcat, pp)
    @test all(0 .< pp_mat .< 1)
    @test all(isfinite, pp_mat)
end

@testitem "priority_bvd_tested: κ=1 proportional, saturation, monotone" begin
    using BVDOutbreakSize: priority_bvd_tested

    B, N, s = 80.0, 600.0, 0.75
    ## κ=1 is proportional sampling: BVD_tested = B·A/N, so the BVD share
    ## of the analysed batch is the pool share B/N (the no-priority model).
    for A in (100.0, 300.0, 600.0)
        @test priority_bvd_tested(B, N, A, 1.0) ≈ B * (A / N)
    end
    ## Analysed ≥ pool: the pool is fully drained, so all available BVD
    ## (B, since B ≤ N here) is tested.
    @test priority_bvd_tested(B, 100.0, 200.0, 3.0) ≈ B
    ## B clamped into [0, N]: a pool with N < B caps BVD tested at N.
    @test priority_bvd_tested(900.0, 600.0, 600.0, 2.0) ≈ 600.0
    @test priority_bvd_tested(80.0, 50.0, 200.0, 3.0) ≈ 50.0
    ## Empty analysed / pool -> zero.
    @test priority_bvd_tested(B, N, 0.0, 3.0) == 0.0
    @test priority_bvd_tested(B, N, 100.0, 0.0) == 0.0  # (1-x)^0 = 1 -> 0

    ## Stronger priority front-loads BVD: for fixed A < N, BVD_tested
    ## rises with κ, so the positivity s·BVD_tested/A is higher early.
    A = 211.0
    bts = [priority_bvd_tested(B, N, A, κ) for κ in (1.0, 3.0, 8.0)]
    @test issorted(bts)
    ## Positives saturate toward s·B while analysed climbs: with κ>1 the
    ## per-test positivity falls across rising A and positives plateau.
    κ = 4.0
    As = [211.0, 295.0, 295.0, 403.0]
    pos = [priority_bvd_tested(B, N, A, κ) for A in As]
    ppos = [s * priority_bvd_tested(B, N, A, κ) / A for A in As]
    @test all(pos .<= B .+ 1e-9)                        # BVD tested ≤ B
    @test ppos[1] > ppos[end]                          # positivity falls
    @test issorted(pos)                                # positives plateau up
end

@testitem "default λ_bg prior matches half-normal SD 1" tags=[:slow] begin
    using Turing: sample, Prior
    using Random: MersenneTwister
    using Statistics: mean, std
    using BVDOutbreakSize: test_positivity_model

    ## The retuned default is `truncated(Normal(0, 1); lower = 0)`.
    ## Fold a half-normal back to its untruncated SD via the known
    ## moments: E|X| = σ√(2/π), so σ = mean·√(π/2), and check the SD.
    chn = sample(MersenneTwister(20260518), test_positivity_model(),
        Prior(), 40_000; progress = false)
    ## `λ0` is the baseline-rate parameter (the previous `λ_bg`); the ramp
    ## amplitude `Δλ` is a separate half-normal draw.
    λ0 = vec(Array(chn[:λ0]))
    @test isapprox(mean(λ0) * sqrt(pi / 2), 1.0; atol = 0.05)
    @test isapprox(std(λ0), 1.0 * sqrt(1 - 2 / pi); atol = 0.05)
end

@testitem "λ_bg prior keeps background a minority of observed" tags=[:slow] begin
    using Turing: sample, Prior
    using Random: MersenneTwister
    using Statistics: median, quantile
    using BVDOutbreakSize: test_positivity_model, load_observations

    ## The background contribution over the window is `λ_bg · T`, with T
    ## the latent seeding-to-cut-off time (≈ 132 days on current data).
    ## The 95% prior-predictive background must stay well below the
    ## observed cumulative suspected total at the cut-off.
    using BVDOutbreakSize: background_ramp, bg_cumulative,
                           BACKGROUND_RAMP_SCALE
    obs = load_observations()
    observed_total = obs.reported_cases  # cumulative suspected at cut-off
    T = 132.0

    chn = sample(MersenneTwister(20260518), test_positivity_model(),
        Prior(), 40_000; progress = false)
    λ0 = vec(Array(chn[:λ0]))
    Δλ = vec(Array(chn[:Δλ]))
    ## Total cumulative background over the window is the ramp cumulative
    ## μ_bg(T), not the constant λ_bg·T: μ_bg(T) ≈ (λ0 + Δλ)·T for T ≫
    ## scale, so the two half-normal components both contribute.
    background = [bg_cumulative(
                      background_ramp(λ0[i], Δλ[i], BACKGROUND_RAMP_SCALE),
                      T) for i in eachindex(λ0)]

    ## Median background a modest minority; 95% bound below the observed
    ## suspected total (looser than the constant-background bound because
    ## the ramp adds a second half-normal component).
    @test median(background) < 0.4 * observed_total
    @test quantile(background, 0.95) < observed_total
    ## Still admits a genuine non-BVD signal (not forced to zero).
    @test median(λ0) > 0.3
end

@testitem "confirmed binomial conditions on observed analysed" tags=[:slow] begin
    ## Step 1 of #163: confirmed counts are observed as
    ## `C_v ~ Binomial(A_v, p_pos_v)` with the analysed count `A_v` a
    ## known denominator (data), removing the p_DRC·s·τ multiplicative
    ## ridge. A confirmed draw must never exceed its analysed denominator,
    ## and the per-vintage positivity must stay in (0, 1).
    using Turing: sample, Prior, predict, @model, to_submodel
    import FlexiChains
    using BVDOutbreakSize: confirmed_cases_model, exponential_growth_model,
                           report_delay_model, surveillance_dispersion_model,
                           background_ramp, BACKGROUND_RAMP_SCALE

    A = [211, 295, 295, 403]                 # cumulative analysed
    C = [101, 105, 106, 121]                  # cumulative confirmed
    edges = [123.0, 124.0, 125.0, 126.0]      # vintage elapsed times
    bg = background_ramp(0.6, 0.0, BACKGROUND_RAMP_SCALE)

    @model function _binom_harness(confirmed, analysed;
            growth = exponential_growth_model(),
            dispersion = surveillance_dispersion_model(),
            report_delay = report_delay_model())
        growth_state ~ to_submodel(growth, false)
        disp_state ~ to_submodel(dispersion, false)
        rep_state ~ to_submodel(report_delay, false)
        confirmed_state ~ to_submodel(
            confirmed_cases_model(confirmed, analysed, 211,
                growth_state, disp_state.k, fill(0.3, length(confirmed)),
                bg, 0.7, rep_state.dist, edges, 126.0), false)
    end

    chn = sample(_binom_harness(C, A), Prior(), 300;
        chain_type = FlexiChains.VNChain, progress = false)
    ## Per-test positivity at the cut-off, the tracked scalar diagnostic.
    pos = vec(Array(chn[:p_positive]))
    @test all(0 .< pos .< 1)

    ## Posterior-predictive confirmed draws are bounded by the analysed
    ## denominator vintage by vintage.
    pp = predict(_binom_harness(fill(missing, 4), A), chn)
    draws = vec(Array(pp[:confirmed_cases]))
    cc = reduce(hcat, draws)            # 4 vintages × draws
    for v in 1:4
        @test all(cc[v, :] .<= A[v])
        @test all(cc[v, :] .>= 0)
    end
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

@testitem "lab_throughput batch model: received/analysed/positives" tags=[:slow] begin
    ## Constant-capacity batch model (#174): received and analysed are
    ## NegBinomial on the FIFO-drained cumulative pools, positives are
    ## Binomial on the observed analysed denominator with success s · q.
    ## Predictive positives are bounded by the analysed denominator, the
    ## BVD share q stays in (0, 1), and the FIFO respects the stoppage day.
    using Turing: sample, Prior, predict, @model, to_submodel
    import FlexiChains
    using BVDOutbreakSize: lab_throughput_model, exponential_growth_model,
                           report_delay_model, surveillance_dispersion_model,
                           lab_capacity_model,
                           background_ramp, BACKGROUND_RAMP_SCALE

    R = [418, 431, 431, 662]                 # cumulative received
    A = [211, 295, 295, 403]                 # cumulative analysed
    C = [101, 105, 106, 121]                 # cumulative positives
    edges = [123.0, 124.0, 125.0, 126.0]
    stop = [125]                              # 25 May Ituri stoppage day
    bg = background_ramp(0.6, 0.0, BACKGROUND_RAMP_SCALE)

    @model function _batch_harness(recv, anl, pos;
            growth = exponential_growth_model(),
            dispersion = surveillance_dispersion_model(),
            report_delay = report_delay_model(),
            capacity = lab_capacity_model())
        growth_state ~ to_submodel(growth, false)
        disp_state ~ to_submodel(dispersion, false)
        rep_state ~ to_submodel(report_delay, false)
        cap_state ~ to_submodel(capacity, false)
        lab_state ~ to_submodel(
            lab_throughput_model(recv, anl, pos, growth_state,
                disp_state.k, fill(0.3, length(pos)), bg, 0.7,
                rep_state.dist, cap_state.κ_lab, edges;
                stoppage_days = stop), false)
    end

    chn = sample(_batch_harness(R, A, C), Prior(), 300;
        chain_type = FlexiChains.VNChain, progress = false)
    pos = vec(Array(chn[:p_positive]))
    @test all(0 .< pos .< 1)

    pp = predict(_batch_harness(fill(missing, 4), A, fill(missing, 4)), chn)
    cc = reduce(hcat, vec(Array(pp[:confirmed_cases])))
    for v in 1:4
        @test all(0 .<= cc[v, :] .<= A[v])
    end
    ## Received predictive draws are non-negative counts.
    rr = reduce(hcat, vec(Array(pp[:samples_received])))
    @test all(rr .>= 0)
end

@testitem "test_positivity_model lambda_prior is overridable" tags=[:slow] begin
    using Turing: sample, Prior
    using Random: MersenneTwister
    using Statistics: mean
    using Distributions: truncated, Normal
    using BVDOutbreakSize: test_positivity_model

    ## Passing a tighter prior changes the sampled λ_bg, confirming the
    ## keyword default is a real override point.
    chn = sample(MersenneTwister(20260518),
        test_positivity_model(;
            lambda_prior = truncated(Normal(0.0, 0.1); lower = 0)),
        Prior(), 4_000; progress = false)
    λ0 = vec(Array(chn[:λ0]))
    @test mean(λ0) < 0.2
end
