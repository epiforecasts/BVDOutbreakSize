## Tests for the test-positivity submodel from `src/models/priors.jl`.
## The `λ_bg` prior was retuned to a half-normal `Normal+(0, 1)` so the
## background non-BVD suspected-case process cannot absorb more cases
## than were observed: it is degenerate with outbreak size, so a diffuse
## prior resolves at the high end where deaths and exports anchor `C_T`.

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
    λ_bg = vec(Array(chn[:λ_bg]))
    @test isapprox(mean(λ_bg) * sqrt(pi / 2), 1.0; atol = 0.05)
    @test isapprox(std(λ_bg), 1.0 * sqrt(1 - 2 / pi); atol = 0.05)
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
    obs = load_observations()
    observed_total = obs.reported_cases  # cumulative suspected at cut-off
    T = 132.0

    chn = sample(MersenneTwister(20260518), test_positivity_model(),
        Prior(), 40_000; progress = false)
    λ_bg = vec(Array(chn[:λ_bg]))
    background = λ_bg .* T

    ## Median background a modest minority; 95% bound clearly below the
    ## observed suspected total.
    @test median(background) < 0.2 * observed_total
    @test quantile(background, 0.95) < 0.5 * observed_total
    ## Still admits a genuine non-BVD signal (not forced to zero).
    @test median(λ_bg) > 0.3
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
                           report_delay_model, surveillance_dispersion_model

    A = [211, 295, 295, 403]                 # cumulative analysed
    C = [101, 105, 106, 121]                  # cumulative confirmed
    edges = [123.0, 124.0, 125.0, 126.0]      # vintage elapsed times

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
                0.6, 0.7, rep_state.dist, edges, 126.0), false)
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
                           lab_capacity_model

    R = [418, 431, 431, 662]                 # cumulative received
    A = [211, 295, 295, 403]                 # cumulative analysed
    C = [101, 105, 106, 121]                 # cumulative positives
    edges = [123.0, 124.0, 125.0, 126.0]
    stop = [125]                              # 25 May Ituri stoppage day

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
                disp_state.k, fill(0.3, length(pos)), 0.6, 0.7,
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
    λ_bg = vec(Array(chn[:λ_bg]))
    @test mean(λ_bg) < 0.2
end
