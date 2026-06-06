## Tests for the test-positivity submodel from `src/models/priors.jl`.
## The `λ_bg` prior was retuned to a truncated Normal with its CENTRE
## above zero (`truncated(Normal(0.5, 0.3); lower = 0)`) so the background
## non-BVD suspected-case process pins at a sane non-zero rate rather than
## collapsing to zero, while still staying a minority of observed cases.

@testitem "default λ_bg prior is above-zero truncated Normal" tags=[:slow] begin
    using Turing: sample, Prior
    using Random: MersenneTwister
    using Statistics: mean, median
    using Distributions: truncated, Normal
    import Distributions
    using BVDOutbreakSize: test_positivity_model

    ## The retuned default is `truncated(Normal(0.5, 0.3); lower = 0)`,
    ## a half-normal with its centre shifted above zero. Compare the
    ## sampled mean against the analytic truncated-Normal mean and check
    ## the centre is genuinely above zero.
    ref = truncated(Normal(0.5, 0.3); lower = 0)
    chn = sample(MersenneTwister(20260518), test_positivity_model(),
        Prior(), 40_000; progress = false)
    λ_bg = vec(Array(chn[:λ_bg]))
    @test isapprox(mean(λ_bg), Distributions.mean(ref); atol = 0.03)
    @test median(λ_bg) > 0.3
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
