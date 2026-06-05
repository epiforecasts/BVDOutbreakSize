## Tests for the test-positivity submodel from `src/models/priors.jl`.
## The `λ_bg` prior is a diffuse half-normal `Normal+(0, 5)`. The
## background non-BVD suspected-case rate is no longer pinned down by a
## tight prior; the composition link ties suspected positivity to the
## background through the laboratory positivity data, so identification
## comes from the joint fit rather than from a prior that forces the
## background to a minority share. The prior is therefore deliberately
## permissive and admits a background comparable to the observed total.

@testitem "default λ_bg prior matches half-normal SD 5" tags=[:slow] begin
    using Turing: sample, Prior
    using Random: MersenneTwister
    using Statistics: mean, std
    using BVDOutbreakSize: test_positivity_model

    ## The default is `truncated(Normal(0, 5); lower = 0)`. Fold a
    ## half-normal back to its untruncated SD via the known moments:
    ## E|X| = σ√(2/π), so σ = mean·√(π/2), and check the SD.
    chn = sample(MersenneTwister(20260518), test_positivity_model(),
        Prior(), 40_000; progress = false)
    λ_bg = vec(Array(chn[:λ_bg]))
    @test isapprox(mean(λ_bg) * sqrt(pi / 2), 5.0; atol = 0.15)
    @test isapprox(std(λ_bg), 5.0 * sqrt(1 - 2 / pi); atol = 0.15)
end

@testitem "λ_bg prior is permissive, identified downstream" tags=[:slow] begin
    using Turing: sample, Prior
    using Random: MersenneTwister
    using Statistics: median, quantile
    using BVDOutbreakSize: test_positivity_model, load_observations

    ## The background contribution over the window is `λ_bg · T`, with T
    ## the latent seeding-to-cut-off time (≈ 132 days on current data).
    ## The diffuse prior deliberately admits a background reaching the
    ## observed scale; the composition link, not the prior, is what
    ## constrains it in the joint fit.
    obs = load_observations()
    observed_total = obs.reported_cases  # cumulative suspected at cut-off
    T = 132.0

    chn = sample(MersenneTwister(20260518), test_positivity_model(),
        Prior(), 40_000; progress = false)
    λ_bg = vec(Array(chn[:λ_bg]))
    background = λ_bg .* T

    ## A genuine non-BVD signal (not forced to zero) and a prior broad
    ## enough that the 95% prior-predictive background can reach the
    ## observed suspected total — identification is deferred to the joint.
    @test median(λ_bg) > 0.3
    @test quantile(background, 0.95) > 0.5 * observed_total
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
