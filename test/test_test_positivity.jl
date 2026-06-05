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
