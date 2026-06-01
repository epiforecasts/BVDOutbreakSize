## Tests for the infection layer: the incubation-period prior and the
## onset rescale that maps the latent cumulative infections onto the
## cumulative symptom onsets the downstream delays act on. Under
## exponential growth the infection-to-onset convolution is the exact
## constant factor `mgf(incubation, −r)`, so `cumulative_cases =
## onset_fraction · cumulative_infections` holds per draw.

@testitem "incubation_model reproduces the MacNeil 2010 mean CI" tags=[:slow] begin
    using Turing: sample, Prior
    using Random: MersenneTwister
    import Statistics
    using BVDOutbreakSize: incubation_model

    chn = sample(MersenneTwister(20260518), incubation_model(),
        Prior(), 20_000; progress = false)
    mean_inc = vec(Array(chn[:mean_inc]))
    cv = vec(Array(chn[:cv_inc]))
    α = vec(Array(chn[:α_inc]))
    θ = vec(Array(chn[:θ_inc]))

    @test all(isfinite, α) && all(>(0), α)
    @test all(isfinite, θ) && all(>(0), θ)
    @test all(>(0), cv)
    ## The prior on the mean carries MacNeil et al. (2010)'s reported
    ## uncertainty: mean 6.3 d, 95% CI 5.2-7.3. Assert the centre and the
    ## CI endpoints are reproduced (generous tolerance for sampling).
    @test isapprox(Statistics.mean(mean_inc), 6.3; atol = 0.1)
    @test isapprox(Statistics.quantile(mean_inc, 0.025), 5.2; atol = 0.3)
    @test isapprox(Statistics.quantile(mean_inc, 0.975), 7.3; atol = 0.3)
    ## Shape/scale are recovered consistently from mean and CV.
    @test all(@. isapprox(α * θ, mean_inc; rtol = 1e-8))
end

@testitem "onset_rescale equals the incubation mgf at −r" begin
    using Distributions: Gamma, mgf
    using BVDOutbreakSize: onset_rescale

    d = Gamma(3.0, 2.1)
    for r in (0.0, 0.02, 0.05, 0.1)
        @test onset_rescale(d, r) ≈ mgf(d, -r)
        ## Onsets lag infections, so the rescale shrinks the trajectory.
        @test 0 < onset_rescale(d, r) <= 1
    end
end

@testitem "cases_only_model exposes infections then cases" tags=[:slow] begin
    using Turing: sample, Prior
    import FlexiChains
    using Distributions: Gamma, mgf
    using BVDOutbreakSize: cases_only_model

    chn = sample(cases_only_model(missing), Prior(), 200;
        chain_type = FlexiChains.VNChain, progress = false)

    infections = vec(Array(chn[:cumulative_infections]))
    cases = vec(Array(chn[:cumulative_cases]))
    r = vec(Array(chn[:r]))
    α = vec(Array(chn[:α_inc]))
    θ = vec(Array(chn[:θ_inc]))

    @test all(isfinite, infections) && all(>(0), infections)
    @test all(isfinite, cases) && all(>(0), cases)
    ## Symptomatic cases are the infections rescaled by the incubation
    ## mgf, so they never exceed infections.
    @test all(cases .<= infections .* (1 + 1e-8))
    expected = infections .* mgf.(Gamma.(α, θ), -r)
    @test all(@. isapprox(cases, expected; rtol = 1e-6))
end

@testitem "onset-driven composers expose infections and cases" tags=[:slow] begin
    using Turing: sample, Prior
    import FlexiChains
    using Distributions: Gamma, mgf
    using BVDOutbreakSize: deaths_only_model, exports_deaths_only_model,
                           confirmed_only_model

    ## Each onset-driven single-stream composer should expose both the
    ## latent infections and the incubation-rescaled cases, with
    ## cases = infections · mgf(incubation, −r) per draw.
    for mdl in (deaths_only_model(missing),
        exports_deaths_only_model(Int[]),
        confirmed_only_model(missing))
        chn = sample(mdl, Prior(), 100;
            chain_type = FlexiChains.VNChain, progress = false)
        infections = vec(Array(chn[:cumulative_infections]))
        cases = vec(Array(chn[:cumulative_cases]))
        os = vec(Array(chn[:onset_fraction]))
        r = vec(Array(chn[:r]))
        α = vec(Array(chn[:α_inc]))
        θ = vec(Array(chn[:θ_inc]))

        @test all(isfinite, infections) && all(>(0), infections)
        @test all(0 .< os .<= 1)
        @test all(@. isapprox(os, mgf(Gamma(α, θ), -r); rtol = 1e-6))
        @test all(@. isapprox(cases, infections * os; rtol = 1e-6))
    end
end

@testitem "imperial_only_model has cases equal to infections" tags=[:slow] begin
    using Turing: sample, Prior
    import FlexiChains
    using BVDOutbreakSize: imperial_only_model

    ## No incubation layer: cases and infections coincide (onset_fraction 1).
    chn = sample(imperial_only_model(2, 136), Prior(), 100;
        chain_type = FlexiChains.VNChain, progress = false)
    infections = vec(Array(chn[:cumulative_infections]))
    cases = vec(Array(chn[:cumulative_cases]))

    @test all(isfinite, infections) && all(>(0), infections)
    @test infections == cases
end
