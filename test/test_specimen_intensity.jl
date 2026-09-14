## Specimens analysed per suspect sampled. `τ_test` is a probability and the
## receipt kernel conserves mass, so it alone bounds the modelled analysed
## volume below the modelled suspect inflow.

@testitem "specimen intensity lifts the analysed-volume ceiling" begin
    ## Without a factor above one the volume cannot exceed the carried inflow,
    ## whatever `τ_test` is.
    using BVDOutbreakSize: convolve_delay, lab_delay_model
    using Turing.DynamicPPL: fix

    pmf = fix(lab_delay_model(), (; delay_mean = 3.0, delay_sd = 1.5))().pmf
    @test sum(pmf) ≈ 1.0
    carried = convolve_delay(fill(100.0, 60), pmf)
    τ_test = 0.95

    @test sum(τ_test .* carried) < sum(carried)
    @test sum(1.6 .* τ_test .* carried) > sum(carried)
end

@testitem "specimen_intensity_model is centred on no effect" begin
    using BVDOutbreakSize: specimen_intensity_model
    using Turing: returned
    using Random: MersenneTwister
    using Statistics: median

    m = specimen_intensity_model()
    κ = [returned(m, rand(MersenneTwister(s), m)).κ for s in 1:400]
    @test all(>(0), κ)
    ## `LogNormal(0, 0.25)`: median 1, 90% range about 0.66 to 1.51.
    @test 0.9 < median(κ) < 1.1
    @test 0.5 < minimum(κ) && maximum(κ) < 2.5
end

@testitem "a unit intensity reproduces the unmultiplied volume exactly" begin
    ## The nesting claim, tested rather than asserted: at `κ = 1` the analysed
    ## volume must match the `nothing` path element-wise.
    using BVDOutbreakSize: confirmed_cases_model, specimen_intensity_model
    using Turing: returned
    using Turing.DynamicPPL: fix
    using Random: MersenneTwister

    onsets = [10.0 + 0.5i for i in 1:80]
    hist = (; days = Int[60, 70], counts = Int[120, 180])
    lab = (; days = Int[60, 70], counts = Int[400, 500])

    build(si) = confirmed_cases_model(hist, 180, onsets, 0.1, 0.6,
        fill(4.0, 80), 0.8, copy(onsets); lab_history = lab,
        specimen_intensity = si)

    unit = fix(specimen_intensity_model(), (; κ = 1.0))
    off = returned(build(nothing), rand(MersenneTwister(3), build(nothing)))
    on = returned(build(unit), rand(MersenneTwister(3), build(unit)))

    @test off.analysed_daily ≈ on.analysed_daily
    @test off.κ_test === nothing
    @test on.κ_test ≈ 1.0
end

@testitem "a doubled intensity doubles the analysed volume" begin
    ## The factor enters the volume and nothing else: doubling `κ` doubles
    ## `analysed_daily` exactly, leaving the suspect pipeline untouched.
    using BVDOutbreakSize: confirmed_cases_model, specimen_intensity_model
    using Turing: returned
    using Turing.DynamicPPL: fix
    using Random: MersenneTwister

    onsets = [10.0 + 0.5i for i in 1:80]
    hist = (; days = Int[60, 70], counts = Int[120, 180])
    lab = (; days = Int[60, 70], counts = Int[400, 500])

    build(κ) = confirmed_cases_model(hist, 180, onsets, 0.1, 0.6,
        fill(4.0, 80), 0.8, copy(onsets); lab_history = lab,
        specimen_intensity = fix(specimen_intensity_model(), (; κ = κ)))

    one_x = returned(build(1.0), rand(MersenneTwister(5), build(1.0)))
    two_x = returned(build(2.0), rand(MersenneTwister(5), build(2.0)))
    @test two_x.analysed_daily ≈ 2 .* one_x.analysed_daily
    @test two_x.bg_daily ≈ one_x.bg_daily
end
