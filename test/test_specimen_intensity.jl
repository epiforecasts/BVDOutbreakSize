## Specimens analysed per suspect sampled. `confirmed_cases_model` built the
## analysed volume as `τ_test · convolve_delay(suspects, receipt_pmf)`, and
## with `τ_test` a probability that caps modelled analysed below modelled
## suspects — a ceiling the reported data cross.

@testitem "specimen_intensity_model is flat when the trend is zero" begin
    using BVDOutbreakSize: specimen_intensity_model
    using Turing.DynamicPPL: fix

    m = fix(specimen_intensity_model(30; ref_day = 10), (; κ0 = 2.0, β_κ = 0.0))
    st = m()
    @test length(st.κ) == 30
    @test all(≈(2.0), st.κ)
end

@testitem "specimen_intensity_model trends from its reference day" begin
    using BVDOutbreakSize: specimen_intensity_model
    using Turing.DynamicPPL: fix

    ## `κ(t) = κ0 · exp(β_κ (t − ref) / 30)`, so the level is exact at `ref`
    ## and one 30-day step multiplies by `exp(β_κ)`.
    m = fix(specimen_intensity_model(70; ref_day = 40),
        (; κ0 = 1.5, β_κ = log(2.0)))
    st = m()
    @test st.κ[40] ≈ 1.5
    @test st.κ[70] ≈ 1.5 * 2.0
    @test st.κ[10] ≈ 1.5 / 2.0
    ## Monotone increasing for a positive trend.
    @test issorted(st.κ)
end

@testitem "specimen intensity lifts the analysed-volume ceiling" begin
    ## The structural point. With `κ ≡ 1` the modelled analysed volume can
    ## never exceed the modelled suspect inflow, because the receipt PMF sums
    ## to one and `convolve_delay` only loses mass. A multiplier above one
    ## removes that cap, which is what the observed ratio (0.93 June, 1.01
    ## July, 1.50 over 1-5 August) requires.
    using BVDOutbreakSize: convolve_delay, lab_delay_model
    using Turing.DynamicPPL: fix

    pmf = fix(lab_delay_model(), (; delay_mean = 3.0, delay_sd = 1.5))().pmf
    @test sum(pmf) ≈ 1.0
    suspects = fill(100.0, 60)
    carried = convolve_delay(suspects, pmf)
    τ_test = 0.95

    ## Old behaviour: bounded above by the suspect inflow, for any τ_test < 1.
    @test sum(τ_test .* carried) < sum(suspects)

    ## With a multiplier above 1/τ_test the volume exceeds the carried
    ## inflow, which is the structural claim; the data reach 1.50 monthly.
    @test sum(1.6 .* τ_test .* carried) > sum(carried)
end

@testitem "a unit intensity reproduces the unmultiplied volume exactly" begin
    ## The nesting claim, tested rather than asserted: with `κ0 = 1` and
    ## `β_κ = 0` the analysed volume must match the `nothing` path
    ## element-wise, so `specimen_intensity = false` and an unmoved
    ## posterior mean the same thing.
    using BVDOutbreakSize: confirmed_cases_model, specimen_intensity_model
    using Turing: returned
    using Turing.DynamicPPL: fix
    using Random: MersenneTwister

    onsets = [10.0 + 0.5i for i in 1:80]
    bg = fill(4.0, 80)
    hist = (; days = Int[60, 70], counts = Int[120, 180])
    lab = (; days = Int[60, 70], counts = Int[400, 500])

    build(si) = confirmed_cases_model(hist, 180, onsets, 0.1, 0.6, bg, 0.8,
        copy(onsets); lab_history = lab, specimen_intensity = si)

    unit = fix(specimen_intensity_model(80), (; κ0 = 1.0, β_κ = 0.0))
    rng() = MersenneTwister(3)
    off = returned(build(nothing), rand(rng(), build(nothing)))
    on = returned(build(unit), rand(rng(), build(unit)))
    @test off.analysed_daily ≈ on.analysed_daily
    @test off.κ_daily === nothing
    @test all(≈(1.0), on.κ_daily)
end
