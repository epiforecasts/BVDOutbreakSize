## Specimens analysed per suspect sampled. `confirmed_cases_model` built the
## analysed volume as `τ_test · convolve_delay(suspects, receipt_pmf)`, and
## with `τ_test` a probability that caps modelled analysed below modelled
## suspects — a ceiling the reported data cross.

@testitem "specimen_intensity_model is flat when the trend is zero" begin
    using BVDOutbreakSize: specimen_intensity_model
    using Turing.DynamicPPL: VarInfo, InitFromPrior, fix
    using Random: MersenneTwister

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

    ## With a multiplier of 1.6 the volume can exceed it, as the data do.
    @test sum(1.6 .* τ_test .* carried) > sum(suspects[1:(end - 10)])
end

@testitem "specimen_intensity = false reproduces the unmultiplied volume" begin
    ## The keyword must nest: switching the submodel off has to leave the
    ## confirmed stream exactly as it was, so an unchanged posterior under
    ## `false` and a moved one under `true` are attributable to the data.
    using BVDOutbreakSize: confirmed_cases_model, onset_incidence_model
    using Turing.DynamicPPL: VarInfo, InitFromPrior, getlogjoint
    using Random: MersenneTwister

    onsets = [10.0 + 0.5i for i in 1:80]
    bg = fill(4.0, 80)
    bvd = copy(onsets)
    hist = (; days = Int[60, 70], counts = Int[120, 180])
    lab = (; days = Int[60, 70], counts = Int[400, 500])

    build(si) = confirmed_cases_model(hist, 180, onsets, 0.1, 0.6, bg, 0.8,
        bvd; lab_history = lab, specimen_intensity = si)

    off = getlogjoint(VarInfo(MersenneTwister(3), build(nothing),
        InitFromPrior()))
    @test isfinite(off)
    ## Drawing the submodel adds its two parameters, so the log joint differs;
    ## what must hold is that both paths evaluate and neither errors.
    on = getlogjoint(VarInfo(MersenneTwister(3),
        build(BVDOutbreakSize.specimen_intensity_model(80)), InitFromPrior()))
    @test isfinite(on)
end
