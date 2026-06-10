@testitem "lab pipeline daily likelihood runs and is finite" begin
    using BVDOutbreakSize: confirmed_only_model
    using Turing: logjoint
    using Random: MersenneTwister

    ## The confirmed-cases likelihood is exercised end to end through
    ## confirmed_only_model, which draws the shared report kernel, fits the
    ## received-specimen stream and scores the confirmed positives as a
    ## Binomial of the observed analysed denominator.
    m = confirmed_only_model(40, 8;
        confirmed_history = (; days = [20, 40], counts = [3, 8]),
        lab_history = (; days = [20, 40], counts = [5, 9]),
        tests_received_history = (; days = [20, 40], counts = [6, 11]))
    lp = logjoint(m, rand(MersenneTwister(1), m))
    @test isfinite(lp)
end

@testitem "confirmed_positivity_windows anchors dark days on 24h analysed" begin
    using BVDOutbreakSize: confirmed_positivity_windows

    ## Last cumulative laboratory date is day 20, so days 21-23 are dark
    ## late windows. A published 24h analysed count on day 22 flags that
    ## window for the observed-denominator Binomial; the others stay dark.
    confirmed = (; days = [5, 20, 21, 22, 23], counts = [10, 50, 60, 75, 90])
    lab = (; days = [10, 20], counts = [100, 300])
    daily = (; days = [22], counts = [40])
    w = confirmed_positivity_windows(confirmed, lab, daily)
    @test w.late_days == [21, 22, 23]
    @test w.late_increments == [10, 15, 15]
    @test w.late_analysed == [0, 40, 0]

    ## With no daily series every late day is dark (no anchor), and the
    ## field is present and aligned with `late_days`.
    w0 = confirmed_positivity_windows(confirmed, lab)
    @test length(w0.late_analysed) == length(w0.late_days)
    @test all(==(0), w0.late_analysed)

    ## A daily count on or before the last laboratory date is ignored: the
    ## cumulative series already covers it.
    w1 = confirmed_positivity_windows(confirmed, lab,
        (; days = [15, 22], counts = [99, 40]))
    @test w1.late_analysed == [0, 40, 0]
end

@testitem "lab pipeline daily anchor scores dark windows as Binomial" begin
    using BVDOutbreakSize: confirmed_only_model
    using Turing: logjoint
    using Random: MersenneTwister

    ## Post-cutoff dark vintages (days 30, 35, 40 after the last lab date
    ## 20) where day 35 publishes a 24h analysed count: it is scored as a
    ## Binomial of that observed denominator, the others as modelled-volume
    ## dark windows.
    m = confirmed_only_model(40, 20;
        confirmed_history =
        (; days = [20, 30, 35, 40], counts = [5, 9, 14, 20]),
        lab_history = (; days = [10, 20], counts = [12, 28]),
        lab_daily_history = (; days = [35], counts = [30]),
        tests_received_history = (; days = [20, 40], counts = [15, 33]))
    lp = logjoint(m, rand(MersenneTwister(1), m))
    @test isfinite(lp)

    ## When the confirmed increment exceeds the observed denominator it is
    ## clamped into the Binomial support, so the likelihood stays finite.
    m2 = confirmed_only_model(40, 20;
        confirmed_history =
        (; days = [20, 30, 35, 40], counts = [5, 9, 14, 20]),
        lab_history = (; days = [10, 20], counts = [12, 28]),
        lab_daily_history = (; days = [35], counts = [3]),
        tests_received_history = (; days = [20, 40], counts = [15, 33]))
    @test isfinite(logjoint(m2, rand(MersenneTwister(2), m2)))
end
