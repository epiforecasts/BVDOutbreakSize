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
