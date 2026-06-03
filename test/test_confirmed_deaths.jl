## Tests for the confirmed-case window alignment, the per-window
## positivity random effect and the confirmed-death thinning stream.

@testitem "confirmed_positivity_windows merges stalled/inconsistent windows" begin
    using BVDOutbreakSize: confirmed_positivity_windows

    ## The 28 May data: analysed cumulative stalls 24-25 May (flat at 295)
    ## and the confirmed series is a superset of the analysed dates. Walking
    ## the analysed vintages and merging the zero-denominator window
    ## reproduces the per-window positives [101, 4, 16, 4, 85].
    confirmed = (; days = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11],
        counts = [33, 51, 57, 79, 83, 101, 105, 106, 121, 125, 210])
    lab = (; days = [6, 7, 8, 9, 10, 11],
        counts = [211, 295, 295, 403, 648, 755])
    w = confirmed_positivity_windows(confirmed, lab)
    @test w.positives == [101, 4, 16, 4, 85]
    @test w.analysed == [211, 84, 108, 245, 107]
    @test all(0 .<= w.positives .<= w.analysed)
    @test sum(w.positives) == confirmed.counts[end]
end

@testitem "confirmed_positivity_windows is empty without both histories" begin
    using BVDOutbreakSize: confirmed_positivity_windows
    empty = (; days = Int[], counts = Int[])
    lab = (; days = [3, 5], counts = [10, 20])
    @test isempty(confirmed_positivity_windows(empty, lab).analysed)
    conf = (; days = [3, 5], counts = [2, 6])
    @test isempty(confirmed_positivity_windows(conf, empty).analysed)
end

@testitem "confirmed_positivity_model draws per-window probabilities" begin
    using BVDOutbreakSize: confirmed_positivity_model
    using Random: seed!
    seed!(1)
    s = confirmed_positivity_model(5)()
    @test length(s.p_pos) == 5
    @test all(0 .<= s.p_pos .<= 1)
    @test s.σ_q >= 0
end

@testitem "confirmed_deaths_model gives a bounded confirmation probability" begin
    using BVDOutbreakSize: confirmed_deaths_model
    using Turing: returned
    using Random: MersenneTwister

    ## 246 suspected deaths, 17 confirmed. The composition q_susp is built
    ## from a flat unit BVD report series and a small background, and the
    ## confirmation probability is the odds-enriched composition.
    bvd = fill(1.0, 40)
    m = confirmed_deaths_model(17, 246, 250.0, bvd, 0.3, 0.5)
    st = returned(m, rand(MersenneTwister(2), m))
    @test 0 < st.q_susp < 1
    @test 0 < st.p_death_conf < 1
    @test st.m_death > 0
    @test st.expected_confirmed_deaths >= 0
end

@testitem "confirmed_deaths_only_model conditions and stays finite" begin
    using BVDOutbreakSize: confirmed_deaths_only_model
    using Turing.DynamicPPL: logjoint
    using Random: MersenneTwister
    m = confirmed_deaths_only_model(40, 17, 246;
        deaths_history = (; days = [20, 40], counts = [120, 246]))
    draw = rand(MersenneTwister(1), m)
    @test isfinite(logjoint(m, draw))
end
