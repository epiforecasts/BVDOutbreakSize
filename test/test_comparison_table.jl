## Tests for comparison_table: one row per scenario with a
## narrowest_CrI value drawn from the documented label set.

@testitem "comparison_table covers every REPORT_SCENARIOS entry" begin
    using DataFrames: DataFrame, nrow
    using Random: MersenneTwister
    using BVDOutbreakSize: comparison_table, REPORT_SCENARIOS

    _CRI_LABELS = Set(["30%", "60%", "90%", "outside 90%"])

    rng = MersenneTwister(2)
    # Draws span the published scenario range so every label can in
    # principle be assigned, but we only check membership.
    C = randn(rng, 2_000) .* 200 .+ 500

    df = comparison_table(C)
    @test df isa DataFrame
    @test names(df) ==
          ["Scenario", "Reported cases", "Narrowest interval"]
    @test nrow(df) == length(REPORT_SCENARIOS)
    @test df[!, "Scenario"] == [label for (label, _) in REPORT_SCENARIOS]
    @test df[!, "Reported cases"] == [val for (_, val) in REPORT_SCENARIOS]

    for v in df[!, "Narrowest interval"]
        @test v in _CRI_LABELS
    end
end

@testitem "comparison_table handles a far-away posterior" begin
    using Random: MersenneTwister
    using BVDOutbreakSize: comparison_table
    # All draws cluster well above every reported scenario, so the
    # narrowest label is always "outside 90%".
    rng = MersenneTwister(3)
    C = randn(rng, 1_000) .* 10 .+ 10_000
    df = comparison_table(C)
    @test all(==("outside 90%"), df[!, "Narrowest interval"])
end
