## Tests for streams_table: one row per stream with the documented
## column order.

@testitem "streams_table returns expected columns and row count" begin
    using DataFrames: DataFrame, nrow
    using Random: MersenneTwister
    using BVDOutbreakSize: streams_table
    rng = MersenneTwister(1)
    a = randn(rng, 500) .* 50 .+ 400
    b = randn(rng, 500) .* 80 .+ 600

    df = streams_table("fit A" => a, "fit B" => b)

    @test df isa DataFrame
    @test names(df) ==
        [
        "Stream", "Lower 90%", "Lower 60%", "Lower 30%",
        "Upper 30%", "Upper 60%", "Upper 90%",
    ]
    @test nrow(df) == 2
    @test df[!, "Stream"] == ["fit A", "fit B"]

    for r in eachrow(df)
        @test r["Lower 90%"] <= r["Lower 60%"] <= r["Lower 30%"] <=
            r["Upper 30%"] <= r["Upper 60%"] <= r["Upper 90%"]
    end
end
