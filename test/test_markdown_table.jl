## Tests for markdown_table: renders a DataFrame as a GitHub-flavoured
## markdown table, with integer-valued floats printed without a trailing
## `.0` so count columns read as whole numbers.

@testitem "markdown_table renders a DataFrame as markdown" begin
    using DataFrames: DataFrame
    using BVDOutbreakSize: markdown_table

    df = DataFrame("Quantity" => ["Cumulative infections", "Case-fatality"],
        "Lower 90%" => [1234.0, 0.12],
        "Upper 90%" => [5678.0, 0.34])
    md = markdown_table(df)
    lines = split(strip(md), "\n")

    @test lines[1] == "| Quantity | Lower 90% | Upper 90% |"
    @test lines[2] == "| --- | --- | --- |"
    ## Integer-valued counts drop the `.0`; genuine decimals are kept.
    @test lines[3] == "| Cumulative infections | 1234 | 5678 |"
    @test lines[4] == "| Case-fatality | 0.12 | 0.34 |"
    ## Header, separator and one row per data row.
    @test length(lines) == 4
end
