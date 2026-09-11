## Tests for markdown_table and MarkdownTable: a DataFrame rendered as a
## GitHub-flavoured markdown table, with integer-valued floats printed
## without a trailing `.0` so count columns read as whole numbers, and a
## display wrapper that a Literate page shows as markdown rather than as
## html.

@testitem "markdown_table renders a DataFrame as markdown" begin
    using DataFrames: DataFrame
    using BVDOutbreakSize: markdown_table

    df = DataFrame("Quantity" => ["Cumulative infections", "Case-fatality"],
        "Lower 90%" => [1234.0, 0.12],
        "Upper 90%" => [5678.0, 0.34])
    md = markdown_table(df)
    lines = split(strip(md), "\n")

    @test lines[1] == "| Quantity | Lower 90% | Upper 90% |"
    ## Numeric columns are right-aligned, label columns left.
    @test lines[2] == "| --- | ---: | ---: |"
    ## Integer-valued counts drop the `.0`; genuine decimals are kept.
    @test lines[3] == "| Cumulative infections | 1234 | 5678 |"
    @test lines[4] == "| Case-fatality | 0.12 | 0.34 |"
    ## Header, separator and one row per data row.
    @test length(lines) == 4
end

@testitem "markdown_table keeps a table readable" begin
    using DataFrames: DataFrame
    using BVDOutbreakSize: markdown_table

    ## A raw Float64 would print its full expansion and stretch the column.
    ## At or above one the magnitude is kept and the tail dropped; below one
    ## the leading digits are kept so a small score does not become zero.
    df = DataFrame("crps" => [231.13920000000002, 0.00012345678])
    lines = split(strip(markdown_table(df)), "\n")
    @test lines[3] == "| 231.139 |"
    @test lines[4] == "| 0.000123 |"

    ## A `|` inside a cell is escaped so it cannot split the row.
    piped = DataFrame("Stream" => ["cases | deaths"])
    @test split(strip(markdown_table(piped)), "\n")[3] ==
          "| cases \\| deaths |"

    ## `Bool` reads as a label rather than a quantity, so it stays left.
    @test split(strip(markdown_table(DataFrame("ok" => [true]))), "\n")[2] ==
          "| --- |"
end

@testitem "MarkdownTable shows as markdown and not as html" begin
    using DataFrames: DataFrame
    using BVDOutbreakSize: MarkdownTable, markdown_table
    using Markdown: Markdown

    df = DataFrame("Stream" => ["confirmed cases"], "crps" => [1.5])
    t = MarkdownTable(df)

    ## Literate takes `text/html` before `text/markdown`, so the wrapper
    ## must offer markdown and refuse html for the table to render as one.
    @test showable(MIME("text/markdown"), t)
    @test !showable(MIME("text/html"), t)
    @test sprint(show, MIME("text/markdown"), t) == markdown_table(df)
    @test sprint(print, t) == markdown_table(df)

    ## A page may pass a placeholder message instead of a table when an
    ## optional analysis is not run in a given build.
    msg = MarkdownTable(Markdown.md"_Not shown in this build._")
    @test showable(MIME("text/markdown"), msg)
    @test occursin("Not shown in this build.",
        sprint(show, MIME("text/markdown"), msg))
end
