## Tests for the summary dashboard's health-zone table (zone_headline).

@testitem "zone_headline tabulates the zones with most recent cases" begin
    using BVDOutbreakSize
    using DataFrames: DataFrame
    using Markdown

    ## Four zones in the `zone_estimates` frame's columns: Aru has no case
    ## over the past two weeks and Mambasa is below the reproduction-number
    ## floor.
    est = DataFrame(
        label = ["Bunia", "Aru", "Beni", "Mambasa"],
        patch = ["Ituri", "Ituri", "Nord-Kivu", "Ituri"],
        cases = [120, 40, 90, 3],
        cases_last_14 = [12, 0, 30, 2],
        R_T_lower = [0.61, 0.2, 1.05, missing],
        R_T_upper = [1.44, 0.9, 1.87, missing],
        p_rt_above_one = [0.4, 0.01, 0.998, missing],
        forecast_lower = [2.0, 0.0, 10.4, 0.0],
        forecast_upper = [15.0, 1.0, 41.6, 3.0],
        share_lower = [0.301, 0.1, 0.52, 0.0],
        share_upper = [0.452, 0.2, 0.71, 0.02],
        walking = [1, 1, 1, 0]
    )
    md = zone_headline(est; top = 10)
    rows = filter(startswith("| "), split(md, "\n"))
    ## A header, a separator and one row per zone with a recent case,
    ## most recent cases first.
    @test length(rows) == 5
    @test occursin("| Zone | Province | Cases, past two weeks |", rows[1])
    @test occursin("Share of province infections (%)", rows[1])
    @test occursin("R at the cut-off | P(R > 1)", rows[1])
    @test occursin("R modelled separately", rows[1])
    @test startswith(rows[3], "| Beni | Nord-Kivu | 30 | 90 | 52–71 |")
    @test occursin("| 1.05–1.87 | over 99% | 10–42 | yes |", rows[3])
    @test startswith(rows[4], "| Bunia | Ituri | 12 | 120 | 30–45 |")
    @test occursin("| 0.61–1.44 | 40% |", rows[4])
    ## A zone below the floor carries no reproduction number.
    @test startswith(rows[5], "| Mambasa | Ituri | 2 | 3 |")
    @test occursin("| - | - | 0–3 | no |", rows[5])
    @test !occursin("Aru", md)
    @test occursin("3 of 4 zones", md)
    ## The summary page parses the markdown, so it must parse as a table
    ## and a paragraph.
    content = Markdown.parse(md).content
    @test content[1] isa Markdown.Table
    @test content[2] isa Markdown.Paragraph

    ## `top` caps the rows.
    one = zone_headline(est; top = 1)
    @test length(filter(startswith("| "), split(one, "\n"))) == 3

    ## With no recent case anywhere there is a sentence and no table.
    quiet = zone_headline(
        DataFrame(est[est.label .== "Aru", :]); top = 10
    )
    @test !occursin("|", quiet)
    @test occursin("No zone", quiet)
end
