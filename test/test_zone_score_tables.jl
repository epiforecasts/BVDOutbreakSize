## Tests for reading the per-zone release scores back on the report pages
## (zone_score_key, zone_score_rows).

@testitem "zone_score_key reads the zone key off a scored label" begin
    using BVDOutbreakSize

    @test zone_score_key("confirmed cases {ituri.bunia}") == "ituri.bunia"
    @test isnothing(zone_score_key("confirmed cases [ituri]"))
    @test isnothing(zone_score_key("confirmed cases"))
end

@testitem "zone_score_rows keeps the named zones in order, relabelled" begin
    using BVDOutbreakSize
    using DataFrames: DataFrame, nrow, names

    tbl = DataFrame(
        stream = [
            "confirmed cases {ituri.bunia}", "confirmed cases {ituri.aru}",
            "confirmed cases {nord_kivu.beni}",
            "confirmed cases {ituri.bunia}",
        ],
        fit = ["joint", "joint", "joint", "baseline"],
        crps = [1.0, 2.0, 3.0, 4.0]
    )
    out = zone_score_rows(
        tbl, ["nord_kivu.beni", "ituri.bunia"], ["Beni", "Bunia"]
    )
    ## Rows follow the order of the keys, and a zone not named is dropped.
    @test out.stream == ["Beni", "Bunia", "Bunia"]
    @test out.crps == [3.0, 1.0, 4.0]
    @test nrow(tbl) == 4

    ## An empty table keeps its columns.
    none = zone_score_rows(tbl[1:0, :], ["ituri.bunia"], ["Bunia"])
    @test nrow(none) == 0
    @test names(none) == names(tbl)

    ## Keys and labels must pair up.
    @test_throws ArgumentError zone_score_rows(
        tbl, ["ituri.bunia", "ituri.aru"], ["Bunia"]
    )
end
