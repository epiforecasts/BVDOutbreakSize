## Tests for the Tableau 2 health-zone parser
## (`scripts/scan_zone_tableau2.jl`).
##
## Filter target for a scoped run (this worktree only, not every sibling
## worktree `@run_package_tests` would otherwise discover):
##   target = joinpath(pwd(), "test", "test_zone_tableau2_parser.jl")
##   @run_package_tests filter = ti -> string(ti.filename) == target

@testitem "zone parser reads a name set one space from its cases cell" begin
    using BVDOutbreakSize: BVDOutbreakSize
    include(
        joinpath(
            pkgdir(BVDOutbreakSize), "scripts",
            "scan_zone_tableau2.jl"
        )
    )

    caption = "Tableau 2. Répartition des cas et décès confirmés " *
        "par province et zone au 7 novembre 2026"
    ## `sep` sits between each name and its cases cell; the text layer
    ## of a report usually gives two or more spaces and sometimes one.
    table(sep) = join(
        [
            caption,
            "Province / Zone de santé    Nombre cumulatif    Cas (n)    Décès (n)    Létalité",
            "Ituri$(sep)1260    470     37,3%",
            "Bunia$(sep)41      17      41,5%",
            "Source: INSP 2026",
            "Page 3",
            "Rwampara$(sep)1219    444     36,4%",
            "Nord-Kivu$(sep)164    81    49,4%",
            "Musienene$(sep)164     81      49,4%",
            "Sud Ubangi$(sep)2      1       50,0%",
            "Bulu$(sep)2            1       50,0%",
            "Total$(sep)1426    552     38,7%",
        ], "\n"
    )

    wide = scan_tableau2(table("      "))
    narrow = scan_tableau2(table(" "))
    @test wide[4] == ""
    @test narrow[4] == ""
    @test narrow[1] == wide[1]
    @test narrow[2] == wide[2]
    @test narrow[3] == wide[3]
    @test wide[1]["ituri"]["rwampara"] == (1219, 444)
    @test wide[1]["nord_kivu"]["musienene"] == (164, 81)
    @test wide[1]["sud_ubangi"]["bulu"] == (2, 1)
    @test wide[2]["sud_ubangi"] == (2, 1)
    @test wide[3] == (1426, 552)
    @test narrow[1]["ituri"]["bunia"] == (41, 17)
    @test narrow[1]["ituri"]["rwampara"] == (1219, 444)
    @test narrow[2]["ituri"] == (1260, 470)
    @test narrow[2]["sud_ubangi"] == (2, 1)
    @test narrow[3] == (1426, 552)
    ## A footnote or page line ending in a number adds no zone.
    @test sort(collect(keys(narrow[1]["ituri"]))) == ["bunia", "rwampara"]
end

@testitem "no zone or province name in the manifest ends in a digit" begin
    using BVDOutbreakSize: BVDOutbreakSize
    using TOML
    manifest = TOML.parsefile(
        joinpath(pkgdir(BVDOutbreakSize), "data", "observations.toml")
    )
    keys_ = String[]
    for zones in values(manifest["zone_confirmed_history"])
        zones isa AbstractDict && append!(keys_, keys(zones))
    end
    @test !isempty(keys_)
    @test all(k -> !occursin(r"[0-9]$", k), keys_)
end
