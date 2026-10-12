## Tests for the care-section locator in `scripts/scan_province_care.jl`.
##
## The section is found by its title and ends at the next heading outside
## its own number. From SitRep 146 the report nests it one level too deep
## ("2.4.3. Prise en charge holistique" under PCI/EDS) and then numbers
## its subsections as siblings ("2.4.4. Principales actions", "2.4.5.
## Défis"), so a locator that only follows the title's own number stops
## at the first subheading and reads no province at all.
##
## Filter target for a scoped run (this worktree only):
##   root = joinpath(pwd(), "test")
##   TestItemRunner.run_tests(root; filter = ti -> ti.filename ==
##       joinpath(root, "test_province_care_parser.jl"))

@testitem "care section locator reads sibling-numbered subsections" begin
    using BVDOutbreakSize: BVDOutbreakSize
    include(
        joinpath(
            pkgdir(BVDOutbreakSize), "scripts",
            "scan_province_care.jl"
        )
    )

    ## SitRep 146's layout, shortened.
    sibling = [
        "2.4.2. Défis",
        "   • À la Tshopo, les intrants PCI restent insuffisants.",
        "2.4.3. Prise en charge holistique",
        "2.4.4. Principales actions",
        "    • En Ituri, 45 nouvelles admissions ont été enregistrées " *
            "contre 76 sorties, dont 15 guérisons. 360 patients",
        "       demeuraient en isolement à la fin de la journée.",
        "    •   À la Tshopo, 14 patients demeuraient en isolement.",
        "2.4.5. Défis",
        "    • En Ituri, la faible capacité d'accueil demeure.",
        "2.5. Vaccination",
    ]
    got = care_lines(sibling)
    @test any(l -> occursin("En Ituri, 45 nouvelles", l), got)
    @test any(l -> occursin("Tshopo, 14 patients", l), got)
    @test !any(l -> occursin("capacité d'accueil", l), got)
    @test !any(l -> occursin("Vaccination", l), got)

    ## The layout every vintage from 109 to 145 prints, which must read as
    ## before: subsections numbered under the section, ending at Défis.
    nested = [
        "2.5. Prise en charge holistique",
        "2.5.1. Principales actions",
        "    • En Ituri, 86 nouvelles admissions ont été enregistrées.",
        "2.5.2. Défis",
        "    • En Ituri, la faible capacité d'accueil demeure.",
        "2.6. Vaccination",
    ]
    got = care_lines(nested)
    @test any(l -> occursin("En Ituri, 86 nouvelles", l), got)
    @test !any(l -> occursin("capacité d'accueil", l), got)

    ## A heading with any other title still ends the section.
    other = [
        "2.4.3. Prise en charge holistique",
        "2.4.4. Principales actions",
        "    • En Ituri, 45 nouvelles admissions.",
        "2.4.5. Vaccination",
        "    • Au Bas-Uélé, 1 412 personnes vaccinées.",
    ]
    got = care_lines(other)
    @test any(l -> occursin("En Ituri, 45", l), got)
    @test !any(l -> occursin("vaccinées", l), got)
end
