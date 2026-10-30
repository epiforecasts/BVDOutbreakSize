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

@testitem "care parser reads the count before its split" begin
    using BVDOutbreakSize: BVDOutbreakSize
    include(
        joinpath(
            pkgdir(BVDOutbreakSize), "scripts",
            "scan_province_care.jl"
        )
    )

    ## Folded text of SitRep 150 and 151 paragraphs, shortened.
    ## SitRep 150, Ituri: no verb after the count.
    f = parse_prose_entry(
        "76 nouvelles admissions ont ete enregistrees contre 71 sorties, " *
            "dont 16 guerisons. ainsi, 337 patients, soit un taux " *
            "d'occupation de 33,2 %, dont 41 malades critiques",
        "ituri"
    )
    @test f[:patients_isolated] == 337
    @test !haskey(f, :beds)

    ## SitRep 150, Nord-Kivu: an unopened bracket before the verb, then
    ## the normed-structure count and its beds.
    f = parse_prose_entry(
        "79 nouvelles admissions ont ete enregistrees contre 109 sorties, " *
            "dont 10 guerisons. ainsi, 471 patients dont 67 enfants de " *
            "moins de 5 ans) sont hospitalises : 330 dans les cte/ct pour " *
            "374 lits, soit un taux d'occupation de 88,2 %",
        "nord_kivu"
    )
    @test f[:patients_isolated] == 471
    @test f[:patients_isolated_normed] == 330
    @test f[:beds] == 374

    ## SitRep 151: the split in brackets before the verb.
    f = parse_prose_entry(
        "6 nouvelles admissions ont ete enregistrees. ainsi, 45 patients " *
            "(17 confirmes et 28 suspects) sont hospitalises pour 169 lits " *
            "montes, soit un taux d'occupation de 26,6 %",
        "haut_uele"
    )
    @test f[:patients_isolated] == 45
    @test f[:beds] == 169

    f = parse_prose_entry(
        "11 nouveaux suspects ont ete admis contre 3 sorties de non-cas. " *
            "ainsi, 14 patients (3 confirmes et 11 suspects) restent en " *
            "isolement",
        "bas_uele"
    )
    @test f[:patients_isolated] == 14

    ## The wording every vintage from 109 prints still reads as before.
    f = parse_prose_entry(
        "89 nouvelles admissions contre 107 sorties dont 18 gueries. " *
            "ainsi, 442 patients sont hospitalises pour 1015 lits, soit un " *
            "taux d'occupation de 43,5 %",
        "ituri"
    )
    @test f[:patients_isolated] == 442
    @test f[:beds] == 1015
end

@testitem "care entries continue across a line opening on a bracket" begin
    using BVDOutbreakSize: BVDOutbreakSize
    include(
        joinpath(
            pkgdir(BVDOutbreakSize), "scripts",
            "scan_province_care.jl"
        )
    )

    ## SitRep 153's Haut-Uele paragraph wraps onto a line that opens on a
    ## bracket, and the count and beds sit after it.
    lines = [
        "    •    Au Haut-Uélé, 16 nouvelles admissions (7 confirmés, dont 5 à Wamba, et 9 suspects, dont 1 PPL) ont",
        "         été enregistrées contre 18 sorties : 5 guérisons (4 à Wamba et 1 à Isiro), 10 non-cas et 3 décès de confirmés",
        "         (Isiro, Boma Mangbetu et Rungu). Ainsi, 43 patients (16 confirmés et 27 suspects, dont 5 enfants de moins",
        "         de 5 ans) sont hospitalisés pour 169 lits montés, soit un taux d’occupation de 25,4 %",
        "    •    À la Tshopo, 4 nouvelles admissions ont été enregistrées.",
    ]
    got, order = province_entries(lines)
    @test order == ["haut_uele", "tshopo"]
    f = parse_prose_entry(got["haut_uele"], "haut_uele")
    @test f[:patients_isolated] == 43
    @test f[:beds] == 169

    ## A quote cut at 300 bytes keeps a character that straddles the cut.
    q = repeat("a", 298) * "’" * "occupation"
    @test trimmed_quote(q) == repeat("a", 298) * "’"
end

@testitem "care parser reads beds from the CTE/CT rate fraction" begin
    using BVDOutbreakSize: BVDOutbreakSize
    include(
        joinpath(
            pkgdir(BVDOutbreakSize), "scripts",
            "scan_province_care.jl"
        )
    )

    ## Folded text of SitRep 156 to 159 paragraphs, shortened.
    ## SitRep 156, Nord-Kivu: the beds are the denominator of the CTE/CT
    ## rate, with no "lits", and the space before "de" is missing.
    f = parse_prose_entry(
        "50 nouvelles admissions ont ete enregistrees contre 62 sorties : " *
            "1 guerison (masereka). ainsi, 471 patients sont hospitalises, " *
            "dont 2 malades critiques, pour un taux d'occupation des " *
            "cte/ctde 85,1 % : (349/410), de 92,7 % pour les suspects " *
            "(178/192) et de 78,4 % pour les confirmes (171/218)",
        "nord_kivu"
    )
    @test f[:patients_isolated] == 471
    @test f[:patients_isolated_normed] == 349
    @test f[:beds] == 410

    ## SitRep 159, Nord-Kivu: "global", a curly apostrophe and an aside.
    f = parse_prose_entry(
        "87 nouvelles admissions ont ete enregistrees contre 106 sorties, " *
            "dont 6 gueris. ainsi, 488 patients sont hospitalises parmi " *
            "eux 59 enfants de moins de 5 ans et 3 malades critiques au " *
            "cte kitatumba et 148 hors cte/ct, pour un taux d’occupation " *
            "global des cte/ct de 82,9 % (340/410)",
        "nord_kivu"
    )
    @test f[:patients_isolated] == 488
    @test f[:patients_isolated_normed] == 340
    @test f[:beds] == 410

    ## SitRep 159, Bas-Uele: the count follows "portant a", and the three
    ## confirmed beds are not read as the province's beds.
    f = parse_prose_entry(
        "14 nouvelles admissions de suspects ont ete enregistrees et " *
            "isolees dans les ess contre 1 sortie (non-cas a buta), " *
            "portant a 60 le nombre de patients en isolement ; les 3 lits " *
            "pour confirmes sont occupes (100,0 %)",
        "bas_uele"
    )
    @test f[:patients_isolated] == 60
    @test !haskey(f, :beds)
end

@testitem "care parser reads patients pris en charge CTE/CT" begin
    using BVDOutbreakSize: BVDOutbreakSize
    include(
        joinpath(
            pkgdir(BVDOutbreakSize), "scripts",
            "scan_province_care.jl"
        )
    )

    ## SitRep 164, Nord-Kivu: the count with "pris en charge cte/ct" and a
    ## rate with no fraction or beds.
    f = parse_prose_entry(
        "495 nouvelles admissions ont ete enregistrees contre 82 sorties, " *
            "dont 4 patients gueris. ainsi, 497 patients sont pris en charge " *
            "cte/ct, pour un taux d’occupation global des cte/ct de 87,8 %.",
        "nord_kivu"
    )
    @test f[:patients_isolated] == 497
    @test !haskey(f, :beds)
end
