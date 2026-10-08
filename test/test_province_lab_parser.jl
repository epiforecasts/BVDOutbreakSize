## Tests for the per-province laboratory parser
## (`scripts/scan_province_lab.jl`).
##
## The parser reads free prose that INSP rewords between vintages, and a
## wording it cannot read drops that whole vintage. That failed quietly
## once: from SitRep 119 the brief format began writing a wholly-negative
## batch as "reveles" rather than "revenus", and six vintages (119-124)
## went missing from `province_lab_daily_history` before anyone noticed
## (issue #695). Every phrasing below is taken verbatim from a report, so a
## future rewording that breaks one of them fails here rather than
## silently shortening the series.
##
## Filter target for a scoped run (this worktree only, not every sibling
## worktree `@run_package_tests` would otherwise discover):
##   target = joinpath(pwd(), "test", "test_province_lab_parser.jl")
##   @run_package_tests filter = ti -> string(ti.filename) == target

@testitem "province lab parser reads every phrasing in the archive" begin
    using BVDOutbreakSize: BVDOutbreakSize
    include(
        joinpath(
            pkgdir(BVDOutbreakSize), "scripts",
            "scan_province_lab.jl"
        )
    )

    ## `(entry, province, expected)`, each entry the text a bullet leaves
    ## after `province_entries` strips the province name and colon, folded
    ## to the unaccented lower case the parser works in.
    cases = [
        ## The two ordinary shapes, which have always parsed.
        (
            "57 nouveaux resultats positifs (40 vivants et 17 deces) sur 203 " *
                "nouveaux echantillons recus et analyses (positivite : 28,1%).",
            "ituri", (203, 57),
        ),
        (
            "21 nouveaux resultats positifs (13 vivants et 8 deces) sur 151 " *
                "echantillons analyses (positivite de 13,9%).",
            "nord_kivu", (151, 21),
        ),

        ## Wholly negative. "revenus" was the only spelling the parser
        ## read; "reveles" is what SitReps 119-124 actually print.
        (
            "3 echantillons recus et testes (3 vivants), tous sont revenus " *
                "negatifs.",
            "tshopo", (3, 0),
        ),
        (
            "3 echantillons recus et testes (3 vivants), tous se sont " *
                "reveles negatifs.",
            "bas_uele",
            (3, 0),
        ),
        (
            "3 echantillons recus et testes (2 vivants et 1 deces), tous se " *
                "sont reveles negatifs.",
            "bas_uele",
            (3, 0),
        ),
        (
            "13 echantillons recus et analyses, tous reveles negatifs.",
            "haut_uele", (13, 0),
        ),

        ## Wholly positive: the numerator is the whole denominator.
        (
            "3 echantillons recus et testes (vivants) tous se sont reveles " *
                "positifs.",
            "tshopo", (3, 3),
        ),
        (
            "1 echantillon recu et teste (1 vivant), s'est revele positif.",
            "bas_uele", (1, 1),
        ),

        ## A definite singular denominator, printed without a digit
        ## (SitRep 119, Sud Ubangi). This one used to return (0, 0), which
        ## broke the national partition without failing loudly.
        (
            "1 nouveau resultat positif (1 deces) sur l'echantillon analyse " *
                "(positivite : 100,0%)",
            "sud_ubangi",
            (1, 1),
        ),

        ## SitRep 124's Nord-Kivu bullet leads with its positive count.
        ## Admitted only because the printed positivity reproduces the
        ## ratio: 21/160 = 13,1%.
        (
            "21 echantillons recus et testes (14 vivants et 7 deces) sur 160 " *
                "echantillons analyses (positivite de 13,1%).",
            "nord_kivu", (160, 21),
        ),

        ## SitRep 130's Tshopo bullet drops "resultats positifs" and names
        ## no samples before its parenthetical. Admitted by the printed
        ## positivity in the same way: 3/4 = 75%.
        (
            "3 nouveaux (2 vivants et 1 deces) sur 4 echantillons recus et " *
                "analyses au laboratoire d'isiro (positivite : 75%).",
            "tshopo", (4, 3),
        ),

        ## No completed analysis: no denominator, so no numerator either.
        (
            "3 echantillons collectes dont 2 sont expedies a kinshasa (inrb).",
            "sud_ubangi", (0, 0),
        ),
    ]

    for (entry, province, expected) in cases
        @test parse_province_entry(entry, province) == expected
    end
end

@testitem "province lab parser refuses to guess a leading count" begin
    using BVDOutbreakSize: BVDOutbreakSize
    include(
        joinpath(
            pkgdir(BVDOutbreakSize), "scripts",
            "scan_province_lab.jl"
        )
    )

    ## The leading-count reading is admitted by the printed positivity
    ## alone. Strip the rate and the bullet must go back to unparsed, so
    ## the vintage is reported rather than read the wrong way round.
    norate = "21 echantillons recus et testes (14 vivants et 7 deces) sur " *
        "160 echantillons analyses."
    @test parse_province_entry(norate, "nord_kivu") == :unparsed

    ## A rate that does not match the ratio must not admit it either:
    ## 21/160 is 13,1%, not 44,0%.
    wrong = "21 echantillons recus et testes (14 vivants et 7 deces) sur " *
        "160 echantillons analyses (positivite de 44,0%)."
    @test parse_province_entry(wrong, "nord_kivu") == :unparsed

    ## The same holds for SitRep 130's "N nouveaux (" lead.
    norate_nouveaux = "3 nouveaux (2 vivants et 1 deces) sur 4 " *
        "echantillons recus et analyses au laboratoire d'isiro."
    @test parse_province_entry(norate_nouveaux, "tshopo") == :unparsed
end

@testitem "committed province lab block partitions the national totals" begin
    using BVDOutbreakSize: BVDOutbreakSize
    using TOML: parsefile

    obs = parsefile(
        joinpath(
            pkgdir(BVDOutbreakSize), "data",
            "observations.toml"
        )
    )
    lab = obs["province_lab_daily_history"]
    nat = obs["tests_analysed_daily_history"]
    national = Dict(zip(nat["dates"], nat["values"]))
    provinces = (
        "ituri", "nord_kivu", "sud_kivu", "haut_uele", "tshopo",
        "bas_uele", "sud_ubangi",
    )

    ## The per-province analysed counts are an exact partition of the
    ## national series, which is what admits a vintage into this block.
    for (i, date) in enumerate(lab["dates"])
        @test haskey(national, date)
        @test sum(lab[p * "_analysed"][i] for p in provinces) ==
            national[date]
        ## A positive needs a sample to have come from.
        for p in provinces
            @test lab[p * "_positive"][i] <= lab[p * "_analysed"][i]
        end
    end

    ## The window closed by issue #695: every date from SitRep 119 onward
    ## is present, not just the ones whose wording happened to parse.
    for date in (
            "2026-09-10", "2026-09-11", "2026-09-12", "2026-09-13",
            "2026-09-14", "2026-09-15",
        )
        @test date in lab["dates"]
    end
end
