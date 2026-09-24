## Tests for the onset-figure audit (`scripts/audit_onset_curve.jl` and the
## file it writes, `data/onset_curve_figures.csv`).
##
## Two layers: consistency of the committed figures CSV with the scanned
## CSV and the printed-n parser, which need nothing but the repository; and
## a read of two real figures, which needs `data/sitrep_pdfs` (git-ignored,
## fetched by `scripts/download_sitreps.jl`) and skips itself when that is
## absent.
##
## The audit uses no package code, so these items locate the repository
## from this file and run without loading BVDOutbreakSize.
##
## Filter target for a scoped run (this worktree only):
##   target = joinpath(pwd(), "test", "test_onset_audit.jl")
##   @run_package_tests filter = ti -> string(ti.filename) == target

@testitem "onset figure audit CSV matches the scanned CSV" default_imports = false begin
    using Test
    using Dates: Date
    root = normpath(joinpath(@__DIR__, ".."))
    include(joinpath(root, "scripts", "audit_onset_curve.jl"))

    scanned = read_scanned(joinpath(root, "data", "onset_curve_scanned.csv"))
    figures = read_figures(joinpath(root, "data", "onset_curve_figures.csv"))
    reprint_of = reprints(scanned)
    last_tick = Dict(sr => tick for (sr, _, tick) in CONFIG)

    ## One row per scanned vintage, in file order.
    @test [r["sitrep"] for r in figures] == scanned.order
    for r in figures
        sr = r["sitrep"]
        @test Date(r["report_date"]) == scanned.report_date[sr]
        ## The total is derived from the scanned block, so it must not
        ## drift from it.
        total = parse(Int, r["digitised_total"])
        @test total == block_total(scanned.blocks[sr])
        ## Reprints follow the loader's rule: exact equality over the
        ## block, pointing at the earliest vintage that carries it.
        @test r["reprint_of"] == reprint_of[sr]
        @test r["last_tick"] == string(last_tick[sr])
        if isempty(r["printed_n"])
            @test isempty(r["gap_pct"])
        else
            n = parse(Int, r["printed_n"])
            @test parse(Float64, r["gap_pct"]) ≈ 100 * (total - n) / n atol = 1.0e-3
        end
    end
end

@testitem "printed n parses from the figure title and source line" default_imports = false begin
    using Test
    root = normpath(joinpath(@__DIR__, ".."))
    include(joinpath(root, "scripts", "audit_onset_curve.jl"))

    ## As tesseract reads them off SitRep 118's title strip and source
    ## strip.
    title = "Nombre des cas confirmes par date de debut des symptomes " *
        "(n = 5 263)"
    source = "Source : DHIS2 Tracker - Riposte MVE, septembre 2026. " *
        "Donnees incluses :\ncas confirmes avec statut vital et date de " *
        "debut de symptomes renseignes, n = 5 263."
    @test parse_printed_n(title) == [5263]
    @test parse_printed_n(source) == [5263]
    @test parse_printed_n(title * "\n" * source) == [5263, 5263]
    ## Thousands separated by a thin space, a no-break space, a full stop
    ## or nothing at all.
    @test parse_printed_n("symptômes (n = 2 308)") == [2308]
    @test parse_printed_n("symptômes (n = 2 308)") == [2308]
    @test parse_printed_n("symptômes (n = 2.308)") == [2308]
    @test parse_printed_n("symptomes (n=2308)") == [2308]
    ## The report text prints n for provinces and deaths on lines that do
    ## not name the onset basis; those are not read.
    @test isempty(
        parse_printed_n(
            "Environ 90,5 % des cas (n=1 772) et 85,4 % des décès (n=608)"
        )
    )
    ## The SitRep 059-062 figures carry no title, so nothing parses.
    @test isempty(
        parse_printed_n("60 -: Premier résultat positif\ndu laboratoire")
    )
end

@testitem "audit reads the committed figure metadata off the SitRep PDFs" default_imports = false begin
    using Test
    root = normpath(joinpath(@__DIR__, ".."))
    pdf_dir = joinpath(root, "data", "sitrep_pdfs")

    if !isdir(pdf_dir)
        @info "onset figure audit check skipped: $pdf_dir absent. Fetch " *
            "the reports with scripts/download_sitreps.jl to run it."
        @test true
    else
        include(joinpath(root, "scripts", "audit_onset_curve.jl"))
        rows = Dict(
            r["sitrep"] => r for r in read_figures(
                    joinpath(root, "data", "onset_curve_figures.csv")
                )
        )
        ## 108 is the larger render that needs the near-grey fallback;
        ## 118 is a plateau vintage whose title and source strips agree.
        for sr in ("108", "118")
            pdf = joinpath(pdf_dir, "SitRep_MVE_$(sr)_2026.pdf")
            isfile(pdf) || continue
            fig = onset_figure(pdf)
            @test fig !== nothing
            fig === nothing && continue
            r = rows[sr]
            @test string(fig.page) == r["page"]
            @test string(size(fig.R, 2)) == r["width"]
            @test string(size(fig.R, 1)) == r["height"]
            isempty(fig.md5) || @test fig.md5 == r["image_md5"]
            cal = calibration(fig.R, fig.G, fig.B, get(Y_AXIS_STEP, sr, 20))
            @test parse(Float64, r["pixels_per_count"]) ≈
                cal.pixels_per_count atol = 1.0e-3
            @test parse(Float64, r["pixels_per_day"]) ≈
                cal.pixels_per_day atol = 1.0e-3
            if Sys.which("tesseract") !== nothing && !isempty(r["printed_n"])
                n, _, _ = printed_n(pdf, fig.page, fig.R, fig.G, fig.B)
                @test n == parse(Int, r["printed_n"])
            end
        end
    end
end
