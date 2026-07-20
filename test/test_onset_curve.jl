## Structural check on the digitised symptom-onset dataset
## (data/onset_curve_scanned.csv). That file is produced on the Python side
## by scripts/digitize_onset_curve.py and is NOT read by the model, so this
## guards its shape and internal consistency from the Julia test suite: a bad
## regeneration (wrong columns, a split that does not sum, an onset after the
## report date, or a vintage total that shrank) is caught in CI.

@testitem "onset_curve_scanned.csv is well-formed" begin
    using BVDOutbreakSize: BVDOutbreakSize
    using Dates: Date

    path = joinpath(pkgdir(BVDOutbreakSize), "data",
        "onset_curve_scanned.csv")
    @test isfile(path)

    lines = filter(!isempty, readlines(path))
    header = split(lines[1], ',')
    @test header == ["sitrep", "report_date", "onset_date",
        "confirmed_alive", "confirmed_dead", "confirmed_total"]

    rows = [split(l, ',') for l in lines[2:end]]
    sitrep = [r[1] for r in rows]
    report_date = [Date(r[2]) for r in rows]
    onset_date = [Date(r[3]) for r in rows]
    alive = [parse(Int, r[4]) for r in rows]
    dead = [parse(Int, r[5]) for r in rows]
    total = [parse(Int, r[6]) for r in rows]

    ## counts are non-negative and the outcome split sums to the total
    @test all(≥(0), alive)
    @test all(≥(0), dead)
    @test alive .+ dead == total

    ## an onset date cannot post-date the report it was scanned from
    @test all(onset_date .≤ report_date)

    ## the five scanned analytique vintages are all present
    @test Set(sitrep) == Set(["059", "060", "061", "062", "064"])

    ## per-vintage confirmed totals are non-decreasing with report date:
    ## later reports backfill recent onsets, so no vintage total can shrink
    vintages = sort(unique(collect(zip(sitrep, report_date))); by = last)
    totals = [sum(total[sitrep .== s]) for (s, _) in vintages]
    @test first(totals) > 0
    @test issorted(totals)
end
