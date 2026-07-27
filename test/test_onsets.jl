## Tests for the symptom-onset reporting-triangle stream: the loader
## (`src/onset_curve.jl`: parsing, dedup, cut-off filtering, increment
## construction), the reporting-delay hazard's pure functions
## (`src/models/observations.jl`), and the composers that fit it
## (`onsets_only_model`, `bvd_joint`, `src/models/joint.jl`).
##
## Filter target for a scoped run (this worktree only, not every sibling
## worktree `@run_package_tests` would otherwise discover):
##   target = joinpath(pwd(), "test", "test_onsets.jl")
##   @run_package_tests filter = ti -> string(ti.filename) == target

## --- Loader: parsing and dedup -------------------------------------------

@testitem "_dedup_onset_blocks collapses byte-identical reprints" begin
    using BVDOutbreakSize: _read_onset_curve_blocks, _dedup_onset_blocks
    using Dates: Date

    dir = mktempdir()
    path = joinpath(dir, "onset.csv")
    write(path, """
    sitrep,report_date,onset_date,confirmed_alive,confirmed_dead,confirmed_total
    001,2026-03-01,2026-02-25,2,0,2
    001,2026-03-01,2026-02-26,1,0,1
    002,2026-03-02,2026-02-25,2,0,2
    002,2026-03-02,2026-02-26,1,0,1
    003,2026-03-04,2026-02-25,3,0,3
    003,2026-03-04,2026-02-26,2,0,2
    003,2026-03-04,2026-02-27,1,0,1
    004,2026-03-06,2026-02-25,4,0,4
    004,2026-03-06,2026-02-26,2,0,2
    004,2026-03-06,2026-02-27,2,0,2
    """)
    blocks = _dedup_onset_blocks(_read_onset_curve_blocks(path))
    ## 001/002 reprint the same figure: collapse to 001, the earlier report
    ## date. 003 and 004 are genuinely new content and both survive.
    @test [b.sitrep for b in blocks] == ["001", "003", "004"]
    @test [b.report_date for b in blocks] ==
          [Date("2026-03-01"), Date("2026-03-04"), Date("2026-03-06")]
end

@testitem "_dedup_onset_blocks keeps equal-total, differently-split blocks distinct" begin
    ## The real data's 059/060 pair: equal cumulative totals but different
    ## per-date splits are NOT a reprint and must both survive.
    using BVDOutbreakSize: _read_onset_curve_blocks, _dedup_onset_blocks

    dir = mktempdir()
    path = joinpath(dir, "onset.csv")
    write(path, """
    sitrep,report_date,onset_date,confirmed_alive,confirmed_dead,confirmed_total
    010,2026-04-01,2026-03-30,2,0,2
    010,2026-04-01,2026-03-31,3,0,3
    011,2026-04-02,2026-03-30,3,0,3
    011,2026-04-02,2026-03-31,2,0,2
    """)
    blocks = _dedup_onset_blocks(_read_onset_curve_blocks(path))
    @test [b.sitrep for b in blocks] == ["010", "011"]
end

@testitem "load_onset_curve degrades gracefully for a missing file" begin
    using BVDOutbreakSize: load_onset_curve
    using Dates: Date

    h = load_onset_curve("/does/not/exist/onset_curve_scanned.csv";
        cutoff = Date("2026-03-01"), seeding = Date("2026-01-01"))
    @test h.onset_days == Int[]
    @test h.report_days == Int[]
    @test h.prev_report_days == Int[]
    @test h.increments == Int[]
    @test ismissing(h.last_total)
end

@testitem "load_onset_curve cut-off filtering recovers a later vintage with no code change" begin
    using BVDOutbreakSize: load_onset_curve
    using Dates: Date

    dir = mktempdir()
    path = joinpath(dir, "onset.csv")
    write(path, """
    sitrep,report_date,onset_date,confirmed_alive,confirmed_dead,confirmed_total
    001,2026-03-01,2026-02-25,2,0,2
    001,2026-03-01,2026-02-26,1,0,1
    003,2026-03-04,2026-02-25,3,0,3
    003,2026-03-04,2026-02-26,2,0,2
    003,2026-03-04,2026-02-27,1,0,1
    004,2026-03-06,2026-02-25,4,0,4
    004,2026-03-06,2026-02-26,2,0,2
    004,2026-03-06,2026-02-27,2,0,2
    """)
    seeding = Date("2026-01-01")

    ## Cut-off before 004's report date: 001 and 003 survive, so two
    ## scored snapshots (001 against the empty predecessor, then 003
    ## against 001).
    early = load_onset_curve(path; cutoff = Date("2026-03-04"), seeding)
    @test !isempty(early.onset_days)
    @test length(unique(early.report_days)) == 2

    ## Advancing the cut-off past 004's report date recovers it: the same
    ## loader call, no code change, exactly the self-correcting behaviour
    ## the manifest `as_of_date` advance relies on.
    late = load_onset_curve(path; cutoff = Date("2026-03-06"), seeding)
    @test length(unique(late.report_days)) == 3
    @test length(late.onset_days) > length(early.onset_days)
end

@testitem "load_onset_curve: missing onset date within a block's own extent reads as zero" begin
    using BVDOutbreakSize: load_onset_curve
    using Dates: Date, date2epochdays

    dir = mktempdir()
    path = joinpath(dir, "onset.csv")
    ## Both blocks print 03-01..03-03; the second omits the middle row
    ## (03-02), a zero-height bar the digitisation drops rather than a date
    ## the figure never covered.
    write(path, """
    sitrep,report_date,onset_date,confirmed_alive,confirmed_dead,confirmed_total
    001,2026-03-05,2026-03-01,1,0,1
    001,2026-03-05,2026-03-02,1,0,1
    001,2026-03-05,2026-03-03,1,0,1
    002,2026-03-07,2026-03-01,2,0,2
    002,2026-03-07,2026-03-03,1,0,1
    """)
    seeding = Date("2026-01-01")
    h = load_onset_curve(path; cutoff = Date("2026-03-07"), seeding,
        max_delay = 10)
    ## Grid index of 2026-03-02 (seeding 2026-01-01).
    u = Int(date2epochdays(Date("2026-03-02")) - date2epochdays(seeding)) + 1
    R2 = Int(date2epochdays(Date("2026-03-07")) - date2epochdays(seeding)) + 1
    idxs = findall(
        i -> h.onset_days[i] == u && h.report_days[i] == R2,
        eachindex(h.onset_days))
    @test !isempty(idxs)
    ## The omitted row reads as a true zero, so the cell scores 0 - 1 = -1,
    ## not a dropped cell and not an error.
    @test all(i -> h.increments[i] == -1, idxs)
end

@testitem "load_onset_curve: hand-built 3-snapshot triangle matches exact expected cells" begin
    using BVDOutbreakSize: load_onset_curve
    using Dates: Date

    dir = mktempdir()
    path = joinpath(dir, "onset.csv")
    ## A three-snapshot triangle whose printed extents run close to each
    ## report date, the shape the real digitised figures take. Grid days
    ## with seeding 2026-01-01: Feb 1 = 32, ..., Feb 14 = 45.
    ##   001 report Feb 10 (day 41), extent Feb 1-8   (32-39)
    ##   002 report Feb 12 (day 43), extent Feb 1-10  (32-41)
    ##   003 report Feb 14 (day 45), extent Feb 1-12  (32-43)
    write(path, """
    sitrep,report_date,onset_date,confirmed_alive,confirmed_dead,confirmed_total
    001,2026-02-10,2026-02-01,5,0,5
    001,2026-02-10,2026-02-02,3,0,3
    001,2026-02-10,2026-02-03,2,0,2
    001,2026-02-10,2026-02-04,4,0,4
    001,2026-02-10,2026-02-05,3,0,3
    001,2026-02-10,2026-02-06,2,0,2
    001,2026-02-10,2026-02-07,1,0,1
    001,2026-02-10,2026-02-08,1,0,1
    002,2026-02-12,2026-02-01,6,0,6
    002,2026-02-12,2026-02-02,3,0,3
    002,2026-02-12,2026-02-03,2,0,2
    002,2026-02-12,2026-02-04,5,0,5
    002,2026-02-12,2026-02-05,4,0,4
    002,2026-02-12,2026-02-06,3,0,3
    002,2026-02-12,2026-02-07,2,0,2
    002,2026-02-12,2026-02-08,2,0,2
    002,2026-02-12,2026-02-09,1,0,1
    002,2026-02-12,2026-02-10,1,0,1
    003,2026-02-14,2026-02-01,7,0,7
    003,2026-02-14,2026-02-02,3,0,3
    003,2026-02-14,2026-02-03,2,0,2
    003,2026-02-14,2026-02-04,5,0,5
    003,2026-02-14,2026-02-05,4,0,4
    003,2026-02-14,2026-02-06,2,0,2
    003,2026-02-14,2026-02-07,3,0,3
    003,2026-02-14,2026-02-08,2,0,2
    003,2026-02-14,2026-02-09,2,0,2
    003,2026-02-14,2026-02-10,2,0,2
    003,2026-02-14,2026-02-11,1,0,1
    003,2026-02-14,2026-02-12,1,0,1
    """)
    seeding = Date("2026-01-01")   # Jan 1 = grid day 1
    h = load_onset_curve(path; cutoff = Date("2026-02-14"), seeding,
        max_delay = 10)

    ## V1 (report day 41, virtual empty predecessor): window is its own
    ## extent 32:39 intersected with the trailing 10-day horizon 32:41.
    ## V2 (report day 43): both extents cover 32:39, horizon starts at 34.
    ## V3 (report day 45): both extents cover 32:41, horizon starts at 36.
    exp_onset = vcat(32:39, 34:39, 36:41)
    exp_report = vcat(fill(41, 8), fill(43, 6), fill(45, 6))
    exp_prev = vcat(fill(0, 8), fill(41, 6), fill(43, 6))
    exp_inc = vcat(
        [5, 3, 2, 4, 3, 2, 1, 1],      # V1 levels vs the virtual empty
        [0, 1, 1, 1, 1, 1],            # V2 vs V1
        [0, -1, 1, 0, 1, 1])           # V3 vs V2 (Feb 6 revised 3 -> 2)

    @test h.onset_days == exp_onset
    @test h.report_days == exp_report
    @test h.prev_report_days == exp_prev
    @test h.increments == exp_inc
    ## The deliberately negative correction (Feb 6 revised 3 -> 2) survives
    ## into the scored increments rather than being clamped.
    @test any(<(0), h.increments)
    @test h.last_total == 34   # V3's cumulative total
end

@testitem "load_onset_curve: horizon window excludes settled onset dates" begin
    using BVDOutbreakSize: load_onset_curve
    using Dates: Date

    dir = mktempdir()
    path = joinpath(dir, "onset.csv")
    ## A three-snapshot triangle whose printed extents run close to each
    ## report date, the shape the real digitised figures take. Grid days
    ## with seeding 2026-01-01: Feb 1 = 32, ..., Feb 14 = 45.
    ##   001 report Feb 10 (day 41), extent Feb 1-8   (32-39)
    ##   002 report Feb 12 (day 43), extent Feb 1-10  (32-41)
    ##   003 report Feb 14 (day 45), extent Feb 1-12  (32-43)
    write(path, """
    sitrep,report_date,onset_date,confirmed_alive,confirmed_dead,confirmed_total
    001,2026-02-10,2026-02-01,5,0,5
    001,2026-02-10,2026-02-02,3,0,3
    001,2026-02-10,2026-02-03,2,0,2
    001,2026-02-10,2026-02-04,4,0,4
    001,2026-02-10,2026-02-05,3,0,3
    001,2026-02-10,2026-02-06,2,0,2
    001,2026-02-10,2026-02-07,1,0,1
    001,2026-02-10,2026-02-08,1,0,1
    002,2026-02-12,2026-02-01,6,0,6
    002,2026-02-12,2026-02-02,3,0,3
    002,2026-02-12,2026-02-03,2,0,2
    002,2026-02-12,2026-02-04,5,0,5
    002,2026-02-12,2026-02-05,4,0,4
    002,2026-02-12,2026-02-06,3,0,3
    002,2026-02-12,2026-02-07,2,0,2
    002,2026-02-12,2026-02-08,2,0,2
    002,2026-02-12,2026-02-09,1,0,1
    002,2026-02-12,2026-02-10,1,0,1
    003,2026-02-14,2026-02-01,7,0,7
    003,2026-02-14,2026-02-02,3,0,3
    003,2026-02-14,2026-02-03,2,0,2
    003,2026-02-14,2026-02-04,5,0,5
    003,2026-02-14,2026-02-05,4,0,4
    003,2026-02-14,2026-02-06,2,0,2
    003,2026-02-14,2026-02-07,3,0,3
    003,2026-02-14,2026-02-08,2,0,2
    003,2026-02-14,2026-02-09,2,0,2
    003,2026-02-14,2026-02-10,2,0,2
    003,2026-02-14,2026-02-11,1,0,1
    003,2026-02-14,2026-02-12,1,0,1
    """)
    seeding = Date("2026-01-01")
    h = load_onset_curve(path; cutoff = Date("2026-02-14"), seeding,
        max_delay = 10)
    ## Feb 1-2 (days 32-33) sit 11-12 days before the 002 report day (43):
    ## older than the horizon (10), so their revisions are not re-scored in
    ## the 002 window even though both figures print them.
    v2_days = h.onset_days[h.report_days .== 43]
    @test 32 ∉ v2_days
    @test 33 ∉ v2_days
    @test minimum(v2_days) == 34
end

@testitem "load_onset_curve: onset dates past a figure's printed extent are dropped, not zeroed" begin
    ## The published figures stop their x axis short of the report date. A
    ## date the figure never covered carries no observation, so its cell is
    ## dropped; reading it as zero would assert nothing had been reported
    ## yet and bias the fitted hazard towards slow reporting at exactly the
    ## delays the axis gap spans.
    using BVDOutbreakSize: load_onset_curve
    using Dates: Date, date2epochdays

    dir = mktempdir()
    path = joinpath(dir, "onset.csv")
    ## A three-snapshot triangle whose printed extents run close to each
    ## report date, the shape the real digitised figures take. Grid days
    ## with seeding 2026-01-01: Feb 1 = 32, ..., Feb 14 = 45.
    ##   001 report Feb 10 (day 41), extent Feb 1-8   (32-39)
    ##   002 report Feb 12 (day 43), extent Feb 1-10  (32-41)
    ##   003 report Feb 14 (day 45), extent Feb 1-12  (32-43)
    write(path, """
    sitrep,report_date,onset_date,confirmed_alive,confirmed_dead,confirmed_total
    001,2026-02-10,2026-02-01,5,0,5
    001,2026-02-10,2026-02-02,3,0,3
    001,2026-02-10,2026-02-03,2,0,2
    001,2026-02-10,2026-02-04,4,0,4
    001,2026-02-10,2026-02-05,3,0,3
    001,2026-02-10,2026-02-06,2,0,2
    001,2026-02-10,2026-02-07,1,0,1
    001,2026-02-10,2026-02-08,1,0,1
    002,2026-02-12,2026-02-01,6,0,6
    002,2026-02-12,2026-02-02,3,0,3
    002,2026-02-12,2026-02-03,2,0,2
    002,2026-02-12,2026-02-04,5,0,5
    002,2026-02-12,2026-02-05,4,0,4
    002,2026-02-12,2026-02-06,3,0,3
    002,2026-02-12,2026-02-07,2,0,2
    002,2026-02-12,2026-02-08,2,0,2
    002,2026-02-12,2026-02-09,1,0,1
    002,2026-02-12,2026-02-10,1,0,1
    003,2026-02-14,2026-02-01,7,0,7
    003,2026-02-14,2026-02-02,3,0,3
    003,2026-02-14,2026-02-03,2,0,2
    003,2026-02-14,2026-02-04,5,0,5
    003,2026-02-14,2026-02-05,4,0,4
    003,2026-02-14,2026-02-06,2,0,2
    003,2026-02-14,2026-02-07,3,0,3
    003,2026-02-14,2026-02-08,2,0,2
    003,2026-02-14,2026-02-09,2,0,2
    003,2026-02-14,2026-02-10,2,0,2
    003,2026-02-14,2026-02-11,1,0,1
    003,2026-02-14,2026-02-12,1,0,1
    """)
    seeding = Date("2026-01-01")
    h = load_onset_curve(path; cutoff = Date("2026-02-14"), seeding,
        max_delay = 10)
    _day(s) = Int(date2epochdays(Date(s)) - date2epochdays(seeding)) + 1

    ## 001's axis ends at Feb 8 (day 39), two days before its report day, so
    ## Feb 9 and Feb 10 are never scored against it.
    v1 = h.onset_days[h.report_days .== _day("2026-02-10")]
    @test maximum(v1) == _day("2026-02-08")
    ## The 002 correction cells stop at Feb 8 too: a correction needs both
    ## figures to print the date, and 001 is the narrower of the pair.
    v2 = h.onset_days[h.report_days .== _day("2026-02-12")]
    @test maximum(v2) == _day("2026-02-08")
    ## No cell anywhere reaches delay 0 or 1 for this triangle.
    @test minimum(h.report_days .- h.onset_days) == 2
end

## --- Hazard / CDF pure functions ------------------------------------------

@testitem "onset_report_cdf: truncation, range and monotonicity" begin
    using BVDOutbreakSize: onset_report_cdf
    using StatsFuns: logit

    logit_h0 = fill(logit(0.1), 28)
    γ = zeros(60)
    u = 5
    vals = [onset_report_cdf(δ, logit_h0, γ, u, 1) for δ in (-5):(28 + 5)]

    ## δ < 0 is exact right truncation: F = 0.
    @test all(==(0.0), vals[1:5])
    ## F always lies in [0, 1].
    @test all(v -> 0 <= v <= 1, vals)
    ## Monotone non-decreasing in δ for fixed hazards.
    @test issorted(vals)
    ## Saturates (constant) once δ >= D - 1 = 27, since the hazard has no
    ## support beyond the tracked delay window.
    i27 = findfirst(==(27), (-5):(28 + 5))
    i30 = findfirst(==(30), (-5):(28 + 5))
    @test vals[i27] == vals[i30]
end

@testitem "onset_report_cdf_extrapolated agrees with onset_report_cdf in-range and stays safe out-of-range" begin
    using BVDOutbreakSize: onset_report_cdf, onset_report_cdf_extrapolated
    using StatsFuns: logit

    logit_h0 = fill(logit(0.15), 28)
    γ = collect(range(-0.5, 0.5; length = 30))
    grid_start = 10

    ## In-range: identical to `onset_report_cdf`.
    for u in grid_start:(grid_start + 5), δ in 0:10

        @test onset_report_cdf_extrapolated(δ, logit_h0, γ, u, grid_start) ≈
              onset_report_cdf(δ, logit_h0, γ, u, grid_start)
    end

    ## Out-of-range on both sides: finite, in [0, 1], no bounds error, and
    ## δ < 0 is still exact right truncation.
    for u in (-20, 0, 1, 5000)
        @test onset_report_cdf_extrapolated(-1, logit_h0, γ, u, grid_start) ==
              0.0
        for δ in (0, 10, 27, 40)
            v = onset_report_cdf_extrapolated(δ, logit_h0, γ, u, grid_start)
            @test isfinite(v)
            @test 0 <= v <= 1
        end
    end
end

@testitem "onset_report_moments: a later snapshot sees more of the same onset date" begin
    ## The literal right-truncation proof: for the SAME onset date and fixed
    ## hazards, the modelled current-cumulative level is non-decreasing as
    ## the report day moves later.
    using BVDOutbreakSize: onset_report_moments
    using StatsFuns: logit

    onsets = fill(10.0, 60)
    logit_h0 = fill(logit(0.1), 28)
    γ = zeros(60)
    u = 20
    early = onset_report_moments(onsets, logit_h0, γ, 1, [u], [u + 3], [0])
    late = onset_report_moments(onsets, logit_h0, γ, 1, [u], [u + 10], [0])
    @test late.level_cur[1] >= early.level_cur[1]
end

@testitem "onset_report_expected_total / onset_report_ascertainment stay in bounds" begin
    using BVDOutbreakSize: onset_report_expected_total,
                           onset_report_ascertainment
    using StatsFuns: logit

    onsets = fill(5.0, 100)
    logit_h0 = fill(logit(0.1), 28)
    for (grid_start, grid_end) in ((1, 100), (1, 10), (5, 40), (50, 55))
        γ = zeros(max(grid_end - grid_start + 1, 1))
        total = onset_report_expected_total(onsets, logit_h0, γ, grid_start,
            grid_end)
        @test isfinite(total)
        @test total >= 0
        asc = onset_report_ascertainment(logit_h0, γ, grid_start, grid_end)
        @test all(a -> 0 <= a <= 1, asc)
        ## Degenerate short grid (grid_end - grid_start + 1 < D): no `u` has
        ## a fully computable asymptote, so an empty vector is returned
        ## rather than reading `γ` out of bounds.
        if grid_end - grid_start + 1 < length(logit_h0)
            @test isempty(asc)
        end
    end
end

@testitem "onset_report_expected_total covers onset dates before grid_start" begin
    ## Regression test: `expected_onset_reported_T` must sum the FULL
    ## `1:n` onset series like every other stream's `expected_*_T`, not
    ## just the triangle's own `grid_start:grid_end` window (see
    ## `onset_report_cdf_extrapolated`).
    using BVDOutbreakSize: onset_report_expected_total, onset_report_cdf
    using StatsFuns: logit

    n = 100
    onsets = fill(5.0, n)
    logit_h0 = fill(logit(0.2), 28)
    grid_start = 60
    grid_end = 90
    γ = zeros(grid_end - grid_start + 1)

    ## A window-restricted sum (the pre-fix behaviour) only over
    ## `grid_start:grid_end`, for comparison.
    restricted = sum(
        onsets[u] * onset_report_cdf(grid_end - u, logit_h0, γ, u, grid_start)
    for u in grid_start:grid_end)
    total = onset_report_expected_total(onsets, logit_h0, γ, grid_start,
        grid_end)

    ## The onset dates before `grid_start` (days 1:59) are old enough by
    ## `grid_end` that they sit at the hazard's asymptote and each
    ## contribute onsets[u] * F(u, D-1) > 0, so the full total must exceed
    ## the window-restricted sum.
    @test total > restricted
    @test isfinite(total)

    ## The extrapolated contribution for a day well before `grid_start`
    ## should match the flat asymptote computed at the earliest known
    ## calendar day (γ held at γ[1] = 0 here).
    F_edge = onset_report_cdf(length(logit_h0) - 1, logit_h0, γ, grid_start,
        grid_start)
    @test total ≈ restricted + (grid_start - 1) * 5.0 * F_edge
end

@testitem "load_onset_curve: a vintage narrower than an earlier one drops the uncovered date" begin
    ## A date one vintage of a pair does not print carries no observation
    ## from that vintage, so the pair cannot form a correction there. The
    ## cell is dropped rather than read as a zero, which would fabricate a
    ## negative correction against the earlier vintage's positive count.
    using BVDOutbreakSize: load_onset_curve
    using Dates: Date, date2epochdays

    dir = mktempdir()
    path = joinpath(dir, "onset.csv")
    ## Block 001 covers 03-01..03-03. Block 002 covers only 03-02..03-03,
    ## so 03-01 sits outside its printed extent.
    write(path, """
    sitrep,report_date,onset_date,confirmed_alive,confirmed_dead,confirmed_total
    001,2026-03-05,2026-03-01,4,0,4
    001,2026-03-05,2026-03-02,2,0,2
    001,2026-03-05,2026-03-03,1,0,1
    002,2026-03-07,2026-03-02,3,0,3
    002,2026-03-07,2026-03-03,2,0,2
    """)
    seeding = Date("2026-01-01")
    h = load_onset_curve(path; cutoff = Date("2026-03-07"), seeding,
        max_delay = 10)
    _day(x) = Int(date2epochdays(Date(x)) - date2epochdays(seeding)) + 1
    u = _day("2026-03-01")
    R2 = _day("2026-03-07")
    ## No cell for 03-01 in the 001-versus-002 pair.
    @test isempty(findall(i -> h.onset_days[i] == u && h.report_days[i] == R2,
        eachindex(h.onset_days)))
    ## The first pair still scores 03-01 as a level against the empty
    ## predecessor, since block 001 does print it.
    R1 = _day("2026-03-05")
    first_idx = findall(
        i -> h.onset_days[i] == u && h.report_days[i] == R1,
        eachindex(h.onset_days))
    @test length(first_idx) == 1
    @test h.increments[first_idx[1]] == 4
    ## Nothing anywhere fabricates the -4 the dropped cell would have given.
    @test minimum(h.increments) >= 0
end

@testitem "onset_report_scales matches the measured-error formula and grows with magnitude" begin
    using BVDOutbreakSize: onset_report_scales

    level_cur = [0.0, 100.0, 40.0]
    level_prev = [0.0, 80.0, 0.0]
    means = level_cur .- level_prev
    ## Cells 1 and 3 have a virtual (empty) predecessor and so score a
    ## level; cell 2 is a correction between two real snapshots.
    prev_idx = [0, 5, 0]
    s = onset_report_scales(means, level_cur, level_prev, prev_idx;
        pixel_sd = 2.1, scan_frac = 0.04)
    @test s[1] ≈ sqrt(2.1^2 * 1)
    @test s[2] ≈ sqrt(20.0 + 2.1^2 * 2 + 0.04^2 * (100.0^2 + 80.0^2))
    ## A level cell carries the counting variation of the cases it reports,
    ## which for a bar of 40 dominates the ≈2.1-case reading error.
    @test s[3] ≈ sqrt(40.0 + 2.1^2 * 1 + 0.04^2 * 40.0^2)
    @test s[3] > sqrt(40.0)
    ## The scale grows with the modelled magnitude.
    @test s[2] > s[1]
end

@testitem "onset_report_scales floors the counting term at a non-negative mean" begin
    ## `means` is non-negative by construction (F is monotone in δ), but a
    ## degenerate call must not take the square root of a negative variance.
    using BVDOutbreakSize: onset_report_scales

    s = onset_report_scales([-5.0], [1.0], [6.0], [3])
    @test isfinite(s[1])
    @test s[1] > 0
end

@testitem "safe_studentt stays valid under extreme scale/df" begin
    using BVDOutbreakSize: safe_studentt
    using Distributions: mean, std, logpdf

    for (σ, ν) in ((0.0, 4.0), (-1.0, 4.0), (NaN, 4.0), (Inf, 4.0),
        (1.0, 0.0), (1.0, -2.0), (1.0, NaN))
        d = safe_studentt(3.0, σ, ν)
        ## Mean and variance both exist: a degenerate degrees-of-freedom
        ## argument falls back to 4, not to the Cauchy at the domain edge.
        @test isfinite(mean(d))
        @test isfinite(std(d))
        @test isfinite(logpdf(d, 5.0))
    end
    ## A well-posed call passes through as the requested location and scale.
    d = safe_studentt(3.0, 2.0, 4.0)
    @test mean(d) ≈ 3.0
    @test std(d) ≈ 2.0 * sqrt(4 / 2)
    ## A caller-chosen heavy tail is respected rather than overridden.
    @test !isfinite(std(safe_studentt(0.0, 1.0, 1.5)))
end

## --- Model level ------------------------------------------------------

@testitem "onsets_only_model: default empty history is a no-op" begin
    using BVDOutbreakSize: onsets_only_model
    using Turing: Prior, sample

    chn = sample(onsets_only_model(30), Prior(), 5; progress = false)
    et = vec(Array(chn[:expected_onset_reported_T]))
    @test length(et) == 5
    @test all(isfinite, et)
end

@testitem "onset_reporting_model conditions on the increments rather than sampling them" begin
    ## Regression test for a silent-failure mode specific to DynamicPPL:
    ## observe-versus-assume is decided by whether the tilde's symbol is one
    ## of the enclosing model's argument names, so reading the observations
    ## out of a container into a local variable and writing `local[i] ~ d`
    ## turns every cell into a latent parameter, drops the likelihood, and
    ## still samples, differentiates and runs NUTS without complaint. Two
    ## checks: no cell appears as a random variable, and the log-density
    ## actually responds to the data.
    using BVDOutbreakSize: onsets_only_model
    using Turing: DynamicPPL
    using LogDensityProblems: logdensity
    using Random: seed!

    base = (; onset_days = [10, 11, 12, 13, 10, 11, 12, 13, 14],
        report_days = [15, 15, 15, 15, 20, 20, 20, 20, 20],
        prev_report_days = [0, 0, 0, 0, 15, 15, 15, 15, 0])
    oc = (; base..., increments = [2, 3, 1, 0, 1, 2, 3, 4, 5])
    seed!(20260727)
    model = onsets_only_model(40; onset_curve_history = oc)
    vi = DynamicPPL.VarInfo(model)
    names = string.(collect(keys(vi)))
    @test !any(n -> occursin("increments", n), names)

    ## Same latent draw, different data: the log-joint must move. Direction
    ## is not asserted, since at an arbitrary prior draw either data set can
    ## sit closer to the modelled increments.
    other = (; base..., increments = [90, 80, 70, 60, 50, 40, 30, 20, 10])
    model_other = onsets_only_model(40; onset_curve_history = other)
    θ = collect(vi[:])
    lp = logdensity(
        DynamicPPL.LogDensityFunction(model, DynamicPPL.getlogjoint, vi), θ)
    lp_other = logdensity(
        DynamicPPL.LogDensityFunction(model_other, DynamicPPL.getlogjoint,
            vi), θ)
    @test isfinite(lp)
    @test isfinite(lp_other)
    @test lp != lp_other
end

@testitem "onsets_only_model: prior predictive on a synthetic triangle is finite" begin
    using BVDOutbreakSize: onsets_only_model
    using Turing: Prior, sample

    oc = (; onset_days = [10, 11, 12, 13, 10, 11, 12, 13, 14],
        report_days = [15, 15, 15, 15, 20, 20, 20, 20, 20],
        prev_report_days = [0, 0, 0, 0, 15, 15, 15, 15, 0],
        increments = [2, 3, 1, 0, 1, 2, 3, 4, 5])
    chn = sample(onsets_only_model(40; onset_curve_history = oc), Prior(),
        20; progress = false)
    et = vec(Array(chn[:expected_onset_reported_T]))
    @test length(et) == 20
    @test all(isfinite, et)
    @test all(>=(0), et)
end

@testitem "AD gradient: onsets_only_model differentiates (Mooncake)" tags = [
    :ad] begin
    using Turing: DynamicPPL
    using LogDensityProblems: logdensity_and_gradient
    using Random: seed!
    using BVDOutbreakSize: onsets_only_model, default_adtype

    oc = (; onset_days = [10, 11, 12, 13, 10, 11, 12, 13, 14],
        report_days = [15, 15, 15, 15, 20, 20, 20, 20, 20],
        prev_report_days = [0, 0, 0, 0, 15, 15, 15, 15, 0],
        increments = [2, 3, 1, 0, 1, 2, 3, 4, 5])

    seed!(20260518)
    model = onsets_only_model(40; onset_curve_history = oc)
    vi = DynamicPPL.link(DynamicPPL.VarInfo(model), model)
    x0 = collect(vi[:])
    ldf = DynamicPPL.LogDensityFunction(
        model, DynamicPPL.getlogjoint, vi; adtype = default_adtype())
    logp, grad = logdensity_and_gradient(ldf, x0)
    @test isfinite(logp)
    @test length(grad) == length(x0)
    @test all(isfinite, grad)
    @test any(!iszero, grad)
end

@testitem "onsets_only_model fits under a short NUTS run" tags = [:slow] begin
    using BVDOutbreakSize: onsets_only_model, nuts_sample

    oc = (; onset_days = [10, 11, 12, 13, 10, 11, 12, 13, 14],
        report_days = [15, 15, 15, 15, 20, 20, 20, 20, 20],
        prev_report_days = [0, 0, 0, 0, 15, 15, 15, 15, 0],
        increments = [2, 3, 1, 0, 1, 2, 3, 4, 5])
    chn = nuts_sample(onsets_only_model(40; onset_curve_history = oc);
        samples = 25, chains = 1, progress = false)
    et = vec(Array(chn[:expected_onset_reported_T]))
    @test length(et) == 25
    @test all(isfinite, et)
end

@testitem "bvd_joint fits under a short NUTS run with the onset stream wired in" tags = [
    :slow] begin
    using BVDOutbreakSize: bvd_joint, nuts_sample

    n = 40
    dh = (; days = [13, 18, 40], counts = [10, 14, 18])
    rh = (; days = [13, 18, 40], counts = [340, 516, 905])
    ch = (; days = [13, 18, 40], counts = [9, 17, 27])
    oc = (; onset_days = [20, 21, 22, 23, 20, 21, 22, 23, 24],
        report_days = [25, 25, 25, 25, 30, 30, 30, 30, 30],
        prev_report_days = [0, 0, 0, 0, 25, 25, 25, 25, 0],
        increments = [3, 2, 1, 0, 1, 2, 1, 3, 2])
    chn = nuts_sample(
        bvd_joint(n, 2, 18, 905, 0, 27, 50;
            confirmed_deaths = 5,
            deaths_history = dh,
            reported_history = rh,
            confirmed_history = ch,
            lab_history = (; days = [18, 40], counts = [30, 50]),
            onset_curve_history = oc,
            breakpoint = 30);
        samples = 12, chains = 1, progress = false)
    C_T = vec(Array(chn[:C_T]))
    et = vec(Array(chn[:expected_onset_reported_T]))
    @test length(C_T) == 12
    @test all(isfinite, C_T)
    @test all(C_T .> 0)
    @test length(et) == 12
    @test all(isfinite, et)
end
