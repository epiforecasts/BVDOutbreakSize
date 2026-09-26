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
    write(
        path, """
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
        """
    )
    blocks = _dedup_onset_blocks(_read_onset_curve_blocks(path))
    ## 001/002 reprint the same figure: collapse to 001, the earlier report
    ## date. 003 and 004 are genuinely new content and both survive.
    @test [b.sitrep for b in blocks] == ["001", "003", "004"]
    @test [b.report_date for b in blocks] ==
        [Date("2026-03-01"), Date("2026-03-04"), Date("2026-03-06")]
end

@testitem "_dedup_onset_blocks keeps distinct equal-total split blocks" begin
    ## The real data's 059/060 pair: equal cumulative totals but different
    ## per-date splits are not a reprint and must both survive.
    using BVDOutbreakSize: _read_onset_curve_blocks, _dedup_onset_blocks

    dir = mktempdir()
    path = joinpath(dir, "onset.csv")
    write(
        path, """
        sitrep,report_date,onset_date,confirmed_alive,confirmed_dead,confirmed_total
        010,2026-04-01,2026-03-30,2,0,2
        010,2026-04-01,2026-03-31,3,0,3
        011,2026-04-02,2026-03-30,3,0,3
        011,2026-04-02,2026-03-31,2,0,2
        """
    )
    blocks = _dedup_onset_blocks(_read_onset_curve_blocks(path))
    @test [b.sitrep for b in blocks] == ["010", "011"]
end

@testitem "load_onset_curve degrades gracefully for a missing file" begin
    using BVDOutbreakSize: load_onset_curve
    using Dates: Date

    h = load_onset_curve(
        "/does/not/exist/onset_curve_scanned.csv";
        cutoff = Date("2026-03-01"), seeding = Date("2026-01-01")
    )
    @test h.onset_days == Int[]
    @test h.report_days == Int[]
    @test h.prev_report_days == Int[]
    @test h.increments == Int[]
    @test ismissing(h.last_total)
end

@testitem "load_onset_curve cut-off filtering recovers a later vintage" begin
    using BVDOutbreakSize: load_onset_curve
    using Dates: Date

    dir = mktempdir()
    path = joinpath(dir, "onset.csv")
    write(
        path, """
        sitrep,report_date,onset_date,confirmed_alive,confirmed_dead,confirmed_total
        001,2026-03-01,2026-02-25,2,0,2
        001,2026-03-01,2026-02-26,1,0,1
        003,2026-03-04,2026-02-25,3,0,3
        003,2026-03-04,2026-02-26,2,0,2
        003,2026-03-04,2026-02-27,1,0,1
        004,2026-03-06,2026-02-25,4,0,4
        004,2026-03-06,2026-02-26,2,0,2
        004,2026-03-06,2026-02-27,2,0,2
        """
    )
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

@testitem "load_onset_curve: a date missing within a block reads as zero" begin
    using BVDOutbreakSize: load_onset_curve
    using Dates: Date, date2epochdays

    dir = mktempdir()
    path = joinpath(dir, "onset.csv")
    ## Both blocks print 03-01..03-03; the second omits the middle row
    ## (03-02), a zero-height bar the digitisation drops rather than a date
    ## the figure never covered.
    write(
        path, """
        sitrep,report_date,onset_date,confirmed_alive,confirmed_dead,confirmed_total
        001,2026-03-05,2026-03-01,1,0,1
        001,2026-03-05,2026-03-02,1,0,1
        001,2026-03-05,2026-03-03,1,0,1
        002,2026-03-07,2026-03-01,2,0,2
        002,2026-03-07,2026-03-03,1,0,1
        """
    )
    seeding = Date("2026-01-01")
    h = load_onset_curve(path; cutoff = Date("2026-03-07"), seeding)
    ## Grid index of 2026-03-02 (seeding 2026-01-01).
    u = Int(date2epochdays(Date("2026-03-02")) - date2epochdays(seeding)) + 1
    R2 = Int(date2epochdays(Date("2026-03-07")) - date2epochdays(seeding)) + 1
    idxs = findall(
        i -> h.onset_days[i] == u && h.report_days[i] == R2,
        eachindex(h.onset_days)
    )
    @test !isempty(idxs)
    ## The omitted row reads as a true zero, so the cell scores 0 - 1 = -1,
    ## not a dropped cell and not an error.
    @test all(i -> h.increments[i] == -1, idxs)
end

@testitem "load_onset_curve: 3-snapshot triangle matches expected cells" begin
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
    write(
        path, """
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
        """
    )
    seeding = Date("2026-01-01")   # Jan 1 = grid day 1
    h = load_onset_curve(path; cutoff = Date("2026-02-14"), seeding)

    ## Each vintage corrects the dates an earlier one printed and scores a
    ## level for each date it prints first.
    ## V1 (report day 41): levels over its extent, 32:39.
    ## V2 (report day 43): corrections 32:39, first-print levels 40:41.
    ## V3 (report day 45): corrections 32:41, first-print levels 42:43.
    exp_onset = vcat(32:39, 32:41, 32:43)
    exp_report = vcat(fill(41, 8), fill(43, 10), fill(45, 12))
    exp_prev = vcat(fill(0, 8), fill(41, 8), [0, 0], fill(43, 10), [0, 0])
    exp_inc = vcat(
        [5, 3, 2, 4, 3, 2, 1, 1],          # V1 levels
        [1, 0, 0, 1, 1, 1, 1, 1], [1, 1],  # V2 corrections, then levels
        [1, 0, 0, 0, 0, -1, 1, 0, 1, 1],   # V3 corrections (Feb 6 3 -> 2)
        [1, 1]                             # V3 levels
    )

    @test h.onset_days == exp_onset
    @test h.report_days == exp_report
    @test h.prev_report_days == exp_prev
    @test h.increments == exp_inc
    ## The deliberately negative correction (Feb 6 revised 3 -> 2) survives
    ## into the scored increments rather than being clamped.
    @test any(<(0), h.increments)
    @test h.last_total == 34   # V3's cumulative total
end

@testitem "load_onset_curve: out-of-extent dates are unobserved, not zeroed" begin
    ## The published figures stop their x axis short of the report date. A
    ## date the figure never covered carries no observation there; reading
    ## it as zero would assert nothing had been reported yet and bias the
    ## fitted hazard towards slow reporting at exactly the delays the axis
    ## gap spans. Its first cell is the level where a figure first prints
    ## it.
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
    write(
        path, """
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
        """
    )
    seeding = Date("2026-01-01")
    h = load_onset_curve(path; cutoff = Date("2026-02-14"), seeding)
    _day(s) = Int(date2epochdays(Date(s)) - date2epochdays(seeding)) + 1

    ## 001's axis ends at Feb 8 (day 39), two days before its report day, so
    ## Feb 9 and Feb 10 are never scored against it.
    v1 = h.onset_days[h.report_days .== _day("2026-02-10")]
    @test maximum(v1) == _day("2026-02-08")
    ## 002 corrects only the dates 001 printed, and scores Feb 9 and Feb 10
    ## as first-print levels.
    v2 = h.report_days .== _day("2026-02-12")
    corr2 = h.onset_days[v2 .& (h.prev_report_days .> 0)]
    @test maximum(corr2) == _day("2026-02-08")
    lev2 = h.onset_days[v2 .& (h.prev_report_days .== 0)]
    @test lev2 == [_day("2026-02-09"), _day("2026-02-10")]
    ## No cell anywhere reaches delay 0 or 1 for this triangle.
    @test minimum(h.report_days .- h.onset_days) == 2
end

@testitem "load_onset_curve: one level per onset date, corrections inside the support" begin
    ## A synthetic 3-vintage triangle with printed extents wider than the
    ## delay support. Each onset date scores one level where it is first
    ## printed and then corrections only while the delay is inside the
    ## support.
    using BVDOutbreakSize: load_onset_curve, ONSET_REPORT_MAX_DELAY
    using Dates: Date, Day

    dir = mktempdir()
    path = joinpath(dir, "onset.csv")
    seeding = Date("2026-01-01")
    lines = [
        "sitrep,report_date,onset_date,confirmed_alive,confirmed_dead," *
            "confirmed_total",
    ]
    for u in 1:58
        push!(lines, "001,2026-03-01,$(seeding + Day(u - 1)),0,0,1")
    end
    for u in 1:63
        push!(lines, "002,2026-03-06,$(seeding + Day(u - 1)),0,0,2")
    end
    for u in 1:68
        push!(lines, "003,2026-03-11,$(seeding + Day(u - 1)),0,0,3")
    end
    write(path, join(lines, "\n"))

    D = ONSET_REPORT_MAX_DELAY
    h = load_onset_curve(path; cutoff = Date("2026-03-11"), seeding)
    R1, R2, R3 = 60, 65, 70
    levels = h.prev_report_days .== 0
    ## One level per printed onset date, at the vintage that printed it
    ## first.
    @test sort(h.onset_days[levels]) == collect(1:68)
    first_print = Dict(zip(h.onset_days[levels], h.report_days[levels]))
    @test all(first_print[u] == R1 for u in 1:58)
    @test all(first_print[u] == R2 for u in 59:63)
    @test all(first_print[u] == R3 for u in 64:68)
    ## Corrections only inside the delay support, each against the
    ## previous vintage.
    corr = .!levels
    @test all(h.report_days[corr] .- h.onset_days[corr] .< D)
    @test sort(h.onset_days[corr .& (h.report_days .== R2)]) ==
        collect((R2 - D + 1):58)
    @test sort(h.onset_days[corr .& (h.report_days .== R3)]) ==
        collect((R3 - D + 1):63)
    @test all(==(1), h.increments[corr])
    ## Every level is the count printed at first print.
    @test all(
        h.increments[i] == (
            h.report_days[i] == R1 ? 1 :
                h.report_days[i] == R2 ? 2 : 3
        )
            for i in findall(levels)
    )
end

@testitem "onset_report_moments: cells with both reads beyond the support score zero" begin
    ## `onset_report_moments` is a pure per-cell function: a cell whose two
    ## reads are both beyond the delay support has both reads saturated at
    ## the same asymptote, so its mean is exactly 0. The loader never builds
    ## such a cell, so the cells are set here directly.
    using BVDOutbreakSize: onset_report_moments, ONSET_REPORT_MAX_DELAY

    D = ONSET_REPORT_MAX_DELAY
    h = (;
        onset_days = vcat(1:58, 1:63),
        report_days = vcat(fill(60, 58), fill(65, 63)),
        prev_report_days = vcat(fill(0, 58), fill(60, 63)),
    )
    grid_start = minimum(h.onset_days)
    grid_end = maximum(h.report_days)
    onsets = fill(50.0, grid_end)
    logit_h0 = fill(log(0.15 / 0.85), D)
    γ = zeros(grid_end - grid_start + 1)
    alpha = fill(0.6, grid_end - grid_start + 1)

    full = onset_report_moments(
        onsets, logit_h0, γ, grid_start, alpha,
        h.onset_days, h.report_days, h.prev_report_days
    )

    ## Added cells whose current and previous reads (a real predecessor,
    ## not the virtual-empty sentinel `0`) are both beyond the delay
    ## support have both reads saturated at the same asymptote, so their
    ## mean is exactly zero.
    added_idx = findall(
        i -> h.prev_report_days[i] > 0 &&
            h.report_days[i] - h.onset_days[i] > D - 1 &&
            h.prev_report_days[i] - h.onset_days[i] > D - 1,
        eachindex(h.onset_days)
    )
    @test !isempty(added_idx)
    @test all(==(0.0), full.means[added_idx])
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

@testitem "onset_report_cdf_extrapolated matches in-range, safe outside" begin
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

@testitem "onset_report_moments: a later snapshot sees more of one date" begin
    ## The literal right-truncation proof: for the same onset date and fixed
    ## hazards, the modelled current-cumulative level is non-decreasing as
    ## the report day moves later.
    using BVDOutbreakSize: onset_report_moments
    using StatsFuns: logit

    onsets = fill(10.0, 60)
    logit_h0 = fill(logit(0.1), 28)
    γ = zeros(60)
    alpha = fill(0.8, 60)
    u = 20
    early = onset_report_moments(
        onsets, logit_h0, γ, 1, alpha, [u], [u + 3],
        [0]
    )
    late = onset_report_moments(
        onsets, logit_h0, γ, 1, alpha, [u], [u + 10],
        [0]
    )
    @test late.level_cur[1] >= early.level_cur[1]
end

@testitem "onset_report_G reaches one at D-1, monotone under underflow" begin
    using BVDOutbreakSize: onset_report_G
    using Random: seed!

    seed!(20260803)
    grid_start = 1
    u = 10
    for _ in 1:20
        D = 28
        logit_h0 = randn(D) .* 1.5 .- 1.0
        γ = randn(60) .* 0.3
        vals = [
            onset_report_G(δ, logit_h0, γ, u, grid_start)
                for δ in (-3):(D + 5)
        ]
        @test all(==(0.0), vals[1:3])
        @test issorted(vals)
        @test all(v -> -1.0e-8 <= v <= 1 + 1.0e-8, vals)
        g_D1 = onset_report_G(D - 1, logit_h0, γ, u, grid_start)
        @test g_D1 ≈ 1.0 atol = 1.0e-8
    end

    ## Every hazard underflows to ≈ 0: numerator and denominator both
    ## underflow together, so the safe-rate guard must return a finite
    ## ratio rather than `0 / 0 = NaN`.
    logit_h0 = fill(-40.0, 28)
    γ = zeros(40)
    for δ in (0, 5, 27)
        @test isfinite(onset_report_G(δ, logit_h0, γ, 5, 1))
    end
end

@testitem "onset_report_F reaches alpha at D-1 and is zero for delta < 0" begin
    using BVDOutbreakSize: onset_report_F
    using Random: seed!

    seed!(20260803)
    D = 28
    grid_start = 1
    u = 7
    for _ in 1:10
        logit_h0 = randn(D) .* 1.2 .- 0.8
        γ = randn(40) .* 0.2
        α = rand()
        f_D1 = onset_report_F(D - 1, logit_h0, γ, u, grid_start, α)
        @test f_D1 ≈ α atol = 1.0e-8
        @test onset_report_F(-1, logit_h0, γ, u, grid_start, α) == 0.0
    end
end

@testitem "onset_report_ascertainment has a finite gradient at 0 or 1" begin
    using BVDOutbreakSize: onset_report_ascertainment
    using ForwardDiff: gradient

    ## `logit(0)` is `-Inf`, whose forward value survives but whose
    ## gradient is `NaN`, and a `NaN` there would spread to the whole
    ## log-density rather than to this stream alone.
    f(x) = onset_report_ascertainment([x[1]], 0.0, [0.0])[1]
    for anchor in (0.0, 1.0e-300, 0.15, 1.0)
        α = f([anchor])
        @test isfinite(α)
        @test 0.0 < α < 1.0
        @test all(isfinite, gradient(f, [anchor]))
    end

    ## An anchor inside the guarded range passes through untouched when the
    ## offset and the walk are both zero, so the guard costs nothing there.
    @test onset_report_ascertainment([0.15], 0.0, [0.0])[1] ≈ 0.15
end

@testitem "onset_confirmation_anchor averages ahead over the pipeline delay" begin
    using BVDOutbreakSize: onset_confirmation_anchor

    a = [0.1, 0.2, 0.3, 0.4, 0.5]
    ## A cohort is confirmed d days after onset with probability pmf[d+1],
    ## so its anchor averages the confirmation chance ahead of its onset.
    pmf = [0.5, 0.5]
    out = onset_confirmation_anchor(a, pmf)
    @test out[1:4] ≈ [0.15, 0.25, 0.35, 0.45]
    ## Past the end of `a` it is held at its last value.
    @test out[5] ≈ 0.5
    ## A constant is returned exactly, and an unnormalised pmf is normalised.
    @test onset_confirmation_anchor(fill(0.23, 6), [2.0, 1.0, 1.0]) ≈ fill(0.23, 6)
    ## No delay is the identity.
    @test onset_confirmation_anchor(a, [1.0]) ≈ a
end

@testitem "onset_report_expected_total stays in bounds" begin
    using BVDOutbreakSize: onset_report_expected_total
    using StatsFuns: logit

    onsets = fill(5.0, 100)
    logit_h0 = fill(logit(0.1), 28)
    for (grid_start, grid_end) in ((1, 100), (1, 10), (5, 40), (50, 55))
        nt = max(grid_end - grid_start + 1, 1)
        γ = zeros(nt)
        alpha = fill(0.3, nt)
        total = onset_report_expected_total(
            onsets, logit_h0, γ, grid_start,
            alpha, grid_end, grid_start
        )
        @test isfinite(total)
        @test total >= 0
    end
end

@testitem "onset_report_expected_total covers dates before grid_start" begin
    ## `expected_onset_reported_T` must sum the full `1:n` onset series
    ## like every other stream's `expected_*_T`, not just the triangle's
    ## own `grid_start:grid_end` window (see `onset_report_G`).
    using BVDOutbreakSize: onset_report_expected_total, onset_report_F
    using StatsFuns: logit

    n = 100
    onsets = fill(5.0, n)
    logit_h0 = fill(logit(0.2), 28)
    grid_start = 60
    grid_end = 90
    γ = zeros(grid_end - grid_start + 1)
    alpha = fill(0.4, grid_end - grid_start + 1)

    ## A sum restricted to just `grid_start:grid_end`, for comparison
    ## against the full total below.
    restricted = sum(
        onsets[u] * onset_report_F(
            grid_end - u, logit_h0, γ, u, grid_start,
            alpha[u - grid_start + 1]
        )
            for u in grid_start:grid_end
    )
    total = onset_report_expected_total(
        onsets, logit_h0, γ, grid_start,
        alpha, grid_end, grid_start
    )

    ## The onset dates before `grid_start` (days 1:59) are old enough by
    ## `grid_end` that they sit at the ascertainment level and each
    ## contribute onsets[u] * F(u, D-1) > 0, so the full total must exceed
    ## the window-restricted sum.
    @test total > restricted
    @test isfinite(total)

    ## The extrapolated contribution for a day well before `grid_start`
    ## should match the flat asymptote computed at the earliest known
    ## calendar day and ascertainment level (both held flat at index 1).
    F_edge = onset_report_F(
        length(logit_h0) - 1, logit_h0, γ, grid_start,
        grid_start, alpha[1]
    )
    @test total ≈ restricted + (grid_start - 1) * 5.0 * F_edge
end

@testitem "load_onset_curve: a narrower vintage drops the uncovered date" begin
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
    write(
        path, """
        sitrep,report_date,onset_date,confirmed_alive,confirmed_dead,confirmed_total
        001,2026-03-05,2026-03-01,4,0,4
        001,2026-03-05,2026-03-02,2,0,2
        001,2026-03-05,2026-03-03,1,0,1
        002,2026-03-07,2026-03-02,3,0,3
        002,2026-03-07,2026-03-03,2,0,2
        """
    )
    seeding = Date("2026-01-01")
    h = load_onset_curve(path; cutoff = Date("2026-03-07"), seeding)
    _day(x) = Int(date2epochdays(Date(x)) - date2epochdays(seeding)) + 1
    u = _day("2026-03-01")
    R2 = _day("2026-03-07")
    ## No cell for 03-01 in the 001-versus-002 pair.
    @test isempty(
        findall(
            i -> h.onset_days[i] == u && h.report_days[i] == R2,
            eachindex(h.onset_days)
        )
    )
    ## The first pair still scores 03-01 as a level against the empty
    ## predecessor, since block 001 does print it.
    R1 = _day("2026-03-05")
    first_idx = findall(
        i -> h.onset_days[i] == u && h.report_days[i] == R1,
        eachindex(h.onset_days)
    )
    @test length(first_idx) == 1
    @test h.increments[first_idx[1]] == 4
    ## Nothing anywhere fabricates the -4 the dropped cell would have given.
    @test minimum(h.increments) >= 0
end

@testitem "load_onset_curve: an onset past its figure's report is not a print" begin
    ## A figure can print an onset date later than its own report date. That
    ## bar has a negative delay, where `onset_report_cdf` is zero, so it is
    ## not a print of that date: the date's first print is the next figure,
    ## which scores it as a level.
    using BVDOutbreakSize: load_onset_curve
    using Dates: Date, date2epochdays

    dir = mktempdir()
    path = joinpath(dir, "onset.csv")
    ## Both blocks print 03-06, which postdates block 001's own report date.
    write(
        path, """
        sitrep,report_date,onset_date,confirmed_alive,confirmed_dead,confirmed_total
        001,2026-03-05,2026-03-04,4,0,4
        001,2026-03-05,2026-03-06,3,0,3
        002,2026-03-06,2026-03-04,6,0,6
        002,2026-03-06,2026-03-06,9,0,9
        """
    )
    seeding = Date("2026-01-01")
    h = load_onset_curve(path; cutoff = Date("2026-03-06"), seeding)
    _day(x) = Int(date2epochdays(Date(x)) - date2epochdays(seeding)) + 1
    ## 03-06 scores once, as a level at 002.
    i606 = findall(==(_day("2026-03-06")), h.onset_days)
    @test length(i606) == 1
    @test h.report_days[i606[1]] == _day("2026-03-06")
    @test h.prev_report_days[i606[1]] == 0
    @test h.increments[i606[1]] == 9
    ## The 03-04 correction the pair does support is untouched.
    idx = findall(
        i -> h.onset_days[i] == _day("2026-03-04") &&
            h.report_days[i] == _day("2026-03-06"),
        eachindex(h.onset_days)
    )
    @test length(idx) == 1
    @test h.increments[idx[1]] == 2
    ## No scored correction cell ever carries a negative previous delay.
    corr = findall(!=(0), h.prev_report_days)
    @test all(h.prev_report_days[i] >= h.onset_days[i] for i in corr)
end

@testitem "load_onset_curve: the archive scores no negative previous delay" begin
    ## Guards the committed CSV against the shape above: a vintage whose
    ## printed window runs past its predecessor's report date.
    using BVDOutbreakSize: BVDOutbreakSize, load_onset_curve
    using Dates: Date

    path = joinpath(
        pkgdir(BVDOutbreakSize), "data",
        "onset_curve_scanned.csv"
    )
    h = load_onset_curve(
        path; cutoff = Date("2100-01-01"),
        seeding = Date("2026-01-01")
    )
    corr = findall(!=(0), h.prev_report_days)
    @test !isempty(corr)
    @test all(h.prev_report_days[i] >= h.onset_days[i] for i in corr)
end

@testitem "onset_report_scales: Student-t variance matches count plus read" begin
    using BVDOutbreakSize: onset_report_scales
    using Distributions: TDist, var

    means = [0.0, 20.0, 40.0]
    ## Cells 1 and 3 are levels with one read; cell 2 is a correction
    ## between two snapshots, so it carries two reads.
    reads = [1, 2, 1]
    τ, k, ν = 1.2, 5.0, 4.0
    s = onset_report_scales(means, τ, k, reads, ν)
    target = means .+ means .^ 2 ./ k .+ reads .* τ^2
    ## A Student-t with scale `σ` has variance `σ² ν / (ν - 2)`.
    @test s .^ 2 .* var(TDist(ν)) ≈ target
    ## A level of 40 is dominated by its count variation, not the read.
    @test s[3] > s[1] * 4
    ## A negative mean cannot give a negative variance.
    @test isfinite(only(onset_report_scales([-5.0], τ, k, [2], ν)))
end

@testitem "onset_reporting_model: scale is count variation plus a read SD per read" begin
    ## A two-vintage triangle: the first vintage's cells and the date the
    ## second prints first are levels (one read); the rest are corrections
    ## (two reads).
    using BVDOutbreakSize: onset_reporting_model
    using Turing: DynamicPPL
    using Random: seed!

    oc = (;
        onset_days = [10, 11, 12, 13, 10, 11, 12, 13, 14],
        report_days = [15, 15, 15, 15, 20, 20, 20, 20, 20],
        prev_report_days = [0, 0, 0, 0, 15, 15, 15, 15, 0],
        increments = [2, 3, 1, 0, 1, 2, 3, 4, 5],
    )
    model = onset_reporting_model(oc, fill(30.0, 25))
    names = string.(collect(keys(DynamicPPL.VarInfo(model))))
    @test "τ" in names
    @test "inv_sqrt_k" in names
    seed!(20260924)
    out = model()
    @test out.τ > 0
    @test out.k > 0
    reads = [p == 0 ? 1 : 2 for p in oc.prev_report_days]
    μ = max.(out.modelled, 0)
    @test out.scales ≈
        sqrt.((out.ν - 2) / out.ν .* (μ .+ μ .^ 2 ./ out.k .+ reads .* out.τ^2))
end

@testitem "safe_studentt stays valid under extreme scale/df" begin
    using BVDOutbreakSize: safe_studentt
    using Distributions: mean, std, logpdf

    for (σ, ν) in (
            (0.0, 4.0), (-1.0, 4.0), (NaN, 4.0), (Inf, 4.0),
            (1.0, 0.0), (1.0, -2.0), (1.0, NaN),
        )
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

@testitem "onset_reporting_model conditions on increments, not sampled" begin
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

    base = (;
        onset_days = [10, 11, 12, 13, 10, 11, 12, 13, 14],
        report_days = [15, 15, 15, 15, 20, 20, 20, 20, 20],
        prev_report_days = [0, 0, 0, 0, 15, 15, 15, 15, 0],
    )
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
        DynamicPPL.LogDensityFunction(model, DynamicPPL.getlogjoint, vi), θ
    )
    lp_other = logdensity(
        DynamicPPL.LogDensityFunction(
            model_other, DynamicPPL.getlogjoint,
            vi
        ), θ
    )
    @test isfinite(lp)
    @test isfinite(lp_other)
    @test lp != lp_other
end

@testitem "onsets_only_model: prior predictive on a triangle is finite" begin
    using BVDOutbreakSize: onsets_only_model
    using Turing: Prior, sample

    oc = (;
        onset_days = [10, 11, 12, 13, 10, 11, 12, 13, 14],
        report_days = [15, 15, 15, 15, 20, 20, 20, 20, 20],
        prev_report_days = [0, 0, 0, 0, 15, 15, 15, 15, 0],
        increments = [2, 3, 1, 0, 1, 2, 3, 4, 5],
    )
    chn = sample(
        onsets_only_model(40; onset_curve_history = oc), Prior(),
        20; progress = false
    )
    et = vec(Array(chn[:expected_onset_reported_T]))
    @test length(et) == 20
    @test all(isfinite, et)
    @test all(>=(0), et)
end

@testitem "onsets_only_model: unanchored ascertainment median near 0.15" begin
    ## No confirmed pipeline to anchor on, so the ascertainment level's
    ## prior is exactly `logistic(logit(0.15) + β)` with `β ~ Normal(0,
    ## 0.75)` (the onset-axis walk is near-flat a priori): median ≈ 0.15,
    ## 90% interval roughly [0.05, 0.38].
    using BVDOutbreakSize: onsets_only_model
    using Turing: Prior, sample
    using Statistics: median, quantile

    oc = (;
        onset_days = [10, 11, 12, 13, 10, 11, 12, 13, 14],
        report_days = [15, 15, 15, 15, 20, 20, 20, 20, 20],
        prev_report_days = [0, 0, 0, 0, 15, 15, 15, 15, 0],
        increments = [2, 3, 1, 0, 1, 2, 3, 4, 5],
    )
    chn = sample(
        onsets_only_model(40; onset_curve_history = oc), Prior(),
        2000; progress = false
    )
    flat = reduce(vcat, vec(collect(chn[:onset_ascertainment])))
    @test isapprox(median(flat), 0.15; atol = 0.05)
    @test 0.02 < quantile(flat, 0.05) < 0.1
    @test 0.25 < quantile(flat, 0.95) < 0.55
end

@testitem "AD gradient: onsets_only_model differentiates (Mooncake)" tags = [
    :ad,
] begin
    using Turing: DynamicPPL
    using LogDensityProblems: logdensity_and_gradient
    using Random: seed!
    using BVDOutbreakSize: onsets_only_model, default_adtype

    oc = (;
        onset_days = [10, 11, 12, 13, 10, 11, 12, 13, 14],
        report_days = [15, 15, 15, 15, 20, 20, 20, 20, 20],
        prev_report_days = [0, 0, 0, 0, 15, 15, 15, 15, 0],
        increments = [2, 3, 1, 0, 1, 2, 3, 4, 5],
    )

    seed!(20260518)
    model = onsets_only_model(40; onset_curve_history = oc)
    vi = DynamicPPL.link(DynamicPPL.VarInfo(model), model)
    x0 = collect(vi[:])
    ldf = DynamicPPL.LogDensityFunction(
        model, DynamicPPL.getlogjoint, vi; adtype = default_adtype()
    )
    logp, grad = logdensity_and_gradient(ldf, x0)
    @test isfinite(logp)
    @test length(grad) == length(x0)
    @test all(isfinite, grad)
    @test any(!iszero, grad)
end

@testitem "bvd_joint: the onset stream is wired in" begin
    using BVDOutbreakSize: bvd_joint
    using Turing: sample, Prior

    n = 40
    dh = (; days = [13, 18, 40], counts = [10, 14, 18])
    rh = (; days = [13, 18, 40], counts = [340, 516, 905])
    ch = (; days = [13, 18, 40], counts = [9, 17, 27])
    oc = (;
        onset_days = [20, 21, 22, 23, 20, 21, 22, 23, 24],
        report_days = [25, 25, 25, 25, 30, 30, 30, 30, 30],
        prev_report_days = [0, 0, 0, 0, 25, 25, 25, 25, 0],
        increments = [3, 2, 1, 0, 1, 2, 1, 3, 2],
    )
    chn = sample(
        bvd_joint(
            n, 2, 18, 905, 0, 27, 50;
            confirmed_deaths = 5,
            deaths_history = dh,
            reported_history = rh,
            confirmed_history = ch,
            lab_history = (; days = [18, 40], counts = [30, 50]),
            onset_curve_history = oc,
            breakpoint = 30
        ),
        Prior(), 12; progress = false
    )
    C_T = vec(Array(chn[:C_T]))
    et = vec(Array(chn[:expected_onset_reported_T]))
    @test length(C_T) == 12
    @test all(isfinite, C_T)
    @test all(C_T .> 0)
    @test length(et) == 12
    @test all(isfinite, et)
end

## --- Per-vintage totals ------------------------------------------------

@testitem "load_onset_curve reports each vintage's printed total" begin
    using BVDOutbreakSize: load_onset_curve
    using Dates: Date

    dir = mktempdir()
    path = joinpath(dir, "onset.csv")
    ## Three distinct vintages. The third reads a smaller total than the
    ## second, which late reporting cannot produce and the per-scan level
    ## error can: the totals are recorded as read rather than made
    ## monotone.
    write(
        path, """
        sitrep,report_date,onset_date,confirmed_alive,confirmed_dead,confirmed_total
        001,2026-03-05,2026-03-01,2,0,2
        001,2026-03-05,2026-03-02,1,0,1
        002,2026-03-07,2026-03-01,4,0,4
        002,2026-03-07,2026-03-02,3,0,3
        002,2026-03-07,2026-03-03,2,0,2
        003,2026-03-09,2026-03-01,4,0,4
        003,2026-03-09,2026-03-02,2,0,2
        003,2026-03-09,2026-03-03,2,0,2
        """
    )
    oc = load_onset_curve(
        path; cutoff = Date("2026-03-20"),
        seeding = Date("2026-03-01")
    )
    ## Seeding day is grid day 1, so 5/7/9 March are grid days 5/7/9.
    @test oc.total_days == [5, 7, 9]
    @test oc.total_counts == [3, 9, 8]
    @test oc.last_total == 8
end

@testitem "load_onset_curve totals honour the cut-off and the no-op path" begin
    using BVDOutbreakSize: load_onset_curve
    using Dates: Date

    dir = mktempdir()
    path = joinpath(dir, "onset.csv")
    write(
        path, """
        sitrep,report_date,onset_date,confirmed_alive,confirmed_dead,confirmed_total
        001,2026-03-05,2026-03-01,2,0,2
        002,2026-03-09,2026-03-01,5,0,5
        """
    )
    ## The 9 March vintage is past the cut-off, so neither its cells nor its
    ## total survive.
    oc = load_onset_curve(
        path; cutoff = Date("2026-03-06"),
        seeding = Date("2026-03-01")
    )
    @test oc.total_days == [5]
    @test oc.total_counts == [2]
    @test oc.last_total == 2

    noop = load_onset_curve(
        joinpath(dir, "absent.csv");
        cutoff = Date("2026-03-06"), seeding = Date("2026-03-01")
    )
    @test isempty(noop.total_days)
    @test isempty(noop.total_counts)
    @test ismissing(noop.last_total)
end

## --- Hazard reconstruction and the onset nowcast/forecast ---------------

@testitem "fitted_onset_hazard reads the model's own hazard" begin
    ## The hazard is the fitted model's own state at each draw, so feeding it
    ## back into `onset_report_expected_total` with the chain's onset
    ## trajectory reproduces the `expected_onset_reported_T` the model
    ## tracked for that draw, and `alpha` is the tracked ascertainment.
    using BVDOutbreakSize: onsets_only_model, fitted_onset_hazard,
        onset_report_expected_total
    using Turing: Prior, sample
    import FlexiChains

    oc = (;
        onset_days = [10, 11, 12, 13, 10, 11, 12, 13, 14],
        report_days = [15, 15, 15, 15, 20, 20, 20, 20, 20],
        prev_report_days = [0, 0, 0, 0, 15, 15, 15, 15, 0],
        increments = [2, 3, 1, 0, 1, 2, 3, 4, 5],
    )
    n = 40
    m = onsets_only_model(n; onset_curve_history = oc)
    chn = sample(
        m, Prior(), 20; chain_type = FlexiChains.VNChain, progress = false
    )
    grid_start = minimum(oc.onset_days)
    grid_end = maximum(oc.report_days)
    hz = fitted_onset_hazard(m, chn)
    daily = [
        (v = collect(t); vcat(v[1], diff(v)))
            for t in vec(collect(chn[:cumulative_onsets]))
    ]
    et = vec(Array(chn[:expected_onset_reported_T]))
    @test length(hz.logit_h0) == 20
    @test all(length(g) == grid_end - grid_start + 1 for g in hz.γ)
    @test hz.alpha == [collect(a) for a in vec(collect(chn[:onset_ascertainment]))]
    rebuilt = [
        onset_report_expected_total(
            daily[i], hz.logit_h0[i], hz.γ[i], grid_start, hz.alpha[i], n,
            grid_start
        ) for i in 1:20
    ]
    @test rebuilt ≈ et
end

@testitem "forecast_archive carries the onset reporting increment" begin
    using DataFrames: DataFrame
    using BVDOutbreakSize: forecast_archive
    using Dates: Date

    fc = DataFrame(
        onset_reports_new = [10, 12, 14, 9],
        confirmed_new = [3, 4, 5, 6]
    )
    arc = forecast_archive([(7, fc)]; made_date = Date("2026-07-25"))
    @test "onset reports" in arc.stream
    rows = arc[arc.stream .== "onset reports", :]
    @test size(rows, 1) == 4
    @test rows.value == [10.0, 12.0, 14.0, 9.0]
    @test all(rows.target_date .== Date("2026-08-01"))
end

@testitem "summed Student-t likelihood matches one term per cell" begin
    using BVDOutbreakSize: safe_studentt, StudentTVector,
        onset_increments_model
    using Distributions: logpdf
    using Random: Xoshiro
    using Turing: DynamicPPL
    using Turing.DynamicPPL: @varname

    loglik(m) = DynamicPPL.loglikelihood(m, DynamicPPL.VarInfo(Xoshiro(1), m))
    ## A floored scale (zero, non-finite) and negative increments included.
    μ = [12.0, -3.5, 40.5, 7.2, 0.0, 9.0]
    σ = [3.0, 2.5, 0.0, 6.1, NaN, 1.0]
    x = [10, -6, 44, 3, 0, 20]
    per_term(ν) = sum(logpdf(safe_studentt(μ[i], σ[i], ν), x[i]) for i in 1:6)
    ## A defaulted `ν` falls back to 4 in both.
    for ν in (4.0, 1.5, -1.0)
        @test logpdf(StudentTVector(μ, σ, ν), x) ≈ per_term(ν)
        @test loglik(onset_increments_model(μ, σ, x, ν)) ≈ per_term(ν)
    end
    @test logpdf(StudentTVector(μ, σ, -1.0), x) ==
        logpdf(StudentTVector(μ, σ, 4.0), x)
    @test loglik(onset_increments_model(Float64[], Float64[], Int[], 4.0)) == 0

    ## A `missing` vector samples as one variable under the whole-vector
    ## key the predictive path reads. It is unconstrained, so linking it
    ## leaves the log density unchanged.
    m = onset_increments_model(μ, σ, missing, 4.0)
    vi = DynamicPPL.VarInfo(Xoshiro(1), m)
    @test Set(keys(vi)) == Set([@varname(increments)])
    @test DynamicPPL.getlogjoint(DynamicPPL.link(vi, m)) ≈
        DynamicPPL.getlogjoint(vi)
    draw = m(Xoshiro(2)).increments
    @test draw isa Vector{Float64} && length(draw) == 6
end

@testitem "load_onset_curve: the archive scores each onset date once" begin
    ## The seeding day sits before the genetic TMRCA bound and the earliest
    ## digitised onset date is in late April, so every scored cell sits
    ## inside the 1-based grid. Every printed onset date has exactly one
    ## level, its corrections sit inside the delay support, and a date's
    ## cells sum to the last print they reach.
    using BVDOutbreakSize: BVDOutbreakSize, load_observations,
        ONSET_REPORT_MAX_DELAY
    using Dates: Day

    obs = load_observations()
    h = obs.onset_curve_history
    @test !isempty(h.onset_days)
    @test minimum(h.onset_days) >= 1

    path = joinpath(
        pkgdir(BVDOutbreakSize), "data",
        "onset_curve_scanned.csv"
    )
    blocks = filter(
        b -> b.report_date <= obs.cutoff,
        BVDOutbreakSize._dedup_onset_blocks(
            BVDOutbreakSize._read_onset_curve_blocks(path)
        )
    )
    _date(u) = obs.cutoff - Day(obs.n - u)
    levels = h.prev_report_days .== 0
    @test allunique(h.onset_days[levels])
    @test Set(h.onset_days) == Set(h.onset_days[levels])
    corr = .!levels
    @test all(h.report_days[corr] .- h.onset_days[corr] .< ONSET_REPORT_MAX_DELAY)
    by_report = Dict(b.report_date => b for b in blocks)
    for u in unique(h.onset_days)
        idx = findall(==(u), h.onset_days)
        last_R = maximum(h.report_days[idx])
        printed = get(by_report[_date(last_R)].onsets, _date(u), 0)
        @test sum(h.increments[idx]) == printed
    end
end

@testitem "onset_reporting_model: the calendar walk and alpha start on their own grids" begin
    ## Report days starting more than D days after the first onset day, so
    ## the hazard's own grid start (`onset_hazard_grid_start`) diverges from
    ## `minimum(onset_days)`, the grid `alpha` is indexed from. The model's
    ## own state carries both origins, `fitted_onset_hazard` returns `γ` and
    ## `alpha` at their own lengths, and the tracked cut-off total is the
    ## 7-argument `onset_report_expected_total` on those two grids.
    using BVDOutbreakSize: onsets_only_model, fitted_onset_hazard,
        onset_hazard_grid_start, onset_report_expected_total,
        ONSET_REPORT_MAX_DELAY
    using Turing: Prior, sample, returned
    import FlexiChains

    oc = (;
        onset_days = [1, 2, 3, 4, 1, 2, 3, 4, 5],
        report_days = [50, 50, 50, 50, 55, 55, 55, 55, 55],
        prev_report_days = [0, 0, 0, 0, 50, 50, 50, 50, 0],
        increments = [2, 3, 1, 0, 1, 2, 3, 4, 5],
    )
    n = 80
    D = ONSET_REPORT_MAX_DELAY
    alpha_grid_start = minimum(oc.onset_days)
    grid_end = maximum(oc.report_days)
    grid_start = onset_hazard_grid_start(oc.onset_days, oc.report_days; D)
    @test grid_start > alpha_grid_start

    m = onsets_only_model(n; onset_curve_history = oc, breakpoint = 30)
    chn = sample(
        m, Prior(), 5; chain_type = FlexiChains.VNChain, progress = false
    )
    states = [r.onset_report_state for r in vec(returned(m, chn))]
    @test all(st.grid_start == grid_start for st in states)
    @test all(st.alpha_grid_start == alpha_grid_start for st in states)

    hz = fitted_onset_hazard(m, chn)
    @test all(length(g) == grid_end - grid_start + 1 for g in hz.γ)
    @test all(length(a) == grid_end - alpha_grid_start + 1 for a in hz.alpha)

    daily = [
        (v = collect(t); vcat(v[1], diff(v)))
            for t in vec(collect(chn[:cumulative_onsets]))
    ]
    et = vec(Array(chn[:expected_onset_reported_T]))
    rebuilt = [
        onset_report_expected_total(
            daily[i], hz.logit_h0[i], hz.γ[i], grid_start, hz.alpha[i], n,
            alpha_grid_start
        ) for i in 1:5
    ]
    @test rebuilt ≈ et
end
