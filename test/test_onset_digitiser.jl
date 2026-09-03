## Tests for the onset-curve digitiser and the file it writes
## (`scripts/digitize_onset_curve.jl`, its Python twin, and
## `data/onset_curve_scanned.csv`).
##
## Three layers, in increasing order of what they need to run: invariants of
## the committed CSV and the L1 date-alignment sweep, which need nothing but
## the repository; synthetic-figure tests of `digitize`, which build a bar
## chart with known counts in memory; and an end-to-end re-digitisation of
## the real figures, which needs `data/sitrep_pdfs` (git-ignored, fetched by
## `scripts/download_sitreps.jl`) and skips itself when that is absent.
##
## Filter target for a scoped run (this worktree only, not every sibling
## worktree `@run_package_tests` would otherwise discover):
##   target = joinpath(pwd(), "test", "test_onset_digitiser.jl")
##   @run_package_tests filter = ti -> string(ti.filename) == target

## --- Committed CSV: schema and structural invariants ----------------------

@testitem "onset curve CSV rows are internally consistent" begin
    using Dates: Date, Day

    path = joinpath(pkgdir(BVDOutbreakSize), "data",
        "onset_curve_scanned.csv")
    lines = readlines(path)
    @test lines[1] ==
          "sitrep,report_date,onset_date,confirmed_alive," *
          "confirmed_dead,confirmed_total"

    ## Violations are collected rather than asserted row by row: the file
    ## carries four thousand rows, and a failure should name the offending
    ## ones instead of burying them in four thousand passes.
    malformed = String[]
    mismatched = String[]
    negative = String[]
    duplicated = String[]
    late = String[]
    seen = Set{Tuple{String, String}}()
    order = String[]
    for line in lines[2:end]
        isempty(strip(line)) && continue
        f = split(line, ',')
        if length(f) != 6
            push!(malformed, line)
            continue
        end
        alive = parse(Int, f[4])
        dead = parse(Int, f[5])
        total = parse(Int, f[6])
        ## The total column is derived, so it must never drift from its
        ## parts.
        total == alive + dead || push!(mismatched, line)
        (alive >= 0 && dead >= 0) || push!(negative, line)
        key = (String(f[1]), String(f[3]))
        ## One row per vintage per onset date.
        key in seen && push!(duplicated, line)
        push!(seen, key)
        (isempty(order) || order[end] != f[1]) && push!(order, String(f[1]))
        ## An onset date can never sit later than the axis of the report
        ## that draws it, and that axis runs at most a day past the
        ## rapportage date.
        Date(f[3]) <= Date(f[2]) + Day(1) || push!(late, line)
    end
    @test isempty(malformed)
    @test isempty(mismatched)
    @test isempty(negative)
    @test isempty(duplicated)
    @test isempty(late)
    ## Blocks are contiguous and in ascending vintage order.
    @test length(unique(order)) == length(order)
    @test issorted(order)
end

@testitem "SitRep 098 stays out of the digitised onset curve" begin
    ## 098 is the only vintage whose onset figure INSP embedded losslessly
    ## rather than as JPEG. The digitiser's colour masks are fixed
    ## thresholds, so JPEG edge blur costs every other vintage a fringe of
    ## every bar while 098 keeps it, and 098 reads about 7% higher on the
    ## same underlying data. The stream is consumed as between-vintage
    ## increments, so splicing in one vintage measured on a different bias
    ## scale would post a 284-case fall on onset dates that cannot fall.
    ## The evidence is in data/README.md. Issue #636 tracks the same blur
    ## costing every JPEG vintage a smaller, unmeasured slice of every bar.
    path = joinpath(pkgdir(BVDOutbreakSize), "data",
        "onset_curve_scanned.csv")
    rows = filter(!isempty, strip.(readlines(path)[2:end]))
    ids = Set(String(split(l, ',')[1]) for l in rows)
    @test !("098" in ids)
    ## Its neighbours are digitised, so this is a deliberate exclusion and
    ## not a hole in the downloaded reports.
    @test "097" in ids && "099" in ids && "100" in ids
end

@testitem "onset curve reprints collapse to identical blocks" begin
    include(joinpath(@__DIR__, "onset_digitiser_helpers.jl"))

    csv = _read_onset_csv(joinpath(pkgdir(BVDOutbreakSize), "data",
        "onset_curve_scanned.csv"))
    ## Reprints are collapsed by exact value equality over the digitised
    ## block, never by a list of vintage ids, so a pair that reprints must
    ## agree on every onset date. 108/109 is the newest such pair: INSP
    ## reprinted 108's figure unchanged, printed n included.
    @test csv.onsets["108"] == csv.onsets["109"]
    @test csv.onsets["104"] == csv.onsets["105"]
    @test csv.onsets["106"] == csv.onsets["107"]
    ## The vintage before a reprint pair is genuinely fresh content.
    @test csv.onsets["107"] != csv.onsets["108"]
end

## --- Committed CSV: L1 date alignment -------------------------------------

@testitem "onset curve L1 alignment lands on shift 0 bar the documented pairs" begin
    using Dates: Day
    include(joinpath(@__DIR__, "onset_digitiser_helpers.jl"))

    csv = _read_onset_csv(joinpath(pkgdir(BVDOutbreakSize), "data",
        "onset_curve_scanned.csv"))

    ## Consecutive pairs where the check does not land on shift 0. Every one
    ## is investigated and recorded in data/README.md's onset-curve section:
    ## the rightmost axis tick was read off the rendered figure directly,
    ## and a verified direct read outweighs this heuristic when the two
    ## conflict. They disagree in alternating directions rather than showing
    ## a systematic offset, and each is bracketed by pairs that do land on
    ## 0, so none of them can be a misread tick.
    documented = Dict("093" => "094", "096" => "097", "099" => "100",
        "102" => "103")

    unexpected = Tuple{String, String, Int, Int}[]
    resolved = String[]
    for i in 2:length(csv.order)
        p, q = csv.order[i - 1], csv.order[i]
        a, b = csv.onsets[p], csv.onsets[q]
        cut = csv.report_date[q] - Day(21)
        full = [_l1_shift(a, b, s) for s in (-1, 0, 1)]
        stable = [_l1_shift(a, b, s; cut = cut) for s in (-1, 0, 1)]
        best_full = (-1, 0, 1)[argmin(full)]
        best_stable = (-1, 0, 1)[argmin(stable)]
        clean = best_full == 0 && best_stable == 0
        if get(documented, p, nothing) == q
            clean && push!(resolved, "$p->$q")
        elseif !clean
            push!(unexpected, (p, q, best_full, best_stable))
        end
    end
    ## A new pair failing the check means that vintage's rightmost tick
    ## needs reading off the figure and recording before its block is
    ## accepted.
    @test isempty(unexpected)
    ## A documented pair that starts passing means the note describing it is
    ## stale.
    @test isempty(resolved)
end

## --- digitize on a synthetic figure ---------------------------------------

@testitem "digitize recovers a synthetic chart's counts and dates" begin
    using Dates: Date, Day
    include(joinpath(@__DIR__, "onset_digitiser_helpers.jl"))
    include(joinpath(pkgdir(BVDOutbreakSize), "scripts",
        "digitize_onset_curve.jl"))

    bars = [(12, 5) for _ in 1:15]
    R, G, B = _synthetic_chart(bars)
    last_tick = Date(2026, 8, 24)
    drawn = [r for r in digitize(R, G, B, last_tick, 20) if r[2] + r[3] > 0]

    @test length(drawn) == length(bars)
    ## The rightmost bar sits on the rightmost weekly tick, and the rest
    ## step back a day each.
    @test [r[1] for r in drawn] ==
          [last_tick - Day(length(bars) - j) for j in eachindex(bars)]
    ## Bar heights come back through the y-axis scale unchanged, split by
    ## colour into the alive and dead segments.
    @test all(r -> (r[2], r[3]) == (12, 5), drawn)
end

@testitem "digitize reads the count scale off the y-axis gridline step" begin
    using Dates: Date
    include(joinpath(@__DIR__, "onset_digitiser_helpers.jl"))
    include(joinpath(pkgdir(BVDOutbreakSize), "scripts",
        "digitize_onset_curve.jl"))

    ## The same drawn pixels under two printed gridline increments. The
    ## geometry cannot tell them apart, which is why Y_AXIS_STEP is a
    ## per-vintage override read off the figure rather than inferred: the
    ## 0/25/50/75 grid the figures switched to at SitRep 087 carries a
    ## quarter more counts per pixel than the 0/20/40/60 grid before it.
    R, G, B = _synthetic_chart([(8, 4) for _ in 1:10])
    nz(rows) = [r for r in rows if r[2] + r[3] > 0]
    twenty = nz(digitize(R, G, B, Date(2026, 8, 24), 20))
    twentyfive = nz(digitize(R, G, B, Date(2026, 8, 24), 25))
    @test !isempty(twenty) && length(twenty) == length(twentyfive)
    @test all(r -> (r[2], r[3]) == (8, 4), twenty)
    @test all(r -> (r[2], r[3]) == (10, 5), twentyfive)
end

## --- End to end against the real figures ----------------------------------

@testitem "digitiser reproduces the committed onset CSV from the SitRep PDFs" begin
    using Dates: Date

    root = pkgdir(BVDOutbreakSize)
    pdf_dir = joinpath(root, "data", "sitrep_pdfs")

    if !isdir(pdf_dir)
        @info "onset re-digitisation skipped: $pdf_dir absent. Fetch the " *
              "reports with scripts/download_sitreps.jl to run it."
        @test true
    else
        include(joinpath(root, "scripts", "digitize_onset_curve.jl"))
        want = Dict{String, Vector{Tuple{Date, Int, Int}}}()
        for line in readlines(joinpath(root, "data",
            "onset_curve_scanned.csv"))[2:end]
            isempty(strip(line)) && continue
            f = split(line, ',')
            push!(get!(want, String(f[1]), Tuple{Date, Int, Int}[]),
                (Date(f[3]), parse(Int, f[4]), parse(Int, f[5])))
        end
        ## Three vintages, one per detection path the digitiser has had to
        ## grow: 094 exercises the y-axis label strip, 106 the daily bar
        ## windows, and 108 the near-grey fallback its larger render needs.
        wanted = (("094", Date(2026, 8, 17)), ("106", Date(2026, 8, 24)),
            ("108", Date(2026, 8, 31)))
        available = [(sr, tick, joinpath(pdf_dir,
                         "SitRep_MVE_$(sr)_2026.pdf"))
                     for (sr, tick) in wanted
                     if isfile(joinpath(pdf_dir,
            "SitRep_MVE_$(sr)_2026.pdf"))]
        isempty(available) &&
            @info "no onset SitRep PDFs matched; nothing re-digitised"
        for (sr, last_tick, pdf) in available
            img = onset_image(pdf)
            @test img !== nothing
            img === nothing && continue
            @test digitize(img..., last_tick, get(Y_AXIS_STEP, sr, 20)) ==
                  want[sr]
        end
    end
end

@testitem "the Python port reproduces the committed onset CSV" begin
    include(joinpath(@__DIR__, "onset_digitiser_helpers.jl"))

    ## The port is the script the automated data-updater runs, so a drift in
    ## it rewrites committed rows without anyone reading Julia output. The
    ## check has to go through the real figures. A synthetic chart cannot
    ## substitute: the daily bar height is a 75th percentile over a window
    ## of four or five columns, which is deliberately robust to losing one
    ## column at the edge, so on flat drawn bars the two implementations
    ## agree even when their windows differ. Only the partial-height
    ## columns that JPEG anti-aliasing leaves on a real bar edge make the
    ## percentile move, which is why the divergence showed up on SitRep 106
    ## and not on any drawn figure.
    root = pkgdir(BVDOutbreakSize)
    pdf_dir = joinpath(root, "data", "sitrep_pdfs")
    runner = _python_runner()

    if !isdir(pdf_dir) || runner === nothing
        @info "Python parity check skipped" pdfs = isdir(pdf_dir) runner
        @test true
    else
        out = joinpath(mktempdir(), "onset.csv")
        script = joinpath(root, "scripts", "digitize_onset_curve.py")
        run(pipeline(`$runner $script $pdf_dir $out`;
            stdout = devnull, stderr = devnull))
        committed = joinpath(root, "data", "onset_curve_scanned.csv")
        @test read(out, String) == read(committed, String)
    end
end
