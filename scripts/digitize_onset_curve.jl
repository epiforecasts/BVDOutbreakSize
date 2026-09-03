#!/usr/bin/env julia
#
# Digitise the "courbe epidemique par date de debut des symptomes (liste
# lineaire DHIS2)" figure that the INSP analytique-format SitReps carry from
# SitRep 059 onward. That figure is the only published source for confirmed
# cases by symptom-onset date; it is a raster bar chart with no data table,
# so the daily counts are recovered from the figure pixels.
#
# This is the Julia reference implementation of the digitiser. It has no
# image-library dependency: poppler's `pdfimages` renders the embedded
# figure to an uncompressed PPM (P6, its default for RGB) that Base Julia
# parses directly. The
# equivalent `scripts/digitize_onset_curve.py` exists for the automated
# data-updater, which has Python (not Julia) access; both produce the same
# `data/onset_curve_scanned.csv`.
#
# Method (per figure, all self-calibrated from the image):
#   * baseline (count 0) = the widest dark horizontal row in the lower panel;
#   * count scale = the y-axis tick marks (0/20/40/60), evenly spaced, giving
#     pixels-per-count = tick-spacing / 20;
#   * date scale = the weekly x-axis tick marks, anchored on the rightmost
#     tick (whose date is in CONFIG, read off the axis) stepping back 7 days;
#   * each daily bar height = the 75th-percentile column in a one-day window,
#     flooded up from the baseline counting light-blue (Vivant) and crimson
#     (Decede) pixels, bridging the anti-alias gap between the stacked
#     segments but stopping at the wide white gap up to the floating label /
#     dashed line above the bar.
#
# Accuracy: the error is a few percent in either direction, per scan, and it is
# independent between vintages. Against the printed `n` it ranges from -3.0%
# (SitRep 069/070/071: 2260 vs n=2 329) to +1.6% (SitRep 068: 2344 vs n=2 308),
# with SitRep 064 at -2.2% (2018 vs n=2 064) and 072 at +0.4% (2531 vs n=2 521).
# Individual bars carry roughly +/-1-2 cases of pixel noise. Part of the
# shortfall sits in the faded bars of the `donnees potentiellement incompletes`
# band, whose lightened fill falls outside the colour masks, but that mechanism
# is one-sided and does not explain the overshoots, so treat the sign as
# unknown.
#
# The consequence that matters: the scans do not preserve a property the
# underlying data has. Late reporting only ever adds cases, so an onset date's
# count must be non-decreasing across vintages, yet on onset dates more than
# three weeks before the earliest report date in the file (12 July, so onsets
# before 21 June) the scanned totals move both ways between consecutive
# snapshots - 064 -> 065 falls by a net 36 cases across 34 of 54 such days,
# and every other consecutive pair falls somewhere too. So a
# between-vintage increment of a few cases is at or below the noise floor, and
# anything built on those increments (a reporting-delay estimate, say) has to
# account for it.
#
# The values are approximate and are not fitted by the model; they are captured
# for later use (see data/README.md and #488).
#
# Dependencies: poppler (`pdfimages`, `pdftotext`, `pdfinfo`) on PATH. No
# Julia packages beyond stdlib.
#
# Usage:
#   julia scripts/digitize_onset_curve.jl [pdf_dir] [out_csv]
# Defaults: pdf_dir = data/sitrep_pdfs, out_csv = data/onset_curve_scanned.csv
# Download the PDFs first with scripts/download_sitreps.jl.

using Dates: Date, Day, value
using Statistics: median, quantile
using Printf: @printf

# Per-vintage anchors. `report_date` is the SitRep rapportage date;
# `last_tick` is the date of the rightmost weekly x-axis tick, read off the
# figure (the axis range differs between vintages). To add a new vintage,
# append its SitRep number, rapportage date and last x-axis tick date.
const CONFIG = [
    ("059", Date(2026, 7, 12), Date(2026, 7, 12)),
    ("060", Date(2026, 7, 13), Date(2026, 7, 12)),
    ("061", Date(2026, 7, 14), Date(2026, 7, 12)),
    ("062", Date(2026, 7, 15), Date(2026, 7, 12)),
    ("064", Date(2026, 7, 17), Date(2026, 7, 15)),
    ("065", Date(2026, 7, 18), Date(2026, 7, 15)),
    ("066", Date(2026, 7, 19), Date(2026, 7, 15)),
    ("067", Date(2026, 7, 20), Date(2026, 7, 15)),
    ("068", Date(2026, 7, 21), Date(2026, 7, 22)),
    ("069", Date(2026, 7, 22), Date(2026, 7, 22)),
    ("070", Date(2026, 7, 23), Date(2026, 7, 22)),
    ("071", Date(2026, 7, 24), Date(2026, 7, 22)),
    ("072", Date(2026, 7, 25), Date(2026, 7, 22)),
    ("073", Date(2026, 7, 26), Date(2026, 7, 22)),
    ("074", Date(2026, 7, 27), Date(2026, 7, 22)),
    ("077", Date(2026, 7, 30), Date(2026, 7, 29)),
    ("078", Date(2026, 7, 31), Date(2026, 7, 29)),
    ("079", Date(2026, 8, 1), Date(2026, 7, 29)),
    ("080", Date(2026, 8, 2), Date(2026, 7, 29)),
    ("081", Date(2026, 8, 3), Date(2026, 7, 29)),
    ("082", Date(2026, 8, 4), Date(2026, 8, 5)),
    ("083", Date(2026, 8, 5), Date(2026, 8, 5)),
    ("087", Date(2026, 8, 9), Date(2026, 8, 5)),
    ("088", Date(2026, 8, 10), Date(2026, 8, 5)),
    ("089", Date(2026, 8, 11), Date(2026, 8, 5)),
    ("090", Date(2026, 8, 12), Date(2026, 8, 5)),
    ("091", Date(2026, 8, 13), Date(2026, 8, 5)),
    ("092", Date(2026, 8, 14), Date(2026, 8, 10)),
    ("093", Date(2026, 8, 15), Date(2026, 8, 10)),
    ("094", Date(2026, 8, 16), Date(2026, 8, 17)),
    ("095", Date(2026, 8, 17), Date(2026, 8, 17)),
    ("096", Date(2026, 8, 18), Date(2026, 8, 17)),
    ("097", Date(2026, 8, 19), Date(2026, 8, 17)),
    # "098" is deliberately absent. It is the only vintage whose figure
    # INSP embedded losslessly rather than as JPEG (1267x789, against
    # ~830x510 for every neighbour). The colour masks below are fixed
    # thresholds, and JPEG blur pushes a fringe of every bar's edge pixels
    # outside them, so every other vintage loses a slice of every bar and
    # 098 keeps it. It digitises to 4280 against a printed n of 4 140
    # (+3.4%, above the +1.6% ceiling measured over 059-083), and on onset
    # dates more than four weeks old - where late reporting can only add a
    # case or two - it sits ~280 cases above BOTH neighbours (stable-region
    # totals: 097 2875, 098 3145, 099 2867), so 098 -> 099 would post a
    # 284-case fall the underlying data cannot have. The neighbour scatter
    # is ~25 cases, an order of magnitude smaller.
    #
    # The render size is not the cause: resampling the neighbours up to
    # 1267x789 without softening an edge moves them by 1-3%, and resampling
    # 098 down to a neighbour's size leaves it at 4287. A fix would have to
    # make the masks encoding-insensitive, which reclassifies pixels on
    # every vintage and so cannot leave 059-100 unchanged. See
    # data/README.md.
    ("099", Date(2026, 8, 21), Date(2026, 8, 17)),
    ("100", Date(2026, 8, 22), Date(2026, 8, 17)),
    ("101", Date(2026, 8, 23), Date(2026, 8, 24)),
    ("102", Date(2026, 8, 24), Date(2026, 8, 24)),
    ("103", Date(2026, 8, 25), Date(2026, 8, 24)),
    ("104", Date(2026, 8, 26), Date(2026, 8, 24)),
    ("105", Date(2026, 8, 27), Date(2026, 8, 24)),
    ("106", Date(2026, 8, 28), Date(2026, 8, 24)),
    ("107", Date(2026, 8, 29), Date(2026, 8, 24)),
    ("108", Date(2026, 8, 30), Date(2026, 8, 31)),
    ("109", Date(2026, 8, 31), Date(2026, 8, 31))
]

# Every figure through SitRep 083 draws its y-axis on a 0/20/40/60/80 grid,
# which `digitize` assumed as a hard-coded divisor. From SitRep 087 the
# brief-format figure switched to a 0/25/50/75 grid (confirmed by reading
# the printed tick labels directly - the pixel geometry is otherwise
# indistinguishable, so this cannot be self-calibrated any more than
# `last_tick` can). Applying the old /20 divisor to a 25-count grid
# undercounts every bar by a scale-dependent amount and was caught only
# because it made stable, weeks-old onset dates fall (SitRep 083's 15 May
# read 26; the same date misread through the old divisor came out as 8).
# Override per vintage here; anything absent keeps the historical 20.
const Y_AXIS_STEP = Dict(
    "087" => 25,
    "088" => 25,
    "089" => 25,
    "090" => 25,
    "091" => 25,
    "092" => 25,
    "093" => 25,
    "094" => 25,
    "095" => 25,
    "096" => 25,
    "097" => 25,
    "099" => 25,
    "100" => 25,
    "101" => 25,
    "102" => 25,
    "103" => 25,
    "104" => 25,
    "105" => 25,
    "106" => 25,
    "107" => 25,
    "108" => 25,
    "109" => 25
)

# --- PPM (P6) reader ------------------------------------------------------
# Returns R, G, B as Int matrices indexed [row, col].
function read_ppm(path)
    bytes = read(path)
    @assert bytes[1] == UInt8('P') && bytes[2] == UInt8('6') "not a P6 PPM"
    i = 3
    vals = Int[]                          # width, height, maxval
    while length(vals) < 3
        while isspace(Char(bytes[i]))     # skip whitespace
            i += 1
        end
        if bytes[i] == UInt8('#')         # skip a comment line
            while bytes[i] != UInt8('\n')
                i += 1
            end
            continue
        end
        n = 0
        while !isspace(Char(bytes[i]))
            n = n * 10 + (bytes[i] - UInt8('0'))
            i += 1
        end
        push!(vals, n)
    end
    w, h, _ = vals
    i += 1                                # single whitespace before the data
    R = Matrix{Int}(undef, h, w)
    G = Matrix{Int}(undef, h, w)
    B = Matrix{Int}(undef, h, w)
    p = i
    @inbounds for r in 1:h, c in 1:w

        R[r, c] = bytes[p]
        G[r, c] = bytes[p + 1]
        B[r, c] = bytes[p + 2]
        p += 3
    end
    return R, G, B
end

function masks(R, G, B)
    (
        blue = (B .> 150) .& (G .> 150) .& (R .< 210) .& (B .>= R .+ 15),
        red = (R .> 120) .& (R .>= G .+ 50) .& (R .>= B .+ 50),
        dark = (R .< 120) .& (G .< 120) .& (B .< 120)
    )
end

# The onset figure is a blue-dominant daily bar chart with no orange (the
# age/sex pyramids use orange; the notification-week chart uses a darker
# steel blue and prints value labels). The tests are pixel fractions, not
# counts, because INSP re-renders the figure at whatever size the layout
# needs and an absolute threshold silently flips as the size moves: the
# blue floor already had to be lowered once when the figure shrank to
# 1009x583, and SitRep 072's larger 1277x799 rendering then pushed the
# crimson/pink-band anti-aliasing to 631 orange pixels, past a 500 cut.
#
# Measured over SitReps 059-080, on the caption page and its immediate
# neighbours (the page-fallback in onset_image widens the search there,
# which brings the provincial case map into the candidate pool - it is
# blue-heavy too, from the lake/river fill and legend swatches):
#   blue fraction    onset 0.066-0.125   province map 0.045-0.047
#   orange fraction  onset <= 0.0007     age/sex pyramids >= 0.053
#   red fraction     onset >= 0.046      (a floor, not a discriminator:
#                                        the notification-week chart is
#                                        also crimson-heavy)
# The map's blue fraction sits clear below every onset chart seen so far, so
# 0.055 (roughly the midpoint of the two clusters) discriminates with margin
# on both sides without needing a caption-text match.
function is_onset_curve(R, G, B)
    m = masks(R, G, B)
    orange = (R .> 200) .& (G .> 110) .& (G .< 195) .& (B .< 90)
    npx = length(R)
    return sum(m.blue) / npx > 0.055 && sum(orange) / npx < 0.01 &&
           sum(m.red) / npx > 0.01
end

function longest_run(col)
    best = cur = 0
    for v in col
        cur = v ? cur + 1 : 0
        best = max(best, cur)
    end
    return best
end

# Cluster nearly-adjacent indices, returning the mean of each cluster.
function cluster(idx; gap = 3)
    out = Int[]
    cl = Int[]
    for i in idx
        if !isempty(cl) && i - cl[end] <= gap
            push!(cl, i)
        else
            isempty(cl) || push!(out, floor(Int, sum(cl) / length(cl)))
            cl = [i]
        end
    end
    isempty(cl) || push!(out, floor(Int, sum(cl) / length(cl)))
    return out
end

# The 0/20/40/60 y-axis tick rows, read from the label strip just left of
# the vertical axis line (searched from column 30 to skip the left image
# border). Candidate strips are scored by the longest dark vertical run (the
# axis line itself), but only among strips whose rows form a plausible axis:
# at least three clusters, evenly spaced, with the last one (the 0 tick) on
# the baseline. Taking the longest run alone is not enough - in SitRep 067 a
# glyph stroke outruns the real axis line and yields a scale that halves
# every count.
function baseline_row(R, G, B, H)
    # The count-0 baseline is the plot's bottom border: a solid line running
    # almost the full chart width. Score rows by their longest contiguous run
    # under a near-gray threshold (<180); a run-length ranking under that
    # threshold correctly finds the border in every vintage, including
    # tighter-anti-aliased renders, unlike a per-row pixel sum.
    line = (R .< 180) .& (G .< 180) .& (B .< 180)
    lo = floor(Int, H * 0.4) + 1
    best_row, best_run = lo, 0
    for r in lo:H
        run = longest_run(@view line[r, :])
        if run > best_run
            best_run, best_row = run, r
        end
    end
    best_run < 100 && error("no baseline row found")
    return best_row
end

function y_axis_ticks(dark, base, H, W)
    best = nothing
    for x in 30:floor(Int, W * 0.13)
        seg = vec(sum(dark[1:min(base + 3, H), max(1, x - 10):(x - 1)];
            dims = 2))
        yt = cluster([y for y in 1:length(seg) if seg[y] >= 3])
        length(yt) < 3 && continue
        abs(yt[end] - base) > 3 && continue
        d = diff(yt)
        (minimum(d) <= 5 || maximum(d) > 1.15 * minimum(d)) && continue
        rank = (longest_run(@view dark[1:base, x]), -x)  # tie-break leftmost
        if best === nothing || rank > best[1]
            best = (rank, yt)
        end
    end
    best === nothing && error("no y-axis tick strip found")
    return best[2]
end

function digitize(R, G, B, last_tick::Date, y_step::Int = 20)
    H, W = size(R)
    m = masks(R, G, B)
    blue, red, dark = m.blue, m.red, m.dark
    base = baseline_row(R, G, B, H)       # count-0 baseline row
    # count scale from the y-axis ticks (0/20/40/60 through SitRep 083;
    # 0/25/50/75 from SitRep 087 - see Y_AXIS_STEP)
    yt = y_axis_ticks(dark, base, H, W)
    ppc = median(diff(yt)) / float(y_step) # pixels per count
    ytop, y0 = yt[1], yt[end]
    # x scale from the weekly tick marks below the baseline. The tick marks
    # are only a few pixels tall and shrink with the embedded figure
    # resolution (5-6 dark rows in the 1257x698 SitRep 064 rendering, 4 in
    # SitRep 066's 1275x623, 3 in SitRep 069/070's 1009x583), so step the
    # cut down until a full weekly row of ticks resolves instead of fixing
    # it at 4 and losing the axis entirely on the smaller figures. They sit
    # just below the baseline (2-6 rows) and, on the faint JPEG-compressed
    # figures (SitRep 081), can be only 1px tall, so cut must come all the
    # way down to 1 to resolve them; the window stops at base+6 so a wide
    # low-cut scan cannot pick up the x-axis date labels further down. Step
    # down through the cuts and keep the most complete regular weekly tick
    # row (the true axis has a fixed number of weekly ticks, so a too-strict
    # cut silently drops every other tick rather than failing).
    function weekly_ticks(mask)
        band = vec(sum(mask[(base + 2):min(base + 6, H), :]; dims = 1))
        best_n = 0
        found = Int[]
        for cut in (4, 3, 2, 1)
            cand = cluster([x for x in 1:W if band[x] >= cut])
            length(cand) >= 8 || continue
            if length(cand) > best_n
                best_n = length(cand)
                found = cand
            end
        end
        return found
    end
    best = weekly_ticks(dark)
    if isempty(best)
        # Only fall back to the <180 near-gray mask when the strict mask
        # finds nothing, so every already-committed vintage (059-107) keeps
        # digitising under the original threshold.
        line = (R .< 180) .& (G .< 180) .& (B .< 180)
        best = weekly_ticks(line)
    end
    isempty(best) && error("no x-axis weekly tick row found")
    xt = best
    ppd = median(diff(xt)) / 7.0          # pixels per day
    lastx = xt[end]                       # rightmost tick is always real
    # per-column stacked bar height, flooded up from the baseline
    bc = zeros(Int, W)
    rc = zeros(Int, W)
    for x in 1:W
        r = y0 - 1
        miss = b = rr = 0
        while r > ytop - 2 && r >= 1
            if blue[r, x]
                b += 1
                miss = 0
            elseif red[r, x]
                rr += 1
                miss = 0
            else
                miss += 1
                miss > 6 && break
            end
            r -= 1
        end
        bc[x] = b
        rc[x] = rr
    end
    tot = bc .+ rc
    nz = findall(>(2), tot)
    barmin, barmax = minimum(nz), maximum(nz)
    rows = Tuple{Date, Int, Int}[]
    for off in -105:3
        cx = lastx + off * ppd
        (cx < barmin - ppd || cx > barmax + ppd) && continue
        lo = round(Int, cx - ppd * 0.45)
        hi = round(Int, cx + ppd * 0.45)
        cols = max(1, lo):min(W, hi)
        bvals = Float64.(bc[cols])
        rvals = Float64.(rc[cols])
        maximum(bvals .+ rvals) < 1 && continue
        alive = round(Int, quantile(bvals, 0.75) / ppc)
        dead = round(Int, quantile(rvals, 0.75) / ppc)
        push!(rows, (last_tick + Day(off), alive, dead))
    end
    # drop trailing zero rows and isolated tiny strays past the curve tail
    while !isempty(rows) && rows[end][2] + rows[end][3] == 0
        pop!(rows)
    end
    while length(rows) >= 2
        gap = value(rows[end][1] - rows[end - 1][1])
        if gap > 1 && rows[end][2] + rows[end][3] <= 2
            pop!(rows)
        else
            break
        end
    end
    return rows
end

function _onset_page(pdf)
    # The onset figure usually sits on the page whose text carries its
    # caption.
    npages = parse(Int, match(r"Pages:\s*(\d+)",
        read(`pdfinfo $pdf`, String)).captures[1])
    for p in 1:npages
        txt = lowercase(read(`pdftotext -layout -f $p -l $p $pdf -`, String))
        if occursin("date de debut des symptom", txt) ||
           occursin("date de début des symptôm", txt)
            return p, npages
        end
    end
    return nothing, npages
end

function _best_onset_image(pdf, page, wd)
    # pdfimages writes PPM (P6) for RGB images by default; no format flag
    run(`pdfimages -f $page -l $page $pdf $(joinpath(wd, "p"))`)
    best = nothing
    for name in sort(readdir(wd))
        endswith(name, ".ppm") || continue
        R, G, B = read_ppm(joinpath(wd, name))
        if is_onset_curve(R, G, B) &&
           (best === nothing || length(R) > length(best[1]))
            best = (R, G, B)
        end
        rm(joinpath(wd, name))
    end
    return best
end

# Extract the onset-curve figure from a SitRep PDF as R, G, B matrices.
function onset_image(pdf)
    page, npages = _onset_page(pdf)
    page === nothing && return nothing
    best = mktempdir(wd -> _best_onset_image(pdf, page, wd))
    best !== nothing && return best
    # SitRep 080 embeds the chart on page 5 under a mislabelled caption ("par
    # semaine de notification") while the matching "date de debut des
    # symptomes" caption text sits on page 6 with no image of its own, so the
    # caption-text page lookup lands one page short of the real figure. Widen
    # to the immediate neighbours only (not the whole document): the map and
    # other embedded figures elsewhere in the report are large enough, and
    # blue enough in places (lakes, legends), to satisfy is_onset_curve too,
    # so a document-wide scan silently grabs the wrong image.
    for q in (page - 1, page + 1)
        if 1 <= q <= npages
            best = mktempdir(wd -> _best_onset_image(pdf, q, wd))
            best !== nothing && return best
        end
    end
    return nothing
end

function main(pdf_dir = "data/sitrep_pdfs",
        out_csv = "data/onset_curve_scanned.csv")
    open(out_csv, "w") do io
        println(io, "sitrep,report_date,onset_date,confirmed_alive," *
                    "confirmed_dead,confirmed_total")
        for (sr, report_date, last_tick) in CONFIG
            pdf = joinpath(pdf_dir, "SitRep_MVE_$(sr)_2026.pdf")
            if !isfile(pdf)
                @warn "skip $sr: $pdf not found"
                continue
            end
            img = onset_image(pdf)
            if img === nothing
                @warn "skip $sr: no onset curve found"
                continue
            end
            rows = digitize(img..., last_tick, get(Y_AXIS_STEP, sr, 20))
            total = sum(a + d for (_, a, d) in rows)
            @printf("SitRep %s (%s): %d onset days, total %d confirmed\n",
                sr, report_date, length(rows), total)
            for (onset, alive, dead) in rows
                println(io, join((sr, report_date, onset, alive, dead,
                        alive + dead), ","))
            end
        end
    end
    println("wrote $out_csv")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(get(ARGS, 1, "data/sitrep_pdfs"),
        get(ARGS, 2, "data/onset_curve_scanned.csv"))
end
