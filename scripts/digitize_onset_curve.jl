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
#   * count scale = the y-axis tick marks (0/20/40/60 or 0/25/50/75), evenly
#     spaced, giving pixels-per-count = tick-spacing / y_step. The ticks
#     come from the strict dark mask, or from the near-gray mask when that
#     reads a grid finer by more than the strip tolerance;
#   * date scale = the weekly x-axis tick marks. Candidate tick rows come
#     from a strict and a near-gray mask at several cuts, and the one whose
#     regular chain from the rightmost tick is longest wins. Pixels per day
#     is the least-squares slope over that chain and each day is anchored on
#     the nearest chain tick at or before it, so a day 150 days back
#     does not drift with the rounding of one spacing. The rightmost tick's
#     date is in CONFIG, read off the axis;
#   * each daily bar = the bar's own pixel columns: the interval between
#     the two consecutive outline columns about a day apart whose midpoint
#     is nearest the day's grid position, else the window
#     [cx - ppd/2, cx + ppd/2] clipped to the nearest outline column on each
#     side. Outline columns are those mostly dark over their run. Every
#     column is read as the run of non-page pixels up from the baseline,
#     bridging up to three page pixels when bar colour resumes (JPEG ringing
#     between the stacked segments) and skipping neutral pixels on the tick
#     rows and tick columns (gridlines); a day whose columns then neither
#     show a bar nor leave it empty, or that reads below the two equal
#     outline columns around it beside a washed tick column, is read again
#     without the tick-column skip, which washed fill on a tick column can
#     break. The run's top is
#     its highest pixel darker than an anti-alias, which is the bar's
#     outline; chroma-washed
#     fill inside the run is crossed on the way up. The bar height is the
#     height at least two of its interior columns agree on to within a
#     pixel, or the tallest
#     interior column when none do (a 3-4 px bar has one saturated column
#     and one or two washed ones that read low). Inside the pink
#     incomplete-data band, where the faded bars can lose their outlines,
#     a day whose columns disagree reads the level run of columns centred
#     nearest it instead. Half a pixel of
#     outline is subtracted before dividing by pixels-per-count. The dead
#     segment is the count of crimson pixels in the chosen column. A column
#     under the red dashed first-positive-result line (evenly spaced
#     crimson dashes above its run) is read to its highest dark pixel
#     instead, never above its run, with the dead segment the crimson run
#     below that.
#
# The Python port (scripts/digitize_onset_curve.py) must match every step
# above pixel for pixel: the same pixel classes (page, neutral, light,
# crimson, dark, saturated), the same run rule with its gap of three and gridline skip,
# the same outline thresholds (dark fraction 0.25 strict with a floor of
# six pixels, 0.1 soft with a floor of five, both only over runs taller
# than four pixels, no saturated pixel in a tenth of the run, and three
# bridged gaps), the same gap rule, the
# same interval choice, the same mode-or-maximum rule with its support of
# two, and the same
# rounding. All of it is expressible on numpy arrays with cumulative sums
# and per-column loops.
#
# Accuracy: against the printed `n` in every figure that carries one
# (SitReps 064-130, read by OCR and checked by eye), the digitised total is
# within 2.1% everywhere and within 0.5% on 43 of the 60 vintages. The
# largest gaps are SitRep 126 at -2.1% (5650 against n = 5 771), SitRep 129
# at -1.9% (5784 against n = 5 896) and SitRep 118 at -1.3% (5192 against
# n = 5 263). Individual daily bars carry pixel rounding of about +/-1 case
# at the small September renders (2.8 px per count) and less before. The
# faded bars inside the `donnees potentiellement incompletes` band are read
# like any other.
#
# Late reporting only ever adds cases, so an onset date's count must be
# non-decreasing across vintages. On onset dates more than three weeks
# before the earlier vintage's report date, consecutive distinct snapshots
# differ by 0.40 cases per day on average (L1) and fall on 14% of such days,
# almost always by a single case, so a between-vintage increment of one is
# at the noise floor and anything larger is signal.
#
# The values are approximate and are not fitted by the model; they are captured
# for later use (see data/README.md and #488).
#
# Dependencies: poppler (`pdfimages`, `pdftotext`, `pdfinfo`) on PATH. No
# Julia packages beyond stdlib.
#
# Incremental by default. Digitising a vintage means extracting the figure
# and walking it pixel by pixel, and a data update adds one or two vintages
# to a file that already holds every earlier one. So a run reuses the rows
# `out_csv` already carries and opens the PDF only for the CONFIG vintages
# missing from it. The rows are written back in CONFIG order either way, so
# an incremental run and a full one produce the same file.
#
# A change to the digitiser itself does not invalidate those reused rows, so
# re-run with `--rebuild` after touching the digitising code, which re-reads
# every vintage. Each run prints how many vintages it reused and how many it
# read, so a run that should have re-read everything and did not is visible.
#
# Usage:
#   julia scripts/digitize_onset_curve.jl [pdf_dir] [out_csv] [--rebuild]
# Defaults: pdf_dir = data/sitrep_pdfs, out_csv = data/onset_curve_scanned.csv
# Download the PDFs first with scripts/download_sitreps.jl.

using Dates: Date, Day, value
using Statistics: median
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
    # "098" is deliberately absent. It is the only vintage INSP embedded
    # losslessly rather than as JPEG, so the fixed colour thresholds below
    # keep a fringe of each bar that JPEG blur costs every other vintage,
    # and it reads about 7% high on the same underlying data. Excluding it
    # keeps a vintage on a different bias scale out of the between-vintage
    # increments this file feeds. The evidence, and the controls that rule
    # out the render size, are in data/README.md. Read them before adding
    # it back.
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
    ("109", Date(2026, 8, 31), Date(2026, 8, 31)),
    # "110" is deliberately absent. Its page-4 figure carries the same
    # outer caption as every other vintage ("par date de début des
    # symptômes") but the embedded chart's own internal title and x-axis
    # read "par date de NOTIFICATION" (n = 5 710) - a genuine basis change,
    # confirmed by extracting and viewing the raw embedded image rather
    # than trusting the caption. Digitising it would silently inject a
    # different-basis series into the reporting-triangle stream. See
    # data/README.md and issue #644.
    ("111", Date(2026, 9, 2), Date(2026, 8, 31)),
    ("112", Date(2026, 9, 3), Date(2026, 8, 31)),
    ("113", Date(2026, 9, 4), Date(2026, 8, 31)),
    ("114", Date(2026, 9, 5), Date(2026, 8, 31)),
    ("115", Date(2026, 9, 6), Date(2026, 9, 7)),
    ("116", Date(2026, 9, 7), Date(2026, 9, 7)),
    # 117 and 118 keep 116's 07 September tick and plot a window that
    # retreats behind their own report dates (to 09-03 and 09-04). That
    # only shortens the pair's coverage intersection, which
    # `load_onset_curve` drops as unobserved.
    ("117", Date(2026, 9, 8), Date(2026, 9, 7)),
    ("118", Date(2026, 9, 9), Date(2026, 9, 7)),
    ("119", Date(2026, 9, 10), Date(2026, 9, 7)),
    ("120", Date(2026, 9, 11), Date(2026, 9, 7)),
    ("121", Date(2026, 9, 12), Date(2026, 9, 7)),
    ("122", Date(2026, 9, 13), Date(2026, 9, 14)),
    ("123", Date(2026, 9, 14), Date(2026, 9, 14)),
    ("124", Date(2026, 9, 15), Date(2026, 9, 14)),
    ("125", Date(2026, 9, 16), Date(2026, 9, 14)),
    ("126", Date(2026, 9, 17), Date(2026, 9, 14)),
    ("127", Date(2026, 9, 18), Date(2026, 9, 14)),
    ("128", Date(2026, 9, 19), Date(2026, 9, 14)),
    ("129", Date(2026, 9, 20), Date(2026, 9, 21)),
    ("130", Date(2026, 9, 21), Date(2026, 9, 21)),
    ("131", Date(2026, 9, 22), Date(2026, 9, 21)),
    ("132", Date(2026, 9, 23), Date(2026, 9, 21)),
    ("133", Date(2026, 9, 24), Date(2026, 9, 21)),
    ("134", Date(2026, 9, 25), Date(2026, 9, 21)),
    ("136", Date(2026, 9, 27), Date(2026, 9, 28)),
    ("137", Date(2026, 9, 28), Date(2026, 9, 28)),
    # "138" is left out. Its render loses the outline between the 15 and
    # 16 June bars, so 16 June reads 40 where every other vintage and the
    # dashboard read 28 to 29, and its 14 May reads 12 against 8 (#1061).
    ("139", Date(2026, 9, 30), Date(2026, 9, 28)),
    ("140", Date(2026, 10, 1), Date(2026, 9, 28)),
    ("141", Date(2026, 10, 2), Date(2026, 9, 28)),
    # SitRep 142 is not published. 144 and 145 reprint 143's figure byte
    # for byte (same image md5).
    ("143", Date(2026, 10, 4), Date(2026, 10, 5)),
    ("144", Date(2026, 10, 5), Date(2026, 10, 5)),
    ("145", Date(2026, 10, 6), Date(2026, 10, 5)),
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
    "109" => 25,
    "111" => 25,
    "112" => 25,
    "113" => 25,
    "114" => 25,
    "115" => 25,
    "116" => 25,
    "117" => 25,
    "118" => 25,
    "119" => 25,
    "120" => 25,
    "121" => 25,
    "122" => 25,
    "123" => 25,
    "124" => 25,
    "125" => 25,
    "126" => 25,
    "127" => 25,
    "128" => 25,
    "129" => 25,
    "130" => 25,
    "131" => 25,
    "132" => 25,
    "133" => 25,
    "134" => 25,
    "136" => 25,
    "137" => 25,
    "139" => 25,
    "140" => 25,
    "141" => 25,
    "143" => 25,
    "144" => 25,
    "145" => 25
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
    return (
        blue = (B .> 150) .& (G .> 150) .& (R .< 210) .& (B .>= R .+ 15),
        red = (R .> 120) .& (R .>= G .+ 50) .& (R .>= B .+ 50),
        dark = (R .< 120) .& (G .< 120) .& (B .< 120),
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
        seg = vec(
            sum(
                dark[1:min(base + 3, H), max(1, x - 10):(x - 1)];
                dims = 2
            )
        )
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

# The y-axis tick rows from the strict dark mask, or from the near-gray
# `line` mask when that reads a finer grid. The small renders anti-alias
# the tick marks into the 120-180 range, so the strict mask can find no
# regular strip at all (SitRep 112) or keep every other tick and take a
# title glyph above the plot for the top one (SitRep 133, spacing 144-155
# against the line mask's 72). Only a spacing finer by more than the strip
# tolerance counts, so the strict rows stay wherever both masks read the
# same grid.
function y_tick_rows(dark, line, base, H, W)
    strict = try
        y_axis_ticks(dark, base, H, W)
    catch e
        e isa ErrorException || rethrow()
        nothing
    end
    near = try
        y_axis_ticks(line, base, H, W)
    catch e
        e isa ErrorException || rethrow()
        nothing
    end
    strict === nothing && near === nothing &&
        error("no y-axis tick strip found")
    strict === nothing && return near
    near === nothing && return strict
    finer = 1.15 * median(diff(near)) < median(diff(strict))
    return finer ? near : strict
end

# Pixel classes. Page is white (with JPEG chroma noise) and the pink
# `donnees potentiellement incompletes` band. Neutral is gridline gray.
# Light is the pale, low-saturation pixel a bar's top edge leaves above
# its outline: anti-aliasing, chroma noise, washed fill and gridlines;
# saturated fill and the outline are not light. Dark is the outline.
function pixel_classes(R, G, B)
    lo = min.(R, min.(G, B))
    hi = max.(R, max.(G, B))
    spread = hi .- lo
    page = ((lo .>= 228) .& (spread .<= 25)) .|
        ((R .>= 238) .& (G .>= 200) .& (B .>= 200) .& (R .- G .>= 15))
    neutral = (lo .>= 170) .& (spread .<= 12)
    light = (hi .>= 190) .& (spread .< 60)
    crimson = (R .- max.(G, B)) .>= 25
    darkpx = (R .< 150) .& (G .< 150) .& (B .< 150)
    saturated = spread .>= 30
    return page, neutral, light, crimson, darkpx, saturated
end

# Per-column run of non-page pixels up from the baseline. Gridlines lie on
# the tick rows and tick columns, so a neutral pixel there is page. Up to
# `gap` page pixels are bridged when a non-page pixel follows (JPEG ringing
# between the stacked segments), which still stops under the dashed
# vertical line's 7 px gaps. The run's top is the highest pixel darker than
# `light`: the bar's outline, never a stray anti-alias or noise pixel
# above it, and chroma-washed fill inside the run is crossed on the way up.
# Returns the run height (to the outline), the run's full non-page extent,
# the crimson, dark and saturated counts and the number of page gaps
# bridged per column (one at the segment junction; many up a dashed line).
function column_runs(
        page, neutral, light, crimson, darkpx, saturated, y0, gridrows,
        gridcols; gap = 3
    )
    H, W = size(page)
    grid = falses(H, W)
    for r in gridrows, d in -1:1
        1 <= r + d <= H && (grid[r + d, :] .= true)
    end
    for c in gridcols, d in -1:1
        1 <= c + d <= W && (grid[:, c + d] .= true)
    end
    h = zeros(Int, W)
    hp = zeros(Int, W)
    nr = zeros(Int, W)
    nd = zeros(Int, W)
    ns = zeros(Int, W)
    nb = zeros(Int, W)
    for x in 1:W
        r = y0 - 1
        miss = 0
        top = y0
        last = y0
        rr = dd = ss = bb = 0
        cr = cd = cs = cb = 0
        while r >= 1
            if !page[r, x] && !(grid[r, x] && neutral[r, x])
                miss >= 2 && (bb += 1)
                miss = 0
                last = r
                crimson[r, x] && (rr += 1)
                darkpx[r, x] && (dd += 1)
                saturated[r, x] && (ss += 1)
                if !light[r, x]
                    top = r
                    cr, cd, cs, cb = rr, dd, ss, bb
                end
            else
                miss += 1
                miss > gap && break
            end
            r -= 1
        end
        h[x] = y0 - top
        hp[x] = y0 - last
        nr[x] = cr
        nd[x] = cd
        ns[x] = cs
        nb[x] = cb
    end
    return h, hp, nr, nd, ns, nb
end

# The regular weekly chain ending on the rightmost tick, as (week index, x)
# pairs. Walking left, a spacing of one to three weeks within 8% (at least
# 2.5 px) of the median spacing continues the chain; anything else ends it,
# which drops the y-axis line and label strokes that cluster as ticks at
# the left, and a stray cluster right of the last tick leaves a chain of
# one. The near-grey mask on the September renders splits a tick in two or
# sees a stroke beside it (SitReps 138 and 141), so up to two clusters off
# the grid are stepped over when the next one lands on it.
function tick_chain(xt)
    s = median(diff(xt))
    ks = [0]
    xs = [xt[end]]
    fits(d) = (
        k = round(Int, d / s);
        1 <= k <= 3 && abs(d - k * s) <= max(2.5, 0.08 * s)
    )
    j = length(xt) - 1
    while j >= 1
        # a cluster off the weekly grid (a tick split in two, or a stray
        # stroke beside one) is stepped over when one of the next two
        # clusters lands on the grid
        skip = findfirst(i -> j - i >= 1 && fits(xs[1] - xt[j - i]), 0:2)
        skip === nothing && break
        j -= skip - 1
        d = xs[1] - xt[j]
        pushfirst!(ks, ks[1] - round(Int, d / s))
        pushfirst!(xs, xt[j])
        j -= 1
    end
    return ks, xs
end

# Columns under the red dashed first-positive-result line: at least twelve
# crimson dashes above the column's run, evenly spaced, with white page in
# at least half of the gaps between them (the pink band's text and edge
# have none).
function dash_columns(crimson, white, hp, y0)
    W = size(crimson, 2)
    out = falses(W)
    for x in 1:W
        starts = Int[]
        gapwhite = Bool[]
        run = 0
        seen = false
        for r in 1:(y0 - hp[x] - 2)
            if crimson[r, x]
                run += 1
            else
                if run >= 2
                    push!(starts, r - run)
                    push!(gapwhite, seen)
                    seen = false
                end
                run = 0
                white[r, x] && (seen = true)
            end
        end
        length(starts) >= 12 || continue
        d = diff(starts)
        md = median(d)
        out[x] = count(v -> abs(v - md) <= 2, d) >= 0.8 * length(d) &&
            2 * count(gapwhite[2:end]) >= length(d)
    end
    return out
end

# A bar read through a dash column: the run walked as `column_runs` walks
# it, with the height taken to its highest dark pixel (the bar's outline,
# which the dash's red never is) and the dead segment to the crimson run
# down from there. Returns `(0, 0)` when the run has no dark pixel.
function dash_bar(page, neutral, darkpx, crimson, y0, x, gridrows; gap = 3)
    isgrid(r) = any(g -> abs(r - g) <= 1, gridrows)
    top = 0
    miss = 0
    r = y0 - 1
    while r >= 1
        if !page[r, x] && !(isgrid(r) && neutral[r, x])
            miss = 0
            darkpx[r, x] && (top = r)
        else
            miss += 1
            miss > gap && break
        end
        r -= 1
    end
    top == 0 && return 0, 0
    n = 0
    r = top
    while r < y0 && crimson[r, x]
        n += 1
        r += 1
    end
    return y0 - top, n
end

# The value of `h` over `cols` that most columns agree with to within a
# pixel, and how many do; ties go to the value nearest `cx` by column.
function modal_height(h, cols, cx)
    best = 0
    bestkey = (-1, -Inf)
    for u in unique(h[cols])
        c = count(x -> abs(h[x] - u) <= 1, cols)
        d = minimum(abs(x - cx) for x in cols if h[x] == u)
        key = (c, -d)
        if key > bestkey
            best = u
            bestkey = key
        end
    end
    return best, bestkey[1]
end

# One day's bar over its interior `cols`, as the run height and the
# crimson count of the column read, `:empty` when there is no bar, or
# `nothing` when the columns neither show a bar nor leave the day empty.
function bar_height(h, hp, nr, cols, cx, ppc)
    # a day is empty when half or more of its columns are page from
    # the baseline up (the baseline's own anti-alias apart); a column
    # that is fill all the way but never shows an outline pixel
    # (chroma-washed) abstains rather than reading 0
    gaps = count(x -> h[x] == 0 && hp[x] <= 2 * ppc, cols)
    2 * gaps >= length(cols) && return :empty
    resolved = [x for x in cols if h[x] > 0]
    if isempty(resolved)
        # no outline pixel in any column: read the fill's extent where
        # two columns agree on it, else it is a halo, not a bar
        hb, support = modal_height(hp, cols, cx)
        support >= 2 || return nothing
        jb = cols[findfirst(x -> hp[x] == hb, cols)]
    else
        hb, support = modal_height(h, resolved, cx)
        support >= 2 || (hb = maximum(h[resolved]))
        jb = resolved[findfirst(x -> h[x] == hb, resolved)]
    end
    hb < 1 && return :empty
    return hb, nr[jb]
end

# Whether a column with an outline reads within a pixel of another column.
function two_agree(h, cols)
    return any(
        count(y -> abs(h[y] - h[x]) <= 1, cols) >= 2 for x in cols if h[x] > 0
    )
end

# The level run (two or more consecutive columns, each within a pixel of the
# one before it) that reaches into `cols` with its centre nearest `cx`, as
# its first and last column, or `nothing` when no such run does. A faded bar
# with no outline reads as one such run; the smear between two faded bars
# is a single column.
function nearest_plateau(h, cols, cx)
    W = length(h)
    level(x) = h[x] > 0 && h[x + 1] > 0 && abs(h[x + 1] - h[x]) <= 1
    a = first(cols)
    while a > 1 && level(a - 1)
        a -= 1
    end
    # the best run is kept in plain integers rather than a tuple rebuilt in
    # a short-circuit: Julia 1.10 returned that tuple with its start
    # overwritten by the loop's later start
    found = false
    best_d, best_lo, best_hi = Inf, 0, 0
    x = a
    while x <= last(cols)
        b = x
        while b < W && level(b)
            b += 1
        end
        if b > x && b >= first(cols)
            d = abs((x + b) / 2 - cx)
            if !found || d < best_d
                found = true
                best_d, best_lo, best_hi = d, x, b
            end
        end
        x = b + 1
    end
    return found ? (best_lo, best_hi) : nothing
end

# Columns inside the pink `donnees potentiellement incompletes` band: above
# the column's run, from the top y-axis tick down, the page is at least as
# often pink as white.
function band_columns(R, G, B, white, hp, y0, ytop)
    pink = (R .>= 238) .& (G .>= 200) .& (B .>= 200) .& (R .- G .>= 15)
    W = size(R, 2)
    out = falses(W)
    for x in 1:W
        top = y0 - hp[x] - 2
        top > ytop || continue
        np = count(@view pink[ytop:top, x])
        out[x] = np > 0 && np >= count(@view white[ytop:top, x])
    end
    return out
end

function digitize(R, G, B, last_tick::Date, y_step::Int = 20)
    H, W = size(R)
    m = masks(R, G, B)
    dark = m.dark
    base = baseline_row(R, G, B, H)       # count-0 baseline row
    # count scale from the y-axis ticks (0/20/40/60 through SitRep 083;
    # 0/25/50/75 from SitRep 087 - see Y_AXIS_STEP)
    line = (R .< 180) .& (G .< 180) .& (B .< 180)
    yt = y_tick_rows(dark, line, base, H, W)
    ppc = median(diff(yt)) / float(y_step) # pixels per count
    y0 = yt[end]
    # x scale from the weekly tick marks 2-5 rows below the baseline (on
    # the short September renders the row below that reaches the tops of
    # the date labels, whose strokes merge into the tick clusters). The
    # marks shrink with the render (down to 1 px tall on the JPEG figures)
    # and the strict mask loses some of them on the small renders, so both
    # masks are tried at every cut and the tick row whose regular weekly
    # chain from the rightmost tick is longest wins. A mask that adds a
    # stray cluster right of the last tick starts its chain at length one.
    best_n = 0
    xt = Int[]
    for mask in (dark, line)
        band = vec(sum(mask[(base + 2):min(base + 5, H), :]; dims = 1))
        for cut in (4, 3, 2, 1)
            cand = cluster([x for x in 1:W if band[x] >= cut])
            length(cand) >= 8 || continue
            n = length(tick_chain(cand)[1])
            if n > best_n
                best_n = n
                xt = cand
            end
        end
    end
    isempty(xt) && error("no x-axis weekly tick row found")
    ks, xs = tick_chain(xt)
    n = length(xs)
    n >= 2 || error("weekly tick chain too short")
    w = 7.0 .* ks
    ppd = (n * sum(w .* xs) - sum(w) * sum(xs)) /
        (n * sum(w .^ 2) - sum(w)^2)   # pixels per day
    page, neutral, light, crimson, darkpx, saturated = pixel_classes(R, G, B)
    h, hp, nr, nd, ns, nb = column_runs(
        page, neutral, light, crimson, darkpx, saturated, y0, yt, xs
    )
    h1, hp1, nr1 = column_runs(
        page, neutral, light, crimson, darkpx, saturated, y0, yt, Int[]
    )
    # outline columns are mostly dark over their run (a short bar's top
    # and junction lines are a few dark pixels in every column, so the
    # floor keeps its interior as interior), and a column with
    # no saturated pixel is a gray line (the y-axis, the panel border and
    # their anti-alias) rather than a bar, as is one that bridged three or
    # more page gaps (the dashed first-positive-result line, whose gaps
    # the small renders shrink inside the bridge); an outline drawn across
    # two columns leaves a softer second column that still carries the
    # neighbour's height, dropped when anything else is left
    isborder = (h .> 4) .&
        ((nd .>= max.(0.25 .* h, 6)) .| (ns .< 0.1 .* h) .| (nb .>= 3))
    outline = (h .> 4) .& ((nd .>= max.(0.25 .* h, 6)) .| (nb .>= 3))
    soft = (h .> 4) .& (nd .>= max.(0.1 .* h, 5))
    nz = findall((h .> 2) .& .!isborder)
    # a column under the dashed line reads its bar from the outline; the
    # dash can only have added height. An outline read under two counts is
    # taken as the dash's own edge and the run is kept, so a real bar of one
    # count under the line keeps its dash-inflated height
    white = (min.(R, G, B) .>= 228) .& (max.(R, G, B) .- min.(R, G, B) .<= 25)
    hread = copy(h)
    nread = copy(nr)
    for x in findall(dash_columns(crimson, white, hp, y0))
        hd, nd_x = dash_bar(page, neutral, darkpx, crimson, y0, x, yt)
        hd > 2 * ppc || continue
        hread[x] = min(hd, h[x])
        nread[x] = min(nd_x, hread[x])
    end
    barmin, barmax = minimum(nz), maximum(nz)
    inband = band_columns(R, G, B, white, hp, y0, yt[1])
    rows = Tuple{Date, Int, Int}[]
    for off in (7 * ks[1] - 7):3
        # anchor on the nearest chain tick at or before the day
        j = findlast(k -> 7 * k <= off, ks)
        j === nothing && (j = 1)
        cx = xs[j] + (off - 7 * ks[j]) * ppd
        (cx < barmin - ppd || cx > barmax + ppd) && continue
        lo = max(1, ceil(Int, cx - ppd / 2 + 0.5))
        hi = min(W, floor(Int, cx + ppd / 2 - 0.5))
        # the bar is the interval between two consecutive outline columns
        # about a day apart whose midpoint is nearest cx; when the day grid
        # lands on an outline that picks the right side of it. With no such
        # pair (an outline the render lost) the window is clipped to the
        # nearest outline on each side instead.
        c = round(Int, cx)
        reach = ceil(Int, ppd)
        # past the last tick the bars are faded into the band and carry no
        # saturated pixel, so only a dark outline bounds a bar there
        border = off > 0 ? outline : isborder
        near = [x for x in max(1, c - reach):min(W, c + reach) if border[x]]
        best = nothing
        for i in 1:(length(near) - 1)
            a, b = near[i], near[i + 1]
            abs(b - a - ppd) <= 1.5 || continue
            d = abs((a + b) / 2 - cx)
            (best === nothing || d < best[1]) && (best = (d, a, b))
        end
        if best !== nothing && best[1] <= ppd / 2
            lo, hi = best[2] + 1, best[3] - 1
        else
            bl = findlast(x -> border[x], max(1, c - reach):(c - 1))
            bl === nothing || (lo = max(lo, max(1, c - reach) + bl))
            br = findfirst(x -> border[x], (c + 1):min(W, c + reach))
            br === nothing || (hi = min(hi, c + br - 1))
        end
        cols = [x for x in lo:hi if !border[x]]
        if isempty(cols) && best !== nothing
            # the clipped window holds only outline columns (two adjacent
            # outlines on the day's centre): read the nearest outline pair
            lo, hi = best[2] + 1, best[3] - 1
            cols = [x for x in lo:hi if !border[x]]
        end
        interior = [x for x in cols if !soft[x]]
        isempty(interior) || (cols = interior)
        isempty(cols) && continue
        bar = bar_height(hread, hp, nread, cols, cx, ppc)
        # past the last tick a window can take in a column of the faded bar
        # before it; when no two columns agree, the column nearest the day's
        # centre is read (its fill extent when it shows no outline)
        if off > 0 && bar isa Tuple && length(cols) >= 2 &&
                !two_agree(hread, cols)
            jn = argmin(x -> (abs(x - cx), x), cols)
            hn = hread[jn] > 0 ? hread[jn] : hp[jn]
            bar = (hn, nread[jn])
        end
        # before the last tick the band's faded bars can lose their outline
        # too, and a window then straddles two of them; when no two columns
        # agree, the level run of columns centred nearest the day is read,
        # from its column nearest the day's centre
        if off <= 0 && bar isa Tuple && length(cols) >= 2 &&
                !two_agree(hread, cols) && all(x -> inband[x], cols)
            run = nearest_plateau(hread, cols, cx)
            if run !== nothing
                jn = argmin(x -> (abs(x - cx), x), run[1]:run[2])
                bar = (hread[jn], nread[jn])
            end
        end
        # washed fill on a tick column can read as gridline and break each run
        # at a different row; a day that reads as neither empty nor a bar, or
        # as a bar with no outline above the baseline's anti-alias in any
        # column, is read again without the skip. That read can climb the
        # gridline above the bar, so it is capped at the taller of the
        # outline columns bounding the window. A day between two outline
        # columns of the same height, with a washed column beside the tick
        # that shows fill but no outline, is read again too when it reads
        # below them: its fill broke at a gridline row
        bounds = (first(cols) - 1, last(cols) + 1)
        paired = all(x -> 1 <= x <= W && isborder[x], bounds) &&
            abs(h[bounds[1]] - h[bounds[2]]) <= 1
        washed = any(
            x -> hread[x] == 0 && hp[x] > 2 * ppc &&
                any(t -> abs(x - t) <= 1, xs), cols
        )
        if bar === nothing || (bar isa Tuple && maximum(hread[cols]) <= 1) ||
                (bar isa Tuple && paired && washed && bar[1] < h[bounds[1]] - 1)
            bar = bar_height(h1, hp1, nr1, cols, cx, ppc)
            edge = [
                h[x] for x in (first(cols) - 1, last(cols) + 1)
                    if 1 <= x <= W && isborder[x]
            ]
            if bar isa Tuple && !isempty(edge) && bar[1] > maximum(edge)
                bar = (maximum(edge), min(bar[2], maximum(edge)))
            end
        end
        bar isa Tuple || continue
        hb, hr = bar
        total = round(Int, max(0.0, hb - 0.5) / ppc)
        dead = min(total, round(Int, max(0.0, hr - 0.5) / ppc))
        push!(rows, (last_tick + Day(off), total - dead, dead))
    end
    # drop leading and trailing zero rows (a stray anti-alias column near
    # the y-axis or the band edge reads as a bar of height 0) and isolated
    # tiny strays past the curve tail
    while !isempty(rows) && rows[1][2] + rows[1][3] == 0
        popfirst!(rows)
    end
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
    npages = parse(
        Int, match(
            r"Pages:\s*(\d+)",
            read(`pdfinfo $pdf`, String)
        ).captures[1]
    )
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

const OUT_HEADER = "sitrep,report_date,onset_date,confirmed_alive," *
    "confirmed_dead,confirmed_total"

## Rows `out_csv` already holds, keyed by the SitRep number in its first
## field, so a run can reuse a vintage it has already read rather than open
## the PDF again. An absent, empty or differently-headed file yields nothing
## and every vintage is read, which is what a first run does anyway.
function digitised_rows(out_csv)
    rows = Dict{String, Vector{String}}()
    isfile(out_csv) || return rows
    lines = readlines(out_csv)
    (isempty(lines) || strip(lines[1]) != OUT_HEADER) && return rows
    for line in lines[2:end]
        isempty(strip(line)) && continue
        sr = first(split(line, ','))
        push!(get!(rows, String(sr), String[]), line)
    end
    return rows
end

function main(
        pdf_dir = "data/sitrep_pdfs",
        out_csv = "data/onset_curve_scanned.csv";
        rebuild::Bool = false
    )
    ## Read before the file is opened for writing, which truncates it.
    cached = rebuild ? Dict{String, Vector{String}}() : digitised_rows(out_csv)
    reused = 0
    read_now = 0
    open(out_csv, "w") do io
        println(io, OUT_HEADER)
        for (sr, report_date, last_tick) in CONFIG
            ## Already digitised, so its rows are carried through untouched.
            ## They are written in CONFIG order like any other, so reusing
            ## them cannot reorder the file.
            if haskey(cached, sr)
                foreach(line -> println(io, line), cached[sr])
                reused += 1
                continue
            end
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
            ## An onset date can never sit later than the axis of the
            ## report that draws it, and that axis runs at most a day past
            ## the rapportage date (the date-de-publication lag). The
            ## window above self-calibrates from pixel content and can
            ## read a few stray days past the last labelled tick when the
            ## figure's own "donnees potentiellement incompletes" band
            ## extends that far (SitRep 115); drop those here rather than
            ## loosen the invariant test/test_onset_digitiser.jl checks.
            filter!(r -> r[1] <= report_date + Day(1), rows)
            total = sum(a + d for (_, a, d) in rows)
            @printf(
                "SitRep %s (%s): %d onset days, total %d confirmed\n",
                sr, report_date, length(rows), total
            )
            for (onset, alive, dead) in rows
                println(
                    io, join(
                        (
                            sr, report_date, onset, alive, dead,
                            alive + dead,
                        ), ","
                    )
                )
            end
            read_now += 1
        end
    end
    @printf(
        "wrote %s: %d vintages read, %d reused from the existing file\n",
        out_csv, read_now, reused
    )
    return reused > 0 && println(
        "re-run with --rebuild to re-read every vintage " *
            "after changing the digitiser"
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    args = filter(a -> a != "--rebuild", ARGS)
    main(
        get(args, 1, "data/sitrep_pdfs"),
        get(args, 2, "data/onset_curve_scanned.csv");
        rebuild = "--rebuild" in ARGS
    )
end
