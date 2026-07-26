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
#   * date scale = the weekly x-axis tick marks, anchored on the RIGHTMOST
#     tick (whose date is in CONFIG, read off the axis) stepping back 7 days;
#   * each daily bar height = the 75th-percentile column in a one-day window,
#     flooded up from the baseline counting light-blue (Vivant) and crimson
#     (Decede) pixels, bridging the anti-alias gap between the stacked
#     segments but stopping at the wide white gap up to the floating label /
#     dashed line above the bar.
#
# Accuracy: digitised per-vintage totals run 2-5% below the printed figure `n`
# (SitRep 064: 2018 vs printed n=2 064; SitRep 070: 2260 vs n=2 329); individual
# bars carry roughly +/-1-2 cases of pixel noise. The shortfall sits in the
# faded bars of the `donnees potentiellement incompletes` band, whose lightened
# fill falls outside the colour masks. The values are approximate and are NOT
# fitted by the model; they are captured for later use (see data/README.md).
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
    ("069", Date(2026, 7, 22), Date(2026, 7, 22)),
    ("070", Date(2026, 7, 23), Date(2026, 7, 22))
]

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
# steel blue and prints value labels). The blue floor is well below the
# smallest onset figure seen (SitReps 069/070 embed it at 1009x583, blue
# ~50k, against ~66-95k for the larger 059-067 renderings) and well above
# the largest non-onset image sharing its page (~5k).
function is_onset_curve(R, G, B)
    m = masks(R, G, B)
    orange = (R .> 200) .& (G .> 110) .& (G .< 195) .& (B .< 90)
    return sum(m.blue) > 20000 && sum(orange) < 500 && sum(m.red) > 20000
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

function digitize(R, G, B, last_tick::Date)
    H, W = size(R)
    m = masks(R, G, B)
    blue, red, dark = m.blue, m.red, m.dark
    drow = vec(sum(dark; dims = 2))
    drow[1:floor(Int, H * 0.4)] .= 0
    base = argmax(drow)                   # count-0 baseline row
    # count scale from the 0/20/40/60 y-axis ticks
    yt = y_axis_ticks(dark, base, H, W)
    ppc = median(diff(yt)) / 20.0         # pixels per count
    ytop, y0 = yt[1], yt[end]
    # x scale from the weekly tick marks below the baseline. The tick marks
    # are only a few pixels tall and shrink with the embedded figure
    # resolution (5-6 dark rows in the 1257x698 SitRep 064 rendering, 4 in
    # SitRep 066's 1275x623, 3 in SitRep 069/070's 1009x583), so step the
    # cut down until a full weekly row of ticks resolves instead of fixing
    # it at 4 and losing the axis entirely on the smaller figures.
    band = vec(sum(dark[(base + 2):min(base + 8, H), :]; dims = 1))
    xt = Int[]
    for cut in (4, 3, 2)
        xt = cluster([x for x in 1:W if band[x] >= cut])
        length(xt) >= 8 && break
    end
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

# Extract the onset-curve figure from a SitRep PDF as R, G, B matrices.
function onset_image(pdf)
    npages = parse(Int, match(r"Pages:\s*(\d+)",
        read(`pdfinfo $pdf`, String)).captures[1])
    page = nothing
    for p in 1:npages
        txt = lowercase(read(`pdftotext -layout -f $p -l $p $pdf -`, String))
        if occursin("date de debut des symptom", txt) ||
           occursin("date de début des symptôm", txt)
            page = p
            break
        end
    end
    page === nothing && return nothing
    best = nothing
    mktempdir() do wd
        # pdfimages writes PPM (P6) for RGB images by default; no format flag
        run(`pdfimages -f $page -l $page $pdf $(joinpath(wd, "p"))`)
        for name in sort(readdir(wd))
            endswith(name, ".ppm") || continue
            R, G, B = read_ppm(joinpath(wd, name))
            if is_onset_curve(R, G, B) &&
               (best === nothing || length(R) > length(best[1]))
                best = (R, G, B)
            end
        end
    end
    return best
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
            rows = digitize(img..., last_tick)
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
