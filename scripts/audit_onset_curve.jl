#!/usr/bin/env julia
#
# Audit the digitised symptom-onset curve (data/onset_curve_scanned.csv)
# against the figures it was read from.
#
# For every vintage in the scanned file this writes one row to
# data/onset_curve_figures.csv: the PDF and page the figure came from, the
# embedded image's size and md5, the pixel scales the digitiser calibrated
# (pixels per count and per day), the CONFIG last tick, the `n` the figure
# prints in its title, the digitised total and the gap between the two, and
# the earlier vintage whose digitised block it reprints, if any.
#
# The printed n is burned into the raster (the chart title reads "... des
# symptomes (n = 5 263)" and the source line repeats it), so it is not in
# the PDF text layer for any vintage. The text layer of the figure page is
# tried first; when that finds nothing, the title strip and the source
# strip of the embedded image are upscaled and read with tesseract. The
# title reading is the value and the source reading corroborates it.
# `printed_n_source` records which path produced the value; `note` says
# why one is missing or where the two strips disagreed.
#
# The audit itself (printed, and written to output/onset_curve_audit.md)
# holds the gap table, and, for each consecutive pair of distinct
# snapshots, what happened to the settled bars: onset dates more than 56
# days before the later report date, where late reporting has finished and
# a bar can only rise. It counts the bars that fell, the net change, and the
# L1 distance between the two blocks at day shifts -2..+2. A minimum away
# from 0 flags a date misalignment between the two vintages' axes.
#
# Dependencies: poppler (`pdfimages`, `pdftotext`, `pdfinfo`) on PATH, and
# tesseract for the printed n (without it the column is left empty).
# `md5sum` (or `md5` on macOS) hashes the extracted image. No Julia
# packages beyond stdlib. The digitiser's own helpers are `include`d from
# scripts/digitize_onset_curve.jl so the calibration reported here is the
# one the digitiser used.
#
# Usage:
#   julia scripts/audit_onset_curve.jl [pdf_dir] [scanned_csv] [figures_csv]
#       [audit_md] [--crops=DIR]
# `--crops=DIR` keeps the upscaled title and source strips each n was read
# from, as `<sitrep>_title.ppm` and `<sitrep>_source.ppm`, for a second
# reader. A full run over 64 vintages takes about a minute.
# Defaults: data/sitrep_pdfs, data/onset_curve_scanned.csv,
#           data/onset_curve_figures.csv, output/onset_curve_audit.md

using Dates: Date, Day, value
using Statistics: median
using Printf: @printf, @sprintf

include(joinpath(@__DIR__, "digitize_onset_curve.jl"))

# --- printed n -------------------------------------------------------------

# A count printed as "n = 5 263" (thousands separated by a space, a thin
# space or a full stop) or "n=1772". Only lines mentioning the symptom-onset
# basis are read, because the report text also prints "n = ..." for
# provinces and deaths.
const PRINTED_N_RE = Regex(
    "n\\s*[=:]\\s*(\\d{1,3}(?:[ \u202f\u00a0\u2009.]\\d{3})+|\\d+)", "i"
)

"""
Every printed n found on a line of `text` that mentions the symptom-onset
basis ("sympt"), in order of appearance. The digits are returned with
their separators removed.
"""
function parse_printed_n(text::AbstractString)
    out = Int[]
    for line in split(text, '\n')
        occursin("sympt", lowercase(line)) || continue
        for m in eachmatch(PRINTED_N_RE, line)
            digits = filter(isdigit, m.captures[1])
            isempty(digits) || push!(out, parse(Int, digits))
        end
    end
    return out
end

# Write `rows` of an RGB image to a P6 PPM at `path`, scaled up by an
# integer factor with nearest-neighbour replication. Tesseract reads the
# small chart text far better at 3x than at native size.
function write_ppm_crop(path, R, G, B, rows; scale = 3)
    H, W = length(rows), size(R, 2)
    open(path, "w") do io
        write(io, "P6\n$(W * scale) $(H * scale)\n255\n")
        buf = Vector{UInt8}(undef, 3 * W * scale)
        for r in rows, _ in 1:scale
            p = 1
            for c in 1:W, _ in 1:scale
                buf[p] = UInt8(R[r, c])
                buf[p + 1] = UInt8(G[r, c])
                buf[p + 2] = UInt8(B[r, c])
                p += 3
            end
            write(io, buf)
        end
    end
    return path
end

# OCR of `rows` of the image. With `keep` set, the upscaled crop is left at
# that path so a second reader can check the digits by eye.
function ocr_rows(R, G, B, rows; keep = nothing)
    return mktempdir() do wd
        path = keep === nothing ? joinpath(wd, "crop.ppm") : keep
        img = write_ppm_crop(path, R, G, B, rows)
        read(pipeline(`tesseract $img stdout --psm 6`; stderr = devnull), String)
    end
end

"""
The n printed on the onset figure of `pdf`, as `(n, source, note)`. `n` is
`nothing` when no reading was found; `source` names where it came from.
"""
function printed_n(pdf, page, R, G, B; crop_dir = nothing, sr = "")
    txt = read(`pdftotext -layout -f $page -l $page $pdf -`, String)
    found = parse_printed_n(txt)
    isempty(found) || return (found[1], "text layer", "")
    if Sys.which("tesseract") === nothing
        return (nothing, "", "not in the text layer; tesseract not on PATH")
    end
    H = size(R, 1)
    crop(name) = crop_dir === nothing ? nothing :
        joinpath(crop_dir, "$(sr)_$(name).ppm")
    title = parse_printed_n(
        ocr_rows(R, G, B, 1:round(Int, 0.12H); keep = crop("title"))
    )
    source = parse_printed_n(
        ocr_rows(R, G, B, round(Int, 0.86H):H; keep = crop("source"))
    )
    if isempty(title) && isempty(source)
        return (
            nothing, "", "not in the text layer; OCR of the title " *
                "and source strips found no n",
        )
    end
    ## The title is set in a larger bold face and reads reliably; the
    ## source line is small italic and drops or swaps digits, so it only
    ## corroborates.
    if isempty(title)
        return (source[1], "ocr source", "title strip gave no n")
    end
    n = title[1]
    if !isempty(source) && source[1] != n
        return (n, "ocr title", "source strip read $(source[1])")
    end
    return (n, isempty(source) ? "ocr title" : "ocr title+source", "")
end

# --- figure location and calibration --------------------------------------

function file_md5(path)
    tool = Sys.which("md5sum")
    tool === nothing && (tool = Sys.which("md5"))
    tool === nothing && return ""
    cmd = basename(tool) == "md5" ? `$tool -q $path` : `$tool $path`
    out = read(cmd, String)
    return String(first(split(strip(out))))
end

# The largest onset-curve image among those embedded on `page`, with its
# md5 and size, or `nothing`.
function onset_figure_on_page(pdf, page)
    return mktempdir() do wd
        run(`pdfimages -f $page -l $page $pdf $(joinpath(wd, "p"))`)
        best = nothing
        for name in sort(readdir(wd))
            endswith(name, ".ppm") || continue
            path = joinpath(wd, name)
            R, G, B = read_ppm(path)
            if is_onset_curve(R, G, B) &&
                    (best === nothing || length(R) > length(best.R))
                best = (; R, G, B, md5 = file_md5(path), page)
            end
        end
        best
    end
end

# Same page search as the digitiser's `onset_image`: the caption page, then
# its immediate neighbours.
function onset_figure(pdf)
    page, npages = _onset_page(pdf)
    page === nothing && return nothing
    for q in (page, page - 1, page + 1)
        1 <= q <= npages || continue
        fig = onset_figure_on_page(pdf, q)
        fig === nothing || return fig
    end
    return nothing
end

# The weekly x-axis tick columns, as `digitize` finds them.
function weekly_tick_columns(mask, base, H, W)
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

"""
Least-squares pixels per week through marks known to sit at a weekly
pitch, given their columns `xs` (ascending) and a first guess `week` at
the pitch. Each mark is indexed by rounding its distance from the last one
to whole weeks, so a missed mark leaves a gap of two or more and a stray
well off the grid is dropped. The pitch is refined and the indexing
repeated, so marks far from the last one are not lost to the guess's
rounding. Returns `(pixels_per_week, points, weeks)`, `weeks` being how
many weeks before the last mark the earliest indexed mark sits.
"""
function fit_week(xs, week)
    fit = week
    idx = Int[]
    cols = Float64[]
    for _ in 1:4
        empty!(idx)
        empty!(cols)
        for x in reverse(xs)
            r = (xs[end] - x) / fit
            k = round(Int, r)
            (abs(r - k) > 0.25 || k in idx) && continue
            push!(idx, k)
            push!(cols, x)
        end
        length(idx) < 3 && return (NaN, length(idx), 0)
        kbar, xbar = sum(idx) / length(idx), sum(cols) / length(cols)
        fit = -sum((idx .- kbar) .* (cols .- xbar)) / sum((idx .- kbar) .^ 2)
    end
    return (fit, length(idx), maximum(idx))
end

"""
Column centroids of the first x-axis label line below the baseline. Each
label is centred under its tick, so the centroids run at the weekly pitch
and give a second, usually more complete, set of marks than the tick marks
themselves, which the small renders draw only a pixel tall.
"""
function label_centroids(R, G, B, base)
    H, W = size(R)
    dark = (R .< 150) .& (G .< 150) .& (B .< 150)
    rows = (base + 7):min(base + 40, H)
    on = [r for r in rows if sum(@view dark[r, :]) >= 5]
    isempty(on) && return Float64[]
    r0 = r1 = on[1]
    while r1 < H && sum(@view dark[r1 + 1, :]) >= 5
        r1 += 1
    end
    band = vec(sum(dark[r0:r1, :]; dims = 1))
    out = Float64[]
    cl = Int[]
    for x in 1:W
        band[x] > 0 || continue
        if !isempty(cl) && x - cl[end] <= 4
            push!(cl, x)
        else
            isempty(cl) || push!(out, sum(cl) / length(cl))
            cl = [x]
        end
    end
    isempty(cl) || push!(out, sum(cl) / length(cl))
    return out
end

"""
The pixel scales `digitize` calibrates on this figure: pixels per count
from the y-axis tick spacing over `y_step`, and pixels per day from the
weekly x-axis ticks, with the same near-grey fallbacks the digitiser uses.

`pixels_per_day` is the digitiser's value, the median tick spacing over 7.
Tick columns are integers, so on a figure whose true spacing is not a whole
number of pixels the median lands on one of the two neighbouring integers
and the day grid stepped back from the last tick drifts away from the
drawn bars. `pixels_per_day_fit` is the least-squares pitch through the
weekly marks (`fit_source` says whether the tick marks or the x-axis
labels gave more points, `fit_points` how many), and `drift_days` is how
many days the digitiser's grid is off at `span_days` before the last tick.
"""
function calibration(R, G, B, y_step; span_days = 100)
    H, W = size(R)
    m = masks(R, G, B)
    line = (R .< 180) .& (G .< 180) .& (B .< 180)
    base = baseline_row(R, G, B, H)
    yt = try
        y_axis_ticks(m.dark, base, H, W)
    catch e
        e isa ErrorException || rethrow()
        y_axis_ticks(line, base, H, W)
    end
    xt = weekly_tick_columns(m.dark, base, H, W)
    isempty(xt) && (xt = weekly_tick_columns(line, base, H, W))
    isempty(xt) && error("no x-axis weekly tick row found")
    week = median(diff(xt))
    pixels_per_day = week / 7
    ticks = fit_week(xt, week)
    labels = fit_week(label_centroids(R, G, B, base), week)
    (fit, fit_points, weeks_seen), fit_source = labels[2] > ticks[2] ?
        (labels, "labels") : (ticks, "ticks")
    pixels_per_day_fit = fit / 7
    drift_days = span_days * (pixels_per_day_fit - pixels_per_day) /
        pixels_per_day_fit
    return (;
        pixels_per_count = median(diff(yt)) / y_step,
        pixels_per_day, pixels_per_day_fit, fit_source, fit_points,
        weeks_seen, drift_days,
    )
end

# --- scanned CSV -----------------------------------------------------------

"""
Read the digitised onset CSV into `(; order, report_date, blocks)`: the
SitRep ids in file order, each one's report date, and each one's rows as
`onset_date => (alive, dead)` in a `Dict`.
"""
function read_scanned(path)
    order = String[]
    report_date = Dict{String, Date}()
    blocks = Dict{String, Dict{Date, Tuple{Int, Int}}}()
    for (i, line) in enumerate(eachline(path))
        (i == 1 || isempty(strip(line))) && continue
        f = split(line, ',')
        length(f) == 6 || error("malformed onset row $i: $line")
        sr = String(f[1])
        if !haskey(blocks, sr)
            push!(order, sr)
            report_date[sr] = Date(f[2])
            blocks[sr] = Dict{Date, Tuple{Int, Int}}()
        end
        blocks[sr][Date(f[3])] = (parse(Int, f[4]), parse(Int, f[5]))
    end
    return (; order, report_date, blocks)
end

block_total(block) = sum(a + d for (a, d) in values(block); init = 0)

"""
For each vintage, the earliest earlier vintage whose digitised block is
identical, else `""`. This is the loader's reprint rule: exact equality
over the whole block, not a list of ids.
"""
function reprints(scanned)
    out = Dict{String, String}()
    for (i, sr) in enumerate(scanned.order)
        out[sr] = ""
        for j in 1:(i - 1)
            p = scanned.order[j]
            if scanned.blocks[p] == scanned.blocks[sr]
                out[sr] = isempty(out[p]) ? p : out[p]
                break
            end
        end
    end
    return out
end

# --- figures table ---------------------------------------------------------

const FIGURES_HEADER = join(
    (
        "sitrep", "report_date", "pdf", "page", "width", "height",
        "image_md5", "y_axis_step", "pixels_per_count", "pixels_per_day",
        "pixels_per_day_fit", "fit_source", "fit_points", "drift_days",
        "last_tick",
        "printed_n", "printed_n_source", "digitised_total",
        "gap_pct", "reprint_of", "axis_start", "loop_start", "first_onset",
        "cases_before_first_onset", "note",
    ), ","
)

fmt(x::Nothing) = ""
fmt(x::AbstractFloat) = @sprintf("%.4f", x)
fmt(x) = string(x)

# The digitiser walks days from 105 before the last tick, so a block cannot
# start earlier than that whatever the axis shows. `peak` is the largest
# count any vintage read for each onset date, a floor on what the loop
# start leaves unread.
function figure_row(
        sr, scanned, config, reprint_of, pdf_dir, peak;
        crop_dir = nothing
    )
    report_date = scanned.report_date[sr]
    block = scanned.blocks[sr]
    total = block_total(block)
    last_tick = haskey(config, sr) ? config[sr] : nothing
    first_onset = isempty(block) ? nothing : minimum(keys(block))
    loop_start = last_tick === nothing ? nothing : last_tick - Day(105)
    y_step = get(Y_AXIS_STEP, sr, 20)
    row = Dict{String, Any}(
        "sitrep" => sr, "report_date" => report_date,
        "y_axis_step" => y_step, "last_tick" => last_tick,
        "digitised_total" => total, "reprint_of" => reprint_of,
        "pdf" => "", "page" => nothing, "width" => nothing,
        "height" => nothing, "image_md5" => "",
        "pixels_per_count" => nothing, "pixels_per_day" => nothing,
        "pixels_per_day_fit" => nothing, "fit_source" => "",
        "fit_points" => nothing, "drift_days" => nothing,
        "printed_n" => nothing, "printed_n_source" => "",
        "gap_pct" => nothing, "axis_start" => nothing,
        "loop_start" => loop_start, "first_onset" => first_onset,
        "cases_before_first_onset" => nothing, "note" => "",
    )
    name = "SitRep_MVE_$(sr)_2026.pdf"
    pdf = joinpath(pdf_dir, name)
    if !isfile(pdf)
        row["note"] = "PDF not found"
        return row
    end
    row["pdf"] = name
    fig = onset_figure(pdf)
    if fig === nothing
        row["note"] = "no onset figure found"
        return row
    end
    row["page"] = fig.page
    row["height"], row["width"] = size(fig.R)
    row["image_md5"] = fig.md5
    try
        cal = calibration(fig.R, fig.G, fig.B, y_step)
        row["pixels_per_count"] = cal.pixels_per_count
        row["pixels_per_day"] = cal.pixels_per_day
        row["pixels_per_day_fit"] = cal.pixels_per_day_fit
        row["fit_source"] = cal.fit_source
        row["fit_points"] = cal.fit_points
        row["drift_days"] = cal.drift_days
        if last_tick !== nothing && first_onset !== nothing
            axis_start = last_tick - Day(7 * cal.weeks_seen)
            row["axis_start"] = axis_start
            row["cases_before_first_onset"] = sum(
                (n for (d, n) in peak if axis_start <= d < first_onset);
                init = 0
            )
        end
    catch e
        e isa ErrorException || rethrow()
        row["note"] = "calibration failed: " * e.msg
    end
    n, src, note = printed_n(
        pdf, fig.page, fig.R, fig.G, fig.B; crop_dir, sr
    )
    row["printed_n"] = n
    row["printed_n_source"] = src
    if n !== nothing
        row["gap_pct"] = 100 * (total - n) / n
    end
    isempty(note) || (
        row["note"] = isempty(row["note"]) ? note :
            row["note"] * "; " * note
    )
    return row
end

function write_figures(path, rows)
    cols = split(FIGURES_HEADER, ',')
    open(path, "w") do io
        println(io, FIGURES_HEADER)
        for row in rows
            println(io, join((fmt(row[c]) for c in cols), ','))
        end
    end
    return path
end

function read_figures(path)
    lines = readlines(path)
    cols = split(lines[1], ',')
    return [
        Dict(zip(cols, split(line, ',')))
            for line in lines[2:end] if !isempty(strip(line))
    ]
end

# --- settled-bar audit -----------------------------------------------------

# L1 distance between blocks `a` (earlier) and `b` (later) with `b` moved
# `s` days, over the onset dates of `a` before `cut` that `b` also carries.
# The sign convention is the one test/onset_digitiser_helpers.jl uses: a
# minimum at s = +1 means the later axis reads a day ahead of the earlier.
function l1_shift(a, b, s, cut)
    tot = 0
    n = 0
    for (d, (alive, dead)) in a
        d < cut || continue
        haskey(b, d - Day(s)) || continue
        bl, bd = b[d - Day(s)]
        tot += abs(alive + dead - bl - bd)
        n += 1
    end
    return n == 0 ? nothing : tot
end

const SHIFTS = -2:2

"""
For each consecutive pair of distinct snapshots (reprints collapsed onto
their first vintage), what the settled bars did: onset dates more than
`settle` days before the later report date that both blocks carry.
"""
function settled_pairs(scanned, reprint_of; settle = 56)
    distinct = [sr for sr in scanned.order if isempty(reprint_of[sr])]
    out = []
    for i in 2:length(distinct)
        p, q = distinct[i - 1], distinct[i]
        a, b = scanned.blocks[p], scanned.blocks[q]
        cut = scanned.report_date[q] - Day(settle)
        shared = [d for d in keys(a) if d < cut && haskey(b, d)]
        fell = count(d -> sum(b[d]) < sum(a[d]), shared)
        rose = count(d -> sum(b[d]) > sum(a[d]), shared)
        net = sum(sum(b[d]) - sum(a[d]) for d in shared; init = 0)
        l1 = [l1_shift(a, b, s, cut) for s in SHIFTS]
        present = [s for (s, v) in zip(SHIFTS, l1) if v !== nothing]
        best = isempty(present) ? nothing :
            present[argmin([v for v in l1 if v !== nothing])]
        push!(
            out, (;
                from = p, to = q, cut, bars = length(shared), fell, rose,
                net, l1, best,
            )
        )
    end
    return out
end

# --- report ----------------------------------------------------------------

function md_table(io, header, rows)
    println(io, "| ", join(header, " | "), " |")
    println(io, "|", join(fill(" --- ", length(header)), "|"), "|")
    for r in rows
        println(io, "| ", join(string.(r), " | "), " |")
    end
    return println(io)
end

pct(x) = x === nothing ? "" : @sprintf("%+.1f", x)

function write_audit(io, rows, pairs)
    println(io, "# Onset curve audit\n")
    println(
        io, "Digitised totals against the n each figure prints. ",
        "gap_pct is 100 (digitised - printed) / printed.\n"
    )
    md_table(
        io,
        (
            "sitrep", "report", "printed n", "digitised", "gap %", "source",
            "reprint of", "note",
        ),
        [
            (
                r["sitrep"], r["report_date"], fmt(r["printed_n"]),
                r["digitised_total"], pct(r["gap_pct"]),
                r["printed_n_source"], r["reprint_of"], r["note"],
            )
                for r in rows
        ]
    )
    with_n = [r for r in rows if r["printed_n"] !== nothing]
    println(io, "Vintages with a printed n: $(length(with_n)) of $(length(rows)).")
    missing_n = [r["sitrep"] for r in rows if r["printed_n"] === nothing]
    isempty(missing_n) ||
        println(io, "Without one: ", join(missing_n, ", "), ".")
    println(io)

    println(io, "## Settled bars between consecutive distinct snapshots\n")
    println(
        io, "Settled bars are onset dates more than 56 days before the ",
        "later report date that both blocks carry. A settled bar can only ",
        "rise, so every fall is scan noise. The L1 columns are the ",
        "distance between the two blocks with the later one moved by ",
        "that many days; the minimum should sit at shift 0.\n"
    )
    md_table(
        io,
        (
            "pair", "settled bars", "fell", "rose", "net", "L1 -2", "L1 -1",
            "L1 0", "L1 +1", "L1 +2", "best shift",
        ),
        [
            (
                "$(p.from)->$(p.to)", p.bars, p.fell, p.rose, p.net,
                (v === nothing ? "" : v for v in p.l1)..., fmt(p.best),
            )
                for p in pairs
        ]
    )
    off = [p for p in pairs if p.best !== nothing && p.best != 0]
    println(
        io, "Pairs: $(length(pairs)). With a fall on at least one settled ",
        "bar: $(count(p -> p.fell > 0, pairs)). Net negative: ",
        "$(count(p -> p.net < 0, pairs)). L1 minimum away from shift 0: ",
        "$(length(off))",
        isempty(off) ? "." :
            " (" * join(("$(p.from)->$(p.to) at $(p.best)" for p in off), ", ") * ").",
    )
    println(io)

    println(io, "## Axis coverage\n")
    println(
        io, "axis start is the earliest weekly mark found on the axis; ",
        "loop start is 105 days before the last tick, where the ",
        "digitiser's day loop begins; first onset is the first date in ",
        "the committed block. cases before is the largest count any ",
        "vintage read on the dates between axis start and first onset, a ",
        "floor on what the block leaves unread; gap+ adds them back.\n"
    )
    with_axis = [r for r in rows if r["axis_start"] !== nothing]
    md_table(
        io,
        (
            "sitrep", "axis start", "loop start", "first onset", "cases before",
            "gap %", "gap+ %",
        ),
        [
            (
                r["sitrep"], r["axis_start"], r["loop_start"],
                r["first_onset"], r["cases_before_first_onset"],
                pct(r["gap_pct"]),
                r["printed_n"] === nothing ? "" : pct(
                        100 * (
                            r["digitised_total"] +
                            r["cases_before_first_onset"] - r["printed_n"]
                        ) / r["printed_n"]
                    ),
            )
                for r in with_axis
        ]
    )

    println(io, "## Day-scale drift\n")
    println(
        io, "pixels_per_day is the digitiser's median tick spacing over 7; ",
        "fit is the least-squares pitch through the weekly tick marks or the ",
        "x-axis labels, whichever gives more points. ",
        "drift_days is how far the digitiser's day grid is from the drawn ",
        "bars 100 days before the last tick, in days: past 0.5 the oldest ",
        "bars are read into neighbouring days.\n"
    )
    with_cal = [r for r in rows if r["drift_days"] !== nothing]
    md_table(
        io,
        (
            "sitrep", "width", "pixels/day", "fit", "fit from", "points",
            "drift days",
        ),
        [
            (
                r["sitrep"], r["width"], fmt(r["pixels_per_day"]),
                fmt(r["pixels_per_day_fit"]), r["fit_source"],
                r["fit_points"], @sprintf("%+.2f", r["drift_days"]),
            )
                for r in with_cal
        ]
    )
    println(
        io, "Vintages with |drift| of at least half a day at 100 days: ",
        "$(count(r -> abs(r["drift_days"]) >= 0.5, with_cal)) of ",
        "$(length(with_cal)).\n"
    )

    println(io, "## Worst ten by |gap|\n")
    worst = sort(with_n; by = r -> -abs(r["gap_pct"]))
    md_table(
        io, ("sitrep", "printed n", "digitised", "gap %"),
        [
            (r["sitrep"], r["printed_n"], r["digitised_total"], pct(r["gap_pct"]))
                for r in worst[1:min(10, end)]
        ]
    )
    println(io, "## Worst ten pairs by settled-bar L1 at shift 0\n")
    at0 = [p for p in pairs if p.l1[3] !== nothing]
    worst_l1 = sort(at0; by = p -> -p.l1[3])
    md_table(
        io, ("pair", "settled bars", "L1 0", "fell", "net"),
        [
            ("$(p.from)->$(p.to)", p.bars, p.l1[3], p.fell, p.net)
                for p in worst_l1[1:min(10, end)]
        ]
    )
    return nothing
end

function main(
        pdf_dir = "data/sitrep_pdfs",
        scanned_csv = "data/onset_curve_scanned.csv",
        figures_csv = "data/onset_curve_figures.csv",
        audit_md = "output/onset_curve_audit.md";
        crop_dir = nothing
    )
    scanned = read_scanned(scanned_csv)
    reprint_of = reprints(scanned)
    config = Dict(sr => last_tick for (sr, _, last_tick) in CONFIG)
    peak = Dict{Date, Int}()
    for block in values(scanned.blocks), (d, (a, dd)) in block
        peak[d] = max(get(peak, d, 0), a + dd)
    end
    crop_dir === nothing || mkpath(crop_dir)
    rows = Dict{String, Any}[]
    for sr in scanned.order
        row = figure_row(
            sr, scanned, config, reprint_of[sr], pdf_dir, peak; crop_dir
        )
        push!(rows, row)
        @printf(
            "%s %s n=%-6s digitised=%-5d gap=%-7s %s %s\n", sr,
            row["report_date"], fmt(row["printed_n"]),
            row["digitised_total"], pct(row["gap_pct"]),
            row["printed_n_source"], row["note"]
        )
    end
    write_figures(figures_csv, rows)
    println("wrote $figures_csv")
    pairs = settled_pairs(scanned, reprint_of)
    buf = IOBuffer()
    write_audit(buf, rows, pairs)
    report = String(take!(buf))
    print(report)
    mkpath(dirname(audit_md))
    write(audit_md, report)
    println("wrote $audit_md")
    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    crops = filter(a -> startswith(a, "--crops="), ARGS)
    main(
        filter(a -> !startswith(a, "--crops="), ARGS)...;
        crop_dir = isempty(crops) ? nothing : last(split(crops[end], "=", limit = 2))
    )
end
