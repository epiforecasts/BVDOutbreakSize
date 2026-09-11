## Helpers shared by the `@testitem`s in `test/test_onset_digitiser.jl`.
## Each test item runs in its own module, so this file is `include`d rather
## than relying on file-level definitions.

using Dates: Date, Day

"""
Read the digitised onset CSV at `path` into `(; header, order, report_date,
onsets)`. `order` is the SitRep ids in file order, `report_date` maps an id
to its rapportage date, and `onsets` maps an id to `onset_date =>
confirmed_total`.
"""
function _read_onset_csv(path)
    order = String[]
    report_date = Dict{String, Date}()
    onsets = Dict{String, Dict{Date, Int}}()
    header = ""
    for (i, line) in enumerate(eachline(path))
        if i == 1
            header = line
            continue
        end
        isempty(strip(line)) && continue
        f = split(line, ',')
        length(f) == 6 || error("malformed onset row $i: $line")
        sr = String(f[1])
        if !haskey(onsets, sr)
            push!(order, sr)
            report_date[sr] = Date(f[2])
            onsets[sr] = Dict{Date, Int}()
        end
        onsets[sr][Date(f[3])] = parse(Int, f[6])
    end
    return (; header, order, report_date, onsets)
end

"""
L1 distance between two vintages' onset curves with the later block moved
`s` days, summed over the onset dates the two share, optionally restricted
to onset dates before `cut`. This is the section-4b date-alignment check
`data/README.md` describes, in that section's sign convention: a minimum at
`s = +1` means the later vintage's axis reads one day ahead of the earlier
one's.
"""
function _l1_shift(a, b, s; cut = nothing)
    ds = [d for d in keys(a)
          if haskey(b, d - Day(s)) && (cut === nothing || d < cut)]
    isempty(ds) && return typemax(Int)
    return sum(abs(a[d] - b[d - Day(s)]) for d in ds)
end

"""
Command that runs a PEP 723 Python script with its declared dependencies,
or `nothing` when neither `uv` nor a Python carrying Pillow and numpy is on
`PATH`. `uv` is preferred because it fetches the dependencies itself.
"""
function _python_runner()
    if Sys.which("uv") !== nothing
        return `uv run --no-project`
    end
    py = Sys.which("python3")
    py === nothing && return nothing
    ok = success(pipeline(`$py -c "import numpy, PIL"`;
        stdout = devnull, stderr = devnull))
    return ok ? `$py` : nothing
end

"""
Build a synthetic stacked bar chart as `(R, G, B)` matrices in the form
`digitize` expects: a plot border on the count-zero baseline, a y-axis
label strip on an evenly spaced tick grid, weekly x-axis tick marks, and
one bar per day ending on the rightmost tick. `bars` is a vector of
`(alive, dead)` counts drawn left to right. `tick_px` is the y-axis
gridline spacing and `y_step` the count between gridlines, so
pixels-per-count is `tick_px / y_step`; `week_px` is the x-axis tick
spacing, giving `week_px / 7` pixels per day.
"""
function _synthetic_chart(bars; y_step = 20, tick_px = 80, week_px = 30)
    W, H = 800, 400
    R, G, B = fill(255, H, W), fill(255, H, W), fill(255, H, W)
    function paint!(rows, cols, c)
        R[rows, cols] .= c[1]
        G[rows, cols] .= c[2]
        B[rows, cols] .= c[3]
        return nothing
    end
    grey, ink = (100, 100, 100), (60, 60, 60)
    blue, crimson = (100, 180, 230), (200, 60, 60)
    base, lastx = 300, 693
    ## Plot border: by far the longest contiguous near-grey run in the
    ## image, which is what `baseline_row` ranks rows on.
    paint!(base:base, 40:780, grey)
    ## Vertical axis line, and the label strip just left of it that
    ## `y_axis_ticks` reads the count scale from.
    ticks = [base - k * tick_px for k in ((base - 60) ÷ tick_px):-1:0]
    paint!(ticks[1]:base, 38:38, grey)
    for y in ticks
        paint!(y:y, 30:34, ink)
    end
    ## Weekly x-axis tick marks, in the band just below the baseline.
    for k in 0:19
        x = lastx - k * week_px
        paint!((base + 2):(base + 4), x:x, grey)
    end
    ## One bar per day, each filling its own day slot so the
    ## 75th-percentile window sees a single flat height.
    ppc, ppd, n = tick_px / y_step, week_px / 7, length(bars)
    for (j, (alive, dead)) in enumerate(bars)
        cx = lastx + (j - n) * ppd
        cols = (round(Int, cx) - 2):(round(Int, cx) + 2)
        ha, hd = round(Int, alive * ppc), round(Int, dead * ppc)
        ha > 0 && paint!((base - ha):(base - 1), cols, blue)
        hd > 0 && paint!((base - ha - hd):(base - ha - 1), cols, crimson)
    end
    return R, G, B
end
