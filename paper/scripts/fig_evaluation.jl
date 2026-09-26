# Figure: evaluation. Forecast skill by release, outbreak size by stream,
# and the comparisons with McCabe et al. and Chamla et al.
#
# Run from the worktree root:
#
#     julia --project=. paper/scripts/fig_evaluation.jl
#
# Reads data/forecast_scores.csv (the cross-release scores the scoring
# script writes), the latest results release's stream_estimates.csv,
# frozen_matched_cutoffs.csv, observations.toml and site.zip under
# paper/data/release/ (git-ignored; paper/scripts/release_assets.jl
# downloads them), and the published scenario constants in
# src/constants.jl. Writes paper/figures/fig-evaluation.pdf and
# fig-evaluation.png (180 mm wide, 600 dpi). Colours are Makie's Wong
# palette, as in src/plots.jl and fig_journey.jl.

using CairoMakie
using Dates
using TOML

const REPO = normpath(joinpath(@__DIR__, "..", ".."))
const DATA = joinpath(REPO, "paper", "data")
const RELEASE = joinpath(DATA, "release")
const FIGS = joinpath(REPO, "paper", "figures")
mkpath(FIGS)

# ---------------------------------------------------------------------
# CSV reading (quoted fields, no embedded newlines)
# ---------------------------------------------------------------------
function split_csv_line(l)
    out = String[]
    buf = IOBuffer()
    inq = false
    i = firstindex(l)
    while i <= lastindex(l)
        c = l[i]
        if inq
            if c == '"' && i < lastindex(l) && l[nextind(l, i)] == '"'
                write(buf, '"'); i = nextind(l, i)
            elseif c == '"'
                inq = false
            else
                write(buf, c)
            end
        elseif c == '"'
            inq = true
        elseif c == ','
            push!(out, String(take!(buf)))
        else
            write(buf, c)
        end
        i = nextind(l, i)
    end
    push!(out, String(take!(buf)))
    return out
end
function read_csv(path)
    lines = filter(!isempty, split(read(path, String), '\n'))
    header = Symbol.(split_csv_line(lines[1]))
    return [NamedTuple{Tuple(header)}(Tuple(split_csv_line(l))) for l in lines[2:end]]
end
num(s) = parse(Float64, s)

# ---------------------------------------------------------------------
# Inputs
# ---------------------------------------------------------------------
scores = read_csv(joinpath(REPO, "data", "forecast_scores.csv"))
estimates = read_csv(joinpath(RELEASE, "stream_estimates.csv"))
frozen = read_csv(joinpath(RELEASE, "frozen_matched_cutoffs.csv"))
obs = TOML.parsefile(joinpath(RELEASE, "observations.toml"))
const CUTOFF = Date(obs["as_of_date"])

# The published scenario tables, evaluated from the committed source so
# the figure cannot drift from the constants the report uses.
function source_const(name)
    src = read(joinpath(REPO, "src", "constants.jl"), String)
    m = match(Regex("const $name = (\\[.*?\\n\\])", "s"), src)
    m === nothing && error("no `const $name` in src/constants.jl")
    return eval(Meta.parse(m.captures[1]))
end
mccabe = source_const("REPORT_SCENARIOS_CI")
chamla_central = source_const("CHAMLA_CONFIRMED_CENTRAL")
chamla_w12 = source_const("CHAMLA_CONFIRMED_W12")

# Our confirmed-case projection from the 8 June frozen fit at Chamla's
# report dates, from the table the site's sensitivity page prints. The
# projection is not published as a CSV asset, so it is read off the
# rendered page inside site.zip.
function chamla_projection_rows()
    zip = joinpath(RELEASE, "site.zip")
    isfile(zip) || error("$zip is missing; see paper/scripts/release_assets.jl")
    html = read(`unzip -p $zip BVDOutbreakSize/dev/sensitivity.html`, String)
    tables = collect(eachmatch(r"<table.*?</table>"s, html))
    t = findfirst(m -> occursin("Chamla central (90% PI)", m.match), tables)
    t === nothing && error("no Chamla comparison table in sensitivity.html")
    strip_tags(s) = replace(
        s, r"<[^>]+>" => "", "&ndash;" => "–",
        "&#8211;" => "–", "&#40;" => "(", "&#41;" => ")"
    )
    rows = map(eachmatch(r"<tr.*?</tr>"s, tables[t].match)) do r
        [
            strip_tags(c.captures[1])
                for c in eachmatch(r"<t[hd][^>]*>(.*?)</t[hd]>"s, r.match)
        ]
    end
    header = rows[1]
    idate = findfirst(==("Date"), header)
    iours = findfirst(==("Our projection (90% CrI)"), header)
    return map(rows[2:end]) do r
        m = match(r"(\d+) \((\d+)–(\d+)\)", r[iours])
        d = Date(r[idate] * " 2026", dateformat"d U yyyy")
        (
            date = d, median = num(m.captures[1]), lo = num(m.captures[2]),
            hi = num(m.captures[3]),
        )
    end
end
ours_chamla = chamla_projection_rows()

# ---------------------------------------------------------------------
# Colours and sizes
# ---------------------------------------------------------------------
wong = Makie.wong_colors()
const C_JOINT = wong[6]      # vermillion, our joint model throughout
const C_SINGLE = :grey45
const C_M1 = wong[1]         # McCabe geographic spread
const C_M2 = wong[2]         # McCabe back-calculation from deaths
const C_CHAMLA = wong[1]
const MM = 72 / 25.4
const FS = 7.0
const FS_SMALL = 5.5

# Streams the situation reports still update, in the order the panels
# list them, with the colour each keeps in panel A.
const STREAMS = [
    ("confirmed cases", wong[1]), ("confirmed deaths", wong[2]),
    ("isolation beds", wong[3]), ("onset reports", wong[4]),
    ("recovered", wong[5]),
]
# Single-stream fit ids of stream_estimates.csv, as the paper names them.
const FIT_LABELS = [
    ("confirmed", "confirmed cases"), ("confirmed_deaths", "confirmed deaths"),
    ("cases", "reported cases"), ("deaths", "reported deaths"),
    ("treatment", "isolation beds"), ("onsets", "onset reports"),
    ("exports", "exports"),
]

fig = Figure(;
    size = (180 * MM, 150 * MM), fontsize = FS,
    figure_padding = (4, 8, 4, 4)
)
axkw = (;
    xgridvisible = false, ygridvisible = false,
    xticklabelsize = FS, yticklabelsize = FS, xlabelsize = FS,
    ylabelsize = FS, spinewidth = 0.6, xtickwidth = 0.6, ytickwidth = 0.6,
)
function panel_letter(pos, s)
    return Label(
        pos, s; font = :bold, fontsize = 9, padding = (0, 6, 6, 0),
        halign = :right, valign = :bottom
    )
end

# ---------------------------------------------------------------------
# A: one-week-ahead relative skill against persistence, by release
# ---------------------------------------------------------------------
# One row per (release, stream): the joint's CRPS over the baseline's for
# the seven-day-ahead forecast, the joint's 90% coverage of that target,
# and whether the release is a reconstructed backfill of a tagged model
# version rather than a forecast published at the time.
const SkillRow = @NamedTuple{
    release::String, date::Date, stream::String,
    skill::Float64, covered::Bool, backfill::Bool,
}
skill = let
    h7 = filter(r -> r.horizon == "7", scores)
    rows = SkillRow[]
    for r in filter(r -> r.fit == "joint", h7)
        b = findfirst(
            s -> s.fit == "baseline" && s.release == r.release &&
                s.stream == r.stream, h7
        )
        b === nothing && continue
        push!(
            rows, (
                release = r.release, date = Date(r.made_date),
                stream = r.stream, skill = num(r.crps) / num(h7[b].crps),
                covered = r.coverage_90 == "true",
                backfill = occursin("(backfill)", r.release),
            )
        )
    end
    sort(rows; by = r -> r.date)
end
const D0 = Date(2026, 6, 1)
dnum(d::Date) = Dates.value(d - D0)
a_dates = unique(r.date for r in skill)
a_xlims = (dnum(minimum(a_dates)) - 3, dnum(maximum(a_dates)) + 3)
a_months = month(minimum(a_dates)):month(maximum(a_dates))
a_ticks = [Date(2026, m, d) for m in a_months for d in (1, 15)]
a_xticks = (dnum.(a_ticks), Dates.format.(a_ticks, "d u"))

axa = Axis(
    fig[1, 1:3]; ylabel = "Relative CRPS (joint / persistence)",
    yscale = log10, yticks = (
        [0.1, 0.3, 1, 3, 10, 30],
        ["0.1", "0.3", "1", "3", "10", "30"],
    ),
    xticks = a_xticks, xminorticksvisible = true,
    xminorticks = IntervalsBetween(4), xminortickwidth = 0.4, axkw...
)
hlines!(axa, [1.0]; color = :black, linewidth = 0.6, linestyle = :dash)
for (s, col) in STREAMS
    rs = filter(r -> r.stream == s, skill)
    live = filter(r -> !r.backfill, rs)
    lines!(
        axa, [dnum(r.date) for r in live], [r.skill for r in live];
        color = col, linewidth = 0.8
    )
    for r in rs
        scatter!(
            axa, [dnum(r.date)], [r.skill]; markersize = 5,
            color = r.backfill ? :white : col, strokecolor = col,
            strokewidth = 0.8
        )
    end
end
ylims!(axa, 0.025, 40)
a_elems = vcat(
    [
        MarkerElement(; color = col, marker = :circle, markersize = 5)
            for (_, col) in STREAMS
    ],
    [
        MarkerElement(;
            color = :white, strokecolor = :black,
            strokewidth = 0.8, marker = :circle, markersize = 5
        ),
        MarkerElement(; color = :black, marker = :rect, markersize = 4.5),
        MarkerElement(; color = :black, marker = :xcross, markersize = 4.5),
    ]
)
a_labels = vcat(
    first.(STREAMS), [
        "reconstructed backfill",
        "inside 90% interval", "outside 90% interval",
    ]
)
Legend(
    fig[1, 1:3], a_elems, a_labels; tellwidth = false,
    tellheight = false, halign = :center, valign = :bottom,
    framevisible = false, labelsize = FS, rowgap = 0, colgap = 8,
    padding = (2, 2, 2, 2), margin = (2, 0, 2, 0), nbanks = 4
)

# Coverage strip: one row per stream, a filled square where the 90%
# interval of the joint's one-week-ahead forecast covered the observation
# and a cross where it did not (both keyed in panel A's legend).
axs = Axis(
    fig[2, 1:3]; xlabel = "Release data cut-off (2026)",
    xticks = a_xticks, xminorticksvisible = true,
    xminorticks = IntervalsBetween(4), xminortickwidth = 0.4,
    yticks = (1:length(STREAMS), first.(STREAMS)),
    yticklabelsize = FS_SMALL, yreversed = true, axkw...
)
for (i, (s, col)) in enumerate(STREAMS)
    rs = filter(r -> r.stream == s, skill)
    hit = filter(r -> r.covered, rs)
    miss = filter(r -> !r.covered, rs)
    scatter!(
        axs, [dnum(r.date) for r in hit], fill(i, length(hit));
        marker = :rect, markersize = 4.5, color = col
    )
    scatter!(
        axs, [dnum(r.date) for r in miss], fill(i, length(miss));
        marker = :xcross, markersize = 4.5, color = col
    )
end
ylims!(axs, length(STREAMS) + 0.7, 0.3)
for ax in (axa, axs)
    xlims!(ax, a_xlims...)
end
hidexdecorations!(axa; ticks = false, minorticks = false)
linkxaxes!(axa, axs)

# ---------------------------------------------------------------------
# Nested 30/60/90% interval bars with a median dot (helper for B–D)
# ---------------------------------------------------------------------
function interval_bar!(
        ax, x, r; color, horizontal = true, median = true,
        base = 1.6
    )
    spans = [
        (r.lo90, r.hi90, base), (r.lo60, r.hi60, 2 * base),
        (r.lo30, r.hi30, 3 * base),
    ]
    for (lo, hi, w) in spans
        if horizontal
            lines!(ax, [lo, hi], [x, x]; color, linewidth = w)
        else
            lines!(ax, [x, x], [lo, hi]; color, linewidth = w)
        end
    end
    median || return
    pt = horizontal ? ([r.median], [x]) : ([x], [r.median])
    return scatter!(
        ax, pt...; color = :white, strokecolor = color,
        strokewidth = 1.0, markersize = 5
    )
end
quantiles(r) = (
    median = num(r.median), lo30 = num(r.lo30),
    hi30 = num(r.hi30), lo60 = num(r.lo60), hi60 = num(r.hi60),
    lo90 = num(r.lo90), hi90 = num(r.hi90),
)

# ---------------------------------------------------------------------
# B: cumulative infections at the latest release, by fit
# ---------------------------------------------------------------------
ct = Dict(r.fit => quantiles(r) for r in estimates if r.quantity == "C_T")
b_labels = vcat(["joint"], last.(FIT_LABELS))
b_fits = vcat(["joint"], first.(FIT_LABELS))
axb = Axis(
    fig[3, 1]; xscale = log10,
    xlabel = "Cumulative infections to $(Dates.format(CUTOFF, "d u yyyy"))",
    xticks = (
        [1.0e4, 1.0e5, 1.0e6, 1.0e7, 1.0e8],
        ["10⁴", "10⁵", "10⁶", "10⁷", "10⁸"],
    ),
    xminorticksvisible = true, xminorticks = IntervalsBetween(9),
    xminortickwidth = 0.4,
    yticks = (1:length(b_fits), b_labels), yreversed = true, axkw...
)
for (i, f) in enumerate(b_fits)
    interval_bar!(axb, i, ct[f]; color = f == "joint" ? C_JOINT : C_SINGLE)
end
hspan!(axb, 0.5, 1.5; color = (C_JOINT, 0.08))
ylims!(axb, length(b_fits) + 0.6, 0.4)

# ---------------------------------------------------------------------
# C: McCabe et al. scenarios against the joint frozen at matched cut-offs
# ---------------------------------------------------------------------
# Scenario rows carry (vintage date, label, mean, lower, upper). Within a
# vintage, geographic-spread (M1) and back-calculation (M2) scenarios are
# dodged either side of the date, and the frozen joint's nested
# intervals sit on the date itself. Only frozen cut-offs that are also a
# McCabe vintage are drawn (20 and 27 May; the 18 May report has no
# frozen fit of its own, and the 23 May and 8 June fits match no McCabe
# vintage). The frozen asset publishes the interval bounds only, so the
# joint carries no median marker here.
const DC = Date(2026, 5, 16)
cnum(d::Date) = Dates.value(d - DC)
frozen_ct = Dict(
    Date(replace(r.Stream, "frozen " => "") * " 2026", dateformat"d U yyyy") =>
        (
        lo90 = num(r.var"Lower 90%"), hi90 = num(r.var"Upper 90%"),
        lo60 = num(r.var"Lower 60%"), hi60 = num(r.var"Upper 60%"),
        lo30 = num(r.var"Lower 30%"), hi30 = num(r.var"Upper 30%"),
    )
        for r in frozen if startswith(r.Stream, "frozen")
)
c_dates = sort(unique(Date(first(m)) for m in mccabe))
c_ticks = c_dates
axc = Axis(
    fig[3, 2]; ylabel = "Cumulative cases",
    xlabel = "Estimate cut-off (2026)",
    xticks = (cnum.(c_ticks), Dates.format.(c_ticks, "d u")),
    xticklabelrotation = pi / 4, axkw...
)
for (d, q) in frozen_ct
    d in c_dates || continue
    interval_bar!(
        axc, cnum(d), q; color = C_JOINT, horizontal = false,
        median = false, base = 1.8
    )
end
for d in c_dates
    rows = filter(m -> Date(m[1]) == d, mccabe)
    for (method, col, side) in (("M1", C_M1, -1), ("M2", C_M2, 1))
        rs = filter(m -> startswith(m[2], method), rows)
        n = length(rs)
        n == 0 && continue
        offs = side .* (0.35 .+ 0.9 .* (0:(n - 1)) ./ max(n - 1, 1))
        x = cnum(d) .+ offs
        rangebars!(
            axc, x, [m[4] for m in rs], [m[5] for m in rs];
            color = col, linewidth = 0.6
        )
        scatter!(axc, x, [m[3] for m in rs]; color = col, markersize = 3)
    end
end
xlims!(axc, cnum(Date(2026, 5, 16)), cnum(Date(2026, 5, 29)))
ylims!(axc, 0, 5200)
Legend(
    fig[3, 2],
    [
        MarkerElement(; color = C_M1, marker = :circle, markersize = 3),
        MarkerElement(; color = C_M2, marker = :circle, markersize = 3),
        LineElement(; color = C_JOINT, linewidth = 3),
    ],
    [
        "McCabe et al., geographic spread",
        "McCabe et al., back-calculation",
        "joint model frozen at the cut-off",
    ];
    tellwidth = false, tellheight = false, halign = :left, valign = :top,
    framevisible = false, labelsize = FS_SMALL, rowgap = -2,
    padding = (2, 2, 2, 2), margin = (2, 0, 0, 2)
)

# ---------------------------------------------------------------------
# D: Chamla et al. confirmed-case projection against ours and observed
# ---------------------------------------------------------------------
const DD = Date(2026, 5, 16)
ddnum(d::Date) = Dates.value(d - DD)
d_lo, d_hi = Date(2026, 5, 16), Date(2026, 6, 28)
d_ticks = [
    Date(2026, 5, 18), Date(2026, 6, 1), Date(2026, 6, 15),
    Date(2026, 6, 24),
]
axd = Axis(
    fig[3, 3]; ylabel = "Cumulative confirmed cases",
    xlabel = "Date (2026)",
    xticks = (ddnum.(d_ticks), Dates.format.(d_ticks, "d u")),
    xticklabelrotation = pi / 4, axkw...
)
# Observed confirmed cases, from the release's own observation history.
hist = obs["confirmed_case_history"]
h_dates = Date.(hist["dates"])
h_sel = d_lo .<= h_dates .<= d_hi
lines!(
    axd, ddnum.(h_dates[h_sel]), Float64.(hist["values"][h_sel]);
    color = :black, linewidth = 0.8
)
anchor = Date(2026, 6, 8)
vlines!(
    axd, [ddnum(anchor)]; color = (:black, 0.3), linewidth = 0.5,
    linestyle = :dot
)
# Chamla's central scenario over the comparison window, with the low and
# high scenarios at week 12 (24 June).
cw = filter(c -> Date(c[1]) <= d_hi, chamla_central)
cx = ddnum.(Date.(first.(cw))) .- 0.7
rangebars!(
    axd, cx, [c[3] for c in cw], [c[4] for c in cw];
    color = C_CHAMLA, linewidth = 0.8
)
scatter!(axd, cx, [c[2] for c in cw]; color = C_CHAMLA, markersize = 4)
w12 = filter(c -> !occursin("central", c[1]), chamla_w12)
wx = ddnum(Date(2026, 6, 24)) .- 0.7 .+ [-0.9, 0.9]
rangebars!(
    axd, wx, [c[3] for c in w12], [c[4] for c in w12];
    color = (C_CHAMLA, 0.6), linewidth = 0.8
)
scatter!(
    axd, wx, [c[2] for c in w12]; color = (C_CHAMLA, 0.6),
    markersize = 3
)
# Our projection from the fit frozen at 8 June.
ox = ddnum.(getfield.(ours_chamla, :date)) .+ 0.7
rangebars!(
    axd, ox, getfield.(ours_chamla, :lo), getfield.(ours_chamla, :hi);
    color = C_JOINT, linewidth = 0.8
)
scatter!(
    axd, ox, getfield.(ours_chamla, :median); color = :white,
    strokecolor = C_JOINT, strokewidth = 1.0, markersize = 5
)
xlims!(axd, ddnum(d_lo), ddnum(d_hi))
ylims!(axd, 0, 5200)
Legend(
    fig[3, 3],
    [
        MarkerElement(; color = C_CHAMLA, marker = :circle, markersize = 4),
        MarkerElement(;
            color = (C_CHAMLA, 0.6), marker = :circle,
            markersize = 3
        ),
        MarkerElement(;
            color = :white, strokecolor = C_JOINT,
            strokewidth = 1.0, marker = :circle, markersize = 5
        ),
        LineElement(; color = :black, linewidth = 0.8),
    ],
    [
        "Chamla et al., central", "Chamla et al., low, high",
        "joint model frozen at 8 June", "observed confirmed",
    ];
    tellwidth = false, tellheight = false, halign = :left, valign = :top,
    framevisible = false, labelsize = FS_SMALL, rowgap = -2,
    padding = (2, 2, 2, 2), margin = (2, 0, 0, 2)
)

# ---------------------------------------------------------------------
# Layout
# ---------------------------------------------------------------------
panel_letter(fig[1, 1, TopLeft()], "A")
panel_letter(fig[3, 1, TopLeft()], "B")
panel_letter(fig[3, 2, TopLeft()], "C")
panel_letter(fig[3, 3, TopLeft()], "D")
rowsize!(fig.layout, 1, Auto(1.0))
rowsize!(fig.layout, 2, Fixed(13 * MM))
rowsize!(fig.layout, 3, Auto(1.15))
rowgap!(fig.layout, 1, 0)
rowgap!(fig.layout, 2, 8)
colgap!(fig.layout, 10)

save(joinpath(FIGS, "fig-evaluation.pdf"), fig; pt_per_unit = 1)
save(joinpath(FIGS, "fig-evaluation.png"), fig; px_per_unit = 600 / 72)
println("wrote fig-evaluation.pdf and fig-evaluation.png")

# ---------------------------------------------------------------------
# Numbers for the text
# ---------------------------------------------------------------------
let
    beat(rs) = string(count(r -> r.skill < 1, rs), " of ", length(rs))
    live = filter(r -> !r.backfill, skill)
    println(
        "one-week-ahead joint beat persistence: all pairs ", beat(skill),
        "; live releases ", beat(live), "; backfills ",
        beat(filter(r -> r.backfill, skill))
    )
    for (s, _) in STREAMS
        println("  ", s, ": ", beat(filter(r -> r.stream == s, skill)))
    end
    cov = count(r -> r.covered, skill)
    println(
        "90% coverage of the one-week-ahead joint: ", cov, " of ",
        length(skill)
    )
    singles = [ct[f].median for f in first.(FIT_LABELS)]
    hi = first.(FIT_LABELS)[argmax(singles)]
    lo = first.(FIT_LABELS)[argmin(singles)]
    println(
        "single-stream medians at ", CUTOFF, ": largest ", hi, " ",
        round(Int, maximum(singles)), ", smallest ", lo, " ",
        round(Int, minimum(singles)), ", ratio ",
        round(maximum(singles) / minimum(singles); digits = 2),
        "; joint ", round(Int, ct["joint"].median)
    )
    println("Chamla projection rows: ", ours_chamla)
end
