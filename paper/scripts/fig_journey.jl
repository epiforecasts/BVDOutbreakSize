# Figure: the development journey, four panels on one date axis.
#
# Run from the worktree root:
#
#     julia --project=. paper/scripts/fig_journey.jl
#
# Reads paper/data/{code_size,commits_weekly,release_estimates,
# data_events,fit_cost}.csv and writes paper/figures/fig-journey.pdf
# and fig-journey.png (180 mm wide, 600 dpi). Colours are Makie's Wong
# palette, as in src/plots.jl.

using CairoMakie
using Dates
using Statistics: mean

const REPO = normpath(joinpath(@__DIR__, "..", ".."))
const DATA = joinpath(REPO, "paper", "data")
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
    out
end
function read_csv(name)
    lines = filter(!isempty, split(read(joinpath(DATA, name), String), '\n'))
    header = Symbol.(split_csv_line(lines[1]))
    [NamedTuple{Tuple(header)}(Tuple(split_csv_line(l))) for l in lines[2:end]]
end
num(s) = isempty(s) ? nothing : parse(Float64, s)
int(s) = parse(Int, s)

code = read_csv("code_size.csv")
weekly = read_csv("commits_weekly.csv")
rel = read_csv("release_estimates.csv")
events = read_csv("data_events.csv")
cost = read_csv("fit_cost.csv")

# ---------------------------------------------------------------------
# Date axis: days since 18 May 2026
# ---------------------------------------------------------------------
const D0 = Date(2026, 5, 18)
dnum(d::Date) = Dates.value(d - D0)
dnum(s::AbstractString) = dnum(Date(s))
last_date = maximum(vcat(
    [Date(r.date) for r in code],
    [Date(r.week_start) + Day(7) for r in weekly],
    [Date(r.cutoff) for r in rel if !isempty(r.cutoff)],
    [Date(r.date) for r in events],
    [Date(r.date) for r in cost]))
xlims = (dnum(D0) - 1, dnum(last_date) + 1)
month_starts = [Date(2026, m, 1) for m in 6:month(last_date)]
xticks = (dnum.(month_starts), Dates.format.(month_starts, "d u"))

# ---------------------------------------------------------------------
# Colours and sizes
# ---------------------------------------------------------------------
wong = Makie.wong_colors()
const C_SRC = wong[1]      # blue
const C_TEST = wong[2]     # orange
const C_HUMAN = wong[6]    # vermillion
const C_AGENT = wong[3]    # green
const C_OTHER = :grey60
const C_TAG = (:black, 0.25)
const BANDS = [
    ("closed-form", "v1.0.0", "v1.4.0", (wong[4], 0.12)),
    ("renewal", "v1.4.0", "V2.0.0", (wong[5], 0.12)),
    ("provincial", "V2.0.0", nothing, (wong[7], 0.12)),
]
const MM = 72 / 25.4
const FS = 7.0            # base font size (pt)
const FS_SMALL = 5.0

tagdate = Dict(r.tag => Date(r.date) for r in code)
tagx = [dnum(Date(r.date)) for r in code]

function band_edges(b)
    x0 = b[2] == "v1.0.0" ? xlims[1] : dnum(tagdate[b[2]])
    x1 = b[3] === nothing ? xlims[2] : dnum(tagdate[b[3]])
    x0, x1
end
function decorate!(ax)
    for b in BANDS
        x0, x1 = band_edges(b)
        vspan!(ax, x0, x1; color = b[4])
    end
    vlines!(ax, tagx; color = C_TAG, linewidth = 0.4)
end
# Vertical stagger for labels closer than `gap` days.
function stagger(xs, gap)
    lvl = zeros(Int, length(xs))
    for i in 2:length(xs)
        lvl[i] = xs[i] - xs[i - 1] < gap ? 1 - lvl[i - 1] : 0
    end
    lvl
end

fig = Figure(; size = (180 * MM, 235 * MM), fontsize = FS,
    figure_padding = (4, 8, 4, 4))
axkw = (; xticks, xgridvisible = false, ygridvisible = false,
    xticklabelsize = FS, yticklabelsize = FS, xlabelsize = FS,
    ylabelsize = FS, xminorticksvisible = true,
    xminorticks = IntervalsBetween(4), spinewidth = 0.6,
    xtickwidth = 0.6, ytickwidth = 0.6, xminortickwidth = 0.4)
# Legend drawn inside a panel, top-left or top-right.
function inset_legend!(pos, elements, labels; halign = :left)
    Legend(pos, elements, labels; tellwidth = false, tellheight = false,
        halign, valign = :top, framevisible = false, labelsize = FS,
        patchsize = (8, 6), rowgap = 0, padding = (2, 2, 2, 2),
        margin = (2, 2, 0, 2))
end
panel_letter(row, s) = Label(fig[row, 1, TopLeft()], s; font = :bold,
    fontsize = 9, padding = (0, 6, 0, 0), halign = :right, valign = :bottom)

# ---------------------------------------------------------------------
# Row 1: release-tag labels and model-version band names
# ---------------------------------------------------------------------
axl = Axis(fig[1, 1]; axkw...)
hidedecorations!(axl); hidespines!(axl)
for b in BANDS
    x0, x1 = band_edges(b)
    vspan!(axl, x0, x1; color = b[4])
    text!(axl, x0 + 0.6, 1.0; text = b[1], fontsize = FS,
        align = (:left, :top), font = :italic)
end
vlines!(axl, tagx; color = C_TAG, linewidth = 0.4)
lvl = stagger(tagx, 3)
for (r, x, l) in zip(code, tagx, lvl)
    text!(axl, x, 0.02 + 0.40 * l; text = r.tag, fontsize = FS_SMALL,
        rotation = pi / 2, align = (:left, :center))
end
ylims!(axl, 0, 1)

# ---------------------------------------------------------------------
# A: lines of code and fitted streams
# ---------------------------------------------------------------------
axa = Axis(fig[2, 1]; ylabel = "Lines of code (thousands)", axkw...)
decorate!(axa)
stairs!(axa, tagx, [int(r.loc_src) / 1000 for r in code]; step = :post,
    color = C_SRC, linewidth = 1.2, label = "src/")
stairs!(axa, tagx, [int(r.loc_test) / 1000 for r in code]; step = :post,
    color = C_TEST, linewidth = 1.2, label = "test/")
ylims!(axa, 0, nothing)
axislegend(axa; position = :lt, framevisible = false, labelsize = FS,
    patchsize = (10, 4), rowgap = 0, padding = (2, 2, 2, 2))
axa2 = Axis(fig[2, 1]; yaxisposition = :right, ylabel = "Fitted streams",
    yticks = 0:2:12, ygridvisible = false, xgridvisible = false,
    yticklabelsize = FS, ylabelsize = FS, spinewidth = 0.6,
    ytickwidth = 0.6, rightspinecolor = :grey30)
hidexdecorations!(axa2)
hidespines!(axa2, :l, :t, :b)
stairs!(axa2, tagx, [int(r.n_streams_fitted) for r in code]; step = :post,
    color = :grey30, linewidth = 1.0, linestyle = :dash)
text!(axa2, tagx[end] - 1, int(code[end].n_streams_fitted) - 0.4;
    text = "streams", fontsize = FS_SMALL, color = :grey30,
    align = (:right, :top))
ylims!(axa2, 0, 13)
linkxaxes!(axa, axa2)

# ---------------------------------------------------------------------
# B: commits per week, stacked, and merged PRs
# ---------------------------------------------------------------------
axb = Axis(fig[3, 1]; ylabel = "Commits and merged PRs per week", axkw...)
decorate!(axb)
wx = [dnum(Date(r.week_start)) + 3.5 for r in weekly]
groups = [(:n_commits_human, C_HUMAN, "human"),
    (:n_commits_agent, C_AGENT, "agent"),
    (:n_commits_other, C_OTHER, "other")]
barplot!(axb,
    repeat(wx, length(groups)),
    vcat([[int(getfield(r, g[1])) for r in weekly] for g in groups]...);
    stack = repeat(1:length(groups); inner = length(wx)),
    color = repeat([g[2] for g in groups]; inner = length(wx)),
    width = 6.4, gap = 0, strokewidth = 0)
prs = [int(r.n_prs_merged_human) + int(r.n_prs_merged_agent) +
       int(r.n_prs_merged_other) for r in weekly]
scatterlines!(axb, wx, prs; color = :black, linewidth = 0.8,
    markersize = 3.5, label = "merged PRs")
ylims!(axb, 0, nothing)
inset_legend!(fig[3, 1],
    [[PolyElement(; color = g[2]) for g in groups]...,
     LineElement(; color = :black, linewidth = 0.8)],
    [[g[3] * " commits" for g in groups]..., "merged PRs"];
    halign = :center)

# ---------------------------------------------------------------------
# C: released headline estimate against its data cut-off, data events
# ---------------------------------------------------------------------
axc = Axis(fig[4, 1]; ylabel = "Cumulative infections (log scale)",
    yscale = log10, yticks = [500, 1000, 2000, 5000, 10000, 20000],
    ytickformat = v -> string.(round.(Int, v)), axkw...)
decorate!(axc)
have = filter(r -> !isempty(r.median), rel)
rx = Float64[dnum(r.cutoff) for r in have]
# Releases sharing a cut-off (v1.0.0/v1.1.0, v1.3.0/v1.4.0) are drawn
# 0.7 days either side of it so both intervals show.
for i in 2:length(rx)
    if rx[i] == rx[i - 1]
        rx[i - 1] -= 0.7
        rx[i] += 0.7
    end
end
rangebars!(axc, rx, [num(r.lower90) for r in have],
    [num(r.upper90) for r in have]; color = :black, linewidth = 0.8)
scatter!(axc, rx, [num(r.median) for r in have]; color = :black,
    markersize = 5, strokecolor = :white, strokewidth = 0.5)
ylims!(axc, 300, 30000)
inset_legend!(fig[4, 1],
    [MarkerElement(; color = :black, marker = :circle, markersize = 5),
     LineElement(; color = :black, linewidth = 0.8),
     MarkerElement(; color = C_HUMAN, marker = :dtriangle, markersize = 5),
     MarkerElement(; color = C_AGENT, marker = :dtriangle, markersize = 5)],
    ["release median", "90% interval", "data event, human decision",
     "data event, no human decision"])

# Data events on a strip under panel C, coloured by human decision.
const EVENT_LABELS = Dict(
    "2026-05-27" => "source moves to INSP",
    "2026-05-28" => "suspects reclassified",
    "2026-06-02" => "GeneXpert misses BDBV",
    "2026-06-06" => "Uganda exports lag",
    "2026-06-11" => "partial 24h analysed",
    "2026-06-19" => "DHIS2 isolation step",
    "2026-07-05" => "Tableau 6 absent",
    "2026-07-12" => "onset curve appears",
    "2026-07-22" => "harmonisation step",
    "2026-08-06" => "short brief format",
    "2026-08-12" => "samples total dropped",
    "2026-08-13" => "isolation census drop",
    "2026-08-20" => "onset figure 7% high",
    "2026-09-07" => "onset past report date",
    "2026-09-16" => "suspects resumed",
)
axe = Axis(fig[5, 1]; axkw...)
hidedecorations!(axe); hidespines!(axe, :l, :r, :t)
decorate!(axe)
ex = [dnum(r.date) for r in events]
ecol = [r.human_decision == "yes" ? C_HUMAN : C_AGENT for r in events]
elvl = stagger(ex, 2)
scatter!(axe, ex, fill(0.06, length(ex)); color = ecol,
    marker = :dtriangle, markersize = 5)
for (r, x, c, l) in zip(events, ex, ecol, elvl)
    text!(axe, x, 0.14; text = EVENT_LABELS[r.date], fontsize = FS_SMALL,
        rotation = pi / 2, align = (:left, :center), color = c,
        offset = (3 * l, 0))
end
ylims!(axe, 0, 1)

# ---------------------------------------------------------------------
# D: joint gradient time per PR, and the one wall-clock figure
# ---------------------------------------------------------------------
axd = Axis(fig[6, 1]; ylabel = "Joint gradient (ms, log scale)",
    xlabel = "Date (2026)", yscale = log10, yticks = [1, 2, 5, 10, 20],
    ytickformat = v -> string.(round.(Int, v)), axkw...)
decorate!(axd)
is_prod(r) = startswith(r.model, "production joint")
is_ci(r) = startswith(r.model, "CI benchmark")
series = [("production joint", is_prod, wong[1], :circle),
    ("CI benchmark joint", is_ci, wong[2], :rect)]
for (name, pred, col, mk) in series
    rows = filter(r -> pred(r) && r.joint_gradient_ms != "", cost)
    for (k, pr) in enumerate(unique(r.tag_or_pr for r in rows))
        rr = filter(r -> r.tag_or_pr == pr, rows)
        x = dnum(rr[1].date)
        ys = [num(r.joint_gradient_ms) for r in rr]
        length(ys) == 2 && lines!(axd, [x, x], ys; color = col,
            linewidth = 0.8)
        for r in rr
            y = num(r.joint_gradient_ms)
            scatter!(axd, [x], [y]; marker = mk, markersize = 5,
                color = r.arm == "after" ? col : :white,
                strokecolor = col, strokewidth = 0.8)
        end
        # PRs a day apart: label above, then right, then below.
        pos, al, off = (k % 3 == 1) ? ((x, maximum(ys)), (:center, :bottom), (0, 2)) :
            (k % 3 == 2) ? ((x + 1.0, exp(mean(log.(ys)))), (:left, :center), (0, 0)) :
            ((x, minimum(ys)), (:center, :top), (0, -2))
        text!(axd, pos...; text = pr, fontsize = FS_SMALL, align = al,
            offset = off)
    end
end
wall = only(filter(r -> r.joint_fit_minutes != "", cost))
wy = 4.0   # the wall-clock has no ms value; its height is arbitrary
scatter!(axd, [dnum(wall.date)], [wy]; marker = :diamond, markersize = 6,
    color = :black)
text!(axd, dnum(wall.date) - 1.5, wy;
    text = "$(wall.tag_or_pr): $(wall.joint_fit_minutes) min fit wall-clock",
    fontsize = FS_SMALL, align = (:right, :center))
ylims!(axd, 0.8, 28)
inset_legend!(fig[6, 1],
    [MarkerElement(; color = wong[1], marker = :circle, markersize = 5),
     MarkerElement(; color = wong[2], marker = :rect, markersize = 5),
     MarkerElement(; color = :white, strokecolor = :black,
         strokewidth = 0.8, marker = :circle, markersize = 5),
     MarkerElement(; color = :black, marker = :circle, markersize = 5)],
    ["production joint", "CI benchmark joint", "before PR", "after PR"])

# ---------------------------------------------------------------------
# Layout
# ---------------------------------------------------------------------
for (row, s) in zip((2, 3, 4, 6), ["A", "B", "C", "D"])
    panel_letter(row, s)
end
for ax in (axl, axa, axb, axc, axe, axd)
    xlims!(ax, xlims...)
end
for ax in (axa, axb, axc)
    hidexdecorations!(ax; ticks = false, minorticks = false)
end
linkxaxes!(axl, axa, axb, axc, axe, axd)
rowsize!(fig.layout, 1, Fixed(16 * MM))
rowsize!(fig.layout, 2, Auto(1.0))
rowsize!(fig.layout, 3, Auto(1.0))
rowsize!(fig.layout, 4, Auto(1.3))
rowsize!(fig.layout, 5, Fixed(22 * MM))
rowsize!(fig.layout, 6, Auto(1.0))
rowgap!(fig.layout, 3)
rowgap!(fig.layout, 4, 0)

save(joinpath(FIGS, "fig-journey.pdf"), fig; pt_per_unit = 1)
save(joinpath(FIGS, "fig-journey.png"), fig; px_per_unit = 600 / 72)
println("wrote fig-journey.pdf and fig-journey.png")
