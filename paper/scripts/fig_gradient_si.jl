# Supplementary figure: the recorded joint gradient time per pull request
# on the development date axis, with the one fit wall-clock figure.
#
# Run from the worktree root:
#
#     julia --project=. paper/scripts/fig_gradient_si.jl
#
# Reads paper/data/{code_size,fit_cost}.csv and writes
# paper/figures/fig-gradient-si.pdf and fig-gradient-si.png (180 mm wide,
# 600 dpi), in the style of fig_journey.jl: the same date axis, release-tag
# lines, model-version bands and Wong palette.

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
    return out
end
function read_csv(name)
    lines = filter(!isempty, split(read(joinpath(DATA, name), String), '\n'))
    header = Symbol.(split_csv_line(lines[1]))
    return [NamedTuple{Tuple(header)}(Tuple(split_csv_line(l))) for l in lines[2:end]]
end
num(s) = isempty(s) ? nothing : parse(Float64, s)

code = read_csv("code_size.csv")
cost = read_csv("fit_cost.csv")

# ---------------------------------------------------------------------
# Date axis: days since 18 May 2026, to the last tag or cost row
# ---------------------------------------------------------------------
const D0 = Date(2026, 5, 18)
dnum(d::Date) = Dates.value(d - D0)
dnum(s::AbstractString) = dnum(Date(s))
last_date = maximum(
    vcat([Date(r.date) for r in code], [Date(r.date) for r in cost])
)
xlims = (dnum(D0) - 1, dnum(last_date) + 4)
month_starts = [Date(2026, m, 1) for m in 6:month(last_date)]
xticks = (dnum.(month_starts), Dates.format.(month_starts, "d u"))

# ---------------------------------------------------------------------
# Colours and sizes
# ---------------------------------------------------------------------
wong = Makie.wong_colors()
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
    return x0, x1
end
function decorate!(ax)
    for b in BANDS
        x0, x1 = band_edges(b)
        vspan!(ax, x0, x1; color = b[4])
    end
    return vlines!(ax, tagx; color = C_TAG, linewidth = 0.4)
end
# Vertical stagger for labels closer than `gap` days.
function stagger(xs, gap)
    lvl = zeros(Int, length(xs))
    for i in 2:length(xs)
        lvl[i] = xs[i] - xs[i - 1] < gap ? 1 - lvl[i - 1] : 0
    end
    return lvl
end

fig = Figure(;
    size = (180 * MM, 75 * MM), fontsize = FS,
    figure_padding = (4, 8, 4, 4)
)
axkw = (;
    xticks, xgridvisible = false, ygridvisible = false,
    xticklabelsize = FS, yticklabelsize = FS, xlabelsize = FS,
    ylabelsize = FS, xminorticksvisible = true,
    xminorticks = IntervalsBetween(4), spinewidth = 0.6,
    xtickwidth = 0.6, ytickwidth = 0.6, xminortickwidth = 0.4,
)

# ---------------------------------------------------------------------
# Row 1: release-tag labels and model-version band names
# ---------------------------------------------------------------------
axl = Axis(fig[1, 1]; axkw...)
hidedecorations!(axl); hidespines!(axl)
for b in BANDS
    x0, x1 = band_edges(b)
    vspan!(axl, x0, x1; color = b[4])
    text!(
        axl, x0 + 0.6, 1.0; text = b[1], fontsize = FS,
        align = (:left, :top), font = :italic
    )
end
vlines!(axl, tagx; color = C_TAG, linewidth = 0.4)
lvl = stagger(tagx, 3)
for (r, x, l) in zip(code, tagx, lvl)
    text!(
        axl, x, 0.02 + 0.4 * l; text = r.tag, fontsize = FS_SMALL,
        rotation = pi / 2, align = (:left, :center)
    )
end
ylims!(axl, 0, 1)

# ---------------------------------------------------------------------
# Joint gradient time per PR, and the one wall-clock figure
# ---------------------------------------------------------------------
axd = Axis(
    fig[2, 1]; ylabel = "Joint gradient (ms, log scale)",
    xlabel = "Date (2026)", yscale = log10, yticks = [1, 2, 5, 10, 20],
    ytickformat = v -> string.(round.(Int, v)), axkw...
)
decorate!(axd)
is_prod(r) = startswith(r.model, "production joint")
is_ci(r) = startswith(r.model, "CI benchmark")
series = [
    ("production joint", is_prod, wong[1], :circle),
    ("CI benchmark joint", is_ci, wong[2], :rect),
]
for (name, pred, col, mk) in series
    rows = filter(r -> pred(r) && r.joint_gradient_ms != "", cost)
    for (k, pr) in enumerate(unique(r.tag_or_pr for r in rows))
        rr = filter(r -> r.tag_or_pr == pr, rows)
        x = dnum(rr[1].date)
        ys = [num(r.joint_gradient_ms) for r in rr]
        length(ys) == 2 && lines!(
            axd, [x, x], ys; color = col,
            linewidth = 0.8
        )
        for r in rr
            y = num(r.joint_gradient_ms)
            scatter!(
                axd, [x], [y]; marker = mk, markersize = 5,
                color = r.arm == "after" ? col : :white,
                strokecolor = col, strokewidth = 0.8
            )
        end
        # PRs a day apart: label above, then right, then below.
        pos, al, off = (k % 3 == 1) ? ((x, maximum(ys)), (:center, :bottom), (0, 2)) :
            (k % 3 == 2) ? ((x + 1.0, exp(mean(log.(ys)))), (:left, :center), (0, 0)) :
            ((x, minimum(ys)), (:center, :top), (0, -2))
        text!(
            axd, pos...; text = pr, fontsize = FS_SMALL, align = al,
            offset = off
        )
    end
end
wall = only(filter(r -> r.joint_fit_minutes != "", cost))
wy = 4.0   # the wall-clock has no ms value; its height is arbitrary
scatter!(
    axd, [dnum(wall.date)], [wy]; marker = :diamond, markersize = 6,
    color = :black
)
text!(
    axd, dnum(wall.date) - 1.5, wy;
    text = "$(wall.tag_or_pr): $(wall.joint_fit_minutes) min fit wall-clock",
    fontsize = FS_SMALL, align = (:right, :center)
)
ylims!(axd, 0.8, 28)
Legend(
    fig[2, 1],
    [
        MarkerElement(; color = wong[1], marker = :circle, markersize = 5),
        MarkerElement(; color = wong[2], marker = :rect, markersize = 5),
        MarkerElement(;
            color = :white, strokecolor = :black,
            strokewidth = 0.8, marker = :circle, markersize = 5
        ),
        MarkerElement(; color = :black, marker = :circle, markersize = 5),
    ],
    ["production joint", "CI benchmark joint", "before PR", "after PR"];
    tellwidth = false, tellheight = false, halign = :left, valign = :top,
    framevisible = false, labelsize = FS, patchsize = (8, 6), rowgap = 0,
    padding = (2, 2, 2, 2), margin = (2, 2, 0, 2)
)

# ---------------------------------------------------------------------
# Layout
# ---------------------------------------------------------------------
for ax in (axl, axd)
    xlims!(ax, xlims...)
end
linkxaxes!(axl, axd)
rowsize!(fig.layout, 1, Fixed(16 * MM))
rowsize!(fig.layout, 2, Auto(1.0))
rowgap!(fig.layout, 3)

save(joinpath(FIGS, "fig-gradient-si.pdf"), fig; pt_per_unit = 1)
save(joinpath(FIGS, "fig-gradient-si.png"), fig; px_per_unit = 600 / 72)
println("wrote fig-gradient-si.pdf and fig-gradient-si.png")
