# Figure: where the work started, the first model's release (results-v1.0.0).
#
# Run from the worktree root:
#
#     julia --project=. paper/scripts/fig_first.jl
#
# Reads the assets of the results-v1.0.0 GitHub release from
# paper/data/release/results-v1.0.0/ (git-ignored; downloaded with
# `gh release download results-v1.0.0 -R epiforecasts/BVDOutbreakSize`
# on the first run) and writes paper/figures/fig-first.pdf and
# fig-first.png (180 mm wide, 600 dpi). Colours are Makie's Wong
# palette, sizes and axis style as in fig_journey.jl.
#
# Assets used:
#   imperial_comparison.csv       McCabe et al. headline scenarios with
#                                 their reported 90% intervals, the
#                                 Method 2 reproduction, and the joint
#                                 fit to the report's data and to the
#                                 18 May data (median, 90% interval)
#   scenario_coverage.csv         the fifteen McCabe et al. scenario
#                                 point estimates
#   cumulative_cases_by_stream.csv 30/60/90% intervals of C(T) from the
#                                 four single-stream fits and the joint
#   posterior_draws.csv           400 thinned joint draws of τ, r, m,
#                                 T, CFR, p_drc, p_uganda,
#                                 cumulative_cases
#   observations.toml             the four counts fitted and the
#                                 data cut-off
#
# Nothing is estimated beyond quantiles of the released draws. The
# latent curve in panel C is the model's own C(s) = exp(r s) from
# seeding at T days before the cut-off, evaluated per draw.

using CairoMakie
using Dates
using Statistics

const REPO = normpath(joinpath(@__DIR__, "..", ".."))
const TAG = "results-v1.0.0"
const REL = joinpath(REPO, "paper", "data", "release", TAG)
const FIGS = joinpath(REPO, "paper", "figures")
mkpath(FIGS)

const ASSETS = (
    "imperial_comparison.csv", "scenario_coverage.csv",
    "cumulative_cases_by_stream.csv", "posterior_draws.csv",
    "observations.toml",
)
for name in ASSETS
    isfile(joinpath(REL, name)) && continue
    mkpath(REL)
    run(`gh release download $TAG -R epiforecasts/BVDOutbreakSize
        -p $name -D $REL`)
end

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
    lines = filter(!isempty, split(read(joinpath(REL, name), String), '\n'))
    header = Symbol.(replace.(split_csv_line(lines[1]), r"[ %]" => "_"))
    return [NamedTuple{Tuple(header)}(Tuple(split_csv_line(l))) for l in lines[2:end]]
end
num(s) = parse(Float64, s)

comparison = read_csv("imperial_comparison.csv")
scenarios = read_csv("scenario_coverage.csv")
by_stream = read_csv("cumulative_cases_by_stream.csv")
draws = read_csv("posterior_draws.csv")

# observations.toml: `as_of_date` and the `value` of each [section]
function read_observations(path)
    vals = Dict{String, Any}()
    section = ""
    for l in eachline(path)
        s = strip(l)
        if startswith(s, "[")
            section = strip(s, ['[', ']'])
        elseif startswith(s, "as_of_date")
            vals["as_of_date"] = Date(strip(split(s, "=")[2], [' ', '"']))
        elseif startswith(s, "value") && !isempty(section)
            vals[section] = parse(Float64, strip(split(s, "=")[2]))
        end
    end
    return vals
end
obs = read_observations(joinpath(REL, "observations.toml"))
const CUTOFF = obs["as_of_date"]

# ---------------------------------------------------------------------
# Colours and sizes (as fig_journey.jl)
# ---------------------------------------------------------------------
wong = Makie.wong_colors()
const C_JOINT = wong[1]     # blue: this model, joint fit
const C_PUB = :grey35       # published scenario estimates
const C_STREAM = :grey55    # single-stream fits
const C_OBS = wong[6]       # vermillion: observed counts
const MM = 72 / 25.4
const FS = 7.0
const FS_SMALL = 5.5

fig = Figure(;
    size = (180 * MM, 150 * MM), fontsize = FS,
    figure_padding = (4, 8, 4, 4)
)
axkw = (;
    xgridvisible = false, ygridvisible = false,
    xticklabelsize = FS, yticklabelsize = FS, xlabelsize = FS,
    ylabelsize = FS, spinewidth = 0.6, xtickwidth = 0.6,
    ytickwidth = 0.6,
)
panel_letter(pos, s) = Label(
    pos, s; font = :bold, fontsize = 9,
    padding = (0, 6, 0, 0), halign = :right, valign = :bottom
)

# Nested interval bar: outer (thin) and inner (thick) ranges, optional
# point.
function interval!(ax, y, lo_out, hi_out, lo_in, hi_in, mid, colour)
    lines!(ax, [lo_out, hi_out], [y, y]; color = colour, linewidth = 0.8)
    lo_in === nothing ||
        lines!(ax, [lo_in, hi_in], [y, y]; color = colour, linewidth = 2.6)
    mid === nothing && return
    return scatter!(
        ax, [mid], [y]; color = colour, markersize = 6,
        strokecolor = :white, strokewidth = 0.6
    )
end

# ---------------------------------------------------------------------
# Panel A: McCabe et al. scenarios against the first joint fit
# ---------------------------------------------------------------------
# Rows from imperial_comparison.csv, keyed by their Source label.
row(src) = only(filter(r -> r.Source == src, comparison))
m1 = row("McCabe Method 1 (Ituri, w=15 d)")
m2 = row("McCabe Method 2 (τ=14 d, CFR 30%)")
repro = row("Our Method 2 reproduction")
joint_report = row("Our joint (report data, 16 May)")
joint_now = row("Our joint (current data)")
# The fifteen scenario point estimates, split by method.
m1_points = [
    num(r.Reported_cases) for r in scenarios
        if startswith(r.Scenario, "Method 1")
]
m2_points = [
    num(r.Reported_cases) for r in scenarios
        if startswith(r.Scenario, "Method 2")
]
# 50% interval of the headline joint from the thinned draws.
C_draws = [num(r.cumulative_cases) for r in draws]
q(p) = quantile(C_draws, p)

a_labels = [
    "McCabe et al. Method 1\n(exports, 6 scenarios)",
    "McCabe et al. Method 2\n(deaths, 9 scenarios)",
    "Method 2 reproduction\n(fixed inputs)",
    "Joint fit, report data\n(16 May)",
    "Joint fit, $(Dates.format(CUTOFF, "d U")) data",
]
ya = collect(length(a_labels):-1:1)
axa = Axis(
    fig[1, 1]; axkw...,
    xlabel = "Cumulative infections at cut-off",
    yticks = (ya, a_labels), yticklabelsize = FS_SMALL,
    xticks = 0:500:2500
)
xlims!(axa, -30, 2600)
ylims!(axa, 0.4, length(a_labels) + 0.6)
# McCabe et al. rows: headline interval as a thin bar, every scenario
# point as an open marker, the headline scenario filled.
for (y, hl, pts) in ((ya[1], m1, m1_points), (ya[2], m2, m2_points))
    lines!(
        axa, [num(hl.Lower_90_), num(hl.Upper_90_)], [y, y];
        color = C_PUB, linewidth = 0.8
    )
    scatter!(
        axa, pts, fill(y, length(pts)); color = :white,
        strokecolor = C_PUB, strokewidth = 0.8, markersize = 5
    )
    scatter!(
        axa, [num(hl.Central_estimate)], [y]; color = C_PUB,
        markersize = 6, strokecolor = :white, strokewidth = 0.6
    )
end
for (y, r, colour) in ((ya[3], repro, C_PUB), (ya[4], joint_report, C_JOINT))
    interval!(
        axa, y, num(r.Lower_90_), num(r.Upper_90_),
        nothing, nothing, num(r.Central_estimate), colour
    )
end
interval!(
    axa, ya[5], num(joint_now.Lower_90_), num(joint_now.Upper_90_),
    q(0.25), q(0.75), num(joint_now.Central_estimate), C_JOINT
)
Legend(
    fig[1, 1],
    [
        MarkerElement(; color = C_PUB, marker = :circle, markersize = 6),
        MarkerElement(;
            color = :white, strokecolor = C_PUB,
            strokewidth = 0.8, marker = :circle, markersize = 5
        ),
        LineElement(; color = C_PUB, linewidth = 0.8),
        LineElement(; color = C_JOINT, linewidth = 2.6),
    ],
    [
        "headline scenario", "other scenarios", "90% interval",
        "50% interval",
    ];
    tellwidth = false, tellheight = false, halign = :right,
    valign = :top, framevisible = false, labelsize = FS_SMALL,
    patchsize = (8, 6), rowgap = -1, padding = (2, 2, 2, 2),
    margin = (2, 4, 2, 4)
)

# ---------------------------------------------------------------------
# Panel B: what each stream implied on its own
# ---------------------------------------------------------------------
stream_rows = [
    ("exports (cases)", "Exported cases (Uganda)", "exported_cases"),
    ("exports (deaths)", "Deaths among exports", "exports_deaths"),
    ("deaths (DRC)", "Suspected deaths (DRC)", "total_deaths"),
    ("cases (DRC)", "Suspected cases (DRC)", "reported_cases"),
    ("joint", "All four streams (joint)", nothing),
]
# Row label with the count that stream's fit conditioned on.
function stream_label((key, label, obskey))
    return obskey === nothing ? label :
        "$label\n(fitted count $(round(Int, obs[obskey])))"
end
yb = collect(length(stream_rows):-1:1)
axb = Axis(
    fig[1, 2]; axkw...,
    xlabel = "Cumulative infections at cut-off",
    yticks = (yb, stream_label.(stream_rows)), yticklabelsize = FS_SMALL,
    xticks = 0:1000:5000
)
xlims!(axb, -60, 5200)
ylims!(axb, 0.4, length(stream_rows) + 0.6)
for (y, (key, _, obskey)) in zip(yb, stream_rows)
    r = only(filter(r -> r.Stream == key, by_stream))
    colour = obskey === nothing ? C_JOINT : C_STREAM
    lines!(
        axb, [num(r.Lower_90_), num(r.Upper_90_)], [y, y];
        color = colour, linewidth = 0.8
    )
    lines!(
        axb, [num(r.Lower_60_), num(r.Upper_60_)], [y, y];
        color = colour, linewidth = 1.8
    )
    lines!(
        axb, [num(r.Lower_30_), num(r.Upper_30_)], [y, y];
        color = colour, linewidth = 3.2
    )
end
Legend(
    fig[1, 2],
    [
        LineElement(; color = C_STREAM, linewidth = 0.8),
        LineElement(; color = C_STREAM, linewidth = 1.8),
        LineElement(; color = C_STREAM, linewidth = 3.2),
    ],
    ["90% interval", "60% interval", "30% interval"];
    tellwidth = false, tellheight = false, halign = :right,
    valign = :bottom, framevisible = false, labelsize = FS_SMALL,
    patchsize = (8, 6), rowgap = -1, padding = (2, 2, 2, 2),
    margin = (2, 4, 4, 2)
)

# ---------------------------------------------------------------------
# Panel C: the latent constant-growth curve and the seeding date
# ---------------------------------------------------------------------
r_draws = [num(r.r) for r in draws]
T_draws = [num(r.T) for r in draws]
dnum(d::Date) = Dates.value(d - CUTOFF)      # days relative to the cut-off
seed_days = -T_draws                         # seeding day per draw
t0 = floor(Int, quantile(seed_days, 0.01))
tgrid = t0:0
curve(i, t) = t < seed_days[i] ? 0.0 : exp(r_draws[i] * (t - seed_days[i]))
Cmat = [curve(i, t) for t in tgrid, i in eachindex(r_draws)]
qs = [
    quantile(view(Cmat, k, :), p) for k in axes(Cmat, 1),
        p in (0.05, 0.25, 0.5, 0.75, 0.95)
]
first_of_month = [Date(y, m, 1) for y in 2025:2026, m in 1:12]
months = filter(d -> t0 <= dnum(d) <= 0, vec(first_of_month))
cticks = (dnum.(months), Dates.format.(months, "u yyyy"))

axc = Axis(
    fig[2, 1:2]; axkw..., ylabel = "Cumulative infections",
    xticks = cticks, xticklabelsvisible = false, xticksvisible = false,
    xminorticksvisible = false
)
tx = collect(Float64, tgrid)
band!(axc, tx, qs[:, 1], qs[:, 5]; color = (C_JOINT, 0.15))
band!(axc, tx, qs[:, 2], qs[:, 4]; color = (C_JOINT, 0.35))
lines!(axc, tx, qs[:, 3]; color = C_JOINT, linewidth = 1.2)
vlines!(axc, [0.0]; color = :grey30, linewidth = 0.6, linestyle = :dash)
# The two DRC counts the model fitted, at the cut-off.
for (key, label, dy) in (
        ("reported_cases", "suspected cases", 0),
        ("total_deaths", "suspected deaths", 0),
    )
    scatter!(
        axc, [0.0], [obs[key]]; color = C_OBS, markersize = 6,
        strokecolor = :white, strokewidth = 0.6
    )
    text!(
        axc, -4, obs[key] + dy; text = label, fontsize = FS_SMALL,
        color = C_OBS, align = (:right, :center)
    )
end
text!(
    axc, -4, 2550; text = "data cut-off\n$(Dates.format(CUTOFF, "d U yyyy"))",
    fontsize = FS_SMALL, color = :grey30, align = (:right, :top)
)
xlims!(axc, t0, 8)
ylims!(axc, -40, 2650)
Legend(
    fig[2, 1:2],
    [
        LineElement(; color = C_JOINT, linewidth = 1.2),
        PolyElement(; color = (C_JOINT, 0.35)),
        PolyElement(; color = (C_JOINT, 0.15)),
        MarkerElement(; color = C_OBS, marker = :circle, markersize = 6),
    ],
    [
        "posterior median", "50% interval", "90% interval",
        "fitted DRC count",
    ];
    tellwidth = false, tellheight = false, halign = :left,
    valign = :top, framevisible = false, labelsize = FS_SMALL,
    patchsize = (8, 6), rowgap = -1, padding = (2, 2, 2, 2),
    margin = (4, 2, 2, 4)
)

# Seeding-date marginal on the same date axis.
axs = Axis(
    fig[3, 1:2]; axkw..., xlabel = "Date",
    ylabel = "Seeding\ndensity", xticks = cticks, yticksvisible = false,
    yticklabelsvisible = false, leftspinevisible = false,
    rightspinevisible = false, topspinevisible = false
)
density!(
    axs, seed_days; color = (C_JOINT, 0.35), strokecolor = C_JOINT,
    strokewidth = 0.8, bandwidth = 7.0
)
vlines!(axs, [0.0]; color = :grey30, linewidth = 0.6, linestyle = :dash)
xlims!(axs, t0, 8)
linkxaxes!(axc, axs)
rowsize!(fig.layout, 3, Relative(0.13))
rowgap!(fig.layout, 2, 2)
colsize!(fig.layout, 1, Relative(0.5))

panel_letter(fig[1, 1, TopLeft()], "A")
panel_letter(fig[1, 2, TopLeft()], "B")
panel_letter(fig[2, 1, TopLeft()], "C")

save(joinpath(FIGS, "fig-first.pdf"), fig; pt_per_unit = 1)
save(joinpath(FIGS, "fig-first.png"), fig; px_per_unit = 600 / 72)

# Numbers the caption and README quote, printed for checking.
println("cut-off: ", CUTOFF)
println(
    "draws: ", length(C_draws), "  median C_T: ", round(Int, q(0.5)),
    "  50%: ", round.(Int, (q(0.25), q(0.75))),
    "  90%: ", round.(Int, (q(0.05), q(0.95)))
)
seed_q = [
    CUTOFF - Day(round(Int, quantile(T_draws, p)))
        for p in (0.95, 0.5, 0.05)
]
println("seeding date 5/50/95%: ", seed_q)
