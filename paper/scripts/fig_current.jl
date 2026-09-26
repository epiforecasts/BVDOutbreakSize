# Figure 2 of the paper: the current model's outputs at one release, drawn
# from the release assets alone into paper/figures/fig-current.{pdf,png}.
#
#     julia --project=paper paper/scripts/fig_current.jl [tag]
#
# Panels:
#   A  national daily symptom onsets (30/60/90% bands) over the whole
#      window with the cut-off marked, and the reported confirmed-case
#      increments per day on a strip beneath;
#   B  reproduction number at the cut-off by province and nationally
#      (30/60/90% bars, line at one);
#   C  ascertainment by stream: posterior 30% and 90% intervals, the
#      median where the assets carry it, and the prior median;
#   D  one-week-ahead forecasts of new confirmed cases and deaths against
#      the weekly increments reported, for the validation origin a week
#      before the cut-off and for the cut-off itself.
#
# The assets carry no daily infection series (only symptom onsets, without a
# median column), no provincial reproduction-number series (only the
# cut-off value) and no median for the fraction tested, so those are not
# drawn.

include(joinpath(@__DIR__, "release_assets.jl"))

using CairoMakie
using Distributions: Beta

const TAG = isempty(ARGS) ? DEFAULT_TAG : ARGS[1]

ensure_assets(TAG)
meta = release_meta(TAG)

summary = read_asset("posterior_summary.csv")
draws = read_asset("posterior_draws.csv")
streams = read_asset("stream_estimates.csv")
forecast = read_asset("forecast.csv")
validation = read_asset("forecast_validation.csv")
onsets = read_asset("onsets_over_time.csv")
obs = TOML.parsefile(asset_path("observations.toml"))
cutoff = Date(obs["as_of_date"])

## Colours follow src/plots.jl and src/maps.jl: onsets seagreen, confirmed
## counts steelblue, the four patches as in the health-zone maps.
const ONSET_COLOUR = :seagreen
const CONFIRMED_COLOUR = :steelblue
const FORECAST_COLOUR = :firebrick
const PATCH_COLOURS = [:firebrick, :steelblue, :seagreen, :darkorange]
const PROVINCES = ["Ituri", "Nord-Kivu", "Haut-Uele", "Other provinces"]
const CUTOFF_LINE = (;
    color = (:black, 0.6), linestyle = :dash, linewidth = 0.8,
)

## Dates on a numeric axis: days since the epoch, ticked as "d Mon".
xday(d::Date) = Float64(Dates.value(d))
function dateticks(from::Date, to::Date; step = Month(1))
    t0 = Date(year(from), month(from), 1)
    t0 < from && (t0 += step)
    ticks = collect(t0:step:to)
    return (xday.(ticks), Dates.format.(ticks, "d u"))
end

## Nested 90/60/30% interval bars along x at row `y`; `v` holds the six
## bounds in `Lower 90%, Lower 60%, Lower 30%, Upper 30%, Upper 60%,
## Upper 90%` order.
function interval_bars!(ax, y, v, colour)
    for (lo, hi, lw) in ((1, 6, 1), (2, 5, 2.5), (3, 4, 4.5))
        rangebars!(
            ax, [y], [v[lo]], [v[hi]]; direction = :x, color = colour,
            linewidth = lw
        )
    end
    return nothing
end

## Open circle for a posterior median.
median_dot!(ax, x, y, colour) = scatter!(
    ax, [x], [y]; color = :white, strokecolor = colour, strokewidth = 1,
    markersize = 6
)

const SUMMARY_COLUMNS = (
    "Lower 90%", "Lower 60%", "Lower 30%", "Upper 30%", "Upper 60%",
    "Upper 90%",
)

## 180 mm wide at 1 pt per unit.
const WIDTH_PT = 180 / 25.4 * 72
fig = Figure(; size = (WIDTH_PT, 1.1 * WIDTH_PT), fontsize = 7.5)

panel_label!(pos, letter) = Label(
    pos, letter; font = :bold, fontsize = 10, halign = :left,
    padding = (0, 4, 0, 0), tellwidth = false, tellheight = false
)

## --- A: national symptom onsets and confirmed increments --------------------

ga = fig[1, 1:2] = GridLayout()
first_day = Date(2026, 4, 1)
keep = onsets.date .>= first_day
ox = xday.(onsets.date[keep])
axa = Axis(
    ga[1, 1]; ylabel = "Symptom onsets per day",
    xticks = dateticks(first_day, cutoff), xticklabelsvisible = false,
    xticksvisible = false
)
for (level, alpha) in (("90", 0.15), ("60", 0.28), ("30", 0.42))
    band!(
        axa, ox, Float64.(onsets[keep, "new_onsets_lower_$level"]),
        Float64.(onsets[keep, "new_onsets_upper_$level"]);
        color = (ONSET_COLOUR, alpha)
    )
end
vlines!(axa, [xday(cutoff)]; CUTOFF_LINE...)
ylims!(axa, 0, nothing)

cdates, cvals = history(obs, "confirmed_case_history")
rate = diff(cvals) ./ Dates.value.(diff(cdates))
axc = Axis(
    ga[2, 1]; xlabel = "Date (2026)", ylabel = "Confirmed cases\nper day",
    xticks = dateticks(first_day, cutoff)
)
scatter!(
    axc, xday.(cdates[2:end]), rate; color = CONFIRMED_COLOUR,
    markersize = 3.5, strokewidth = 0
)
vlines!(axc, [xday(cutoff)]; CUTOFF_LINE...)
ylims!(axc, 0, nothing)
linkxaxes!(axa, axc)
xlims!(axa, xday(first_day), xday(cutoff + Day(7)))
rowsize!(ga, 2, Relative(0.32))
rowgap!(ga, 4)
panel_label!(ga[1, 1, TopLeft()], "A")

## --- B: reproduction number at the cut-off by province -----------------

detail = Dict(headed_tables(site_page("estimates/province.html"), "h4"))
axb = Axis(
    fig[2, 1]; xlabel = "Reproduction number at the cut-off",
    yticks = (1:5, [PROVINCES; "National"]), yreversed = true
)
vlines!(axb, [1.0]; color = (:black, 0.5), linestyle = :dash, linewidth = 0.8)
for (i, label) in enumerate(PROVINCES)
    r = table_row([detail[label]], "Reproduction number")
    interval_bars!(axb, i, parse.(Float64, r[2:7]), PATCH_COLOURS[i])
end
rt_summary = summary_row(summary, "R_T")
interval_bars!(axb, 5, [rt_summary[k] for k in SUMMARY_COLUMNS], :black)
median_dot!(axb, stream_row(streams, "joint", "R_T").median, 5, :black)
xlims!(axb, 0.6, nothing)
panel_label!(fig[2, 1, TopLeft()], "B")

## --- C: ascertainment by stream ----------------------------------------
## Prior medians from src/models/priors.jl at the release commit, as in
## paper_numbers.jl.

priors_src = read(
    `git -C $(repo_dir()) show $(meta.commit):src/models/priors.jl`, String
)
asc_prior = parse(
    Float64,
    match(r"mu_prior = Normal\(logit\(([0-9.]+)\)", priors_src).captures[1]
)
mb = match(r"fraction_tested_prior = Beta\(([0-9.]+), ([0-9.]+)\)", priors_src)
tau_prior = median(
    Beta(parse(Float64, mb.captures[1]), parse(Float64, mb.captures[2]))
)
onset_row = table_row(
    html_tables(site_page("estimates/national.html")),
    "median modelled ascertainment"
)
## 90% and 30% bounds of one posterior_summary.csv quantity.
bounds(q) = [
    summary_row(summary, q)[k]
        for k in ("Lower 90%", "Lower 30%", "Upper 30%", "Upper 90%")
]
## Label, the four bounds, posterior median (or nothing) and prior median
## (or nothing) per stream.
asc_rows = [
    ("Suspected cases", bounds("p_drc"), median(draws.p_drc), asc_prior),
    ("Exported cases", bounds("p_uganda"), median(draws.p_uganda), asc_prior),
    ("Suspects tested", bounds("tau_test"), nothing, tau_prior),
    (
        "Onset reports",
        parse.(Float64, onset_row[[2, 4, 5, 7]]),
        nothing, nothing,
    ),
]
axcc = Axis(
    fig[2, 2]; xlabel = "Ascertainment (proportion)",
    yticks = (1:length(asc_rows), first.(asc_rows)), yreversed = true
)
prior_style = (;
    marker = :diamond, color = :white, strokecolor = FORECAST_COLOUR,
    strokewidth = 1.2, markersize = 8,
)
for (i, (_, b, m, p)) in enumerate(asc_rows)
    rangebars!(
        axcc, [i], [b[1]], [b[4]]; direction = :x, color = :black,
        linewidth = 1
    )
    rangebars!(
        axcc, [i], [b[2]], [b[3]]; direction = :x, color = :black,
        linewidth = 4
    )
    m === nothing || median_dot!(axcc, m, i, :black)
    p === nothing || scatter!(axcc, [p], [i]; prior_style...)
end
xlims!(axcc, 0, 1)
Legend(
    fig[2, 2],
    [
        MarkerElement(; prior_style...),
        MarkerElement(;
            marker = :circle, color = :white, strokecolor = :black,
            strokewidth = 1, markersize = 6
        ),
    ],
    ["Prior median", "Posterior median"];
    tellwidth = false, tellheight = false, halign = :left, valign = :top,
    framevisible = false, padding = (2, 2, 2, 2), rowgap = 0
)
panel_label!(fig[2, 2, TopLeft()], "C")

## --- D: one-week-ahead forecasts against the reported increments -------

gd = fig[3, 1:2] = GridLayout()
ddates, dvals = history(obs, "confirmed_death_history")
vmade = only(unique(validation.made_date))
fx_from = cutoff - Day(70)
fx_to = cutoff + Day(10)
obs_handle = fc_handle = nothing
for (j, (stream, ylabel, dates, values)) in enumerate(
        (
            (
                "confirmed cases", "New confirmed\ncases in 7 days",
                cdates, cvals,
            ),
            (
                "confirmed deaths", "New confirmed\ndeaths in 7 days",
                ddates, dvals,
            ),
        )
    )
    ax = Axis(
        gd[1, j]; xlabel = "Target date (2026)", ylabel,
        xticks = dateticks(fx_from, fx_to; step = Week(2))
    )
    ## Weekly increments reported, each on the week's last day.
    tx = Float64[]
    ty = Float64[]
    for k in 0:10
        target = cutoff - Day(7k)
        inc = increment_over(dates, values, target - Day(7), target)
        inc === missing && continue
        push!(tx, xday(target))
        push!(ty, inc)
    end
    for (df, made) in ((validation, vmade), (forecast, cutoff))
        v = forecast_draws(df, stream; made, horizon = 7)
        x = xday(made + Day(7))
        lo90, hi90 = q90(v)
        lo30, hi30 = q30(v)
        rangebars!(
            ax, [x], [lo90], [hi90]; color = FORECAST_COLOUR, linewidth = 1
        )
        rangebars!(
            ax, [x], [lo30], [hi30]; color = FORECAST_COLOUR, linewidth = 4
        )
        global fc_handle = median_dot!(ax, x, median(v), FORECAST_COLOUR)
    end
    global obs_handle = scatter!(ax, tx, ty; color = :black, markersize = 5)
    vlines!(ax, [xday(cutoff)]; CUTOFF_LINE...)
    xlims!(ax, xday(fx_from), xday(fx_to))
    ylims!(ax, 0, nothing)
end
Legend(
    gd[2, 1:2], [obs_handle, fc_handle],
    [
        "Reported in the 7 days to the target date",
        "Forecast median with 30% and 90% intervals",
    ];
    orientation = :horizontal, framevisible = false, padding = (0, 0, 0, 0)
)
rowgap!(gd, 2)
panel_label!(gd[1, 1, TopLeft()], "D")

rowsize!(fig.layout, 1, Relative(0.34))
rowsize!(fig.layout, 3, Relative(0.3))
rowgap!(fig.layout, 14)
colgap!(fig.layout, 10)

outdir = joinpath(paper_dir(), "figures")
mkpath(outdir)
save(joinpath(outdir, "fig-current.pdf"), fig; pt_per_unit = 1)
save(joinpath(outdir, "fig-current.png"), fig; px_per_unit = 600 / 72)
println("wrote fig-current.pdf and fig-current.png to ", outdir)
