# All package figures: posterior densities of `C_T`, posterior- and
# prior-predictive panel grids, pair plots, point-and-interval
# comparison, CFR prior, start-date and no-onward-transmission
# densities, and the one-week-ahead forecast figures.

"""
Overlaid posterior densities of `C_T` from one or more fits, built
through AlgebraOfGraphics. The 15 published scenario point estimates
are drawn as faint dashed Makie `vlines` on top of the AoG figure.
"""
function plot_cumulative_cases(
        streams::Pair{String, <:AbstractVector}...;
        scenarios = REPORT_SCENARIOS,
        xmax::Union{Nothing, Real} = nothing,
        xlabel::AbstractString = "Cumulative infections",
        title::AbstractString = "Outbreak size estimated by each data stream")
    upper = isnothing(xmax) ?
            1.05 * maximum(quantile(s.second, 0.995) for s in streams) :
            xmax
    df = @chain DataFrame(
        stream = String[], C_T = Float64[]
    ) begin
        let df = _
            for (label, draws) in streams
                for x in draws
                    0 < x < upper * 1.05 && push!(df, (label, float(x)))
                end
            end
            df
        end
    end

    spec = AoG.data(df) *
           AoG.mapping(:C_T => xlabel,
               color = :stream => "Data stream") *
           AoG.AlgebraOfGraphics.density() *
           AoG.subvisual(:line, linewidth = 2)
    fg = AoG.draw(spec;
        axis = (; ylabel = "Posterior density",
            title = title,
            limits = ((0, upper), nothing)),
        figure = (; size = (760, 420))
    )

    scenario_xs = Float64[val for (_, val) in scenarios if val < upper]
    isempty(scenario_xs) || vlines!(fg.figure.content[1], scenario_xs;
        color = (:grey, 0.4), linestyle = :dash)
    return fg
end

"""
Headline 3x2 cumulative figure. Rows are cumulative infections, cumulative
symptom onsets and cumulative deaths, all modelled BVD-only latent renewal
quantities (the deaths row is the BVD death series, excluding the non-BVD
background, so it stays as smooth as the infection and onset rows). The
left column is the modelled expected cumulative trajectory over the grid as
50% and 90% ribbons; the right column is the posterior
density of the current cut-off cumulative. The chain must carry the vector
deterministics `cumulative_infections`, `cumulative_onsets` and
`cumulative_expected_deaths` (one per draw). `seeding` is the calendar date
of grid day 1, so day `d` is `seeding + (d - 1)`. No observed data is
overlaid: each row is a latent quantity that sits upstream of ascertainment,
confirmation and reporting delays, so the observed counts are not on the
same scale.
"""
function plot_cumulative_trajectories(chn;
        n::Integer, seeding::Date)
    epoch = date2epochdays(seeding)
    x = Float64[epoch + (d - 1) for d in 1:n]

    ## Each trajectory deterministic is an iter×chain matrix of per-draw
    ## vectors; flatten to one vector of per-draw trajectories.
    function _trajectories(key)
        mat = chn[key]
        return [collect(v) for v in vec(collect(mat))]
    end
    function _ribbon(trajs)
        q(d, pr) = quantile(Float64[t[d] for t in trajs], pr)
        lo90 = [q(d, 0.05) for d in 1:n]
        hi90 = [q(d, 0.95) for d in 1:n]
        lo50 = [q(d, 0.25) for d in 1:n]
        hi50 = [q(d, 0.75) for d in 1:n]
        return lo90, hi90, lo50, hi50
    end

    rows = (
        (:cumulative_infections, "infections", :steelblue),
        (:cumulative_onsets, "symptom onsets", :seagreen),
        (:cumulative_expected_deaths, "deaths", :firebrick)
    )

    fig = Figure(; size = (940, 1020))
    for (i, (key, name, colour)) in enumerate(rows)
        trajs = _trajectories(key)
        lo90, hi90, lo50, hi50 = _ribbon(trajs)
        ax = Axis(fig[i, 1];
            xlabel = "Date", ylabel = "Cumulative $name",
            title = "Modelled cumulative $name over time",
            xticklabelrotation = pi / 6)
        band!(ax, x, lo90, hi90; color = (colour, 0.15))
        band!(ax, x, lo50, hi50; color = (colour, 0.30))
        loax = floor(Int, minimum(x))
        hiax = ceil(Int, maximum(x))
        ax.xticks = collect(loax:14:hiax)
        ax.xtickformat = vals -> [string(epochdays2date(round(Int, v)))
                                  for v in vals]

        finals = Float64[t[n] for t in trajs]
        axd = Axis(fig[i, 2];
            xlabel = "Cumulative $name at the cut-off",
            ylabel = "Posterior density",
            title = "Current cumulative $name")
        density!(axd, finals; color = (colour, 0.5),
            strokecolor = colour, strokewidth = 2)
    end
    return fig
end

"""
Overlaid cumulative-infection trajectories, one per single-stream fit, each
projected out to the cut-off on day `n` even when that stream's data stops
earlier. Each stream is drawn as 50% and 90% credible ribbons only, no
median line, matching the ribbon style of [`plot_rt`](@ref) and
[`plot_cumulative_trajectories`](@ref). A dotted vertical rule on each
stream's colour marks the date that stream's data stops reporting, so the
projection beyond the data reads apart from the fitted span.

Each `stream` is a `NamedTuple` `(; label, trajs, last_day, colour)`, where
`trajs` is a vector of per-draw cumulative-infection vectors of length `n`
(one per posterior draw) and `last_day` the 1-based grid day that stream's
data last reports (or `nothing` to omit the rule). `seeding` is the
calendar date of grid day 1, so day `d` is `seeding + (d - 1)`.
"""
function plot_stream_trajectories(streams::AbstractVector;
        n::Integer, seeding::Date,
        title::AbstractString =
        "Outbreak size projected to the cut-off by each data stream")
    epoch = date2epochdays(seeding)
    x = Float64[epoch + (d - 1) for d in 1:n]

    fig = Figure(; size = (900, 480))
    ax = Axis(fig[1, 1];
        xlabel = "Date", ylabel = "Cumulative infections",
        title = title, xticklabelrotation = pi / 6)

    handles = Any[]
    labels = String[]
    ymax = 0.0
    for s in streams
        trajs = s.trajs
        q(d, pr) = quantile(Float64[t[d] for t in trajs], pr)
        lo90 = [q(d, 0.05) for d in 1:n]
        hi90 = [q(d, 0.95) for d in 1:n]
        lo50 = [q(d, 0.25) for d in 1:n]
        hi50 = [q(d, 0.75) for d in 1:n]
        ymax = max(ymax, maximum(hi90))
        colour = s.colour
        band!(ax, x, lo90, hi90; color = (colour, 0.12))
        h = band!(ax, x, lo50, hi50; color = (colour, 0.30))
        push!(handles, h)
        push!(labels, s.label)
        ## Dotted rule in the stream's colour where its data stops reporting.
        ld = get(s, :last_day, nothing)
        ld === nothing || vlines!(ax, [Float64(epoch + ld - 1)];
            color = (colour, 0.8), linestyle = :dot, linewidth = 2)
    end

    loax = floor(Int, minimum(x))
    hiax = ceil(Int, maximum(x))
    ax.xticks = collect(loax:14:hiax)
    ax.xtickformat = vals -> [string(epochdays2date(round(Int, v)))
                              for v in vals]
    CairoMakie.ylims!(ax, 0, ymax * 1.05)
    CairoMakie.axislegend(ax, handles, labels; position = :lt,
        framevisible = true)
    return fig
end

"""
Overlaid posterior densities of an arbitrary scalar quantity from one
or more fits, built through AlgebraOfGraphics. Pass each fit as
`"label" => draws`; `xlabel` and `title` set the axis text.
"""
function plot_density_overlay(
        streams::Pair{String, <:AbstractVector}...;
        xlabel::AbstractString = "Value",
        title::AbstractString = "Posterior density")
    df = @chain DataFrame(stream = String[], value = Float64[]) begin
        let df = _
            for (label, draws) in streams
                for x in draws
                    push!(df, (label, float(x)))
                end
            end
            df
        end
    end

    spec = AoG.data(df) *
           AoG.mapping(:value => xlabel, color = :stream => "Fit") *
           AoG.AlgebraOfGraphics.density() *
           AoG.subvisual(:line, linewidth = 2)
    return AoG.draw(spec;
        axis = (; ylabel = "Posterior density", title = title),
        figure = (; size = (760, 420))
    )
end

_panel_pos(pos::Integer) = (1, pos)
_panel_pos(pos::Tuple) = pos

# Makie 0.24 (pulled in by AlgebraOfGraphics 0.12) computes data
# limits by calling `isfinite` elementwise, which has no method for
# integer vectors. Predictions of vector-valued observations (e.g.
# per-vintage `total_deaths`) arrive as a `Vector{Vector{Int}}`, so we
# flatten any nesting and convert to `Float64` before plotting. Scalar
# integer draws are floated for the same reason.
_pp_floats(pp::AbstractVector{<:Real}) = float.(pp)
function _pp_floats(pp::AbstractVector{<:AbstractVector})
    return Float64[float(x) for v in pp for x in v]
end
_pp_floats(pp) = Float64[float(x) for x in Iterators.flatten(pp)]

# Observed markers go through the same limit machinery, so float the
# scalar (or each element of a vector-valued observation).
_obs_floats(obs::Real) = Float64[float(obs)]
_obs_floats(obs::AbstractVector{<:Real}) = float.(obs)

function _panel_exports!(fig, pos, pp, obs; predictive_label = "Posterior")
    r, c = _panel_pos(pos)
    ppf = _pp_floats(pp)
    upper = max(20, ceil(Int, quantile(ppf, 0.99)))
    ax = Axis(fig[r, c];
        xlabel = "Replicated exported cases",
        ylabel = "$(predictive_label) predictive frequency",
        title = "Exports (cases)",
        limits = ((0, upper), nothing)
    )
    hist!(ax, ppf; bins = 0:1:upper, color = (:steelblue, 0.7))
    vlines!(ax, _obs_floats(obs); color = :red, linewidth = 2)
    return ax
end

function _panel_exports_deaths!(fig, pos, pp, obs;
        predictive_label = "Posterior")
    r, c = _panel_pos(pos)
    ppf = _pp_floats(pp)
    upper = max(3, ceil(Int, quantile(ppf, 0.995)))
    ax = Axis(fig[r, c];
        xlabel = "Replicated deaths among exports",
        ylabel = "$(predictive_label) predictive frequency",
        title = "Exports (deaths)",
        limits = ((0, upper), nothing)
    )
    hist!(ax, ppf; bins = 0:1:upper, color = (:rebeccapurple, 0.7))
    vlines!(ax, _obs_floats(obs); color = :red, linewidth = 2)
    return ax
end

function _panel_confirmed_deaths!(fig, pos, pp, obs;
        predictive_label = "Posterior")
    r, c = _panel_pos(pos)
    ppf = _pp_floats(pp)
    upper = max(3, ceil(Int, quantile(ppf, 0.995)))
    obs === nothing || (upper = max(upper, ceil(Int, 1.1 * obs)))
    ax = Axis(fig[r, c];
        xlabel = "Replicated confirmed deaths",
        ylabel = "$(predictive_label) predictive frequency",
        title = "Confirmed deaths (DRC)",
        limits = ((0, upper), nothing)
    )
    hist!(ax, ppf; bins = 0:1:upper, color = (:darkorange3, 0.7))
    obs === nothing || vlines!(ax, _obs_floats(obs);
        color = :red, linewidth = 2)
    return ax
end

function _panel_deaths!(fig, pos, pp, obs; predictive_label = "Posterior")
    r, c = _panel_pos(pos)
    ppf = _pp_floats(pp)
    upper = max(1.0, quantile(ppf, 0.995))
    ax = Axis(fig[r, c];
        xlabel = "Replicated deaths",
        ylabel = "$(predictive_label) predictive frequency",
        title = "Deaths (DRC)",
        limits = ((0, upper), nothing)
    )
    hist!(ax, ppf; bins = range(0, upper; length = 40),
        color = (:firebrick, 0.7))
    vlines!(ax, _obs_floats(obs); color = :red, linewidth = 2)
    return ax
end

function _panel_confirmed!(fig, pos, pp, obs;
        predictive_label = "Posterior")
    r, c = _panel_pos(pos)
    ppf = _pp_floats(pp)
    upper = max(1.0, quantile(ppf, 0.995))
    if obs !== nothing
        upper = max(upper, 1.05 * maximum(_obs_floats(obs)))
    end
    ax = Axis(fig[r, c];
        xlabel = "Replicated confirmed cases",
        ylabel = "$(predictive_label) predictive frequency",
        title = "Confirmed cases (DRC)",
        limits = ((0, upper), nothing)
    )
    hist!(ax, ppf; bins = range(0, upper; length = 40),
        color = (:goldenrod, 0.7))
    if obs !== nothing
        vlines!(ax, _obs_floats(obs); color = :red, linewidth = 2)
    end
    return ax
end

function _panel_tests!(fig, pos, pp, obs;
        predictive_label = "Posterior")
    r, c = _panel_pos(pos)
    ppf = _pp_floats(pp)
    upper = max(1.0, quantile(ppf, 0.995))
    if obs !== nothing
        upper = max(upper, 1.05 * maximum(_obs_floats(obs)))
    end
    ax = Axis(fig[r, c];
        xlabel = "Replicated tests analysed",
        ylabel = "$(predictive_label) predictive frequency",
        title = "Tests analysed (DRC)",
        limits = ((0, upper), nothing)
    )
    hist!(ax, ppf; bins = range(0, upper; length = 40),
        color = (:teal, 0.7))
    if obs !== nothing
        vlines!(ax, _obs_floats(obs); color = :red, linewidth = 2)
    end
    return ax
end

function _panel_cases!(fig, pos, pp, obs; predictive_label = "Posterior")
    r, c = _panel_pos(pos)
    ppf = _pp_floats(pp)
    upper = max(1.0, quantile(ppf, 0.995))
    ax = Axis(fig[r, c];
        xlabel = "Replicated reported cases",
        ylabel = "$(predictive_label) predictive frequency",
        title = "Reported cases (DRC)",
        limits = ((0, upper), nothing)
    )
    hist!(ax, ppf; bins = range(0, upper; length = 40),
        color = (:seagreen, 0.7))
    if obs !== nothing
        vlines!(ax, _obs_floats(obs); color = :red, linewidth = 2)
    end
    return ax
end

"""
Posterior predictive histogram with one panel per supplied data
stream. Pass `pp_exports`/`pp_deaths` as `nothing` to suppress
either of the first two panels, and supply `pp_cases` and/or
`pp_exports_deaths` to add the reported-cases and deaths-among-exports
panels. Observed values are drawn as red `vlines`. With four streams
the panels are laid out as a 2×2 grid (exports cases, exports deaths,
DRC deaths, DRC reported cases); fewer streams are placed in a single
row.
"""
function plot_posterior_predictive(
        pp_exports::Union{Nothing, AbstractVector},
        pp_deaths::Union{Nothing, AbstractVector},
        obs_exports::Union{Nothing, Real},
        obs_deaths::Union{Nothing, Real};
        pp_cases::Union{Nothing, AbstractVector} = nothing,
        obs_cases::Union{Nothing, Real} = nothing,
        pp_exports_deaths::Union{Nothing, AbstractVector} = nothing,
        obs_exports_deaths::Union{Nothing, Real} = nothing,
        pp_confirmed_deaths::Union{Nothing, AbstractVector} = nothing,
        obs_confirmed_deaths::Union{Nothing, Real} = nothing,
        pp_confirmed::Union{Nothing, AbstractVector} = nothing,
        obs_confirmed::Union{Nothing, Real} = nothing,
        pp_tests::Union{Nothing, AbstractVector} = nothing,
        obs_tests::Union{Nothing, Real} = nothing,
        predictive_label::AbstractString = "Posterior")
    panels = Tuple{Symbol, Any, Any}[]
    pp_exports === nothing ||
        push!(panels, (:exports, pp_exports, obs_exports))
    pp_exports_deaths === nothing ||
        push!(panels, (:exports_deaths, pp_exports_deaths,
            obs_exports_deaths))
    pp_confirmed_deaths === nothing ||
        push!(panels, (:confirmed_deaths, pp_confirmed_deaths,
            obs_confirmed_deaths))
    pp_deaths === nothing ||
        push!(panels, (:deaths, pp_deaths, obs_deaths))
    pp_cases === nothing ||
        push!(panels, (:cases, pp_cases, obs_cases))
    pp_tests === nothing ||
        push!(panels, (:tests, pp_tests, obs_tests))
    pp_confirmed === nothing ||
        push!(panels, (:confirmed, pp_confirmed, obs_confirmed))

    isempty(panels) && error(
        "plot_posterior_predictive needs at least one stream")

    ncols = length(panels) >= 4 ? 3 : length(panels)
    ncols = min(ncols, length(panels))
    nrows = cld(length(panels), ncols)
    fig = Figure(; size = (450 * ncols, 380 * nrows))
    for (i, (kind, pp, obs)) in enumerate(panels)
        pos = (cld(i, ncols), mod1(i, ncols))
        if kind === :exports
            _panel_exports!(fig, pos, pp, obs; predictive_label)
        elseif kind === :exports_deaths
            _panel_exports_deaths!(fig, pos, pp, obs; predictive_label)
        elseif kind === :confirmed_deaths
            _panel_confirmed_deaths!(fig, pos, pp, obs; predictive_label)
        elseif kind === :deaths
            _panel_deaths!(fig, pos, pp, obs; predictive_label)
        elseif kind === :cases
            _panel_cases!(fig, pos, pp, obs; predictive_label)
        elseif kind === :tests
            _panel_tests!(fig, pos, pp, obs; predictive_label)
        else
            _panel_confirmed!(fig, pos, pp, obs; predictive_label)
        end
    end
    return fig
end

## Panel painter for each stream key used by the comparison grid.
const _GRID_PANELS = (
    (:exports, _panel_exports!),
    (:exports_deaths, _panel_exports_deaths!),
    (:deaths, _panel_deaths!),
    (:cases, _panel_cases!),
    (:tests, _panel_tests!),
    (:confirmed, _panel_confirmed!)
)

"""
Two-row comparison of posterior-predictive distributions, one column
per stream. Top row: replicates from the per-stream fits. Bottom row:
replicates from the joint fit, conditioning on all observed streams.
Observed values shown as red vertical lines.

Each `NamedTuple` carries a subset of `(; exports, exports_deaths,
deaths, cases, tests, confirmed)`; columns are drawn in that canonical
order for whichever streams are present in `individual` (the
`confirmed`/`tests` columns appear only when the laboratory pipeline is
included). Each panel is a histogram of replicated counts; rows share
the same x-axis (the stream's count) so the per-stream and joint
predictives are directly comparable.
"""
function plot_posterior_predictive_grid(;
        individual::NamedTuple,
        joint::NamedTuple,
        observed::NamedTuple
)
    streams = [(key, painter)
               for (key, painter) in _GRID_PANELS
               if hasproperty(individual, key)]
    ncols = length(streams)
    fig = Figure(; size = (400 * ncols, 640))
    rows = ((:individual, individual, "per-stream fit"),
        (:joint, joint, "joint fit"))
    for (i, (_, pp, label)) in enumerate(rows)
        for (j, (key, painter)) in enumerate(streams)
            painter(fig, (i, j), getproperty(pp, key),
                getproperty(observed, key); predictive_label = label)
        end
    end
    return fig
end

"""
Prior predictive variant of `plot_posterior_predictive`, with the
panel labels switched to "Prior".
"""
function plot_prior_predictive(
        pp_exports::Union{Nothing, AbstractVector},
        pp_deaths::Union{Nothing, AbstractVector},
        obs_exports::Union{Nothing, Real},
        obs_deaths::Union{Nothing, Real};
        pp_cases::Union{Nothing, AbstractVector} = nothing,
        obs_cases::Union{Nothing, Real} = nothing,
        pp_confirmed::Union{Nothing, AbstractVector} = nothing,
        obs_confirmed::Union{Nothing, Real} = nothing,
        pp_tests::Union{Nothing, AbstractVector} = nothing,
        obs_tests::Union{Nothing, Real} = nothing)
    return plot_posterior_predictive(
        pp_exports, pp_deaths, obs_exports, obs_deaths;
        pp_cases, obs_cases, pp_confirmed, obs_confirmed,
        pp_tests, obs_tests, predictive_label = "Prior")
end

"""
PairPlots.jl corner plot over the named posterior parameters,
thinned by `thin`. Pass `prior` (another chain holding the same
parameters) to overlay the prior as a second series with a legend,
so the data's contribution to each marginal is visible.

`labels` is an optional map from the raw chain symbol to a clean display
name (e.g. `Symbol("rt_state.sigma_rw") => "Rt step size"`), applied to the
axis labels only; the model's variable names are unchanged. Symbols absent
from the map keep their raw name.
"""
function plot_pair(chn, params::AbstractVector{Symbol};
        thin::Integer = 2, prior = nothing,
        labels::AbstractDict = Dict{Symbol, String}())
    _name(p) = Symbol(get(labels, p, string(p)))
    _table(c) = DataFrame(
        NamedTuple(_name(p) => _draws(c, p) for p in params))[1:thin:end, :]
    post = _table(chn)
    prior === nothing && return PairPlots.pairplot(post)
    colours = CairoMakie.Makie.wong_colors()
    return PairPlots.pairplot(
        PairPlots.Series(post; label = "Posterior", color = colours[1]),
        PairPlots.Series(_table(prior); label = "Prior",
            color = colours[2])
    )
end

"""
Posterior correlation heatmap over the named scalar quantities `params`
(tracked deterministics or sampled parameters). Each cell is the Pearson
correlation of two quantities' posterior draws, drawn on a symmetric
red–blue scale with the value printed in the cell. Unlike the block-grouped
[`plot_pair`](@ref) corners — which keep the infection, delay and
surveillance parameters in separate plots — this surfaces the whole joint
identifiability structure in one panel, including the CROSS-block
degeneracies those corners split apart: the size–ascertainment seesaw
(`C_T` vs `p_drc`), the weaker size–fatality tilt (`C_T` vs `CFR`), and the
pooled `p_drc`–`p_uganda` link. `labels` maps a parameter symbol to a
display string, interpreted as LaTeX math (without the `\$` delimiters), so
`"p_\\mathrm{drc}"` renders with a subscript and `"\\lambda_\\mathrm{bg}"` as
`λ_bg`; a parameter absent from `labels` falls back to its symbol name.
Returns the `Figure`.
"""
function plot_correlation_heatmap(chn, params::AbstractVector{Symbol};
        labels::AbstractDict = Dict{Symbol, String}())
    ## Render tick labels as LaTeX so subscripts (R_T, p_drc, λ_bg) typeset
    ## properly; callers pass plain LaTeX math strings.
    name(p) = CairoMakie.Makie.latexstring(get(labels, p, string(p)))
    mat = reduce(hcat, (_draws(chn, p) for p in params))
    R = cor(mat)
    n = length(params)
    labs = [name(p) for p in params]
    fig = Figure(; size = (78 * n + 180, 78 * n + 140))
    ax = Axis(fig[1, 1]; xticks = (1:n, labs), yticks = (1:n, labs),
        xticklabelrotation = pi / 4, title = "Posterior correlation",
        aspect = 1)
    hm = CairoMakie.heatmap!(ax, 1:n, 1:n, R; colormap = :RdBu,
        colorrange = (-1, 1))
    for i in 1:n, j in 1:n

        CairoMakie.text!(ax, i, j;
            text = string(round(R[i, j]; digits = 2)),
            align = (:center, :center), fontsize = 10,
            color = abs(R[i, j]) > 0.6 ? :white : :black)
    end
    CairoMakie.Colorbar(fig[1, 2], hm)
    return fig
end

"""
Pairs plot of the per-stream modelled totals: a [`PairPlots`](https://sefffal.github.io/PairPlots.jl/)
corner of the per-draw modelled stream totals (`modelled`, a `NamedTuple`
of one per-draw vector per stream, each summed to that stream's own
observed support) with the observed totals (`observed`, a `NamedTuple` of
scalars) drawn as crosshair reference lines. The diagonals show how much
predictive density sits above or below each observed value (the
over/under-shoot), and the off-diagonals whether those shoots move together
across draws — the posterior-predictive view of data-stream conflict, the
counterpart to the parameter-space [`plot_correlation_heatmap`](@ref).
Returns the `Figure`.
"""
function plot_stream_pairs(modelled::NamedTuple, observed::NamedTuple)
    return PairPlots.pairplot(modelled,
        PairPlots.Truth(observed; label = "observed"))
end

"""
Estimate-evolution plot: how the outbreak-size estimate moves as the
data cut-off advances, drawn against the calendar date.

`released` is a vector of `(cutoff_date, median, lo30, hi30, lo60, hi60,
lo90, hi90)` tuples, one per published project release, drawn in blue.
`renewal` is the same tuple shape, the current renewal model re-fit
frozen at each release date, drawn in red: the current method evaluated
at a past cut-off. Each release and each frozen re-fit is its own
independent fit, so both are drawn as discrete per-date estimates — a
median marker with nested 30/60/90% vertical interval bars — rather than
a connected ribbon. Marks sharing a date (an integral and a renewal
release at one cut-off, or a release and its frozen re-fit) are dodged
horizontally so each reads as a separate estimate.

`trajectory` is the current-data, current-model cumulative-infection
trajectory over the day grid, a `(dates, lo30, hi30, lo60, hi60, lo90,
hi90)` tuple where `dates` is the calendar date of each grid day and the
remaining entries are the per-day credible bounds. It is a single fit
shown over time, so it is drawn in a third colour as one continuous
time-varying ribbon on the same calendar axis as the discrete marks.

Release dates are marked with dotted vertical rules.

`xlabel`/`ylabel`/`title` set the axis text; `released_label`,
`renewal_label` and `trajectory_label` name the three series.
"""
function plot_estimate_evolution(
        released::AbstractVector;
        renewal::AbstractVector = NamedTuple[],
        trajectory::Union{Nothing, Tuple} = nothing,
        xlabel::AbstractString = "Date",
        ylabel::AbstractString = "Cumulative infections",
        title::AbstractString = "Outbreak-size estimate as data accrued",
        released_label::AbstractString =
        "Released estimates (per project release)",
        renewal_label::AbstractString =
        "Current model re-fit frozen at each release date",
        trajectory_label::AbstractString =
        "Current model, current data",
        refline::Union{Nothing, Real} = nothing)
    ## Calendar dates → numeric day-offsets so the x-axis is to scale,
    ## then relabel the ticks with the dates. The released marks, the
    ## frozen renewal marks and the current-model trajectory all share
    ## this single calendar mapping, so they line up on the same axis.
    rdates = [Date(String(r[1])) for r in released]
    ndates = [Date(String(p[1])) for p in renewal]
    tdates = isnothing(trajectory) ? Date[] :
             [d isa Date ? d : Date(String(d)) for d in trajectory[1]]
    tickdates = sort(unique(vcat(rdates, ndates)))
    alldates = sort(unique(vcat(rdates, ndates, tdates)))
    ref = minimum(alldates)
    _x(d) = Float64((d - ref).value)

    _hi(tuples) = isempty(tuples) ? 0.0 : maximum(float(t[8]) for t in tuples)
    upper = max(_hi(released), _hi(renewal))
    isnothing(trajectory) ||
        (upper = max(upper, maximum(float.(trajectory[7]))))

    xlo = _x(ref) - 1
    xhi = _x(maximum(alldates)) + 1
    fig = Figure(; size = (860, 480))
    ax = Axis(fig[1, 1];
        xlabel = xlabel, ylabel = ylabel, title = title,
        xticks = (_x.(tickdates), [string(d) for d in tickdates]),
        xticklabelrotation = pi / 4,
        limits = ((xlo, xhi), (0, upper * 1.08)))

    ## Each release and each frozen re-fit is its own fit, so collect them
    ## as discrete marks and dodge any that share a date so they read
    ## apart rather than overplotting.
    marks = NamedTuple[]
    for (d, t) in zip(rdates, released)
        push!(marks, (; date = d, colour = :steelblue, t = t))
    end
    for (d, t) in zip(ndates, renewal)
        push!(marks, (; date = d, colour = :firebrick, t = t))
    end
    dodge = 0.7
    bydate = Dict{Date, Vector{Int}}()
    for (i, m) in enumerate(marks)
        push!(get!(bydate, m.date, Int[]), i)
    end
    markx = zeros(Float64, length(marks))
    for (d, idxs) in bydate, (k, i) in enumerate(idxs)

        markx[i] = _x(d) + (k - (length(idxs) + 1) / 2) * dodge
    end

    ## Vertical interval bars at one x: 90% (thin), 60%, 30% (thick).
    function _bars!(xs, los, his, colour, lw, alpha)
        bx = Float64[]
        by = Float64[]
        for (x, lo, hi) in zip(xs, los, his)
            append!(bx, (x, x))
            append!(by, (lo, hi))
        end
        return linesegments!(ax, bx, by;
            color = (colour, alpha), linewidth = lw)
    end

    ## One discrete series: nested 30/60/90% bars topped by a median dot.
    function _series!(colour)
        sel = [i for i in eachindex(marks) if marks[i].colour == colour]
        isempty(sel) && return nothing
        xs = markx[sel]
        ts = [marks[i].t for i in sel]
        _bars!(xs, [float(t[7]) for t in ts], [float(t[8]) for t in ts],
            colour, 1.4, 0.45)
        _bars!(xs, [float(t[5]) for t in ts], [float(t[6]) for t in ts],
            colour, 3.2, 0.55)
        _bars!(xs, [float(t[3]) for t in ts], [float(t[4]) for t in ts],
            colour, 6.5, 0.70)
        return scatter!(ax, xs, [float(t[2]) for t in ts];
            color = colour, markersize = 9,
            strokecolor = :white, strokewidth = 1)
    end

    handles = Any[]
    labels = String[]
    ## Current-data estimate as the cumulative-infection trajectory over the
    ## grid: a single fit, so one continuous ribbon on the same calendar
    ## axis. `trajectory` is `(dates, lo30, hi30, lo60, hi60, lo90, hi90)`.
    if !isnothing(trajectory)
        cc = :seagreen
        xs = _x.(tdates)
        ord = sortperm(xs)
        xs = xs[ord]
        band!(ax, xs, float.(trajectory[6])[ord], float.(trajectory[7])[ord];
            color = (cc, 0.10))
        band!(ax, xs, float.(trajectory[4])[ord], float.(trajectory[5])[ord];
            color = (cc, 0.16))
        th = band!(ax, xs, float.(trajectory[2])[ord],
            float.(trajectory[3])[ord]; color = (cc, 0.24))
        push!(handles, th)
        push!(labels, trajectory_label * " (30/60/90% band)")
    end
    renewal_mark = _series!(:firebrick)
    if !isnothing(renewal_mark)
        push!(handles, renewal_mark)
        push!(labels, renewal_label * " (median, 30/60/90% bars)")
    end
    released_mark = _series!(:steelblue)
    if !isnothing(released_mark)
        push!(handles, released_mark)
        push!(labels, released_label * " (median, 30/60/90% bars)")
    end

    ## Dotted vertical rule at each release date.
    isempty(rdates) || vlines!(ax, _x.(rdates);
        color = (:grey, 0.55), linestyle = :dot, linewidth = 1)

    ## Optional horizontal reference line, e.g. Rt = 1 for a reproduction
    ## number, drawn faint so it reads behind the estimates.
    isnothing(refline) ||
        hlines!(ax, [float(refline)];
            color = (:black, 0.4), linestyle = :dash, linewidth = 1)

    CairoMakie.axislegend(ax, handles, labels; position = :lt,
        framevisible = true)
    return fig
end

"""
Forecasts-versus-now overlay: one panel per observed stream, each showing the
forecasts made at every release (median with the 90% predictive interval as a
vertical bar, coloured by horizon) against the value observed since (black
points). `overlay` is the `data/forecast_overlay.csv` table restricted to the
`"ours"` model rows, with columns `stream`, `target_date`, `horizon`,
`observed`, `median`, `lo90` and `hi90`. Returns a figure carrying a short
note in place of the panels when no forecasts have been scored yet.
"""
function plot_forecast_overlay(overlay::DataFrame)
    streams = unique(overlay.stream)
    ## Scores accrue only once a release stores a forecast, so an empty
    ## table is the expected early state and says so rather than
    ## returning a blank panel.
    if isempty(streams)
        fig = Figure(; size = (860, 160))
        CairoMakie.Label(fig[1, 1],
            "No forecasts scored yet. No release carries a stored forecast.";
            tellwidth = false, tellheight = false, color = (:black, 0.55))
        return fig
    end
    fig = Figure(; size = (860, 240 * length(streams) + 50))
    horizons = sort(unique(overlay.horizon))
    palette = [:steelblue, :seagreen, :goldenrod, :firebrick]
    hcol = Dict(h => palette[mod1(i, length(palette))]
    for (i, h) in enumerate(horizons))
    ref = minimum(Date.(string.(overlay.target_date)))
    _x(d) = Float64((Date(string(d)) - ref).value)
    handles = Any[]
    labels = String[]
    for (row, s) in enumerate(streams)
        sub = overlay[overlay.stream .== s, :]
        ax = Axis(fig[row, 1]; title = s, xlabel = "Target date",
            ylabel = "Count")
        ## Observed value at each target date, once per date.
        od = sort(unique([(Date(string(t)), float(o))
                          for (t, o) in zip(sub.target_date, sub.observed)]))
        oh = scatter!(ax, [_x(t) for (t, _) in od], [o for (_, o) in od];
            color = :black, markersize = 7)
        if row == 1
            push!(handles, oh)
            push!(labels, "observed")
        end
        for h in horizons
            hs = sub[sub.horizon .== h, :]
            isempty(hs) && continue
            xs = [_x(t) for t in hs.target_date]
            bx = Float64[]
            by = Float64[]
            for (x, lo, hi) in zip(xs, hs.lo90, hs.hi90)
                append!(bx, (x, x))
                append!(by, (float(lo), float(hi)))
            end
            linesegments!(ax, bx, by; color = (hcol[h], 0.4), linewidth = 2)
            mh = scatter!(ax, xs, float.(hs.median);
                color = hcol[h], markersize = 6)
            if row == 1
                push!(handles, mh)
                push!(labels, "$(h)-day forecast")
            end
        end
        ts = sort(unique(Date.(string.(sub.target_date))))
        ax.xticks = (_x.(ts), string.(ts))
        ax.xticklabelrotation = pi / 4
    end
    ## Legend below the panels, since inside the first axis it covers the
    ## data and the tick labels once several horizons are scored.
    isempty(handles) ||
        CairoMakie.Legend(fig[length(streams) + 1, 1], handles, labels;
            orientation = :horizontal, framevisible = true,
            tellheight = true, tellwidth = false)
    return fig
end

"""
Horizontal point-and-interval comparison of cumulative-case estimates
from several sources. `rows` is a vector of
`(label, central, lower, upper)` tuples, drawn top to bottom with the
central estimate as a point and `[lower, upper]` as a bar. A row whose
lower and upper match its central is a deterministic point estimate and is
drawn as a bare marker with no bar.

`groups` is an optional vector of group keys, one per row, matched against
`group_colours` (a vector of `key => colour` pairs) to colour each row's
marker and bar and build a legend, so several sources read apart at a
glance. Without `groups` the rows share a single colour.
"""
function plot_estimate_comparison(
        rows::AbstractVector;
        xlabel::AbstractString = "Cumulative cases",
        xmax::Union{Nothing, Real} = nothing,
        groups::Union{Nothing, AbstractVector} = nothing,
        group_colours::AbstractVector = Pair[])
    n = length(rows)
    labels = [String(r[1]) for r in rows]
    central = [float(r[2]) for r in rows]
    lo = [float(r[3]) for r in rows]
    hi = [float(r[4]) for r in rows]
    top = isnothing(xmax) ? maximum(hi) * 1.08 : xmax

    cmap = Dict(group_colours)
    _colour(i) = isnothing(groups) ? :steelblue :
                 get(cmap, groups[i], :steelblue)

    fig = Figure(; size = (840, 120 + 46n))
    ax = Axis(fig[1, 1];
        xlabel = xlabel,
        yticks = (collect(1:n), reverse(labels)),
        limits = ((0, top), (0.5, n + 0.5))
    )
    for i in 1:n
        y = n - i + 1
        col = _colour(i)
        ## A deterministic point estimate has no interval, so draw a bare
        ## marker; otherwise draw the bar with the central point on top.
        if hi[i] > lo[i]
            lines!(ax, [lo[i], hi[i]], [y, y];
                color = (col, 0.8), linewidth = 3)
        end
        scatter!(ax, [central[i]], [y]; color = col, markersize = 12)
    end
    if !isnothing(groups) && !isempty(group_colours)
        handles = [CairoMakie.MarkerElement(; color = c, marker = :circle,
                       markersize = 12) for (_, c) in group_colours]
        glabels = [String(k) for (k, _) in group_colours]
        CairoMakie.axislegend(ax, handles, glabels; position = :rb,
            framevisible = true)
    end
    return fig
end

"""
Calendar time-series comparison of cumulative-count projections against the
data observed since. `external` is another group's published projection and
`ours` is our own forward projection, each drawn as a central line with a
shaded `[lower, upper]` band so the two fans read directly against each other;
`observed` is the data observed so far, drawn as a marked line. Each is a
vector of `(date, ...)` tuples with `date` an ISO string: `external` and `ours`
are `(date, central, lower, upper)`, `observed` is `(date, value)`. The dates
share one calendar x-axis, so the two projections read directly against what
the outbreak actually did. Used to set our forward projection beside the
Chamla et al. [chamla2026](@cite) confirmed-case projection.
"""
function plot_projection_comparison(;
        external::AbstractVector,
        ours::AbstractVector,
        observed::AbstractVector,
        external_label::AbstractString = "External projection",
        ours_label::AbstractString = "Our projection",
        observed_label::AbstractString = "Observed",
        external_colour = :steelblue,
        ours_colour = :firebrick,
        ylabel::AbstractString = "Cumulative confirmed cases",
        title::AbstractString = "Projected versus observed cumulative cases")
    _x(d) = Float64(date2epochdays(Date(String(d))))

    fig = Figure(; size = (820, 460))
    ax = Axis(fig[1, 1]; xlabel = "Date", ylabel = ylabel, title = title)

    ## External projection: shaded 90% band, central line and markers.
    ex_x = [_x(r[1]) for r in external]
    ex_m = [float(r[2]) for r in external]
    ex_lo = [float(r[3]) for r in external]
    ex_hi = [float(r[4]) for r in external]
    ord = sortperm(ex_x)
    band!(ax, ex_x[ord], ex_lo[ord], ex_hi[ord];
        color = (external_colour, 0.15))
    lines!(ax, ex_x[ord], ex_m[ord]; color = external_colour, linewidth = 2)
    ex_h = scatter!(ax, ex_x, ex_m; color = external_colour, markersize = 8)

    ## Observed data so far: a marked black line.
    ob_x = [_x(r[1]) for r in observed]
    ob_y = [float(r[2]) for r in observed]
    obord = sortperm(ob_x)
    lines!(ax, ob_x[obord], ob_y[obord]; color = :black, linewidth = 1.5)
    ob_h = scatter!(ax, ob_x, ob_y; color = :black, markersize = 7)

    ## Our forward projection: shaded 90% band, central line and markers, the
    ## same ribbon form as the external projection so the two fans read against
    ## each other.
    our_x = [_x(r[1]) for r in ours]
    our_m = [float(r[2]) for r in ours]
    our_lo = [float(r[3]) for r in ours]
    our_hi = [float(r[4]) for r in ours]
    oord = sortperm(our_x)
    band!(ax, our_x[oord], our_lo[oord], our_hi[oord];
        color = (ours_colour, 0.15))
    lines!(ax, our_x[oord], our_m[oord]; color = ours_colour, linewidth = 2)
    our_h = scatter!(ax, our_x, our_m;
        color = ours_colour, markersize = 8, marker = :diamond)

    allx = vcat(ex_x, ob_x, our_x)
    lo, hi = minimum(allx), maximum(allx)
    ax.xticks = collect(lo:14:hi)
    ax.xtickformat = vals -> [string(epochdays2date(round(Int, v)))
                              for v in vals]
    CairoMakie.axislegend(ax, [ex_h, our_h, ob_h],
        [external_label, ours_label, observed_label];
        position = :lt, framevisible = true)
    return fig
end

"""
Faceted point-and-interval comparison of published scenario estimates, one
panel per date of estimation so the figure spreads sideways instead of into one
tall column. `scenarios` is `REPORT_SCENARIOS_CI`-shaped
`(date, label, central, lower, upper)` with `label` of the form
`"M1|M2 <family>, <swept level>"` (for example `"M2 τ=14 d, CFR 26%"`). Within a
panel each `(method, family)` is one row, and the swept nuisance level (the CFR,
window or doubling time) is dodged onto that single line so every scenario keeps
its own interval while the sweep no longer adds rows. Method sets the colour.
`ours` maps a date string to our matched `(median, lower, upper)` estimate, drawn
as a grey reference band with a dashed median in that date's panel.
`date_titles` are `date => title` pairs giving each panel its heading.
"""
function plot_scenario_comparison(scenarios::AbstractVector;
        ours::AbstractDict = Dict{String, Any}(),
        date_titles::AbstractVector = ["2026-05-18" => "18 May report",
            "2026-05-20" => "20 May update",
            "2026-05-27" => "27 May (Lancet)"],
        method_names = Dict("M1" => "geographic", "M2" => "back-calc"),
        method_colours = Dict("M1" => :steelblue, "M2" => :darkorange),
        xlabel::AbstractString = "Cumulative cases",
        title::AbstractString = "Published scenarios versus our estimate")
    title_of = Dict(date_titles)
    dates = sort(unique(String[String(s[1]) for s in scenarios]))

    ## Parse "M2 τ=14 d, CFR 26%" → (method, family, level).
    function parts(label)
        head, level = split(String(label), ", "; limit = 2)
        method, family = split(head, " "; limit = 2)
        return (String(method), String(family), String(level))
    end

    ## First pass: group each date's scenarios into ordered (method, family)
    ## rows, so a shared row count keeps the panels' rows the same height.
    function group(date)
        fams = Tuple{String, String}[]
        members = Dict{Tuple{String, String}, Vector{Any}}()
        for s in scenarios
            String(s[1]) == date || continue
            m, fam, lvl = parts(s[2])
            key = (m, fam)
            haskey(members, key) || (push!(fams, key); members[key] = [])
            push!(members[key], (lvl, float(s[3]), float(s[4]), float(s[5])))
        end
        return fams, members
    end
    grouped = Dict(d => group(d) for d in dates)
    ## Reserve the geographic rows at the top and the back-calculation rows at
    ## the bottom of every panel, sized to the busiest panel, so the two method
    ## blocks line up across dates even where a panel has no geographic row.
    countm(d, m) = count(f -> f[1] == m, first(grouped[d]))
    maxgeo = maximum(countm(d, "M1") for d in dates)
    maxbc = maximum(countm(d, "M2") for d in dates)
    maxrow = maxgeo + maxbc

    xmax = 1.05 * max(maximum(float(s[5]) for s in scenarios),
        maximum((float(v[3]) for v in values(ours)); init = 0.0))

    fig = Figure(; size = (340 * length(dates) + 140, 110 + 70 * maxrow))
    for (j, d) in enumerate(dates)
        fams, members = grouped[d]
        geo = [f for f in fams if f[1] == "M1"]
        bc = [f for f in fams if f[1] == "M2"]
        ## Geographic families fill the top block (`maxrow` down); back-calc
        ## families fill the bottom block (`maxbc` down), so both align in y.
        ypos = Dict{Tuple{String, String}, Int}()
        for (i, f) in enumerate(geo)
            ypos[f] = maxrow - i + 1
        end
        for (i, f) in enumerate(bc)
            ypos[f] = maxbc - i + 1
        end
        yvals = [ypos[f] for f in fams]
        ylabels = [string(get(method_names, m, m), " · ", fam)
                   for (m, fam) in fams]
        ax = Axis(fig[1, j];
            title = get(title_of, d, d), xlabel = xlabel,
            yticks = (yvals, ylabels),
            limits = ((0, xmax), (0.4, maxrow + 0.6)))

        ## Our matched estimate for this vintage: a reference band + median.
        if haskey(ours, d)
            med, lo, hi = ours[d]
            vspan!(ax, float(lo), float(hi); color = (:grey, 0.18))
            vlines!(ax, [float(med)];
                color = :black, linestyle = :dash, linewidth = 1.5)
        end

        ## Each family row carries its swept levels dodged around the row centre,
        ## each a point with its interval, coloured by method.
        for key in fams
            y = ypos[key]
            col = get(method_colours, key[1], :grey)
            ms = members[key]
            k = length(ms)
            offs = k == 1 ? [0.0] : collect(LinRange(0.26, -0.26, k))
            for (t, (_, c, lo, hi)) in enumerate(ms)
                yy = y + offs[t]
                lines!(ax, [lo, hi], [yy, yy];
                    color = (col, 0.85), linewidth = 2.5)
                scatter!(ax, [c], [yy]; color = col, markersize = 9)
            end
        end
    end

    CairoMakie.Label(fig[0, 1:length(dates)], title;
        fontsize = 16, font = :bold)
    handles = [
        CairoMakie.MarkerElement(; color = method_colours["M1"],
            marker = :circle, markersize = 11),
        CairoMakie.MarkerElement(; color = method_colours["M2"],
            marker = :circle, markersize = 11),
        CairoMakie.PolyElement(; color = (:grey, 0.4))]
    labels = [get(method_names, "M1", "M1") * " spread",
        get(method_names, "M2", "M2") * " spread", "our estimate (90%)"]
    CairoMakie.Legend(fig[2, 1:length(dates)], handles, labels;
        orientation = :horizontal, framevisible = false)
    return fig
end

"""
Density of a prior over the case-fatality ratio (CFR) on `[0, 1]`,
plotted on the sub-range `[0, 0.7]`. The CDC central estimate of
55/169 ≈ 0.33 is drawn as a solid vertical rule, and the report's 26%
and 40% scenario bounds as dashed rules, so the prior can be read
against the published CFR scenarios.
"""
function plot_cfr_prior(prior::Distribution)
    colours = CairoMakie.Makie.wong_colors()
    xs = range(0.0, 0.7; length = 400)
    ys = pdf.(Ref(prior), xs)

    fig = Figure(; size = (760, 420))
    ax = Axis(fig[1, 1];
        xlabel = "Case-fatality ratio (CFR)",
        ylabel = "Prior density",
        title = "Prior over the case-fatality ratio",
        limits = ((0, 0.7), nothing)
    )
    lines!(ax, xs, ys; color = colours[1], linewidth = 2)
    vlines!(ax, [55 / 169]; color = :firebrick, linewidth = 2)
    vlines!(ax, [0.26, 0.40];
        color = (:grey, 0.6), linestyle = :dash, linewidth = 2)
    return fig
end

"""
Posterior densities of the delay-corrected confirmed case-fatality ratio
and the structural (infection-based) CFR from a
[`delay_corrected_confirmed_cfr`](@ref) result `res`, with the naive observed
confirmed ratio drawn as a solid vertical rule and the median uncorrected
modelled confirmed ratio as a dashed rule. The gap between the naive rule and
the corrected density shows the real-time delay debiasing; the gap to the
structural density shows the residual case/death ascertainment difference.
Plotted on the CFR percentage scale.
"""
function plot_confirmed_cfr(res)
    colours = CairoMakie.Makie.wong_colors()
    corrected = 100 .* filter(isfinite, res.corrected)
    structural = 100 .* filter(isfinite, res.structural)
    naive = 100 * res.naive_observed
    modelled_naive = 100 * quantile(filter(isfinite, res.modelled_naive), 0.5)

    hi = max(maximum(corrected), maximum(structural), naive) * 1.05
    fig = Figure(; size = (760, 420))
    ax = Axis(fig[1, 1];
        xlabel = "Case-fatality ratio (%)",
        ylabel = "Posterior density",
        title = "Delay-corrected confirmed CFR versus the structural CFR",
        limits = ((0, hi), nothing)
    )
    h_corr = density!(ax, corrected; color = (colours[1], 0.5),
        strokecolor = colours[1], strokewidth = 2)
    h_struct = density!(ax, structural; color = (colours[2], 0.4),
        strokecolor = colours[2], strokewidth = 2)
    h_naive = vlines!(ax, [naive]; color = :firebrick, linewidth = 2)
    h_mod = vlines!(ax, [modelled_naive];
        color = (:grey, 0.7), linestyle = :dash, linewidth = 2)
    CairoMakie.axislegend(ax,
        [h_corr, h_struct, h_naive, h_mod],
        ["Delay-corrected confirmed CFR", "Structural CFR",
            "Naive observed confirmed ratio",
            "Uncorrected modelled confirmed ratio (median)"];
        position = :rt, framevisible = true)
    return fig
end

"""
One-row, two-panel figure summarising when the outbreak began. The
left panel is the posterior density of the outbreak start date, the
calendar date of the import that started the outbreak, obtained by
rescaling the outbreak age `T` (origin to cut-off) to a calendar date
(`as_of_date` minus `T`). The right panel is the joint
`(doubling_time, T)` posterior pair plot: shorter doubling times
correspond to faster early growth, which reaches the same epidemic
size in less time (smaller `T`).
"""
function plot_start_date_pair(chn;
        as_of_date::AbstractString, thin::Integer = 2)
    T_draws = _draws(chn, :T)
    cutoff_days = date2epochdays(Date(as_of_date))
    start_days = cutoff_days .- T_draws

    fig = Figure(; size = (1100, 460))
    ax = Axis(fig[1, 1];
        xlabel = "Outbreak start date",
        ylabel = "Posterior density",
        title = "Estimated outbreak start date",
        xticklabelrotation = π / 6
    )
    density!(ax, start_days; color = (:steelblue, 0.5),
        strokecolor = :steelblue, strokewidth = 2)
    ## Date ticks every four weeks across the posterior range, so the
    ## start date stays readable rather than relying on the default
    ## locator or crowding the axis as the range widens.
    lo = floor(Int, minimum(start_days))
    hi = ceil(Int, maximum(start_days))
    ax.xticks = collect(lo:28:hi)
    ax.xtickformat = vals -> [string(epochdays2date(round(Int, v))) for v in vals]

    dt_draws = _draws(chn, :doubling_time)
    ## Clip extreme doubling times (near-zero growth) to keep the pair
    ## plot readable; a finite cap at 200 days covers the credible range.
    dt_clipped = clamp.(dt_draws, -200.0, 200.0)
    pair_df = DataFrame(doubling_time = dt_clipped, T = T_draws)
    PairPlots.pairplot(fig[1, 2], pair_df[1:thin:end, :])
    return fig
end

"""
Reconstruct each posterior draw's daily reproduction-number trajectory
`Rt` from the sampled weekly random-walk parameters, returning a
`ndraws × n` matrix masked to each draw's established window (`missing`
before `rt_start`). The saved chain stores only the cut-off `R_T`, so each
draw's daily `Rt` is rebuilt by mirroring [`rt_walk_model`](@ref): weekly
knots ([`knot_days`](@ref)) from `rt_walk_start` follow a non-centred
Gaussian walk (`rt_state.log_R0` plus the cumulative sum of
`rt_state.sigma_rw .* rt_state.z`), linearly interpolated to the day grid
([`interpolate_knots`](@ref)) and shifted by the sampled
`rt_state.intervention_effect` along a logistic ramp
([`sigmoid_ramp`](@ref)) centred at the outbreak-response `breakpoint`.
Shared by [`plot_rt`](@ref) and [`plot_rt_streams`](@ref).
"""
function reconstruct_rt(chn; n::Integer, breakpoint::Real,
        rt_start::Integer = 1, rt_walk_start::Integer = rt_start,
        week::Integer = 7, ramp::Real = 14.0)
    log_R0 = _draws(chn, Symbol("rt_state.log_R0"))
    sigma = _draws(chn, Symbol("rt_state.sigma_rw"))
    effect = _draws(chn, Symbol("rt_state.intervention_effect"))
    ## `rt_state.z` is vector-valued: one standard-normal innovation vector
    ## per draw. Pull each draw's full vector from the chain slice.
    zmat = chn[Symbol("rt_state.z")]
    zrows = [collect(z) for z in vec(collect(zmat))]

    ## The knot grid is built from the model's WALK start `rt_walk_start`
    ## (the first situation report, the breakpoint grid day), which is
    ## decoupled from `rt_start` (the renewal/established-window start used for
    ## the mask below). The innovation vector length is fixed by that walk
    ## start; if the caller passes a mismatching one the knot grid here will
    ## not match, so fail with a clear message rather than a downstream bounds
    ## error.
    days = knot_days(n; week, start = rt_walk_start)
    nb = length(days)
    if !isempty(zrows) && length(zrows[1]) != nb - 1
        error("reconstruct_rt: rt_walk_start = $rt_walk_start gives " *
              "$(nb - 1) random-walk steps but the chain has " *
              "$(length(zrows[1])); pass the same walk start the model used " *
              "(the breakpoint grid day, n - who_first_sitrep_days).")
    end
    ramp_shape = sigmoid_ramp(n, breakpoint; ramp)
    ndraws = length(log_R0)

    ## Per-draw daily Rt, masked to the draw's own established window
    ## (cumulative infections ≥ 1, i.e. grid day ≥ n - round(T)).
    rt = Matrix{Union{Missing, Float64}}(missing, ndraws, n)
    for i in 1:ndraws
        z = zrows[i]
        steps = sigma[i] .* z[1:(nb - 1)]
        log_R = log_R0[i] .+ vcat(0.0, cumsum(steps))
        walk = interpolate_knots(log_R, days, n)
        ## Days before the renewal start clamp to the established R0 (the
        ## walk base); they are filled by the analytic cryptic exponential in
        ## the model and are not plotted (masked from `rt_start` onward).
        log_Rt = walk .+ effect[i] .* ramp_shape
        start = clamp(rt_start, 1, n)
        for d in start:n
            rt[i, d] = exp(log_Rt[d])
        end
    end
    return rt
end

## Per-day quantile `pr` of an established-window Rt matrix, skipping the
## masked (pre-renewal) days; `missing` where a day has no established draws.
function _rt_quantile(rt::AbstractMatrix, d::Integer, pr::Real)
    col = collect(skipmissing(@view rt[:, d]))
    return isempty(col) ? missing : quantile(col, pr)
end

"""
Reconstruct the daily reproduction-number trajectory `Rt` per posterior
draw from the sampled weekly random-walk parameters and plot it over the
established-outbreak window. The saved chain stores only the cut-off
`R_T`, so each draw's daily `Rt` is rebuilt by mirroring
[`rt_walk_model`](@ref): weekly knots ([`knot_days`](@ref)) follow a
non-centred Gaussian walk (`rt_state.log_R0` plus the cumulative sum of
`rt_state.sigma_rw .* rt_state.z`), linearly interpolated to the day grid
([`interpolate_knots`](@ref)) and shifted by the sampled
`rt_state.intervention_effect` along a logistic ramp
([`sigmoid_ramp`](@ref)) centred at the outbreak-response `breakpoint`.

The estimated window runs from `rt_start` (the renewal start, where the
random walk begins) to the cut-off; only that period is drawn, the median
with 50% and 90% ribbons, and
about `n_traj` thinned sampled trajectories are overlaid thin and faint to
show the per-draw spread. The intervention breakpoint, the end of the
intervention scale-up (`breakpoint + ramp`) as a dotted rule, and the cut-off
are marked. `seeding` is the calendar date of grid day 1 (so day `d` is
`seeding + (d - 1)`).
"""
function plot_rt(chn; n::Integer, breakpoint::Real,
        as_of_date::AbstractString, seeding::Date,
        rt_start::Integer = 1, rt_walk_start::Integer = rt_start,
        week::Integer = 7, ramp::Real = 14.0,
        n_traj::Integer = 100)
    rt = reconstruct_rt(chn; n, breakpoint, rt_start, rt_walk_start, week, ramp)
    ndraws = size(rt, 1)

    ## Median and ribbons over established draws only (skip masked days).
    q(d, pr) = _rt_quantile(rt, d, pr)
    med = [q(d, 0.5) for d in 1:n]
    lo90 = [q(d, 0.05) for d in 1:n]
    hi90 = [q(d, 0.95) for d in 1:n]
    lo50 = [q(d, 0.25) for d in 1:n]
    hi50 = [q(d, 0.75) for d in 1:n]
    est = findall(!ismissing, med)

    ## 30% inner band for the third credible level alongside the 50/90.
    lo30 = [q(d, 0.35) for d in 1:n]
    hi30 = [q(d, 0.65) for d in 1:n]

    epoch = date2epochdays(seeding)
    x = [epoch + (d - 1) for d in 1:n]
    xe = x[est]
    ## Cap the y-axis just above the upper 90% credible band so the focus
    ## stays on the Rt trajectory: a 20% headroom keeps the band off the top
    ## edge without leaving a wide empty strip, and the handful of higher
    ## sampled trajectories above the cap are simply clipped.
    hi90_est = Float64[hi90[d] for d in est if !ismissing(hi90[d])]
    ytop = isempty(hi90_est) ? 4.0 :
           max(1.2, 1.2 * maximum(hi90_est))
    fig = Figure(; size = (900, 440))
    ax = Axis(fig[1, 1]; xlabel = "Date", ylabel = "Reproduction number Rt",
        title = "Estimated Rt over the established outbreak",
        limits = (nothing, (0.0, ytop)),
        xticklabelrotation = pi / 6)
    ## Thin sampled trajectories over the estimated window, faint, so the
    ## per-draw spread reads alongside the ribbons.
    if n_traj > 0 && !isempty(est)
        step = max(1, fld(ndraws, n_traj))
        for i in 1:step:ndraws
            yi = Float64[rt[i, d] for d in est]
            lines!(ax, xe, yi; color = (:purple, 0.15), linewidth = 0.5)
        end
    end
    ## 30/60/90% credible ribbons over the established window, no median line.
    band!(ax, xe, Float64[lo90[d] for d in est], Float64[hi90[d] for d in est];
        color = (:purple, 0.15))
    band!(ax, xe, Float64[lo50[d] for d in est], Float64[hi50[d] for d in est];
        color = (:purple, 0.28))
    band!(ax, xe, Float64[lo30[d] for d in est], Float64[hi30[d] for d in est];
        color = (:purple, 0.42))
    ## Horizontal grey dashed line at the no-growth threshold Rt = 1.
    hlines!(ax, [1.0]; color = (:grey, 0.8), linestyle = :dash, linewidth = 2)
    vlines!(ax, [Float64(epoch + breakpoint - 1)];
        color = :firebrick, linestyle = :dash, linewidth = 2)
    ## End of the intervention scale-up (breakpoint plus the ramp length).
    vlines!(ax, [Float64(epoch + breakpoint - 1 + ramp)];
        color = :firebrick, linestyle = :dot, linewidth = 2)
    vlines!(ax, [Float64(date2epochdays(Date(as_of_date)))];
        color = :grey, linestyle = :dash)
    ## Limit the x-axis to the estimated window so only the period being
    ## estimated is shown.
    lo = isempty(xe) ? floor(Int, minimum(x)) : floor(Int, minimum(xe))
    hi = ceil(Int, maximum(x))
    CairoMakie.xlims!(ax, lo, hi)
    CairoMakie.ylims!(ax, 0, ytop)
    ## Weekly date ticks across the estimated window.
    ax.xticks = collect(lo:7:hi)
    ax.xtickformat = vals -> [string(epochdays2date(round(Int, v)))
                              for v in vals]
    return fig
end

## Per-fit credible bands over the shared display window (`ds` onward), with
## the fit's own established mask applied through `reconstruct_rt`. Returns a
## NamedTuple of the 30/60/90% lower/upper quantiles and the established days
## `est` to draw, mirroring the band style used across the report.
function _rt_bands(chn; n, breakpoint, rt_start, rt_walk_start, week, ramp, ds)
    rt = reconstruct_rt(chn; n, breakpoint, rt_start, rt_walk_start, week, ramp)
    q(pr) = [_rt_quantile(rt, d, pr) for d in 1:n]
    med = q(0.5)
    est = findall(d -> d >= ds && !ismissing(med[d]), 1:n)
    return (; lo90 = q(0.05), hi90 = q(0.95), lo60 = q(0.20), hi60 = q(0.80),
        lo30 = q(0.35), hi30 = q(0.65), est)
end

## Draw nested 30/60/90% credible ribbons (no median line) for one fit's
## bands in `colour`, matching the ribbon style of [`plot_rt`](@ref).
function _draw_rt_bands!(ax, x, b, colour;
        alphas = (0.15, 0.28, 0.42))
    isempty(b.est) && return
    xe = x[b.est]
    band!(ax, xe, Float64[b.lo90[d] for d in b.est],
        Float64[b.hi90[d] for d in b.est]; color = (colour, alphas[1]))
    band!(ax, xe, Float64[b.lo60[d] for d in b.est],
        Float64[b.hi60[d] for d in b.est]; color = (colour, alphas[2]))
    band!(ax, xe, Float64[b.lo30[d] for d in b.est],
        Float64[b.hi30[d] for d in b.est]; color = (colour, alphas[3]))
    return
end

"""
Faceted implied reproduction number, one panel per single-stream fit, each
with the joint fit overlaid as the reference. Each fit's daily `Rt` is
reconstructed from its own sampled random walk exactly as in
[`plot_rt`](@ref) (see [`reconstruct_rt`](@ref)), so the figure shows what
reproduction number each data stream implies on its own against the
all-streams-together joint estimate.

Every panel draws 30/60/90% credible ribbons with no median line, matching
the band style used across the report: the joint in grey first, then the
stream in its colour on top (the panel title is in that colour). One panel
per stream rather than a single overlay, so the wide per-stream ribbons stay
legible. The y-axis is shared across panels and capped from the panels'
90% bands so a weakly-informed stream (the confirmed-only fit) does not
stretch the scale; ribbon above the cap is clipped.

Each `stream` is a `NamedTuple`
`(; label, chn, rt_start, rt_walk_start, colour)` and `joint` is
`(; label, chn, rt_start, rt_walk_start)`, where `chn` is that fit's chain
and `rt_start`/`rt_walk_start` are the renewal start and random-walk start
that fit used (the per-stream fits walk from day 1, the joint from the
breakpoint lead). `display_start` is the shared grid day the panels draw from
(the joint renewal start), so every stream reads over the same window.
`seeding` is the calendar date of grid day 1, so day `d` is
`seeding + (d - 1)`. The intervention breakpoint, the end of the scale-up
(`breakpoint + ramp`, dotted) and the cut-off are marked as in
[`plot_rt`](@ref).
"""
function plot_rt_streams(streams::AbstractVector;
        joint, n::Integer, breakpoint::Real,
        as_of_date::AbstractString, seeding::Date,
        display_start::Integer = 1, week::Integer = 7, ramp::Real = 14.0,
        ncols::Integer = 2, joint_colour = :grey25)
    epoch = date2epochdays(seeding)
    x = Float64[epoch + (d - 1) for d in 1:n]
    ds = clamp(display_start, 1, n)

    ## Reconstruct the joint once and each stream's bands; the joint is the
    ## shared reference drawn behind every stream.
    bj = _rt_bands(joint.chn; n, breakpoint, rt_start = joint.rt_start,
        rt_walk_start = joint.rt_walk_start, week, ramp, ds)
    sbands = Tuple{Any, NamedTuple}[]
    for s in streams
        b = _rt_bands(s.chn; n, breakpoint, rt_start = s.rt_start,
            rt_walk_start = s.rt_walk_start, week, ramp, ds)
        push!(sbands, (s, b))
    end

    ## Shared y-cap from the panels' typical 90% upper band (the median over
    ## days, robust to a single spiky day), so the weakly-informed streams do
    ## not stretch the axis while the Rt = 1 line stays visible.
    function panel_top(b)
        v = Float64[b.hi90[d] for d in b.est if !ismissing(b.hi90[d])]
        return isempty(v) ? 0.0 : quantile(v, 0.5)
    end
    tops = Float64[panel_top(bj); [panel_top(b) for (_, b) in sbands]]
    ytop = max(2.5, ceil(1.3 * maximum(tops) * 2) / 2)

    lo = floor(Int, x[ds])
    hi = ceil(Int, maximum(x))
    nrows = cld(length(sbands), ncols)
    fig = Figure(; size = (480 * ncols, 300 * nrows + 70))

    for (i, (s, b)) in enumerate(sbands)
        r = cld(i, ncols)
        c = i - (r - 1) * ncols
        ax = Axis(fig[r, c]; xlabel = "Date", ylabel = "Rt",
            title = s.label, titlecolor = s.colour,
            xticklabelrotation = pi / 6)
        ## Joint reference behind the stream, both as 30/60/90% ribbons.
        _draw_rt_bands!(ax, x, bj, joint_colour;
            alphas = (0.10, 0.16, 0.22))
        _draw_rt_bands!(ax, x, b, s.colour)
        hlines!(ax, [1.0]; color = (:grey, 0.8), linestyle = :dash,
            linewidth = 2)
        vlines!(ax, [Float64(epoch + breakpoint - 1)];
            color = :firebrick, linestyle = :dash, linewidth = 2)
        vlines!(ax, [Float64(epoch + breakpoint - 1 + ramp)];
            color = :firebrick, linestyle = :dot, linewidth = 2)
        vlines!(ax, [Float64(date2epochdays(Date(as_of_date)))];
            color = :grey, linestyle = :dash)
        CairoMakie.xlims!(ax, lo, hi)
        CairoMakie.ylims!(ax, 0, ytop)
        ax.xticks = collect(lo:14:hi)
        ax.xtickformat = vals -> [string(epochdays2date(round(Int, v)))
                                  for v in vals]
    end

    CairoMakie.Label(fig[nrows + 1, 1:ncols],
        "Bands are 30/60/90% credible intervals. Grey is the joint fit " *
        "(the same in every panel); the coloured band is the single-stream " *
        "fit named in the panel title.";
        fontsize = 12, padding = (0, 0, 0, 6))
    CairoMakie.Label(fig[0, 1:ncols],
        "Implied Rt by data stream, with the joint fit overlaid";
        fontsize = 16, font = :bold)
    return fig
end

"""
Two-panel density of the no-onward-transmission counterfactual from
[`predict_no_onward_deaths`](@ref). The left panel shows the *still
expected* deaths (`:delta_deaths`, the future deaths in cases already
infected by `T`, net of the `obs_deaths` already observed). The right
panel shows the *projected total* (`:total_projected = obs_deaths +
delta_deaths`) with a dashed black rule at `obs_deaths`. Both are
lower bounds: they assume every onward transmission stops at time `T`.
"""
function plot_no_onward_deaths(df::DataFrame; obs_deaths::Real)
    fig = Figure(; size = (980, 420))

    ax1 = Axis(fig[1, 1];
        xlabel = "Still expected deaths (beyond those already observed)",
        ylabel = "Posterior density",
        title = "Still expected (future)")
    density!(ax1, df.delta_deaths; color = (:firebrick, 0.5),
        strokecolor = :firebrick, strokewidth = 2)

    ax2 = Axis(fig[1, 2];
        xlabel = "Projected total deaths (no onward transmission)",
        ylabel = "Posterior density",
        title = "Projected total")
    density!(ax2, df.total_projected; color = (:firebrick, 0.5),
        strokecolor = :firebrick, strokewidth = 2)
    vlines!(ax2, [float(obs_deaths)];
        color = :black, linestyle = :dash, linewidth = 2)

    return fig
end

## Shared panel painter for the forecast figures: a histogram of the
## forecast new-count draws with its 90% predictive interval shaded.
function _forecast_count_panel!(fig, pos, v, title, colour)
    r, c = pos
    ## 98th-percentile x-axis cap so a heavy forecast tail does not squash the
    ## readable bulk of the histogram.
    upper = max(1.0, quantile(v, 0.98))
    lo = quantile(v, 0.05)
    hi = quantile(v, 0.95)
    ax = Axis(fig[r, c];
        xlabel = title, ylabel = "Predictive frequency",
        title = "One week ahead", limits = ((0, upper), nothing))
    vspan!(ax, lo, hi; color = (colour, 0.15))
    hist!(ax, v; bins = range(0, upper; length = 30), color = (colour, 0.7))
    return ax
end

"""
One-week-ahead forecast of the unobserved (latent) quantities from
[`forecast_reported`](@ref): new infections, new symptom onsets and new
deaths over the horizon, with the reproduction number left to keep evolving
over the horizon. Each count panel histograms the projected new latent
count with its 90% predictive interval shaded; the reproduction-number
panel shows the posterior of the end-of-horizon forecast `R_t` with the
no-growth line at one marked.
These are the latent counterparts of the observed-stream forecast in
[`plot_forecast`](@ref).
"""
function plot_forecast_latent(fc::DataFrame)
    count_cols = [
        (:infections_new, "New infections (DRC)", :steelblue),
        (:onsets_new, "New symptom onsets (DRC)", :seagreen),
        (:deaths_latent_new, "New deaths (DRC)", :firebrick)
    ]
    npanels = length(count_cols) + 1
    ncols = 2
    nrows = cld(npanels, ncols)
    fig = Figure(; size = (400 * ncols, 360 * nrows))
    for (i, (col, title, colour)) in enumerate(count_cols)
        pos = (cld(i, ncols), mod1(i, ncols))
        _forecast_count_panel!(fig, pos, fc[!, col], title, colour)
    end
    ## Forecast reproduction number panel (a value, not a count).
    i = npanels
    r, c = cld(i, ncols), mod1(i, ncols)
    rt = fc[!, :rt_forecast]
    ax = Axis(fig[r, c];
        xlabel = "Forecast reproduction number (DRC)",
        ylabel = "Posterior density", title = "One week ahead")
    density!(ax, rt; color = (:purple, 0.5),
        strokecolor = :purple, strokewidth = 2)
    vlines!(ax, [1.0]; color = :black, linestyle = :dash, linewidth = 2)
    return fig
end

"""
One-week-ahead forecast of the observed confirmed quantities from
[`forecast_reported`](@ref): new laboratory-confirmed cases and new
confirmed deaths over the horizon. Each panel histograms the projected new
count with its 90% predictive interval shaded. The panels are drawn only
when the forecast carries the laboratory streams. The suspected reported
streams are no longer reported, so they are not forecast here. The latent
counterparts are shown by [`plot_forecast_latent`](@ref).
"""
function plot_forecast(fc::DataFrame)
    count_cols = Tuple{Symbol, String, Symbol}[]
    :confirmed_new in propertynames(fc) && push!(count_cols,
        (:confirmed_new, "New confirmed cases (DRC)", :goldenrod))
    :confirmed_deaths_new in propertynames(fc) && push!(count_cols,
        (:confirmed_deaths_new, "New confirmed deaths (DRC)", :darkorange3))
    npanels = length(count_cols)
    ncols = min(npanels, 2)
    nrows = cld(npanels, ncols)
    fig = Figure(; size = (400 * ncols, 360 * nrows))
    for (i, (col, title, colour)) in enumerate(count_cols)
        pos = (cld(i, ncols), mod1(i, ncols))
        _forecast_count_panel!(fig, pos, fc[!, col], title, colour)
    end
    return fig
end

"""
One-week-ahead forecast of the daily isolation/treatment flows from
[`forecast_reported`](@ref): the projected new admissions, in-care deaths and
rule-outs a day at the horizon. Each panel histograms the projected daily
count with its 90% predictive interval shaded, drawn only when the forecast
carries the flow streams (`admissions_fc`, `incare_deaths_fc`, `ruleouts_fc`).
These are the daily-flow counterparts of the bed-stock forecast in
[`plot_forecast_beds`](@ref).
"""
function plot_forecast_flows(fc::DataFrame)
    count_cols = Tuple{Symbol, String, Symbol}[]
    :admissions_fc in propertynames(fc) && push!(count_cols,
        (:admissions_fc, "New isolation admissions (DRC)", :steelblue))
    :incare_deaths_fc in propertynames(fc) && push!(count_cols,
        (:incare_deaths_fc, "New in-care deaths (DRC)", :firebrick))
    :ruleouts_fc in propertynames(fc) && push!(count_cols,
        (:ruleouts_fc, "New rule-outs (DRC)", :seagreen))
    npanels = length(count_cols)
    npanels == 0 && return Figure()
    ncols = min(npanels, 2)
    nrows = cld(npanels, ncols)
    fig = Figure(; size = (400 * ncols, 360 * nrows))
    for (i, (col, title, colour)) in enumerate(count_cols)
        pos = (cld(i, ncols), mod1(i, ncols))
        _forecast_count_panel!(fig, pos, fc[!, col], title, colour)
    end
    return fig
end

"""
One-week-ahead isolation/treatment-bed forecast from
[`forecast_reported`](@ref): the projected bed DEMAND (the need a week ahead,
under unconstrained supply) against the supply-limited occupancy (the beds
actually filled), and the shortfall between them. The left panel overlays the
two predictive distributions, so the gap between the need and the
supply-limited occupancy is the unmet demand; the right panel histograms the
shortfall directly. Drawn only when the forecast carries the bed streams
(`bed_demand` and `isolation_level`).

Because the model carries a SINGLE national bed capacity it cannot represent
LOCAL saturation — on 13 June Ituri was at 93.9% occupancy while Sud-Kivu was
at 21.9% — so the national shortfall understates the local unmet need (beds
free in one province cannot serve patients who need them in another).
"""
function plot_forecast_beds(fc::DataFrame)
    (:bed_demand in propertynames(fc) &&
     :isolation_level in propertynames(fc)) || return Figure()
    demand = float.(fc[!, :bed_demand])
    occ = float.(fc[!, :isolation_level])
    shortfall = max.(demand .- occ, 0.0)
    fig = Figure(; size = (800, 360))
    ## Cap the x-axis at the 98th percentile of demand: the unconstrained
    ## bed-demand projection is heavy-tailed (it grows with the reproduction
    ## number over the horizon), so its long upper tail otherwise squashes the
    ## readable bulk of both densities. Occupancy is capped at capacity, so it
    ## sits below this bound.
    upper = max(1.0, quantile(demand, 0.98))
    ax1 = Axis(fig[1, 1];
        xlabel = "Isolation beds a week ahead (DRC)",
        ylabel = "Predictive density", title = "Need vs supply-limited use",
        limits = ((0, upper), nothing))
    density!(ax1, demand; color = (:darkorange, 0.35),
        strokecolor = :darkorange, strokewidth = 2, label = "Demand (need)")
    density!(ax1, occ; color = (:steelblue, 0.35),
        strokecolor = :steelblue, strokewidth = 2,
        label = "Occupancy (supply-limited)")
    CairoMakie.axislegend(ax1; position = :rt, framevisible = false)
    _forecast_count_panel!(fig, (1, 2), shortfall, "Bed shortfall (DRC)",
        :firebrick)
    return fig
end

"""
Validate a [`forecast_reported`](@ref) bed projection against the beds
actually occupied a week later. Histograms the projected supply-limited
isolation-bed occupancy with the 90% predictive interval shaded and the
`isolation` count observed at the target date drawn as a dashed black rule,
so last week's bed forecast is scored against what the beds held. Drawn only
when the forecast carries `isolation_level`.

At a one-week-back freeze the bed capacity has no implied-capacity anchor
(the reported occupancy rate starts only on 9 June), so the projected
occupancy rides the capacity random walk back to the freeze date and the
interval is wide — the bed forecast is the weakest of the validated streams.
"""
function plot_forecast_beds_vs_truth(fc::DataFrame;
        isolation::Union{Real, Missing})
    (isolation !== missing && :isolation_level in propertynames(fc)) ||
        return Figure()
    v = float.(fc[!, :isolation_level])
    lo = quantile(v, 0.05)
    hi = quantile(v, 0.95)
    upper = max(1.0, quantile(v, 0.995), float(isolation) * 1.05)
    fig = Figure(; size = (440, 360))
    ax = Axis(fig[1, 1];
        xlabel = "Isolation beds occupied at the target date (DRC)",
        ylabel = "Predictive frequency", title = "Forecast vs observed",
        limits = ((0, upper), nothing))
    vspan!(ax, lo, hi; color = (:steelblue, 0.15))
    hist!(ax, v; bins = range(0, upper; length = 30),
        color = (:steelblue, 0.7))
    vlines!(ax, [float(isolation)]; color = :black, linestyle = :dash,
        linewidth = 2)
    return fig
end

"""
Confirmed-stream validation figure for a [`forecast_reported`](@ref)
projection, laid out as a two-row grid. The top row shows the cumulative
forecast distribution per confirmed stream (DRC confirmed cases and
confirmed deaths); the bottom row shows the new counts forecast over the
horizon. Each panel is a histogram with the 90% predictive interval shaded
and the later-observed count drawn as a dashed black rule, so the forecast
distribution is scored against the count that was actually observed.
`confirmed` and `confirmed_deaths` are the observed cumulative counts;
`baseline_*` are the counts at the forecast origin, so the observed new
count is the cumulative truth minus the baseline. The suspected reported
streams are no longer reported, so they are not scored here. The latent
counterparts are scored distribution-versus-distribution by
[`plot_forecast_vs_truth_latent`](@ref).
"""
function plot_forecast_vs_truth(fc::DataFrame;
        confirmed::Real, confirmed_deaths::Real,
        baseline_confirmed::Real = 0,
        baseline_confirmed_deaths::Real = 0)
    streams = Vector{Tuple{Symbol, Symbol, String, Symbol, Float64, Float64}}()
    :confirmed_cum in propertynames(fc) &&
        push!(streams,
            (:confirmed_cum, :confirmed_new, "confirmed cases (DRC)",
                :goldenrod, float(confirmed),
                float(confirmed) - float(baseline_confirmed)))
    :confirmed_deaths_cum in propertynames(fc) &&
        push!(streams,
            (:confirmed_deaths_cum, :confirmed_deaths_new,
                "confirmed deaths (DRC)", :darkorange3,
                float(confirmed_deaths),
                float(confirmed_deaths) - float(baseline_confirmed_deaths)))
    ncols = length(streams)
    fig = Figure(; size = (370 * ncols, 680))
    function panel!(row, col, v, obs, title, colour)
        lo = quantile(v, 0.05)
        hi = quantile(v, 0.95)
        upper = max(1.0, quantile(v, 0.995), obs * 1.05)
        ax = Axis(fig[row, col];
            xlabel = title, ylabel = "Predictive frequency",
            limits = ((0, upper), nothing))
        vspan!(ax, lo, hi; color = (colour, 0.15))
        hist!(ax, v; bins = range(0, upper; length = 30),
            color = (colour, 0.7))
        vlines!(ax, [obs]; color = :black, linestyle = :dash, linewidth = 2)
    end
    for (j, (ccol, ncol, name, colour, obs_cum, obs_new)) in enumerate(streams)
        panel!(1, j, fc[!, ccol], obs_cum, "Cumulative $name", colour)
        panel!(2, j, fc[!, ncol], max(obs_new, 0.0), "New $name", colour)
    end
    return fig
end

"""
Latent-quantity validation figure: for each unobserved quantity (new
infections, new symptom onsets, new deaths over the past week) overlay the
distribution the frozen (last-week) fit forecast against the distribution
the current fit now estimates for the same window. Both are latent, so the
comparison is density versus density rather than density versus a single
observed count. `fc` is the frozen forecast from [`forecast_reported`](@ref)
(its `*_new` latent columns); `now` is a `NamedTuple` carrying the current
fit's draws of the same new-count quantities,
`(; infections_new, onsets_new, deaths_latent_new)`.
"""
function plot_forecast_vs_truth_latent(fc::DataFrame; now::NamedTuple)
    panels = [
        (:infections_new, "New infections (DRC)", :steelblue),
        (:onsets_new, "New symptom onsets (DRC)", :seagreen),
        (:deaths_latent_new, "New deaths (DRC)", :firebrick)
    ]
    ncols = length(panels)
    fig = Figure(; size = (370 * ncols, 380))
    local frozen_h, now_h
    for (j, (col, title, colour)) in enumerate(panels)
        vf = float.(fc[!, col])
        vn = float.(getproperty(now, col))
        upper = max(1.0, quantile(vf, 0.99), quantile(vn, 0.99))
        ax = Axis(fig[1, j];
            xlabel = title, ylabel = "Posterior density",
            limits = ((0, upper), nothing))
        frozen_h = density!(ax, vf; color = (colour, 0.25),
            strokecolor = colour, strokewidth = 2)
        now_h = density!(ax, vn; color = (:grey, 0.0),
            strokecolor = :black, strokewidth = 2, linestyle = :dash)
    end
    CairoMakie.Legend(fig[1, ncols + 1], [frozen_h, now_h],
        ["Forecast last week", "Estimated now"])
    return fig
end

"""
Per-vintage conditional one-step-ahead posterior-predictive for the DRC
streams. For each `panel` the predicted cumulative count at vintage `v`
conditions on the *observed* cumulative at the previous vintage and adds
only the posterior-predictive between-vintage increment,
``\\hat{y}_v = y_{v-1} + \\Delta_v`` with ``y_0 = 0``. The increment
``\\Delta_v`` is the per-bin replicate draw the model already samples in
predictive mode, so each step carries the full posterior uncertainty of
the *new* increment while grounding on what was actually observed at the
preceding sitrep. This is the "filtered" one-step-ahead predictive: it
isolates the model's per-step increment prediction and does not let
errors compound across the series, unlike a running sum of the modelled
increments. The result is summarised by vintage as shaded 30/60/90%
predictive ribbons; the observed cumulative counts are overlaid as
points. Each `panel` is a `NamedTuple`
`(; title, dates, replicates, observed)`, where `replicates` is a vector
of per-draw increment vectors (one entry per vintage, oldest first) and
`observed` the matching observed cumulative counts used as the
conditioning baselines. `colour` is optional per panel. A panel may set
`cumulative = false` to plot standalone per-day counts instead of a
cumulative series (e.g. the post-26 May daily new-suspect inflow, or the 24h
analysed volume): there is no previous-vintage baseline, so each replicate is
plotted as its own daily count against the observed daily count, and the
y-axis reads "Daily count".

`max_date` (an ISO date string or `Date`) truncates every panel to the
vintages on or before that date, so streams that keep reporting past the
others (the laboratory streams run to the cut-off while the suspected
streams freeze earlier) are cut back to the shared last date. Without
this the confirmed panel runs further along the date axis than the
suspected panel and reads as though it overtakes it, when the two are
simply shown to different end dates.
"""
function plot_vintage_conditional_ppc(
        panels::AbstractVector; xlabel = "Sitrep date",
        max_date::Union{Nothing, Date, AbstractString} = nothing)
    cap = isnothing(max_date) ? nothing :
          (max_date isa Date ? max_date : Date(String(max_date)))
    npanels = length(panels)
    ## Cap the grid at four columns so a large stream set lays out over
    ## several rows rather than one very wide strip that the page downscales
    ## into tiny panels; with the treatment-centre flows there are ~14 streams,
    ## giving a readable 4-column, ~4-row grid.
    ncols = min(npanels, 4)
    nrows = cld(npanels, ncols)
    fig = Figure(; size = (460 * ncols, 420 * nrows))
    for (j, p) in enumerate(panels)
        row, col = cld(j, ncols), mod1(j, ncols)
        ## Drop vintages past the shared cap so every panel ends on the
        ## same date; the replicates and observed counts are truncated to
        ## match, keeping the conditional baselines aligned.
        keep = isnothing(cap) ? eachindex(p.dates) :
               [i for i in eachindex(p.dates) if Date(p.dates[i]) <= cap]
        dates = p.dates[keep]
        observed = p.observed[keep]
        replicates = [collect(r)[keep] for r in vec(collect(p.replicates))]
        n = length(dates)
        colour = get(p, :colour, :steelblue)
        ## A panel with `cumulative = false` (the post-26 May daily
        ## new-suspect inflow, or the 24h analysed volume) carries standalone
        ## per-day counts rather than a cumulative series, so there is no
        ## previous-vintage baseline to condition on: each replicate is its
        ## own per-day count directly.
        cumulative = get(p, :cumulative, true)
        ## Observed cumulative at the previous vintage is the conditioning
        ## baseline for each step (`y_0 = 0`); `obs_prev[v]` is `y_{v-1}`.
        obs_cum = float.(observed)
        obs_prev = cumulative ?
                   [v == 1 ? 0.0 : obs_cum[v - 1] for v in 1:n] : zeros(n)
        ## `replicates` is already flattened to one vector of per-draw
        ## increment vectors and truncated to the kept vintages above.
        ## Each draw's conditional cumulative at vintage `v` is the
        ## observed previous cumulative plus the drawn increment `Δ_v` (the
        ## baseline is zero for a non-cumulative panel).
        cond = [obs_prev .+ r for r in replicates]
        q(i, pr) = quantile([c[i] for c in cond], pr)
        lo90 = [q(i, 0.05) for i in 1:n]
        hi90 = [q(i, 0.95) for i in 1:n]
        lo60 = [q(i, 0.20) for i in 1:n]
        hi60 = [q(i, 0.80) for i in 1:n]
        lo30 = [q(i, 0.35) for i in 1:n]
        hi30 = [q(i, 0.65) for i in 1:n]
        x = collect(1:n)
        ## Truncate the y-axis to a robust ceiling driven by the observed
        ## counts and the 60% band, so a heavy upper tail in any one stream
        ## (the 90% band can run far above the data) does not flatten the
        ## visible detail; the band clips at the axis limit.
        yupper = 1.6 * max(isempty(obs_cum) ? 1.0 : maximum(obs_cum),
            isempty(hi60) ? 1.0 : maximum(hi60), 1.0)
        ax = Axis(fig[row, col]; title = p.title, xlabel = xlabel,
            ylabel = cumulative ? (col == 1 ? "Cumulative count" : "") :
                     "Daily count",
            xticks = (x, string.(dates)),
            xticklabelrotation = pi / 4, xticklabelsize = 9,
            limits = (nothing, (0, yupper)))
        ## 30/60/90% credible ribbons over the situation-report dates.
        band!(ax, x, lo90, hi90; color = (colour, 0.15))
        band!(ax, x, lo60, hi60; color = (colour, 0.28))
        band!(ax, x, lo30, hi30; color = (colour, 0.42))
        scatter!(ax, x, float.(observed); color = :black, markersize = 9)
    end
    return fig
end

"""
Per-vintage incidence posterior-predictive check: the same panels as
[`plot_vintage_conditional_ppc`](@ref) but plotting the count between
consecutive vintages rather than the running cumulative, so trends (a rise
or a slowdown) read directly off the height of each bar-like step instead of
the slope of a near-straight cumulative line. For a cumulative panel the
observed incidence is the between-vintage increment (the first vintage is
its own baseline); for a non-cumulative panel (a standalone per-day count
such as the 24h analysed volume or the daily new-suspect inflow) it is the
count itself. The replicates are already per-vintage increments, so they are
the modelled incidence directly and are summarised as 30/60/90% credible
ribbons with the observed incidence overlaid. `panels` and `max_date` match
[`plot_vintage_conditional_ppc`](@ref).
"""
function plot_vintage_incidence_ppc(
        panels::AbstractVector; xlabel = "Sitrep date",
        max_date::Union{Nothing, Date, AbstractString} = nothing)
    cap = isnothing(max_date) ? nothing :
          (max_date isa Date ? max_date : Date(String(max_date)))
    npanels = length(panels)
    ## Cap the grid at four columns so a large stream set lays out over
    ## several rows rather than one very wide strip that the page downscales
    ## into tiny panels; with the treatment-centre flows there are ~14 streams,
    ## giving a readable 4-column, ~4-row grid.
    ncols = min(npanels, 4)
    nrows = cld(npanels, ncols)
    fig = Figure(; size = (460 * ncols, 420 * nrows))
    for (j, p) in enumerate(panels)
        row, col = cld(j, ncols), mod1(j, ncols)
        keep = isnothing(cap) ? eachindex(p.dates) :
               [i for i in eachindex(p.dates) if Date(p.dates[i]) <= cap]
        dates = p.dates[keep]
        observed = p.observed[keep]
        replicates = [collect(r)[keep] for r in vec(collect(p.replicates))]
        n = length(dates)
        colour = get(p, :colour, :steelblue)
        cumulative = get(p, :cumulative, true)
        ## Observed per-vintage incidence: the increment for a cumulative
        ## series (the first vintage is its own baseline at zero), or the
        ## standalone count for a non-cumulative panel.
        obs_cum = float.(observed)
        obs_inc = cumulative ?
                  [v == 1 ? obs_cum[v] : obs_cum[v] - obs_cum[v - 1]
                   for v in 1:n] : obs_cum
        ## The replicates are already per-vintage increments (per-day counts
        ## for a non-cumulative panel), so they are the modelled incidence.
        q(i, pr) = quantile([r[i] for r in replicates], pr)
        lo90 = [q(i, 0.05) for i in 1:n]
        hi90 = [q(i, 0.95) for i in 1:n]
        lo60 = [q(i, 0.20) for i in 1:n]
        hi60 = [q(i, 0.80) for i in 1:n]
        lo30 = [q(i, 0.35) for i in 1:n]
        hi30 = [q(i, 0.65) for i in 1:n]
        x = collect(1:n)
        yupper = 1.6 * max(isempty(obs_inc) ? 1.0 : maximum(obs_inc),
            isempty(hi60) ? 1.0 : maximum(hi60), 1.0)
        ax = Axis(fig[row, col]; title = p.title, xlabel = xlabel,
            ylabel = col == 1 ? "New per vintage" : "",
            xticks = (x, string.(dates)),
            xticklabelrotation = pi / 4, xticklabelsize = 9,
            limits = (nothing, (0, yupper)))
        band!(ax, x, lo90, hi90; color = (colour, 0.15))
        band!(ax, x, lo60, hi60; color = (colour, 0.28))
        band!(ax, x, lo30, hi30; color = (colour, 0.42))
        scatter!(ax, x, float.(obs_inc); color = :black, markersize = 9)
    end
    return fig
end

"""
Per-stream calibration of the one-step-ahead conditional posterior
predictive, plotting the table from [`stream_calibration`](@ref). Pass the
table that function returns (its prettified columns `Stream`, `50% coverage`,
`90% coverage`, `Bias`). The figure has two panels sharing a categorical
y-axis of streams: the left panel marks each stream's empirical 50% and 90%
coverage with vertical dashed reference lines at the nominal 0.5 and 0.9, so a
well-calibrated stream sits on its line and a marker to the left of its line
flags under-coverage; the right panel marks the mean forecast `Bias` (negative
= the stream is under-predicted, positive = over-predicted) with a dashed line
at zero. Streams are read off the shared row labels, so all ~14 stay legible
without splitting the figure into columns.
"""
function plot_stream_calibration(tbl::DataFrame)
    ## Read the prettified columns the table carries; oldest-first order is
    ## kept but reversed for the y-axis so the first stream reads at the top.
    streams = string.(tbl[!, "Stream"])
    cov50 = float.(tbl[!, "50% coverage"])
    cov90 = float.(tbl[!, "90% coverage"])
    bias = float.(tbl[!, "Bias"])
    n = length(streams)
    ## Categorical y positions, top-to-bottom in table order.
    y = collect(n:-1:1)
    height = max(360, 60 + 26 * n)
    fig = Figure(; size = (980, height))

    ax1 = Axis(fig[1, 1];
        xlabel = "Empirical coverage", title = "Interval coverage",
        yticks = (y, streams), limits = ((0, 1), nothing))
    ## Nominal reference lines: a marker on its line is well calibrated.
    vlines!(ax1, [0.5]; color = (:steelblue, 0.6), linestyle = :dash,
        linewidth = 2)
    vlines!(ax1, [0.9]; color = (:seagreen, 0.6), linestyle = :dash,
        linewidth = 2)
    h50 = scatter!(ax1, cov50, y; color = :steelblue, markersize = 11)
    h90 = scatter!(ax1, cov90, y; color = :seagreen, markersize = 11,
        marker = :diamond)
    CairoMakie.axislegend(ax1, [h50, h90], ["50% interval", "90% interval"];
        position = :lt, framevisible = false)

    ## Bias panel: zero is unbiased; sign flags over/under-prediction.
    bmax = max(1.0, maximum(abs.(bias)) * 1.1)
    ax2 = Axis(fig[1, 2];
        xlabel = "Mean forecast bias", title = "Forecast bias",
        yticks = (y, fill("", n)), limits = ((-bmax, bmax), nothing))
    vlines!(ax2, [0.0]; color = :black, linestyle = :dash, linewidth = 2)
    scatter!(ax2, bias, y; color = :firebrick, markersize = 11)
    return fig
end
