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
symptom onsets and cumulative deaths, all modelled latent quantities. The
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
Estimate-evolution plot: how the outbreak-size estimate moves as the
data cut-off advances, drawn against the calendar date.

`released` is a vector of `(cutoff_date, median, lo30, hi30, lo60, hi60,
lo90, hi90)` tuples, one per published project release, drawn in blue with
nested 30/60/90% credible-interval ribbons along the date axis. `renewal`
is a vector of the same tuple shape, the current renewal model re-fit
frozen at each release date, drawn in red: the current method evaluated at
each past cut-off, so the series rises as the outbreak grows. Both are
summarised by ribbons only, no central line or marker.

`trajectory` is the current-data, current-model cumulative-infection
trajectory over the day grid, a `(dates, lo30, hi30, lo60, hi60, lo90,
hi90)` tuple where `dates` is the calendar date of each grid day and the
remaining entries are the per-day credible bounds. It is drawn in a third
colour as a time-varying ribbon, so the latest estimate rises across the
period on the same calendar axis as the release points and lines up with
them.

Release dates are marked with dotted vertical rules so each release reads
against the rising trajectory.

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
        "Current model, current data")
    ## Calendar dates → numeric day-offsets so the x-axis is to scale,
    ## then relabel the ticks with the dates. The released points, the
    ## frozen renewal points and the current-model trajectory all share
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

    ## One per-vintage series: nested 30/60/90% credible-interval ribbons
    ## along the date axis, with no central line or marker.
    function _series!(dates, tuples, colour)
        xs = _x.(dates)
        ord = sortperm(xs)
        xs = xs[ord]
        band!(ax, xs, [float(t[7]) for t in tuples][ord],
            [float(t[8]) for t in tuples][ord]; color = (colour, 0.12))
        band!(ax, xs, [float(t[5]) for t in tuples][ord],
            [float(t[6]) for t in tuples][ord]; color = (colour, 0.20))
        return band!(ax, xs, [float(t[3]) for t in tuples][ord],
            [float(t[4]) for t in tuples][ord]; color = (colour, 0.30))
    end

    handles = Any[]
    labels = String[]
    ## Current-data estimate as the cumulative-infection trajectory over the
    ## grid, a time-varying ribbon on the same calendar axis as the points.
    ## `trajectory` is `(dates, lo30, hi30, lo60, hi60, lo90, hi90)`.
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
    if !isempty(renewal)
        renewal_band = _series!(ndates, renewal, :firebrick)
        push!(handles, renewal_band)
        push!(labels, renewal_label * " (30/60/90% bands)")
    end
    released_band = _series!(rdates, released, :steelblue)
    push!(handles, released_band)
    push!(labels, released_label * " (30/60/90% bands)")

    ## Dotted vertical rule at each release date.
    isempty(rdates) || vlines!(ax, _x.(rdates);
        color = (:grey, 0.55), linestyle = :dot, linewidth = 1)

    CairoMakie.axislegend(ax, handles, labels; position = :lt,
        framevisible = true)
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
        rt_start::Integer = 1, week::Integer = 7, ramp::Real = 14.0,
        n_traj::Integer = 100)
    log_R0 = _draws(chn, Symbol("rt_state.log_R0"))
    sigma = _draws(chn, Symbol("rt_state.sigma_rw"))
    effect = _draws(chn, Symbol("rt_state.intervention_effect"))
    T_draws = _draws(chn, :T)
    ## `rt_state.z` is vector-valued: one standard-normal innovation vector
    ## per draw. Pull each draw's full vector from the chain slice.
    zmat = chn[Symbol("rt_state.z")]
    zrows = [collect(z) for z in vec(collect(zmat))]

    days = knot_days(n; week, start = rt_start)
    nb = length(days)
    ## The innovation vector length is fixed by the model's own `rt_start`
    ## (the renewal start, `n - tmrca_days + RENEWAL_START_LEAD`). If the
    ## caller passes a different `rt_start` the knot grid here will not match,
    ## so fail with a clear message rather than a downstream bounds error.
    if !isempty(zrows) && length(zrows[1]) != nb - 1
        error("plot_rt: rt_start = $rt_start gives $(nb - 1) random-walk " *
              "steps but the chain has $(length(zrows[1])); pass the same " *
              "renewal start the model used " *
              "(n - tmrca_days + RENEWAL_START_LEAD).")
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
        ## the model and are not plotted (masked below from `rt_start`
        ## onward).
        log_Rt = walk .+ effect[i] .* ramp_shape
        ## Only the established R_t, from the genetic bound onward; the
        ## pre-bound seeding window is not plotted.
        start = clamp(rt_start, 1, n)
        for d in start:n
            rt[i, d] = exp(log_Rt[d])
        end
    end

    ## Median and ribbons over established draws only (skip masked days).
    function q(d, pr)
        col = collect(skipmissing(@view rt[:, d]))
        isempty(col) ? missing : quantile(col, pr)
    end
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
    ## Cap the y-axis a little above the upper 90% credible band so a handful
    ## of high sampled trajectories do not stretch the scale; the thin
    ## trajectories above the cap are simply clipped.
    ytop = isempty(est) ? 6.0 :
           1.2 * maximum(Float64[hi90[d] for d in est])
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
    ## Cap the y-axis at roughly twice the 90% upper bound (rounded up to a
    ## tidy step), so stray trajectories do not stretch the axis into
    ## whitespace while the Rt = 1 line stays visible.
    hi90_est = [hi90[d] for d in est if !ismissing(hi90[d])]
    ytop = isempty(hi90_est) ? 4.0 :
           max(1.2, ceil(2 * maximum(hi90_est) * 2) / 2)
    CairoMakie.ylims!(ax, 0, ytop)
    ## Weekly date ticks across the estimated window.
    ax.xticks = collect(lo:7:hi)
    ax.xtickformat = vals -> [string(epochdays2date(round(Int, v)))
                              for v in vals]
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
    upper = max(1.0, quantile(v, 0.995))
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
One-week-ahead forecast of the observed (reported) quantities from
[`forecast_reported`](@ref): new reported cases, new confirmed cases and
new confirmed deaths over the horizon. Each panel histograms the projected
new count with its 90% predictive interval shaded. The confirmed panels are
drawn only when the forecast carries the laboratory streams. The latent
counterparts are shown by [`plot_forecast_latent`](@ref).
"""
function plot_forecast(fc::DataFrame)
    count_cols = Tuple{Symbol, String, Symbol}[
    (
        :cases_new, "New reported cases (DRC)", :seagreen)]
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
Observed-stream validation figure for a [`forecast_reported`](@ref)
projection, laid out as a two-row grid. The top row shows the cumulative
forecast distribution per stream (DRC reported cases, DRC deaths, and the
laboratory streams when present); the bottom row shows the new counts
forecast over the horizon. Each panel is a histogram with the 90%
predictive interval shaded and the later-observed count drawn as a dashed
black rule, so the forecast distribution is scored against the count that
was actually observed. `cases`, `deaths` and `exports` are the observed
cumulative counts; `baseline_*` are the counts at the forecast origin, so
the observed new count is the cumulative truth minus the baseline. The
latent counterparts are scored distribution-versus-distribution by
[`plot_forecast_vs_truth_latent`](@ref).
"""
function plot_forecast_vs_truth(fc::DataFrame;
        cases::Real, deaths::Real,
        exports::Union{Real, Missing} = missing,
        confirmed::Union{Real, Missing} = missing,
        tests::Union{Real, Missing} = missing,
        baseline_cases::Real = 0, baseline_deaths::Real = 0,
        baseline_exports::Real = 0,
        baseline_confirmed::Real = 0,
        baseline_tests::Real = 0)
    streams = Vector{Tuple{Symbol, Symbol, String, Symbol, Float64, Float64}}([
        (:cases_cum, :cases_new, "reported cases (DRC)", :steelblue,
            float(cases), float(cases) - float(baseline_cases)),
        (:deaths_cum, :deaths_new, "deaths (DRC)", :firebrick,
            float(deaths), float(deaths) - float(baseline_deaths))
    ])
    exports !== missing && :exports_cum in propertynames(fc) &&
        push!(streams,
            (:exports_cum, :exports_new, "exports (Uganda)", :seagreen,
                float(exports), float(exports) - float(baseline_exports)))
    tests !== missing && :tests_cum in propertynames(fc) &&
        push!(streams,
            (:tests_cum, :tests_new, "tests analysed (DRC)", :teal,
                float(tests), float(tests) - float(baseline_tests)))
    confirmed !== missing && :confirmed_cum in propertynames(fc) &&
        push!(streams,
            (:confirmed_cum, :confirmed_new, "confirmed cases (DRC)",
                :goldenrod, float(confirmed),
                float(confirmed) - float(baseline_confirmed)))
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
conditioning baselines. `colour` is optional per panel.

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
    nrows = npanels > 1 ? 2 : 1
    ncols = cld(npanels, nrows)
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
        ## Observed cumulative at the previous vintage is the conditioning
        ## baseline for each step (`y_0 = 0`); `obs_prev[v]` is `y_{v-1}`.
        obs_cum = float.(observed)
        obs_prev = [v == 1 ? 0.0 : obs_cum[v - 1] for v in 1:n]
        ## `replicates` is already flattened to one vector of per-draw
        ## increment vectors and truncated to the kept vintages above.
        ## Each draw's conditional cumulative at vintage `v` is the
        ## observed previous cumulative plus the drawn increment `Δ_v`.
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
            ylabel = col == 1 ? "Cumulative count" : "",
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
