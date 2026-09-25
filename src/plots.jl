# All package figures: posterior densities of `C_T`, posterior- and
# prior-predictive panel grids, pair plots, point-and-interval
# comparison, CFR prior, start-date and no-onward-transmission
# densities, and the one-week-ahead forecast figures.

"""
Kernel density for a quantity that cannot fall below `lower`, with the axis
clipped to the side of the bound the quantity can reach. A Gaussian KDE puts
mass past the smallest draw, so a bounded quantity otherwise shows a tail on
the impossible side.

The axis is clipped rather than the estimator narrowed with `boundary`. That
tabulation drops draws landing on the first grid point while still
normalising by the full sample size, so it discards the mass piled against
the bound.
"""
function _bounded_density!(ax, x; lower::Real, kwargs...)
    h = density!(ax, x; kwargs...)
    ## Only the impossible side is clipped. Pinning the upper limit too would
    ## cut the curve off at the largest draw.
    CairoMakie.xlims!(ax, float(lower), nothing)
    return h
end

"""
Overlaid posterior densities of `C_T` from one or more fits, built through
AlgebraOfGraphics. The published scenario point estimates are drawn as faint
dashed `vlines` on top.
"""
function plot_cumulative_cases(
        streams::Pair{String, <:AbstractVector}...;
        scenarios = REPORT_SCENARIOS,
        xmax::Union{Nothing, Real} = nothing,
        xlabel::AbstractString = "Cumulative infections",
        title::AbstractString = "Outbreak size estimated by each data stream"
    )
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
        AoG.mapping(
        :C_T => xlabel,
        color = :stream => "Data stream"
    ) *
        AoG.AlgebraOfGraphics.density() *
        AoG.subvisual(:line, linewidth = 2)
    fg = AoG.draw(
        spec;
        axis = (;
            ylabel = "Posterior density",
            title = title,
            limits = ((0, upper), nothing),
        ),
        figure = (; size = (760, 420))
    )

    scenario_xs = Float64[val for (_, val) in scenarios if val < upper]
    isempty(scenario_xs) || vlines!(
        fg.figure.content[1], scenario_xs;
        color = (:grey, 0.4), linestyle = :dash
    )
    return fg
end

"""
Headline 3x2 cumulative figure. Rows are cumulative infections, cumulative
symptom onsets and cumulative deaths, all modelled BVD-only latent renewal
quantities. The deaths row excludes the non-BVD background, so it stays as
smooth as the infection and onset rows. The left column is the modelled
expected cumulative trajectory over the grid as 30%, 60% and 90% ribbons with
no median line. The right column is the posterior density of the cut-off
cumulative.

The chain must carry the vector deterministics `cumulative_infections`,
`cumulative_onsets` and `cumulative_expected_deaths`, one per draw.
`seeding` is the calendar date of grid day 1, so day `d` is
`seeding + (d - 1)`. No observed data is overlaid. Each row is a latent
quantity upstream of ascertainment, confirmation and reporting delays, so the
observed counts are not on the same scale.
"""
function plot_cumulative_trajectories(
        chn;
        n::Integer, seeding::Date
    )
    epoch = date2epochdays(seeding)
    x = Float64[epoch + (d - 1) for d in 1:n]

    ## Each trajectory deterministic is an iter×chain matrix of per-draw
    ## vectors, so flatten to one vector of per-draw trajectories.
    function _trajectories(key)
        mat = chn[key]
        return [collect(v) for v in vec(collect(mat))]
    end
    function _ribbon(trajs)
        q(d, pr) = quantile(Float64[t[d] for t in trajs], pr)
        lo90 = [q(d, 0.05) for d in 1:n]
        hi90 = [q(d, 0.95) for d in 1:n]
        lo60 = [q(d, 0.2) for d in 1:n]
        hi60 = [q(d, 0.8) for d in 1:n]
        lo30 = [q(d, 0.35) for d in 1:n]
        hi30 = [q(d, 0.65) for d in 1:n]
        return lo90, hi90, lo60, hi60, lo30, hi30
    end

    rows = (
        (:cumulative_infections, "infections", :steelblue),
        (:cumulative_onsets, "symptom onsets", :seagreen),
        (:cumulative_expected_deaths, "deaths", :firebrick),
    )

    fig = Figure(; size = (940, 1020))
    for (i, (key, name, colour)) in enumerate(rows)
        trajs = _trajectories(key)
        lo90, hi90, lo60, hi60, lo30, hi30 = _ribbon(trajs)
        ax = Axis(
            fig[i, 1];
            xlabel = "Date", ylabel = "Cumulative $name",
            title = "Modelled cumulative $name over time",
            xticklabelrotation = pi / 6
        )
        band!(ax, x, lo90, hi90; color = (colour, 0.15))
        band!(ax, x, lo60, hi60; color = (colour, 0.28))
        band!(ax, x, lo30, hi30; color = (colour, 0.42))
        loax = floor(Int, minimum(x))
        hiax = ceil(Int, maximum(x))
        ax.xticks = collect(loax:14:hiax)
        ax.xtickformat = vals -> [
            string(epochdays2date(round(Int, v)))
                for v in vals
        ]

        ## Deliberately unbounded. These counts sit in the thousands, far
        ## from zero, so anchoring the axis there would squash the posterior
        ## into a spike against the right-hand edge.
        finals = Float64[t[n] for t in trajs]
        axd = Axis(
            fig[i, 2];
            xlabel = "Cumulative $name at the cut-off",
            ylabel = "Posterior density",
            title = "Current cumulative $name"
        )
        density!(
            axd, finals; color = (colour, 0.5),
            strokecolor = colour, strokewidth = 2
        )
    end
    return fig
end

## Open triangles marking where a series ran past a cropped axis, one per
## x-position given. Drawn a little inside the crop, since a marker placed
## exactly on the limit is cut in half by the plot area's clipping.
##
## Every cropped axis in this file marks its overflow through here, so the
## three figures that crop cannot drift apart in marker shape, size or
## placement.
function _mark_overflow!(ax, xs::AbstractVector, cap::Real, colour)
    isempty(xs) && return nothing
    return scatter!(
        ax, Float64.(collect(xs)), fill(0.97 * float(cap), length(xs));
        color = colour, marker = :utriangle, markersize = 10
    )
end

"""
Overlaid cumulative-infection trajectories, one per single-stream fit, each
projected out to the cut-off on day `n` even when that stream's data stops
earlier. Each stream is drawn as 30%, 60% and 90% credible ribbons with no
median line. A dotted vertical rule in each stream's colour marks the date
that stream's data stops reporting, so the projection beyond the data reads
apart from the fitted span.

Each `stream` is a `NamedTuple` `(; label, trajs, last_day, colour)`, where
`trajs` is a vector of per-draw cumulative-infection vectors of length `n`
(one per posterior draw) and `last_day` the 1-based grid day that stream's
data last reports (or `nothing` to omit the rule). `seeding` is the
calendar date of grid day 1, so day `d` is `seeding + (d - 1)`.

`ymax` crops the count axis. A stream whose data barely bounds the infection
count runs to a 90% upper on the scale of the source population, which on a
free axis flattens every other stream onto the baseline. Pass a multiple of
a reference fit's upper bound, as the cut-off density figure does, and the
streams that run past it are clipped and marked with an open triangle at the
day each one first leaves the axis rather than dropped. The default `nothing`
sizes the axis to the widest stream.
"""
function plot_stream_trajectories(
        streams::AbstractVector;
        n::Integer, seeding::Date,
        ymax::Union{Nothing, Real} = nothing,
        title::AbstractString =
            "Outbreak size projected to the cut-off by each data stream"
    )
    epoch = date2epochdays(seeding)
    x = Float64[epoch + (d - 1) for d in 1:n]

    fig = Figure(; size = (900, 480))
    ax = Axis(
        fig[1, 1];
        xlabel = "Date", ylabel = "Cumulative infections",
        title = title, xticklabelrotation = pi / 6
    )

    handles = Any[]
    labels = String[]
    cap = isnothing(ymax) ? nothing : float(ymax)
    ## The widest stream's 90% upper, which sizes the axis when no crop is
    ## asked for.
    datamax = 0.0
    for s in streams
        trajs = s.trajs
        q(d, pr) = quantile(Float64[t[d] for t in trajs], pr)
        lo90 = [q(d, 0.05) for d in 1:n]
        hi90 = [q(d, 0.95) for d in 1:n]
        lo60 = [q(d, 0.2) for d in 1:n]
        hi60 = [q(d, 0.8) for d in 1:n]
        lo30 = [q(d, 0.35) for d in 1:n]
        hi30 = [q(d, 0.65) for d in 1:n]
        datamax = max(datamax, maximum(hi90))
        colour = s.colour
        band!(ax, x, lo90, hi90; color = (colour, 0.12))
        band!(ax, x, lo60, hi60; color = (colour, 0.28))
        h = band!(ax, x, lo30, hi30; color = (colour, 0.42))
        push!(handles, h)
        push!(labels, s.label)
        ## Where a cropped axis is asked for, mark the day this stream's 90%
        ## upper first leaves it, so a stream running off the top reads as
        ## running off rather than as one that simply stops.
        if !isnothing(cap)
            over = findfirst(v -> v > cap, hi90)
            isnothing(over) || _mark_overflow!(ax, [x[over]], cap, colour)
        end
        ## Dotted rule in the stream's colour where its data stops reporting.
        ld = get(s, :last_day, nothing)
        ld === nothing || vlines!(
            ax, [Float64(epoch + ld - 1)];
            color = (colour, 0.8), linestyle = :dot, linewidth = 2
        )
    end

    loax = floor(Int, minimum(x))
    hiax = ceil(Int, maximum(x))
    ax.xticks = collect(loax:14:hiax)
    ax.xtickformat = vals -> [
        string(epochdays2date(round(Int, v)))
            for v in vals
    ]
    CairoMakie.ylims!(ax, 0, isnothing(cap) ? datamax * 1.05 : cap)
    CairoMakie.axislegend(
        ax, handles, labels; position = :lt,
        framevisible = true
    )
    return fig
end

"""
Overlaid posterior densities of an arbitrary scalar quantity from one
or more fits, built through AlgebraOfGraphics. Pass each fit as
`"label" => draws`; `xlabel` and `title` set the axis text.

`lower` clips the axis for a quantity that cannot fall below it, such as a
count or a duration. The kernel density spreads mass past the smallest draw,
so without it the curve runs onto the impossible side of the bound. The
estimate itself is left alone, as on `_bounded_density!`.
"""
function plot_density_overlay(
        streams::Pair{String, <:AbstractVector}...;
        xlabel::AbstractString = "Value",
        title::AbstractString = "Posterior density",
        lower::Union{Nothing, Real} = nothing
    )
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
    ax = isnothing(lower) ?
        (; ylabel = "Posterior density", title = title) :
        (;
            ylabel = "Posterior density", title = title,
            limits = ((float(lower), nothing), nothing),
        )
    return AoG.draw(spec; axis = ax, figure = (; size = (760, 420)))
end

_panel_pos(pos::Integer) = (1, pos)
_panel_pos(pos::Tuple) = pos

# Makie 0.24 computes data limits by calling `isfinite` elementwise, which has
# no method for integer vectors. Predictions of vector-valued observations
# arrive as a `Vector{Vector{Int}}`, so flatten any nesting and convert to
# `Float64` before plotting.
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
    ax = Axis(
        fig[r, c];
        xlabel = "Replicated exported cases",
        ylabel = "$(predictive_label) predictive frequency",
        title = "Exports (cases)",
        limits = ((0, upper), nothing)
    )
    hist!(ax, ppf; bins = 0:1:upper, color = (:steelblue, 0.7))
    vlines!(ax, _obs_floats(obs); color = :red, linewidth = 2)
    return ax
end

function _panel_exports_deaths!(
        fig, pos, pp, obs;
        predictive_label = "Posterior"
    )
    r, c = _panel_pos(pos)
    ppf = _pp_floats(pp)
    upper = max(3, ceil(Int, quantile(ppf, 0.995)))
    ax = Axis(
        fig[r, c];
        xlabel = "Replicated deaths among exports",
        ylabel = "$(predictive_label) predictive frequency",
        title = "Exports (deaths)",
        limits = ((0, upper), nothing)
    )
    hist!(ax, ppf; bins = 0:1:upper, color = (:rebeccapurple, 0.7))
    vlines!(ax, _obs_floats(obs); color = :red, linewidth = 2)
    return ax
end

function _panel_confirmed_deaths!(
        fig, pos, pp, obs;
        predictive_label = "Posterior"
    )
    r, c = _panel_pos(pos)
    ppf = _pp_floats(pp)
    upper = max(3, ceil(Int, quantile(ppf, 0.995)))
    obs === nothing || (upper = max(upper, ceil(Int, 1.1 * obs)))
    ax = Axis(
        fig[r, c];
        xlabel = "Replicated confirmed deaths",
        ylabel = "$(predictive_label) predictive frequency",
        title = "Confirmed deaths (DRC)",
        limits = ((0, upper), nothing)
    )
    hist!(ax, ppf; bins = 0:1:upper, color = (:darkorange3, 0.7))
    obs === nothing || vlines!(
        ax, _obs_floats(obs);
        color = :red, linewidth = 2
    )
    return ax
end

function _panel_deaths!(fig, pos, pp, obs; predictive_label = "Posterior")
    r, c = _panel_pos(pos)
    ppf = _pp_floats(pp)
    upper = max(1.0, quantile(ppf, 0.995))
    ax = Axis(
        fig[r, c];
        xlabel = "Replicated deaths",
        ylabel = "$(predictive_label) predictive frequency",
        title = "Deaths (DRC)",
        limits = ((0, upper), nothing)
    )
    hist!(
        ax, ppf; bins = range(0, upper; length = 40),
        color = (:firebrick, 0.7)
    )
    vlines!(ax, _obs_floats(obs); color = :red, linewidth = 2)
    return ax
end

function _panel_confirmed!(
        fig, pos, pp, obs;
        predictive_label = "Posterior"
    )
    r, c = _panel_pos(pos)
    ppf = _pp_floats(pp)
    upper = max(1.0, quantile(ppf, 0.995))
    if obs !== nothing
        upper = max(upper, 1.05 * maximum(_obs_floats(obs)))
    end
    ax = Axis(
        fig[r, c];
        xlabel = "Replicated confirmed cases",
        ylabel = "$(predictive_label) predictive frequency",
        title = "Confirmed cases (DRC)",
        limits = ((0, upper), nothing)
    )
    hist!(
        ax, ppf; bins = range(0, upper; length = 40),
        color = (:goldenrod, 0.7)
    )
    if obs !== nothing
        vlines!(ax, _obs_floats(obs); color = :red, linewidth = 2)
    end
    return ax
end

function _panel_tests!(
        fig, pos, pp, obs;
        predictive_label = "Posterior"
    )
    r, c = _panel_pos(pos)
    ppf = _pp_floats(pp)
    upper = max(1.0, quantile(ppf, 0.995))
    if obs !== nothing
        upper = max(upper, 1.05 * maximum(_obs_floats(obs)))
    end
    ax = Axis(
        fig[r, c];
        xlabel = "Replicated tests analysed",
        ylabel = "$(predictive_label) predictive frequency",
        title = "Tests analysed (DRC)",
        limits = ((0, upper), nothing)
    )
    hist!(
        ax, ppf; bins = range(0, upper; length = 40),
        color = (:teal, 0.7)
    )
    if obs !== nothing
        vlines!(ax, _obs_floats(obs); color = :red, linewidth = 2)
    end
    return ax
end

function _panel_cases!(fig, pos, pp, obs; predictive_label = "Posterior")
    r, c = _panel_pos(pos)
    ppf = _pp_floats(pp)
    upper = max(1.0, quantile(ppf, 0.995))
    ax = Axis(
        fig[r, c];
        xlabel = "Replicated reported cases",
        ylabel = "$(predictive_label) predictive frequency",
        title = "Reported cases (DRC)",
        limits = ((0, upper), nothing)
    )
    hist!(
        ax, ppf; bins = range(0, upper; length = 40),
        color = (:seagreen, 0.7)
    )
    if obs !== nothing
        vlines!(ax, _obs_floats(obs); color = :red, linewidth = 2)
    end
    return ax
end

"""
Posterior predictive histogram with one panel per supplied data stream. Pass
`pp_exports`/`pp_deaths` as `nothing` to suppress either of the first two
panels, and supply `pp_cases` and/or `pp_exports_deaths` to add the
reported-cases and deaths-among-exports panels. Observed values are drawn as
red `vlines`. Four or more streams lay out over three columns, fewer in a
single row.
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
        predictive_label::AbstractString = "Posterior"
    )
    panels = Tuple{Symbol, Any, Any}[]
    pp_exports === nothing ||
        push!(panels, (:exports, pp_exports, obs_exports))
    pp_exports_deaths === nothing ||
        push!(
        panels, (
            :exports_deaths, pp_exports_deaths,
            obs_exports_deaths,
        )
    )
    pp_confirmed_deaths === nothing ||
        push!(
        panels, (
            :confirmed_deaths, pp_confirmed_deaths,
            obs_confirmed_deaths,
        )
    )
    pp_deaths === nothing ||
        push!(panels, (:deaths, pp_deaths, obs_deaths))
    pp_cases === nothing ||
        push!(panels, (:cases, pp_cases, obs_cases))
    pp_tests === nothing ||
        push!(panels, (:tests, pp_tests, obs_tests))
    pp_confirmed === nothing ||
        push!(panels, (:confirmed, pp_confirmed, obs_confirmed))

    isempty(panels) && error(
        "plot_posterior_predictive needs at least one stream"
    )

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
    (:confirmed, _panel_confirmed!),
)

"""
Two-row comparison of posterior-predictive distributions, one column per
stream. The top row holds replicates from the per-stream fits, the bottom row
replicates from the joint fit conditioning on every observed stream. Observed
values are drawn as red vertical lines.

Each `NamedTuple` carries a subset of `(; exports, exports_deaths, deaths,
cases, tests, confirmed)`. Columns are drawn in that order for whichever
streams are present in `individual`. Each panel is a histogram of replicated
counts, and the two rows share one x-axis so they read against each other.
"""
function plot_posterior_predictive_grid(;
        individual::NamedTuple,
        joint::NamedTuple,
        observed::NamedTuple
    )
    streams = [
        (key, painter)
            for (key, painter) in _GRID_PANELS
            if hasproperty(individual, key)
    ]
    ncols = length(streams)
    fig = Figure(; size = (400 * ncols, 640))
    rows = (
        (:individual, individual, "per-stream fit"),
        (:joint, joint, "joint fit"),
    )
    for (i, (_, pp, label)) in enumerate(rows)
        for (j, (key, painter)) in enumerate(streams)
            painter(
                fig, (i, j), getproperty(pp, key),
                getproperty(observed, key); predictive_label = label
            )
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
        obs_tests::Union{Nothing, Real} = nothing
    )
    return plot_posterior_predictive(
        pp_exports, pp_deaths, obs_exports, obs_deaths;
        pp_cases, obs_cases, pp_confirmed, obs_confirmed,
        pp_tests, obs_tests, predictive_label = "Prior"
    )
end

"""
PairPlots.jl corner plot over the named posterior parameters, thinned by
`thin`. Pass `prior`, another chain holding the same parameters, to overlay
the prior as a second series with a legend, so the data's contribution to
each marginal is visible. Draws non-finite in any parameter are left out
of their series with a warning naming the parameter.

`labels` maps a raw chain symbol to a display name (e.g.
`Symbol("rt_state.sigma_rw") => "Rt step size"`), applied to the axis labels
only. Symbols absent from the map keep their raw name.

`patch` selects one entry of vector-valued deterministics such as
`R_T_patch` or `province_ascertainment`, so the corner plot shows one
province. Every parameter in `params` is then read as a per-patch vector.

`plot_pair(draws::NamedTuple; ...)` takes one draw vector per named quantity
instead of a chain, for quantities a chain holds only inside a vector
deterministic, such as one province's entry of `C_T_patch`. `prior` is then a
`NamedTuple` with the same names.
"""
function plot_pair(
        chn, params::AbstractVector{Symbol};
        thin::Integer = 2, prior = nothing,
        labels::AbstractDict = Dict{Symbol, String}(),
        patch::Union{Nothing, Integer} = nothing
    )
    _col(c, p) = patch === nothing ? _draws(c, p) :
        [v[patch] for v in _draw_vectors(c, p)]
    _named(c) = NamedTuple(p => _col(c, p) for p in params)
    return plot_pair(
        _named(chn); thin,
        prior = prior === nothing ? nothing : _named(prior), labels
    )
end

function plot_pair(
        draws::NamedTuple;
        thin::Integer = 2, prior::Union{Nothing, NamedTuple} = nothing,
        labels::AbstractDict = Dict{Symbol, String}()
    )
    _name(p) = Symbol(get(labels, p, string(p)))
    _table(d) = DataFrame(
        NamedTuple(_name(p) => v for (p, v) in pairs(d))
    )[_finite_rows(d), :][1:thin:end, :]
    post = _table(draws)
    prior === nothing && return PairPlots.pairplot(post)
    colours = CairoMakie.Makie.wong_colors()
    return PairPlots.pairplot(
        PairPlots.Series(post; label = "Posterior", color = colours[1]),
        PairPlots.Series(
            _table(prior); label = "Prior",
            color = colours[2]
        )
    )
end

## The draws finite in every quantity of `draws`. A pool the outbreak
## exhausts gives `R_T = 0` and so `r = -Inf`, which the density and
## histogram panels cannot bin. Each quantity with such draws is named.
function _finite_rows(draws::NamedTuple)
    keep = trues(length(first(draws)))
    for (p, v) in pairs(draws)
        bad = .!isfinite.(v)
        any(bad) || continue
        @warn "plot_pair drops $(count(bad)) non-finite draw(s) of $p"
        keep .&= .!bad
    end
    return keep
end

"""
Posterior correlation heatmap over the named scalar quantities `params`,
tracked deterministics or sampled parameters. Each cell is the Pearson
correlation of two quantities' posterior draws, drawn on a symmetric
red–blue scale with the value printed in the cell. The whole joint
identifiability structure sits in one panel, including the cross-block
degeneracies the [`plot_pair`](@ref) corners split apart: the
size–ascertainment seesaw (`C_T` vs `p_drc`), the weaker size–fatality tilt
(`C_T` vs `CFR`), and the pooled `p_drc`–`p_uganda` link.

`labels` maps a parameter symbol to a display string, read as LaTeX math
without the `\$` delimiters, so `"p_\\mathrm{drc}"` renders with a subscript.
A parameter absent from `labels` falls back to its symbol name. Returns the
`Figure`.

`plot_correlation_heatmap(draws::NamedTuple; labels)` takes one draw vector
per named quantity instead of a chain, for quantities a chain holds only
inside a vector deterministic, such as one province's entry of `C_T_patch`.
"""
function plot_correlation_heatmap(
        chn, params::AbstractVector{Symbol};
        labels::AbstractDict = Dict{Symbol, String}()
    )
    return plot_correlation_heatmap(
        NamedTuple(p => _draws(chn, p) for p in params); labels
    )
end

function plot_correlation_heatmap(
        draws::NamedTuple;
        labels::AbstractDict = Dict{Symbol, String}()
    )
    ## Render tick labels as LaTeX so subscripts (R_T, p_drc, λ_bg) typeset
    ## properly. Callers pass plain LaTeX math strings.
    name(p) = CairoMakie.Makie.latexstring(get(labels, p, string(p)))
    params = collect(keys(draws))
    mat = reduce(hcat, (float.(v) for v in values(draws)))
    R = cor(mat)
    n = length(params)
    labs = [name(p) for p in params]
    fig = Figure(; size = (78 * n + 180, 78 * n + 140))
    ax = Axis(
        fig[1, 1]; xticks = (1:n, labs), yticks = (1:n, labs),
        xticklabelrotation = pi / 4, title = "Posterior correlation",
        aspect = 1
    )
    hm = CairoMakie.heatmap!(
        ax, 1:n, 1:n, R; colormap = :RdBu,
        colorrange = (-1, 1)
    )
    for i in 1:n, j in 1:n

        CairoMakie.text!(
            ax, i, j;
            text = string(round(R[i, j]; digits = 2)),
            align = (:center, :center), fontsize = 10,
            color = abs(R[i, j]) > 0.6 ? :white : :black
        )
    end
    CairoMakie.Colorbar(fig[1, 2], hm)
    return fig
end

"""
Pairs plot of the per-stream modelled totals. `modelled` is a `NamedTuple` of
one per-draw vector per stream, each summed to that stream's own observed
support, and `observed` a `NamedTuple` of scalars drawn as crosshair
reference lines.

The diagonals show how much predictive density sits above or below each
observed value, and the off-diagonals whether those shoots move together
across draws. This is the posterior-predictive view of data-stream conflict,
the counterpart to the parameter-space
[`plot_correlation_heatmap`](@ref). Returns the `Figure`.
"""
function plot_stream_pairs(modelled::NamedTuple, observed::NamedTuple)
    return PairPlots.pairplot(
        modelled,
        PairPlots.Truth(observed; label = "observed")
    )
end

## Calendar dates carried by one panel's three series.
function _evolution_dates(released, renewal, trajectory)
    rdates = [Date(String(r[1])) for r in released]
    ndates = [Date(String(p[1])) for p in renewal]
    tdates = isnothing(trajectory) ? Date[] :
        [d isa Date ? d : Date(String(d)) for d in trajectory[1]]
    return rdates, ndates, tdates
end

## Top of one panel's series, the 90% upper over every discrete mark and
## over the trajectory's upper band.
function _evolution_upper(released, renewal, trajectory)
    _hi(ts) = isempty(ts) ? 0.0 : maximum(float(t[8]) for t in ts)
    upper = max(_hi(released), _hi(renewal))
    isnothing(trajectory) ||
        (upper = max(upper, maximum(float.(trajectory[7]))))
    return upper
end

## Draw one estimate-evolution panel into `ax`, shared by the single-axis and
## faceted plots. `_x` maps a calendar date to the axis' numeric day-offset
## and is passed in so faceted panels share one mapping. `cap` is the panel's
## value-axis crop: every interval, median and band is clamped to it and a
## series that runs past it is marked with an open triangle, so a wide fit
## reads as running off the axis rather than being dropped. Returns the
## legend handles and labels for the series actually drawn.
function _evolution_panel!(
        ax, _x, released, renewal, trajectory,
        refline, labels::NamedTuple;
        cap::Union{Nothing, Real} = nothing
    )
    _clamp(v) = isnothing(cap) ? float(v) : min(float(v), float(cap))
    _over(v) = !isnothing(cap) && float(v) > float(cap)
    overflow_x = Dict{Symbol, Vector{Float64}}()
    rdates, ndates, tdates = _evolution_dates(released, renewal, trajectory)

    ## Each release and each frozen re-fit is its own fit, so collect them as
    ## discrete marks and dodge any that share a date.
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
            append!(by, (_clamp(lo), _clamp(hi)))
        end
        return linesegments!(
            ax, bx, by;
            color = (colour, alpha), linewidth = lw
        )
    end

    ## One discrete series: nested 30/60/90% bars topped by a median dot.
    function _series!(colour)
        sel = [i for i in eachindex(marks) if marks[i].colour == colour]
        isempty(sel) && return nothing
        xs = markx[sel]
        ts = [marks[i].t for i in sel]
        _bars!(
            xs, [float(t[7]) for t in ts], [float(t[8]) for t in ts],
            colour, 1.4, 0.45
        )
        _bars!(
            xs, [float(t[5]) for t in ts], [float(t[6]) for t in ts],
            colour, 3.2, 0.55
        )
        _bars!(
            xs, [float(t[3]) for t in ts], [float(t[4]) for t in ts],
            colour, 6.5, 0.7
        )
        ## The 90% upper and the median are what a reader takes off the
        ## panel, so either running past the crop is what the marker notes.
        for (x, t) in zip(xs, ts)
            (_over(t[8]) || _over(t[2])) &&
                push!(get!(overflow_x, colour, Float64[]), x)
        end
        return scatter!(
            ax, xs, [_clamp(t[2]) for t in ts];
            color = colour, markersize = 9,
            strokecolor = :white, strokewidth = 1
        )
    end

    handles = Any[]
    llabels = String[]
    ## Current-data estimate as the cumulative-infection trajectory over the
    ## grid, one continuous ribbon. `trajectory` is
    ## `(dates, lo30, hi30, lo60, hi60, lo90, hi90)`.
    if !isnothing(trajectory)
        cc = :seagreen
        xs = _x.(tdates)
        ord = sortperm(xs)
        xs = xs[ord]
        lo30 = float.(trajectory[2])[ord]
        hi30 = float.(trajectory[3])[ord]
        lo60 = float.(trajectory[4])[ord]
        hi60 = float.(trajectory[5])[ord]
        lo90 = float.(trajectory[6])[ord]
        hi90 = float.(trajectory[7])[ord]
        ## A trajectory whose dates collapse to a single x-position would
        ## draw as an invisible zero-width band, so widen it across the span
        ## of the discrete marks.
        if length(unique(xs)) <= 1 && !(isempty(rdates) && isempty(ndates))
            markx = _x.(vcat(rdates, ndates))
            wlo, whi = minimum(markx), maximum(markx)
            wlo == whi && (wlo -= 0.5; whi += 0.5)
            xs = [wlo, whi]
            lo30 = fill(lo30[1], 2)
            hi30 = fill(hi30[1], 2)
            lo60 = fill(lo60[1], 2)
            hi60 = fill(hi60[1], 2)
            lo90 = fill(lo90[1], 2)
            hi90 = fill(hi90[1], 2)
        end
        any(_over, hi90) &&
            push!(
            get!(overflow_x, cc, Float64[]),
            xs[findfirst(_over, hi90)]
        )
        band!(ax, xs, _clamp.(lo90), _clamp.(hi90); color = (cc, 0.1))
        band!(ax, xs, _clamp.(lo60), _clamp.(hi60); color = (cc, 0.16))
        th = band!(ax, xs, _clamp.(lo30), _clamp.(hi30); color = (cc, 0.24))
        push!(handles, th)
        push!(llabels, labels.trajectory * " (30/60/90% band)")
    end
    renewal_mark = _series!(:firebrick)
    if !isnothing(renewal_mark)
        push!(handles, renewal_mark)
        push!(llabels, labels.renewal * " (median, 30/60/90% bars)")
    end
    released_mark = _series!(:steelblue)
    if !isnothing(released_mark)
        push!(handles, released_mark)
        push!(llabels, labels.released * " (median, 30/60/90% bars)")
    end

    ## Dotted vertical rule at each release date.
    isempty(rdates) || vlines!(
        ax, _x.(rdates);
        color = (:grey, 0.55), linestyle = :dot, linewidth = 1
    )

    ## Optional horizontal reference line, e.g. Rt = 1 for a reproduction
    ## number, drawn faint so it reads behind the estimates. A reference
    ## line past the crop is not drawn at all rather than clamped onto it,
    ## since a rule sitting at the top of the axis reads as a reference at
    ## that value.
    (isnothing(refline) || _over(refline)) ||
        hlines!(
        ax, [float(refline)];
        color = (:black, 0.4), linestyle = :dash, linewidth = 1
    )

    for (colour, xs) in overflow_x
        _mark_overflow!(ax, xs, cap, colour)
    end
    return handles, llabels
end

"""
Estimate-evolution plot: how the outbreak-size estimate moves as the data
cut-off advances, drawn against the calendar date.

`released` is a vector of `(cutoff_date, median, lo30, hi30, lo60, hi60,
lo90, hi90)` tuples, one per published project release, drawn in blue.
`renewal` is the same tuple shape for the current renewal model re-fit frozen
at each release date, drawn in red. Each release and each frozen re-fit is
its own fit, so both are drawn as discrete per-date estimates, a median
marker with nested 30/60/90% vertical interval bars. Marks sharing a date are
dodged horizontally so each reads as a separate estimate.

`trajectory` is the current-data, current-model cumulative-infection
trajectory over the day grid, a `(dates, lo30, hi30, lo60, hi60, lo90, hi90)`
tuple where `dates` is the calendar date of each grid day. It is a single fit
shown over time, so it is drawn in a third colour as one continuous ribbon on
the same calendar axis. When its dates collapse to a single day the ribbon is
widened across the span of the discrete marks, so a flat reference still
reads as a band.

Release dates are marked with dotted vertical rules.

`xlabel`, `ylabel` and `title` set the axis text. `released_label`,
`renewal_label` and `trajectory_label` name the three series.

`ymax` fixes the value axis rather than sizing it to the widest interval.
An interval, median or band past it is clamped and the series marked with
an open triangle where it leaves the axis, so a quantity whose estimates
sit in a narrow range is not flattened by one wide tail.
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
        refline::Union{Nothing, Real} = nothing,
        ymax::Union{Nothing, Real} = nothing
    )
    ## Calendar dates → numeric day-offsets so the x-axis is to scale, then
    ## relabel the ticks with the dates. All three series share this one
    ## mapping, so they line up on the same axis.
    rdates, ndates, tdates = _evolution_dates(released, renewal, trajectory)
    tickdates = sort(unique(vcat(rdates, ndates)))
    alldates = sort(unique(vcat(rdates, ndates, tdates)))
    ref = minimum(alldates)
    _x(d) = Float64((d - ref).value)

    upper = isnothing(ymax) ?
        _evolution_upper(released, renewal, trajectory) * 1.08 : float(ymax)
    xlo = _x(ref) - 1
    xhi = _x(maximum(alldates)) + 1
    fig = Figure(; size = (860, 480))
    ax = Axis(
        fig[1, 1];
        xlabel = xlabel, ylabel = ylabel, title = title,
        xticks = (_x.(tickdates), [string(d) for d in tickdates]),
        xticklabelrotation = pi / 4,
        limits = ((xlo, xhi), (0, upper))
    )

    handles, labels = _evolution_panel!(
        ax, _x, released, renewal,
        trajectory, refline,
        (;
            released = released_label, renewal = renewal_label,
            trajectory = trajectory_label,
        );
        cap = ymax
    )

    CairoMakie.axislegend(
        ax, handles, labels; position = :lt,
        framevisible = true
    )
    return fig
end

"""
Faceted estimate evolution, one panel per group, for comparing how the same
quantity moved across releases under each of several fits.

`groups` is a vector of `label => released` pairs, where `released` is the
`(cutoff_date, median, lo30, hi30, lo60, hi60, lo90, hi90)` tuple vector
`plot_estimate_evolution` takes. Panels share one calendar mapping and are
laid out over `ncols` columns. Groups with no estimates are dropped rather
than drawn as an empty panel. `refline` draws a faint horizontal rule in
every panel, e.g. Rt = 1.

`trajectories` optionally maps a group's label to that group's own
current-data, current-model trajectory, a `(dates, lo30, hi30, lo60, hi60,
lo90, hi90)` tuple in the shape `plot_estimate_evolution` takes, drawn with
the same band styling. A group absent from `trajectories` still draws its
released points with no band.

`shared_yrange` (default `true`) draws every panel on one y range, sized to
the largest of them, so groups on a comparable scale read against each other.
Set it to `false` when the groups span very different scales, so a wide-scale
group does not squash every other panel's band. Each panel then uses its own
range, floored at `1.0`.

`ymax` overrides both, fixing every panel's value axis at it. An interval,
median or band past it is clamped and marked with an open triangle, as in
[`plot_estimate_evolution`](@ref).

Returns a figure carrying `empty_note` in place of the panels when no group
has any estimate.
"""
function plot_evolution_by_group(
        groups::AbstractVector;
        trajectories::AbstractDict = Dict{Any, Any}(),
        xlabel::AbstractString = "Date",
        ylabel::AbstractString = "Estimate",
        title::AbstractString = "",
        released_label::AbstractString = "Released estimate (per release)",
        trajectory_label::AbstractString = "Current model, current data",
        refline::Union{Nothing, Real} = nothing,
        ncols::Int = 2,
        shared_yrange::Bool = true,
        ymax::Union{Nothing, Real} = nothing,
        empty_note::AbstractString = "No per-dataset estimates yet."
    )
    ## A group with no estimates is dropped, so the panels show only fits
    ## that exist.
    shown = [g for g in groups if !isempty(last(g))]
    if isempty(shown)
        fig = Figure(; size = (860, 160))
        CairoMakie.Label(
            fig[1, 1], empty_note;
            tellwidth = false, tellheight = false, color = (:black, 0.55)
        )
        return fig
    end
    _traj(g) = get(trajectories, first(g), nothing)

    ## One calendar mapping and one y range across every panel, so a
    ## dataset's estimates read against the others. Trajectory dates are
    ## folded in, so a band reaching back past any release point still fits.
    function _group_dates(g)
        rd, _, td = _evolution_dates(last(g), NamedTuple[], _traj(g))
        return vcat(rd, td)
    end
    alldates = sort(unique(reduce(vcat, [_group_dates(g) for g in shown])))
    ref = minimum(alldates)
    _x(d) = Float64((d - ref).value)
    shared_upper = max(
        1.0,
        maximum(
            _evolution_upper(last(g), NamedTuple[], _traj(g))
                for g in shown
        )
    )
    xlo = _x(ref) - 1
    xhi = _x(maximum(alldates)) + 1

    ## Facet panels are narrow, so label a subset of the release dates.
    step = max(1, cld(length(alldates), 5))
    tickdates = alldates[1:step:end]

    nrows = cld(length(shown), ncols)
    fig = Figure(; size = (460 * ncols, 250 * nrows + 90))
    ## Collect each series' legend handle wherever it first appears. Panels
    ## draw different subsets of the series, and the legend must cover every
    ## one drawn in any panel.
    handle_map = Dict{String, Any}()
    order = String[]
    for (i, g) in enumerate(shown)
        r, c = fldmod1(i, ncols)
        panel_upper = if !isnothing(ymax)
            float(ymax)
        elseif shared_yrange
            shared_upper * 1.08
        else
            max(1.0, _evolution_upper(last(g), NamedTuple[], _traj(g))) * 1.08
        end
        ax = Axis(
            fig[r, c]; title = string(first(g)),
            xlabel = r == nrows ? xlabel : "",
            ylabel = c == 1 ? ylabel : "",
            xticks = (_x.(tickdates), [string(d) for d in tickdates]),
            xticklabelrotation = pi / 4,
            limits = ((xlo, xhi), (0, panel_upper))
        )
        h, l = _evolution_panel!(
            ax, _x, last(g), NamedTuple[], _traj(g),
            refline, (;
                released = released_label, renewal = "",
                trajectory = trajectory_label,
            );
            cap = ymax
        )
        for (hh, ll) in zip(h, l)
            haskey(handle_map, ll) && continue
            handle_map[ll] = hh
            push!(order, ll)
        end
    end
    handles = [handle_map[l] for l in order]
    labels = order

    isempty(title) || CairoMakie.Label(
        fig[0, 1:ncols], title;
        font = :bold, tellwidth = false
    )
    isempty(handles) ||
        CairoMakie.Legend(
        fig[nrows + 1, 1:ncols], handles, labels;
        orientation = :horizontal, framevisible = true,
        tellheight = true, tellwidth = false
    )
    return fig
end

"""
Forecasts-versus-now overlay: a grid of panels, one row per observed stream
and one column per forecast horizon, each showing the forecasts made at every
release (median with the 90% predictive interval as a vertical bar) against
the value observed since (black points). `overlay` is the
`data/forecast_overlay.csv` table with columns `stream`, `made_date`,
`horizon`, `fit`, `observed`, `median`, `lo90` and `hi90`.

The x-axis is the date the forecast was made, not the target date. For an
incident stream the observed value is the new count over the forecast's own
`(made_date, target_date]` window, so it depends on the made date. Keying on
the target date would stack several different observed windows at one x with
no way to pair each forecast to its own truth.

Forecasts are coloured by fit role (`baseline`, `individual`, `joint`) and
dodged by a small fraction of the made-date spacing so the series read apart.
A stream that carries only some roles is drawn with those alone, and the
legend covers every role drawn in any panel.

Each panel is cropped to zero and about three times the larger of its own
observed values and forecast medians, so one very wide predictive interval
elsewhere cannot squash every other panel to a line near the bottom. An
interval or median past that crop is clamped and marked with an open triangle
at the top of the axis. Every made date gets its own x tick, thinning to
about a dozen for a busier release history.

Returns a figure carrying `empty_message` in place of the panels when no
forecasts have been scored yet.
"""
function plot_forecast_overlay(
        overlay::DataFrame;
        empty_message::AbstractString =
            "No forecasts scored yet. No release carries a stored forecast."
    )
    streams = unique(overlay.stream)
    ## An empty table is the expected early state, so say so rather than
    ## returning a blank panel.
    if isempty(streams)
        fig = Figure(; size = (860, 160))
        CairoMakie.Label(
            fig[1, 1], empty_message;
            tellwidth = false, tellheight = false, color = (:black, 0.55)
        )
        return fig
    end
    horizons = sort(unique(overlay.horizon))
    role_order = ["baseline", "individual", "joint"]
    role_colour = Dict(
        "baseline" => :goldenrod, "individual" => :steelblue,
        "joint" => :firebrick
    )

    ## One made-date axis shared across every cell, with the made dates
    ## evenly spaced rather than placed to calendar scale. The releases are a
    ## handful of discrete cut-offs, and to scale the ones a few days apart
    ## land on top of each other with their rotated labels overprinting.
    alldates = sort(unique(Date.(string.(overlay.made_date))))
    slot_of = Dict(d => Float64(i) for (i, d) in enumerate(alldates))
    _x(d) = slot_of[Date(string(d))]
    spacing = 1.0
    ## A fraction of the made-date spacing, so the observed point and the
    ## three fit roles read apart at each made date.
    dodge = 0.1 * spacing
    ## Every made date gets its own tick, thinned to about sixteen once the
    ## release history is long enough to crush the labels.
    step = length(alldates) <= 16 ? 1 : cld(length(alldates), 16)
    tickdates = alldates[1:step:end]

    ## Fixed x-slot per series, so the observed truth and the three fit roles
    ## each sit at their own offset from the made date, and the slot is stable
    ## across cells whichever roles are present.
    slot = Dict(
        "observed" => -1.5, "baseline" => -0.5, "individual" => 0.5,
        "joint" => 1.5
    )

    ncols = length(horizons)
    nrows = length(streams)
    fig = Figure(; size = (240 * ncols + 80, 200 * nrows + 90))

    ## A stream may lack a role, e.g. recovered has no individual fit, so
    ## collect each role's legend handle wherever it first appears.
    obs_handle = nothing
    role_handles = Dict{String, Any}()
    for (r, s) in enumerate(streams), (c, h) in enumerate(horizons)

        cell = overlay[(overlay.stream .== s) .& (overlay.horizon .== h), :]
        ## Crop to zero and about three times the larger of this panel's own
        ## observed values and forecast medians. The bare floor keeps an
        ## empty or all-zero panel off a zero-height axis.
        cap = isempty(cell) ? 1.0 :
            max(1.0, 3.0 * max(maximum(cell.observed), maximum(cell.median)))
        ax = Axis(
            fig[r, c];
            title = r == 1 ? "$(h)-day ahead" : "",
            xlabel = r == nrows ? "Forecast made" : "",
            ylabel = c == 1 ? string(s) : "",
            xticks = (_x.(tickdates), [string(d) for d in tickdates]),
            xticklabelrotation = pi / 4,
            limits = ((0.5, length(alldates) + 0.5), (0, cap))
        )
        ## Every panel shares the one made-date axis, so only the bottom row
        ## carries the rotated date labels.
        r == nrows || CairoMakie.hidexdecorations!(
            ax;
            ticklabels = true, ticks = false, grid = false,
            label = false, minorgrid = false, minorticks = false
        )
        ## A stream and horizon with nothing scored says so rather than
        ## presenting empty axes. Its value ticks go too, there being no
        ## scale to read.
        if isempty(cell)
            CairoMakie.hideydecorations!(ax)
            CairoMakie.text!(
                ax, (length(alldates) + 1) / 2, cap / 2;
                text = "not scored", align = (:center, :center),
                color = (:black, 0.4), fontsize = 10
            )
            continue
        end
        ## Observed new count over each forecast's own window, one value per
        ## made date within this horizon.
        od = sort(
            unique(
                [
                    (Date(string(m)), float(o))
                        for (m, o) in zip(cell.made_date, cell.observed)
                ]
            )
        )
        oh = scatter!(
            ax,
            [_x(m) + slot["observed"] * dodge for (m, _) in od],
            [o for (_, o) in od]; color = :black, markersize = 7
        )
        isnothing(obs_handle) && (obs_handle = oh)
        for role in role_order
            rs = select_fit_role(cell, role)
            isempty(rs) && continue
            col = role_colour[role]
            off = slot[role] * dodge
            xs = [_x(m) + off for m in rs.made_date]
            ## Clamp any bound past the crop and note where the interval or
            ## median ran off the top, so a wide forecast reads as running
            ## off the axis rather than being dropped.
            bx = Float64[]
            by = Float64[]
            meds = Float64[]
            overflow_x = Float64[]
            for (x, lo, hi, med) in zip(xs, rs.lo90, rs.hi90, rs.median)
                append!(bx, (x, x))
                append!(by, (min(float(lo), cap), min(float(hi), cap)))
                push!(meds, min(float(med), cap))
                (float(hi) > cap || float(med) > cap) && push!(overflow_x, x)
            end
            linesegments!(ax, bx, by; color = (col, 0.45), linewidth = 2)
            mh = scatter!(ax, xs, meds; color = col, markersize = 6)
            _mark_overflow!(ax, overflow_x, cap, col)
            get!(role_handles, role, mh)
        end
    end

    handles = Any[]
    labels = String[]
    if !isnothing(obs_handle)
        push!(handles, obs_handle)
        push!(labels, "observed")
    end
    for role in role_order
        haskey(role_handles, role) || continue
        push!(handles, role_handles[role])
        push!(labels, role)
    end
    isempty(handles) ||
        CairoMakie.Legend(
        fig[nrows + 1, 1:ncols], handles, labels;
        orientation = :horizontal, framevisible = true,
        tellheight = true, tellwidth = false
    )
    return fig
end

## Ratio ticks for a logarithmic skill axis, labelled as the ratios
## themselves. Makie places log ticks on round powers of ten, which over a
## range of a few leaves only fractional exponents to label.
const _SKILL_TICKS = [
    0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0,
    20.0, 50.0, 100.0, 200.0, 500.0,
]
const _skill_tick_labels = [
    t >= 1 ? string(round(Int, t)) : string(t)
        for t in _SKILL_TICKS
]

"""
By-horizon relative-skill figure: one panel per stream, plotting relative
skill against the persistence baseline (`rel_to_baseline`, or
`log_rel_to_baseline` on the log scale) against the forecast horizon, with
one series per fit role (`individual`, `joint`). `scores` is a
[`forecast_score_by_horizon`](@ref)-shaped table, baseline rows already
excluded, carrying `stream`, `horizon`, `fit` and the column named by
`value_col`.

The skill axis is log-scaled, so a fit twice as good and a fit twice as bad
sit the same distance from the reference line at one, drawn as a dashed
horizontal rule. A series below the line beats the baseline on average at
that horizon. A `(stream, horizon, fit)` cell whose skill is missing or
non-finite is absent from its series rather than drawn as a break, and a
stream with no individual fit is drawn with the joint series alone.

`empty_message` is shown in place of the panels when `scores` has no rows.
The caller sets it, since an empty table can mean either that no release
carries a stored forecast or that no stored forecast's target has resolved.
"""
function plot_forecast_relative_skill(
        scores::DataFrame;
        value_col::Symbol = :rel_to_baseline,
        ylabel::AbstractString = "Relative skill (log scale, 1 = baseline)",
        title::AbstractString =
            "Relative skill against the baseline, by horizon",
        ncols::Integer = 3,
        empty_message::AbstractString =
            "No scored forecasts yet. No release carries a stored forecast."
    )
    streams = sort(unique(scores.stream))
    if isempty(streams)
        fig = Figure(; size = (860, 160))
        CairoMakie.Label(
            fig[1, 1], empty_message;
            tellwidth = false, tellheight = false, color = (:black, 0.55)
        )
        return fig
    end
    role_order = ["individual", "joint"]
    role_colour = Dict("individual" => :steelblue, "joint" => :firebrick)
    horizons = sort(unique(scores.horizon))

    usedcols = min(ncols, length(streams))
    nrows = cld(length(streams), usedcols)
    fig = Figure(; size = (320 * usedcols, 260 * nrows + 90))

    role_handles = Dict{String, Any}()
    for (i, s) in enumerate(streams)
        r, c = fldmod1(i, usedcols)
        cell = scores[scores.stream .== s, :]
        ## Skill is a ratio, so the axis is logarithmic to put a factor of
        ## two better and a factor of two worse the same distance from one.
        ax = Axis(
            fig[r, c]; title = string(s),
            xlabel = r == nrows ? "Forecast horizon (days)" : "",
            ylabel = c == 1 ? ylabel : "",
            xticks = horizons, yscale = log10,
            yticks = (_SKILL_TICKS, _skill_tick_labels)
        )
        hlines!(
            ax, [1.0]; color = (:grey, 0.6), linestyle = :dash,
            linewidth = 2
        )
        for role in role_order
            rs = select_fit_role(cell, role)
            isempty(rs) && continue
            keep = [!ismissing(v) && isfinite(v) for v in rs[!, value_col]]
            any(keep) || continue
            xs = Float64.(rs.horizon[keep])
            ys = Float64.(rs[keep, value_col])
            ord = sortperm(xs)
            h = scatterlines!(
                ax, xs[ord], ys[ord];
                color = role_colour[role], markersize = 8, linewidth = 2
            )
            get!(role_handles, role, h)
        end
    end

    handles = Any[]
    labels = String[]
    for role in role_order
        haskey(role_handles, role) || continue
        push!(handles, role_handles[role])
        push!(labels, role)
    end
    isempty(title) || CairoMakie.Label(
        fig[0, 1:usedcols], title;
        font = :bold, tellwidth = false
    )
    isempty(handles) ||
        CairoMakie.Legend(
        fig[nrows + 1, 1:usedcols], handles, labels;
        orientation = :horizontal, framevisible = true,
        tellheight = true, tellwidth = false
    )
    return fig
end

"""
By-vintage relative-skill figure: one panel per stream, plotting relative
skill against the persistence baseline against the release that made the
forecast, with one series per fit role. `scores` is a
[`forecast_score_by_vintage`](@ref)-shaped table, carrying `stream`,
`release`, `release_date`, `fit` and the column named by `value_col`.

Releases sit at evenly spaced slots in `release_date` order, not to
calendar scale, and are labelled with the date each was cut. The skill
axis is log-scaled about a reference line at one, as in
[`plot_forecast_relative_skill`](@ref).

A cell whose skill is missing or non-finite is absent from its series.
`empty_message` replaces the panels when `scores` has no rows.
"""
function plot_forecast_skill_by_vintage(
        scores::DataFrame;
        value_col::Symbol = :rel_to_baseline,
        ylabel::AbstractString = "Relative skill (log scale, 1 = baseline)",
        title::AbstractString =
            "Relative skill against the baseline, by release",
        ncols::Integer = 3,
        empty_message::AbstractString =
            "No cut-off has been forecast by more than one release yet."
    )
    ## One shared slot per release across every panel, so a stream missing a
    ## release leaves a gap rather than shifting against the others.
    keys_of(tbl) = tbl.release
    order_of(tbl) = Dict(
        tbl.release[i] => (tbl.release_date[i], tbl.release[i])
            for i in 1:size(tbl, 1)
    )
    return _plot_skill_by_slot(
        scores; keys_of = keys_of, order_of = order_of,
        label_of = k -> string(first(k)),
        xlabel = "Release cut from", value_col = value_col,
        ylabel = ylabel, title = title, ncols = ncols,
        empty_message = empty_message
    )
end

"""
By-cut-off relative-skill figure: one panel per stream, plotting relative
skill against the persistence baseline against the date the forecast was
made, with one series per fit role. `scores` is a
[`forecast_score_by_release`](@ref)-shaped table, carrying `stream`,
`made_date`, `fit` and the column named by `value_col`.

The made dates sit at evenly spaced slots in date order, not to calendar
scale, since the cut-offs are a handful of discrete dates and the ones a
few days apart would otherwise overprint. The skill axis is log-scaled
about a reference line at one, as in
[`plot_forecast_relative_skill`](@ref).

This is the figure for the by-cut-off (or by-release) score tables, which
carry one row per stream and cut-off and are too long to read as numbers.

A cell whose skill is missing or non-finite is absent from its series.
`empty_message` replaces the panels when `scores` has no rows.
"""
function plot_forecast_skill_by_cutoff(
        scores::DataFrame;
        value_col::Symbol = :rel_to_baseline,
        ylabel::AbstractString = "Relative skill (log scale, 1 = baseline)",
        xlabel::AbstractString = "Forecast made",
        title::AbstractString =
            "Relative skill against the baseline, by cut-off",
        ncols::Integer = 3,
        empty_message::AbstractString =
            "No scored forecasts yet. No release carries a stored forecast."
    )
    keys_of(tbl) = [Date(string(d)) for d in tbl.made_date]
    order_of(tbl) = Dict(
        Date(string(d)) => (Date(string(d)),) for d in tbl.made_date
    )
    return _plot_skill_by_slot(
        scores; keys_of = keys_of, order_of = order_of,
        label_of = k -> string(first(k)),
        xlabel = xlabel, value_col = value_col,
        ylabel = ylabel, title = title, ncols = ncols,
        empty_message = empty_message
    )
end

## Shared body of the two by-slot skill figures. `keys_of` takes a table to
## the per-row key each point sits at, `order_of` to a key => sort-tuple
## mapping (the tuple's first entry is the tick label's source) and
## `label_of` to a tick label. Keeping one body means the release and
## cut-off figures cannot drift apart in their axis, colours or legend.
function _plot_skill_by_slot(
        scores::DataFrame;
        keys_of, order_of, label_of, xlabel::AbstractString,
        value_col::Symbol, ylabel::AbstractString,
        title::AbstractString, ncols::Integer,
        empty_message::AbstractString
    )
    streams = sort(unique(scores.stream))
    if isempty(streams)
        fig = Figure(; size = (860, 160))
        CairoMakie.Label(
            fig[1, 1], empty_message;
            tellwidth = false, tellheight = false, color = (:black, 0.55)
        )
        return fig
    end
    role_order = ["individual", "joint"]
    role_colour = Dict("individual" => :steelblue, "joint" => :firebrick)

    ordering = order_of(scores)
    slots = sort(collect(keys(ordering)); by = k -> ordering[k])
    slot = Dict(k => Float64(i) for (i, k) in enumerate(slots))
    ## Thinned to about eight labels once the history is long enough.
    step = length(slots) <= 8 ? 1 : cld(length(slots), 8)
    ticks = 1:step:length(slots)
    ticklabels = [label_of(ordering[slots[i]]) for i in ticks]

    usedcols = min(ncols, length(streams))
    nrows = cld(length(streams), usedcols)
    fig = Figure(; size = (340 * usedcols, 260 * nrows + 90))

    role_handles = Dict{String, Any}()
    for (i, s) in enumerate(streams)
        r, c = fldmod1(i, usedcols)
        cell = scores[scores.stream .== s, :]
        ax = Axis(
            fig[r, c]; title = string(s),
            xlabel = r == nrows ? xlabel : "",
            ylabel = c == 1 ? ylabel : "",
            xticks = (Float64.(collect(ticks)), ticklabels),
            xticklabelrotation = pi / 4,
            yscale = log10,
            yticks = (_SKILL_TICKS, _skill_tick_labels),
            limits = ((0.5, length(slots) + 0.5), nothing)
        )
        hlines!(
            ax, [1.0]; color = (:grey, 0.6), linestyle = :dash,
            linewidth = 2
        )
        for role in role_order
            rs = select_fit_role(cell, role)
            isempty(rs) && continue
            keep = [!ismissing(v) && isfinite(v) for v in rs[!, value_col]]
            any(keep) || continue
            xs = [slot[k] for k in keys_of(rs)[keep]]
            ys = Float64.(rs[keep, value_col])
            ord = sortperm(xs)
            h = scatterlines!(
                ax, xs[ord], ys[ord];
                color = role_colour[role], markersize = 8, linewidth = 2
            )
            get!(role_handles, role, h)
        end
    end

    handles = Any[]
    labels = String[]
    for role in role_order
        haskey(role_handles, role) || continue
        push!(handles, role_handles[role])
        push!(labels, role)
    end
    isempty(title) || CairoMakie.Label(
        fig[0, 1:usedcols], title;
        font = :bold, tellwidth = false
    )
    isempty(handles) ||
        CairoMakie.Legend(
        fig[nrows + 1, 1:usedcols], handles, labels;
        orientation = :horizontal, framevisible = true,
        tellheight = true, tellwidth = false
    )
    return fig
end

## Colour per CRPS component, shared by the by-horizon decomposition figure
## and its legend so the two cannot disagree.
const _CRPS_PARTS = (
    (:dispersion, "dispersion", :slategray3),
    (:overprediction, "overprediction", :goldenrod),
    (:underprediction, "underprediction", :firebrick),
)

"""
By-horizon CRPS figure: one panel per stream, the mean CRPS split into the
three parts it decomposes into (dispersion, overprediction and
underprediction) and drawn as one stacked bar per horizon, dodged by fit
role where a stream carries more than one.

`scores` is a [`forecast_score_by_horizon`](@ref)-shaped table carrying
`stream`, `horizon`, `fit` and the three component columns. The bar's full
height is the mean CRPS, so a panel shows both how the error grows with the
horizon and whether it is made of width, of forecasting too high or of
forecasting too low. This is the figure for the by-horizon score table,
which is too long to read as numbers.

Each panel takes its own y range, since a stream's CRPS is on the scale of
its own counts. `empty_message` replaces the panels when `scores` has no
rows.
"""
function plot_forecast_crps_by_horizon(
        scores::DataFrame;
        ylabel::AbstractString = "Mean CRPS",
        title::AbstractString = "CRPS decomposition by horizon",
        ncols::Integer = 3,
        empty_message::AbstractString =
            "No scored forecasts yet. No release carries a stored forecast."
    )
    streams = sort(unique(scores.stream))
    if isempty(streams)
        fig = Figure(; size = (860, 160))
        CairoMakie.Label(
            fig[1, 1], empty_message;
            tellwidth = false, tellheight = false, color = (:black, 0.55)
        )
        return fig
    end
    role_order = ["individual", "joint"]
    horizons = sort(unique(scores.horizon))
    hslot = Dict(h => i for (i, h) in enumerate(horizons))

    usedcols = min(ncols, length(streams))
    nrows = cld(length(streams), usedcols)
    fig = Figure(; size = (320 * usedcols, 260 * nrows + 110))

    ## The dodge slot is taken over every role the table carries, not over
    ## the roles one panel happens to draw, so a role keeps its position
    ## across panels and the note under the legend reads true of all of
    ## them. A panel missing a role leaves its slot empty.
    drawn_roles = [
        role for role in role_order
            if !isempty(select_fit_role(scores, role))
    ]
    dodge_of = Dict(role => i for (i, role) in enumerate(drawn_roles))
    for (i, s) in enumerate(streams)
        r, c = fldmod1(i, usedcols)
        cell = scores[scores.stream .== s, :]
        ax = Axis(
            fig[r, c]; title = string(s),
            xlabel = r == nrows ? "Forecast horizon (days)" : "",
            ylabel = c == 1 ? ylabel : "",
            xticks = (
                Float64.(1:length(horizons)),
                [string(h) for h in horizons],
            )
        )
        ## One stacked bar per (horizon, role), for the roles this stream
        ## carries.
        roles = [
            role for role in drawn_roles
                if !isempty(select_fit_role(cell, role))
        ]
        isempty(roles) && continue
        xs = Float64[]
        ys = Float64[]
        stack = Int[]
        dodge = Int[]
        colours = Symbol[]
        for role in roles
            di = dodge_of[role]
            rs = select_fit_role(cell, role)
            for row in eachrow(rs)
                for (pi, (col, _, colour)) in enumerate(_CRPS_PARTS)
                    v = row[col]
                    (ismissing(v) || !isfinite(v) || v <= 0) && continue
                    push!(xs, Float64(hslot[row.horizon]))
                    push!(ys, Float64(v))
                    push!(stack, pi)
                    push!(dodge, di)
                    push!(colours, colour)
                end
            end
        end
        isempty(xs) && continue
        CairoMakie.barplot!(
            ax, xs, ys; stack = stack, dodge = dodge,
            color = colours, n_dodge = length(drawn_roles),
            gap = 0.25, dodge_gap = 0.06
        )
    end

    handles = Any[
        CairoMakie.PolyElement(; color = colour)
            for (_, _, colour) in _CRPS_PARTS
    ]
    labels = String[name for (_, name, _) in _CRPS_PARTS]
    isempty(title) || CairoMakie.Label(
        fig[0, 1:usedcols], title;
        font = :bold, tellwidth = false
    )
    CairoMakie.Legend(
        fig[nrows + 1, 1:usedcols], handles, labels;
        orientation = :horizontal, framevisible = true,
        tellheight = true, tellwidth = false
    )
    ## The roles share the component colours, so a panel drawing more than
    ## one names them by the order its bars are dodged in rather than by
    ## colour. A note under the legend rather than an entry in it, since an
    ## entry that long crowds the colour keys off a narrow figure.
    length(drawn_roles) > 1 && CairoMakie.Label(
        fig[nrows + 2, 1:usedcols],
        "bars at each horizon, left to right: " * join(drawn_roles, ", ");
        tellwidth = false, color = (:black, 0.6)
    )
    return fig
end

"""
Horizontal point-and-interval comparison of cumulative-case estimates from
several sources. `rows` is a vector of `(label, central, lower, upper)`
tuples, drawn top to bottom with the central estimate as a point and
`[lower, upper]` as a bar. A row whose lower and upper match its central is a
deterministic point estimate and is drawn as a bare marker with no bar.

`groups` is an optional vector of group keys, one per row, matched against
`group_colours` (a vector of `key => colour` pairs) to colour each row's
marker and bar and build a legend. Without `groups` the rows share a single
colour.
"""
function plot_estimate_comparison(
        rows::AbstractVector;
        xlabel::AbstractString = "Cumulative cases",
        xmax::Union{Nothing, Real} = nothing,
        groups::Union{Nothing, AbstractVector} = nothing,
        group_colours::AbstractVector = Pair[]
    )
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
    ax = Axis(
        fig[1, 1];
        xlabel = xlabel,
        yticks = (collect(1:n), reverse(labels)),
        limits = ((0, top), (0.5, n + 0.5))
    )
    for i in 1:n
        y = n - i + 1
        col = _colour(i)
        ## A deterministic point estimate has no interval, so it gets a bare
        ## marker.
        if hi[i] > lo[i]
            lines!(
                ax, [lo[i], hi[i]], [y, y];
                color = (col, 0.8), linewidth = 3
            )
        end
        scatter!(ax, [central[i]], [y]; color = col, markersize = 12)
    end
    if !isnothing(groups) && !isempty(group_colours)
        handles = [
            CairoMakie.MarkerElement(;
                color = c, marker = :circle,
                markersize = 12
            ) for (_, c) in group_colours
        ]
        glabels = [String(k) for (k, _) in group_colours]
        CairoMakie.axislegend(
            ax, handles, glabels; position = :rb,
            framevisible = true
        )
    end
    return fig
end

"""
Calendar time-series comparison of cumulative-count projections against the
data observed since. `external` is another group's published projection and
`ours` is our own forward projection, each drawn as a central line with a
shaded `[lower, upper]` band. `observed` is the data observed so far, drawn
as a marked line.

Each is a vector of `(date, ...)` tuples with `date` an ISO string.
`external` and `ours` are `(date, central, lower, upper)`, `observed` is
`(date, value)`. The dates share one calendar x-axis, so the projections read
against what the outbreak did. Used to set our forward projection beside the
[chamla2026](@citet) confirmed-case projection.
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
        title::AbstractString = "Projected versus observed cumulative cases"
    )
    _x(d) = Float64(date2epochdays(Date(String(d))))

    fig = Figure(; size = (820, 460))
    ax = Axis(fig[1, 1]; xlabel = "Date", ylabel = ylabel, title = title)

    ## External projection: shaded 90% band, central line and markers.
    ex_x = [_x(r[1]) for r in external]
    ex_m = [float(r[2]) for r in external]
    ex_lo = [float(r[3]) for r in external]
    ex_hi = [float(r[4]) for r in external]
    ord = sortperm(ex_x)
    band!(
        ax, ex_x[ord], ex_lo[ord], ex_hi[ord];
        color = (external_colour, 0.15)
    )
    lines!(ax, ex_x[ord], ex_m[ord]; color = external_colour, linewidth = 2)
    ex_h = scatter!(ax, ex_x, ex_m; color = external_colour, markersize = 8)

    ## Observed data so far: a marked black line.
    ob_x = [_x(r[1]) for r in observed]
    ob_y = [float(r[2]) for r in observed]
    obord = sortperm(ob_x)
    lines!(ax, ob_x[obord], ob_y[obord]; color = :black, linewidth = 1.5)
    ob_h = scatter!(ax, ob_x, ob_y; color = :black, markersize = 7)

    ## Our forward projection: shaded 90% band, central line and markers, in
    ## the same ribbon form as the external projection.
    our_x = [_x(r[1]) for r in ours]
    our_m = [float(r[2]) for r in ours]
    our_lo = [float(r[3]) for r in ours]
    our_hi = [float(r[4]) for r in ours]
    oord = sortperm(our_x)
    band!(
        ax, our_x[oord], our_lo[oord], our_hi[oord];
        color = (ours_colour, 0.15)
    )
    lines!(ax, our_x[oord], our_m[oord]; color = ours_colour, linewidth = 2)
    our_h = scatter!(
        ax, our_x, our_m;
        color = ours_colour, markersize = 8, marker = :diamond
    )

    allx = vcat(ex_x, ob_x, our_x)
    lo, hi = minimum(allx), maximum(allx)
    ax.xticks = collect(lo:14:hi)
    ax.xtickformat = vals -> [
        string(epochdays2date(round(Int, v)))
            for v in vals
    ]
    CairoMakie.axislegend(
        ax, [ex_h, our_h, ob_h],
        [external_label, ours_label, observed_label];
        position = :lt, framevisible = true
    )
    return fig
end

"""
Faceted point-and-interval comparison of published scenario estimates, one
panel per date of estimation. `scenarios` is `REPORT_SCENARIOS_CI`-shaped
`(date, label, central, lower, upper)` with `label` of the form
`"M1|M2 <family>, <swept level>"`, for example `"M2 τ=14 d, CFR 26%"`.

Within a panel each `(method, family)` is one row, and the swept nuisance
level (the CFR, window or doubling time) is dodged onto that line so every
scenario keeps its own interval while the sweep adds no rows. Method sets the
colour. `ours` maps a date string to our matched `(median, lower, upper)`
estimate, drawn as a grey reference band with a dashed median in that date's
panel. `date_titles` are `date => title` pairs giving each panel its
heading.
"""
function plot_scenario_comparison(
        scenarios::AbstractVector;
        ours::AbstractDict = Dict{String, Any}(),
        date_titles::AbstractVector = [
            "2026-05-18" => "18 May report",
            "2026-05-20" => "20 May update",
            "2026-05-27" => "27 May (Lancet)",
        ],
        method_names = Dict("M1" => "geographic", "M2" => "back-calc"),
        method_colours = Dict("M1" => :steelblue, "M2" => :darkorange),
        xlabel::AbstractString = "Cumulative cases",
        title::AbstractString = "Published scenarios versus our estimate"
    )
    title_of = Dict(date_titles)
    dates = sort(unique(String[String(s[1]) for s in scenarios]))

    ## Parse "M2 τ=14 d, CFR 26%" → (method, family, level).
    function parts(label)
        head, level = split(String(label), ", "; limit = 2)
        method, family = split(head, " "; limit = 2)
        return (String(method), String(family), String(level))
    end

    ## Group each date's scenarios into ordered (method, family) rows, so a
    ## shared row count keeps the panels' rows the same height.
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
    ## blocks line up across dates.
    countm(d, m) = count(f -> f[1] == m, first(grouped[d]))
    maxgeo = maximum(countm(d, "M1") for d in dates)
    maxbc = maximum(countm(d, "M2") for d in dates)
    maxrow = maxgeo + maxbc

    xmax = 1.05 * max(
        maximum(float(s[5]) for s in scenarios),
        maximum((float(v[3]) for v in values(ours)); init = 0.0)
    )

    fig = Figure(; size = (340 * length(dates) + 140, 110 + 70 * maxrow))
    for (j, d) in enumerate(dates)
        fams, members = grouped[d]
        geo = [f for f in fams if f[1] == "M1"]
        bc = [f for f in fams if f[1] == "M2"]
        ## Geographic families fill the top block (`maxrow` down). Back-calc
        ## families fill the bottom block (`maxbc` down), so both align in y.
        ypos = Dict{Tuple{String, String}, Int}()
        for (i, f) in enumerate(geo)
            ypos[f] = maxrow - i + 1
        end
        for (i, f) in enumerate(bc)
            ypos[f] = maxbc - i + 1
        end
        yvals = [ypos[f] for f in fams]
        ylabels = [
            string(get(method_names, m, m), " · ", fam)
                for (m, fam) in fams
        ]
        ax = Axis(
            fig[1, j];
            title = get(title_of, d, d), xlabel = xlabel,
            yticks = (yvals, ylabels),
            limits = ((0, xmax), (0.4, maxrow + 0.6))
        )

        ## Our matched estimate for this vintage: a reference band + median.
        if haskey(ours, d)
            med, lo, hi = ours[d]
            vspan!(ax, float(lo), float(hi); color = (:grey, 0.18))
            vlines!(
                ax, [float(med)];
                color = :black, linestyle = :dash, linewidth = 1.5
            )
        end

        ## Each family row carries its swept levels dodged around the row
        ## centre, each a point with its interval, coloured by method.
        for key in fams
            y = ypos[key]
            col = get(method_colours, key[1], :grey)
            ms = members[key]
            k = length(ms)
            offs = k == 1 ? [0.0] : collect(LinRange(0.26, -0.26, k))
            for (t, (_, c, lo, hi)) in enumerate(ms)
                yy = y + offs[t]
                lines!(
                    ax, [lo, hi], [yy, yy];
                    color = (col, 0.85), linewidth = 2.5
                )
                scatter!(ax, [c], [yy]; color = col, markersize = 9)
            end
        end
    end

    CairoMakie.Label(
        fig[0, 1:length(dates)], title;
        fontsize = 16, font = :bold
    )
    handles = [
        CairoMakie.MarkerElement(;
            color = method_colours["M1"],
            marker = :circle, markersize = 11
        ),
        CairoMakie.MarkerElement(;
            color = method_colours["M2"],
            marker = :circle, markersize = 11
        ),
        CairoMakie.PolyElement(; color = (:grey, 0.4)),
    ]
    labels = [
        get(method_names, "M1", "M1") * " spread",
        get(method_names, "M2", "M2") * " spread", "our estimate (90%)",
    ]
    CairoMakie.Legend(
        fig[2, 1:length(dates)], handles, labels;
        orientation = :horizontal, framevisible = false
    )
    return fig
end

"""
Density of a prior over the case-fatality ratio (CFR) on `[0, 1]`, plotted
on the sub-range `[0, 0.7]`. The CDC central estimate of 55/169 ≈ 0.33 is
drawn as a solid vertical rule and the report's 26% and 40% scenario bounds
as dashed rules, so the prior reads against the published CFR scenarios.
"""
function plot_cfr_prior(prior::Distribution)
    colours = CairoMakie.Makie.wong_colors()
    xs = range(0.0, 0.7; length = 400)
    ys = pdf.(Ref(prior), xs)

    fig = Figure(; size = (760, 420))
    ax = Axis(
        fig[1, 1];
        xlabel = "Case-fatality ratio (CFR)",
        ylabel = "Prior density",
        title = "Prior over the case-fatality ratio",
        limits = ((0, 0.7), nothing)
    )
    lines!(ax, xs, ys; color = colours[1], linewidth = 2)
    vlines!(ax, [55 / 169]; color = :firebrick, linewidth = 2)
    vlines!(
        ax, [0.26, 0.4];
        color = (:grey, 0.6), linestyle = :dash, linewidth = 2
    )
    return fig
end

"""
Posterior densities of the delay-corrected confirmed case-fatality ratio and
the structural (infection-based) CFR from a
[`delay_corrected_confirmed_cfr`](@ref) result `res`, on the CFR percentage
scale. The naive observed confirmed ratio is a solid vertical rule and the
median uncorrected modelled confirmed ratio a dashed one. The gap from the
naive rule to the corrected density is the real-time delay debiasing, and the
gap to the structural density the residual case/death ascertainment
difference.
"""
function plot_confirmed_cfr(res)
    colours = CairoMakie.Makie.wong_colors()
    corrected = 100 .* filter(isfinite, res.corrected)
    structural = 100 .* filter(isfinite, res.structural)
    naive = 100 * res.naive_observed
    modelled_naive = 100 * quantile(filter(isfinite, res.modelled_naive), 0.5)

    hi = max(maximum(corrected), maximum(structural), naive) * 1.05
    fig = Figure(; size = (760, 420))
    ax = Axis(
        fig[1, 1];
        xlabel = "Case-fatality ratio (%)",
        ylabel = "Posterior density",
        title = "Delay-corrected confirmed CFR versus the structural CFR",
        limits = ((0, hi), nothing)
    )
    h_corr = density!(
        ax, corrected; color = (colours[1], 0.5),
        strokecolor = colours[1], strokewidth = 2
    )
    h_struct = density!(
        ax, structural; color = (colours[2], 0.4),
        strokecolor = colours[2], strokewidth = 2
    )
    h_naive = vlines!(ax, [naive]; color = :firebrick, linewidth = 2)
    h_mod = vlines!(
        ax, [modelled_naive];
        color = (:grey, 0.7), linestyle = :dash, linewidth = 2
    )
    CairoMakie.axislegend(
        ax,
        [h_corr, h_struct, h_naive, h_mod],
        [
            "Delay-corrected confirmed CFR", "Structural CFR",
            "Naive observed confirmed ratio",
            "Uncorrected modelled confirmed ratio (median)",
        ];
        position = :rt, framevisible = true
    )
    return fig
end

"""
One-row, two-panel figure summarising when the outbreak began. The left
panel is the posterior density of the outbreak start date, the calendar date
of the import that started the outbreak, from the outbreak age `T` as
`as_of_date` minus `T`. The right panel is the joint `(doubling_time, T)`
posterior pair plot. Shorter doubling times mean faster early growth, which
reaches the same epidemic size in less time.
"""
function plot_start_date_pair(
        chn;
        as_of_date::AbstractString, thin::Integer = 2
    )
    T_draws = _draws(chn, :T)
    cutoff_days = date2epochdays(Date(as_of_date))
    start_days = cutoff_days .- T_draws

    fig = Figure(; size = (1100, 460))
    ax = Axis(
        fig[1, 1];
        xlabel = "Outbreak start date",
        ylabel = "Posterior density",
        title = "Estimated outbreak start date",
        xticklabelrotation = π / 6
    )
    density!(
        ax, start_days; color = (:steelblue, 0.5),
        strokecolor = :steelblue, strokewidth = 2
    )
    ## Date ticks every four weeks across the posterior range, so the axis
    ## does not crowd as the range widens.
    lo = floor(Int, minimum(start_days))
    hi = ceil(Int, maximum(start_days))
    ax.xticks = collect(lo:28:hi)
    ax.xtickformat = vals -> [string(epochdays2date(round(Int, v))) for v in vals]

    dt_draws = _draws(chn, :doubling_time)
    ## Clip extreme doubling times (near-zero growth) to keep the pair plot
    ## readable. A cap at 200 days covers the credible range.
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
Each day is then scaled by the draw's susceptible fraction at the end of
the day before, the chain's `susceptible_fraction`, so the trajectory is
net of depletion ([`adjusted_rt`](@ref)). A chain without that series
predates depletion and is returned unscaled.
Shared by [`plot_rt`](@ref) and [`plot_rt_streams`](@ref).
"""
function reconstruct_rt(
        chn; n::Integer, breakpoint::Real,
        rt_start::Integer = 1, rt_walk_start::Integer = rt_start,
        week::Integer = 7, ramp::Real = RT_INTERVENTION_RAMP
    )
    rt = _reconstruct_rt_walk(
        chn; n, breakpoint, rt_start, rt_walk_start, week, ramp
    )
    _has_key(chn, :susceptible_fraction) || return rt
    fractions = _draw_vectors(chn, :susceptible_fraction)
    for i in axes(rt, 1)
        _deplete_rt!(view(rt, i, :), fractions[i])
    end
    return rt
end

## Scale an established-window Rt row by the susceptible fraction at the end
## of each day before, as `adjusted_rt` does in the model.
function _deplete_rt!(row::AbstractVector, fraction::AbstractVector)
    for d in 2:length(row)
        ismissing(row[d]) && continue
        row[d] *= fraction[d - 1]
    end
    return row
end

## The walk's daily Rt per draw, before depletion.
function _reconstruct_rt_walk(
        chn; n::Integer, breakpoint::Real,
        rt_start::Integer = 1, rt_walk_start::Integer = rt_start,
        week::Integer = 7, ramp::Real = RT_INTERVENTION_RAMP
    )
    log_R0 = _draws(chn, Symbol("rt_state.log_R0"))
    sigma = _draws(chn, Symbol("rt_state.sigma_rw"))
    effect = _draws(chn, Symbol("rt_state.intervention_effect"))
    ## `rt_state.z` is vector-valued: one standard-normal innovation vector
    ## per draw. Pull each draw's full vector from the chain slice.
    zmat = chn[Symbol("rt_state.z")]
    zrows = [collect(z) for z in vec(collect(zmat))]

    ## The knot grid is built from the model's walk start `rt_walk_start`
    ## (the breakpoint grid day), which is decoupled from `rt_start` (the
    ## established-window start used for the mask below). The innovation
    ## vector length is fixed by that walk start, so a mismatching one fails
    ## here rather than as a downstream bounds error.
    days = knot_days(n; week, start = rt_walk_start)
    nb = length(days)
    if !isempty(zrows) && length(zrows[1]) != nb - 1
        error(
            "reconstruct_rt: rt_walk_start = $rt_walk_start gives " *
                "$(nb - 1) random-walk steps but the chain has " *
                "$(length(zrows[1])); pass the same walk start the model used " *
                "(the breakpoint grid day, n - who_first_sitrep_days)."
        )
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
        ## Days before the renewal start clamp to the established R0, the
        ## walk base. The model fills them with the analytic cryptic
        ## exponential and they are masked out here.
        log_Rt = walk .+ effect[i] .* ramp_shape
        start = clamp(rt_start, 1, n)
        for d in start:n
            rt[i, d] = exp(log_Rt[d])
        end
    end
    return rt
end

"""
    fitted_onset_hazard(model, chn) -> NamedTuple

Every posterior draw's fitted symptom-onset reporting hazard, read off the
fitted `model`'s own onset-reporting state at each draw of `chn`
(`returned`): `logit_h0` (the baseline delay hazard), `γ` (the report-date
calendar walk) and `alpha` (the ascertainment level over the triangle's
onset dates), one vector per draw, as the composers return them. The
model is the one the chain was fitted with, such as a fit spec's `model()`
in the report.
"""
function fitted_onset_hazard(model::Model, chn)
    states = [r.onset_report_state for r in vec(returned(model, chn))]
    return (;
        logit_h0 = [collect(Float64, st.logit_h0) for st in states],
        γ = [collect(Float64, st.γ) for st in states],
        alpha = [collect(Float64, st.alpha) for st in states],
    )
end

"""
    onset_nowcast_draws(days, observed, delays, onsets, hazard;
                        grid_start, target_delays)

[`onset_nowcast`](@ref) per posterior draw, one `ndraws`-long vector per
onset day in `days`. `observed[k]` is the count a digitised figure prints
for `days[k]` and `delays[k]` is that figure's own reporting delay, so a
snapshot's cells are nowcast from the delay that snapshot had run to.

`target_delays` is the delay to nowcast to, one per day; the default
`nothing` targets each day's eventual total. Pass the delay of the figure
the prediction will be compared against to keep the two like for like.

`onsets` holds each draw's daily onsets indexed by grid day, the `diff` of
the chain's `cumulative_onsets`. `hazard` is
[`fitted_onset_hazard`](@ref)'s `(; logit_h0, γ, alpha)`, with `alpha`
indexed from `grid_start` and held flat outside the fitted grid. The two are
paired draw by draw and must come from one fit. Summarised by
[`plot_onset_nowcast_grid`](@ref).
"""
function onset_nowcast_draws(
        days::AbstractVector{<:Integer},
        observed::AbstractVector{<:Real},
        delays::AbstractVector{<:Integer},
        onsets::AbstractVector{<:AbstractVector{<:Real}},
        hazard::NamedTuple; grid_start::Integer,
        target_delays::Union{Nothing, AbstractVector{<:Integer}} = nothing
    )
    n = length(days)
    if length(observed) != n || length(delays) != n ||
            (!isnothing(target_delays) && length(target_delays) != n)
        error(
            "onset_nowcast_draws: `days`, `observed`, `delays` and any " *
                "`target_delays` must have the same length, got $n, " *
                "$(length(observed)), $(length(delays)) and " *
                "$(isnothing(target_delays) ? "none" : length(target_delays))."
        )
    end
    nd = length(onsets)
    if length(hazard.alpha) != nd || length(hazard.logit_h0) != nd ||
            length(hazard.γ) != nd
        error(
            "onset_nowcast_draws: `onsets` and `hazard` must come from " *
                "the same fit, got $nd onset draws against " *
                "$(length(hazard.alpha)) ascertainment, " *
                "$(length(hazard.logit_h0)) baseline-hazard and " *
                "$(length(hazard.γ)) calendar-walk draws."
        )
    end
    ## The onset series is indexed by grid day, so name an out-of-range day
    ## rather than raising a bare `BoundsError` inside the draw loop.
    ndays = nd == 0 ? 0 : length(first(onsets))
    for d in days
        (1 <= d <= ndays) ||
            error(
            "onset_nowcast_draws: day $d is outside the onset " *
                "series, which runs 1:$ndays."
        )
    end
    out = Vector{Vector{Float64}}(undef, n)
    for k in 1:n
        u = Int(days[k])
        δ = Int(delays[k])
        y = float(observed[k])
        until = isnothing(target_delays) ? nothing : Int(target_delays[k])
        out[k] = [
            begin
                a = hazard.alpha[i]
                α = a[clamp(u - Int(grid_start) + 1, 1, length(a))]
                onset_nowcast(
                    y, onsets[i][u], δ, hazard.logit_h0[i],
                    hazard.γ[i], u, grid_start, α; until
                )
            end
                for i in eachindex(onsets)
        ]
    end
    return out
end

"""
    plot_onset_nowcast_grid(panels; kwargs...)

Nowcast of the symptom-onset reporting triangle, one panel per digitised
snapshot: what that snapshot's own figure implied for the onset dates it
printed, against what the figures print for them now.

Each `panel` is a `NamedTuple` of `title` (the snapshot's report date),
`dates` (the onset dates, one per x position), `observed` (that snapshot's
own counts, grey crosses), `nowcast` (per-draw predictions of the count the
latest figure prints, drawn as 30/60/90% ribbons with a median line) and
`latest` (that count, black points).

The panel is read on whether the ribbon covers the black points. Pass
`nowcast` as a predictive rather than the latent count. Two scans of one bar
disagree by the read error the stream estimates, so on the onset dates where
reporting had already finished the latent quantity is the grey cross exactly
and would be scored against a reading it cannot match. Build it from
[`onset_nowcast_draws`](@ref) at the latest figure's own delay, then through
the stream's bar measurement error.

A panel whose series disagree in length raises, and an empty `panels`
returns a blank figure.
"""
function plot_onset_nowcast_grid(
        panels::AbstractVector;
        ncol::Integer = 3, colour = :steelblue,
        title = "Symptom-onset reporting triangle: nowcast vs digitised"
    )
    isempty(panels) && return Figure()
    ncols = min(length(panels), Int(ncol))
    nrows = cld(length(panels), ncols)
    fig = Figure(; size = (360 * ncols, 280 * nrows + 60))
    for (j, p) in enumerate(panels)
        n = length(p.dates)
        if length(p.observed) != n || length(p.nowcast) != n ||
                length(p.latest) != n
            error(
                "plot_onset_nowcast_grid: panel $(p.title) must carry " *
                    "one `observed`, `nowcast` and `latest` entry per onset " *
                    "date, got $(length(p.observed)), $(length(p.nowcast)) " *
                    "and $(length(p.latest)) for $n dates."
            )
        end
        row, col = fldmod1(j, ncols)
        x = collect(1:n)
        ## Sorted once per cell, else the seven ribbon and median quantiles
        ## below sort every draw vector seven times.
        sorted = [sort(float.(d)) for d in p.nowcast]
        q(pr) = [quantile(d, pr; sorted = true) for d in sorted]
        hi60 = q(0.8)
        ## The 90% tail on the newest onset dates runs well past the counts
        ## the panel is read on, so the axis is set by the 60% ribbon.
        yupper = 1.6 * max(
            1.0,
            isempty(p.latest) ? 1.0 : maximum(float.(p.latest)),
            isempty(hi60) ? 1.0 : maximum(hi60)
        )
        ax = Axis(
            fig[row, col]; title = string(p.title),
            xlabel = row == nrows ? "onset date" : "",
            ylabel = col == 1 ? "cases at this onset date" : "",
            xticks = _vintage_ticks(p.dates),
            xticklabelrotation = pi / 4, xticklabelsize = 11,
            limits = (nothing, (0, yupper))
        )
        band!(ax, x, q(0.05), q(0.95); color = (colour, 0.15))
        band!(ax, x, q(0.2), hi60; color = (colour, 0.28))
        band!(ax, x, q(0.35), q(0.65); color = (colour, 0.42))
        lines!(ax, x, q(0.5); color = colour, linewidth = 2)
        scatter!(
            ax, x, float.(p.observed); color = (:grey40, 0.9),
            marker = :cross, markersize = 8
        )
        scatter!(ax, x, float.(p.latest); color = :black, markersize = 6)
    end
    CairoMakie.Label(
        fig[0, 1:ncols], title; font = :bold,
        tellwidth = false
    )
    CairoMakie.Legend(
        fig[nrows + 1, 1:ncols],
        [
            CairoMakie.MarkerElement(color = (:grey40, 0.9), marker = :cross),
            CairoMakie.LineElement(color = colour),
            CairoMakie.PolyElement(color = (colour, 0.28)),
            CairoMakie.MarkerElement(color = :black, marker = :circle),
        ],
        [
            "digitised at this snapshot", "nowcast median",
            "nowcast 30/60/90%", "digitised to date",
        ];
        orientation = :horizontal, tellwidth = false
    )
    return fig
end

## Per-day quantile `pr` of an established-window Rt matrix, skipping the
## masked (pre-renewal) days: `missing` where a day has no established draws.
function _rt_quantile(rt::AbstractMatrix, d::Integer, pr::Real)
    col = collect(skipmissing(@view rt[:, d]))
    return isempty(col) ? missing : quantile(col, pr)
end

"""
Daily reproduction-number trajectory `Rt` per posterior draw, rebuilt by
[`reconstruct_rt`](@ref) and plotted over the established-outbreak window.

That window runs from `rt_start`, the renewal start where the random walk
begins, to the cut-off. Only that period is drawn, as 30%, 60% and 90%
credible ribbons with no median line, with about `n_traj` thinned sampled
trajectories overlaid faint to show the per-draw spread. The intervention
breakpoint, the end of the scale-up (`breakpoint + ramp`, dotted) and the
cut-off are marked. `seeding` is the calendar date of grid day 1, so day `d`
is `seeding + (d - 1)`.
"""
function plot_rt(
        chn; n::Integer, breakpoint::Real,
        as_of_date::AbstractString, seeding::Date,
        rt_start::Integer = 1, rt_walk_start::Integer = rt_start,
        week::Integer = 7, ramp::Real = RT_INTERVENTION_RAMP,
        n_traj::Integer = 100
    )
    rt = reconstruct_rt(chn; n, breakpoint, rt_start, rt_walk_start, week, ramp)
    ndraws = size(rt, 1)

    ## Median and ribbons over established draws only (skip masked days).
    q(d, pr) = _rt_quantile(rt, d, pr)
    med = [q(d, 0.5) for d in 1:n]
    lo90 = [q(d, 0.05) for d in 1:n]
    hi90 = [q(d, 0.95) for d in 1:n]
    lo60 = [q(d, 0.2) for d in 1:n]
    hi60 = [q(d, 0.8) for d in 1:n]
    est = findall(!ismissing, med)

    lo30 = [q(d, 0.35) for d in 1:n]
    hi30 = [q(d, 0.65) for d in 1:n]

    epoch = date2epochdays(seeding)
    x = [epoch + (d - 1) for d in 1:n]
    xe = x[est]
    ## Cap the y-axis just above the upper 90% credible band, with 20%
    ## headroom. The handful of sampled trajectories above the cap are
    ## clipped.
    hi90_est = Float64[hi90[d] for d in est if !ismissing(hi90[d])]
    ytop = isempty(hi90_est) ? 4.0 :
        max(1.2, 1.2 * maximum(hi90_est))
    fig = Figure(; size = (900, 440))
    ax = Axis(
        fig[1, 1]; xlabel = "Date", ylabel = "Reproduction number Rt",
        title = "Estimated Rt over the established outbreak",
        limits = (nothing, (0.0, ytop)),
        xticklabelrotation = pi / 6
    )
    ## Thin sampled trajectories, faint, so the per-draw spread reads
    ## alongside the ribbons.
    if n_traj > 0 && !isempty(est)
        step = max(1, fld(ndraws, n_traj))
        for i in 1:step:ndraws
            yi = Float64[rt[i, d] for d in est]
            lines!(ax, xe, yi; color = (:purple, 0.15), linewidth = 0.5)
        end
    end
    band!(
        ax, xe, Float64[lo90[d] for d in est], Float64[hi90[d] for d in est];
        color = (:purple, 0.15)
    )
    band!(
        ax, xe, Float64[lo60[d] for d in est], Float64[hi60[d] for d in est];
        color = (:purple, 0.28)
    )
    band!(
        ax, xe, Float64[lo30[d] for d in est], Float64[hi30[d] for d in est];
        color = (:purple, 0.42)
    )
    ## No-growth threshold at Rt = 1.
    hlines!(ax, [1.0]; color = (:grey, 0.8), linestyle = :dash, linewidth = 2)
    vlines!(
        ax, [Float64(epoch + breakpoint - 1)];
        color = :firebrick, linestyle = :dash, linewidth = 2
    )
    ## End of the intervention scale-up.
    vlines!(
        ax, [Float64(epoch + breakpoint - 1 + ramp)];
        color = :firebrick, linestyle = :dot, linewidth = 2
    )
    vlines!(
        ax, [Float64(date2epochdays(Date(as_of_date)))];
        color = :grey, linestyle = :dash
    )
    ## Limit the x-axis to the estimated window.
    lo = isempty(xe) ? floor(Int, minimum(x)) : floor(Int, minimum(xe))
    hi = ceil(Int, maximum(x))
    CairoMakie.xlims!(ax, lo, hi)
    CairoMakie.ylims!(ax, 0, ytop)
    ax.xticks = collect(lo:7:hi)
    ax.xtickformat = vals -> [
        string(epochdays2date(round(Int, v)))
            for v in vals
    ]
    return fig
end

## Per-fit credible bands over the shared display window (`ds` onward), with
## the fit's own established mask applied through `reconstruct_rt`. Returns a
## NamedTuple of the 30/60/90% lower/upper quantiles and the established days
## `est` to draw.
function _rt_bands(chn; n, breakpoint, rt_start, rt_walk_start, week, ramp, ds)
    rt = reconstruct_rt(chn; n, breakpoint, rt_start, rt_walk_start, week, ramp)
    return _rt_bands_matrix(rt; n, ds)
end

## The band quantiles for an already-reconstructed `ndraws × n` Rt matrix.
## Split out from `_rt_bands` so the per-province trajectories, built from the
## national walk plus a deviation, summarise through the same code.
function _rt_bands_matrix(rt::AbstractMatrix; n, ds)
    q(pr) = [_rt_quantile(rt, d, pr) for d in 1:n]
    med = q(0.5)
    est = findall(d -> d >= ds && !ismissing(med[d]), 1:n)
    return (;
        lo90 = q(0.05), hi90 = q(0.95), lo60 = q(0.2), hi60 = q(0.8),
        lo30 = q(0.35), hi30 = q(0.65), est,
    )
end

## Draw nested 30/60/90% credible ribbons (no median line) for one fit's
## bands in `colour`.
function _draw_rt_bands!(
        ax, x, b, colour;
        alphas = (0.15, 0.28, 0.42)
    )
    isempty(b.est) && return
    xe = x[b.est]
    band!(
        ax, xe, Float64[b.lo90[d] for d in b.est],
        Float64[b.hi90[d] for d in b.est]; color = (colour, alphas[1])
    )
    band!(
        ax, xe, Float64[b.lo60[d] for d in b.est],
        Float64[b.hi60[d] for d in b.est]; color = (colour, alphas[2])
    )
    band!(
        ax, xe, Float64[b.lo30[d] for d in b.est],
        Float64[b.hi30[d] for d in b.est]; color = (colour, alphas[3])
    )
    return
end

"""
Faceted implied reproduction number, one panel per single-stream fit, each
with the joint fit overlaid as the reference. Each fit's daily `Rt` comes
from [`reconstruct_rt`](@ref) on its own sampled random walk, so the figure
shows what reproduction number each data stream implies on its own against
the all-streams-together joint estimate.

Every panel draws 30/60/90% credible ribbons with no median line, the joint
in grey first and then the stream in its colour on top, with the panel title
in that colour. The y-axis is shared across panels and capped from the
panels' 90% bands, so a weakly-informed stream does not stretch the scale.
Ribbon above the cap is clipped.

Each `stream` is a `NamedTuple`
`(; label, chn, rt_start, rt_walk_start, colour)` and `joint` is
`(; label, chn, rt_start, rt_walk_start)`, where `chn` is that fit's chain
and `rt_start`/`rt_walk_start` are the renewal start and random-walk start
that fit used. The per-stream fits walk from day 1, the joint from the
breakpoint lead. `display_start` is the shared grid day the panels draw from,
so every stream reads over the same window. `seeding` is the calendar date of
grid day 1, so day `d` is `seeding + (d - 1)`. The intervention breakpoint,
the end of the scale-up (`breakpoint + ramp`, dotted) and the cut-off are
marked as in [`plot_rt`](@ref).

`title`, `reference_label` and `panel_label` name what the figure is showing.
They default to the per-stream reading, and the sensitivity page passes its
own to set two model structures against each other instead.
"""
function plot_rt_streams(
        streams::AbstractVector;
        joint, n::Integer, breakpoint::Real,
        as_of_date::AbstractString, seeding::Date,
        display_start::Integer = 1, week::Integer = 7,
        ramp::Real = RT_INTERVENTION_RAMP,
        ncols::Integer = 2, joint_colour = :grey25,
        title::AbstractString =
            "Implied Rt by data stream, with the joint fit overlaid",
        reference_label::AbstractString = "the joint fit",
        panel_label::AbstractString =
            "the single-stream fit named in the panel title"
    )
    epoch = date2epochdays(seeding)
    x = Float64[epoch + (d - 1) for d in 1:n]
    ds = clamp(display_start, 1, n)

    ## The joint is the shared reference drawn behind every stream, so
    ## reconstruct it once.
    bj = _rt_bands(
        joint.chn; n, breakpoint, rt_start = joint.rt_start,
        rt_walk_start = joint.rt_walk_start, week, ramp, ds
    )
    sbands = Tuple{Any, NamedTuple}[]
    for s in streams
        b = _rt_bands(
            s.chn; n, breakpoint, rt_start = s.rt_start,
            rt_walk_start = s.rt_walk_start, week, ramp, ds
        )
        push!(sbands, (s, b))
    end

    ## Shared y-cap from the panels' typical 90% upper band, a median over
    ## days rather than a maximum, so one spiky day does not stretch the axis
    ## and the Rt = 1 line stays visible.
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
        ax = Axis(
            fig[r, c]; xlabel = "Date", ylabel = "Rt",
            title = s.label, titlecolor = s.colour,
            xticklabelrotation = pi / 6
        )
        _draw_rt_bands!(
            ax, x, bj, joint_colour;
            alphas = (0.1, 0.16, 0.22)
        )
        _draw_rt_bands!(ax, x, b, s.colour)
        hlines!(
            ax, [1.0]; color = (:grey, 0.8), linestyle = :dash,
            linewidth = 2
        )
        vlines!(
            ax, [Float64(epoch + breakpoint - 1)];
            color = :firebrick, linestyle = :dash, linewidth = 2
        )
        vlines!(
            ax, [Float64(epoch + breakpoint - 1 + ramp)];
            color = :firebrick, linestyle = :dot, linewidth = 2
        )
        vlines!(
            ax, [Float64(date2epochdays(Date(as_of_date)))];
            color = :grey, linestyle = :dash
        )
        CairoMakie.xlims!(ax, lo, hi)
        CairoMakie.ylims!(ax, 0, ytop)
        ax.xticks = collect(lo:14:hi)
        ax.xtickformat = vals -> [
            string(epochdays2date(round(Int, v)))
                for v in vals
        ]
    end

    CairoMakie.Label(
        fig[nrows + 1, 1:ncols],
        "Bands are 30/60/90% credible intervals. Grey is " *
            reference_label * ", the same in every panel. The coloured band " *
            "is " * panel_label * ".";
        fontsize = 12, padding = (0, 0, 0, 6)
    )
    CairoMakie.Label(fig[0, 1:ncols], title; fontsize = 16, font = :bold)
    return fig
end

"""
Reconstruct each posterior draw's daily reproduction number for every
province, returning a vector of `ndraws × n` matrices, one per patch, each
masked to the draw's established window exactly as [`reconstruct_rt`](@ref)
masks the national trajectory.

The chain stores the provincial `Rt` only at the cut-off (`R_T_patch`), so
the trajectory is rebuilt by mirroring the model: the central trend from
[`reconstruct_rt`](@ref), times `exp(δ_p(t))` with `δ_p` the sum-to-zero
deviation interpolated from the weekly knots the chain carries as
`delta_knots` ([`interpolate_knots`](@ref)). `delta_knots` is the
`(n_patches × n_knots)` deviation matrix flattened column-major.

Each province runs its own renewal at its own `Rt` and nothing rescales it,
so `μ(t) · exp(δ_p(t))` is what the model used. Scaled by the province's
susceptible fraction the day before (the chain's
`susceptible_fraction_patch`), it is what `R_T_patch` reports. A chain
without that series predates depletion and is returned unscaled.
The national reproduction number is not `μ` but the value implied by the
summed infections, which is why [`plot_rt_patches`](@ref) draws it from the
chain's own national trajectory rather than from these.

A chain carrying no usable deviations is an error here rather than a silent
national trajectory repeated per panel.
"""
function reconstruct_patch_rt(
        chn; n::Integer, breakpoint::Real,
        n_patches::Integer = length(PROVINCE_NAMES),
        rt_start::Integer = 1, rt_walk_start::Integer = rt_start,
        week::Integer = 7, ramp::Real = RT_INTERVENTION_RAMP
    )
    national = _reconstruct_rt_walk(
        chn; n, breakpoint, rt_start, rt_walk_start,
        week, ramp
    )
    days = knot_days(n; week, start = rt_walk_start)
    nb = length(days)
    knots = try
        [collect(v) for v in vec(collect(chn[:delta_knots]))]
    catch
        error(
            "reconstruct_patch_rt: the chain is missing `delta_knots`, " *
                "so the provincial Rt trajectories cannot be rebuilt. It was " *
                "sampled either with `n_patches = 1` or before that was " *
                "surfaced; refit with the patch structure on."
        )
    end
    ndraws = size(national, 1)
    expected = n_patches * nb
    isempty(knots) || length(knots[1]) == expected ||
        error(
        "reconstruct_patch_rt: `delta_knots` holds $(length(knots[1])) " *
            "entries but $n_patches patches by $nb knots is $expected; " *
            "pass the same `n_patches` and `rt_walk_start` the model used."
    )
    out = [
        Matrix{Union{Missing, Float64}}(missing, ndraws, n)
            for _ in 1:n_patches
    ]
    for i in 1:ndraws
        δ_knots = reshape(knots[i], n_patches, nb)
        for p in 1:n_patches
            δ_daily = interpolate_knots(δ_knots[p, :], days, n)
            for d in 1:n
                ismissing(national[i, d]) && continue
                out[p][i, d] = national[i, d] * exp(δ_daily[d])
            end
        end
    end
    _has_key(chn, :susceptible_fraction_patch) || return out
    fractions = _draw_vectors(chn, :susceptible_fraction_patch)
    for i in 1:ndraws
        f = reshape(fractions[i], n_patches, n)
        for p in 1:n_patches
            _deplete_rt!(view(out[p], i, :), view(f, p, :))
        end
    end
    return out
end

"""
Faceted reproduction number by province, one panel per patch, each with the
national trajectory overlaid in grey as the shared reference. Provincial `Rt`
is rebuilt by [`reconstruct_patch_rt`](@ref). The national reference is
[`reconstruct_rt`](@ref), the same trajectory [`plot_rt`](@ref) draws.

Every panel draws 30/60/90% credible ribbons with no median line, on a shared
y-axis so the provinces are compared rather than each rescaled to its own
range. The intervention breakpoint (dashed), the end of the scale-up (dotted)
and the cut-off are marked as in [`plot_rt`](@ref).

The grey reference is the reproduction number implied by the summed
provinces, the incidence-weighted mean of the panels rather than any one
province or the central trend they pool toward. A panel tracking the grey
band says that province moves with the country. Separation between panels is
the spatial signal, and `region_drift_sd`, the sd of each province's
weekly deviation innovation, measures it.
"""
function plot_rt_patches(
        chn; n::Integer, breakpoint::Real,
        as_of_date::AbstractString, seeding::Date,
        n_patches::Integer = length(PROVINCE_NAMES),
        patch_labels::AbstractVector = PROVINCE_LABELS,
        rt_start::Integer = 1, rt_walk_start::Integer = rt_start,
        display_start::Integer = rt_start,
        week::Integer = 7, ramp::Real = RT_INTERVENTION_RAMP,
        ncols::Integer = 3,
        colours = [:firebrick, :steelblue, :seagreen],
        national_colour = :grey25
    )
    np = min(n_patches, length(patch_labels))
    epoch = date2epochdays(seeding)
    x = Float64[epoch + (d - 1) for d in 1:n]
    ds = clamp(display_start, 1, n)

    patch_rt = reconstruct_patch_rt(
        chn; n, breakpoint, n_patches = np,
        rt_start, rt_walk_start, week, ramp
    )
    bands = [_rt_bands_matrix(patch_rt[p]; n, ds) for p in 1:np]
    bn = _rt_bands(chn; n, breakpoint, rt_start, rt_walk_start, week, ramp, ds)

    ## Shared y-cap from the panels' typical 90% upper band, as in
    ## `plot_rt_streams`. A median over days rather than a maximum, so one
    ## spiky day in the weakest-informed province does not flatten the rest.
    function panel_top(b)
        v = Float64[b.hi90[d] for d in b.est if !ismissing(b.hi90[d])]
        return isempty(v) ? 0.0 : quantile(v, 0.5)
    end
    tops = Float64[panel_top(bn); [panel_top(b) for b in bands]]
    ytop = max(2.5, ceil(1.3 * maximum(tops) * 2) / 2)

    lo = floor(Int, x[ds])
    hi = ceil(Int, maximum(x))
    nrows = cld(np, ncols)
    fig = Figure(; size = (480 * min(np, ncols), 320 * nrows + 70))
    for p in 1:np
        r = cld(p, ncols)
        c = p - (r - 1) * ncols
        colour = colours[mod1(p, length(colours))]
        ax = Axis(
            fig[r, c]; xlabel = "Date", ylabel = "Rt",
            title = patch_labels[p], titlecolor = colour,
            xticklabelrotation = pi / 6
        )
        _draw_rt_bands!(
            ax, x, bn, national_colour;
            alphas = (0.1, 0.16, 0.22)
        )
        _draw_rt_bands!(ax, x, bands[p], colour)
        hlines!(
            ax, [1.0]; color = (:grey, 0.8), linestyle = :dash,
            linewidth = 2
        )
        vlines!(
            ax, [Float64(epoch + breakpoint - 1)];
            color = :firebrick, linestyle = :dash, linewidth = 2
        )
        vlines!(
            ax, [Float64(epoch + breakpoint - 1 + ramp)];
            color = :firebrick, linestyle = :dot, linewidth = 2
        )
        vlines!(
            ax, [Float64(date2epochdays(Date(as_of_date)))];
            color = :grey, linestyle = :dash
        )
        CairoMakie.xlims!(ax, lo, hi)
        CairoMakie.ylims!(ax, 0, ytop)
        ax.xticks = collect(lo:14:hi)
        ax.xtickformat = vals -> [
            string(epochdays2date(round(Int, v)))
                for v in vals
        ]
    end
    CairoMakie.Label(
        fig[nrows + 1, 1:min(np, ncols)],
        "Bands are 30/60/90% credible intervals. Grey is the national " *
            "trajectory, the same in every panel; the coloured band is the " *
            "province named in the panel title.";
        fontsize = 12, padding = (0, 0, 0, 6)
    )
    CairoMakie.Label(
        fig[0, 1:min(np, ncols)],
        "Reproduction number by province";
        fontsize = 16, font = :bold
    )
    return fig
end

## Per-patch daily series from a flattened `(n_patches x n)` vector
## deterministic (`infections_patch`, `importation_patch`). The matrix reaches
## the chain flattened column-major, so day `t` of patch `p` sits at
## `(t - 1) * np + p`.
function _patch_daily(chn, sym::Symbol, np::Integer, n::Integer)
    vs = try
        [collect(v) for v in vec(collect(chn[sym]))]
    catch
        error(
            "plot: the chain carries no `$(sym)`, so the per-province " *
                "trajectories cannot be drawn. It was sampled either with " *
                "`n_patches = 1` or before that was surfaced."
        )
    end
    length(first(vs)) == np * n || error(
        "plot: `$(sym)` holds $(length(first(vs))) entries but $np patches " *
            "by $n days is $(np * n)."
    )
    return [
        [Float64[v[(t - 1) * np + p] for t in 1:n] for v in vs]
            for p in 1:np
    ]
end

## 30/60/90% ribbons of a set of per-draw daily series.
function _traj_bands(trajs, n::Integer)
    q(d, pr) = quantile(Float64[t[d] for t in trajs], pr)
    return (
        lo90 = [q(d, 0.05) for d in 1:n], hi90 = [q(d, 0.95) for d in 1:n],
        lo60 = [q(d, 0.2) for d in 1:n], hi60 = [q(d, 0.8) for d in 1:n],
        lo30 = [q(d, 0.35) for d in 1:n], hi30 = [q(d, 0.65) for d in 1:n],
    )
end

function _draw_traj_bands!(ax, x, b, colour)
    band!(ax, x, b.lo90, b.hi90; color = (colour, 0.15))
    band!(ax, x, b.lo60, b.hi60; color = (colour, 0.28))
    band!(ax, x, b.lo30, b.hi30; color = (colour, 0.42))
    return ax
end

"""
Modelled infections by province over time, one column per province, with the
daily infections on the top row and the cumulative total on the bottom. Every
panel is 30/60/90% credible ribbons with no median line.

Each panel carries its own y-axis. The provinces differ by orders of
magnitude, so a shared axis would flatten every province but the epicentre
into the floor. The cross-province comparison belongs in
[`patch_overview_table`](@ref). Reads the `infections_patch` deterministic,
the daily per-province infection matrix flattened column-major.
"""
function plot_infections_patches(
        chn; n::Integer, seeding::Date,
        n_patches::Integer = length(PROVINCE_NAMES),
        patch_labels::AbstractVector = PROVINCE_LABELS,
        colours = [:firebrick, :steelblue, :seagreen]
    )
    np = min(n_patches, length(patch_labels))
    epoch = date2epochdays(seeding)
    x = Float64[epoch + (d - 1) for d in 1:n]
    daily = _patch_daily(chn, :infections_patch, np, n)
    cumul = [[cumsum(t) for t in daily[p]] for p in 1:np]
    lo = floor(Int, minimum(x))
    hi = ceil(Int, maximum(x))
    fig = Figure(; size = (460 * np, 700))
    for p in 1:np
        colour = colours[mod1(p, length(colours))]
        for (r, (trajs, lab)) in enumerate(
                (
                    (daily[p], "Daily infections"),
                    (cumul[p], "Cumulative infections"),
                )
            )
            ax = Axis(
                fig[r, p]; xlabel = "Date", ylabel = lab,
                title = r == 1 ? patch_labels[p] : "",
                titlecolor = colour, xticklabelrotation = pi / 6
            )
            _draw_traj_bands!(ax, x, _traj_bands(trajs, n), colour)
            CairoMakie.xlims!(ax, lo, hi)
            ax.xticks = collect(lo:28:hi)
            ax.xtickformat = vals -> [
                string(epochdays2date(round(Int, v)))
                    for v in vals
            ]
        end
    end
    CairoMakie.Label(
        fig[3, 1:np],
        "Bands are 30/60/90% credible intervals. Each panel has its own " *
            "y-axis, so panels are read for shape and timing rather than " *
            "compared by height.";
        fontsize = 12, padding = (0, 0, 0, 6)
    )
    CairoMakie.Label(
        fig[0, 1:np], "Modelled infections by province";
        fontsize = 16, font = :bold
    )
    return fig
end

"""
Imported infections by province over time, one panel per province: the daily
infections a province received from the others through the importation
kernel, as 30/60/90% credible ribbons.

Every arrival is debited from its origin the same day, so these are
transmission relocated rather than transmission added. The national total
still moves with the coupling, because the destination then grows at its own
reproduction number. The intensity is weakly identified against the secondary
provinces' seeds, since both raise a secondary province's early incidence, so
the level is read as the coupling the data tolerate rather than as a measured
flow. Reads the `importation_patch` deterministic.
"""
function plot_imports_patches(
        chn; n::Integer, seeding::Date,
        n_patches::Integer = length(PROVINCE_NAMES),
        patch_labels::AbstractVector = PROVINCE_LABELS,
        colours = [:firebrick, :steelblue, :seagreen]
    )
    np = min(n_patches, length(patch_labels))
    epoch = date2epochdays(seeding)
    x = Float64[epoch + (d - 1) for d in 1:n]
    imports = _patch_daily(chn, :importation_patch, np, n)
    lo = floor(Int, minimum(x))
    hi = ceil(Int, maximum(x))
    fig = Figure(; size = (460 * np, 400))
    for p in 1:np
        colour = colours[mod1(p, length(colours))]
        ax = Axis(
            fig[1, p]; xlabel = "Date",
            ylabel = "Imported infections per day",
            title = patch_labels[p], titlecolor = colour,
            xticklabelrotation = pi / 6
        )
        _draw_traj_bands!(ax, x, _traj_bands(imports[p], n), colour)
        CairoMakie.xlims!(ax, lo, hi)
        ax.xticks = collect(lo:28:hi)
        ax.xtickformat = vals -> [
            string(epochdays2date(round(Int, v)))
                for v in vals
        ]
    end
    CairoMakie.Label(
        fig[2, 1:np],
        "Bands are 30/60/90% credible intervals. Each panel has its own " *
            "y-axis. Importation relocates transmission between provinces " *
            "rather than adding it.";
        fontsize = 12, padding = (0, 0, 0, 6)
    )
    CairoMakie.Label(
        fig[0, 1:np], "Imported infections by province";
        fontsize = 16, font = :bold
    )
    return fig
end

## Nested 30/60/90% vertical interval bars at one x, topped by a median dot.
## Widths run the other way from the alphas so the nesting reads at a glance:
## the 90% is the thin outer bar and the 30% the thick inner one.
function _draw_patch_interval!(ax, x, draws, colour)
    q = posterior_summary(draws)
    lines!(ax, [x, x], [q.lo90, q.hi90]; color = (colour, 0.35), linewidth = 2)
    lines!(ax, [x, x], [q.lo60, q.hi60]; color = (colour, 0.55), linewidth = 6)
    lines!(
        ax, [x, x], [q.lo30, q.hi30]; color = (colour, 0.85),
        linewidth = 11
    )
    return scatter!(ax, [x], [median(draws)]; color = :black, markersize = 8)
end

"""
Per-province posterior summary as a figure: one panel per quantity, the
provinces side by side on a shared axis, each drawn as a median dot over
nested 30/60/90% credible bars.

This is the figure form of [`patch_summary_table`](@ref) and reads the same
chain deterministics in the same order: the cut-off cumulative infections
`C_T`, the cut-off reproduction number `R_T`, the daily infections at the
cut-off, and the log-Rt deviation `δ` from the common national trend, plus
the contrast against the primary patch, the deviation-walk scale and the
relative case ascertainment wherever the chain carries them.

Each panel carries its own y-axis, the quantities having different units and
differing by orders of magnitude. Panels whose quantity has a meaningful
reference value are drawn with it as a dashed rule, at one for the
reproduction number and the relative ascertainment and at zero for the log-Rt
deviations.

The deviations are sum-to-zero contrasts around the national trend (see
[`patch_rt_model`](@ref)), so `δ` is read relative to the national average
across provinces rather than to any one patch. Ascertainment and the
reproduction number must be read together. The case composition identifies
only their product, and it is the per-province deaths that tilt the balance
between them.
"""
function plot_patch_summary(
        chn, n_patches::Integer = length(PROVINCE_NAMES);
        patch_labels::AbstractVector = PROVINCE_LABELS,
        colours = [:firebrick, :steelblue, :seagreen],
        ncols::Integer = 4,
        title::AbstractString = "Per-province posterior summary"
    )
    required = [:C_T_patch, :R_T_patch, :infections_T_patch, :delta_patch]
    absent = filter(q -> !_has_key(chn, q), required)
    isempty(absent) || error(
        "chain is missing the per-patch deterministics $(absent); it was " *
            "not sampled from `bvd_joint`."
    )
    np = min(n_patches, length(patch_labels))
    ## Quantity, panel label, and the reference value worth a rule, in the
    ## order `patch_summary_table` reports them. The optional ones are absent
    ## from a chain fitted without the matching model piece.
    panels = Tuple{Symbol, String, Union{Nothing, Float64}}[
        (:C_T_patch, "Cumulative infections", nothing),
        (:R_T_patch, "Reproduction number", 1.0),
        (:infections_T_patch, "Daily infections at cut-off", nothing),
        (:delta_patch, "log-Rt deviation from trend", 0.0),
    ]
    optional = [
        (:log_rt_contrast, "log-Rt vs primary patch", 0.0),
        (:region_drift_sd, "Rt deviation innovation sd", nothing),
        (:province_ascertainment, "Relative case ascertainment", 1.0),
    ]
    for o in optional
        _has_key(chn, first(o)) && push!(panels, o)
    end
    nc = min(ncols, length(panels))
    nr = cld(length(panels), nc)
    fig = Figure(; size = (420 * nc, 340 * nr))
    xs = Float64.(1:np)
    for (k, (sym, label, reference)) in enumerate(panels)
        r, c = cld(k, nc), mod1(k, nc)
        ax = Axis(
            fig[r, c]; ylabel = label, title = label,
            xticks = (xs, String.(patch_labels[1:np])),
            xticklabelrotation = pi / 6
        )
        reference === nothing || hlines!(
            ax, [reference]; color = :black,
            linestyle = :dash, linewidth = 1
        )
        draws = _per_patch(chn, sym, np)
        for p in 1:np
            _draw_patch_interval!(
                ax, xs[p], draws[p],
                colours[mod1(p, length(colours))]
            )
        end
        ## A single province would otherwise sit on the axis edge.
        CairoMakie.xlims!(ax, 0.5, np + 0.5)
    end
    CairoMakie.Label(
        fig[nr + 1, 1:nc],
        "Bars are 30/60/90% credible intervals, thickest for the 30%, with " *
            "the median as a dot. Each panel has its own y-axis. Dashed rules " *
            "mark the reference value: one for the reproduction number and the " *
            "relative ascertainment, zero for the log-Rt deviations.";
        fontsize = 12, padding = (0, 0, 0, 6)
    )
    CairoMakie.Label(fig[0, 1:nc], title; fontsize = 16, font = :bold)
    return fig
end

## Quantile bands over per-draw trajectories that may be undefined at a
## point. `_traj_bands` throws on a NaN, and a vintage whose observed total
## is zero has no share to allocate, so it has no predictive share either.
## Those points stay NaN, which Makie draws as a gap in the band.
function _traj_bands_missing(trajs, n::Integer)
    function q(d, pr)
        v = Float64[t[d] for t in trajs]
        any(isnan, v) && return NaN
        return quantile(v, pr)
    end
    return (
        lo90 = [q(d, 0.05) for d in 1:n], hi90 = [q(d, 0.95) for d in 1:n],
        lo60 = [q(d, 0.2) for d in 1:n], hi60 = [q(d, 0.8) for d in 1:n],
        lo30 = [q(d, 0.35) for d in 1:n], hi30 = [q(d, 0.65) for d in 1:n],
    )
end

## Predictive band behind an expected-value ribbon: the same 30/60/90%
## nesting in grey, with the 90% edges outlined so the outer extent stays
## legible where the coloured ribbon sits inside it.
function _draw_pred_bands!(ax, x, b)
    band!(ax, x, b.lo90, b.hi90; color = (:grey30, 0.12))
    band!(ax, x, b.lo60, b.hi60; color = (:grey30, 0.2))
    band!(ax, x, b.lo30, b.hi30; color = (:grey30, 0.28))
    lines!(
        ax, x, b.lo90; color = (:grey20, 0.7), linestyle = :dash,
        linewidth = 1
    )
    lines!(
        ax, x, b.hi90; color = (:grey20, 0.7), linestyle = :dash,
        linewidth = 1
    )
    return ax
end

## Chain keys carrying a composition's overdispersion, most specific first.
## The case composition exposes it as `province_composition_rho` and the
## death composition as `province_death_composition_rho`. A chain carrying
## neither still has the submodel's own sampled `rho` under its prefix.
function _composition_rho_keys(share_key::Symbol)
    return share_key === :province_death_shares ?
        [
            :province_death_composition_rho,
            Symbol("death_composition_state.ρ"),
        ] :
        [:province_composition_rho, Symbol("composition_state.ρ")]
end

## Per-draw composition overdispersion, or `nothing` when the chain carries
## none of the candidate keys, which then plots without the predictive band.
function _composition_rho_draws(chn, keys, nd::Integer)
    for k in keys
        _has_key(chn, k) || continue
        v = Float64[x for x in vec(collect(chn[k]))]
        length(v) == nd && return v
    end
    return nothing
end

## Predictive shares for every patch at every vintage, one trajectory per
## posterior draw. Each draw's expected shares `m` and overdispersion
## `rho[d]` are pushed back through the stick-breaking allocation
## [`province_composition_model`](@ref) scores (see `_composition_counts`).
## The returned shares carry the composition's extra-Multinomial scatter as
## well as the posterior width of the expected share.
##
## The vintage's observed total is the trial count, matching the fitted
## likelihood, so the check is on the split alone. A vintage with no observed
## cases has no split to predict and stays NaN.
function _composition_predictive(
        ms, rho, totals, nv::Integer;
        seed::Integer = 20_240
    )
    counts = _composition_counts(ms, rho, fill(totals, length(ms)); seed)
    return [
        [
            [totals[i] > 0 ? c[i] / totals[i] : NaN for i in 1:nv]
                for c in cp
        ]
            for cp in counts
    ]
end

## Predictive counts for every patch at every vintage, one trajectory per
## posterior draw, allocating `totals[d][i]` for draw `d` at vintage `i`.
## Each draw's split is drawn from [`composition_split_model`](@ref), the
## stick-breaking the fitted composition scores and the forecast splits
## with, so there is one province split. A vintage with a total of zero
## allocates zero to every patch.
##
## The seed is fixed so a rebuilt report redraws the same band rather than
## moving it by the Monte Carlo error of the simulation.
function _composition_counts(ms, rho, totals; seed::Integer = 20_240)
    rng = MersenneTwister(seed)
    np, nv = size(first(ms))
    out = [[zeros(Int, nv) for _ in ms] for _ in 1:np]
    for (d, m) in enumerate(ms)
        split = composition_split_model(
            missing, m, [max(Int(t), 0) for t in totals[d]], rho[d]
        )
        x = first(
            init!!(rng, split, VarInfo(), InitFromPrior())
        ).obs_increments
        for p in 1:np, i in 1:nv
            out[p][d][i] = x[p, i]
        end
    end
    return out
end

"""
Posterior predictive check on a per-province composition: the modelled share
of each province at every spatial vintage, with the observed share drawn
over it as black points.

Each panel carries two bands. The grey predictive band is what the observed
points are drawn from. Every posterior draw's expected shares and composition
overdispersion are pushed back through the stick-breaking allocation
[`province_composition_model`](@ref) scores, at that vintage's observed
total, so the band carries the composition's extra-Multinomial scatter on top
of the posterior width. The coloured ribbon inside it is the expected share
alone, the modelled centre the points scatter around.

`share_key` is the chain's share deterministic (`province_shares` for the
confirmed cases, `province_death_shares` for the confirmed deaths),
`obs_increments` the matching `(n_patches x n_vintages)` observed increments
and `days` their grid days, both from [`province_increment_matrix`](@ref).

Shares rather than counts is what the model scores. The per-province totals
are an exact partition of the national totals the confirmed streams already
carry, so only the split is new information (see
[`province_composition_model`](@ref)).

Every panel starts at zero but takes its own upper limit. The shares differ
by orders of magnitude, so a common axis leaves every province but the
epicentre pinned to the floor.

`rho_key` names the chain's overdispersion for this composition and defaults
to the one matching `share_key`. A chain carrying neither that deterministic
nor the submodel's own draw is drawn with the expected-share ribbon only.
"""
function plot_province_composition_ppc(
        chn; share_key::Symbol,
        obs_increments::AbstractMatrix, days::AbstractVector{<:Integer},
        seeding::Date,
        n_patches::Integer = length(PROVINCE_NAMES),
        patch_labels::AbstractVector = PROVINCE_LABELS,
        colours = [:firebrick, :steelblue, :seagreen],
        rho_key::Union{Nothing, Symbol} = nothing,
        title::AbstractString = "Province share, modelled against observed"
    )
    np = min(n_patches, length(patch_labels))
    ms = [collect(v) for v in vec(collect(chn[share_key]))]
    nv = size(first(ms), 2)
    nv == length(days) || error(
        "plot_province_composition_ppc: `$(share_key)` holds $nv vintages " *
            "but `days` holds $(length(days))."
    )
    epoch = date2epochdays(seeding)
    x = Float64[epoch + (d - 1) for d in days]
    totals = [sum(@view obs_increments[:, i]) for i in 1:nv]
    ## The composition's overdispersion sets how far the observed shares are
    ## expected to scatter from the modelled ones, so the predictive band
    ## needs it.
    rho_keys = rho_key === nothing ? _composition_rho_keys(share_key) :
        [rho_key]
    rho = _composition_rho_draws(chn, rho_keys, length(ms))
    preds = rho === nothing ? nothing :
        _composition_predictive(ms, rho, totals, nv)
    fig = Figure(; size = (460 * np, 400))
    for p in 1:np
        colour = colours[mod1(p, length(colours))]
        trajs = [Float64[m[p, i] for i in 1:nv] for m in ms]
        b = _traj_bands(trajs, nv)
        ax = Axis(
            fig[1, p]; xlabel = "Vintage", ylabel = "Share of total",
            title = patch_labels[p], titlecolor = colour,
            xticklabelrotation = pi / 6
        )
        ## Predictive first, so the narrower expected-share ribbon stays
        ## readable on top of it.
        preds === nothing ||
            _draw_pred_bands!(ax, x, _traj_bands_missing(preds[p], nv))
        _draw_traj_bands!(ax, x, b, colour)
        obs = [
            totals[i] > 0 ? obs_increments[p, i] / totals[i] : NaN
                for i in 1:nv
        ]
        CairoMakie.scatter!(ax, x, obs; color = :black, markersize = 8)
        CairoMakie.ylims!(ax, 0, nothing)
        loax = floor(Int, minimum(x))
        hiax = ceil(Int, maximum(x))
        ax.xticks = collect(loax:7:hiax)
        ax.xtickformat = vals -> [
            string(epochdays2date(round(Int, v)))
                for v in vals
        ]
    end
    caption = preds === nothing ?
        "Bands are 30/60/90% credible intervals on the expected " *
        "share. Black points are the observed share at each vintage. " *
        "Each panel starts at zero and takes its own upper limit." :
        "Grey band is the 30/60/90% posterior predictive interval on " *
        "the observed share, dashed at its 90% edges. The coloured " *
        "ribbon inside it is the same intervals on the expected " *
        "share. Black points are the observed share at each vintage " *
        "and should fall inside the grey band. Each panel starts at " *
        "zero and takes its own upper limit."
    CairoMakie.Label(
        fig[2, 1:np], caption;
        fontsize = 12, padding = (0, 0, 0, 6)
    )
    CairoMakie.Label(fig[0, 1:np], title; fontsize = 16, font = :bold)
    return fig
end

"""
Two-panel density of the no-onward-transmission counterfactual from
[`predict_no_onward_deaths`](@ref). The left panel shows the *still expected*
deaths (`:delta_deaths`, the future deaths in cases already infected by `T`,
net of the `obs_deaths` already observed). The right panel shows the
*projected total* (`:total_projected = obs_deaths + delta_deaths`), whose
axis starts at `obs_deaths`. Both are lower bounds, assuming every onward
transmission stops at time `T`.
"""
function plot_no_onward_deaths(df::DataFrame; obs_deaths::Real)
    fig = Figure(; size = (980, 420))

    ## `delta_deaths` is clamped at zero in `predict_no_onward_deaths`, so the
    ## projected total cannot fall below the deaths already observed.
    ax1 = Axis(
        fig[1, 1];
        xlabel = "Still expected deaths (beyond those already observed)",
        ylabel = "Posterior density",
        title = "Still expected (future)"
    )
    _bounded_density!(
        ax1, df.delta_deaths; lower = 0,
        color = (:firebrick, 0.5), strokecolor = :firebrick, strokewidth = 2
    )

    ## The axis starts at the deaths already observed, so the left spine is
    ## that reference and a rule drawn on it would be invisible.
    ax2 = Axis(
        fig[1, 2];
        xlabel = "Projected total deaths, from the $(obs_deaths) observed",
        ylabel = "Posterior density",
        title = "Projected total"
    )
    _bounded_density!(
        ax2, df.total_projected; lower = obs_deaths,
        color = (:firebrick, 0.5), strokecolor = :firebrick, strokewidth = 2
    )

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
    ax = Axis(
        fig[r, c];
        xlabel = title, ylabel = "Predictive frequency",
        title = "One week ahead", limits = ((0, upper), nothing)
    )
    vspan!(ax, lo, hi; color = (colour, 0.15))
    hist!(ax, v; bins = range(0, upper; length = 30), color = (colour, 0.7))
    return ax
end

"""
One-week-ahead forecast of the unobserved (latent) quantities from
[`forecast_reported`](@ref): new infections, new symptom onsets and new
deaths over the horizon, with the reproduction number left to keep evolving.
Each count panel histograms the projected new latent count with its 90%
predictive interval shaded. The reproduction-number panel shows the posterior
of the end-of-horizon forecast `R_t` with the no-growth line at one marked.
These are the latent counterparts of the observed-stream forecast in
[`plot_forecast`](@ref).
"""
function plot_forecast_latent(fc::DataFrame)
    count_cols = [
        (:infections_new, "New infections (DRC)", :steelblue),
        (:onsets_new, "New symptom onsets (DRC)", :seagreen),
        (:deaths_latent_new, "New deaths (DRC)", :firebrick),
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
    ax = Axis(
        fig[r, c];
        xlabel = "Forecast reproduction number (DRC)",
        ylabel = "Posterior density", title = "One week ahead"
    )
    ## A reproduction number cannot be negative, and the horizon walk leaves
    ## draws close enough to zero for the kernel to spill past zero.
    _bounded_density!(
        ax, rt; lower = 0, color = (:purple, 0.5),
        strokecolor = :purple, strokewidth = 2
    )
    vlines!(ax, [1.0]; color = :black, linestyle = :dash, linewidth = 2)
    return fig
end

## Panel colours of the confirmed streams, keyed by forecast column, shared
## by the national and per-province forecast figures.
const _CONFIRMED_FORECAST_COLOURS = (
    confirmed_new = :goldenrod, confirmed_deaths_new = :darkorange3,
)

"""
One-week-ahead forecast of the observed count streams from
[`forecast_reported`](@ref): the new count each stream adds over the horizon.
Panels cover reported cases, suspected deaths, laboratory-confirmed cases,
confirmed deaths and recovered, each drawn only when the forecast carries
that stream's `*_new` column. Passing a column subset restricts the panels
further, which is how the report draws the stopped streams on their own. Each
panel histograms the projected new count with its 90% predictive interval
shaded. The latent counterparts are in [`plot_forecast_latent`](@ref).
"""
function plot_forecast(fc::DataFrame)
    count_cols = Tuple{Symbol, String, Symbol}[]
    for (col, title, colour) in (
            (:cases_new, "New reported cases (DRC)", :steelblue),
            (:deaths_new, "New suspected deaths (DRC)", :firebrick),
            (
                :confirmed_new, "New confirmed cases (DRC)",
                _CONFIRMED_FORECAST_COLOURS.confirmed_new,
            ),
            (
                :confirmed_deaths_new, "New confirmed deaths (DRC)",
                _CONFIRMED_FORECAST_COLOURS.confirmed_deaths_new,
            ),
            (:recovered_new, "New recovered among confirmed (DRC)", :seagreen),
        )
        col in propertynames(fc) || continue
        push!(count_cols, (col, title, colour))
    end
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
    :admissions_fc in propertynames(fc) && push!(
        count_cols,
        (:admissions_fc, "New isolation admissions (DRC)", :steelblue)
    )
    :incare_deaths_fc in propertynames(fc) && push!(
        count_cols,
        (:incare_deaths_fc, "New in-care deaths (DRC)", :firebrick)
    )
    :ruleouts_fc in propertynames(fc) && push!(
        count_cols,
        (:ruleouts_fc, "New rule-outs (DRC)", :seagreen)
    )
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
One-week-ahead forecast split by province, for the two streams the spatial
tables report: the new confirmed cases and confirmed deaths expected in each
province over the week to `T + 7`. One panel per stream, the provinces side
by side on a shared axis, each drawn as a median dot over nested 30/60/90%
credible bars, in the style of [`plot_patch_summary`](@ref).

This is the figure form of [`province_forecast_table`](@ref), and the figure
the release archive [`province_forecast_archive`](@ref) carries the draws
behind.

`fc` is a [`forecast_provinces`](@ref) frame. A national
[`forecast_reported`](@ref) result is replaced by the one-week province
forecast read from the posterior-predictive draws `pp`.

Panels are drawn only for the streams `fc` carries, so a forecast without the
confirmed deaths column shows the cases panel alone, and a forecast carrying
neither returns an empty figure.
"""
function plot_province_forecast(
        pp, fc::DataFrame;
        n_patches::Integer = length(PROVINCE_NAMES),
        patch_labels::AbstractVector = PROVINCE_LABELS,
        colours = [:firebrick, :steelblue, :seagreen],
        title::AbstractString = "One-week-ahead forecast by province"
    )
    np = min(n_patches, length(patch_labels))
    entries = _province_forecast_draws(pp, fc, np, patch_labels)
    isempty(entries) && return Figure()
    ## One panel per stream, each holding every province's interval on the
    ## shared province axis.
    labels = unique(first.(entries))
    nc = length(labels)
    fig = Figure(; size = (420 * nc, 380))
    xs = Float64.(1:np)
    for (k, label) in enumerate(labels)
        sel = [e for e in entries if e[1] == label]
        ax = Axis(
            fig[1, k]; ylabel = "Forecast count over the week",
            title = "New $(label) by T+7",
            xticks = (xs, String.(patch_labels[1:np])),
            xticklabelrotation = pi / 6
        )
        for (p, e) in enumerate(sel)
            _draw_patch_interval!(
                ax, xs[p], e[3],
                colours[mod1(p, length(colours))]
            )
        end
        ## A single province would otherwise sit on the axis edge.
        CairoMakie.xlims!(ax, 0.5, np + 0.5)
        ## A count cannot be negative and the panel is read against zero, so
        ## the axis starts there rather than at the smallest lower bound.
        CairoMakie.ylims!(ax, 0, nothing)
    end
    ## Two panels is a narrower figure than the per-province summary grid, so
    ## the caption wraps to the layout width.
    CairoMakie.Label(
        fig[2, 1:nc],
        "Bars are 30/60/90% credible intervals, thickest for the 30%, with " *
            "the median as a dot. Each province is forecast by the fitted " *
            "patch model, and the provinces add up to the national forecast.";
        fontsize = 12, word_wrap = true, padding = (0, 0, 0, 6)
    )
    CairoMakie.Label(fig[0, 1:nc], title; fontsize = 16, font = :bold)
    return fig
end

"""
One-week-ahead forecast for a single province, the per-province counterpart
of [`plot_forecast`](@ref): the new confirmed cases and confirmed deaths
expected in patch `province` over the week to `T + 7`, one histogram panel
per stream with its 90% predictive interval shaded.

The draws are the ones [`plot_province_forecast`](@ref) summarises, from a
[`forecast_provinces`](@ref) frame. A national [`forecast_reported`](@ref)
result is replaced by the one-week province forecast read from the draws
`pp`.

`observed` optionally gives a recent observed week per stream, keyed by the
forecast column (`confirmed_new`, `confirmed_deaths_new`), for example from
[`province_recent_counts`](@ref). Each is drawn as a dashed rule, and the
axis widens to hold it. Panels are drawn only for the streams `fc` carries.
"""
function plot_province_forecast_detail(
        pp, fc::DataFrame;
        province::Integer,
        n_patches::Integer = length(PROVINCE_NAMES),
        patch_labels::AbstractVector = PROVINCE_LABELS,
        observed::NamedTuple = (;)
    )
    np = min(n_patches, length(patch_labels))
    1 <= province <= np || throw(
        ArgumentError("province must be in 1:$np; got $province")
    )
    label = patch_labels[province]
    entries = [
        e for e in _province_forecast_draws(pp, fc, np, patch_labels)
            if e[2] == label
    ]
    isempty(entries) && return Figure()
    cols = Dict(
        label => col for (col, label) in _PROVINCE_FORECAST_STREAMS
    )
    ncols = length(entries)
    fig = Figure(; size = (400 * ncols, 360))
    for (i, (stream, _, draws)) in enumerate(entries)
        col = cols[stream]
        ax = _forecast_count_panel!(
            fig, (1, i), draws, "New $(stream) ($(label))",
            _CONFIRMED_FORECAST_COLOURS[col]
        )
        haskey(observed, col) || continue
        o = float(observed[col])
        vlines!(ax, [o]; color = :black, linestyle = :dash, linewidth = 2)
        CairoMakie.xlims!(
            ax, 0, max(1.0, quantile(draws, 0.98), 1.05 * o)
        )
    end
    return fig
end

"""
One-week-ahead isolation/treatment-bed forecast from
[`forecast_reported`](@ref): the projected bed demand (the need a week ahead,
under unconstrained supply) against the occupancy the situation reports would
print, and the demand above the beds available. The left panel overlays the
two predictive distributions. The right panel histograms the shortfall, which
is the need above the capacity rather than the gap between the two panels'
densities, since the occupancy carries the fitted reporting-basis offset and
the shortfall does not. Drawn only when the forecast carries the bed streams
(`bed_demand` and `isolation_level`).

The model carries a single national bed capacity, so it cannot represent
local saturation. On 13 June Ituri was at 93.9% occupancy while Sud-Kivu was
at 21.9%, and beds free in one province cannot serve patients in another, so
the national shortfall understates the local unmet need.
"""
function plot_forecast_beds(fc::DataFrame)
    (
        :bed_demand in propertynames(fc) &&
            :isolation_level in propertynames(fc)
    ) || return Figure()
    demand = float.(fc[!, :bed_demand])
    occ = float.(fc[!, :isolation_level])
    ## The shortfall is the need above the beds available. Without that
    ## column the occupancy is the demand capped at the capacity, so their
    ## difference is the same quantity.
    shortfall = :bed_shortfall in propertynames(fc) ?
        float.(fc[!, :bed_shortfall]) : max.(demand .- occ, 0.0)
    fig = Figure(; size = (800, 360))
    ## Cap the x-axis at the 98th percentile of demand. The unconstrained
    ## projection is heavy-tailed, so its long upper tail otherwise squashes
    ## the readable bulk of both densities. Occupancy is capped at capacity,
    ## so it sits below this bound.
    upper = max(1.0, quantile(demand, 0.98))
    ax1 = Axis(
        fig[1, 1];
        xlabel = "Isolation beds a week ahead (DRC)",
        ylabel = "Predictive density", title = "Need vs supply-limited use",
        limits = ((0, upper), nothing)
    )
    density!(
        ax1, demand; color = (:darkorange, 0.35),
        strokecolor = :darkorange, strokewidth = 2, label = "Demand (need)"
    )
    density!(
        ax1, occ; color = (:steelblue, 0.35),
        strokecolor = :steelblue, strokewidth = 2,
        label = "Occupancy (supply-limited)"
    )
    CairoMakie.axislegend(ax1; position = :rt, framevisible = false)
    _forecast_count_panel!(
        fig, (1, 2), shortfall, "Bed shortfall (DRC)",
        :firebrick
    )
    return fig
end

"""
Validate a [`forecast_reported`](@ref) bed projection against the beds
actually occupied a week later. Histograms the projected reported
isolation-bed occupancy with the 90% predictive interval shaded and the
`isolation` count observed at the target date drawn as a dashed black rule,
so last week's bed forecast is scored against what the beds held. Drawn only
when the forecast carries `isolation_level`.

`individual`, when given, is a second predictive sample: the frozen
individual (treatment-only) model's own forecast draws at the same cut-off,
from [`forecast_stream`](@ref), overlaid as a dotted step outline on
the joint's own bins so both bed
forecasts are visible against the observed occupancy.

At a one-week-back freeze the bed capacity has no implied-capacity anchor,
the reported occupancy rate starting only on 9 June, so the projected
occupancy rides the capacity random walk back to the freeze date and the
interval is wide.
"""
function plot_forecast_beds_vs_truth(
        fc::DataFrame;
        isolation::Union{Real, Missing},
        individual::Union{Nothing, AbstractVector} = nothing
    )
    (isolation !== missing && :isolation_level in propertynames(fc)) ||
        return Figure()
    v = float.(fc[!, :isolation_level])
    indiv = isnothing(individual) ? nothing : float.(individual)
    lo = quantile(v, 0.05)
    hi = quantile(v, 0.95)
    upper = max(
        1.0, quantile(v, 0.995), float(isolation) * 1.05,
        isnothing(indiv) || isempty(indiv) ? 0.0 :
            quantile(indiv, 0.995)
    )
    fig = Figure(; size = (440, 360))
    ax = Axis(
        fig[1, 1];
        xlabel = "Isolation beds occupied at the target date (DRC)",
        ylabel = "Predictive frequency", title = "Forecast vs observed",
        limits = ((0, upper), nothing)
    )
    vspan!(ax, lo, hi; color = (:steelblue, 0.15))
    bins = range(0, upper; length = 30)
    joint_h = hist!(ax, v; bins = bins, color = (:steelblue, 0.7))
    handles = Any[joint_h]
    labels = String["joint"]
    ## The same bins and the same draw count as the joint's histogram, so
    ## the two read on one scale. See `plot_forecast_vs_truth` for why a
    ## kernel density cannot be drawn against a count axis.
    if !isnothing(indiv) && !isempty(indiv) && length(unique(indiv)) > 1
        indiv_h = CairoMakie.stephist!(
            ax, indiv; bins = bins,
            weights = fill(length(v) / length(indiv), length(indiv)),
            color = :black, linewidth = 2, linestyle = :dot
        )
        push!(handles, indiv_h)
        push!(labels, "individual")
    end
    vlines!(
        ax, [float(isolation)]; color = :black, linestyle = :dash,
        linewidth = 2
    )
    length(handles) > 1 && CairoMakie.axislegend(
        ax, handles, labels;
        position = :rt, framevisible = false
    )
    return fig
end

"""
Validation figure for a [`forecast_reported`](@ref) projection, laid out as a
two-row grid with one column per scored stream: the top row shows the
cumulative forecast distribution, the bottom row the new count forecast over
the horizon. Each panel is a histogram with the 90% predictive interval shaded
and the later-observed count drawn as a dashed black rule, so the forecast
distribution is scored against the count that was actually observed.

Streams are drawn in the order reported cases, suspected deaths,
laboratory-confirmed cases, confirmed deaths and recovered, each shown only
when the forecast carries that stream's `*_cum`/`*_new` columns and an
observed cumulative count is supplied for it.

`observed` maps a stream's cumulative column (`:confirmed_cum`, `:cases_cum`,
…) to its observed cumulative count at the target date. A stream absent from
`observed` is skipped, which is how the caller withholds a stream whose
reporting does not cover the target date (see [`stream_reporting`](@ref)).
`baseline` maps the same columns to the cumulative count at the forecast
origin (default `0`), and `breaks` to each stream's retrospective
harmonisation correction over the forecast window (see
[`confirmed_break_correction`](@ref), default `0`), so the observed new count
is `max(observed − baseline − breaks, 0)`.

`individual` maps a stream's new-count column (`:confirmed_new`,
`:cases_new`, …) to that stream's own frozen single-stream model's forecast
draws of the new count, from [`forecast_stream`](@ref). A stream present in
`individual` gets a second, dotted outline over the joint's histogram on
both its panels, on the same bins and reweighted to the joint's draw count,
the cumulative panel from `baseline + individual` new-count draws. A stream
absent from it draws the joint alone. The latent counterparts are scored
distribution-versus-distribution by
[`plot_forecast_vs_truth_latent`](@ref).
"""
function plot_forecast_vs_truth(
        fc::DataFrame;
        observed::NamedTuple, baseline::NamedTuple = NamedTuple(),
        breaks::NamedTuple = NamedTuple(),
        individual::NamedTuple = NamedTuple()
    )
    specs = (
        (:cases_cum, :cases_new, "reported cases (DRC)", :steelblue),
        (:deaths_cum, :deaths_new, "suspected deaths (DRC)", :firebrick),
        (:confirmed_cum, :confirmed_new, "confirmed cases (DRC)", :goldenrod),
        (
            :confirmed_deaths_cum, :confirmed_deaths_new,
            "confirmed deaths (DRC)", :darkorange3,
        ),
        (
            :recovered_cum, :recovered_new,
            "recovered among confirmed (DRC)", :seagreen,
        ),
    )
    streams = Vector{
        Tuple{
            Symbol, Symbol, String, Symbol, Float64, Float64, Float64,
            Union{Nothing, Vector{Float64}},
        },
    }()
    for (cumcol, newcol, name, colour) in specs
        (cumcol in propertynames(fc) && haskey(observed, cumcol)) || continue
        ## Both truths are what was notified across the week, so a
        ## retrospective harmonisation sitting in the reported cumulative
        ## comes out of each. The projected cumulative carries none.
        brk = float(get(breaks, cumcol, 0))
        obs = float(observed[cumcol]) - brk
        base = float(get(baseline, cumcol, 0))
        indiv_new = haskey(individual, newcol) ?
            Float64.(individual[newcol]) : nothing
        push!(
            streams,
            (
                cumcol, newcol, name, colour, obs, obs - base, base,
                indiv_new,
            )
        )
    end
    ncols = length(streams)
    ncols == 0 && return Figure()
    fig = Figure(; size = (370 * ncols, 680))
    any_indiv = false
    function panel!(row, col, v, obs, title, colour, indiv)
        lo = quantile(v, 0.05)
        hi = quantile(v, 0.95)
        upper = max(
            1.0, quantile(v, 0.995), obs * 1.05,
            isnothing(indiv) || isempty(indiv) ? 0.0 :
                quantile(indiv, 0.995)
        )
        ax = Axis(
            fig[row, col];
            xlabel = title, ylabel = "Predictive frequency",
            limits = ((0, upper), nothing)
        )
        vspan!(ax, lo, hi; color = (colour, 0.15))
        bins = range(0, upper; length = 30)
        hist!(ax, v; bins = bins, color = (colour, 0.7))
        ## A dotted outline over the same bins as the histogram, so the two
        ## fits' forecasts read apart on one scale. A kernel density would be
        ## drawn on the density scale against an axis counting draws, which
        ## puts it flat on the floor of the panel whatever it says. The
        ## individual fit's draws are reweighted to the joint's count, so the
        ## two outlines are comparable even when the chains differ in length.
        if !isnothing(indiv) && !isempty(indiv) && length(unique(indiv)) > 1
            CairoMakie.stephist!(
                ax, indiv; bins = bins,
                weights = fill(length(v) / length(indiv), length(indiv)),
                color = :black, linewidth = 2, linestyle = :dot
            )
            any_indiv = true
        end
        return vlines!(ax, [obs]; color = :black, linestyle = :dash, linewidth = 2)
    end
    for (j, entry) in enumerate(streams)
        ccol, ncol, name, colour, obs_cum, obs_new, origin,
            indiv_new = entry
        ## The individual fit forecasts new counts from the frozen origin, so
        ## its cumulative overlay is anchored there.
        indiv_cum = isnothing(indiv_new) ? nothing : indiv_new .+ origin
        panel!(
            1, j, fc[!, ccol], obs_cum, "Cumulative $name", colour,
            indiv_cum
        )
        panel!(
            2, j, fc[!, ncol], max(obs_new, 0.0), "New $name", colour,
            indiv_new
        )
    end
    if any_indiv
        joint_marker = CairoMakie.PolyElement(; color = (:grey, 0.7))
        indiv_marker = CairoMakie.LineElement(;
            color = :black,
            linestyle = :dot, linewidth = 2
        )
        CairoMakie.Legend(
            fig[0, 1:ncols], [joint_marker, indiv_marker],
            ["joint", "individual"]; orientation = :horizontal,
            framevisible = false, tellwidth = false
        )
    end
    return fig
end

"""
Latent-quantity validation figure. For each unobserved quantity (new
infections, new symptom onsets, new deaths over the past week) the
distribution the frozen last-week fit forecast is overlaid against the
distribution the current fit now estimates for the same window. Both are
latent, so the comparison is density versus density rather than density
versus a single observed count.

`fc` is the frozen forecast from [`forecast_reported`](@ref), its `*_new`
latent columns. `now` is a `NamedTuple` carrying the current fit's draws of
the same quantities, `(; infections_new, onsets_new, deaths_latent_new)`.
"""
function plot_forecast_vs_truth_latent(fc::DataFrame; now::NamedTuple)
    panels = [
        (:infections_new, "New infections (DRC)", :steelblue),
        (:onsets_new, "New symptom onsets (DRC)", :seagreen),
        (:deaths_latent_new, "New deaths (DRC)", :firebrick),
    ]
    ncols = length(panels)
    fig = Figure(; size = (370 * ncols, 380))
    local frozen_h, now_h
    for (j, (col, title, colour)) in enumerate(panels)
        vf = float.(fc[!, col])
        vn = float.(getproperty(now, col))
        upper = max(1.0, quantile(vf, 0.99), quantile(vn, 0.99))
        ax = Axis(
            fig[1, j];
            xlabel = title, ylabel = "Posterior density",
            limits = ((0, upper), nothing)
        )
        frozen_h = density!(
            ax, vf; color = (colour, 0.25),
            strokecolor = colour, strokewidth = 2
        )
        now_h = density!(
            ax, vn; color = (:grey, 0.0),
            strokecolor = :black, strokewidth = 2, linestyle = :dash
        )
    end
    CairoMakie.Legend(
        fig[1, ncols + 1], [frozen_h, now_h],
        ["Forecast last week", "Estimated now"]
    )
    return fig
end

## Date ticks for a vintage-indexed axis, at about one label a week. `dates`
## are the panel's vintage labels, one per x position, so the returned
## positions index into them. Labels are picked by calendar date rather than
## by index stride, the vintages being irregular. Selection walks back from
## the latest vintage, so the last vintage is always labelled. Once a weekly
## cadence would need more than `max_ticks` labels the step widens in whole
## weeks. A series too short to carry three ticks keeps every vintage.
function _vintage_ticks(
        dates::AbstractVector; step_days::Integer = 7,
        max_ticks::Integer = 18
    )
    n = length(dates)
    n == 0 && return (Int[], String[])
    ds = [d isa Date ? d : Date(String(d)) for d in dates]
    weeks = (ds[end] - ds[begin]).value ÷ step_days + 1
    step = step_days * max(1, cld(weeks, max_ticks))
    keep = [n]
    for i in (n - 1):-1:1
        (ds[last(keep)] - ds[i]).value >= step && push!(keep, i)
    end
    reverse!(keep)
    length(keep) < 3 && (keep = collect(1:n))
    return (keep, [string(ds[i]) for i in keep])
end

"""
Per-vintage conditional one-step-ahead posterior-predictive for the DRC
streams. For each `panel` the predicted cumulative count at vintage `v`
conditions on the *observed* cumulative at the previous vintage and adds only
the posterior-predictive between-vintage increment,
``\\hat{y}_v = y_{v-1} + \\Delta_v`` with ``y_0 = 0``. The increment
``\\Delta_v`` is the per-bin replicate draw the model already samples in
predictive mode, so each step carries the full posterior uncertainty of the
*new* increment while grounding on what was observed at the preceding sitrep.
This is the filtered one-step-ahead predictive, so errors do not compound
across the series as they would under a running sum of the modelled
increments. The result is summarised by vintage as shaded 30/60/90%
predictive ribbons with the observed cumulative counts overlaid as points.

Each `panel` is a `NamedTuple` `(; title, dates, replicates, observed)`,
where `replicates` is a vector of per-draw increment vectors, one entry per
vintage oldest first, and `observed` the matching observed cumulative counts
used as the conditioning baselines. `colour` is optional per panel. A panel
may set `cumulative = false` to plot standalone per-day counts instead of a
cumulative series, such as the daily new-suspect inflow or the 24h analysed
volume. There is no previous-vintage baseline then, so each replicate is
plotted as its own daily count and the y-axis reads "Daily count". A panel
may also set `ylabel` to name its own y-axis, which is how an occupancy
census is kept from reading as either a running total or a count of new
events.

`max_date` (an ISO date string or `Date`) truncates every panel to the
vintages on or before that date, so streams that keep reporting past the
others are cut back to the shared last date. Without it the confirmed panel
runs further along the date axis than the suspected panel and reads as though
it overtakes it, when the two are simply shown to different end dates.

The date axis carries about one label a week and always labels the last
vintage.
"""
function plot_vintage_conditional_ppc(
        panels::AbstractVector; xlabel = "Sitrep date",
        max_date::Union{Nothing, Date, AbstractString} = nothing
    )
    cap = isnothing(max_date) ? nothing :
        (max_date isa Date ? max_date : Date(String(max_date)))
    ## An empty panel set has no grid to lay out, so return a blank figure
    ## rather than dividing by zero.
    isempty(panels) && return Figure()
    npanels = length(panels)
    ## Cap the grid at four columns so a large stream set lays out over
    ## several rows rather than one wide strip the page downscales into tiny
    ## panels.
    ncols = min(npanels, 4)
    nrows = cld(npanels, ncols)
    fig = Figure(; size = (460 * ncols, 420 * nrows))
    for (j, p) in enumerate(panels)
        row, col = cld(j, ncols), mod1(j, ncols)
        ## Drop vintages past the shared cap so every panel ends on the same
        ## date. The replicates and observed counts are truncated to match,
        ## keeping the conditional baselines aligned.
        keep = isnothing(cap) ? eachindex(p.dates) :
            [i for i in eachindex(p.dates) if Date(p.dates[i]) <= cap]
        dates = p.dates[keep]
        observed = p.observed[keep]
        replicates = [collect(r)[keep] for r in vec(collect(p.replicates))]
        n = length(dates)
        colour = get(p, :colour, :steelblue)
        ## A `cumulative = false` panel carries standalone per-day counts, so
        ## there is no previous-vintage baseline to condition on and each
        ## replicate is its own per-day count.
        cumulative = get(p, :cumulative, true)
        ## Observed cumulative at the previous vintage is the conditioning
        ## baseline for each step (`y_0 = 0`). `obs_prev[v]` is `y_{v-1}`.
        obs_cum = float.(observed)
        obs_prev = cumulative ?
            [v == 1 ? 0.0 : obs_cum[v - 1] for v in 1:n] : zeros(n)
        ## `ylabel` lets a panel that is neither a running total nor a
        ## per-day flow name its own axis.
        ylabel = get(
            p, :ylabel,
            cumulative ? (col == 1 ? "Cumulative count" : "") : "Daily count"
        )
        ## Each draw's conditional cumulative at vintage `v` is the observed
        ## previous cumulative plus the drawn increment `Δ_v`, with a zero
        ## baseline for a non-cumulative panel.
        cond = [obs_prev .+ r for r in replicates]
        q(i, pr) = quantile([c[i] for c in cond], pr)
        lo90 = [q(i, 0.05) for i in 1:n]
        hi90 = [q(i, 0.95) for i in 1:n]
        lo60 = [q(i, 0.2) for i in 1:n]
        hi60 = [q(i, 0.8) for i in 1:n]
        lo30 = [q(i, 0.35) for i in 1:n]
        hi30 = [q(i, 0.65) for i in 1:n]
        x = collect(1:n)
        ## Truncate the y-axis to a ceiling driven by the observed counts and
        ## the 60% band, so a heavy upper tail does not flatten the visible
        ## detail. The band clips at the axis limit.
        yupper = 1.6 * max(
            isempty(obs_cum) ? 1.0 : maximum(obs_cum),
            isempty(hi60) ? 1.0 : maximum(hi60), 1.0
        )
        ax = Axis(
            fig[row, col]; title = p.title, xlabel = xlabel,
            ylabel = ylabel,
            xticks = _vintage_ticks(dates),
            xticklabelrotation = pi / 4, xticklabelsize = 11,
            limits = (nothing, (0, yupper))
        )
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
consecutive vintages rather than the running cumulative, so a rise or a
slowdown reads off the height of each step instead of the slope of a
near-straight cumulative line.

For a cumulative panel the observed incidence is the between-vintage
increment, the first vintage being its own baseline. For a non-cumulative
panel it is the count itself. The replicates are already per-vintage
increments, so they are the modelled incidence directly and are summarised as
30/60/90% credible ribbons with the observed incidence overlaid. `panels`,
`max_date` and the weekly date axis match
[`plot_vintage_conditional_ppc`](@ref), including the optional per-panel
`ylabel`.
"""
function plot_vintage_incidence_ppc(
        panels::AbstractVector; xlabel = "Sitrep date",
        max_date::Union{Nothing, Date, AbstractString} = nothing
    )
    cap = isnothing(max_date) ? nothing :
        (max_date isa Date ? max_date : Date(String(max_date)))
    ## An empty panel set has no grid to lay out, so return a blank figure
    ## rather than dividing by zero.
    isempty(panels) && return Figure()
    npanels = length(panels)
    ## Cap the grid at four columns so a large stream set lays out over
    ## several rows rather than one wide strip the page downscales into tiny
    ## panels.
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
        ## series, the first vintage being its own baseline at zero, or the
        ## standalone count for a non-cumulative panel.
        obs_cum = float.(observed)
        obs_inc = cumulative ?
            [
                v == 1 ? obs_cum[v] : obs_cum[v] - obs_cum[v - 1]
                for v in 1:n
            ] : obs_cum
        ## As in the conditional view, `ylabel` overrides the default.
        ylabel = get(p, :ylabel, col == 1 ? "New per vintage" : "")
        ## The replicates are already per-vintage increments (per-day counts
        ## for a non-cumulative panel), so they are the modelled incidence.
        q(i, pr) = quantile([r[i] for r in replicates], pr)
        lo90 = [q(i, 0.05) for i in 1:n]
        hi90 = [q(i, 0.95) for i in 1:n]
        lo60 = [q(i, 0.2) for i in 1:n]
        hi60 = [q(i, 0.8) for i in 1:n]
        lo30 = [q(i, 0.35) for i in 1:n]
        hi30 = [q(i, 0.65) for i in 1:n]
        x = collect(1:n)
        yupper = 1.6 * max(
            isempty(obs_inc) ? 1.0 : maximum(obs_inc),
            isempty(hi60) ? 1.0 : maximum(hi60), 1.0
        )
        ax = Axis(
            fig[row, col]; title = p.title, xlabel = xlabel,
            ylabel = ylabel,
            xticks = _vintage_ticks(dates),
            xticklabelrotation = pi / 4, xticklabelsize = 11,
            limits = (nothing, (0, yupper))
        )
        band!(ax, x, lo90, hi90; color = (colour, 0.15))
        band!(ax, x, lo60, hi60; color = (colour, 0.28))
        band!(ax, x, lo30, hi30; color = (colour, 0.42))
        scatter!(ax, x, float.(obs_inc); color = :black, markersize = 9)
    end
    return fig
end

"""
Per-stream calibration of the one-step-ahead conditional posterior
predictive, plotting the table from [`stream_calibration`](@ref). Pass that
table as returned, with its prettified columns `Stream`, `50% coverage`,
`90% coverage` and `Bias`.

Two panels share a categorical y-axis of streams. The left panel marks each
stream's empirical 50% and 90% coverage against vertical dashed reference
lines at the nominal 0.5 and 0.9, so a well-calibrated stream sits on its
line and a marker to the left of it flags under-coverage. The right panel
marks the mean forecast `Bias`, negative for under-prediction and positive
for over-prediction, with a dashed line at zero.
"""
function plot_stream_calibration(tbl::DataFrame)
    ## Table order is kept but reversed for the y-axis, so the first stream
    ## reads at the top.
    streams = string.(tbl[!, "Stream"])
    cov50 = float.(tbl[!, "50% coverage"])
    cov90 = float.(tbl[!, "90% coverage"])
    bias = float.(tbl[!, "Bias"])
    n = length(streams)
    ## Categorical y positions, top-to-bottom in table order.
    y = collect(n:-1:1)
    height = max(360, 60 + 26 * n)
    fig = Figure(; size = (980, height))

    ax1 = Axis(
        fig[1, 1];
        xlabel = "Empirical coverage", title = "Interval coverage",
        yticks = (y, streams), limits = ((0, 1), nothing)
    )
    ## Nominal reference lines. A marker on its line is well calibrated.
    vlines!(
        ax1, [0.5]; color = (:steelblue, 0.6), linestyle = :dash,
        linewidth = 2
    )
    vlines!(
        ax1, [0.9]; color = (:seagreen, 0.6), linestyle = :dash,
        linewidth = 2
    )
    h50 = scatter!(ax1, cov50, y; color = :steelblue, markersize = 11)
    h90 = scatter!(
        ax1, cov90, y; color = :seagreen, markersize = 11,
        marker = :diamond
    )
    CairoMakie.axislegend(
        ax1, [h50, h90], ["50% interval", "90% interval"];
        position = :lt, framevisible = false
    )

    ## Zero is unbiased. The sign flags over- or under-prediction.
    bmax = max(1.0, maximum(abs.(bias)) * 1.1)
    ax2 = Axis(
        fig[1, 2];
        xlabel = "Mean forecast bias", title = "Forecast bias",
        yticks = (y, fill("", n)), limits = ((-bmax, bmax), nothing)
    )
    vlines!(ax2, [0.0]; color = :black, linestyle = :dash, linewidth = 2)
    scatter!(ax2, bias, y; color = :firebrick, markersize = 11)
    return fig
end

## --- Per-parameter fit diagnostics --------------------------------------

"""
Cumulative share of parameters at or below each R-hat, one line per fit.
Pass each fit as `"label" => chain`, or as `"label" => frame` where the frame
is one [`parameter_diagnostics`](@ref) has already produced.

A line that climbs to one just past the left edge is a fit where a handful
of parameters are bad and the rest are fine. A line that stays low across
the axis is a fit where most of the model has not converged. The axis runs
to the worst fit's `clip` quantile so one extreme parameter cannot stretch
it, and anything beyond that is drawn at the right edge.
"""
function plot_rhat_spread(
        fits::Pair{String}...; xmax = nothing,
        clip::Real = 0.995, thresholds = (1.01, 1.1),
        title::AbstractString = "Spread of R-hat across parameters"
    )
    series = [
        (f.first, sort(filter(isfinite, _as_diagnostics(f.second).rhat)))
            for f in fits
    ]
    series = [s for s in series if !isempty(s[2])]
    if isempty(series)
        fig = Figure(; size = (860, 160))
        CairoMakie.Label(
            fig[1, 1], "No fit carries R-hat diagnostics.";
            tellwidth = false, tellheight = false, color = (:black, 0.55)
        )
        return fig
    end
    hi = isnothing(xmax) ?
        max(
            maximum(quantile(v, clip) for (_, v) in series),
            maximum(thresholds) + 0.01
        ) : float(xmax)
    colours = CairoMakie.Makie.wong_colors()
    fig = Figure(; size = (860, 420))
    ax = Axis(
        fig[1, 1]; title = title, xlabel = "R-hat",
        ylabel = "Share of parameters at or below"
    )
    vlines!(
        ax, collect(thresholds); color = (:grey, 0.6),
        linestyle = :dash, linewidth = 1.5
    )
    handles = Any[]
    labels = String[]
    for (i, (label, v)) in enumerate(series)
        x = clamp.(v, 1.0, hi)
        y = collect(1:length(x)) ./ length(x)
        h = lines!(
            ax, x, y; color = colours[mod1(i, length(colours))],
            linewidth = 2
        )
        push!(handles, h)
        push!(labels, String(label))
    end
    CairoMakie.xlims!(ax, 1.0, hi)
    CairoMakie.ylims!(ax, 0.0, 1.02)
    CairoMakie.Legend(
        fig[2, 1], handles, labels;
        orientation = :horizontal, framevisible = true,
        tellheight = true, tellwidth = false, nbanks = 2
    )
    return fig
end

# Vector-valued parameters with the lowest bulk effective sample size, most
# degraded first, keeping only those with at least `min_elements` entries so
# each panel is a series rather than a handful of points. A deterministic
# copy of a sampled vector carries identical diagnostics, so a parameter
# whose diagnostics repeat one already picked is dropped rather than drawn
# twice.
function _worst_vector_parameters(
        df::DataFrame, n::Integer;
        min_elements::Integer = 8
    )
    groups = [
        g for g in unique(df.parameter)
            if count(==(g), df.parameter) >= min_elements
    ]
    isempty(groups) && return String[]
    mins = [_min_finite(df.ess_bulk[df.parameter .== g]) for g in groups]
    ord = sortperm(replace(mins, NaN => Inf))
    picked = String[]
    seen = Set{Vector{Float64}}()
    for g in groups[ord]
        key = sort(df.ess_bulk[df.parameter .== g])
        key in seen && continue
        push!(seen, key)
        push!(picked, String(g))
        length(picked) == n && break
    end
    return picked
end

"""
Bulk effective sample size against element index for the vector-valued
parameters of a fit that mix worst, one panel each. `fit` is a chain or a
frame [`parameter_diagnostics`](@ref) has already produced. Points are
coloured by whether the element's R-hat exceeds `rhat_threshold`.

The element index runs in model order, so for a random walk or a daily
latent series a higher index is later in the outbreak. Bad mixing piled up
at one end of a panel is a problem confined to that stretch of the window.
Bad mixing spread evenly across a panel is a problem with the whole walk.
"""
function plot_parameter_index_diagnostics(
        fit;
        groups::Union{Nothing, AbstractVector} = nothing,
        n_groups::Integer = 3, min_elements::Integer = 8,
        rhat_threshold::Real = 1.1, ess_threshold::Real = 100,
        labels = Dict{Symbol, String}(),
        title::AbstractString =
            "Mixing along the worst vector-valued parameters"
    )
    df = _as_diagnostics(fit)
    picked = isnothing(groups) ?
        _worst_vector_parameters(
            df, n_groups;
            min_elements = min_elements
        ) : [String(g) for g in groups]
    if isempty(picked)
        fig = Figure(; size = (860, 160))
        CairoMakie.Label(
            fig[1, 1],
            "No vector-valued parameter carries diagnostics.";
            tellwidth = false, tellheight = false, color = (:black, 0.55)
        )
        return fig
    end
    fig = Figure(; size = (860, 230 * length(picked) + 110))
    ok_handle = nothing
    bad_handle = nothing
    for (i, g) in enumerate(picked)
        sub = df[df.parameter .== g, :]
        sub = sub[sortperm(sub.index), :]
        ax = Axis(
            fig[i, 1]; title = String(get(labels, Symbol(g), g)),
            xlabel = i == length(picked) ? "Element index" : "",
            ylabel = "Bulk effective sample size", yscale = log10
        )
        hlines!(
            ax, [float(ess_threshold)]; color = (:grey, 0.6),
            linestyle = :dash, linewidth = 1.5
        )
        y = [isnan(v) ? 1.0 : max(v, 1.0) for v in sub.ess_bulk]
        bad = sub.rhat .> rhat_threshold
        if any(.!bad)
            h = scatter!(
                ax, sub.index[.!bad], y[.!bad];
                color = :steelblue, markersize = 6
            )
            ok_handle = something(ok_handle, h)
        end
        if any(bad)
            h = scatter!(
                ax, sub.index[bad], y[bad];
                color = :firebrick, markersize = 6
            )
            bad_handle = something(bad_handle, h)
        end
    end
    CairoMakie.Label(fig[0, 1], title; font = :bold, tellwidth = false)
    handles = Any[]
    legend_labels = String[]
    if !isnothing(ok_handle)
        push!(handles, ok_handle)
        push!(legend_labels, "R-hat at most $(rhat_threshold)")
    end
    if !isnothing(bad_handle)
        push!(handles, bad_handle)
        push!(legend_labels, "R-hat above $(rhat_threshold)")
    end
    isempty(handles) ||
        CairoMakie.Legend(
        fig[length(picked) + 1, 1], handles, legend_labels;
        orientation = :horizontal, framevisible = true,
        tellheight = true, tellwidth = false
    )
    return fig
end

"""
Where the divergent transitions sit against the posterior, one panel per
parameter in `params`. Each panel draws the full posterior as a density and
the divergent draws as ticks along the axis, with the middle 90% of the
divergent draws shaded.

Ticks spread under the whole density are divergences scattered through the
posterior, which points at the sampler settings. Ticks piled into one shaded
stretch are divergences confined to one region, which points at the geometry
there.
"""
function plot_divergence_locations(
        chn, params::AbstractVector;
        labels = Dict{Symbol, String}(), ncols::Integer = 3,
        title::AbstractString =
            "Divergent draws against the full posterior"
    )
    flag = _divergent_flags(chn)
    if !any(flag)
        fig = Figure(; size = (860, 160))
        CairoMakie.Label(
            fig[1, 1], "No divergent transitions to place.";
            tellwidth = false, tellheight = false, color = (:black, 0.55)
        )
        return fig
    end
    usedcols = min(ncols, length(params))
    nrows = cld(length(params), usedcols)
    fig = Figure(; size = (300 * usedcols, 230 * nrows + 100))
    for (i, p) in enumerate(params)
        r, c = fldmod1(i, usedcols)
        x = Float64.(vec(collect(chn[p])))
        ax = Axis(
            fig[r, c];
            title = String(get(labels, Symbol(p), string(p))),
            ylabel = c == 1 ? "Density" : ""
        )
        xd = x[flag]
        vspan!(
            ax, quantile(xd, 0.05), quantile(xd, 0.95);
            color = (:firebrick, 0.12)
        )
        density!(
            ax, x; color = (:steelblue, 0.35),
            strokecolor = :steelblue, strokewidth = 1.5
        )
        scatter!(
            ax, xd, fill(0.0, length(xd)); color = (:firebrick, 0.6),
            marker = :vline, markersize = 10
        )
    end
    CairoMakie.Label(
        fig[0, 1:usedcols], title; font = :bold,
        tellwidth = false
    )
    return fig
end

"""
Bulk effective sample size in a reference fit against the same parameter's
bulk effective sample size in the fits it is compared with, from the frame
[`diagnostic_contrast`](@ref) returns. Both axes are logarithmic and the
dashed line is equality.

A point on the line is a parameter the reference fit handles as well as the
comparison does. A point far below it is a parameter that mixes on its own
and stops mixing in the reference, so the cause is what the reference adds
rather than the parameter.
"""
function plot_diagnostic_contrast(
        df::DataFrame;
        xlabel::AbstractString = "Bulk effective sample size, comparison fit",
        ylabel::AbstractString = "Bulk effective sample size, reference fit",
        title::AbstractString = "Mixing in the reference against each fit"
    )
    if isempty(df)
        fig = Figure(; size = (860, 160))
        CairoMakie.Label(
            fig[1, 1], "No parameter is shared between fits.";
            tellwidth = false, tellheight = false, color = (:black, 0.55)
        )
        return fig
    end
    colours = CairoMakie.Makie.wong_colors()
    fig = Figure(; size = (860, 460))
    ax = Axis(
        fig[1, 1]; title = title, xlabel = xlabel, ylabel = ylabel,
        xscale = log10, yscale = log10
    )
    ## Both axes cover the same span so the equality line runs corner to
    ## corner and the distance below it reads the same on either axis.
    lo = max(
        1.0, 0.8 * min(
            minimum(df.ess_bulk),
            minimum(df.ess_bulk_reference)
        )
    )
    hi = 1.25 * max(
        maximum(df.ess_bulk), maximum(df.ess_bulk_reference),
        lo + 1
    )
    lines!(
        ax, [lo, hi], [lo, hi]; color = (:grey, 0.7), linestyle = :dash,
        linewidth = 1.5
    )
    handles = Any[]
    labels = String[]
    for (i, f) in enumerate(unique(df.fit))
        cell = df[df.fit .== f, :]
        h = scatter!(
            ax, clamp.(cell.ess_bulk, lo, hi),
            clamp.(cell.ess_bulk_reference, lo, hi);
            color = (colours[mod1(i, length(colours))], 0.6),
            markersize = 7
        )
        push!(handles, h)
        push!(labels, String(f))
    end
    CairoMakie.xlims!(ax, lo, hi)
    CairoMakie.ylims!(ax, lo, hi)
    CairoMakie.Legend(
        fig[2, 1], handles, labels;
        orientation = :horizontal, framevisible = true,
        tellheight = true, tellwidth = false, nbanks = 2
    )
    return fig
end

## --- Parameter recovery -------------------------------------------------

"""
Parameter recovery over every seed (see `scripts/recovery.jl`). `params` is
the stacked [`recovery_table`](@ref) of each seed with a `seed` column,
`draws` maps each seed to a frame of its thinned posterior draws and `prior`
is a frame of prior draws, one column per quantity in both.

The top panel puts every quantity and seed on one scale relative to the
truth: each seed's posterior median with its 50% and 90% intervals divided
by that seed's true value, on a log axis about a line at one. A quantity in
`difference` is shown as `exp(value - truth)` instead, the ratio of daily
growth factors for a growth rate. Below, one panel per quantity on its
natural scale holds the prior draws in grey, each seed's posterior in its
colour and each seed's true value as a dashed line; quantities in `log_x`
are drawn on a log axis. Panels are laid out `ncols` to a row in the order
of `quantities`, titled by `panel_labels` (by default `labels`), with
`row_labels` down the left when given. The axes cover the posteriors and
the truths, so a wide prior shows as a low grey floor. Returns the
`Figure`.
"""
function plot_recovery(
        params::DataFrame, draws::AbstractDict, prior::DataFrame;
        quantities::AbstractVector{<:AbstractString} = unique(params.quantity),
        labels::AbstractDict = Dict{String, String}(),
        panel_labels::AbstractDict = labels,
        row_labels::AbstractVector{<:AbstractString} = String[],
        log_x = String[], difference = String[], ncols::Integer = 4
    )
    label(q) = get(labels, q, q)
    seeds = sort(unique(params.seed))
    colours = CairoMakie.Makie.wong_colors()
    colour(s) = colours[mod1(findfirst(==(s), seeds), length(colours))]
    nq = length(quantities)
    nrows = cld(nq, ncols)
    fig = Figure(
        ;
        size = (
            240 * ncols + 120,
            100 + (8 + 12 * length(seeds)) * nq + 210 * nrows,
        )
    )

    ## Every quantity against its truth, top to bottom in `quantities` order.
    ax = Axis(
        fig[1, 1:ncols];
        xlabel = "Posterior relative to the truth (log scale)",
        yticks = (collect(nq:-1:1), label.(quantities)), xscale = log10,
        xticks = [0.1, 0.2, 0.5, 1, 2, 5, 10],
        xtickformat = vs -> [_recovery_tick(v) for v in vs],
        title = "Recovery across seeds"
    )
    vlines!(ax, [1.0]; color = :black, linestyle = :dash, linewidth = 1.5)
    dodge = length(seeds) > 1 ?
        range(-0.3, 0.3; length = length(seeds)) : [0.0]
    for (k, s) in enumerate(seeds)
        tab = params[params.seed .== s, :]
        for (i, q) in enumerate(quantities)
            j = findfirst(==(q), tab.quantity)
            j === nothing && continue
            r = tab[j, :]
            rel(x) = q in difference ? exp(x - r.truth) : x / r.truth
            y = nq - i + 1 + dodge[k]
            linesegments!(
                ax, [rel(r.lower_90), rel(r.upper_90)], [y, y];
                color = colour(s), linewidth = 1.5
            )
            linesegments!(
                ax, [rel(r.lower_50), rel(r.upper_50)], [y, y];
                color = colour(s), linewidth = 4
            )
            scatter!(ax, [rel(r.median)], [y]; color = colour(s))
        end
    end
    CairoMakie.Legend(
        fig[1, ncols + 1],
        [
            CairoMakie.LineElement(; color = colour(s), linewidth = 3)
                for s in seeds
        ],
        ["seed $s" for s in seeds];
        framevisible = false
    )

    ## Each quantity on its own scale, the prior behind the posteriors.
    for (i, q) in enumerate(quantities)
        r, c = fldmod1(i, ncols)
        logged = q in log_x
        tf(x) = logged ? log10.(filter(>(0), x)) : x
        axq = Axis(
            fig[r + 1, c]; title = get(panel_labels, q, q),
            yticklabelsvisible = false,
            yticksvisible = false,
            xtickformat = logged ?
                (vs -> [_recovery_tick(10^v) for v in vs]) :
                (vs -> [_recovery_tick(v) for v in vs])
        )
        tab = params[params.quantity .== q, :]
        post = [
            tf(Float64.(draws[s][!, q])) for s in seeds
                if haskey(draws, s) && q in names(draws[s])
        ]
        truths = tf(Float64.(tab.truth))
        span = vcat(truths, [quantile(p, [0.005, 0.995]) for p in post]...)
        if q in names(prior)
            pv = tf(filter(isfinite, Float64.(prior[!, q])))
            length(pv) > 1 &&
                density!(axq, pv; color = (:grey, 0.35), strokewidth = 0)
        end
        for s in seeds
            haskey(draws, s) && q in names(draws[s]) || continue
            v = tf(Float64.(draws[s][!, q]))
            density!(
                axq, v; color = (colour(s), 0.15), strokecolor = colour(s),
                strokewidth = 1.5
            )
        end
        for (t, s) in zip(truths, tab.seed)
            vlines!(
                axq, [t]; color = colour(s), linestyle = :dash, linewidth = 1.5
            )
        end
        lo, hi = extrema(span)
        pad = hi > lo ? 0.05 * (hi - lo) : max(abs(lo), 1.0) * 0.1
        CairoMakie.xlims!(axq, lo - pad, hi + pad)
        if logged
            ## Round values on the log axis, spaced to the decades shown, or
            ## finer when the posteriors span less than a factor of five.
            decades = hi - lo
            ms = decades > 2 ? (1,) : decades > 1 ? (1, 3) : (1, 2, 5)
            steps = [m * 10.0^k for k in -4:9 for m in ms]
            ticks = filter(v -> lo - pad <= log10(v) <= hi + pad, steps)
            if length(ticks) < 3
                steps = [m * 10.0^k for k in -4:9 for m in 1:9]
                ticks = filter(v -> lo - pad <= log10(v) <= hi + pad, steps)
                ticks = ticks[1:cld(length(ticks), 4):end]
            end
            axq.xticks = log10.(ticks)
        end
    end
    for (r, t) in enumerate(row_labels)
        CairoMakie.Label(
            fig[r + 1, 0], t; rotation = pi / 2, font = :bold,
            tellheight = false
        )
    end
    return fig
end

## Tick labels for the recovery panels: thousands separated, small values to
## two significant figures.
function _recovery_tick(v::Real)
    a = abs(v)
    a >= 1000 && return replace(
        string(round(Int, v)), r"(?<=\d)(?=(\d{3})+$)" => ","
    )
    a >= 10 && return string(round(Int, v))
    return replace(string(round(v; sigdigits = 2)), r"\.0$" => "")
end
