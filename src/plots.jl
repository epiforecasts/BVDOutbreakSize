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
        xmax::Union{Nothing, Real} = nothing)
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
           AoG.mapping(:C_T => "Cumulative cases C_T",
               color = :stream => "Data stream") *
           AoG.AlgebraOfGraphics.density() *
           AoG.subvisual(:line, linewidth = 2)
    fg = AoG.draw(spec;
        axis = (; ylabel = "Posterior density",
            title = "Posterior C_T by data stream",
            limits = ((0, upper), nothing)),
        figure = (; size = (760, 420))
    )

    scenario_xs = Float64[val for (_, val) in scenarios if val < upper]
    isempty(scenario_xs) || vlines!(fg.figure.content[1], scenario_xs;
        color = (:grey, 0.4), linestyle = :dash)
    return fg
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
"""
function plot_pair(chn, params::AbstractVector{Symbol};
        thin::Integer = 2, prior = nothing)
    _table(c) = DataFrame(
        NamedTuple(p => _draws(c, p) for p in params))[1:thin:end, :]
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
Horizontal point-and-interval comparison of cumulative-case estimates
from several sources. `rows` is a vector of
`(label, central, lower, upper)` tuples, drawn top to bottom with the
central estimate as a point and `[lower, upper]` as a bar. Use it to
place model posteriors next to published point estimates and their
intervals.
"""
function plot_estimate_comparison(
        rows::AbstractVector;
        xlabel::AbstractString = "Cumulative cases C(T)",
        xmax::Union{Nothing, Real} = nothing)
    n = length(rows)
    labels = [String(r[1]) for r in rows]
    central = [float(r[2]) for r in rows]
    lo = [float(r[3]) for r in rows]
    hi = [float(r[4]) for r in rows]
    top = isnothing(xmax) ? maximum(hi) * 1.08 : xmax

    fig = Figure(; size = (840, 120 + 46n))
    ax = Axis(fig[1, 1];
        xlabel = xlabel,
        yticks = (collect(1:n), reverse(labels)),
        limits = ((0, top), (0.5, n + 0.5))
    )
    for i in 1:n
        y = n - i + 1
        lines!(ax, [lo[i], hi[i]], [y, y];
            color = (:steelblue, 0.8), linewidth = 3)
        scatter!(ax, [central[i]], [y];
            color = :firebrick, markersize = 12)
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
left panel is the posterior density of the outbreak start date,
obtained by rescaling the days-since-seeding `T` to a calendar date
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
        title = "Implied start of sustained transmission",
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
([`sigmoid_ramp`](@ref)) anchored at the WHO-response `breakpoint`.

Each draw's `Rt` is masked to days at or after its own establishment day
(`n - round(T)`, where cumulative infections cross one), so the median
and 50/90% ribbons cover only the established epidemic; the
pre-establishment seeding window is prior-driven and is left unshaded and
annotated. The cut-off and the intervention breakpoint are marked.
`seeding` is the calendar date of grid day 1 (so day `d` is
`seeding + (d - 1)`).
"""
function plot_rt(chn; n::Integer, breakpoint::Real,
        as_of_date::AbstractString, seeding::Date,
        rt_start::Integer = 1, week::Integer = 7, ramp::Real = 14.0)
    log_R0 = _draws(chn, Symbol("rt_state.log_R0"))
    log_R0_seed = _draws(chn, Symbol("rt_state.log_R0_seed"))
    sigma = _draws(chn, Symbol("rt_state.sigma_rw"))
    effect = _draws(chn, Symbol("rt_state.intervention_effect"))
    T_draws = _draws(chn, :T)
    ## `rt_state.z` is vector-valued: one standard-normal innovation vector
    ## per draw. Pull each draw's full vector from the chain slice.
    zmat = chn[Symbol("rt_state.z")]
    zrows = [collect(z) for z in vec(collect(zmat))]

    days = knot_days(n; week, start = rt_start)
    nb = length(days)
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
        ## Pre-bound days take the low seeding R_t, not the genetic walk base.
        log_Rt = [d < rt_start ? log_R0_seed[i] : walk[d] for d in 1:n] .+
                 effect[i] .* ramp_shape
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
    ## Median establishment day for the unshaded seeding-window annotation.
    est_start = isempty(est) ? n : minimum(est)

    epoch = date2epochdays(seeding)
    x = [epoch + (d - 1) for d in 1:n]
    xe = x[est]
    fig = Figure(; size = (900, 440))
    ax = Axis(fig[1, 1]; xlabel = "Date", ylabel = "Reproduction number Rt",
        title = "Estimated Rt over the established outbreak",
        xticklabelrotation = pi / 6)
    ## Shade the pre-establishment seeding window (prior-driven).
    est_start > 1 && vspan!(ax, x[1], x[est_start];
        color = (:grey, 0.12))
    band!(ax, xe, Float64[lo90[d] for d in est], Float64[hi90[d] for d in est];
        color = (:purple, 0.15))
    band!(ax, xe, Float64[lo50[d] for d in est], Float64[hi50[d] for d in est];
        color = (:purple, 0.30))
    lines!(ax, xe, Float64[med[d] for d in est];
        color = :purple, linewidth = 2)
    hlines!(ax, [1.0]; color = :black, linestyle = :dot)
    vlines!(ax, [Float64(epoch + breakpoint - 1)];
        color = :firebrick, linestyle = :dash, linewidth = 2)
    vlines!(ax, [Float64(date2epochdays(Date(as_of_date)))];
        color = :grey, linestyle = :dash)
    lo = floor(Int, minimum(x))
    hi = ceil(Int, maximum(x))
    ax.xticks = collect(lo:14:hi)
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

"""
One-week-ahead forecast figure from [`forecast_reported`](@ref): new
confirmed cases, new confirmed deaths, new latent infections, and the
reproduction number carried over the horizon. The count panels are
predictive-frequency histograms (new confirmed cases and deaths); the
infections panel histograms the projected new latent infections; the
reproduction-number panel shows the posterior of the forecast `R_t` with
the no-growth line at one marked. Confirmed panels are drawn only when the
forecast carries the laboratory streams.
"""
function plot_forecast(fc::DataFrame)
    count_cols = Tuple{Symbol, String, Symbol}[]
    :confirmed_new in propertynames(fc) && push!(count_cols,
        (:confirmed_new, "New confirmed cases (DRC)", :goldenrod))
    :confirmed_deaths_new in propertynames(fc) && push!(count_cols,
        (:confirmed_deaths_new, "New confirmed deaths (DRC)", :darkorange3))
    push!(count_cols, (:infections_new, "New infections (DRC)", :steelblue))
    npanels = length(count_cols) + 1
    ncols = min(npanels, 2)
    nrows = cld(npanels, ncols)
    fig = Figure(; size = (400 * ncols, 360 * nrows))
    for (i, (col, title, colour)) in enumerate(count_cols)
        v = fc[!, col]
        upper = max(1.0, quantile(v, 0.995))
        r, c = cld(i, ncols), mod1(i, ncols)
        ax = Axis(fig[r, c];
            xlabel = title, ylabel = "Predictive frequency",
            title = "One week ahead", limits = ((0, upper), nothing))
        hist!(ax, v; bins = range(0, upper; length = 30),
            color = (colour, 0.7))
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
Validation figure for a [`forecast_reported`](@ref) projection, laid out
as a 2×3 grid. The top row shows the cumulative forecast distribution per
stream (DRC reported cases, DRC deaths, Uganda exports); the bottom row
shows the new counts forecast over the horizon, mirroring the
one-week-ahead forecast. Each panel is a histogram with the 90%
predictive interval shaded and the later-observed count drawn as a dashed
black rule. `cases`, `deaths` and `exports` are the observed cumulative
counts; `baseline_*` are the counts at the forecast origin, so the
observed new count is the cumulative truth minus the baseline.
"""
function plot_forecast_vs_truth(fc::DataFrame;
        cases::Real, deaths::Real, exports::Real,
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
            float(deaths), float(deaths) - float(baseline_deaths)),
        (:exports_cum, :exports_new, "exports (Uganda)", :seagreen,
            float(exports), float(exports) - float(baseline_exports))
    ])
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
Per-vintage conditional one-step-ahead posterior-predictive for the DRC
streams. For each `panel` the predicted cumulative count at vintage `v`
conditions on the *observed* cumulative at the previous vintage and adds
only the posterior-predictive between-vintage increment,
``\\hat{y}_v = y_{v-1} + \\Delta_v`` with ``y_0 = 0``. The increment
``\\Delta_v`` is the per-bin replicate draw the model already samples in
predictive mode, so each step carries the full posterior uncertainty of
the *new* increment while anchoring on what was actually observed at the
preceding sitrep. This is the "filtered" one-step-ahead predictive: it
isolates the model's per-step increment prediction and does not let
errors compound across the series, unlike a running sum of the modelled
increments. The result is summarised by vintage as a median line with
shaded 50% and 90% predictive ribbons; the observed cumulative counts
are overlaid as points. Each `panel` is a `NamedTuple`
`(; title, dates, replicates, observed)`, where `replicates` is a vector
of per-draw increment vectors (one entry per vintage, oldest first) and
`observed` the matching observed cumulative counts used as the
conditioning baselines. `colour` is optional per panel.
"""
function plot_vintage_conditional_ppc(
        panels::AbstractVector; xlabel = "Sitrep date")
    npanels = length(panels)
    nrows = npanels > 1 ? 2 : 1
    ncols = cld(npanels, nrows)
    fig = Figure(; size = (460 * ncols, 420 * nrows))
    for (j, p) in enumerate(panels)
        row, col = cld(j, ncols), mod1(j, ncols)
        n = length(p.dates)
        colour = get(p, :colour, :steelblue)
        ## Observed cumulative at the previous vintage is the conditioning
        ## baseline for each step (`y_0 = 0`); `obs_prev[v]` is `y_{v-1}`.
        obs_cum = float.(p.observed)
        obs_prev = [v == 1 ? 0.0 : obs_cum[v - 1] for v in 1:n]
        ## `replicates` may arrive as a draws×chains matrix of per-bin
        ## vectors (FlexiChains slice); flatten to one vector of draws.
        ## Each draw's conditional cumulative at vintage `v` is the
        ## observed previous cumulative plus the drawn increment `Δ_v`.
        cond = [obs_prev .+ collect(r) for r in vec(collect(p.replicates))]
        q(i, pr) = quantile([c[i] for c in cond], pr)
        med = [q(i, 0.5) for i in 1:n]
        lo90 = [q(i, 0.05) for i in 1:n]
        hi90 = [q(i, 0.95) for i in 1:n]
        lo50 = [q(i, 0.25) for i in 1:n]
        hi50 = [q(i, 0.75) for i in 1:n]
        x = collect(1:n)
        ## Truncate the y-axis to a robust ceiling driven by the observed
        ## counts and the 50% band, so a heavy upper tail in any one stream
        ## (the 90% band can run far above the data) does not flatten the
        ## visible detail; the band clips at the axis limit.
        yupper = 1.6 * max(maximum(obs_cum), maximum(hi50), 1.0)
        ax = Axis(fig[row, col]; title = p.title, xlabel = xlabel,
            ylabel = col == 1 ? "Cumulative count" : "",
            xticks = (x, string.(p.dates)),
            xticklabelrotation = pi / 4, xticklabelsize = 9,
            limits = (nothing, (0, yupper)))
        band!(ax, x, lo90, hi90; color = (colour, 0.15))
        band!(ax, x, lo50, hi50; color = (colour, 0.30))
        lines!(ax, x, med; color = colour, linewidth = 2)
        scatter!(ax, x, float.(p.observed); color = :black, markersize = 9)
    end
    return fig
end
