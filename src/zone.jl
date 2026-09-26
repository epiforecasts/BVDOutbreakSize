# Post-processing of a health-zone chain (`bvd_zone`): daily trajectories
# rebuilt from the stored knots, the zone forecast drawn from the model with
# `predict` and its archive, the overview and forecast tables, scoring
# against the observed split, and the fit diagnostics. Every draw carries
# its own draw of the shared quantity, so the province model's uncertainty
# is already inside each zone draw and nothing here pairs the two stages
# after the fact. Everything but the forecast reads the chain and the
# `zone_fit_inputs` it was fitted with.

## One draw's state from the chain: the deviation knots `(n_zones ×
## n_knots)`, the initial shares, the AR retention, the draw of the shared
## quantity and the patch trajectory it implies, and, with mixing, the
## per-patch mixing fractions. A chain with no kept meld cell carries an
## empty `η` and every draw reads the province model's mean curve.
function _zone_states(chn, inputs; week::Integer = inputs.week)
    zd = inputs.model_data
    nz = length(inputs.zone_keys)
    np = length(inputs.patch_names)
    K = length(inputs.knots)
    knots = _draw_vectors(chn, :delta_knots_zone)
    starts = _draw_vectors(chn, :share_start_zone)
    halflife = _draws(chn, :region_halflife_zone)
    eps_ = _has_key(chn, :mixing_epsilon_zone) ?
        _draw_vectors(chn, :mixing_epsilon_zone) : nothing
    eta = _has_key(chn, :parent_eta_zone) ?
        _draw_vectors(chn, :parent_eta_zone) : nothing
    ndraws = length(knots)
    isempty(knots) || length(knots[1]) == nz * K ||
        error(
        "_zone_states: `delta_knots_zone` holds $(length(knots[1])) " *
            "entries but $nz zones by $K knots is $(nz * K); the chain was " *
            "fitted to different inputs."
    )
    eta === nothing || ndraws == 0 || length(eta[1]) in (0, zd.meld_d) ||
        error(
        "_zone_states: `parent_eta_zone` holds $(length(eta[1])) " *
            "entries but the shared quantity has $(zd.meld_d); the chain " *
            "was fitted to different inputs."
    )
    scale(i) = (eta === nothing || isempty(eta[i])) ? nothing :
        zone_parent_scale(
            zd.meld_weights, zd.meld_L,
            Float64.(eta[i]), np, zd.n
        )
    return [
        (;
            δ_knots = reshape(Float64.(knots[i]), nz, K),
            w0 = Float64.(starts[i]),
            φ = exp2(-week / halflife[i]),
            ε = eps_ === nothing ? nothing : Float64.(eps_[i]),
            η = eta === nothing ? Float64[] : Float64.(eta[i]),
            def = zone_deformation(zd, scale(i)),
        )
            for i in 1:ndraws
    ]
end

## Per-zone `(ndraws × n)` matrices from a per-draw builder returning an
## `(n_days × n_zones)` matrix over `t0 … n`, with the days before `t0`
## filled by `pre(state, z)`.
function _zone_per_zone_matrices(states, inputs, build, pre)
    nz = length(inputs.zone_keys)
    n = inputs.n
    t0 = inputs.t0
    ndraws = length(states)
    out = [Matrix{Float64}(undef, ndraws, n) for _ in 1:nz]
    for (i, st) in enumerate(states)
        m = build(st)
        for z in 1:nz
            v = pre(st, z)
            for t in 1:(t0 - 1)
                out[z][i, t] = v
            end
            for t in t0:n
                out[z][i, t] = m[t - t0 + 1, z]
            end
        end
    end
    return out
end

"""
$(TYPEDSIGNATURES)

Each zone's daily share of its patch's infections, rebuilt for every draw
by re-running the share renewal ([`zone_forward`](@ref)) from the knots
and the fixed inputs. Returns one `(ndraws × n)` matrix per zone, in the
order of `inputs.zone_keys`; the days before the grid start hold the
draw's initial share, as the model does.
"""
function reconstruct_zone_shares(chn, inputs)
    zd = inputs.model_data
    states = _zone_states(chn, inputs)
    return _zone_per_zone_matrices(
        states, inputs,
        st -> zone_forward(zd, st.δ_knots, st.w0, st.ε, st.def).shares,
        (st, z) -> st.w0[z]
    )
end

"""
$(TYPEDSIGNATURES)

How the fitted draw of the shared quantity compares with the province
model's own posterior on it. The zone stage places no prior of its own on
that quantity beyond the province model's, so if the zone data say nothing
about it the whitened draw should read back as its prior, a standard
normal in every component. A component whose posterior mean is far from
zero or whose spread is well under one says the zone structure is pulling
the province model's trajectory in that week.

Returns `(; table, mean_norm_sq, dimension)`. The table has one row per
kept patch-week cell: the patch, the window's midpoint day and date, the
posterior mean and standard deviation of that component of the whitened
draw, and the province model's own posterior standard deviation of the log
weekly infections there. `mean_norm_sq` is the mean of `‖η‖²`, which is
`dimension` when the fit reproduces the prior.
"""
function zone_meld_check(chn, inputs)
    meld = inputs.meld
    out = DataFrame(
        patch = String[], midpoint = Int[], date = Date[],
        eta_mean = Float64[], eta_sd = Float64[], parent_sd = Float64[]
    )
    empty_result = (; table = out, mean_norm_sq = NaN, dimension = meld.d)
    _has_key(chn, :parent_eta_zone) || return empty_result
    draws = _draw_vectors(chn, :parent_eta_zone)
    (isempty(draws) || isempty(draws[1])) && return empty_result
    E = reduce(hcat, [Float64.(v) for v in draws])
    parent_sd = sqrt.(vec(sum(abs2, meld.L; dims = 2)))
    for c in 1:meld.d
        mid = meld.midpoints[meld.cells_week[c]]
        push!(
            out,
            (
                inputs.patch_labels[meld.cells_patch[c]], mid,
                inputs.seeding + Day(mid - 1),
                mean(view(E, c, :)), std(view(E, c, :)), parent_sd[c],
            )
        )
    end
    return (;
        table = out, mean_norm_sq = mean(vec(sum(abs2, E; dims = 1))),
        dimension = meld.d,
    )
end

## `I_z / Λ_z` on the grid days from one draw's daily shares `(n_days ×
## n_zones)` and patch infections `I_p` `(n_patches × n)`, the days before
## `t0` at the initial share, with the first grid day on which each zone's
## cumulative infections reach `rt_floor` (`n_days + 1` if never).
function _zone_rt_daily(
        shares::AbstractMatrix, w0::AbstractVector,
        I_p::AbstractMatrix, g::AbstractVector,
        patch_of_zone::AbstractVector{<:Integer}, t0::Integer,
        rt_floor::Real
    )
    nd, nz = size(shares)
    L = length(g)
    out = zeros(Float64, nd, nz)
    first_day = fill(nd + 1, nz)
    for z in 1:nz
        p = patch_of_zone[z]
        cum = w0[z] * sum(@view I_p[p, 1:(t0 - 1)])
        for j in 1:nd
            t = t0 + j - 1
            Λ = 0.0
            for s in 1:min(L, t - 1)
                u = t - s
                w = u < t0 ? w0[z] : shares[u - t0 + 1, z]
                Λ += g[s] * I_p[p, u] * w
            end
            I = I_p[p, t] * shares[j, z]
            cum += I
            cum >= rt_floor && first_day[z] > nd && (first_day[z] = j)
            out[j, z] = I / max(Λ, eps())
        end
    end
    return out, first_day
end

"""
$(TYPEDSIGNATURES)

Each zone's implied daily reproduction number `I_z(t) / Λ_z(t)`, rebuilt
per draw from the shares and that draw's own patch trajectory, so the
interval carries the province model's uncertainty as the fit saw it.
A zone is reported from the day its cumulative infections
reach `rt_floor` (the fit's `inputs.rt_floor` by default) in the median
draw, and holds `NaN` in every draw before that day, so the reported draws
are never a selected subset. Returns one `(ndraws × n)` matrix per zone.
"""
function reconstruct_zone_rt(
        chn, inputs;
        rt_floor::Real = inputs.rt_floor
    )
    zd = inputs.model_data
    states = _zone_states(chn, inputs)
    ndraws = length(states)
    nz = length(inputs.zone_keys)
    nd = zd.n - zd.t0 + 1
    out = [fill(NaN, ndraws, inputs.n) for _ in 1:nz]
    first_day = zeros(Int, ndraws, nz)
    for (i, st) in enumerate(states)
        shares = zone_forward(zd, st.δ_knots, st.w0, st.ε, st.def).shares
        r, first_day[i, :] = _zone_rt_daily(
            shares, st.w0, st.def.I_bar, zd.g,
            inputs.patch_of_zone, zd.t0, rt_floor
        )
        for z in 1:nz, j in 1:nd

            out[z][i, zd.t0 + j - 1] = r[j, z]
        end
    end
    for z in 1:nz
        start = ceil(Int, median(view(first_day, :, z)))
        out[z][:, 1:min(zd.t0 + start - 2, inputs.n)] .= NaN
    end
    return out
end

"""
$(TYPEDSIGNATURES)

Daily zone infections with the patch uncertainty carried: each draw's
share of its own draw of the patch trajectory,
`I_z = I_p^{(j)}(t) w_z^{(j)}(t)`, so the two stages are one object rather
than paired after the fact. Returns one `(ndraws × n)` matrix per zone.
"""
function zone_infections(chn, inputs)
    zd = inputs.model_data
    states = _zone_states(chn, inputs)
    shares = _zone_per_zone_matrices(
        states, inputs,
        st -> zone_forward(zd, st.δ_knots, st.w0, st.ε, st.def).shares,
        (st, z) -> st.w0[z]
    )
    ndraws = length(states)
    out = [similar(s) for s in shares]
    for z in eachindex(shares)
        p = inputs.patch_of_zone[z]
        for i in 1:ndraws
            I = states[i].def.I_bar
            for t in 1:inputs.n
                out[z][i, t] = I[p, t] * shares[z][i, t]
            end
        end
    end
    return out
end

"""
$(TYPEDSIGNATURES)

Posterior-predictive draws of the zone forecast: `predict` on
[`bvd_zone`](@ref) run past the cut-off ([`forecast_draws`](@ref)) over the
zone chain `chn`. `inputs` must carry the parent's forecast
([`zone_fit_inputs`](@ref) with `parent_forecast`), whose horizon the
forecast runs to. Each draw keeps its fitted parameters, draws the future
shared quantity from the parent's posterior conditional on the fitted one,
fresh deviation innovations for the future knots and a parent draw of the
patch totals, and splits each total over the zones by the fitted
composition.
"""
function zone_forecast(chn, inputs; seed::Integer = 20260518)
    zf = inputs.model_data.forecast
    zf === nothing && throw(
        ArgumentError(
            "zone_forecast: the inputs carry no parent forecast; build them " *
                "with `zone_fit_inputs(...; parent_forecast)`."
        )
    )
    return forecast_draws(
        bvd_zone(inputs.model_data), chn; horizon = zf.horizon, seed
    )
end

"""
$(TYPEDSIGNATURES)

The zone forecast draws `fc` ([`zone_forecast`](@ref)) per zone and patch:
the new confirmed cases over the horizon, one draw vector per zone in the
order of `inputs.zone_keys`, the paired parent patch totals, one draw
vector per patch, each draw's zone shares `(ndraws × n_zones)` and each
draw's composition concentration. Returns `(; zones, patches, shares,
kappa)`.
"""
function zone_forecast_draws(fc, inputs)
    nz = length(inputs.zone_keys)
    np = length(inputs.patch_names)
    counts = _draw_vectors(fc, :forecast_zone_confirmed)
    totals = _draw_vectors(fc, :forecast_patch_confirmed)
    shares = _draw_vectors(fc, :forecast_zone_share)
    return (;
        zones = [Float64[c[z] for c in counts] for z in 1:nz],
        patches = [Float64[t[p] for t in totals] for p in 1:np],
        shares = Float64[s[z] for s in shares, z in 1:nz],
        kappa = Float64.(_draws(fc, :forecast_zone_concentration)),
    )
end

"""
$(TYPEDSIGNATURES)

Probability that each zone reports at least `K` new confirmed cases over
the horizon, for every `K` in `thresholds`. Under the fitted composition a
zone's count given its patch total `N` is Beta-binomial, the marginal of
the Dirichlet-multinomial with concentration `κ π`, so over the forecast
draws `i` ([`zone_forecast_draws`](@ref)), each with its patch total `N_i`,
zone share `π_i` and concentration `κ_i`:

```math
P(y_z ≥ K) = \\frac{1}{n}\\sum_i \\Bigl[1 − F_{\\mathrm{BB}(N_i,\\, κ_i π_{z,i},\\, κ_i (1 − π_{z,i}))}(K − 1)\\Bigr].
```

Returns an `(n_zones × length(thresholds))` matrix.
"""
function zone_forecast_probabilities(
        fc, inputs;
        thresholds = (1, 5, 10, 20),
        draws = zone_forecast_draws(fc, inputs)
    )
    nz = length(inputs.zone_keys)
    out = zeros(Float64, nz, length(thresholds))
    nd = length(draws.kappa)
    for (p, zs) in enumerate(inputs.patch_ranges), z in zs
        for (k, K) in enumerate(thresholds)
            acc = 0.0
            for i in 1:nd
                N = round(Int, draws.patches[p][i])
                N < K && continue
                κ = draws.kappa[i]
                α = max(κ * draws.shares[i, z], eps())
                β = max(κ - α, eps())
                acc += 1 - cdf(BetaBinomial(N, α, β), K - 1)
            end
            out[z, k] = acc / nd
        end
    end
    return out
end

"""
$(TYPEDSIGNATURES)

Confirmed cases allocated to each zone over the `window` days ending at the
fit's cut-off, from the fitted count matrix and vintage days, in the order
of `inputs.zone_keys`.
"""
zone_recent_cases(inputs; window::Integer = 7) =
    _zone_recent_counts(inputs; window)

"""
$(TYPEDSIGNATURES)

Date of the last vintage on which each zone's allocated confirmed count
rose, or `missing` for a zone with no allocated case, in the order of
`inputs.zone_keys`.
"""
function zone_last_case_dates(inputs)
    return [
        let vs = findall(>(0), @view inputs.counts[z, :])
            isempty(vs) ? missing : inputs.dates[vs[end]]
        end
            for z in eachindex(inputs.zone_keys)
    ]
end

"""
$(TYPEDSIGNATURES)

Long-format archive of the zone forecast draws `fc` ([`zone_forecast`](@ref))
made from the cut-off `made_date`, in the
[`province_forecast_archive`](@ref) schema plus a `zone` column. `province`
is the patch key and `zone` the manifest's dotted `province.zone` key, and
each value is one draw of the zone's new confirmed cases over the forecast
horizon, under the `confirmed cases` stream label. `thin` keeps every
`thin`-th draw.
"""
function zone_forecast_archive(
        fc, inputs; made_date::Date, thin::Integer = 1
    )
    out = DataFrame(
        made_date = Date[], horizon = Int[], target_date = Date[],
        province = String[], zone = String[], stream = String[],
        draw = Int[], value = Float64[]
    )
    h = inputs.model_data.forecast.horizon
    target = made_date + Day(h)
    draws = zone_forecast_draws(fc, inputs)
    for z in eachindex(inputs.zone_keys)
        vals = draws.zones[z]
        prov = inputs.patch_names[inputs.patch_of_zone[z]]
        for (d, i) in enumerate(1:thin:length(vals))
            push!(
                out, (
                    made_date, h, target, prov, inputs.zone_keys[z],
                    "confirmed cases", d, vals[i],
                )
            )
        end
    end
    return out
end

## Per-zone draws of a vector deterministic, one vector per zone.
function _zone_draws(chn, key::Symbol, nz::Integer)
    vs = _draw_vectors(chn, key)
    return [[Float64(v[z]) for v in vs] for z in 1:nz]
end

"""
$(TYPEDSIGNATURES)

One row per health zone: its patch, the confirmed cases to date, the share
of its patch's infections at the cut-off with a 90% interval, the implied
reproduction number at the cut-off with its 90% interval and the posterior
probability that it exceeds one, the log-transmission deviation at the
cut-off with a 90% interval and whether the zone carries a time-varying
deviation. The
reproduction number is also given numerically as `rt_median`, `rt_lo90`
and `rt_hi90` (`NaN` below the reporting floor), the probability unrounded
as `p_rt_above_one`, and the patch as its index `patch_index`. It is
rebuilt by [`reconstruct_zone_rt`](@ref). Walking zones sort first, in
decreasing order of
the probability that the reproduction number exceeds one with those below
the reporting floor last among them, then the level-only zones in the
same order.
"""
function zone_overview_table(chn, inputs; digits::Integer = 2)
    nz = length(inputs.zone_keys)
    rt = reconstruct_zone_rt(chn, inputs)
    R = [m[:, inputs.n] for m in rt]
    share = _zone_draws(chn, :share_T_zone, nz)
    δ = _zone_draws(chn, :delta_T_zone, nz)
    rows = NamedTuple[]
    for z in 1:nz
        r = filter(isfinite, R[z])
        p_above = isempty(r) ? NaN : mean(r .> 1)
        push!(
            rows,
            (
                zone = inputs.zone_labels[z],
                patch = inputs.patch_labels[inputs.patch_of_zone[z]],
                cases = inputs.cumulative[z],
                share = _median_ci(100 .* share[z]; digits = 1),
                R_T = isempty(r) ? "" : _median_ci(r; digits),
                rt_median = isempty(r) ? NaN : median(r),
                rt_lo90 = isempty(r) ? NaN : quantile(r, 0.05),
                rt_hi90 = isempty(r) ? NaN : quantile(r, 0.95),
                p_R_above_1 = isnan(p_above) ? NaN : round(p_above; digits),
                p_rt_above_one = p_above,
                patch_index = inputs.patch_of_zone[z],
                delta_T = _median_ci(δ[z]; digits),
                walking = inputs.walking[z],
            )
        )
    end
    order = sortperm(
        [
            (r.walking ? 2.0 : 0.0) +
                (isnan(r.p_R_above_1) ? -1.0 : r.p_R_above_1)
                for r in rows
        ]; rev = true
    )
    return DataFrame(rows[order])
end

"""
$(TYPEDSIGNATURES)

Per-zone new confirmed cases over the forecast horizon from the zone
forecast draws `fc` ([`zone_forecast`](@ref)), as the median and the
90/60/30% intervals the other forecast tables report, then the probability
of at least `K` cases for each `K` in `thresholds`
([`zone_forecast_probabilities`](@ref)) as `p_ge_K`, one row per zone with
a patch-total row per patch (its threshold columns `missing`).
"""
function zone_forecast_table(
        fc, inputs; digits::Integer = 0, thresholds = (1, 5, 10, 20)
    )
    draws = zone_forecast_draws(fc, inputs)
    probs = zone_forecast_probabilities(fc, inputs; thresholds, draws)
    pcols = [Symbol("p_ge_", K) for K in thresholds]
    row(label, patch, v, pr) = begin
        s = posterior_summary(v)
        merge(
            (
                zone = label, patch = patch,
                median = round(median(v); digits),
                lower_90 = round(s.lo90; digits),
                lower_60 = round(s.lo60; digits),
                lower_30 = round(s.lo30; digits),
                upper_30 = round(s.hi30; digits),
                upper_60 = round(s.hi60; digits),
                upper_90 = round(s.hi90; digits),
            ),
            NamedTuple{Tuple(pcols)}(Tuple(pr))
        )
    end
    rows = NamedTuple[]
    for (p, zs) in enumerate(inputs.patch_ranges)
        isempty(zs) && continue
        label = inputs.patch_labels[p]
        for z in zs
            push!(
                rows, row(
                    inputs.zone_labels[z], label, draws.zones[z],
                    round.(probs[z, :]; digits = 2)
                )
            )
        end
        push!(
            rows, row(
                "Patch total", label, draws.patches[p],
                fill(missing, length(thresholds))
            )
        )
    end
    return DataFrame(rows)
end

"""
$(TYPEDSIGNATURES)

The observed new confirmed cases per zone over `(made_date, made_date +
horizon]`, from the zone histories in `obs` (a later vintage than the fit's),
in the order of `inputs.zone_keys`: the cumulative at the vintage on the
target date minus the cumulative at the last vintage on or before
`made_date`, clamped at zero. A zone whose history has no vintage on the
target date itself is `missing`, since the nearest vintage either side
would count a different window; [`zone_forecast_vs_truth`](@ref) shows
such a zone unscored and [`zone_forecast_scores`](@ref) skips its patch.
"""
function zone_forecast_truth(
        obs, inputs; made_date::Date,
        horizon::Integer = 7
    )
    hist = obs.zone_confirmed_history
    target = made_date + Day(horizon)
    dates(h) = [obs.seeding + Day(d - 1) for d in h.days]
    at(h, date) = begin
        idx = findlast(<=(date), dates(h))
        idx === nothing ? 0 : Int(h.counts[idx])
    end
    on(h, date) = begin
        idx = findfirst(==(date), dates(h))
        idx === nothing ? nothing : Int(h.counts[idx])
    end
    out = Vector{Union{Missing, Int}}(undef, length(inputs.zone_keys))
    for z in eachindex(inputs.zone_keys)
        prov = inputs.zone_province[z]
        zone = inputs.zone_names[z]
        h = haskey(hist, prov) && haskey(hist[prov], zone) ?
            hist[prov][zone] : nothing
        final = h === nothing ? nothing : on(h, target)
        out[z] = final === nothing ? missing :
            max(final - at(h, made_date), 0)
    end
    return out
end

"""
$(TYPEDSIGNATURES)

Per-zone forecast against what was observed: the median, 50% and 90%
intervals of the zone forecast draws `fc` ([`zone_forecast`](@ref)), the
observed count from `truth` ([`zone_forecast_truth`](@ref)) and whether it
fell inside the interval, one row per zone plus a patch-total row.
"""
function zone_forecast_vs_truth(
        fc, inputs; truth::AbstractVector, digits::Integer = 0
    )
    draws = zone_forecast_draws(fc, inputs)
    row(label, patch, v, obs) = begin
        s = posterior_summary(v)
        (
            zone = label, patch = patch,
            central_estimate = round(median(v); digits),
            lower_90 = round(s.lo90; digits),
            lower_50 = round(quantile(v, 0.25); digits),
            upper_50 = round(quantile(v, 0.75); digits),
            upper_90 = round(s.hi90; digits),
            observed = obs,
            within_90 = ismissing(obs) ? missing : s.lo90 <= obs <= s.hi90,
        )
    end
    rows = NamedTuple[]
    for (p, zs) in enumerate(inputs.patch_ranges)
        isempty(zs) && continue
        label = inputs.patch_labels[p]
        for z in zs
            push!(
                rows, row(
                    inputs.zone_labels[z], label, draws.zones[z],
                    truth[z]
                )
            )
        end
        tot = any(ismissing, truth[zs]) ? missing : sum(truth[zs])
        push!(rows, row("Patch total", label, draws.patches[p], tot))
    end
    return DataFrame(rows)
end

## Multinomial log mass of the split `y` with probabilities `π`.
function _multinomial_logpmf(y::AbstractVector{<:Integer}, π::AbstractVector)
    N = sum(y)
    lp = loggamma(N + 1.0)
    for i in eachindex(y)
        lp -= loggamma(y[i] + 1.0)
        y[i] == 0 && continue
        lp += y[i] * log(max(π[i], eps()))
    end
    return lp
end

## Log of the mean over draws of the multinomial mass, each draw's split
## probabilities in a row of `Π`.
function _mixture_multinomial_logscore(
        y::AbstractVector{<:Integer},
        Π::AbstractMatrix
    )
    lps = [_multinomial_logpmf(y, view(Π, i, :)) for i in 1:size(Π, 1)]
    m = maximum(lps)
    return m + log(mean(exp.(lps .- m)))
end

## Observed increments per zone over the `window` days ending at the fit's
## cut-off, from the fitted count matrix and vintage days.
function _zone_recent_counts(inputs; window::Integer = 7)
    nz = length(inputs.zone_keys)
    out = zeros(Int, nz)
    for (v, d) in enumerate(inputs.days)
        d > inputs.n - window || continue
        out .+= inputs.counts[:, v]
    end
    return out
end

"""
$(TYPEDSIGNATURES)

Scores of the zone forecast draws `fc` ([`zone_forecast`](@ref)) against
the observed split `truth` ([`zone_forecast_truth`](@ref)). Per patch, the multinomial log score of
the observed zone split of the patch's weekly total under the zone model
(the log of the draw-averaged multinomial mass at the projected shares),
against two nulls: share persistence, the observed cumulative zone shares
at the cut-off, and naive persistence, the zone split of the last
`window` days before the cut-off. Both nulls carry a pseudo-count of 0.5
per zone. Patches whose observed total is zero or whose truth is
incomplete are skipped. The patch total is scored alongside the split by
[`score_draws`](@ref), the zone model's total against a naive persistence
total of the last `window` days with negative-binomial noise at dispersion
`k_baseline`. Returns a `DataFrame` with `patch`, `method`, `observed` and
`log_score`, then the [`score_draws`](@ref) columns `crps`, `log_crps`,
`dispersion`, `overprediction`, `underprediction`, `coverage_50`,
`coverage_90`, `bias` and `n`. Share persistence forecasts no total, so
its score columns are `missing`.
"""
function zone_forecast_scores(
        fc, inputs;
        truth::AbstractVector, window::Integer = 7, k_baseline::Real = 5.0,
        rng::AbstractRNG = MersenneTwister(20260518)
    )
    draws = zone_forecast_draws(fc, inputs)
    recent = _zone_recent_counts(inputs; window)
    rows = NamedTuple[]
    for (p, zs) in enumerate(inputs.patch_ranges)
        isempty(zs) && continue
        any(ismissing, truth[zs]) && continue
        y = Int[truth[z] for z in zs]
        N = sum(y)
        label = inputs.patch_labels[p]
        if N > 0
            Π = draws.shares[:, zs]
            cum = Float64[inputs.cumulative[z] + 0.5 for z in zs]
            last_ = Float64[recent[z] + 0.5 for z in zs]
            row(method, log_score, total) = merge(
                (; patch = label, method, observed = N, log_score),
                total === nothing ? _zone_no_total_scores() :
                    score_draws(N, total)
            )
            push!(
                rows, row(
                    "zone model", _mixture_multinomial_logscore(y, Π),
                    draws.patches[p]
                )
            )
            ## Share persistence forecasts the split alone, so it has no
            ## total to score.
            push!(
                rows, row(
                    "share persistence",
                    _multinomial_logpmf(y, cum ./ sum(cum)), nothing
                )
            )
            naive_tot = sum(recent[zs])
            naive = _zone_nb_draws(
                rng, k_baseline,
                max(float(naive_tot), 0.5), length(draws.patches[p])
            )
            push!(
                rows, row(
                    "naive persistence",
                    _multinomial_logpmf(y, last_ ./ sum(last_)), naive
                )
            )
        end
    end
    return DataFrame(rows)
end

## The [`score_draws`](@ref) columns for a rule that forecasts no total,
## so every row of the table carries the same columns.
function _zone_no_total_scores()
    keys_ = keys(score_draws(0, Float64[]))
    return NamedTuple{keys_}(ntuple(_ -> missing, length(keys_)))
end

## Negative-binomial draws with mean `μ` and dispersion `k`, `m` of them.
function _zone_nb_draws(rng::AbstractRNG, k::Real, μ::Real, m::Integer)
    p = k / (k + μ)
    d = NegativeBinomial(k, p)
    return Float64[rand(rng, d) for _ in 1:m]
end

## The two compositions `bvd_zone` scores, as the names the render reads
## them by: the `model_data` field holding the observed counts, the
## per-draw relative multiplier the allocated rates carry, the chain's own
## concentration and the word the figures label the stream with.
const _ZONE_COMPOSITIONS = (
    cases = (;
        counts = :counts, multiplier = :zone_ascertainment_relative,
        rho = :composition_rho_zone, word = "cases",
    ),
    deaths = (;
        counts = :death_counts, multiplier = :zone_severity_relative,
        rho = :composition_rho_death_zone, word = "deaths",
    ),
)

## The composition entry of `stream`, `:cases` or `:deaths`.
function _zone_composition(stream::Symbol)
    haskey(_ZONE_COMPOSITIONS, stream) || error(
        "unknown zone composition '$stream'; it is `:cases` or `:deaths`."
    )
    return _ZONE_COMPOSITIONS[stream]
end

## Observed zone counts of a composition, `(n_zones × n_vintages)`.
function _zone_observed_counts(inputs, stream::Symbol)
    return getproperty(inputs.model_data, _zone_composition(stream).counts)
end

## One draw's allocated rates per zone and vintage. The deaths run the same
## infections through the infection-to-confirmed-death delay rather than
## the case delay, as the model does.
function _zone_allocated_increments(zd, st, stream::Symbol)
    fw = zone_forward(zd, st.δ_knots, st.w0, st.ε, st.def)
    stream === :cases && return fw.increments
    daily = zd.death_matrix * fw.infections .+
        st.def.death_pre_rows .* transpose(st.w0)
    return zone_report_increments(
        daily, st.w0, zd.patch_ranges, zd.death_days, zd.t0,
        st.def.death_pre_cum
    )
end

## Per-draw relative multiplier of a composition, one vector per draw, or
## `nothing` for a chain that does not carry it. It is centred within a
## patch rather than constant, so it does not cancel in the share.
function _zone_multiplier_draws(chn, stream::Symbol)
    key = _zone_composition(stream).multiplier
    _has_key(chn, key) || return nothing
    return _draw_vectors(chn, key)
end

## Modelled share of each zone in its patch's allocated total per vintage
## for every draw of a chain, `(ndraws × n_zones × n_vintages)`. `stream`
## picks the composition, `:cases` or `:deaths`.
function _zone_modelled_shares(chn, inputs; stream::Symbol = :cases)
    zd = inputs.model_data
    states = _zone_states(chn, inputs)
    nz = length(inputs.zone_keys)
    nv = length(inputs.days)
    mult = _zone_multiplier_draws(chn, stream)
    out = zeros(Float64, length(states), nz, nv)
    for (i, st) in enumerate(states)
        inc = _zone_allocated_increments(zd, st, stream)
        rate(z, v) = safe_rate(inc[z, v]) *
            (mult === nothing ? 1.0 : Float64(mult[i][z]))
        for zs in inputs.patch_ranges, v in 1:nv

            isempty(zs) && continue
            tot = max(sum(rate(z, v) for z in zs), eps())
            for z in zs
                out[i, z, v] = rate(z, v) / tot
            end
        end
    end
    return out
end

## Zone totals a composition's panels rank by: the manifest cumulative for
## the cases and the summed allocated counts for the deaths.
function _zone_stream_totals(inputs, stream::Symbol)
    stream === :cases && return inputs.cumulative
    return vec(sum(_zone_observed_counts(inputs, stream); dims = 2))
end

## The `top` zones of one patch, largest `score` first.
function _zone_largest_in(zs, score, top::Integer)
    order = sort(collect(zs); by = z -> -score[z])
    return order[1:min(top, length(order))]
end

## The `top` zones by `score` in each patch, patch by patch.
function _zone_largest(inputs, top::Integer, score = inputs.cumulative)
    out = Int[]
    for zs in inputs.patch_ranges
        append!(out, _zone_largest_in(zs, score, top))
    end
    return out
end

"""
$(TYPEDSIGNATURES)

Composition posterior predictive check: the observed and modelled share of
each patch's allocated counts per vintage, for the `top` zones with most
of them in each patch. The modelled share is the median
and 90% interval over draws of the expected increments normalised within
the patch. With `prior_chain` (draws from the model prior, `sample(model,
Prior(), n)`) the same summaries of the prior predictive are added as
`prior_lower_90`, `prior_median` and `prior_upper_90`. `stream` picks the
composition, `:cases` or `:deaths`. Returns a long
`DataFrame` with
columns `patch`, `zone`, `date`, `day`, `observed`, `observed_share`,
`total`, `lower_90`, `median` and `upper_90`; vintages whose allocated
patch total is zero are omitted.
"""
function zone_composition_ppc(
        chn, inputs; top::Integer = 10,
        prior_chain = nothing, stream::Symbol = :cases
    )
    post = _zone_modelled_shares(chn, inputs; stream)
    prior = prior_chain === nothing ? nothing :
        _zone_modelled_shares(prior_chain, inputs; stream)
    counts = _zone_observed_counts(inputs, stream)
    largest = _zone_largest(inputs, top, _zone_stream_totals(inputs, stream))
    rows = NamedTuple[]
    for z in largest, v in eachindex(inputs.days)

        p = inputs.patch_of_zone[z]
        zs = inputs.patch_ranges[p]
        N = sum(@view counts[zs, v])
        N > 0 || continue
        m = view(post, :, z, v)
        row = (
            patch = inputs.patch_labels[p], zone = inputs.zone_labels[z],
            date = inputs.dates[v], day = inputs.days[v],
            observed = counts[z, v],
            observed_share = counts[z, v] / N, total = N,
            lower_90 = quantile(m, 0.05), median = median(m),
            upper_90 = quantile(m, 0.95),
        )
        if prior !== nothing
            q = view(prior, :, z, v)
            row = merge(
                row,
                (;
                    prior_lower_90 = quantile(q, 0.05),
                    prior_median = median(q),
                    prior_upper_90 = quantile(q, 0.95),
                )
            )
        end
        push!(rows, row)
    end
    return DataFrame(rows)
end

"""
$(TYPEDSIGNATURES)

One figure per patch of the composition check: for the `top` zones with
most of the stream's counts, the share band of `chn` (5–95% and 25–75%,
named `label` in
the title), the prior predictive band from `prior_chain` when given, and
the observed share per vintage as points. `stream` picks the composition,
`:cases` or `:deaths`. Returns a vector of figures in
patch order, one per patch with zones.
"""
function plot_zone_composition_ppc(
        chn, inputs; prior_chain = nothing,
        top::Integer = 6, ncols::Integer = 3,
        label::AbstractString = "posterior predictive",
        stream::Symbol = :cases
    )
    word = _zone_composition(stream).word
    post = _zone_modelled_shares(chn, inputs; stream)
    prior = prior_chain === nothing ? nothing :
        _zone_modelled_shares(prior_chain, inputs; stream)
    counts = _zone_observed_counts(inputs, stream)
    score = _zone_stream_totals(inputs, stream)
    x = Float64[date2epochdays(d) for d in inputs.dates]
    figs = Figure[]
    for (p, zs) in enumerate(inputs.patch_ranges)
        isempty(zs) && continue
        zones = _zone_largest_in(zs, score, top)
        nrows = cld(length(zones), ncols)
        fig = Figure(; size = (360 * ncols, 240 * nrows + 60))
        for (k, z) in enumerate(zones)
            r = cld(k, ncols)
            c = k - (r - 1) * ncols
            ax = Axis(
                fig[r, c]; title = inputs.zone_labels[z],
                xlabel = "Date", ylabel = "Share of patch " * word,
                xticklabelrotation = pi / 6
            )
            _share_bands!(ax, x, prior, z, :grey50)
            _share_bands!(ax, x, post, z, :steelblue)
            obs = [
                let N = sum(@view counts[zs, v])
                    N > 0 ? counts[z, v] / N : NaN
                end
                    for v in eachindex(inputs.days)
            ]
            keep = .!isnan.(obs)
            scatter!(ax, x[keep], obs[keep]; color = :black, markersize = 5)
        end
        _date_ticks!(fig, nrows, ncols, length(zones), x)
        CairoMakie.Label(
            fig[0, 1:ncols],
            inputs.patch_labels[p] * ": observed (points) and " * label *
                " (blue)" * (prior === nothing ? "" : ", prior predictive (grey)") *
                " zone " * word * " shares"; fontsize = 14, font = :bold
        )
        push!(figs, fig)
    end
    return figs
end

## 5–95% and 25–75% bands of `shares[:, z, :]` over the vintages.
function _share_bands!(ax, x, shares, z, colour)
    shares === nothing && return nothing
    q(f) = [f(view(shares, :, z, v)) for v in 1:size(shares, 3)]
    band!(
        ax, x, q(v -> quantile(v, 0.05)), q(v -> quantile(v, 0.95));
        color = (colour, 0.2)
    )
    band!(
        ax, x, q(v -> quantile(v, 0.25)), q(v -> quantile(v, 0.75));
        color = (colour, 0.35)
    )
    lines!(ax, x, q(median); color = colour)
    return nothing
end

## Date tick labels on every axis of a grid of `n` panels.
function _date_ticks!(fig, nrows, ncols, n, x)
    lo, hi = floor(Int, minimum(x)), ceil(Int, maximum(x))
    ticks = collect(lo:14:hi)
    labels = [string(epochdays2date(t)) for t in ticks]
    for k in 1:n
        r = cld(k, ncols)
        c = k - (r - 1) * ncols
        ax = fig.content[k]
        ax isa Axis || continue
        ax.xticks = (ticks, labels)
    end
    return nothing
end

"""
$(TYPEDSIGNATURES)

Per-draw composition draws for the posterior predictive check, one
`(ndraws × n_vintages)` matrix per zone in chain order: `expected`, the
modelled share of the patch's allocated total; `predictive`, counts drawn
per draw and vintage from the Dirichlet-multinomial on the observed
allocated total at that draw's `ρ`; and `predictive_share`, those counts
over the total. `observed` is the `(n_zones × n_vintages)` observed share.
With `cumulative = true` counts accumulate over the vintages and the
shares are of the patch's cumulative allocated total. Shares are `NaN`
where the total is zero.

`stream` picks the composition: `:cases` for the per-zone confirmed cases
and `:deaths` for the per-zone confirmed deaths, each with its own
concentration and its own observed counts. Both are scored on one vintage
grid.
"""
function zone_composition_draws(
        chn, inputs; stream::Symbol = :cases, cumulative::Bool = false,
        rng::AbstractRNG = MersenneTwister(20260518)
    )
    comp = _zone_composition(stream)
    shares = _zone_modelled_shares(chn, inputs; stream)
    counts = _zone_observed_counts(inputs, stream)
    ρs = _draws(chn, comp.rho)
    ndraws, nz, nv = size(shares)
    expected = [fill(NaN, ndraws, nv) for _ in 1:nz]
    predictive = [zeros(Int, ndraws, nv) for _ in 1:nz]
    predictive_share = [fill(NaN, ndraws, nv) for _ in 1:nz]
    observed = fill(NaN, nz, nv)
    for zs in inputs.patch_ranges
        isempty(zs) && continue
        m = length(zs)
        N = [sum(@view counts[zs, v]) for v in 1:nv]
        tot = cumulative ? cumsum(N) : N
        obs = Float64.(counts[zs, :])
        cumulative && cumsum!(obs, obs; dims = 2)
        for (k, z) in enumerate(zs), v in 1:nv

            tot[v] > 0 && (observed[z, v] = obs[k, v] / tot[v])
        end
        for i in 1:ndraws
            κ = _zone_kappa(ρs[i])
            y = zeros(Int, m, nv)
            e = zeros(Float64, m, nv)
            for v in 1:nv
                N[v] > 0 || continue
                π = shares[i, zs, v]
                ## Floored as the likelihood floors its rates: a zone whose
                ## share underflows to zero would otherwise fail the draw.
                y[:, v] = rand(
                    rng, DirichletMultinomial(N[v], max.(κ .* π, eps()))
                )
                e[:, v] = N[v] .* π
            end
            if cumulative
                cumsum!(y, y; dims = 2)
                cumsum!(e, e; dims = 2)
            end
            for (k, z) in enumerate(zs), v in 1:nv

                tot[v] > 0 || continue
                expected[z][i, v] = e[k, v] / tot[v]
                predictive[z][i, v] = y[k, v]
                predictive_share[z][i, v] = y[k, v] / tot[v]
            end
        end
    end
    return (; expected, predictive, predictive_share, observed)
end

"""
$(TYPEDSIGNATURES)

Calibration of the composition per patch, in the columns of
[`stream_calibration`](@ref): over every zone and vintage with a positive
allocated patch total, the mean [`bias_sample`](@ref) of the observed
count against the predictive counts of [`zone_composition_draws`](@ref)
and the fraction of cells inside the central 50% and 90% predictive
intervals. One row per patch and an `All zones` row. `stream` picks the
composition, `:cases` or `:deaths`.
"""
function zone_composition_calibration(
        chn, inputs; stream::Symbol = :cases,
        rng::AbstractRNG = MersenneTwister(20260518)
    )
    draws = zone_composition_draws(chn, inputs; stream, rng)
    counts = _zone_observed_counts(inputs, stream)
    nv = length(inputs.days)
    cells(zs) = [
        (z, v) for v in 1:nv for z in zs
            if sum(@view counts[zs, v]) > 0
    ]
    function row(label, cs)
        n = length(cs)
        obs(c) = counts[c[1], c[2]]
        pred(c) = view(draws.predictive[c[1]], :, c[2])
        bias = [bias_sample(obs(c), pred(c)) for c in cs]
        cov50 = [_covered(obs(c), pred(c), 0.5) for c in cs]
        cov90 = [_covered(obs(c), pred(c), 0.9) for c in cs]
        return (
            stream = label, n = n,
            bias = round(n == 0 ? NaN : mean(bias); digits = 2),
            coverage_50 = round(n == 0 ? NaN : mean(cov50); digits = 2),
            coverage_90 = round(n == 0 ? NaN : mean(cov90); digits = 2),
        )
    end
    patches = [
        (p, zs) for (p, zs) in enumerate(inputs.patch_ranges)
            if !isempty(zs)
    ]
    rows = [row(inputs.patch_labels[p], cells(zs)) for (p, zs) in patches]
    push!(
        rows,
        row("All zones", reduce(vcat, cells(zs) for (_, zs) in patches))
    )
    return _prettify(DataFrame(rows))
end

## Per-element diagnostic values of a vector deterministic from a
## FlexiChains summary, as a length-`nz` vector; `NaN` where absent.
function _zone_summary_vector(summary, key::Symbol, nz::Integer)
    out = fill(NaN, nz)
    name = string(key)
    for p in FlexiChains.parameters(summary)
        s = string(p)
        v = summary[p]
        if s == name
            vals = collect(skipmissing(vec(collect(v))))
            for (z, x) in enumerate(vals)
                z <= nz && (out[z] = Float64(x))
            end
        elseif startswith(s, name * "[")
            z = parse(Int, s[(length(name) + 2):(end - 1)])
            z <= nz && !ismissing(v) && (out[z] = Float64(v))
        end
    end
    return out
end

"""
$(TYPEDSIGNATURES)

Split R-hat and bulk and tail effective sample sizes of the cut-off
quantities of every zone: the reproduction number `R_T_zone`, the share
`share_T_zone` and the deviation `delta_T_zone`. One row per zone, in
chain order, with the zone's label when `inputs` is given. A quantity
that is constant, or undefined in any draw (a zone below the reporting
floor), shows `NaN`. A level-only zone's `R_T_zone` varies across draws
only at round-off, so its R-hat is meaningless; read the `R_T` columns for
the walking zones (the `walking` column when `inputs` is given).
"""
function zone_diagnostics_table(chn, inputs = nothing)
    nz = if inputs === nothing
        length(first(_draw_vectors(chn, :share_T_zone)))
    else
        length(inputs.zone_keys)
    end
    rhat = FlexiChains.rhat(chn)
    bulk = FlexiChains.ess(chn; kind = :bulk)
    tail = FlexiChains.ess(chn; kind = :tail)
    df = DataFrame(
        zone = inputs === nothing ? string.(1:nz) :
            inputs.zone_labels
    )
    inputs === nothing || (df[!, "walking"] = collect(inputs.walking))
    defined = _zone_rt_defined(chn, nz)
    for key in (:R_T_zone, :share_T_zone, :delta_T_zone)
        stem = replace(string(key), "_zone" => "")
        r = _zone_summary_vector(rhat, key, nz)
        b = _zone_summary_vector(bulk, key, nz)
        t = _zone_summary_vector(tail, key, nz)
        if key === :R_T_zone
            r[.!defined] .= NaN
            b[.!defined] .= NaN
            t[.!defined] .= NaN
        end
        df[!, "rhat_$(stem)"] = r
        df[!, "ess_bulk_$(stem)"] = b
        df[!, "ess_tail_$(stem)"] = t
    end
    return df
end

## Whether each zone's `R_T_zone` is finite in every draw. A zone undefined
## in some draw gets a finite but meaningless R-hat from the rank
## normalisation, so it is reported as having none.
function _zone_rt_defined(chn, nz::Integer)
    _has_key(chn, :R_T_zone) || return trues(nz)
    out = trues(nz)
    for v in _draw_vectors(chn, :R_T_zone), z in 1:min(nz, length(v))
        isfinite(v[z]) || (out[z] = false)
    end
    return out
end

"""
$(TYPEDSIGNATURES)

Sampler-level diagnostics of a zone chain: the number of divergent
transitions, the fraction of iterations that hit the tree-depth cap
`max_depth`, the energy Bayesian fraction of missing information (E-BFMI)
per chain and the adapted step size per chain, plus the worst R-hat and
smallest effective sample sizes over the stored quantities as
[`fit_diagnostics`](@ref) computes them, with `R_T_zone` left out of that
pool (see [`zone_diagnostics_table`](@ref)). With `inputs` the walking
zones' `R_T_zone` R-hat is reported separately as
`max_rhat_R_T_walking`.
"""
function zone_sampler_diagnostics(
        chn, inputs = nothing;
        max_depth::Integer = 8
    )
    depth = _zone_stat(chn, :tree_depth)
    energy = _zone_stat(chn, :hamiltonian_energy)
    steps = _zone_stat(chn, :step_size)
    per_chain(m, f) = m === nothing ? Float64[] :
        [f(view(m, :, c)) for c in 1:size(m, 2)]
    ebfmi(E) = sum(abs2, diff(E)) / max(sum(abs2, E .- mean(E)), floatmin())
    exclude = (_DIAGNOSTIC_EXCLUDE..., "R_T_zone")
    finite(v) = filter(isfinite, v)
    rhats = finite(_scalar_stats(FlexiChains.rhat(chn); exclude))
    bulk = finite(_scalar_stats(FlexiChains.ess(chn; kind = :bulk); exclude))
    tail = finite(_scalar_stats(FlexiChains.ess(chn; kind = :tail); exclude))
    base = (
        max_rhat = isempty(rhats) ? NaN : maximum(rhats),
        min_ess_bulk = isempty(bulk) ? NaN : minimum(bulk),
        min_ess_tail = isempty(tail) ? NaN : minimum(tail),
        n_divergent = _num_divergences(chn),
    )
    walking_rt = if inputs === nothing
        NaN
    else
        nz = length(inputs.zone_keys)
        r = _zone_summary_vector(FlexiChains.rhat(chn), :R_T_zone, nz)
        keep = collect(inputs.walking) .& _zone_rt_defined(chn, nz)
        v = finite(r[keep])
        isempty(v) ? NaN : maximum(v)
    end
    return (;
        base..., max_rhat_R_T_walking = walking_rt,
        depth_cap_fraction = per_chain(depth, x -> mean(x .>= max_depth)),
        ebfmi = per_chain(energy, ebfmi),
        step_size = per_chain(steps, x -> x[end]),
    )
end
