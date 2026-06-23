# One-week-ahead posterior-predictive forecast. Continues the fitted
# renewal trajectory `horizon` days past the cut-off, letting the
# reproduction number keep evolving over the horizon (continuing the
# reconstructed terminal drift of the weekly walk) rather than freezing the
# cut-off growth rate, then scales the fitted expected counts for each
# stream and replicates an integer count per draw so the intervals carry
# both parameter and observation uncertainty.

function _nb_rand(rng, k, μ)
    μs = max(μ, eps(typeof(μ)))
    p = clamp(k / (k + μs), eps(typeof(k)), one(k) - eps(typeof(k)))
    return rand(rng, NegativeBinomial(k, p))
end

## New latent count over the horizon, continuing the cut-off daily value
## `daily_T` at the constant daily growth `r`: the geometric sum
## `daily_T Σ_{d=1}^{h} e^{r d}`, which tends to `daily_T · h` as `r → 0`.
function _geometric_new(daily_T, r, horizon)
    er = exp(r)
    h = float(horizon)
    return abs(er - 1) < 1e-8 ? daily_T * h :
           daily_T * er * (exp(r * h) - 1) / (er - 1)
end

## New latent count over the horizon under a per-day evolving growth rate
## `rs[d]` (one entry per horizon day): the cut-off daily value `daily_T`
## carried forward as `Σ_{d=1}^{h} daily_T · Π_{j≤d} e^{rs[j]}`. With a
## constant rate this reduces to `_geometric_new`.
function _evolving_new(daily_T, rs)
    total = zero(float(daily_T))
    fac = one(float(daily_T))
    @inbounds for r in rs
        fac *= exp(r)
        total += daily_T * fac
    end
    return total
end

## Per-draw growth-rate path over the horizon. The reproduction number keeps
## evolving: the terminal daily drift of the reconstructed log-Rt walk is
## continued forward, the future Rt converted back to a daily growth rate
## through the per-draw generation interval (Euler–Lotka). Returns a
## `Vector{Vector{Float64}}` (one length-`horizon` rate path per draw) and
## the matching terminal forecast reproduction numbers, or `nothing` when
## the chain does not carry the walk and generation-interval parameters
## (single-stream fits), so the caller can fall back to the constant rate.
function _evolving_rates(chn, horizon::Integer)
    has(k) =
        try
            chn[k]
            true
        catch
            false
        end
    (has(Symbol("rt_state.z")) && has(Symbol("gi_state.α")) &&
     has(Symbol("gi_state.θ")) && has(:R_T)) || return nothing

    sigma = _draws(chn, Symbol("rt_state.sigma_rw"))
    R_T = _draws(chn, :R_T)
    α = _draws(chn, Symbol("gi_state.α"))
    θ = _draws(chn, Symbol("gi_state.θ"))
    zmat = chn[Symbol("rt_state.z")]
    zrows = [collect(z) for z in vec(collect(zmat))]
    nd = length(R_T)

    paths = Vector{Vector{Float64}}(undef, nd)
    rt_term = Vector{Float64}(undef, nd)
    @inbounds for i in 1:nd
        ## Terminal weekly drift of the non-centred walk: the last sampled
        ## innovation scaled by its step SD, converted to a per-day log-Rt
        ## slope (a week is seven days). This continues the recent trend of
        ## the walk forward rather than holding Rt flat.
        z = zrows[i]
        last_step = isempty(z) ? 0.0 : sigma[i] * z[end]
        daily_slope = last_step / 7
        ## Generation-interval PMF for this draw (lag-1 indexed), matching
        ## the model's discretisation, so Rt maps to a growth rate the same
        ## way the fit does.
        g = _gi_pmf(α[i], θ[i])
        rs = Vector{Float64}(undef, horizon)
        log_rt = log(max(R_T[i], 1e-6))
        for d in 1:horizon
            log_rt += daily_slope
            rt_d = exp(log_rt)
            rs[d] = euler_lotka_r(rt_d, g)
        end
        paths[i] = rs
        rt_term[i] = exp(log_rt)
    end
    return (; paths, rt_term)
end

## Generation-interval PMF (lag-1 indexed, renormalised) for a Gamma with
## shape `α` and scale `θ`, mirroring `generation_interval_model`.
function _gi_pmf(α, θ; nmax::Integer = 30)
    pmf = discretise_censored(Gamma(α, θ), nmax)
    return pmf[2:end] ./ sum(pmf[2:end])
end

## Per-draw daily latent value at the cut-off from a cumulative-trajectory
## deterministic (the last increment of the cumulative vector). The chain
## stores each trajectory as an iter×chain matrix of per-draw vectors.
## Returns `nothing` when the trajectory is not carried by the chain, so the
## caller can fall back to a scalar daily proxy.
function _daily_at_cutoff(chn, key)
    mat = try
        chn[key]
    catch
        return nothing
    end
    trajs = [collect(v) for v in vec(collect(mat))]
    return Float64[length(t) < 2 ? t[end] : t[end] - t[end - 1]
                   for t in trajs]
end

"""
One-week-ahead (default `horizon = 7` days) posterior-predictive
forecast. For each draw, continue the reproduction number over the horizon
(letting it keep evolving by carrying the walk's terminal drift forward and
mapping back to a per-day growth rate), scale the fitted expected counts by
the resulting horizon growth factor, then replicate the cumulative counts.
Returns a `DataFrame` with one row per draw and columns:

- `:cases_cum`, `:deaths_cum` — replicated cumulative suspected reported
  cases and deaths by the cut-off plus the horizon.
- `:confirmed_cum`, `:confirmed_deaths_cum` — laboratory-confirmed case
  and confirmed-death counterparts, present when `obs_confirmed` and
  `obs_confirmed_deaths` are supplied.
- `:cases_new`, … `:confirmed_deaths_new` — new counts over the coming
  week (`*_cum` minus the corresponding observed count at the cut-off,
  floored at zero).
- `:bed_demand`, `:isolation_level` — the projected isolation/treatment-bed
  DEMAND (need under unconstrained supply) and the supply-limited occupancy
  it produces against the bed capacity, both at the horizon, present when the
  chain carries `expected_bed_demand_T` and `bed_capacity`. The demand grows
  by the horizon factor like the inflow; the occupancy is that demand capped
  at the capacity, `min(demand, C)` (matching the fitted occupancy), so
  `bed_demand − isolation_level` is the projected bed shortfall. Replicated
  with the isolation stream's own dispersion.
- `:recovered_cum`, `:recovered_new` — cumulative recovered-among-confirmed
  by the horizon (and new this week when `obs_recovered` is supplied),
  present when the chain carries `expected_recovered_T`.
- `:infections_new`, `:onsets_new`, `:deaths_latent_new` — new latent
  infections, symptom onsets and deaths projected over the horizon under the
  evolving growth rate (deterministic per-draw quantities, so they carry
  parameter uncertainty only). These are the unobserved counterparts of the
  observed-stream forecasts above.
- `:rt_forecast` — the reproduction number at the end of the horizon, the
  walk's terminal drift continued forward (it evolves rather than freezing
  at the cut-off `R_T`).

Reads `:r`, `:expected_reports_T`, `:expected_deaths_T`,
`:expected_infections_T`, `:R_T`, `:k`, the reproduction-number walk and
generation-interval parameters (to let the reproduction number keep
evolving over the horizon), the cumulative-onset and cumulative-death
trajectories (for the latent onset and death forecasts) and (for the
laboratory streams) `:expected_confirmed_T` and
`:expected_confirmed_deaths_T` from `chn`. When the walk parameters are
not carried (single-stream fits) the cut-off growth rate is held constant
instead. Exports are not forecast: cross-border travel is unlikely to
continue at its baseline rate, so the forward travel rate the export model
relies on no longer holds. The reproduction number is allowed to keep
evolving over the horizon, but no further interventions and no saturation
are imposed.
"""
function forecast_reported(chn;
        horizon::Real = 7,
        obs_cases::Real,
        obs_deaths::Real,
        obs_confirmed::Union{Real, Missing} = missing,
        obs_confirmed_deaths::Union{Real, Missing} = missing,
        obs_recovered::Union{Real, Missing} = missing,
        seed::Integer = 20260520)
    r = _draws(chn, :r)
    cases_T = _draws(chn, :expected_reports_T)
    deaths_T = _draws(chn, :expected_deaths_T)
    infections_T = _draws(chn, :expected_infections_T)
    ## Daily latent onset and death incidence at the cut-off, from the last
    ## increment of the cumulative-trajectory deterministics, grown over the
    ## horizon the same way as the latent infections. When the chain does not
    ## carry the trajectories fall back to the cut-off cumulative reported
    ## cases and deaths as a daily-incidence proxy.
    onsets_T = something(
        _daily_at_cutoff(chn, :cumulative_onsets), cases_T)
    deaths_daily_T = something(
        _daily_at_cutoff(chn, :cumulative_expected_deaths), deaths_T)
    ## Per-draw growth-rate path over the horizon, letting the reproduction
    ## number keep evolving (the walk's terminal drift continued forward).
    ## When the chain does not carry the walk/generation-interval parameters
    ## fall back to the constant cut-off rate and the held terminal R_T.
    evolving = _evolving_rates(chn, Int(horizon))
    rt_forecast = isnothing(evolving) ? _draws(chn, :R_T) : evolving.rt_term
    k = _draws(chn, :k)
    has_conf = obs_confirmed !== missing
    has_conf_deaths = obs_confirmed_deaths !== missing
    conf_T = has_conf ? _draws(chn, :expected_confirmed_T) : nothing
    conf_deaths_T = has_conf_deaths ?
                    _draws(chn, :expected_confirmed_deaths_T) : nothing
    ## Isolation beds and recovered-among-confirmed, forecast when the chain
    ## carries them and using each stream's OWN dispersion. The isolation
    ## stream projects the latent bed DEMAND under unconstrained supply (the
    ## cut-off demand grown by the horizon factor) — the need a week ahead —
    ## and the supply-limited occupancy that demand produces against the bed
    ## capacity, so the forecast quantifies both the need and the shortfall.
    _has(key) =
        try
            chn[key]
            true
        catch
            false
        end
    has_iso = _has(:expected_bed_demand_T) && _has(:bed_capacity) &&
              _has(:isolation_dispersion)
    has_rec = _has(:expected_recovered_T) && _has(:recovered_dispersion)
    demand_T = has_iso ? _draws(chn, :expected_bed_demand_T) : nothing
    cap = has_iso ? _draws(chn, :bed_capacity) : nothing
    k_iso = has_iso ? _draws(chn, :isolation_dispersion) : nothing
    rec_T = has_rec ? _draws(chn, :expected_recovered_T) : nothing
    k_rec = has_rec ? _draws(chn, :recovered_dispersion) : nothing

    rng = MersenneTwister(seed)
    n = length(r)
    cases_cum = Vector{Int}(undef, n)
    deaths_cum = Vector{Int}(undef, n)
    infections_new = Vector{Float64}(undef, n)
    onsets_new = Vector{Float64}(undef, n)
    deaths_latent_new = Vector{Float64}(undef, n)
    confirmed_cum = has_conf ? Vector{Int}(undef, n) : nothing
    confirmed_deaths_cum = has_conf_deaths ? Vector{Int}(undef, n) : nothing
    bed_demand = has_iso ? Vector{Int}(undef, n) : nothing
    isolation_level = has_iso ? Vector{Int}(undef, n) : nothing
    recovered_cum = has_rec ? Vector{Int}(undef, n) : nothing

    @inbounds for i in 1:n
        ## Growth factor over the whole horizon: the product of the daily
        ## evolving factors when the reproduction number keeps evolving,
        ## else the constant cut-off rate compounded over the horizon.
        rs = isnothing(evolving) ? nothing : evolving.paths[i]
        grow = isnothing(rs) ? exp(r[i] * horizon) :
               prod(exp, rs)
        cases_cum[i] = _nb_rand(rng, k[i], cases_T[i] * grow)
        deaths_cum[i] = _nb_rand(rng, k[i], deaths_T[i] * grow)
        ## New latent infections over the horizon, continuing the cut-off
        ## daily infections `I_T` forward. With the reproduction number
        ## evolving the daily rate drifts across the horizon; otherwise it
        ## is the constant-rate geometric sum. Latent, so carries parameter
        ## (not observation) uncertainty across draws.
        if isnothing(rs)
            infections_new[i] = _geometric_new(infections_T[i], r[i], horizon)
            onsets_new[i] = _geometric_new(onsets_T[i], r[i], horizon)
            deaths_latent_new[i] = _geometric_new(
                deaths_daily_T[i], r[i], horizon)
        else
            infections_new[i] = _evolving_new(infections_T[i], rs)
            onsets_new[i] = _evolving_new(onsets_T[i], rs)
            deaths_latent_new[i] = _evolving_new(deaths_daily_T[i], rs)
        end
        has_conf && (confirmed_cum[i] = _nb_rand(rng, k[i], conf_T[i] * grow))
        ## Confirmed deaths grow with the suspected-death signal but are
        ## bounded by it (a thinning), so cap the replicate at the forecast
        ## cumulative suspected deaths.
        if has_conf_deaths
            μ = conf_deaths_T[i] * grow
            confirmed_deaths_cum[i] = min(deaths_cum[i],
                _nb_rand(rng, k[i], μ))
        end
        ## Projected bed demand (need under unconstrained supply) and the
        ## supply-limited occupancy it produces against the bed capacity, plus
        ## cumulative recovered. The demand replicate carries the dispersion;
        ## the occupancy is that same replicate capped at the capacity,
        ## `min(demand, C)`, matching the fitted occupancy and the censored
        ## likelihood, so per draw the occupancy never exceeds the demand or
        ## the capacity.
        if has_iso
            d = _nb_rand(rng, k_iso[i], demand_T[i] * grow)
            bed_demand[i] = d
            isolation_level[i] = min(d, round(Int, cap[i]))
        end
        has_rec && (recovered_cum[i] = _nb_rand(rng, k_rec[i],
            rec_T[i] * grow))
    end

    _new(cum, obs) = max.(cum .- round(Int, obs), 0)
    df = DataFrame(
        cases_cum = cases_cum,
        deaths_cum = deaths_cum,
        cases_new = _new(cases_cum, obs_cases),
        deaths_new = _new(deaths_cum, obs_deaths),
        infections_new = infections_new,
        onsets_new = onsets_new,
        deaths_latent_new = deaths_latent_new,
        rt_forecast = rt_forecast
    )
    if has_conf
        df.confirmed_cum = confirmed_cum
        df.confirmed_new = _new(confirmed_cum, obs_confirmed)
    end
    if has_conf_deaths
        df.confirmed_deaths_cum = confirmed_deaths_cum
        df.confirmed_deaths_new = _new(confirmed_deaths_cum,
            obs_confirmed_deaths)
    end
    if has_iso
        df.bed_demand = bed_demand
        df.isolation_level = isolation_level
    end
    if has_rec
        df.recovered_cum = recovered_cum
        obs_recovered !== missing &&
            (df.recovered_new = _new(recovered_cum, obs_recovered))
    end
    return df
end

"""
Summarise a [`forecast_reported`](@ref) result into a `DataFrame` with
one row per confirmed stream (laboratory-confirmed cases and confirmed
deaths) and quantity (cumulative total by the cut-off plus the horizon, or
new this week), reporting the equal-tailed 30/60/90% credible interval
endpoints (`lower_90 … upper_90`) used by the other summary tables. The
suspected reported-case and suspected-death streams are no longer reported,
so they are not shown as forecast targets.
"""
function forecast_table(fc::DataFrame; digits::Integer = 0)
    _row(label,
        quantity,
        draws) = begin
        s = posterior_summary(draws)
        (stream = label, quantity = quantity,
            lower_90 = round(s.lo90; digits), lower_60 = round(s.lo60; digits),
            lower_30 = round(s.lo30; digits), upper_30 = round(s.hi30; digits),
            upper_60 = round(s.hi60; digits), upper_90 = round(s.hi90; digits))
    end
    streams = Tuple{String, Symbol, Symbol}[]
    :confirmed_cum in propertynames(fc) && push!(streams,
        ("DRC confirmed cases", :confirmed_cum, :confirmed_new))
    :confirmed_deaths_cum in propertynames(fc) && push!(streams,
        ("DRC confirmed deaths", :confirmed_deaths_cum,
            :confirmed_deaths_new))
    rows = NamedTuple[]
    for (label, cum, new) in streams
        push!(rows, _row(label, "cumulative by T+7", fc[!, cum]))
        push!(rows, _row(label, "new this week", fc[!, new]))
    end
    ## Isolation beds: the projected bed DEMAND (need under unconstrained
    ## supply) and the supply-limited occupancy that demand produces, both
    ## levels at the horizon. The gap between them is the bed shortfall.
    :bed_demand in propertynames(fc) && push!(rows,
        _row("DRC isolation beds", "demand at T+7", fc[!, :bed_demand]))
    :isolation_level in propertynames(fc) && push!(rows,
        _row("DRC isolation beds", "occupancy at T+7", fc[!, :isolation_level]))
    if :recovered_cum in propertynames(fc)
        push!(rows,
            _row("DRC recovered", "cumulative by T+7", fc[!, :recovered_cum]))
        :recovered_new in propertynames(fc) && push!(rows,
            _row("DRC recovered", "new this week", fc[!, :recovered_new]))
    end
    return _prettify(DataFrame(rows))
end

"""
Validate a [`forecast_reported`](@ref) projection against the counts that
were later observed. `confirmed` and `confirmed_deaths` are the observed
cumulative DRC laboratory-confirmed cases and confirmed deaths at the
forecast target date. When `isolation` (the observed bed occupancy at the
target date) is supplied and the forecast carries the beds, the projected
supply-limited occupancy is scored against it too — a last-week's-forecast
versus what the beds actually held. Returns a `DataFrame` with one row per
scored stream giving the observed count, the equal-tailed 30/60/90%
predictive intervals (the same endpoints as the other summary tables), and
whether the observed count falls inside the 90% interval. The suspected
reported streams are no longer reported, so they are not scored.

Note that at a one-week-back freeze the bed capacity is weakly informed (the
reported occupancy rate starts only on 9 June), so the projected bed
occupancy rides the capacity random walk back to the freeze date and its
interval is wide.
"""
function forecast_vs_truth(fc::DataFrame;
        confirmed::Real, confirmed_deaths::Real,
        isolation::Union{Real, Missing} = missing,
        digits::Integer = 0)
    _row(label,
        draws,
        obs) = begin
        s = posterior_summary(draws)
        lo = round(s.lo90; digits)
        hi = round(s.hi90; digits)
        (stream = label, observed = round(obs; digits),
            lower_90 = lo, lower_60 = round(s.lo60; digits),
            lower_30 = round(s.lo30; digits), upper_30 = round(s.hi30; digits),
            upper_60 = round(s.hi60; digits), upper_90 = hi,
            within_90 = lo <= obs <= hi ? "yes" : "no")
    end
    rows = NamedTuple[]
    :confirmed_cum in propertynames(fc) &&
        push!(rows, _row("DRC confirmed cases", fc[!, :confirmed_cum],
            confirmed))
    :confirmed_deaths_cum in propertynames(fc) &&
        push!(rows, _row("DRC confirmed deaths", fc[!, :confirmed_deaths_cum],
            confirmed_deaths))
    isolation !== missing && :isolation_level in propertynames(fc) &&
        push!(rows, _row("DRC isolation beds", fc[!, :isolation_level],
            isolation))
    return _prettify(DataFrame(rows))
end

"""
Roll the one-week-ahead forecast across an observed cumulative
trajectory. `targets` is a vector of `(label, horizon_days,
observed_cumulative)` triples: for each, the fitted current growth rate
`r` is projected `horizon_days` past the cut-off and the predicted
cumulative reported cases compared against `observed_cumulative`. Returns
a `DataFrame` with one row per target giving the horizon, the observed
count, the equal-tailed 30/60/90% predictive intervals, and whether the
observed count falls inside the 90% interval. Unlike
[`forecast_vs_truth`](@ref), which scores only the endpoint, this scores
the whole observed trajectory across the horizon. Reads `:r`,
`:expected_reports_T` and `:k` from `chn`.
"""
function forecast_vs_truth_trajectory(
        chn; targets::AbstractVector,
        seed::Integer = 20260520)
    r = _draws(chn, :r)
    cases_T = _draws(chn, :expected_reports_T)
    k = _draws(chn, :k)
    rng = MersenneTwister(seed)
    rows = NamedTuple[]
    for (label, horizon, obs) in targets
        grow = exp.(r .* horizon)
        cases_cum = [_nb_rand(rng, k[i], cases_T[i] * grow[i])
                     for i in eachindex(r)]
        s = posterior_summary(cases_cum)
        lo = round(s.lo90)
        hi = round(s.hi90)
        push!(rows,
            (label = label, horizon_days = horizon,
                observed = round(obs),
                lower_90 = lo, lower_60 = round(s.lo60),
                lower_30 = round(s.lo30), upper_30 = round(s.hi30),
                upper_60 = round(s.hi60), upper_90 = hi,
                within_90 = lo <= obs <= hi ? "yes" : "no"))
    end
    return _prettify(DataFrame(rows))
end
