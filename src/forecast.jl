# One-week-ahead posterior-predictive forecast. Continues the fitted
# renewal trajectory `horizon` days past the cut-off, letting the
# reproduction number keep evolving over the horizon (continuing the
# reconstructed terminal drift of the weekly walk) rather than freezing the
# cut-off growth rate. Cumulative streams add the NEW counts projected over
# the horizon to the cut-off cumulative; rate and prevalence streams (the
# bed demand and daily treatment flows) scale by the horizon growth factor.
# Each is replicated as an integer count per draw so the intervals carry
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

## Daily incidence at the cut-off implied by a cumulative total `C` under
## exponential growth at rate `r` over an outbreak of age `T` days. Inverts
## `C = ∫₀ᵀ i_T e^{r(t-T)} dt`, giving `i_T = C r / (1 - e^{-rT})` (the
## current daily rate, larger than the average `C/T` while growing, smaller
## while declining), with the `r → 0` limit `C/T`. Always non-negative, even
## for `r < 0`. Used as the daily-incidence proxy for observed streams whose
## cumulative trajectory the chain does not carry, so the horizon projection
## can mirror the latent streams instead of scaling the cumulative stock.
function _approx_daily(C, r, T)
    (T <= 0 || !isfinite(T)) && return max(float(C), 0.0) * max(float(r), 0.0)
    rt = r * T
    daily = abs(rt) < 1e-6 ? C / T : C * r / (1 - exp(-rt))
    return max(float(daily), 0.0)
end

## Per-draw growth-rate path over the horizon. The reproduction number keeps
## evolving: the terminal daily drift of the reconstructed log-Rt walk is
## continued forward, the future Rt converted back to a daily growth rate
## through the per-draw generation interval (Euler–Lotka). Returns a
## `Vector{Vector{Float64}}` (one length-`horizon` rate path per draw) and
## the matching terminal forecast reproduction numbers, or `nothing` when
## the chain does not carry the walk and generation-interval parameters, so
## the caller can fall back to the constant rate. Single-stream fits carry
## the walk and generation interval but not the cut-off `R_T` (the joint
## alone exposes it un-prefixed); pass their reconstructed `R_T` draws as
## `R_T` so they take the same evolving path as the joint.
function _evolving_rates(chn, horizon::Integer;
        R_T::Union{Nothing, AbstractVector} = nothing)
    has(k) =
        try
            chn[k]
            true
        catch
            false
        end
    (has(Symbol("rt_state.z")) && has(Symbol("gi_state.α")) &&
     has(Symbol("gi_state.θ")) &&
     (!isnothing(R_T) || has(:R_T))) || return nothing

    sigma = _draws(chn, Symbol("rt_state.sigma_rw"))
    R_T = isnothing(R_T) ? _draws(chn, :R_T) : R_T
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
mapping back to a per-day growth rate), project the new counts each stream
adds over the horizon from its cut-off daily incidence, add them to the
cut-off cumulative and replicate. Returns a `DataFrame` with one row per
draw and columns:

- `:cases_cum`, `:deaths_cum` — replicated cumulative suspected reported
  cases and deaths by the cut-off plus the horizon: the observed cut-off
  cumulative plus the replicated new counts over the horizon, so the
  cumulative never falls below the cut-off even when the reproduction number
  is below one. The cut-off daily incidence projected forward is the last
  increment of the stream's cumulative trajectory when the chain carries it
  (`:cumulative_confirmed` for the laboratory cases), otherwise inferred
  from the cut-off cumulative total under exponential growth.
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
- `:admissions_fc`, `:incare_deaths_fc`, `:ruleouts_fc` — the projected
  one-week-ahead daily isolation/treatment flows (new admissions, in-care
  deaths and rule-outs), present when the chain carries
  `expected_admissions_T`, `expected_incare_deaths_T` and
  `expected_ruleouts_T`. Each cut-off daily flow grows by the horizon factor
  and is replicated through the isolation stream's dispersion, mirroring the
  bed-demand projection.
- `:recovered_cum`, `:recovered_new` — cumulative recovered-among-confirmed
  by the horizon (the cut-off cumulative plus the projected new this week)
  and the new this week, present when the chain carries
  `expected_recovered_T`. The cut-off cumulative is `obs_recovered` when
  supplied, otherwise the fitted `expected_recovered_T`.
- `:infections_new`, `:onsets_new`, `:deaths_latent_new` — new latent
  infections, symptom onsets and deaths projected over the horizon under the
  evolving growth rate (deterministic per-draw quantities, so they carry
  parameter uncertainty only). These are the unobserved counterparts of the
  observed-stream forecasts above.
- `:rt_forecast` — the reproduction number at the end of the horizon, the
  walk's terminal drift continued forward (it evolves rather than freezing
  at the cut-off `R_T`).

Reads `:r`, `:T`, `:expected_reports_T`, `:expected_deaths_T`,
`:expected_infections_T`, `:R_T`, `:k`, the reproduction-number walk and
generation-interval parameters (to let the reproduction number keep
evolving over the horizon), the cumulative-onset and cumulative-death
trajectories (for the latent onset and death forecasts) and (for the
laboratory streams) `:expected_confirmed_T`, `:cumulative_confirmed` and
`:expected_confirmed_deaths_T` from `chn`. The suspected-case and
suspected-death daily incidence at the cut-off is taken from
`:cumulative_reports` / `:cumulative_deaths_total` when carried, otherwise
inferred from the cut-off cumulative total and `:T`. When the walk
parameters are not carried (single-stream fits) the cut-off growth rate is
held constant instead. Exports are not forecast: cross-border travel is unlikely to
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
    ## The three daily treatment flows (admissions, in-care deaths, rule-outs)
    ## share the isolation dispersion, so they need it carried too.
    has_flows = _has(:expected_admissions_T) && _has(:expected_incare_deaths_T) &&
                _has(:expected_ruleouts_T) && _has(:isolation_dispersion)
    demand_T = has_iso ? _draws(chn, :expected_bed_demand_T) : nothing
    cap = has_iso ? _draws(chn, :bed_capacity) : nothing
    k_iso = has_iso ? _draws(chn, :isolation_dispersion) : nothing
    rec_T = has_rec ? _draws(chn, :expected_recovered_T) : nothing
    k_rec = has_rec ? _draws(chn, :recovered_dispersion) : nothing
    admit_T = has_flows ? _draws(chn, :expected_admissions_T) : nothing
    incare_deaths_T = has_flows ? _draws(chn, :expected_incare_deaths_T) :
                      nothing
    ruleout_T = has_flows ? _draws(chn, :expected_ruleouts_T) : nothing
    k_flow = has_flows ? _draws(chn, :isolation_dispersion) : nothing

    ## Daily incidence at the cut-off for each OBSERVED cumulative stream,
    ## projected forward and ADDED to the cut-off cumulative (rather than
    ## scaling the cumulative stock by the horizon growth factor, which would
    ## shrink the cumulative whenever the growth rate is negative). Prefer the
    ## last increment of the stream's cumulative trajectory when the chain
    ## carries it (laboratory-confirmed cases do), otherwise infer the daily
    ## rate from the cumulative total under exponential growth at the cut-off
    ## rate over the outbreak age `T`.
    T_age = _has(:T) ? _draws(chn, :T) : fill(Inf, length(r))
    _approx(C) = _approx_daily.(C, r, T_age)
    cases_daily = something(_daily_at_cutoff(chn, :cumulative_reports),
        _approx(cases_T))
    deaths_daily_obs = something(_daily_at_cutoff(chn, :cumulative_deaths_total),
        _approx(deaths_T))
    conf_daily = has_conf ?
                 something(_daily_at_cutoff(chn, :cumulative_confirmed),
        _approx(conf_T)) : nothing
    conf_deaths_daily = has_conf_deaths ? _approx(conf_deaths_T) : nothing
    rec_daily = has_rec ? _approx(rec_T) : nothing

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
    admissions_fc = has_flows ? Vector{Int}(undef, n) : nothing
    incare_deaths_fc = has_flows ? Vector{Int}(undef, n) : nothing
    ruleouts_fc = has_flows ? Vector{Int}(undef, n) : nothing

    @inbounds for i in 1:n
        ## Growth factor over the whole horizon: the product of the daily
        ## evolving factors when the reproduction number keeps evolving,
        ## else the constant cut-off rate compounded over the horizon.
        rs = isnothing(evolving) ? nothing : evolving.paths[i]
        grow = isnothing(rs) ? exp(r[i] * horizon) :
               prod(exp, rs)
        ## New count over the horizon for a stream continuing forward from its
        ## cut-off daily rate: the constant-rate geometric sum, or the sum
        ## under the evolving per-day rate path. Shared by the latent streams
        ## and the observed cumulative streams.
        _new_h(daily) = isnothing(rs) ? _geometric_new(daily, r[i], horizon) :
                        _evolving_new(daily, rs)
        ## Observed cumulative streams: project NEW counts over the horizon and
        ## add them to the cut-off cumulative, so the projected cumulative
        ## never falls below the cut-off (replicating the cumulative stock
        ## scaled by `grow` would shrink it whenever the growth rate is
        ## negative). The dispersion applies to the projected new counts.
        cases_cum[i] = round(Int, obs_cases) +
                       _nb_rand(rng, k[i], _new_h(cases_daily[i]))
        deaths_cum[i] = round(Int, obs_deaths) +
                        _nb_rand(rng, k[i], _new_h(deaths_daily_obs[i]))
        ## New latent infections over the horizon, continuing the cut-off
        ## daily infections `I_T` forward. With the reproduction number
        ## evolving the daily rate drifts across the horizon; otherwise it
        ## is the constant-rate geometric sum. Latent, so carries parameter
        ## (not observation) uncertainty across draws.
        infections_new[i] = _new_h(infections_T[i])
        onsets_new[i] = _new_h(onsets_T[i])
        deaths_latent_new[i] = _new_h(deaths_daily_T[i])
        has_conf && (confirmed_cum[i] = round(Int, obs_confirmed) +
                            _nb_rand(rng, k[i], _new_h(conf_daily[i])))
        ## Confirmed deaths grow with the suspected-death signal but are
        ## bounded by it (a thinning), so cap the replicate at the forecast
        ## cumulative suspected deaths.
        if has_conf_deaths
            confirmed_deaths_cum[i] = min(deaths_cum[i],
                round(Int, obs_confirmed_deaths) +
                _nb_rand(rng, k[i], _new_h(conf_deaths_daily[i])))
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
        if has_rec
            base_rec = obs_recovered === missing ? round(Int, rec_T[i]) :
                       round(Int, obs_recovered)
            recovered_cum[i] = base_rec +
                               _nb_rand(rng, k_rec[i], _new_h(rec_daily[i]))
        end
        ## One-week-ahead daily treatment flows: each cut-off daily rate grown
        ## by the horizon factor and replicated through the isolation
        ## dispersion, mirroring the bed-demand projection.
        if has_flows
            admissions_fc[i] = _nb_rand(rng, k_flow[i], admit_T[i] * grow)
            incare_deaths_fc[i] = _nb_rand(rng, k_flow[i],
                incare_deaths_T[i] * grow)
            ruleouts_fc[i] = _nb_rand(rng, k_flow[i], ruleout_T[i] * grow)
        end
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
    if has_flows
        df.admissions_fc = admissions_fc
        df.incare_deaths_fc = incare_deaths_fc
        df.ruleouts_fc = ruleouts_fc
    end
    if has_rec
        df.recovered_cum = recovered_cum
        obs_recovered !== missing &&
            (df.recovered_new = _new(recovered_cum, obs_recovered))
    end
    return df
end

"""
    forecast_archive(fcs; made_date, thin = 1) -> DataFrame

Long-format archive of one or more [`forecast_reported`](@ref) results made
from a single cut-off, for the observed streams scored across releases. `fcs`
is an iterable of `(horizon, fc)` pairs and `made_date` is the cut-off `Date`
the forecasts were made from. Returns one row per `(stream, horizon, draw)`
with columns `made_date`, `horizon`, `target_date` (`made_date` plus the
horizon), `stream`, `draw` and `value`.

Only the incident and level quantities are archived: `confirmed cases` and
`confirmed deaths` new over the horizon, `recovered` new over the horizon and
the supply-limited `isolation beds` occupancy. The cumulative totals are not
archived because they are revised across data vintages, so scoring is done on
the incident and level quantities instead. Streams a forecast does not carry
(a single-stream fit, say) are skipped. `thin` keeps every `thin`-th draw so
the archive stays compact when it is saved as a release asset.
"""
function forecast_archive(fcs; made_date::Date, thin::Integer = 1)
    ## Incident (new-over-horizon) and level quantities only, never cumulative.
    streams = (
        (:confirmed_new, "confirmed cases"),
        (:confirmed_deaths_new, "confirmed deaths"),
        (:recovered_new, "recovered"),
        (:isolation_level, "isolation beds"))
    out = DataFrame(made_date = Date[], horizon = Int[], target_date = Date[],
        stream = String[], draw = Int[], value = Float64[])
    for (horizon, fc) in fcs
        h = Int(horizon)
        target = made_date + Day(h)
        for (col, label) in streams
            col in propertynames(fc) || continue
            vals = fc[!, col]
            for (d, i) in enumerate(1:thin:length(vals))
                push!(out,
                    (made_date, h, target, label, d, Float64(vals[i])))
            end
        end
    end
    return out
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
    ## One-week-ahead daily isolation/treatment flows (a single-day rate at
    ## the horizon, not a cumulative total).
    :admissions_fc in propertynames(fc) && push!(rows,
        _row("DRC isolation admissions", "daily at T+7", fc[!, :admissions_fc]))
    :incare_deaths_fc in propertynames(fc) && push!(rows,
        _row("DRC in-care deaths", "daily at T+7", fc[!, :incare_deaths_fc]))
    :ruleouts_fc in propertynames(fc) && push!(rows,
        _row("DRC isolation rule-outs", "daily at T+7", fc[!, :ruleouts_fc]))
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

## Does the chain carry `key`? `FlexiChains` throws rather than returning
## `nothing` for an absent key, so the lookup is probed.
function _has_key(chn, key)
    try
        chn[key]
        return true
    catch
        return false
    end
end

## Draws of the first key in `candidates` the chain carries, or `nothing`
## when it carries none. Single-stream fits bind each stream's submodel
## under a composer-level name (`cases_state.expected_reports`), and the
## joint attaches the same submodels prefixed, so it carries those nested
## names too alongside its own un-prefixed aliases (`expected_reports_T`).
## One ordered candidate list therefore serves both fit kinds.
function _resolve_draws(chn, candidates)
    for key in candidates
        _has_key(chn, key) && return _draws(chn, key)
    end
    return nothing
end

## Per-stream chain-key resolution for [`forecast_stream`](@ref). The joint
## and the single-stream fits do not name these consistently, so the mapping
## is spelled out per stream rather than derived by munging the stream name.
##
## - `expected`: the stream's cut-off expected count. The nested name comes
##   first because both fit kinds carry it; the un-prefixed joint alias is
##   the fallback. The names differ in whether they carry a `_T` suffix
##   (`deaths_state.expected_deaths_T` and `exports_state.expected_exports_T`
##   do, `cases_state.expected_reports` does not).
## - `dispersion`: the joint's un-prefixed per-stream alias comes FIRST here,
##   inverting the order above, because the joint's nested
##   `dispersion_state.k` is the six-element partially-pooled VECTOR, not
##   this stream's scalar. Only the single-stream fits, which sample one
##   scalar `k`, may fall through to it.
## - `trajectory`: the stream's cumulative trajectory, whose last increment
##   is the cut-off daily rate. Only the joint exposes any (as un-prefixed
##   `:=` aliases); an empty list means the daily rate is always inferred
##   from the cut-off cumulative instead.
## - `kind`: `:cumulative` streams report a total accrued to the cut-off, so
##   the forecast is the NEW count over the horizon; `:level` streams report
##   a prevalence, so the forecast is the level at the horizon.
## - `noise`: the observation model the replicate is drawn through.
##   Exports are Poisson ([`exports_model`](@ref)) and carry no dispersion;
##   every other stream is negative-binomial.
const _STREAM_SPEC = Dict{Symbol, NamedTuple}(
    :reported_cases => (
        expected = [Symbol("cases_state.expected_reports"),
            :expected_reports_T],
        dispersion = [:k_cases, Symbol("dispersion_state.k")],
        trajectory = [:cumulative_reports],
        kind = :cumulative, noise = :nb),
    :suspected_deaths => (
        expected = [Symbol("deaths_state.expected_deaths_T"),
            :expected_deaths_T],
        dispersion = [:k_deaths, Symbol("dispersion_state.k")],
        trajectory = [:cumulative_deaths_total],
        kind = :cumulative, noise = :nb),
    :confirmed_cases => (
        expected = [Symbol("confirmed_state.expected_confirmed"),
            :expected_confirmed_T],
        dispersion = [:k_confirmed, Symbol("dispersion_state.k")],
        trajectory = [:cumulative_confirmed],
        kind = :cumulative, noise = :nb),
    :confirmed_deaths => (
        expected = [
            Symbol("confirmed_deaths_state.expected_confirmed_deaths"),
            :expected_confirmed_deaths_T],
        dispersion = [:k_confirmed_deaths, Symbol("dispersion_state.k")],
        trajectory = Symbol[],
        kind = :cumulative, noise = :nb),
    :exports => (
        expected = [Symbol("exports_state.expected_exports_T"),
            :expected_exports_T],
        dispersion = Symbol[],
        trajectory = Symbol[],
        kind = :cumulative, noise = :poisson),
    :isolation_beds => (
        expected = [Symbol("treatment_state.expected_bed_demand"),
            :expected_bed_demand_T],
        dispersion = [:isolation_dispersion,
            Symbol("treatment_state.disp_state.k")],
        trajectory = Symbol[],
        kind = :level, noise = :nb)
)

## Cut-off bed capacity draws. The joint aliases it un-prefixed
## (`bed_capacity := treatment_state.capacity`), but `capacity` is a return
## field of [`treatment_flow_model`](@ref) rather than a `:=`, so a
## standalone treatment fit carries no such key. It does expose
## `bed_utilisation := occ_T / C_T`, so the capacity is recovered exactly by
## dividing the cut-off occupancy through it.
function _bed_capacity(chn)
    _has_key(chn, :bed_capacity) && return _draws(chn, :bed_capacity)
    occ = Symbol("treatment_state.expected_isolation")
    util = Symbol("treatment_state.bed_utilisation")
    (_has_key(chn, occ) && _has_key(chn, util)) || return nothing
    return _draws(chn, occ) ./ _draws(chn, util)
end

## Cut-off reproduction-number draws. The joint exposes `R_T :=
## infection_state.Rt[n]`; single-stream composers do not (the alias lives
## in `bvd_joint`, not in the shared `_latent` submodel), so `R_T` is
## rebuilt from the walk parameters every chain carries. That needs the grid
## length `n` and the intervention `breakpoint`, which are data rather than
## chain contents, so it is only possible when the caller supplies them.
function _cutoff_rt(chn; n, breakpoint, rt_start, rt_walk_start)
    _has_key(chn, :R_T) && return _draws(chn, :R_T)
    (isnothing(n) || isnothing(breakpoint)) && return nothing
    rt = reconstruct_rt(chn; n = n, breakpoint = breakpoint,
        rt_start = rt_start, rt_walk_start = rt_walk_start)
    return Float64[rt[i, n] for i in axes(rt, 1)]
end

## Cut-off growth-rate draws. The joint exposes `r := infection_state.r`;
## single-stream fits do not, so it is derived the way
## [`infection_model`](@ref) itself derives it — forward Euler–Lotka from
## the cut-off reproduction number through the draw's own generation
## interval, which keeps `r < 0` exactly when `R_T < 1`. The chain's
## `growth_state.r` is NOT this quantity: that is the cryptic clock rate the
## joint exposes as `r0`, which is the early growth of the seeded phase
## rather than the growth at the cut-off.
function _cutoff_r(chn, R_T)
    _has_key(chn, :r) && return _draws(chn, :r)
    isnothing(R_T) && return nothing
    α = _draws(chn, Symbol("gi_state.α"))
    θ = _draws(chn, Symbol("gi_state.θ"))
    return Float64[euler_lotka_r(R_T[i], _gi_pmf(α[i], θ[i]))
                   for i in eachindex(R_T)]
end

## Total outbreak-age draws at the cut-off. The joint exposes `T :=
## infection_state.T`; single-stream fits carry only `growth_state.T`, the
## CRYPTIC duration, to which [`infection_model`](@ref) adds the observed
## span `n - rt_start` to reach the total. Falls back to an infinite age
## (the `_approx_daily` limit) when neither is recoverable, matching
## [`forecast_reported`](@ref).
function _outbreak_age(chn, ndraws; n, rt_start)
    _has_key(chn, :T) && return _draws(chn, :T)
    key = Symbol("growth_state.T")
    (_has_key(chn, key) && !isnothing(n)) || return fill(Inf, ndraws)
    return _draws(chn, key) .+ (n - clamp(Int(rt_start), 1, n))
end

"""
    forecast_stream(chn, stream; horizon, obs_value, n, breakpoint) -> Vector

Posterior-predictive forecast of a SINGLE observed stream from either a
single-stream fit or the joint, returning one replicate per draw. Where
[`forecast_reported`](@ref) projects every stream the joint carries at once,
this projects one named stream and resolves its chain keys for both fit
kinds, so each single-stream fit can be scored on forecasting its own
dataset against the joint and against a baseline.

`stream` is one of `:reported_cases`, `:suspected_deaths`,
`:confirmed_cases`, `:confirmed_deaths`, `:isolation_beds` and `:exports`.
The incident streams (everything but `:isolation_beds`) return the NEW count
accrued over the horizon, matching [`forecast_archive`](@ref)'s convention;
`:isolation_beds` returns the supply-limited occupancy LEVEL at the horizon
(the projected demand replicate capped at the bed capacity, `min(demand, C)`,
as the fitted occupancy is).

`obs_value` is the stream's observed count at the cut-off: the cumulative
total for the incident streams (the base the projected new counts accrue on
top of, mirroring [`forecast_reported`](@ref)), or the observed occupancy for
`:isolation_beds`. Since the incident streams return the increment rather
than the cumulative, and the level stream is projected from the fitted stock,
it anchors the reported quantity rather than changing it.

The reproduction number keeps evolving over the horizon exactly as in
[`forecast_reported`](@ref). The joint carries the cut-off `R_T`, `r` and `T`
as un-prefixed deterministics; a single-stream fit carries none of them, so
they are rebuilt from the reproduction-number walk ([`reconstruct_rt`](@ref))
and the generation interval, which every chain does carry. That needs the
grid length `n` and the intervention `breakpoint` — data rather than chain
contents. Pass both for a single-stream chain (with `rt_start` /
`rt_walk_start` if the fit used non-default ones); without them a
single-stream chain cannot be projected and an error says so. They are
unnecessary for the joint, which is read directly.

Exports are forecast here for completeness of the per-stream comparison, but
the caveat [`forecast_reported`](@ref) documents still applies: the export
model relies on a forward travel rate that cross-border movement is unlikely
to keep once the outbreak is known, so the export projection assumes a
baseline travel rate that is unlikely to hold.
"""
function forecast_stream(chn, stream::Symbol;
        horizon::Real = 7,
        obs_value::Real,
        n::Union{Nothing, Integer} = nothing,
        breakpoint::Union{Nothing, Real} = nothing,
        rt_start::Integer = 1,
        rt_walk_start::Integer = rt_start,
        seed::Integer = 20260520)
    spec = get(_STREAM_SPEC, stream, nothing)
    isnothing(spec) && throw(ArgumentError(
        "forecast_stream: unknown stream `:$stream`; expected one of " *
        join(sort!([":$s" for s in keys(_STREAM_SPEC)]), ", ")))

    expected_T = _resolve_draws(chn, spec.expected)
    isnothing(expected_T) && throw(ArgumentError(
        "forecast_stream: chain carries no expected count for `:$stream` " *
        "(tried " * join(["`$k`" for k in spec.expected], ", ") *
        "); the fit does not include this stream."))
    nd = length(expected_T)

    R_T = _cutoff_rt(chn; n = n, breakpoint = breakpoint,
        rt_start = rt_start, rt_walk_start = rt_walk_start)
    r = _cutoff_r(chn, R_T)
    isnothing(r) && throw(ArgumentError(
        "forecast_stream: chain carries neither `r` nor `R_T`, so the " *
        "cut-off growth rate cannot be recovered; pass `n` and " *
        "`breakpoint` to rebuild them from the walk (a single-stream fit " *
        "exposes neither)."))

    ## Per-draw growth-rate path over the horizon, letting the reproduction
    ## number keep evolving. Falls back to the constant cut-off rate when the
    ## chain does not carry the walk and generation interval.
    evolving = _evolving_rates(chn, Int(horizon); R_T = R_T)

    k = nothing
    if spec.noise === :nb
        k = _resolve_draws(chn, spec.dispersion)
        isnothing(k) && throw(ArgumentError(
            "forecast_stream: chain carries no dispersion for `:$stream` " *
            "(tried " * join(["`$key`" for key in spec.dispersion], ", ") *
            ")."))
        ## The joint's nested `dispersion_state.k` is the pooled six-element
        ## vector, so a flattened resolution would silently mismatch the
        ## draws. Every candidate above should be per-draw scalar; fail loudly
        ## if one is not rather than replicating against the wrong `k`.
        length(k) == nd || throw(ArgumentError(
            "forecast_stream: dispersion for `:$stream` resolved to " *
            "$(length(k)) values for $nd draws; expected a scalar per draw."))
    end

    rng = MersenneTwister(seed)
    out = Vector{Int}(undef, nd)

    ## Replicate a projected mean through the stream's observation model.
    _replicate(i, μ) = spec.noise === :poisson ?
                       rand(rng, Poisson(safe_rate(float(μ)))) :
                       _nb_rand(rng, k[i], μ)

    if spec.kind === :level
        ## Prevalence: the cut-off bed demand grown by the horizon factor,
        ## replicated, then capped at the bed capacity (held at its cut-off
        ## value, as `forecast_reported` holds it).
        cap = _bed_capacity(chn)
        isnothing(cap) && throw(ArgumentError(
            "forecast_stream: chain carries no bed capacity (tried " *
            "`bed_capacity` and `treatment_state.expected_isolation` / " *
            "`treatment_state.bed_utilisation`)."))
        @inbounds for i in 1:nd
            rs = isnothing(evolving) ? nothing : evolving.paths[i]
            grow = isnothing(rs) ? exp(r[i] * horizon) : prod(exp, rs)
            demand = _replicate(i, expected_T[i] * grow)
            out[i] = min(demand, round(Int, cap[i]))
        end
        return out
    end

    ## Cumulative: project the NEW count the stream adds over the horizon
    ## from its cut-off daily rate. That rate is the last increment of the
    ## stream's cumulative trajectory when the chain carries it, otherwise
    ## inferred from the cut-off cumulative under exponential growth over the
    ## outbreak age.
    T_age = _outbreak_age(chn, nd; n = n, rt_start = rt_start)
    daily = something(_daily_at_cutoff_any(chn, spec.trajectory),
        _approx_daily.(expected_T, r, T_age))
    @inbounds for i in 1:nd
        rs = isnothing(evolving) ? nothing : evolving.paths[i]
        new_h = isnothing(rs) ? _geometric_new(daily[i], r[i], horizon) :
                _evolving_new(daily[i], rs)
        ## The projected cumulative is `obs_value + replicate` and the new
        ## count over the horizon is that minus `obs_value`, floored at zero,
        ## exactly as `forecast_reported` forms its `*_new` columns.
        cum = round(Int, obs_value) + _replicate(i, new_h)
        out[i] = max(cum - round(Int, obs_value), 0)
    end
    return out
end

## Cut-off daily rate from the first cumulative trajectory in `candidates`
## the chain carries, or `nothing` when it carries none.
function _daily_at_cutoff_any(chn, candidates)
    for key in candidates
        d = _daily_at_cutoff(chn, key)
        isnothing(d) || return d
    end
    return nothing
end
