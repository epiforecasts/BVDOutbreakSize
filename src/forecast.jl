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
- `:confirmed_ward`, `:suspect_ward`, `:confirmed_occupancy`,
  `:suspect_occupancy`, `:confirmed_share` — the confirmed/suspect ward split
  of the beds, present when the chain also carries
  `expected_confirmed_incare_T`. The projected total demand and occupancy are
  partitioned by the per-draw cut-off confirmed share `confirmed_share`
  (`expected_confirmed_incare_T / expected_bed_demand_T`, held flat over the
  horizon), so `confirmed_ward + suspect_ward ≡ bed_demand` and
  `confirmed_occupancy + suspect_occupancy ≡ isolation_level` by construction.
  The split is purely additive: the total `bed_demand` is unchanged.
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
    ## Confirmed/suspect ward split of the isolation beds: partition the total
    ## bed demand and occupancy by the cut-off confirmed SHARE
    ## `expected_confirmed_incare_T / expected_bed_demand_T`. Purely additive —
    ## the total (`bed_demand`) is untouched — so it needs only the confirmed
    ## in-care prevalence the chain already exposes, no model change.
    has_split = has_iso && _has(:expected_confirmed_incare_T)
    has_rec = _has(:expected_recovered_T) && _has(:recovered_dispersion)
    ## The three daily treatment flows (admissions, in-care deaths, rule-outs)
    ## share the isolation dispersion, so they need it carried too.
    has_flows = _has(:expected_admissions_T) && _has(:expected_incare_deaths_T) &&
                _has(:expected_ruleouts_T) && _has(:isolation_dispersion)
    demand_T = has_iso ? _draws(chn, :expected_bed_demand_T) : nothing
    cap = has_iso ? _draws(chn, :bed_capacity) : nothing
    k_iso = has_iso ? _draws(chn, :isolation_dispersion) : nothing
    conf_incare_T = has_split ? _draws(chn, :expected_confirmed_incare_T) :
                    nothing
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
    confirmed_ward = has_split ? Vector{Int}(undef, n) : nothing
    suspect_ward = has_split ? Vector{Int}(undef, n) : nothing
    confirmed_occ = has_split ? Vector{Int}(undef, n) : nothing
    suspect_occ = has_split ? Vector{Int}(undef, n) : nothing
    confirmed_share = has_split ? Vector{Float64}(undef, n) : nothing
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
        ## Confirmed/suspect ward split. Hold the confirmed SHARE flat at its
        ## cut-off value over the horizon: the confirmed and suspect wards move
        ## together with the total, so partitioning by a single per-draw share
        ## preserves the identity `confirmed_ward + suspect_ward ≡ bed_demand`.
        ## (A drift-continued share is a possible future refinement, but flat
        ## vs drift agree to within ~5 beds at T+7, swamped by the total's own
        ## uncertainty, and flat needs no model change since the share is
        ## already exposed on the chain.)
        if has_split
            s_i = clamp(conf_incare_T[i] / max(demand_T[i], eps()), 0.0, 1.0)
            confirmed_share[i] = s_i
            confirmed_ward[i] = round(Int, s_i * bed_demand[i])
            suspect_ward[i] = bed_demand[i] - confirmed_ward[i]
            ## The occupancy split assumes proportional capping: the confirmed
            ## and suspect wards fill the capped beds in the same ratio as
            ## demand. If confirmed cases were prioritised once demand exceeds
            ## capacity this would understate the confirmed occupancy; the
            ## demand split above carries no such assumption.
            confirmed_occ[i] = round(Int, s_i * isolation_level[i])
            suspect_occ[i] = isolation_level[i] - confirmed_occ[i]
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
    if has_split
        df.confirmed_ward = confirmed_ward
        df.suspect_ward = suspect_ward
        df.confirmed_occupancy = confirmed_occ
        df.suspect_occupancy = suspect_occ
        df.confirmed_share = confirmed_share
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
Summarise a [`forecast_reported`](@ref) result into a `DataFrame` with
one row per confirmed stream (laboratory-confirmed cases and confirmed
deaths) and quantity (cumulative total by the cut-off plus the horizon, or
new this week), reporting the equal-tailed 30/60/90% credible interval
endpoints (`lower_90 … upper_90`) used by the other summary tables. The
suspected reported-case and suspected-death streams are no longer reported,
so they are not shown as forecast targets. When the forecast carries the
confirmed/suspect bed split, the beds are reported as three types — the total,
the isolation (suspected) beds and the treatment (confirmed) beds — each with
its demand and occupancy at the horizon. The "Patients en isolement" census is
one occupancy pool split by confirmation status; the report label map is
suspected → isolation beds and confirmed → treatment beds, the total being
their sum.
"""
function forecast_table(fc::DataFrame; digits::Integer = 0)
    _row(label,
        quantity,
        draws) = begin
        s = posterior_summary(draws)
        (stream = label, quantity = quantity,
            lower_90 = round(s.lo90; digits),
            lower_60 = round(s.lo60; digits),
            lower_30 = round(s.lo30; digits),
            upper_30 = round(s.hi30; digits),
            upper_60 = round(s.hi60; digits),
            upper_90 = round(s.hi90; digits))
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
    ## Beds: the projected DEMAND (the need under unconstrained supply) and the
    ## supply-limited OCCUPANCY it produces, both as levels at the horizon. The
    ## gap between demand and occupancy is the bed shortfall. The census is one
    ## occupancy pool split by confirmation status; when the split is present
    ## the total is reported alongside the isolation (suspected) and treatment
    ## (confirmed) beds. Column names stay by confirmation status (unambiguous
    ## facts); only the display labels map suspected → isolation, confirmed →
    ## treatment.
    bed_streams = Tuple{String, Symbol, Symbol}[]
    (:bed_demand in propertynames(fc) &&
     :isolation_level in propertynames(fc)) && push!(bed_streams,
        ("DRC beds (total)", :bed_demand, :isolation_level))
    (:suspect_ward in propertynames(fc) &&
     :suspect_occupancy in propertynames(fc)) && push!(bed_streams,
        ("DRC isolation beds (suspected)", :suspect_ward, :suspect_occupancy))
    (:confirmed_ward in propertynames(fc) &&
     :confirmed_occupancy in propertynames(fc)) && push!(bed_streams,
        ("DRC treatment beds (confirmed)", :confirmed_ward,
            :confirmed_occupancy))
    for (label, dem, occ) in bed_streams
        push!(rows, _row(label, "demand at T+7", fc[!, dem]))
        push!(rows, _row(label, "occupancy at T+7", fc[!, occ]))
    end
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
        push!(rows,
            _row("DRC confirmed deaths", fc[!, :confirmed_deaths_cum],
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
