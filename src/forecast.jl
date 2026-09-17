# One-week-ahead posterior-predictive forecast. Continues the fitted renewal
# trajectory `horizon` days past the cut-off, letting the reproduction number
# keep evolving over the horizon. Cumulative streams add the projected new
# counts to the cut-off cumulative. Rate and prevalence streams scale by the
# horizon growth factor. Each is replicated as an integer count per draw, so
# the intervals carry both parameter and observation uncertainty.

## Floor on the reproduction number passed to `euler_lotka_r`. Its Newton
## solve overflows to a non-finite growth rate below R ≈ 1e-3, which a draw
## forecasting a steep decline can reach over the horizon. Below 1e-2 the
## regime has already collapsed and the forecast need not resolve it.
const _RT_EULER_FLOOR = 1e-2

function _nb_rand(rng, k, μ)
    μs = max(μ, eps(typeof(μ)))
    p = clamp(k / (k + μs), eps(typeof(k)), one(k) - eps(typeof(k)))
    return rand(rng, NegativeBinomial(k, p))
end

## Per-day projected latent means over the horizon: the cut-off daily value
## `daily_T` carried forward as `daily_T · Π_{j≤d} e^{rs[j]}`, one entry per
## horizon day, under the per-day rate path `rs`. A chain that does not carry
## the walk holds the cut-off rate `r` constant, which reduces the product to
## the geometric `daily_T · e^{r d}`.
function _daily_means(daily_T, rs, r, horizon::Integer)
    means = Vector{Float64}(undef, horizon)
    fac = 1.0
    @inbounds for d in 1:horizon
        fac *= exp(isnothing(rs) ? r : rs[d])
        means[d] = daily_T * fac
    end
    return means
end

## Replicated new count over the horizon: one draw per projected day, summed.
## The dispersion is fitted against single-day counts and between-vintage
## increments spanning a day or three, so it belongs at that resolution
## rather than on the horizon total.
function _nb_new(rng, k, means)
    total = 0
    @inbounds for m in means
        total += _nb_rand(rng, k, m)
    end
    return total
end

## Daily incidence at the cut-off implied by a cumulative total `C` under
## exponential growth at rate `r` over an outbreak of age `T` days. Inverts
## `C = ∫₀ᵀ i_T e^{r(t-T)} dt`, giving `i_T = C r / (1 - e^{-rT})`, with the
## `r → 0` limit `C/T`. Always non-negative, even for `r < 0`. Used as the
## daily-incidence proxy for observed streams whose cumulative trajectory the
## chain does not carry.
function _approx_daily(C, r, T)
    (T <= 0 || !isfinite(T)) && return max(float(C), 0.0) * max(float(r), 0.0)
    rt = r * T
    daily = abs(rt) < 1e-6 ? C / T : C * r / (1 - exp(-rt))
    return max(float(daily), 0.0)
end

## Per-draw growth-rate path over the horizon. The fitted weekly walk is
## continued past the cut-off with fresh innovations at its own step scale,
## and the future Rt converted back to a daily growth rate through the
## per-draw generation interval (Euler–Lotka). Returns one length-`horizon`
## rate path per draw and the matching terminal reproduction numbers, or
## `nothing` when the chain carries no walk and generation-interval
## parameters, so the caller can fall back to the constant rate.
## Single-stream fits carry the walk but not the cut-off `R_T`, which only
## the joint exposes un-prefixed. Pass their reconstructed `R_T` draws as
## `R_T` so they take the same evolving path as the joint.
function _evolving_rates(chn, horizon::Integer;
        R_T::Union{Nothing, AbstractVector} = nothing,
        rng::AbstractRNG = MersenneTwister(20260520),
        week::Integer = 7)
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
    nd = length(R_T)
    nknots = cld(horizon, week)

    paths = Vector{Vector{Float64}}(undef, nd)
    rt_term = Vector{Float64}(undef, nd)
    @inbounds for i in 1:nd
        ## The fitted walk sits on weekly knots with innovations of SD
        ## `sigma_rw`, and daily log-Rt is the linear interpolation between
        ## them ([`rt_walk_model`](@ref)). The horizon draws fresh weekly
        ## innovations from the cut-off level and interpolates across the days
        ## between, so the log-Rt spread grows as the square root of the
        ## horizon, as the fitted walk's does.
        innov = sigma[i] .* randn(rng, nknots)
        cum_innov = cumsum(innov)
        ## Generation-interval PMF for this draw (lag-1 indexed), matching the
        ## model's discretisation, so Rt maps to a growth rate as the fit does.
        g = _gi_pmf(α[i], θ[i])
        rs = Vector{Float64}(undef, horizon)
        log_R = log(max(R_T[i], 1e-6))
        log_rt = log_R
        for d in 1:horizon
            weeks = d / week
            j = floor(Int, weeks)
            whole = j == 0 ? 0.0 : cum_innov[j]
            part = j < nknots ? (weeks - j) * innov[j + 1] : 0.0
            log_rt = log_R + whole + part
            rt_d = exp(log_rt)
            rs[d] = euler_lotka_r(max(rt_d, _RT_EULER_FLOOR), g)
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
One-week-ahead (default `horizon = 7` days) posterior-predictive forecast.
For each draw, continue the reproduction number over the horizon, project the
new counts each stream adds from its cut-off daily incidence, add them to the
cut-off cumulative and replicate. Each count stream is replicated day by day
through its own fitted dispersion and the daily replicates summed, matching
the resolution the dispersion was fitted at. Returns a `DataFrame` with one
row per draw and columns:

- `:cases_cum`, `:deaths_cum`: replicated cumulative suspected reported cases
  and deaths by the cut-off plus the horizon, the observed cut-off cumulative
  plus the replicated new counts. The cumulative never falls below the cut-off
  even when the reproduction number is below one. The cut-off daily incidence
  is the last increment of the stream's cumulative trajectory when the chain
  carries it, otherwise inferred from the cut-off cumulative total under
  exponential growth.
- `:confirmed_cum`, `:confirmed_deaths_cum`: laboratory-confirmed case and
  confirmed-death counterparts, present when `obs_confirmed` and
  `obs_confirmed_deaths` are supplied. `:confirmed_deaths_cum` is not bounded
  by `:deaths_cum`. The subset relationship is imposed per day on the latent
  pool ([`confirmed_deaths_model`](@ref)), and either reported cumulative can
  stall while the other accrues.
- `:cases_new`, … `:confirmed_deaths_new`: new counts over the coming week
  (`*_cum` minus the corresponding observed count at the cut-off, floored at
  zero).
- `:bed_demand`, `:isolation_level`, `:bed_shortfall`: the projected
  isolation/treatment-bed demand (need under unconstrained supply), the
  occupancy that demand produces as it is reported, and the unmet demand, all
  at the horizon, present when the chain carries `expected_bed_demand_T` and
  `bed_capacity`. The demand grows by the horizon factor like the inflow. The
  occupancy is that demand shifted by the fitted reclassification offset Δ
  standing at the cut-off and capped at the capacity, `min(demand + Δ, C)`,
  the quantity the occupancy likelihood scores against the reported series
  (see [`cumulative_occupancy_offset`](@ref)). The shortfall is
  `max(demand − C, 0)`, taken against the capacity rather than as
  `bed_demand − isolation_level`, since Δ is a change of reporting basis and
  not a bed. Replicated with the isolation stream's own dispersion.
- `:admissions_fc`, `:incare_deaths_fc`, `:ruleouts_fc`: the projected
  one-week-ahead daily isolation/treatment flows (new admissions, in-care
  deaths and rule-outs), present when the chain carries
  `expected_admissions_T`, `expected_incare_deaths_T` and
  `expected_ruleouts_T`. Each cut-off daily flow grows by the horizon factor
  and is replicated through the isolation stream's dispersion.
- `:recovered_cum`, `:recovered_new`: cumulative recovered-among-confirmed by
  the horizon and the new this week, present when the chain carries
  `expected_recovered_T`. The cut-off cumulative is `obs_recovered` when
  supplied, otherwise the fitted `expected_recovered_T`.
- `:infections_new`, `:onsets_new`, `:deaths_latent_new`: new latent
  infections, symptom onsets and deaths projected over the horizon
  (deterministic per-draw quantities, so they carry parameter uncertainty
  only).
- `:rt_forecast`: the reproduction number at the end of the horizon, the
  fitted weekly walk continued past the cut-off.
- the [`forecast_onsets`](@ref) columns (`onsets_to_date`,
  `onsets_unreported`, `onset_reports_backfill`, `onset_reports_new`, …),
  present when `onset_grid_start` and `onset_grid_end` are supplied and the
  chain carries the fitted reporting hazard. Of them only `onset_reports_new`
  is archived and scored, since it is the one the triangle gives an
  observation for. The model cut-off grid day goes in as `grid_n` rather than
  `n`, which is already the draw count here.

Reads `:r`, `:T`, `:expected_reports_T`, `:expected_deaths_T`,
`:expected_infections_T`, `:R_T`, the per-stream dispersions (`:k_cases`,
`:k_deaths`, `:k_confirmed`, `:k_confirmed_deaths`, falling back to the
population-mean `:k`), the reproduction-number walk and generation-interval
parameters, the cumulative-onset and cumulative-death trajectories and (for
the laboratory streams) `:expected_confirmed_T` and
`:expected_confirmed_deaths_T` from `chn`. Each observed stream's daily
incidence at the cut-off is the last increment of its own cumulative
trajectory (`:cumulative_reports`, `:cumulative_deaths_total`,
`:cumulative_confirmed`, `:cumulative_confirmed_deaths`,
`:cumulative_recovered`), otherwise inferred from the cut-off cumulative
total and `:T`. When the walk parameters are not carried (single-stream fits)
the cut-off growth rate is held constant instead.

Exports are not forecast, since cross-border travel is unlikely to continue
at the baseline rate the export model relies on. The reproduction number is
allowed to keep evolving over the horizon, but no further interventions and
no saturation are imposed.
"""
function forecast_reported(chn;
        horizon::Real = 7,
        obs_cases::Real,
        obs_deaths::Real,
        obs_confirmed::Union{Real, Missing} = missing,
        obs_confirmed_deaths::Union{Real, Missing} = missing,
        obs_recovered::Union{Real, Missing} = missing,
        onset_grid_start::Union{Nothing, Integer} = nothing,
        onset_grid_end::Union{Nothing, Integer} = nothing,
        grid_n::Union{Nothing, Integer} = nothing,
        seed::Integer = 20260520)
    r = _draws(chn, :r)
    cases_T = _draws(chn, :expected_reports_T)
    deaths_T = _draws(chn, :expected_deaths_T)
    infections_T = _draws(chn, :expected_infections_T)
    ## Daily latent onset and death incidence at the cut-off, from the last
    ## increment of the cumulative-trajectory deterministics. A chain that
    ## does not carry the trajectories falls back to the cut-off cumulative
    ## reported cases and deaths as a daily-incidence proxy.
    onsets_T = something(
        _daily_at_cutoff(chn, :cumulative_onsets), cases_T)
    deaths_daily_T = something(
        _daily_at_cutoff(chn, :cumulative_expected_deaths), deaths_T)
    ## The rate paths and the observation replicates below draw from one
    ## seeded stream, so the whole forecast is reproducible.
    nh = Int(horizon)
    rng = MersenneTwister(seed)
    evolving = _evolving_rates(chn, nh; rng = rng)
    rt_forecast = isnothing(evolving) ? _draws(chn, :R_T) : evolving.rt_term
    ## Per-stream dispersions, resolved from the same keys `forecast_stream`
    ## uses so the two forecasters agree. `k` is the population-mean
    ## dispersion and the fallback for a chain carrying no per-stream alias.
    k = _draws(chn, :k)
    _stream_k(stream) = begin
        d = _resolve_draws(chn, _STREAM_SPEC[stream].dispersion)
        (isnothing(d) || length(d) != length(k)) ? k : d
    end
    k_cases = _stream_k(:reported_cases)
    k_deaths = _stream_k(:suspected_deaths)
    k_conf = _stream_k(:confirmed_cases)
    k_conf_deaths = _stream_k(:confirmed_deaths)
    has_conf = obs_confirmed !== missing
    has_conf_deaths = obs_confirmed_deaths !== missing
    conf_T = has_conf ? _draws(chn, :expected_confirmed_T) : nothing
    conf_deaths_T = has_conf_deaths ?
                    _draws(chn, :expected_confirmed_deaths_T) : nothing
    ## Isolation beds and recovered-among-confirmed, forecast when the chain
    ## carries them and using each stream's own dispersion.
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
    occ_offset = has_iso ? _occupancy_offset(chn, length(r)) : nothing
    k_iso = has_iso ? _draws(chn, :isolation_dispersion) : nothing
    rec_T = has_rec ? _draws(chn, :expected_recovered_T) : nothing
    k_rec = has_rec ? _draws(chn, :recovered_dispersion) : nothing
    admit_T = has_flows ? _draws(chn, :expected_admissions_T) : nothing
    incare_deaths_T = has_flows ? _draws(chn, :expected_incare_deaths_T) :
                      nothing
    ruleout_T = has_flows ? _draws(chn, :expected_ruleouts_T) : nothing
    k_flow = has_flows ? _draws(chn, :isolation_dispersion) : nothing

    ## Daily incidence at the cut-off for each observed cumulative stream.
    ## Prefer the last increment of the stream's cumulative trajectory when
    ## the chain carries it, otherwise infer the daily rate from the
    ## cumulative total under exponential growth at the cut-off rate over the
    ## outbreak age `T`.
    T_age = _has(:T) ? _draws(chn, :T) : fill(Inf, length(r))
    _approx(C) = _approx_daily.(C, r, T_age)
    cases_daily = something(_daily_at_cutoff(chn, :cumulative_reports),
        _approx(cases_T))
    deaths_daily_obs = something(_daily_at_cutoff(chn, :cumulative_deaths_total),
        _approx(deaths_T))
    conf_daily = has_conf ?
                 something(_daily_at_cutoff(chn, :cumulative_confirmed),
        _approx(conf_T)) : nothing
    conf_deaths_daily = has_conf_deaths ?
                        something(
        _daily_at_cutoff(chn, :cumulative_confirmed_deaths),
        _approx(conf_deaths_T)) : nothing
    rec_daily = has_rec ?
                something(_daily_at_cutoff(chn, :cumulative_recovered),
        _approx(rec_T)) : nothing

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
    bed_shortfall = has_iso ? Vector{Int}(undef, n) : nothing
    recovered_cum = has_rec ? Vector{Int}(undef, n) : nothing
    admissions_fc = has_flows ? Vector{Int}(undef, n) : nothing
    incare_deaths_fc = has_flows ? Vector{Int}(undef, n) : nothing
    ruleouts_fc = has_flows ? Vector{Int}(undef, n) : nothing

    @inbounds for i in 1:n
        ## Growth factor over the whole horizon: the product of the daily
        ## evolving factors, else the constant cut-off rate compounded.
        rs = isnothing(evolving) ? nothing : evolving.paths[i]
        grow = isnothing(rs) ? exp(r[i] * horizon) :
               prod(exp, rs)
        _means(daily) = _daily_means(daily, rs, r[i], nh)
        _new_h(daily) = sum(_means(daily))
        ## Observed cumulative streams: project new counts over the horizon
        ## and add them to the cut-off cumulative, so the projected cumulative
        ## never falls below the cut-off. Scaling the cumulative stock by
        ## `grow` would shrink it whenever the growth rate is negative.
        cases_cum[i] = round(Int, obs_cases) +
                       _nb_new(rng, k_cases[i], _means(cases_daily[i]))
        deaths_cum[i] = round(Int, obs_deaths) +
                        _nb_new(rng, k_deaths[i], _means(deaths_daily_obs[i]))
        ## Latent streams, so these carry parameter uncertainty across draws
        ## and no observation uncertainty.
        infections_new[i] = _new_h(infections_T[i])
        onsets_new[i] = _new_h(onsets_T[i])
        deaths_latent_new[i] = _new_h(deaths_daily_T[i])
        has_conf && (confirmed_cum[i] = round(Int, obs_confirmed) +
                            _nb_new(rng, k_conf[i], _means(conf_daily[i])))
        ## Confirmed deaths carry no cap against the suspected-death
        ## cumulative. The subset relationship is imposed inside the model,
        ## per day and on the latent pool, by capping the tested death volume
        ## at `susp_death[t]` before positivity scales it
        ## ([`confirmed_deaths_model`](@ref)), so a forecast projected from
        ## that trajectory carries the constraint with it.
        if has_conf_deaths
            confirmed_deaths_cum[i] = round(Int, obs_confirmed_deaths) +
                                      _nb_new(rng, k_conf_deaths[i],
                _means(conf_deaths_daily[i]))
        end
        ## The demand replicate carries the dispersion. The occupancy is that
        ## same replicate shifted by the fitted reclassification offset and
        ## capped at the capacity, `min(demand + Δ, C)`, which the censored
        ## occupancy likelihood scores against the reported series. The
        ## shortfall is taken against the capacity, so a change of reporting
        ## basis is not counted as unmet demand.
        if has_iso
            d = _nb_rand(rng, k_iso[i], demand_T[i] * grow)
            bed_demand[i] = d
            isolation_level[i] = _reported_occupancy(d, occ_offset[i], cap[i])
            bed_shortfall[i] = max(d - round(Int, cap[i]), 0)
        end
        if has_rec
            base_rec = obs_recovered === missing ? round(Int, rec_T[i]) :
                       round(Int, obs_recovered)
            recovered_cum[i] = base_rec +
                               _nb_new(rng, k_rec[i], _means(rec_daily[i]))
        end
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
        df.bed_shortfall = bed_shortfall
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
    ## Symptom-onset nowcast and forecast, attached when the caller supplies
    ## the fitted triangle's grid (data, not chain contents) and the chain
    ## carries the reporting hazard.
    if !isnothing(onset_grid_start) && !isnothing(onset_grid_end) &&
       !isnothing(grid_n) &&
       _has_key(chn, Symbol("onset_report_state.σ_mult"))
        onset_fc = forecast_onsets(chn; grid_start = onset_grid_start,
            grid_end = onset_grid_end, n = grid_n, horizon = horizon,
            seed = seed)
        ## `onsets_new` is already a column here, built from the same cut-off
        ## daily onset rate and the same evolving growth path. The existing
        ## one is kept rather than overwritten, so a divergence between the
        ## two would show rather than being hidden.
        for col in propertynames(onset_fc)
            col in propertynames(df) && continue
            df[!, col] = onset_fc[!, col]
        end
    end
    return df
end

## Per-draw daily symptom-onset series over `1:n`, the increments of the
## chain's `cumulative_onsets` trajectory (the model stores the running sum,
## not the daily series). Returns `nothing` when the chain does not carry
## the trajectory, so the caller can say the fit has no onset stream.
function _onset_daily_series(chn)
    mat = try
        chn[:cumulative_onsets]
    catch
        return nothing
    end
    return [(v = collect(t); vcat(v[1], diff(v)))
            for t in vec(collect(mat))]
end

"""
    forecast_onsets(chn; grid_start, grid_end, n, horizon = 7, ...)
        -> DataFrame

Nowcast and forecast of symptom onsets and of the reporting triangle that
observes them, one row per posterior draw. It separates the two things a
reporting triangle can tell apart: cases that have already had their onset
but have not yet been reported, and cases whose onset has not yet happened.

The reported total as of grid day `a` is
[`onset_report_expected_total`](@ref)'s `Σ_u onsets[u] · F(u, a - u)`, with
`F(u, δ) = α(u) · G(u, δ)` ([`onset_report_F`](@ref)) the fitted delay
hazard's cumulative reported proportion and `α(u)` the fitted ascertainment
level. Splitting that sum at the cut-off `n` and differencing it between `n`
and `n + h` gives the columns below. Onsets after the cut-off are projected
by carrying the cut-off daily onset rate forward under the same evolving
growth-rate path every other forecast uses (see
[`forecast_reported`](@ref)). The calendar-time effect `γ` and the
ascertainment level `α` are both held flat at their last fitted value over
the horizon ([`onset_report_F`](@ref) already does this, see
[`onset_report_G`](@ref)), the same assumption a nowcast makes.

Columns, all per draw:

  - `onsets_to_date`: symptom onsets that have already happened by the
    cut-off, `Σ_{u ≤ n} onsets[u]`. The nowcast, free of reporting delay
    and of the ascertainment level `α`.
  - `onset_reports_to_date`: of those, the number the triangle should
    already have printed, `Σ_{u ≤ n} onsets[u] · F(u, n - u)`. The
    quantity the digitised total is an observation of.
  - `onsets_unreported`: their difference, onsets that have happened but
    are not yet reported. `F` does not tend to one, so this carries both
    the reporting backlog and the cases ascertainment will never pick up.
  - `onsets_new`: onsets over the horizon, `Σ_{n < u ≤ n+h} onsets[u]`,
    what has not yet happened.
  - `onset_reports_backfill`: reports arriving over the horizon for
    onsets on or before the cut-off, the part of `onsets_unreported` the
    triangle should fill in this week.
  - `onset_reports_future`: reports arriving over the horizon for onsets
    after the cut-off.
  - `onset_reports_new`: their sum, the replicated new reported count
    over the horizon, and the scored forecast. It is an increment rather
    than a level because the digitised total is revised between vintages
    (see [`forecast_archive`](@ref)). Floored at zero to match the
    `max(new, 0)` convention the observed truth is built with.
  - `onset_reports_cum`: `obs_value + onset_reports_new` when
    `obs_value` (the triangle's own cumulative total at the cut-off) is
    supplied, so the projection can be plotted on the observed scale.

`onset_reports_new` is replicated through the stream's own observation model
rather than a negative binomial. The increment is scored under
[`onset_report_scale`](@ref)'s three-term scale (counting variation,
pixel-reading noise, per-scan level error), inflated by the fitted `σ_mult`,
and perturbed by a Student-t with the same `ν` the likelihood uses. The scale
is applied once to the whole projected total, so its scan-level term takes
the fitted shared per-scan coefficient `σ_scan` rather than the whole
measured per-bar error. At a total of a couple of thousand cases the per-scan
level term dominates by an order of magnitude, so the interval on a weekly
increment is mostly digitisation error rather than epidemic uncertainty. The
forecast is therefore more useful as a check that the fitted delay and
ascertainment reproduce the next vintage than as a case-count prediction.

`grid_start` and `grid_end` are the fitted triangle's own onset/report-day
grid (see [`reconstruct_onset_hazard`](@ref)). `n` is the model cut-off
grid day and `breakpoint` the intervention breakpoint, needed only when
the chain is a single-stream fit that does not carry `R_T`/`r`
(see [`forecast_stream`](@ref)).
"""
function forecast_onsets(chn;
        grid_start::Integer,
        grid_end::Integer,
        n::Integer,
        horizon::Real = 7,
        obs_value::Union{Real, Missing} = missing,
        week::Integer = 7,
        ν::Real = 4.0,
        pixel_sd::Real = 2.1,
        scan_frac::Real = 0.04,
        breakpoint::Union{Nothing, Real} = nothing,
        rt_start::Integer = 1,
        rt_walk_start::Integer = rt_start,
        seed::Integer = 20260520)
    daily = _onset_daily_series(chn)
    isnothing(daily) && throw(ArgumentError(
        "forecast_onsets: chain carries no `cumulative_onsets` trajectory, " *
        "so the latent onset series cannot be recovered. This fit has no " *
        "onset stream."))
    hazard = reconstruct_onset_hazard(chn; grid_start, grid_end, week)
    σ_mult = _draws(chn, Symbol("onset_report_state.σ_mult"))
    ## A chain that does not sample the onset-report overdispersion falls
    ## back to no quadratic term, which is that fit's own likelihood.
    k_onset = _has_key(chn, Symbol("onset_report_state.k_onset")) ?
              _draws(chn, Symbol("onset_report_state.k_onset")) : nothing
    ## Per-scan level error for the projected snapshot. The likelihood carries
    ## this on the modelled level rather than in the per-cell scale (see
    ## `onset_reporting_model`), but a projected total is scored against a
    ## scan that has not happened yet, so its own level error has to enter the
    ## scale here. A chain that does not sample it falls back to `scan_frac`.
    σ_scan = _has_key(chn, Symbol("onset_report_state.σ_scan")) ?
             _draws(chn, Symbol("onset_report_state.σ_scan")) : nothing

    R_T = _cutoff_rt(chn; n = n, breakpoint = breakpoint,
        rt_start = rt_start, rt_walk_start = rt_walk_start)
    r = _cutoff_r(chn, R_T)
    isnothing(r) && throw(ArgumentError(
        "forecast_onsets: chain carries neither `r` nor `R_T`, so the " *
        "cut-off growth rate cannot be recovered. Pass `n` and " *
        "`breakpoint` to rebuild them from the walk."))
    h = Int(horizon)
    rng = MersenneTwister(seed)
    evolving = _evolving_rates(chn, h; R_T = R_T, rng = rng)

    nd = length(daily)
    onsets_to_date = Vector{Float64}(undef, nd)
    reports_to_date = Vector{Float64}(undef, nd)
    onsets_new = Vector{Float64}(undef, nd)
    backfill = Vector{Float64}(undef, nd)
    future = Vector{Float64}(undef, nd)
    reports_new = Vector{Int}(undef, nd)

    @inbounds for i in 1:nd
        o = daily[i]
        lh0 = hazard.logit_h0[i]
        γ = hazard.γ[i]
        α = hazard.alpha[i]
        na = length(α)
        ## Future daily onsets, built exactly as `forecast_reported` builds
        ## its `onsets_new`.
        rs = isnothing(evolving) ? nothing : evolving.paths[i]
        daily_T = o[min(n, length(o))]
        fut = _daily_means(daily_T, rs, r[i], h)

        ## Reported totals at the cut-off and at the horizon, split by whether
        ## the onset date is on or after the cut-off. `F` is evaluated with
        ## the calendar walk and the ascertainment level both held flat
        ## outside their fitted support.
        past_now = 0.0
        past_then = 0.0
        ge = min(n, length(o))
        for u in 1:ge
            αu = α[clamp(u - grid_start + 1, 1, na)]
            f_now = onset_report_F(n - u, lh0, γ, u, grid_start, αu)
            f_then = onset_report_F(n + h - u, lh0, γ, u, grid_start, αu)
            past_now += o[u] * f_now
            past_then += o[u] * f_then
        end
        fut_then = 0.0
        for d in 1:h
            u = n + d
            αu = α[clamp(u - grid_start + 1, 1, na)]
            fut_then += fut[d] *
                        onset_report_F(n + h - u, lh0, γ, u, grid_start, αu)
        end

        onsets_to_date[i] = sum(@view o[1:ge])
        reports_to_date[i] = past_now
        onsets_new[i] = sum(fut)
        backfill[i] = past_then - past_now
        future[i] = fut_then

        ## Replicate the projected increment through the same three-term
        ## observation scale a scored correction cell carries: two reads
        ## (`reads = 2`, this is a difference of two vintages), the levels
        ## being the reported totals at the two ends of the horizon.
        μ = backfill[i] + future[i]
        base = onset_report_scale(μ, past_then, past_now, 2;
            pixel_sd, scan_sd = isnothing(σ_scan) ? scan_frac : σ_scan[i])
        σ = σ_mult[i] * (isnothing(k_onset) ? base :
             sqrt(base^2 + μ^2 / max(k_onset[i], eps(Float64))))
        reports_new[i] = max(round(Int, μ + σ * rand(rng, TDist(ν))), 0)
    end

    df = DataFrame(
        onsets_to_date = onsets_to_date,
        onset_reports_to_date = reports_to_date,
        onsets_unreported = onsets_to_date .- reports_to_date,
        onsets_new = onsets_new,
        onset_reports_backfill = backfill,
        onset_reports_future = future,
        onset_reports_new = reports_new)
    obs_value === missing ||
        (df.onset_reports_cum = round(Int, obs_value) .+ reports_new)
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
`confirmed deaths` new over the horizon, `recovered` new over the horizon,
`onset reports` new over the horizon, and the supply-limited `isolation beds`
occupancy. The cumulative totals are revised across data vintages, so they
are not archived. The reporting triangle's own total is revised by more than
the rest, since the ≈4% per-scan level error alone moves it by tens of cases
and the printed total falls between consecutive vintages (2531 to 2523, and
2018 to 1996), which late reporting cannot produce.

When the forecast carries the confirmed/suspect ward split
(`confirmed_occupancy` / `suspect_occupancy`, the occupancy partitioned by
the cut-off confirmed share), the two ward occupancy levels are archived too
as `treatment beds` (confirmed) and `isolation beds (suspected)`, each scored
against its own Tableau 6 occupancy sub-stock. The total `isolation beds`
occupancy stays as its own stream. Nothing produces those two columns yet:
[`forecast_reported`](@ref) projects the total occupancy alone, so the two
ward entries are reachable only from a forecast a caller has partitioned
itself, which today is the tests. Streams a forecast does not carry are
skipped, so this costs a real forecast nothing. `thin` keeps every
`thin`-th draw so the archive stays compact when it is saved as a release
asset.
"""
function forecast_archive(fcs; made_date::Date, thin::Integer = 1)
    streams = (
        (:confirmed_new, "confirmed cases"),
        (:confirmed_deaths_new, "confirmed deaths"),
        (:recovered_new, "recovered"),
        (:isolation_level, "isolation beds"),
        (:confirmed_occupancy, "treatment beds"),
        (:suspect_occupancy, "isolation beds (suspected)"),
        (:onset_reports_new, "onset reports"))
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
    province_forecast_archive(chn, fcs; made_date, thin = 1) -> DataFrame

Long-format archive of the per-province split of one or more
[`forecast_reported`](@ref) results made from a single cut-off, in the
[`forecast_archive`](@ref) schema plus a `province` column. `chn` is the
patch chain the forecasts were made from, `fcs` an iterable of
`(horizon, fc)` pairs and `made_date` the cut-off `Date`. Returns one row
per `(province, stream, horizon, draw)` with columns `made_date`,
`horizon`, `target_date`, `province`, `stream`, `draw` and `value`.

The two incident streams the spatial tables report are archived, under the
same `confirmed cases` and `confirmed deaths` labels `forecast_archive`
gives the national streams each row is a share of. `province` is the patch
key from [`PROVINCE_NAMES`](@ref), which [`PROVINCE_MEMBERS`](@ref) maps to
the source provinces a patch pools, so a scorer can build each patch's truth
from the per-province histories in the same release's `observations.toml`.

Each province's value is the national draw times that province's modelled
share at the most recent spatial vintage, multiplied draw by draw so the
archived draws carry the correlation between the two factors. The share is
held over the horizon rather than projected forward, so a province whose
share is moving is archived as though it were not. `thin` keeps every
`thin`-th draw so the archive stays compact as a release asset.
"""
function province_forecast_archive(chn, fcs; made_date::Date,
        n_patches::Integer = length(PROVINCE_NAMES),
        patch_labels::AbstractVector = PROVINCE_NAMES,
        thin::Integer = 1)
    np = min(n_patches, length(patch_labels))
    out = DataFrame(made_date = Date[], horizon = Int[], target_date = Date[],
        province = String[], stream = String[], draw = Int[],
        value = Float64[])
    for (horizon, fc) in fcs
        h = Int(horizon)
        target = made_date + Day(h)
        for (label, province, vals) in _province_forecast_draws(
            chn, fc, np, patch_labels)
            for (d, i) in enumerate(1:thin:length(vals))
                push!(out,
                    (made_date, h, target, province, label, d,
                        Float64(vals[i])))
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
    ## Both are levels at the horizon. The gap between them is the shortfall.
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
    onset_forecast_table(fc; digits = 0) -> DataFrame

Summarise a [`forecast_onsets`](@ref) result into the report's usual
30/60/90% credible-interval table, one row per quantity, ordered so the
nowcast is read before the forecast. The `Quantity` labels name what each
row is a statement about rather than which column it came from, since
"onsets" and "onsets reported" are different numbers:

  - `symptom onsets to date` and `of those, reported by T`: the nowcast
    pair. Their difference is the next row.
  - `onsets not yet reported at T`: already happened, not yet in the
    figure (reporting backlog and never-ascertained cases together, see
    [`forecast_onsets`](@ref)).
  - `reports this week of onsets before T`: the backlog the coming week
    should clear.
  - `reports this week of onsets after T`: reports of cases that have not
    yet had their symptom onset.
  - `new onset reports this week`: the replicated total of the two, the
    quantity scored against the triangle's own increment.
  - `new symptom onsets this week`: the latent forecast, unobserved.
"""
function onset_forecast_table(fc::DataFrame; digits::Integer = 0)
    _row(quantity, draws) = begin
        s = posterior_summary(draws)
        (quantity = quantity,
            lower_90 = round(s.lo90; digits), lower_60 = round(s.lo60; digits),
            lower_30 = round(s.lo30; digits), upper_30 = round(s.hi30; digits),
            upper_60 = round(s.hi60; digits), upper_90 = round(s.hi90; digits))
    end
    rows = [
        _row("symptom onsets to date", fc.onsets_to_date),
        _row("of those, reported by T", fc.onset_reports_to_date),
        _row("onsets not yet reported at T", fc.onsets_unreported),
        _row("reports this week of onsets before T",
            fc.onset_reports_backfill),
        _row("reports this week of onsets after T", fc.onset_reports_future),
        _row("new onset reports this week", fc.onset_reports_new),
        _row("new symptom onsets this week", fc.onsets_new)]
    return _prettify(DataFrame(rows))
end

"""
Validate a [`forecast_reported`](@ref) projection against the counts that
were later observed. `observed` is a `NamedTuple` mapping each stream's
cumulative column (`:confirmed_cum`, `:cases_cum`, …) to its observed
cumulative count at the forecast target date. `baseline` maps the same
columns to the cumulative count at the forecast origin (default `0`). A
stream is scored only when its cumulative column is in the forecast and a key
for it is in `observed`. The scored streams are the reported cases, suspected
deaths, laboratory-confirmed cases, confirmed deaths and recovered (the same
set [`plot_forecast_vs_truth`](@ref) draws), so the call site can pass the
same `observed`/`baseline` NamedTuples it builds for the plot.

Each scored stream gets two rows, mirroring the plot's two panels and the
`Quantity` split of [`forecast_table`](@ref): a `cumulative by T+7` row
scoring the projected cumulative against `observed`, and a `new this week`
row scoring the projected new count against
`max(observed − breaks − baseline, 0)`. `breaks` is keyed like `observed`
and carries each stream's retrospective harmonisation correction over the
forecast window (see [`confirmed_break_correction`](@ref)), defaulting to
zero. It comes out of both truths, because it lands in the reported
cumulative without having been notified in that week and a projection cannot
contain it.

When `isolation` (the observed bed occupancy at the target date) is supplied
and the forecast carries the beds, the projected supply-limited occupancy is
scored against it too as a single level row. Returns a `DataFrame` with the
observed count, the equal-tailed 30/60/90% predictive intervals (the same
endpoints as the other summary tables), and whether the observed count falls
inside the 90% interval.

Note that at a one-week-back freeze the bed capacity is weakly informed (the
reported occupancy rate starts only on 9 June), so the projected bed
occupancy rides the capacity random walk back to the freeze date and its
interval is wide.
"""
function forecast_vs_truth(fc::DataFrame;
        observed::NamedTuple, baseline::NamedTuple = NamedTuple(),
        breaks::NamedTuple = NamedTuple(),
        isolation::Union{Real, Missing} = missing,
        digits::Integer = 0)
    _row(label, quantity, draws, obs) = begin
        s = posterior_summary(draws)
        lo = round(s.lo90; digits)
        hi = round(s.hi90; digits)
        (stream = label, quantity = quantity, observed = round(obs; digits),
            lower_90 = lo, lower_60 = round(s.lo60; digits),
            lower_30 = round(s.lo30; digits), upper_30 = round(s.hi30; digits),
            upper_60 = round(s.hi60; digits), upper_90 = hi,
            within_90 = lo <= obs <= hi ? "yes" : "no")
    end
    specs = (
        (:cases_cum, :cases_new, "DRC reported cases"),
        (:deaths_cum, :deaths_new, "DRC suspected deaths"),
        (:confirmed_cum, :confirmed_new, "DRC confirmed cases"),
        (:confirmed_deaths_cum, :confirmed_deaths_new, "DRC confirmed deaths"),
        (:recovered_cum, :recovered_new, "DRC recovered")
    )
    rows = NamedTuple[]
    for (cumcol, newcol, label) in specs
        (cumcol in propertynames(fc) && haskey(observed, cumcol)) || continue
        brk = float(get(breaks, cumcol, 0))
        obs_cum = float(observed[cumcol]) - brk
        push!(rows, _row(label, "cumulative by T+7", fc[!, cumcol], obs_cum))
        newcol in propertynames(fc) || continue
        obs_new = max(obs_cum - float(get(baseline, cumcol, 0)), 0.0)
        push!(rows, _row(label, "new this week", fc[!, newcol], obs_new))
    end
    isolation !== missing && :isolation_level in propertynames(fc) &&
        push!(rows, _row("DRC isolation beds", "occupancy at T+7",
            fc[!, :isolation_level], isolation))
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
## when it carries none. Single-stream fits bind each stream's submodel under
## a composer-level name (`cases_state.expected_reports`), and the joint
## carries those nested names alongside its own un-prefixed aliases
## (`expected_reports_T`), so one ordered list serves both fit kinds.
function _resolve_draws(chn, candidates)
    for key in candidates
        _has_key(chn, key) && return _draws(chn, key)
    end
    return nothing
end

## Per-stream chain-key resolution for [`forecast_stream`](@ref). The joint
## and the single-stream fits do not name these consistently, so the mapping
## is spelled out per stream rather than derived from the stream name.
##
## - `expected`: the stream's cut-off expected count. The nested name comes
##   first because both fit kinds carry it, with the un-prefixed joint alias
##   as the fallback. The names differ in whether they carry a `_T` suffix
##   (`deaths_state.expected_deaths_T` and `exports_state.expected_exports_T`
##   do, `cases_state.expected_reports` does not). `:confirmed_cases` has no
##   nested name at all, since [`confirmed_cases_model`](@ref) keeps its
##   derived quantities on plain `=` rather than `:=` (a `:=` there builds a
##   tracking closure over the boxed `p_pos` that Enzyme cannot
##   differentiate through). Both [`bvd_joint`](@ref) and
##   [`confirmed_only_model`](@ref) alias the cut-off count un-prefixed as
##   `expected_confirmed_T` instead.
## - `dispersion`: each fit's own dispersion, read as it stands. The two
##   candidate sets are disjoint by fit kind. The joint carries only the
##   per-stream aliases (`k_cases`, `isolation_dispersion`, …), since
##   `pooled_dispersion_model` returns its pooled `k` vector rather than
##   `:=`-exposing it and `treatment_flow_model` skips its own dispersion
##   when the joint injects `k_external`. The single-stream fits carry only
##   the nested scalar they sample. The joint alias is listed first so a
##   pooled vector cannot be reached ahead of a stream's scalar, and
##   `forecast_stream` guards on the resolved length as a backstop.
## - `trajectory`: the stream's cumulative trajectory, whose last increment
##   is the cut-off daily rate. The joint and the stream's own single-stream
##   composer expose it under the same un-prefixed `:=` name, so one entry
##   serves both. An empty list means neither fit carries one and the daily
##   rate is inferred from the cut-off cumulative instead.
## - `kind`: `:cumulative` streams report a total accrued to the cut-off, so
##   the forecast is the new count over the horizon. `:level` streams report
##   a prevalence, so the forecast is the level at the horizon.
## - `noise`: the observation model the replicate is drawn through.
##   Exports are Poisson ([`exports_model`](@ref)) and carry no dispersion.
##   Every other stream is negative-binomial.
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
        expected = [:expected_confirmed_T],
        dispersion = [:k_confirmed, Symbol("dispersion_state.k")],
        trajectory = [:cumulative_confirmed],
        kind = :cumulative, noise = :nb),
    :confirmed_deaths => (
        expected = [
            Symbol("confirmed_deaths_state.expected_confirmed_deaths"),
            :expected_confirmed_deaths_T],
        dispersion = [:k_confirmed_deaths, Symbol("dispersion_state.k")],
        trajectory = [:cumulative_confirmed_deaths],
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
        kind = :level, noise = :nb),
    ## The onset stream does not fit this shape. Its forecast is a difference
    ## of two reported totals under the fitted delay hazard rather than a
    ## cut-off expectation grown by a horizon factor, and its noise is the
    ## digitisation scale rather than a count dispersion. The entry exists so
    ## `:onset_reports` is a recognised stream name. The `kind = :onset`
    ## branch in `forecast_stream` hands it to [`forecast_onsets`](@ref)
    ## rather than reading these fields.
    :onset_reports => (
        expected = [:expected_onset_reported_T],
        dispersion = Symbol[],
        trajectory = [:cumulative_onsets],
        kind = :onset, noise = :onset)
)

## Cut-off bed capacity draws. The joint aliases it un-prefixed
## (`bed_capacity := treatment_state.capacity`), but `capacity` is a return
## field of [`treatment_flow_model`](@ref) rather than a `:=`, so a
## standalone treatment fit carries no such key. Its capacity walk exposes
## only its parameters (`cap_state.C0`, `cap_state.z`, `cap_state.σ_cap`).
##
## The capacity is instead recovered from the two deterministics the submodel
## does expose, `expected_isolation := safe_rate(occ_T)` and
## `bed_utilisation := safe_rate(occ_T) / safe_rate(C_T)`. Their ratio is
## `safe_rate(C_T)` exactly, since the occupancy cancels including its
## [`safe_rate`](@ref) flooring. The cancellation holds because `safe_rate`
## floors at `eps`, so the denominator is positive even at the near-zero
## early occupancy where the division would otherwise be ill-posed. A
## non-finite or non-positive result means that invariant has been broken
## upstream. The isolation forecast is supply-limited, so a wrong capacity
## would silently unbound the projected level and this fails loudly instead.
function _bed_capacity(chn)
    _has_key(chn, :bed_capacity) && return _draws(chn, :bed_capacity)
    occ = Symbol("treatment_state.expected_isolation")
    util = Symbol("treatment_state.bed_utilisation")
    (_has_key(chn, occ) && _has_key(chn, util)) || return nothing
    cap = _draws(chn, occ) ./ _draws(chn, util)
    all(c -> isfinite(c) && c > 0, cap) || throw(ArgumentError(
        "forecast_stream: bed capacity recovered from " *
        "`treatment_state.expected_isolation` / " *
        "`treatment_state.bed_utilisation` is not finite and positive for " *
        "every draw; the occupancy did not cancel as expected, so the " *
        "supply limit cannot be trusted."))
    return cap
end

## Chain keys carrying the cumulative occupancy reclassification offset
## standing at the cut-off, `Δ(n)` in [`cumulative_occupancy_offset`](@ref).
## `treatment_flow_model` exposes it as the `occupancy_break` deterministic,
## which both the joint and a standalone treatment fit carry under the
## submodel prefix.
const _OCCUPANCY_OFFSET_KEYS = [Symbol("treatment_state.occupancy_break"),
    :occupancy_break]

## Per-draw cut-off occupancy reclassification offset, or a zero vector for a
## chain that carries none. A fit with no declared occupancy break day pins
## the offset at zero, so both cases mean no reclassification to carry and
## the occupancy is projected on the demand scale.
function _occupancy_offset(chn, nd::Integer)
    d = _resolve_draws(chn, _OCCUPANCY_OFFSET_KEYS)
    (isnothing(d) || length(d) != nd) && return zeros(Float64, nd)
    return d
end

## The reported occupancy a demand replicate implies: the demand shifted by
## the cut-off reclassification offset and capped at the bed capacity, the
## same `min(demand + Δ, C)` the censored occupancy likelihood scores against
## the reported series. Floored at zero, since a reclassification that
## revises the reported basis down cannot take the reported stock below it.
function _reported_occupancy(demand::Integer, offset::Real, capacity::Real)
    return clamp(demand + round(Int, offset), 0, round(Int, capacity))
end

## Cut-off reproduction-number draws. The joint exposes
## `R_T := infection_state.Rt[n]`. Single-stream composers do not, since the
## alias lives in `bvd_joint` rather than the shared `_latent` submodel, so
## `R_T` is rebuilt from the walk parameters every chain carries. That needs
## the grid length `n` and the intervention `breakpoint`, which are data
## rather than chain contents, so it is only possible when the caller
## supplies them.
##
## The intervention ramp is passed explicitly to match the model's own
## default (`sigmoid_ramp` / `rt_walk_model`) rather than taking
## `reconstruct_rt`'s lighter one. The cut-off `Rt[n]` is ramp-sensitive, and
## at the model's own ramp the reconstruction reproduces the joint's `R_T`
## exactly.
function _cutoff_rt(chn; n, breakpoint, rt_start, rt_walk_start)
    _has_key(chn, :R_T) && return _draws(chn, :R_T)
    (isnothing(n) || isnothing(breakpoint)) && return nothing
    rt = reconstruct_rt(chn; n = n, breakpoint = breakpoint,
        rt_start = rt_start, rt_walk_start = rt_walk_start,
        ramp = RT_INTERVENTION_RAMP)
    return Float64[rt[i, n] for i in axes(rt, 1)]
end

## Cut-off growth-rate draws. The joint exposes `r := infection_state.r`.
## Single-stream fits do not, so it is derived the way
## [`infection_model`](@ref) itself derives it, forward Euler–Lotka from the
## cut-off reproduction number through the draw's own generation interval,
## which keeps `r < 0` exactly when `R_T < 1`. The chain's `growth_state.r`
## is not this quantity. That is the cryptic clock rate the joint exposes as
## `r0`, the early growth of the seeded phase rather than the growth at the
## cut-off.
##
## The generation interval is discretised at the same truncation
## `infection_model` uses (`cdf_nmax` of the prior centre, not `_gi_pmf`'s
## lighter default), so the rate reconstructed here matches the `r` the joint
## reads straight off its chain.
function _cutoff_r(chn, R_T)
    _has_key(chn, :r) && return _draws(chn, :r)
    isnothing(R_T) && return nothing
    α = _draws(chn, Symbol("gi_state.α"))
    θ = _draws(chn, Symbol("gi_state.θ"))
    nmax = cdf_nmax(Gamma(2.71, 5.65))
    return Float64[euler_lotka_r(
                       max(R_T[i], _RT_EULER_FLOOR),
                       _gi_pmf(α[i], θ[i]; nmax = nmax))
                   for i in eachindex(R_T)]
end

## Total outbreak-age draws at the cut-off. The joint exposes `T :=
## infection_state.T`. Single-stream fits carry only `growth_state.T`, the
## cryptic duration, to which [`infection_model`](@ref) adds the observed
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

Posterior-predictive forecast of a single observed stream from either a
single-stream fit or the joint, returning one replicate per draw. Where
[`forecast_reported`](@ref) projects every stream the joint carries at once,
this projects one named stream and resolves its chain keys for both fit
kinds, so each single-stream fit can be scored on forecasting its own
dataset against the joint and against a baseline.

`stream` is one of `:reported_cases`, `:suspected_deaths`,
`:confirmed_cases`, `:confirmed_deaths`, `:isolation_beds`, `:exports` and
`:onset_reports`. The incident streams (everything but `:isolation_beds`)
return the new count accrued over the horizon, matching
[`forecast_archive`](@ref)'s convention. `:isolation_beds` returns the
reported occupancy level at the horizon, the projected demand replicate
shifted by the fitted reclassification offset standing at the cut-off and
capped at the bed capacity, `min(demand + Δ, C)` (see
[`cumulative_occupancy_offset`](@ref)).

`:onset_reports` is incident like the rest but projected differently. It is
the new reported count the digitised triangle should add over the horizon,
which [`forecast_onsets`](@ref) builds by differencing two reported totals
under the fitted delay hazard rather than by growing a cut-off expectation.
That needs the fitted triangle's own onset/report-day grid, so
`onset_grid_start` and `onset_grid_end` must be passed for this stream. Of
[`forecast_onsets`](@ref)'s quantities this returns only the scored one.

`obs_value` is the stream's observed count at the cut-off, the cumulative
total for the incident streams or the observed occupancy for
`:isolation_beds`. No returned quantity depends on it: the incident streams
return the increment their cut-off cumulative cancels out of, and the level
stream is projected from the fitted stock. It is kept so a caller states the
base each stream is being projected on top of, and so this signature reads
against [`forecast_reported`](@ref)'s `obs_*` arguments, which do set that
function's cumulative columns.

The reproduction number keeps evolving over the horizon exactly as in
[`forecast_reported`](@ref). The joint carries the cut-off `R_T`, `r` and `T`
as un-prefixed deterministics. A single-stream fit carries none of them, so
they are rebuilt from the reproduction-number walk ([`reconstruct_rt`](@ref))
and the generation interval, which every chain does carry. That needs the
grid length `n` and the intervention `breakpoint` (data rather than chain
contents). Pass both for a single-stream chain, with `rt_start` /
`rt_walk_start` if the fit used non-default ones. Without them a
single-stream chain cannot be projected and an error says so. They are
unnecessary for the joint, which is read directly.

Exports are forecast here for completeness of the per-stream comparison, but
the caveat [`forecast_reported`](@ref) documents still applies. The export
projection assumes a baseline travel rate that cross-border movement is
unlikely to keep once the outbreak is known.
"""
function forecast_stream(chn, stream::Symbol;
        horizon::Real = 7,
        obs_value::Real,
        n::Union{Nothing, Integer} = nothing,
        breakpoint::Union{Nothing, Real} = nothing,
        rt_start::Integer = 1,
        rt_walk_start::Integer = rt_start,
        onset_grid_start::Union{Nothing, Integer} = nothing,
        onset_grid_end::Union{Nothing, Integer} = nothing,
        seed::Integer = 20260520)
    spec = get(_STREAM_SPEC, stream, nothing)
    isnothing(spec) && throw(ArgumentError(
        "forecast_stream: unknown stream `:$stream`; expected one of " *
        join(sort!([":$s" for s in keys(_STREAM_SPEC)]), ", ")))

    ## The onset stream is projected through the fitted reporting hazard
    ## rather than a growth factor, so it is handed straight to
    ## `forecast_onsets` before any of the spec-driven resolution below.
    if spec.kind === :onset
        (isnothing(onset_grid_start) || isnothing(onset_grid_end)) &&
            throw(ArgumentError(
                "forecast_stream: `:onset_reports` needs the fitted " *
                "triangle's own grid. Pass `onset_grid_start` and " *
                "`onset_grid_end` (the minimum onset day and maximum " *
                "report day of the scored cells)."))
        isnothing(n) && throw(ArgumentError(
            "forecast_stream: `:onset_reports` needs the grid length `n` " *
            "to know which onset dates are in the past."))
        fc = forecast_onsets(chn; grid_start = onset_grid_start,
            grid_end = onset_grid_end, n = n, horizon = horizon,
            breakpoint = breakpoint, rt_start = rt_start,
            rt_walk_start = rt_walk_start, seed = seed)
        return fc.onset_reports_new
    end

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

    ## The rate paths and the replicates below draw from one seeded stream,
    ## the same way `forecast_reported` does, so the two stay comparable.
    nh = Int(horizon)
    rng = MersenneTwister(seed)
    evolving = _evolving_rates(chn, nh; R_T = R_T, rng = rng)

    k = nothing
    if spec.noise === :nb
        k = _resolve_draws(chn, spec.dispersion)
        isnothing(k) && throw(ArgumentError(
            "forecast_stream: chain carries no dispersion for `:$stream` " *
            "(tried " * join(["`$key`" for key in spec.dispersion], ", ") *
            ")."))
        ## The joint's nested `dispersion_state.k` is the pooled six-element
        ## vector, so a flattened resolution would silently mismatch the
        ## draws. Every candidate above should be per-draw scalar.
        length(k) == nd || throw(ArgumentError(
            "forecast_stream: dispersion for `:$stream` resolved to " *
            "$(length(k)) values for $nd draws; expected a scalar per draw."))
    end

    out = Vector{Int}(undef, nd)

    ## Replicate a projected mean through the stream's observation model.
    _replicate(i, μ) = spec.noise === :poisson ?
                       rand(rng, Poisson(safe_rate(float(μ)))) :
                       _nb_rand(rng, k[i], μ)
    ## Incident streams are replicated day by day and summed. The dispersion
    ## is fitted at daily and few-day resolution, so a single draw on the
    ## summed mean would carry the overdispersion term at the wrong scale.
    _replicate_new(i, means) = sum(_replicate(i, m) for m in means)

    if spec.kind === :level
        ## The capacity is held at its cut-off value, as `forecast_reported`
        ## holds it, so the projected level is on the reported scale the
        ## truth is read on.
        cap = _bed_capacity(chn)
        isnothing(cap) && throw(ArgumentError(
            "forecast_stream: chain carries no bed capacity (tried " *
            "`bed_capacity` and `treatment_state.expected_isolation` / " *
            "`treatment_state.bed_utilisation`)."))
        occ_offset = _occupancy_offset(chn, nd)
        @inbounds for i in 1:nd
            rs = isnothing(evolving) ? nothing : evolving.paths[i]
            grow = isnothing(rs) ? exp(r[i] * horizon) : prod(exp, rs)
            demand = _replicate(i, expected_T[i] * grow)
            out[i] = _reported_occupancy(demand, occ_offset[i], cap[i])
        end
        return out
    end

    ## Cumulative: project the new count the stream adds over the horizon from
    ## its cut-off daily rate. That rate is the last increment of the stream's
    ## cumulative trajectory when the chain carries it, otherwise inferred
    ## from the cut-off cumulative under exponential growth over the outbreak
    ## age.
    T_age = _outbreak_age(chn, nd; n = n, rt_start = rt_start)
    daily = something(_daily_at_cutoff_any(chn, spec.trajectory),
        _approx_daily.(expected_T, r, T_age))
    @inbounds for i in 1:nd
        rs = isnothing(evolving) ? nothing : evolving.paths[i]
        means = _daily_means(daily[i], rs, r[i], nh)
        ## The new count over the horizon, as `forecast_reported` forms its
        ## `*_new` columns: `obs_value + replicate` less `obs_value`, floored
        ## at zero. `_replicate_new` sums non-negative count draws, so the
        ## cut-off cumulative cancels exactly and the floor never binds. The
        ## increment is taken directly rather than round-tripped through a
        ## cumulative that is added and removed again.
        out[i] = _replicate_new(i, means)
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
