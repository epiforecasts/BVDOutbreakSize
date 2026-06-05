# One-week-ahead posterior-predictive forecast. Continues the fitted
# renewal trajectory `horizon` days past the cut-off at the current growth
# rate `r` (the daily growth of infections at the cut-off, held constant
# over the short horizon) and scales the fitted expected counts for each
# stream, then replicates an integer count per draw so the intervals carry
# both parameter and observation uncertainty.

function _nb_rand(rng, k, μ)
    μs = max(μ, eps(typeof(μ)))
    p = clamp(k / (k + μs), eps(typeof(k)), one(k) - eps(typeof(k)))
    return rand(rng, NegativeBinomial(k, p))
end

"""
One-week-ahead (default `horizon = 7` days) posterior-predictive
forecast. For each draw, continue the current growth rate `r` over the
horizon and scale the fitted expected counts by `exp(r · horizon)`, then
replicate the cumulative counts. Returns a `DataFrame` with one row per
draw and columns:

- `:cases_cum`, `:deaths_cum` — replicated cumulative suspected reported
  cases and deaths by the cut-off plus the horizon.
- `:confirmed_cum`, `:confirmed_deaths_cum` — laboratory-confirmed case
  and confirmed-death counterparts, present when `obs_confirmed` and
  `obs_confirmed_deaths` are supplied.
- `:cases_new`, … `:confirmed_deaths_new` — new counts over the coming
  week (`*_cum` minus the corresponding observed count at the cut-off,
  floored at zero).
- `:infections_new` — new latent infections projected over the horizon at
  the constant cut-off growth rate (a deterministic per-draw quantity, so
  it carries parameter uncertainty only).
- `:rt_forecast` — the reproduction number held over the horizon (the
  terminal `R_T`).

Reads `:r`, `:expected_reports_T`, `:expected_deaths_T`,
`:expected_infections_T`, `:R_T`, `:k` and (for the laboratory streams)
`:expected_confirmed_T` and `:expected_confirmed_deaths_T` from `chn`.
Exports are not forecast:
cross-border travel is unlikely to continue at its baseline rate, so the
forward travel rate the export model relies on no longer holds. Assumes
the current growth rate continues unchanged over the horizon (no further
interventions, no saturation).
"""
function forecast_reported(chn;
        horizon::Real = 7,
        obs_cases::Real,
        obs_deaths::Real,
        obs_confirmed::Union{Real, Missing} = missing,
        obs_confirmed_deaths::Union{Real, Missing} = missing,
        seed::Integer = 20260520)
    r = _draws(chn, :r)
    cases_T = _draws(chn, :expected_reports_T)
    deaths_T = _draws(chn, :expected_deaths_T)
    infections_T = _draws(chn, :expected_infections_T)
    rt_forecast = _draws(chn, :R_T)
    k = _draws(chn, :k)
    has_conf = obs_confirmed !== missing
    has_conf_deaths = obs_confirmed_deaths !== missing
    conf_T = has_conf ? _draws(chn, :expected_confirmed_T) : nothing
    conf_deaths_T = has_conf_deaths ?
                    _draws(chn, :expected_confirmed_deaths_T) : nothing

    rng = MersenneTwister(seed)
    n = length(r)
    cases_cum = Vector{Int}(undef, n)
    deaths_cum = Vector{Int}(undef, n)
    infections_new = Vector{Float64}(undef, n)
    confirmed_cum = has_conf ? Vector{Int}(undef, n) : nothing
    confirmed_deaths_cum = has_conf_deaths ? Vector{Int}(undef, n) : nothing

    @inbounds for i in 1:n
        grow = exp(r[i] * horizon)
        cases_cum[i] = _nb_rand(rng, k[i], cases_T[i] * grow)
        deaths_cum[i] = _nb_rand(rng, k[i], deaths_T[i] * grow)
        ## New latent infections over the horizon, continuing the cut-off
        ## daily infections `I_T` at the constant daily growth `r`: the
        ## geometric sum `I_T Σ_{d=1}^{h} e^{r d}`, which tends to `I_T · h`
        ## as `r → 0`. Latent, so carries parameter (not observation)
        ## uncertainty across draws.
        er = exp(r[i])
        h = float(horizon)
        infections_new[i] = abs(er - 1) < 1e-8 ? infections_T[i] * h :
                            infections_T[i] * er * (exp(r[i] * h) - 1) /
                            (er - 1)
        has_conf && (confirmed_cum[i] = _nb_rand(rng, k[i], conf_T[i] * grow))
        ## Confirmed deaths grow with the suspected-death signal but are
        ## bounded by it (a thinning), so cap the replicate at the forecast
        ## cumulative suspected deaths.
        if has_conf_deaths
            μ = conf_deaths_T[i] * grow
            confirmed_deaths_cum[i] = min(deaths_cum[i],
                _nb_rand(rng, k[i], μ))
        end
    end

    _new(cum, obs) = max.(cum .- round(Int, obs), 0)
    df = DataFrame(
        cases_cum = cases_cum,
        deaths_cum = deaths_cum,
        cases_new = _new(cases_cum, obs_cases),
        deaths_new = _new(deaths_cum, obs_deaths),
        infections_new = infections_new,
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
    return df
end

"""
Summarise a [`forecast_reported`](@ref) result into a `DataFrame` with
one row per stream (suspected cases and deaths, plus confirmed cases and
confirmed deaths when present) and quantity (cumulative total by the
cut-off plus the horizon, or new this week), reporting the equal-tailed
30/60/90% credible interval endpoints (`lower_90 … upper_90`) used by the
other summary tables.
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
    streams = Tuple{String, Symbol, Symbol}[
    (
        "DRC reported cases", :cases_cum, :cases_new),
    (
        "DRC deaths", :deaths_cum, :deaths_new)]
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
    return _prettify(DataFrame(rows))
end

"""
Validate a [`forecast_reported`](@ref) projection against the counts that
were later observed. `cases` and `deaths` are the observed cumulative DRC
suspected reported cases and deaths at the forecast target date;
`confirmed` and `confirmed_deaths` add the laboratory-confirmed streams
when the forecast carries them. Returns a `DataFrame` with one row per
stream giving the observed count, the equal-tailed 30/60/90% predictive
intervals (the same endpoints as the other summary tables), and whether
the observed count falls inside the 90% interval.
"""
function forecast_vs_truth(fc::DataFrame;
        cases::Real, deaths::Real,
        confirmed::Union{Real, Missing} = missing,
        confirmed_deaths::Union{Real, Missing} = missing,
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
    rows = [
        _row("DRC reported cases", fc[!, :cases_cum], cases),
        _row("DRC deaths", fc[!, :deaths_cum], deaths)
    ]
    confirmed !== missing && :confirmed_cum in propertynames(fc) &&
        push!(rows, _row("DRC confirmed cases", fc[!, :confirmed_cum],
            confirmed))
    confirmed_deaths !== missing &&
        :confirmed_deaths_cum in propertynames(fc) &&
        push!(rows, _row("DRC confirmed deaths", fc[!, :confirmed_deaths_cum],
            confirmed_deaths))
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
