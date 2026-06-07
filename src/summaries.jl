# Posterior summary tables and fit diagnostics: per-parameter
# credible-interval rows, R-hat/ESS/divergence diagnostics, per-stream
# comparisons, and the published-scenario lookup.

_draws(chn, name::Symbol) = vec(Array(chn[name]))

"""
Return `(lo90, lo60, lo30, hi30, hi60, hi90)` equal-tailed credible
interval endpoints from a vector of draws.
"""
function posterior_summary(xs)
    return (
        lo90 = quantile(xs, 0.05),
        lo60 = quantile(xs, 0.20),
        lo30 = quantile(xs, 0.35),
        hi30 = quantile(xs, 0.65),
        hi60 = quantile(xs, 0.80),
        hi90 = quantile(xs, 0.95)
    )
end

## Human-readable headers for the displayed summary tables. Internal
## column keys stay machine-friendly; this maps them to nice labels at
## the point each table is returned.
const _PRETTY_COLS = Dict(
    "quantity" => "Quantity",
    "stream" => "Stream",
    "scenario" => "Scenario",
    "central_estimate" => "Central estimate",
    "reported_cases" => "Reported cases",
    "narrowest_interval" => "Narrowest interval",
    "observed" => "Observed",
    "within_90" => "Within 90% PI",
    "date" => "Date",
    "horizon_days" => "Horizon (days)",
    "lower_90" => "Lower 90%", "lower_60" => "Lower 60%",
    "lower_30" => "Lower 30%", "upper_30" => "Upper 30%",
    "upper_60" => "Upper 60%", "upper_90" => "Upper 90%"
)

_prettify(df::DataFrame) = rename(df, [n => get(_PRETTY_COLS, n, n) for n in names(df)])

"""
`DataFrame` with one row per posterior parameter and the columns
`Quantity, Lower 90%, Lower 60%, Lower 30%, Upper 30%, Upper 60%,
Upper 90%` giving the lower and upper endpoints of the equal-tailed
30%, 60% and 90% credible intervals.
"""
function summary_table(chn, params::AbstractVector{Symbol};
        digits::Integer = 2)
    df = @chain DataFrame(
        quantity = String[],
        lower_90 = Float64[], lower_60 = Float64[],
        lower_30 = Float64[], upper_30 = Float64[],
        upper_60 = Float64[], upper_90 = Float64[]
    ) begin
        let df = _
            for p in params
                s = posterior_summary(_draws(chn, p))
                push!(df,
                    (string(p),
                        round(s.lo90; digits), round(s.lo60; digits),
                        round(s.lo30; digits), round(s.hi30; digits),
                        round(s.hi60; digits), round(s.hi90; digits)))
            end
            df
        end
    end
    return _prettify(df)
end

## --- Fit diagnostics ----------------------------------------------------

# Derived daily latent trajectories carried only for the figures. Their
# early cryptic-phase entries are near-degenerate deterministic functions
# of the seed, so their R-hat / ESS would dominate the headline fit
# summary without reflecting genuine sampler mixing. Excluded from the
# convergence pool; the sampled vectors (random-walk innovations) stay in.
const _DIAGNOSTIC_EXCLUDE = (
    "cumulative_infections", "cumulative_onsets", "cumulative_expected_deaths")

# Flat vector of a scalar diagnostic (R-hat or ESS), one entry per scalar
# parameter in a FlexiChains summary; vector-valued sampled parameters
# contribute one entry per element. Excluded derived trajectories skipped.
function _scalar_stats(summary; exclude = _DIAGNOSTIC_EXCLUDE)
    out = Float64[]
    for p in FlexiChains.parameters(summary)
        ## Vector deterministics surface as indexed scalars
        ## (`cumulative_infections[1]`, ...), so match on the name prefix.
        any(b -> startswith(string(p), b), exclude) && continue
        v = summary[p]
        if v isa Number
            ismissing(v) && continue
            push!(out, Float64(v))
        else
            for x in skipmissing(vec(collect(v)))
                push!(out, Float64(x))
            end
        end
    end
    return out
end

# Number of divergent NUTS transitions recorded in the chain.
function _num_divergences(chn)
    for e in FlexiChains.extras(chn)
        e.name === :numerical_error || continue
        return Int(sum(skipmissing(vec(chn[e]))))
    end
    return 0
end

"""
NUTS fit-quality summary for one chain: the worst (maximum) R-hat and
the smallest bulk effective sample size across parameters, and the
number of divergent transitions.
"""
function fit_diagnostics(chn)
    ## Drop non-finite entries: a fixed or degenerate quantity has an
    ## undefined R-hat / ESS (NaN) that would otherwise mask the worst
    ## genuine value across the sampled parameters.
    rhats = filter(isfinite, _scalar_stats(FlexiChains.rhat(chn)))
    esses = filter(isfinite, _scalar_stats(FlexiChains.ess(chn; kind = :bulk)))
    return (max_rhat = isempty(rhats) ? NaN : maximum(rhats),
        min_ess_bulk = isempty(esses) ? NaN : minimum(esses),
        n_divergent = _num_divergences(chn))
end

"""
`DataFrame` of fit-quality diagnostics with one row per fit. Pass each
fit as `"label" => chain`. Columns `:fit, :max_rhat, :min_ess_bulk,
:divergences`.
"""
function diagnostics_table(fits::Pair{String}...)
    rows = map(fits) do (label, chn)
        d = fit_diagnostics(chn)
        (fit = label,
            max_rhat = round(d.max_rhat; digits = 3),
            min_ess_bulk = round(d.min_ess_bulk; digits = 0),
            divergences = d.n_divergent)
    end
    return DataFrame(rows)
end

"""
Side-by-side credible intervals for `C_T` from several fits. Pass
each fit as `"label" => draws_vector`.
"""
function streams_table(streams::Pair{String, <:AbstractVector}...;
        digits::Integer = 0)
    rows = map(streams) do (label, draws)
        s = posterior_summary(draws)
        (stream = label,
            lower_90 = round(s.lo90; digits), lower_60 = round(s.lo60; digits),
            lower_30 = round(s.lo30; digits), upper_30 = round(s.hi30; digits),
            upper_60 = round(s.hi60; digits), upper_90 = round(s.hi90; digits))
    end
    return _prettify(DataFrame(rows))
end

"""
For each published `C_T` scenario, the narrowest joint posterior
credible interval (30, 60 or 90%) that contains it, or "outside
90%".
"""
function comparison_table(C_draws::AbstractVector;
        scenarios = REPORT_SCENARIOS)
    s = posterior_summary(C_draws)
    rows = map(scenarios) do (label, val)
        crI = if s.lo30 <= val <= s.hi30
            "30%"
        elseif s.lo60 <= val <= s.hi60
            "60%"
        elseif s.lo90 <= val <= s.hi90
            "90%"
        else
            "outside 90%"
        end
        (scenario = label, reported_cases = val,
            narrowest_interval = crI)
    end
    return _prettify(DataFrame(rows))
end
