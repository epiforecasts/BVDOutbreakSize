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
## column keys stay machine-friendly. This maps them to nice labels at
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
    "upper_60" => "Upper 60%", "upper_90" => "Upper 90%",
    "n" => "Vintages", "bias" => "Bias",
    "coverage_50" => "50% coverage", "coverage_90" => "90% coverage"
)

_prettify(df::DataFrame) = rename(df, [n => get(_PRETTY_COLS, n, n) for n in names(df)])

"""
`DataFrame` with one row per posterior parameter and the columns
`Quantity, Lower 90%, Lower 60%, Lower 30%, Upper 30%, Upper 60%,
Upper 90%` giving the lower and upper endpoints of the equal-tailed
30%, 60% and 90% credible intervals.

`labels` is an optional map from the raw chain symbol to a clean display
name (e.g. `Symbol("rt_state.sigma_rw") => "Rt step size"`), applied to the
`Quantity` column only. The model's variable names are unchanged. Symbols
absent from the map keep their raw name.
"""
function summary_table(chn, params::AbstractVector{Symbol};
        digits::Integer = 2,
        labels::AbstractDict = Dict{Symbol, String}())
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
                    (get(labels, p, string(p)),
                        round(s.lo90; digits), round(s.lo60; digits),
                        round(s.lo30; digits), round(s.hi30; digits),
                        round(s.hi60; digits), round(s.hi90; digits)))
            end
            df
        end
    end
    return _prettify(df)
end

## --- Markdown rendering --------------------------------------------------

# Format one cell for a markdown table: integer-valued floats print without
# a trailing `.0` so count columns read as whole numbers, everything else
# prints with its default string form.
_md_cell(x::Real) = isinteger(x) ? string(Integer(x)) : string(x)
_md_cell(x) = string(x)

"""
GitHub-flavoured markdown table for a `DataFrame`, one header row from the
column names and one body row per data row. Used to persist a rendered
summary table to disk so a static documentation page can embed it without
re-running the fit.
"""
function markdown_table(df::DataFrame)
    cols = names(df)
    header = "| " * join(cols, " | ") * " |"
    sep = "| " * join(fill("---", length(cols)), " | ") * " |"
    rows = map(eachrow(df)) do row
        "| " * join((_md_cell(row[c]) for c in cols), " | ") * " |"
    end
    return join(vcat(header, sep, rows), "\n") * "\n"
end

## --- Fit diagnostics ----------------------------------------------------

# Derived daily latent trajectories carried only for the figures. Their
# early cryptic-phase entries are near-degenerate deterministic functions
# of the seed, so their R-hat / ESS would dominate the headline fit
# summary without reflecting genuine sampler mixing. Excluded from the
# convergence pool. The sampled vectors (random-walk innovations) stay in.
const _DIAGNOSTIC_EXCLUDE = (
    "cumulative_infections", "cumulative_onsets", "cumulative_expected_deaths")

# Flat vector of a scalar diagnostic (R-hat or ESS), one entry per scalar
# parameter in a FlexiChains summary. Vector-valued sampled parameters
# contribute one entry per element. Trajectories in `exclude` are skipped.
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
Per-date posterior summary of the latent symptom-onset trajectory, the
"symptomatic cases" curve plotted to show the outbreak over time. One row
per grid day from `seeding` (grid day 1) to the cut-off (grid day `n`),
giving the equal-tailed 30%, 60% and 90% credible intervals (the same
intervals as [`summary_table`](@ref) and [`streams_table`](@ref)) of both
the daily new symptom onsets and the cumulative symptom onsets to that date.
The chain must carry the vector deterministic `cumulative_onsets` (one
trajectory per draw), as the joint fit does.

Columns: `date`, then for each of `new_onsets` and `cumulative_onsets` the
six endpoints `_lower_90, _lower_60, _lower_30, _upper_30, _upper_60,
_upper_90`.
"""
function onsets_over_time(chn; n::Integer, seeding::Date)
    ## Per-draw cumulative-onset trajectories (the plotted ribbons), then the
    ## per-draw daily new onsets: the first grid day carries the seed
    ## cumulative, later days the day-on-day increment.
    cumulative = [collect(v) for v in vec(collect(chn[:cumulative_onsets]))]
    daily = [vcat(c[1], diff(c)) for c in cumulative]
    ## Prefix a `posterior_summary` NamedTuple's fields with the quantity name
    ## so each day's row carries both onset series side by side.
    function _bounds(prefix, xs)
        s = posterior_summary(xs)
        cols = (Symbol(prefix, "_lower_90"), Symbol(prefix, "_lower_60"),
            Symbol(prefix, "_lower_30"), Symbol(prefix, "_upper_30"),
            Symbol(prefix, "_upper_60"), Symbol(prefix, "_upper_90"))
        return NamedTuple{cols}((s.lo90, s.lo60, s.lo30, s.hi30, s.hi60, s.hi90))
    end
    function _row(d)
        new = Float64[v[d] for v in daily]
        cum = Float64[c[d] for c in cumulative]
        merge((date = seeding + Day(d - 1),),
            _bounds("new_onsets", new), _bounds("cumulative_onsets", cum))
    end
    return DataFrame([_row(d) for d in 1:n])
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

## --- Posterior-predictive calibration -----------------------------------

"""
Sample-based forecast bias for a single observation against a predictive
sample, following `scoringutils::bias_sample`. With `n_lt` and `n_eq` the
counts of predictive draws strictly below and exactly equal to the
observation, the bias is

``1 - 2\\,(n_{lt} + n_{eq}/2)/m``

over the `m` draws. It lies in ``[-1, 1]``: negative when the predictive
distribution sits below the observation (under-prediction), positive when
it sits above (over-prediction), and zero when the observation falls at the
predictive median. The equal-count term handles ties in count data so a
mass of draws exactly at the observation does not bias the score.
"""
function bias_sample(observed::Real, predicted::AbstractVector{<:Real})
    m = length(predicted)
    m == 0 && return NaN
    n_lt = count(<(observed), predicted)
    n_eq = count(==(observed), predicted)
    return 1 - 2 * (n_lt + n_eq / 2) / m
end

# Whether `observed` lies inside the central equal-tailed `level` predictive
# interval of the sample (e.g. `level = 0.9` → the 5–95% interval).
function _covered(observed::Real, predicted::AbstractVector{<:Real}, level::Real)
    lo = quantile(predicted, (1 - level) / 2)
    hi = quantile(predicted, (1 + level) / 2)
    return lo <= observed <= hi
end

# Per-vintage conditional predictive samples for one PPC panel, mirroring
# `plot_vintage_conditional_ppc`: each cumulative-stream draw at vintage `v`
# is the observed previous cumulative plus the drawn increment (baseline
# zero for a `cumulative = false` daily panel), and the matching observed
# value is the cumulative (or daily) count at that vintage. Returns
# `(samples, observed)` with `samples[v]` the draw vector at vintage `v`.
function _panel_conditional(panel)
    observed = float.(panel.observed)
    n = length(observed)
    cumulative = get(panel, :cumulative, true)
    obs_prev = cumulative ?
               [v == 1 ? 0.0 : observed[v - 1] for v in 1:n] : zeros(n)
    replicates = [collect(r) for r in vec(collect(panel.replicates))]
    samples = [[obs_prev[v] + r[v] for r in replicates] for v in 1:n]
    return (samples = samples, observed = observed)
end

"""
Per-stream posterior-predictive calibration for the per-vintage
one-step-ahead conditional checks. Pass the same `panels` given to
[`plot_vintage_conditional_ppc`](@ref): each is a `NamedTuple`
`(; title, observed, replicates, …)` with optional `cumulative`. For each
stream the conditional predictive at every vintage is scored against the
observed count, and the per-vintage scores are averaged into one row.

Columns: `stream`, the number of scored vintages `n`, the mean forecast
`bias` (see [`bias_sample`](@ref): negative means the stream is
under-predicted, positive means over-predicted), and the empirical
`coverage_50`/`coverage_90`: the fraction of vintages whose observed count
falls inside the central 50% and 90% predictive intervals. A well-calibrated
stream has bias near zero and coverage near its nominal level. Departures
flag the streams the joint fit reproduces less well.
"""
function stream_calibration(panels::AbstractVector)
    rows = map(panels) do panel
        c = _panel_conditional(panel)
        n = length(c.observed)
        biases = [bias_sample(c.observed[v], c.samples[v]) for v in 1:n]
        cov50 = [_covered(c.observed[v], c.samples[v], 0.5) for v in 1:n]
        cov90 = [_covered(c.observed[v], c.samples[v], 0.9) for v in 1:n]
        (stream = panel.title, n = n,
            bias = round(n == 0 ? NaN : mean(biases); digits = 2),
            coverage_50 = round(n == 0 ? NaN : mean(cov50); digits = 2),
            coverage_90 = round(n == 0 ? NaN : mean(cov90); digits = 2))
    end
    return _prettify(DataFrame(rows))
end

## Whether a chain carries a given key. Chain types throw on a missing key
## rather than returning a sentinel, so presence has to be probed.
function _has_key(chn, key::Symbol)
    try
        chn[key]
        return true
    catch
        return false
    end
end

"""
Per-patch outbreak summary for the patch model: one row per province, with
the cut-off cumulative infections `C_T`, the cut-off reproduction number
`R_T`, the daily infections at the cut-off, and the log-Rt deviation `δ`
from the common national trend. Each is reported as the same 90/60/30%
credible intervals [`summary_table`](@ref) uses.

The deviations are SUM-TO-ZERO contrasts around the national trend (see
[`patch_rt_model`](@ref)), so `δ` is read relative to the national average
across provinces, not relative to any one patch: a negative `δ` means that
province transmits below the national trend. Every patch, including the
primary, carries its own deviation, and they sum to zero in every draw.
For the log-Rt of a province relative to ITURI specifically, read the
chain's `log_rt_contrast` instead.

Expects a chain from [`bvd_joint`](@ref), which stores the per-patch
quantities as vector deterministics (`C_T_patch`, `R_T_patch`,
`infections_T_patch`, `delta_patch`), one entry per patch.
"""
function patch_summary_table(chn, n_patches::Integer = length(PROVINCE_NAMES);
        digits::Integer = 2,
        patch_labels::AbstractVector = ["Ituri", "Nord-Kivu", "Sud-Kivu"])
    required = [:C_T_patch, :R_T_patch, :infections_T_patch, :delta_patch]
    absent = filter(p -> !_has_key(chn, p), required)
    isempty(absent) || error(
        "chain is missing the per-patch deterministics $(absent); it was " *
        "not sampled from `bvd_joint`.")
    np = min(n_patches, length(patch_labels))
    ## `_draw_vectors` gives one vector per draw; transpose to per-patch
    ## draw vectors so each patch can be summarised independently.
    per_patch(sym) =
        let vs = _draw_vectors(chn, sym)
            [[v[p] for v in vs] for p in 1:np]
        end
    C_T = per_patch(:C_T_patch)
    R_T = per_patch(:R_T_patch)
    inf_T = per_patch(:infections_T_patch)
    δ = per_patch(:delta_patch)
    ## Case ascertainment and the Rt contrast against the primary patch are the
    ## two quantities that must be read TOGETHER: the case composition
    ## identifies only their product, and it is the per-province deaths that
    ## tilt the balance between them. Reporting one without the other invites
    ## a low provincial Rt to be read as epidemiology when it is case-finding.
    asc = _has_key(chn, :province_ascertainment) ?
          per_patch(:province_ascertainment) : nothing
    contrast = _has_key(chn, :log_rt_contrast) ?
               per_patch(:log_rt_contrast) : nothing
    df = DataFrame(
        patch = String[],
        quantity = String[],
        lower_90 = Float64[], lower_60 = Float64[], lower_30 = Float64[],
        upper_30 = Float64[], upper_60 = Float64[], upper_90 = Float64[]
    )
    for p in 1:np
        rows = Any[("Cumulative infections", C_T[p], 0),
            ("Reproduction number", R_T[p], digits),
            ("Daily infections at cut-off", inf_T[p], 0),
            ("log-Rt deviation from trend", δ[p], digits)]
        contrast === nothing ||
            push!(rows, ("log-Rt vs primary patch", contrast[p], digits))
        asc === nothing ||
            push!(rows, ("Relative case ascertainment", asc[p], digits))
        for (label, draws, dg) in rows
            s = posterior_summary(draws)
            push!(df,
                (patch_labels[p], label,
                    round(s.lo90; digits = dg), round(s.lo60; digits = dg),
                    round(s.lo30; digits = dg), round(s.hi30; digits = dg),
                    round(s.hi60; digits = dg), round(s.hi90; digits = dg)))
        end
    end
    return _prettify(rename(df, [:patch => "Patch", :quantity => "Quantity"]))
end
