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
        lo60 = quantile(xs, 0.2),
        lo30 = quantile(xs, 0.35),
        hi30 = quantile(xs, 0.65),
        hi60 = quantile(xs, 0.8),
        hi90 = quantile(xs, 0.95),
    )
end

## Quantities reported as a transform of another parameter rather than from
## their own draws, because their own quantiles do not bound them.
const _DERIVED_FROM = Dict{Symbol, Tuple{Symbol, Function}}(
    :doubling_time => (:r, doubling_time)
)

## Interval endpoints for one reported quantity. A caller may pass a chain
## carrying a derived quantity without the parameter it is derived from, so
## fall back to its own draws there.
function _summary_for(chn, p::Symbol)
    haskey(_DERIVED_FROM, p) ||
        return posterior_summary(_scalar_draws(chn, p))
    src, f = _DERIVED_FROM[p]
    src_draws = try
        _draws(chn, src)
    catch
        return posterior_summary(_scalar_draws(chn, p))
    end
    return map(f, posterior_summary(src_draws))
end

## Draws of a scalar quantity. A vector-valued deterministic reaches here as
## a vector of vectors, and `quantile` then fails with a cryptic
## `isfinite(::Vector)` MethodError, so check and name the real problem.
function _scalar_draws(chn, p::Symbol)
    d = _draws(chn, p)
    eltype(d) <: Number || error(
        "summary_table: `$(p)` is vector-valued, not a scalar. " *
            "Per-patch quantities go in `patch_summary_table`; daily " *
            "trajectories are read with `_draw_vectors`."
    )
    return d
end

## Display headers for the summary tables. Internal column keys stay
## machine-friendly and are mapped here as each table is returned.
const _PRETTY_COLS = Dict(
    "quantity" => "Quantity",
    "province" => "Province",
    "stream" => "Stream",
    "scenario" => "Scenario",
    "estimate" => "Estimate",
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
Upper 90%` giving the endpoints of the equal-tailed 30%, 60% and 90%
credible intervals.

`doubling_time` is reported as the image of `r`'s interval rather than from
its own draws, because it is unbounded at zero growth. Its row runs from the
fastest decline through the zero-growth pole to the fastest growth, and is
not sorted by value.

`labels` maps a raw chain symbol to a display name (e.g.
`Symbol("rt_state.sigma_rw") => "Rt step size"`), applied to the `Quantity`
column only. Symbols absent from the map keep their raw name.
"""
function summary_table(
        chn, params::AbstractVector{Symbol};
        digits::Integer = 2,
        labels::AbstractDict = Dict{Symbol, String}()
    )
    df = @chain DataFrame(
        quantity = String[],
        lower_90 = Float64[], lower_60 = Float64[],
        lower_30 = Float64[], upper_30 = Float64[],
        upper_60 = Float64[], upper_90 = Float64[]
    ) begin
        let df = _
            for p in params
                s = _summary_for(chn, p)
                push!(
                    df,
                    (
                        get(labels, p, string(p)),
                        round(s.lo90; digits), round(s.lo60; digits),
                        round(s.lo30; digits), round(s.hi30; digits),
                        round(s.hi60; digits), round(s.hi90; digits),
                    )
                )
            end
            df
        end
    end
    return _prettify(df)
end

## --- Markdown rendering --------------------------------------------------

# Trim a float to a readable column width. Values at or above one keep their
# magnitude and lose the tail. Values below one keep their leading digits, so
# a small score does not round away to zero.
_md_round(x::AbstractFloat) = abs(x) < 1 ? round(x; sigdigits = 3) :
    round(x; digits = 3)

# Format one cell for a markdown table. Integer-valued floats print without a
# trailing `.0`, and a literal `|` is escaped so it cannot split the row.
_md_cell(x::Real) = isinteger(x) ? string(Integer(x)) : string(x)
_md_cell(x::AbstractFloat) = isinteger(x) ? string(Integer(x)) :
    string(_md_round(x))
_md_cell(x) = replace(string(x), "|" => "\\|")

# Right-align numeric columns and left-align everything else. `Bool` is a
# `Real` but reads as a label rather than a quantity, so it stays left.
function _md_align(col)
    return eltype(col) <: Union{Missing, Bool} ? "---" :
        eltype(col) <: Union{Missing, Real} ? "---:" : "---"
end

"""
GitHub-flavoured markdown table for a `DataFrame`, one header row from the
column names and one body row per data row. Used to persist a rendered table
to disk so a static page can embed it without re-running the fit, and by
[`MarkdownTable`](@ref) to put a table on a Literate page.
"""
function markdown_table(df::DataFrame)
    cols = names(df)
    header = "| " * join(cols, " | ") * " |"
    sep = "| " * join((_md_align(df[!, c]) for c in cols), " | ") * " |"
    rows = map(eachrow(df)) do row
        "| " * join((_md_cell(row[c]) for c in cols), " | ") * " |"
    end
    return join(vcat(header, sep, rows), "\n") * "\n"
end

"""
A table for a Literate page to render as a table rather than as a block of
printed output. Display it as the last expression of a chunk.

Literate picks a chunk's output format from what the value is showable as,
taking `text/html` before `text/markdown`. A bare `DataFrame` is
html-showable, so it goes out as a `@raw html` block, and Documenter compiles
a regex from each raw block's own text, which fails once the block passes
PCRE's ~64KB compiled-pattern limit. Several of the scoring tables cross it.
This wrapper is showable as markdown and not as html, so the table goes out
as an ordinary markdown table.

A `DataFrame` is rendered by [`markdown_table`](@ref). Anything else is taken
through its own markdown rendering, or its plain-text one where it has none,
so a page can pass either a table or a placeholder message.
"""
struct MarkdownTable
    text::String
    ## Inner constructor so Julia generates no outer one: its automatic
    ## `MarkdownTable(::Any)` converting constructor would collide with the
    ## fallback method below.
    MarkdownTable(text::AbstractString) = new(text)
end

MarkdownTable(df::DataFrame) = MarkdownTable(markdown_table(df))
function MarkdownTable(x)
    return MarkdownTable(
        showable(MIME("text/markdown"), x) ?
            sprint(show, MIME("text/markdown"), x) :
            sprint(show, MIME("text/plain"), x)
    )
end

Base.show(io::IO, t::MarkdownTable) = print(io, t.text)
Base.show(io::IO, ::MIME"text/markdown", t::MarkdownTable) = print(io, t.text)

## --- Fit diagnostics ----------------------------------------------------

# Derived daily latent trajectories carried only for the figures. Their early
# cryptic-phase entries are near-degenerate functions of the seed, so their
# R-hat / ESS would dominate the headline fit summary without reflecting
# sampler mixing. The sampled random-walk innovations stay in.
const _DIAGNOSTIC_EXCLUDE = (
    "cumulative_infections", "cumulative_onsets", "cumulative_expected_deaths",
)

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

# The sampler's per-draw divergence flag, or `nothing` where the chain
# carries no sampler extras. One lookup behind both counts below, so they
# cannot come to be taken over different sets of draws.
function _numerical_error_flags(chn)
    for e in FlexiChains.extras(chn)
        e.name === :numerical_error || continue
        return vec(chn[e])
    end
    return nothing
end

# Divergent transitions across every chain.
function _num_divergences(chn)
    flags = _numerical_error_flags(chn)
    flags === nothing && return 0
    return Int(sum(skipmissing(flags)))
end

# Post-warmup draws across every chain. A divergence count means little
# without it: six divergences in 3200 draws and six in twelve are not the
# same fit. Draws whose flag is missing are left out, as they are from the
# count above, so the two are over the same set and their ratio is the
# divergence rate among the draws whose outcome is known.
function _num_draws(chn)
    flags = _numerical_error_flags(chn)
    flags === nothing && return 0
    return count(!ismissing, flags)
end

"""
NUTS fit-quality summary for one chain: the worst (maximum) R-hat, the
smallest bulk and tail effective sample sizes across parameters, the number
of divergent transitions and the number of post-warmup draws they are out of.
"""
function fit_diagnostics(chn)
    ## Drop non-finite entries: a fixed or degenerate quantity has an
    ## undefined R-hat / ESS (NaN) that would otherwise mask the worst
    ## genuine value across the sampled parameters.
    rhats = filter(isfinite, _scalar_stats(FlexiChains.rhat(chn)))
    bulk = filter(isfinite, _scalar_stats(FlexiChains.ess(chn; kind = :bulk)))
    tail = filter(isfinite, _scalar_stats(FlexiChains.ess(chn; kind = :tail)))
    return (
        max_rhat = isempty(rhats) ? NaN : maximum(rhats),
        min_ess_bulk = isempty(bulk) ? NaN : minimum(bulk),
        min_ess_tail = isempty(tail) ? NaN : minimum(tail),
        n_divergent = _num_divergences(chn),
        n_draws = _num_draws(chn),
    )
end

"""
`DataFrame` of fit-quality diagnostics with one row per fit. Pass each
fit as `"label" => chain`. Columns `:fit, :max_rhat, :min_ess_bulk,
:divergences`.
"""
function diagnostics_table(fits::Pair{String}...)
    rows = map(fits) do (label, chn)
        d = fit_diagnostics(chn)
        (
            fit = label,
            max_rhat = round(d.max_rhat; digits = 3),
            min_ess_bulk = round(d.min_ess_bulk; digits = 0),
            divergences = d.n_divergent,
        )
    end
    return DataFrame(rows)
end

"""
Side-by-side credible intervals for `C_T` from several fits. Pass
each fit as `"label" => draws_vector`.
"""
function streams_table(
        streams::Pair{String, <:AbstractVector}...;
        digits::Integer = 0
    )
    rows = map(streams) do (label, draws)
        s = posterior_summary(draws)
        (
            stream = label,
            lower_90 = round(s.lo90; digits), lower_60 = round(s.lo60; digits),
            lower_30 = round(s.lo30; digits), upper_30 = round(s.hi30; digits),
            upper_60 = round(s.hi60; digits), upper_90 = round(s.hi90; digits),
        )
    end
    return _prettify(DataFrame(rows))
end

"""
Per-date posterior summary of the latent symptom-onset trajectory, one row
per grid day from `seeding` (grid day 1) to the cut-off (grid day `n`).
Each row gives the equal-tailed 30%, 60% and 90% credible intervals of both
the daily new symptom onsets and the cumulative symptom onsets to that date.
The chain must carry the vector deterministic `cumulative_onsets`, one
trajectory per draw, as the joint fit does.

Columns: `date`, then for each of `new_onsets` and `cumulative_onsets` the
six endpoints `_lower_90, _lower_60, _lower_30, _upper_30, _upper_60,
_upper_90`.
"""
function onsets_over_time(chn; n::Integer, seeding::Date)
    ## Per-draw cumulative trajectories, then the per-draw daily new onsets.
    ## The first grid day carries the seed cumulative, later days the
    ## day-on-day increment.
    cumulative = [collect(v) for v in vec(collect(chn[:cumulative_onsets]))]
    daily = [vcat(c[1], diff(c)) for c in cumulative]
    ## Prefix the interval fields with the quantity name so each day's row
    ## carries both onset series side by side.
    function _bounds(prefix, xs)
        s = posterior_summary(xs)
        cols = (
            Symbol(prefix, "_lower_90"), Symbol(prefix, "_lower_60"),
            Symbol(prefix, "_lower_30"), Symbol(prefix, "_upper_30"),
            Symbol(prefix, "_upper_60"), Symbol(prefix, "_upper_90"),
        )
        return NamedTuple{cols}((s.lo90, s.lo60, s.lo30, s.hi30, s.hi60, s.hi90))
    end
    function _row(d)
        new = Float64[v[d] for v in daily]
        cum = Float64[c[d] for c in cumulative]
        return merge(
            (date = seeding + Day(d - 1),),
            _bounds("new_onsets", new), _bounds("cumulative_onsets", cum)
        )
    end
    return DataFrame([_row(d) for d in 1:n])
end

"""
For each published `C_T` scenario, the narrowest joint posterior
credible interval (30, 60 or 90%) that contains it, or "outside
90%".
"""
function comparison_table(
        C_draws::AbstractVector;
        scenarios = REPORT_SCENARIOS
    )
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
        (
            scenario = label, reported_cases = val,
            narrowest_interval = crI,
        )
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

over the `m` draws. It lies in ``[-1, 1]``. Negative means the predictive
distribution sits below the observation (under-prediction), positive that it
sits above, and zero that the observation falls at the predictive median. The
equal-count term handles ties in count data, so a mass of draws exactly at
the observation does not bias the score.
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
# `plot_vintage_conditional_ppc`. Each cumulative-stream draw at vintage `v`
# is the observed previous cumulative plus the drawn increment, with a zero
# baseline for a `cumulative = false` daily panel. Returns
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
`bias` from [`bias_sample`](@ref), and the empirical
`coverage_50`/`coverage_90`, the fraction of vintages whose observed count
falls inside the central 50% and 90% predictive intervals. A well-calibrated
stream has bias near zero and coverage near its nominal level.
"""
function stream_calibration(panels::AbstractVector)
    rows = map(panels) do panel
        c = _panel_conditional(panel)
        n = length(c.observed)
        biases = [bias_sample(c.observed[v], c.samples[v]) for v in 1:n]
        cov50 = [_covered(c.observed[v], c.samples[v], 0.5) for v in 1:n]
        cov90 = [_covered(c.observed[v], c.samples[v], 0.9) for v in 1:n]
        (
            stream = panel.title, n = n,
            bias = round(n == 0 ? NaN : mean(biases); digits = 2),
            coverage_50 = round(n == 0 ? NaN : mean(cov50); digits = 2),
            coverage_90 = round(n == 0 ? NaN : mean(cov90); digits = 2),
        )
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

## Per-patch draw vectors for a vector deterministic: `_draw_vectors` gives one
## vector per draw, so transpose to one draw vector per patch.
function _per_patch(chn, sym::Symbol, np::Integer)
    vs = _draw_vectors(chn, sym)
    return [[v[p] for v in vs] for p in 1:np]
end

## Median and 90% credible interval as one cell, `median (lower–upper)`. The
## cross-province overview puts several quantities side by side, so it trades
## the six-column interval layout for one column per quantity.
function _median_ci(draws; digits::Integer = 2)
    fmt(x) = digits <= 0 ? string(round(Int, x)) : string(round(x; digits))
    return string(
        fmt(median(draws)), " (", fmt(quantile(draws, 0.05)), "–",
        fmt(quantile(draws, 0.95)), ")"
    )
end

"""
Cross-province overview for the patch model: one row per province and one
column per quantity, each a median with a 90% credible interval. Reads the
reproduction number at the cut-off, cumulative infections, the province's
share of national infections, and its case ascertainment relative to the
national average. The per-province detail, with the full 30/60/90% intervals
and the deviation parameters, is in [`patch_summary_table`](@ref).

The infection share is computed per draw before summarising, so its interval
carries the correlation between provinces rather than dividing two
independently summarised numbers. Ascertainment and the reproduction number
must be read together. The case composition identifies only their product,
and it is the per-province deaths that tilt the balance between them.
"""
function patch_overview_table(
        chn, n_patches::Integer = length(PROVINCE_NAMES);
        digits::Integer = 2,
        patch_labels::AbstractVector = PROVINCE_LABELS
    )
    required = [:C_T_patch, :R_T_patch]
    absent = filter(p -> !_has_key(chn, p), required)
    isempty(absent) || error(
        "chain is missing the per-patch deterministics $(absent); it was " *
            "not sampled from `bvd_joint`."
    )
    np = min(n_patches, length(patch_labels))
    C_T = _per_patch(chn, :C_T_patch, np)
    R_T = _per_patch(chn, :R_T_patch, np)
    ## Share per draw, so the interval reflects that the provinces' shares are
    ## constrained to sum to one rather than varying independently.
    totals = sum(C_T)
    share = [100 .* C_T[p] ./ totals for p in 1:np]
    asc = _has_key(chn, :province_ascertainment) ?
        _per_patch(chn, :province_ascertainment, np) : nothing
    df = DataFrame(
        "Province" => String[],
        "Reproduction number" => String[],
        "Cumulative infections" => String[],
        "Share of infections (%)" => String[]
    )
    asc === nothing || (df[!, "Relative ascertainment"] = String[])
    for p in 1:np
        row = Any[
            patch_labels[p], _median_ci(R_T[p]; digits),
            _median_ci(C_T[p]; digits = 0), _median_ci(share[p]; digits = 1),
        ]
        asc === nothing || push!(row, _median_ci(asc[p]; digits))
        push!(df, row)
    end
    return df
end

## Equal-tailed 30%, 60% and 90% intervals as one phrase, `30% a–b, 60% c–d,
## 90% e–f`. `digits = 0` writes whole numbers.
function _interval_text(s; digits::Integer = 2, unit::AbstractString = "")
    fmt(x) = digits <= 0 ? string(round(Int, x)) : string(round(x; digits))
    return join(
        (
            string(
                lvl, "% ",
                fmt(getproperty(s, Symbol("lo", lvl))), "–",
                fmt(getproperty(s, Symbol("hi", lvl))), unit
            )
                for lvl in (30, 60, 90)
        ), ", "
    )
end

"""
Markdown bullets summarising each province of the patch model, one bullet per
province. Each gives the cumulative infections and the reproduction number at
the cut-off, and, when the chain carries them, the case-fatality ratio
(`CFR_patch`, as a percentage) and the relative case ascertainment
(`province_ascertainment`). Every quantity is written as its equal-tailed
30%, 60% and 90% credible intervals.

The reproduction number and the relative ascertainment are identified only
as a product by the case composition, so the two are given together.
"""
function patch_headline(
        chn, n_patches::Integer = length(PROVINCE_NAMES);
        patch_labels::AbstractVector = PROVINCE_LABELS
    )
    required = [:C_T_patch, :R_T_patch]
    absent = filter(p -> !_has_key(chn, p), required)
    isempty(absent) || error(
        "chain is missing the per-patch deterministics $(absent); it was " *
            "not sampled from `bvd_joint`."
    )
    np = min(n_patches, length(patch_labels))
    per_patch(sym) = _has_key(chn, sym) ? _per_patch(chn, sym, np) : nothing
    C_T = per_patch(:C_T_patch)
    R_T = per_patch(:R_T_patch)
    cfr = per_patch(:CFR_patch)
    asc = per_patch(:province_ascertainment)
    bullets = map(1:np) do p
        lines = [
            "- **$(patch_labels[p]):** " *
                _interval_text(posterior_summary(C_T[p]); digits = 0) *
                " infections to date.",
            "  The reproduction number at the cut-off is " *
                _interval_text(posterior_summary(R_T[p])) * ".",
        ]
        cfr === nothing || push!(
            lines,
            "  The case-fatality ratio is " *
                _interval_text(
                posterior_summary(100 .* cfr[p]); digits = 1, unit = "%"
            ) * "."
        )
        asc === nothing || push!(
            lines,
            "  Case ascertainment relative to the national average is " *
                _interval_text(posterior_summary(asc[p])) * "."
        )
        join(lines, "\n")
    end
    return join(bullets, "\n") * "\n"
end

"""
Per-patch outbreak summary for the patch model: one row per province, with
the cut-off cumulative infections `C_T`, the cut-off reproduction number
`R_T`, the daily infections at the cut-off, and the log-Rt deviation `δ`
from the common national trend. Each is reported as the same 90/60/30%
credible intervals [`summary_table`](@ref) uses.

Pass `patch` to restrict the table to a single province, by index or by
label. The `Patch` column is then dropped, since it would repeat one value.

The deviations are sum-to-zero contrasts around the national trend (see
[`patch_rt_model`](@ref)), so `δ` is read relative to the national average
across provinces, not relative to any one patch: a negative `δ` means that
province transmits below the national trend. Every patch, including the
primary, carries its own deviation, and they sum to zero in every draw.
For the log-Rt of a province relative to Ituri specifically, read the
chain's `log_rt_contrast` instead.

Expects a chain from [`bvd_joint`](@ref), which stores the per-patch
quantities as vector deterministics (`C_T_patch`, `R_T_patch`,
`infections_T_patch`, `delta_patch`), one entry per patch.
"""
function patch_summary_table(
        chn, n_patches::Integer = length(PROVINCE_NAMES);
        digits::Integer = 2,
        patch::Union{Nothing, Integer, AbstractString} = nothing,
        patch_labels::AbstractVector = PROVINCE_LABELS
    )
    required = [:C_T_patch, :R_T_patch, :infections_T_patch, :delta_patch]
    absent = filter(p -> !_has_key(chn, p), required)
    isempty(absent) || error(
        "chain is missing the per-patch deterministics $(absent); it was " *
            "not sampled from `bvd_joint`."
    )
    np = min(n_patches, length(patch_labels))
    ## Which patches to report. A label is matched against `patch_labels`, so
    ## the caller names the province rather than tracking its index.
    selected = if patch === nothing
        1:np
    elseif patch isa Integer
        1 <= patch <= np || error(
            "patch = $patch is out of range; the chain has $np patches."
        )
        patch:patch
    else
        i = findfirst(==(patch), patch_labels[1:np])
        i === nothing && error(
            "patch = \"$patch\" is not one of $(patch_labels[1:np])."
        )
        i:i
    end
    per_patch(sym) = _per_patch(chn, sym, np)
    C_T = per_patch(:C_T_patch)
    R_T = per_patch(:R_T_patch)
    inf_T = per_patch(:infections_T_patch)
    δ = per_patch(:delta_patch)
    ## Absent on a chain fitted without the per-province compositions.
    asc = _has_key(chn, :province_ascertainment) ?
        per_patch(:province_ascertainment) : nothing
    ## The deviation-walk scale is per patch, so it belongs here rather than
    ## with the scalar hyperparameters. Near zero means that province's Rt
    ## tracks the national trend.
    drift = _has_key(chn, :region_drift_sd) ?
        per_patch(:region_drift_sd) : nothing
    contrast = _has_key(chn, :log_rt_contrast) ?
        per_patch(:log_rt_contrast) : nothing
    df = DataFrame(
        patch = String[],
        quantity = String[],
        lower_90 = Float64[], lower_60 = Float64[], lower_30 = Float64[],
        upper_30 = Float64[], upper_60 = Float64[], upper_90 = Float64[]
    )
    for p in selected
        rows = Any[
            ("Cumulative infections", C_T[p], 0),
            ("Reproduction number", R_T[p], digits),
            ("Daily infections at cut-off", inf_T[p], 0),
            ("log-Rt deviation from trend", δ[p], digits),
        ]
        contrast === nothing ||
            push!(rows, ("log-Rt vs primary patch", contrast[p], digits))
        drift === nothing ||
            push!(rows, ("Rt deviation drift", drift[p], 3))
        asc === nothing ||
            push!(rows, ("Relative case ascertainment", asc[p], digits))
        for (label, draws, dg) in rows
            s = posterior_summary(draws)
            push!(
                df,
                (
                    patch_labels[p], label,
                    round(s.lo90; digits = dg), round(s.lo60; digits = dg),
                    round(s.lo30; digits = dg), round(s.hi30; digits = dg),
                    round(s.hi60; digits = dg), round(s.hi90; digits = dg),
                )
            )
        end
    end
    ## A single-province table would repeat one patch name down every row, so
    ## drop the column. The province belongs in the surrounding heading.
    return patch === nothing ? _prettify(rename(df, :patch => "Patch")) :
        _prettify(select(df, Not(:patch)))
end

## Per-patch draws of the final column of a `(n_patches × n_vintages)` matrix
## deterministic (`province_shares`, `province_death_shares`): the modelled
## split at the most recent vintage.
function _per_patch_last_share(chn, sym::Symbol, np::Integer)
    ms = [collect(v) for v in vec(collect(chn[sym]))]
    nv = size(first(ms), 2)
    return [[m[p, nv] for m in ms] for p in 1:np]
end

"""
Per-province case-fatality ratios, set against the national ones the
[`confirmed_cfr_table`](@ref) reports. `res` is a
[`delay_corrected_confirmed_cfr`](@ref) result and `chn` the patch chain.

Three quantities per province, each a median with a 90% credible interval
except the observed ratio, which is a count:

- the naive observed confirmed ratio, that province's reported confirmed
  deaths over its reported confirmed cases;
- the delay-corrected confirmed ratio, the national corrected ratio scaled
  by the province's relative lethality and death confirmation over its
  relative case ascertainment, which is what varies by province once the
  delays are corrected for;
- the structural, infection-based ratio, read from `CFR_patch`, the national
  ratio times that province's sum-to-zero lethality contrast.

The death composition identifies only the product of the lethality and the
death-confirmation contrasts, so their split is set by their priors. The
lethality prior is the looser of the two, so a provincial excess of deaths
over cases is read first as lethality. A chain carrying no `CFR_patch` falls
back to the national ratio in the structural column.

`province_cases` and `province_deaths` are the observed per-province
confirmed case and death totals over the fitted window, in the order of
`patch_labels`.
"""
function province_cfr_table(
        chn, res;
        province_cases::AbstractVector, province_deaths::AbstractVector,
        n_patches::Integer = length(PROVINCE_NAMES),
        patch_labels::AbstractVector = PROVINCE_LABELS,
        digits::Integer = 1
    )
    np = min(n_patches, length(patch_labels))
    _has_key(chn, :province_ascertainment) || error(
        "chain carries no `province_ascertainment`; it was not sampled " *
            "from `bvd_joint` with the per-province compositions on."
    )
    case_asc = _per_patch(chn, :province_ascertainment, np)
    death_asc = _has_key(chn, :province_death_ascertainment) ?
        _per_patch(chn, :province_death_ascertainment, np) :
        [ones(length(case_asc[1])) for _ in 1:np]
    ## Per-province lethality contrast, and the per-province structural ratio
    ## it implies. Both absent on a chain that reports the national ratio in
    ## every row.
    sev = _has_key(chn, :province_cfr_relative) ?
        _per_patch(chn, :province_cfr_relative, np) :
        [ones(length(case_asc[1])) for _ in 1:np]
    cfr_patch = _has_key(chn, :CFR_patch) ?
        _per_patch(chn, :CFR_patch, np) : nothing
    ## Mask rather than filter, so the corrected draws and the per-province
    ## scaling below stay aligned draw for draw.
    mask = isfinite.(res.corrected)
    corrected = res.corrected[mask]
    structural = filter(isfinite, res.structural)
    pct(x) = round(100 * x; digits)
    cell(v) = string(
        pct(quantile(v, 0.5)), "% (",
        pct(quantile(v, 0.05)), "–", pct(quantile(v, 0.95)), "%)"
    )
    df = DataFrame(
        "Province" => String[],
        "Naive observed confirmed ratio" => String[],
        "Delay-corrected confirmed CFR" => String[],
        "Structural (infection-based) CFR" => String[]
    )
    for p in 1:np
        naive = province_cases[p] > 0 ?
            string(pct(province_deaths[p] / province_cases[p]), "%") : "—"
        ## The delay correction is national, so the province enters only
        ## through the ratio of its relative death confirmation to its
        ## relative case ascertainment. Both are sum-to-zero on the log scale,
        ## so the corrected ratios sit around the national one.
        scale = ((sev[p] .* death_asc[p]) ./ case_asc[p])[mask]
        struc_p = cfr_patch === nothing ? structural :
            filter(isfinite, cfr_patch[p])
        push!(
            df, (
                patch_labels[p], naive,
                cell(corrected .* scale), cell(struc_p),
            )
        )
    end
    return df
end

"""
    province_recent_counts(province_history, province_names, n_patches;
                           window = 7)

New counts in each patch over the most recent `window` days of the
per-province cumulative histories, read off the same increments
[`province_increment_matrix`](@ref) builds, so a downward revision counts as
no new cases.

The span runs from the latest vintage at least `window` days before the last
one to the last vintage, so it can be longer than `window` when the tables
skip days. Returns `(; start_day, last_day, counts)` with the span's grid day
indices and one count per patch, or `nothing` when there is no history or it
is shorter than the window.
"""
function province_recent_counts(
        province_history,
        province_names::AbstractVector, n_patches::Integer;
        window::Integer = 7
    )
    inc = province_increment_matrix(province_history, province_names, n_patches)
    isempty(inc.days) && return nothing
    last_day = inc.days[end]
    ref = findlast(<=(last_day - window), inc.days)
    ref === nothing && return nothing
    counts = vec(sum(inc.increments[:, (ref + 1):end]; dims = 2))
    return (; start_day = inc.days[ref], last_day, counts)
end

## Streams the per-province split covers, as the forecast column and the
## stream label. The label is the one `forecast_archive` gives the national
## stream, so a province row names the stream it is a share of, and the
## tables and figure build their headings from it.
const _PROVINCE_FORECAST_STREAMS = (
    (:confirmed_new, "confirmed cases"),
    (:confirmed_deaths_new, "confirmed deaths"),
)

## Per-province forecast draws from one [`forecast_reported`](@ref) result:
## the national draw times that province's modelled share at the most recent
## spatial vintage, multiplied draw by draw so the split carries the
## correlation between the two factors. Returns one
## `(stream_label, province, draws)` entry per province and per stream the
## forecast carries, provinces outer. Shared by the province forecast table,
## figure and release archive, so all three read one split.
function _province_forecast_draws(
        chn, fc, np::Integer,
        patch_labels::AbstractVector
    )
    _has_key(chn, :province_shares) || error(
        "chain carries no `province_shares`; it was not sampled from " *
            "`bvd_joint` with the per-province compositions on."
    )
    case_share = _per_patch_last_share(chn, :province_shares, np)
    death_share = _has_key(chn, :province_death_shares) ?
        _per_patch_last_share(chn, :province_death_shares, np) :
        case_share
    cols = propertynames(fc)
    out = Tuple{String, String, Vector{Float64}}[]
    for p in 1:np, (col, label) in _PROVINCE_FORECAST_STREAMS

        col in cols || continue
        v = float.(fc[!, col])
        share = col === :confirmed_new ? case_share : death_share
        push!(out, (label, patch_labels[p], v .* share[p][1:length(v)]))
    end
    return out
end

"""
Per-province split of the one-week-ahead national forecast `fc` from
[`forecast_reported`](@ref): the new confirmed cases and confirmed deaths
expected in each province over the week to `T + 7`, as the same 90/60/30%
credible intervals [`forecast_table`](@ref) reports nationally. The same
content is drawn by [`plot_province_forecast`](@ref) and archived for
scoring by [`province_forecast_archive`](@ref).

Each province's count is the national draw times that province's modelled
share at the most recent spatial vintage. The split is held at its current
value over the week rather than projected forward. The provincial
compositions are fitted only where the spatial tables report, so a province
whose share is moving is not tracked past the last vintage.
"""
function province_forecast_table(
        chn, fc;
        n_patches::Integer = length(PROVINCE_NAMES),
        patch_labels::AbstractVector = PROVINCE_LABELS,
        digits::Integer = 0
    )
    np = min(n_patches, length(patch_labels))
    rows = NamedTuple[]
    for (label, province, draws) in _province_forecast_draws(
            chn, fc, np, patch_labels
        )
        s = posterior_summary(draws)
        push!(
            rows,
            (
                province = province, quantity = "New $(label) by T+7",
                lower_90 = round(s.lo90; digits),
                lower_60 = round(s.lo60; digits),
                lower_30 = round(s.lo30; digits),
                upper_30 = round(s.hi30; digits),
                upper_60 = round(s.hi60; digits),
                upper_90 = round(s.hi90; digits),
            )
        )
    end
    return _prettify(DataFrame(rows))
end

"""
Per-province forecast against what was observed. `fc` is a
[`forecast_reported`](@ref) result from a frozen patch fit, `chn` that same
chain, and `observed` and `baseline` the per-province cumulative counts at
the target date and at the forecast origin, so the truth is their difference.

Each province's forecast is the national draw times that province's modelled
share at the frozen fit's most recent spatial vintage. The two factors are
multiplied draw by draw, so the interval carries their correlation. The share
itself is held over the horizon, which the width does not express. A province
whose share is moving is scored as though it were not.

Reports the 90% predictive interval, the observed count, and whether the
observation fell inside the interval, one row per province and stream. No
central estimate is reported.
"""
function province_forecast_vs_truth(
        chn, fc;
        observed::AbstractVector, baseline::AbstractVector,
        death_observed::Union{Nothing, AbstractVector} = nothing,
        death_baseline::Union{Nothing, AbstractVector} = nothing,
        n_patches::Integer = length(PROVINCE_NAMES),
        patch_labels::AbstractVector = PROVINCE_LABELS,
        digits::Integer = 0
    )
    np = min(n_patches, length(patch_labels))
    _has_key(chn, :province_shares) || error(
        "chain carries no `province_shares`; the frozen fit was not run " *
            "with the per-province compositions on."
    )
    case_share = _per_patch_last_share(chn, :province_shares, np)
    death_share = _has_key(chn, :province_death_shares) ?
        _per_patch_last_share(chn, :province_death_shares, np) :
        case_share
    cols = propertynames(fc)
    rows = NamedTuple[]
    function add!(stream, p, draws, truth)
        s = posterior_summary(draws)
        return push!(
            rows,
            (
                province = patch_labels[p], stream = stream,
                lower_90 = round(s.lo90; digits),
                upper_90 = round(s.hi90; digits),
                observed = truth,
                within_90 = s.lo90 <= truth <= s.hi90,
            )
        )
    end
    for p in 1:np
        if :confirmed_new in cols
            v = fc[!, :confirmed_new]
            add!(
                "Confirmed cases", p, v .* case_share[p][1:length(v)],
                observed[p] - baseline[p]
            )
        end
        if :confirmed_deaths_new in cols && death_observed !== nothing
            v = fc[!, :confirmed_deaths_new]
            add!(
                "Confirmed deaths", p, v .* death_share[p][1:length(v)],
                death_observed[p] - death_baseline[p]
            )
        end
    end
    return _prettify(DataFrame(rows))
end
