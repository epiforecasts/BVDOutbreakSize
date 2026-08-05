# Delay-corrected confirmed case-fatality ratio. The naive confirmed CFR,
# cumulative confirmed deaths over cumulative confirmed cases, is biased low
# in real time because recently-confirmed cases have not yet had time to die.
# We debias it by shrinking the denominator to the confirmed cases expected
# to have had a fatal outcome confirmed by the cut-off, using the residual
# delay between a confirmed case and its confirmed death (the onset-to-death-
# confirmation lag minus the onset-to-confirmation lag). The correction
# (Nishiura et al. 2009) is computed per posterior draw on the model's own
# confirmed-case trajectory and sampled delays, so it propagates the joint
# uncertainty. It sits alongside the structural (infection-based) CFR the
# model already reports, which is harder to identify because case and death
# ascertainment differ.

# Probability a case confirmed `δ` days before the cut-off has, if fatal,
# had its death confirmed by the cut-off: `P(X_d − X_c ≤ δ)`. The onset-to-
# death-confirmation lag `X_d ~ Kd` and the onset-to-confirmation lag
# `X_c ~ Kc` are assumed independent, computed as `Σ_xc Kc[xc] · Fd(xc + δ)`
# from the onset-to-death-confirmation CDF `Fd` (cumulative `Kd`). `δ` may be
# negative (a case confirmed after the cut-off horizon contributes no resolved
# outcome), in which case `Fd` is read at a negative lag and returns zero.
function _residual_outcome_cdf(Kc::AbstractVector, Fd::AbstractVector, δ::Integer)
    Ld = length(Fd)
    acc = 0.0
    @inbounds for j in eachindex(Kc)
        m = (j - 1) + δ                  # death-confirmation lag threshold
        f = m < 0 ? 0.0 : (m >= Ld ? Fd[Ld] : Fd[m + 1])
        acc += Kc[j] * f
    end
    return acc
end

"""
Delay-corrected confirmed case-fatality ratio for one posterior draw.

`c_daily` is the modelled daily confirmed-case incidence over the day grid
(day `length(c_daily)` is the cut-off), `Kc` the onset-to-confirmation delay
PMF (lag 0 at index 1), `Kd` the onset-to-death-confirmation delay PMF, and
`deaths_total` the cumulative confirmed deaths at the cut-off. The corrected
denominator is `Σ_t c_daily[t] · P(outcome resolved by the cut-off | confirmed
at t)`, smaller than the raw cumulative confirmed cases, so the corrected
ratio `deaths_total / denominator` lifts the naive ratio toward the eventual
confirmed CFR. Returns `NaN` when no confirmed cases have had time to resolve.
"""
function delay_corrected_cfr(c_daily::AbstractVector, Kc::AbstractVector,
        Kd::AbstractVector, deaths_total::Real)
    T = length(c_daily)
    (T == 0 || isempty(Kc) || isempty(Kd)) && return NaN
    Fd = cumsum(Kd)
    denom = 0.0
    @inbounds for t in 1:T
        denom += float(c_daily[t]) * _residual_outcome_cdf(Kc, Fd, T - t)
    end
    return denom > 0 ? float(deaths_total) / denom : NaN
end

# Per-draw vectors stored by a vector deterministic (`cumulative_confirmed`,
# the delay PMFs): the chain holds an iter×chain matrix of per-draw vectors.
function _draw_vectors(chn, key::Symbol)
    mat = chn[key]
    return [collect(v) for v in vec(collect(mat))]
end

# Daily series from a cumulative trajectory: the first day carries the
# cut-in level, each later day the increment.
function _to_daily(cum::AbstractVector)
    n = length(cum)
    n == 0 && return Float64[]
    out = Vector{Float64}(undef, n)
    out[1] = float(cum[1])
    @inbounds for t in 2:n
        out[t] = float(cum[t]) - float(cum[t - 1])
    end
    return out
end

"""
Delay-corrected confirmed case-fatality ratio across the posterior, read off
a joint `bvd_joint` chain. For each draw the modelled daily confirmed-case
incidence (rescaled so its total matches the scored expected confirmed cases),
the onset-to-confirmation and onset-to-death-confirmation delay PMFs, and the
cumulative confirmed deaths give the corrected ratio
([`delay_corrected_cfr`](@ref)).

`obs_confirmed` and `obs_confirmed_deaths` are the observed cumulative
laboratory-confirmed cases and confirmed deaths at the cut-off, used for the
naive observed confirmed ratio.

Returns a `NamedTuple` with the per-draw vectors `corrected` (delay-corrected
confirmed CFR), `modelled_naive` (uncorrected modelled confirmed deaths over
confirmed cases) and `structural` (the model's infection/onset-level `CFR`),
and the scalar `naive_observed = obs_confirmed_deaths / obs_confirmed`.
"""
function delay_corrected_confirmed_cfr(chn;
        obs_confirmed::Real, obs_confirmed_deaths::Real)
    cum_conf = _draw_vectors(chn, :cumulative_confirmed)
    Kc = _draw_vectors(chn, :onset_to_confirmation_pmf)
    Kd = _draw_vectors(chn, :onset_to_death_confirmation_pmf)
    deaths_T = _draws(chn, :expected_confirmed_deaths_T)
    cases_T = _draws(chn, :expected_confirmed_T)
    structural = _draws(chn, :CFR)

    nd = length(deaths_T)
    corrected = Vector{Float64}(undef, nd)
    modelled_naive = Vector{Float64}(undef, nd)
    @inbounds for i in 1:nd
        c_daily = _to_daily(cum_conf[i])
        ## Rescale the modelled confirmed-case incidence so its total equals
        ## the scored expected confirmed cases `cases_T`. The daily series is
        ## built on the modelled analysed volume, but the observed-denominator
        ## windows score against the observed analysed counts, so the two can
        ## differ by a level. Rescaling keeps the uncorrected limit equal to
        ## the modelled confirmed CFR `deaths_T / cases_T`.
        s = sum(c_daily)
        if s > 0
            c_daily .*= cases_T[i] / s
        end
        corrected[i] = delay_corrected_cfr(c_daily, Kc[i], Kd[i], deaths_T[i])
        modelled_naive[i] = cases_T[i] > 0 ? deaths_T[i] / cases_T[i] : NaN
    end
    naive_observed = obs_confirmed > 0 ?
                     float(obs_confirmed_deaths) / float(obs_confirmed) : NaN
    return (; corrected, modelled_naive, structural, naive_observed)
end

"""
One-row `DataFrame` summarising the confirmed-CFR comparison from a
[`delay_corrected_confirmed_cfr`](@ref) result `res`: the median and
equal-tailed 90% credible interval of the delay-corrected confirmed CFR and
the structural (infection-based) CFR, the median uncorrected modelled
confirmed ratio, and the naive observed confirmed ratio. Percentages rounded
to `digits` decimal places.
"""
function confirmed_cfr_table(res; digits::Integer = 1)
    pct(x) = round(100 * x; digits)
    med(v) = pct(quantile(filter(isfinite, v), 0.5))
    function ci(v)
        s = posterior_summary(filter(isfinite, v))
        return string(pct(s.lo90), "–", pct(s.hi90), "%")
    end
    return DataFrame(
        quantity = ["Delay-corrected confirmed CFR",
            "Structural (infection-based) CFR",
            "Uncorrected modelled confirmed ratio",
            "Naive observed confirmed ratio"],
        central_estimate = [string(med(res.corrected), "%"),
            string(med(res.structural), "%"),
            string(med(res.modelled_naive), "%"),
            string(pct(res.naive_observed), "%")],
        narrowest_interval = [ci(res.corrected), ci(res.structural),
            ci(res.modelled_naive), "—"]
    ) |> _prettify
end
