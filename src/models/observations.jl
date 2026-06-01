# Observation submodels: each ties one data stream to the shared daily
# onset incidence from the generating infection process. A stream
# convolves the onsets with its own sampled onset-to-event delay, scales
# by the relevant ascertainment / CFR / positivity / travel factor, and
# reads the cumulative expected count off the daily series at each vintage
# day, fitting the between-vintage increments (and the cut-off total).
# Convolutions replace the integral model's continuous onset-to-event
# integrals while preserving the v1.3.0 per-vintage time-series
# likelihoods.

"""
NaN / Inf-safe `NegativeBinomial` constructor parameterised by mean `μ`
and dispersion `k`, with clamping on the success probability so extreme
NUTS proposals during warmup do not trip the distribution domain check.
Shared by the count-stream observation submodels.
"""
function safe_nbinomial(k, μ)
    p_raw = k / (k + max(μ, eps(typeof(μ))))
    p = isfinite(p_raw) ?
        clamp(p_raw, eps(typeof(k)), one(k) - eps(typeof(k))) :
        eps(typeof(k))
    return NegativeBinomial(k, p)
end

"""
Cumulative of a daily series `daily` read off at the vintage day indices
`days` (1-based into the grid). Each entry is the running sum up to and
including that day, giving the modelled cumulative count at each vintage.
Days beyond the series length are clamped to the final cumulative. Pure
and AD-transparent; the output element type follows `daily`.
"""
function cumulative_at(daily::AbstractVector, days::AbstractVector{<:Integer})
    cum = cumsum(daily)
    n = length(cum)
    return [cum[clamp(d, 1, n)] for d in days]
end

"""
Resolve a stream's per-vintage observation into the vintage day indices
and the observed between-vintage increment vector to score, given the
dated cumulative `history` `(; days, counts)`, the cut-off total `total`
and the grid length `n` (day `n` is the cut-off). When the history is
non-empty it already ends at the cut-off (the last vintage count equals
`total`), so the increments are differenced from the history alone and the
separate cut-off total is not scored again. When the history is empty but
a cut-off `total` is supplied (e.g. the tests-analysed stream, whose dated
vintage history is absent), the cut-off becomes a single vintage point at
day `n` so the total is still scored as one increment. When both are empty
or `total` is `missing`, the increments are `missing` so the caller is a
predictive generator. Returns `(; days, obs_increments)` with
`obs_increments` either an `Int` vector or `missing`.
"""
function vintage_obs(history, total::Union{Missing, Integer}, n::Integer)
    if !isempty(history.counts)
        days = history.days
        obs_increments = diff(vcat(zero(eltype(history.counts)),
            collect(history.counts)))
        return (; days, obs_increments)
    elseif !ismissing(total)
        return (; days = [Int(n)], obs_increments = [Int(total)])
    else
        return (; days = Int[], obs_increments = missing)
    end
end

"""
Per-vintage cumulative-history likelihood for one stream, expressed as a
proper vector likelihood. Given the modelled cumulative `expected_cum` at
each vintage day, scores the between-vintage increments with
NegativeBinomials sharing the dispersion `k` (one per vintage). The
modelled increment for vintage `i` is `expected_cum[i] - expected_cum[i-1]`
(the first is `expected_cum[1]`).

The observed increments are passed as `obs_increments` and scored with a
loop of scalar `~` so the stream is real observed data that `predict`
replicates. The increment variable is named `increments`, so under a
prefixed submodel attachment the predict keys are `<prefix>.increments[i]`.
When `obs_increments` is `missing` the increments are sampled, making the
submodel a predictive generator. When it is empty (zero vintages) the
likelihood is a no-op. Returns the modelled and (when present) observed
increments for reuse.
"""
@model function vintage_increments_model(expected_cum::AbstractVector,
        increments::Union{Missing, AbstractVector{<:Integer}},
        k::Real)
    modelled = diff(vcat(zero(eltype(expected_cum)), expected_cum))
    n = length(modelled)
    ## `increments` is the model argument scored on the LHS of `~`, so a
    ## supplied vector is observed data DynamicPPL conditions on and a
    ## `missing` argument is sampled (the predictive-generator path). The
    ## indexed `increments[i]` keeps the predict keys
    ## (`<prefix>.increments[i]`) replicable.
    if ismissing(increments)
        increments = Vector{Union{Missing, Int}}(missing, n)
    end
    for i in 1:n
        increments[i] ~ safe_nbinomial(k, safe_rate(modelled[i]))
    end
    return (; modelled, increments)
end

"""
DRC suspected-deaths likelihood, per-vintage time series. Convolves the
daily onsets with the sampled onset-to-death delay, scales by the CFR, and
reads the modelled cumulative deaths at each vintage day off the daily
series, fitting the between-vintage increments as observed `~` data with a
NegativeBinomial sharing the surveillance dispersion `k`
([`surveillance_dispersion_model`](@ref)). The death history ends at the
cut-off, so the cut-off total is the final increment and is not scored
separately. Samples the onset-to-death delay
and the CFR via injected submodels. The onset-to-death prior is centred on
the Bayesian BDBV line-list reanalysis (mean 11.2 d, SD 5.4 d; the
`bdbv-linelist-analysis` submodule), the same source the integral model
used. Returns the cut-off expected count, the daily death series, the
onset-to-death PMF and the CFR for reuse by [`exports_deaths_model`](@ref).
"""
@model function deaths_model(
        deaths_history,
        total_deaths::Union{Missing, Integer},
        onsets::AbstractVector, k::Real;
        cfr = cfr_model(),
        onset_to_death = censored_delay_model(60;
            mean_prior = truncated(Normal(11.2, 2.0); lower = 1),
            sd_prior = truncated(Normal(5.4, 1.5); lower = 1)))
    cfr_state ~ to_submodel(cfr)
    od_state ~ to_submodel(onset_to_death)
    CFR = cfr_state.CFR
    deaths_daily = CFR .* convolve_delay(onsets, od_state.pmf)

    n = length(deaths_daily)
    vobs = vintage_obs(deaths_history, total_deaths, n)
    expected_cum = cumulative_at(deaths_daily, vobs.days)
    death_increments ~ to_submodel(
        vintage_increments_model(expected_cum, vobs.obs_increments, k))

    raw_total = sum(deaths_daily)
    expected_deaths_T := safe_rate(raw_total)

    return (; CFR, od_pmf = od_state.pmf, deaths_daily, expected_deaths_T)
end

"""
DRC suspected (reported) cases likelihood, per-vintage time series. The
expected daily suspected cases are a BVD-driven onset-to-report
convolution scaled by the DRC ascertainment fraction `p_drc`, plus an
additive non-BVD background rate `λ_bg` per day (so a suspected case need
not be a true BVD infection). Reads the modelled cumulative suspected
cases at each vintage day off the daily series and fits the increments
with a NegativeBinomial sharing `k`. The background and testing fraction
are sampled by an injected [`test_positivity_model`](@ref), and the
onset-to-report delay is injected, defaulting to a weakly-informative
prior on the onset-to-notification delay (mean 4.5 d, SD 3.6 d),
consistent with Ebola surveillance reporting delays.

Returns the onset-to-report PMF and the BVD onset-to-report daily series
(unit ascertainment, no background) so [`confirmed_cases_model`](@ref) can
reuse the same report kernel, the sampled background rate and testing
fraction, and the implied per-suspected positivity (the BVD share of the
expected suspected total) as a derived quantity for comparison with the
sitrep.
"""
@model function reported_cases_model(
        reported_history,
        reported_cases::Union{Missing, Integer},
        onsets::AbstractVector, k::Real, p_drc::Real;
        positivity = test_positivity_model(),
        onset_to_report = censored_delay_model(30;
            mean_prior = truncated(Normal(4.5, 1.5); lower = 1),
            sd_prior = truncated(Normal(3.6, 1.2); lower = 1)))
    pos_state ~ to_submodel(positivity)
    report_state ~ to_submodel(onset_to_report)
    λ_bg = pos_state.λ_bg
    τ_test = pos_state.τ_test
    report_pmf = report_state.pmf

    ## Unit-ascertainment BVD onset-to-report daily series, reused by the
    ## confirmed stream. Suspected daily cases add the p_drc-scaled BVD
    ## signal and the constant non-BVD background.
    bvd_reports_daily = convolve_delay(onsets, report_pmf)
    reports_daily = p_drc .* bvd_reports_daily .+ λ_bg

    n = length(reports_daily)
    vobs = vintage_obs(reported_history, reported_cases, n)
    expected_cum = cumulative_at(reports_daily, vobs.days)
    reported_increments ~ to_submodel(
        vintage_increments_model(expected_cum, vobs.obs_increments, k))

    raw_total = sum(reports_daily)
    expected_reports := safe_rate(raw_total)

    ## Implied per-suspected positivity at the cut-off: BVD share of the
    ## expected suspected total.
    bvd_total = p_drc * sum(bvd_reports_daily)
    positivity := safe_rate(bvd_total) / expected_reports

    return (; p_drc, λ_bg, τ_test, report_pmf, bvd_reports_daily,
        reports_daily, expected_reports, positivity)
end

"""
Laboratory pipeline likelihood, per-vintage time series. Two coupled
streams driven by the shared renewal onsets:

- Confirmed cases. The BVD onset-to-report daily series (the unit
  ascertainment `bvd_reports_daily` shared from
  [`reported_cases_model`](@ref)) is convolved with the sampled
  report-to-confirmation (lab-turnaround) delay, then scaled by the DRC
  ascertainment `p_drc`, the PCR sensitivity `s_test` and the testing
  fraction `τ_test`. The modelled cumulative confirmed cases are read off
  at each `confirmed_history` vintage day and the increments fitted as
  observed `~` data with a NegativeBinomial sharing `k`. The confirmed
  history ends at the cut-off, so the cut-off total is the final increment
  and is not scored separately.
- Tests analysed. The tested daily volume is `τ_test` times the sum of the
  `p_drc`-scaled BVD lab series and the non-BVD background `λ_bg` carried
  through the same lab delay. Its cumulative is read off at the
  `lab_history` vintage days (so a lab stream that lags the case cut-off is
  right-truncated by reading the running cumulative at its own days) and
  the increments fitted with a NegativeBinomial sharing `k`. When the
  dated `lab_history` is absent, the cut-off `tests_analysed` total is
  scored as a single vintage point at day `n`. Confirmed
  counts are not re-observed conditional on the tested volume, so the two
  lab streams are not double-counted.

Samples the lab-turnaround delay and the PCR sensitivity via injected
submodels; the testing fraction `τ_test` and background rate `λ_bg` come
from [`reported_cases_model`](@ref) so the suspected and laboratory
streams share them. Exposes the per-test positivity
`s_test · BVD_tested / tested` and the expected confirmed and tested
totals at the cut-off as derived quantities.
"""
@model function confirmed_cases_model(
        confirmed_history,
        confirmed_cases::Union{Missing, Integer},
        onsets::AbstractVector, k::Real, p_drc::Real,
        λ_bg::Real, τ_test::Real, bvd_reports_daily::AbstractVector;
        lab_history = (; days = Int[], counts = Int[]),
        tests_analysed::Union{Missing, Integer} = missing,
        sensitivity = test_sensitivity_model(),
        report_to_confirm = lab_delay_model())
    sens_state ~ to_submodel(sensitivity)
    lab_state ~ to_submodel(report_to_confirm)
    s_test = sens_state.s_test
    lab_pmf = lab_state.pmf

    ## BVD report signal carried through the lab-turnaround delay (unit
    ## ascertainment), and the non-BVD background carried likewise.
    bvd_lab_daily = convolve_delay(bvd_reports_daily, lab_pmf)
    bg_lab_daily = convolve_delay(fill(λ_bg, length(onsets)), lab_pmf)

    confirmed_daily = (p_drc * s_test * τ_test) .* bvd_lab_daily
    tested_daily = τ_test .* (p_drc .* bvd_lab_daily .+ bg_lab_daily)

    n = length(confirmed_daily)
    cvobs = vintage_obs(confirmed_history, confirmed_cases, n)
    confirmed_cum = cumulative_at(confirmed_daily, cvobs.days)
    confirmed_increments ~ to_submodel(
        vintage_increments_model(confirmed_cum, cvobs.obs_increments, k))

    tvobs = vintage_obs(lab_history, tests_analysed, n)
    tested_cum = cumulative_at(tested_daily, tvobs.days)
    tested_increments ~ to_submodel(
        vintage_increments_model(tested_cum, tvobs.obs_increments, k))

    expected_confirmed := safe_rate(sum(confirmed_daily))
    expected_tested := safe_rate(sum(tested_daily))

    ## Per-test positivity: PCR sensitivity times the BVD share of the
    ## tested pool, for comparison with the sitrep positivity rate.
    bvd_tested = p_drc * sum(bvd_lab_daily)
    p_positive := s_test * safe_rate(bvd_tested) /
                  safe_rate(p_drc * sum(bvd_lab_daily) + sum(bg_lab_daily))

    return (; s_test, τ_test, λ_bg, lab_pmf, confirmed_daily, tested_daily,
        expected_confirmed, expected_tested, p_positive)
end

"""
Uganda exports likelihood (geographic spread). Builds the export onset
incidence `p_uganda · q · onsets` (with `q = daily_travellers /
source_population` the per-capita travel rate) and convolves it with the
sampled onset-to-detection delay, replacing the integral model's
detection-window survival term with a convolved set of delays. Sums to the
expected detected exports by the cut-off, fitted with Poisson (Uganda's
stream is small). Samples the traveller volume and the onset-to-detection
delay via injected submodels. The onset-to-detection prior is centred on
the Ebola onset-to-hospitalisation delay (mean 5.0 d, SD 4.7 d; WHO Ebola
Response Team 2014, NEJM), used here as the delay from symptom onset to
detection at a point of entry abroad. Returns the expected count, the
detection-timed series and the export onsets for reuse by
[`exports_deaths_model`](@ref).
"""
@model function exports_model(
        exported_cases::Union{Missing, Integer},
        onsets::AbstractVector, p_uganda::Real;
        source_population::Real = ITURI_POPULATION,
        traveller = traveller_volume_model(),
        onset_to_detection = censored_delay_model(30;
            mean_prior = truncated(Normal(5.0, 2.0); lower = 1),
            sd_prior = truncated(Normal(4.7, 1.5); lower = 1)))
    travel_state ~ to_submodel(traveller)
    daily_travellers = travel_state.daily_travellers
    q = daily_travellers / source_population

    export_onsets = p_uganda .* q .* onsets
    detect_state ~ to_submodel(onset_to_detection)
    detect_daily = convolve_delay(export_onsets, detect_state.pmf)

    raw_exports = sum(detect_daily)
    expected_exports_T := safe_rate(raw_exports)
    exported_cases ~ Poisson(expected_exports_T)

    return (; p_uganda, daily_travellers, export_onsets,
        expected_exports = expected_exports_T)
end

"""
Deaths-among-detected-exports likelihood. Convolves the export onsets
(timed from onset, the same staging as detection) with the onset-to-death
PMF shared from [`deaths_model`](@ref), scales by the CFR, and sums to the
expected cumulative export deaths by the cut-off, fitted with Poisson.
"""
@model function exports_deaths_model(
        exports_deaths::Union{Missing, Integer},
        export_onsets::AbstractVector, CFR::Real,
        od_pmf::AbstractVector)
    series = CFR .* convolve_delay(export_onsets, od_pmf)

    raw = sum(series)
    expected_exports_deaths_T := safe_rate(raw)
    exports_deaths ~ Poisson(expected_exports_deaths_T)

    return (; expected_exports_deaths_T, series)
end
