# Observation submodels: each ties one data stream to the shared daily
# onset incidence from the generating infection process. A stream
# convolves the onsets with its own sampled onset-to-event delay, scales
# by the relevant ascertainment / CFR / positivity / travel factor, and
# bins the daily series into per-vintage increments (the first bin is the
# cumulative up to the first vintage, the rest inter-vintage sums), fitting
# them against the observed increments. Convolutions replace the integral
# model's continuous onset-to-event integrals while preserving the v1.3.0
# per-vintage time-series likelihoods.

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
Modelled between-vintage increments of a daily series `daily`, summed
directly into the bins delimited by the vintage day indices `days` (1-based
into the grid, ascending). The first increment is the cumulative count up
to the first vintage day (`sum(daily[1:days[1]])`); each later increment is
the inter-vintage sum `sum(daily[days[i-1]+1:days[i]])`, so only the first
bin needs the cumulative. This avoids the cumulative-then-difference round
trip — the daily series is binned once into the quantities the likelihood
scores. Day indices are clamped to the grid. Pure and AD-transparent; the
output element type follows `daily`.
"""
function bin_increments(daily::AbstractVector,
        days::AbstractVector{<:Integer})
    n = length(daily)
    out = Vector{eltype(daily)}(undef, length(days))
    prev = 0
    @inbounds for (i, d) in enumerate(days)
        hi = clamp(Int(d), 0, n)
        lo = clamp(prev, 0, n)
        out[i] = hi > lo ? sum(@view daily[(lo + 1):hi]) : zero(eltype(daily))
        prev = hi
    end
    return out
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
day `n` so the total is still scored as one increment. A history with
`days` but empty `counts` keeps the vintage day grid while leaving the
increments `missing`, the posterior-predictive generator path that
`predict` resamples. When both are empty or `total` is `missing`, the
increments are `missing` over zero days. Returns `(; days, obs_increments)`
with `obs_increments` either an `Int` vector or `missing`.
"""
function vintage_obs(history, total::Union{Missing, Integer}, n::Integer)
    if !isempty(history.counts)
        days = history.days
        obs_increments = diff(vcat(zero(eltype(history.counts)),
            collect(history.counts)))
        return (; days, obs_increments)
    elseif !isempty(history.days)
        return (; days = collect(history.days), obs_increments = missing)
    elseif !ismissing(total)
        return (; days = [Int(n)], obs_increments = [Int(total)])
    else
        return (; days = Int[], obs_increments = missing)
    end
end

"""
Per-vintage increment likelihood for one stream, expressed as a proper
vector likelihood. Given the modelled per-vintage increments `modelled`
(see [`bin_increments`](@ref)), scores them against the observed
increments with NegativeBinomials sharing the dispersion `k` (one per
vintage).

The observed increments are passed as `increments` and scored with a loop
of scalar `~` so the stream is real observed data that `predict`
replicates. The increment variable is named `increments`, so under a
prefixed submodel attachment the predict keys are `<prefix>.increments[i]`.
When `increments` is `missing` the increments are sampled, making the
submodel a predictive generator. When it is empty (zero vintages) the
likelihood is a no-op. Returns the modelled and (when present) observed
increments for reuse.
"""
@model function vintage_increments_model(modelled::AbstractVector,
        increments::Union{Missing, AbstractVector{<:Integer}},
        k::Real)
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
    modelled_increments = bin_increments(deaths_daily, vobs.days)
    death_increments ~ to_submodel(
        vintage_increments_model(modelled_increments, vobs.obs_increments, k))

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
    modelled_increments = bin_increments(reports_daily, vobs.days)
    reported_increments ~ to_submodel(
        vintage_increments_model(modelled_increments, vobs.obs_increments, k))

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
Per-window confirmed-positives likelihood, expressed as a vector
likelihood scored against the observed analysed denominators. Given the
per-window analysed counts `analysed` and positivities `p_pos`, scores the
observed `positives` with one `Binomial(analysed[i], p_pos[i])` per window.
`positives` is the model argument on the LHS of `~`, so a supplied vector
is observed data DynamicPPL conditions on (mirroring
[`vintage_increments_model`](@ref)) and a `missing` argument is sampled,
making the submodel a predictive generator. Under a prefixed submodel
attachment the predict keys are `<prefix>.positives[i]`. Returns the
observed (or sampled) positives.
"""
@model function confirmed_positives_model(
        positives::Union{Missing, AbstractVector{<:Integer}},
        analysed::AbstractVector{<:Integer}, p_pos::AbstractVector)
    nv = length(analysed)
    if ismissing(positives)
        positives = Vector{Union{Missing, Int}}(missing, nv)
    end
    for i in 1:nv
        positives[i] ~ Binomial(analysed[i], p_pos[i])
    end
    return (; positives)
end

"""
Align the confirmed-case counts onto the laboratory windows, splitting
them into two non-overlapping groups so all the confirmed data is used:

- *Observed* windows, where an analysed-specimen increment is available
  (the laboratory series only differences cleanly from its second vintage
  onward; the first cumulative value is the baseline). Each window's
  analysed increment is the `Binomial` denominator and the matching
  confirmed increment the positives. A zero analysed increment (the
  24-25 May analysis stall) or a window whose positives exceed its
  denominator is merged forward, so every observed window has
  `obs_analysed > 0` and `0 ≤ obs_positives ≤ obs_analysed`.
- *Early* windows, the confirmed vintages up to and including the first
  laboratory date (18-23 May), which have no per-vintage analysed
  denominator. Their confirmed increments are returned so the model can
  score them against the modelled laboratory volume, with the per-window
  positivity partially pooled with the observed windows.

The two groups partition the confirmed counts at the first laboratory
date, so no confirmed case is counted twice. Returns
`(; obs_days, obs_positives, obs_analysed, early_days, early_increments)`
of grid day-indices and per-window counts; the observed group is empty
when no laboratory history is present and every confirmed vintage becomes
an early window. Pure integer bookkeeping on the observed data, so it
carries no gradient.
"""
function confirmed_positivity_windows(confirmed_history, lab_history)
    empty = (; obs_days = Int[], obs_positives = Int[], obs_analysed = Int[],
        early_days = Int[], early_increments = Int[])
    isempty(confirmed_history.counts) && return empty
    cdays = confirmed_history.days
    ccounts = confirmed_history.counts
    ## Cumulative confirmed at (or most recently before) a grid day.
    function confirmed_at(day)
        i = searchsortedlast(cdays, day)
        return i == 0 ? 0 : Int(ccounts[i])
    end

    ## No laboratory denominators: every confirmed vintage is an early
    ## window scored through the modelled laboratory volume.
    if isempty(lab_history.counts)
        inc = diff(vcat(zero(eltype(ccounts)), collect(Int.(ccounts))))
        return (; obs_days = Int[], obs_positives = Int[],
            obs_analysed = Int[], early_days = collect(Int.(cdays)),
            early_increments = Int.(inc))
    end

    first_lab_day = Int(lab_history.days[1])
    ## Early windows: confirmed vintages up to the first laboratory date.
    early_days = Int[]
    early_increments = Int[]
    prev = 0
    for (i, d) in enumerate(cdays)
        Int(d) > first_lab_day && break
        push!(early_days, Int(d))
        push!(early_increments, Int(ccounts[i]) - prev)
        prev = Int(ccounts[i])
    end

    ## Observed windows: analysed increments between laboratory vintages
    ## (from the second onward) paired with confirmed increments, merging
    ## zero-denominator stalls forward.
    obs_days = Int[]
    obs_positives = Int[]
    obs_analysed = Int[]
    a_acc = 0
    c_acc = 0
    prev_conf = confirmed_at(first_lab_day)
    for j in 2:length(lab_history.days)
        d = Int(lab_history.days[j])
        a_inc = Int(lab_history.counts[j]) - Int(lab_history.counts[j - 1])
        conf_here = confirmed_at(d)
        a_acc += a_inc
        c_acc += conf_here - prev_conf
        prev_conf = conf_here
        if a_acc > 0 && c_acc <= a_acc
            push!(obs_days, d)
            push!(obs_positives, max(c_acc, 0))
            push!(obs_analysed, a_acc)
            a_acc = 0
            c_acc = 0
        end
    end
    if a_acc > 0
        push!(obs_days, Int(lab_history.days[end]))
        push!(obs_positives, clamp(c_acc, 0, a_acc))
        push!(obs_analysed, a_acc)
    elseif c_acc > 0 && !isempty(obs_analysed)
        obs_positives[end] = clamp(obs_positives[end] + c_acc, 0,
            obs_analysed[end])
    end
    return (; obs_days, obs_positives, obs_analysed, early_days,
        early_increments)
end

"""
Laboratory pipeline likelihood. Two streams driven by the shared renewal
onsets and the suspected-case pipeline from [`reported_cases_model`](@ref):

- Specimens received. The suspected daily pipeline (`p_drc`-scaled BVD
  onset-to-report signal plus the non-BVD background `λ_bg`) is carried
  through a sampled report-to-laboratory receipt delay and thinned by the
  tested fraction `τ_test`, giving the expected daily received-specimen
  volume. Its between-vintage increments are fitted against
  `tests_received_history` as observed `~` data with a NegativeBinomial
  sharing the surveillance dispersion `k`, identifying `τ_test` and the
  receipt delay.
- Confirmed positives. The confirmed counts are scored as a `Binomial` of
  the observed specimens-*analysed* denominator in each laboratory window
  ([`confirmed_positivity_windows`](@ref)), with a partially-pooled
  per-window positivity ([`confirmed_positivity_model`](@ref)).
  Conditioning the positives on the observed analysed denominator (rather
  than a modelled count scaled by `p_drc · s_test · τ_test`) removes the
  multiplicative ascertainment ridge that basin-split the joint, so the
  outbreak size is pinned by the deaths and exports streams while the
  laboratory positivity is free to track the noisy per-vintage data. The
  early confirmed vintages with no per-vintage analysed denominator
  (18-23 May) are scored as NegativeBinomial counts against the modelled
  laboratory volume with the same pooled positivity, so all the confirmed
  data is used and the early per-vintage shape informs the fit.

The tested fraction `τ_test` and background rate `λ_bg` come from
[`reported_cases_model`](@ref) so the suspected and laboratory streams
share them. Exposes the per-window positivity, the expected received and
confirmed totals at the cut-off, and the cut-off positivity as derived
quantities.
"""
@model function confirmed_cases_model(
        confirmed_history,
        confirmed_cases::Union{Missing, Integer},
        onsets::AbstractVector, k::Real, p_drc::Real,
        λ_bg::Real, τ_test::Real, bvd_reports_daily::AbstractVector;
        lab_history = (; days = Int[], counts = Int[]),
        tests_received_history = (; days = Int[], counts = Int[]),
        tests_received::Union{Missing, Integer} = missing,
        receipt = lab_delay_model(),
        positivity = confirmed_positivity_model)
    n = length(onsets)

    ## Received-specimen volume: the suspected pipeline carried through the
    ## receipt delay and thinned by the tested fraction, fit as a count.
    receipt_state ~ to_submodel(receipt)
    suspected_daily = p_drc .* bvd_reports_daily .+ λ_bg
    received_daily = τ_test .* convolve_delay(suspected_daily,
        receipt_state.pmf)
    rvobs = vintage_obs(tests_received_history, tests_received, n)
    received_inc = bin_increments(received_daily, rvobs.days)
    received_increments ~ to_submodel(
        vintage_increments_model(received_inc, rvobs.obs_increments, k))

    ## Confirmed positives in two groups sharing one partially-pooled
    ## positivity: early windows (no observed analysed) scored as counts
    ## against the modelled laboratory volume, observed windows scored as a
    ## Binomial of the observed analysed denominator.
    windows = confirmed_positivity_windows(confirmed_history, lab_history)
    n_early = length(windows.early_days)
    n_obs = length(windows.obs_analysed)
    nv = n_early + n_obs
    pos_state ~ to_submodel(positivity(nv))
    p_pos = pos_state.p_pos
    have_data = !ismissing(confirmed_cases)

    ## Early windows: confirmed increment ~ NegBinomial(positivity ×
    ## modelled analysed volume), the volume binned from the received
    ## series so the early counts inform the fit through partial pooling.
    early_p = p_pos[1:n_early]
    early_volume = bin_increments(received_daily, windows.early_days)
    early_mean = early_p .* early_volume
    early_obs = (have_data && n_early > 0) ? windows.early_increments : missing
    early_increments ~ to_submodel(
        vintage_increments_model(early_mean, early_obs, k))

    ## Observed windows: Binomial of the observed analysed denominator.
    obs_p = p_pos[(n_early + 1):nv]
    obs_positives = (have_data && n_obs > 0) ? collect(windows.obs_positives) :
                    missing
    confirmed_positives ~ to_submodel(
        confirmed_positives_model(obs_positives, windows.obs_analysed, obs_p))

    expected_received := safe_rate(sum(received_daily))
    ## Expected confirmed at the cut-off and the overall positivity, over
    ## both the modelled early volume and the observed analysed windows.
    denom = sum(early_volume) + float(sum(windows.obs_analysed))
    expected_positives = sum(early_mean) +
                         (n_obs > 0 ? sum(obs_p .* windows.obs_analysed) :
                          zero(eltype(p_pos)))
    expected_confirmed := safe_rate(expected_positives)
    p_positive := safe_rate(expected_positives) / safe_rate(denom)

    return (; τ_test, λ_bg, p_pos, windows, received_daily,
        expected_received, expected_confirmed, p_positive)
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

"""
Laboratory-confirmed-deaths likelihood. Confirmed deaths are a thinning of
the observed suspected deaths: the cut-off confirmed-death count is scored
as a `Binomial` of the suspected-death total `total_deaths`, with a
confirmation probability linked to the suspected-case BVD composition
`q_susp` (the BVD share of the expected suspected total, from the
[`reported_cases_model`](@ref) pipeline) and enriched on the odds scale by
the sampled `m_death` ([`confirmed_death_enrichment_model`](@ref)):

```math
p = \\operatorname{logistic}(\\operatorname{logit}(q_\\text{susp}) +
    \\log m_\\text{death}).
```

The odds-scale enrichment keeps `p` in `(0, 1)` without a hard clamp, and
ties the death-confirmation rate to the same composition that drives the
case streams, so a confirmed-death observation informs the background
`λ_bg` and ascertainment `p_drc` rather than introducing a free rate. The
confirmed-death series is flat over the fitted window, so a single cut-off
binomial captures its information; the suspected-death total is the
denominator. Returns the enrichment, the composition, the confirmation
probability and the expected confirmed-death count.
"""
@model function confirmed_deaths_model(
        confirmed_deaths::Union{Missing, Integer},
        total_deaths::Union{Missing, Integer},
        expected_deaths::Real,
        bvd_reports_daily::AbstractVector, p_drc::Real, λ_bg::Real;
        enrichment = confirmed_death_enrichment_model())
    enr_state ~ to_submodel(enrichment)
    m_death = enr_state.m_death

    n = length(bvd_reports_daily)
    bvd_total = p_drc * sum(bvd_reports_daily)
    q_susp := safe_rate(bvd_total) / safe_rate(bvd_total + λ_bg * n)
    qc = clamp(q_susp, eps(typeof(q_susp)), one(q_susp) - eps(typeof(q_susp)))
    p_death_conf := logistic(logit(qc) + log(m_death))

    ## Suspected-death denominator: the observed total when conditioned, or
    ## the expected suspected deaths on the predictive-generator path.
    n_deaths = ismissing(total_deaths) ?
               max(round(Int, safe_rate(expected_deaths)), 0) :
               Int(total_deaths)
    confirmed_deaths ~ Binomial(n_deaths, p_death_conf)
    expected_confirmed_deaths := p_death_conf * n_deaths

    return (; m_death, q_susp, p_death_conf, expected_confirmed_deaths)
end
