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
Expand a per-vintage rate vector `rate` onto a length-`n` daily grid,
assigning each day the rate of the vintage window it falls in. The
windows are delimited by the ascending day indices `days` (1-based into
the grid); day `t ≤ days[1]` takes `rate[1]`, a day in `(days[i-1],
days[i]]` takes `rate[i]`, and any day beyond the last vintage takes the
last rate (a flat carry-forward of the final window). When `days` is
empty the whole grid takes `rate[1]` if present, else zero, so a scalar
background is recovered. Pure and AD-transparent; the element type
follows `rate`. Used to turn the per-vintage background random effect
([`background_re_model`](@ref)) into the additive daily background the
suspected-case and suspected-death streams consume.
"""
function expand_vintage_rate(rate::AbstractVector,
        days::AbstractVector{<:Integer}, n::Integer)
    T = eltype(rate)
    out = Vector{T}(undef, n)
    if isempty(days) || isempty(rate)
        fill!(out, isempty(rate) ? zero(T) : rate[1])
        return out
    end
    nv = length(rate)
    prev = 0
    @inbounds for (i, d) in enumerate(days)
        i > nv && break
        hi = clamp(Int(d), 0, n)
        for t in (prev + 1):hi
            out[t] = rate[i]
        end
        prev = hi
    end
    ## Carry the final window's rate forward over any tail days beyond the
    ## last vintage edge.
    last_rate = rate[min(length(days), nv)]
    @inbounds for t in (prev + 1):n
        out[t] = last_rate
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
        death_background = nothing,
        background_re = nothing,
        onset_to_death = censored_delay_model(40;
            mean_prior = truncated(Normal(11.2, 2.0); lower = 1),
            sd_prior = truncated(Normal(5.4, 1.5); lower = 1)))
    cfr_state ~ to_submodel(cfr)
    od_state ~ to_submodel(onset_to_death)
    CFR = cfr_state.CFR
    bvd_deaths_daily = CFR .* convolve_delay(onsets, od_state.pmf)

    n = length(bvd_deaths_daily)
    vobs = vintage_obs(deaths_history, total_deaths, n)

    ## Daily non-BVD background deaths. The renewal default has no
    ## background (`death_background === nothing`), so the deaths stream is
    ## pure BVD. With a scalar `death_background` the background is constant
    ## over the grid; with a per-vintage `background_re` it is the
    ## regularised time-varying random effect expanded to a daily series,
    ## sharing the suspected-case background's structure.
    if background_re !== nothing
        bg_state ~ to_submodel(background_re(length(vobs.days)))
        λ_bg_death = bg_state.λ_mu
        bg_death_sigma = bg_state.σ_bg
        bg_death_daily = expand_vintage_rate(bg_state.λ, vobs.days, n)
    elseif death_background !== nothing
        dbg_state ~ to_submodel(death_background)
        λ_bg_death = dbg_state.λ_bg_death
        bg_death_sigma = zero(λ_bg_death)
        bg_death_daily = fill(λ_bg_death, n)
    else
        λ_bg_death = zero(CFR)
        bg_death_sigma = zero(CFR)
        bg_death_daily = fill(zero(CFR), n)
    end

    deaths_daily = bvd_deaths_daily .+ bg_death_daily

    modelled_increments = bin_increments(deaths_daily, vobs.days)
    death_increments ~ to_submodel(
        vintage_increments_model(modelled_increments, vobs.obs_increments, k))

    raw_total = sum(deaths_daily)
    expected_deaths_T := safe_rate(raw_total)
    bg_death_total = sum(bg_death_daily)

    return (; CFR, od_pmf = od_state.pmf, deaths_daily, bvd_deaths_daily,
        expected_deaths_T, λ_bg_death, bg_death_sigma, bg_death_daily,
        bg_death_total)
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
        background_re = nothing,
        onset_to_report = censored_delay_model(30;
            mean_prior = truncated(Normal(4.5, 1.5); lower = 1),
            sd_prior = truncated(Normal(3.6, 1.2); lower = 1)))
    pos_state ~ to_submodel(positivity)
    report_state ~ to_submodel(onset_to_report)
    λ_bg = pos_state.λ_bg
    τ_test = pos_state.τ_test
    report_pmf = report_state.pmf

    ## Unit-ascertainment BVD onset-to-report daily series, reused by the
    ## confirmed stream.
    bvd_reports_daily = convolve_delay(onsets, report_pmf)

    n = length(bvd_reports_daily)
    vobs = vintage_obs(reported_history, reported_cases, n)

    ## Daily non-BVD background. With `background_re === nothing` this is
    ## the constant scalar `λ_bg` over the grid (the renewal default). When
    ## a per-vintage random effect is injected, the scalar baseline is
    ## perturbed window-by-window and expanded to a daily series; the
    ## baseline `λ_bg` from `positivity` is overridden by the random
    ## effect's `λ_mu`, and the per-vintage rates are tightly pooled toward
    ## it so the background stays a regularised minority of suspected cases.
    if background_re === nothing
        λ_bg_base = λ_bg
        bg_sigma = zero(λ_bg)
        bg_daily = fill(λ_bg, n)
    else
        bg_state ~ to_submodel(background_re(length(vobs.days)))
        λ_bg_base = bg_state.λ_mu
        bg_sigma = bg_state.σ_bg
        bg_daily = expand_vintage_rate(bg_state.λ, vobs.days, n)
    end

    ## Suspected daily cases add the p_drc-scaled BVD signal and the
    ## non-BVD background.
    reports_daily = p_drc .* bvd_reports_daily .+ bg_daily

    modelled_increments = bin_increments(reports_daily, vobs.days)
    reported_increments ~ to_submodel(
        vintage_increments_model(modelled_increments, vobs.obs_increments, k))

    raw_total = sum(reports_daily)
    expected_reports := safe_rate(raw_total)

    ## Implied per-suspected positivity at the cut-off: BVD share of the
    ## expected suspected total.
    bvd_total = p_drc * sum(bvd_reports_daily)
    positivity := safe_rate(bvd_total) / expected_reports

    ## Cumulative background suspected cases over the grid, exposed for
    ## comparison with the observed suspected total.
    bg_total = sum(bg_daily)

    return (; p_drc, λ_bg = λ_bg_base, τ_test, report_pmf, bvd_reports_daily,
        reports_daily, expected_reports, positivity, bg_daily, bg_sigma,
        bg_total)
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
        bg_daily::AbstractVector, τ_test::Real,
        bvd_reports_daily::AbstractVector;
        lab_history = (; days = Int[], counts = Int[]),
        tests_received_history = (; days = Int[], counts = Int[]),
        tests_received::Union{Missing, Integer} = missing,
        receipt = lab_delay_model(),
        positivity = confirmed_positivity_model,
        positivity_link::Symbol = :free,
        severity_enrichment = severity_enrichment_model())
    n = length(onsets)

    ## Received-specimen volume: the suspected pipeline carried through the
    ## receipt delay and thinned by the tested fraction, fit as a count.
    ## `bg_daily` is the per-day non-BVD background from the suspected-case
    ## stream (a constant series when the scalar background is used, the
    ## per-vintage random effect when it is on).
    receipt_state ~ to_submodel(receipt)
    suspected_daily = p_drc .* bvd_reports_daily .+ bg_daily
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
    have_data = !ismissing(confirmed_cases)

    ## Per-window tested BVD share `p_pos`. Two links:
    ## `:free` (default) — a free partially-pooled per-window random effect
    ## ([`confirmed_positivity_model`](@ref)), decoupled from `λ_bg`.
    ## `:composition` — the tested share is the suspect-pool composition
    ## `φ_v = (p_drc·BVD)_v / ((p_drc·BVD)_v + λ_bg_v)` over each laboratory
    ## window, upsampled by a decaying severity enrichment δ0 (see
    ## [`severity_enrichment_model`](@ref)), so the lab positivity identifies
    ## the background `λ_bg` rather than absorbing it into a free curve.
    window_days = vcat(windows.early_days, windows.obs_days)
    if positivity_link === :composition
        enrich_state ~ to_submodel(severity_enrichment, false)
        δ0 = enrich_state.δ0
        decay_scale = enrich_state.decay_scale
        ## Suspect-pool composition over each window from the SAME suspected-
        ## stream series that drive `reported_cases_model`.
        bvd_window = bin_increments(p_drc .* bvd_reports_daily, window_days)
        bg_window = bin_increments(bg_daily, window_days)
        Tt = eltype(bvd_window)
        ## Testing clock: cumulative modelled analysed volume at each window.
        vol_window = bin_increments(received_daily, window_days)
        c_window = cumsum(vol_window)
        lo = convert(Tt, 1e-8)
        hi = one(Tt) - lo
        ## Floor the decay scale so a near-zero `decay_scale` draw cannot make
        ## the clock ratio `0/0` (NaN) and break the downstream Binomial.
        dscale = max(convert(Tt, decay_scale), one(Tt))
        p_pos = map(eachindex(window_days)) do i
            ## Pool composition φ = (p_drc·BVD) / ((p_drc·BVD) + λ_bg) over the
            ## window, guarded against a zero/negative denominator.
            num = bvd_window[i]
            den = bvd_window[i] + bg_window[i]
            ratio = num / (den + lo)
            φ = clamp(isfinite(ratio) ? ratio : convert(Tt, 0.5), lo, hi)
            δ_i = convert(Tt, δ0) * exp(-c_window[i] / dscale)
            q = logistic(logit(φ) + δ_i)
            ## Final guard: clamp into (0,1) and replace any non-finite value
            ## with the pool composition so the confirmed Binomial always sees
            ## a valid probability even under an AD perturbation.
            qf = isfinite(q) ? q : φ
            clamp(qf, lo, hi)
        end
    else
        pos_state ~ to_submodel(positivity(nv))
        p_pos = pos_state.p_pos
    end

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

    return (; τ_test, bg_daily, p_pos, windows, received_daily,
        expected_received, expected_confirmed, p_positive)
end

"""
Uganda exports likelihood (geographic spread). The exports stream is
travel-gated, so the at-risk clock starts at infection: a traveller
moves and is exported during incubation (pre-symptomatic) and stays at
risk of being exported and detected abroad only until the
infection→detection delay has elapsed. The expected detected exports by
the cut-off therefore accumulate the per-capita travel rate `q =
daily_travellers / source_population` over the at-risk PERSON-TIME, not
over single-day onset events:

```math
\\mathbb{E}[\\text{exports}(T)] = p_\\text{uganda}\\, q
    \\sum_{s=1}^{T} \\text{prevalence}(s),
\\qquad
\\text{prevalence}(s) = C(s) - \\text{detected}(s),
```

with `C(s)` the cumulative infections and `detected(s)` the cumulative
infections that have already completed the infection→detection delay by
day `s`. The infection→detection delay is the sampled onset-to-detection
delay convolved with the shared incubation PMF (so incubation sits inside
it, keyed to infection like `C(s)`); `detected` is the running sum of
`convolve_delay(infections, f_det)`. Summing the daily at-risk
prevalence is the discrete person-time integral, the discrete analogue of
the integral model's at-risk person-time export integral; the earlier
onset-incidence form summed `q · onsets`, charging each case only a
single day of travel risk and so under-counting exports by roughly the
mean at-risk dwell time. Fitted with Poisson (Uganda's stream is small).
Samples the traveller volume and the onset-to-detection delay via
injected submodels. The onset-to-detection prior is centred on the Ebola
onset-to-hospitalisation delay (mean 5.0 d, SD 4.7 d; WHO Ebola Response
Team 2014, NEJM), the delay from symptom onset to detection at a point of
entry abroad. Returns the expected count, the per-capita travel rate,
the daily at-risk prevalence and the daily expected export incidence for
reuse by [`exports_deaths_model`](@ref).
"""
@model function exports_model(
        exported_cases::Union{Missing, Integer},
        infections::AbstractVector, p_uganda::Real;
        incubation_pmf::AbstractVector,
        source_population::Real = ITURI_POPULATION,
        traveller = traveller_volume_model(),
        onset_to_detection = censored_delay_model(30;
            mean_prior = truncated(Normal(5.0, 2.0); lower = 1),
            sd_prior = truncated(Normal(4.7, 1.5); lower = 1)))
    travel_state ~ to_submodel(traveller)
    daily_travellers = travel_state.daily_travellers
    q = daily_travellers / source_population

    detect_state ~ to_submodel(onset_to_detection)
    ## Infection→detection delay: onset→detection convolved with the shared
    ## incubation PMF, so the survival clock runs from infection.
    f_det = convolve_pmf(incubation_pmf, detect_state.pmf)
    detected_daily = convolve_delay(infections, f_det)
    ## At-risk prevalence (person-days): infected but not yet detected.
    prevalence = cumsum(infections) .- cumsum(detected_daily)
    export_prevalence = p_uganda .* q .* prevalence

    raw_exports = sum(export_prevalence)
    expected_exports_T := safe_rate(raw_exports)
    exported_cases ~ Poisson(expected_exports_T)

    return (; p_uganda, daily_travellers, q, prevalence,
        export_prevalence,
        expected_exports = expected_exports_T)
end

"""
Deaths-among-detected-exports likelihood. Deaths accrue among the at-risk
export person-time of [`exports_model`](@ref): each day's at-risk
prevalence is weighted by the infection→death CDF at its remaining age to
the cut-off and summed, scaled by the CFR, the discrete analogue of the
integral model's `∫ C(s)·S_det(T−s)·F_death(T−s) ds`. The onset-to-death
PMF shared from [`deaths_model`](@ref) is convolved with the incubation
PMF to give the infection→death distribution, keyed to infection like the
prevalence. Fitted with Poisson.
"""
@model function exports_deaths_model(
        exports_deaths::Union{Missing, Integer},
        export_prevalence::AbstractVector, CFR::Real,
        od_pmf::AbstractVector, incubation_pmf::AbstractVector)
    n = length(export_prevalence)
    ## Infection→death CDF by age (age 0 = same day).
    fd_pmf = convolve_pmf(incubation_pmf, od_pmf)
    death_cdf = cumsum(fd_pmf)
    ## Weight each day's at-risk prevalence by the death CDF at its
    ## remaining age to the cut-off (day n).
    acc = zero(eltype(export_prevalence))
    @inbounds for s in 1:n
        age = n - s
        w = age + 1 <= length(death_cdf) ? death_cdf[age + 1] :
            (isempty(death_cdf) ? zero(eltype(death_cdf)) : death_cdf[end])
        acc += export_prevalence[s] * w
    end
    raw = CFR * acc

    expected_exports_deaths_T := safe_rate(raw)
    exports_deaths ~ Poisson(expected_exports_deaths_T)

    return (; expected_exports_deaths_T)
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
        bvd_reports_daily::AbstractVector, p_drc::Real,
        bg_daily::AbstractVector;
        enrichment = confirmed_death_enrichment_model())
    enr_state ~ to_submodel(enrichment)
    m_death = enr_state.m_death

    bvd_total = p_drc * sum(bvd_reports_daily)
    q_susp := safe_rate(bvd_total) / safe_rate(bvd_total + sum(bg_daily))
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
