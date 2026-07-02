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
    ## Guard the dispersion `r`, not just `p`: a fit can drive `k` to `0`
    ## (or underflow / non-finite), and `NegativeBinomial(0, p)` throws
    ## `DomainError: r > 0`, which aborts the whole fit inside a Mooncake
    ## gradient rather than rejecting the step. Floor it at `eps` so the
    ## likelihood stays defined everywhere the sampler can reach.
    r = (isfinite(k) && k > zero(k)) ? k : eps(typeof(k))
    p_raw = r / (r + max(μ, eps(typeof(μ))))
    p = isfinite(p_raw) ?
        clamp(p_raw, eps(typeof(r)), one(r) - eps(typeof(r))) :
        eps(typeof(r))
    return NegativeBinomial(r, p)
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
    ## Concrete zero for empty bins. `zero(first(daily))` dispatches on the
    ## runtime element type, so it works even when the container has widened
    ## to `Vector{Any}` in predict mode (`zero(eltype(daily))` would call
    ## `zero(::Type{Any})` and error).
    zfill = isempty(daily) ? zero(eltype(daily)) : zero(@inbounds daily[begin])
    prev = 0
    @inbounds for (i, d) in enumerate(days)
        hi = clamp(Int(d), 0, n)
        lo = clamp(prev, 0, n)
        out[i] = hi > lo ? sum(@view daily[(lo + 1):hi]) : zfill
        prev = hi
    end
    return out
end

"""
Zero a daily series `v` before grid day `start` (1-based), modelling a
process that did not exist before `start`. Used to gate the laboratory
analysed-specimen CAPACITY before testing began: no specimens are analysed
before the testing system exists, so the modelled analysed volume (and the
confirmed counts derived from it) must not accrue over the pre-surveillance
cryptic phase. The suspected-case and suspected-death streams are NOT gated —
those counts did accumulate over the cryptic phase, so their first per-vintage
bin legitimately rolls from grid day 1. `start ≤ 1` returns `v` unchanged.
Pure and AD-transparent; the element type follows `v`. Callers concretise any
`Vector{Any}` (predict mode) before gating, so `zero(eltype(v))` is well
defined.
"""
function gate_before(v::AbstractVector, start::Integer)
    start <= 1 && return v
    T = eltype(v)
    out = Vector{T}(undef, length(v))
    @inbounds for i in eachindex(v)
        out[i] = i < start ? zero(T) : v[i]
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
Group a sorted event-day list `event_days` (grid day-indices, ascending,
one entry per dated event, possibly with repeats) into its unique days and
the per-day occupancy count, clamped to `[1, n]`. Returns `(days, counts)`
with `days` the unique detection/death days and `counts[i]` the number of
events on `days[i]`. Used by [`exports_model`](@ref) and the export-deaths
likelihood to place one Poisson term per detection/death day with an
observed count equal to that day's occupancy, so simultaneous imports on
one day share a single edge. `event_days` must be non-empty.
"""
function dated_event_bins(event_days::AbstractVector{<:Integer}, n::Integer)
    clamped = [clamp(Int(d), 1, Int(n)) for d in event_days]
    days = unique(clamped)
    counts = Int[count(==(d), clamped) for d in days]
    return days, counts
end

"""
Per-day Poisson likelihood for a dated event series. Scores the observed
per-day counts `obs` against the modelled per-day means `means` with one
Poisson term each, NaN/Inf-safe via [`safe_rate`](@ref). When `obs` is
`missing` the counts are sampled, the predictive-generator path; the
indexed `counts[i]` keeps the predict keys (`<prefix>.counts[i]`)
replicable. Used by [`exports_model`](@ref) and the export-deaths
likelihood for the dated Uganda export series.
"""
@model function dated_poisson_model(means::AbstractVector,
        counts::Union{Missing, AbstractVector{<:Integer}})
    n = length(means)
    ## `counts` is the model argument scored on the LHS of `~`, so a supplied
    ## vector is observed data DynamicPPL conditions on and a `missing`
    ## argument is sampled (the predictive-generator path). The indexed
    ## `counts[i]` keeps the predict keys (`<prefix>.counts[i]`) replicable.
    if ismissing(counts)
        counts = Vector{Union{Missing, Int}}(missing, n)
    end
    for i in 1:n
        counts[i] ~ Poisson(safe_rate(means[i]))
    end
    return (; means, counts)
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
Right-censored occupancy likelihood. Each day's observed bed count `obs[i]`
is a NegBinomial around the latent bed DEMAND `means[i]`, right-censored at
the effective capacity `ceilings[i]`: below the ceiling the count reflects
demand directly, and once demand reaches the ceiling the count is censored
there. Unlike a smooth saturation of the mean, the censored tail probability
still depends on the demand above the ceiling, so demand stays identified
when beds are full rather than the occupancy going flat in demand. A
`missing` `obs` samples (the predictive path). Shares the surveillance
dispersion `k`.
"""
@model function censored_occupancy_model(means::AbstractVector,
        ceilings::AbstractVector,
        obs::Union{Missing, AbstractVector{<:Integer}}, k::Real)
    n = length(means)
    if ismissing(obs)
        obs = Vector{Union{Missing, Int}}(missing, n)
    end
    for i in 1:n
        obs[i] ~ censored(safe_nbinomial(k, safe_rate(means[i]));
            upper = safe_rate(ceilings[i]))
    end
    return (; means, ceilings, obs)
end

"""
    censoring_cap(iso_days, iso_obs, capacity_history)

Per-report-day right-censoring bound for the isolation-occupancy likelihood,
built from DATA only (no sampled parameter), so the censored NegativeBinomial
never sits against a moving `-Inf` wall. Each isolation day takes the nearest
recorded implied-capacity value (`capacity_history`, the occupancy / reported
occupancy-rate series), floored at that day's observed occupancy so the bound
is never below the count (a count above the cap has zero probability under the
censored NB, the discontinuity that drives the divergences). With the cap
fixed, the censored likelihood is smooth in the sampled demand: below the cap
the count identifies demand directly, and at the cap it contributes the
one-sided `demand ≥ capacity` tail. A `capacity_history` with no counts
(empty, or days-only as the predictive generator supplies) returns a large
finite cap far above any bed count, so the censoring is a no-op and the
likelihood is the plain NB (a literal `Inf` cannot be used because
[`safe_rate`](@ref) maps non-finite values to `eps`); a `missing` `iso_obs`
skips the occupancy floor.

The latent bed demand itself is left UNCENSORED: the cap enters only as the
observation bound, so demand above capacity is carried by the renewal /
length-of-stay demand model and its priors, not by the bound. Demand above a
saturated capacity is only partially identified from occupancy (occupancy can
say "demand was at least the beds filled", not how much more), so the bed
shortfall is a prior/model-informed quantity, not pinned by the occupancy data
— see [`treatment_flow_model`](@ref).
"""
function censoring_cap(iso_days, iso_obs, capacity_history)
    cdays = Int.(capacity_history.days)
    ccounts = Float64.(capacity_history.counts)
    ## Large finite no-op cap (bed counts are in the hundreds), used when there
    ## is no recorded capacity. Not `Inf`: `safe_rate` maps non-finite to eps.
    nocap = 1.0e6
    ## Censor only where recorded capacity COUNTS exist. The predictive
    ## generator passes the capacity history days-only (counts emptied), so the
    ## guard is on the counts, not the days — otherwise the nearest-day lookup
    ## would index an empty counts vector.
    have_cap = !isempty(ccounts)
    m = length(iso_days)
    cap = Vector{Float64}(undef, m)
    @inbounds for (i, d) in enumerate(iso_days)
        di = Int(d)
        c = if !have_cap
            nocap
        else
            ## nearest recorded capacity by day distance (data lookup, no AD)
            j = argmin(abs.(cdays .- di))
            ccounts[j]
        end
        o = iso_obs === missing ? 0.0 : Float64(iso_obs[i])
        cap[i] = max(c, o)
    end
    return cap
end

"""
DRC suspected-deaths likelihood, per-vintage time series. Convolves the
daily onsets with the sampled onset-to-death delay, scales by the CFR and
the death ascertainment `p_death`, and reads the modelled cumulative deaths
at each vintage day off the daily series, fitting the between-vintage
increments as observed `~` data with a NegativeBinomial sharing the
surveillance dispersion `k` ([`surveillance_dispersion_model`](@ref)). The
death history ends at the cut-off, so the cut-off total is the final
increment and is not scored separately. Samples the onset-to-death delay,
the CFR and the death ascertainment via injected submodels. The onset-to-
death prior is the convolution of the two atomic line-list components
(onset-to-admission and admission-to-death; the `bdbv-linelist-analysis`
submodule), implied mean ≈ 12.8 d, SD ≈ 7.0 d, the same source the integral
model used.

The expected BVD suspected deaths are `p_death · CFR` of the onset-to-death-
convolved infections: a fatal BVD infection is counted as a suspected death
only when ascertained, the death analogue of the suspected-case
ascertainment `p_drc` (see [`death_ascertainment_model`](@ref)). The default
death ascertainment prior is informative and high (centre ≈ 0.9), reflecting
that a death is more reliably reported than a living suspect.

Non-BVD background suspected deaths are added on top. The joint passes the
per-day non-BVD suspected-case background `case_bg_daily` and a background CFR
submodel ([`background_cfr_model`](@ref)): the background deaths are `cfr_bg`
times that case background, lagged by the onset-to-death delay so a background
death follows its background case the way the BVD deaths follow the onsets.
The death background therefore inherits a level and time profile from the
identified case background without a second free, outbreak-size-degenerate
rate. The `death_background` ([`death_background_model`](@ref)) scalar, the
per-vintage `background_re` and the pure-BVD stream are sensitivity fallbacks.

An optional `suspected_daily_deaths_history` adds the post-26 May daily new
suspected deaths ("cas suspects du jour N (M deces)"): per-day counts of
suspected deaths scored against the modelled daily suspected-death series
`deaths_daily` at each report day (a single-day mean, not a between-vintage
increment) with NegativeBinomials sharing `k`. This is the deaths analogue of
the suspected-case daily inflow ([`reported_cases_model`](@ref)): it continues
the suspected-death signal where the cumulative `deaths_history` stops once
INSP stopped publishing a national suspected-death total. The inflow is a
genuine per-day count that never falls, so it fits directly; its days fall
strictly after the cumulative series ends, so the two suspected-death
likelihoods cover disjoint days. Empty by default.

Returns the cut-off expected count, the daily death series (total and the
BVD and background components), the onset-to-death PMF, the CFR, the death
ascertainment and the background CFR for reuse by
[`confirmed_deaths_model`](@ref) and [`exports_deaths_model`](@ref).
"""
@model function deaths_model(
        deaths_history,
        total_deaths::Union{Missing, Integer},
        onsets::AbstractVector, k::Real;
        suspected_daily_deaths_history = (; days = Int[], counts = Int[]),
        cfr = cfr_model(),
        ascertainment = death_ascertainment_model(),
        case_bg_daily = nothing,
        background_cfr = background_cfr_model(),
        death_background = nothing,
        background_re = nothing,
        ## nmax covers 98% of the convolved onset->death sum (the two atomic
        ## Gammas moment-matched to a single Gamma only for the truncation).
        onset_to_death = onset_to_death_model(cdf_nmax(Gamma(3.33, 3.83));
            oa_alpha_prior = truncated(Normal(1.178, 0.285); lower = 0.01),
            oa_theta_prior = truncated(Normal(3.694, 1.198); lower = 0.1),
            ad_alpha_prior = truncated(Normal(2.151, 0.604); lower = 0.01),
            ad_theta_prior = truncated(Normal(3.906, 1.381); lower = 0.1)))
    cfr_state ~ to_submodel(cfr)
    od_state ~ to_submodel(onset_to_death)
    asc_state ~ to_submodel(ascertainment)
    CFR = cfr_state.CFR
    p_death = asc_state.p_death
    bvd_deaths_daily = (p_death * CFR) .* convolve_delay(onsets, od_state.pmf)

    n = length(bvd_deaths_daily)
    vobs = vintage_obs(deaths_history, total_deaths, n)

    ## Daily non-BVD background deaths: the background CFR `cfr_bg` applied to
    ## the non-BVD suspected-case background `case_bg_daily`, lagged by the
    ## onset-to-death delay so a background death follows its background case
    ## the way the BVD deaths follow the onsets. `λ_bg_death` is the mean daily
    ## background death rate. The `background_re`, `death_background` and
    ## pure-BVD branches are sensitivity fallbacks.
    if case_bg_daily !== nothing
        bgcfr_state ~ to_submodel(background_cfr)
        cfr_bg = bgcfr_state.cfr_bg
        bg_death_daily = cfr_bg .* convolve_delay(case_bg_daily, od_state.pmf)
        λ_bg_death = sum(bg_death_daily) / n
        bg_death_sigma = zero(CFR)
    elseif background_re !== nothing
        bg_state ~ to_submodel(background_re(n))
        cfr_bg = zero(CFR)
        λ_bg_death = bg_state.λ_mu
        bg_death_sigma = bg_state.σ_bg
        bg_death_daily = bg_state.λ
    elseif death_background !== nothing
        dbg_state ~ to_submodel(death_background)
        cfr_bg = zero(CFR)
        λ_bg_death = dbg_state.λ_bg_death
        bg_death_sigma = zero(λ_bg_death)
        bg_death_daily = fill(λ_bg_death, n)
    else
        cfr_bg = zero(CFR)
        λ_bg_death = zero(CFR)
        bg_death_sigma = zero(CFR)
        bg_death_daily = fill(zero(CFR), n)
    end

    deaths_daily = bvd_deaths_daily .+ bg_death_daily

    modelled_increments = bin_increments(deaths_daily, vobs.days)
    death_increments ~ to_submodel(
        vintage_increments_model(modelled_increments, vobs.obs_increments, k))

    ## Daily new suspected deaths ("cas suspects du jour N (M deces)"): per-day
    ## counts scored against the modelled daily suspected-death series at each
    ## report day. The mean for day `d` is the single-day `deaths_daily[d]`
    ## (clamped into the grid), NOT a between-vintage increment — this is a
    ## genuine daily count, so it never differences a falling cumulative. Empty
    ## by default; a `missing` count vector samples (the predictive path). The
    ## deaths analogue of the suspected-case daily inflow.
    sdd_days = suspected_daily_deaths_history.days
    sdd_modelled = [deaths_daily[clamp(Int(d), 1, n)] for d in sdd_days]
    sdd_obs = isempty(suspected_daily_deaths_history.counts) ? missing :
              collect(Int.(suspected_daily_deaths_history.counts))
    suspected_daily_deaths ~ to_submodel(
        vintage_increments_model(sdd_modelled, sdd_obs, k))

    raw_total = sum(deaths_daily)
    expected_deaths_T := safe_rate(raw_total)
    bg_death_total = sum(bg_death_daily)

    return (; CFR, p_death, cfr_bg, od_pmf = od_state.pmf, deaths_daily,
        bvd_deaths_daily, expected_deaths_T, λ_bg_death, bg_death_sigma,
        bg_death_daily, bg_death_total)
end

"""
DRC suspected (reported) cases likelihood, per-vintage time series. The
expected daily suspected cases are a BVD-driven onset-to-report
convolution scaled by the DRC ascertainment fraction `p_drc`, plus an
additive non-BVD background rate `λ_bg` per day (so a suspected case need
not be a true BVD infection). Reads the modelled cumulative suspected
cases at each vintage day off the daily series and fits the increments
with a NegativeBinomial sharing `k`.

An optional `suspected_daily_history` adds the post-26 May daily new-suspect
inflow ("nouveaux cas suspects du jour"): per-day counts of newly reported
suspects scored against the modelled daily suspected series `reports_daily`
at each report day (a single-day mean, not a between-vintage increment) with
NegativeBinomials sharing `k`. This continues the suspected signal where the
cumulative `reported_history` stops, once INSP began reclassifying suspects
downward and the cumulative total fell. The inflow is a genuine per-day
incidence that never falls, so it fits directly; it shares the suspect
pipeline and dispersion with the cumulative stream and is empty by default.
Its days fall strictly after the cumulative series ends, so the two suspected
likelihoods cover disjoint days.

The background and testing fraction
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
        suspected_daily_history = (; days = Int[], counts = Int[]),
        positivity = test_positivity_model(),
        background_re = nothing,
        ## Onset to a suspected case being detected/reported, from the
        ## line-list onset→admission delay (d_oa, ~4 d): a case enters
        ## surveillance when first formally seen, so one delay serves both the
        ## suspect-case and export streams. The line-list onset→notification
        ## delay (~20 d) is NOT used — we assume it reflects a longer pathway
        ## (likely confirmation and administrative processing), though what it
        ## captures is uncertain.
        onset_to_report = gamma_delay_model(cdf_nmax(Gamma(1.178, 3.694));
            alpha_prior = truncated(Normal(1.178, 0.285); lower = 0.01),
            theta_prior = truncated(Normal(3.694, 1.198); lower = 0.1)))
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
    ## `background_re` is injected it is the smooth daily random-walk
    ## background ([`background_walk_model`](@ref)): a length-`n` daily series
    ## that is zero before the surveillance onset and follows a tight lognormal
    ## random walk after it, so the baseline `λ_bg` from `positivity` is
    ## overridden by the walk's level `λ_mu` and the background varies smoothly
    ## day-to-day rather than in per-vintage steps.
    if background_re === nothing
        λ_bg_base = λ_bg
        bg_sigma = zero(λ_bg)
        bg_daily = fill(λ_bg, n)
    else
        bg_state ~ to_submodel(background_re(n))
        λ_bg_base = bg_state.λ_mu
        bg_sigma = bg_state.σ_bg
        bg_daily = bg_state.λ
    end

    ## Suspected daily cases add the p_drc-scaled BVD signal and the
    ## non-BVD background.
    reports_daily = p_drc .* bvd_reports_daily .+ bg_daily

    modelled_increments = bin_increments(reports_daily, vobs.days)
    reported_increments ~ to_submodel(
        vintage_increments_model(modelled_increments, vobs.obs_increments, k))

    ## Daily new-suspect inflow ("nouveaux cas suspects du jour"): per-day
    ## counts scored against the modelled daily series at each report day. The
    ## mean for day `d` is the single-day `reports_daily[d]` (clamped into the
    ## grid), NOT a between-vintage increment — this is a genuine daily
    ## incidence, so it never differences a falling cumulative. Empty by
    ## default; a `missing` count vector samples (the predictive path).
    sd_days = suspected_daily_history.days
    sd_modelled = [reports_daily[clamp(Int(d), 1, n)] for d in sd_days]
    sd_obs = isempty(suspected_daily_history.counts) ? missing :
             collect(Int.(suspected_daily_history.counts))
    suspected_daily ~ to_submodel(
        vintage_increments_model(sd_modelled, sd_obs, k))

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
Late confirmed-vintage likelihood, one per-window confirmed increment over
the days after the last cumulative laboratory date. Each window scores its
increment two ways depending on whether a 24h analysed denominator was
published for that day (`analysed[i] > 0`): an anchored day is a
`Binomial(analysed[i], p_pos[i])` of the observed denominator (like an
observed window), a unanchored day a `NegativeBinomial` of the modelled
laboratory volume `modelled[i]` sharing the surveillance dispersion `k`.
Keeping both in one submodel preserves a single per-window predict-key
vector (`<prefix>.increments[i]`) over all late days in order, so the
posterior-predictive trajectory reconstructs without interleaving.

`increments` is the model argument on the LHS of `~`: a supplied vector is
observed data DynamicPPL conditions on and a `missing` argument is sampled
(the predictive-generator path). A `Vector{Union{Missing, Int}}` with some
entries `missing` scores only the present ones, used to observe the
anchored days while leaving the unanchored days latent under the
no-extrapolation probe.
"""
@model function late_confirmed_model(
        increments::Union{Missing, AbstractVector},
        modelled::AbstractVector, analysed::AbstractVector{<:Integer},
        p_pos::AbstractVector, k::Real)
    n = length(modelled)
    if ismissing(increments)
        increments = Vector{Union{Missing, Int}}(missing, n)
    end
    for i in 1:n
        if analysed[i] > 0
            increments[i] ~ Binomial(analysed[i], p_pos[i])
        else
            increments[i] ~ safe_nbinomial(k, safe_rate(modelled[i]))
        end
    end
    return (; modelled, increments)
end

"""
Align the confirmed-case counts onto the laboratory windows, splitting
them into three non-overlapping groups so all the confirmed data is used:

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
- *Late* windows, the confirmed vintages strictly after the last
  laboratory date (the days after national testing stops, INSP's
  confirmed-only format with no national analysed-specimen denominators).
  Like the early windows, their
  confirmed increments are scored against the modelled laboratory volume
  with the pooled positivity. Their day grid carries a `late_start` day
  (the last laboratory day) so the model bins the modelled volume over
  each late window's *own* day range, pinned at `late_start`, rather
  than from day 0 (which would double-count the observed window volume).
  A late day that publishes a 24h analysed count (`lab_daily_history`,
  the post-cutoff daily denominators the unanchored windows otherwise lack) is
  flagged in `late_analysed` so the model can score it as a Binomial of
  that observed denominator instead, anchoring its positivity to data.

The three groups partition the confirmed counts at the first and last
laboratory dates, so no confirmed case is counted twice. Returns
`(; obs_days, obs_positives, obs_analysed, early_days, early_increments,
late_days, late_increments, late_analysed, late_start)` of grid
day-indices and per-window counts; `late_analysed[i]` is the observed 24h
analysed denominator for late day `i` (0 when none was published). The
observed and late groups are empty when no laboratory history is present
and every confirmed vintage becomes an early window. Pure integer
bookkeeping on the observed data, so it carries no gradient.
"""
function confirmed_positivity_windows(confirmed_history, lab_history,
        lab_daily_history = (; days = Int[], counts = Int[]))
    empty = (; obs_days = Int[], obs_positives = Int[], obs_analysed = Int[],
        early_days = Int[], early_increments = Int[], early_start = 0,
        late_days = Int[], late_increments = Int[], late_analysed = Int[],
        late_start = 0)
    isempty(confirmed_history.counts) && return empty
    cdays = confirmed_history.days
    ccounts = confirmed_history.counts
    ## Cumulative confirmed at (or most recently before) a grid day.
    function confirmed_at(day)
        i = searchsortedlast(cdays, day)
        return i == 0 ? 0 : Int(ccounts[i])
    end

    ## The first confirmed vintage is the BASELINE (the initial cumulative
    ## level the surveillance system was at when reporting began); it is not
    ## scored. The vintaging "starts with the data" from there, so no window's
    ## modelled volume rolls over the pre-surveillance cryptic phase. This
    ## matches how the observed and late laboratory windows already treat their
    ## first value as a baseline. `early_start` is that first confirmed day; the
    ## early windows pin their modelled-volume accumulation at it.
    early_start = Int(cdays[1])

    ## No laboratory denominators: every confirmed vintage AFTER the baseline is
    ## an early window scored through the modelled laboratory volume, binned
    ## from `early_start`.
    if isempty(lab_history.counts)
        ed = Int[Int(d) for d in cdays[2:end]]
        inc = Int[Int(ccounts[i]) - Int(ccounts[i - 1]) for i in 2:length(ccounts)]
        return (; obs_days = Int[], obs_positives = Int[],
            obs_analysed = Int[], early_days = ed, early_increments = inc,
            early_start, late_days = Int[], late_increments = Int[],
            late_analysed = Int[], late_start = 0)
    end

    first_lab_day = Int(lab_history.days[1])
    ## Early windows: confirmed vintages AFTER the first confirmed (baseline) up
    ## to and including the first laboratory date, scored against the modelled
    ## laboratory volume binned over each window's own range from `early_start`.
    early_days = Int[]
    early_increments = Int[]
    prev = Int(ccounts[1])
    for (i, d) in enumerate(cdays)
        i == 1 && continue
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

    ## Late windows: confirmed vintages strictly after the last laboratory
    ## date, which carry no analysed-specimen denominator. Their increments
    ## are scored against the modelled laboratory volume binned over each
    ## window's own day range, so the running edge must start at the last
    ## laboratory day (`late_start`) — the `bin_increments` running `prev`
    ## starts at 0, so the model pins it at `late_start` to avoid
    ## re-counting the observed-window volume already scored above.
    last_lab_day = Int(lab_history.days[end])
    late_days = Int[]
    late_increments = Int[]
    prev_conf_late = confirmed_at(last_lab_day)
    for (i, d) in enumerate(cdays)
        Int(d) <= last_lab_day && continue
        push!(late_days, Int(d))
        push!(late_increments, Int(ccounts[i]) - prev_conf_late)
        prev_conf_late = Int(ccounts[i])
    end

    ## Observed 24h analysed denominators on the late days: a published
    ## daily analysed count (`lab_daily_history`) on a late day becomes that
    ## window's Binomial denominator, 0 leaves it a modelled-volume
    ## window. Only days strictly after the last cumulative laboratory date
    ## are anchored (the cumulative series already covers the rest). Built
    ## with a plain loop (no Dict/filtered generator) so the AD backend can
    ## compile a rule for this pure-integer bookkeeping on every Julia
    ## version.
    late_analysed = zeros(Int, length(late_days))
    for (d, c) in zip(lab_daily_history.days, lab_daily_history.counts)
        di = Int(d)
        di > last_lab_day || continue
        pos = findfirst(==(di), late_days)
        pos === nothing || (late_analysed[pos] = Int(c))
    end

    return (; obs_days, obs_positives, obs_analysed, early_days,
        early_increments, early_start, late_days, late_increments,
        late_analysed, late_start = last_lab_day)
end

"""
Laboratory pipeline likelihood. Two streams driven by the shared renewal
onsets and the suspected-case pipeline from [`reported_cases_model`](@ref):

- Specimens analysed. The suspected daily pipeline (`p_drc`-scaled BVD
  onset-to-report signal plus the non-BVD background `λ_bg`) is carried
  through a sampled report-to-analysed delay and thinned by the tested
  fraction `τ_test`, giving the expected daily analysed-specimen volume.
  Its between-vintage increments are fitted against `lab_history`
  (specimens analysed) as observed `~` data with a NegativeBinomial
  sharing the surveillance dispersion `k`, identifying `τ_test` and the
  delay. This single suspected->analysed volume is also the denominator
  proxy reused in the early and unanchored late windows below, so the volume
  that is fitted is the same quantity used as the denominator where no
  analysed count is observed. The specimens-received series is not
  modelled: analysed is the throughput that actually produces confirmed
  cases, and received overshoots it by the laboratory backlog.
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

The per-window positivity has two links, set by `positivity_link`. The
default `:composition` ties the tested BVD share to the suspect-pool
composition `φ_v = (p_drc·BVD)_v / ((p_drc·BVD)_v + λ_bg_v)`, severity-
upsampled by a decaying enrichment δ0, then mapped to the tested-positive
probability through the assay sensitivity and specificity,
`p = s · q + (1 − spec)(1 − q)`. The false-positive term `(1 − spec)(1 − q)`
makes the confirmed counts respond to the non-BVD share `1 − q`, so the
laboratory positivity identifies the background `λ_bg` rather than
absorbing it into a free curve — the structural link that lets the lab data
pin `λ_bg`. The alternative `:free` link uses a free partially-pooled
per-window random effect ([`confirmed_positivity_model`](@ref)) decoupled
from `λ_bg`; it leaves `λ_bg` weakly identified and is kept for sensitivity
analysis.

The observed-window positives are conditioned on the observed analysed
denominator, so the Binomial conditioning that removes the multiplicative
ascertainment ridge is preserved; only the early and unanchored late windows use the
modelled analysed volume as the denominator.

The tested fraction `τ_test` and background rate `λ_bg` come from
[`reported_cases_model`](@ref) so the suspected and laboratory streams
share them. Exposes the per-window positivity, the expected analysed and
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
        lab_daily_history = (; days = Int[], counts = Int[]),
        tests_analysed::Union{Missing, Integer} = missing,
        receipt = lab_delay_model(),
        positivity = confirmed_positivity_model,
        positivity_link::Symbol = :composition,
        severity_enrichment = severity_enrichment_model(),
        sensitivity = test_sensitivity_model(),
        specificity = test_specificity_model(),
        ## When false, the early/late windows (confirmed vintages with NO
        ## observed analysed denominator) are not scored — only the
        ## observed-denominator Binomial windows contribute, so confirmed
        ## informs positivity without extrapolating a denominator from
        ## incidence. Used to probe the no-test-data extrapolation.
        fit_unanchored::Bool = true)
    n = length(onsets)
    ## `missing` cut-off scalar means generator mode: observed increments are
    ## left missing so `predict` resamples them.
    have_data = !ismissing(confirmed_cases)

    ## Laboratory capacity onset. No specimens are analysed before testing
    ## existed, so the modelled analysed volume is gated to zero before the
    ## first confirmed-case vintage (the earliest evidence of testing; the
    ## first laboratory date is the fallback). Modelling a pre-testing analysed
    ## volume would invent capacity that did not exist AND roll it into the
    ## first laboratory and early-confirmed bins, vastly over-predicting the
    ## early confirmed counts. The suspected-case pipeline feeding the volume is
    ## NOT gated — suspected cases did accumulate over the cryptic phase.
    cap_start = !isempty(confirmed_history.days) ?
                clamp(Int(confirmed_history.days[1]), 1, n) :
                (!isempty(lab_history.days) ?
                 clamp(Int(lab_history.days[1]), 1, n) : 1)

    ## Analysed-specimen volume: the suspected pipeline carried through the
    ## report-to-analysed delay and thinned by the tested fraction, fit to the
    ## analysed series and reused as the denominator in the early and unanchored late
    ## windows below. `bg_daily` is the per-day non-BVD background.
    receipt_state ~ to_submodel(receipt)
    suspected_daily = p_drc .* bvd_reports_daily .+ bg_daily
    analysed_daily = τ_test .* convolve_delay(suspected_daily,
        receipt_state.pmf)
    ## In predict mode (no AD) the daily series can infer as `Vector{Any}`
    ## on some Julia versions, which then makes `reduce_empty` / `zero(Any)`
    ## fail on the empty derived window vectors below. Concretise to the
    ## working scalar type; this runs only when the element type has widened,
    ## so the AD/fit path (concrete eltype) is left untouched.
    if eltype(analysed_daily) === Any
        analysed_daily = convert(Vector{typeof(τ_test)}, analysed_daily)
    end
    ## Gate the capacity to the testing window: zero before `cap_start`.
    analysed_daily = gate_before(analysed_daily, cap_start)
    rvobs = vintage_obs(lab_history, tests_analysed, n)
    analysed_inc = bin_increments(analysed_daily, rvobs.days)
    ## Generator mode leaves the volume increments missing so `predict`
    ## resamples them, like the early/late windows below.
    vol_obs = have_data ? rvobs.obs_increments : missing
    analysed_increments ~ to_submodel(
        vintage_increments_model(analysed_inc, vol_obs, k))

    ## Post-cutoff 24h analysed volume. After the national cumulative analysed
    ## series stops, INSP publishes a 24h analysed count on some days; the
    ## modelled daily analysed volume on each such day is scored against that
    ## count, so the post-cutoff testing throughput is fitted from the same
    ## stream rather than only used as a confirmed denominator. Same
    ## `have_data` gate so `predict` resamples it.
    daily_days = [clamp(Int(d), 1, n) for d in lab_daily_history.days]
    daily_modelled = isempty(daily_days) ? similar(analysed_daily, 0) :
                     [analysed_daily[d] for d in daily_days]
    daily_obs = have_data ? lab_daily_history.counts : missing
    analysed_daily_increments ~ to_submodel(
        vintage_increments_model(daily_modelled, daily_obs, k))

    ## Confirmed positives in three groups sharing one partially-pooled
    ## positivity: early windows (before the first lab date, no observed
    ## analysed) scored as counts against the modelled laboratory volume,
    ## observed windows scored as a Binomial of the observed analysed
    ## denominator, and late windows (after the last lab date, INSP's
    ## confirmed-only format) scored as counts against the modelled volume
    ## like the early windows.
    windows = confirmed_positivity_windows(confirmed_history, lab_history,
        lab_daily_history)
    n_early = length(windows.early_days)
    n_obs = length(windows.obs_analysed)
    n_late = length(windows.late_days)
    nv = n_early + n_obs + n_late

    ## Per-window tested BVD share `p_pos`. Two links:
    ## `:free` — a free partially-pooled per-window random effect
    ## ([`confirmed_positivity_model`](@ref)), decoupled from `λ_bg`.
    ## `:composition` (default) — the tested share is the suspect-pool
    ## `φ_v = (p_drc·BVD)_v / ((p_drc·BVD)_v + λ_bg_v)` over each laboratory
    ## window, upsampled by a decaying severity enrichment δ0 (see
    ## [`severity_enrichment_model`](@ref)), so the lab positivity identifies
    ## the background `λ_bg` rather than absorbing it into a free curve.
    window_days = vcat(windows.early_days, windows.obs_days,
        windows.late_days)
    if positivity_link === :composition
        enrich_state ~ to_submodel(severity_enrichment, false)
        δ0 = enrich_state.δ0
        decay_scale = enrich_state.decay_scale
        ## PCR sensitivity and specificity. The tested-positive probability
        ## is `p = s · q + (1 − spec)(1 − q)` with `q` the tested BVD share:
        ## the false-positive term `(1 − spec)(1 − q)` makes the confirmed
        ## counts respond to the non-BVD share `1 − q`, so the laboratory
        ## data identify the background `λ_bg` through the composition `φ`
        ## rather than the BVD signal alone. Without it the confirmed
        ## positivity tracks only `q`, leaving `λ_bg` weakly identified.
        sens_state ~ to_submodel(sensitivity, false)
        spec_state ~ to_submodel(specificity, false)
        s_test = sens_state.s_test
        spec = spec_state.spec
        ## Suspect-pool composition over each window, carried through the
        ## report-to-analysed delay so it reflects the composition of the
        ## specimens actually analysed in the window, consistent with the
        ## modelled volume `analysed_daily`. The `τ_test` factor cancels in the
        ## ratio φ, so it is omitted here.
        analysed_bvd_daily = convolve_delay(p_drc .* bvd_reports_daily,
            receipt_state.pmf)
        analysed_bg_daily = convolve_delay(bg_daily, receipt_state.pmf)
        if eltype(analysed_bvd_daily) === Any
            analysed_bvd_daily = convert(Vector{typeof(τ_test)},
                analysed_bvd_daily)
            analysed_bg_daily = convert(Vector{typeof(τ_test)},
                analysed_bg_daily)
        end
        ## Gate the tested composition to the testing window too, so the
        ## composition clock and the per-window BVD share start at the testing
        ## onset rather than rolling the cryptic phase.
        analysed_bvd_daily = gate_before(analysed_bvd_daily, cap_start)
        analysed_bg_daily = gate_before(analysed_bg_daily, cap_start)
        bvd_window = bin_increments(analysed_bvd_daily, window_days)
        bg_window = bin_increments(analysed_bg_daily, window_days)
        Tt = eltype(bvd_window)
        ## Testing clock: cumulative modelled analysed volume at each window.
        vol_window = bin_increments(analysed_daily, window_days)
        c_window = cumsum(vol_window)
        lo = convert(Tt, 1e-8)
        hi = one(Tt) - lo
        ## Floor the decay scale so a near-zero `decay_scale` draw cannot make
        ## the clock ratio `0/0` (NaN) and break the downstream Binomial.
        dscale = max(convert(Tt, decay_scale), one(Tt))
        s_t = convert(Tt, s_test)
        sp_t = convert(Tt, spec)
        p_pos = map(eachindex(window_days)) do i
            ## Pool composition φ = (p_drc·BVD) / ((p_drc·BVD) + λ_bg) over the
            ## window, guarded against a zero/negative denominator.
            num = bvd_window[i]
            den = bvd_window[i] + bg_window[i]
            ratio = num / (den + lo)
            φ = clamp(isfinite(ratio) ? ratio : convert(Tt, 0.5), lo, hi)
            δ_i = convert(Tt, δ0) * exp(-c_window[i] / dscale)
            ## Severity-enriched tested BVD share, then the assay
            ## sensitivity/specificity transform to the tested-positive
            ## probability so the false-positive term identifies `λ_bg`.
            q = logistic(logit(φ) + δ_i)
            qf = isfinite(q) ? q : φ
            qe = clamp(qf, lo, hi)
            p = s_t * qe + (one(Tt) - sp_t) * (one(Tt) - qe)
            ## Final guard: clamp into (0,1) and replace any non-finite value
            ## with the composition so the confirmed Binomial always sees a
            ## valid probability even under an AD perturbation.
            clamp(isfinite(p) ? p : φ, lo, hi)
        end
    else
        pos_state ~ to_submodel(positivity(nv))
        p_pos = pos_state.p_pos
    end

    ## Early windows: confirmed increment ~ NegBinomial(positivity ×
    ## modelled analysed volume), the volume binned over each window's OWN day
    ## range pinned at `early_start` (the first confirmed vintage, the testing-
    ## onset baseline), so the first early increment is scored from the data
    ## start rather than rolling the (now-gated) pre-testing volume. Mirrors the
    ## late-window pinning at `late_start`.
    early_p = p_pos[1:n_early]
    early_volume = n_early > 0 ?
                   bin_increments(analysed_daily,
        vcat(windows.early_start, windows.early_days))[2:end] :
                   similar(analysed_daily, 0)
    early_mean = early_p .* early_volume
    early_obs = (have_data && n_early > 0 && fit_unanchored) ?
                windows.early_increments : missing
    early_increments ~ to_submodel(
        vintage_increments_model(early_mean, early_obs, k))

    ## Observed windows: Binomial of the observed analysed denominator.
    obs_p = p_pos[(n_early + 1):(n_early + n_obs)]
    obs_positives = (have_data && n_obs > 0) ? collect(windows.obs_positives) :
                    missing
    confirmed_positives ~ to_submodel(
        confirmed_positives_model(obs_positives, windows.obs_analysed, obs_p))

    ## Late windows: confirmed-only vintages after the last laboratory date.
    ## A day that publishes a 24h analysed count (`late_analysed > 0`) is
    ## scored as a Binomial of that observed denominator — like an observed
    ## window, anchoring its positivity to data — and each remaining unanchored day
    ## as NegBinomial(positivity × modelled volume). The modelled volume is
    ## binned over each late window's own day range, with the running edge
    ## PINNED at the last laboratory day (`late_start`): `bin_increments`
    ## runs its running `prev` from day 0, so prepending `late_start` to the
    ## late day edges and dropping the synthetic first bin starts the
    ## accumulation at `late_start`, avoiding double-counting the
    ## observed-window volume.
    late_p = p_pos[(n_early + n_obs + 1):nv]
    if n_late > 0
        late_edges = vcat(windows.late_start, windows.late_days)
        late_volume = bin_increments(analysed_daily, late_edges)[2:end]
    else
        late_volume = similar(analysed_daily, 0)
    end
    late_mean = late_p .* late_volume
    ## Observed late increments: anchored days (24h denominator) carry the
    ## confirmed increment clamped into the Binomial support and are always
    ## scored; unanchored days are scored only when `fit_unanchored` (the
    ## no-extrapolation probe leaves them latent). A per-entry
    ## `missing`/value vector lets the one submodel observe each accordingly.
    if have_data && n_late > 0
        late_obs = Vector{Union{Missing, Int}}(undef, n_late)
        for i in 1:n_late
            a = windows.late_analysed[i]
            if a > 0
                late_obs[i] = clamp(windows.late_increments[i], 0, a)
            elseif fit_unanchored
                late_obs[i] = windows.late_increments[i]
            else
                late_obs[i] = missing
            end
        end
    else
        late_obs = missing
    end
    late_increments ~ to_submodel(
        late_confirmed_model(late_obs, late_mean, windows.late_analysed,
        late_p, k))

    expected_analysed := safe_rate(sum(analysed_daily))
    ## Expected confirmed at the cut-off and the overall positivity, over the
    ## modelled early volume, the observed cumulative analysed windows and the
    ## late windows (anchored days contribute `p · analysed`, unanchored days the
    ## modelled `p · volume`). The window vectors are empty when a vintage has
    ## no such window, and in predict mode their element type can widen to
    ## `Any`, so each sum is given a concrete `init` to skip `reduce_empty`'s
    ## `zero(Any)`. The init is taken from the scalar `τ_test` (always
    ## concrete), NOT from `eltype(p_pos)`, which can widen to `Any`.
    z = zero(τ_test)
    amask = windows.late_analysed .> 0
    late_den_a = float.(windows.late_analysed)
    late_den = n_late > 0 ? ifelse.(amask, late_den_a, late_volume) :
               similar(late_volume, 0)
    late_expected = n_late > 0 ? ifelse.(amask, late_p .* late_den_a, late_mean) :
                    similar(late_mean, 0)
    denom = sum(early_volume; init = z) + float(sum(windows.obs_analysed)) +
            sum(late_den; init = z)
    expected_positives = sum(early_mean; init = z) +
                         sum(late_expected; init = z) +
                         (n_obs > 0 ? sum(obs_p .* windows.obs_analysed) : z)
    expected_confirmed := safe_rate(expected_positives)
    p_positive := safe_rate(expected_positives) / safe_rate(denom)

    ## Modelled daily confirmed-case incidence: the per-window tested-positive
    ## probability expanded onto the daily grid times the modelled analysed
    ## volume. Exposed so the composer's cumulative-confirmed trajectory and a
    ## survivors-among-confirmed (recovered) stream can reuse one consistent
    ## daily series. In predict / check-model mode `p_pos` can widen to
    ## `Vector{Any}`, so pin it to the analysed volume's (always-concrete)
    ## element type before expanding, matching the other guards here.
    p_pos_daily = p_pos
    if eltype(p_pos_daily) === Any
        p_pos_daily = convert(Vector{eltype(analysed_daily)}, p_pos_daily)
    end
    confirmed_daily = expand_vintage_rate(p_pos_daily, window_days, n) .*
                      analysed_daily

    return (; τ_test, bg_daily, p_pos, windows, analysed_daily,
        confirmed_daily,
        receipt_pmf = receipt_state.pmf,
        expected_analysed, expected_confirmed, p_positive)
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
mean at-risk dwell time. Samples the traveller volume and the
onset-to-detection delay via injected submodels. The onset-to-detection
prior is centred on the Ebola onset-to-hospitalisation delay (mean 5.0 d,
SD 4.7 d; WHO Ebola Response Team 2014, NEJM), the delay from symptom
onset to detection at a point of entry abroad.

## Dated per-day likelihood

The observed Uganda imports are a dated series — `export_case_days` gives
the grid day-index of each detection (sorted ascending) — fitted as an
inhomogeneous Poisson process rather than a single cumulative count. The
cumulative export intensity is `Λ(t) = sum(export_prevalence[1:t])`, so
the per-day expected export increment `Λ(t) − Λ(t−1)` is just the day-`t`
at-risk person-time `export_prevalence[t]`. The likelihood places one
Poisson term per import at its detection day, the increment between
consecutive detection-day edges (via [`bin_increments`](@ref) on the daily
prevalence series), with two extra terms:

  - `pre_detection_exports ~ Poisson(Λ(d₁−1))` observed at zero, the
    first-detection timing bound: no export is expected before the
    earliest detection day `d₁`. The first import's increment is measured
    from `Λ(d₁−1)`, so the pre-detection weight and the import increments
    partition `Λ(t_last)` and the model still conditions on the same total
    as a cumulative single-Poisson would.
  - `last_offset` truncation: the travel-gated export clock stops at the
    last observed import `t_last = day of the most recent detection`. Days
    after the last import carry no informative zero (cross-border movement
    shifts over the outbreak and the most recent days are right-truncated
    by reporting lag), so prevalence past `t_last` does not accrue. With an
    empty `export_case_days` the model falls back to the cut-off cumulative
    Poisson `exported_cases ~ Poisson(Λ(T))`.

Returns the expected cumulative count at `t_last`, the per-capita travel
rate and the daily at-risk prevalence for reuse by
[`exports_deaths_model`](@ref).
"""
@model function exports_model(
        exported_cases::Union{Missing, Integer},
        infections::AbstractVector, p_uganda::Real;
        export_case_days::AbstractVector{<:Integer} = Int[],
        pre_detection_exports::Union{Missing, Integer} = 0,
        incubation_pmf::AbstractVector,
        source_population::Real = ITURI_POPULATION,
        traveller = traveller_volume_model(),
        ## Export detection abroad uses the same line-list onset→admission
        ## delay (d_oa) as the suspect-case report: a case is detected at a
        ## point of entry when first formally seen, ~4 days after onset.
        onset_to_detection = gamma_delay_model(cdf_nmax(Gamma(1.178, 3.694));
            alpha_prior = truncated(Normal(1.178, 0.285); lower = 0.01),
            theta_prior = truncated(Normal(3.694, 1.198); lower = 0.1)))
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
    n = length(export_prevalence)

    if isempty(export_case_days)
        ## No dated series: cumulative single-total Poisson at the cut-off.
        raw_exports = sum(export_prevalence)
        expected_exports_T := safe_rate(raw_exports)
        exported_cases ~ Poisson(expected_exports_T)
    else
        ## Dated per-day Poisson. The export clock stops at the last import
        ## `t_last` (the `last_offset` truncation); prevalence past it does
        ## not accrue. `d₁` is the earliest detection day.
        days, counts = dated_event_bins(export_case_days, n)
        d₁ = days[1]
        ## Pre-detection survival weight Λ(d₁−1): the cumulative export
        ## intensity up to the day before the earliest detection.
        pre = d₁ > 1 ? sum(@view export_prevalence[1:(d₁ - 1)]) :
              zero(@inbounds export_prevalence[begin])
        pre_detection_exports ~ Poisson(safe_rate(pre))
        ## Per-day-edge increments between consecutive detection days; the
        ## first is measured from the pre-detection weight `pre`, so the
        ## pre-detection term and the increments partition Λ(t_last).
        raw_inc = bin_increments(export_prevalence, days)
        μ_day = [i == 1 ? raw_inc[1] - pre : raw_inc[i]
                 for i in eachindex(raw_inc)]
        obs = ismissing(exported_cases) ? missing : counts
        export_obs ~ to_submodel(dated_poisson_model(μ_day, obs))
        ## Reported expected count is the cumulative intensity to `t_last`.
        expected_exports_T := safe_rate(pre + sum(μ_day))
    end

    ## Travel-scaled at-risk prevalence WITHOUT the export-case ascertainment
    ## `p_uganda`: a death among an exported case would be reported whether or
    ## not the case itself was ascertained as an import, so the export-death
    ## model accrues over the travelled person-time `q · prevalence`, not the
    ## ascertained `export_prevalence = p_uganda · q · prevalence`.
    travelled_prevalence = q .* prevalence
    return (; p_uganda, daily_travellers, q, prevalence,
        export_prevalence, travelled_prevalence,
        expected_exports = expected_exports_T)
end

"""
Deaths-among-exported-cases likelihood, dated per-day. Deaths accrue among
the TRAVELLED at-risk person-time `q · prevalence` from
[`exports_model`](@ref), BEFORE the export-case ascertainment `p_uganda`: a
death among an exported case would be reported whether or not the case was
ascertained as an import, so the death model is not thinned by `p_uganda`
(unlike the export-case count). The day-`t` expected export-death increment
is the CFR-scaled convolution of that travelled prevalence with the
infection→death PMF,

```math
\\mathrm{d}\\Lambda_\\text{d}(t) = \\mathrm{CFR}
    \\sum_{s\\le t} q\\,\\text{prevalence}(s)\\, f_\\text{d}(t-s),
```

so the cumulative export-death intensity `Λ_d(t)` is its running sum, the
discrete analogue of the integral model's `∫ C(s)·S_det·F_death ds`. The
onset-to-death PMF `od_pmf` shared from [`deaths_model`](@ref) is convolved
with the incubation PMF to give the infection→death distribution `f_d`,
keyed to infection like the prevalence.

The observed Uganda export deaths are a dated series `export_death_days`
(grid day-indices, ascending), fitted as an inhomogeneous Poisson exactly
as [`exports_model`](@ref) fits the imports: one Poisson term per death day
(the increment of `Λ_d` between consecutive death-day edges), a
`pre_death_exports ~ Poisson(Λ_d(δ₁−1))` zero term bounding the first death
day `δ₁`, and the `last_offset` truncation that stops the clock at the last
observed death day. With an empty `export_death_days` the model falls back
to the cut-off cumulative Poisson `exports_deaths ~ Poisson(Λ_d(n))`.
"""
@model function exports_deaths_model(
        exports_deaths::Union{Missing, Integer},
        travelled_prevalence::AbstractVector, CFR::Real,
        od_pmf::AbstractVector, incubation_pmf::AbstractVector;
        export_death_days::AbstractVector{<:Integer} = Int[],
        pre_death_exports::Union{Missing, Integer} = 0)
    n = length(travelled_prevalence)
    ## Infection→death PMF by age (age 0 = same day).
    fd_pmf = convolve_pmf(incubation_pmf, od_pmf)
    ## Per-day expected export-death increment: CFR-scaled convolution of
    ## the daily at-risk prevalence with the infection→death PMF. Its
    ## running sum is the cumulative export-death intensity `Λ_d`.
    death_daily = CFR .* convolve_delay(travelled_prevalence, fd_pmf)

    if isempty(export_death_days)
        ## No dated series: cumulative single-total Poisson at the cut-off.
        expected_exports_deaths_T := safe_rate(sum(death_daily))
        exports_deaths ~ Poisson(expected_exports_deaths_T)
    else
        ## Dated per-day Poisson; the clock stops at the last death day.
        days, counts = dated_event_bins(export_death_days, n)
        δ₁ = days[1]
        pre = δ₁ > 1 ? sum(@view death_daily[1:(δ₁ - 1)]) :
              zero(@inbounds death_daily[begin])
        pre_death_exports ~ Poisson(safe_rate(pre))
        raw_inc = bin_increments(death_daily, days)
        μ_day = [i == 1 ? raw_inc[1] - pre : raw_inc[i]
                 for i in eachindex(raw_inc)]
        obs = ismissing(exports_deaths) ? missing : counts
        death_obs ~ to_submodel(dated_poisson_model(μ_day, obs))
        expected_exports_deaths_T := safe_rate(pre + sum(μ_day))
    end

    return (; expected_exports_deaths_T)
end

"""
Laboratory-confirmed-deaths likelihood, the death analogue of the
confirmed-case laboratory pipeline ([`confirmed_cases_model`](@ref)). The
confirmed cases score a modelled analysed-specimen volume (the suspected-case
pipeline carried to laboratory receipt and thinned by the testing fraction)
times a composition-linked assay positivity. The death side has no published
analysed denominator, so the confirmed-death increments are `NegBinomial(k)`
counts of the modelled death volume, the same modelled-volume route the early
and post-lab confirmed-case windows use.

Three pieces:

  - Death analysed volume. Deaths are tested out of the same laboratory as
    cases, so the death volume tracks the modelled case analysed volume
    `case_analysed_daily` at the per-day suspected death-to-case ratio, times a
    testing-intensity scaling, `v = scaling · case_analysed_daily ·
    susp_death / susp_case`, with `susp_death` and `susp_case` the suspected
    deaths and cases carried to receipt by the confirmed cases' delay. The case
    volume carries the laboratory capacity onset, so the death volume inherits
    it. The scaling ([`death_testing_scaling_model`](@ref)) is the per-suspect
    testing-intensity difference between deaths and living suspects; the
    death-to-case ratio carries the suspect-pool severity and the
    suspected-death level. The death-only composer has no case stream and falls
    back to a death testing fraction ([`death_testing_fraction_model`](@ref)).
  - Death-pool composition. The BVD share of the suspected deaths at receipt,
    `q_death = bvd_death / (bvd_death + bg_death)` per day, from the death
    series' own BVD and background components. The death background, tied to
    the case background by `cfr_bg` (see [`deaths_model`](@ref) and
    [`background_cfr_model`](@ref)), keeps `q_death` below one.
  - Assay positivity. `p = s · q_death + (1 − spec)(1 − q_death)` with PCR
    sensitivity `s` ([`test_sensitivity_model`](@ref)) and specificity `spec`
    ([`test_specificity_model`](@ref)), the same form as the confirmed-case
    positivity, drawn from the same priors as separate death-stream parameters.

Returns the cut-off realised death testing fraction `τ_death`, the testing
scaling, the cut-off death-pool composition `q_death`, the confirmation
positivity and the expected confirmed-death count.
"""
@model function confirmed_deaths_model(
        confirmed_deaths::Union{Missing, Integer},
        total_deaths::Union{Missing, Integer},
        deaths_daily::AbstractVector,
        bvd_deaths_daily::AbstractVector,
        bg_death_daily::AbstractVector, k::Real;
        confirmed_deaths_history = (; days = Int[], counts = Int[]),
        receipt_pmf::AbstractVector = [1.0],
        capacity_start::Integer = 0,
        case_analysed_daily = nothing,
        case_suspected_daily = nothing,
        scaling = death_testing_scaling_model(),
        testing = death_testing_fraction_model(),
        sensitivity = test_sensitivity_model(),
        specificity = test_specificity_model())
    sens_state ~ to_submodel(sensitivity)
    spec_state ~ to_submodel(specificity)
    s = sens_state.s_test
    spec = spec_state.spec
    n = length(deaths_daily)

    ## Suspected deaths carried to laboratory receipt by the same
    ## report-to-receipt delay the confirmed cases use, with the BVD component.
    susp_death = convolve_delay(deaths_daily, receipt_pmf)
    bvd_death = convolve_delay(bvd_deaths_daily, receipt_pmf)
    ## In predict or check-model mode the series can widen to `Vector{Any}`,
    ## which trips `zero(Any)` downstream; pin to the sampled scalar type,
    ## leaving the fit path (concrete dual eltype) untouched.
    if eltype(susp_death) === Any
        susp_death = convert(Vector{typeof(s)}, susp_death)
        bvd_death = convert(Vector{typeof(s)}, bvd_death)
    end

    ## Death-pool BVD composition per day, q = bvd / (bvd + bg), and the assay
    ## tested-positive probability p = s·q + (1 − spec)(1 − q). The false-
    ## positive term lets the confirmed deaths respond to the non-BVD share,
    ## the same link the confirmed cases use.
    lo = eps(typeof(s))
    hi = one(s) - lo
    q_death_daily = map(eachindex(susp_death)) do t
        den = susp_death[t]
        ratio = den > lo ? bvd_death[t] / den : one(s)
        clamp(isfinite(ratio) ? ratio : one(s), lo, hi)
    end
    p_pos_daily = s .* q_death_daily .+ (one(s) - spec) .*
                                        (one(s) .- q_death_daily)

    ## Death analysed volume. Deaths are tested out of the same laboratory as
    ## cases, so the death volume tracks the modelled case analysed volume at
    ## the per-day suspected death-to-case ratio, times a testing-intensity
    ## scaling, v = scaling · analysed_case · susp_death / susp_case. The case
    ## volume already carries the laboratory capacity onset, so the death volume
    ## inherits it. The death-only composer has no case stream and falls back to
    ## a death testing fraction of the suspected deaths, gated at the onset.
    if case_analysed_daily !== nothing
        scale_state ~ to_submodel(scaling)
        sc = scale_state.scaling
        susp_case = convolve_delay(case_suspected_daily, receipt_pmf)
        death_volume = map(eachindex(susp_death)) do t
            den = susp_case[t]
            den > lo ? sc * case_analysed_daily[t] * susp_death[t] / den :
            zero(sc)
        end
        τ_death = susp_death[n] > lo ? death_volume[n] / susp_death[n] : zero(sc)
    else
        test_state ~ to_submodel(testing)
        τ_death = test_state.τ_death
        sc = one(τ_death)
        death_volume = τ_death .* gate_before(susp_death, capacity_start)
    end

    confirmed_death_daily = p_pos_daily .* death_volume
    vobs = vintage_obs(confirmed_deaths_history, confirmed_deaths, n)
    modelled_inc = bin_increments(confirmed_death_daily, vobs.days)
    cdeath_increments ~ to_submodel(
        vintage_increments_model(modelled_inc, vobs.obs_increments, k))

    expected_confirmed_deaths := safe_rate(sum(confirmed_death_daily))
    ## Cut-off death-pool composition and confirmation positivity, surfaced as
    ## `death_composition` and `death_confirmation`.
    q_death := q_death_daily[n]
    p_death_conf := p_pos_daily[n]

    return (; τ_death, scaling = sc, s_test = s, spec, q_death, p_death_conf,
        expected_confirmed_deaths)
end

"""
DRC treatment-centre patient-flow likelihood: the supply-limited
isolation/treatment-bed occupancy stream together with the in-care outcome
flows. The "Patients en isolement" figure is the daily count of occupied
beds, modelled as a latent bed demand right-censored at the bed capacity (bed
demand can outstrip supply, with occupancy catching up as capacity expands).

Admissions are the reported suspects ([`reported_cases_model`](@ref)) carried
through a short suspected→admission delay (triage, transport, bed-wait) and
split into a BVD/background mixture: BVD demand `p_iso_bvd · p_drc ·
bvd_reports` admitted at the severity-skewed rate `p_iso_bvd`
([`isolation_severity_model`](@ref)), and non-BVD demand `p_iso · bg_daily` at
the base rate ([`isolation_admission_model`](@ref)). The bed length-of-stay is
an outcome mixture — not competing risks: an admission is assigned by a
probability to an outcome and then leaves on that outcome's length-of-stay
distribution. The BVD bed stay is the mixture of an admission→death stay
(weight `CFR_iso`, the in-care fatality) and a longer admission→recovery stay
(weight `1 − CFR_iso`); the non-BVD stay is a rule-out stay. Occupancy is the
admissions carried through the corresponding survival
([`convolve_survival`](@ref)) and summed.

The in-care fatality `CFR_iso = 1 − p_recover` is the infection CFR adjusted
for the admitted population by a sampled log-odds modifier `β_iso`
([`recovery_probability_model`](@ref)). It weights the death versus recovery
branch and is reported alongside `β_iso`; it is a conditional-on-admission
(in-care) fatality, not a causal treatment effect.

The occupancy likelihood is a NegativeBinomial around the latent demand,
right-censored at the fixed implied-capacity bound ([`censoring_cap`](@ref),
[`censored_occupancy_model`](@ref)) — a fixed data bound, not a sampled
ceiling, so the censored likelihood has no moving `-Inf` wall to diverge
against. The capacity walk `C(t)` ([`bed_capacity_walk_model`](@ref)) carries
the implied-capacity likelihood and the forecast cap, and `occupancy =
min(demand, C)` is the supply-capped stock the derived quantities report.

The daily in-care outcome flows are scored against the Tableau 6 patient-
movement counts where available, each an optional NegativeBinomial stream that
is a no-op when its history is empty: in-care deaths (`CFR_iso ·` BVD
admissions through the death stay), rule-out discharges (non-BVD admissions
through the rule-out stay), admissions (the BVD + non-BVD inflow) and
absconded (a small fraction of occupancy). The recovered-among-confirmed
stream ("cumul guéris") is the recovery branch (`(1 − CFR_iso) ·` BVD
admissions through the recovery stay), its cumulative increments scored as
observed data — the Tableau 6 daily recoveries equal these increments. The
longer occupancy / recovered / capacity stock streams carry the model over
the window before the daily flows begin.

Exposes the cut-off occupancy, bed demand, their difference (the bed
shortfall), the utilisation, the BVD share of demand, the in-care fatality
`CFR_iso` and modifier `β_iso`, the overall length-of-stay (the outcome-
mixture mean) and the daily series for forecasting and replication.
"""
@model function treatment_flow_model(
        isolation_history,
        bvd_reports_daily::AbstractVector,
        bg_daily::AbstractVector,
        p_drc::Real,
        CFR::Real;
        capacity_history = (; days = Int[], counts = Int[]),
        admissions_history = (; days = Int[], counts = Int[]),
        deaths_history = (; days = Int[], counts = Int[]),
        ruleout_history = (; days = Int[], counts = Int[]),
        absconded_history = (; days = Int[], counts = Int[]),
        admission = isolation_admission_model(),
        severity = isolation_severity_model(),
        capacity = bed_capacity_walk_model,
        dispersion = surveillance_dispersion_model(),
        ## Occupancy / flow dispersion can be injected from the joint composer's
        ## pooled set (`k_external`); standalone it samples its own.
        k_external::Union{Nothing, Real} = nothing,
        ## In-care fatality modifier prior: β_iso on the infection CFR.
        cfr_modifier_prior = Normal(0.0, 0.5),
        ## Small abscond / loss-to-follow-up fraction of occupancy per day.
        abscond_prior = truncated(Normal(0.01, 0.01); lower = 0),
        ## Short suspected→admission delay (report → reaching a bed: triage,
        ## transport, bed-wait), distinct from the report→lab receipt delay.
        admission_delay = censored_delay_model(
            cdf_nmax(lognormal_meansd(2.0, 1.5); q = 0.99);
            mean_prior = truncated(Normal(2.0, 1.0); lower = 0.1),
            sd_prior = truncated(Normal(1.5, 1.0); lower = 0.3)),
        ## Outcome-mixture BVD bed stay: admission→death (the admission→death
        ## atomic delay the onset→death convolution also uses, mean ≈ 8.4 d) and
        ## the longer admission→recovery stay (mean ≈ 14 d). Built to a common
        ## nmax so the two PMFs align for the elementwise mixture.
        death_los = gamma_delay_model(
            cdf_nmax(lognormal_meansd(14.0, 8.0); q = 0.99);
            alpha_prior = truncated(Normal(2.151, 0.604); lower = 0.01),
            theta_prior = truncated(Normal(3.906, 1.381); lower = 0.1)),
        recovery_los = censored_delay_model(
            cdf_nmax(lognormal_meansd(14.0, 8.0); q = 0.99);
            mean_prior = truncated(Normal(14.0, 5.0); lower = 1),
            sd_prior = truncated(Normal(8.0, 4.0); lower = 1)),
        ## Non-BVD rule-out stay (report→receipt turnaround plus sign-off).
        ruleout_los = censored_delay_model(
            cdf_nmax(lognormal_meansd(4.5, 4.0); q = 0.99);
            mean_prior = truncated(Normal(4.5, 2.0); lower = 1),
            sd_prior = truncated(Normal(4.0, 1.5); lower = 1)))
    adm_state ~ to_submodel(admission)
    p_iso = adm_state.p_iso
    sev_state ~ to_submodel(severity)
    ## BVD suspects are admitted at a higher rate than non-BVD rule-outs,
    ## skewed up from `p_iso` by the severity log-odds `δ_iso`.
    p_iso_bvd = logistic(logit(p_iso) + sev_state.δ_iso)
    if k_external === nothing
        disp_state ~ to_submodel(dispersion)
        k = disp_state.k
    else
        k = k_external
    end
    n = length(bvd_reports_daily)
    ## In-care fatality CFR_iso = logistic(logit(CFR) + β_iso): a log-odds
    ## modifier on the infection CFR, identified by the in-care death flow
    ## (Tableau 6 décédés) relative to admissions and occupancy. β_iso < 0 means
    ## treatment lowers the in-care fatality below the infection CFR. This is a
    ## conditional-on-admission (in-care) fatality, not a causal treatment
    ## effect. The recovered-among-confirmed ("cumul guéris") stream is modelled
    ## separately off the confirmed cases ([`recovered_model`](@ref)) — guéris is
    ## the confirmed-and-discharged subset, not all in-care recoveries.
    β_iso ~ cfr_modifier_prior
    CFR_iso = logistic(logit(CFR) + β_iso)
    ## Time-varying bed capacity `C(t)` (a random walk), started at the first
    ## day with occupancy or capacity data.
    cap_obs_days = vcat(Int.(isolation_history.days),
        Int.(capacity_history.days))
    cap_start = isempty(cap_obs_days) ? 1 : minimum(cap_obs_days)
    cap_state ~ to_submodel(capacity(n; start = cap_start))
    C = cap_state.C
    C_T = isempty(C) ? zero(eltype(C)) : C[end]
    adm_delay_state ~ to_submodel(admission_delay)
    death_los_state ~ to_submodel(death_los)
    recovery_los_state ~ to_submodel(recovery_los)
    ruleout_los_state ~ to_submodel(ruleout_los)
    abscond_frac ~ abscond_prior

    ## Admission inflow carried through the short suspected→admission delay,
    ## split into BVD (admitted at `p_iso_bvd`) and non-BVD (`p_iso`) demand.
    bvd_adm = convolve_delay(p_iso_bvd .* p_drc .* bvd_reports_daily,
        adm_delay_state.pmf)
    bg_adm = convolve_delay(p_iso .* bg_daily, adm_delay_state.pmf)

    ## Outcome-mixture (not competing risks) BVD bed stay: a fraction `CFR_iso`
    ## of BVD admissions leave on the admission→death stay, the rest on the
    ## admission→recovery stay. The two PMFs share an nmax, so the mixture is
    ## elementwise.
    dpmf = death_los_state.pmf
    rpmf = recovery_los_state.pmf
    f_bvd = CFR_iso .* dpmf .+ (one(CFR_iso) - CFR_iso) .* rpmf

    ## Latent bed demand: admissions carried through the survival of their
    ## length-of-stay (occupancy = admissions still in a bed).
    bvd_demand = convolve_survival(bvd_adm, f_bvd)
    bg_demand = convolve_survival(bg_adm, ruleout_los_state.pmf)
    demand = bvd_demand .+ bg_demand
    ## Predict / check-model mode can widen the series to `Vector{Any}`, which
    ## `min.(demand, C)` cannot broadcast over; pin to the capacity's concrete
    ## element type, leaving the AD/fit path untouched.
    if eltype(demand) === Any
        demand = convert(Vector{eltype(C)}, demand)
    end
    occupancy = min.(demand, C)

    ## Occupancy likelihood: each day's bed count is a NegativeBinomial around
    ## the latent demand, right-censored at the fixed implied-capacity bound.
    ## Day-to-day reporting noise is absorbed by the NegativeBinomial
    ## dispersion rather than by fitted reclassification steps.
    iso_days = isolation_history.days
    iso_obs = isempty(isolation_history.counts) ? missing :
              collect(Int.(isolation_history.counts))
    iso_demand = [demand[clamp(Int(d), 1, n)] for d in iso_days]
    iso_means = iso_demand
    iso_ceil = censoring_cap(iso_days, iso_obs, capacity_history)
    isolation ~ to_submodel(
        censored_occupancy_model(iso_means, iso_ceil, iso_obs, k))

    ## Capacity likelihood: the implied bed count is a noisy observation of C(t).
    cap_days = capacity_history.days
    cap_modelled = [C[clamp(Int(d), 1, n)] for d in cap_days]
    cap_obs = isempty(capacity_history.counts) ? missing :
              collect(Int.(capacity_history.counts))
    bed_capacity ~ to_submodel(
        vintage_increments_model(cap_modelled, cap_obs, k))

    ## Modelled daily in-care outcome flows. The death and recovery flows sum
    ## to the total BVD discharge (consistent with the mixture survival).
    deaths_daily = convolve_delay(CFR_iso .* bvd_adm, dpmf)
    recover_daily = convolve_delay((one(CFR_iso) - CFR_iso) .* bvd_adm, rpmf)
    ruleout_daily = convolve_delay(bg_adm, ruleout_los_state.pmf)
    admit_daily = bvd_adm .+ bg_adm
    abscond_daily = abscond_frac .* occupancy

    ## Optional daily Tableau 6 flow likelihoods. Each is a no-op when its
    ## history is empty (no days → no scored terms).
    dth_days = deaths_history.days
    dth_obs = isempty(deaths_history.counts) ? missing :
              collect(Int.(deaths_history.counts))
    incare_deaths ~ to_submodel(vintage_increments_model(
        [deaths_daily[clamp(Int(d), 1, n)] for d in dth_days], dth_obs, k))
    ro_days = ruleout_history.days
    ro_obs = isempty(ruleout_history.counts) ? missing :
             collect(Int.(ruleout_history.counts))
    ruleouts ~ to_submodel(vintage_increments_model(
        [ruleout_daily[clamp(Int(d), 1, n)] for d in ro_days], ro_obs, k))
    adm_h_days = admissions_history.days
    adm_h_obs = isempty(admissions_history.counts) ? missing :
                collect(Int.(admissions_history.counts))
    admissions ~ to_submodel(vintage_increments_model(
        [admit_daily[clamp(Int(d), 1, n)] for d in adm_h_days], adm_h_obs, k))
    ab_days = absconded_history.days
    ab_obs = isempty(absconded_history.counts) ? missing :
             collect(Int.(absconded_history.counts))
    absconded ~ to_submodel(vintage_increments_model(
        [abscond_daily[clamp(Int(d), 1, n)] for d in ab_days], ab_obs, k))

    ## Cut-off reported quantities.
    z0 = zero(eltype(C))
    occ_T = isempty(occupancy) ? z0 : occupancy[end]
    dem_T = isempty(demand) ? z0 : demand[end]
    bvd_dem_T = isempty(bvd_demand) ? z0 : bvd_demand[end]
    overall_los = CFR_iso * death_los_state.mean +
                  (one(CFR_iso) - CFR_iso) * recovery_los_state.mean
    expected_isolation := safe_rate(occ_T)
    expected_bed_demand := safe_rate(dem_T)
    ## Cut-off daily flows: the end-of-grid value of each modelled daily
    ## series, the one-week-ahead forecast base for admissions, in-care
    ## deaths and rule-outs.
    expected_admissions := safe_rate(isempty(admit_daily) ? z0 : admit_daily[end])
    expected_incare_deaths := safe_rate(isempty(deaths_daily) ? z0 :
                                        deaths_daily[end])
    expected_ruleouts := safe_rate(isempty(ruleout_daily) ? z0 :
                                   ruleout_daily[end])
    bed_shortfall := safe_rate(max(dem_T - occ_T, z0))
    bed_utilisation := safe_rate(occ_T) / safe_rate(C_T)
    isolation_bvd_share := safe_rate(bvd_dem_T) / safe_rate(dem_T)
    isolation_severity := sev_state.δ_iso
    isolation_bvd_admission := p_iso_bvd
    incare_cfr := CFR_iso
    incare_cfr_modifier := β_iso
    treatment_overall_los := overall_los

    return (; p_iso, p_iso_bvd, δ_iso = sev_state.δ_iso,
        CFR_iso, β_iso, capacity = C_T,
        death_los_mean = death_los_state.mean,
        recovery_los_mean = recovery_los_state.mean,
        ruleout_los_mean = ruleout_los_state.mean,
        admission_delay_mean = adm_delay_state.mean,
        overall_los, abscond_frac, k_isolation = k,
        demand, occupancy, isolation,
        deaths_daily, recover_daily, ruleout_daily, admit_daily,
        expected_isolation = safe_rate(occ_T),
        expected_bed_demand = safe_rate(dem_T),
        expected_admissions = safe_rate(isempty(admit_daily) ? z0 :
                                        admit_daily[end]),
        expected_incare_deaths = safe_rate(isempty(deaths_daily) ? z0 :
                                           deaths_daily[end]),
        expected_ruleouts = safe_rate(isempty(ruleout_daily) ? z0 :
                                      ruleout_daily[end]))
end

"""
DRC recovered-among-confirmed likelihood ("cumul guéris"), an incidence
(scaled-convolution) stream — the renewal analogue of the convolution-and-
scaling secondary-observation model of EpiNow2 [epinow2](@cite). Recoveries
are the survivors among laboratory-confirmed cases: the modelled daily
confirmed-case incidence `confirmed_daily` (from
[`confirmed_cases_model`](@ref)) is scaled by the
recovery probability `p_recover` (the confirmed-case survival fraction,
[`recovery_probability_model`](@ref)) and convolved with a sampled
confirmation-to-recovery delay,

```math
\\text{recovered}_t = p_\\text{recover} \\sum_{s \\ge 0}
    \\text{confirmed}_{t-s}\\, f_{\\text{rec},s}.
```

The recovery fraction `p_recover` is grounded on the case-fatality ratio
`CFR` (a recovered case is one that did not die), adjusted for the confirmed
population by a sampled log-odds offset (see
[`recovery_probability_model`](@ref)). A case is taken to be confirmed
BEFORE it is recorded as recovered: the report's "cumul guéris" counts
recoveries among confirmed cases, so the recovery follows confirmation by
the confirmation-to-recovery delay. In principle a positive result could
return after a patient has already recovered and been discharged, but we
assume the reported total reflects confirmed cases carefully recorded as
recovered, so the confirmation-then-recovery ordering holds.

The cumulative recovered series ends at the cut-off, so its between-vintage
increments are fitted as observed `~` data with a NegativeBinomial whose
dispersion is sampled here, NOT shared with the other streams (the recovered
signal has its own observation noise). The convolution right-censors
recoveries that have not yet resolved by the cut-off, so a small observed
recovered count is consistent with a high eventual survival fraction and a
long recovery delay. Empty by default; a `missing` cut-off total leaves the
increments missing (the predictive-generator path). Returns the recovery
probability, the recovery-delay mean, the dispersion, the daily recovered
series and the cut-off total.
"""
@model function recovered_model(
        recovered_history,
        recovered_total::Union{Missing, Integer},
        confirmed_daily::AbstractVector, CFR::Real;
        recovery = recovery_probability_model,
        dispersion = surveillance_dispersion_model(),
        ## Confirmation-to-recovery (discharge) delay; an Ebola survivor is
        ## discharged a couple of weeks after confirmation, so the default is
        ## a mean ~14 d stay before recovery is recorded.
        confirmation_to_recovery = censored_delay_model(
            cdf_nmax(lognormal_meansd(14.0, 8.0); q = 0.99);
            mean_prior = truncated(Normal(14.0, 5.0); lower = 1),
            sd_prior = truncated(Normal(8.0, 4.0); lower = 1)),
        ## Dispersion can be injected from the joint composer's pooled set
        ## (`k_external`); standalone it samples its own from `dispersion`.
        k_external::Union{Nothing, Real} = nothing)
    ## Recovery fraction grounded on the CFR complement (see
    ## `recovery_probability_model`), adjusted for the confirmed population.
    rec_state ~ to_submodel(recovery(CFR))
    p_recover = rec_state.p_recover
    if k_external === nothing
        disp_state ~ to_submodel(dispersion)
        k = disp_state.k
    else
        k = k_external
    end
    delay_state ~ to_submodel(confirmation_to_recovery)

    ## Survivors among confirmed cases, lagged by the confirmation-to-recovery
    ## delay: a scaled convolution of the modelled daily confirmed incidence.
    recovered_daily = p_recover .* convolve_delay(confirmed_daily,
        delay_state.pmf)

    n = length(confirmed_daily)
    vobs = vintage_obs(recovered_history, recovered_total, n)
    modelled_inc = bin_increments(recovered_daily, vobs.days)
    recovered_increments ~ to_submodel(
        vintage_increments_model(modelled_inc, vobs.obs_increments, k))

    expected_recovered := safe_rate(sum(recovered_daily))

    return (; p_recover, recovery_delay_mean = delay_state.mean,
        k_recovered = k, recovered_daily, expected_recovered)
end
