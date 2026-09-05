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
NaN / Inf-safe overdispersed-`Binomial` (`BetaBinomial`) constructor
parameterised by the trial count `n`, the mean positive probability `p`
and an intra-window overdispersion `ρ ∈ (0, 1)`. The `BetaBinomial(n, α, β)`
has mean `n·p` and variance `n·p·(1 − p)·(1 + (n − 1)·ρ)` with `α = s·p`,
`β = s·(1 − p)` and concentration `s = (1 − ρ)/ρ`, so `ρ → 0` (`s → ∞`)
recovers the plain `Binomial(n, p)` and larger `ρ` inflates the variance
above the Binomial. This adds the extra-Binomial variation the confirmed
positives carry (day-to-day laboratory batching and within-window
positivity heterogeneity the pooled per-window `p` does not capture),
mirroring how the count streams use the overdispersed `NegativeBinomial`
rather than a Poisson. `ρ` is floored away from `0` (capping `s`) and `p`
clamped into `(0, 1)` so the distribution stays defined under extreme NUTS
proposals during warmup. Shared by the confirmed-positives windows.
"""
function safe_betabinomial(n::Integer, p, ρ)
    T = float(promote_type(typeof(p), typeof(ρ)))
    lo = eps(T)
    ## Clamp the positivity into the open unit interval.
    pc = isfinite(p) ? clamp(T(p), lo, one(T) - lo) : one(T) / 2
    ## Floor ρ at 1e-6 (cap the concentration `s` at ≈1e6) so a near-zero
    ## draw stays a well-conditioned near-Binomial rather than an infinite
    ## `s`, and cap it below 1 so `s` stays positive.
    ρc = isfinite(ρ) ? clamp(T(ρ), T(1e-6), one(T) - T(1e-6)) : T(1e-6)
    s = (one(T) - ρc) / ρc
    α = max(s * pc, lo)
    β = max(s * (one(T) - pc), lo)
    return BetaBinomial(n, α, β)
end

"""
Modelled between-vintage increments of a daily series `daily`, summed
directly into the bins delimited by the vintage day indices `days` (1-based
into the grid, ascending). The first increment is the cumulative count up
to the first vintage day (`sum(daily[1:days[1]])`). Each later increment is
the inter-vintage sum `sum(daily[days[i-1]+1:days[i]])`, so only the first
bin needs the cumulative. This avoids the cumulative-then-difference round
trip — the daily series is binned once into the quantities the likelihood
scores. Day indices are clamped to the grid. Pure and AD-transparent. The
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
analysed-specimen capacity before testing began: no specimens are
analysed before the testing system exists, so the modelled analysed
volume (and the confirmed counts derived from it) must not accrue over
the pre-surveillance cryptic phase. The suspected-case and suspected-death
streams are not gated: those counts did accumulate over the cryptic
phase, so their first per-vintage bin legitimately rolls from grid day 1.
`start ≤ 1` returns `v` unchanged. Pure and AD-transparent. The element
type follows `v`. Callers concretise any `Vector{Any}` (predict mode)
before gating, so `zero(eltype(v))` is well defined.
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
the grid). Day `t ≤ days[1]` takes `rate[1]`, a day in `(days[i-1],
days[i]]` takes `rate[i]`, and any day beyond the last vintage takes the
last rate (a flat carry-forward of the final window). When `days` is
empty the whole grid takes `rate[1]` if present, else zero, so a scalar
background is recovered. Pure and AD-transparent. The element type
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
`missing` the counts are sampled, the predictive-generator path. The
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
is a NegBinomial around the latent bed demand `means[i]`, right-censored at
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

Per-report-day right-censoring bound for the isolation-occupancy
likelihood, built from data only (no sampled parameter), so the censored
NegativeBinomial never sits against a moving `-Inf` wall. Each isolation
day takes the nearest recorded implied-capacity value (`capacity_history`,
the occupancy / reported occupancy-rate series), floored at that day's
observed occupancy so the bound is never below the count (a count above
the cap has zero probability under the censored NB, the discontinuity
that drives the divergences). With the cap fixed, the censored likelihood
is smooth in the sampled demand: below the cap the count identifies
demand directly, and at the cap it contributes the one-sided
`demand ≥ capacity` tail. A `capacity_history` with no counts (empty, or
days-only as the predictive generator supplies) returns a large finite
cap far above any bed count, so the censoring is a no-op and the
likelihood is the plain NB (a literal `Inf` cannot be used because
[`safe_rate`](@ref) maps non-finite values to `eps`). A `missing`
`iso_obs` skips the occupancy floor.

The latent bed demand itself is left uncensored: the cap enters only as
the observation bound, so demand above capacity is carried by the
renewal / length-of-stay demand model and its priors, not by the bound.
Demand above a saturated capacity is only partially identified from
occupancy (occupancy can say "demand was at least the beds filled", not
how much more), so the bed shortfall is a prior/model-informed quantity,
not pinned by the occupancy data. See [`treatment_flow_model`](@ref).
"""
function censoring_cap(iso_days, iso_obs, capacity_history)
    cdays = Int.(capacity_history.days)
    ccounts = Float64.(capacity_history.counts)
    ## Large finite no-op cap (bed counts are in the hundreds), used when there
    ## is no recorded capacity. Not `Inf`: `safe_rate` maps non-finite to eps.
    nocap = 1.0e6
    ## Censor only where recorded capacity counts exist. The predictive
    ## generator passes the capacity history days-only (counts emptied), so
    ## the guard is on the counts, not the days: otherwise the nearest-day
    ## lookup would index an empty counts vector.
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
    cumulative_occupancy_offset(iso_days, break_days, b)

Per-isolation-day cumulative occupancy reclassification offset `Δ` from the
fitted break steps `b` on the manually specified `break_days` (grid day-
indices). `Δ(t)` sums `b_j` over the break days at or before `iso_days[t]`,
zero before the first break day. Added to the modelled occupancy mean in the
censored-occupancy likelihood so the modelled occupancy tracks a between-report
measurement-basis discontinuity in the observed series (e.g. a Tableau 6-sum →
page-1 headline transition) and the residual stays smooth. Each `b_j` is
sampled in [`treatment_flow_model`](@ref), so `Δ` is a fitted quantity: the fit
partitions each step into reporting artifact vs real demand change. Returns a
length-`length(iso_days)` vector of `eltype(b)` (zeros when there are no break
days, a no-op).
"""
function cumulative_occupancy_offset(iso_days, break_days, b)
    m = length(iso_days)
    T = isempty(b) ? Float64 : eltype(b)
    Δ = zeros(T, m)
    isempty(break_days) && return Δ
    idays = Int.(iso_days)
    bdays = Int.(break_days)
    ## Cumulative sum of the fitted steps up to and including each iso day.
    @inbounds for t in 1:m
        s = zero(T)
        for (j, bd) in enumerate(bdays)
            bd <= idays[t] && (s += b[j])
        end
        Δ[t] = s
    end
    return Δ
end

"""
    break_step_centres(window_days, increments, break_days, gross)

Break days that land on a scored window, paired with the data-derived
centre for each one's fitted step. The centre is
`observed increment − gross`: the part of that vintage's step the report
itself attributes to integrating a harmonised provincial base rather
than to the day's own notifications, where `gross` is the printed 24h
new-confirmed count. On 22 July 2026 that is `369 − 97 = 272` cases and
`236 − 62 = 174` deaths.

Centring on published data rather than zero is what makes the step
identifiable: a zero-centred step with a wide prior has to discover the
harmonisation magnitude from a single observation, whereas the report
states it. The sampled deviation around the centre then carries only the
residual uncertainty about how much of the discrepancy is truly
retrospective.

Break days not matching a window are dropped so no inert step is
sampled. A missing/short `gross` is treated as a gross of zero, so the
centre becomes the whole increment: the entire vintage step is
attributed to the artefact. This is a deliberately conservative fallback
rather than a neutral one: the de-anchor drops the day's positivity
denominator, and a step centred near zero over a de-anchored day is the
pathological configuration (measured on `confirmed_only_model`: 94
divergences and a min bulk ESS of 15, against 20 and 522 with no break
day, because nothing absorbs the backlog and the fit books it as
incidence). Supply the printed count so the split is data-derived.
Omitting it errs towards artefact, never towards incidence.

`increments` may itself be `missing`, the generator path where no
cumulative counts are supplied to difference: the day still gets a step,
so a predictive keeps the fitted chain's dimensions, but there is no
published discrepancy to centre it on and the centre falls back to zero.
Returns `(days, centres)`, both possibly empty.
"""
function break_step_centres(window_days, increments, break_days, gross)
    days = Int[]
    centres = Float64[]
    isempty(break_days) && return (days, centres)
    wdays = Int.(window_days)
    for (j, bd) in enumerate(Int.(break_days))
        pos = findfirst(==(bd), wdays)
        pos === nothing && continue
        push!(days, bd)
        if ismissing(increments)
            push!(centres, 0.0)
        else
            g = j <= length(gross) ? Int(gross[j]) : 0
            push!(centres, Float64(Int(increments[pos]) - g))
        end
    end
    return (days, centres)
end

"""
    confirmed_break_offset(late_days, break_days, b)

Per-late-window confirmed harmonisation offset from the fitted break
steps `b` on the manually specified `break_days` (grid day-indices).
Unlike [`cumulative_occupancy_offset`](@ref) this is not cumulative: the
confirmed likelihood fits between-vintage increments, so a one-off
retrospective base integration inflates exactly one increment. Every
later vintage differences against the new, higher base and is
unaffected. Each entry is therefore the step for that window alone, zero
on every non-break window.

Added to the modelled confirmed mean of the break window in
[`confirmed_cases_model`](@ref) and [`confirmed_deaths_model`](@ref), so the
fit can explain an increment that is mostly reattached backlog without
inflating the latent incidence that drives `Rt`. Each `b_j` is sampled there
around the published discrepancy ([`break_step_centres`](@ref)), symmetric
because a harmonisation can revise down as well as up. Returns a
length-`length(late_days)` vector of `eltype(b)` (zeros when there are no break
days, a no-op).
"""
function confirmed_break_offset(late_days, break_days, b)
    m = length(late_days)
    T = isempty(b) ? Float64 : eltype(b)
    Δ = zeros(T, m)
    isempty(break_days) && return Δ
    ldays = Int.(late_days)
    bdays = Int.(break_days)
    @inbounds for (j, bd) in enumerate(bdays)
        for t in 1:m
            ldays[t] == bd && (Δ[t] += b[j])
        end
    end
    return Δ
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
(onset-to-admission and admission-to-death, the `bdbv-linelist-analysis`
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
genuine per-day count that never falls, so it fits directly. Its days fall
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
    ## (clamped into the grid), not a between-vintage increment: this is a
    ## genuine daily count, so it never differences a falling cumulative. Empty
    ## by default. A `missing` count vector samples (the predictive path). The
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
incidence that never falls, so it fits directly. It shares the suspect
pipeline and dispersion with the cumulative stream and is empty by default.
Its days fall strictly after the cumulative series ends, so the two suspected
likelihoods cover disjoint days.

The background and testing fraction are sampled by an injected
[`test_positivity_model`](@ref), and the onset-to-report delay is
injected, defaulting to a weakly-informative prior on the
onset-to-notification delay (mean 4.5 d, SD 3.6 d), consistent with
Ebola surveillance reporting delays.

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
        ## delay (~20 d) is not used: it is assumed to reflect a longer
        ## pathway (likely confirmation and administrative processing),
        ## though what it captures is uncertain.
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
    ## grid), not a between-vintage increment: this is a genuine daily
    ## incidence, so it never differences a falling cumulative. Empty by
    ## default. A `missing` count vector samples (the predictive path).
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

    return (; p_drc, λ_bg = λ_bg_base, τ_test, report_pmf,
        report_mean = report_state.mean, report_sd = report_state.sd,
        bvd_reports_daily,
        reports_daily, expected_reports, positivity, bg_daily, bg_sigma,
        bg_total)
end

"""
Per-window confirmed-positives likelihood, expressed as a vector
likelihood scored against the observed analysed denominators. Given the
per-window analysed counts `analysed`, positivities `p_pos` and the
intra-window overdispersion `ρ`, scores the observed `positives` with one
`BetaBinomial(analysed[i], p_pos[i], ρ)` per window (see
[`safe_betabinomial`](@ref)). The overdispersion adds the extra-Binomial
variation the confirmed counts carry. `ρ → 0` recovers the plain
`Binomial(analysed[i], p_pos[i])`.
`positives` is the model argument on the LHS of `~`, so a supplied vector
is observed data DynamicPPL conditions on (mirroring
[`vintage_increments_model`](@ref)) and a `missing` argument is sampled,
making the submodel a predictive generator. Under a prefixed submodel
attachment the predict keys are `<prefix>.positives[i]`. Returns the
observed (or sampled) positives.
"""
@model function confirmed_positives_model(
        positives::Union{Missing, AbstractVector{<:Integer}},
        analysed::AbstractVector{<:Integer}, p_pos::AbstractVector,
        ρ::Real = 0.0)
    nv = length(analysed)
    if ismissing(positives)
        positives = Vector{Union{Missing, Int}}(missing, nv)
    end
    for i in 1:nv
        positives[i] ~ safe_betabinomial(analysed[i], p_pos[i], ρ)
    end
    return (; positives)
end

"""
Late confirmed-vintage likelihood, one per-window confirmed increment over
the days after the last cumulative laboratory date. Each window scores its
increment two ways depending on whether a 24h analysed denominator was
published for that day (`analysed[i] > 0`): an anchored day is an
overdispersed `BetaBinomial(analysed[i], p_pos[i], ρ)` of the observed
denominator (like an observed window, see [`safe_betabinomial`](@ref)), a
unanchored day a `NegativeBinomial` of the modelled laboratory volume
`modelled[i]` sharing the surveillance dispersion `k`.
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
        p_pos::AbstractVector, k::Real, ρ::Real = 0.0)
    n = length(modelled)
    if ismissing(increments)
        increments = Vector{Union{Missing, Int}}(missing, n)
    end
    for i in 1:n
        if analysed[i] > 0
            increments[i] ~ safe_betabinomial(analysed[i], p_pos[i], ρ)
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
  onward. The first cumulative value is the baseline). Each window's
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
  A day listed in `confirmed_break_days` is the exception: its denominator is
  dropped (`late_analysed` set to 0) so it stays a modelled-volume window.
  On such a day INSP has retrospectively integrated a harmonised provincial
  base, so most of the confirmed increment is reattached backlog rather than
  positives out of that day's specimens, and pairing the two would score a
  wildly inflated positivity (22 July 2026: 369 of 414 analysed, 89% implied
  against the 23% actually reported). The increment is still scored. The
  backlog is absorbed by a fitted level step in
  [`confirmed_cases_model`](@ref) (see [`confirmed_break_offset`](@ref)).

The three groups partition the confirmed counts at the first and last
laboratory dates, so no confirmed case is counted twice. Returns
`(; obs_days, obs_positives, obs_analysed, early_days, early_increments,
late_days, late_increments, late_analysed, late_start)` of grid
day-indices and per-window counts. `late_analysed[i]` is the observed 24h
analysed denominator for late day `i` (0 when none was published). The
observed and late groups are empty when no laboratory history is present
and every confirmed vintage becomes an early window. Pure integer
bookkeeping on the observed data, so it carries no gradient.
"""
function confirmed_positivity_windows(confirmed_history, lab_history,
        lab_daily_history = (; days = Int[], counts = Int[]),
        confirmed_break_days = Int[])
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

    ## The first confirmed vintage is the baseline (the initial cumulative
    ## level the surveillance system was at when reporting began). It is not
    ## scored. The vintaging "starts with the data" from there, so no window's
    ## modelled volume rolls over the pre-surveillance cryptic phase. This
    ## matches how the observed and late laboratory windows already treat their
    ## first value as a baseline. `early_start` is that first confirmed day.
    ## The early windows pin their modelled-volume accumulation at it.
    early_start = Int(cdays[1])

    ## No laboratory denominators: every confirmed vintage after the baseline
    ## is an early window scored through the modelled laboratory volume,
    ## binned from `early_start`.
    if isempty(lab_history.counts)
        ed = Int[Int(d) for d in cdays[2:end]]
        inc = Int[Int(ccounts[i]) - Int(ccounts[i - 1])
                  for i in 2:length(ccounts)]
        return (; obs_days = Int[], obs_positives = Int[],
            obs_analysed = Int[], early_days = ed, early_increments = inc,
            early_start, late_days = Int[], late_increments = Int[],
            late_analysed = Int[], late_start = 0)
    end

    first_lab_day = Int(lab_history.days[1])
    ## Early windows: confirmed vintages after the first confirmed (baseline)
    ## up to and including the first laboratory date, scored against the
    ## modelled laboratory volume binned over each window's own range from
    ## `early_start`.
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

    ## Retrospective harmonisation days are de-anchored: their published 24h
    ## analysed count is dropped, so the window falls to the modelled-volume
    ## NegativeBinomial branch instead of entering the BetaBinomial as
    ## same-day positives. On such a day most of the confirmed increment is a
    ## provincial base integration, not specimens that tested positive out of
    ## that day's denominator, so pairing the two would score a wildly
    ## inflated positivity (22 July 2026: 369 of 414 = 89%, against the 23%
    ## actually reported) and drag `λ_bg`, the composition link and the
    ## ascertainment with it. The increment itself is still scored, with a
    ## fitted level step absorbing the backlog (see `confirmed_cases_model`).
    ## Plain integer bookkeeping in a loop, like the block above, so the AD
    ## backends can compile a rule for it on every Julia version.
    for d in confirmed_break_days
        di = Int(d)
        pos = findfirst(==(di), late_days)
        pos === nothing || (late_analysed[pos] = 0)
    end

    return (; obs_days, obs_positives, obs_analysed, early_days,
        early_increments, early_start, late_days, late_increments,
        late_analysed, late_start = last_lab_day)
end

## Per-window tested-positive probability `p_pos` under the `:composition`
## link, extracted from `confirmed_cases_model` as a plain function so its
## working variables are function arguments rather than closure captures.
## Inside the `@model` this lived in a `map(...) do i` closure nested in an
## `if positivity_link === :composition` branch. The conditionally-scoped
## captures were boxed in a `Base.RefValue`, which Enzyme's reverse mode
## cannot differentiate through (Mooncake tolerates the box). As a
## top-level function the captures become parameters, and the explicit
## indexed loop (not `map(do i)`) creates no anonymous-closure shadow for
## Enzyme to construct. The result is bit-identical to the in-model `map`
## under Mooncake.
function composition_positivity(window_days, bvd_window, bg_window,
        c_window, δ0, dscale, s_test, spec, lo, hi)
    Tt = eltype(bvd_window)
    s_t = convert(Tt, s_test)
    sp_t = convert(Tt, spec)
    nw = length(window_days)
    p_pos = Vector{Tt}(undef, nw)
    @inbounds for i in 1:nw
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
        ## with the composition so the confirmed BetaBinomial always sees
        ## a valid probability even under an AD perturbation.
        p_pos[i] = clamp(isfinite(p) ? p : φ, lo, hi)
    end
    return p_pos
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
- Confirmed positives. The confirmed counts are scored as an overdispersed
  `BetaBinomial` of the observed specimens-*analysed* denominator in each
  laboratory window ([`confirmed_positivity_windows`](@ref),
  [`safe_betabinomial`](@ref)), with a partially-pooled per-window
  positivity ([`confirmed_positivity_model`](@ref)) and a shared
  intra-window overdispersion ([`confirmed_overdispersion_model`](@ref)) so
  the confirmed intervals capture the extra-Binomial day-to-day laboratory
  variation rather than the far-too-tight plain Binomial spread.
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
from `λ_bg`. It leaves `λ_bg` weakly identified and is kept for sensitivity
analysis.

The observed-window positives are conditioned on the observed analysed
denominator, so the Binomial conditioning that removes the multiplicative
ascertainment ridge is preserved. Only the early and unanchored late
windows use the modelled analysed volume as the denominator.

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
        overdispersion = confirmed_overdispersion_model(),
        ## When false, the early/late windows (confirmed vintages with no
        ## observed analysed denominator) are not scored: only the
        ## observed-denominator Binomial windows contribute, so confirmed
        ## informs positivity without extrapolating a denominator from
        ## incidence. Used to probe the no-test-data extrapolation.
        fit_unanchored::Bool = true,
        ## Opt-in retrospective harmonisation-break days (grid day-indices):
        ## days whose cumulative confirmed step is mostly a provincial base
        ## integration rather than 24h notifications. De-anchored from the
        ## positivity denominator and given a fitted level step in the
        ## modelled mean. Empty (the default) is a no-op.
        confirmed_break_days::AbstractVector{<:Integer} = Int[],
        ## Printed 24h new-confirmed counts on each break day. The step is
        ## centred on `observed increment − gross`, the part of the vintage
        ## step the report attributes to base integration, so the magnitude
        ## comes from published data rather than a prior guess. Empty or
        ## all-zero centres on the whole increment, attributing all of it to
        ## the artefact: conservative, not neutral.
        confirmed_break_gross::AbstractVector{<:Integer} = Int[],
        ## Residual uncertainty about how much of that discrepancy is truly
        ## retrospective rather than coincident same-day incidence, not the
        ## harmonisation magnitude, which `confirmed_break_gross` supplies.
        ## Zero pins the step at the published discrepancy and samples no
        ## parameter, taking the printed count as exact.
        confirmed_break_sd::Real = 25.0)
    n = length(onsets)
    ## `missing` cut-off scalar means generator mode: observed increments are
    ## left missing so `predict` resamples them.
    have_data = !ismissing(confirmed_cases)

    ## Intra-window overdispersion for the confirmed positives. The observed
    ## and anchored-late windows are scored as an overdispersed BetaBinomial
    ## of the observed analysed denominator (`safe_betabinomial`), so the
    ## confirmed intervals are not the far-too-tight plain Binomial. Sampled
    ## once and shared across all confirmed windows regardless of the
    ## positivity link.
    od_state ~ to_submodel(overdispersion, false)
    ρ_conf = od_state.ρ

    ## Laboratory capacity onset. No specimens are analysed before testing
    ## existed, so the modelled analysed volume is gated to zero before the
    ## first confirmed-case vintage (the earliest evidence of testing. The
    ## first laboratory date is the fallback). Modelling a pre-testing
    ## analysed volume would invent capacity that did not exist and roll it
    ## into the first laboratory and early-confirmed bins, vastly
    ## over-predicting the early confirmed counts. The suspected-case
    ## pipeline feeding the volume is not gated: suspected cases did
    ## accumulate over the cryptic phase.
    cap_start = !isempty(confirmed_history.days) ?
                clamp(Int(confirmed_history.days[1]), 1, n) :
                (!isempty(lab_history.days) ?
                 clamp(Int(lab_history.days[1]), 1, n) : 1)

    ## Analysed-specimen volume: the suspected pipeline carried through the
    ## report-to-analysed delay and thinned by the tested fraction, fit to
    ## the analysed series and reused as the denominator in the early and
    ## unanchored late windows below. `bg_daily` is the per-day non-BVD
    ## background.
    receipt_state ~ to_submodel(receipt)
    suspected_daily = p_drc .* bvd_reports_daily .+ bg_daily
    analysed_daily = τ_test .* convolve_delay(suspected_daily,
        receipt_state.pmf)
    ## In predict mode (no AD) the daily series can infer as `Vector{Any}`
    ## on some Julia versions, which then makes `reduce_empty` / `zero(Any)`
    ## fail on the empty derived window vectors below. Concretise to the
    ## working scalar type. This runs only when the element type has widened,
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
    ## series stops, INSP publishes a 24h analysed count on some days. The
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
        lab_daily_history, confirmed_break_days)
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
        ## The per-window pool composition φ, severity upshift, assay
        ## sensitivity/specificity transform and (0,1) guards live in
        ## `composition_positivity` (a plain function with an explicit loop,
        ## so its working variables are not boxed closure captures under
        ## Enzyme).
        dscale = max(convert(Tt, decay_scale), one(Tt))
        p_pos = composition_positivity(window_days, bvd_window, bg_window,
            c_window, δ0, dscale, s_test, spec, lo, hi)
    else
        pos_state ~ to_submodel(positivity(nv))
        p_pos = pos_state.p_pos
    end

    ## Early windows: confirmed increment ~ NegBinomial(positivity ×
    ## modelled analysed volume), the volume binned over each window's own
    ## day range pinned at `early_start` (the first confirmed vintage, the
    ## testing-onset baseline), so the first early increment is scored from
    ## the data start rather than rolling the (now-gated) pre-testing volume.
    ## Mirrors the late-window pinning at `late_start`.
    early_p = p_pos[1:n_early]
    ## `early_volume` is bound unconditionally to one array allocation site
    ## (`bin_increments(...)[2:end]`) rather than a two-branch ternary: the
    ## ternary compiles to a pointer that is one of two array allocations
    ## depending on a branch, which Enzyme's `nodecayed_phis!` LLVM pass
    ## cannot trace (an `EnzymeInternalError`). With no late/early window
    ## days the edge is the singleton `[start]`, so
    ## `bin_increments(...)[2:end]` is empty, the same value a
    ## `similar(analysed_daily, 0)` branch would give. Mooncake is
    ## unaffected.
    early_edges = n_early > 0 ?
                  vcat(windows.early_start, windows.early_days) :
                  [windows.early_start]
    early_volume = bin_increments(analysed_daily, early_edges)[2:end]
    early_mean = early_p .* early_volume
    early_obs = (have_data && n_early > 0 && fit_unanchored) ?
                windows.early_increments : missing
    early_increments ~ to_submodel(
        vintage_increments_model(early_mean, early_obs, k))

    ## Observed windows: overdispersed BetaBinomial of the observed analysed
    ## denominator (`ρ_conf` the intra-window overdispersion).
    obs_p = p_pos[(n_early + 1):(n_early + n_obs)]
    obs_positives = (have_data && n_obs > 0) ? collect(windows.obs_positives) :
                    missing
    confirmed_positives ~ to_submodel(
        confirmed_positives_model(obs_positives, windows.obs_analysed, obs_p,
        ρ_conf))

    ## Late windows: confirmed-only vintages after the last laboratory date.
    ## A day that publishes a 24h analysed count (`late_analysed > 0`) is
    ## scored as an overdispersed BetaBinomial of that observed denominator,
    ## like an observed window, anchoring its positivity to data. Each
    ## remaining unanchored day is NegBinomial(positivity × modelled
    ## volume). The modelled volume is binned over each late window's own
    ## day range, with the running edge pinned at the last laboratory day
    ## (`late_start`): `bin_increments` runs its running `prev` from day 0,
    ## so prepending `late_start` to the late day edges and dropping the
    ## synthetic first bin starts the accumulation at `late_start`, avoiding
    ## double-counting the observed-window volume.
    late_p = p_pos[(n_early + n_obs + 1):nv]
    ## Unconditional single-allocation binding, avoiding the same
    ## two-allocation pointer-PHI Enzyme's `nodecayed_phis!` pass cannot
    ## trace (see `early_volume` above). Empty when `n_late = 0`.
    late_edges = n_late > 0 ? vcat(windows.late_start, windows.late_days) :
                 [windows.late_start]
    late_volume = bin_increments(analysed_daily, late_edges)[2:end]
    ## Opt-in retrospective harmonisation step. On a declared break day the
    ## observed increment is mostly a provincial base integration, so a single
    ## level step is fitted into that window's modelled mean and the fit
    ## partitions the increment into reporting artefact vs real incidence with
    ## uncertainty. Sampled non-centred around the published discrepancy
    ## (`break_step_centres`: the vintage increment minus the printed 24h
    ## count), with `confirmed_break_sd` the residual uncertainty about how much
    ## of that discrepancy is truly retrospective. Symmetric, because that
    ## residual can fall either side. Only break days that land
    ## on a late window can move the likelihood. Others are dropped so no inert
    ## step is sampled. Same conditional-`b`-then-pure-function shape as the
    ## occupancy break (`cumulative_occupancy_offset`), which keeps the
    ## allocation single-site for the AD backends. Empty → Δ = 0, a no-op.
    late_day_ints = Int.(windows.late_days)
    conf_brk_days, conf_brk_centre = break_step_centres(late_day_ints,
        windows.late_increments, confirmed_break_days, confirmed_break_gross)
    if isempty(conf_brk_days)
        cb = Float64[]
    elseif iszero(confirmed_break_sd)
        ## `confirmed_break_sd = 0` is the deterministic correction: the report
        ## publishes the 24h count, so the artefact size is taken as exact and
        ## the increment is fitted against the gross count with no sampled
        ## parameter. Two fewer parameters, and no ridge between a step and the
        ## ascertainment it trades off against.
        cb = conf_brk_centre
    else
        confirmed_step ~ product_distribution(
            fill(Normal(0, 1), length(conf_brk_days)))
        cb = conf_brk_centre .+ confirmed_break_sd .* confirmed_step
    end
    late_break_offset = confirmed_break_offset(windows.late_days,
        conf_brk_days, cb)
    late_mean = late_p .* late_volume .+ late_break_offset
    ## Observed late increments: anchored days (24h denominator) carry the
    ## confirmed increment clamped into the Binomial support and are always
    ## scored. Unanchored days are scored only when `fit_unanchored` (the
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
        late_p, k, ρ_conf))

    ## Plain `=`, not `:=`: these derived quantities are surfaced through the
    ## returned NamedTuple and re-tracked at the joint level (`joint.jl`
    ## `expected_confirmed_T`/`expected_analysed_T`/`test_positivity`), so the
    ## submodel-level `:=` tracking is redundant. `:=` also builds a
    ## DynamicPPL tracking closure that captures `p_pos` (assigned in both the
    ## `:composition` and `else` positivity branches, so boxed in a
    ## `Core.Box`), and Enzyme's `nodecayed_phis!` LLVM pass cannot
    ## differentiate through the resulting boxed pointer-PHI.
    expected_analysed = safe_rate(sum(analysed_daily))
    ## Expected confirmed at the cut-off and the overall positivity, over the
    ## modelled early volume, the observed cumulative analysed windows and the
    ## late windows (anchored days contribute `p · analysed`, unanchored days
    ## the modelled `p · volume`). The window vectors are empty when a
    ## vintage has no such window, and in predict mode their element type
    ## can widen to `Any`, so each sum is given a concrete `init` to skip
    ## `reduce_empty`'s `zero(Any)`. The init is taken from the scalar
    ## `τ_test` (always concrete), not from `eltype(p_pos)`, which can widen
    ## to `Any`.
    z = zero(τ_test)
    amask = windows.late_analysed .> 0
    late_den_a = float.(windows.late_analysed)
    ## Unconditional broadcasts (single allocation site each): with
    ## `n_late = 0` every operand (`amask`, `late_den_a`, `late_volume`,
    ## `late_mean`) is empty, so the `ifelse.` result is empty too, the same
    ## value a `? … : similar(…, 0)` ternary would produce, without the
    ## two-allocation pointer-PHI that Enzyme's `nodecayed_phis!` pass
    ## cannot trace. Mooncake is unaffected.
    late_den = ifelse.(amask, late_den_a, late_volume)
    late_expected = ifelse.(amask, late_p .* late_den_a, late_mean)
    denom = sum(early_volume; init = z) + float(sum(windows.obs_analysed)) +
            sum(late_den; init = z)
    expected_positives = sum(early_mean; init = z) +
                         sum(late_expected; init = z) +
                         (n_obs > 0 ? sum(obs_p .* windows.obs_analysed) : z)
    expected_confirmed = safe_rate(expected_positives)
    p_positive = safe_rate(expected_positives) / safe_rate(denom)

    ## Modelled daily confirmed-case incidence: per-window tested-positive
    ## probability expanded onto the daily grid times the modelled analysed
    ## volume. In predict mode `p_pos` can widen to `Vector{Any}`, so pin it to
    ## the analysed volume's element type before expanding.
    p_pos_daily = p_pos
    if eltype(p_pos_daily) === Any
        p_pos_daily = convert(Vector{eltype(analysed_daily)}, p_pos_daily)
    end
    ## Per-day positivity `p_pos_grid`, exposed with `τ_test` so the treatment
    ## model can form the in-care confirmation hazard `τ_test · p_pos_grid[t]`.
    p_pos_grid = expand_vintage_rate(p_pos_daily, window_days, n)
    confirmed_daily = p_pos_grid .* analysed_daily

    return (; τ_test, bg_daily, p_pos, p_pos_grid, windows, analysed_daily,
        confirmed_daily,
        receipt_pmf = receipt_state.pmf,
        receipt_mean = receipt_state.mean, receipt_sd = receipt_state.sd,
        expected_analysed, expected_confirmed, p_positive)
end

"""
Uganda exports likelihood (geographic spread). The exports stream is
travel-gated, so the at-risk clock starts at infection: a traveller
moves and is exported during incubation (pre-symptomatic) and stays at
risk of being exported and detected abroad only until the
infection→detection delay has elapsed. The expected detected exports by
the cut-off therefore accumulate the per-capita travel rate `q =
daily_travellers / source_population` over the at-risk person-time, not
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
it, keyed to infection like `C(s)`). `detected` is the running sum of
`convolve_delay(infections, f_det)`. Summing the daily at-risk prevalence
is the discrete person-time integral. Summing `q · onsets` instead would
charge each case only a single day of travel risk, under-counting exports
by roughly the mean at-risk dwell time. Samples the traveller volume and
the onset-to-detection delay via injected submodels. The onset-to-detection
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
        ## `t_last` (the `last_offset` truncation). Prevalence past it does
        ## not accrue. `d₁` is the earliest detection day.
        days, counts = dated_event_bins(export_case_days, n)
        d₁ = days[1]
        ## Pre-detection survival weight Λ(d₁−1): the cumulative export
        ## intensity up to the day before the earliest detection.
        pre = d₁ > 1 ? sum(@view export_prevalence[1:(d₁ - 1)]) :
              zero(@inbounds export_prevalence[begin])
        pre_detection_exports ~ Poisson(safe_rate(pre))
        ## Per-day-edge increments between consecutive detection days. The
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

    ## Travel-scaled at-risk prevalence without the export-case ascertainment
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
the travelled at-risk person-time `q · prevalence` from
[`exports_model`](@ref), before the export-case ascertainment `p_uganda`: a
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
        ## Dated per-day Poisson. The clock stops at the last death day.
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
    testing-intensity difference between deaths and living suspects. The
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
        ## Opt-in retrospective harmonisation-break days, as in
        ## `confirmed_cases_model`: the 22 July 2026 base integration steps the
        ## confirmed-death cumulative by +236 against a printed 24h count of
        ## +62, and this stream fits between-vintage increments too, so the
        ## backlog would otherwise be read as one day of confirmed deaths.
        ## Empty (the default) is a no-op.
        confirmed_break_days::AbstractVector{<:Integer} = Int[],
        confirmed_break_gross::AbstractVector{<:Integer} = Int[],
        ## As on the cases path: the residual uncertainty around the published
        ## discrepancy, with zero pinning the step and sampling no parameter.
        confirmed_break_sd::Real = 25.0,
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
    ## which trips `zero(Any)` downstream. Pin to the sampled scalar type,
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
            v = den > lo ? sc * case_analysed_daily[t] * susp_death[t] / den :
                zero(sc)
            ## Cap the volume at the suspected-death pool so confirmed deaths
            ## stay a subset of suspected and the realised τ_death ≤ 1.
            min(v, susp_death[t])
        end
        τ_death = susp_death[n] > lo ?
                  death_volume[n] / susp_death[n] : zero(sc)
    else
        test_state ~ to_submodel(testing)
        τ_death = test_state.τ_death
        sc = one(τ_death)
        death_volume = τ_death .* gate_before(susp_death, capacity_start)
    end

    confirmed_death_daily = p_pos_daily .* death_volume
    vobs = vintage_obs(confirmed_deaths_history, confirmed_deaths, n)
    modelled_inc = bin_increments(confirmed_death_daily, vobs.days)
    ## Retrospective harmonisation step, mirroring `confirmed_cases_model`:
    ## centred on the published discrepancy (vintage increment minus the
    ## printed 24h death count) and added to that window's modelled mean, so
    ## the increment likelihood does not attribute reattached deaths to the
    ## day. Sampled whenever a break day lands on a vintage, as on the cases
    ## path, so a posterior predictive carries the same dimensions as the
    ## fitted chain and replicates the break instead of leaving the vintage as
    ## an outlier.
    cd_brk_days, cd_brk_centre = break_step_centres(vobs.days,
        vobs.obs_increments, confirmed_break_days, confirmed_break_gross)
    if isempty(cd_brk_days)
        cdb = Float64[]
    elseif iszero(confirmed_break_sd)
        ## Deterministic correction, as on the cases path: the printed 24h death
        ## count is taken as exact and no step parameter is sampled.
        cdb = cd_brk_centre
    else
        cdeath_step ~ product_distribution(
            fill(Normal(0, 1), length(cd_brk_days)))
        cdb = cd_brk_centre .+ confirmed_break_sd .* cdeath_step
    end
    modelled_inc = modelled_inc .+
                   confirmed_break_offset(vobs.days, cd_brk_days, cdb)
    ## The cut-off scalar is the generator gate, as in `confirmed_cases_model`:
    ## nulling it leaves `predict` to resample the increments while the dated
    ## history still supplies the vintage grid and the break-step centres. That
    ## keeps the published discrepancy available to a predictive, which
    ## differencing an emptied history cannot do.
    cdeath_obs = ismissing(confirmed_deaths) ? missing : vobs.obs_increments
    cdeath_increments ~ to_submodel(
        vintage_increments_model(modelled_inc, cdeath_obs, k))

    expected_confirmed_deaths := safe_rate(sum(confirmed_death_daily))
    ## Cut-off death-pool composition and confirmation positivity, surfaced as
    ## `death_composition` and `death_confirmation`.
    q_death := q_death_daily[n]
    p_death_conf := p_pos_daily[n]

    return (; τ_death, scaling = sc, s_test = s, spec, q_death, p_death_conf,
        confirmed_death_daily, expected_confirmed_deaths)
end

"""
    accumulate_occupancy(A_bvd, A_bg, deaths, recover, ruleout, κ, conf_hazard)

Treatment-centre occupancy built as a forward day-by-day running balance of
admission, discharge and abscond events, carved into a confirmed and a suspect
sub-stock by a confirmation overlay.

The clinical layer is label-independent and carried as two sub-stocks that sum
to the total demand by construction. The occupied BVD true-case stock `O_bvd`
accumulates `A_bvd` minus the BVD discharges (`deaths + recover`) and its share
of the abscond outflow. The occupied non-case stock `O_bg` accumulates `A_bg`
minus the rule-out exits (`ruleout`) and its share of the abscond outflow:

```math
O_\\text{bvd}(t) = O_\\text{bvd}(t-1) + A_\\text{bvd}(t)
    - \\text{deaths}(t) - \\text{recover}(t) - a_\\text{bvd}(t),
\\qquad
O_\\text{bg}(t) = O_\\text{bg}(t-1) + A_\\text{bg}(t)
    - \\text{ruleout}(t) - a_\\text{bg}(t),
```

and the total demand is their sum `D(t) = O_bvd(t) + O_bg(t)` exactly. The
abscond outflow `κ · O_susp(t-1)` drains the suspect (unconfirmed) pool only —
confirmed patients do not abscond — and is split proportionally between its
unconfirmed-BVD portion `O_bvd(t-1) − O_conf(t-1)` and the non-case portion
`O_bg(t-1)`, the two parts of the suspect pool `O_susp(t-1)`:

```math
a_\\text{bvd}(t) = κ\\,O_\\text{susp}(t-1)\\,
    \\frac{O_\\text{bvd}(t-1) - O_\\text{conf}(t-1)}{O_\\text{susp}(t-1)},
\\qquad
a_\\text{bg}(t) = κ\\,O_\\text{susp}(t-1)\\,
    \\frac{O_\\text{bg}(t-1)}{O_\\text{susp}(t-1)}.
```

The confirmation overlay relabels occupied true cases at the daily hazard
`conf_hazard(t) = τ_test · p_pos(t)` borrowed from the lab pipeline, so the
confirmed sub-stock is the confirmed subset of `O_bvd`, never an extra
compartment, draining at the confirmed share of the BVD discharges:

```math
O_\\text{conf}(t) = O_\\text{conf}(t-1)
    + \\text{conf\\_hazard}(t)\\,(O_\\text{bvd}(t-1) - O_\\text{conf}(t-1))
    - (\\text{deaths}(t) + \\text{recover}(t))\\,
        \\frac{O_\\text{conf}(t-1)}{O_\\text{bvd}(t-1)},
```

so `O_conf ≤ O_bvd ≤ D`. The suspect sub-stock is the remainder
`O_susp(t) = D(t) − O_conf(t) = (O_bvd(t) − O_conf(t)) + O_bg(t) ≥ 0`, the
not-yet-confirmed BVD occupancy plus the non-case rule-out occupancy, so
`O_conf + O_susp = D` closes without a clamp. Each sub-stock is floored at zero,
`O_conf` is clamped into `[0, O_bvd]`, and the abscond denominator `O_susp(t-1)`
is floored with `eps` so a zero suspect pool gives a zero split. Returns
`(; demand, O_bvd, O_conf, O_susp, abscond)`, each a length-`n` vector.
"""
function accumulate_occupancy(A_bvd::AbstractVector, A_bg::AbstractVector,
        deaths::AbstractVector, recover::AbstractVector,
        ruleout::AbstractVector, κ::Real, conf_hazard::AbstractVector)
    n = length(A_bvd)
    T = promote_type(eltype(A_bvd), eltype(A_bg), eltype(deaths),
        eltype(recover), eltype(ruleout), typeof(κ), eltype(conf_hazard))
    demand = Vector{T}(undef, n)
    O_bvd = Vector{T}(undef, n)
    O_conf = Vector{T}(undef, n)
    O_susp = Vector{T}(undef, n)
    abscond = Vector{T}(undef, n)
    z = zero(T)
    Obvd_prev = z
    Obg_prev = z
    Oconf_prev = z
    Osusp_prev = z
    ## Floor for the abscond denominator so a zero suspect pool gives a zero
    ## split rather than 0/0.
    ε = eps(T)
    @inbounds for t in 1:n
        bvd_out = deaths[t] + recover[t]
        ab = κ * Osusp_prev
        unconf = max(Obvd_prev - Oconf_prev, z)
        denom = max(Osusp_prev, ε)
        ab_bvd = ab * (unconf / denom)
        ab_bg = ab * (Obg_prev / denom)
        Obvd_t = max(Obvd_prev + A_bvd[t] - bvd_out - ab_bvd, z)
        Obg_t = max(Obg_prev + A_bg[t] - ruleout[t] - ab_bg, z)
        Dt = Obvd_t + Obg_t
        conf_in = conf_hazard[t] * unconf
        share = Obvd_prev > z ? Oconf_prev / Obvd_prev : z
        Oconf_t = clamp(Oconf_prev + conf_in - bvd_out * share, z, Obvd_t)
        Osusp_t = max(Dt - Oconf_t, z)
        demand[t] = Dt
        O_bvd[t] = Obvd_t
        O_conf[t] = Oconf_t
        O_susp[t] = Osusp_t
        abscond[t] = ab
        Obvd_prev = Obvd_t
        Obg_prev = Obg_t
        Oconf_prev = Oconf_t
        Osusp_prev = Osusp_t
    end
    return (; demand, O_bvd, O_conf, O_susp, abscond)
end

"""
    clinical_stay_survival(death_pmf, recover_pmf, cfr)

Clinical-stay survival of the BVD discharge mixture, indexed from cohort age 0,
matching the discharge balance [`accumulate_occupancy`](@ref) runs. An admitted
true case leaves by death (weight `cfr`, the in-care fatality `CFR_iso`) on the
admission→death stay or by recovery (weight `1 − cfr`) on the admission→recovery
stay, so the cumulative discharge fraction by age `d` is `G(d) = Σ_{j=0}^{d}
(cfr · death_pmf[j+1] + (1 − cfr) · recover_pmf[j+1])` and the survival is its
complement

```math
S_\\text{clin}(d) = 1 - G(d) = \\Pr(\\text{still in a bed after day } d),
```

so `S_clin[d+1]` is the probability the case has not yet been discharged by the
end of day `d`. This is the per-cohort weight the running-balance BVD stock
carries: with no absconds `O_bvd(t) = Σ_{u ≤ t} A_bvd(u) · S_clin(t − u)`. The
discharge-complement `P(stay > d)`, rather than the inclusive `P(stay ≥ d)` of
[`convolve_survival`](@ref), keeps the two-clock confirmed sub-stock `≤ O_bvd`
by construction. The two PMFs need not share a length. Each contributes zero
beyond its support, so the result takes the longer length.
"""
function clinical_stay_survival(death_pmf::AbstractVector,
        recover_pmf::AbstractVector, cfr::Real)
    L = max(length(death_pmf), length(recover_pmf))
    T = promote_type(eltype(death_pmf), eltype(recover_pmf), typeof(cfr))
    S = Vector{T}(undef, L)
    one_T = one(T)
    cum = zero(T)
    @inbounds for d in 1:L
        pd = d <= length(death_pmf) ? death_pmf[d] : zero(eltype(death_pmf))
        pr = d <= length(recover_pmf) ? recover_pmf[d] :
             zero(eltype(recover_pmf))
        cum += cfr * pd + (one_T - cfr) * pr
        S[d] = one_T - cum
    end
    return S
end

"""
    two_clock_confirmed(A_bvd, conf_hazard, S_clin)

Cohort-tracked confirmed-in-care sub-stock. A true-case admission cohort
admitted on day `u` (`A_bvd[u]`) contributes to the confirmed stock on day `t`
only once it has been confirmed and while it has not yet clinically departed:

```math
O_\\text{conf}(t) = \\sum_{u \\le t} A_\\text{bvd}(u)\\,
    \\text{CDF}_\\text{conf}(u, t)\\, S_\\text{clin}(t - u),
```

a product of two clocks running from admission. The confirmation clock is
the per-cohort cumulative confirmation probability under the time-varying
daily hazard `conf_hazard(t) = τ_test · p_pos(t)` borrowed from the lab
pipeline, with exposure from the day after admission:

```math
\\text{CDF}_\\text{conf}(u, t) = 1 - \\prod_{j = u+1}^{t}\\bigl(1 -
    \\text{conf\\_hazard}(j)\\bigr),
```

so a just-admitted cohort (`u = t`) has `CDF_conf = 0`. The clinical clock
is the BVD stay survival `S_clin(d) = P(stay > d)`
([`clinical_stay_survival`](@ref)), so
`Σ_{u ≤ t} A_bvd(u) · S_clin(t − u)` reconstructs the abscond-free BVD
stock. The confirmed-and-present cohort is a subset of the present
cohort, so `O_conf ≤ O_bvd` holds by construction, without the
proportional split's mean-field approximation: a true case that dies
before its test returns is never counted as confirmed. Returns the
length-`n` confirmed-in-care sub-stock. The caller forms the suspect
sub-stock as the demand remainder `D − O_conf`.
"""
function two_clock_confirmed(A_bvd::AbstractVector, conf_hazard::AbstractVector,
        S_clin::AbstractVector)
    n = length(A_bvd)
    T = promote_type(eltype(A_bvd), eltype(conf_hazard), eltype(S_clin))
    O_conf = Vector{T}(undef, n)
    L = length(S_clin)
    one_T = one(T)
    @inbounds for t in 1:n
        acc = zero(T)
        ## `prod_unconf` = probability cohort `u` is still unconfirmed at day
        ## `t`, `∏_{j=u+1}^{t}(1 − hazard_j)`, extended one factor per step as
        ## `u` walks back from `t`.
        prod_unconf = one_T
        for u in t:-1:1
            d = t - u
            s = d < L ? S_clin[d + 1] : zero(T)
            cdf_conf = one_T - prod_unconf
            acc += A_bvd[u] * cdf_conf * s
            prod_unconf *= (one_T - conf_hazard[u])
        end
        O_conf[t] = acc
    end
    return O_conf
end

"""
    admission_headroom(adm_days, capacity_history, isolation_history)

Per-admission-day right-censoring bound for the admissions likelihood, built
from data alone, the admissions analogue of [`censoring_cap`](@ref). Each
admission day takes the recorded free-bed headroom `C_fix(t) − occupancy(t-1)`:
the nearest recorded implied-capacity value (`capacity_history`) less the
previous day's observed occupancy (`isolation_history`, forward-filled). Days
before any occupancy record and days with no recorded capacity get a large
no-op bound. The bound is kept strictly above each day's observed admissions
(floored at `obs + 0.5`), so a day whose admissions exceeded the headroom is
scored as a plain count rather than on the censoring boundary. Returns a
length-`length(adm_days)` `Float64` vector.
"""
function admission_headroom(adm_days, adm_obs, capacity_history,
        isolation_history)
    cdays = Int.(capacity_history.days)
    ccounts = Float64.(capacity_history.counts)
    idays = Int.(isolation_history.days)
    icounts = Float64.(isolation_history.counts)
    nocap = 1.0e6
    have_cap = !isempty(ccounts)
    have_occ = !isempty(icounts)
    m = length(adm_days)
    head = Vector{Float64}(undef, m)
    @inbounds for (i, d) in enumerate(adm_days)
        di = Int(d)
        cap = if have_cap
            j = argmin(abs.(cdays .- di))
            ccounts[j]
        else
            nocap
        end
        ## Previous day's observed occupancy, the most recent record strictly
        ## before this admission day. No prior record ⇒ a large no-op headroom.
        prev_occ = if have_occ
            prior = findall(<(di), idays)
            isempty(prior) ? -nocap : icounts[prior[argmax(idays[prior])]]
        else
            -nocap
        end
        h = cap - prev_occ
        o = adm_obs === missing ? 0.0 : Float64(adm_obs[i])
        head[i] = max(h, o + 0.5)
    end
    return head
end

"""
DRC treatment-centre patient-flow likelihood. Occupancy is built as a running
balance of latent admission, discharge and abscond events rather than a
convolution. The "Patients en isolement" figure is the daily occupied-bed
count, fitted as latent bed demand right-censored at the implied-capacity
bound.

Admissions are the reported suspects ([`reported_cases_model`](@ref))
carried through a short suspected→admission delay and split into a BVD
true-case inflow `A_bvd = p_iso_bvd · p_drc · bvd_reports` at the
severity-skewed rate `p_iso_bvd` ([`isolation_severity_model`](@ref)) and
a non-BVD inflow `A_bg = p_iso · bg_daily` at the base rate
([`isolation_admission_model`](@ref)).

The clinical layer sets the total occupancy and discharge flows and is
label-independent, so a true case can die before confirmation. A BVD true case
leaves by death (weight `CFR_iso`) on the admission→death stay or recovery
(weight `1 − CFR_iso`) on the admission→recovery stay. A non-BVD admission
rules out on the rule-out stay. Absconds drain the suspect pool at
`κ · O_susp(t−1)`. The total latent demand is the running balance

```math
D(t) = D(t-1) + A(t) - \\text{deaths}(t) - \\text{recoveries}(t)
       - \\text{rule-outs}(t) - \\text{absconds}(t),
```

with `A = A_bvd + A_bg`. The in-care deaths flow is suspect and confirmed
combined (the Tableau 6 `décédés` row), scored directly and not gated by
confirmation. The in-care fatality `CFR_iso = logistic(logit(CFR) + β_iso)`
adjusts the infection CFR to the admitted population by a sampled modifier
`β_iso`. It is conditional on admission, not a causal treatment effect.

A label overlay carves the occupied stock into confirmed and suspect
sub-stocks without removing anyone from the total. Confirmation relabels
an occupied true case at the daily hazard `ρ · τ_test · p_pos`, the
community hazard `τ_test · p_pos` borrowed from the lab pipeline
([`confirmed_cases_model`](@ref)) scaled by a sampled in-care modifier
`ρ = exp(γ_conf)`. In-care confirmation runs slower than the community
rate (`ρ < 1`), the confirmed/suspect census identifying the modifier.
The confirmed sub-stock is the confirmed subset of the occupied BVD
stock and the suspect sub-stock the remainder
`O_susp(t) = D(t) − O_conf(t)`, so a case that dies before confirmation
is a suspect death in the combined deaths flow.

Capacity enters only as a fixed, data-derived censoring bound. The
latent demand is uncapped. The occupancy likelihood is a
NegativeBinomial around the demand, right-censored at the recorded
implied-capacity series ([`censoring_cap`](@ref),
[`censored_occupancy_model`](@ref)). Admissions are censored at the
recorded free-bed headroom `C_fix(t) − occupancy(t-1)`. The capacity
walk `C(t)` ([`bed_capacity_walk_model`](@ref)) carries the
implied-capacity likelihood, and unmet demand
`D(t) − censored occupancy` is a returned diagnostic.

The daily in-care outcome flows (deaths, rule-outs, admissions,
absconds) are each an optional NegativeBinomial stream scored against
the Tableau 6 patient-movement counts, a no-op when their history is
empty. The two sub-stocks are scored against the Tableau 6 census
breakdown (`dont confirmés` / `dont suspects`) where published, in place
of the total on those days. An opt-in occupancy offset Δ(t) on the
supplied `occupancy_break_days` ([`cumulative_occupancy_offset`](@ref))
absorbs a between-report measurement-basis discontinuity in the
isolation series. Empty (the default) is a no-op.

Exposes the cut-off occupancy, bed demand and shortfall, the utilisation, the
BVD share of demand, `CFR_iso` and `β_iso`, the length-of-stay, the two
sub-stock prevalences and their stay means, the in-care fraction and the daily
series for forecasting and replication.
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
        ## Tableau 6 occupancy split (`dont confirmés` / `dont suspects`): two
        ## census series scored in place of the total occupancy on the days they
        ## are present, only when the confirmation hazard is non-zero. Empty by
        ## default, leaving the total-occupancy likelihood alone.
        confirmed_incare_history = (; days = Int[], counts = Int[]),
        suspect_incare_history = (; days = Int[], counts = Int[]),
        ## Daily in-care confirmation hazard `τ_test · p_pos` borrowed from
        ## the lab pipeline ([`confirmed_cases_model`](@ref)). `nothing`
        ## (standalone, no lab stream) gives a zero hazard, so the confirmed
        ## sub-stock stays empty and the suspect sub-stock carries the whole
        ## occupancy.
        conf_hazard_daily::Union{Nothing, AbstractVector} = nothing,
        admission = isolation_admission_model(),
        severity = isolation_severity_model(),
        capacity = bed_capacity_walk_model,
        dispersion = surveillance_dispersion_model(),
        ## Occupancy / flow dispersion can be injected from the joint composer's
        ## pooled set (`k_external`). Standalone it samples its own.
        k_external::Union{Nothing, Real} = nothing,
        ## In-care fatality modifier prior: β_iso on the infection CFR.
        cfr_modifier_prior = Normal(0.0, 0.5),
        ## Small abscond / loss-to-follow-up fraction of occupancy per day.
        abscond_prior = truncated(Normal(0.01, 0.01); lower = 0),
        ## In-care confirmation-rate modifier prior (log scale): γ_conf scales
        ## the borrowed community hazard to the effective in-care rate
        ## ρ = exp(γ_conf). Centred on zero (ρ = 1) so the census sets the
        ## split.
        incare_confirm_log_prior = Normal(0.0, 0.5),
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
            sd_prior = truncated(Normal(4.0, 1.5); lower = 1)),
        ## Opt-in occupancy reclassification-break days (grid indices). A level
        ## step is fitted into the modelled total at each, absorbing a
        ## measurement-basis discontinuity in the isolation series. Empty (the
        ## default) is a no-op. See `cumulative_occupancy_offset`.
        occupancy_break_days::AbstractVector{<:Integer} = Int[],
        ## Prior sd of each occupancy break step (beds), centred on zero.
        occupancy_break_sd::Real = 25.0)
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
    ## (Tableau 6 décédés) relative to admissions and occupancy. β_iso < 0
    ## means treatment lowers the in-care fatality below the infection CFR.
    ## This is a conditional-on-admission (in-care) fatality, not a causal
    ## treatment effect. The recovered-among-confirmed ("cumul guéris")
    ## stream is modelled separately off the confirmed cases
    ## ([`recovered_model`](@ref)): guéris is the confirmed-and-discharged
    ## subset, not all in-care recoveries.
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
    ## Abscond (loss-to-follow-up) drains the suspect pool at a small daily
    ## fraction κ of the previous-day suspect occupancy.
    abscond_frac ~ abscond_prior
    κ = abscond_frac

    ## Opt-in cumulative occupancy reclassification-break offset Δ(t) added to
    ## the modelled total occupancy, absorbing a between-report measurement-
    ## basis discontinuity in the observed isolation series (the `au-lit-J-1`
    ## versus `Fin-J` reclassification) on the manually supplied
    ## `occupancy_break_days`. A level step is fitted at each break day, sampled
    ## non-centred and centred on zero, so the fit partitions it into reporting
    ## artifact vs real demand. Applied additively to the modelled total only.
    ## Demand (the diagnostic) stays the un-offset latent stock. Empty (the
    ## default) → no sampled step, Δ = 0, a no-op. Only break days on or before
    ## an observed occupancy day can move the likelihood. Later ones are dropped
    ## so no inert step is sampled. See `cumulative_occupancy_offset`.
    iso_last = isempty(isolation_history.days) ? 0 :
               maximum(Int.(isolation_history.days))
    brk_days = [Int(d) for d in occupancy_break_days if Int(d) <= iso_last]
    if isempty(brk_days)
        b = Float64[]
    else
        occupancy_step ~ product_distribution(
            fill(Normal(0, 1), length(brk_days)))
        b = occupancy_break_sd .* occupancy_step
    end
    ## Cumulative reclassification offset Δ over the 1:n grid.
    break_grid_days = brk_days
    occ_break_offset = cumulative_occupancy_offset(1:n, break_grid_days, b)

    ## Admission inflow through the suspected→admission delay, split into BVD
    ## true-case (`p_iso_bvd`) and non-BVD (`p_iso`) inflows. Uncapped latent
    ## demand. Capacity enters only as a censoring bound below.
    A_bvd = convolve_delay(p_iso_bvd .* p_drc .* bvd_reports_daily,
        adm_delay_state.pmf)
    A_bg = convolve_delay(p_iso .* bg_daily, adm_delay_state.pmf)
    if eltype(A_bvd) === Any
        A_bvd = convert(Vector{eltype(C)}, A_bvd)
        A_bg = convert(Vector{eltype(C)}, A_bg)
    end

    ## Label-independent clinical discharge events. Deaths and recoveries split
    ## `A_bvd` by `CFR_iso`. Rule-outs discharge `A_bg`.
    dpmf = death_los_state.pmf
    rpmf = recovery_los_state.pmf
    deaths_daily = convolve_delay(CFR_iso .* A_bvd, dpmf)
    recover_daily = convolve_delay((one(CFR_iso) - CFR_iso) .* A_bvd, rpmf)
    ruleout_daily = convolve_delay(A_bg, ruleout_los_state.pmf)
    admit_daily = A_bvd .+ A_bg

    ## Community confirmation hazard `τ_test · p_pos` borrowed from the lab
    ## pipeline. `nothing` (standalone) gives a zero hazard.
    borrowed_hazard = if conf_hazard_daily === nothing
        zeros(eltype(A_bvd), n)
    elseif eltype(conf_hazard_daily) === Any
        convert(Vector{eltype(C)}, conf_hazard_daily)
    else
        conf_hazard_daily
    end

    ## The split census is identified only when the borrowed hazard is
    ## non-zero (a lab stream supplies it). With a structural zero the
    ## confirmed sub-stock is empty, so the split likelihood no-ops and
    ## those days stay on the total.
    split_active = any(>(zero(eltype(borrowed_hazard))), borrowed_hazard)

    ## In-care confirmation-rate modifier ρ = exp(γ_conf) on the borrowed
    ## hazard. Sampled only when the hazard is non-zero, so no unidentified
    ## dimension is added when the split is absent (matching the
    ## `k_external` pattern).
    if split_active
        incare_confirm_log ~ incare_confirm_log_prior
    else
        incare_confirm_log = zero(eltype(borrowed_hazard))
    end
    ρ_conf = exp(incare_confirm_log)
    conf_hazard = ρ_conf .* borrowed_hazard

    ## Forward running-balance occupancy: total demand and the BVD/non-case
    ## stocks. The scored abscond flow is recomputed below off the two-clock
    ## suspect stock.
    acc = accumulate_occupancy(A_bvd, A_bg, deaths_daily, recover_daily,
        ruleout_daily, κ, conf_hazard)
    demand = acc.demand
    O_bvd = acc.O_bvd

    ## Two-clock confirmed-in-care sub-stock: cohort-tracked
    ## confirmed-and-present prevalence, exact in the fast-death tail where
    ## the running balance's proportional drain is only mean-field. Demand
    ## and `O_bvd` stay as `accumulate_occupancy` built them. `O_conf` (and
    ## `O_susp = D − O_conf`) is replaced.
    S_clin = clinical_stay_survival(dpmf, rpmf, CFR_iso)
    O_conf_raw = two_clock_confirmed(A_bvd, conf_hazard, S_clin)
    ## `O_conf ≤ O_bvd` holds by construction. Clamp into `[0, O_bvd]` as a
    ## guard under any prior draw. The suspect sub-stock is the demand
    ## remainder.
    O_conf = map((c, b) -> clamp(c, zero(eltype(demand)), b), O_conf_raw, O_bvd)
    O_susp = map((d, c) -> max(d - c, zero(eltype(demand))), demand, O_conf)
    ## Abscond outflow off the two-clock suspect stock, `κ · O_susp(t-1)`.
    ## Day 1 has no prior stock.
    abscond_daily = [t == 1 ? zero(eltype(demand)) : κ * O_susp[t - 1]
                     for t in 1:n]
    if eltype(demand) === Any
        demand = convert(Vector{eltype(C)}, demand)
        O_conf = convert(Vector{eltype(C)}, O_conf)
        O_susp = convert(Vector{eltype(C)}, O_susp)
        abscond_daily = convert(Vector{eltype(C)}, abscond_daily)
    end

    ## Add the reclassification offset Δ(t) to the modelled total only.
    ## Demand (the diagnostic) stays the un-offset latent stock.
    occ_offset = eltype(occ_break_offset) === Any ?
                 convert(Vector{eltype(C)}, occ_break_offset) : occ_break_offset
    ## Broadcasts, not `map(1:n) do t`: the `do`-block builds an anonymous
    ## closure over `demand`/`occ_offset`/`O_conf` whose reverse-mode shadow
    ## Enzyme cannot construct. An elementwise broadcast creates no closure
    ## and is bit-identical under Mooncake.
    occ_obs_total = demand .+ occ_offset

    ## Confirmed and suspect census means, summing to the offset total.
    conf_split = copy(O_conf)
    susp_split = max.(occ_obs_total .- conf_split, zero(eltype(demand)))

    ## Occupancy likelihood: NegativeBinomial around the latent demand,
    ## right-censored at the implied-capacity bound. Days with a published split
    ## are scored as the two sub-stocks below, not the total, so the total and
    ## its parts are never both scored on one day. When the overlay is
    ## unidentified the split is unscored and those days stay on the total.
    split_days = split_active ? Set(Int.(confirmed_incare_history.days)) :
                 Set{Int}()
    iso_all_days = isolation_history.days
    iso_keep = [!(Int(d) in split_days) for d in iso_all_days]
    iso_days = iso_all_days[iso_keep]
    iso_obs = isempty(isolation_history.counts) ? missing :
              collect(Int.(isolation_history.counts))[iso_keep]
    iso_means = [occ_obs_total[clamp(Int(d), 1, n)] for d in iso_days]
    iso_ceil = censoring_cap(iso_days, iso_obs, capacity_history)
    isolation ~ to_submodel(
        censored_occupancy_model(iso_means, iso_ceil, iso_obs, k))

    ## Capacity likelihood: the implied bed count is a noisy observation of
    ## C(t).
    cap_days = capacity_history.days
    cap_modelled = [C[clamp(Int(d), 1, n)] for d in cap_days]
    cap_obs = isempty(capacity_history.counts) ? missing :
              collect(Int.(capacity_history.counts))
    bed_capacity ~ to_submodel(
        vintage_increments_model(cap_modelled, cap_obs, k))

    ## Split likelihoods: on days with a published split, score the two
    ## sub-stock censuses instead of the total. Guarded by `split_active`,
    ## so they no-op when the hazard is structurally zero.
    ci_days = split_active ? confirmed_incare_history.days : Int[]
    ci_obs = (!split_active || isempty(confirmed_incare_history.counts)) ?
             missing : collect(Int.(confirmed_incare_history.counts))
    confirmed_incare_obs ~ to_submodel(vintage_increments_model(
        [conf_split[clamp(Int(d), 1, n)] for d in ci_days], ci_obs, k))
    si_days = split_active ? suspect_incare_history.days : Int[]
    si_obs = (!split_active || isempty(suspect_incare_history.counts)) ?
             missing : collect(Int.(suspect_incare_history.counts))
    suspect_incare_obs ~ to_submodel(vintage_increments_model(
        [max(susp_split[clamp(Int(d), 1, n)], zero(eltype(susp_split)))
         for d in si_days], si_obs, k))
    ## Confirmed-in-care deaths, attributed as the death flow times the
    ## confirmed share of the BVD stock. Exposed for the report, not
    ## separately scored.
    bvd_stock = acc.O_bvd
    conf_share = [bvd_stock[t] > zero(eltype(demand)) ?
                  O_conf[t] / bvd_stock[t] : zero(eltype(demand)) for t in 1:n]
    confirmed_incare_deaths_daily = deaths_daily .* conf_share

    ## Optional daily Tableau 6 flow likelihoods, each a no-op on empty history.
    ## Admissions are right-censored at the recorded free-bed headroom.
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
    adm_means = [admit_daily[clamp(Int(d), 1, n)] for d in adm_h_days]
    adm_ceil = admission_headroom(adm_h_days, adm_h_obs, capacity_history,
        isolation_history)
    admissions ~ to_submodel(
        censored_occupancy_model(adm_means, adm_ceil, adm_h_obs, k))
    ab_days = absconded_history.days
    ab_obs = isempty(absconded_history.counts) ? missing :
             collect(Int.(absconded_history.counts))
    absconded ~ to_submodel(vintage_increments_model(
        [abscond_daily[clamp(Int(d), 1, n)] for d in ab_days], ab_obs, k))

    ## Cut-off reported quantities. Occupancy is the censored stock — capacity
    ## bounds it at the cut-off, so the reported occupancy is the demand capped
    ## at the fixed bed count. Bed demand is the uncapped latent stock.
    z0 = zero(eltype(C))
    dem_T = isempty(demand) ? z0 : demand[end]
    occ_T = min(dem_T, C_T)
    overall_los = CFR_iso * death_los_state.mean +
                  (one(CFR_iso) - CFR_iso) * recovery_los_state.mean
    expected_isolation := safe_rate(occ_T)
    expected_bed_demand := safe_rate(dem_T)
    ## Cut-off daily flows: the end-of-grid value of each modelled daily
    ## series, the one-week-ahead forecast base for admissions, in-care
    ## deaths and rule-outs.
    expected_admissions := safe_rate(isempty(admit_daily) ? z0 :
                                     admit_daily[end])
    expected_incare_deaths := safe_rate(isempty(deaths_daily) ? z0 :
                                        deaths_daily[end])
    expected_ruleouts := safe_rate(isempty(ruleout_daily) ? z0 :
                                   ruleout_daily[end])
    ## Unmet demand: the uncapped demand above the censored occupancy.
    bed_shortfall := safe_rate(max(dem_T - occ_T, z0))
    bed_utilisation := safe_rate(occ_T) / safe_rate(C_T)
    isolation_severity := sev_state.δ_iso
    isolation_bvd_admission := p_iso_bvd
    incare_cfr := CFR_iso
    incare_cfr_modifier := β_iso
    treatment_overall_los := overall_los
    ## Cut-off occupancy split: the two sub-stock prevalences at the grid end,
    ## and the confirmed share of the occupied stock.
    conf_incare_T = isempty(conf_split) ? z0 : conf_split[end]
    susp_incare_T = isempty(susp_split) ? z0 : max(susp_split[end], z0)
    expected_confirmed_incare := safe_rate(conf_incare_T)
    expected_suspect_incare := safe_rate(susp_incare_T)
    incare_confirmed_share := safe_rate(conf_incare_T) / safe_rate(dem_T)
    ## In-care confirmation-rate modifier ρ (raw, can exceed one, so reported
    ## directly). ρ < 1 means occupied suspects are confirmed slower than the
    ## borrowed community hazard, held for repeated exclusion testing.
    incare_confirm_modifier := ρ_conf
    ## Cut-off cumulative occupancy reclassification offset (fitted, can be
    ## negative, so reported raw rather than through `safe_rate`). Reports how
    ## much of the observed reclassification the model absorbed as a reporting
    ## artifact, the rest carried by real demand. The per-day steps `b` and the
    ## grid offset `occ_break_offset` are returned for the report.
    occupancy_break := isempty(occ_break_offset) ? z0 : last(occ_break_offset)

    return (; p_iso, p_iso_bvd, δ_iso = sev_state.δ_iso,
        CFR_iso, β_iso, capacity = C_T,
        death_los_mean = death_los_state.mean,
        recovery_los_mean = recovery_los_state.mean,
        ruleout_los_mean = ruleout_los_state.mean,
        admission_delay_mean = adm_delay_state.mean,
        overall_los, abscond_frac, k_isolation = k,
        demand, occupancy = min.(demand, C), isolation,
        deaths_daily, recover_daily, ruleout_daily, admit_daily,
        abscond_daily,
        break_steps = b, break_offset = occ_break_offset,
        break_grid_days,
        occupancy_break = isempty(occ_break_offset) ? z0 :
                          last(occ_break_offset),
        confirmed_incare = conf_split, suspect_incare = susp_split,
        confirmed_incare_deaths_daily, incare_confirm_modifier = ρ_conf,
        expected_confirmed_incare = safe_rate(conf_incare_T),
        expected_suspect_incare = safe_rate(susp_incare_T),
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
(scaled-convolution) stream: the renewal analogue of the convolution-and-
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
before it is recorded as recovered: the report's "cumul guéris" counts
recoveries among confirmed cases, so the recovery follows confirmation by
the confirmation-to-recovery delay. In principle a positive result could
return after a patient has already recovered and been discharged. The
reported total is assumed to reflect confirmed cases carefully recorded
as recovered, so the confirmation-then-recovery ordering holds.

The cumulative recovered series ends at the cut-off, so its between-vintage
increments are fitted as observed `~` data with a NegativeBinomial whose
dispersion is sampled here, not shared with the other streams (the recovered
signal has its own observation noise). The convolution right-censors
recoveries that have not yet resolved by the cut-off, so a small observed
recovered count is consistent with a high eventual survival fraction and a
long recovery delay. Empty by default. A `missing` cut-off total leaves the
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
        ## Confirmation-to-recovery (discharge) delay. An Ebola survivor is
        ## discharged a couple of weeks after confirmation, so the default is
        ## a mean ~14 d stay before recovery is recorded.
        confirmation_to_recovery = censored_delay_model(
            cdf_nmax(lognormal_meansd(14.0, 8.0); q = 0.99);
            mean_prior = truncated(Normal(14.0, 5.0); lower = 1),
            sd_prior = truncated(Normal(8.0, 4.0); lower = 1)),
        ## Dispersion can be injected from the joint composer's pooled set
        ## (`k_external`). Standalone it samples its own from `dispersion`.
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

# Symptom-onset reporting-triangle observation model. The digitised
# reporting triangle (`data/onset_curve_scanned.csv`, loaded by
# `load_onset_curve`, `src/onset_curve.jl`) is the only direct observation
# of the shared latent onset series: every other stream sees onsets only
# after a further convolution (suspected via onset-to-report, deaths via
# onset-to-death, laboratory via onset-to-report ⊕ receipt). A discrete
# reporting-delay hazard, nonparametric over the delay and drifting over
# calendar time, is fitted to the between-vintage increments of the
# triangle (not its levels, see `onset_reporting_model`), so ascertainment
# and delay speed are read off the same fitted hazard rather than a
# separate multiplicative factor.

"""
    safe_studentt(μ, σ, ν)

NaN / Inf-safe location-scale Student-t distribution `μ + σ · Tν`, built
from `Distributions.TDist(ν)` via the affine-combination operators. `σ` is
floored away from zero and non-finite values, mirroring
[`safe_nbinomial`](@ref)'s domain-guard idiom so a wild NUTS warmup proposal
cannot throw a `DomainError` and abort a Mooncake gradient. A non-positive
or non-finite `ν` falls back to `4`, the caller's own default, rather than
to the smallest value `TDist` accepts: `TDist(1)` is Cauchy, so a floor at
the domain edge would quietly turn a bad degrees-of-freedom argument into a
likelihood with no mean or variance, which is a worse failure than the
`DomainError` it replaces.
Used by [`onset_reporting_model`](@ref) to score the reporting-triangle
increments, which are frequently negative (a later scan reads fewer cases
at some onset date than an earlier one, from digitisation noise rather than
a real reporting reversal) under a heavy-tailed likelihood a count
distribution cannot represent.
"""
function safe_studentt(μ::Real, σ::Real, ν::Real)
    σc = (isfinite(σ) && σ > zero(σ)) ? σ : eps(typeof(float(σ)))
    νc = (isfinite(ν) && ν > zero(ν)) ? ν : oftype(float(ν), 4)
    return μ + σc * TDist(νc)
end

"""
    onset_increments_model(means, sds, increments, ν)

Heavy-tailed likelihood for the reporting-triangle increment cells: cell `i`
is Student-t about the modelled increment `means[i]` with scale `sds[i]` and
fixed degrees of freedom `ν` (see [`safe_studentt`](@ref)).

`increments` is a positional model argument, not a field read out of a
container inside the model body, and that is load bearing rather than
stylistic: DynamicPPL decides observe-versus-assume by checking whether the
tilde's symbol appears in the enclosing model's argument names
(`DynamicPPL.inargnames`), so a local variable on the left of `~` is treated
as a latent quantity, silently turning every observation into a sampled
parameter and dropping the likelihood. Same argument shape as
[`vintage_increments_model`](@ref). A `missing` argument samples
instead (the predictive-generator path).
"""
@model function onset_increments_model(means::AbstractVector,
        sds::AbstractVector,
        increments::Union{Missing, AbstractVector{<:Real}}, ν::Real)
    n = length(means)
    if ismissing(increments)
        increments = Vector{Union{Missing, Float64}}(missing, n)
    end
    for i in 1:n
        increments[i] ~ safe_studentt(means[i], sds[i], ν)
    end
    return (; means, sds, increments)
end

"""
    onset_report_cdf(δ, logit_h0, γ, u, grid_start)

Un-normalised reported proportion of onset date `u`'s eventual symptom-onset
cases reported within `δ` days, under the discrete hazard `h(d, t) =
logistic(logit_h0[d+1] + γ[t - grid_start + 1])` (delay `d`, report calendar
day `t = u + d`, indexed into the calendar-time random-walk vector `γ` which
starts at grid day `grid_start`):

```math
\\text{cdf}(u, \\delta) = \\begin{cases} 0 & \\delta < 0 \\\\
    1 - \\prod_{j=0}^{\\min(\\delta, D-1)} \\bigl(1 - h(j, u + j)\\bigr)
    & \\delta \\ge 0 \\end{cases}, \\qquad D = \\text{length}(logit\\_h0).
```

This is the delay shape only: it is not the reported proportion `F` the
model scores, which also carries ascertainment (see [`onset_report_G`](@ref)
and [`onset_report_F`](@ref)). `δ < 0` returns exactly `0`, the right-
truncation case, and is the building block both `G` and the delay-weighted
[`onset_report_anchor`](@ref) rely on. `δ` is capped at `D-1`, so the
returned value is constant for every `δ >= D-1`.

Pure, top-level and allocation-free (a single indexed `@inbounds` loop, no
closures), so it AD-transparently composes into the vectorised moment and
total functions below without the Enzyme boxing failure mode a `map(...)
do` or a closure over `logit_h0`/`γ` inside an `@model` body can hit (see
`two_clock_confirmed` for the same pattern).
"""
function onset_report_cdf(δ::Integer, logit_h0::AbstractVector,
        γ::AbstractVector, u::Integer, grid_start::Integer)
    T = promote_type(eltype(logit_h0), eltype(γ))
    δ < 0 && return zero(T)
    D = length(logit_h0)
    δc = min(Int(δ), D - 1)
    surv = one(T)
    @inbounds for j in 0:δc
        gi = u + j - grid_start + 1
        h = logistic(logit_h0[j + 1] + γ[gi])
        surv *= (one(T) - h)
    end
    return one(T) - surv
end

"""
    onset_report_G(δ, logit_h0, γ, u, grid_start)

Normalised delay CDF `G(u, δ) = cdf(u, δ) / cdf(u, D-1)`
([`onset_report_cdf_extrapolated`](@ref) supplying both terms), so
`G(u, D-1) = 1` by construction and `G` is a proper delay distribution
rather than an asymptote that drifts with the hazard level. Built on the
extrapolated cdf (calendar index clamped rather than assumed in-range)
because the denominator always reaches `D - 1` days ahead of `u`, which for
an onset date within `D - 1` days of `grid_end` runs past `γ`'s fitted
support even though every scored `δ` itself stays in range. The calendar
effect is held flat at the walk's nearest edge there, same as everywhere
else this package extrapolates a fitted walk. Agrees exactly with the
in-range calculation whenever every index touched is already in-range. The
denominator is guarded with [`safe_rate`](@ref) so it is never divided by
zero. `G(u, δ) = 0` for `δ < 0`, since the numerator already returns `0`
there.

`G(u, D-1) = 1` holds except in one corner. The numerator never exceeds
the denominator, so the guard can only shrink the ratio, and if every
`h(j, ·)` is small enough that each `1 - h` rounds to `1` then both terms
underflow to exactly `0` and `G` is `0` rather than `1`. That needs all
`D` baseline hazards below about `1e-16` at once, a joint excursion the
prior puts roughly 30 standard deviations away, and it degrades quietly
rather than producing `NaN`. [`onset_report_ascertainment`](@ref) clamps
against the one consequence that would spread, an anchor of exactly zero.
Pure, top-level, allocation-free.
"""
function onset_report_G(δ::Integer, logit_h0::AbstractVector,
        γ::AbstractVector, u::Integer, grid_start::Integer)
    D = length(logit_h0)
    num = onset_report_cdf_extrapolated(δ, logit_h0, γ, u, grid_start)
    den = onset_report_cdf_extrapolated(D - 1, logit_h0, γ, u, grid_start)
    return num / safe_rate(den)
end

"""
    onset_report_F(δ, logit_h0, γ, u, grid_start, α)

Cumulative reported proportion `F(u, δ) = α * G(u, δ)` of onset date `u`'s
eventual symptom-onset cases reported within `δ` days, with `α` the
ascertainment level for onset date `u` ([`onset_report_ascertainment`](@ref))
and `G` the normalised delay CDF ([`onset_report_G`](@ref)). `F(u, D-1) = α`
exactly, so the hazard's delay shape and its ascertainment level are two
separate factors rather than one asymptote. Pure, top-level, allocation-free.
"""
function onset_report_F(δ::Integer, logit_h0::AbstractVector,
        γ::AbstractVector, u::Integer, grid_start::Integer, α::Real)
    return α * onset_report_G(δ, logit_h0, γ, u, grid_start)
end

"""
    onset_nowcast(observed, onsets_u, δ, logit_h0, γ, u, grid_start, α)

Eventual reported count for onset date `u`, given the count `observed`
already printed for it at reporting delay `δ`:

```math
N(u) = y(u, \\delta) + n_u \\, \\bigl(\\alpha - F(u, \\delta)\\bigr),
```

with `n_u` the modelled symptom onsets on that date, `α` its ascertainment
level and `F` the cumulative reported proportion ([`onset_report_F`](@ref)).
The first term is what the latest figure carries and the second is what the
model expects still to arrive, so the estimate is anchored on the data
rather than on the fitted curve alone.

This is the nowcast the report plots, and it behaves differently from the
unconditional expectation `n_u α` at the two ends of the delay axis. At a
long delay `F(u, δ)` has reached `α`, the correction term is zero and the
nowcast is the observed count exactly, with no posterior spread left. At a
short delay most of `α` is still outstanding, so the correction carries the
posterior's full uncertainty about `n_u` and the delay curve. The interval
therefore widens towards the most recent onset dates and closes on the data
going back from them, which the unconditional expectation does neither.

`F(u, δ) ≤ α` for every `δ` ([`onset_report_G`](@ref) is a normalised CDF),
so the correction is never negative and the nowcast never falls below the
count already reported. Pure, top-level, allocation-free.
"""
function onset_nowcast(observed::Real, onsets_u::Real, δ::Integer,
        logit_h0::AbstractVector, γ::AbstractVector, u::Integer,
        grid_start::Integer, α::Real)
    outstanding = α - onset_report_F(δ, logit_h0, γ, u, grid_start, α)
    return observed + onsets_u * outstanding
end

"""
    onset_report_anchor(logit_h0, γ, u, grid_start, a)

Delay-weighted average of the calendar-indexed daily ascertainment series
`a` over onset date `u`'s reporting window, `anchor(u) = Σ_d g(u, d) *
a[clamp(u + d, 1, length(a))]`, with `g(u, d) = G(u, d) - G(u, d-1)`
([`onset_report_G`](@ref)) the normalised delay PMF, `d = 0 … D-1`. Since
`Σ_d g(u, d) = 1`, `anchor(u)` is a genuine weighted average of `a` and so
lies within `[minimum(a), maximum(a)]`, and a constant `a` gives back that
constant exactly. Both hold wherever `G` reaches one, so they inherit the
single underflow corner [`onset_report_G`](@ref) documents, where the
weights sum to zero and the anchor with them. `a` is indexed on the
calendar/report axis and clamped at
both ends, which is what reconciles the onset-indexed `anchor(u)` with a
report-indexed series such as the confirmed pipeline's daily ascertainment
(see [`onset_reporting_model`](@ref)). Pure, top-level, allocation-free.
"""
function onset_report_anchor(logit_h0::AbstractVector, γ::AbstractVector,
        u::Integer, grid_start::Integer, a::AbstractVector)
    D = length(logit_h0)
    T = promote_type(eltype(logit_h0), eltype(γ), eltype(a))
    na = length(a)
    ng = length(γ)
    ## `G(u, d)` shares one survival product across every `d`, and its
    ## denominator does not depend on `d`, so the whole weighted sum runs in
    ## a single pass rather than rebuilding `G` from scratch per delay.
    invden = inv(safe_rate(onset_report_cdf_extrapolated(
        D - 1, logit_h0, γ, u, grid_start)))
    surv = one(T)
    g_prev = zero(T)
    acc = zero(T)
    @inbounds for d in 0:(D - 1)
        gi = clamp(u + d - Int(grid_start) + 1, 1, ng)
        surv *= (one(T) - logistic(logit_h0[d + 1] + γ[gi]))
        g_cur = (one(T) - surv) * invden
        acc += (g_cur - g_prev) * a[clamp(u + d, 1, na)]
        g_prev = g_cur
    end
    return acc
end

"""
    onset_report_anchor_series(logit_h0, γ, grid_start, grid_end, a)

[`onset_report_anchor`](@ref) evaluated for every onset date `u in
grid_start:grid_end`, the onset-date grid the ascertainment walk spans (see
[`onset_reporting_model`](@ref)). Returns an empty vector when `grid_end <
grid_start`. Pure, top-level, single indexed loop.
"""
function onset_report_anchor_series(logit_h0::AbstractVector,
        γ::AbstractVector, grid_start::Integer, grid_end::Integer,
        a::AbstractVector)
    lo = Int(grid_start)
    hi = Int(grid_end)
    T = promote_type(eltype(logit_h0), eltype(γ), eltype(a))
    hi < lo && return T[]
    out = Vector{T}(undef, hi - lo + 1)
    @inbounds for (k, u) in enumerate(lo:hi)
        out[k] = onset_report_anchor(logit_h0, γ, u, grid_start, a)
    end
    return out
end

"""
    onset_report_moments(onsets, logit_h0, γ, grid_start, alpha, onset_idx,
        cur_report_idx, prev_report_idx)

Per-cell modelled increment mean and the two modelled cumulative levels it
is the difference of, for a batch of `(onset day, current report day,
previous report day)` triples (the increment cells [`load_onset_curve`](@ref)
builds). For cell `i`,

```math
\\ell_{\\text{cur},i} = \\text{onsets}[u_i] \\cdot
    F(u_i, \\delta_{\\text{cur},i}), \\qquad
\\ell_{\\text{prev},i} = \\text{onsets}[u_i] \\cdot
    F(u_i, \\delta_{\\text{prev},i}), \\qquad
\\text{mean}_i = \\ell_{\\text{cur},i} - \\ell_{\\text{prev},i},
```

with `δ = report_idx - u_i` and [`onset_report_F`](@ref) supplying `F`,
`alpha` indexed at onset date `u_i` (clamped into `1:length(alpha)`, the
`grid_start:grid_end` ascertainment grid) for its `α` argument (so
`δ_prev < 0` at the sentinel `prev_report_idx = 0`, the virtual empty
predecessor for the very first scored vintage, contributes `ℓ_prev = 0` with
no special-casing). An `onset_idx` outside `1:length(onsets)` contributes a
zero rate rather than indexing out of bounds. Returns `(; means, level_cur,
level_prev)`, each a length-`length(onset_idx)` vector. `level_cur`/
`level_prev` are reused by [`onset_report_scales`](@ref) to grow the
observation scale with the modelled (not observed) magnitude. Pure,
top-level, single indexed loop (see [`onset_report_cdf`](@ref) for the
AD-safety rationale).
"""
function onset_report_moments(onsets::AbstractVector,
        logit_h0::AbstractVector, γ::AbstractVector, grid_start::Integer,
        alpha::AbstractVector,
        onset_idx::AbstractVector{<:Integer},
        cur_report_idx::AbstractVector{<:Integer},
        prev_report_idx::AbstractVector{<:Integer})
    m = length(onset_idx)
    T = promote_type(eltype(onsets), eltype(logit_h0), eltype(γ),
        eltype(alpha))
    means = Vector{T}(undef, m)
    level_cur = Vector{T}(undef, m)
    level_prev = Vector{T}(undef, m)
    n = length(onsets)
    na = length(alpha)
    @inbounds for i in 1:m
        u = onset_idx[i]
        onset_rate = (u >= 1 && u <= n) ? onsets[u] : zero(T)
        α = alpha[clamp(u - Int(grid_start) + 1, 1, na)]
        δ_cur = cur_report_idx[i] - u
        δ_prev = prev_report_idx[i] - u
        F_cur = onset_report_F(δ_cur, logit_h0, γ, u, grid_start, α)
        F_prev = onset_report_F(δ_prev, logit_h0, γ, u, grid_start, α)
        level_cur[i] = onset_rate * F_cur
        level_prev[i] = onset_rate * F_prev
        means[i] = level_cur[i] - level_prev[i]
    end
    return (; means, level_cur, level_prev)
end

"""
    onset_vintage_indices(report_idx, prev_report_idx)

Map each increment cell onto the vintage its two reads come from. The
vintages are the sorted distinct report days of the scored cells, so
vintage `s` is the `s`-th surviving snapshot [`load_onset_curve`](@ref)
kept. Returns `(; vintage_idx, prev_vintage_idx, n_vintages)`, the first
two length-matched to `report_idx`. The sentinel `prev_report_idx = 0`
(the virtual empty predecessor of the very first scored vintage, see
[`load_onset_curve`](@ref)) maps to `prev_vintage_idx = 0`, which
[`onset_scan_adjust`](@ref) reads as "no scan to dilate".

Derived here rather than carried on the history so every existing
`(; onset_days, report_days, prev_report_days, increments)` history keeps
working unchanged. It is integer work on the observation grid with no
sampled quantity in it, so it costs one pass per model evaluation and
contributes nothing to the gradient.
"""
function onset_vintage_indices(report_idx::AbstractVector{<:Integer},
        prev_report_idx::AbstractVector{<:Integer})
    days = sort(unique(report_idx))
    m = length(report_idx)
    vintage_idx = Vector{Int}(undef, m)
    prev_vintage_idx = Vector{Int}(undef, m)
    @inbounds for i in 1:m
        vintage_idx[i] = searchsortedfirst(days, report_idx[i])
        p = prev_report_idx[i]
        j = p > 0 ? searchsortedfirst(days, p) : 0
        prev_vintage_idx[i] = (j >= 1 && j <= length(days) && days[j] == p) ?
                              j : 0
    end
    return (; vintage_idx, prev_vintage_idx, n_vintages = length(days))
end

"""
    onset_scan_adjust(level_cur, level_prev, scan_level, vintage_idx,
        prev_vintage_idx)

Apply each vintage's own scan level to the modelled cumulative levels the
increment cells difference. Cell `i` reads its current level off scan
`vintage_idx[i]` and its previous level off scan `prev_vintage_idx[i]`, so

```math
\\ell_{\\text{cur},i} \\mapsto \\ell_{\\text{cur},i} \\, c_{s_i}, \\qquad
\\ell_{\\text{prev},i} \\mapsto \\ell_{\\text{prev},i} \\, c_{p_i},
\\qquad \\text{mean}_i = \\ell_{\\text{cur},i} c_{s_i} -
    \\ell_{\\text{prev},i} c_{p_i},
```

with `scan_level[s] = 1 + σ_scan · z_scan[s]` the level a whole published
figure was read at ([`onset_reporting_model`](@ref)). One multiplier per
scan, not per cell: the ~28 bars digitised off one figure share a single
level error, so a scan that reads high moves that snapshot's whole row of
cells together. An index outside `1:length(scan_level)` (the sentinel `0`
for the virtual empty predecessor) contributes a multiplier of exactly
one, so the first vintage's level cells carry only their own scan's
error. Returns the same `(; means, level_cur, level_prev)` shape
[`onset_report_moments`](@ref) does, ready for
[`onset_report_scales`](@ref). Pure, top-level, single indexed loop.
"""
function onset_scan_adjust(level_cur::AbstractVector,
        level_prev::AbstractVector, scan_level::AbstractVector,
        vintage_idx::AbstractVector{<:Integer},
        prev_vintage_idx::AbstractVector{<:Integer})
    m = length(level_cur)
    T = promote_type(eltype(level_cur), eltype(level_prev),
        eltype(scan_level))
    means = Vector{T}(undef, m)
    lc = Vector{T}(undef, m)
    lp = Vector{T}(undef, m)
    ns = length(scan_level)
    @inbounds for i in 1:m
        s = vintage_idx[i]
        p = prev_vintage_idx[i]
        cs = (s >= 1 && s <= ns) ? scan_level[s] : one(T)
        cp = (p >= 1 && p <= ns) ? scan_level[p] : one(T)
        lc[i] = level_cur[i] * cs
        lp[i] = level_prev[i] * cp
        means[i] = lc[i] - lp[i]
    end
    return (; means, level_cur = lc, level_prev = lp)
end

"""
    onset_report_scales(means, level_cur, level_prev, prev_report_idx;
        pixel_sd = 2.1, scan_sd = 0.0)

Per-cell observation scale for the reporting-triangle increment likelihood,
the square root of a variance built from three sources.

  - Counting variation of the cases the cell actually reports. The cell is a
    count of newly reported cases, so it carries its own sampling variation
    of about its mean, `means[i]`, on top of any reading error. This term is
    what makes the scale correct for the very first snapshot's cells, which
    are differenced against an empty predecessor and so score a level rather
    than a correction (see [`load_onset_curve`](@ref)): a level of 40 cases
    has counting variation of about `sqrt(40) ≈ 6`, far larger than the
    reading error below, and scoring it on reading error alone would let 28
    level cells dominate the joint likelihood.
  - Pixel-reading noise, roughly constant per bar read (`pixel_sd`, ≈2.1
    cases). An increment differences two independent reads, so its variance
    doubles. The first snapshot's level cells read only one bar, so theirs
    does not.
  - An optional multiplicative level error `scan_sd` on each read's own
    cumulative level, off by default. A bar's height is read in pixels and
    converted with the figure's own axis scale, so the multiplicative part
    of the digitisation error belongs to the whole scan and not to the bar:
    the likelihood carries it on the modelled level instead
    ([`onset_scan_adjust`](@ref)) and leaves this at zero. Callers that
    score a quantity the scan level has not already been applied to (a
    single digitised bar, or a whole projected snapshot total) pass the
    fitted per-scan level SD here.

```math
\\sigma_i = \\sqrt{\\max(\\mu_i, 0) + \\text{pixel\\_sd}^2 \\cdot r_i +
    \\text{scan\\_sd}^2 \\cdot
    (\\ell_{\\text{cur},i}^2 + \\ell_{\\text{prev},i}^2)},
\\qquad r_i = \\begin{cases} 1 & \\text{prev\\_report\\_idx}_i = 0
    \\ \\text{(virtual first pair)} \\\\ 2 & \\text{otherwise} \\end{cases},
```

with `μ_i = means[i]` the modelled increment. The counting term cancels for
a genuine correction between two snapshots only to the extent that the two
reads share the same realised cases: the newly reported cases in between are
a fresh count, and `μ_i` is exactly their expected number, so the same
formula covers both cell kinds without a branch.

Every magnitude entering the scale is modelled (`means`, `level_cur`,
`level_prev` from [`onset_report_moments`](@ref)), never the raw observed
count: feeding the likelihood's own noisy observation back into its variance
would bias towards overconfidence on cells that happen to undershoot. The
caller ([`onset_reporting_model`](@ref)) applies a sampled multiplicative
slack on top, so this fixed, measurement-derived formula is correctable by
the data rather than treated as exact, and the slack also stands in for
count overdispersion beyond the Poisson-like term above. Pure, top-level,
single indexed loop.

One approximation is left: the counting term is Poisson-like, with no
separate overdispersion parameter. `σ_mult` is the lever that can absorb
that shortfall, so it is the diagnostic for it, but read it in the right
direction: its prior (`onset_reporting_model`'s `slack_prior`) is bounded
below at 1 and unbounded above, so there is no upper edge for a posterior
to press against. The signature that says the term needs its own parameter
is `σ_mult` mass well above 1. A posterior sitting on the lower bound says
the opposite: the fit would like a tighter likelihood than the measurement
floor allows. A missing overdispersion parameter grows the shortfall with
the cell mean (a uniform multiplier inflates every cell but cannot change
the mean-variance shape), tested against the empirical/modelled residual
ratio across bins of `means` versus the `sqrt(ν/(ν-2))` a Student-t
likelihood implies. Two fits have now run that test and both put the
quadratic term at zero, which is why there is no negative-binomial term
here.

The shared part of the scan error is a different failure and needs a
different test. It leaves individual cells alone and shows up only once
the cells of one snapshot are summed, so it is tested on per-snapshot
rather than per-cell coverage: summing independent errors can reproduce
the right total spread while putting it in the wrong place, and central
coverage collapses while the 90% interval still looks nominal. That is
what the per-cell-only scale did, and why the shared component now sits
on the modelled level instead. Both tests are in the report's
symptom-onset reporting-delay section.
"""
function onset_report_scales(means::AbstractVector,
        level_cur::AbstractVector,
        level_prev::AbstractVector,
        prev_report_idx::AbstractVector{<:Integer};
        pixel_sd::Real = 2.1, scan_sd::Real = 0.0)
    m = length(level_cur)
    T = promote_type(eltype(means), eltype(level_cur), eltype(level_prev),
        typeof(float(pixel_sd)), typeof(float(scan_sd)))
    out = Vector{T}(undef, m)
    @inbounds for i in 1:m
        r = prev_report_idx[i] > 0 ? 2 : 1
        out[i] = onset_report_scale(means[i], level_cur[i], level_prev[i], r;
            pixel_sd, scan_sd)
    end
    return out
end

"""
    onset_report_scale(μ, level_cur, level_prev, reads;
        pixel_sd = 2.1, scan_sd = 0.0)

Scalar form of [`onset_report_scales`](@ref)'s per-cell formula, for one
increment mean `μ` between two modelled cumulative levels `level_cur`
and `level_prev` read off `reads` bars (`1` for a level differenced
against an empty predecessor, `2` for a genuine correction). The vector
method calls this, so the two cannot drift apart. The forecast
([`forecast_onsets`](@ref)) calls it directly to give a projected
reporting increment the same three-term observation scale the likelihood
gives a scored cell. See [`onset_report_scales`](@ref) for what each term
means and which of them a fit can correct.
"""
function onset_report_scale(μ::Real, level_cur::Real, level_prev::Real,
        reads::Integer; pixel_sd::Real = 2.1, scan_sd::Real = 0.0)
    T = promote_type(typeof(float(μ)), typeof(float(level_cur)),
        typeof(float(level_prev)), typeof(float(pixel_sd)),
        typeof(float(scan_sd)))
    return sqrt(max(μ, zero(T)) + pixel_sd^2 * reads +
                scan_sd^2 * (level_cur^2 + level_prev^2))
end

"""
    onset_report_cdf_extrapolated(δ, logit_h0, γ, u, grid_start)

Like [`onset_report_cdf`](@ref), but the calendar-time index into `γ` is
clamped to `[1, length(γ)]` rather than assumed in-range, so an onset date
`u` outside `[grid_start, grid_start + length(γ) - 1]` (a calendar day the
fitted `γ` walk has no support for) still returns a well-defined value: the
calendar effect at that unseen report day is held flat at the walk's
nearest known edge (`γ[1]` if the whole delay window falls before
`grid_start`, `γ[end]` if after). Agrees exactly with
[`onset_report_cdf`](@ref) whenever every index touched is already
in-range, so it is a strict extension, not a different formula for the
in-range case.

Used only by [`onset_report_expected_total`](@ref) to extend the cut-off
total to onset dates the digitised triangle itself never covers (days
`1:grid_start-1`, see that function's docstring): the increment likelihood
in [`onset_reporting_model`](@ref) never needs it, because every scored
onset date is inside `[grid_start, grid_end]` by construction
(`grid_start = minimum(onset_days)`). Pure, top-level, allocation-free.
"""
function onset_report_cdf_extrapolated(δ::Integer, logit_h0::AbstractVector,
        γ::AbstractVector, u::Integer, grid_start::Integer)
    T = promote_type(eltype(logit_h0), eltype(γ))
    δ < 0 && return zero(T)
    D = length(logit_h0)
    δc = min(Int(δ), D - 1)
    ng = length(γ)
    surv = one(T)
    @inbounds for j in 0:δc
        gi = clamp(u + j - grid_start + 1, 1, ng)
        h = logistic(logit_h0[j + 1] + γ[gi])
        surv *= (one(T) - h)
    end
    return one(T) - surv
end

"""
    onset_report_expected_total(onsets, logit_h0, γ, grid_start, alpha, as_of)

Expected reported symptom-onset total as of grid day `as_of`,
`Σ_u onsets[u] · F(u, as_of - u)` for `u` in `1:as_of` (clamped to
`1:length(onsets)`), the onset-stream analogue of every other stream's
`expected_*_T` cut-off deterministic. Callers pass the model cut-off `n` for
`as_of`, so the total is directly comparable to
`expected_confirmed_T`/`expected_deaths_T`, which sum their full `1:n`
series. Passing the triangle's own last report day instead would give a
total anchored a few days earlier than every sibling.

`γ` and `alpha` both span only the digitised triangle's own grid (see
[`onset_reporting_model`](@ref)), which starts after grid day 1 and ends at
or before `as_of`, so both the oldest and the most recent terms need a
calendar effect and an ascertainment level the fit has no estimate for.
[`onset_report_F`](@ref) already holds both flat at their nearest fitted
edge (see [`onset_report_G`](@ref)), so no separate extrapolated form is
needed here.

Safe for any `as_of` and any `γ`/`alpha` length, including the degenerate
`length(γ) < D` case, because both indices are clamped rather than assumed
in range. Pure, top-level, single indexed loop.
"""
function onset_report_expected_total(onsets::AbstractVector,
        logit_h0::AbstractVector, γ::AbstractVector,
        grid_start::Integer, alpha::AbstractVector, as_of::Integer)
    T = promote_type(eltype(onsets), eltype(logit_h0), eltype(γ),
        eltype(alpha))
    total = zero(T)
    n = length(onsets)
    na = length(alpha)
    ge = min(Int(as_of), n)
    @inbounds for u in 1:ge
        δ = as_of - u
        α = alpha[clamp(u - Int(grid_start) + 1, 1, na)]
        total += onsets[u] * onset_report_F(δ, logit_h0, γ, u, grid_start, α)
    end
    return total
end

"""
    onset_report_ascertainment(anchor_series, β, ω)

Ascertainment level `α(u) = logistic(logit(anchor(u)) + β + ω(u))` for
every onset date `u` the ascertainment walk spans, `anchor_series[k]` being
the delay-weighted anchor at onset date `grid_start + k - 1`
([`onset_report_anchor_series`](@ref)), `β` a sampled logit-scale offset and
`ω` the calendar-time random walk over the same onset dates
([`onset_ascertainment_model`](@ref)). `α` is a level rather than an
asymptote needing `D` days of walk to become observable, so it is returned
over the walk's full span with no restriction.

The anchor is clamped into `(0, 1)` before the logit, the same guard
`composition_positivity` applies to its own probability. An anchor
of exactly `0` or `1` is not reachable through the confirmed pipeline,
whose positivity is already clamped, but it is reachable if every hazard
underflows so that [`onset_report_anchor`](@ref)'s weights sum to zero.
`logit(0)` is `-Inf`, and while the forward value stays finite the
gradient is `NaN`, which would poison the whole log-density rather than
this stream's part of it. Pure, top-level, elementwise broadcast.
"""
function onset_report_ascertainment(anchor_series::AbstractVector,
        β::Real, ω::AbstractVector)
    T = promote_type(eltype(anchor_series), typeof(β), eltype(ω))
    lo = convert(T, 1e-8)
    hi = one(T) - lo
    safe = clamp.(ifelse.(isfinite.(anchor_series), anchor_series, lo),
        lo, hi)
    return logistic.(logit.(safe) .+ β .+ ω)
end

"""
Discrete symptom-onset reporting-delay hazard, nonparametric over the delay
and drifting over calendar time. Two non-centred random effects, matching
the repo's established idioms exactly:

  - a baseline logit hazard over the delay dimension `d = 0 … D-1`
    ([`confirmed_positivity_model`](@ref)'s per-vintage positivity random
    effect, reindexed to delay instead of vintage):
    ```math
    \\eta_0 \\sim \\text{baseline\\_prior}, \\quad
    \\sigma_{h0} \\sim \\text{pooling\\_prior}, \\quad
    z_{h0,d} \\sim \\mathcal N(0,1), \\quad
    \\text{logit\\_h0}(d) = \\eta_0 + \\sigma_{h0} z_{h0,d};
    ```
  - a calendar-time random walk on report date, weekly knots linearly
    interpolated to the daily grid ([`rt_walk_model`](@ref)'s non-centred
    cumulative-sum walk, same construction):
    ```math
    \\sigma_\\gamma \\sim \\text{walk\\_sigma\\_prior}, \\quad
    z_{\\gamma,k} \\sim \\mathcal N(0,1), \\quad
    \\gamma_{\\text{knot},1} = 0, \\
    \\gamma_{\\text{knot},k+1} = \\gamma_{\\text{knot},k} +
        \\sigma_\\gamma z_{\\gamma,k}.
    ```

The walk is indexed on the grid `[grid_start, grid_end]`, the report-date
axis, not the onset/infection-date axis that [`rt_walk_model`](@ref)
already carries a smooth calendar-time multiplicative effect on: indexing
this walk on onset date would make it unidentified against the
reproduction-number walk (both would be free to explain the same
onset-date-indexed rise and fall). Report date lags onset date by the
delay itself (mean ≈6 d) plus the ≈6 d incubation period the onset series
is already convolved from, so the two walks act on genuinely different,
if overlapping, calendar windows. See `onset_reporting_model`'s docstring
for what the vintage structure does and does not separate here.

The default `baseline_prior = Normal(logit(0.13), 0.7)` targets a median
onset-to-report delay of roughly 5 days under a constant-hazard
approximation (`(1 - 0.13)^5 ≈ 0.5`), close to the ≈6-day median a
line-list re-analysis found. The wide SD (0.7 logit-scale) leaves the
per-delay random effect free to depart substantially, so the fitted delay
should be read with the same honesty that analysis needed (its 7-day
reporting fraction carried a 43-68% 95% interval, wide because the
digitisation noise is comparable in size to the between-vintage increments
the estimate rests on).

`walk_sigma_prior` is a half-normal with SD 0.3, sized so a drift the
scored triangle can show is reachable rather than extreme. A 14% shift in
a snapshot's printed level needs a calendar shift of roughly `0.32` at
delay 8 and `0.58` at delay 12 on the logit hazard, holding the baseline
at its prior median. The walk is a cumulative sum of `σ_γ .* z_γ`, so
after the six weekly knots the scored window spans `γ` has SD
`E[σ_γ]·√6` = `0.586`, putting that shift at 0.5-1σ. It stays a
proper half-normal concentrated at zero, so a flat reporting profile is
the default the data has to argue away from.

The reason not to widen it further is the reproduction-number walk rather
than noise absorption. The two act on overlapping calendar windows and
both are least constrained over the final fortnight, so a calendar walk
with too much freedom can start explaining recent onset-date structure
that belongs to `R_t`. Any change here should report `R_t` over the final
fortnight and `C_T` either side, and treat a material move as a reason to
stop.

Returns `(; logit_h0, γ, grid_start, η0, σ_h0, σ_γ)`, with `γ` length
`max(grid_end - grid_start + 1, 1)`, indexed from `grid_start`.
"""
@model function onset_report_hazard_model(grid_start::Integer,
        grid_end::Integer;
        D::Integer = ONSET_REPORT_MAX_DELAY,
        baseline_prior = Normal(logit(0.13), 0.7),
        pooling_prior = truncated(Normal(0.0, 1.0); lower = 0),
        walk_sigma_prior = truncated(Normal(0.0, 0.3); lower = 0),
        week::Integer = 7)
    ## Non-centred logit random effect over the delay dimension, one value
    ## per delay day `d = 0 … D-1`.
    η0 ~ baseline_prior
    σ_h0 ~ pooling_prior
    z_h0 ~ product_distribution(fill(Normal(0, 1), D))
    logit_h0 = η0 .+ σ_h0 .* z_h0

    ## Non-centred cumulative random walk on weekly knots over the report-
    ## date grid `[grid_start, grid_end]`, linearly interpolated to the
    ## daily grid (see `knot_days`/`interpolate_knots`, `renewal.jl`). The
    ## local day count `nt` is floored at 1 so an empty/degenerate grid
    ## (the opt-in no-op path) still returns a well-formed length-1 `γ`.
    nt = max(Int(grid_end) - Int(grid_start) + 1, 1)
    days = knot_days(nt; week, start = 1)
    nb = length(days)
    σ_γ ~ walk_sigma_prior
    z_γ ~ product_distribution(fill(Normal(0, 1), max(nb - 1, 1)))
    steps = σ_γ .* z_γ[1:max(nb - 1, 0)]
    γ_knots = vcat(zero(σ_γ), cumsum(steps))
    γ = interpolate_knots(γ_knots, days, nt)

    return (; logit_h0, γ, grid_start = Int(grid_start), η0, σ_h0, σ_γ)
end

"""
Ascertainment level over the onset-date grid `[grid_start, grid_end]`: a
logit-scale offset and slow random walk on top of a delay-weighted anchor
series `anchor_series` ([`onset_report_anchor_series`](@ref)),

```math
\\beta \\sim \\text{beta\\_prior}, \\quad
\\sigma_a \\sim \\text{pooling\\_prior}, \\quad
z_{a,k} \\sim \\mathcal N(0,1), \\quad
\\omega_{\\text{knot},1} = 0, \\
\\omega_{\\text{knot},k+1} = \\omega_{\\text{knot},k} + \\sigma_a z_{a,k},
```

with `alpha` given by [`onset_report_ascertainment`](@ref). `beta_prior =
Normal(0, 0.75)` is a zero-centred logit-scale departure from the anchor,
so the triangle's ascertainment defaults to the anchor series' own level
and can depart by roughly a factor of two either way. `pooling_prior =
truncated(Normal(0, 0.1); lower = 0)` is deliberately tight: `omega` shares
the onset-date axis with [`rt_walk_model`](@ref), so a flat ascertainment
level is the default the data has to argue away from. Weekly knots,
linearly interpolated ([`knot_days`](@ref)/[`interpolate_knots`](@ref)),
first knot pinned at zero, mirroring [`onset_report_hazard_model`](@ref)'s
calendar walk exactly.

Returns `(; alpha, β, σ_a, z_a, ω)`, `alpha` and `ω` length `nt =
max(grid_end - grid_start + 1, 1)`.
"""
@model function onset_ascertainment_model(anchor_series::AbstractVector,
        grid_start::Integer, grid_end::Integer;
        beta_prior = Normal(0.0, 0.75),
        pooling_prior = truncated(Normal(0.0, 0.1); lower = 0),
        week::Integer = 7)
    nt = max(Int(grid_end) - Int(grid_start) + 1, 1)
    days = knot_days(nt; week, start = 1)
    nb = length(days)
    β ~ beta_prior
    σ_a ~ pooling_prior
    z_a ~ product_distribution(fill(Normal(0, 1), max(nb - 1, 1)))
    steps = σ_a .* z_a[1:max(nb - 1, 0)]
    ω_knots = vcat(zero(σ_a), cumsum(steps))
    ω = interpolate_knots(ω_knots, days, nt)
    alpha = onset_report_ascertainment(anchor_series, β, ω)
    return (; alpha, β, σ_a, z_a, ω)
end

"""
Symptom-onset reporting-triangle observation model: fits the between-vintage
increments of the digitised reporting triangle
([`load_onset_curve`](@ref)) against the shared latent daily onset series
`onsets`, through a nonparametric delay hazard
([`onset_report_hazard_model`](@ref)) and an explicit ascertainment level
([`onset_ascertainment_model`](@ref)). See [`onset_report_F`](@ref) for how
the two combine and [`onset_report_cdf`](@ref) for how right truncation
enters.

Close in spirit to how the R package `baselinenowcast` treats a reporting
triangle: a triangle of between-vintage increments, not a single total
column. It differs from how [`reported_cases_model`](@ref) and every other
stream here is scored (and from EpiNow2, which nowcasts a single evolving
total): the same right-truncation mechanism (`F(u, δ) = 0` for `δ < 0`) is
applied to corrections between snapshots rather than to a level, so a case
already scored in an earlier snapshot is never counted again.

**What the triangle separates.** Four time-varying multiplicative objects
now act on or near the same latent series: the reproduction-number walk
([`rt_walk_model`](@ref)) on the infection/onset axis, this stream's
calendar walk `γ` on the report axis, the delay shape's baseline
`logit_h0`, and the ascertainment walk `ω` on the onset axis. A change in
the onset series moves a column of scored cells. A change in `γ` moves a
row. A change in `logit_h0` moves a diagonal band. Column, row and band are
distinguishable once there is more than one snapshot, which is the
structural reason this stream is worth fitting, and the reason `γ` is
indexed on the report day rather than the onset day (indexing it on onset
date would make it unidentified against `rt_walk_model`).

Ascertainment is anchored on the confirmed pipeline's own daily
ascertainment (`p_drc · τ_test · p_pos_grid`, delay-weighted onto the onset
axis by [`onset_report_anchor`](@ref)) rather than left free: `bvd_joint`
passes that series in as `anchor`, `onsets_only_model` falls back to a
constant `0.15` anchor (no confirmed pipeline to borrow from). `β` and `ω`
let the fitted level depart from that anchor, with `ω`'s tight prior making
a flat departure the default the data has to argue away from (see
[`onset_ascertainment_model`](@ref)).

Three things stay genuinely weak. First, `logit_h0` and `alpha` are pinned
by levels, not by corrections: corrections constrain only differences of
`F`, so what breaks the tie is the first snapshot's cells (differenced
against an empty predecessor, see [`load_onset_curve`](@ref)) together with
the onset series being pinned by the other streams. A single-stream
[`onsets_only_model`](@ref) fit has neither, so its ascertainment and `C_T`
stay close to prior-driven: worth running as a comparison, not an estimate.
Second, the hazard at the shortest delays is barely observed (published
figures stop short of the report date), so those hazards rest on partial
pooling to `η0`, not data. Third, a falling ascertainment and a slowing
delay shape both suppress recent bars within one snapshot. Right truncation
self-corrects for the delay explanation but not a genuine ascertainment
fall, so the vintage structure is what separates them, and only as well as
the (under-ten) snapshots covering the affected dates allow.

What has not been checked is whether `γ`, `ω` and `rt_walk_model` pull on
each other in practice: all three act on overlapping calendar windows and
are least constrained over the final fortnight. Any change here should
report `R_t` over the final fortnight and `C_T` either side, treating a
material move as a reason to stop rather than a result.

The alive/dead split the raw triangle carries is not modelled separately:
only `confirmed_total` is fitted, since the confirmed-death stream already
carries that split from other data.

`onset_curve_history` is the [`load_onset_curve`](@ref) return shape
`(; onset_days, report_days, prev_report_days, increments)`. The default
empty history makes every loop here a no-op, the degrade-gracefully path
for a missing input file. `increments` may be `missing` to sample instead
of condition (the predictive-generator path).

The observation scale ([`onset_report_scales`](@ref)) is built from
counting variation, the measured per-bar pixel noise (`pixel_sd` ≈2.1
cases/bar), and a sampled multiplicative slack
`σ_mult ~ slack_prior` bounded below at 1 (each scale term is a lower bound
on the truth, so a fitted scale below them would let a couple of hundred
cells outvote every other stream). A short onsets-only run pulling the
slack to the bound is expected there. In the joint fit a `σ_mult` posterior
well above 1 says the scale is missing a term.

**The scan error belongs to the figure, not to the bar.** A bar's height
is read in pixels and converted with the axis scale that scan calibrated,
so the digitisation error splits by construction: an absolute per-bar
term (the ≈2.1-case pixel noise) and a multiplicative term that is one
number for the whole figure. `data/README.md` measures the second
directly, as each vintage's digitised total against the total the figure
prints: -5.0% to +1.6% over the audited vintages, one value per scan.
Two mechanisms feed it. The axis calibration is read per figure, and the
digitiser's colour masks are fixed thresholds against a render blur that
varies with the embedded image's size, so a blurrier scan loses more of
every bar's edge. Neither explains SitRep 088, whose digitised total falls
184 below SitRep 087 on the same render size, the same edge softness and
the same plotted window, so the level moves for reasons the observable
covariates do not cover and is sampled free rather than regressed on them.

Scoring that multiplicative term as independent per-cell noise, which is
what this stream first did, reproduces the right total spread per
snapshot and puts it in the wrong place. A snapshot's cells can then only
miss in uncorrelated directions, so the net correction the report's panel
plots is far too tightly predicted: 1 of 11 snapshots inside a nominal
50% interval, against an aggregate variance ratio of 1.07 and 90%
coverage at nominal. It also made a single cell's scale about twice as
wide as its residuals wanted.

The per-scan level is therefore a sampled level `1 + σ_scan · z_scan[s]`
on the modelled cumulative levels each cell differences
([`onset_scan_adjust`](@ref)), and the per-cell scale keeps counting and
pixel noise alone. `σ_scan ~ scan_sd_prior` is estimated rather than
fixed at the audited spread, which rests on a handful of vintages, but
its prior is a half-normal centred to put that spread (an SD of about
2.5%) in its bulk. The upper bound is the audit's own reach: no vintage
has ever read more than 5% away from its printed total, so a level error
past 8% is excluded by a data check rather than by taste.

The shared level also gives the fit somewhere to put a snapshot that
reprints nothing new. A vintage whose figure reads at the same level as
its predecessor scores increments of about zero across its whole row, and
with no per-scan level the only way to fit that row is to drive the
reporting hazard towards zero at the delays it covers, which then
propagates into every other snapshot's delay shape. With one, the row is
explained as two scans read at the same level, and the hazard is left to
the snapshots that do carry new reports.

The likelihood is Student-t with fixed degrees of freedom `ν` (default 4):
with only a few hundred cells, `ν` is weakly identified, and the repo's own
`lab_delay_model` docstring makes the same argument against sampling it.
The heavy tail lets the (frequently negative) measured increments score as
large-but-plausible residuals rather than breaking a count likelihood.

Returns `(; increments, modelled, unscanned, scan_level, logit_h0, γ,
grid_start, grid_end, alpha, σ_mult, σ_scan, η0, σ_h0, σ_γ, β, σ_a)` with
`modelled` the per-cell increment means the likelihood scores (each
vintage's scan level applied), `unscanned` the same means before it (the
epidemiological signal alone), `scan_level` the per-vintage multipliers,
`grid_end` the report-date grid day the calendar walk was built up to
(`max(report_days)`, or `grid_start` when the history is empty), and the
hyperparameters re-exposed at this level for the pairs-plot summary.
"""
@model function onset_reporting_model(
        onset_curve_history, onsets::AbstractVector;
        hazard = onset_report_hazard_model,
        ascertainment = onset_ascertainment_model,
        anchor::AbstractVector = [0.15],
        D::Integer = ONSET_REPORT_MAX_DELAY,
        pixel_sd::Real = 2.1,
        scan_sd_prior = truncated(Normal(0.0, 0.03);
            lower = 0.0, upper = 0.08),
        slack_prior = truncated(Normal(1.0, 0.5); lower = 1.0),
        ν::Real = 4.0)
    onset_days = onset_curve_history.onset_days
    report_days = onset_curve_history.report_days
    prev_report_days = onset_curve_history.prev_report_days
    m = length(onset_days)
    ## Report-date grid the calendar walk spans: the union of every onset
    ## and report day a scored cell can touch. Falls back to a degenerate
    ## length-1 grid `[1, 1]` when the history is empty (the no-op path),
    ## which `onset_report_hazard_model` handles via its own `nt` floor.
    grid_start = m > 0 ? minimum(onset_days) : 1
    grid_end = m > 0 ? max(maximum(report_days), grid_start) : 1
    ## Unprefixed (`false`): the hazard model has no `:=` deterministics to
    ## collide with, and hoisting its sampled variables (`η0`, `σ_h0`,
    ## `σ_γ`, …) straight into this frame means they surface as a single
    ## `onset_report_state.η0` etc at the composer level, rather than the
    ## double-nested `onset_report_state.hazard_state.η0` a prefixed
    ## attachment would give: the flat form the pairs-plot summary uses.
    hazard_state ~ to_submodel(hazard(grid_start, grid_end; D), false)

    ## Delay-weighted anchor series over the onset-date grid, built from the
    ## fitted hazard and the caller-supplied calendar-indexed daily
    ## ascertainment `anchor` (the confirmed pipeline's own series, or the
    ## length-1 constant default). Attached unprefixed for the same reason
    ## `hazard_state` is: the ascertainment hyperparameters (`β`, `σ_a`, …)
    ## surface as flat `onset_report_state.β` etc.
    anchor_series = onset_report_anchor_series(hazard_state.logit_h0,
        hazard_state.γ, grid_start, grid_end, anchor)
    asc_state ~ to_submodel(
        ascertainment(anchor_series, grid_start, grid_end), false)
    alpha = asc_state.alpha
    σ_mult ~ slack_prior

    ## Per-scan level error: one non-centred multiplier per surviving
    ## vintage, the level that scan's whole figure was read at. It sits on
    ## the modelled level rather than in the per-cell scale because a
    ## figure's axis calibration is one number for the whole figure, so the
    ## scale below carries counting and pixel noise alone.
    vintages = onset_vintage_indices(report_days, prev_report_days)
    σ_scan ~ scan_sd_prior
    z_scan ~ product_distribution(fill(Normal(0, 1),
        max(vintages.n_vintages, 1)))
    scan_level = one(σ_scan) .+ σ_scan .* z_scan

    moments = onset_report_moments(onsets, hazard_state.logit_h0,
        hazard_state.γ, hazard_state.grid_start, alpha, onset_days,
        report_days, prev_report_days)
    scanned = onset_scan_adjust(moments.level_cur, moments.level_prev,
        scan_level, vintages.vintage_idx, vintages.prev_vintage_idx)
    scales = onset_report_scales(scanned.means, scanned.level_cur,
        scanned.level_prev, prev_report_days; pixel_sd)

    ## Scored in a dedicated submodel so `increments` is a model argument
    ## on the left of `~`: DynamicPPL decides observe-versus-assume from
    ## whether the tilde's symbol is in the enclosing model's argument
    ## names, so pulling the observations out of `onset_curve_history` into
    ## a local here would make every cell a latent variable and drop the
    ## likelihood silently. Attached unprefixed so the cells keep the flat
    ## `increments[i]` names the predictive path indexes.
    increments_state ~ to_submodel(
        onset_increments_model(scanned.means, σ_mult .* scales,
            onset_curve_history.increments, ν), false)
    increments = increments_state.increments

    return (; increments, modelled = scanned.means,
        unscanned = moments.means, scan_level,
        logit_h0 = hazard_state.logit_h0, γ = hazard_state.γ,
        grid_start = hazard_state.grid_start, grid_end, alpha, σ_mult,
        σ_scan, η0 = hazard_state.η0, σ_h0 = hazard_state.σ_h0,
        σ_γ = hazard_state.σ_γ, β = asc_state.β, σ_a = asc_state.σ_a)
end
