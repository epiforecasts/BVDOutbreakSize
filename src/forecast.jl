# One-week-ahead posterior-predictive forecast. Continues the fitted
# exponential growth `h` days past the cut-off `T` and projects the four
# trusted quantities over the coming week: new infections, new true BVD
# deaths, new laboratory-confirmed cases (DRC) and new laboratory-confirmed
# deaths (DRC). The untrusted suspected/reported cases, suspected deaths and
# tests-analysed streams are dropped. Each draw produces a replicated
# integer count so the intervals include both parameter and observation
# uncertainty.

# Deaths convolution evaluated at horizon `Th = T + h`, sharing the
# package-wide `delay_convolution` integrator. `onset_fraction` maps the
# latent infection trajectory onto onsets (the incubation mgf), matching
# `deaths_model`.
function _forecast_deaths_mean(r, Th, α, θ, CFR;
        onset_fraction = 1.0, alg = DEATH_INTEGRAL_ALG)
    return onset_fraction * delay_convolution(CFR, r, Th, Gamma(α, θ); alg)
end

# Reported (suspected) cases at horizon `Th`: the truth-anchored BVD
# contribution plus the constant-rate non-BVD background. The BVD
# contribution reuses `delay_convolution` as a delay-convolved cumulative
# integrator at unit ascertainment to compute
# `∫₀^{Th} exp(r·s) · f_rep(Th-s) ds`, scaled by `onset_fraction` (the
# incubation mgf) since the latent trajectory is infections; the background
# contribution is the constant cumulative `μ_bg(Th) = λ_bg · Th`, non-BVD
# and unscaled.
function _forecast_cases_mean(r, Th, α_rep, θ_rep, p_drc, λ_bg;
        onset_fraction = 1.0, alg = DEATH_INTEGRAL_ALG)
    conv = delay_convolution(one(p_drc), r, Th, Gamma(α_rep, θ_rep); alg)
    return onset_fraction * p_drc * conv + max(λ_bg * Th, zero(λ_bg))
end

# Cumulative samples received at horizon `Th`: the forwarded fraction
# `τ_forward` of the suspect backlog `N_susp = μ_BVD + μ_bg` convolved with
# the receipt-delay kernel `f_receipt` (specimens reach the lab after a
# transport delay), matching the received stream of `confirmed_cases_model`.
function _forecast_received(r, Th, α_rep, θ_rep, p_drc, λ_bg, τ_forward,
        f_receipt; onset_fraction = 1.0)
    f_rep = Gamma(α_rep, θ_rep)
    N_susp = let r = r, f_rep = f_rep, p_drc = p_drc, onset_fraction = onset_fraction,
        λ_bg = λ_bg

        u -> u <= zero(u) ? zero(u) :
             onset_fraction * p_drc *
             delay_convolution(one(p_drc), r, u, f_rep) +
             max(λ_bg * u, zero(λ_bg))
    end
    N_recv = delay_convolution(N_susp, Th, f_receipt)
    return τ_forward * max(N_recv, zero(N_recv))
end

# Composition-linked per-test positivity at horizon `Th`:
# `p_pos = s·q + (1−spec)·(1−q)`, with the tested BVD share `q` the
# suspect-pool composition `φ = μ_BVD/(μ_BVD+μ_bg)` at the horizon. Far out
# the cumulative analysed volume is large, so the severity enrichment
# `δ0·exp(−c/decay)` has fully relaxed and the tested share equals `φ`.
function _forecast_positivity(r, Th, α_rep, θ_rep, p_drc, λ_bg, s_test,
        spec_test; onset_fraction = 1.0)
    f_rep = Gamma(α_rep, θ_rep)
    μ_bvd = Th <= zero(Th) ? zero(Th) :
            onset_fraction * p_drc *
            delay_convolution(one(p_drc), r, Th, f_rep)
    μ_bg = max(λ_bg * Th, zero(λ_bg))
    denom = μ_bvd + μ_bg
    q = denom > zero(denom) ? clamp(μ_bvd / denom, zero(denom), one(denom)) :
        zero(denom)
    return clamp(s_test * q + (one(q) - spec_test) * (one(q) - q),
        zero(q), one(q))
end

# Cumulative confirmed deaths at horizon `Th`: a fraction `τ_death` of the
# suspect-death backlog increment over the horizon, weighted by the
# death-specimen positivity `p_pos_death = s·q_death + (1−spec)(1−q_death)`,
# with `q_death = μ_BVD_death / N_death_susp` the BVD share of the
# suspect-death pool. The suspect-death backlog is the modelled BVD-death
# trajectory (`CFR·p_deaths·os` convolved with the onset-to-death delay)
# plus the constant-rate `λ_bg_death` background, matching
# `confirmed_deaths_model` and the joint's `nsusp_death_edges`. Returns the
# forwarded positive increment over `(T, Th]` to add to the observed
# cumulative confirmed deaths.
function _forecast_confirmed_deaths_increment(r, T, Th, α, θ, CFR, p_deaths,
        λ_bg_death, τ_death, s, spec; onset_fraction = 1.0,
        alg = DEATH_INTEGRAL_ALG)
    f_death = Gamma(α, θ)
    bvd_death(u) = u <= zero(u) ? zero(u) :
                   onset_fraction * p_deaths *
                   delay_convolution(CFR, r, u, f_death; alg)
    nsusp(u) = max(bvd_death(u), zero(u)) +
               max(λ_bg_death * u, zero(λ_bg_death))
    ΔN = max(nsusp(Th) - nsusp(T), zero(Th))
    ## BVD share of the suspect-death pool at the horizon sets positivity.
    denom = nsusp(Th)
    q_d = denom > zero(denom) ?
          clamp(max(bvd_death(Th), zero(Th)) / denom, zero(Th), one(Th)) :
          zero(Th)
    p_pos = clamp(s * q_d + (one(q_d) - spec) * (one(q_d) - q_d),
        zero(q_d), one(q_d))
    return max(τ_death * p_pos * ΔN, zero(ΔN))
end

# Capacity-limited analysed throughput over the horizon: the lab processes
# a fraction `1 − exp(−κ·h/backlog)` of the received backlog not yet
# analysed (`backlog = received − analysed_cutoff`), with daily capacity `κ`
# held at its cut-off value over the `h`-day horizon.
function _forecast_analysed_increment(received, analysed_cutoff, κ, horizon)
    backlog = max(received - analysed_cutoff, eps(typeof(received)))
    cap = max(κ * horizon, zero(received))
    return backlog * (one(received) - exp(-cap / backlog))
end

function _nb_rand(rng, k, μ)
    μs = max(μ, eps(typeof(μ)))
    p = clamp(k / (k + μs), eps(typeof(k)), one(k) - eps(typeof(k)))
    return rand(rng, NegativeBinomial(k, p))
end

"""
One-week-ahead (default `horizon = 7` days) posterior-predictive
forecast of the four trusted quantities. For each draw, continue
exponential growth to `T + horizon` and project the new counts over the
coming week, returning a `DataFrame` with one row per draw and columns:

- `:infections_new` — new latent BVD infections over the week, the
  continued-growth increment `C(T) · (exp(r·h) − 1)` of the cumulative
  infection trajectory `C(T) = 2^m` (chain symbol `:cumulative_infections`).
- `:bvd_deaths_new` — new true BVD deaths over the week, the latent
  all-BVD-deaths mean increment `μ_bvd(T+h) − μ_bvd(T)` with
  `μ_bvd(s) = os · CFR · ∫₀^s exp(r·u) f_death(s−u) du` (no `p_deaths`, no
  background, no observed-count subtraction; these are latent toll counts).
- `:confirmed_new` — new laboratory-confirmed cases over the week from the
  queue of [`confirmed_cases_model`](@ref): the received backlog
  (`τ_forward · N_susp` convolved with the receipt delay) is analysed at
  the cut-off capacity `κ` over the horizon,
  `ΔA = backlog·(1 − exp(−κ·h/backlog))`, and the new positives are `ΔA`
  times the composition-linked positivity. Present when the chain carries
  the confirmed-stream parameters (`:s_test`, `:spec_test`, `:τ_forward`,
  `:α_recv`, `:θ_recv`, `:capacity_cutoff`) and `obs_confirmed` and
  `obs_analysed` are supplied; otherwise absent.
- `:confirmed_deaths_new` — new laboratory-confirmed deaths over the week,
  the forwarded positive increment of the suspect-death backlog, when the
  chain additionally carries `:τ_death`, `:p_deaths`, `:λ_bg_death`
  (sharing the case-lab `:s_test` / `:spec_test`) and `obs_confirmed_deaths`
  is supplied; otherwise absent.

`obs_confirmed` and `obs_analysed` anchor the confirmed-case queue at the
cut-off (the lab backlog already received and analysed by `T`);
`obs_confirmed_deaths` anchors the confirmed-death queue. Exports remain
available behind `forecast_exports = true`, which adds `:exports_new`
(`p_uganda · q`, `q = daily_travellers / source_population`); the report
default is `false` (cross-border travel disrupted, so the forward travel
rate no longer holds). Assumes growth continues unchanged over the horizon
(no interventions, no saturation).

`report_onset_offset` is accepted for call-site compatibility with the
fit; the composition-linked horizon positivity does not use it.
The default `nothing` keeps the seeding-anchored clock (`t_report = 0`).
"""
function forecast_reported(chn;
        horizon::Real = 7,
        daily_travellers::Real,
        source_population::Real,
        obs_confirmed::Union{Real, Missing} = missing,
        obs_confirmed_deaths::Union{Real, Missing} = missing,
        obs_analysed::Union{Real, Missing} = missing,
        obs_exports::Union{Real, Missing} = missing,
        forecast_exports::Bool = false,
        seed::Integer = 20260520,
        report_onset_offset::Union{Nothing, Real} = nothing,
        alg = DEATH_INTEGRAL_ALG)
    r = _draws(chn, :r)
    T = _draws(chn, :T)
    CFR = _draws(chn, :CFR)
    α = _draws(chn, :α)
    θ = _draws(chn, :θ)
    ## Latent cumulative infections `C(T) = 2^m` (a `:=` on the model);
    ## continued growth over the horizon gives the new-infections increment.
    cumulative_infections = _draws(chn, :cumulative_infections)
    ## Uganda exports use either the McCabe rectangular detection window
    ## `w` or the explicit onset-to-detection delay convolution. The delay
    ## reuses the DRC onset-to-report delay `Gamma(α_rep, θ_rep)`, so there
    ## is no separate export-delay parameter; the chain is identified as
    ## window- vs delay-fitted by the presence of `w`.
    has_window = haskey_chain(chn, :w)
    w = has_window ? _draws(chn, :w) : nothing
    pr = _draws(chn, :p_drc)
    pu = _draws(chn, :p_uganda)
    k = _draws(chn, :k)
    α_rep = _draws(chn, :α_rep)
    θ_rep = _draws(chn, :θ_rep)
    ## Constant non-BVD background rate.
    λ_bg = _draws(chn, :λ_bg)
    ## Incubation draws map the latent infection trajectory onto onsets
    ## (the `onset_fraction = mgf(incubation, −r)` of the observation
    ## models). Absent on chains predating the infection layer, where
    ## `onset_fraction = 1` recovers the previous behaviour.
    has_incubation = all(haskey_chain(chn, n) for n in (:α_inc, :θ_inc))
    α_inc = has_incubation ? _draws(chn, :α_inc) : nothing
    θ_inc = has_incubation ? _draws(chn, :θ_inc) : nothing
    ## Severe-first confirmed-stream draws (PCR sensitivity / specificity and
    ## the q-curve shape) live on the joint chain only; their absence drops
    ## the confirmed-cases columns.
    ## The forwarding fraction is `τ_forward_out` on the production queue
    ## chain (where `τ_forward` is fixed, not sampled) and `τ_forward` on the
    ## prior/test chains; resolve whichever is present.
    fwd_key = _forward_key(chn)
    has_lab = fwd_key !== nothing &&
              all(haskey_chain(chn, n)
              for n in (:s_test, :spec_test, :α_recv, :θ_recv,
                  :capacity_cutoff)) &&
              obs_confirmed !== missing && obs_analysed !== missing
    s_test = has_lab ? _draws(chn, :s_test) : nothing
    spec_test = has_lab ? _draws(chn, :spec_test) : nothing
    τ_forward = has_lab ? _draws(chn, fwd_key) : nothing
    α_recv = has_lab ? _draws(chn, :α_recv) : nothing
    θ_recv = has_lab ? _draws(chn, :θ_recv) : nothing
    capacity = has_lab ? _draws(chn, :capacity_cutoff) : nothing

    ## Confirmed-death draws: the death-specimen forwarding fraction
    ## `τ_death`, the deaths drift `p_deaths` and background `λ_bg_death`,
    ## sharing the case-lab sensitivity `s_test` / specificity `spec_test`
    ## and the deaths CFR / onset-to-death delay (`CFR`, `α`, `θ`). Present
    ## on the joint chain only; their absence (or a missing
    ## `obs_confirmed_deaths`) drops the confirmed-death columns.
    has_lab_deaths = has_lab &&
                     all(haskey_chain(chn, n)
                     for n in (:τ_death, :p_deaths, :λ_bg_death)) &&
                     obs_confirmed_deaths !== missing
    τ_death = has_lab_deaths ? _draws(chn, :τ_death) : nothing
    p_deaths = has_lab_deaths ? _draws(chn, :p_deaths) : nothing
    λ_bg_death = has_lab_deaths ? _draws(chn, :λ_bg_death) : nothing

    rng = MersenneTwister(seed)
    n = length(r)
    q = daily_travellers / source_population
    infections_new = Vector{Int}(undef, n)
    bvd_deaths_new = Vector{Int}(undef, n)
    exports_new = forecast_exports ? Vector{Int}(undef, n) : nothing
    confirmed_new = has_lab ? Vector{Int}(undef, n) : nothing
    confirmed_deaths_new = has_lab_deaths ? Vector{Int}(undef, n) : nothing

    @inbounds for i in 1:n
        Th = T[i] + horizon
        os = has_incubation ?
             onset_rescale(Gamma(α_inc[i], θ_inc[i]), r[i]) : 1.0
        ## New latent infections over the week: the continued-growth
        ## increment of the cumulative infection trajectory C(T) = 2^m.
        infections_new[i] = max(
            round(Int,
                cumulative_infections[i] * (exp(r[i] * horizon) - 1)), 0)
        ## New true BVD deaths over the week: the latent all-BVD-deaths mean
        ## increment μ_bvd(T+h) − μ_bvd(T) with μ_bvd(s) = os · CFR ·
        ## ∫₀^s exp(r·u) f_death(s−u) du (no p_deaths, no background, no
        ## observed-count subtraction; these are latent toll counts).
        μ_bvd_Th = _forecast_deaths_mean(r[i], Th, α[i], θ[i], CFR[i];
            onset_fraction = os, alg)
        μ_bvd_T = _forecast_deaths_mean(r[i], T[i], α[i], θ[i], CFR[i];
            onset_fraction = os, alg)
        bvd_deaths_new[i] = _nb_rand(rng, k[i],
            max(μ_bvd_Th - μ_bvd_T, zero(μ_bvd_Th)))
        ## Uganda exports. With the McCabe window: p_uganda · q ·
        ## ∫_{T+h−w}^{T+h} C(s) ds (closed form). With the delay
        ## mechanism: the at-risk export window runs from infection, so
        ## the detection delay is incubation ⊕ onset-to-report (combined
        ## to one Gamma); the incubation period lives in the delay, so the
        ## person-time integral is not separately rescaled by `os`. Skipped
        ## when `forecast_exports = false` (e.g. cross-border travel
        ## disrupted, so the forward travel rate no longer holds).
        if forecast_exports
            if has_window
                lo = max(Th - w[i], zero(Th))
                μ_exports = pu[i] * q * (exp(r[i] * Th) - exp(r[i] * lo)) /
                            r[i]
            else
                f_det = has_incubation ?
                        combined_delay(Gamma(α_inc[i], θ_inc[i]),
                    Gamma(α_rep[i], θ_rep[i])) : Gamma(α_rep[i], θ_rep[i])
                μ_exports = expected_exports_delay(r[i], pu[i], q, Th, f_det;
                    alg = alg)
            end
            exports_new[i] = rand(rng,
                Poisson(max(μ_exports, eps(μ_exports))))
        end
        if has_lab
            ## Received backlog (receipt-delayed) and the capacity-limited
            ## analysed increment over the horizon; the new confirmed
            ## positives this week are the horizon positivity times that
            ## increment (the queue increment over `(T, T+h]`).
            f_receipt = Gamma(α_recv[i], θ_recv[i])
            received = _forecast_received(r[i], Th, α_rep[i], θ_rep[i],
                pr[i], λ_bg[i], τ_forward[i], f_receipt; onset_fraction = os)
            Δanalysed = _forecast_analysed_increment(received, obs_analysed,
                capacity[i], horizon)
            p_pos = _forecast_positivity(r[i], Th, α_rep[i], θ_rep[i],
                pr[i], λ_bg[i], s_test[i], spec_test[i]; onset_fraction = os)
            confirmed_new[i] = _nb_rand(rng, k[i], p_pos * Δanalysed)
        end
        if has_lab_deaths
            ## New confirmed deaths this week: the forwarded positive
            ## increment of the suspect-death backlog over `(T, T+h]`,
            ## sharing the case-lab sensitivity / specificity.
            μ_cdeath = _forecast_confirmed_deaths_increment(r[i], T[i], Th,
                α[i], θ[i], CFR[i], p_deaths[i], λ_bg_death[i], τ_death[i],
                s_test[i], spec_test[i]; onset_fraction = os, alg)
            confirmed_deaths_new[i] = _nb_rand(rng, k[i], μ_cdeath)
        end
    end

    df = DataFrame(
        infections_new = infections_new,
        bvd_deaths_new = bvd_deaths_new
    )
    if forecast_exports
        df.exports_new = exports_new
    end
    if has_lab
        df.confirmed_new = confirmed_new
    end
    if has_lab_deaths
        df.confirmed_deaths_new = confirmed_deaths_new
    end
    return df
end

## Best-effort presence check for a chain key across the FlexiChains /
## MCMCChains containers in use. Avoids loading the FlexiChains type
## just to dispatch.
function haskey_chain(chn, name::Symbol)
    try
        chn[name]
        return true
    catch
        return false
    end
end

## Resolve the suspect-case forwarding fraction symbol. The production queue
## chain pins `τ_forward` (unsampled) and exposes the derived
## `τ_forward_out`; prior/test chains sample `τ_forward` directly. Returns
## the present symbol, preferring the sampled `τ_forward`, or `nothing`.
function _forward_key(chn)
    haskey_chain(chn, :τ_forward) && return :τ_forward
    haskey_chain(chn, :τ_forward_out) && return :τ_forward_out
    return nothing
end

"""
Summarise a [`forecast_reported`](@ref) result into a `DataFrame` with
one row per projected quantity (new infections, new BVD deaths, and the
new confirmed cases / confirmed deaths when present), reporting the same
equal-tailed 30/60/90% credible interval endpoints (`lower_90 …
upper_90`) as the other summary tables. Every quantity is the new count
over the coming week.
"""
function forecast_table(fc::DataFrame; digits::Integer = 0)
    _row(label,
        draws) = begin
        s = posterior_summary(draws)
        (quantity = label,
            lower_90 = round(s.lo90; digits), lower_60 = round(s.lo60; digits),
            lower_30 = round(s.lo30; digits), upper_30 = round(s.hi30; digits),
            upper_60 = round(s.hi60; digits), upper_90 = round(s.hi90; digits))
    end
    rows = NamedTuple[
    _row("New infections", fc[!, :infections_new]),
    _row(
        "New BVD deaths (DRC)", fc[!, :bvd_deaths_new])]
    :exports_new in propertynames(fc) &&
        push!(rows, _row("New exports (Uganda)", fc[!, :exports_new]))
    :confirmed_new in propertynames(fc) &&
        push!(rows, _row("New confirmed cases (DRC)", fc[!, :confirmed_new]))
    :confirmed_deaths_new in propertynames(fc) &&
        push!(rows,
            _row("New confirmed deaths (DRC)", fc[!, :confirmed_deaths_new]))
    return _prettify(DataFrame(rows))
end

"""
Validate a [`forecast_reported`](@ref) projection against the counts that
were later observed, for the two directly observable forecast targets:
laboratory-confirmed cases and laboratory-confirmed deaths. Fit the joint
at an earlier cut-off, forecast forward to the current cut-off, and pass
the now-known truth here to score the prediction.

`fc` is a [`forecast_reported`](@ref) result carrying the per-draw new
counts over the forecast horizon (`:confirmed_new`, `:confirmed_deaths_new`).
`baseline_confirmed` / `baseline_confirmed_deaths` are the cut-off
cumulative counts the forecast started from (the confirmed cases / deaths
observed by the earlier as-of date); the predicted cumulative at the
target date is `baseline + new`. `confirmed` / `confirmed_deaths` are the
cumulative counts actually observed at the target date.

Returns a `DataFrame` with one row per (quantity, horizon-view): the
cumulative and the new-over-horizon prediction for each of confirmed cases
and confirmed deaths, giving the observed count, the equal-tailed
30/60/90% predictive intervals (the same endpoints as the other summary
tables) and whether the observed count falls inside the 90% interval.
Latent infections and all-BVD deaths are not scored here: neither is
directly observed, so neither can be validated against data.
"""
function forecast_vs_truth(fc::DataFrame;
        confirmed::Real, confirmed_deaths::Real,
        baseline_confirmed::Real, baseline_confirmed_deaths::Real,
        digits::Integer = 0)
    _row(label,
        draws,
        obs) = begin
        s = posterior_summary(draws)
        lo = round(s.lo90; digits)
        hi = round(s.hi90; digits)
        (quantity = label, observed = round(obs; digits),
            lower_90 = lo, lower_60 = round(s.lo60; digits),
            lower_30 = round(s.lo30; digits), upper_30 = round(s.hi30; digits),
            upper_60 = round(s.hi60; digits), upper_90 = hi,
            within_90 = lo <= obs <= hi ? "yes" : "no")
    end
    :confirmed_new in propertynames(fc) ||
        error("forecast_vs_truth needs a :confirmed_new column; pass a " *
              "forecast_reported result fitted with the laboratory streams")
    :confirmed_deaths_new in propertynames(fc) ||
        error("forecast_vs_truth needs a :confirmed_deaths_new column; " *
              "pass a forecast_reported result fitted with the death-lab " *
              "stream")
    new_cases = fc[!, :confirmed_new]
    new_deaths = fc[!, :confirmed_deaths_new]
    cum_cases = new_cases .+ baseline_confirmed
    cum_deaths = new_deaths .+ baseline_confirmed_deaths
    rows = NamedTuple[
    _row("Confirmed cases (DRC), cumulative", cum_cases, confirmed),
    _row(
        "Confirmed cases (DRC), new", new_cases,
        confirmed - baseline_confirmed),
    _row(
        "Confirmed deaths (DRC), cumulative", cum_deaths,
        confirmed_deaths),
    _row(
        "Confirmed deaths (DRC), new", new_deaths,
        confirmed_deaths - baseline_confirmed_deaths)
]
    return _prettify(DataFrame(rows))
end

"""
Validate a fitted chain's confirmed-stream projection against the full
observed daily trajectory rather than a single endpoint.

Given the per-vintage observed cumulative confirmed series (`dates`, with
`confirmed` and `confirmed_deaths` the cumulative confirmed cases /
deaths at each sitrep date) and the fit's data cut-off `snapshot_date`,
project the chain forward to every vintage date that falls *after* the
cut-off and compare the predicted cumulative count against the observed
one. This scores the whole forecast horizon, not just its endpoint.

`baseline_confirmed` / `baseline_confirmed_deaths` are the fit's cut-off
cumulative confirmed counts (the forecast origin), and
`baseline_analysed` the cut-off samples-analysed total that anchors the
laboratory queue. The remaining keywords match [`forecast_reported`](@ref).

Returns a `DataFrame` with one row per (quantity, date): the horizon in
days, the observed count, the equal-tailed 30/60/90% predictive intervals
and whether the observed count falls within the 90% interval.
"""
function forecast_vs_truth_trajectory(chn;
        dates::AbstractVector,
        confirmed::AbstractVector,
        confirmed_deaths::AbstractVector,
        snapshot_date,
        baseline_confirmed::Real,
        baseline_confirmed_deaths::Real,
        baseline_analysed::Real,
        daily_travellers::Real = 0,
        source_population::Real = 1,
        seed::Integer = 20260520,
        report_onset_offset::Union{Nothing, Real} = nothing,
        digits::Integer = 0,
        alg = DEATH_INTEGRAL_ALG)
    length(dates) == length(confirmed) == length(confirmed_deaths) ||
        error("dates, confirmed and confirmed_deaths must be equal length " *
              "(got $(length(dates)), $(length(confirmed)), " *
              "$(length(confirmed_deaths)))")
    snap = date2epochdays(Date(snapshot_date))
    rows = NamedTuple[]
    for i in eachindex(dates)
        ## Only vintages strictly after the fit cut-off are a forecast.
        h = date2epochdays(Date(dates[i])) - snap
        h > 0 || continue
        fc = forecast_reported(chn; horizon = h,
            daily_travellers, source_population,
            obs_confirmed = baseline_confirmed,
            obs_confirmed_deaths = baseline_confirmed_deaths,
            obs_analysed = baseline_analysed,
            forecast_exports = false, seed, report_onset_offset, alg)
        for (label, col, base, obs) in (
            ("Confirmed cases (DRC)", :confirmed_new,
                baseline_confirmed, confirmed[i]),
            ("Confirmed deaths (DRC)", :confirmed_deaths_new,
                baseline_confirmed_deaths, confirmed_deaths[i]))
            col in propertynames(fc) || continue
            s = posterior_summary(fc[!, col] .+ base)
            lo = round(s.lo90; digits)
            hi = round(s.hi90; digits)
            push!(rows,
                (quantity = label, date = string(dates[i]),
                    horizon_days = h, observed = round(obs; digits),
                    lower_90 = lo, lower_60 = round(s.lo60; digits),
                    lower_30 = round(s.lo30; digits),
                    upper_30 = round(s.hi30; digits),
                    upper_60 = round(s.hi60; digits), upper_90 = hi,
                    within_90 = lo <= obs <= hi ? "yes" : "no"))
        end
    end
    return _prettify(DataFrame(rows))
end
