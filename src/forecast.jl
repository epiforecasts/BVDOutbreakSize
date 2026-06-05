# One-week-ahead posterior-predictive forecast. Continues the fitted
# exponential growth `h` days past the cut-off `T` and applies the same
# observation models to forecast the cumulative reported cases (DRC),
# deaths (DRC) and exports (Uganda) by `T + h`, plus the new counts
# expected over the coming week. Each draw produces a replicated integer
# count so the intervals include both parameter and observation
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
forecast. For each draw, continue exponential growth to `T + horizon`
and apply the observation models, returning a `DataFrame` with one row
per draw and columns:

- `:cases_cum`, `:deaths_cum`, `:exports_cum` — replicated cumulative
  counts reported by `T + horizon`.
- `:cases_new`, `:deaths_new`, `:exports_new` — new counts over the
  coming week (`*_cum` minus the corresponding observed count at `T`,
  floored at zero).
- `:confirmed_cum`, `:confirmed_new` — laboratory-confirmed counterparts
  when the chain carries the confirmed-stream parameters
  (`:s_test`, `:spec_test`, `:τ_forward`, `:α_recv`, `:θ_recv`,
  `:capacity_cutoff`) and `obs_confirmed` and `obs_analysed` are supplied.
  Otherwise these columns are absent.
- `:confirmed_deaths_cum`, `:confirmed_deaths_new` — laboratory-confirmed
  deaths when the chain additionally carries the confirmed-death
  parameters (`:τ_death`, `:p_deaths`, `:λ_bg_death`, sharing the case-lab
  `:s_test` / `:spec_test`) and `obs_confirmed_deaths` is supplied.
  Otherwise these columns are absent.

DRC reported cases follow the additive expectation
`p_drc · ∫₀^{T+h} exp(r·s) · f_rep(T+h-s) ds + λ_bg·(T+h)`, with
`f_rep = Gamma(α_rep, θ_rep)` for the BVD-driven contribution and the
constant-rate non-BVD background. Laboratory-confirmed cases continue the
queue of [`confirmed_cases_model`](@ref): the received backlog
(`τ_forward · N_susp` convolved with the receipt delay) is analysed at the
cut-off capacity `κ` over the horizon, `ΔA = backlog·(1 − exp(−κ·h/backlog))`,
and the new positives are `ΔA` times the composition-linked positivity,
added to the observed cumulative confirmed. Exports use `p_uganda · q` with
`q = daily_travellers / source_population`; pass `forecast_exports = false`
to drop the export columns (e.g. when cross-border travel is disrupted so
the forward travel rate no longer holds). Assumes growth continues
unchanged over the horizon (no interventions, no saturation).

`report_onset_offset` is accepted for call-site compatibility with the
fit; the composition-linked horizon positivity does not use it.
The default `nothing` keeps the seeding-anchored clock (`t_report = 0`).
"""
function forecast_reported(chn;
        horizon::Real = 7,
        daily_travellers::Real,
        source_population::Real,
        obs_cases::Real,
        obs_deaths::Real,
        obs_exports::Union{Real, Missing} = missing,
        obs_confirmed::Union{Real, Missing} = missing,
        obs_confirmed_deaths::Union{Real, Missing} = missing,
        obs_tests::Union{Real, Missing} = missing,
        obs_analysed::Union{Real, Missing} = missing,
        forecast_exports::Bool = true,
        seed::Integer = 20260520,
        report_onset_offset::Union{Nothing, Real} = nothing,
        alg = DEATH_INTEGRAL_ALG)
    r = _draws(chn, :r)
    T = _draws(chn, :T)
    CFR = _draws(chn, :CFR)
    α = _draws(chn, :α)
    θ = _draws(chn, :θ)
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
    has_lab = all(haskey_chain(chn, n)
              for n in (:s_test, :spec_test, :τ_forward, :α_recv, :θ_recv,
                  :capacity_cutoff)) &&
              obs_confirmed !== missing && obs_analysed !== missing
    s_test = has_lab ? _draws(chn, :s_test) : nothing
    spec_test = has_lab ? _draws(chn, :spec_test) : nothing
    τ_forward = has_lab ? _draws(chn, :τ_forward) : nothing
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

    has_tests = has_lab && obs_tests !== missing
    rng = MersenneTwister(seed)
    n = length(r)
    q = daily_travellers / source_population
    cases_cum = Vector{Int}(undef, n)
    deaths_cum = Vector{Int}(undef, n)
    exports_cum = forecast_exports ? Vector{Int}(undef, n) : nothing
    confirmed_cum = has_lab ? Vector{Int}(undef, n) : nothing
    confirmed_deaths_cum = has_lab_deaths ? Vector{Int}(undef, n) : nothing
    tests_cum = has_tests ? Vector{Int}(undef, n) : nothing

    @inbounds for i in 1:n
        Th = T[i] + horizon
        os = has_incubation ?
             onset_rescale(Gamma(α_inc[i], θ_inc[i]), r[i]) : 1.0
        ## DRC reported cases: onset_fraction · p_drc · ∫₀^{T+h} exp(r·s) ·
        ## f_rep(T+h-s) ds + λ_bg·(T+h) (constant-rate background).
        μ_cases = _forecast_cases_mean(r[i], Th, α_rep[i], θ_rep[i],
            pr[i], λ_bg[i]; onset_fraction = os, alg)
        cases_cum[i] = _nb_rand(rng, k[i], μ_cases)
        ## DRC deaths: onset_fraction · CFR · ∫_0^{T+h} exp(r·s) f(T+h−s) ds.
        μ_deaths = _forecast_deaths_mean(r[i], Th, α[i], θ[i], CFR[i];
            onset_fraction = os, alg)
        deaths_cum[i] = _nb_rand(rng, k[i], μ_deaths)
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
            exports_cum[i] = rand(rng,
                Poisson(max(μ_exports, eps(μ_exports))))
        end
        if has_lab
            ## Received backlog (receipt-delayed) and the capacity-limited
            ## analysed increment over the horizon; new confirmed positives
            ## are the horizon positivity times that increment, added to the
            ## observed cumulative confirmed at the cut-off.
            f_receipt = Gamma(α_recv[i], θ_recv[i])
            received = _forecast_received(r[i], Th, α_rep[i], θ_rep[i],
                pr[i], λ_bg[i], τ_forward[i], f_receipt; onset_fraction = os)
            Δanalysed = _forecast_analysed_increment(received, obs_analysed,
                capacity[i], horizon)
            p_pos = _forecast_positivity(r[i], Th, α_rep[i], θ_rep[i],
                pr[i], λ_bg[i], s_test[i], spec_test[i]; onset_fraction = os)
            confirmed_cum[i] = round(Int, obs_confirmed) +
                               _nb_rand(rng, k[i], p_pos * Δanalysed)
            if has_tests
                tests_cum[i] = _nb_rand(rng, k[i], received)
            end
        end
        if has_lab_deaths
            ## Forwarded positive increment of the suspect-death backlog
            ## over the horizon, sharing the case-lab sensitivity /
            ## specificity, added to the observed cumulative confirmed
            ## deaths at the cut-off.
            μ_cdeath = _forecast_confirmed_deaths_increment(r[i], T[i], Th,
                α[i], θ[i], CFR[i], p_deaths[i], λ_bg_death[i], τ_death[i],
                s_test[i], spec_test[i]; onset_fraction = os, alg)
            confirmed_deaths_cum[i] = round(Int, obs_confirmed_deaths) +
                                      _nb_rand(rng, k[i], μ_cdeath)
        end
    end

    _new(cum, obs) = max.(cum .- round(Int, obs), 0)
    df = DataFrame(
        cases_cum = cases_cum,
        deaths_cum = deaths_cum,
        cases_new = _new(cases_cum, obs_cases),
        deaths_new = _new(deaths_cum, obs_deaths)
    )
    if forecast_exports
        df.exports_cum = exports_cum
        df.exports_new = _new(exports_cum, obs_exports)
    end
    if has_lab
        df.confirmed_cum = confirmed_cum
        df.confirmed_new = _new(confirmed_cum, obs_confirmed)
    end
    if has_lab_deaths
        df.confirmed_deaths_cum = confirmed_deaths_cum
        df.confirmed_deaths_new = _new(confirmed_deaths_cum,
            obs_confirmed_deaths)
    end
    if has_tests
        df.tests_cum = tests_cum
        df.tests_new = _new(tests_cum, obs_tests)
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

"""
Summarise a [`forecast_reported`](@ref) result into a `DataFrame` with
one row per stream (cases, deaths, exports) and quantity (cumulative
total by `T + 7`, or new this week), reporting the same equal-tailed
30/60/90% credible interval endpoints (`lower_90 … upper_90`) as the
other summary tables.
"""
function forecast_table(fc::DataFrame; digits::Integer = 0)
    _row(label,
        quantity,
        draws) = begin
        s = posterior_summary(draws)
        (stream = label, quantity = quantity,
            lower_90 = round(s.lo90; digits), lower_60 = round(s.lo60; digits),
            lower_30 = round(s.lo30; digits), upper_30 = round(s.hi30; digits),
            upper_60 = round(s.hi60; digits), upper_90 = round(s.hi90; digits))
    end
    rows = NamedTuple[]
    streams = [
        ("DRC reported cases", :cases_cum, :cases_new),
        ("DRC deaths", :deaths_cum, :deaths_new)]
    :exports_cum in propertynames(fc) &&
        push!(streams, ("Uganda exports", :exports_cum, :exports_new))
    :tests_cum in propertynames(fc) &&
        push!(streams, ("DRC tests analysed", :tests_cum, :tests_new))
    :confirmed_cum in propertynames(fc) &&
        push!(streams, ("DRC confirmed cases",
            :confirmed_cum, :confirmed_new))
    :confirmed_deaths_cum in propertynames(fc) &&
        push!(streams, ("DRC confirmed deaths",
            :confirmed_deaths_cum, :confirmed_deaths_new))
    for (label, cum, new) in streams
        push!(rows, _row(label, "cumulative by T+7", fc[!, cum]))
        push!(rows, _row(label, "new this week", fc[!, new]))
    end
    return _prettify(DataFrame(rows))
end

"""
Validate a [`forecast_reported`](@ref) projection against the counts
that were later observed. `cases`, `deaths` and `exports` are the
observed cumulative DRC reported cases, DRC deaths and Uganda exports at
the forecast target date. Returns a `DataFrame` with one row per stream
giving the observed count, the equal-tailed 30/60/90% predictive
intervals (the same endpoints as the other summary tables), and whether
the observed count falls inside the 90% interval. Use it to forecast
from an earlier data snapshot and check the now-known truth against the
projection.
"""
function forecast_vs_truth(fc::DataFrame;
        cases::Real, deaths::Real, exports::Real,
        confirmed::Union{Real, Missing} = missing,
        tests::Union{Real, Missing} = missing,
        digits::Integer = 0)
    _row(label,
        draws,
        obs) = begin
        s = posterior_summary(draws)
        lo = round(s.lo90; digits)
        hi = round(s.hi90; digits)
        (stream = label, observed = round(obs; digits),
            lower_90 = lo, lower_60 = round(s.lo60; digits),
            lower_30 = round(s.lo30; digits), upper_30 = round(s.hi30; digits),
            upper_60 = round(s.hi60; digits), upper_90 = hi,
            within_90 = lo <= obs <= hi ? "yes" : "no")
    end
    rows = NamedTuple[
    _row("DRC reported cases", fc[!, :cases_cum], cases),
    _row(
        "DRC deaths", fc[!, :deaths_cum], deaths),
    _row(
        "Uganda exports", fc[!, :exports_cum], exports)
]
    tests !== missing && :tests_cum in propertynames(fc) &&
        push!(rows,
            _row("DRC tests analysed", fc[!, :tests_cum], tests))
    confirmed !== missing && :confirmed_cum in propertynames(fc) &&
        push!(rows,
            _row("DRC confirmed cases", fc[!, :confirmed_cum], confirmed))
    return _prettify(DataFrame(rows))
end

"""
Validate a fitted chain's projection against the full observed daily
trajectory rather than a single endpoint.

Given the per-vintage observed cumulative series (`dates`, with `cases`
and `deaths` the cumulative DRC reported cases and suspected deaths at
each date) and the fit's data cut-off `snapshot_date`, project the chain
forward to every vintage date that falls *after* the cut-off and compare
the predicted cumulative count against the observed one. This scores the
whole forecast horizon, not just its endpoint.

Returns a `DataFrame` with one row per (stream, date): the horizon in
days, the observed count, the equal-tailed 30/60/90% predictive
intervals (the same endpoints as the other tables) and whether the
observed count falls within the 90% interval. Streams with no
per-vintage series (Uganda exports) stay with [`forecast_vs_truth`](@ref).

`baseline_cases` / `baseline_deaths` are the fit's cut-off counts; they
only feed the unused `*_new` columns of the per-horizon projection.
"""
function forecast_vs_truth_trajectory(chn;
        dates::AbstractVector,
        cases::AbstractVector,
        deaths::AbstractVector,
        snapshot_date,
        daily_travellers::Real,
        source_population::Real,
        baseline_cases::Real = 0,
        baseline_deaths::Real = 0,
        seed::Integer = 20260520,
        report_onset_offset::Union{Nothing, Real} = nothing,
        digits::Integer = 0,
        alg = DEATH_INTEGRAL_ALG)
    length(dates) == length(cases) == length(deaths) ||
        error("dates, cases and deaths must be equal length (got " *
              "$(length(dates)), $(length(cases)), $(length(deaths)))")
    snap = date2epochdays(Date(snapshot_date))
    rows = NamedTuple[]
    for i in eachindex(dates)
        ## Only vintages strictly after the fit cut-off are a forecast.
        h = date2epochdays(Date(dates[i])) - snap
        h > 0 || continue
        fc = forecast_reported(chn; horizon = h,
            daily_travellers, source_population,
            obs_cases = baseline_cases, obs_deaths = baseline_deaths,
            obs_exports = 0, seed, report_onset_offset, alg)
        for (label,
            col,
            obs) in (
            ("DRC reported cases", :cases_cum, cases[i]),
            ("DRC deaths", :deaths_cum, deaths[i]))
            s = posterior_summary(fc[!, col])
            lo = round(s.lo90; digits)
            hi = round(s.hi90; digits)
            push!(rows,
                (stream = label, date = string(dates[i]),
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
