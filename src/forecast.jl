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
# contribution plus the non-BVD background. The BVD contribution reuses
# `delay_convolution` as a delay-convolved cumulative integrator at unit
# ascertainment to compute `∫₀^{Th} exp(r·s) · f_rep(Th-s) ds`, scaled by
# `onset_fraction` (the incubation mgf) since the latent trajectory is
# infections; the background contribution is the ramp cumulative
# `μ_bg(Th)` (`λ0 · Th` when `Δλ = 0`), non-BVD and unscaled.
function _forecast_cases_mean(r, Th, α_rep, θ_rep, p_drc, bg::BackgroundRamp;
        onset_fraction = 1.0, alg = DEATH_INTEGRAL_ALG)
    conv = delay_convolution(one(p_drc), r, Th, Gamma(α_rep, θ_rep); alg)
    ## Background contribution is the ramp cumulative μ_bg(Th); with Δλ = 0
    ## this is the old constant `λ_bg · Th`.
    return onset_fraction * p_drc * conv + bg_cumulative(bg, Th)
end

function _forecast_confirmed_mean(r, Th, α_rep, θ_rep, α_lab, θ_lab,
        p_drc, s_test, τ_test; onset_fraction = 1.0, alg = DEATH_INTEGRAL_ALG)
    d_rep = Gamma(α_rep, θ_rep)
    d_lab = Gamma(α_lab, θ_lab)
    ## Inner BVD-reported trajectory in closed form; only the outer lab
    ## convolution is a quadrature. Omitting `alg` selects the Gamma
    ## analytic `delay_convolution` method.
    bvd_reported_at = let r = r, p_drc = p_drc, d_rep = d_rep
        u -> p_drc * delay_convolution(one(p_drc), r, u, d_rep)
    end
    return onset_fraction * s_test * τ_test *
           delay_convolution(bvd_reported_at, Th, d_lab; alg)
end

# Cumulative tests analysed at horizon `Th`: τ_test gates both the BVD
# contribution (lab-completed BVD samples, scaled by `onset_fraction`) and
# the non-BVD background (saturating-ramp rate λ_bg(u) convolved against
# F_lab, unscaled).
function _forecast_tests_mean(r, Th, α_rep, θ_rep, α_lab, θ_lab,
        p_drc, bg::BackgroundRamp, τ_test;
        onset_fraction = 1.0, alg = DEATH_INTEGRAL_ALG)
    d_rep = Gamma(α_rep, θ_rep)
    d_lab = Gamma(α_lab, θ_lab)
    ## Inner BVD-reported trajectory in closed form; only the outer lab
    ## convolution is a quadrature. Omitting `alg` selects the Gamma
    ## analytic `delay_convolution` method.
    bvd_reported_at = let r = r, p_drc = p_drc, d_rep = d_rep
        u -> p_drc * delay_convolution(one(p_drc), r, u, d_rep)
    end
    bvd_tested = onset_fraction *
                 delay_convolution(bvd_reported_at, Th, d_lab; alg)
    ## Background tested volume at `Th`: ∫₀^Th λ_bg(u) F_lab(Th - u) du, the
    ## ramp arrivals convolved against the lab-delay CDF. With Δλ = 0 this
    ## reduces to the constant `λ0 · ∫₀^Th F_lab(v) dv`.
    bg_tested = bg_tested_integral(bg, α_lab, θ_lab, Th; alg)
    return τ_test * (bvd_tested + bg_tested)
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
  when the chain carries the lab-turnaround delay (`:α_lab`, `:θ_lab`)
  and `obs_confirmed` is supplied. Otherwise these columns are absent.

DRC reported cases follow the additive expectation
`p_drc · ∫₀^{T+h} exp(r·s) · f_rep(T+h-s) ds + μ_bg(T+h)`, with
`f_rep = Gamma(α_rep, θ_rep)` for the BVD-driven contribution and
`μ_bg` the saturating-ramp non-BVD background cumulative (`λ0 · t` when
`Δλ = 0`). Laboratory-confirmed cases
follow `s_test · p_drc · ∫₀^{T+h} exp(r·s) · f_conf(T+h-s) ds`, with
`f_conf` the moment-matched Gamma of `f_rep ∗ f_lab` and `s_test` the
PCR sensitivity. Exports use `p_uganda · q` with
`q = daily_travellers / source_population`. Assumes growth continues
unchanged over the horizon (no interventions, no saturation).

`report_onset_offset` anchors the background ramp to reporting onset
(`t_report = T − report_onset_offset`), matching the fit; pass the same
value used in [`bvd_joint`](@ref) (see [`report_onset_offset`](@ref)).
The default `nothing` keeps the seeding-anchored clock (`t_report = 0`).
"""
function forecast_reported(chn;
        horizon::Real = 7,
        daily_travellers::Real,
        source_population::Real,
        obs_cases::Real,
        obs_deaths::Real,
        obs_exports::Real,
        obs_confirmed::Union{Real, Missing} = missing,
        obs_tests::Union{Real, Missing} = missing,
        seed::Integer = 20260520,
        report_onset_offset::Union{Nothing, Real} = nothing,
        alg = DEATH_INTEGRAL_ALG)
    r = _draws(chn, :r)
    T = _draws(chn, :T)
    CFR = _draws(chn, :CFR)
    α = _draws(chn, :α)
    θ = _draws(chn, :θ)
    w = _draws(chn, :w)
    pr = _draws(chn, :p_drc)
    pu = _draws(chn, :p_uganda)
    k = _draws(chn, :k)
    α_rep = _draws(chn, :α_rep)
    θ_rep = _draws(chn, :θ_rep)
    ## Background ramp draws. New chains carry the baseline `λ0` and ramp
    ## amplitude `Δλ`; chains predating the ramp carry only the constant
    ## `λ_bg`, recovered here as `λ0` with `Δλ = 0` (so `bg_*` reduce to the
    ## old constant-rate forms).
    if haskey_chain(chn, :λ0)
        λ0 = _draws(chn, :λ0)
        Δλ = haskey_chain(chn, :Δλ) ? _draws(chn, :Δλ) : zeros(length(λ0))
    else
        λ0 = _draws(chn, :λ_bg)
        Δλ = zeros(length(λ0))
    end
    ## Incubation draws map the latent infection trajectory onto onsets
    ## (the `onset_fraction = mgf(incubation, −r)` of the observation
    ## models). Absent on chains predating the infection layer, where
    ## `onset_fraction = 1` recovers the previous behaviour.
    has_incubation = all(haskey_chain(chn, n) for n in (:α_inc, :θ_inc))
    α_inc = has_incubation ? _draws(chn, :α_inc) : nothing
    θ_inc = has_incubation ? _draws(chn, :θ_inc) : nothing
    ## Lab-turnaround and PCR sensitivity draws live on the joint chain
    ## only; their absence drops the confirmed-cases columns.
    has_lab = all(haskey_chain(chn, n) for n in (:α_lab, :θ_lab, :s_test)) &&
              obs_confirmed !== missing
    α_lab = has_lab ? _draws(chn, :α_lab) : nothing
    θ_lab = has_lab ? _draws(chn, :θ_lab) : nothing
    s_test = has_lab ? _draws(chn, :s_test) : nothing
    ## τ_test is optional on the chain: when absent (older fits or
    ## synthetic chains predating the testing-rate extension) fall back
    ## to τ = 1, matching the previous "all suspected get tested"
    ## assumption.
    τ_test_draws = (has_lab && haskey_chain(chn, :τ_test)) ?
                   _draws(chn, :τ_test) :
                   (has_lab ? ones(length(r)) : nothing)

    has_tests = has_lab && obs_tests !== missing
    rng = MersenneTwister(seed)
    n = length(r)
    q = daily_travellers / source_population
    cases_cum = Vector{Int}(undef, n)
    deaths_cum = Vector{Int}(undef, n)
    exports_cum = Vector{Int}(undef, n)
    confirmed_cum = has_lab ? Vector{Int}(undef, n) : nothing
    tests_cum = has_tests ? Vector{Int}(undef, n) : nothing

    @inbounds for i in 1:n
        Th = T[i] + horizon
        os = has_incubation ?
             onset_rescale(Gamma(α_inc[i], θ_inc[i]), r[i]) : 1.0
        ## DRC reported cases: onset_fraction · p_drc · ∫₀^{T+h} exp(r·s) ·
        ## f_rep(T+h-s) ds + μ_bg(T+h). The background ramp is anchored to
        ## reporting onset, `t_report = T − report_onset_offset` (clamped at
        ## 0); `nothing` keeps the seeding-anchored clock (`t_report = 0`).
        t_report_i = report_onset_offset === nothing ? zero(T[i]) :
                     max(T[i] - report_onset_offset, zero(T[i]))
        bg_i = background_ramp(λ0[i], Δλ[i], BACKGROUND_RAMP_SCALE,
            t_report_i)
        μ_cases = _forecast_cases_mean(r[i], Th, α_rep[i], θ_rep[i],
            pr[i], bg_i; onset_fraction = os, alg)
        cases_cum[i] = _nb_rand(rng, k[i], μ_cases)
        ## DRC deaths: onset_fraction · CFR · ∫_0^{T+h} exp(r·s) f(T+h−s) ds.
        μ_deaths = _forecast_deaths_mean(r[i], Th, α[i], θ[i], CFR[i];
            onset_fraction = os, alg)
        deaths_cum[i] = _nb_rand(rng, k[i], μ_deaths)
        ## Uganda exports: p_uganda · q · ∫_{T+h−w}^{T+h} C(s) ds (closed
        ## form for exponential growth). Exports see infections directly —
        ## the detection window absorbs the infection-to-detection delay —
        ## so they are not rescaled by the incubation period.
        lo = max(Th - w[i], zero(Th))
        μ_exports = pu[i] * q * (exp(r[i] * Th) - exp(r[i] * lo)) / r[i]
        exports_cum[i] = rand(rng, Poisson(max(μ_exports, eps(μ_exports))))
        if has_lab
            μ_confirmed = _forecast_confirmed_mean(r[i], Th,
                α_rep[i], θ_rep[i], α_lab[i], θ_lab[i], pr[i],
                s_test[i], τ_test_draws[i]; onset_fraction = os, alg)
            confirmed_cum[i] = _nb_rand(rng, k[i], μ_confirmed)
            if has_tests
                μ_tests = _forecast_tests_mean(r[i], Th,
                    α_rep[i], θ_rep[i], α_lab[i], θ_lab[i],
                    pr[i], bg_i, τ_test_draws[i]; onset_fraction = os, alg)
                tests_cum[i] = _nb_rand(rng, k[i], μ_tests)
            end
        end
    end

    _new(cum, obs) = max.(cum .- round(Int, obs), 0)
    df = DataFrame(
        cases_cum = cases_cum,
        deaths_cum = deaths_cum,
        exports_cum = exports_cum,
        cases_new = _new(cases_cum, obs_cases),
        deaths_new = _new(deaths_cum, obs_deaths),
        exports_new = _new(exports_cum, obs_exports)
    )
    if has_lab
        df.confirmed_cum = confirmed_cum
        df.confirmed_new = _new(confirmed_cum, obs_confirmed)
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
        ("DRC deaths", :deaths_cum, :deaths_new),
        ("Uganda exports", :exports_cum, :exports_new)]
    :tests_cum in propertynames(fc) &&
        push!(streams, ("DRC tests analysed", :tests_cum, :tests_new))
    :confirmed_cum in propertynames(fc) &&
        push!(streams, ("DRC confirmed cases",
            :confirmed_cum, :confirmed_new))
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
