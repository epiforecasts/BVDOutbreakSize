# Joint composer models: build the full generative model for each analysis
# by running the generating infection process once, staging it to daily
# onset incidence, and routing the shared onsets into the relevant
# observation submodels. Single-stream composers condition on one stream
# each. [`bvd_joint`](@ref) conditions on all streams plus the optional
# genetic seeding bound. Any count passed as `missing` is dropped, so the
# composers double as prior- and posterior-predictive generators. With a
# `forecast` horizon each composer runs past its cut-off and draws its
# streams' future counts, which is how every forecast is made.
#
# Submodels whose `:=` deterministics are re-exposed at composer level are
# attached with a prefixed `to_submodel(x)`. Attaching them with
# `to_submodel(x, false)` re-introduces the nested `:=` names at the parent
# and trips Turing's MustNotOverwriteError on the duplicates.

## Run the generating infection process and onset staging, returning the
## infection state and the daily onsets shared by every stream. With a
## `forecast` horizon both run past the cut-off `n`, and the future
## infections, onsets and reproduction number are tracked.
@model function _latent(
        n::Integer, breakpoint, infection, onset_incidence;
        rt_start::Integer = 1, rt_walk_start::Integer = rt_start,
        forecast::Union{Nothing, ForecastHorizon} = nothing
    )
    fkw = forecast === nothing ? (;) : (; forecast)
    infection_state ~ to_submodel(
        infection(n; breakpoint, rt_start, rt_walk_start, fkw...), false
    )
    onset_state ~ to_submodel(
        onset_incidence(infection_state.infections), false
    )
    cumulative_infections := infection_state.cumulative
    C_T := infection_state.C_T
    cumulative_onsets := cumsum(onset_state.onsets)
    ## The susceptible fraction reaches only `:=` quantities, so it is built
    ## only when they are recorded.
    if _reporting(__varinfo__)
        fraction = pool_fraction(
            infection_state.cumulative, infection_state.population
        )
        susceptible_fraction := fraction[1:n]
        ## The forecast reproduction number is net of depletion.
        if forecast !== nothing
            fd = forecast_days(n, forecast)
            forecast_rt := adjusted_rt(infection_state.Rt, fraction, fd)
        end
    end
    if forecast !== nothing
        fd = forecast_days(n, forecast)
        forecast_infections := infection_state.infections[fd]
        forecast_onsets := onset_state.onsets[fd]
    end
    return (;
        infection_state, onsets = onset_state.onsets,
        incubation_pmf = onset_state.incubation_pmf,
    )
end

## Keywords handing one stream its simulated observations, keyed by the
## stream's name in the joint (`simulate_recovery`). Nothing is passed for a
## fit to real data.
_sim_kw(::Nothing, ::Symbol) = (;)
function _sim_kw(simulated_data, state::Symbol)
    return haskey(simulated_data, state) ?
        (; simulated = simulated_data[state]) : (;)
end

"""
Grid days a `forecast` covers past the cut-off `n`, `n + 1` to
`n + horizon`.
"""
forecast_days(n::Integer, forecast) = (n + 1):(n + horizon_days(forecast))

## Future Uganda exports: the per-day Poisson the dated export series is
## scored with, at the modelled export intensity on each forecast day.
_forecast_exports(state, fd) = dated_poisson_model(
    state.export_prevalence[fd], missing
)

## National forecast counts summed over each future week, the totals the
## province split divides. `edges` is the cut-off followed by the future
## vintage days. A week past `typemax(Int)` saturates there.
function _weekly_totals(daily, n, edges)
    return [
        Int(_saturating_sum(view(daily, (edges[j] - n + 1):(edges[j + 1] - n))))
            for j in 1:(length(edges) - 1)
    ]
end

## Future counts of one stream: its modelled daily series on the forecast
## days, drawn one day at a time through the stream's own negative binomial.
## The counts are `missing`, so `predict` draws them.
_forecast_counts(daily, fd, k) = vintage_increments_model(daily[fd], missing, k)

"""
Future observations of the isolation and treatment flows, from a fitted
[`treatment_flow_model`](@ref) state run past the cut-off. Each is drawn
through the likelihood the fitted days are scored with. The occupancy is
the censored negative binomial around the modelled demand plus the
reclassification offset. Its cap is the modelled bed capacity, which the
capacity walk carries past the cut-off, floored at the last fitted cap
([`censoring_cap`](@ref)). Admissions are censored at the free-bed headroom,
that cap less the previous day's occupancy: the last recorded occupancy on
the first day, then the drawn one ([`admission_headroom`](@ref)). In-care
deaths and rule-outs are daily negative binomials. The latent bed demand and
capacity are tracked alongside. With no recorded capacity or occupancy the
fitted days are uncensored, and so are the future ones.
"""
@model function treatment_forecast_model(
        state, fd, capacity_history, isolation_history, k
    )
    nocap = 1.0e6
    have_cap = !isempty(capacity_history.counts)
    have_occ = !isempty(isolation_history.counts)
    last_cap = have_occ ?
        only(
            censoring_cap(
                isolation_history.days[end:end],
                isolation_history.counts[end:end], capacity_history
            )
        ) : 0.0
    ceilings = have_cap ?
        [max(float(c), last_cap) for c in state.C[fd]] :
        fill(nocap, length(fd))
    occupancy = state.occupancy_mean[fd]
    forecast_isolation ~ to_submodel(
        censored_occupancy_model(occupancy, ceilings, missing, k)
    )
    occ = forecast_isolation.obs
    head = if have_cap && have_occ
        prev = vcat(
            float(isolation_history.counts[end]),
            [float(occ[j]) for j in 1:(length(fd) - 1)]
        )
        max.(ceilings .- prev, 0.5)
    else
        fill(nocap, length(fd))
    end
    admissions = state.admit_daily[fd]
    forecast_admissions ~ to_submodel(
        censored_occupancy_model(admissions, head, missing, k)
    )
    forecast_incare_deaths ~ to_submodel(
        _forecast_counts(state.deaths_daily, fd, k)
    )
    forecast_ruleouts ~ to_submodel(
        _forecast_counts(state.ruleout_daily, fd, k)
    )
    forecast_bed_demand := state.demand[fd]
    forecast_bed_capacity := state.C[fd]
    return (;
        isolation = occupancy, admissions, occupancy = occ,
        incare_deaths = state.deaths_daily[fd],
        ruleouts = state.ruleout_daily[fd],
    )
end

## Reported onset total as of `as_of` under a fitted onset-reporting state.
_onset_total(onsets, state, as_of) = onset_report_expected_total(
    onsets, state.logit_h0, state.γ, state.grid_start, state.alpha, as_of
)

"""
Future vintages of the symptom-onset reporting triangle. Each future
vintage `v` in `vintages` prints a total the fitted reporting hazard puts at
[`onset_report_expected_total`](@ref) as of `v`, over the onsets the model
runs past the cut-off `n`. Its increment over the total at the cut-off is
drawn with the Student-t the scored cells take, at the scale
[`onset_report_scale`](@ref) gives a correction read off two bars with the
fitted read SD `τ`. The increment is drawn on the whole figure, not per
onset date, since the figure's total is the quantity the forecast is
scored on. The split of each
increment into reports of onsets up to the cut-off (`backfill`) and after it
(`future`) is tracked, with the total the triangle should already have
printed by the cut-off.
"""
@model function onset_forecast_model(onsets, state, n, vintages)
    past = onsets[1:n]
    now = _onset_total(onsets, state, n)
    then = [_onset_total(onsets, state, v) for v in vintages]
    then_past = [_onset_total(past, state, v) for v in vintages]
    means = then .- now
    sds = [onset_report_scale(means[j], state.τ, 2) for j in eachindex(vintages)]
    forecast_onset_reports ~ to_submodel(
        onset_increments_model(means, sds, missing, state.ν)
    )
    forecast_onset_reports_to_date := now
    forecast_onset_reports_backfill := then_past .- now
    forecast_onset_reports_future := then .- then_past
    return (; means, sds)
end

## Cumulative confirmed-case trajectory on the observed scale, shared by the
## joint and the confirmed-only composer. The first confirmed vintage is the
## initial condition and is not scored, so the reconstruction counts only the
## fitted increments. Adding that first count back from the testing onset
## makes the trajectory comparable to the observed total.
function _cumulative_confirmed(confirmed_daily, confirmed_history, n::Integer)
    base = isempty(confirmed_history.counts) ? 0 :
        Int(confirmed_history.counts[1])
    cap = isempty(confirmed_history.days) ? 1 :
        clamp(Int(confirmed_history.days[1]), 1, n)
    return cumsum(confirmed_daily) .+
        [t >= cap ? base : 0 for t in eachindex(confirmed_daily)]
end

## Every composer exposes its stream's cumulative trajectory under the same
## un-prefixed `:=` name, so one reader serves the joint and the
## single-stream fits.

"""
Exports-only composer (geographic-spread analogue). Runs the infection
process and onset staging, samples ascertainment, then conditions on the
exports likelihood only. See [`exports_model`](@ref).
"""
@model function exports_only_model(
        n::Integer, exported_cases::Union{Missing, Integer};
        export_case_days::AbstractVector{<:Integer} = Int[],
        breakpoint::Union{Missing, Real} = missing,
        source_population::Real = ITURI_POPULATION,
        infection = infection_model,
        onset_incidence = onset_incidence_model,
        exports = exports_model,
        ascertainment = pooled_ascertainment_model()
    )
    latent ~ to_submodel(
        _latent(n, breakpoint, infection, onset_incidence), false
    )
    asc_state ~ to_submodel(ascertainment)
    exports_state ~ to_submodel(
        exports(
            exported_cases, latent.infection_state.infections,
            asc_state.p_uganda; export_case_days,
            incubation_pmf = latent.incubation_pmf,
            source_population
        )
    )
end

"""
Deaths-only composer (back-calculation analogue). Runs the infection
process and onset staging, samples dispersion, then conditions on the
deaths likelihood only. See [`deaths_model`](@ref).
"""
@model function deaths_only_model(
        n::Integer, total_deaths::Union{Missing, Integer};
        deaths_history = (; days = Int[], counts = Int[]),
        suspected_daily_deaths_history = (; days = Int[], counts = Int[]),
        breakpoint::Union{Missing, Real} = missing,
        infection = infection_model,
        onset_incidence = onset_incidence_model,
        deaths = deaths_model,
        dispersion = surveillance_dispersion_model(),
        forecast::Union{Nothing, ForecastHorizon} = nothing
    )
    ckw = forecast === nothing ? (;) : (; cutoff = n)
    latent ~ to_submodel(
        _latent(n, breakpoint, infection, onset_incidence; forecast), false
    )
    dispersion_state ~ to_submodel(dispersion)
    deaths_state ~ to_submodel(
        deaths(
            deaths_history, total_deaths, latent.onsets,
            dispersion_state.k; suspected_daily_deaths_history, ckw...
        )
    )
    cumulative_deaths_total := cumsum(deaths_state.deaths_daily)
    if forecast !== nothing
        forecast_deaths ~ to_submodel(
            _forecast_counts(
                deaths_state.deaths_daily, forecast_days(n, forecast),
                dispersion_state.k
            )
        )
    end
end

"""
Cases-only composer (reported-cases ascertainment). Runs the infection
process and onset staging, samples dispersion and pooled ascertainment,
then conditions on the reported-cases likelihood. See
[`reported_cases_model`](@ref).
"""
@model function cases_only_model(
        n::Integer, reported_cases::Union{Missing, Integer};
        reported_history = (; days = Int[], counts = Int[]),
        suspected_daily_history = (; days = Int[], counts = Int[]),
        breakpoint::Union{Missing, Real} = missing,
        infection = infection_model,
        onset_incidence = onset_incidence_model,
        cases = reported_cases_model,
        dispersion = surveillance_dispersion_model(),
        ascertainment = pooled_ascertainment_model(),
        forecast::Union{Nothing, ForecastHorizon} = nothing
    )
    ckw = forecast === nothing ? (;) : (; cutoff = n)
    latent ~ to_submodel(
        _latent(n, breakpoint, infection, onset_incidence; forecast), false
    )
    dispersion_state ~ to_submodel(dispersion)
    asc_state ~ to_submodel(ascertainment)
    cases_state ~ to_submodel(
        cases(
            reported_history, reported_cases, latent.onsets,
            dispersion_state.k, asc_state.p_drc; suspected_daily_history,
            ckw...
        )
    )
    cumulative_reports := cumsum(cases_state.reports_daily)
    if forecast !== nothing
        forecast_reports ~ to_submodel(
            _forecast_counts(
                cases_state.reports_daily, forecast_days(n, forecast),
                dispersion_state.k
            )
        )
    end
end

"""
Confirmed-cases-only composer (laboratory pipeline in isolation). Runs
the infection process and onset staging, samples dispersion and pooled
ascertainment, then runs the suspected-case stream in predictive mode (to
draw the shared background rate, testing fraction and onset-to-report
kernel) and conditions on the laboratory pipeline alone. That is the
confirmed positives, a Binomial of the observed analysed denominator in
`lab_history`, and the modelled analysed-specimen volume. See
[`confirmed_cases_model`](@ref) and [`reported_cases_model`](@ref).

Exposes the cut-off expected confirmed count as `expected_confirmed_T`,
the same un-prefixed name [`bvd_joint`](@ref) uses. With `forecast` it draws
the confirmed counts past the cut-off as `forecast_confirmed`, the future
draws [`forecast_stream`](@ref) reads.
"""
@model function confirmed_only_model(
        n::Integer, confirmed_cases::Union{Missing, Integer};
        confirmed_history = (; days = Int[], counts = Int[]),
        lab_history = (; days = Int[], counts = Int[]),
        lab_daily_history = (; days = Int[], counts = Int[]),
        tests_analysed::Union{Missing, Integer} = missing,
        breakpoint::Union{Missing, Real} = missing,
        confirmed_break_days::AbstractVector{<:Integer} = Int[],
        confirmed_break_gross_cases::AbstractVector{<:Integer} = Int[],
        confirmed_break_sd::Real = 25.0,
        infection = infection_model,
        onset_incidence = onset_incidence_model,
        cases = reported_cases_model,
        confirmed = confirmed_cases_model,
        dispersion = surveillance_dispersion_model(),
        ascertainment = pooled_ascertainment_model(),
        forecast::Union{Nothing, ForecastHorizon} = nothing
    )
    ckw = forecast === nothing ? (;) : (; cutoff = n)
    latent ~ to_submodel(
        _latent(n, breakpoint, infection, onset_incidence; forecast), false
    )
    dispersion_state ~ to_submodel(dispersion)
    asc_state ~ to_submodel(ascertainment)
    k = dispersion_state.k
    p_drc = asc_state.p_drc
    cases_state ~ to_submodel(
        cases(
            (; days = Int[], counts = Int[]), missing, latent.onsets,
            k, p_drc; ckw...
        )
    )

    confirmed_state ~ to_submodel(
        confirmed(
            confirmed_history, confirmed_cases, latent.onsets, k,
            p_drc, cases_state.bg_daily, cases_state.τ_test,
            cases_state.bvd_reports_daily;
            lab_history, lab_daily_history,
            tests_analysed, confirmed_break_days,
            confirmed_break_gross = confirmed_break_gross_cases,
            confirmed_break_sd, ckw...
        )
    )

    expected_confirmed_T := confirmed_state.expected_confirmed
    cumulative_confirmed := _cumulative_confirmed(
        confirmed_state.confirmed_daily, confirmed_history, n
    )
    ## A future day has no published analysed count, so its confirmed
    ## cases are the negative binomial the late windows take.
    if forecast !== nothing
        forecast_confirmed ~ to_submodel(
            _forecast_counts(
                confirmed_state.confirmed_daily,
                forecast_days(n, forecast), k
            )
        )
    end
end

"""
Isolation-occupancy-only composer (treatment-bed prevalence in isolation).
Runs the infection process and onset staging, samples dispersion and pooled
ascertainment, then runs the suspected-case stream in predictive mode (to
draw the shared background rate, testing fraction and onset-to-report
kernel) and conditions on the isolation/treatment-bed occupancy alone. See
[`treatment_flow_model`](@ref) and [`reported_cases_model`](@ref).
"""
@model function treatment_only_model(
        n::Integer;
        isolation_history = (; days = Int[], counts = Int[]),
        bed_capacity_history = (; days = Int[], counts = Int[]),
        treatment_admissions_history = (; days = Int[], counts = Int[]),
        treatment_deaths_history = (; days = Int[], counts = Int[]),
        treatment_ruleout_history = (; days = Int[], counts = Int[]),
        treatment_absconded_history = (; days = Int[], counts = Int[]),
        treatment_confirmed_incare_history = (; days = Int[], counts = Int[]),
        treatment_suspect_incare_history = (; days = Int[], counts = Int[]),
        confirmed_history = (; days = Int[], counts = Int[]),
        confirmed_cases::Union{Missing, Integer} = missing,
        lab_history = (; days = Int[], counts = Int[]),
        lab_daily_history = (; days = Int[], counts = Int[]),
        tests_analysed::Union{Missing, Integer} = missing,
        occupancy_break_days::AbstractVector{<:Integer} = Int[],
        confirmed_break_days::AbstractVector{<:Integer} = Int[],
        confirmed_break_gross_cases::AbstractVector{<:Integer} = Int[],
        confirmed_break_sd::Real = 25.0,
        breakpoint::Union{Missing, Real} = missing,
        infection = infection_model,
        onset_incidence = onset_incidence_model,
        cases = reported_cases_model,
        confirmed = confirmed_cases_model,
        treatment = treatment_flow_model,
        cfr = cfr_model(),
        dispersion = surveillance_dispersion_model(),
        ascertainment = pooled_ascertainment_model(),
        ## Built once, with the model, and passed to `treatment`.
        treatment_defaults = treatment_flow_defaults(),
        forecast::Union{Nothing, ForecastHorizon} = nothing
    )
    ckw = forecast === nothing ? (;) : (; cutoff = n)
    latent ~ to_submodel(
        _latent(n, breakpoint, infection, onset_incidence; forecast), false
    )
    dispersion_state ~ to_submodel(dispersion)
    asc_state ~ to_submodel(ascertainment)
    cfr_state ~ to_submodel(cfr)
    k = dispersion_state.k
    p_drc = asc_state.p_drc
    cases_state ~ to_submodel(
        cases(
            (; days = Int[], counts = Int[]), missing, latent.onsets,
            k, p_drc; ckw...
        )
    )
    ## Confirmed-case lab pipeline, run so the treatment model can borrow the
    ## daily testing intensity and positivity for the in-care confirmation
    ## overlay.
    confirmed_state ~ to_submodel(
        confirmed(
            confirmed_history, confirmed_cases, latent.onsets, k,
            p_drc, cases_state.bg_daily, cases_state.τ_test,
            cases_state.bvd_reports_daily;
            lab_history, lab_daily_history, tests_analysed,
            confirmed_break_days,
            confirmed_break_gross = confirmed_break_gross_cases,
            confirmed_break_sd, ckw...
        )
    )
    ## In-care confirmation hazard `τ_test · p_pos` on the daily grid.
    conf_hazard_daily = confirmed_state.τ_test .* confirmed_state.p_pos_grid
    treatment_state ~ to_submodel(
        treatment(
            isolation_history, cases_state.bvd_reports_daily,
            cases_state.bg_daily, p_drc, cfr_state.CFR;
            capacity_history = bed_capacity_history,
            admissions_history = treatment_admissions_history,
            deaths_history = treatment_deaths_history,
            ruleout_history = treatment_ruleout_history,
            absconded_history = treatment_absconded_history,
            confirmed_incare_history = treatment_confirmed_incare_history,
            suspect_incare_history = treatment_suspect_incare_history,
            occupancy_break_days = occupancy_break_days,
            conf_hazard_daily = conf_hazard_daily,
            defaults = treatment_defaults, ckw...
        )
    )
    if forecast !== nothing
        treatment_forecast ~ to_submodel(
            treatment_forecast_model(
                treatment_state, forecast_days(n, forecast),
                bed_capacity_history, isolation_history,
                treatment_state.k_isolation
            ), false
        )
    end
end

"""
Onsets-only composer (the direct-observation analogue). Runs the infection
process and onset staging, then conditions on the symptom-onset reporting-
triangle likelihood alone. See [`onset_reporting_model`](@ref) for the
delay hazard, calendar-time drift, right-truncation and ascertainment
maths.

This stream needs no dispersion submodel of its own. `onset_report`
samples its own ascertainment level and is the only injected submodel. No
confirmed pipeline is available to anchor ascertainment on, so it falls
back to its constant `0.15` anchor.

Exposes the cut-off expected onset-reported count as
`expected_onset_reported_T`, the un-prefixed name [`bvd_joint`](@ref) uses,
and the modelled ascertainment as `onset_ascertainment`. With `forecast` it
draws future vintages of the triangle ([`onset_forecast_model`](@ref)).
Those future draws and the shared `cumulative_onsets` trajectory from
`_latent` are what [`forecast_onsets`](@ref) reads, so this fit nowcasts and
forecasts the onset stream, scored on the reported increment rather than on
the digitised level.
"""
@model function onsets_only_model(
        n::Integer;
        onset_curve_history = (;
            onset_days = Int[], report_days = Int[],
            prev_report_days = Int[], increments = Int[],
        ),
        breakpoint::Union{Missing, Real} = missing,
        infection = infection_model,
        onset_incidence = onset_incidence_model,
        onset_report = onset_reporting_model,
        forecast::Union{Nothing, ForecastHorizon} = nothing
    )
    latent ~ to_submodel(
        _latent(n, breakpoint, infection, onset_incidence; forecast), false
    )
    onset_report_state ~ to_submodel(
        onset_report(onset_curve_history, latent.onsets)
    )
    ## Reported only, so built only when `:=` values are recorded.
    if _reporting(__varinfo__)
        expected_onset_reported_T := onset_report_expected_total(
            latent.onsets,
            onset_report_state.logit_h0, onset_report_state.γ,
            onset_report_state.grid_start, onset_report_state.alpha, n
        )
    end
    onset_ascertainment := onset_report_state.alpha
    if forecast !== nothing && !isempty(onset_curve_history.onset_days)
        onset_forecast ~ to_submodel(
            onset_forecast_model(
                latent.onsets, onset_report_state, n,
                future_knot_days(n, horizon_days(forecast))
            ), false
        )
    end
    return (; onsets = latent.onsets, onset_report_state)
end

"""
Confirmed-deaths-only composer. Runs the infection process and onset
staging, samples dispersion and pooled ascertainment, runs the reported-
cases stream (in predictive mode, to supply the non-BVD background the death
background is scaled from) and the suspected-deaths stream, then conditions
on the confirmed-death likelihood alone. See
[`confirmed_deaths_model`](@ref).
"""
@model function confirmed_deaths_only_model(
        n::Integer, confirmed_deaths::Union{Missing, Integer},
        total_deaths::Union{Missing, Integer} = missing;
        deaths_history = (; days = Int[], counts = Int[]),
        confirmed_deaths_history = (; days = Int[], counts = Int[]),
        confirmed_break_days::AbstractVector{<:Integer} = Int[],
        confirmed_break_gross_deaths::AbstractVector{<:Integer} = Int[],
        confirmed_break_sd::Real = 25.0,
        breakpoint::Union{Missing, Real} = missing,
        infection = infection_model,
        onset_incidence = onset_incidence_model,
        deaths = deaths_model,
        cases = reported_cases_model,
        confirmed_deaths_stream = confirmed_deaths_model,
        dispersion = surveillance_dispersion_model(),
        ascertainment = pooled_ascertainment_model(),
        forecast::Union{Nothing, ForecastHorizon} = nothing
    )
    ckw = forecast === nothing ? (;) : (; cutoff = n)
    latent ~ to_submodel(
        _latent(n, breakpoint, infection, onset_incidence; forecast), false
    )
    dispersion_state ~ to_submodel(dispersion)
    asc_state ~ to_submodel(ascertainment)
    k = dispersion_state.k
    p_drc = asc_state.p_drc
    cases_state ~ to_submodel(
        cases(
            (; days = Int[], counts = Int[]), missing, latent.onsets,
            k, p_drc; ckw...
        )
    )
    deaths_state ~ to_submodel(
        deaths(
            deaths_history, total_deaths, latent.onsets, k;
            case_bg_daily = cases_state.bg_daily, ckw...
        )
    )
    confirmed_deaths_state ~ to_submodel(
        confirmed_deaths_stream(
            confirmed_deaths, total_deaths,
            deaths_state.deaths_daily, deaths_state.bvd_deaths_daily,
            deaths_state.bg_death_daily, k;
            confirmed_deaths_history, confirmed_break_days,
            confirmed_break_gross = confirmed_break_gross_deaths,
            confirmed_break_sd, ckw...
        )
    )
    cumulative_confirmed_deaths := cumsum(
        confirmed_deaths_state.confirmed_death_daily
    )
    if forecast !== nothing
        forecast_confirmed_deaths ~ to_submodel(
            _forecast_counts(
                confirmed_deaths_state.confirmed_death_daily,
                forecast_days(n, forecast), k
            )
        )
    end
end

"""
Deaths-among-exports-only composer. Runs the infection process and onset
staging, samples ascertainment and the deaths submodel (for the CFR and
onset-to-death delay), then conditions on the export-deaths likelihood.
See [`exports_deaths_model`](@ref).
"""
@model function exports_deaths_only_model(
        n::Integer, exports_deaths::Union{Missing, Integer};
        export_death_days::AbstractVector{<:Integer} = Int[],
        breakpoint::Union{Missing, Real} = missing,
        source_population::Real = ITURI_POPULATION,
        infection = infection_model,
        onset_incidence = onset_incidence_model,
        deaths = deaths_model,
        exports = exports_model,
        dispersion = surveillance_dispersion_model(),
        ascertainment = pooled_ascertainment_model()
    )
    latent ~ to_submodel(
        _latent(n, breakpoint, infection, onset_incidence), false
    )
    dispersion_state ~ to_submodel(dispersion)
    asc_state ~ to_submodel(ascertainment)
    deaths_state ~ to_submodel(
        deaths(
            (; days = Int[], counts = Int[]), missing, latent.onsets,
            dispersion_state.k
        )
    )
    exports_state ~ to_submodel(
        exports(
            missing, latent.infection_state.infections,
            asc_state.p_uganda; incubation_pmf = latent.incubation_pmf,
            source_population
        )
    )
    exports_deaths_state ~ to_submodel(
        exports_deaths_model(
            exports_deaths,
            exports_state.travelled_prevalence, deaths_state.CFR,
            deaths_state.od_pmf, latent.incubation_pmf; export_death_days
        )
    )
end

"""
Exports-joint composer. The Uganda export cases and deaths fit together as
one geographic-spread stream. Runs the infection process and onset staging,
samples ascertainment and the deaths submodel (for the CFR and
onset-to-death delay), then conditions on both the export-case and
export-death likelihoods over the one travel-gated at-risk prevalence, so
the two inform the outbreak size jointly. Either count may be `missing` to
drop it. See [`exports_model`](@ref) and [`exports_deaths_model`](@ref).
"""
@model function exports_joint_only_model(
        n::Integer, exported_cases::Union{Missing, Integer},
        exports_deaths::Union{Missing, Integer};
        export_case_days::AbstractVector{<:Integer} = Int[],
        export_death_days::AbstractVector{<:Integer} = Int[],
        breakpoint::Union{Missing, Real} = missing,
        source_population::Real = ITURI_POPULATION,
        infection = infection_model,
        onset_incidence = onset_incidence_model,
        deaths = deaths_model,
        exports = exports_model,
        dispersion = surveillance_dispersion_model(),
        ascertainment = pooled_ascertainment_model(),
        forecast::Union{Nothing, ForecastHorizon} = nothing
    )
    ckw = forecast === nothing ? (;) : (; cutoff = n)
    latent ~ to_submodel(
        _latent(n, breakpoint, infection, onset_incidence; forecast), false
    )
    dispersion_state ~ to_submodel(dispersion)
    asc_state ~ to_submodel(ascertainment)
    deaths_state ~ to_submodel(
        deaths(
            (; days = Int[], counts = Int[]), missing, latent.onsets,
            dispersion_state.k; ckw...
        )
    )
    exports_state ~ to_submodel(
        exports(
            exported_cases, latent.infection_state.infections,
            asc_state.p_uganda; export_case_days,
            incubation_pmf = latent.incubation_pmf, source_population,
            ckw...
        )
    )
    exports_deaths_state ~ to_submodel(
        exports_deaths_model(
            exports_deaths,
            exports_state.travelled_prevalence, deaths_state.CFR,
            deaths_state.od_pmf, latent.incubation_pmf; export_death_days,
            ckw...
        )
    )
    if forecast !== nothing
        forecast_exports ~ to_submodel(
            _forecast_exports(exports_state, forecast_days(n, forecast))
        )
    end
end

## --- Patch (meta-population) latent process ----------------------------

## Run the patch renewal and expose the per-patch state alongside the
## national aggregates the observation submodels consume. Mirrors
## [`_latent`](@ref) and is attached unprefixed, so the shared
## latent-trajectory deterministics surface bare under the same names a
## single-population chain carries.
@model function _patch_latent(
        n::Integer, n_patches::Integer,
        breakpoint, patch_infection;
        rt_start::Integer = 1,
        rt_walk_start::Integer = rt_start,
        importation_kernel::AbstractMatrix = province_importation_kernel(
            PROVINCE_POPULATIONS[1:min(n_patches, end)]
        ),
        forecast::Union{Nothing, ForecastHorizon} = nothing
    )
    fkw = forecast === nothing ? (;) : (; forecast)
    patch_state ~ to_submodel(
        patch_infection(
            n, n_patches;
            breakpoint, rt_start, rt_walk_start,
            importation_kernel, fkw...
        ), false
    )
    ## Summed over the patches with one matrix-vector product.
    onsets_total = transpose(patch_state.onsets_matrix) *
        ones(eltype(patch_state.onsets_matrix), n_patches)
    cumulative_infections := patch_state.cumulative_total
    C_T := patch_state.C_T
    if _reporting(__varinfo__)
        cumulative_onsets := cumsum(onsets_total)
        susceptible_fraction := pool_fraction(
            view(patch_state.cumulative_total, 1:n),
            sum(patch_state.populations)
        )
    end
    ## Past the cut-off the national reproduction number is read off the
    ## summed infections, as `R_T` is at the cut-off.
    if forecast !== nothing
        fd = forecast_days(n, forecast)
        forecast_infections := patch_state.infections_total[fd]
        forecast_onsets := onsets_total[fd]
        forecast_rt := _implied_rt_path(
            patch_state.infections_total, patch_state.g, fd
        )
    end
    return (; patch_state, onsets_total)
end

## Each patch's share of its population left susceptible each day, one row
## per patch ([`pool_fraction`](@ref)). Only `:=` quantities read it.
function _patch_fractions(state)
    I = state.infections_matrix
    return reduce(
        vcat,
        [
            pool_fraction(
                cumsum(view(I, p, :)), state.populations[p]
            )' for p in axes(I, 1)
        ]
    )
end

## Patch `p`'s reproduction number net of its own depletion on `days`, from
## the fractions `fr` of `_patch_fractions`.
_patch_adjusted_rt(state, fr, p, days) = adjusted_rt(
    view(state.Rt_matrix, p, :), view(fr, p, :), days
)

## The implied national reproduction number on each forecast day. A function
## rather than a comprehension in the model body, which would capture
## model-scope locals and box them.
_implied_rt_path(infections, g, days) =
    [implied_national_Rt_at(infections, g, t) for t in days]

"""
Modelled per-province confirmed-death increments, binned to the vintages
of the Tableau 1 spatial tables. Each patch's onsets are pushed through the
onset-to-death delay and then the report-to-receipt delay, so `kernel` is
the convolution of the two, giving the deaths confirmed by each vintage.

This is the term that makes the provincial split identifiable. The
case-fatality ratio and the death-confirmation probability are properties
of the virus and of a national laboratory pipeline, not of a province, so
they cancel out of the normalised shares. What remains weights each patch
by its delay-convolved incidence alone, free of case ascertainment.

The delay convolution also separates ascertainment from epidemic phase. A
fast-growing province has proportionally fewer deaths to date than a flat
one at the same true CFR, because its recent cases have not died yet.
Predicting each province's deaths to date from its own incidence curve
accounts for that censoring, so only the residual is ascertainment.
"""
function _patch_death_increments(
        onsets_matrix::AbstractMatrix,
        kernel::AbstractVector,
        province_days::AbstractVector{<:Integer}
    )
    np = size(onsets_matrix, 1)
    nv = length(province_days)
    first_daily = convolve_delay(onsets_matrix[1, :], kernel)
    out = Matrix{eltype(first_daily)}(undef, np, nv)
    @inbounds out[1, :] = bin_increments(first_daily, province_days)
    @inbounds for p in 2:np
        daily = convolve_delay(onsets_matrix[p, :], kernel)
        out[p, :] = bin_increments(daily, province_days)
    end
    return out
end

"""
Modelled per-province confirmed increments, binned to the vintages of the
per-province spatial tables. Each patch's onsets are pushed through the
same onset-to-confirmation kernel and test sensitivity as the national
confirmed stream (`kernel` is the onset-to-report pmf convolved with the
report-to-receipt pmf; the laboratory pipeline is national, only the
incidence feeding it is provincial) and then binned onto the shared
vintage days.

`province_days` are the (shared) grid-day indices of the spatial-table
vintages. Returns an `(n_patches × n_vintages)` matrix.

!!! note "Ascertainment and incidence are confounded by province"
    This matrix reaches [`province_composition_model`](@ref) only through
    normalised shares, so every factor common to all provinces cancels
    (`s_test` here, and ascertainment, background and positivity in the
    national confirmed stream). That keeps the composition from re-scoring
    the national total.

    It does not fix the case-finding probability across provinces. The
    composition weights each patch by `asc_p * lambda_p`, its relative
    ascertainment times its modelled incidence, and the data identify only
    that product. The per-province laboratory series shows very
    differently-selected testing pools (Ituri 31.8% positivity against
    Nord-Kivu's 5.5%), so `asc_p` is sampled and partially pooled toward
    equality rather than fixed, which widens the per-patch `Rt` contrast
    and `C_T` split to their honest width.

    Read `log_rt_contrast` and `province_ascertainment` together. Neither
    is interpretable alone. The national headline does not depend on the
    split.
"""
function _patch_confirmed_increments(
        onsets_matrix::AbstractMatrix,
        kernel::AbstractVector, s_test::Real,
        province_days::AbstractVector{<:Integer}
    )
    return _patch_confirmed_increments(
        _patch_carried(onsets_matrix, kernel), s_test, province_days
    )
end

## The same binning from an already carried `(n_patches × n)` matrix.
function _patch_confirmed_increments(
        carried::AbstractMatrix, s_test::Real,
        province_days::AbstractVector{<:Integer}
    )
    np = size(carried, 1)
    nv = length(province_days)
    T = promote_type(eltype(carried), typeof(float(s_test)))
    out = Matrix{T}(undef, np, nv)
    @inbounds for p in 1:np
        out[p, :] = bin_increments(
            s_test .* carried[p, :], province_days
        )
    end
    return out
end

"""
Each patch's onsets carried through a delay kernel, an `(n_patches × n)`
matrix. The confirmed composition and the laboratory composition both carry
the onsets through the onset-to-confirmation kernel, so [`bvd_joint`](@ref)
carries them once and hands the matrix to both.
"""
function _patch_carried(onsets_matrix::AbstractMatrix, kernel::AbstractVector)
    np, n = size(onsets_matrix)
    T = promote_type(eltype(onsets_matrix), eltype(kernel))
    out = Matrix{T}(undef, np, n)
    @inbounds for p in 1:np
        out[p, :] = convolve_delay(onsets_matrix[p, :], kernel)
    end
    return out
end

## Province rows are `nothing` or a `(; days, patches, counts)` tuple.
_has_province_rows(rows) = rows !== nothing && !isempty(rows.days)

"""
Per-patch BVD reports, each patch's onsets through the shared
onset-to-report delay, the rows summing to the national `bvd_reports_daily`
of [`reported_cases_model`](@ref). Returns an `(n_patches × n)` matrix.
"""
function _patch_reports(onsets_matrix::AbstractMatrix, report_pmf::AbstractVector)
    np, n = size(onsets_matrix)
    T = promote_type(eltype(onsets_matrix), eltype(report_pmf))
    out = Matrix{T}(undef, np, n)
    @inbounds for p in 1:np
        out[p, :] = convolve_delay(onsets_matrix[p, :], report_pmf)
    end
    return out
end

"""
Modelled per-patch analysed-specimen volume summed over the laboratory
bins (`lab_bins[i]` is the bin of printed day `lab_days[i]`), up to
the national factors the composition normalises away. `carried` is each
patch's onsets through the onset-to-confirmation kernel (report ⊕ receipt,
[`_patch_carried`](@ref)); the national BVD volume `p_drc` times the
patch sum is split by ascertainment-weighted incidence (`asc`, the case
composition's relative ascertainment), so the BVD-to-background proportion
stays the national one, and joined by the patch's share `w[p]` of the
non-BVD
background carried to receipt (`bg_carried`). Summed over the patches this
is the national suspect pipeline [`confirmed_cases_model`](@ref) scales
into the analysed volume, so the shares are exact. `lab_days` are grid day
indices. Returns an `(n_patches × n_days)` matrix.
"""
function _patch_analysed_increments(
        carried::AbstractMatrix,
        p_drc::Real, asc::AbstractVector, bg_carried::AbstractVector,
        w::AbstractVector, lab_days::AbstractVector{<:Integer},
        lab_bins::AbstractVector{<:Integer} = 1:length(lab_days)
    )
    np, n = size(carried)
    nb = isempty(lab_bins) ? 0 : maximum(lab_bins)
    T = promote_type(
        eltype(carried), typeof(float(p_drc)), eltype(asc), eltype(bg_carried),
        eltype(w)
    )
    out = zeros(T, np, nb)
    @inbounds for (i, day) in enumerate(lab_days)
        d = clamp(Int(day), 1, n)
        ## The national BVD volume that day, split by ascertainment-weighted
        ## incidence, so the BVD-to-background proportion stays national.
        total = zero(T)
        weighted = zero(T)
        for q in 1:np
            total += carried[q, d]
            weighted += asc[q] * carried[q, d]
        end
        bvd_nat = p_drc * total
        for p in 1:np
            share = weighted > zero(T) ? asc[p] * carried[p, d] / weighted :
                one(T) / np
            out[p, lab_bins[i]] += bvd_nat * share + w[p] * bg_carried[d]
        end
    end
    return out
end

"""
Joint composer over all data streams. Runs the generating infection
process once on a daily grid of length `n` (day `n` is the cut-off),
stages it to daily onset incidence, then conditions on the DRC suspected
cases, deaths and the laboratory pipeline, the confirmed deaths, the
Uganda exports and deaths-among-exports, the digitised symptom-onset
reporting triangle ([`onset_reporting_model`](@ref), the only direct
observation of the shared onset series, a no-op when `onset_curve_history`
is empty), and the optional genetic seeding bound on the outbreak age.
Each stream argument may be `missing` to drop it, so the model doubles as
a prior- and posterior-predictive generator.

The confirmed-case stream shares its onset-to-report kernel with the
suspected-case stream. A single analysed-specimen volume is fit through a
report-to-analysed delay and the tested fraction. The confirmed positives
are scored as a Binomial of the observed specimens-analysed denominator
(`lab_history`) with a partially-pooled per-window positivity, so they do
not pass through the multiplicative ascertainment ridge. After the national
cumulative analysed series stops, the reporting format gives a 24h analysed
count on some days (`lab_daily_history`). Those are fitted as per-day
analysed volumes and also anchor that day's confirmed positives. Windows
with no published denominator use the modelled analysed volume as the
denominator, with the positivity (hence `λ_bg`) carried over from the
windows that do have data (see [`confirmed_cases_model`](@ref)).

The optional `suspected_daily_history` adds the post-26 May daily
new-suspect inflow ("nouveaux cas suspects du jour"), scored against the
modelled daily suspected series on the report days where the frozen
cumulative suspected stream stops, disjoint from it. The optional
`suspected_daily_deaths_history` adds the deaths analogue, the daily new
suspected deaths ("cas suspects du jour N (M deces)").

The confirmed deaths mirror the confirmed-case laboratory pipeline. A
death "analysed" volume (the suspected deaths carried to laboratory
receipt and scaled by the death testing intensity `tau_death`, specimens
per suspected death, which may exceed one) is scored
through a death-pool composition positivity
`p = s·q_death + (1−spec)(1−q_death)`, with `q_death` the BVD share of
the suspected deaths (see [`confirmed_deaths_model`](@ref)). The
suspected deaths carry a death ascertainment `p_death` and a non-BVD
background tied to the case background by a background CFR `cfr_bg`
(see [`deaths_model`](@ref)).

The optional `isolation_history` adds the daily isolation/treatment-bed
occupancy ("Patients en isolement"), a prevalence stream fitted as the
suspect inflow (BVD treatment stay plus non-BVD rule-out stay) carried
through a length-of-stay survival into a daily stock (see
[`treatment_flow_model`](@ref)). The optional `recovered_history` adds
the recovered-among-confirmed stream ("cumul guéris"), survivors among
the modelled daily confirmed cases scaled by the recovery probability
and lagged by a confirmation-to-recovery delay (see
[`recovered_model`](@ref)).

`breakpoint` is the intervention day passed to the reproduction-number
walk (e.g. the first WHO situation report). `genetic` injects the genetic
seeding submodel when `tmrca_days` is given. The analysed volume is scaled
by specimens analysed per suspect ([`specimen_intensity_model`](@ref)),
which `τ_test` alone, being a probability, cannot exceed one of.

Tracked deterministics: `C_T` (cumulative infections by the cut-off), the
established reproduction number `R0` (= the first `R_t`), `r` and
`doubling_time` (current growth), `r0` (the `R0`-implied cryptic growth
rate), `T` (outbreak age), `R_T` (current reproduction number), the
per-stream expected counts, the testing fraction `tau_test`, the specimens
analysed per suspect (`specimens_per_suspect`), the background rate
`lambda_bg`, the death ascertainment `death_ascertainment`, the background
CFR `background_cfr`, the death testing intensity `tau_death`, the implied
per-suspected (`suspected_positivity`) and per-test (`test_positivity`)
positivities, and the death-pool BVD composition (`death_composition`) and
death-confirmation positivity (`death_confirmation`).

With `n_patches > 1` the latent process is the meta-population renewal of
[`patch_infection_model`](@ref), one renewal equation per province coupled
by importation, with every national stream above fitted against the summed
provinces. The default `n_patches = 1` collapses it onto the
single-population model, since the sum-to-zero deviations vanish, there is
nothing to import between, and no per-province likelihood is scored.

The spatial information enters through two composition terms. The
per-province confirmed cases and confirmed deaths in the situation
reports' spatial tables are exact partitions of the national totals, so
scoring them with their own count likelihoods would put the same data into
the joint density twice. [`province_composition_model`](@ref) scores each
split conditional on the national total instead, adding only the spatial
signal the national series does not carry. Pass the reshaped data as
`province_increments` with `province_days`, and `province_death_increments`
with `province_death_days`, both built by
[`province_increment_matrix`](@ref). Omit them and the composition terms
are skipped. The case composition identifies only the product of a
province's incidence and its case-finding. The death composition
identifies the incidence split, since the case-fatality ratio and the
death-confirmation probability are national, and the case composition then
identifies the relative case ascertainment as the residual.

A third composition scores the per-province analysed-specimen volume
conditional on the national total in each laboratory bin. The modelled
split is each patch's BVD suspects (its onsets through the
onset-to-confirmation kernel, thinned by `p_drc`) plus its share of the
non-BVD background, the share a simplex
([`background_split_model`](@ref)) carries and this term identifies.
The BVD suspects carry the case composition's
relative ascertainment, so the two compositions agree on how many of a
patch's cases reach the laboratory. Pass `province_lab_increments` with
`province_lab_days` and `province_lab_bins`, built by
[`province_lab_increment_matrix`](@ref), which sums the printed days into
calendar weeks for the production fit.
The testing fraction stays national and the composition samples no
ascertainment contrast of its own. The per-province positives are not
fitted, being the differencing of the confirmed counts already scored.

Uganda exports are driven by the provinces in proportion to sampled
relative export weights, with Ituri the reference at weight one (see
[`province_export_pressure_model`](@ref)).

Per-patch quantities are surfaced as vector deterministics, one entry per
patch, for [`patch_summary_table`](@ref) and
[`patch_overview_table`](@ref): `C_T_patch`, `R_T_patch`,
`infections_T_patch`, `delta_patch`, `log_rt_contrast` and the deviation
knots `delta_knots` from which [`reconstruct_patch_rt`](@ref) rebuilds the
provincial trajectories. `R_T` is the force-of-infection-weighted
reproduction number implied by the summed patch infections.

With `forecast` a [`ForecastHorizon`](@ref) (see [`with_horizon`](@ref))
the model runs past the cut-off `n` and draws each stream's future counts
through its own likelihood, as `missing` observations `predict` generates:
`forecast_reports`, `forecast_deaths`, `forecast_confirmed`,
`forecast_confirmed_deaths`, `forecast_recovered`, the treatment flows
([`treatment_forecast_model`](@ref)), `forecast_exports` and the onset
triangle's future vintages ([`onset_forecast_model`](@ref)). With province
data, each future week's national confirmed cases and deaths are split
across the provinces by the fitted compositions
([`composition_split_model`](@ref)), as `forecast_province_confirmed` and
`forecast_province_deaths`, and the national isolation and bed forecasts
likewise as `forecast_province_isolation` and `forecast_province_beds`. Every quantity up to the cut-off, and the
density there, is the fitted model's.
"""
@model function bvd_joint(
        n::Integer,
        exported_cases::Union{Missing, Integer},
        total_deaths::Union{Missing, Integer},
        reported_cases::Union{Missing, Integer} = missing,
        exports_deaths::Union{Missing, Integer} = missing,
        confirmed_cases::Union{Missing, Integer} = missing,
        tests_analysed::Union{Missing, Integer} = missing;
        n_patches::Integer = 1,
        importation_kernel::AbstractMatrix = province_importation_kernel(
            PROVINCE_POPULATIONS[1:min(n_patches, end)]
        ),
        confirmed_deaths::Union{Missing, Integer} = missing,
        recovered_cases::Union{Missing, Integer} = missing,
        deaths_history = (; days = Int[], counts = Int[]),
        reported_history = (; days = Int[], counts = Int[]),
        confirmed_history = (; days = Int[], counts = Int[]),
        confirmed_deaths_history = (; days = Int[], counts = Int[]),
        lab_history = (; days = Int[], counts = Int[]),
        lab_daily_history = (; days = Int[], counts = Int[]),
        suspected_daily_history = (; days = Int[], counts = Int[]),
        suspected_daily_deaths_history = (; days = Int[], counts = Int[]),
        isolation_history = (; days = Int[], counts = Int[]),
        bed_capacity_history = (; days = Int[], counts = Int[]),
        recovered_history = (; days = Int[], counts = Int[]),
        treatment_admissions_history = (; days = Int[], counts = Int[]),
        treatment_deaths_history = (; days = Int[], counts = Int[]),
        treatment_ruleout_history = (; days = Int[], counts = Int[]),
        treatment_absconded_history = (; days = Int[], counts = Int[]),
        treatment_confirmed_incare_history = (; days = Int[], counts = Int[]),
        treatment_suspect_incare_history = (; days = Int[], counts = Int[]),
        occupancy_break_days::AbstractVector{<:Integer} = Int[],
        confirmed_break_days::AbstractVector{<:Integer} = Int[],
        confirmed_break_gross_cases::AbstractVector{<:Integer} = Int[],
        confirmed_break_gross_deaths::AbstractVector{<:Integer} = Int[],
        confirmed_break_sd::Real = 25.0,
        export_case_days::AbstractVector{<:Integer} = Int[],
        export_death_days::AbstractVector{<:Integer} = Int[],
        onset_curve_history = (;
            onset_days = Int[], report_days = Int[],
            prev_report_days = Int[], increments = Int[],
        ),
        breakpoint::Union{Missing, Real} = missing,
        source_population::Real = ITURI_POPULATION,
        patch_infection = patch_infection_model,
        composition = province_composition_model,
        province_increments::Union{
            Missing, AbstractMatrix{<:Integer},
        } = missing,
        province_days::AbstractVector{<:Integer} = Int[],
        province_lab_increments::Union{
            Missing, AbstractMatrix{<:Integer},
        } = missing,
        province_lab_days::AbstractVector{<:Integer} = Int[],
        province_lab_bins::AbstractVector{<:Integer} = 1:length(province_lab_days),
        lab_composition = province_composition_model,
        background_split = background_split_model,
        province_isolation = nothing,
        province_capacity = nothing,
        province_death_increments::Union{
            Missing, AbstractMatrix{<:Integer},
        } = missing,
        province_death_days::AbstractVector{<:Integer} = Int[],
        death_composition = province_composition_model,
        death_ascertainment_sd_prior = truncated(
            Normal(0, 0.1); lower = 0
        ),
        province_cfr_sd_prior = truncated(Normal(0, 0.1); lower = 0),
        export_pressure = province_export_pressure_model,
        exports = exports_model,
        deaths = deaths_model,
        cases = reported_cases_model,
        confirmed = confirmed_cases_model,
        confirmed_deaths_stream = confirmed_deaths_model,
        treatment = treatment_flow_model,
        recovered = recovered_model,
        onset_report = onset_reporting_model,
        dispersion = pooled_dispersion_model,
        ascertainment = pooled_ascertainment_model(),
        background_pooling = nothing,
        genetic = nothing,
        onset_to_sample = nejm_onset_to_sample(),
        tmrca_days::Union{Missing, Real} = missing,
        tmrca_days_sd::Real = 16.0,
        renewal_start_lead::Integer = RENEWAL_START_LEAD,
        rt_walk_lead::Integer = RT_WALK_LEAD,
        ## Days the suspected-case background starts before the first
        ## reported case: the support of the default report-to-receipt
        ## kernel, so the convolution into the analysed volume is fully
        ## formed by that report.
        background_onset_lead::Integer = cdf_nmax(lognormal_meansd(4.5, 4.0)),
        ## Built once, with the model, and passed to `treatment`.
        treatment_defaults = treatment_flow_defaults(),
        forecast::Union{Nothing, ForecastHorizon} = nothing,
        simulated_data = nothing
    )

    if n_patches == 1 && (
            !isempty(province_days) || !isempty(province_death_days) ||
                !isempty(province_lab_days) ||
                _has_province_rows(province_isolation) ||
                _has_province_rows(province_capacity)
        )
        error(
            "per-province data was supplied but n_patches = 1. The " *
                "spatial structure would be silently dropped. Pass " *
                "n_patches = $(length(PROVINCE_NAMES)) (or the number of " *
                "patches the data covers)."
        )
    end

    rt_start = ismissing(tmrca_days) ? 1 :
        clamp(n - round(Int, tmrca_days) + renewal_start_lead, 1, n)

    rt_walk_start = ismissing(breakpoint) ? rt_start :
        clamp(round(Int, breakpoint) - rt_walk_lead, rt_start, n)
    latent ~ to_submodel(
        _patch_latent(
            n, n_patches, breakpoint, patch_infection;
            rt_start, rt_walk_start, importation_kernel, forecast
        ), false
    )
    patch_state = latent.patch_state
    onsets = latent.onsets_total
    ## Each stream reads its cut-off quantities at day `n` whatever length
    ## the latent grid runs to.
    ckw = forecast === nothing ? (;) : (; cutoff = n)

    dispersion_state ~ to_submodel(dispersion(6))
    asc_state ~ to_submodel(ascertainment)
    kv = dispersion_state.k
    k_cases = kv[1]
    k_deaths = kv[2]
    k_confirmed = kv[3]
    k_confirmed_deaths = kv[4]
    k_isolation = kv[5]
    k_recovered = kv[6]
    p_drc = asc_state.p_drc
    p_uganda = asc_state.p_uganda

    bg_onset = isempty(reported_history.days) ? 1 :
        clamp(Int(reported_history.days[1]) - background_onset_lead, 1, n)

    ## `nothing` holds the non-BVD background at the constant rate the
    ## testing submodel samples. An injected pooling submodel gives it a
    ## smooth daily random walk instead, whose scale is partially pooled.
    ## Injected rather than switched on a flag so the unused arm is a
    ## `Nothing` the compiler folds away, not a second branch: a `Bool`
    ## reaches the model as a value, so both arms are inferred and the
    ## resulting `Union`-typed argument specialises the whole
    ## suspected-case submodel twice.
    case_bg_re = if background_pooling === nothing
        nothing
    else
        bg_pool ~ to_submodel(background_pooling())
        σ_rw_shared = bg_pool.σ_bg
        (nn; kw...) -> background_walk_model(
            nn, σ_rw_shared; onset = bg_onset, kw...
        )
    end

    ## Cases first so the suspected-case background `bg_daily` is available
    ## to the deaths stream, which scales it by `cfr_bg`, and to the
    ## laboratory pipeline.
    cases_state ~ to_submodel(
        cases(
            reported_history, reported_cases, onsets, k_cases, p_drc;
            suspected_daily_history, background_re = case_bg_re, ckw...,
            _sim_kw(simulated_data, :cases_state)...
        )
    )
    deaths_state ~ to_submodel(
        deaths(
            deaths_history, total_deaths, onsets, k_deaths;
            suspected_daily_deaths_history,
            case_bg_daily = cases_state.bg_daily, ckw...,
            _sim_kw(simulated_data, :deaths_state)...
        )
    )
    confirmed_state ~ to_submodel(
        confirmed(
            confirmed_history, confirmed_cases, onsets, k_confirmed,
            p_drc, cases_state.bg_daily, cases_state.τ_test,
            cases_state.bvd_reports_daily;
            lab_history, lab_daily_history,
            tests_analysed, confirmed_break_days,
            confirmed_break_gross = confirmed_break_gross_cases,
            confirmed_break_sd,
            specimen_intensity = specimen_intensity_model(), ckw...,
            _sim_kw(simulated_data, :confirmed_state)...
        )
    )

    ## Split of the non-BVD suspected background across the patches, shared
    ## by the laboratory composition below and the isolation stream. One
    ## patch takes the whole background and samples nothing.
    bg_split_state ~ to_submodel(background_split(n_patches))
    province_background_split := bg_split_state.w

    ## The anchor stops at the cut-off, so a longer grid leaves the fitted
    ## ascertainment where it was.
    onset_anchor_daily = p_drc .* confirmed_state.τ_test .*
        upto(confirmed_state.p_pos_grid, n)
    onset_report_state ~ to_submodel(
        onset_report(
            onset_curve_history, onsets;
            anchor = onset_anchor_daily
        )
    )

    confirmed_deaths_state ~ to_submodel(
        confirmed_deaths_stream(
            confirmed_deaths, total_deaths,
            deaths_state.deaths_daily, deaths_state.bvd_deaths_daily,
            deaths_state.bg_death_daily, k_confirmed_deaths;
            confirmed_deaths_history, receipt_pmf = confirmed_state.receipt_pmf,
            confirmed_break_days,
            confirmed_break_gross = confirmed_break_gross_deaths,
            confirmed_break_sd,
            case_analysed_daily = confirmed_state.analysed_daily,
            case_suspected_daily = cases_state.reports_daily, ckw...,
            _sim_kw(simulated_data, :confirmed_deaths_state)...
        )
    )

    ## The confirmed and laboratory compositions both carry the patch
    ## onsets through the onset-to-confirmation kernel, once here.
    if !isempty(province_days) || !isempty(province_lab_days)
        confirmed_kernel = convolve_pmf(
            cases_state.report_pmf, confirmed_state.receipt_pmf
        )
        confirmed_carried = _patch_carried(
            patch_state.onsets_matrix, confirmed_kernel
        )
    end

    if !isempty(province_days)
        modelled_prov = _patch_confirmed_increments(
            confirmed_carried, confirmed_state.s_test, province_days
        )
        composition_state ~ to_submodel(
            composition(province_increments, modelled_prov)
        )
        province_shares := composition_state.shares
        province_composition_rho := composition_state.rho
        ## Relative province case ascertainment, the probability an
        ## infection there becomes a confirmed case, partially pooled and
        ## sum-to-zero on the log scale. On its own the case composition
        ## identifies only the product of ascertainment and incidence. The
        ## death composition below separates them.
        province_ascertainment := composition_state.province_ascertainment
        province_ascertainment_sd := composition_state.ascertainment_sd
    end

    ## Relative case ascertainment by patch, shared by every stream that is
    ## driven by reported suspects: the laboratory composition and the
    ## per-patch bed demand split the national BVD volume by
    ## ascertainment-weighted incidence. Ones when the case composition is
    ## not scored.
    patch_asc = isempty(province_days) ? ones(n_patches) :
        composition_state.province_ascertainment

    conf_hazard_daily = confirmed_state.τ_test .* confirmed_state.p_pos_grid
    ## Per-patch BVD reports for the province occupancy split, the rows
    ## summing to the national series the flows are built on.
    bvd_reports_matrix = n_patches > 1 ?
        _patch_reports(patch_state.onsets_matrix, cases_state.report_pmf) :
        nothing
    treatment_state ~ to_submodel(
        treatment(
            isolation_history, cases_state.bvd_reports_daily,
            cases_state.bg_daily, p_drc, deaths_state.CFR;
            bvd_reports_matrix,
            background_split = bg_split_state.w,
            patch_ascertainment = patch_asc,
            province_isolation, province_capacity,
            capacity_history = bed_capacity_history,
            admissions_history = treatment_admissions_history,
            deaths_history = treatment_deaths_history,
            ruleout_history = treatment_ruleout_history,
            absconded_history = treatment_absconded_history,
            confirmed_incare_history = treatment_confirmed_incare_history,
            suspect_incare_history = treatment_suspect_incare_history,
            occupancy_break_days = occupancy_break_days,
            conf_hazard_daily = conf_hazard_daily,
            k_external = k_isolation,
            defaults = treatment_defaults, ckw...,
            _sim_kw(simulated_data, :treatment_state)...
        )
    )

    recovered_state ~ to_submodel(
        recovered(
            recovered_history, recovered_cases,
            confirmed_state.confirmed_daily, deaths_state.CFR;
            k_external = k_recovered, ckw...,
            _sim_kw(simulated_data, :recovered_state)...
        )
    )

    export_pressure_state ~ to_submodel(export_pressure(n_patches))
    export_weight := export_pressure_state.weights
    export_pressure_sd := export_pressure_state.pooling_sd
    ## The patches' infections weighted by export propensity, one
    ## matrix-vector product, over the grid past the cut-off when forecasting.
    export_infections = transpose(patch_state.infections_matrix) *
        export_pressure_state.weights
    exports_state ~ to_submodel(
        exports(
            exported_cases, export_infections, p_uganda;
            export_case_days, incubation_pmf = patch_state.incubation_pmf,
            source_population, ckw...,
            _sim_kw(simulated_data, :exports_state)...
        )
    )
    exports_deaths_state ~ to_submodel(
        exports_deaths_model(
            exports_deaths,
            exports_state.travelled_prevalence, deaths_state.CFR,
            deaths_state.od_pmf, patch_state.incubation_pmf; export_death_days,
            ckw...,
            _sim_kw(simulated_data, :exports_deaths_state)...
        )
    )

    if genetic !== nothing
        genetic_state ~ to_submodel(
            genetic(patch_state.T, tmrca_days; tmrca_days_sd), false
        )
    end


    if !isempty(province_lab_days)
        modelled_lab = _patch_analysed_increments(
            confirmed_carried, p_drc, patch_asc,
            convolve_delay(cases_state.bg_daily, confirmed_state.receipt_pmf),
            bg_split_state.w, province_lab_days, province_lab_bins
        )
        ## No ascertainment contrast of its own: the split is carried by the
        ## background shares and the case composition's ascertainment, and
        ## the testing fraction is national.
        lab_composition_state ~ to_submodel(
            lab_composition(
                province_lab_increments, modelled_lab;
                ascertainment_sd_prior = nothing
            )
        )
        province_lab_shares := lab_composition_state.shares
        province_lab_composition_rho := lab_composition_state.rho
    end

    if !isempty(province_death_days)
        death_kernel = convolve_pmf(
            deaths_state.od_pmf, confirmed_state.receipt_pmf
        )
        modelled_deaths_prov = _patch_death_increments(
            patch_state.onsets_matrix, death_kernel, province_death_days
        )
        death_composition_state ~ to_submodel(
            death_composition(
                province_death_increments,
                modelled_deaths_prov;
                ascertainment_sd_prior = death_ascertainment_sd_prior,
                severity_sd_prior = province_cfr_sd_prior
            )
        )
        province_death_shares := death_composition_state.shares
        province_death_composition_rho := death_composition_state.rho
        province_death_ascertainment := death_composition_state.province_ascertainment
        province_death_ascertainment_sd := death_composition_state.ascertainment_sd
        ## Per-province case-fatality ratio, the national ratio times that
        ## province's sum-to-zero contrast, so the provinces are reported on
        ## the same scale as the national quantity they pool toward.
        province_cfr_relative := death_composition_state.province_severity
        CFR_patch := deaths_state.CFR .*
            death_composition_state.province_severity
        province_cfr_sd := death_composition_state.severity_sd
    end

    ## The cumulative series and the combined delay PMFs below are reported
    ## only, so they are built only when `:=` values are recorded.
    if _reporting(__varinfo__)
        cumulative_expected_deaths := cumsum(deaths_state.bvd_deaths_daily)
        cumulative_confirmed := _cumulative_confirmed(
            confirmed_state.confirmed_daily, confirmed_history, n
        )
        ## Each of the remaining count streams sums to its own cut-off
        ## expected total, so none needs the baseline re-add the confirmed
        ## path takes.
        cumulative_reports := cumsum(cases_state.reports_daily)
        cumulative_deaths_total := cumsum(deaths_state.deaths_daily)
        cumulative_confirmed_deaths := cumsum(
            confirmed_deaths_state.confirmed_death_daily
        )
        cumulative_recovered := cumsum(recovered_state.recovered_daily)
        onset_to_confirmation_pmf := convolve_pmf(
            cases_state.report_pmf, confirmed_state.receipt_pmf
        )
        onset_to_death_confirmation_pmf := convolve_pmf(
            deaths_state.od_pmf, confirmed_state.receipt_pmf
        )
    end

    onset_to_sample_mean := cases_state.report_mean +
        confirmed_state.receipt_mean
    onset_to_sample_sd := sqrt(
        cases_state.report_sd^2 +
            confirmed_state.receipt_sd^2
    )
    if onset_to_sample !== nothing
        @addlogprob! onset_to_sample_logweight(
            cases_state.report_mean,
            cases_state.report_sd, confirmed_state.receipt_mean,
            confirmed_state.receipt_sd, onset_to_sample
        )
    end
    R0 := patch_state.R0
    r := patch_state.r
    r0 := patch_state.r0
    doubling_time := patch_state.doubling_time
    T := patch_state.T

    R_T := patch_state.R_T
    expected_infections_T := @inbounds(patch_state.infections_total[n])
    CFR := deaths_state.CFR
    ## Per-patch quantities, as vector deterministics (one entry per patch).
    if n_patches > 1 && size(treatment_state.capacity_patch, 1) == n_patches
        ## Cut-off beds, demand and censored occupancy by patch, from the
        ## province splits of the isolation stream. Present only when a
        ## province split scored them.
        province_bed_capacity := treatment_state.capacity_patch[:, n]
        province_bed_demand := treatment_state.demand_patch[:, n]
        province_expected_isolation := min.(
            treatment_state.demand_patch[:, n],
            treatment_state.capacity_patch[:, n]
        )
        province_bed_utilisation := min.(
            treatment_state.demand_patch[:, n],
            treatment_state.capacity_patch[:, n]
        ) ./ treatment_state.capacity_patch[:, n]
        province_bed_shortfall := max.(
            treatment_state.demand_patch[:, n] .-
                treatment_state.capacity_patch[:, n],
            0.0
        )
        province_capacity_share := treatment_state.capacity_shares
        ## Daily share of the national bed demand by patch, the modelled
        ## centre of the occupancy split.
        province_occupancy_share := treatment_state.demand_patch ./
            sum(treatment_state.demand_patch; dims = 1)
        province_occupancy_split_rho := treatment_state.occupancy_split_rho
        province_capacity_split_rho := treatment_state.capacity_split_rho
    end
    C_T_patch := patch_state.C_T_patch
    ## Each province's reproduction number net of its own depletion, built
    ## only when recorded.
    if _reporting(__varinfo__)
        fr = _patch_fractions(patch_state)
        susceptible_fraction_patch := vec(fr[:, 1:n])
        R_T_patch := [
            only(_patch_adjusted_rt(patch_state, fr, p, n:n))
                for p in 1:n_patches
        ]
    end
    infections_T_patch := [
        @inbounds(patch_state.infections_matrix[p, n])
            for p in 1:n_patches
    ]

    infections_patch := vec(patch_state.infections_matrix)
    importation_patch := vec(patch_state.importation_matrix)

    delta_patch := [@inbounds(patch_state.δ_patch[p, n]) for p in 1:n_patches]
    delta_patch_start := [
        @inbounds(patch_state.δ_patch[p, rt_walk_start])
            for p in 1:n_patches
    ]
    ## The deviation at every weekly knot, flattened column-major from the
    ## `(n_patches × n_knots)` matrix, so the provincial Rt trajectory can be
    ## rebuilt for plotting ([`reconstruct_patch_rt`](@ref)).
    delta_knots := vec(patch_state.δ_knots)

    region_sd := patch_state.σ_level
    region_drift_sd := patch_state.σ_δ

    region_halflife := patch_state.δ_halflife

    region_corr_primary_secondary := n_patches > 1 ?
        @inbounds(patch_state.Ω[1, 2]) :
        one(eltype(patch_state.Ω))
    log_rt_contrast := [
        @inbounds(
            patch_state.δ_patch[p, n] -
                patch_state.δ_patch[1, n]
        )
            for p in 1:n_patches
    ]
    ## Population-level dispersion (`k`, the headline scalar) plus the
    ## partially-pooled per-stream dispersions and the pooling SD.
    k := dispersion_state.k_pop
    k_cases := kv[1]
    k_deaths := kv[2]
    k_confirmed := kv[3]
    k_confirmed_deaths := kv[4]
    dispersion_sd := dispersion_state.τ
    p_drc := asc_state.p_drc
    p_uganda := asc_state.p_uganda
    expected_deaths_T := deaths_state.expected_deaths_T
    expected_reports_T := cases_state.expected_reports
    expected_confirmed_T := confirmed_state.expected_confirmed
    expected_analysed_T := confirmed_state.expected_analysed
    _ecd = confirmed_deaths_state.expected_confirmed_deaths
    expected_confirmed_deaths_T := _ecd
    expected_exports_T := exports_state.expected_exports
    expected_exports_deaths_T := exports_deaths_state.expected_exports_deaths_T
    ## Cut-off expected onset-reported total and the modelled per-onset-date
    ## ascertainment level, off the same fitted hazard and ascertainment
    ## walk. The total is reported only, so it is built only when `:=`
    ## values are recorded.
    if _reporting(__varinfo__)
        expected_onset_reported_T := onset_report_expected_total(
            onsets,
            onset_report_state.logit_h0, onset_report_state.γ,
            onset_report_state.grid_start, onset_report_state.alpha, n
        )
    end
    onset_ascertainment := onset_report_state.alpha
    expected_isolation_T := treatment_state.expected_isolation
    expected_bed_demand_T := treatment_state.expected_bed_demand
    bed_shortfall_T := safe_rate(
        treatment_state.expected_bed_demand -
            treatment_state.expected_isolation
    )
    ## Cut-off occupancy split, the confirmed-in-care and suspect-in-care
    ## sub-stock prevalences carved from the occupied true-case stock by the
    ## confirmation overlay.
    expected_confirmed_incare_T := treatment_state.expected_confirmed_incare
    expected_suspect_incare_T := treatment_state.expected_suspect_incare
    ## Cut-off daily treatment flows.
    expected_admissions_T := treatment_state.expected_admissions
    expected_incare_deaths_T := treatment_state.expected_incare_deaths
    expected_ruleouts_T := treatment_state.expected_ruleouts
    bed_capacity := treatment_state.capacity
    isolation_admission := treatment_state.p_iso
    isolation_bvd_admission := treatment_state.p_iso_bvd
    isolation_severity := treatment_state.δ_iso
    ## BVD bed stay outcome mixture: `isolation_bvd_los_mean` is the mixture
    ## mean (overall length-of-stay), with the death and recovery branch means
    ## surfaced separately.
    isolation_bvd_los_mean := treatment_state.overall_los
    isolation_death_los_mean := treatment_state.death_los_mean
    isolation_recovery_los_mean := treatment_state.recovery_los_mean
    isolation_ruleout_los_mean := treatment_state.ruleout_los_mean
    isolation_admission_delay_mean := treatment_state.admission_delay_mean
    isolation_dispersion := treatment_state.k_isolation
    ## In-care fatality CFR_iso, a modifier on the infection CFR, and the
    ## abscond fraction.
    incare_cfr := treatment_state.CFR_iso
    incare_cfr_modifier := treatment_state.β_iso
    abscond_fraction := treatment_state.abscond_frac
    ## In-care confirmation-rate modifier ρ on the borrowed community
    ## confirmation hazard, identified by the confirmed/suspected-in-care split.
    incare_confirm_modifier := treatment_state.incare_confirm_modifier
    expected_recovered_T := recovered_state.expected_recovered
    recovery_probability := recovered_state.p_recover
    recovery_delay_mean := recovered_state.recovery_delay_mean
    recovered_dispersion := recovered_state.k_recovered
    tau_test := cases_state.τ_test
    ## Specimens analysed per suspect sampled. `1.0` when the factor is off.
    specimens_per_suspect := confirmed_state.κ_test === nothing ? 1.0 :
        confirmed_state.κ_test
    lambda_bg := cases_state.λ_bg
    bg_sigma := cases_state.bg_sigma
    background_total := cases_state.bg_total
    death_ascertainment := deaths_state.p_death
    background_cfr := deaths_state.cfr_bg
    lambda_bg_death := deaths_state.λ_bg_death
    bg_death_sigma := deaths_state.bg_death_sigma
    background_death_total := deaths_state.bg_death_total
    tau_death := confirmed_deaths_state.τ_death
    death_testing_scaling := confirmed_deaths_state.scaling
    suspected_positivity := cases_state.positivity
    test_positivity := confirmed_state.p_positive
    death_composition := confirmed_deaths_state.q_death
    death_confirmation := confirmed_deaths_state.p_death_conf

    ## Past the cut-off each stream draws its future counts through its own
    ## likelihood, as `missing` observations `predict` generates.
    forecast_means = nothing
    if forecast !== nothing
        fd = forecast_days(n, forecast)
        vintages = future_knot_days(n, horizon_days(forecast))
        forecast_reports ~ to_submodel(
            _forecast_counts(cases_state.reports_daily, fd, k_cases)
        )
        forecast_deaths ~ to_submodel(
            _forecast_counts(deaths_state.deaths_daily, fd, k_deaths)
        )
        ## A future day has no published analysed count, so its confirmed
        ## cases are the negative binomial the late windows take.
        forecast_confirmed ~ to_submodel(
            _forecast_counts(confirmed_state.confirmed_daily, fd, k_confirmed)
        )
        forecast_confirmed_deaths ~ to_submodel(
            _forecast_counts(
                confirmed_deaths_state.confirmed_death_daily, fd,
                k_confirmed_deaths
            )
        )
        forecast_recovered ~ to_submodel(
            _forecast_counts(recovered_state.recovered_daily, fd, k_recovered)
        )
        treatment_forecast ~ to_submodel(
            treatment_forecast_model(
                treatment_state, fd, bed_capacity_history,
                isolation_history, k_isolation
            ), false
        )
        forecast_exports ~ to_submodel(_forecast_exports(exports_state, fd))
        forecast_latent_deaths := deaths_state.bvd_deaths_daily[fd]
        onset_means = nothing
        if !isempty(onset_curve_history.onset_days)
            onset_forecast ~ to_submodel(
                onset_forecast_model(
                    onsets, onset_report_state, n, vintages
                ), false
            )
            onset_means = onset_forecast.means
        end
        ## Each province's share of the national forecast, one split a
        ## week, by the fitted composition: the same delays, ascertainment
        ## and overdispersion the province tables are fitted with, so the
        ## provinces add up to the national counts drawn above.
        forecast_infections_patch := vec(patch_state.infections_matrix[:, fd])
        if _reporting(__varinfo__)
            frf = _patch_fractions(patch_state)
            forecast_rt_patch := vec(
                permutedims(
                    reduce(
                        hcat,
                        [
                            _patch_adjusted_rt(patch_state, frf, p, fd)
                                for p in 1:n_patches
                        ]
                    )
                )
            )
        end
        edges = vcat(n, vintages)
        if !isempty(province_days)
            province_future = _patch_confirmed_increments(
                confirmed_carried, confirmed_state.s_test, edges
            )[:, 2:end]
            forecast_province_split ~ to_submodel(
                composition_split_model(
                    missing,
                    composition_shares(
                        composition_state.province_ascertainment .*
                            composition_state.province_severity,
                        province_future
                    ),
                    _weekly_totals(forecast_confirmed.increments, n, edges),
                    composition_state.rho
                )
            )
            forecast_province_confirmed := vec(
                forecast_province_split.obs_increments
            )
        end
        if !isempty(province_death_days)
            province_deaths_future = _patch_death_increments(
                patch_state.onsets_matrix, death_kernel, edges
            )[:, 2:end]
            forecast_province_death_split ~ to_submodel(
                composition_split_model(
                    missing,
                    composition_shares(
                        death_composition_state.province_ascertainment .*
                            death_composition_state.province_severity,
                        province_deaths_future
                    ),
                    _weekly_totals(
                        forecast_confirmed_deaths.increments, n, edges
                    ),
                    death_composition_state.rho
                )
            )
            forecast_province_deaths := vec(
                forecast_province_death_split.obs_increments
            )
        end
        ## Each province's occupancy at each future vintage, the national
        ## occupancy split by the fitted occupancy split over the per-patch
        ## demand, and its beds, its static share of the national capacity.
        if _has_province_rows(province_isolation)
            vj = vintages .- n
            forecast_province_isolation_split ~ to_submodel(
                composition_split_model(
                    missing,
                    treatment_state.demand_patch[:, vintages] ./
                        reshape(
                        vec(sum(treatment_state.demand_patch[:, vintages]; dims = 1)),
                        1, :
                    ),
                    ## A draw at the bed ceiling is censored, so not whole;
                    ## `forecast_reported` rounds it the same way.
                    round.(Int, treatment_forecast.occupancy[vj]),
                    treatment_state.occupancy_split_rho
                )
            )
            forecast_province_isolation := vec(
                forecast_province_isolation_split.obs_increments
            )
        end
        if length(treatment_state.capacity_shares) == n_patches
            forecast_province_beds := vec(
                treatment_state.capacity_shares .*
                    reshape(treatment_state.capacity_series[vintages], 1, :)
            )
        end
        forecast_means = (;
            reports = forecast_reports.modelled,
            deaths = forecast_deaths.modelled,
            confirmed = forecast_confirmed.modelled,
            confirmed_deaths = forecast_confirmed_deaths.modelled,
            recovered = forecast_recovered.modelled,
            treatment_forecast.isolation, treatment_forecast.admissions,
            onset_reports = onset_means,
        )
    end
    return (;
        patch_state, onsets, cases_state, deaths_state, confirmed_state,
        confirmed_deaths_state, treatment_state, recovered_state,
        exports_state, onset_report_state, forecast_means,
    )
end
