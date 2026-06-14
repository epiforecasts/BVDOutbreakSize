# Joint composer models: build the full generative model for each
# analysis by running the generating infection process once, staging it to
# daily onset incidence, and routing the shared onsets into the relevant
# observation submodels. Single-stream composers condition on one stream
# each; [`bvd_joint`](@ref) conditions on all streams plus the optional
# genetic seeding bound. Any count passed as `missing` is dropped, so the
# composers double as prior- and posterior-predictive generators.
#
# Submodels whose `:=` deterministics are re-exposed at composer level are
# attached with a prefixed `to_submodel(x)` (no `false`); attaching them
# with `to_submodel(x, false)` would re-introduce the nested `:=` names at
# the parent and trip Turing's MustNotOverwriteError on the duplicates.

## Run the generating infection process and onset staging, returning the
## infection state and the daily onsets shared by every stream.
@model function _latent(n::Integer, breakpoint, infection, onset_incidence;
        rt_start::Integer = 1)
    infection_state ~ to_submodel(infection(n; breakpoint, rt_start), false)
    onset_state ~ to_submodel(
        onset_incidence(infection_state.infections), false)
    return (; infection_state, onsets = onset_state.onsets,
        incubation_pmf = onset_state.incubation_pmf)
end

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
        ascertainment = pooled_ascertainment_model())
    latent ~ to_submodel(
        _latent(n, breakpoint, infection, onset_incidence), false)
    asc_state ~ to_submodel(ascertainment)
    exports_state ~ to_submodel(
        exports(exported_cases, latent.infection_state.infections,
        asc_state.p_uganda; export_case_days,
        incubation_pmf = latent.incubation_pmf,
        source_population))
    cumulative_infections := cumsum(latent.infection_state.infections)
    C_T := latent.infection_state.C_T
end

"""
Deaths-only composer (back-calculation analogue). Runs the infection
process and onset staging, samples dispersion, then conditions on the
deaths likelihood only. See [`deaths_model`](@ref).
"""
@model function deaths_only_model(
        n::Integer, total_deaths::Union{Missing, Integer};
        deaths_history = (; days = Int[], counts = Int[]),
        breakpoint::Union{Missing, Real} = missing,
        infection = infection_model,
        onset_incidence = onset_incidence_model,
        deaths = deaths_model,
        dispersion = surveillance_dispersion_model())
    latent ~ to_submodel(
        _latent(n, breakpoint, infection, onset_incidence), false)
    dispersion_state ~ to_submodel(dispersion)
    deaths_state ~ to_submodel(
        deaths(deaths_history, total_deaths, latent.onsets,
        dispersion_state.k))
    cumulative_infections := cumsum(latent.infection_state.infections)
    C_T := latent.infection_state.C_T
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
        ascertainment = pooled_ascertainment_model())
    latent ~ to_submodel(
        _latent(n, breakpoint, infection, onset_incidence), false)
    dispersion_state ~ to_submodel(dispersion)
    asc_state ~ to_submodel(ascertainment)
    cases_state ~ to_submodel(
        cases(reported_history, reported_cases, latent.onsets,
        dispersion_state.k, asc_state.p_drc; suspected_daily_history))
    cumulative_infections := cumsum(latent.infection_state.infections)
    C_T := latent.infection_state.C_T
end

"""
Confirmed-cases-only composer (laboratory pipeline in isolation). Runs
the infection process and onset staging, samples dispersion and pooled
ascertainment, then runs the suspected-case stream (in predictive mode,
to draw the shared background rate, testing fraction and onset-to-report
kernel) and conditions on the laboratory pipeline alone: the confirmed
positives (a Binomial of the observed analysed denominator in
`lab_history`) and the modelled analysed-specimen volume. See
[`confirmed_cases_model`](@ref) and [`reported_cases_model`](@ref).
"""
@model function confirmed_only_model(
        n::Integer, confirmed_cases::Union{Missing, Integer};
        confirmed_history = (; days = Int[], counts = Int[]),
        lab_history = (; days = Int[], counts = Int[]),
        lab_daily_history = (; days = Int[], counts = Int[]),
        tests_analysed::Union{Missing, Integer} = missing,
        breakpoint::Union{Missing, Real} = missing,
        infection = infection_model,
        onset_incidence = onset_incidence_model,
        cases = reported_cases_model,
        confirmed = confirmed_cases_model,
        dispersion = surveillance_dispersion_model(),
        ascertainment = pooled_ascertainment_model(),
        confirmed_positivity_link::Symbol = :composition)
    latent ~ to_submodel(
        _latent(n, breakpoint, infection, onset_incidence), false)
    dispersion_state ~ to_submodel(dispersion)
    asc_state ~ to_submodel(ascertainment)
    k = dispersion_state.k
    p_drc = asc_state.p_drc
    cases_state ~ to_submodel(
        cases((; days = Int[], counts = Int[]), missing, latent.onsets,
        k, p_drc))
    confirmed_state ~ to_submodel(
        confirmed(confirmed_history, confirmed_cases, latent.onsets, k,
        p_drc, cases_state.bg_daily, cases_state.τ_test,
        cases_state.bvd_reports_daily;
        lab_history, lab_daily_history,
        tests_analysed,
        positivity_link = confirmed_positivity_link))
    cumulative_infections := cumsum(latent.infection_state.infections)
    C_T := latent.infection_state.C_T
end

"""
Isolation-occupancy-only composer (treatment-bed prevalence in isolation).
Runs the infection process and onset staging, samples dispersion and pooled
ascertainment, then runs the suspected-case stream in predictive mode (to
draw the shared background rate, testing fraction and onset-to-report
kernel) and conditions on the isolation/treatment-bed occupancy alone. See
[`treatment_admission_model`](@ref) and [`reported_cases_model`](@ref).
"""
@model function treatment_only_model(
        n::Integer;
        isolation_history = (; days = Int[], counts = Int[]),
        breakpoint::Union{Missing, Real} = missing,
        infection = infection_model,
        onset_incidence = onset_incidence_model,
        cases = reported_cases_model,
        treatment = treatment_admission_model,
        receipt = lab_delay_model(),
        dispersion = surveillance_dispersion_model(),
        ascertainment = pooled_ascertainment_model())
    latent ~ to_submodel(
        _latent(n, breakpoint, infection, onset_incidence), false)
    dispersion_state ~ to_submodel(dispersion)
    asc_state ~ to_submodel(ascertainment)
    k = dispersion_state.k
    p_drc = asc_state.p_drc
    cases_state ~ to_submodel(
        cases((; days = Int[], counts = Int[]), missing, latent.onsets,
        k, p_drc))
    ## Report-to-receipt laboratory delay for the non-BVD rule-out stay (the
    ## confirmed pipeline that supplies it in the joint is not run here).
    receipt_state ~ to_submodel(receipt)
    treatment_state ~ to_submodel(
        treatment(isolation_history, cases_state.bvd_reports_daily,
        cases_state.bg_daily, p_drc, receipt_state.pmf))
    cumulative_infections := cumsum(latent.infection_state.infections)
    C_T := latent.infection_state.C_T
end

"""
Confirmed-deaths-only composer. Runs the infection process and onset
staging, samples dispersion and pooled ascertainment, runs the suspected
deaths and reported-cases streams (the latter in predictive mode for the
shared background rate and onset-to-report kernel), then conditions on the
confirmed-death thinning alone. See [`confirmed_deaths_model`](@ref).
"""
@model function confirmed_deaths_only_model(
        n::Integer, confirmed_deaths::Union{Missing, Integer},
        total_deaths::Union{Missing, Integer} = missing;
        deaths_history = (; days = Int[], counts = Int[]),
        confirmed_deaths_history = (; days = Int[], counts = Int[]),
        breakpoint::Union{Missing, Real} = missing,
        infection = infection_model,
        onset_incidence = onset_incidence_model,
        deaths = deaths_model,
        cases = reported_cases_model,
        confirmed_deaths_stream = confirmed_deaths_model,
        dispersion = surveillance_dispersion_model(),
        ascertainment = pooled_ascertainment_model())
    latent ~ to_submodel(
        _latent(n, breakpoint, infection, onset_incidence), false)
    dispersion_state ~ to_submodel(dispersion)
    asc_state ~ to_submodel(ascertainment)
    k = dispersion_state.k
    p_drc = asc_state.p_drc
    deaths_state ~ to_submodel(
        deaths(deaths_history, total_deaths, latent.onsets, k))
    cases_state ~ to_submodel(
        cases((; days = Int[], counts = Int[]), missing, latent.onsets,
        k, p_drc))
    confirmed_deaths_state ~ to_submodel(
        confirmed_deaths_stream(confirmed_deaths, total_deaths,
        deaths_state.deaths_daily, cases_state.bvd_reports_daily,
        p_drc, cases_state.bg_daily, k;
        confirmed_deaths_history))
    cumulative_infections := cumsum(latent.infection_state.infections)
    C_T := latent.infection_state.C_T
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
        ascertainment = pooled_ascertainment_model())
    latent ~ to_submodel(
        _latent(n, breakpoint, infection, onset_incidence), false)
    dispersion_state ~ to_submodel(dispersion)
    asc_state ~ to_submodel(ascertainment)
    deaths_state ~ to_submodel(
        deaths((; days = Int[], counts = Int[]), missing, latent.onsets,
        dispersion_state.k))
    exports_state ~ to_submodel(
        exports(missing, latent.infection_state.infections,
        asc_state.p_uganda; incubation_pmf = latent.incubation_pmf,
        source_population))
    exports_deaths_state ~ to_submodel(
        exports_deaths_model(exports_deaths,
        exports_state.travelled_prevalence, deaths_state.CFR,
        deaths_state.od_pmf, latent.incubation_pmf; export_death_days))
    cumulative_infections := cumsum(latent.infection_state.infections)
    C_T := latent.infection_state.C_T
end

"""
Joint composer over all data streams. Runs the generating infection
process once on a daily grid of length `n` (day `n` is the cut-off),
stages it to daily onset incidence, then conditions on the DRC suspected
cases, deaths and the laboratory pipeline (the analysed-specimen volume as
a per-vintage time series and the confirmed positives as a Binomial of the
observed analysed denominator), the confirmed deaths, the Uganda exports
and deaths-among-exports, and the optional genetic seeding bound on the
outbreak age. Each stream argument may be `missing` to drop it, so the
model doubles as a prior- and posterior-predictive generator.

onset-to-report kernel with the suspected-case stream. A single
analysed-specimen volume is fit through a report-to-analysed delay and the
tested fraction; the confirmed positives are scored as a Binomial of the
observed specimens-analysed denominator (`lab_history`) with a
partially-pooled per-window positivity, so the confirmed counts no longer
pass through the multiplicative ascertainment ridge. After the national
cumulative analysed series stops, the reporting format gives a 24h analysed
count on some days (`lab_daily_history`); these are fitted as per-day
analysed volumes and also anchor that day's confirmed positives as a
Binomial of the observed denominator. The early and unanchored windows
(days with no published denominator) use the modelled analysed volume as
the denominator, with the positivity (hence `λ_bg`) carried over from the
windows that do have data (see [`confirmed_cases_model`](@ref)). The
optional `suspected_daily_history` adds the post-26 May daily new-suspect
inflow ("nouveaux cas suspects du jour"), scored against the modelled daily
suspected series at each report day where the frozen cumulative suspected
stream stops, on days disjoint from it. The confirmed deaths are a thinning of
the suspected deaths whose confirmation probability is the suspected-case
BVD composition enriched on the odds scale (`confirmed_deaths`,
`total_deaths` the denominator). The optional `isolation_history` adds the
daily isolation/treatment-bed occupancy ("Patients en isolement"), a
prevalence stream fitted as the suspect inflow (BVD treatment stay plus
non-BVD rule-out stay) carried through a length-of-stay survival into a
daily stock (see [`treatment_admission_model`](@ref)). The optional
`recovered_history` adds the recovered-among-confirmed stream ("cumul
guéris"), survivors among the modelled daily confirmed cases scaled by the
recovery probability and lagged by a confirmation-to-recovery delay (see
[`recovered_model`](@ref)).

`breakpoint` is the intervention day passed to the reproduction-number
walk (e.g. the first WHO situation report); `genetic` injects the genetic
seeding submodel when `tmrca_days` is given. Tracked deterministics:
`C_T` (cumulative infections by the cut-off), the single established
reproduction number `R0` (= the first `R_t`), `r` and `doubling_time`
(current growth), `r0` (the `R0`-implied cryptic growth rate), `T`
(outbreak age),
`R_T` (current reproduction number), the per-stream expected counts, the
testing fraction `tau_test`, the background rate `lambda_bg`, the
confirmed-death enrichment `m_death`, the implied per-suspected
(`suspected_positivity`) and per-test (`test_positivity`) positivities,
and the suspected BVD composition (`death_composition`) and
death-confirmation probability (`death_confirmation`).
"""
@model function bvd_joint(
        n::Integer,
        exported_cases::Union{Missing, Integer},
        total_deaths::Union{Missing, Integer},
        reported_cases::Union{Missing, Integer} = missing,
        exports_deaths::Union{Missing, Integer} = missing,
        confirmed_cases::Union{Missing, Integer} = missing,
        tests_analysed::Union{Missing, Integer} = missing;
        confirmed_deaths::Union{Missing, Integer} = missing,
        recovered_cases::Union{Missing, Integer} = missing,
        deaths_history = (; days = Int[], counts = Int[]),
        reported_history = (; days = Int[], counts = Int[]),
        confirmed_history = (; days = Int[], counts = Int[]),
        confirmed_deaths_history = (; days = Int[], counts = Int[]),
        lab_history = (; days = Int[], counts = Int[]),
        lab_daily_history = (; days = Int[], counts = Int[]),
        suspected_daily_history = (; days = Int[], counts = Int[]),
        isolation_history = (; days = Int[], counts = Int[]),
        recovered_history = (; days = Int[], counts = Int[]),
        export_case_days::AbstractVector{<:Integer} = Int[],
        export_death_days::AbstractVector{<:Integer} = Int[],
        breakpoint::Union{Missing, Real} = missing,
        source_population::Real = ITURI_POPULATION,
        infection = infection_model,
        onset_incidence = onset_incidence_model,
        exports = exports_model,
        deaths = deaths_model,
        cases = reported_cases_model,
        confirmed = confirmed_cases_model,
        confirmed_deaths_stream = confirmed_deaths_model,
        treatment = treatment_admission_model,
        recovered = recovered_model,
        dispersion = surveillance_dispersion_model(),
        ascertainment = pooled_ascertainment_model(),
        background_re::Bool = false,
        confirmed_positivity_link::Symbol = :composition,
        genetic = nothing,
        tmrca_days::Union{Missing, Real} = missing,
        tmrca_days_sd::Real = 15.0,
        renewal_start_lead::Integer = RENEWAL_START_LEAD)
    ## Fix R_t at R0 over the pre-establishment seeding window, letting the
    ## random walk vary R_t only from the renewal start onward — before that
    ## point the outbreak dynamics are unidentified, so a free walk there only
    ## adds unsupported drift. The renewal start sits `renewal_start_lead`
    ## days AFTER the genetic TMRCA day (`n - tmrca_days + lead`), past the
    ## TMRCA's uncertainty where sustained transmission is confident. The lead
    ## keeps the observed span `τ_obs = n − renewal_start` strictly shorter
    ## than `tmrca_days`, so the genetic bound on the total age
    ## `T = m·τ + τ_obs` stays informative (it bounds the cryptic duration
    ## `m·τ` from below).
    rt_start = ismissing(tmrca_days) ? 1 :
               clamp(n - round(Int, tmrca_days) + renewal_start_lead, 1, n)
    latent ~ to_submodel(
        _latent(n, breakpoint, infection, onset_incidence; rt_start), false)
    infection_state = latent.infection_state
    onsets = latent.onsets

    dispersion_state ~ to_submodel(dispersion)
    asc_state ~ to_submodel(ascertainment)
    k = dispersion_state.k
    p_drc = asc_state.p_drc
    p_uganda = asc_state.p_uganda

    ## Per-vintage background random effect, with ONE pooling SD `σ_bg`
    ## shared between the suspected-case and suspected-death streams (the
    ## same non-BVD reporting environment drives both). With `background_re
    ## = false` (the renewal default) the case stream keeps its scalar
    ## `λ_bg` and the death stream has no background. Each stream still
    ## samples its own baseline; the death baseline prior is tighter (deaths
    ## are far fewer than suspected cases).
    if background_re
        bg_pool ~ to_submodel(background_pooling_model())
        σ_bg_shared = bg_pool.σ_bg
        case_bg_re = nv -> background_re_model(nv, σ_bg_shared)
        death_bg_re = nv -> background_re_model(nv, σ_bg_shared;
            baseline_prior = truncated(Normal(0.0, 0.25); lower = 0))
    else
        case_bg_re = nothing
        death_bg_re = nothing
    end

    deaths_state ~ to_submodel(
        deaths(deaths_history, total_deaths, onsets, k;
        background_re = death_bg_re))
    cases_state ~ to_submodel(
        cases(reported_history, reported_cases, onsets, k, p_drc;
        suspected_daily_history, background_re = case_bg_re))
    confirmed_state ~ to_submodel(
        confirmed(confirmed_history, confirmed_cases, onsets, k, p_drc,
        cases_state.bg_daily, cases_state.τ_test,
        cases_state.bvd_reports_daily;
        lab_history, lab_daily_history,
        tests_analysed,
        positivity_link = confirmed_positivity_link))
    confirmed_deaths_state ~ to_submodel(
        confirmed_deaths_stream(confirmed_deaths, total_deaths,
        deaths_state.deaths_daily, cases_state.bvd_reports_daily,
        p_drc, cases_state.bg_daily, k;
        confirmed_deaths_history, receipt_pmf = confirmed_state.receipt_pmf))
    ## Isolation/treatment-bed occupancy: the suspect inflow (BVD treatment
    ## stay plus non-BVD rule-out stay) carried through a length-of-stay
    ## survival into a daily stock (see [`treatment_admission_model`](@ref)).
    treatment_state ~ to_submodel(
        treatment(isolation_history, cases_state.bvd_reports_daily,
        cases_state.bg_daily, p_drc, confirmed_state.receipt_pmf))
    ## Recovered among confirmed ("cumul guéris"): survivors among the modelled
    ## daily confirmed cases, scaled by the recovery probability and lagged by
    ## a confirmation-to-recovery delay (see [`recovered_model`](@ref)).
    recovered_state ~ to_submodel(
        recovered(recovered_history, recovered_cases,
        confirmed_state.confirmed_daily))
    exports_state ~ to_submodel(
        exports(exported_cases, infection_state.infections, p_uganda;
        export_case_days, incubation_pmf = latent.incubation_pmf,
        source_population))
    exports_deaths_state ~ to_submodel(
        exports_deaths_model(exports_deaths,
        exports_state.travelled_prevalence, deaths_state.CFR,
        deaths_state.od_pmf, latent.incubation_pmf; export_death_days))

    if genetic !== nothing
        genetic_state ~ to_submodel(
            genetic(infection_state.T, tmrca_days; tmrca_days_sd), false)
    end

    ## Daily cumulative trajectories for the headline 3x2 figure: the
    ## modelled expected cumulative infections, symptom onsets and deaths
    ## over the grid. Exposed as vector deterministics so the ribbon panels
    ## reconstruct from the chain without re-running the renewal.
    cumulative_infections := cumsum(infection_state.infections)
    cumulative_onsets := cumsum(onsets)
    cumulative_expected_deaths := cumsum(deaths_state.deaths_daily)
    ## Modelled daily laboratory-confirmed cases (from `confirmed_cases_model`:
    ## the per-window tested-positive probability expanded onto the daily grid
    ## and applied to the modelled analysed volume), so the cumulative
    ## trajectory carries the confirmed-case timing for the delay-corrected
    ## confirmed-CFR reconstruction. The onset-to-confirmation kernel
    ## (onset-to-report ⊕ receipt) and the onset-to-death-confirmation kernel
    ## (onset-to-death ⊕ receipt, the death carrying the same report-to-receipt
    ## laboratory delay) are exposed alongside so the residual delay between a
    ## confirmed case and its confirmed death can be rebuilt per draw off the
    ## chain.
    cumulative_confirmed := cumsum(confirmed_state.confirmed_daily)
    onset_to_confirmation_pmf := convolve_pmf(cases_state.report_pmf, confirmed_state.receipt_pmf)
    onset_to_death_confirmation_pmf := convolve_pmf(deaths_state.od_pmf, confirmed_state.receipt_pmf)
    C_T := infection_state.C_T
    R0 := infection_state.R0
    r := infection_state.r
    r0 := infection_state.r0
    doubling_time := infection_state.doubling_time
    T := infection_state.T
    R_T := infection_state.Rt[n]
    expected_infections_T := infection_state.infections[n]
    CFR := deaths_state.CFR
    k := dispersion_state.k
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
    expected_isolation_T := treatment_state.expected_isolation
    isolation_admission := treatment_state.p_iso
    isolation_bvd_los_mean := treatment_state.bvd_los_mean
    isolation_dispersion := treatment_state.k_isolation
    expected_recovered_T := recovered_state.expected_recovered
    recovery_probability := recovered_state.p_recover
    recovery_delay_mean := recovered_state.recovery_delay_mean
    recovered_dispersion := recovered_state.k_recovered
    tau_test := cases_state.τ_test
    lambda_bg := cases_state.λ_bg
    bg_sigma := cases_state.bg_sigma
    background_total := cases_state.bg_total
    lambda_bg_death := deaths_state.λ_bg_death
    bg_death_sigma := deaths_state.bg_death_sigma
    background_death_total := deaths_state.bg_death_total
    m_death := confirmed_deaths_state.m_death
    suspected_positivity := cases_state.positivity
    test_positivity := confirmed_state.p_positive
    death_composition := confirmed_deaths_state.q_susp
    death_confirmation := confirmed_deaths_state.p_death_conf
end
