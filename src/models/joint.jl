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
        rt_start::Integer = 1, rt_walk_start::Integer = rt_start)
    infection_state ~ to_submodel(
        infection(n; breakpoint, rt_start, rt_walk_start), false)
    onset_state ~ to_submodel(
        onset_incidence(infection_state.infections), false)
    ## Shared latent-trajectory deterministics, exposed once here rather than
    ## repeated in every composer: `_latent` is attached unprefixed
    ## (`to_submodel(_latent(...), false)`) by all composers and by
    ## [`bvd_joint`](@ref), so these surface bare in each. `cumulative` is the
    ## renewal `cumsum` already computed by [`infection_model`](@ref) (no
    ## recompute), and `C_T` its cut-off value (`cumulative[n]`).
    ## `cumulative_onsets` is here rather than in [`bvd_joint`](@ref) so a
    ## single-stream fit carries the latent onset trajectory too: the onset
    ## nowcast ([`forecast_onsets`](@ref)) reads its increments, and the
    ## onsets-only fit is exactly the one that has no other route to them.
    cumulative_infections := infection_state.cumulative
    C_T := infection_state.C_T
    cumulative_onsets := cumsum(onset_state.onsets)
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
        dispersion = surveillance_dispersion_model())
    latent ~ to_submodel(
        _latent(n, breakpoint, infection, onset_incidence), false)
    dispersion_state ~ to_submodel(dispersion)
    deaths_state ~ to_submodel(
        deaths(deaths_history, total_deaths, latent.onsets,
        dispersion_state.k; suspected_daily_deaths_history))
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

Exposes the cut-off expected confirmed count as `expected_confirmed_T`,
the same un-prefixed name [`bvd_joint`](@ref) uses, so the confirmed
stream can be forecast from this fit ([`forecast_stream`](@ref)).
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
        tests_analysed, confirmed_break_days,
        confirmed_break_gross = confirmed_break_gross_cases,
        confirmed_break_sd,
        positivity_link = confirmed_positivity_link))
    ## Cut-off expected confirmed count, aliased under the same un-prefixed
    ## name [`bvd_joint`](@ref) uses so both fit kinds carry one key and the
    ## confirmed stream can be forecast from this fit
    ## ([`forecast_stream`](@ref)). Aliased here, at the composer level,
    ## rather than `:=`-tracked inside [`confirmed_cases_model`](@ref): that
    ## submodel keeps its derived quantities on plain `=` because a `:=`
    ## there would build a DynamicPPL tracking closure over the boxed
    ## `p_pos`, which Enzyme's `nodecayed_phis!` pass cannot differentiate
    ## through. The alias below closes over the submodel's returned
    ## NamedTuple, which is assigned once and not boxed.
    expected_confirmed_T := confirmed_state.expected_confirmed
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
        ascertainment = pooled_ascertainment_model())
    latent ~ to_submodel(
        _latent(n, breakpoint, infection, onset_incidence), false)
    dispersion_state ~ to_submodel(dispersion)
    asc_state ~ to_submodel(ascertainment)
    cfr_state ~ to_submodel(cfr)
    k = dispersion_state.k
    p_drc = asc_state.p_drc
    cases_state ~ to_submodel(
        cases((; days = Int[], counts = Int[]), missing, latent.onsets,
        k, p_drc))
    ## Confirmed-case lab pipeline, run so the treatment model can borrow the
    ## daily testing intensity and positivity for the in-care confirmation
    ## overlay.
    confirmed_state ~ to_submodel(
        confirmed(confirmed_history, confirmed_cases, latent.onsets, k,
        p_drc, cases_state.bg_daily, cases_state.τ_test,
        cases_state.bvd_reports_daily;
        lab_history, lab_daily_history, tests_analysed,
        confirmed_break_days,
        confirmed_break_gross = confirmed_break_gross_cases,
        confirmed_break_sd))
    ## In-care confirmation hazard `τ_test · p_pos` on the daily grid.
    conf_hazard_daily = confirmed_state.τ_test .* confirmed_state.p_pos_grid
    treatment_state ~ to_submodel(
        treatment(isolation_history, cases_state.bvd_reports_daily,
        cases_state.bg_daily, p_drc, cfr_state.CFR;
        capacity_history = bed_capacity_history,
        admissions_history = treatment_admissions_history,
        deaths_history = treatment_deaths_history,
        ruleout_history = treatment_ruleout_history,
        absconded_history = treatment_absconded_history,
        confirmed_incare_history = treatment_confirmed_incare_history,
        suspect_incare_history = treatment_suspect_incare_history,
        occupancy_break_days = occupancy_break_days,
        conf_hazard_daily = conf_hazard_daily))
end

"""
Onsets-only composer (the direct-observation analogue). Runs the infection
process and onset staging, then conditions on the symptom-onset reporting-
triangle likelihood alone. See [`onset_reporting_model`](@ref) for the
delay hazard, calendar-time drift, right-truncation and ascertainment
maths.

Unlike every other single-stream composer, this stream needs no dispersion
submodel of its own; `onset_report` (the reporting-triangle model, which
samples its own ascertainment level) is the only injected submodel,
mirroring [`deaths_only_model`](@ref)'s shape. No confirmed pipeline is
available to anchor ascertainment on, so `onset_report` falls back to its
constant `0.15` anchor.

Exposes the cut-off expected onset-reported count as
`expected_onset_reported_T`, the same un-prefixed name [`bvd_joint`](@ref)
uses, so the two carry one key and can be compared directly (see
`expected_confirmed_T` on [`confirmed_only_model`](@ref) for the
precedent), and the modelled ascertainment as `onset_ascertainment`. With
the shared `cumulative_onsets` trajectory from `_latent` that is everything
[`forecast_onsets`](@ref) needs, so this fit nowcasts and forecasts the
onset stream like any other, scored on the reported increment rather than
on the digitised level (see [`forecast_stream`](@ref)).
"""
@model function onsets_only_model(n::Integer;
        onset_curve_history = (; onset_days = Int[], report_days = Int[],
            prev_report_days = Int[], increments = Int[]),
        breakpoint::Union{Missing, Real} = missing,
        infection = infection_model,
        onset_incidence = onset_incidence_model,
        onset_report = onset_reporting_model)
    latent ~ to_submodel(
        _latent(n, breakpoint, infection, onset_incidence), false)
    onset_report_state ~ to_submodel(
        onset_report(onset_curve_history, latent.onsets))
    expected_onset_reported_T := onset_report_expected_total(
        latent.onsets, onset_report_state.logit_h0, onset_report_state.γ,
        onset_report_state.grid_start, onset_report_state.alpha, n)
    onset_ascertainment := onset_report_state.alpha
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
    deaths_state ~ to_submodel(
        deaths(deaths_history, total_deaths, latent.onsets, k;
        case_bg_daily = cases_state.bg_daily))
    confirmed_deaths_state ~ to_submodel(
        confirmed_deaths_stream(confirmed_deaths, total_deaths,
        deaths_state.deaths_daily, deaths_state.bvd_deaths_daily,
        deaths_state.bg_death_daily, k;
        confirmed_deaths_history, confirmed_break_days,
        confirmed_break_gross = confirmed_break_gross_deaths,
        confirmed_break_sd))
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
end

"""
Exports-joint composer: the Uganda export cases and deaths fit together
as one geographic-spread stream. Runs the infection process and onset
staging, samples ascertainment and the deaths submodel (for the CFR and
onset-to-death delay), then conditions on both the export-case and
export-death likelihoods over the one travel-gated at-risk prevalence, so
the imports and the import deaths inform the outbreak size jointly rather
than as two separate single-stream fits. Either count may be `missing` to
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
        ascertainment = pooled_ascertainment_model())
    latent ~ to_submodel(
        _latent(n, breakpoint, infection, onset_incidence), false)
    dispersion_state ~ to_submodel(dispersion)
    asc_state ~ to_submodel(ascertainment)
    deaths_state ~ to_submodel(
        deaths((; days = Int[], counts = Int[]), missing, latent.onsets,
        dispersion_state.k))
    exports_state ~ to_submodel(
        exports(exported_cases, latent.infection_state.infections,
        asc_state.p_uganda; export_case_days,
        incubation_pmf = latent.incubation_pmf, source_population))
    exports_deaths_state ~ to_submodel(
        exports_deaths_model(exports_deaths,
        exports_state.travelled_prevalence, deaths_state.CFR,
        deaths_state.od_pmf, latent.incubation_pmf; export_death_days))
end

"""
Joint composer over all data streams. Runs the generating infection
process once on a daily grid of length `n` (day `n` is the cut-off),
stages it to daily onset incidence, then conditions on the DRC suspected
cases, deaths and the laboratory pipeline (the analysed-specimen volume as
a per-vintage time series and the confirmed positives as a Binomial of the
observed analysed denominator), the confirmed deaths, the Uganda exports
and deaths-among-exports, the digitised symptom-onset reporting triangle
(the only direct observation of the shared onset series, see
[`onset_reporting_model`](@ref); not opt-in, degrades to a no-op when
`onset_curve_history` is empty), and the optional genetic seeding bound on
the outbreak age. Each stream argument may be `missing` to drop it, so the
model doubles as a prior- and posterior-predictive generator.

The confirmed-case stream shares its onset-to-report kernel with the
suspected-case stream. A single analysed-specimen volume is fit through a
report-to-analysed delay and the tested fraction; the confirmed positives
are scored as a Binomial of the observed specimens-analysed denominator
(`lab_history`) with a partially-pooled per-window positivity, so the
confirmed counts do not pass through the multiplicative ascertainment
ridge. After the national cumulative analysed series stops, the
reporting format gives a 24h analysed count on some days
(`lab_daily_history`); these are fitted as per-day analysed volumes and
also anchor that day's confirmed positives as a Binomial of the observed
denominator. The early and unanchored windows (days with no published
denominator) use the modelled analysed volume as the denominator, with
the positivity (hence `λ_bg`) carried over from the windows that do have
data (see [`confirmed_cases_model`](@ref)).

The optional `suspected_daily_history` adds the post-26 May daily
new-suspect inflow ("nouveaux cas suspects du jour"), scored against the
modelled daily suspected series at each report day where the frozen
cumulative suspected stream stops, on days disjoint from it. The optional
`suspected_daily_deaths_history` adds the deaths analogue, the post-26
May daily new suspected deaths ("cas suspects du jour N (M deces)"),
scored against the modelled daily suspected-death series where the
frozen cumulative suspected-death stream stops.

The confirmed deaths mirror the confirmed-case laboratory pipeline: a
death "analysed" volume (the suspected deaths carried to laboratory
receipt and thinned by the death testing fraction `tau_death`) scored
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
walk (e.g. the first WHO situation report); `genetic` injects the genetic
seeding submodel when `tmrca_days` is given. Tracked deterministics:
`C_T` (cumulative infections by the cut-off), the established
reproduction number `R0` (= the first `R_t`), `r` and `doubling_time`
(current growth), `r0` (the `R0`-implied cryptic growth rate), `T`
(outbreak age), `R_T` (current reproduction number), the per-stream
expected counts, the testing fraction `tau_test`, the background rate
`lambda_bg`, the death ascertainment `death_ascertainment`, the
background CFR `background_cfr`, the death testing fraction `tau_death`,
the implied per-suspected (`suspected_positivity`) and per-test
(`test_positivity`) positivities, and the death-pool BVD composition
(`death_composition`) and death-confirmation positivity
(`death_confirmation`).
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
        onset_curve_history = (; onset_days = Int[], report_days = Int[],
            prev_report_days = Int[], increments = Int[]),
        breakpoint::Union{Missing, Real} = missing,
        source_population::Real = ITURI_POPULATION,
        infection = infection_model,
        onset_incidence = onset_incidence_model,
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
        background_re::Bool = false,
        confirmed_positivity_link::Symbol = :composition,
        genetic = nothing,
        onset_to_sample = nejm_onset_to_sample(),
        tmrca_days::Union{Missing, Real} = missing,
        tmrca_days_sd::Real = 16.0,
        renewal_start_lead::Integer = RENEWAL_START_LEAD,
        rt_walk_lead::Integer = RT_WALK_LEAD)
    ## The renewal start sits `renewal_start_lead` days after the genetic
    ## TMRCA day (`n - tmrca_days + lead`), past the TMRCA's uncertainty
    ## where sustained transmission is confident. The lead keeps the
    ## observed span `τ_obs = n − renewal_start` strictly shorter than
    ## `tmrca_days`, so the genetic bound on the total age
    ## `T = m·τ + τ_obs` stays informative (it bounds the cryptic duration
    ## `m·τ` from below). The renewal seeds and grows from here.
    rt_start = ismissing(tmrca_days) ? 1 :
               clamp(n - round(Int, tmrca_days) + renewal_start_lead, 1, n)
    ## Start the random walk `rt_walk_lead` days (a month by default) before
    ## the first situation report (`breakpoint`) rather than exactly at it,
    ## so R_t is free to move over the weeks of transmission leading up to
    ## the first report instead of being held flat at R0 right to it (the
    ## response decline can begin before the outbreak is first reported).
    ## The start is floored at the renewal start so the walk never precedes
    ## the seeded trajectory. With no breakpoint the walk falls back to the
    ## renewal start.
    rt_walk_start = ismissing(breakpoint) ? rt_start :
                    clamp(round(Int, breakpoint) - rt_walk_lead, rt_start, n)
    latent ~ to_submodel(
        _latent(n, breakpoint, infection, onset_incidence;
            rt_start, rt_walk_start), false)
    infection_state = latent.infection_state
    onsets = latent.onsets

    ## Partially-pooled per-stream dispersions: every count stream draws its
    ## own negative-binomial dispersion from a shared population rather than
    ## sharing one global `k`, so a stream's noise is not pulled around by
    ## whichever stream dominates the likelihood while the sparse streams
    ## still borrow strength. Order: 1 suspected cases, 2 suspected deaths,
    ## 3 confirmed cases, 4 confirmed deaths, 5 isolation occupancy,
    ## 6 recovered. The isolation and recovered dispersions are injected into
    ## their submodels (which sample their own only when run standalone).
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

    ## Non-BVD background as a smooth daily lognormal random walk over the
    ## surveillance window ([`background_walk_model`](@ref)), with the tight
    ## innovation SD `σ_rw` driving the suspected-case stream. The background
    ## is gated to zero before the surveillance onset (a report-to-receipt
    ## lead before the first suspected-case report): it does not exist
    ## before surveillance began. The tight innovation SD keeps it fairly
    ## constant, which regularises the background/outbreak-size degeneracy:
    ## a per-vintage step random effect's multiplicative blow-up can open a
    ## second posterior mode that breaks convergence. The suspected-death
    ## background is not a separate random effect: it is tied to the case
    ## background by a background CFR (`cfr_bg · case_bg_daily`, see
    ## [`deaths_model`](@ref)), so it inherits the case background's smooth,
    ## gated, ramped level and time-variation rather than competing as a
    ## second free, outbreak-size-degenerate rate. With `background_re =
    ## false` (the renewal default) the case stream keeps its scalar `λ_bg`.
    ## The pooling SD is sampled only when the random effect is active (the
    ## tilde must stay gated), but `σ_rw_shared` and `bg_onset` are bound
    ## unconditionally to plain values so the closure below captures
    ## non-conditional, write-once variables. Conditionally-scoped captures
    ## get boxed in a `Base.RefValue`, which Enzyme's reverse mode cannot
    ## differentiate through (Mooncake tolerates it); binding them in both
    ## paths keeps the capture type-stable and the closure un-boxed.
    if background_re
        bg_pool ~ to_submodel(background_pooling_model())
        σ_rw_shared = bg_pool.σ_bg
    else
        σ_rw_shared = 0.0
    end
    ## Onset of the suspected pool's non-BVD background: a report-to-receipt
    ## lead before the first suspected-case report, not exactly at it. The
    ## suspects in the first report were already in the pipeline, and the
    ## background feeds the laboratory analysed volume through the
    ## report-to-receipt convolution, so it must begin early enough for that
    ## convolution to be fully formed by the first report. The lead is the
    ## max lag of the report-to-receipt kernel (its truncation `nmax`, the
    ## default `lab_delay_model` support), not its mean, so no tail
    ## contribution is cut off at the onset. Bound unconditionally (unused
    ## when the effect is off).
    bg_lead = cdf_nmax(lognormal_meansd(4.5, 4.0))
    bg_onset = isempty(reported_history.days) ? 1 :
               clamp(Int(reported_history.days[1]) - bg_lead, 1, n)
    ## The closure is built unconditionally (one concrete closure type, not
    ## closure-or-`Nothing`); `background_re` then selects the closure or the
    ## `nothing` sentinel. Identical behaviour to the gated form: with
    ## `background_re = false` the closure is never passed, so the unused
    ## `σ_rw_shared = 0` never enters the log-density.
    make_case_bg = nn -> background_walk_model(nn, σ_rw_shared;
        onset = bg_onset)
    case_bg_re = background_re ? make_case_bg : nothing

    ## Cases first so the suspected-case background `bg_daily` is available to
    ## the deaths stream (which scales it by `cfr_bg` for the death background)
    ## and to the laboratory pipeline.
    cases_state ~ to_submodel(
        cases(reported_history, reported_cases, onsets, k_cases, p_drc;
        suspected_daily_history, background_re = case_bg_re))
    deaths_state ~ to_submodel(
        deaths(deaths_history, total_deaths, onsets, k_deaths;
        suspected_daily_deaths_history, case_bg_daily = cases_state.bg_daily))
    confirmed_state ~ to_submodel(
        confirmed(confirmed_history, confirmed_cases, onsets, k_confirmed,
        p_drc, cases_state.bg_daily, cases_state.τ_test,
        cases_state.bvd_reports_daily;
        lab_history, lab_daily_history,
        tests_analysed, confirmed_break_days,
        confirmed_break_gross = confirmed_break_gross_cases,
        confirmed_break_sd,
        positivity_link = confirmed_positivity_link))
    ## Symptom-onset reporting-triangle stream
    ## ([`onset_reporting_model`](@ref)): the only direct observation of the
    ## shared latent onset series `onsets`, every other stream sees it only
    ## after a further onset-to-event convolution. Runs after
    ## `confirmed_state` so its daily ascertainment
    ## (`p_drc · τ_test · p_pos_grid`) is available to anchor this stream's
    ## own ascertainment level on.
    onset_anchor_daily = p_drc .* confirmed_state.τ_test .*
                         confirmed_state.p_pos_grid
    onset_report_state ~ to_submodel(
        onset_report(onset_curve_history, onsets;
        anchor = onset_anchor_daily))
    ## Confirmed deaths mirror the confirmed-case lab pipeline: the death
    ## analysed volume scales the modelled case analysed volume
    ## (`confirmed_state.analysed_daily`) at the per-day suspected
    ## death-to-case ratio, scored through a death-pool composition positivity
    ## from the death series' own BVD and background components. The case
    ## volume carries the laboratory capacity onset, so the death volume
    ## inherits it and no deaths are confirmed before testing began.
    confirmed_deaths_state ~ to_submodel(
        confirmed_deaths_stream(confirmed_deaths, total_deaths,
        deaths_state.deaths_daily, deaths_state.bvd_deaths_daily,
        deaths_state.bg_death_daily, k_confirmed_deaths;
        confirmed_deaths_history, receipt_pmf = confirmed_state.receipt_pmf,
        confirmed_break_days,
        confirmed_break_gross = confirmed_break_gross_deaths,
        confirmed_break_sd,
        case_analysed_daily = confirmed_state.analysed_daily,
        case_suspected_daily = cases_state.reports_daily))
    ## Treatment-centre patient flow ([`treatment_flow_model`](@ref)):
    ## occupancy plus the in-care outcome flows, with the in-care fatality
    ## CFR_iso identified by the in-care death flow. The occupancy split
    ## borrows the in-care confirmation hazard `τ_test · p_pos` from the
    ## confirmed pipeline to carve the occupied true-case stock into
    ## confirmed and suspect sub-stocks, scored against the Tableau 6
    ## `dont confirmés` / `dont suspects` census. The known DHIS2
    ## harmonisation days carry the overnight total reporting break.
    conf_hazard_daily = confirmed_state.τ_test .* confirmed_state.p_pos_grid
    treatment_state ~ to_submodel(
        treatment(isolation_history, cases_state.bvd_reports_daily,
        cases_state.bg_daily, p_drc, deaths_state.CFR;
        capacity_history = bed_capacity_history,
        admissions_history = treatment_admissions_history,
        deaths_history = treatment_deaths_history,
        ruleout_history = treatment_ruleout_history,
        absconded_history = treatment_absconded_history,
        confirmed_incare_history = treatment_confirmed_incare_history,
        suspect_incare_history = treatment_suspect_incare_history,
        occupancy_break_days = occupancy_break_days,
        conf_hazard_daily = conf_hazard_daily,
        k_external = k_isolation))
    ## Recovered among confirmed ("cumul guéris"): survivors among the
    ## modelled daily confirmed cases (the confirmed-and-discharged subset,
    ## not all in-care recoveries), with a recovery fraction grounded on the
    ## CFR and lagged by a confirmation-to-recovery delay (see
    ## [`recovered_model`](@ref)).
    recovered_state ~ to_submodel(
        recovered(recovered_history, recovered_cases,
        confirmed_state.confirmed_daily, deaths_state.CFR;
        k_external = k_recovered))
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
    ## reconstruct from the chain without re-running the renewal. All three
    ## are BVD-only latent renewal quantities: deaths uses the BVD death
    ## series (onsets convolved with the onset-to-death delay), not the
    ## fitted total, so it stays smooth like infections and onsets. The
    ## additive non-BVD background is a daily random walk and belongs to the
    ## observation side, not this latent trajectory. `cumulative_infections`,
    ## `cumulative_onsets` and `C_T` are exposed once by the shared `_latent`
    ## submodel above.
    cumulative_expected_deaths := cumsum(deaths_state.bvd_deaths_daily)
    ## Modelled daily laboratory-confirmed cases (from
    ## `confirmed_cases_model`: the per-window tested-positive probability
    ## applied to the modelled, testing-onset-gated analysed volume), so the
    ## cumulative trajectory carries the confirmed-case timing for the
    ## delay-corrected confirmed-CFR reconstruction. The onset-to-confirmation
    ## kernel (onset-to-report ⊕ receipt) and the onset-to-death-confirmation
    ## kernel (onset-to-death ⊕ receipt) are exposed alongside so the residual
    ## delay between a confirmed case and its confirmed death can be rebuilt
    ## per draw off the chain.
    ## Re-add the testing-onset baseline: the laboratory capacity is gated to
    ## zero before testing began and the first confirmed vintage is treated as
    ## the initial condition (a baseline the early windows do not score), so the
    ## reconstructed cumulative counts only the fitted increments. Adding the
    ## first observed confirmed count back from the testing onset onward makes
    ## the trajectory comparable to the observed confirmed total (and keeps the
    ## delay-corrected confirmed-CFR denominator on the right level).
    _conf_inc_cum = cumsum(confirmed_state.confirmed_daily)
    _conf_base = isempty(confirmed_history.counts) ? 0 :
                 Int(confirmed_history.counts[1])
    _conf_cap = isempty(confirmed_history.days) ? 1 :
                clamp(Int(confirmed_history.days[1]), 1, n)
    _conf_base_vec = [t >= _conf_cap ? _conf_base : 0 for t in 1:n]
    cumulative_confirmed := _conf_inc_cum .+ _conf_base_vec
    onset_to_confirmation_pmf := convolve_pmf(
        cases_state.report_pmf, confirmed_state.receipt_pmf)
    onset_to_death_confirmation_pmf := convolve_pmf(
        deaths_state.od_pmf, confirmed_state.receipt_pmf)
    ## External onset-to-sample constraint on the confirmed sampling delay
    ## (grounded on the NEJM DRC 2026 cohort by default, see
    ## [`nejm_onset_to_sample`](@ref)). The onset→report and report→receipt
    ## legs convolve to the confirmed onset-to-sample delay, so its
    ## continuous mean is the sum of the two legs' means and its continuous
    ## SD the root-sum of their variances; both are exposed here. The
    ## cohort's reported (continuous) mean
    ## and median are fitted to these as soft Normal observations
    ## ([`onset_to_sample_logweight`](@ref)), grounding the otherwise-
    ## unidentified receipt (lab-turnaround) leg without touching either prior.
    ## The term only exists on the confirmed report⊕receipt path, so single-
    ## stream and isolation composers carry none; passing `nothing` drops it.
    onset_to_sample_mean := cases_state.report_mean +
                            confirmed_state.receipt_mean
    onset_to_sample_sd := sqrt(cases_state.report_sd^2 +
                               confirmed_state.receipt_sd^2)
    if onset_to_sample !== nothing
        @addlogprob! onset_to_sample_logweight(cases_state.report_mean,
            cases_state.report_sd, confirmed_state.receipt_mean,
            confirmed_state.receipt_sd, onset_to_sample)
    end
    R0 := infection_state.R0
    r := infection_state.r
    r0 := infection_state.r0
    doubling_time := infection_state.doubling_time
    T := infection_state.T
    R_T := infection_state.Rt[n]
    expected_infections_T := infection_state.infections[n]
    CFR := deaths_state.CFR
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
    ## walk; see [`onset_reporting_model`](@ref) for what the vintage
    ## structure does and does not separate here.
    expected_onset_reported_T := onset_report_expected_total(
        onsets, onset_report_state.logit_h0, onset_report_state.γ,
        onset_report_state.grid_start, onset_report_state.alpha, n)
    onset_ascertainment := onset_report_state.alpha
    expected_isolation_T := treatment_state.expected_isolation
    expected_bed_demand_T := treatment_state.expected_bed_demand
    bed_shortfall_T := safe_rate(treatment_state.expected_bed_demand -
                                 treatment_state.expected_isolation)
    ## Cut-off occupancy split: the confirmed-in-care and suspect-in-care
    ## sub-stock prevalences carved from the occupied true-case stock by the
    ## confirmation overlay.
    expected_confirmed_incare_T := treatment_state.expected_confirmed_incare
    expected_suspect_incare_T := treatment_state.expected_suspect_incare
    ## Cut-off daily treatment flows surfaced for the one-week-ahead forecast.
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
    ## In-care fatality CFR_iso (a modifier on the infection CFR) and the
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
end

## --- Patch (multi-population) joint models ------------------------------

"""
Patch latent process: run [`patch_infection_model`](@ref) and expose
per-patch infections, onsets, and the total (summed-across-patches)
infection trajectory that feeds the national observation submodels.

Returns `(; patch_state, onsets_total, cumulative_total, C_T_total)` where
`C_T_total` is the national total (sum of patch cut-off cumulatives).
"""
@model function _patch_latent(n::Integer, n_patches::Integer,
        breakpoint, patch_infection;
        rt_start::Integer = 1,
        rt_walk_start::Integer = rt_start,
        importation_kernel::AbstractMatrix = zeros(n_patches, n_patches))
    patch_state ~ to_submodel(
        patch_infection(n, n_patches;
            breakpoint, rt_start, rt_walk_start,
            importation_kernel), false)
    ## Sum the per-patch onsets to get the national total.
    onsets_total = vec(sum(patch_state.onsets_matrix; dims = 1))
    cumulative_total = vec(sum(patch_state.cumulative_matrix; dims = 1))
    C_T_total_out := patch_state.C_T_total
    return (; patch_state, onsets_total, cumulative_total, C_T_total = C_T_total_out)
end

"""
Joint composer over all data streams using the PATCH (meta-population)
latent process. Runs [`patch_infection_model`](@ref) with the given
number of patches and importation kernel, sums the patch onsets into a
national trajectory, then conditions on all the existing national-level
observation submodels: DRC suspected cases, deaths, confirmed cases and
deaths, laboratory pipeline, treatment flows, and Uganda exports.

Per-province observation models are NOT yet wired; only national-level
streams are fitted, but the per-patch state (`patch_state`) is surfaced
in the return tuple for diagnostics and downstream extension. See also
[`bvd_joint`](@ref) for the single-patch analogue and
[`patch_infection_model`](@ref) for the patch latent process.
"""
@model function bvd_patch_joint(
        n::Integer, n_patches::Integer,
        exported_cases::Union{Missing, Integer},
        total_deaths::Union{Missing, Integer},
        reported_cases::Union{Missing, Integer} = missing,
        exports_deaths::Union{Missing, Integer} = missing,
        confirmed_cases::Union{Missing, Integer} = missing,
        tests_analysed::Union{Missing, Integer} = missing;
        importation_kernel::AbstractMatrix = zeros(n_patches, n_patches),
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
        occupancy_break_days::AbstractVector{<:Integer} = Int[],
        export_case_days::AbstractVector{<:Integer} = Int[],
        export_death_days::AbstractVector{<:Integer} = Int[],
        breakpoint::Union{Missing, Real} = missing,
        source_population::Real = ITURI_POPULATION,
        patch_infection = patch_infection_model,
        exports = exports_model,
        deaths = deaths_model,
        cases = reported_cases_model,
        confirmed = confirmed_cases_model,
        confirmed_deaths_stream = confirmed_deaths_model,
        treatment = treatment_flow_model,
        recovered = recovered_model,
        dispersion = pooled_dispersion_model,
        ascertainment = pooled_ascertainment_model(),
        background_re::Bool = false,
        confirmed_positivity_link::Symbol = :composition,
        genetic = nothing,
        onset_to_sample = nejm_onset_to_sample(),
        tmrca_days::Union{Missing, Real} = missing,
        tmrca_days_sd::Real = 15.0,
        renewal_start_lead::Integer = RENEWAL_START_LEAD,
        rt_walk_lead::Integer = RT_WALK_LEAD)
    ## Renewal-start and walk-start: same logic as [`bvd_joint`](@ref).
    rt_start = ismissing(tmrca_days) ? 1 :
               clamp(n - round(Int, tmrca_days) + renewal_start_lead, 1, n)
    rt_walk_start = ismissing(breakpoint) ? rt_start :
                    clamp(round(Int, breakpoint) - rt_walk_lead, rt_start, n)
    latent ~ to_submodel(
        _patch_latent(n, n_patches, breakpoint, patch_infection;
            rt_start, rt_walk_start, importation_kernel), false)
    patch_state = latent.patch_state
    onsets = latent.onsets_total
    ## Partially-pooled per-stream dispersions (same layout as bvd_joint).
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
    ## Background (safe `bg_onset` with guard against empty histories,
    ## matching the pattern from bvd_joint). The lead is the MAX lag of
    ## the report-to-receipt kernel (its truncation `nmax`), using 7 as
    ## an approximation for `cdf_nmax(lognormal_meansd(4.5, 4.0))`.
    case_bg_re = background_re ?
    begin
        bg_pool ~ to_submodel(background_pooling_model())
        σ_rw_shared = bg_pool.σ_bg
        bg_lead = 7
        bg_onset = isempty(reported_history.days) ? 1 :
                   clamp(Int(reported_history.days[1]) - bg_lead, 1, n)
        nn -> background_walk_model(nn, σ_rw_shared;
            onset = bg_onset)
    end : nothing
    ## --- Observation submodels (all national-level, same as bvd_joint) ---
    ## 1. Reported (suspected) cases.
    cases_state ~ to_submodel(cases(reported_history, reported_cases,
        onsets, k_cases, p_drc;
        suspected_daily_history, background_re = case_bg_re))
    ## 2. Deaths (suspected).
    deaths_state ~ to_submodel(deaths(deaths_history, total_deaths, onsets, k_deaths;
        suspected_daily_deaths_history, case_bg_daily = cases_state.bg_daily))
    ## 3. Confirmed cases (laboratory pipeline).
    confirmed_state ~ to_submodel(confirmed(confirmed_history,
        confirmed_cases, onsets, k_confirmed, p_drc,
        cases_state.bg_daily, cases_state.τ_test,
        cases_state.bvd_reports_daily;
        lab_history, lab_daily_history, tests_analysed,
        positivity_link = confirmed_positivity_link))
    ## 4. Confirmed deaths.
    confirmed_deaths_state ~ to_submodel(
        confirmed_deaths_stream(confirmed_deaths, total_deaths,
        deaths_state.deaths_daily, deaths_state.bvd_deaths_daily,
        deaths_state.bg_death_daily, k_confirmed_deaths;
        confirmed_deaths_history, receipt_pmf = confirmed_state.receipt_pmf,
        case_analysed_daily = confirmed_state.analysed_daily,
        case_suspected_daily = cases_state.reports_daily))
    ## 5. Uganda exports (cases and deaths).
    exports_state ~ to_submodel(exports(exported_cases,
        patch_state.infections_matrix[1, :], p_uganda;
        export_case_days, incubation_pmf = patch_state.incubation_pmf,
        source_population))
    exports_deaths_state ~ to_submodel(exports_deaths_model(
        exports_deaths, exports_state.travelled_prevalence,
        deaths_state.CFR, deaths_state.od_pmf,
        patch_state.incubation_pmf;
        export_death_days))
    ## Exports from Ituri (patch 1) only, matching the single-patch model.
    ## 6. Treatment flows (isolation, bed capacity, LOS). Uses the summed
    ##    (national) suspected-case inflow, matching the national-level data.
    treatment_state ~ to_submodel(treatment(isolation_history,
        cases_state.bvd_reports_daily, cases_state.bg_daily, p_drc,
        deaths_state.CFR;
        capacity_history = bed_capacity_history,
        admissions_history = treatment_admissions_history,
        deaths_history = treatment_deaths_history,
        ruleout_history = treatment_ruleout_history,
        absconded_history = treatment_absconded_history,
        occupancy_break_days = occupancy_break_days,
        k_external = k_isolation))
    ## 7. Recovered (among confirmed).
    recovered_state ~ to_submodel(recovered(recovered_history,
        recovered_cases, confirmed_state.confirmed_daily,
        deaths_state.CFR))
    ## --- Deterministics surfaced for reporting --------------------------
    ## Patch-level parameters are already traced by the submodel; we surface
    ## only those that are NOT already `:=` or `~` inside the nested models.
    R0 := patch_state.R0
    r := patch_state.r
    ## Derived national Rt at the cut-off (exposed here with a distinct name
    ## since the inner `:=` is hidden by to_submodel(..., false)).
    implied_Rt_national_cutoff := patch_state.implied_Rt_national[n]
    ## Per-patch outbreak summaries at the cut-off.
    C_T_patch_1 := patch_state.C_T_patch[1]
    C_T_patch_2 := patch_state.C_T_patch[2]
    C_T_patch_3 := patch_state.C_T_patch[3]
    R_T_patch_1 := patch_state.Rt_matrix[1, n]
    R_T_patch_2 := patch_state.Rt_matrix[2, n]
    R_T_patch_3 := patch_state.Rt_matrix[3, n]
    infections_T_patch_1 := patch_state.infections_matrix[1, n]
    infections_T_patch_2 := patch_state.infections_matrix[2, n]
    infections_T_patch_3 := patch_state.infections_matrix[3, n]
    k := dispersion_state.k_pop
    k_cases := kv[1]
    k_deaths := kv[2]
    k_confirmed := kv[3]
    k_confirmed_deaths := kv[4]
    dispersion_sd := dispersion_state.τ
    p_drc := asc_state.p_drc
    p_uganda := asc_state.p_uganda
    CFR := deaths_state.CFR
end
