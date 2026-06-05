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
        asc_state.p_uganda; incubation_pmf = latent.incubation_pmf,
        source_population))
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
        dispersion_state.k, asc_state.p_drc))
    C_T := latent.infection_state.C_T
end

"""
Confirmed-cases-only composer (laboratory pipeline in isolation). Runs
the infection process and onset staging, samples dispersion and pooled
ascertainment, then runs the suspected-case stream (in predictive mode,
to draw the shared background rate, testing fraction and onset-to-report
kernel) and conditions on the laboratory pipeline alone: the confirmed
positives (a Binomial of the observed analysed denominator in
`lab_history`) and, when present, the received-specimen stream. See
[`confirmed_cases_model`](@ref) and [`reported_cases_model`](@ref).
"""
@model function confirmed_only_model(
        n::Integer, confirmed_cases::Union{Missing, Integer};
        confirmed_history = (; days = Int[], counts = Int[]),
        lab_history = (; days = Int[], counts = Int[]),
        tests_received_history = (; days = Int[], counts = Int[]),
        breakpoint::Union{Missing, Real} = missing,
        infection = infection_model,
        onset_incidence = onset_incidence_model,
        cases = reported_cases_model,
        confirmed = confirmed_cases_model,
        dispersion = surveillance_dispersion_model(),
        ascertainment = pooled_ascertainment_model(),
        confirmed_positivity_link::Symbol = :free)
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
        lab_history, tests_received_history,
        positivity_link = confirmed_positivity_link))
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
        exports_state.export_prevalence, deaths_state.CFR,
        deaths_state.od_pmf, latent.incubation_pmf))
    C_T := latent.infection_state.C_T
end

"""
Joint composer over all data streams. Runs the generating infection
process once on a daily grid of length `n` (day `n` is the cut-off),
stages it to daily onset incidence, then conditions on the DRC suspected
cases, deaths and the laboratory pipeline (the received-specimen volume as
a per-vintage time series and the confirmed positives as a Binomial of the
observed analysed denominator), the confirmed deaths, the Uganda exports
and deaths-among-exports, and the optional genetic seeding bound on the
outbreak age. Each stream argument may be `missing` to drop it, so the
model doubles as a prior- and posterior-predictive generator.

The laboratory pipeline shares the testing fraction, background rate and
onset-to-report kernel with the suspected-case stream. The received
specimens are fit through a receipt delay and the tested fraction; the
confirmed positives are scored as a Binomial of the observed
specimens-analysed denominator (`lab_history`) with a partially-pooled
per-window positivity, so the confirmed counts no longer pass through the
multiplicative ascertainment ridge. The confirmed deaths are a thinning of
the suspected deaths whose confirmation probability is the suspected-case
BVD composition enriched on the odds scale (`confirmed_deaths`,
`total_deaths` the denominator).

`breakpoint` is the intervention day passed to the reproduction-number
walk (e.g. the first WHO situation report); `genetic` injects the genetic
seeding submodel when `tmrca_days` is given. Tracked deterministics:
`C_T` (cumulative infections by the cut-off), `r` and `doubling_time`
(current growth), `r0` (implied initial growth), `T` (outbreak age),
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
        deaths_history = (; days = Int[], counts = Int[]),
        reported_history = (; days = Int[], counts = Int[]),
        confirmed_history = (; days = Int[], counts = Int[]),
        confirmed_deaths_history = (; days = Int[], counts = Int[]),
        lab_history = (; days = Int[], counts = Int[]),
        tests_received_history = (; days = Int[], counts = Int[]),
        breakpoint::Union{Missing, Real} = missing,
        source_population::Real = ITURI_POPULATION,
        infection = infection_model,
        onset_incidence = onset_incidence_model,
        exports = exports_model,
        deaths = deaths_model,
        cases = reported_cases_model,
        confirmed = confirmed_cases_model,
        confirmed_deaths_stream = confirmed_deaths_model,
        dispersion = surveillance_dispersion_model(),
        ascertainment = pooled_ascertainment_model(),
        background_re::Bool = false,
        confirmed_positivity_link::Symbol = :free,
        genetic = nothing,
        tmrca_days::Union{Missing, Real} = missing,
        tmrca_days_sd::Real = 15.0)
    ## Fix R_t at R0 over the pre-establishment seeding window, letting the
    ## random walk vary R_t only from the conservative genetic TMRCA bound
    ## (`n - tmrca_days`) onward — before that point the outbreak dynamics
    ## are unidentified, so a free walk there only adds unsupported drift.
    rt_start = ismissing(tmrca_days) ? 1 :
               clamp(n - round(Int, tmrca_days), 1, n)
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
        background_re = case_bg_re))
    confirmed_state ~ to_submodel(
        confirmed(confirmed_history, confirmed_cases, onsets, k, p_drc,
        cases_state.bg_daily, cases_state.τ_test,
        cases_state.bvd_reports_daily;
        lab_history, tests_received_history,
        positivity_link = confirmed_positivity_link))
    confirmed_deaths_state ~ to_submodel(
        confirmed_deaths_stream(confirmed_deaths, total_deaths,
        deaths_state.deaths_daily, cases_state.bvd_reports_daily,
        p_drc, cases_state.bg_daily, k;
        confirmed_deaths_history))
    exports_state ~ to_submodel(
        exports(exported_cases, infection_state.infections, p_uganda;
        incubation_pmf = latent.incubation_pmf, source_population))
    exports_deaths_state ~ to_submodel(
        exports_deaths_model(exports_deaths,
        exports_state.export_prevalence, deaths_state.CFR,
        deaths_state.od_pmf, latent.incubation_pmf))

    if genetic !== nothing
        genetic_state ~ to_submodel(
            genetic(infection_state.T, tmrca_days; tmrca_days_sd), false)
    end

    C_T := infection_state.C_T
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
    expected_received_T := confirmed_state.expected_received
    _ecd = confirmed_deaths_state.expected_confirmed_deaths
    expected_confirmed_deaths_T := _ecd
    expected_exports_T := exports_state.expected_exports
    expected_exports_deaths_T := exports_deaths_state.expected_exports_deaths_T
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
