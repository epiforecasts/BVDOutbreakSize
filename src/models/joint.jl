# Joint composer models: build the full generative model for each
# analysis by sampling the shared priors once and routing them into
# the relevant observation submodels. Single-stream composers are
# provided for the four count-based streams; [`bvd_joint`](@ref)
# conditions on all streams simultaneously.

"""
Exports-only composer (Method 1 analogue). Samples growth and
ascertainment, then conditions on the exports likelihood only. See
[`exports_model`](@ref).
"""
@model function exports_only_model(
        exported_cases::Union{Missing, Integer};
        growth = exponential_growth_model(),
        exports = exports_model,
        ascertainment = pooled_ascertainment_model())
    growth_state ~ to_submodel(growth, false)
    asc_state ~ to_submodel(ascertainment, false)

    exports_state ~ to_submodel(
        exports(exported_cases, growth_state, asc_state.p_uganda), false)

    cumulative_cases := growth_state.C_T
end

"""
Deaths-only composer (Method 2 analogue). Samples growth and
dispersion, then conditions on the deaths likelihood only. See
[`deaths_model`](@ref).
"""
@model function deaths_only_model(
        total_deaths::Union{Missing, Integer};
        growth = exponential_growth_model(),
        deaths = deaths_model,
        dispersion = surveillance_dispersion_model())
    growth_state ~ to_submodel(growth, false)
    dispersion_state ~ to_submodel(dispersion, false)
    k = dispersion_state.k

    ## Single cumulative total at the cut-off: a length-1 vintage vector
    ## whose only edge is `T`, so the per-vintage deaths likelihood
    ## reduces to the cumulative single-total NegBinomial (Method 2).
    deaths_vec = Union{Missing, Int}[total_deaths]
    deaths_state ~ to_submodel(
        deaths(deaths_vec, growth_state, k, [growth_state.T]), false)

    cumulative_cases := growth_state.C_T
end

"""
Cases-only composer (ascertainment extension). Samples growth,
dispersion and pooled ascertainment, then conditions on the
reported-cases likelihood. See [`reported_cases_model`](@ref).
"""
@model function cases_only_model(
        reported_cases::Union{Missing, Integer};
        growth = exponential_growth_model(),
        reported_cases_submodel = reported_cases_model,
        dispersion = surveillance_dispersion_model(),
        ascertainment = pooled_ascertainment_model())
    growth_state ~ to_submodel(growth, false)
    dispersion_state ~ to_submodel(dispersion, false)
    asc_state ~ to_submodel(ascertainment, false)
    k = dispersion_state.k

    ## Single cumulative total at the cut-off: a length-1 vintage vector
    ## with the pooled scalar ascertainment and the single edge `T`, so
    ## the per-vintage reported likelihood reduces to the cumulative
    ## single-total NegBinomial.
    reported_vec = Union{Missing, Int}[reported_cases]
    reported_state ~ to_submodel(
        reported_cases_submodel(reported_vec, growth_state, k,
            [asc_state.p_drc], [growth_state.T]), false)

    cumulative_cases := growth_state.C_T
end

"""
Confirmed-cases-only composer (laboratory pipeline in isolation).
Samples growth, dispersion, pooled ascertainment, the background /
testing-fraction prior and the report and lab delays, then conditions
on the laboratory pipeline alone. The single cumulative confirmed count
(and optional `tests_analysed`) is wrapped into the length-1 per-vintage
form at the cut-off, so it reduces to the cumulative laboratory
likelihood. See [`confirmed_cases_model`](@ref).
"""
@model function confirmed_only_model(
        confirmed_cases::Union{Missing, Integer},
        tests_analysed::Union{Missing, Integer} = missing;
        growth = exponential_growth_model(),
        confirmed = confirmed_cases_model,
        dispersion = surveillance_dispersion_model(),
        ascertainment = pooled_ascertainment_model(),
        test_positivity = test_positivity_model(),
        report_delay = report_delay_model(),
        lab_delay = lab_delay_model(),
        test_sensitivity = test_sensitivity_model())
    growth_state ~ to_submodel(growth, false)
    dispersion_state ~ to_submodel(dispersion, false)
    asc_state ~ to_submodel(ascertainment, false)
    report_state ~ to_submodel(report_delay, false)
    test_positivity_state ~ to_submodel(test_positivity, false)
    k = dispersion_state.k
    λ_bg = test_positivity_state.λ_bg
    τ_test = test_positivity_state.τ_test
    f_rep = report_state.dist
    T = growth_state.T

    confirmed_vec = Union{Missing, Int}[confirmed_cases]
    confirmed_state ~ to_submodel(
        confirmed(confirmed_vec, tests_analysed, growth_state, k,
            [asc_state.p_drc], λ_bg, τ_test, f_rep, [T], T;
            lab_delay = lab_delay,
            test_sensitivity = test_sensitivity), false)

    cumulative_cases := growth_state.C_T
end

"""
Deaths-among-exports-only composer. Samples growth, delay, CFR,
detection window, traveller volume and ascertainment, then conditions
on the dated export-deaths likelihood. See
[`exports_deaths_model`](@ref).
"""
@model function exports_deaths_only_model(
        export_deaths_daily::AbstractVector;
        growth = exponential_growth_model(),
        delay = delay_model(),
        cfr = cfr_model(),
        window = detection_window_model(),
        traveller = traveller_volume_model(),
        exports_deaths_model = exports_deaths_model,
        ascertainment = pooled_ascertainment_model(),
        source_population::Real = ITURI_POPULATION,
        pre_start_deaths::Union{Missing, Integer} = 0)
    growth_state ~ to_submodel(growth, false)
    delay_state ~ to_submodel(delay, false)
    cfr_state ~ to_submodel(cfr, false)
    window_state ~ to_submodel(window, false)
    asc_state ~ to_submodel(ascertainment, false)

    travel_state ~ to_submodel(traveller, false)
    daily_travellers = travel_state.daily_travellers

    exports_deaths_state ~ to_submodel(
        exports_deaths_model(export_deaths_daily, growth_state,
            cfr_state.CFR, delay_state.dist, asc_state.p_uganda;
            pre_start_deaths = pre_start_deaths,
            window = window_state.w,
            daily_travellers = daily_travellers,
            source_population = source_population),
        false)

    cumulative_cases := growth_state.C_T
end

"""
Joint composer over all data streams. Conditions on exports, the DRC
suspected-deaths, reported (suspected) and laboratory-confirmed case
streams, the dated deaths-among-exports series, and optional
export-detection timing and genetic seeding bound.

The DRC deaths, reported and confirmed streams are fitted per
sitrep vintage: `total_deaths`, `reported_cases` and `confirmed_cases`
are cumulative-count vectors (index 1 = oldest) and the model fits the
between-vintage increments. Each carries its own offsets vector
(`death_offsets`, `reported_offsets`, `confirmed_offsets`; days before
the cut-off), converted to elapsed times `T - offset` internally so the
bin edges track the latent `T`. A length-1 vector with offset `0`
reduces a stream to its cumulative single-total likelihood, recovering
the McCabe et al. Method 2 configuration. Pass a vector of `missing`
entries (with matching offsets) to drop a stream while keeping the model
usable as a prior- and posterior-predictive generator.

DRC ascertainment is a single fixed fraction `p_drc` (the pooled scalar
from [`pooled_ascertainment_model`](@ref)) applied to every reported and
confirmed vintage bin, shared between the two streams.

`tests_analysed` is a single cumulative testing-volume count observed at
its own elapsed time `tests_offset` before the cut-off, so it stays
robust if lab reporting lags or stops before the case cut-off. Per-test
positivity is exposed as a derived quantity rather than re-observing the
confirmed counts conditional on tests. Pass `tests_analysed = missing`
to drop it.

`deaths_ascertainment` samples a multiplicative drift factor `p_deaths`
on the expected-deaths trajectory (see
[`deaths_ascertainment_model`](@ref)); pass `p_deaths_fixed = 1.0` to
disable the factor entirely.
"""
@model function bvd_joint(
        exported_cases::Union{Missing, Integer},
        total_deaths::AbstractVector,
        reported_cases::AbstractVector,
        export_deaths_daily::AbstractVector = Int[];
        reported_offsets::AbstractVector,
        death_offsets::AbstractVector = reported_offsets,
        confirmed_cases::AbstractVector = Union{Missing, Int}[],
        confirmed_offsets::AbstractVector = reported_offsets,
        tests_analysed::Union{Missing, Integer} = missing,
        tests_offset::Real = 0,
        growth = exponential_growth_model(),
        exports = exports_model,
        deaths = deaths_model,
        reported_cases_submodel = reported_cases_model,
        confirmed = confirmed_cases_model,
        exports_deaths_model = exports_deaths_model,
        exports_detection_timing = exports_detection_timing_model,
        dispersion = surveillance_dispersion_model(),
        ascertainment = pooled_ascertainment_model(),
        deaths_ascertainment = deaths_ascertainment_model(),
        p_deaths_fixed::Union{Nothing, Real} = nothing,
        test_positivity = test_positivity_model(),
        report_delay = report_delay_model(),
        lab_delay = lab_delay_model(),
        test_sensitivity = test_sensitivity_model(),
        genetic = nothing,
        source_population::Real = ITURI_POPULATION,
        pre_start_deaths::Union{Missing, Integer} = 0,
        pre_detection_exports::Union{Missing, Integer} = 0,
        first_export_detection_delta::Union{Missing, Real} = missing)
    growth_state ~ to_submodel(growth, false)
    if genetic !== nothing
        genetic_state ~ to_submodel(genetic(growth_state.T), false)
    end
    dispersion_state ~ to_submodel(dispersion, false)
    asc_state ~ to_submodel(ascertainment, false)
    k = dispersion_state.k
    p_uganda = asc_state.p_uganda
    if p_deaths_fixed === nothing
        deaths_asc_state ~ to_submodel(deaths_ascertainment, false)
        p_deaths = deaths_asc_state.p_deaths
    else
        p_deaths = p_deaths_fixed
    end
    T = growth_state.T

    exports_state ~ to_submodel(
        exports(exported_cases, growth_state, p_uganda), false)

    death_edges = [T - δ for δ in death_offsets]
    deaths_state ~ to_submodel(
        deaths(total_deaths, growth_state, k, death_edges;
            p_deaths = p_deaths), false)

    ## Fixed per-stream DRC ascertainment: every reported and confirmed
    ## vintage bin uses the pooled scalar `p_drc`, shared between the two
    ## streams.
    n_rep = length(reported_offsets)
    n_conf = length(confirmed_offsets)
    p_drc_per_bin = fill(asc_state.p_drc, max(n_rep, n_conf))

    reported_edges = [T - δ for δ in reported_offsets]
    reported_state ~ to_submodel(
        reported_cases_submodel(reported_cases, growth_state, k,
            p_drc_per_bin[1:n_rep], reported_edges;
            report_delay = report_delay,
            test_positivity = test_positivity), false)

    if !isempty(confirmed_cases)
        confirmed_edges = [T - δ for δ in confirmed_offsets]
        tests_edge = T - tests_offset
        confirmed_state ~ to_submodel(
            confirmed(confirmed_cases, tests_analysed, growth_state, k,
                p_drc_per_bin[1:n_conf], reported_state.λ_bg,
                reported_state.τ_test, reported_state.report_delay_dist,
                confirmed_edges, tests_edge;
                lab_delay = lab_delay,
                test_sensitivity = test_sensitivity), false)
    end

    exports_deaths_state ~ to_submodel(
        exports_deaths_model(export_deaths_daily, growth_state,
            deaths_state.CFR, deaths_state.delay_dist, p_uganda;
            pre_start_deaths = pre_start_deaths,
            window = exports_state.w,
            daily_travellers = exports_state.daily_travellers,
            source_population = source_population),
        false)
    detection_timing_state ~ to_submodel(
        exports_detection_timing(growth_state, p_uganda;
            delta = first_export_detection_delta,
            pre_detection_exports = pre_detection_exports,
            window = exports_state.w,
            daily_travellers = exports_state.daily_travellers,
            source_population = source_population),
        false)

    cumulative_cases := growth_state.C_T
end

"""
McCabe et al. reimplementation composer: exports and deaths only,
mirroring the joint configuration in the Imperial report (no
reported-cases or deaths-among-exports likelihood). Passing
`missing` for `exported_cases` reduces to a pure Method 2 (deaths-
only) fit.
"""
@model function imperial_only_model(
        exported_cases::Union{Missing, Integer},
        total_deaths::Union{Missing, Integer};
        growth = exponential_growth_model(),
        exports = exports_model,
        deaths = deaths_model,
        dispersion = surveillance_dispersion_model(),
        ascertainment = pooled_ascertainment_model())
    growth_state ~ to_submodel(growth, false)
    dispersion_state ~ to_submodel(dispersion, false)
    asc_state ~ to_submodel(ascertainment, false)
    k = dispersion_state.k
    p_uganda = asc_state.p_uganda

    if !ismissing(exported_cases)
        exports_state ~ to_submodel(
            exports(exported_cases, growth_state, p_uganda), false)
    end
    ## Single cumulative deaths total at the cut-off (Method 2), kept
    ## deliberately as one observation to mirror the McCabe et al.
    ## configuration rather than the per-vintage fit.
    deaths_vec = Union{Missing, Int}[total_deaths]
    deaths_state ~ to_submodel(
        deaths(deaths_vec, growth_state, k, [growth_state.T]), false)

    cumulative_cases := growth_state.C_T
end
