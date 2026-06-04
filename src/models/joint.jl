# Joint composer models: build the full generative model for each
# analysis by sampling the shared priors once and routing them into
# the relevant observation submodels. Single-stream composers are
# provided for the four count-based streams; [`bvd_joint`](@ref)
# conditions on all streams simultaneously.

"""
Exports-only composer (Method 1 analogue). Samples growth and
ascertainment, then conditions on the exports likelihood only. Defaults
to the explicit infection→detection delay mechanism
[`exports_delay_model`](@ref).

The exports stream is travel-gated: a case is at risk of being exported
from infection (the traveller moves during incubation, pre-symptomatic)
and detected abroad after the full infection→detection delay. That delay
is the incubation period (see [`incubation_model`](@ref)) convolved with
the DRC onset-to-report delay (see [`report_delay_model`](@ref)),
moment-matched to one Gamma via [`combined_delay`](@ref): entering
surveillance is taken to be the same process abroad as in the DRC. With
no reported-cases stream here both delays are sampled from their priors;
in [`bvd_joint`](@ref) the same draws are shared with the reported and
incubation streams, which pin them. The incubation period sits inside
`f_det`, so the exports likelihood carries no separate incubation
rescale. `cumulative_cases = C_T · onset_fraction` is still exposed as a
prior-predictive case-equivalent of the latent infections.
"""
@model function exports_only_model(
        exported_cases::Union{Missing, Integer};
        growth = exponential_growth_model(),
        exports = exports_delay_model,
        report_delay = report_delay_model(),
        ascertainment = pooled_ascertainment_model(),
        incubation = incubation_model())
    growth_state ~ to_submodel(growth, false)
    asc_state ~ to_submodel(ascertainment, false)
    report_state ~ to_submodel(report_delay, false)
    incubation_state ~ to_submodel(incubation, false)
    os = onset_rescale(incubation_state.dist, growth_state.r)

    ## Exports are travel-gated, so the at-risk clock runs from infection:
    ## the detection delay is incubation ⊕ onset-to-report, moment-matched
    ## to one Gamma. No separate incubation rescale (it lives in `f_det`).
    f_det = combined_delay(incubation_state.dist, report_state.dist)
    exports_state ~ to_submodel(
        exports(exported_cases, growth_state, asc_state.p_uganda, f_det),
        false)

    onset_fraction := os
    cumulative_infections := growth_state.C_T
    cumulative_cases := growth_state.C_T * os
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
        dispersion = surveillance_dispersion_model(),
        incubation = incubation_model())
    growth_state ~ to_submodel(growth, false)
    dispersion_state ~ to_submodel(dispersion, false)
    incubation_state ~ to_submodel(incubation, false)
    k = dispersion_state.k
    os = onset_rescale(incubation_state.dist, growth_state.r)

    ## Single cumulative total at the cut-off: a length-1 vintage vector
    ## whose only edge is `T`, so the per-vintage deaths likelihood
    ## reduces to the cumulative single-total NegBinomial (Method 2).
    deaths_vec = Union{Missing, Int}[total_deaths]
    deaths_state ~ to_submodel(
        deaths(deaths_vec, growth_state, k, [growth_state.T];
            onset_fraction = os), false)

    onset_fraction := os
    cumulative_infections := growth_state.C_T
    cumulative_cases := growth_state.C_T * os
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
        ascertainment = pooled_ascertainment_model(),
        incubation = incubation_model())
    growth_state ~ to_submodel(growth, false)
    dispersion_state ~ to_submodel(dispersion, false)
    asc_state ~ to_submodel(ascertainment, false)
    incubation_state ~ to_submodel(incubation, false)
    k = dispersion_state.k
    os = onset_rescale(incubation_state.dist, growth_state.r)

    ## Single cumulative total at the cut-off: a length-1 vintage vector
    ## with the pooled scalar ascertainment and the single edge `T`, so
    ## the per-vintage reported likelihood reduces to the cumulative
    ## single-total NegBinomial.
    reported_vec = Union{Missing, Int}[reported_cases]
    reported_state ~ to_submodel(
        reported_cases_submodel(reported_vec, growth_state, k,
            [asc_state.p_drc], [growth_state.T]; onset_fraction = os), false)

    onset_fraction := os
    cumulative_infections := growth_state.C_T
    cumulative_cases := growth_state.C_T * os
end

"""
Confirmed-deaths-only composer (laboratory-confirmed deaths in
isolation). Samples the growth, CFR, onset-to-death delay and incubation
machinery the deaths stream uses to build the modelled BVD-death
trajectory, the non-BVD background `λ_bg_death`, and the shared PCR
sensitivity and specificity, then conditions a single cumulative
confirmed-death total (`Cumul décès parmi les confirmés`) through the
lab/positivity process of [`confirmed_deaths_model`](@ref): the BVD share
of the suspect-death pool (BVD plus background) sets the death-specimen
positivity `s·q_death + (1−spec)(1−q_death)` and a forwarded fraction
`τ_death` of the suspect-death backlog gives the expected confirmed count.
`τ_death` is the free death-specimen forwarding rate; `s` and `spec` are
the confirmed-case sensitivity / specificity. `confirmed_deaths` defaults
to `missing` (posterior-predictive generator); pass an integer to
condition. `susp_deaths` is retained for the public signature but no longer
feeds the likelihood (the suspect-death backlog is modelled, not the
observed suspected count).
"""
@model function confirmed_deaths_only_model(
        susp_deaths::Integer,
        confirmed_deaths::Union{Missing, Integer} = missing;
        growth = exponential_growth_model(),
        confirmed_deaths_submodel = confirmed_deaths_model,
        delay = delay_model(),
        cfr = cfr_model(),
        dispersion = surveillance_dispersion_model(),
        death_background = death_background_model(),
        test_sensitivity = test_sensitivity_model(),
        test_specificity = test_specificity_model(),
        incubation = incubation_model())
    growth_state ~ to_submodel(growth, false)
    delay_state ~ to_submodel(delay, false)
    cfr_state ~ to_submodel(cfr, false)
    dispersion_state ~ to_submodel(dispersion, false)
    death_bg_state ~ to_submodel(death_background, false)
    sensitivity_state ~ to_submodel(test_sensitivity, false)
    specificity_state ~ to_submodel(test_specificity, false)
    incubation_state ~ to_submodel(incubation, false)
    os = onset_rescale(incubation_state.dist, growth_state.r)
    λ_bg_death = death_bg_state.λ_bg_death
    edges = [growth_state.T]
    bvd_death_edges = bvd_death_trajectory(growth_state.r, cfr_state.CFR,
        delay_state.dist, edges; onset_fraction = os)
    nsusp_death_edges = [max(bvd_death_edges[i], zero(eltype(bvd_death_edges))) +
                         max(λ_bg_death * edges[i],
                             zero(λ_bg_death * edges[i]))
                         for i in eachindex(edges)]
    confirmed_vec = Union{Missing, Int}[confirmed_deaths]
    confirmed_deaths_state ~ to_submodel(
        confirmed_deaths_submodel(confirmed_vec, bvd_death_edges,
            nsusp_death_edges, sensitivity_state.s_test,
            specificity_state.spec_test, dispersion_state.k), false)
    cumulative_infections := growth_state.C_T
    cumulative_cases := growth_state.C_T * os
end

"""
Confirmed-cases-only composer (laboratory pipeline in isolation).
Samples growth, dispersion, pooled ascertainment, the background /
testing-fraction prior and the report and lab delays, then conditions
on the laboratory pipeline alone. The single cumulative confirmed count
(and optional `tests_analysed`) is wrapped into the length-1 per-vintage
form at the cut-off, so it reduces to the cumulative laboratory
likelihood. See [`confirmed_cases_model`](@ref).

`tests_received` (cumulative `Cumul échantillons reçus` at the cut-off)
conditions the forwarded fraction `τ_forward` via the received-count
NegBinomial. Pass it when fitting: leaving it `missing` turns the
received count into a sampled discrete latent, which HMC / NUTS cannot
handle (the gradient-based sampler rejects discrete unknowns). The
`missing` default is for predictive generation, where the received count
is drawn from the posterior rather than conditioned on.
"""
@model function confirmed_only_model(
        confirmed_cases::Union{Missing, Integer},
        tests_analysed::Union{Missing, Integer} = missing,
        tests_received::Union{Missing, Integer} = missing;
        growth = exponential_growth_model(),
        confirmed = confirmed_cases_model,
        dispersion = surveillance_dispersion_model(),
        ascertainment = pooled_ascertainment_model(),
        test_positivity = test_positivity_model(),
        report_delay = report_delay_model(),
        test_sensitivity = test_sensitivity_model(),
        test_specificity = test_specificity_model(),
        test_selection = test_selection_model(),
        incubation = incubation_model(),
        report_onset_offset::Union{Nothing, Real} = nothing)
    growth_state ~ to_submodel(growth, false)
    dispersion_state ~ to_submodel(dispersion, false)
    asc_state ~ to_submodel(ascertainment, false)
    report_state ~ to_submodel(report_delay, false)
    test_positivity_state ~ to_submodel(test_positivity, false)
    incubation_state ~ to_submodel(incubation, false)
    k = dispersion_state.k
    λ_bg = test_positivity_state.λ_bg
    f_rep = report_state.dist
    T = growth_state.T
    os = onset_rescale(incubation_state.dist, growth_state.r)

    confirmed_vec = Union{Missing, Int}[confirmed_cases]
    ## Single-vintage analysed denominator for the confirmed Binomial: the
    ## cumulative tests-analysed count at the cut-off.
    analysed_vec = Union{Missing, Int}[tests_analysed]
    received_vec = Union{Missing, Int}[tests_received]
    confirmed_state ~ to_submodel(
        confirmed(confirmed_vec, analysed_vec, received_vec, tests_analysed,
            growth_state, k,
            [asc_state.p_drc], λ_bg, test_positivity_state.τ_forward,
            f_rep, [T], T;
            test_sensitivity = test_sensitivity,
            test_specificity = test_specificity,
            test_selection = test_selection,
            report_onset_offset = report_onset_offset,
            onset_fraction = os), false)

    onset_fraction := os
    cumulative_infections := growth_state.C_T
    cumulative_cases := growth_state.C_T * os
end

"""
Deaths-among-exports-only composer. Samples growth, onset-to-death
delay, CFR, the incubation and onset-to-report delays, traveller volume
and ascertainment, then conditions on the dated export-deaths
likelihood. Defaults to the explicit infection-keyed delay mechanism
[`exports_deaths_delay_model`](@ref). Both the export at-risk clock and
the export-death timing run from infection: the infection→detection
delay is incubation ⊕ onset-to-report and the infection→death delay is
incubation ⊕ onset-to-death, each moment-matched to one Gamma via
[`combined_delay`](@ref); see [`exports_deaths_delay_model`](@ref).
"""
@model function exports_deaths_only_model(
        export_deaths_daily::AbstractVector;
        growth = exponential_growth_model(),
        delay = delay_model(),
        cfr = cfr_model(),
        report_delay = report_delay_model(),
        traveller = traveller_volume_model(),
        exports_deaths_model = exports_deaths_delay_model,
        ascertainment = pooled_ascertainment_model(),
        incubation = incubation_model(),
        source_population::Real = ITURI_POPULATION,
        pre_start_deaths::Union{Missing, Integer} = 0)
    growth_state ~ to_submodel(growth, false)
    delay_state ~ to_submodel(delay, false)
    cfr_state ~ to_submodel(cfr, false)
    report_state ~ to_submodel(report_delay, false)
    asc_state ~ to_submodel(ascertainment, false)
    incubation_state ~ to_submodel(incubation, false)
    os = onset_rescale(incubation_state.dist, growth_state.r)

    travel_state ~ to_submodel(traveller, false)
    daily_travellers = travel_state.daily_travellers

    ## Infection→detection delay = incubation ⊕ onset-to-report,
    ## moment-matched to one Gamma; the export at-risk clock runs from
    ## infection, so the survival is 1 at age 0. Export deaths are also
    ## timed from infection: the death delay is incubation ⊕ onset-to-death
    ## (incubation appears in both, a slight accepted double-count of the
    ## shared incubation, better than omitting it on death entirely).
    f_det = combined_delay(incubation_state.dist, report_state.dist)
    f_death = combined_delay(incubation_state.dist, delay_state.dist)
    exports_deaths_state ~ to_submodel(
        exports_deaths_model(export_deaths_daily, growth_state,
            cfr_state.CFR, f_death, asc_state.p_uganda;
            pre_start_deaths = pre_start_deaths,
            f_det = f_det,
            daily_travellers = daily_travellers,
            source_population = source_population),
        false)

    onset_fraction := os
    cumulative_infections := growth_state.C_T
    cumulative_cases := growth_state.C_T * os
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

`samples_analysed` is the per-vintage analysed-count vector aligned with
`confirmed_offsets`; each confirmed vintage is observed as
`C_v ~ Binomial(A_v, p_pos_v)` conditional on its analysed denominator
`A_v` (data, not modelled), which removes the multiplicative ridge of the
old NegBinomial-increment confirmed likelihood (#163). When left empty it
defaults to the cumulative `tests_analysed` for every vintage, recovering
the single-total denominator.

`tests_analysed` is a single cumulative testing-volume count observed at
its own elapsed time `tests_offset` before the cut-off, so it stays
robust if lab reporting lags or stops before the case cut-off. It also
supplies the binomial denominator for any `missing` `samples_analysed`
entry. Per-test positivity is exposed as a derived quantity. Pass
`tests_analysed = missing` to drop the tested-volume NegBinomial.

`confirmed_analysed_impute` toggles the imputed-denominator experiment:
when set to the [`analysed_impute_model`](@ref) submodel, confirmed
vintages whose analysed count is `missing` (the early 18-22 May and late
29-31 May lab windows, where no national analysed total was published)
get a TIGHT partially-pooled log-random-walk denominator anchored to the
observed 23-28 May increments. This funnels and does not converge (see
[`analysed_impute_model`](@ref)), so it is off by default; the default
`nothing` keeps the model on the 23-28 May observed-denominator vintages.

`confirmed_queue` (default `false`) switches the confirmed stream to a
coherent laboratory-throughput queue that fits ALL confirmed vintages,
including the dark early 18-22 May and late 29-31 May windows that lack a
published analysed denominator. The queue drains the received backlog by a
capacity-limited rate `μ_A = backlog·(1 − exp(−κ·Δt/backlog))` for every
window. Observed-denominator windows condition on the real analysed count
(`confirmed ~ Binomial(ΔA_obs, p_pos)` and `ΔA_obs ~ Poisson(μ_A)`); dark
windows use the exact marginal of `Binomial(Poisson(μ_A), p_pos)`, namely
`confirmed ~ Poisson(μ_A · p_pos)`, so the unobserved denominator is
integrated out rather than carried as a free per-vintage latent. This
replaces the [`analysed_impute_model`](@ref) funnel. Pass
`confirmed_epi_exclusion = epi_exclusion_model()` to add the opt-in
epi-exclusion fraction `e ~ Beta(2, 12)` (received asymptotes to
`(1 − e)·N_susp`); the default `nothing` pins `e = 0` (forward fraction 1)
for the headline fit. `τ_forward` is dropped in this path.

`samples_received` is the per-vintage cumulative received-count vector
(`Cumul échantillons reçus`). When supplied it conditions the forwarded
fraction `τ_forward` via `R_v ~ NegBinomial(τ_forward · N_susp,v, k)`,
with `N_susp,v` the cumulative suspect backlog (BVD-suspected plus
background) the confirmed model already uses for the positivity baseline
(see [`confirmed_cases_model`](@ref)). This pins `τ_forward` directly from
received-versus-suspected; left empty it drops the received likelihood.

`confirmed_deaths` is an optional per-vintage laboratory-confirmed-death
increment vector (the sitrep front-page `Cumul décès parmi les
confirmés`, deaths that got confirmed) aligned with
`confirmed_death_offsets` (default `death_offsets`). Each increment is a
genuine lab/positivity process on the death specimens forwarded to the
laboratory (issue #193): the suspect-death backlog at the confirmed-death
edges (the BVD-death CFR-weighted convolution plus the constant-rate
non-BVD background `λ_bg_death`) is forwarded at fraction `τ_death`, and
its BVD share sets the death-specimen positivity `s·q_death +
(1−spec)(1−q_death)`. The PCR sensitivity `s` and specificity `spec` are
imported *shared* from the confirmed-case lab pipeline (see
[`confirmed_deaths_model`](@ref)), so the stream requires the
confirmed-case stream. Left empty (the default) the stream is off and
existing callers are unchanged.

`death_background` samples the constant-rate non-BVD suspected-death
background `λ_bg_death` (see [`death_background_model`](@ref)), the death
analogue of the case `λ_bg`, added to the BVD-driven suspected-death
expectation; pass `λ_bg_death_fixed = 0.0` to disable it. `death_forward`
samples the death-specimen forwarding fraction `τ_death` for the confirmed
deaths (see [`death_forward_model`](@ref)).

`confirmed_q_random_effect` is the per-vintage tested-BVD-share random
effect for the confirmed positivity (see [`confirmed_q_re_model`](@ref)),
on by default because the observed per-window positivity is non-monotone
and a monotone severe-first q-curve cannot match it. Pass `nothing` to
recover the smooth severe-first baseline.

`confirmed_positivity_link` chooses how the tested BVD share is set.
`:free` (default) samples the severe-first selection curve (`q0 → qinf`),
a free description of the tested share. `:composition` instead ties the
tested share to the suspect-pool composition `φ = μ_BVD/(μ_BVD+μ_bg)`
upsampled by the decaying severity enrichment
`confirmed_severity_enrichment` (see [`severity_enrichment_model`](@ref)
and [`confirmed_cases_model`](@ref)), so the confirmed/positivity data
identify the non-BVD background `λ_bg` rather than it being absorbed by a
free curve.

`deaths_ascertainment` samples a multiplicative drift factor `p_deaths`
on the expected-deaths trajectory (see
[`deaths_ascertainment_model`](@ref)); pass `p_deaths_fixed = 1.0` to
disable the factor entirely.

`report_onset_offset` sets the testing-onset clock for the severe-first
BVD-share decay: the tested BVD fraction relaxes from `q0` toward the
count-implied baseline over elapsed time since `t_report = T − offset`
(see [`confirmed_cases_model`](@ref) and [`report_onset_offset`](@ref)).
Pass `report_onset_offset(as_of_date)` (8 days for the 26 May cut-off);
the default `nothing` anchors the clock at seeding (`t_report = 0`).

The Uganda exports streams default to the explicit infection→detection
delay mechanism ([`exports_delay_model`](@ref),
[`exports_deaths_delay_model`](@ref),
[`exports_detection_timing_delay_model`](@ref)), which convolves the
infection trajectory with an infection→detection delay. The exports
stream is travel-gated, so the at-risk clock starts at infection (the
traveller moves during incubation, pre-symptomatic): the delay is the
incubation period (`incubation`) convolved with the DRC onset-to-report
delay `f_rep` (`report_delay`, the same draws the incubation and
reported-cases streams use), moment-matched to one Gamma via
[`combined_delay`](@ref). The dated cases are anchored on their Uganda
admittance/detection dates, so this onset-to-report delay is reinterpreted
as the onset-to-admittance (care-seeking) delay abroad: the export count
is a single datum that cannot identify its own delay, and sharing the
exact DRC draws lets the incubation and reported-cases streams pin it.
This is a limitation — the imports were active care-seekers, so their true
onset-to-admittance delay is plausibly shorter than the DRC
onset-to-report delay (see the analysis limitations). Because the
incubation period sits inside that delay the export likelihoods carry no
separate incubation rescale. The export-death timing is likewise keyed to
infection: its death delay is incubation ⊕ onset-to-death (again via
[`combined_delay`](@ref)), so deaths are not timed one incubation period
too early. Incubation therefore enters both the detection and death
delays, a slight accepted double-count of the shared incubation period.
By default the exports are a single cumulative count `exported_cases` at
the cut-off, with the earliest Uganda detection entering separately as the
one-sided timing bound `first_export_detection_delta`. Passing a non-empty
`exported_cases_daily` (earliest detection day to the cut-off, see
[`load_observations`](@ref)) switches to the time-resolved
[`exports_daily_delay_model`](@ref), which fits each import at its Uganda
report date. Its pre-detection survival term already carries the
first-detection timing bound, so the separate timing term is then disabled
to avoid double-counting the earliest detection. `export_last_offset`
(default 0) stops both travel-gated streams (exports and deaths-among-
exports) `export_last_offset` days before the cut-off — the date of the
most recent reported import to Uganda — since movement patterns likely
shift over the outbreak and the most recent days are right-truncated by
reporting lag.

The McCabe et al. rectangular detection-window configuration is provided
as a separate comparison by
[`imperial_only_model`](@ref) (exports and deaths, window-based); the
window submodels (`exports_model`, `exports_deaths_model`,
`exports_detection_timing_model`) take a window rather than a delay, so
they are not drop-in `bvd_joint` kwargs. Pass a `report_delay` (or
`incubation`) returning a fixed `Gamma` (with no `~` sampling) to pin the
export infection→detection delay rather than learning it.
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
        confirmed_deaths::AbstractVector = Union{Missing, Int}[],
        confirmed_death_offsets::AbstractVector = death_offsets,
        samples_analysed::AbstractVector = Union{Missing, Int}[],
        samples_received::AbstractVector = Union{Missing, Int}[],
        tests_analysed::Union{Missing, Integer} = missing,
        tests_offset::Real = 0,
        exported_cases_daily::AbstractVector = Union{Missing, Int}[],
        export_last_offset::Real = 0,
        growth = exponential_growth_model(),
        exports = exports_delay_model,
        exports_daily = exports_daily_delay_model,
        deaths = deaths_model,
        reported_cases_submodel = reported_cases_model,
        confirmed = confirmed_cases_model,
        confirmed_deaths_submodel = confirmed_deaths_model,
        exports_deaths_model = exports_deaths_delay_model,
        exports_detection_timing = exports_detection_timing_delay_model,
        dispersion = surveillance_dispersion_model(),
        reported_dispersion = nothing,
        ascertainment = pooled_ascertainment_model(),
        deaths_ascertainment = deaths_ascertainment_model(),
        p_deaths_fixed::Union{Nothing, Real} = nothing,
        death_background = death_background_model(),
        λ_bg_death_fixed::Union{Nothing, Real} = nothing,
        death_forward = death_forward_model(),
        test_positivity = test_positivity_model(),
        report_delay = report_delay_model(),
        test_sensitivity = test_sensitivity_model(),
        test_specificity = test_specificity_model(),
        test_selection = test_selection_model(),
        confirmed_severity_enrichment = severity_enrichment_model(),
        confirmed_positivity_link::Symbol = :free,
        confirmed_overdispersion = nothing,
        confirmed_q_random_effect = confirmed_q_re_model,
        confirmed_analysed_impute = nothing,
        confirmed_queue::Bool = false,
        confirmed_epi_exclusion = nothing,
        confirmed_selection_clock::Symbol = :time,
        confirmed_volume_scale::Real = 200.0,
        incubation = incubation_model(),
        genetic = nothing,
        source_population::Real = ITURI_POPULATION,
        pre_start_deaths::Union{Missing, Integer} = 0,
        pre_detection_exports::Union{Missing, Integer} = 0,
        first_export_detection_delta::Union{Missing, Real} = missing,
        report_onset_offset::Union{Nothing, Real} = nothing)
    growth_state ~ to_submodel(growth, false)
    if genetic !== nothing
        genetic_state ~ to_submodel(genetic(growth_state.T), false)
    end
    dispersion_state ~ to_submodel(dispersion, false)
    asc_state ~ to_submodel(ascertainment, false)
    incubation_state ~ to_submodel(incubation, false)
    k = dispersion_state.k
    p_uganda = asc_state.p_uganda
    if p_deaths_fixed === nothing
        deaths_asc_state ~ to_submodel(deaths_ascertainment, false)
        p_deaths = deaths_asc_state.p_deaths
    else
        p_deaths = p_deaths_fixed
    end
    ## Constant-rate non-BVD background for the suspected deaths, the death
    ## analogue of the case `λ_bg`. Sampled here (like `p_deaths`) and passed
    ## into the deaths submodel; `λ_bg_death_fixed = 0.0` disables it,
    ## recovering the pure BVD suspected-death likelihood.
    if λ_bg_death_fixed === nothing
        death_bg_state ~ to_submodel(death_background, false)
        λ_bg_death = death_bg_state.λ_bg_death
    else
        λ_bg_death = λ_bg_death_fixed
    end
    T = growth_state.T
    ## The latent trajectory `exp(r·s)` is cumulative infections; the
    ## incubation mgf at −r maps it onto symptom onsets, the series every
    ## onset-driven stream below conditions on. With the default delay
    ## mechanism the exports and export-detection-timing streams also see
    ## onsets through the onset-to-detection delay convolution, so they
    ## are onset-rescaled by `os` too. The exports streams are sampled
    ## after the reported-cases stream so they can reuse its `f_rep`.
    os = onset_rescale(incubation_state.dist, growth_state.r)

    death_edges = [T - δ for δ in death_offsets]
    deaths_state ~ to_submodel(
        deaths(total_deaths, growth_state, k, death_edges;
            p_deaths = p_deaths, λ_bg_death = λ_bg_death,
            onset_fraction = os), false)

    ## Fixed per-stream DRC ascertainment: every reported and confirmed
    ## vintage bin uses the pooled scalar `p_drc`, shared between the two
    ## streams.
    n_rep = length(reported_offsets)
    n_conf = length(confirmed_offsets)
    p_drc_per_bin = fill(asc_state.p_drc, max(n_rep, n_conf))

    ## Optional separate dispersion for the suspected (reported) stream, so
    ## the reported likelihood can be loosened (down-weighted) without
    ## touching the deaths / received dispersion `k`. When `nothing` the
    ## reported stream shares `k` (the original behaviour).
    if reported_dispersion === nothing
        k_rep = k
    else
        reported_dispersion_state ~ to_submodel(reported_dispersion, false)
        k_rep = reported_dispersion_state.k
    end

    ## In the queue path forwarding is governed by `1 − e` from the
    ## exclusion submodel, so the `τ_forward` dimension is dead. Rebuild
    ## the test-positivity submodel with `sample_forward = false`,
    ## preserving any caller-supplied priors, so it is not sampled.
    test_positivity_eff = if confirmed_queue
        d = test_positivity.defaults
        test_positivity_model(;
            lambda_prior = d.lambda_prior,
            fraction_forwarded_prior = d.fraction_forwarded_prior,
            sample_forward = false)
    else
        test_positivity
    end

    reported_edges = [T - δ for δ in reported_offsets]
    reported_state ~ to_submodel(
        reported_cases_submodel(reported_cases, growth_state, k_rep,
            p_drc_per_bin[1:n_rep], reported_edges;
            report_delay = report_delay,
            test_positivity = test_positivity_eff,
            report_onset_offset = report_onset_offset,
            onset_fraction = os), false)

    if !isempty(confirmed_cases)
        confirmed_edges = [T - δ for δ in confirmed_offsets]
        tests_edge = T - tests_offset
        ## Per-vintage analysed denominators for the confirmed Binomial.
        ## Default to the cumulative `tests_analysed` for every vintage
        ## when no per-vintage denominators are supplied, so a single
        ## confirmed total still conditions on its analysed count.
        analysed_vec = isempty(samples_analysed) ?
                       fill(tests_analysed, n_conf) :
                       samples_analysed
        ## Per-vintage received counts (cumulative `Cumul échantillons
        ## reçus`). When supplied they condition the forwarded fraction
        ## `τ_forward` via a NegBinomial on the suspect backlog; an empty
        ## vector drops the received likelihood.
        received_vec = isempty(samples_received) ?
                       fill(missing, n_conf) : samples_received
        confirmed_state ~ to_submodel(
            confirmed(confirmed_cases, analysed_vec, received_vec,
                tests_analysed, growth_state, k,
                p_drc_per_bin[1:n_conf], reported_state.λ_bg,
                reported_state.τ_forward,
                reported_state.report_delay_dist,
                confirmed_edges, tests_edge;
                test_sensitivity = test_sensitivity,
                test_specificity = test_specificity,
                test_selection = test_selection,
                severity_enrichment = confirmed_severity_enrichment,
                positivity_link = confirmed_positivity_link,
                overdispersion = confirmed_overdispersion,
                q_random_effect = confirmed_q_random_effect,
                analysed_impute = confirmed_analysed_impute,
                confirmed_queue = confirmed_queue,
                epi_exclusion = confirmed_epi_exclusion,
                selection_clock = confirmed_selection_clock,
                volume_scale = confirmed_volume_scale,
                report_onset_offset = report_onset_offset,
                onset_fraction = os), false)
    end

    if !isempty(confirmed_deaths)
        ## Laboratory-confirmed deaths (`Cumul décès parmi les confirmés`):
        ## a genuine lab/positivity process on the post-mortem death
        ## specimens forwarded to the laboratory (issue #193), importing the
        ## SHARED case-lab sensitivity `s` and specificity `spec`. The
        ## suspect-death backlog at the confirmed-death edges is the BVD-death
        ## CFR-weighted convolution (sharing the deaths CFR / delay / growth /
        ## drift) plus the constant-rate non-BVD background `λ_bg_death`; the
        ## BVD share of that pool sets the death-specimen positivity and a
        ## forwarded fraction `τ_death` gives the expected confirmed count.
        isempty(confirmed_cases) &&
            error("confirmed_deaths requires confirmed_cases for the " *
                  "shared sensitivity `s` and specificity `spec`")
        cdeath_edges = [T - δ for δ in confirmed_death_offsets]
        bvd_death_edges = bvd_death_trajectory(growth_state.r,
            deaths_state.CFR, deaths_state.delay_dist, cdeath_edges;
            p_deaths = p_deaths, onset_fraction = os)
        nsusp_death_edges = [max(bvd_death_edges[i], zero(eltype(bvd_death_edges))) +
                             max(λ_bg_death * cdeath_edges[i],
                                 zero(λ_bg_death * cdeath_edges[i]))
                             for i in eachindex(cdeath_edges)]
        confirmed_deaths_state ~ to_submodel(
            confirmed_deaths_submodel(confirmed_deaths, bvd_death_edges,
                nsusp_death_edges, confirmed_state.s_test,
                confirmed_state.spec_test, k;
                death_forward = death_forward), false)
    end

    ## Exports are travel-gated, so the at-risk clock runs from infection:
    ## the detection delay is incubation ⊕ onset-to-report (the same draws
    ## the incubation and reported-cases streams use), moment-matched to
    ## one Gamma. No separate incubation rescale (it lives in `f_det`).
    f_det = combined_delay(
        incubation_state.dist, reported_state.report_delay_dist)
    ## Two export configurations. With a dated detection series the
    ## time-resolved likelihood uses each import's Uganda report date and
    ## its pre-detection survival term already carries the
    ## first-detection timing bound, so the separate timing submodel is
    ## made a no-op (`delta = missing`) to avoid double-counting the
    ## earliest detection. Otherwise the scalar single-total count is fit
    ## at the cut-off with the timing bound as a separate term.
    if isempty(exported_cases_daily)
        exports_state ~ to_submodel(
            exports(exported_cases, growth_state, p_uganda, f_det), false)
        detection_timing_delta = first_export_detection_delta
    else
        exports_state ~ to_submodel(
            exports_daily(exported_cases_daily, growth_state, p_uganda,
                f_det; last_offset = export_last_offset), false)
        detection_timing_delta = missing
    end

    ## Export deaths are also timed from infection: the death delay is
    ## incubation ⊕ onset-to-death, moment-matched to one Gamma, so the
    ## death timing carries the incubation period rather than firing one
    ## incubation period too early. (Detection and death share the same
    ## onset, so incubation appears in both `f_det` and `f_death`; this
    ## slight double-count of the shared incubation is accepted as better
    ## than omitting incubation on death entirely.)
    f_death = combined_delay(incubation_state.dist, deaths_state.delay_dist)
    exports_deaths_state ~ to_submodel(
        exports_deaths_model(export_deaths_daily, growth_state,
            deaths_state.CFR, f_death, p_uganda;
            pre_start_deaths = pre_start_deaths,
            last_offset = export_last_offset,
            f_det = exports_state.f_det,
            daily_travellers = exports_state.daily_travellers,
            source_population = source_population),
        false)
    detection_timing_state ~ to_submodel(
        exports_detection_timing(growth_state, p_uganda;
            delta = detection_timing_delta,
            pre_detection_exports = pre_detection_exports,
            f_det = exports_state.f_det,
            daily_travellers = exports_state.daily_travellers,
            source_population = source_population),
        false)

    onset_fraction := os
    cumulative_infections := growth_state.C_T
    cumulative_cases := growth_state.C_T * os
end

"""
McCabe et al. reimplementation composer: exports and deaths only,
mirroring the joint configuration in the Imperial report (no
reported-cases or deaths-among-exports likelihood). Passing
`missing` for `exported_cases` reduces to a pure Method 2 (deaths-
only) fit.

Unlike the other composers this fits no incubation period: McCabe et al.
model cases directly, so `C_T` is the case trajectory (`onset_fraction = 1`)
and `cumulative_infections` and `cumulative_cases` coincide. This keeps
`cumulative_cases` comparable with the report's published case estimates.

This composer is left entirely McCabe-based as the comparison / revert
path: it retains the rectangular detection-window assumption of
[`exports_model`](@ref) and [`detection_window_model`](@ref), in contrast
with the onset-to-detection delay convolution used by default in
[`bvd_joint`](@ref) and [`exports_only_model`](@ref).
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

    ## No incubation layer: McCabe et al. model cases directly, so the
    ## onset-to-death convolution acts on `C_T` as the case trajectory
    ## (`onset_fraction = 1`). `cumulative_infections` and `cumulative_cases`
    ## therefore coincide here, both equal to `C_T`.
    cumulative_infections := growth_state.C_T
    cumulative_cases := growth_state.C_T
end
