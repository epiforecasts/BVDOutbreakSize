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


## --- Patch (multi-population) joint models ------------------------------

"""
Patch latent process: run [`patch_infection_model`](@ref) and expose the
per-patch state alongside the national aggregates (summed onsets and
infections) that the national observation submodels consume.
"""
@model function _patch_latent(n::Integer, n_patches::Integer,
        breakpoint, patch_infection;
        rt_start::Integer = 1,
        rt_walk_start::Integer = rt_start,
        importation_kernel::AbstractMatrix = province_importation_kernel(
            PROVINCE_POPULATIONS[1:min(n_patches, end)]))
    patch_state ~ to_submodel(
        patch_infection(n, n_patches;
            breakpoint, rt_start, rt_walk_start,
            importation_kernel), false)
    onsets_total = vec(sum(patch_state.onsets_matrix; dims = 1))
    return (; patch_state, onsets_total)
end

"""
Modelled per-province confirmed-DEATH increments, binned to the vintages of
the Tableau 1 spatial tables. Each patch's onsets are pushed through the
onset-to-death delay and then the report-to-receipt delay (so `kernel` is
the convolution of the two), giving the deaths CONFIRMED by each vintage.

This is the term that makes the provincial split identifiable. The
composition of DEATHS across provinces weights each patch by its
delay-convolved incidence and by nothing else: the case-fatality ratio and
the death-confirmation probability are properties of the virus and of a
national laboratory pipeline, not of a province, so they are common factors
and cancel out of the normalised shares. What remains is the provincial
INCIDENCE split — free of case-ascertainment.

Convolving through the onset-to-death delay (rather than comparing raw
death-to-case ratios) is what separates ascertainment from epidemic phase. A
fast-growing province has proportionally fewer deaths *to date* than a flat
one at the same true CFR, simply because its recent cases have not died yet.
Ituri grows faster than Nord-Kivu, so a naive CFR comparison would read that
right-censoring as a difference in case-finding. The delay convolution
predicts each province's deaths-to-date from its own incidence curve, so the
censoring is accounted for and only the residual is ascertainment.
"""
function _patch_death_increments(onsets_matrix::AbstractMatrix,
        kernel::AbstractVector,
        province_days::AbstractVector{<:Integer})
    np = size(onsets_matrix, 1)
    nv = length(province_days)
    first_daily = convolve_delay(vec(@view onsets_matrix[1, :]), kernel)
    out = Matrix{eltype(first_daily)}(undef, np, nv)
    @inbounds out[1, :] = bin_increments(first_daily, province_days)
    @inbounds for p in 2:np
        daily = convolve_delay(vec(@view onsets_matrix[p, :]), kernel)
        out[p, :] = bin_increments(daily, province_days)
    end
    return out
end

"""
Modelled per-province confirmed increments, binned to the vintages of the
per-province spatial tables. Each patch's onsets are pushed through the
SAME report-to-receipt delay and test sensitivity as the national confirmed
stream — the laboratory pipeline is national, only the incidence feeding it
is provincial — and then binned onto the shared vintage days.

`province_days` are the (shared) grid-day indices of the spatial-table
vintages. Returns an `(n_patches × n_vintages)` matrix.

!!! note "Ascertainment and incidence are confounded by province"
    This matrix reaches [`province_composition_model`](@ref) only through
    NORMALISED shares, so every factor common to all provinces (`s_test`
    here; ascertainment, background and positivity in the national confirmed
    stream) cancels. That is what keeps the composition from re-scoring the
    national total.

    What it does NOT do is fix the case-finding probability across provinces.
    The composition weights each patch by `asc_p * lambda_p` — its relative
    ascertainment times its modelled incidence — and the data identify only
    that PRODUCT. A province with fewer infections but better case-finding
    looks exactly like one with more infections and worse case-finding.

    That is a property of the data, not of this code. The per-province
    laboratory series shows the provinces are testing very
    differently-selected pools (Ituri 31.8% positivity against Nord-Kivu's
    5.5%), so assuming equal ascertainment would push the whole difference
    into the provincial `Rt` and report a case-finding artefact as
    epidemiology. `asc_p` is therefore sampled and partially pooled toward
    equality rather than fixed, which widens the per-patch `Rt` contrast and
    `C_T` split to their honest width.

    Read `log_rt_contrast` and `province_ascertainment` together; neither is
    interpretable alone. The NATIONAL headline does not depend on the split
    and is unaffected by the confound.
"""
function _patch_confirmed_increments(onsets_matrix::AbstractMatrix,
        receipt_pmf::AbstractVector, s_test::Real,
        province_days::AbstractVector{<:Integer})
    np = size(onsets_matrix, 1)
    nv = length(province_days)
    first_daily = s_test .* convolve_delay(
        vec(@view onsets_matrix[1, :]), receipt_pmf)
    out = Matrix{eltype(first_daily)}(undef, np, nv)
    @inbounds out[1, :] = bin_increments(first_daily, province_days)
    @inbounds for p in 2:np
        daily = s_test .* convolve_delay(
            vec(@view onsets_matrix[p, :]), receipt_pmf)
        out[p, :] = bin_increments(daily, province_days)
    end
    return out
end

"""
Joint composer over all data streams using the PATCH (meta-population)
latent process.

Runs [`patch_infection_model`](@ref) over `n_patches` provinces, sums the
per-patch onsets into a national trajectory, and conditions on exactly the
national observation submodels that [`bvd_joint`](@ref) uses: DRC suspected
cases, deaths, confirmed cases and deaths, the laboratory pipeline,
treatment flows and Uganda exports. The national streams therefore see the
same latent onset trajectory they would in the single-patch model.

The spatial information enters through one extra term:
[`province_composition_model`](@ref), which scores how the confirmed cases
divide between provinces, conditional on the national total. The
per-province spatial-table counts are an exact partition of the national
confirmed counts, so scoring them with their own count likelihood would put
the same data into the joint density twice; the composition term adds only
the spatial signal that the national series does not carry.

Pass that data as `province_increments` (an `n_patches × n_vintages` matrix
of new-confirmed counts) and `province_days` (the shared vintage day
indices), which [`province_increment_matrix`](@ref) builds from the
per-province histories that [`load_observations`](@ref) returns. Omit them
and the composition term is skipped, leaving a purely national fit over a
patch latent process.

Uganda exports are driven by the primary patch (Ituri) alone, since the
border crossings the export stream describes are Ituri-to-Uganda.

Headline quantities (`C_T`, `R_T`, `r`, `T`, `doubling_time`, `R0`) are
surfaced under the same names as [`bvd_joint`](@ref), so a patch chain
summarises through [`summary_table`](@ref) exactly like a single-patch one.
`R_T` is the incidence-weighted aggregate reproduction number implied by
the summed patch infections. Per-patch quantities (`C_T_patch`, `R_T_patch`,
`delta_patch`) are surfaced as vector deterministics for
[`patch_summary_table`](@ref).
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
            PROVINCE_POPULATIONS[1:min(n_patches, end)]),
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
        confirmed_break_days::AbstractVector{<:Integer} = Int[],
        confirmed_break_gross_cases::AbstractVector{<:Integer} = Int[],
        confirmed_break_gross_deaths::AbstractVector{<:Integer} = Int[],
        confirmed_break_sd::Real = 25.0,
        treatment_confirmed_incare_history = (; days = Int[], counts = Int[]),
        treatment_suspect_incare_history = (; days = Int[], counts = Int[]),
        breakpoint::Union{Missing, Real} = missing,
        source_population::Real = ITURI_POPULATION,
        patch_infection = patch_infection_model,
        composition = province_composition_model,
        province_increments::Union{Missing, AbstractMatrix{<:Integer}} = missing,
        province_days::AbstractVector{<:Integer} = Int[],
        province_death_increments::Union{
            Missing, AbstractMatrix{<:Integer}} = missing,
        province_death_days::AbstractVector{<:Integer} = Int[],
        death_composition = province_composition_model,
        death_ascertainment_sd_prior = truncated(
            Normal(0, 0.1); lower = 0),
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
        genetic = nothing,
        confirmed_positivity_link::Symbol = :composition,
        onset_to_sample = nejm_onset_to_sample(),
        tmrca_days::Union{Missing, Real} = missing,
        tmrca_days_sd::Real = 15.0,
        renewal_start_lead::Integer = RENEWAL_START_LEAD,
        rt_walk_lead::Integer = RT_WALK_LEAD)
    ## Guard against the silent failure mode: per-province data supplied but
    ## `n_patches` left at its default of 1. The compositions would then be
    ## scored against a single patch that holds the entire national total, the
    ## spatial structure would quietly vanish, and the fit would look fine.
    if n_patches == 1 &&
       (!isempty(province_days) || !isempty(province_death_days))
        error("per-province data was supplied but n_patches = 1. The spatial " *
              "structure would be silently dropped. Pass n_patches = " *
              "$(length(PROVINCE_NAMES)) (or the number of patches the data " *
              "covers).")
    end
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
    ## Background random effect, guarded against an empty reported history.
    ## The lead is the max lag of the report-to-receipt kernel.
    case_bg_re = background_re ?
    begin
        bg_pool ~ to_submodel(background_pooling_model())
        σ_rw_shared = bg_pool.σ_bg
        bg_lead = 7
        bg_onset = isempty(reported_history.days) ? 1 :
                   clamp(Int(reported_history.days[1]) - bg_lead, 1, n)
        nn -> background_walk_model(nn, σ_rw_shared; onset = bg_onset)
    end : nothing
    ## --- National observation submodels (identical to bvd_joint) --------
    ## 1. Reported (suspected) cases.
    cases_state ~ to_submodel(cases(reported_history, reported_cases,
        onsets, k_cases, p_drc;
        suspected_daily_history, background_re = case_bg_re))
    ## 2. Deaths (suspected).
    deaths_state ~ to_submodel(deaths(deaths_history, total_deaths, onsets,
        k_deaths;
        suspected_daily_deaths_history, case_bg_daily = cases_state.bg_daily))
    ## 3. Confirmed cases (laboratory pipeline). This scores the national
    ##    confirmed TOTAL; the provincial split is scored separately below.
    confirmed_state ~ to_submodel(confirmed(confirmed_history,
        confirmed_cases, onsets, k_confirmed, p_drc,
        cases_state.bg_daily, cases_state.τ_test,
        cases_state.bvd_reports_daily;
        lab_history, lab_daily_history, tests_analysed,
        break_days = confirmed_break_days,
        break_gross_cases = confirmed_break_gross_cases,
        break_sd = confirmed_break_sd,
        positivity_link = confirmed_positivity_link))
    ## 4. Confirmed deaths.
    confirmed_deaths_state ~ to_submodel(
        confirmed_deaths_stream(confirmed_deaths, total_deaths,
        deaths_state.deaths_daily, deaths_state.bvd_deaths_daily,
        deaths_state.bg_death_daily, k_confirmed_deaths;
        confirmed_deaths_history, receipt_pmf = confirmed_state.receipt_pmf,
        case_analysed_daily = confirmed_state.analysed_daily,
        case_suspected_daily = cases_state.reports_daily,
        break_days = confirmed_break_days,
        break_gross_deaths = confirmed_break_gross_deaths))
    ## 5. Uganda exports, from the primary patch (Ituri) only: the border
    ##    crossings this stream describes are Ituri-to-Uganda.
    exports_state ~ to_submodel(exports(exported_cases,
        vec(patch_state.infections_matrix[1, :]), p_uganda;
        export_case_days, incubation_pmf = patch_state.incubation_pmf,
        source_population))
    exports_deaths_state ~ to_submodel(exports_deaths_model(
        exports_deaths, exports_state.travelled_prevalence,
        deaths_state.CFR, deaths_state.od_pmf,
        patch_state.incubation_pmf;
        export_death_days))
    ## Genetic seeding bound on the total outbreak age, as in the national
    ## composer: the molecular-clock TMRCA pulls the origin to sit at or before
    ## the most recent common ancestor.
    if genetic !== nothing
        genetic_state ~ to_submodel(
            genetic(patch_state.T, tmrca_days; tmrca_days_sd), false)
    end

    ## 6. Treatment flows (isolation, bed capacity, LOS), on the national
    ##    suspected-case inflow, matching the national-level data.
    treatment_state ~ to_submodel(treatment(isolation_history,
        cases_state.bvd_reports_daily, cases_state.bg_daily, p_drc,
        deaths_state.CFR;
        capacity_history = bed_capacity_history,
        admissions_history = treatment_admissions_history,
        deaths_history = treatment_deaths_history,
        ruleout_history = treatment_ruleout_history,
        absconded_history = treatment_absconded_history,
        confirmed_incare_history = treatment_confirmed_incare_history,
        suspect_incare_history = treatment_suspect_incare_history,
        occupancy_break_days = occupancy_break_days,
        k_external = k_isolation))
    ## 7. Recovered (among confirmed).
    recovered_state ~ to_submodel(recovered(recovered_history,
        recovered_cases, confirmed_state.confirmed_daily,
        deaths_state.CFR; k_external = k_recovered))
    ## --- Spatial observation submodel -----------------------------------
    ## 8. Per-province composition of the confirmed cases. Conditional on the
    ##    national total (already scored above), so no observation is counted
    ##    twice. Skipped when no spatial-table data is supplied.
    ##
    ##    `province_increments` and `province_days` are the already-reshaped
    ##    spatial-table data (see [`province_increment_matrix`](@ref)), NOT the
    ##    raw per-province history dict. The reshaping is pure data handling
    ##    with no dependence on any parameter, and it looks provinces up by
    ##    name in a `Dict{String}`; doing that inside the model body puts a
    ##    string comparison (`memcmp`) on the AD tape, which Mooncake has no
    ##    rule for and which aborts the gradient. It must stay hoisted out.
    if !isempty(province_days)
        modelled_prov = _patch_confirmed_increments(
            patch_state.onsets_matrix, confirmed_state.receipt_pmf,
            confirmed_state.s_test, province_days)
        composition_state ~ to_submodel(
            composition(province_increments, modelled_prov))
        province_shares := composition_state.shares
        province_composition_rho := composition_state.rho
        ## Relative province CASE ascertainment: the probability an infection
        ## there becomes a CONFIRMED case, partially pooled and sum-to-zero on
        ## the log scale. On its own the case composition identifies only the
        ## PRODUCT of ascertainment and incidence. The DEATH composition below
        ## is what separates them.
        province_ascertainment := composition_state.province_ascertainment
        province_ascertainment_sd := composition_state.ascertainment_sd
    end
    ## 9. Per-province composition of the confirmed DEATHS. This is the term
    ##    that identifies the provincial split.
    ##
    ##    The case composition weights each patch by `asc_p * lambda_p` and can
    ##    never separate the two: a province with fewer infections but better
    ##    case-finding is observationally identical to one with more infections
    ##    and worse case-finding.
    ##
    ##    Deaths break the tie. They are far harder to miss than cases, and the
    ##    case-fatality ratio and the death-confirmation probability belong to
    ##    the virus and to a national laboratory pipeline, not to a province —
    ##    so they are common factors and cancel out of the normalised death
    ##    shares. What is left weights each patch by its delay-convolved
    ##    INCIDENCE alone, free of case-ascertainment. The deaths therefore pin
    ##    `lambda_p`, and the case composition then identifies `asc_p` as the
    ##    residual.
    ##
    ##    `death_ascertainment_sd_prior` is deliberately TIGHT: the identifying
    ##    assumption is that death ascertainment is near-uniform across
    ##    provinces, which is far weaker and more defensible than assuming case
    ##    ascertainment is. It is not fixed at zero, so the assumption can bend
    ##    where the data insist rather than snapping.
    ##
    ##    The data say this matters: Nord-Kivu holds a steady 8-9% of confirmed
    ##    cases but 14-19% of confirmed deaths at every vintage.
    if !isempty(province_death_days)
        death_kernel = convolve_pmf(
            deaths_state.od_pmf, confirmed_state.receipt_pmf)
        modelled_deaths_prov = _patch_death_increments(
            patch_state.onsets_matrix, death_kernel, province_death_days)
        death_composition_state ~ to_submodel(
            death_composition(province_death_increments, modelled_deaths_prov;
            ascertainment_sd_prior = death_ascertainment_sd_prior))
        province_death_shares := death_composition_state.shares
        province_death_ascertainment := death_composition_state.province_ascertainment
    end
    ## --- Deterministics surfaced for reporting --------------------------
    ## Headline quantities, under the same names as bvd_joint so a patch
    ## chain drops into the existing summary and forecast machinery.
    R0 := patch_state.R0
    r := patch_state.r
    r0 := patch_state.r0
    doubling_time := patch_state.doubling_time
    T := patch_state.T
    C_T := patch_state.C_T
    R_T := patch_state.R_T
    expected_infections_T := @inbounds(patch_state.infections_total[n])
    ## The national cumulative infection trajectory, summed over patches. The
    ## single-population composer surfaced this from its `_latent` submodel
    ## rather than from its own body, so a diff of the two function bodies does
    ## not show it as missing -- which is exactly how it was dropped. The docs
    ## read it off the chain for the headline ribbon panels.
    cumulative_infections := patch_state.cumulative_total
    cumulative_onsets := cumsum(onsets)
    Rt_national_implied := patch_state.Rt_national_implied
    ## Per-patch quantities, as vector deterministics (one entry per patch).
    C_T_patch := patch_state.C_T_patch
    R_T_patch := [@inbounds(patch_state.Rt_matrix[p, n]) for p in 1:n_patches]
    infections_T_patch := [@inbounds(patch_state.infections_matrix[p, n])
                           for p in 1:n_patches]
    ## The per-patch log-Rt modifier at the cut-off (one entry per patch),
    ## and its spread at the START of the walk, so a change in the provincial
    ## Rt gap over the window is visible as the difference between them.
    delta_patch := [@inbounds(patch_state.δ_patch[p, n]) for p in 1:n_patches]
    delta_patch_start := [@inbounds(patch_state.δ_patch[p, rt_walk_start])
                          for p in 1:n_patches]
    ## THE spatial diagnostic: the per-patch scale of the log-Rt deviation
    ## walk. A posterior concentrated near zero says the provinces share one
    ## temporal Rt shape (a fixed ratio between them); pushed away from zero
    ## it is direct evidence that provincial Rt trajectories are separating.
    ## See [`patch_rt_model`](@ref).
    region_sd := patch_state.σ_level
    region_drift_sd := patch_state.σ_δ
    ## `seed_fraction` (each secondary patch's seed as a fraction of the
    ## primary patch's) is sampled inside the latent submodel and already
    ## reaches the chain under that name, so it is NOT re-surfaced here. Read
    ## it alongside `log_rt_contrast`: the seed fraction sets the LEVEL of the
    ## provincial case split, and if it were pinned far below what the data
    ## need, the Rt contrast would silently absorb the difference and the
    ## provincial Rt gap would be an artefact of the seed prior.
    ## Learned cross-patch correlation of the deviation innovations. With
    ## three patches only the Ituri / Nord-Kivu entry carries real
    ## information (Sud-Kivu has no signal), so the rest tracks the LKJ prior.
    ## With a single patch there is no cross-patch correlation to report; the
    ## 1x1 correlation matrix is trivially 1.
    region_corr_primary_secondary := n_patches > 1 ?
                                     @inbounds(patch_state.Ω[1, 2]) :
                                     one(eltype(patch_state.Ω))
    ## The provincial log-Rt CONTRASTS at the cut-off: what the per-province
    ## composition data actually measure. Entry p is log R_p - log R_1, so a
    ## negative value means province p is transmitting less than the primary
    ## patch. Sum-to-zero deviations make these the interpretable quantity
    ## rather than the deviations themselves.
    log_rt_contrast := [@inbounds(patch_state.δ_patch[p, n] -
                                  patch_state.δ_patch[1, n])
                        for p in 1:n_patches]
    ## --- The full bvd_joint deterministic set --------------------------
    ## The patch model is the headline joint, so a patch chain must carry
    ## EVERY quantity a single-patch chain does: `analysis.jl`, the forecast
    ## machinery (`forecast_reported` reads `expected_reports_T`,
    ## `cumulative_confirmed`, ...) and the plots all key off these names. The
    ## observation submodels are identical to bvd_joint's, so the states are
    ## already here; only the latent-derived quantities differ, and those are
    ## surfaced above from the patch state.
    cumulative_expected_deaths := cumsum(deaths_state.bvd_deaths_daily)
    ## Modelled daily laboratory-confirmed cases (from `confirmed_cases_model`:
    ## the per-window tested-positive probability applied to the modelled,
    ## testing-onset-gated analysed volume), so the cumulative trajectory carries
    ## the confirmed-case timing for the delay-corrected confirmed-CFR
    ## reconstruction. The onset-to-confirmation kernel (onset-to-report ⊕
    ## receipt) and the onset-to-death-confirmation kernel (onset-to-death ⊕
    ## receipt) are exposed alongside so the residual delay between a confirmed
    ## case and its confirmed death can be rebuilt per draw off the chain.
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
    onset_to_confirmation_pmf := convolve_pmf(cases_state.report_pmf, confirmed_state.receipt_pmf)
    onset_to_death_confirmation_pmf := convolve_pmf(deaths_state.od_pmf, confirmed_state.receipt_pmf)
    ## External onset-to-sample constraint on the confirmed sampling delay
    ## (grounded on the NEJM DRC 2026 cohort by default, see
    ## [`nejm_onset_to_sample`](@ref)). The onset→report and report→receipt legs
    ## convolve to the confirmed onset-to-sample delay, so its continuous mean is
    ## the sum of the two legs' means and its continuous SD the root-sum of their
    ## variances; both are exposed here. The cohort's reported (continuous) mean
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
    ## Population-level dispersion (`k`, the headline scalar) plus the
    ## partially-pooled per-stream dispersions and the pooling SD.
    expected_deaths_T := deaths_state.expected_deaths_T
    expected_reports_T := cases_state.expected_reports
    expected_confirmed_T := confirmed_state.expected_confirmed
    expected_analysed_T := confirmed_state.expected_analysed
    _ecd = confirmed_deaths_state.expected_confirmed_deaths
    expected_confirmed_deaths_T := _ecd
    expected_exports_T := exports_state.expected_exports
    expected_exports_deaths_T := exports_deaths_state.expected_exports_deaths_T
    expected_isolation_T := treatment_state.expected_isolation
    expected_bed_demand_T := treatment_state.expected_bed_demand
    bed_shortfall_T := safe_rate(treatment_state.expected_bed_demand -
                                 treatment_state.expected_isolation)
    ## Cut-off daily treatment flows surfaced for the one-week-ahead forecast.
    expected_admissions_T := treatment_state.expected_admissions
    expected_incare_deaths_T := treatment_state.expected_incare_deaths
    expected_ruleouts_T := treatment_state.expected_ruleouts
    bed_capacity := treatment_state.capacity
    isolation_admission := treatment_state.p_iso
    isolation_bvd_admission := treatment_state.p_iso_bvd
    isolation_severity := treatment_state.δ_iso
    ## BVD bed stay is now the outcome mixture; `isolation_bvd_los_mean`
    ## reports the mixture mean (overall length-of-stay), with the death and
    ## recovery branch means surfaced separately.
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

    ## Shared observation-model parameters.
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
