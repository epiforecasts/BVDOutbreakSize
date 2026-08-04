# Observation loading from the dated TOML manifest. The manifest stores
# calendar dates and cumulative counts (never grid day-indices), so the
# data stay aligned by date as vintages are added, revised or arrive
# sparsely. The renewal model works on a daily grid, so this loader
# derives the grid length and the per-vintage day-indices from the dates
# at load time: the cut-off is the last grid day, and the seeding day (the
# first grid day) is placed a fixed lead before the genetic TMRCA so the
# grid always contains the inferred seeding time.

## Days the grid extends before the genetic TMRCA date, so the seeding
## crossing has room to be inferred below the molecular-clock bound.
const SEEDING_LEAD_DAYS = 30

"""
Load the BVD observation manifest from `path` (a dated TOML file) and
return a named tuple for the renewal model. Calendar dates are converted
to 1-based grid day-indices (day 1 is the seeding day, day `n` the
cut-off) using the cut-off (`as_of_date`) and a seeding day placed
`seeding_lead` days before the genetic TMRCA date.

Returns the grid length `n`, the `cutoff` and `seeding` dates, and the
per-stream cumulative totals at the cut-off (`reported_cases`,
`total_deaths`, `confirmed_cases`, `confirmed_deaths`, `tests_analysed`,
`exported_cases`, `exports_deaths`).

The dated Uganda export series are grid day-indices (`export_case_days`,
`export_death_days`), each a sorted list of detection/death days on or
before the cut-off, fitted with a per-day Poisson likelihood.

The per-vintage histories are returned as `(; days, counts)` with `days`
the grid day-indices: `reported_history`, `confirmed_history`,
`confirmed_deaths_history`, `deaths_history`, `lab_history` (the
cumulative analysed-specimen series), `lab_daily_history` (post-cutoff
24h analysed counts on the trusted post-cutoff days),
`suspected_daily_history` (post-cutoff daily new-suspect inflow,
"nouveaux cas suspects du jour"), `suspected_daily_deaths_history`
(post-cutoff daily new suspected deaths, "cas suspects du jour N (M
deces)", the deaths analogue), `isolation_history` (post-cutoff daily
isolation/hospitalisation occupancy, "Patients en isolement", a daily
bed count fitted by the length-of-stay submodel), `bed_capacity_history`
(the implied bed capacity, occupancy / reported occupancy rate),
`recovered_history` (the cumulative recovered-among-confirmed series,
"cumul guéris"), `treatment_confirmed_incare_history` and
`treatment_suspect_incare_history` (the Tableau 6 occupancy split, `dont
confirmes (NC+AC)` and `dont suspects`, two prevalence sub-stocks that
sum to the total occupancy), and `tests_received_history`.

Also returned: the genetic TMRCA bound `tmrca_days` (days before the
cut-off), and `who_first_sitrep_days` (days from the first situation
report, the earliest reported-case vintage, to the cut-off). The
intervention breakpoint grid day is `n - who_first_sitrep_days`. A
cut-off scalar with no explicit TOML block is derived from the final
vintage of the matching history.
"""
function load_observations(
        path::AbstractString = joinpath(@__DIR__, "..", "data",
            "observations.toml");
        seeding_lead::Integer = SEEDING_LEAD_DAYS,
        cutoff_date::Union{Nothing, Date, AbstractString} = nothing)
    raw = TOML.parsefile(path)
    _val(k) = raw[k]["value"]
    ## The cut-off is the manifest `as_of_date` unless an earlier
    ## `cutoff_date` is supplied (used to freeze the data to a past
    ## date, see `freeze_observations`). Freezing only ever moves the
    ## cut-off earlier and never invents vintages, so the grid stays
    ## date-aligned with the full-data fit.
    cutoff = isnothing(cutoff_date) ? Date(String(raw["as_of_date"])) :
             (cutoff_date isa Date ? cutoff_date :
              Date(String(cutoff_date)))
    ## Grid day-index (1-based) of a calendar date: seeding day is day 1.
    _gap(d) = Int(date2epochdays(cutoff) - date2epochdays(Date(String(d))))
    tmrca_date = Date(String(raw["genetic_tmrca"]["date"]))
    seeding = tmrca_date - Day(seeding_lead)
    n = Int(date2epochdays(cutoff) - date2epochdays(seeding)) + 1
    _index(d) = n - _gap(d)

    ## A dated cumulative history → grid day-indices and counts, sorted
    ## oldest-first so the model differences consecutive vintages into
    ## daily increments. Empty when the block is absent. Vintages dated
    ## after the cut-off are dropped, so freezing to an earlier date
    ## keeps only the data that was available by then.
    function history(key)
        haskey(raw, key) || return (; days = Int[], counts = Int[])
        block = raw[key]
        keep = [Date(String(d)) <= cutoff for d in block["dates"]]
        idx = Int[_index(d) for d in block["dates"][keep]]
        vals = Int.(block["values"][keep])
        ord = sortperm(idx)
        return (; days = idx[ord], counts = vals[ord])
    end

    ## A dated list of event dates (not a cumulative block) → the grid
    ## day-indices on or before the cut-off, sorted ascending. Used for the
    ## dated Uganda export-case and export-death series, each a list of
    ## detection/death dates fitted with a per-day Poisson likelihood.
    ## Empty when the block is absent or every date falls after the cut-off.
    function event_days(key)
        haskey(raw, key) || return Int[]
        ds = String.(raw[key]["value"])
        keep = [Date(d) <= cutoff for d in ds]
        idx = Int[_index(d) for d in ds[keep]]
        return sort(idx)
    end

    ## Dated Uganda export-case detection days and export-death days
    ## (1-based grid indices). Each detection/death contributes one Poisson
    ## term at its day. Days carry the per-day expected export count
    ## differenced from the at-risk export person-time (see `exports_model`).
    export_case_days = event_days("export_case_dates")
    export_death_days = event_days("export_death_dates")

    ## Manually specified occupancy reclassification-break days (opt-in): grid
    ## days on which the treatment-flow model fits a level step into the
    ## modelled occupancy mean, absorbing a between-report measurement-basis
    ## discontinuity in the observed isolation series (e.g. a Tableau 6-sum →
    ## page-1 headline transition) without bending Rt. A dated list filtered
    ## to the cut-off like the histories. Absent or empty → no break days, a
    ## no-op.
    occupancy_break_days = event_days("occupancy_break_dates")

    ## Manually specified retrospective harmonisation-break days (opt-in) for
    ## the confirmed streams: grid days on which INSP integrated a harmonised
    ## provincial base, so the cumulative confirmed headline steps by far more
    ## than that day's own notifications. On such a day the confirmed submodels
    ## de-anchor the laboratory positivity denominator (the reattached cases
    ## are not same-day positives) and fit a level step into the modelled
    ## confirmed mean, so the increment likelihood does not read the backlog as
    ## one day of incidence. A dated list filtered to the cut-off like the
    ## histories. Absent or empty → no break days, a no-op.
    ## The grid days are returned sorted while the TOML arrays keep the order
    ## they are written in, so the cut-off filter and the sort permutation are
    ## built once here and shared with the gross vectors below. Pairing them
    ## per-vector instead would silently mismatch an out-of-order block.
    _brk = let key = "confirmed_break_dates"
        if haskey(raw, key)
            blk = raw[key]
            ds = String.(blk["value"])
            keep = [Date(String(d)) <= cutoff for d in ds]
            days = Int[_index(d) for d in ds[keep]]
            (; blk, ds, keep, ord = sortperm(days), days)
        else
            (; blk = nothing, ds = String[], keep = Bool[], ord = Int[],
                days = Int[])
        end
    end
    confirmed_break_days = _brk.days[_brk.ord]
    ## The same days as calendar dates, in the same order, so a validation
    ## failure below names the date the user wrote rather than a grid index.
    confirmed_break_date_labels = String.(_brk.ds[_brk.keep])[_brk.ord]

    ## The printed 24h new-confirmed counts on each break day, which make the
    ## step data-derived rather than a guessed prior width. The submodels
    ## centre each step on `observed increment − gross`, the part of the
    ## vintage step the report itself attributes to base integration rather
    ## than to the day's notifications. Filtered and permuted with the dates so
    ## the three stay aligned. Absent counts default to zeros, which makes each
    ## centre the whole increment rather than anything neutral. The entire
    ## step is then attributed to the artefact. That errs towards artefact
    ## rather than towards incidence, which is the safe direction for a
    ## de-anchored day (see `break_step_centres`).
    function break_gross(key)
        _brk.blk === nothing && return Int[]
        haskey(_brk.blk, key) || return zeros(Int, length(_brk.ord))
        vals = Int.(_brk.blk[key])
        length(vals) == length(_brk.ds) || error(
            "confirmed_break_dates: $key has $(length(vals)) entries for " *
            "$(length(_brk.ds)) dates")
        return vals[_brk.keep][_brk.ord]
    end
    confirmed_break_gross_cases = break_gross("gross_cases")
    confirmed_break_gross_deaths = break_gross("gross_deaths")

    reported_history = history("reported_case_history")
    confirmed_history = history("confirmed_case_history")
    confirmed_deaths_history = history("confirmed_death_history")
    deaths_history = history("death_history")

    ## Validate each listed break day, once here rather than inside the models,
    ## where it would re-run on every likelihood evaluation. Both confirmed
    ## streams are checked separately against their own gross vector and their
    ## own increments: a harmonisation can be well behaved on one and not the
    ## other, so the failing stream is named.
    ##
    ## Two configurations are refused, both of which are silent today.
    ##
    ## 1. `gross >= increment`. A day whose printed 24h count already covers its
    ##    whole vintage step is a provincial transfer, not a base integration.
    ##    The diagnostic is the direction: an integration reattaches cases and
    ##    deaths and so adds to both, while a transfer moves both down (SitRep
    ##    065, 18 July 2026, is the worked example: +83 gross against +77 net
    ##    cases and +40 against +37 net deaths, as Haut-Uele's cumulative was
    ##    revised down to 16 cases and 10 deaths). Listing such a day is harmful
    ##    rather than merely useless: `confirmed_cases_model` de-anchors a listed
    ##    day from the laboratory positivity denominator, so it loses its
    ##    `BetaBinomial` term while the step that should absorb the backlog is
    ##    centred at or below zero. Measured on `confirmed_only_model`
    ##    (`scripts/diag_break_attribution.jl`, 500 draws x 2 chains): 94
    ##    divergences and a min bulk ESS of 15, against 20 and 522 with no break
    ##    day declared, and the cut-off infection count inflated 14% as the fit
    ##    books the artefact as incidence. Pinning the step at a published
    ##    discrepancy instead restores 22 divergences and 477 ESS.
    ##
    ## 2. A date matching no vintage in the history. It does nothing at all: no
    ##    step, no de-anchor, no error, and the `gross` check above cannot fire
    ##    because there is no increment to compare. A transposed digit or wrong
    ##    month therefore presents as silence while the user believes a
    ##    harmonisation is being absorbed.
    ##
    ## A gross of zero (the default when the key is absent) is legal but warned:
    ## it makes the centre the whole increment, attributing all of it to the
    ## artefact. That errs towards artefact rather than towards incidence, the
    ## safe direction, but it is not the neutral choice it looks like.
    function check_break_gross(gross, hist, label)
        isempty(confirmed_break_days) && return nothing
        isempty(hist.counts) && return nothing
        inc = diff(vcat(0, collect(hist.counts)))
        hdays = collect(hist.days)
        for (i, d) in enumerate(confirmed_break_days)
            date = confirmed_break_date_labels[i]
            pos = findfirst(==(d), hdays)
            if pos === nothing
                error("confirmed_break_dates: $date (grid day $d) matches no " *
                      "vintage in the $label history, so it would be silently " *
                      "ignored — no step and no de-anchor — while appearing to " *
                      "absorb a harmonisation. Check the date against the " *
                      "history's own vintages; a transposed digit or the wrong " *
                      "month presents exactly like this.")
            end
            g = i <= length(gross) ? gross[i] : 0
            if g >= inc[pos]
                error("confirmed_break_dates: $date has a printed 24h $label " *
                      "count of $g against a net vintage increment of " *
                      "$(inc[pos]), so the gross does not sit below the net. " *
                      "That is a provincial TRANSFER, not a base integration: " *
                      "an integration reattaches records and so ADDS to cases " *
                      "and deaths together, whereas a transfer moves both " *
                      "DOWN (SitRep 065, 18 July 2026: +83 gross vs +77 net " *
                      "cases and +40 vs +37 net deaths). Listing it de-anchors " *
                      "the positivity denominator with no backlog to absorb, " *
                      "which measured 94 divergences and a min bulk ESS of 15 " *
                      "against 20 and 522 undeclared, and inflated the cut-off " *
                      "infection count 14%. Remove $date from " *
                      "[confirmed_break_dates].")
            end
            if g == 0
                @warn "confirmed_break_dates: no printed 24h $label count for " *
                      "$date, so its step is centred on the whole increment " *
                      "and attributes all of it to the harmonisation rather " *
                      "than splitting it. Supply gross_$label to make the " *
                      "split data-derived." increment=inc[pos]
            end
        end
        return nothing
    end
    check_break_gross(confirmed_break_gross_cases, confirmed_history, "cases")
    check_break_gross(confirmed_break_gross_deaths, confirmed_deaths_history,
        "deaths")
    ## The analysed-specimen series is the laboratory denominator. The
    ## received series is recorded for the pipeline view but not fitted.
    lab_history = history("tests_analysed_history")
    tests_received_history = history("tests_received_history")
    ## Post-cutoff 24h analysed counts (daily increments, not cumulative)
    ## on the trusted post-cutoff days: the confirmed model pairs each with
    ## that day's confirmed increment as a Binomial-denominator window.
    lab_daily_history = history("tests_analysed_daily_history")
    ## Post-cutoff daily new-suspect inflow ("nouveaux cas suspects du jour"):
    ## per-day counts (not cumulative) of newly reported suspects, fitted as
    ## a daily incidence against the modelled suspected series where the
    ## frozen cumulative suspected stream stops at 26 May.
    suspected_daily_history = history("suspected_daily_history")
    ## Post-cutoff daily new suspected deaths ("cas suspects du jour N (M
    ## deces)"): per-day counts (not cumulative) of suspected deaths in the
    ## preceding 24h, fitted as a daily incidence against the modelled
    ## suspected-death series where the frozen cumulative suspected-death stream
    ## stops at 26 May (the deaths analogue of `suspected_daily_history`).
    suspected_daily_deaths_history = history("suspected_daily_deaths_history")
    ## Post-cutoff daily isolation/hospitalisation occupancy ("Patients en
    ## isolement"): a per-day count of patients in an isolation/treatment bed,
    ## fitted against the modelled bed count on each report day by the
    ## length-of-stay observation submodel. Begins 1 June where the
    ## all-patients column definition is stable (see the manifest note).
    isolation_history = history("isolation_history")
    ## Implied isolation/treatment-bed capacity (occupancy / reported
    ## occupancy rate) on the days a rate is published: fitted by the
    ## isolation submodel as noisy observations of the national bed capacity
    ## the latent bed demand saturates against.
    bed_capacity_history = history("bed_capacity_history")
    ## Post-cutoff cumulative recovered-among-confirmed ("cumul guéris"):
    ## a cumulative count of laboratory-confirmed cases recorded as recovered,
    ## fitted as survivors among the modelled confirmed cases (a scaled
    ## confirmation-to-recovery convolution) by the recovered observation
    ## submodel. Begins 6 June, where the confirmed-based reports first print
    ## the running total.
    recovered_history = history("recovered_history")
    ## Post-cutoff daily treatment-centre patient-movement flows (Tableau 6
    ## "Mouvement des patients", national): per-day counts (not cumulative) of
    ## admissions, in-care deaths (suspects + confirmed), rule-out discharges
    ## (non-cas) and absconded patients. Optional refinements of the
    ## treatment-flow submodel over their 13-22 June overlap. The longer
    ## occupancy / recovered / capacity stock streams carry the earlier window.
    ## An absent block loads empty and is a no-op in the model.
    treatment_admissions_history = history("treatment_admissions_history")
    treatment_deaths_history = history("treatment_deaths_history")
    treatment_ruleout_history = history("treatment_ruleout_history")
    treatment_absconded_history = history("treatment_absconded_history")
    ## Post-cutoff start-of-day in-bed count (Tableau 6 "Patients au lit
    ## (J-1)", national): the bed stock at the start of each report day. The
    ## treatment-flow submodel differences it against the previous day's
    ## occupancy to identify the between-report DHIS2 reclassification days and
    ## centre a fitted break-step prior, accumulating the sampled steps into a
    ## cumulative offset added to the modelled occupancy mean, so the modelled
    ## occupancy tracks the reclassification without bending Rt while still
    ## partitioning each step into reporting-artifact vs real demand. Optional:
    ## an empty block leaves the offset at zero, a no-op.
    ## Tableau 6 occupancy split (13-23 June): the `dont confirmes (NC+AC)` and
    ## `dont suspects` sub-rows of the `Patients en isolement (Fin J)` stock.
    ## Census (prevalence) sub-stocks, not flows: each is the count of that
    ## class of patient occupying a bed at end-of-day, and the two sum to the
    ## total occupancy (= `isolation_history`) exactly. On the days they are
    ## present the treatment-flow submodel scores the two sub-stocks in place of
    ## the total occupancy (a per-day total-or-split switch). An absent or
    ## empty block falls back to the total-occupancy likelihood and is a no-op.
    treatment_confirmed_incare_history = history("treatment_confirmed_incare_history")
    treatment_suspect_incare_history = history("treatment_suspect_incare_history")
    ## Cut-off scalar from an explicit TOML block, else the final
    ## (most recent) vintage of the matching history. When a `cutoff_date`
    ## freeze is active the explicit TOML scalars (which hold the final,
    ## full-data total) no longer match the truncated history, so the
    ## frozen final vintage is used instead.
    _hist_end(h) = isempty(h.counts) ? missing : h.counts[end]
    frozen = !isnothing(cutoff_date)
    _scalar(k, h) = (frozen || !haskey(raw, k)) ? _hist_end(h) : Int(_val(k))
    ## The first WHO joint situation report is the earliest reported-case
    ## vintage. Days from it to the cut-off set the intervention breakpoint.
    who_first_sitrep_days = isempty(reported_history.days) ? n :
                            n - reported_history.days[1] + 1

    ## Cut-off export scalars. When freezing to an earlier date the dated
    ## series is truncated, so the cumulative export totals are the number
    ## of dated detections/deaths kept (matching the per-day series). The
    ## full-data manifest scalars are used otherwise.
    exported_cases = frozen ? length(export_case_days) :
                     Int(_val("exported_cases"))
    exports_deaths = frozen ? length(export_death_days) :
                     Int(_val("exports_deaths"))

    return (; n, cutoff, seeding,
        exported_cases = exported_cases,
        exports_deaths = exports_deaths,
        export_case_days = export_case_days,
        export_death_days = export_death_days,
        total_deaths = frozen ?
                       _hist_end(deaths_history) : Int(_val("total_deaths")),
        reported_cases = frozen ?
                         _hist_end(reported_history) :
                         Int(_val("reported_cases")),
        confirmed_cases = _scalar("confirmed_cases", confirmed_history),
        confirmed_deaths = _scalar("confirmed_deaths",
            confirmed_deaths_history),
        tests_analysed = _scalar("cumulative_tests_analysed", lab_history),
        reported_history = reported_history,
        confirmed_history = confirmed_history,
        confirmed_deaths_history = confirmed_deaths_history,
        deaths_history = deaths_history,
        lab_history = lab_history,
        lab_daily_history = lab_daily_history,
        suspected_daily_history = suspected_daily_history,
        suspected_daily_deaths_history = suspected_daily_deaths_history,
        isolation_history = isolation_history,
        bed_capacity_history = bed_capacity_history,
        recovered_history = recovered_history,
        recovered_cases = _scalar("recovered_cases", recovered_history),
        treatment_admissions_history = treatment_admissions_history,
        treatment_deaths_history = treatment_deaths_history,
        treatment_ruleout_history = treatment_ruleout_history,
        treatment_absconded_history = treatment_absconded_history,
        treatment_confirmed_incare_history =
        treatment_confirmed_incare_history,
        treatment_suspect_incare_history = treatment_suspect_incare_history,
        occupancy_break_days = occupancy_break_days,
        confirmed_break_days = confirmed_break_days,
        confirmed_break_gross_cases = confirmed_break_gross_cases,
        confirmed_break_gross_deaths = confirmed_break_gross_deaths,
        tests_received_history = tests_received_history,
        tmrca_days = _gap(raw["genetic_tmrca"]["date"]),
        who_first_sitrep_days)
end

"""
    freeze_observations(cutoff_date; path = default manifest)

Load the observation manifest frozen to `cutoff_date`: the cut-off is
moved to `cutoff_date` and every dated history is truncated to the
vintages available by then, so the returned named tuple is what the
renewal model would have seen on that date. The cut-off scalar totals
(`reported_cases`, `total_deaths`, `confirmed_cases`, ...) are taken
from the truncated histories rather than the manifest's full-data
scalars. Use it to re-evaluate the renewal estimate at a past report
date (for example a McCabe et al. situation-report cut-off) for a
like-for-like, matched-in-time comparison.

`cutoff_date` accepts a `Date` or an ISO date string. It must be on or
after the earliest history vintage in the manifest (the renewal DRC
series begins 18 May 2026). An earlier date leaves the suspected
streams empty and is not a meaningful renewal fit.
"""
function freeze_observations(
        cutoff_date::Union{Date, AbstractString};
        path::AbstractString = joinpath(@__DIR__, "..", "data",
            "observations.toml"),
        seeding_lead::Integer = SEEDING_LEAD_DAYS)
    return load_observations(path; seeding_lead, cutoff_date)
end

"""
    m_prior_centre(as_of_date; base_date, m_base, doubling_days)

Centre for the doubling-count prior `m`, based on `m_base` doublings at
`base_date` and advancing by one doubling per `doubling_days` of elapsed
time to `as_of_date`:

```math
m_0 = m_\\text{base} +
    \\frac{\\text{as\\_of} - \\text{base}}{\\text{doubling\\_days}}.
```

The base is McCabe et al.'s first report (18 May 2026; Method 2 central
501 cases ⇒ `m ≈ 9`), advancing at the outbreak-specific doubling time
(`M_PRIOR_DOUBLING_DAYS`, the BEAST X estimate of mbalaplacide2026,
mean 11.7 d), so the prior stays centred on the plausible outbreak
size as the cut-off moves. `C_T = 2^m` is the cumulative infection
count. 9 is a weakly-informative centre of the same order. Passed into
[`exponential_growth_model`](@ref) as the centre of the wide `m` prior.
"""
function m_prior_centre(as_of_date::Union{Date, AbstractString};
        base_date::AbstractString = M_PRIOR_BASE_DATE,
        m_base::Real = M_PRIOR_BASE,
        doubling_days::Real = M_PRIOR_DOUBLING_DAYS)
    as_of = as_of_date isa Date ? as_of_date : Date(String(as_of_date))
    elapsed = date2epochdays(as_of) - date2epochdays(Date(base_date))
    return m_base + elapsed / doubling_days
end
