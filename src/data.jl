# Observation loading from the dated TOML manifest. The manifest stores
# calendar dates and cumulative counts, never grid day-indices, so the data
# stay aligned by date as vintages are added, revised or arrive sparsely.
# The renewal model works on a daily grid, so this loader derives the grid
# length and the per-vintage day-indices from the dates at load time. The
# cut-off is the last grid day and the seeding day is the first.

## Days the grid extends before the genetic TMRCA date, so the seeding
## crossing has room to be inferred below the molecular-clock bound.
const SEEDING_LEAD_DAYS = 30

"""
The `(grid day, correction)` pairs a confirmed stream's listed
harmonisation-break days carry, one per listed day the stream's own history
has a matching vintage for.

On a listed break day (`[confirmed_break_dates]`) the reported cumulative
jumps by more than that day's notifications because INSP reattached
previously unlinked records. The excess is `net - gross`, the vintage's own
step less the printed 24h count. Subtracting it turns a raw cumulative
difference into the count actually notified across the window.

`deaths` selects the confirmed-death stream rather than confirmed cases. A
break day with no matching vintage is left out. A break day on the stream's
first vintage takes the whole cumulative as its step. The correction is
floored at zero so a smaller vintage cannot inflate the truth.
"""
function confirmed_break_steps(obs; deaths::Bool = false)
    out = Tuple{Int, Float64}[]
    days = obs.confirmed_break_days
    isempty(days) && return out
    hist = deaths ? obs.confirmed_deaths_history : obs.confirmed_history
    gross = deaths ? obs.confirmed_break_gross_deaths :
        obs.confirmed_break_gross_cases
    isempty(hist.counts) && return out
    hdays = collect(hist.days)
    counts = collect(hist.counts)
    for (i, d) in enumerate(days)
        pos = findfirst(==(d), hdays)
        pos === nothing && continue
        net = counts[pos] - (pos == 1 ? 0 : counts[pos - 1])
        g = i <= length(gross) ? gross[i] : 0
        push!(out, (d, float(max(net - g, 0))))
    end
    return out
end

"""
Total harmonisation correction carried by a confirmed stream over the grid
days `(from_day, to_day]`, summing [`confirmed_break_steps`](@ref) over
every listed break day in the window.

The window is half open on the left, so a break day on the origin belongs to
the window before it. `deaths` selects the confirmed-death stream rather
than confirmed cases. Returns zero when the window holds no break day.
"""
function confirmed_break_correction(
        obs, from_day::Real, to_day::Real;
        deaths::Bool = false
    )
    total = 0.0
    for (d, c) in confirmed_break_steps(obs; deaths = deaths)
        from_day < d <= to_day && (total += c)
    end
    return total
end

"""
Load the BVD observation manifest from `path` (a dated TOML file) and
return a named tuple for the renewal model. Calendar dates are converted
to 1-based grid day-indices, day 1 the seeding day and day `n` the
cut-off, using the cut-off (`as_of_date`) and a seeding day placed
`seeding_lead` days before the genetic TMRCA date.

Returns the grid length `n`, the `cutoff` and `seeding` dates, and the
per-stream cumulative totals at the cut-off (`reported_cases`,
`total_deaths`, `confirmed_cases`, `confirmed_deaths`, `tests_analysed`,
`exported_cases`, `exports_deaths`). A scalar with no explicit TOML block
is taken from the final vintage of the matching history.

The dated Uganda export series are grid day-indices (`export_case_days`,
`export_death_days`), each a sorted list of detection or death days on or
before the cut-off.

The per-vintage histories are returned as `(; days, counts)` with `days`
the grid day-indices: `reported_history`, `confirmed_history`,
`confirmed_deaths_history`, `deaths_history`, `lab_history` (cumulative
analysed specimens), `lab_daily_history` (24h analysed counts),
`suspected_daily_history` (daily new-suspect inflow),
`suspected_daily_deaths_history` (daily new suspected deaths),
`isolation_history` (daily isolation occupancy), `bed_capacity_history`
(occupancy / reported occupancy rate), `recovered_history` (cumulative
recovered among confirmed), `treatment_confirmed_incare_history` and
`treatment_suspect_incare_history` (the occupancy split into two
prevalence sub-stocks that sum to the total), and
`tests_received_history`.

The digitised symptom-onset reporting triangle is returned as
`onset_curve_history`, the per-vintage increments read from
`onset_curve_path`, by default `onset_curve_scanned.csv` alongside `path`.
A manifest read from elsewhere has no such sibling, so name the triangle
explicitly to keep the stream rather than degrade it to a no-op. See
[`load_onset_curve`](@ref) for the dedup and increment construction. The
same triangle's cumulative confirmed-by-onset total is returned as
`onset_report_history` in the usual shape.

Also returned are the genetic TMRCA bound `tmrca_days` (days before the
cut-off) and `who_first_sitrep_days` (days from the earliest
reported-case vintage to the cut-off). The intervention breakpoint grid
day is `n - who_first_sitrep_days`.
"""
function load_observations(
        path::AbstractString = joinpath(
            @__DIR__, "..", "data",
            "observations.toml"
        );
        seeding_lead::Integer = SEEDING_LEAD_DAYS,
        cutoff_date::Union{Nothing, Date, AbstractString} = nothing,
        onset_curve_path::AbstractString = joinpath(
            dirname(path),
            "onset_curve_scanned.csv"
        )
    )
    raw = TOML.parsefile(path)
    _val(k) = raw[k]["value"]
    ## The cut-off is the manifest `as_of_date` unless an earlier
    ## `cutoff_date` is supplied by `freeze_observations`. Freezing only
    ## moves the cut-off earlier, so the grid stays date-aligned with the
    ## full-data fit.
    cutoff = isnothing(cutoff_date) ? Date(String(raw["as_of_date"])) :
        (
            cutoff_date isa Date ? cutoff_date :
            Date(String(cutoff_date))
        )
    ## Grid day-index (1-based) of a calendar date: seeding day is day 1.
    _gap(d) = Int(date2epochdays(cutoff) - date2epochdays(Date(String(d))))
    tmrca_date = Date(String(raw["genetic_tmrca"]["date"]))
    seeding = tmrca_date - Day(seeding_lead)
    n = Int(date2epochdays(cutoff) - date2epochdays(seeding)) + 1
    _index(d) = n - _gap(d)

    ## A dated cumulative history → grid day-indices and counts, sorted
    ## oldest-first so the model differences consecutive vintages into
    ## daily increments. Vintages after the cut-off are dropped. Empty
    ## when the block is absent.
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
    ## day-indices on or before the cut-off, sorted ascending. Empty when
    ## the block is absent or every date falls after the cut-off.
    function event_days(key)
        haskey(raw, key) || return Int[]
        ds = String.(raw[key]["value"])
        keep = [Date(d) <= cutoff for d in ds]
        idx = Int[_index(d) for d in ds[keep]]
        return sort(idx)
    end

    ## Dated Uganda export detection and death days (1-based grid indices).
    ## Each contributes one Poisson term at its day, against the per-day
    ## expected export count (see `exports_model`).
    export_case_days = event_days("export_case_dates")
    export_death_days = event_days("export_death_dates")

    ## Opt-in occupancy reclassification-break days: grid days on which the
    ## treatment-flow model fits a level step into the modelled occupancy
    ## mean, absorbing a between-report measurement-basis discontinuity in
    ## the observed isolation series without bending Rt. Absent or empty is
    ## a no-op.
    occupancy_break_days = event_days("occupancy_break_dates")

    ## Opt-in retrospective harmonisation-break days for the confirmed
    ## streams: grid days on which INSP integrated a harmonised provincial
    ## base, so the cumulative headline steps by far more than that day's own
    ## notifications. The confirmed submodels then de-anchor the laboratory
    ## positivity denominator, since the reattached cases are not same-day
    ## positives, and fit a level step so the increment likelihood does not
    ## read the backlog as one day of incidence. Absent or empty is a no-op.
    ## The grid days come back sorted while the TOML arrays keep their
    ## written order, so the cut-off filter and the sort permutation are
    ## built once here and shared with the gross vectors below. Pairing them
    ## per-vector would silently mismatch an out-of-order block.
    _brk = let key = "confirmed_break_dates"
        if haskey(raw, key)
            blk = raw[key]
            ds = String.(blk["value"])
            keep = [Date(String(d)) <= cutoff for d in ds]
            days = Int[_index(d) for d in ds[keep]]
            (; blk, ds, keep, ord = sortperm(days), days)
        else
            (;
                blk = nothing, ds = String[], keep = Bool[], ord = Int[],
                days = Int[],
            )
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
    ## than to the day's notifications. Filtered and permuted with the dates
    ## so the three stay aligned. Absent counts default to zeros, which
    ## centres on the whole increment and attributes all of it to the
    ## artefact (see `break_step_centres`).
    function break_gross(key)
        _brk.blk === nothing && return Int[]
        haskey(_brk.blk, key) || return zeros(Int, length(_brk.ord))
        vals = Int.(_brk.blk[key])
        length(vals) == length(_brk.ds) || error(
            "confirmed_break_dates: $key has $(length(vals)) entries for " *
                "$(length(_brk.ds)) dates"
        )
        return vals[_brk.keep][_brk.ord]
    end
    confirmed_break_gross_cases = break_gross("gross_cases")
    confirmed_break_gross_deaths = break_gross("gross_deaths")

    ## Per-province history from a TOML block with one shared `dates` array
    ## and one array per series. Series are keyed by province
    ## (`province_confirmed_history`) or by province-and-measure
    ## (`province_lab_daily_history`, e.g. `ituri_analysed`). Returns a Dict
    ## from series name to the same `(; days, counts)` shape as `history`,
    ## empty when the block is absent. Counts are cumulative or daily
    ## according to the block, and the caller knows which.
    function province_history(key)
        ProvHistory = @NamedTuple{days::Vector{Int}, counts::Vector{Int}}
        !haskey(raw, key) && return Dict{String, ProvHistory}()
        block = raw[key]
        !haskey(block, "dates") && return Dict{String, ProvHistory}()
        provinces = sort!([k for k in keys(block) if k != "dates" && k != "source"])
        result = Dict{String, ProvHistory}()
        keep = [Date(String(d)) <= cutoff for d in block["dates"]]
        idx = Int[_index(d) for d in block["dates"][keep]]
        ord = sortperm(idx)
        for prov in provinces
            vals = Int.(block[prov][keep])
            result[prov] = (; days = idx[ord], counts = vals[ord])
        end
        return result
    end

    ## Per-zone history from a TOML block with one shared `dates` array and
    ## one dotted `province.zone` array per series (`zone_confirmed_history`
    ## and `zone_death_history`). Returns a Dict keyed province → zone, each
    ## zone in the same (; days, counts) shape as history(); `unallocated`
    ## is the report's own row of counts not yet attributed to a zone and
    ## is kept as a zone of that name. Empty when the block is absent.
    function zone_history(key)
        ZoneHistory = @NamedTuple{days::Vector{Int}, counts::Vector{Int}}
        result = Dict{String, Dict{String, ZoneHistory}}()
        !haskey(raw, key) && return result
        block = raw[key]
        !haskey(block, "dates") && return result
        keep = [Date(String(d)) <= cutoff for d in block["dates"]]
        idx = Int[_index(d) for d in block["dates"][keep]]
        ord = sortperm(idx)
        for prov in sort!([k for k in keys(block) if block[k] isa AbstractDict])
            zones = Dict{String, ZoneHistory}()
            for (zone, vals) in block[prov]
                length(vals) == length(block["dates"]) || error(
                    "$key: $prov.$zone has $(length(vals)) entries for " *
                        "$(length(block["dates"])) dates"
                )
                v = Int.(vals[keep])
                zones[String(zone)] = (; days = idx[ord], counts = v[ord])
            end
            result[prov] = zones
        end
        return result
    end

    reported_history = history("reported_case_history")
    confirmed_history = history("confirmed_case_history")
    confirmed_deaths_history = history("confirmed_death_history")
    deaths_history = history("death_history")

    ## Validate each listed break day once here rather than inside the
    ## models, where it would re-run on every likelihood evaluation. Each
    ## confirmed stream is checked against its own gross vector and its own
    ## increments, since a harmonisation can be well behaved on one and not
    ## the other, so the failing stream is named.
    ##
    ## Two configurations are refused, both otherwise silent.
    ##
    ## 1. `gross >= increment`. A day whose printed 24h count already covers
    ##    its whole vintage step is a provincial transfer, not a base
    ##    integration. An integration reattaches cases and deaths and so adds
    ##    to both, while a transfer moves both down. Listing one de-anchors
    ##    the positivity denominator with no backlog to absorb, and the step
    ##    that should take the backlog is centred at or below zero. Measured
    ##    on `confirmed_only_model` at 500 draws x 2 chains: 94 divergences
    ##    and a min bulk ESS of 15, against 20 and 522 with no break day
    ##    declared, and the cut-off infection count inflated 14% as the fit
    ##    books the artefact as incidence. Pinning the step at a published
    ##    discrepancy gives 22 divergences and 477 ESS.
    ##
    ## 2. A date matching no vintage in the history. It does nothing at all,
    ##    and the `gross` check cannot fire because there is no increment to
    ##    compare against, so a transposed digit presents as silence.
    ##
    ## A gross of zero, the default when the key is absent, is legal but
    ## warned. It centres on the whole increment and attributes all of it to
    ## the artefact, which is not the neutral choice it looks like.
    function check_break_gross(gross, hist, label)
        isempty(confirmed_break_days) && return nothing
        isempty(hist.counts) && return nothing
        inc = diff(vcat(0, collect(hist.counts)))
        hdays = collect(hist.days)
        for (i, d) in enumerate(confirmed_break_days)
            date = confirmed_break_date_labels[i]
            pos = findfirst(==(d), hdays)
            if pos === nothing
                error(
                    "confirmed_break_dates: $date (grid day $d) matches no " *
                        "vintage in the $label history, so it would be silently " *
                        "ignored — no step and no de-anchor — while appearing to " *
                        "absorb a harmonisation. Check the date against the " *
                        "history's own vintages; a transposed digit or the wrong " *
                        "month presents exactly like this."
                )
            end
            g = i <= length(gross) ? gross[i] : 0
            if g >= inc[pos]
                error(
                    "confirmed_break_dates: $date has a printed 24h $label " *
                        "count of $g against a net vintage increment of " *
                        "$(inc[pos]), so the gross does not sit below the net. " *
                        "That is a provincial transfer, not a base " *
                        "integration. An integration reattaches records and " *
                        "so adds to cases and deaths together, whereas a " *
                        "transfer moves both down (SitRep 065, 18 July " *
                        "2026: +83 gross vs +77 net " *
                        "cases and +40 vs +37 net deaths). Listing it de-anchors " *
                        "the positivity denominator with no backlog to absorb, " *
                        "which measured 94 divergences and a min bulk ESS of 15 " *
                        "against 20 and 522 undeclared, and inflated the cut-off " *
                        "infection count 14%. Remove $date from " *
                        "[confirmed_break_dates]."
                )
            end
            if g == 0
                @warn "confirmed_break_dates: no printed 24h $label count for " *
                    "$date, so its step is centred on the whole increment " *
                    "and attributes all of it to the harmonisation rather " *
                    "than splitting it. Supply gross_$label to make the " *
                    "split data-derived." increment = inc[pos]
            end
        end
        return nothing
    end
    check_break_gross(confirmed_break_gross_cases, confirmed_history, "cases")
    check_break_gross(
        confirmed_break_gross_deaths, confirmed_deaths_history,
        "deaths"
    )
    ## The analysed-specimen series is the laboratory denominator. The
    ## received series is recorded for the pipeline view but not fitted.
    lab_history = history("tests_analysed_history")
    tests_received_history = history("tests_received_history")
    ## 24h analysed counts (daily increments, not cumulative). The confirmed
    ## model pairs each with that day's confirmed increment as a
    ## Binomial-denominator window.
    lab_daily_history = history("tests_analysed_daily_history")
    ## Daily new-suspect inflow, "nouveaux cas suspects du jour". Per-day
    ## counts, fitted as a daily incidence against the modelled suspected
    ## series where the cumulative suspected stream stops at 26 May.
    suspected_daily_history = history("suspected_daily_history")
    ## The deaths analogue, "cas suspects du jour N (M deces)": per-day
    ## suspected deaths in the preceding 24h.
    suspected_daily_deaths_history = history("suspected_daily_deaths_history")
    ## Daily isolation occupancy, "Patients en isolement". A per-day count of
    ## patients in a bed, fitted by the length-of-stay submodel. Begins 1
    ## June where the all-patients column definition is stable.
    isolation_history = history("isolation_history")
    ## Implied bed capacity (occupancy / reported occupancy rate) on the days
    ## a rate is published, fitted as noisy observations of the national
    ## capacity the latent bed demand saturates against.
    bed_capacity_history = history("bed_capacity_history")
    ## Cumulative recovered among confirmed, "cumul guéris", fitted as
    ## survivors among the modelled confirmed cases. Begins 6 June, where the
    ## reports first print the running total.
    recovered_history = history("recovered_history")
    ## Daily treatment-centre patient-movement flows: admissions, in-care
    ## deaths, rule-out discharges and absconded patients. Optional
    ## refinements of the treatment-flow submodel over their 13-22 June
    ## overlap. The longer stock streams carry the earlier window. An absent
    ## block loads empty and is a no-op.
    treatment_admissions_history = history("treatment_admissions_history")
    treatment_deaths_history = history("treatment_deaths_history")
    treatment_ruleout_history = history("treatment_ruleout_history")
    treatment_absconded_history = history("treatment_absconded_history")
    ## The occupancy split, `dont confirmes (NC+AC)` and `dont suspects`.
    ## Census sub-stocks rather than flows: each counts that class of patient
    ## occupying a bed at end-of-day, and the two sum to `isolation_history`
    ## exactly. On the days both are present the treatment-flow submodel
    ## scores them in place of the total occupancy. An absent or empty block
    ## falls back to the total-occupancy likelihood.
    treatment_confirmed_incare_history = history("treatment_confirmed_incare_history")
    treatment_suspect_incare_history = history("treatment_suspect_incare_history")
    ## Digitised symptom-onset reporting-triangle increments, from a sibling
    ## CSV rather than the TOML manifest, since it is a 2D onset-date by
    ## report-date triangle and not a single dated series. Filtered to the
    ## same `cutoff` as every history above, so a freeze also freezes this
    ## stream.
    onset_curve_history = load_onset_curve(
        onset_curve_path; cutoff, seeding
    )
    ## The same triangle's cumulative confirmed-by-onset total, restated in
    ## the shape every other history uses so it is scored by the same
    ## machinery. The fitted stream never sees this series, since it fits the
    ## increments.
    onset_report_history = (;
        days = onset_curve_history.total_days,
        counts = onset_curve_history.total_counts,
    )
    ## Cut-off scalar from an explicit TOML block, else the final vintage of
    ## the matching history. Under a freeze the TOML scalars hold the
    ## full-data total and no longer match the truncated history, so the
    ## frozen final vintage is used instead.
    _hist_end(h) = isempty(h.counts) ? missing : h.counts[end]
    frozen = !isnothing(cutoff_date)
    _scalar(k, h) = (frozen || !haskey(raw, k)) ? _hist_end(h) : Int(_val(k))
    ## The first WHO joint situation report is the earliest reported-case
    ## vintage. Days from it to the cut-off set the intervention breakpoint.
    who_first_sitrep_days = isempty(reported_history.days) ? n :
        n - reported_history.days[1] + 1

    ## Cut-off export scalars. A freeze truncates the dated series, so the
    ## cumulative totals are the number of dated events kept, matching the
    ## per-day series. Otherwise the manifest scalars are used.
    exported_cases = frozen ? length(export_case_days) :
        Int(_val("exported_cases"))
    exports_deaths = frozen ? length(export_death_days) :
        Int(_val("exports_deaths"))

    return (;
        n, cutoff, seeding,
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
        confirmed_deaths = _scalar(
            "confirmed_deaths",
            confirmed_deaths_history
        ),
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
        onset_curve_history = onset_curve_history,
        onset_report_history = onset_report_history,
        province_confirmed_history = province_history("province_confirmed_history"),
        province_death_history = province_history("province_death_history"),
        province_lab_daily_history = province_history("province_lab_daily_history"),
        zone_confirmed_history = zone_history("zone_confirmed_history"),
        zone_death_history = zone_history("zone_death_history"),
        tmrca_days = _gap(raw["genetic_tmrca"]["date"]),
        who_first_sitrep_days,
    )
end

"""
    patch_members(province_names)

Source provinces each patch pools, in the order of `province_names`.

A name with no [`PROVINCE_MEMBERS`](@ref) entry is its own province, which
keeps the province helpers usable with an arbitrary province list.
[`province_increment_matrix`](@ref) and
[`province_testing_covariate`](@ref) both resolve patches through this, so
how a patch maps to the manifest blocks it covers is written once.
"""
function patch_members(province_names::AbstractVector)
    return [get(PROVINCE_MEMBERS, nm, [nm]) for nm in province_names]
end

"""
    province_increment_matrix(province_history, province_names, n_patches)

Reshape the per-province cumulative histories loaded by
[`load_observations`](@ref) into the `(n_patches × n_vintages)` matrix of
new-confirmed counts that [`province_composition_model`](@ref) scores,
together with the shared vintage day indices.

Every province must be reported on the same vintage days, which the
composition likelihood requires. It allocates each vintage's national
total across the provinces, so a province missing from a vintage would
silently shift cases into the others. A mismatch is an error.

Returns `(; days, increments)`. When no per-province data is supplied,
`days` is empty and the caller skips the composition term.
"""
function province_increment_matrix(
        province_history,
        province_names::AbstractVector, n_patches::Integer
    )
    empty = (; days = Int[], increments = Matrix{Int}(undef, 0, 0))
    isempty(province_history) && return empty
    names = province_names[1:min(n_patches, length(province_names))]
    members = patch_members(names)
    any(ms -> any(m -> !haskey(province_history, m), ms), members) &&
        return empty
    hists = [[province_history[m] for m in ms] for ms in members]
    days = hists[1][1].days
    isempty(days) && return empty
    for (ms, hs) in zip(members, hists), (m, h) in zip(ms, hs)

        h.days == days || error(
            "province `$(m)` is reported on different vintage days to " *
                "`$(first(members)[1])`; the composition likelihood needs " *
                "every province on the same vintages."
        )
    end
    ## Cumulative → per-vintage increments. The first increment is the
    ## cumulative to the first vintage day, matching `bin_increments`, which
    ## bins the modelled daily series from day 1 to `days[1]`.
    ##
    ## A province's cumulative count can fall between vintages when cases are
    ## reclassified, so the raw difference goes negative. The composition
    ## likelihood cannot take that, since it scores a count drawn from a
    ## total. Clamping at zero reads a downward revision as no new cases in
    ## that province this vintage. The composition conditions on the sum of
    ## these increments rather than on the national total, so the clamp stays
    ## self-consistent. A pooled patch is the sum of its members' cumulative
    ## counts differenced once. Summing before differencing keeps a downward
    ## revision in one member from being clamped away while another rises,
    ## which would inflate the patch.
    increments = Matrix{Int}(undef, length(names), length(days))
    for (p, hs) in enumerate(hists)
        pooled = reduce(.+, (collect(h.counts) for h in hs))
        increments[p, :] = max.(diff(vcat(0, pooled)), 0)
    end
    return (; days, increments)
end

"""
    zone_cumulative_falls(zone_history; min_fall = 1)

Every place a named zone's cumulative count falls between consecutive
vintages by more than `min_fall`, as a vector of named tuples
`(; province, zone, day, from, to)`. A cumulative series cannot fall on
its own, so each entry is a revision: the report has moved counts
between zones, provinces or the unallocated row. `min_fall` guards
against a single-unit correction being read as a revision.

[`zone_increment_matrix`](@ref) clamps a negative increment to zero, so
without this the revision is absorbed silently. Report these rather than
let the clamp hide them.
"""
function zone_cumulative_falls(zone_history; min_fall::Integer = 1)
    out = @NamedTuple{
        province::String, zone::String, day::Int,
        from::Int, to::Int,
    }[]
    for (prov, zones) in zone_history, (zone, h) in zones

        zone == "unallocated" && continue
        for i in 2:length(h.days)
            h.counts[i - 1] - h.counts[i] > min_fall || continue
            push!(
                out,
                (;
                    province = String(prov), zone = String(zone),
                    day = h.days[i], from = h.counts[i - 1], to = h.counts[i],
                )
            )
        end
    end
    return sort!(out; by = x -> (x.day, x.province, x.zone))
end

"""
    zone_reattribution_days(zone_history; include_zone_falls = false,
        min_fall = 1)

The vintage days on which a province's `unallocated` cumulative count
falls, keyed by province, from the per-health-zone histories loaded by
[`load_observations`](@ref). A fall means the report has attributed
cases (or deaths) it had carried as unallocated to named zones, so on
that day the zones' cumulative counts rise by more than the province's.
A province whose unallocated row never falls, or that has none, is
absent. [`zone_increment_matrix`](@ref) leaves those vintages out of the
composition.

A revision can also move counts the other way, out of named zones,
leaving the unallocated row flat or rising. That vintage is a revision
just as much, but the unallocated rule does not see it. With
`include_zone_falls` the vintages of [`zone_cumulative_falls`](@ref) are
added, so any vintage on which a named zone loses more than `min_fall`
is left out too. It is off by default: the confirmed-case composition
was fitted under the unallocated rule alone, and widening it there
changes that stream rather than this one.
"""
function zone_reattribution_days(
        zone_history;
        include_zone_falls::Bool = false, min_fall::Integer = 1
    )
    out = Dict{String, Vector{Int}}()
    for (prov, zones) in zone_history
        haskey(zones, "unallocated") || continue
        h = zones["unallocated"]
        falls = [
            h.days[i] for i in 2:length(h.days)
                if h.counts[i] < h.counts[i - 1]
        ]
        isempty(falls) || (out[String(prov)] = falls)
    end
    if include_zone_falls
        for f in zone_cumulative_falls(zone_history; min_fall)
            push!(get!(out, f.province, Int[]), f.day)
        end
        for (prov, days) in out
            out[prov] = sort!(unique!(days))
        end
    end
    return out
end

"""
    zone_increment_matrix(zone_history, patch_names, members = PROVINCE_MEMBERS;
        reattribution = zone_reattribution_days(zone_history))

Reshape the per-health-zone cumulative histories loaded by
[`load_observations`](@ref) into one increment matrix per patch, the
within-patch analogue of [`province_increment_matrix`](@ref).

Each patch pools the source provinces `members` gives it (a name with no
entry is its own province), and its matrix has one row per zone of those
provinces, the `unallocated` rows left out, and one column per vintage.
Every zone must be reported on the same vintage days, which the tables
guarantee by sharing one `dates` array; a mismatch is an error rather
than a silent reshape. Increments are the differences of consecutive
cumulative counts, the first from zero, clamped at zero as for the
provinces since a zone's cumulative can fall when cases are reattributed.

A vintage on which a member province's unallocated count falls
(`reattribution`, the days per province of
[`zone_reattribution_days`](@ref)) is a reattribution of counts the
report had carried as unallocated into named zones. The zone increments
would read them as new cases, so that column is set to zero for the
patch, which drops the cell from the composition, and its days are
returned as `excluded`. The next vintage's increment is still the
difference of the cumulative counts, so nothing is counted twice. Pass
the merged days of more than one history (the death block's falls as
well as the case block's) to exclude the union, or an empty `Dict` to
keep every vintage.

Returns a vector with one named tuple per patch: `patch` (its name),
`zones` (a vector of `(province, zone)` key pairs in row order), `days`,
`increments` (the `(n_zones × n_vintages)` matrix), `totals` (the column
sums, the allocated count each vintage's composition conditions on) and
`excluded` (the zeroed vintage days). A patch none of whose provinces
has zone data gets an empty matrix. An empty `zone_history` returns an
empty vector.
"""
function zone_increment_matrix(
        zone_history, patch_names::AbstractVector,
        members::AbstractDict = PROVINCE_MEMBERS;
        reattribution::AbstractDict = zone_reattribution_days(zone_history)
    )
    out = @NamedTuple{
        patch::String, zones::Vector{Tuple{String, String}},
        days::Vector{Int}, increments::Matrix{Int}, totals::Vector{Int},
        excluded::Vector{Int},
    }[]
    isempty(zone_history) && return out
    days = nothing
    for nm in patch_names
        provs = get(members, nm, [nm])
        keys_ = Tuple{String, String}[]
        for prov in provs
            haskey(zone_history, prov) || continue
            for zone in sort!(collect(keys(zone_history[prov])))
                zone == "unallocated" && continue
                push!(keys_, (prov, zone))
            end
        end
        if isempty(keys_)
            push!(
                out,
                (;
                    patch = String(nm), zones = keys_, days = Int[],
                    increments = Matrix{Int}(undef, 0, 0), totals = Int[],
                    excluded = Int[],
                )
            )
            continue
        end
        for (prov, zone) in keys_
            h = zone_history[prov][zone]
            days === nothing && (days = h.days)
            h.days == days || error(
                "zone `$(prov).$(zone)` is reported on different vintage " *
                    "days to `$(keys_[1][1]).$(keys_[1][2])`; the zone " *
                    "composition needs every zone on the same vintages."
            )
        end
        inc = Matrix{Int}(undef, length(keys_), length(days))
        for (i, (prov, zone)) in enumerate(keys_)
            c = zone_history[prov][zone].counts
            inc[i, :] = max.(diff(vcat(0, c)), 0)
        end
        excluded = sort!(
            unique!(
                Int[
                    d
                        for prov in provs
                        for d in get(reattribution, prov, Int[])
                        if d in days
                ]
            )
        )
        for d in excluded
            inc[:, findfirst(==(d), days)] .= 0
        end
        push!(
            out,
            (;
                patch = String(nm), zones = keys_, days = copy(days),
                increments = inc, totals = vec(sum(inc; dims = 1)),
                excluded,
            )
        )
    end
    return out
end

"""
    load_health_zones(path = data/health_zones.csv)

Read the health-zone metadata table: one named tuple per zone with the
manifest key `zone`, the display `label`, the `province` key, the
WorldPop `population`, the polygon centroid `lat` and `lon` in decimal
degrees, and the DHIS2 `zscode` from the health-zone shapefile. Rows are
in patch order then alphabetical by key. The file is written by
`scripts/build_health_zones.py`; see `data/README.md` for its sources.
"""
function load_health_zones(
        path::AbstractString = joinpath(
            @__DIR__, "..", "data",
            "health_zones.csv"
        )
    )
    Row = @NamedTuple{
        zone::String, label::String, province::String,
        population::Int, lat::Float64, lon::Float64, zscode::String,
    }
    rows = Row[]
    lines = readlines(path)
    header = split(lines[1], ',')
    header == [
        "zone", "label", "province", "population", "lat", "lon",
        "zscode",
    ] || error("unexpected header in $(path): $(header)")
    for line in lines[2:end]
        isempty(strip(line)) && continue
        f = split(line, ',')
        length(f) == 7 || error("expected 7 fields in $(path): $(line)")
        push!(
            rows,
            (;
                zone = String(f[1]), label = String(f[2]),
                province = String(f[3]), population = parse(Int, f[4]),
                lat = parse(Float64, f[5]), lon = parse(Float64, f[6]),
                zscode = String(f[7]),
            )
        )
    end
    return rows
end

"""
    province_testing_covariate(province_lab_daily_history, province_names,
                               populations)

Per-capita laboratory effort in each patch, logged and centred to mean
zero, for the covariate on the prior for relative case ascertainment in
[`province_composition_model`](@ref).

Sums each patch's `<province>_analysed` daily counts over the whole
laboratory window, pooling the source provinces a patch covers (see
[`PROVINCE_MEMBERS`](@ref)), and divides by `populations[p]`. Centring
matches the sum-to-zero ascertainment the composition identifies.

Returns a length-`n_patches` vector of zeros when the laboratory history
is absent, when a patch has no analysed series, or when a patch analysed
nothing over the window. A zero covariate recovers the model without it.
"""
function province_testing_covariate(
        province_lab_daily_history,
        province_names::AbstractVector = PROVINCE_NAMES,
        populations::AbstractVector{<:Real} = PROVINCE_POPULATIONS
    )
    np = length(province_names)
    length(populations) == np || error(
        "province_testing_covariate: $(length(populations)) populations " *
            "for $(np) patches."
    )
    none = zeros(np)
    isempty(province_lab_daily_history) && return none
    members = patch_members(province_names)
    series = [["$(m)_analysed" for m in ms] for ms in members]
    any(ks -> any(k -> !haskey(province_lab_daily_history, k), ks), series) &&
        return none
    analysed = [
        sum(sum(province_lab_daily_history[k].counts) for k in ks)
            for ks in series
    ]
    any(iszero, analysed) && return none
    log_rate = log.(analysed ./ populations)
    return log_rate .- (sum(log_rate) / np)
end

"""
    freeze_observations(cutoff_date; path = default manifest)

Load the observation manifest frozen to `cutoff_date`. Every dated
history is truncated to the vintages available by then, so the returned
named tuple is what the renewal model would have seen on that date. The
cut-off scalar totals are taken from the truncated histories rather than
the manifest's full-data scalars. Use it to re-evaluate the renewal
estimate at a past report date for a matched-in-time comparison.

`cutoff_date` accepts a `Date` or an ISO date string. It must be on or
after the earliest history vintage in the manifest, since the DRC series
begins 18 May 2026. An earlier date leaves the suspected streams empty
and is not a meaningful renewal fit.
"""
function freeze_observations(
        cutoff_date::Union{Date, AbstractString};
        path::AbstractString = joinpath(
            @__DIR__, "..", "data",
            "observations.toml"
        ),
        seeding_lead::Integer = SEEDING_LEAD_DAYS
    )
    return load_observations(path; seeding_lead, cutoff_date)
end

"""
One entry per observation stream the manifest loader carries, in the order
the report presents them. Each entry is a `NamedTuple`:

- `id`: the canonical short identifier used across the package.
- `field`: the field of a loaded observation set holding the stream's
  dated history.
- `label`: the display title the report uses for the stream.
- `score_label`: the label the forecast archive and the release scoring
  table use, or `nothing` for an unscored stream.
- `forecast_prefix`: the stem of the stream's forecast columns (`:cases`
  for `cases_cum` and `cases_new`), or `nothing` for a stream that is not
  forecast.

This is the single list of streams, so a stream is named once rather than
once per consumer. [`stream_id`](@ref) resolves any of the four
vocabularies back to `id`.
"""
const OBSERVATION_STREAMS = (
    (;
        id = :suspected_cases, field = :reported_history,
        label = "Suspected cases", score_label = "reported cases",
        forecast_prefix = :cases,
    ),
    (;
        id = :suspected_deaths, field = :deaths_history,
        label = "Suspected deaths", score_label = "suspected deaths",
        forecast_prefix = :deaths,
    ),
    (;
        id = :suspected_daily, field = :suspected_daily_history,
        label = "New suspects/day", score_label = nothing,
        forecast_prefix = nothing,
    ),
    (;
        id = :suspected_daily_deaths, field = :suspected_daily_deaths_history,
        label = "New suspected deaths/day", score_label = nothing,
        forecast_prefix = nothing,
    ),
    (;
        id = :confirmed_cases, field = :confirmed_history,
        label = "Confirmed cases", score_label = "confirmed cases",
        forecast_prefix = :confirmed,
    ),
    (;
        id = :confirmed_deaths, field = :confirmed_deaths_history,
        label = "Confirmed deaths", score_label = "confirmed deaths",
        forecast_prefix = :confirmed_deaths,
    ),
    (;
        id = :recovered, field = :recovered_history,
        label = "Recovered (confirmed)", score_label = "recovered",
        forecast_prefix = :recovered,
    ),
    (;
        id = :tests_analysed, field = :lab_history,
        label = "Specimens analysed (cumulative)", score_label = nothing,
        forecast_prefix = nothing,
    ),
    (;
        id = :tests_analysed_daily, field = :lab_daily_history,
        label = "Specimens analysed (24h)", score_label = nothing,
        forecast_prefix = nothing,
    ),
    (;
        id = :tests_received, field = :tests_received_history,
        label = "Specimens received", score_label = nothing,
        forecast_prefix = nothing,
    ),
    (;
        id = :isolation_beds, field = :isolation_history,
        label = "Patients in isolation", score_label = "isolation beds",
        forecast_prefix = nothing,
    ),
    (;
        id = :bed_capacity, field = :bed_capacity_history,
        label = "Bed capacity", score_label = nothing,
        forecast_prefix = nothing,
    ),
    (;
        id = :treatment_admissions, field = :treatment_admissions_history,
        label = "Admissions/day", score_label = nothing,
        forecast_prefix = nothing,
    ),
    (;
        id = :treatment_deaths, field = :treatment_deaths_history,
        label = "In-care deaths/day", score_label = nothing,
        forecast_prefix = nothing,
    ),
    (;
        id = :treatment_ruleouts, field = :treatment_ruleout_history,
        label = "Rule-outs/day", score_label = nothing,
        forecast_prefix = nothing,
    ),
    (;
        id = :treatment_absconded, field = :treatment_absconded_history,
        label = "Absconded/day", score_label = nothing,
        forecast_prefix = nothing,
    ),
    (;
        id = :treatment_beds, field = :treatment_confirmed_incare_history,
        label = "Confirmed in care", score_label = "treatment beds",
        forecast_prefix = nothing,
    ),
    (;
        id = :suspect_beds, field = :treatment_suspect_incare_history,
        label = "Suspects in care",
        score_label = "isolation beds (suspected)",
        forecast_prefix = nothing,
    ),
    (;
        id = :onset_reports, field = :onset_report_history,
        label = "Onset reports", score_label = "onset reports",
        forecast_prefix = nothing,
    ),
    (;
        id = :exports, field = :export_case_days,
        label = "Uganda exports", score_label = "exports",
        forecast_prefix = nothing,
    ),
)

"""
Calendar date of a grid day-index, where day `obs.n` is the cut-off and
day 1 is the seeding day. Every dated history stores day-indices on this
grid, so this is how a vintage is read back as a date.
"""
grid_date(obs, day::Integer)::Date = obs.cutoff - Day(obs.n - day)

"""
Resolve `stream` to its canonical identifier in [`OBSERVATION_STREAMS`](@ref).
Accepts a canonical identifier, an observation-set history field name, a
scoring label, or a forecast column name (`:cases_cum` and `:cases_new`
both resolve to `:suspected_cases`). The four vocabularies are disjoint,
so the mapping is unambiguous. Errors on an unknown stream, naming the
identifiers it knows.
"""
function stream_id(stream)::Symbol
    for e in OBSERVATION_STREAMS
        stream === e.id && return e.id
        stream === e.field && return e.id
        !isnothing(e.score_label) && stream == e.score_label && return e.id
        if !isnothing(e.forecast_prefix)
            for suffix in ("_cum", "_new")
                stream === Symbol(e.forecast_prefix, suffix) && return e.id
            end
        end
    end
    known = join(string.(getfield.(OBSERVATION_STREAMS, :id)), ", ")
    return error("unknown stream '$stream'. Known streams: $known.")
end

## Registry entry for a canonical id, already resolved by `stream_id`.
function _stream_entry(id::Symbol)
    for e in OBSERVATION_STREAMS
        e.id === id && return e
    end
    return error("no registry entry for stream id '$id'.")
end

"""
Forecast column names for `stream` as `(; cum, new)`, or `nothing` for a
stream the forecast does not carry. Built from the registry's
`forecast_prefix`, so the column names are derived rather than repeated.
"""
function stream_forecast_columns(stream)
    e = _stream_entry(stream_id(stream))
    isnothing(e.forecast_prefix) && return nothing
    return (;
        cum = Symbol(e.forecast_prefix, "_cum"),
        new = Symbol(e.forecast_prefix, "_new"),
    )
end

"""
Date `stream` was last reported in `obs`, or `missing` when `obs` does
not carry the stream or the stream has no vintages. This is the date of
the stream's last vintage, since past it the series is only ever repeated
at its last reported value rather than genuinely observed.

The Uganda exports are a dated list of detections rather than a series of
vintages, so their last reported date is the later of the last detected
import and the last detected import death.
"""
function stream_last_date(obs, stream)::Union{Date, Missing}
    id = stream_id(stream)
    if id === :exports
        days = Int[]
        for f in (:export_case_days, :export_death_days)
            hasproperty(obs, f) && append!(days, Int.(getproperty(obs, f)))
        end
        isempty(days) && return missing
        return grid_date(obs, maximum(days))
    end
    field = _stream_entry(id).field
    hasproperty(obs, field) || return missing
    h = getproperty(obs, field)
    isempty(h.days) && return missing
    return grid_date(obs, maximum(h.days))
end

"""
Days a stream's last vintage may lag the cut-off and still count as
reporting. One reporting week, which is also the forecast horizon, so a
stream that skips a single situation report is not read as stopped.
"""
const STREAM_REPORTING_GRACE_DAYS = 7

"""
Whether `stream` was still being reported at the cut-off of `obs`, that
is whether its last vintage falls within `grace` days of the cut-off. A
stream `obs` does not carry, or one with no vintages, is not reporting.
"""
function stream_reporting(
        obs, stream;
        grace::Integer = STREAM_REPORTING_GRACE_DAYS
    )::Bool
    d = stream_last_date(obs, stream)
    return !ismissing(d) && (obs.cutoff - d) <= Day(grace)
end

"""
Reporting status of every stream `obs` carries, one row per entry of
[`OBSERVATION_STREAMS`](@ref) and in that order. Columns: the canonical
`stream` identifier, its display `label`, the `last_date` it reported,
whether it is still `reporting` at the cut-off, and `days_since` that
last report.
"""
function stream_report_status(
        obs;
        grace::Integer = STREAM_REPORTING_GRACE_DAYS
    )
    entries = [e for e in OBSERVATION_STREAMS if hasproperty(obs, e.field)]
    last_dates = Union{Date, Missing}[
        stream_last_date(obs, e.id)
            for e in entries
    ]
    ## Columns are built typed rather than from a row vector, so a stream
    ## with no vintages (a `missing` date) cannot widen the whole table to
    ## `Any` and leave the reporting flag unusable as an index.
    return DataFrame(
        stream = Symbol[e.id for e in entries],
        label = String[e.label for e in entries],
        last_date = last_dates,
        reporting = Bool[
            stream_reporting(obs, e.id; grace)
                for e in entries
        ],
        days_since = Union{Int, Missing}[
            ismissing(d) ? missing :
                (obs.cutoff - d).value
                for d in last_dates
        ]
    )
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

For the integral backfill only, where `m` counts doublings over the whole
outbreak and `2^m` is the cut-off cumulative case total, so a base of 9
matches McCabe et al.'s Method 2 central 501 cases.

Not for a renewal fit. There `m` counts the cryptic generations and the
seed is the daily incidence reached over them, so an advancing
outbreak-size centre would give a seed of order half a million per day.
`exponential_growth_model` carries its own default.
"""
function m_prior_centre(
        as_of_date::Union{Date, AbstractString};
        base_date::AbstractString = M_PRIOR_BASE_DATE,
        m_base::Real = M_PRIOR_BASE,
        doubling_days::Real = M_PRIOR_DOUBLING_DAYS
    )
    as_of = as_of_date isa Date ? as_of_date : Date(String(as_of_date))
    elapsed = date2epochdays(as_of) - date2epochdays(Date(base_date))
    return m_base + elapsed / doubling_days
end
