# Observation loader: reads `data/observations.toml` and assembles
# the parsed counts, date offsets and citation strings into a
# `NamedTuple` consumed by the analysis pipeline.

"""
Load the observation block from `data/observations.toml` and return
it as a `NamedTuple`. Each observation in the TOML is a subtable
with `value = …` and `source = "…"`; this function returns both the
parsed numeric values and a parallel `sources::NamedTuple` of
citation strings so they can be printed alongside the data.

Fields returned:

- `exported_cases::Int`
- `exports_deaths::Int`
- `total_deaths::Int`
- `reported_cases::Int` — DRC suspected cumulative case count.
- `confirmed_cases::Union{Int, Missing}` — DRC laboratory-confirmed
  cumulative case count at the cut-off, the truth-anchor on the latent
  eventually-confirmable pool ``C(T)`` (reported counts are an inflated
  view). Derived from the final vintage of `confirmed_case_history` (or
  an explicit `confirmed_cases` scalar block if one is present);
  `missing` when neither is present.
- `reported_case_history::Union{NamedTuple, Missing}` — vintage-by-vintage
  cumulative DRC suspected counts, with fields `dates`, `offsets` (days
  before `as_of_date`, sorted ascending) and `values` (cumulative count
  at each sitrep date). Drives the daily reported-cases likelihood by
  per-day differencing. `missing` when no `reported_case_history` block
  is present.
- `confirmed_case_history::Union{NamedTuple, Missing}` — same layout for
  cumulative DRC laboratory-confirmed counts. Drives the daily
  confirmed-cases likelihood. `missing` when absent.
- `confirmed_deaths::Union{Int, Missing}` — DRC cumulative confirmed
  deaths at the cut-off (front-page `Cumul décès parmi les confirmés`).
  Derived from the final vintage of `confirmed_death_history` (or an
  explicit scalar block if present). Recorded for completeness; not
  currently fitted. `missing` when neither is present.
- `confirmed_death_history::Union{NamedTuple, Missing}` — same layout for
  cumulative DRC confirmed deaths. Recorded for completeness; the model
  does not currently fit a confirmed-death stream. `missing` when absent.
- `death_history::Union{NamedTuple, Missing}` — same layout for
  cumulative DRC suspected deaths. Drives the daily deaths likelihood.
  `missing` when absent.
- `cumulative_tests_analysed::Union{Int, Missing}` — cumulative number
  of suspected-case specimens whose lab processing has completed by the
  cut-off. Paired with `confirmed_cases` it gives a per-test positivity
  observation; right-truncation is handled inside the model by the lab
  delay CDF. Derived from the final vintage of `tests_analysed_history`
  (or an explicit scalar block if present); `missing` when neither is
  present.
- `tests_received_history::Union{NamedTuple, Missing}` — same layout as
  the case histories for cumulative specimens received. Recorded for the
  laboratory-pipeline view; not currently fitted. `missing` when absent.
- `tests_analysed_history::Union{NamedTuple, Missing}` — same layout for
  cumulative specimens analysed. Recorded for the laboratory-pipeline
  view; the model uses only the cut-off `cumulative_tests_analysed`
  total. `missing` when absent.
- `daily_outbound_travellers::Real`
- `daily_outbound_travellers_sd::Real`
- `source_population::Int`
- `genetic_tmrca_days::Union{Real, Missing}` — estimated time to the
  most recent common ancestor (TMRCA) in days before `as_of_date`, a
  soft lower bound on the seeding time `T`; `missing` when no
  `genetic_tmrca` block is present.
- `genetic_tmrca_days_sd::Union{Real, Missing}` — SD (days) on the
  location of that floor; `missing` when absent.
- `genetic_tmrca_alt_days::Union{Real, Missing}` — TMRCA (days before
  `as_of_date`) under the alternative clock rate, for the clock-rate
  sensitivity; `missing` when no `alt_date` is present.
- `genetic_tmrca_alt_days_sd::Union{Real, Missing}` — SD (days) on the
  alternative-clock floor; `missing` when absent.
- `sources::NamedTuple` — citation per loaded field, with the same keys
  as the numeric fields above. Optional fields (`confirmed_cases`,
  `genetic_tmrca`) carry `missing` rather than a citation when absent.

The `as_of_override` keyword re-bases the cut-off to an earlier as-of
date for out-of-sample validation. Every vintage history is truncated to
entries on or before the override date, the dated export series keep only
detections by then, all elapsed-time offsets are recomputed relative to
the override, and the derived cut-off scalars (`confirmed_cases`,
`confirmed_deaths`, `cumulative_tests_analysed`, `reported_cases`,
`total_deaths`, `exported_cases`, `exports_deaths`) are taken from the
truncated data so they stay internally consistent with what had been
reported by that date. Pass a `Date` or an ISO `"YYYY-MM-DD"` string;
`nothing` (the default) loads the file's own `as_of_date` unchanged.
"""
function load_observations(
        path::AbstractString = joinpath(@__DIR__, "..", "data",
            "observations.toml");
        as_of_override::Union{AbstractString, Date, Nothing} = nothing)
    raw = TOML.parsefile(path)
    _val(k) = raw[k]["value"]
    _src(k) = String(raw[k]["source"])
    ## `as_of_override` re-bases the cut-off to an earlier as-of date for
    ## out-of-sample validation: every vintage history is truncated to
    ## entries on or before the cut-off, the dated export series keep only
    ## detections by the cut-off, and the derived cut-off scalars
    ## (confirmed cases / deaths, samples analysed, suspected cases /
    ## deaths, exports) are recomputed from the truncated data so they stay
    ## internally consistent. Dates after the cut-off are dropped, not
    ## interpolated.
    as_of = as_of_override === nothing ? String(raw["as_of_date"]) :
            string(Date(as_of_override))
    cutoff_epoch = date2epochdays(Date(as_of))
    truncating = as_of_override !== nothing
    ## Keep a dated entry only when it falls on or before the cut-off.
    _by_cutoff(d) = date2epochdays(Date(String(d))) <= cutoff_epoch
    _gap(d) = cutoff_epoch - date2epochdays(Date(String(d)))
    ## Days between a recorded event date and the cut-off, used as the
    ## elapsed-time offset for the timing terms. A scalar date gives a
    ## `missing` offset when absent (so its term is a no-op), and an event
    ## dated after the cut-off is treated as not-yet-observed (`missing`).
    _delta(k) = haskey(raw, k) && _by_cutoff(_val(k)) ? _gap(_val(k)) : missing
    ## Daily export-death series, earliest dated death (index 1) to the
    ## cut-off day (offset 0, kept); empty when no dates are present.
    export_deaths_daily = if haskey(raw, "export_death_dates")
        ds = filter(_by_cutoff, _val("export_death_dates"))
        offs = Int[_gap(d) for d in ds]
        isempty(offs) ? Int[] :
        Int[count(==(δ), offs) for δ in maximum(offs):-1:0]
    else
        Int[]
    end
    ## Daily Uganda export-case series, earliest detection (index 1) to
    ## the cut-off day (offset 0, kept); empty when no dates are present.
    exported_cases_daily = if haskey(raw, "export_case_dates")
        ds = filter(_by_cutoff, _val("export_case_dates"))
        offs = Int[_gap(d) for d in ds]
        isempty(offs) ? Int[] :
        Int[count(==(δ), offs) for δ in maximum(offs):-1:0]
    else
        Int[]
    end
    ## Cumulative dated-export counts at the cut-off: the number of dated
    ## detections on or before the cut-off. Used to override the bare
    ## `exported_cases` / `exports_deaths` scalars when truncating, so the
    ## export totals match the truncated daily series.
    _dated_count(k) = haskey(raw, k) ? count(_by_cutoff, _val(k)) : nothing
    exported_cases_cut = _dated_count("export_case_dates")
    exports_deaths_cut = _dated_count("export_death_dates")
    ## Cumulative DRC counts at each sitrep vintage: parsed dates,
    ## the elapsed-time offset before the cut-off (days since the
    ## vintage's date, in ascending elapsed-time order) and the
    ## cumulative count. The daily likelihood differences `values`
    ## between consecutive bin edges, so the vector must be monotone
    ## non-decreasing.
    function _history(k)
        haskey(raw, k) || return missing
        block = raw[k]
        ds = String.(block["dates"])
        vs = Int.(block["values"])
        ## Drop vintages after the cut-off when truncating, so the history
        ## reflects only what had been reported by the as-of date.
        if truncating
            keep = _by_cutoff.(ds)
            ds = ds[keep]
            vs = vs[keep]
            isempty(ds) && return missing
        end
        offs = Int[_gap(d) for d in ds]
        ## Sort by ascending elapsed-time (oldest first), so a `diff`
        ## of `values` matches the natural day-by-day increment.
        ord = sortperm(offs; rev = true)
        return (; dates = ds[ord], offsets = offs[ord], values = vs[ord])
    end
    ## Vintage histories are the single source of truth for the confirmed
    ## and laboratory streams: the cut-off scalar each stream needs is the
    ## final (most recent) vintage value. A bare scalar in the TOML still
    ## takes precedence if one is present, for backward compatibility.
    reported_hist = _history("reported_case_history")
    confirmed_hist = _history("confirmed_case_history")
    confirmed_death_hist = _history("confirmed_death_history")
    death_hist = _history("death_history")
    tests_received_hist = _history("tests_received_history")
    tests_analysed_hist = _history("tests_analysed_history")
    _hist_end(h) = h === missing ? missing : h.values[end]
    ## When truncating, the truncated history is the single source of
    ## truth for every cut-off scalar (a bare TOML scalar would carry the
    ## un-truncated value); otherwise a bare scalar still takes precedence.
    _scalar(k, h) = truncating ? _hist_end(h) :
                    (haskey(raw, k) ? Int(_val(k)) : _hist_end(h))
    _scalar_src(k, hk) = haskey(raw, k) ? _src(k) :
                         (haskey(raw, hk) ? String(raw[hk]["source"]) : missing)
    has_gen = haskey(raw, "genetic_tmrca")
    ## Cut-off scalars for the dated-export and frozen suspected streams.
    ## When truncating, prefer the truncated-data view (dated-detection
    ## count for exports, truncated history-end for the frozen suspected
    ## totals); otherwise the bare TOML scalar.
    exported_cases_val = truncating && exported_cases_cut !== nothing ?
                         exported_cases_cut : Int(_val("exported_cases"))
    exports_deaths_val = truncating && exports_deaths_cut !== nothing ?
                         exports_deaths_cut : Int(_val("exports_deaths"))
    total_deaths_val = truncating && death_hist !== missing ?
                       death_hist.values[end] : Int(_val("total_deaths"))
    reported_cases_val = truncating && reported_hist !== missing ?
                         reported_hist.values[end] : Int(_val("reported_cases"))
    return (;
        as_of_date = as_of,
        exported_cases = exported_cases_val,
        exports_deaths = exports_deaths_val,
        total_deaths = total_deaths_val,
        reported_cases = reported_cases_val,
        confirmed_cases = _scalar("confirmed_cases", confirmed_hist),
        confirmed_deaths = _scalar("confirmed_deaths", confirmed_death_hist),
        reported_case_history = reported_hist,
        confirmed_case_history = confirmed_hist,
        confirmed_death_history = confirmed_death_hist,
        death_history = death_hist,
        cumulative_tests_analysed = _scalar("cumulative_tests_analysed",
            tests_analysed_hist),
        tests_received_history = tests_received_hist,
        tests_analysed_history = tests_analysed_hist,
        daily_outbound_travellers = float(
            _val("daily_outbound_travellers")),
        daily_outbound_travellers_sd = float(
            _val("daily_outbound_travellers_sd")),
        source_population = Int(_val("source_population")),
        export_deaths_daily = export_deaths_daily,
        exported_cases_daily = exported_cases_daily,
        first_export_detection_delta = _delta("first_export_detection_date"),
        genetic_tmrca_days = has_gen ?
                             _gap(raw["genetic_tmrca"]["date"]) : missing,
        genetic_tmrca_days_sd = has_gen ?
                                float(raw["genetic_tmrca"]["days_sd"]) : missing,
        genetic_tmrca_alt_days =
        has_gen && haskey(raw["genetic_tmrca"], "alt_date") ?
        _gap(raw["genetic_tmrca"]["alt_date"]) : missing,
        genetic_tmrca_alt_days_sd =
        has_gen && haskey(raw["genetic_tmrca"], "alt_days_sd") ?
        float(raw["genetic_tmrca"]["alt_days_sd"]) : missing,
        sources = (;
            exported_cases = _src("exported_cases"),
            exported_cases_daily = haskey(raw, "export_case_dates") ?
                                   _src("export_case_dates") : missing,
            exports_deaths = _src("exports_deaths"),
            total_deaths = _src("total_deaths"),
            reported_cases = _src("reported_cases"),
            confirmed_cases = _scalar_src("confirmed_cases",
                "confirmed_case_history"),
            confirmed_deaths = _scalar_src("confirmed_deaths",
                "confirmed_death_history"),
            reported_case_history = haskey(raw, "reported_case_history") ?
                                    String(raw["reported_case_history"]["source"]) :
                                    missing,
            confirmed_case_history = haskey(raw, "confirmed_case_history") ?
                                     String(raw["confirmed_case_history"]["source"]) :
                                     missing,
            confirmed_death_history = haskey(raw, "confirmed_death_history") ?
                                      String(raw["confirmed_death_history"]["source"]) :
                                      missing,
            death_history = haskey(raw, "death_history") ?
                            String(raw["death_history"]["source"]) :
                            missing,
            cumulative_tests_analysed = _scalar_src(
                "cumulative_tests_analysed", "tests_analysed_history"),
            tests_received_history = haskey(raw, "tests_received_history") ?
                                     String(raw["tests_received_history"]["source"]) :
                                     missing,
            tests_analysed_history = haskey(raw, "tests_analysed_history") ?
                                     String(raw["tests_analysed_history"]["source"]) :
                                     missing,
            daily_outbound_travellers = _src("daily_outbound_travellers"),
            daily_outbound_travellers_sd = _src("daily_outbound_travellers_sd"),
            source_population = _src("source_population"),
            genetic_tmrca = has_gen ?
                            String(raw["genetic_tmrca"]["source"]) : missing
        )
    )
end

## Doubling-count prior base: McCabe et al.'s first report (18 May
## 2026), whose Method 2 central scenario of 501 cases implies a doubling
## count `log2(501) ≈ 9`. `C(T) = 2^m` is now the cumulative *infection*
## count, but 9 is kept as a weakly-informative centre of the same
## order: the prior only sets where the sampler starts and the fit is
## likelihood-dominated, so the small incubation rescale between
## infections and cases is not worth carrying here. The centre advances
## by one doubling per `M_PRIOR_DOUBLING_DAYS` of elapsed time to the
## cut-off. `M_PRIOR_DOUBLING_DAYS` is the central doubling time, set to
## 20 days from Cuomo-Dannenburg & Ghafari's molecular-clock reanalysis
## (cuomodannenburg2026; mean 15.2-24.5 days across six clock
## assumptions), so the size and growth priors share one central
## doubling time (see [`exponential_growth_model`](@ref)).
const M_PRIOR_BASE_DATE = "2026-05-18"
const M_PRIOR_DOUBLING_DAYS = 20.0
const M_PRIOR_BASE = 9.0

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
501 cases ⇒ `m ≈ 9`), advancing at the central 20-day doubling time
(`M_PRIOR_DOUBLING_DAYS`, the molecular-clock estimate of
cuomodannenburg2026), so the prior stays centred on the plausible outbreak
size as the cut-off moves — it tracks data refreshes without manual edits,
and a McCabe-date fit recovers the base value. `C(T) = 2^m` is now the
cumulative *infection* count; 9 is kept as a weakly-informative centre of
the same order rather than rescaled to a case-equivalent infection count.
"""
function m_prior_centre(as_of_date::AbstractString;
        base_date::AbstractString = M_PRIOR_BASE_DATE,
        m_base::Real = M_PRIOR_BASE,
        doubling_days::Real = M_PRIOR_DOUBLING_DAYS)
    elapsed = date2epochdays(Date(as_of_date)) -
              date2epochdays(Date(base_date))
    return m_base + elapsed / doubling_days
end

"""
    report_onset_offset(as_of_date; base_date)

Days between surveillance/reporting onset (`base_date`, McCabe et al.'s
first report 18 May 2026) and the cut-off `as_of_date`. Used to anchor
the non-BVD background ramp to reporting onset: the background-clock
elapsed time is `t_report = T - report_onset_offset`, so background
suspects accrue only once case-finding has begun rather than from the
latent seeding time. For the 26 May cut-off this is 8 days. Returns zero
when the cut-off is the base date (the McCabe-date configuration), which
collapses the offset and leaves `t_report = T`.
"""
function report_onset_offset(as_of_date::AbstractString;
        base_date::AbstractString = M_PRIOR_BASE_DATE)
    return date2epochdays(Date(as_of_date)) -
           date2epochdays(Date(base_date))
end
