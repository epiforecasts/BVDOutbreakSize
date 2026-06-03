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
"""
function load_observations(
        path::AbstractString = joinpath(@__DIR__, "..", "data",
        "observations.toml"))
    raw = TOML.parsefile(path)
    _val(k) = raw[k]["value"]
    _src(k) = String(raw[k]["source"])
    as_of = String(raw["as_of_date"])
    _gap(d) = date2epochdays(Date(as_of)) - date2epochdays(Date(String(d)))
    ## Days between a recorded event date and the cut-off, used as the
    ## elapsed-time offset for the timing terms. A scalar date gives a
    ## `missing` offset when absent (so its term is a no-op).
    _delta(k) = haskey(raw, k) ? _gap(_val(k)) : missing
    ## Daily export-death series, earliest dated death (index 1) to the
    ## cut-off day (offset 0, kept); empty when no dates are present.
    export_deaths_daily = if haskey(raw, "export_death_dates")
        offs = Int[_gap(d) for d in _val("export_death_dates")]
        isempty(offs) ? Int[] :
        Int[count(==(δ), offs) for δ in maximum(offs):-1:0]
    else
        Int[]
    end
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
        offs = Int[_gap(d) for d in ds]
        vs = Int.(block["values"])
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
    _scalar(k, h) = haskey(raw, k) ? Int(_val(k)) : _hist_end(h)
    _scalar_src(k, hk) = haskey(raw, k) ? _src(k) :
                         (haskey(raw, hk) ? String(raw[hk]["source"]) : missing)
    has_gen = haskey(raw, "genetic_tmrca")
    return (;
        as_of_date = as_of,
        exported_cases = Int(_val("exported_cases")),
        exports_deaths = Int(_val("exports_deaths")),
        total_deaths = Int(_val("total_deaths")),
        reported_cases = Int(_val("reported_cases")),
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
## by one doubling per `M_PRIOR_DOUBLING_DAYS` (the central 14-day
## doubling time) of elapsed time to the cut-off.
const M_PRIOR_BASE_DATE = "2026-05-18"
const M_PRIOR_DOUBLING_DAYS = 14.0
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
501 cases ⇒ `m ≈ 9`), advancing at the central 14-day doubling time, so
the prior stays centred on the plausible outbreak size as the cut-off
moves — it tracks data refreshes without manual edits, and a McCabe-date
fit recovers the base value. `C(T) = 2^m` is now the cumulative
*infection* count; 9 is kept as a weakly-informative centre of the same
order rather than rescaled to a case-equivalent infection count.
"""
function m_prior_centre(as_of_date::AbstractString;
        base_date::AbstractString = M_PRIOR_BASE_DATE,
        m_base::Real = M_PRIOR_BASE,
        doubling_days::Real = M_PRIOR_DOUBLING_DAYS)
    elapsed = date2epochdays(Date(as_of_date)) -
              date2epochdays(Date(base_date))
    return m_base + elapsed / doubling_days
end
