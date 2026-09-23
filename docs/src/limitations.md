# Limitations

## Overall

- The model is fitted to aggregate counts, with the DRC data national and per-province situation-report totals read from the published PDFs, and the Uganda data three export cases with one death from WHO reports.
  We do not have a line list, information on case definitions or reporting completeness, or access to the response's internal data, so every estimate is a model-based extrapolation under strong assumptions, not a measurement.
- The work is done outside the outbreak response, so changes in reporting practice are inferred from the reports rather than known.
- The code and analysis were drafted by a language model and reviewed by people and by language models ([authors](@ref "Authors, funding and acknowledgements")).
- The report is re-run as new data arrive, so the estimates change between updates, and it has not been peer reviewed.
  Every release is signed off by a person before it is published.

## National

### Data

- Most quantities, including nearly all of the delays, the case-fatality ratio and the laboratory assumptions, rest on weakly-informed priors, often from other outbreaks, and their posteriors largely track those priors.
  The onset-to-report delay is the exception, informed by the [onset curve](@ref "Symptom-onset reporting delay").
- Almost every count is report-dated, so the epidemic's timing is recovered mainly through the assumed delays.
  The digitised onset curve, covering confirmed cases from SitRep 059 onward, is the only series carrying symptom-onset dates.
- The suspected streams are no longer published, so the late window rests on the confirmed, laboratory, treatment-centre and onset data, and the suspected forecasts cannot be checked.
  The laboratory analysed-specimen series covers only part of the window ([laboratory pipeline](@ref "Laboratory pipeline")).
- Later situation reports can revise earlier totals up or down as suspects are reclassified and newly-reporting health zones are added, and we do not model this revision process.
- The onset curve is digitised from a figure that each report redraws, so the counts carry scan error that the model treats as noise at a fixed scale.

### Model

- The streams observe overlapping people but are fitted as conditionally independent given latent incidence, which can understate uncertainty.
  The [outbreak size estimated by each data stream](@ref "Outbreak size estimated by each data stream") checks whether they agree.
- The DRC ascertainment and the testing fraction are each one value over the window, although ascertainment probably rose, so a change in either is read as a change in incidence.
- The model inherits McCabe et al.'s assumptions of a single zoonotic seed, a generation interval from earlier Ebola outbreaks and no depletion of susceptibles, which help constrain the estimated outbreak age and early reproduction number.
  The onset-to-death delay and the [genetic seeding bound](@ref "Genetic bound on outbreak age") do not propagate cross-outbreak or clock uncertainty.
- The intervention ramp in the [reproduction number](@ref "Reproduction number") is fixed in time, with an assumed start and three-week duration, and only its size is estimated.
- Occupancy shows demand only up to the beds filled, so the bed shortfall above capacity is not measured and comes from the [treatment-centre flow](@ref "Treatment-centre flow") model and its priors.

### Evaluation

- The [onset forecast](@ref "Symptom-onset nowcast and forecast") mostly measures scan error, since its interval is dominated by the per-scan error rather than the epidemic.
- The persistence baseline for a frozen re-fit reads a data snapshot that can post-date the forecast by weeks, so it can see later corrections ([details](@ref "Forecast scoring against a persistence baseline")).
- Past forecasts that were not stored are rebuilt from each release's own code but with current dependency versions, so they are not exact.

## Provinces

### Data

- Confirmed cases and deaths are the only province data fitted.
  Occupancy and beds by province are not, and there is no province onset nowcast.
- **The province tables cover part of the window.** The province vintages start on 15 June and stop before the cut-off, and harmonisation backfill is published only nationally.
  Outside those vintages the provinces are informed only through the national streams ([province compositions](@ref province-compositions)).

### Model

- Ituri, Nord-Kivu and Haut-Uele are modelled individually and every other affected province is pooled into a fourth patch, each well mixed, so spread inside a province is not represented.
- There is no mobility or origin-destination data for this outbreak, so the [gravity kernel](@ref "Mixing and importation") coupling the provinces is assumed, not measured.
- Provincial testing enters the prior, not the likelihood, so the deaths do most of the work in separating a province's incidence from its case-finding and the prior strongly influences the rest ([province parameters](@ref "Province parameters against their priors")).
- The treatment-centre model carries one national bed capacity and demand, so it cannot represent local saturation.
  Ituri holds most of the occupied beds, so the national bed shortfall understates local unmet need.

### Evaluation

- No province projection has been scored yet, since only [province forecast](@ref "Province forecast") projections are scored and earlier fixed-share archives are not.
- Each province is projected on its own, without the fitted cross-province correlation and with the pooled patch forecast as one, so the province forecasts need not sum to the national one.
