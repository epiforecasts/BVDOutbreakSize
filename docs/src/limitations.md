# Limitations

## Overall

- The model is fitted to aggregate counts: the DRC situation-report totals read from the published PDFs, and three Uganda export cases with one death from WHO reports.
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
- The cumulative suspected streams and the daily suspected deaths are no longer published, so the late window rests on the daily new suspects and the confirmed, laboratory, treatment-centre and onset data, and the cumulative suspected forecasts cannot be checked.
  The laboratory analysed-specimen series covers only part of the window ([laboratory pipeline](@ref "Laboratory pipeline")).
- Later situation reports can revise earlier totals up or down as suspects are reclassified and newly-reporting health zones are added, and we do not model this revision process.
- The onset curve is digitised from a figure that each report redraws, so the counts carry scan error that the model treats as noise with one fitted scale shared by every vintage.
- The reports also print surveillance-effort indicators, such as contact-tracing follow-up and alert investigation, and community suspect-death alerts that continue after 11 July.
  We record these signals but do not fit them.

### Model

- The streams observe overlapping people but are fitted as conditionally independent given latent incidence, which can understate uncertainty.
  The [outbreak size estimated by each data stream](@ref "Outbreak size estimated by each data stream") checks whether they agree.
- The DRC ascertainment, the testing fraction and the specimens analysed per suspect are each one value over the window, although ascertainment probably rose, so a change in any of them is read as a change in incidence.
- The severity-enrichment priors are not set from data.
  Under the decay-scale prior the enrichment is spent within about the first thousand analysed specimens, so it shapes only the earliest laboratory windows.
- The model inherits McCabe et al.'s assumptions of a single zoonotic seed and a generation interval from earlier Ebola outbreaks, which help constrain the estimated outbreak age and early reproduction number.
  The onset-to-death delay does not propagate cross-outbreak uncertainty, and the [genetic seeding bound](@ref "Genetic bound on outbreak age") does not propagate the choice of tree prior.
- The export model counts every border crossing by an undetected infection as a possible detection in Uganda, so a traveller who returns to the DRC before detection is over-counted.
  A constant share of such round trips is absorbed into the Uganda ascertainment, but a share that changes over the window is not.
- The intervention ramp in the [reproduction number](@ref "Reproduction number") is fixed in time, with an assumed start and three-week duration, and only its size is estimated.
- The share of reported suspects admitted to a treatment-centre bed has a weakly-informative prior that no external evidence sets.
- Occupancy shows demand only up to the beds filled, so the bed shortfall above capacity is not measured and comes from the [treatment-centre flow](@ref "Treatment-centre flow") model and its priors.

### Evaluation

- The interval on the [onset forecast](@ref "Symptom-onset nowcast and forecast") is dominated by per-scan error rather than the epidemic.
- The persistence baseline for a frozen re-fit reads a data snapshot that can post-date the forecast by weeks, so it can see later corrections ([details](@ref "Forecast scoring against a persistence baseline")).
- Past forecasts that were not stored are rebuilt from each release's own code but with current dependency versions, so they are not exact.

## Provinces

### Data

- Province data enter as splits of the national or printed totals of confirmed cases, deaths, analysed specimens, patients in isolation, beds and admissions.
  There is no province onset nowcast.
- The province case and death tables start on 15 June, and harmonisation backfill is published only nationally.
  Before then the provinces are informed only through the national streams and the province laboratory split ([province compositions](@ref province-compositions)).

### Model

- Each affected province is one well-mixed patch, so spread inside a province is not represented.
  Sud-Kivu, Tshopo, Bas-Uele, Sud-Ubangi and Nord-Ubangi hold few confirmed cases, so their estimates are informed mostly by the priors their parameters share with the other patches.
- There is no mobility or origin-destination data for this outbreak, so the [gravity kernel](@ref "Mixing and importation") coupling the provinces is assumed, not measured.
- The province laboratory data enter only as the split of analysed specimens, with a national testing fraction, so the deaths do most of the work in separating a province's incidence from its case-finding and the prior strongly influences the rest ([province parameters](@ref "Province parameters against their priors")).
- The treatment-centre model censors bed demand at the national capacity only, and each province's demand and beds are shares of the national ones.
  A province can be over its beds while another has beds free, so the national bed shortfall understates local unmet need.

## Health zones

### Data

- A zone table can lag the cut-off, since a zone's count moves only when a report prints it.
  The [data currency](@ref zone-data-currency) table lists each zone block that stops before the cut-off.

### Model

- A zone below the walking threshold has no transmission walk of its own, so its reproduction number is from its province ([health-zone model](@ref "Health-zone model")).
- The zone stage is melded one way onto the joint fit, so the zone data do not update or correct the national and province estimates.
- The zone forecast splits each patch's forecast across its zones, so it inherits the patch forecast's level and adds only the split.
