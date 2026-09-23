# Limitations

## Overall

- **Fitted to aggregate counts.** The DRC data are national and per-province situation-report totals, read from the published PDFs, and the Uganda data are three export cases with one death from WHO reports.
  We do not have a line list, information on case definitions or reporting completeness, or access to the response's internal data.
  Every estimate is a model-based extrapolation under strong assumptions, not a measurement.
- **An external view.** The work is done outside the outbreak response.
  Changes in reporting practice are inferred from the reports rather than known.
- **LLM-driven implementation.** The code and analysis were drafted by a language model and reviewed by people and by language models ([authors](@ref "Authors, funding and acknowledgements")).
- **Live and not peer reviewed.** The report is re-run as new data arrive, so the estimates change between updates, and it has not been peer reviewed.
  Every release is signed off by a person before it is published.

## National

### Data

- **Most quantities rest on weakly-informed priors.** Nearly all of the delays, the case-fatality ratio and the laboratory assumptions rest on priors, often from other outbreaks.
  The data do little to move them, so these posteriors largely track their priors.
  The onset-to-report delay is the exception, informed by the [onset curve](@ref "Symptom-onset reporting delay").
- **Almost every count is report-dated.** The digitised onset curve is the only series carrying symptom-onset dates, and it covers confirmed cases from SitRep 059 onward.
  Everything else is a total at the report date, so the epidemic's timing is recovered mainly through the assumed delays.
- **The suspected streams are no longer published.** The late window rests on the confirmed, laboratory, treatment-centre and onset data, and the suspected forecasts cannot be checked.
  The laboratory analysed-specimen series covers only part of the window ([laboratory pipeline](@ref "Laboratory pipeline")).
- **Later sitreps revise earlier figures.** A later situation report can revise an earlier total up or down as suspects are reclassified and newly-reporting health zones are added.
  We do not model this revision process.
- **The onset curve is read from a figure.** Each bar is digitised from a figure that each report redraws, so the counts carry scan error that the model treats as noise at a fixed scale ([#793](https://github.com/epiforecasts/BVDOutbreakSize/pull/793), [#824](https://github.com/epiforecasts/BVDOutbreakSize/issues/824)).

### Model

- **Streams share one case pool.** They are fitted as conditionally independent given latent incidence but observe overlapping people.
  This can understate uncertainty ([#307](https://github.com/epiforecasts/BVDOutbreakSize/issues/307)), which the [outbreak size estimated by each data stream](@ref "Outbreak size estimated by each data stream") checks.
- **Ascertainment and testing are constant.** The DRC ascertainment and the testing fraction are each one value over the window, although ascertainment probably rose, so a change in either is read as a change in incidence ([#400](https://github.com/epiforecasts/BVDOutbreakSize/pull/400), [#546](https://github.com/epiforecasts/BVDOutbreakSize/issues/546)).
- **Inherits McCabe et al.'s epidemiological assumptions.** A single zoonotic seed, a generation interval from earlier Ebola outbreaks, and no depletion of susceptibles.
  The onset-to-death delay and the [genetic seeding bound](@ref "Genetic bound on outbreak age") do not propagate cross-outbreak or clock uncertainty.
  These help constrain the estimated outbreak age and early reproduction number.
- **Intervention ramp is fixed in time.** The ramp in the [reproduction number](@ref "Reproduction number") has an assumed start and three-week duration, and only its size is estimated.
- **The bed shortfall is not measured.** Occupancy shows demand only up to the beds filled, so the shortfall above capacity comes from the [treatment-centre flow](@ref "Treatment-centre flow") model and its priors ([#640](https://github.com/epiforecasts/BVDOutbreakSize/issues/640)).

### Evaluation

- **The onset forecast mostly measures scan error.** The interval on the [onset forecast](@ref "Symptom-onset nowcast and forecast") is dominated by the per-scan error rather than the epidemic.
- **Frozen-fit baselines can see later corrections.** The persistence baseline for a frozen re-fit reads a data snapshot that can post-date the forecast by weeks ([details](@ref "Forecast scoring against a persistence baseline")).
- **Reconstructed forecasts are not exact.** Past forecasts that were not stored are rebuilt from each release's own code, but with current dependency versions ([#622](https://github.com/epiforecasts/BVDOutbreakSize/issues/622)).

## Provinces

### Data

- **Province data are confirmed cases and deaths only.** Occupancy and beds by province are not fitted ([#784](https://github.com/epiforecasts/BVDOutbreakSize/pull/784)), and there is no province onset nowcast.
- **The province tables cover part of the window.** The province vintages start on 15 June and stop before the cut-off, and harmonisation backfill is published only nationally.
  Outside those vintages the provinces are informed only through the national streams ([province compositions](@ref province-compositions)).

### Model

- **Four patches, not the full provincial detail.** Ituri, Nord-Kivu and Haut-Uele are modelled individually and every other affected province is pooled into a fourth patch ([#779](https://github.com/epiforecasts/BVDOutbreakSize/pull/779)).
  Transmission within a patch is well mixed, so spread inside a province is not represented.
- **Importation structure is assumed, not measured.** There is no mobility or origin-destination data for this outbreak, so the [gravity kernel](@ref "Mixing and importation") is a structural assumption.
- **Provincial testing enters the prior, not the likelihood.** The deaths do most of the work in separating a province's incidence from its case-finding, and the prior strongly influences the rest ([province parameters](@ref "Province parameters against their priors"), [#784](https://github.com/epiforecasts/BVDOutbreakSize/pull/784)).
- **Single national bed capacity.** The treatment-centre model carries one national bed capacity and one national demand, so it cannot represent local saturation ([#784](https://github.com/epiforecasts/BVDOutbreakSize/pull/784)).
  Ituri holds most of the occupied beds, so the national bed shortfall understates local unmet need.

### Evaluation

- **No province projection has been scored yet.** Only [province forecast](@ref "Province forecast") projections are scored, and earlier fixed-share archives are not.
- **Province forecasts do not sum to the national one.** Each province is projected on its own, without the fitted cross-province correlation, and the pooled patch is forecast as one ([#831](https://github.com/epiforecasts/BVDOutbreakSize/issues/831)).
