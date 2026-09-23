# Limitations

The limitations that apply to the whole work come first, then the national and province ones by data, model and evaluation, with the most consequential first.
Each says what would remove it and links the issue that tracks it, where there is one.
The detail is in the [Methods](@ref "Methods").

## Overall

- **Public data only.** Every count comes from published situation reports and WHO reports, read from the PDFs.
  There is no line list and no access to the response's internal data, so every estimate is a model-based extrapolation from aggregate totals.
  Access to line-list or internal surveillance data would remove it.
- **An external view.** The work is done outside the outbreak response.
  Changes in reporting practice are inferred from the reports rather than known.
  Working with the teams that collect the data would remove it.
- **Drafted and read by a language model.** The code and analysis were drafted by a language model and reviewed by people and by language models ([authors](@ref "Authors, funding and acknowledgements")).
  Language models also read the situation-report figures, and every release is signed off by a person before it is published.
  Errors can remain where no check reaches, since the checks are a second read of each figure, cross-checks of the confirmed totals against the published national series, and tests of the code.
  An independent replication and data read would remove it.
- **Live and not peer reviewed.** The report is re-run as new data arrive, so the estimates change between updates, and it has not been peer reviewed.
  Peer review of a fixed version would remove it.

## National

### Data

- **Delays rest on priors from other outbreaks.** Almost every count is a report-dated aggregate, so the epidemic's timing comes through delays whose estimates rest mainly on priors from earlier outbreaks.
  The case-fatality ratio and the generation interval rest on priors too, and only the onset-to-report delay is informed by the data, through the [onset curve](@ref "Symptom-onset reporting delay").
  A line list from this outbreak would remove it.
- **The suspected streams are no longer published.** The suspected cases and deaths stopped being published, so the late window rests on the confirmed, laboratory, treatment-centre and onset data.
  Their forecasts cannot be checked against later data.
  Resumed publication of the daily suspected counts would remove it.
- **Reported totals are revised and rescanned.** Later reports revise earlier totals, and the onset curve is digitised from a figure that each report redraws.
  The model takes each total as reported and treats rescanned bars as noise at a fixed scale ([#824](https://github.com/epiforecasts/BVDOutbreakSize/issues/824)).
  A revision model and a per-vintage noise scale ([#793](https://github.com/epiforecasts/BVDOutbreakSize/pull/793)) would reduce it.

### Model

- **The streams share one case pool.** The streams count overlapping people but are fitted as conditionally independent, which can understate the uncertainty ([#307](https://github.com/epiforecasts/BVDOutbreakSize/issues/307)).
  The [outbreak size estimated by each data stream](@ref "Outbreak size estimated by each data stream") checks whether the streams agree.
  A likelihood that models the overlap would remove it.
- **Ascertainment and testing are constant.** The DRC ascertainment and the testing fraction are each one value over the window, so a change in either is read as a change in incidence.
  Time-varying ascertainment ([#400](https://github.com/epiforecasts/BVDOutbreakSize/pull/400), [#335](https://github.com/epiforecasts/BVDOutbreakSize/issues/335)) and testing ([#546](https://github.com/epiforecasts/BVDOutbreakSize/issues/546)) would remove it.
- **The outbreak's start rests on fixed assumptions.** The model assumes a single zoonotic seed, no depletion of susceptibles, a fixed molecular clock for the [genetic bound](@ref "Genetic bound on outbreak age") and a fixed timing for the response ramp.
  These help constrain the estimated outbreak age and early reproduction number.
  Estimates from this outbreak's genomes and dated response milestones would remove it.
- **The bed shortfall is not measured.** Occupancy shows demand only up to the beds filled, so the shortfall above capacity comes from the [treatment-centre flow](@ref "Treatment-centre flow") model and its priors ([#640](https://github.com/epiforecasts/BVDOutbreakSize/issues/640)).
  A count of patients waiting for a bed would remove it.

### Evaluation

- **The onset forecast mostly measures scan error.** The interval on the [onset forecast](@ref "Symptom-onset nowcast and forecast") is dominated by the per-scan error rather than the epidemic.
  It checks the fitted delay more than it predicts cases.
  A published onset table would remove it.
- **Frozen-fit baselines can see later corrections.** The persistence baseline for a frozen re-fit reads a data snapshot that can post-date the forecast by weeks, so a later correction can reach it ([details](@ref "Forecast scoring against a persistence baseline")).
  A snapshot archived at each frozen cut-off would remove it.
- **Reconstructed forecasts are not exact.** Past forecasts that were not stored are rebuilt from each release's own code, but with current dependency versions ([#622](https://github.com/epiforecasts/BVDOutbreakSize/issues/622)).
  Their scores are close to, but not exactly, what each release produced.
  Pinned manifests for every release would remove it.

## Provinces

### Data

- **Province data are confirmed cases and deaths only.** The spatial tables give confirmed cases and deaths, the only province data in the likelihood, and there is no province onset nowcast.
  Occupancy, beds and analysed specimens by province are not fitted ([#784](https://github.com/epiforecasts/BVDOutbreakSize/pull/784)).
  Scoring them as splits of their printed sums would narrow it.
- **The province tables cover part of the window.** The province vintages start on 15 June and stop before the cut-off, and harmonisation backfill is published only nationally.
  Outside those vintages the provinces are informed only through the national streams ([province compositions](@ref province-compositions)).
  Province vintages over the whole window would remove it.

### Model

- **Four patches with assumed importation.** Three provinces are modelled on their own and the rest are pooled into a patch with almost no confirmed cases, each well mixed inside.
  The patches are coupled by an assumed [gravity kernel](@ref "Mixing and importation") in the absence of mobility data.
  A health-zone model ([#779](https://github.com/epiforecasts/BVDOutbreakSize/pull/779), [#705](https://github.com/epiforecasts/BVDOutbreakSize/issues/705)) and mobility data would narrow it.
- **Province case-finding rests on the prior.** A province's confirmed share is the product of its incidence and its case-finding, and testing by province enters only the prior.
  The deaths separate the two, and the prior strongly influences where the estimate sits along the ridge ([province parameters](@ref "Province parameters against their priors")).
  Scoring the province analysed volume ([#784](https://github.com/epiforecasts/BVDOutbreakSize/pull/784)) would narrow it.
- **Beds and occupancy are national only.** One national bed capacity and demand cannot show saturation in one province while another has free beds.
  Ituri holds most occupied beds, so the national shortfall understates local unmet need.
  Occupancy and beds by province ([#784](https://github.com/epiforecasts/BVDOutbreakSize/pull/784)) would remove it.

### Evaluation

- **No province projection has been scored yet.** Only [province forecast](@ref "Province forecast") projections are scored, and earlier fixed-share archives are not.
  The scores fill in as projections become old enough to check.
  More releases will remove it.
- **Province forecasts do not sum to the national one.** Each province is projected on its own, without the fitted cross-province correlation, and the pooled patch is forecast as one.
  The provinces need not add up to the national forecast.
  A projection reconciled with the national one ([#831](https://github.com/epiforecasts/BVDOutbreakSize/issues/831)) would remove it.
