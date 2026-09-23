# Limitations

The limitations are grouped by geography and then by data, model and evaluation, with the most consequential first in each group.
Each says what would remove it, with the issue or pull request that tracks it where there is one.
The detail of each method is in the [Methods](@ref "Methods") section linked from it.

## National

### Data

#### The data identify few of the delays

Most delays, the generation interval and the case-fatality ratio are set by priors from earlier Ebola and Bundibugyo outbreaks, and in places by our own judgement.
The between-report increments inform the change in the reproduction number, but they say little about these quantities on their own, so their posteriors largely track their priors.
The onset curve informs the onset-to-report hazard, with the weak points listed under [symptom-onset reporting delay](@ref "Symptom-onset reporting delay").
A line list from this outbreak, with onset, report, sample and death dates, would remove it.

#### Almost every count is dated by report

The digitised onset curve is the only series dated by symptom onset, and it covers confirmed cases from SitRep 059 onward.
Every other stream is a total at the report date, so the epidemic's timing is recovered mainly through the assumed delays.
A published onset-dated series for the other streams would remove it.

#### The model is fitted to aggregate counts

The DRC data are national and per-province situation-report totals, and the Uganda data are a handful of exported cases and one death.
We have no line list, and no information on case definitions or reporting completeness.
Every estimate is therefore a model-based extrapolation under strong assumptions, not a measurement.
Individual-level data would remove it.

#### The suspected streams are no longer published

The cumulative suspected cases and deaths freeze at 26 May, and the daily suspected counts that replace them are no longer published.
The late window therefore rests on the confirmed, laboratory, treatment-centre and onset data.
Resumed publication of the daily suspected counts would remove it.

#### The laboratory denominator covers part of the window

The analysed-specimen count exists only between the first and last laboratory dates.
Confirmed vintages outside that window are scored against the modelled volume, as the [laboratory pipeline](@ref "Laboratory pipeline") describes.
A published analysed count over the whole window would remove it.

#### Nothing measures the testing of deaths

No death-testing data are published, so the specimens taken per suspected death rest on a tight prior centred on one, as set out under [confirmed deaths](@ref "Confirmed deaths").
A published count of specimens taken from deaths would remove it.

#### Later reports revise earlier totals

A later situation report can revise an earlier total up or down as suspects are reclassified and newly reporting health zones are added.
The model does not represent this revision process.
The forecast scoring removes the harmonisation backfill from the confirmed streams, but the fit sees each revised total as reported.
A revision model fitted to the vintage history would remove it.

#### The onset curve is read from a figure

The onset curve has no data table, so each bar is digitised from a raster figure and every scan rereads the whole figure.
The model carries a fixed pixel-noise scale and a per-scan level error, and absorbs re-dating as noise.
Scan noise has not been constant, and some falls between scans have no identified cause ([#824](https://github.com/epiforecasts/BVDOutbreakSize/issues/824)).
A per-vintage noise scale ([#793](https://github.com/epiforecasts/BVDOutbreakSize/pull/793)) would reduce it, and a published data table would remove it.

### Model

#### The streams share one case pool

The streams are fitted as conditionally independent given the latent incidence, but they count overlapping people.
This can understate the uncertainty ([#307](https://github.com/epiforecasts/BVDOutbreakSize/issues/307)).
The [outbreak size estimated by each data stream](@ref "Outbreak size estimated by each data stream") checks whether each stream on its own implies a consistent size.
A likelihood that models the overlap between streams would remove it.

#### Ascertainment and testing are constant over the window

The DRC ascertainment and the laboratory testing fraction are each one value for the whole window, although ascertainment probably rose as the response grew.
A change in either is read as a change in incidence.
Time-varying ascertainment grounded in contact tracing ([#400](https://github.com/epiforecasts/BVDOutbreakSize/pull/400), [#335](https://github.com/epiforecasts/BVDOutbreakSize/issues/335)) and a time-varying testing fraction ([#546](https://github.com/epiforecasts/BVDOutbreakSize/issues/546)) would remove it.

#### The epidemiological assumptions come from McCabe et al.

The model keeps a single zoonotic seed, a generation interval from earlier Ebola outbreaks and no depletion of susceptibles.
The [generation interval](@ref "Generation interval") and the [onset-to-death delay](@ref "Onset-to-death delay") carry the uncertainty of their sources, but not the difference between those outbreaks and this one.
Delay and generation-interval estimates from this outbreak would remove it.

#### The genetic bound fixes the clock

The [genetic bound on outbreak age](@ref "Genetic bound on outbreak age") fixes the molecular clock at the West African Ebola rate and does not propagate clock uncertainty.
A clock estimated from a longer sampling range of this outbreak's genomes would remove it.

#### The response ramp is fixed in time

The response effect on the [reproduction number](@ref "Reproduction number") follows a logistic ramp centred on the first WHO situation report, with an assumed three-week duration.
Only its size is estimated.
Dated response milestones, or a ramp with sampled timing, would remove it.

#### The export model assumes one-way travel

The [exported cases](@ref "Exported cases") model counts outbound travel only, so round trips before detection are over-counted ([#177](https://github.com/epiforecasts/BVDOutbreakSize/issues/177)).
Exports after the last detection date are not modelled, since the baseline travel rate no longer holds.
Return-travel data, or a model of travel that changes over the outbreak, would remove it.

#### The bed shortfall is not measured

The occupancy shows that demand was at least the beds filled, but not how much more.
The shortfall above a saturated capacity is therefore set by the demand model and its priors, as described under [treatment-centre flow](@ref "Treatment-centre flow"), and sits on the latent demand scale ([#640](https://github.com/epiforecasts/BVDOutbreakSize/issues/640)).
A published count of patients turned away or waiting for a bed would remove it.

### Evaluation

#### The suspected streams and exports cannot be checked

The suspected streams are no longer published, so their forecasts cannot be set against a later observation.
Exports are not forecast at all.
Resumed publication of the suspected counts would remove the first.
The second needs a travel rate that holds over the forecast week.

#### The onset forecast mostly measures scan error

The interval on the onset increment is dominated by the per-scan error rather than epidemic uncertainty, as the [symptom-onset nowcast and forecast](@ref "Symptom-onset nowcast and forecast") explains.
It is a check on the fitted delay and ascertainment more than a case-count prediction.
A published onset table would remove it.

#### Frozen-fit baselines can see later corrections

The frozen re-fits forecast from fixed historical cut-offs, but their persistence baseline reads a data snapshot that can post-date the forecast by weeks.
A correction landing in between is then already in the baseline, as set out under [forecast scoring against a persistence baseline](@ref "Forecast scoring against a persistence baseline").
A snapshot archived at each frozen cut-off would remove it.

#### Reconstructed forecasts are not exact

Past releases that did not store a forecast are reconstructed from their own code and data, but with dependencies resolved at current versions.
The reconstruction is therefore close to, but not exactly, what each release would have produced ([#622](https://github.com/epiforecasts/BVDOutbreakSize/issues/622)).
Pinned manifests for every release would remove it for future releases.
The earliest releases archived no dated vintage record, so they carry no relative skill against the baseline.

## Provinces

### Data

#### Province data are confirmed cases and deaths

The situation reports' spatial tables give per-province confirmed cases and deaths, and those are the only province data in the likelihood.
The per-province analysed specimens enter only the prior on case ascertainment, and occupancy and beds by province are not used.
The onset curve is national, so there is no province nowcast.
Scoring the province analysed volume, occupancy and beds as splits of their printed sums ([#784](https://github.com/epiforecasts/BVDOutbreakSize/pull/784)) would narrow it.

#### The province tables cover part of the window

The per-province vintages start on 15 June and stop before the cut-off.
Before and after them the provinces are informed only through the national streams, as the [province compositions](@ref province-compositions) show.
Province vintages read over the whole window would remove it.

#### Harmonisation backfill has no province breakdown

When a report integrates harmonised records retrospectively, the backfill is published for the country and not by province.
Province forecast windows that hold such a day are left unscored.
A province breakdown of each backfill would remove it.

### Model

#### Four patches, not every province

Ituri, Nord-Kivu and Haut-Uele are modelled on their own, and the other affected provinces are pooled into a fourth patch that holds almost no confirmed cases.
Transmission within a patch is well mixed, so spread inside a province is not represented.
With so few patches, the cross-province correlation of the reproduction number is not identified and tracks its prior, as the [reproduction number by province](@ref "Reproduction number by province") shows.
A health-zone model ([#779](https://github.com/epiforecasts/BVDOutbreakSize/pull/779), [#705](https://github.com/epiforecasts/BVDOutbreakSize/issues/705)) would narrow it.

#### Importation structure is assumed

There is no mobility or origin-destination data for this outbreak, so the gravity kernel under [mixing and importation](@ref "Mixing and importation") is a structural assumption.
Only its intensity is estimated.
Mobility data between provinces would remove it.

#### Province testing informs the prior only

The per-province analysed specimens set the prior on each province's case ascertainment rather than entering the likelihood.
Their positives are the per-province confirmed counts differenced, which the case composition already scores, so fitting them would count the same data twice.
A province's confirmed share is the product of its incidence and its case-finding, so the deaths do most of the work in separating the two, and the position along that ridge is set by the prior ([province parameters against their priors](@ref "Province parameters against their priors")).
Scoring the province analysed volume as a split of the national total ([#784](https://github.com/epiforecasts/BVDOutbreakSize/pull/784)) would narrow it.

#### Beds and occupancy are national only

The treatment-centre model carries one national bed capacity and one national demand, so it cannot represent saturation in one province while another has free beds.
Ituri holds most of the occupied beds, so the national shortfall understates local unmet need.
Occupancy and beds by province ([#784](https://github.com/epiforecasts/BVDOutbreakSize/pull/784)) would remove it.

### Evaluation

#### No province projection has been scored yet

Only [province forecast](@ref "Province forecast") projections are scored.
Earlier releases archived province forecasts that split the national forecast by a fixed share, and those are not scored.
The scores fill in as releases carrying a projection become old enough for their targets to be observed.
Time, and more releases, will remove it.

#### Province forecasts do not sum to the national forecast

Each province is projected on its own, so the provinces need not add up to the national forecast.
The fresh deviations of the reproduction number over the horizon also leave out the fitted cross-province correlation.
Each province's confirmed counts grow with its projected infections from the cut-off, with no delay between infection and report.
A joint projection that reconciles the provinces with the national forecast, and province latent onsets and deaths ([#831](https://github.com/epiforecasts/BVDOutbreakSize/issues/831)), would remove it.

#### The pooled patch is forecast as one

The provinces in the pooled patch are forecast together, so none of them has its own forecast or score.
Splitting the pooled patch, which needs more data from those provinces, would remove it.

## Across the report

### The code was drafted by a language model

The model code, priors and analysis were drafted by a language model from the [mccabe2026](@citet) report and the companion delay reanalysis, then reviewed and revised.
The situation-report figures were also read from the PDFs by a language model, with a second pass to re-read them.
The work has not been independently replicated against the authors' code.
An independent replication, and an independent read of the data, would remove it.
