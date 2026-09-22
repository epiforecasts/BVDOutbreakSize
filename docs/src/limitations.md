# Limitations

The limitations are grouped by the data, the model assumptions and design, and the implementation, with the most consequential first in each group.

**Data and what it can support**

- *Most quantities rest on weakly-informed priors.* Nearly all of the delays, the case-fatality ratio and the laboratory assumptions are set by priors informed at best by a handful of literature sources, often from other outbreaks.
  In places the priors instead reflect our own judgement rather than anything from this outbreak.
  The data do little to move them, so these posteriors largely track their priors.
  We fit the between-report increments, so the trajectory informs the change in the reproduction number over the window.
  It is uninformative about the delays, the surveillance dispersion or the reporting fractions on their own.
- *Almost every count is report-dated.* The digitised onset curve is the only series carrying symptom-onset dates, and it covers confirmed cases from SitRep 059 onward.
  Everything else is a total at the report date, so the epidemic's timing is recovered mainly through the assumed delays.
- *Fitted to aggregate counts.* The DRC data are national and per-province situation-report totals, and the Uganda data are three export cases with one death.
  We do not have a line list or information on case definitions or reporting completeness.
  The laboratory testing series gives partial information on testing capacity, but it is incomplete and stops at the cut-off.
  Every estimate is a model-based extrapolation under strong assumptions, not a measurement.
- *Later sitreps revise earlier figures.* A later situation report can revise an earlier total up or down as suspects are reclassified and newly-reporting health zones are added, and ascertainment probably rose over the window.
  We do not model this revision process.
- *Streams share one case pool.* They are fitted as conditionally independent given latent incidence but observe overlapping people.
  This can understate uncertainty.
  Whether they imply mutually consistent outbreak sizes is assessed on the sensitivity page, which fits each count stream on its own and sets the resulting cumulative intervals against the joint.

**Model assumptions and design**

- *Inherits McCabe et al.'s epidemiological assumptions.* A single zoonotic seed, an assumed generation interval, and no depletion of susceptibles.
  The onset-to-death delay is grounded on Isiro 2012 and the genetic seeding bound on an external clock rate.
  Neither propagates cross-outbreak or clock uncertainty.
- *Importation structure is assumed, not measured.* There is no mobility or origin-destination data for this outbreak, so the gravity kernel is a structural assumption.
  Its intensity is weakly identified against the secondary provinces' seeds, since both raise a province's early incidence.
- *Four patches, not the full provincial detail.* Ituri, Nord-Kivu and Haut-Uele are modelled individually and every other affected province is pooled into a fourth patch, which takes the population-weighted mean of its members' capitals.
  Transmission within a patch is well mixed, so spread inside a province is not represented.
- *Provincial testing enters the prior, not the likelihood.* The alternative was a per-patch laboratory process, fitting each province's analysed volume and positives so that the data set each patch's testing capacity directly.
  It was not taken because those positives are the per-province confirmed counts differenced, which the composition already scores, so they would enter the joint density twice.
- *Intervention ramp is weakly identified.* With only a few sitreps straddling it, the ramp effect and the pre-ramp reproduction number are not well separated.
- *Single national bed capacity.* The treatment-centre model carries one national bed capacity and one national demand, so it cannot represent local saturation.
  On 13 June Ituri was at 93.9% occupancy while Sud-Kivu was at 21.9%.
  The national bed shortfall therefore understates local unmet need.

**Implementation**

- *LLM-driven reimplementation.* The model code, priors and analysis were drafted by a language model from the [mccabe2026](@citet) report and the companion delay reanalysis.
  It was then reviewed and revised.
  It has not been independently replicated against the authors' code.
