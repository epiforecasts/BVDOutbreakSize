# Authors, funding and acknowledgements

## Authors

Sam Abbott, Kath Sherratt, Samuel Brand and Sebastian Funk.

The model code and analysis were drafted by a language model, then reviewed and revised under human oversight.
The named authors are responsible for that oversight.
The [limitations](limitations.md) set out what that means for the work.

## Funding

This work was funded by the National Institute for Health and Care Research (NIHR) Health Protection Research Unit in Health Analytics & Modelling, a partnership between the UK Health Security Agency, Imperial College London and the London School of Hygiene & Tropical Medicine (grant code NIHR207404).
The views expressed are those of the authors and not necessarily those of the NIHR, UK Health Security Agency or the Department of Health and Social Care.

## Data sources

The DRC counts are taken from the situation reports of the Institut National de Santé Publique (INSP), which are the source of every national, provincial and health-zone series here.
The Uganda exports and the deaths among them come from WHO situation reports and Disease Outbreak News.
The INRB-UMIE mirror is used to cross-check the scanned INSP totals, not as a source in its own right.
We thank the INSP and WHO teams who compile and publish these reports.

## Prior work this builds on

This work started as a replication of the [mccabe2026](@citet) report on the size of the outbreak.
It has since grown into a joint model fitted to more of the published data streams.
The [aim and origins](aim.md) page sets out what has changed since the replication.
The onset-to-admission and onset-to-death delays come from a Bayesian reanalysis [bdbv_linelist_analysis_2026](@cite) of the line list from the 2012 Bundibugyo outbreak in Isiro [rosello2015](@cite).

## Citing this work

The repository carries a `CITATION.cff` with the current metadata, and each release is archived with a DOI.

```
Abbott, S., Sherratt, K., Brand, S. and Funk, S.
BVDOutbreakSize: estimating the current size of the 2026 DRC Bundibugyo
virus outbreak. https://doi.org/10.5281/zenodo.20312758
```

The DOI above resolves to the latest release.
The [news](news.md) page records what changed in each one.

## Reuse

The code is MIT licensed.
We support reuse and adaptation, and we welcome feedback and corrections through the [repository](https://github.com/epiforecasts/BVDOutbreakSize).
The [contributing guide](contributing.md) covers how to build and extend it.
