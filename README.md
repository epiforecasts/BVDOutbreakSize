# Estimating the current size of the 2026 DRC Bundibugyo virus outbreak

**Authors:** Sam Abbott, Kath Sherratt, Samuel Brand and Sebastian Funk.

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://epiforecasts.io/BVDOutbreakSize/stable) [![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://epiforecasts.io/BVDOutbreakSize/dev) [![Tests](https://github.com/epiforecasts/BVDOutbreakSize/actions/workflows/test.yml/badge.svg)](https://github.com/epiforecasts/BVDOutbreakSize/actions/workflows/test.yml) [![codecov](https://codecov.io/gh/epiforecasts/BVDOutbreakSize/branch/main/graph/badge.svg)](https://codecov.io/gh/epiforecasts/BVDOutbreakSize) [![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl) [![Code Style: SciML](https://img.shields.io/static/v1?label=code%20style&message=SciML&color=9558b2&labelColor=389826)](https://github.com/SciML/SciMLStyle) [![DOI](https://zenodo.org/badge/1243778099.svg)](https://doi.org/10.5281/zenodo.20312758)

**Last updated:** on each rebuild.
This is a live report, re-run as new data arrive, so the estimates change between updates.

**Data as of:** the most recent situation report.
DRC counts come from the situation reports of the Institut National de Santé Publique (INSP).
Uganda imports come from WHO.
The rendered report fills in the build date and the exact data cut-off automatically.

> **Reduced data streams.**
> The situation reports stopped publishing the daily new-suspected-case count on 6 August 2026, so it is frozen at its last published value.
> The treatment-centre patient-movement table was dropped on 3 August, and although the reports returned to their full format on 12 August with per-province prose covering much of the same ground, that stream stays frozen pending a decision on whether the two are comparable.
> Every other stream continues, so the most recent weeks rest on fewer streams than earlier ones.
> See the inclusion rules in `data/README.md`.

**See:**
[current outbreak size](https://epiforecasts.io/BVDOutbreakSize/stable/analysis#Summary) ·
[one-week-ahead forecast](https://epiforecasts.io/BVDOutbreakSize/stable/analysis#One-week-ahead-forecast-results) ·
[estimate evolution across releases](https://epiforecasts.io/BVDOutbreakSize/stable/sensitivity#Estimate-evolution-across-releases) ·
[comparison with McCabe et al.](https://epiforecasts.io/BVDOutbreakSize/stable/sensitivity#Comparison-with-McCabe-et-al.) ·
[comparison with Chamla et al.](https://epiforecasts.io/BVDOutbreakSize/stable/sensitivity#Comparison-with-Chamla-et-al.) ·
[how the data streams compare](https://epiforecasts.io/BVDOutbreakSize/stable/sensitivity#Outbreak-size-estimated-by-each-data-stream) ·
[limitations](https://epiforecasts.io/BVDOutbreakSize/stable/analysis#Limitations) ·
[full joint results](https://epiforecasts.io/BVDOutbreakSize/stable/analysis#Results).

**Abstract.** An outbreak of Ebola disease caused by Bundibugyo virus (BVD) is ongoing in the Democratic Republic of the Congo (DRC), with cases also detected across the border in Uganda.
This is a real-time joint Bayesian estimate of the current size of that outbreak, refreshed as new data arrive.
Most infections are not yet reported, so the current size has to be inferred from the surveillance data that are available.
The model is a discrete-time renewal process on a daily grid.
It fits in a single posterior the daily counts the INSP situation reports publish, from the suspected and confirmed case and death series through to the laboratory, isolation, treatment-centre and recovery records.
The digitised epidemic curve by symptom-onset date is fitted with them, the only direct observation of onsets.
So are the cases and deaths exported to Uganda, taken from the WHO situation reports and Disease Outbreak News.
The same infection process generates all of them, staged to daily symptom onsets and routed into every observation stream.
A genetic bound on the time to the most recent common ancestor and priors from the McCabe et al. report complete the inputs.
From these it estimates the infections and deaths to date, reported and unreported, the time-varying reproduction number with its growth rate and doubling time, the case-fatality ratio, and the ascertainment of each surveillance system.
Every release projects each DRC stream a week ahead.
The forecasts are scored with the continuous ranked probability score against the data that arrive later and against a persistence baseline.
The aim is a transparent estimate of how far the outbreak has already grown, and of how much each published stream contributes to it.

**Scope.** This work adds an external view of the current situation, drawing on our understanding of real-time infectious disease dynamics and the infection process behind the observed surveillance counts.
We are developing it and encourage feedback, so please get in touch.
We support reuse and adaptation.
Find out more in the [contributing guide](https://epiforecasts.io/BVDOutbreakSize/stable/contributing).

**Use of AI.** The model code and analysis were drafted by a language model, then reviewed and revised under human oversight.
The named authors are responsible for that oversight.

<!-- SHARED:END -->

## Installing the package

To use the model and the bundled outbreak data from your own Julia environment, add the package from the repository:

```julia
using Pkg
Pkg.add(url = "https://github.com/epiforecasts/BVDOutbreakSize")
```

You can then load the model machinery and the data the report is fitted to:

```julia
using BVDOutbreakSize
obs = load_observations()
```

This gives you the exported model components, constants and data loaders, enough to build your own analysis on top of the package.
Reproducing the full report (fitting the models and writing the result tables and plots) can be done in a few ways, described next.

## Running

There are a couple of ways to re-fit the model.

### Re-fit from a clone

```bash
git clone --recurse-submodules https://github.com/epiforecasts/BVDOutbreakSize
cd BVDOutbreakSize
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. scripts/run.jl
```

`scripts/run.jl` fits the models and writes the output CSVs (the analysis literate is also run as part of the docs build).
Running `docs/examples/analysis.jl` directly instead steps through the full narrative.

### Re-fit without cloning

`scripts/reproduce.jl` fetches the package, instantiates its environment, and runs the fit:

```bash
curl -fsSL https://raw.githubusercontent.com/epiforecasts/BVDOutbreakSize/main/scripts/reproduce.jl | julia
```

Outputs land in `./bvd-output`.
Set `BVD_OUTPUT_DIR` to write them elsewhere, or `BVD_REF` to a release tag to reproduce a specific version.
The script clones into a temporary directory and runs from there, so it leaves your own Julia environments untouched.

### Render the docs page

Executes the literate and produces HTML at `docs/build/`:

```bash
julia --project=docs -e 'using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.instantiate()'
julia --project=docs docs/make.jl
```

### Updating the data

The observations live in `data/observations.toml`.
The literate picks up new numbers automatically.
The figures come from the INSP situation reports, transcribed by [INRB-UMIE/BDBV2026-Data](https://github.com/INRB-UMIE/BDBV2026-Data), and are refreshed for each new sitrep with the `scripts/` tooling (`task download-sitreps`, `task confirm-data`):

- The cumulative confirmed-case and confirmed-death series come from the upstream `national_*` daily CSVs, regenerated and spot-confirmed against our own scan by `scripts/confirm_insp_data.jl`.
- Every other stream (suspected totals, laboratory cumulatives, the 24h analysed volume, daily new suspects, isolation occupancy, bed capacity and recoveries) is read directly from the situation-report PDFs (fetched by `scripts/download_sitreps.jl` into `data/sitrep_pdfs/`) and recorded in `data/insp_sitrep_scanned.csv`.

`scripts/confirm_insp_data.jl` cross-checks the scan against the upstream national series and exits non-zero on any disagreement.

## Outputs and releases

Each push to `main` regenerates the model outputs as part of the documentation build and publishes them as a GitHub Release.
The [latest release](https://github.com/epiforecasts/BVDOutbreakSize/releases/latest) bundles: the saved result tables and plots; thinned posterior draws; a copy of the input `observations.toml` that produced them; a `site.zip` snapshot of the rendered report site; and `analysis.html`, a self-contained single-file copy of the report that opens offline ([download the latest](https://github.com/epiforecasts/BVDOutbreakSize/releases/latest/download/analysis.html)).
The same artifacts are written to the repository's `output/` directory on each build.
Browse [all releases](https://github.com/epiforecasts/BVDOutbreakSize/releases) for earlier output bundles.
Major versions of the report are kept as GitHub Releases.

The rendered report is published from the
[`gh-pages` branch](https://github.com/epiforecasts/BVDOutbreakSize/tree/gh-pages), where past and development versions of the analysis page can be found.

## Submodules

- `external/bdbv-linelist-analysis` — Bayesian reanalysis of the 2012 Isiro BDBV line list (Rosello et al. 2015).
  Source of the informative onset-to-death gamma shape and scale priors.

## Citation

If you use or build on this project, please cite the works this repository depends on:

- **This project** — Abbott, S., Brand, S., Funk, S. (2026).
  *BVDOutbreakSize: joint forward-generative Turing model for the
  2026 DRC Bundibugyo outbreak.*
  <https://github.com/epiforecasts/BVDOutbreakSize>.
  DOI: [10.5281/zenodo.20312758](https://doi.org/10.5281/zenodo.20312758).
- **INSP situation reports** that supply the DRC suspected-case and
  suspected-death counts and the sitrep cumulative trajectory —
  Institut National de Santé Publique, Démocratique du Congo (2026).
  *Situation reports on the 17th Ebola epidemic.*
  <https://insp.cd/ebola-17eme-epidemie/>. Transcribed by
  INRB-UMIE/BDBV2026-Data,
  <https://github.com/INRB-UMIE/BDBV2026-Data>.
- **WHO situation reports and Disease Outbreak News** that supply
  the Uganda import-case counts and the dated detection and death
  events for the first Uganda import —
  World Health Organization Regional Office for Africa (2026).
  *Ebola disease caused by Bundibugyo virus outbreak, Democratic
  Republic of the Congo and Uganda — Weekly External Situation
  Report 01.* Data as of 18 May 2026.
  World Health Organization (2026). *Disease Outbreak News:
  Ebola disease caused by Bundibugyo virus — Democratic Republic
  of the Congo and Uganda (DON602, DON603).* Source of the first
  Uganda import hospital-admission date (11 May 2026) and the fatal
  import death date (14 May 2026).
- **McCabe et al. estimates** that this work re-implements and compares
  against, in all three released versions (the two Imperial reports and
  the peer-reviewed Lancet publication) —
  McCabe, R., Ebbarnezh, L., Okware, S., Fotsing, R., Koua, E.,
  Mbaka, P., Lofungola, A., van Elsland, S. L., McMenamin, M.,
  Ferguson, N., le Polain de Waroux, O., Cori, A. (2026).
  *Estimation of the size of the outbreak of Ebola disease caused
  by Bundibugyo virus in DRC.* Imperial College London, 18 May 2026.
  DOI: [10.25560/130007](https://doi.org/10.25560/130007).
  [Report page](https://www.imperial.ac.uk/mrc-global-infectious-disease-analysis/research-themes/preparedness-and-response-to-emerging-threats/report-ebola-18-05-2026/).
  McCabe, R. and others (2026).
  *Estimation of the size of the Ebola outbreak caused by Bundibugyo
  virus in DRC: May 20, 2026 update.* Imperial College London,
  20 May 2026.
  DOI: [10.25560/13005307](https://doi.org/10.25560/13005307).
  [Report PDF](https://www.imperial.ac.uk/media/imperial-college/medicine/mrc-gida/Report-ebola-update-20-05-2026.pdf).
  McCabe, R., Ebbarnezh, L., Okware, S., Fotsing, R., Koua, E.,
  Mbaka, P., Lofungola, A., Ebengo, D. M., Mbala, P. K., Bishola, T. T.,
  Ibolobolo, C. M., Matondo, H. M., Sibo, J.-C. M., van Elsland, S. L.,
  McMenamin, M., Ferguson, N., le Polain de Waroux, O., Cori, A. (2026).
  *Estimation of the Ebola outbreak size in the Democratic Republic
  of the Congo.* The Lancet Infectious Diseases, correspondence,
  online first 9 June 2026; estimates as of 27 May 2026.
  DOI: [10.1016/S1473-3099(26)00299-9](https://doi.org/10.1016/S1473-3099(26)00299-9).
- **Chamla et al. projection** that this work compares against — a
  stochastic SEIRD ensemble that projects cumulative laboratory-confirmed
  cases forward (a confirmed-case count, not an ascertainment-corrected
  total) — Chamla, D., Belizaire, M. R. D., Co, I. F., Jinadu, A.,
  Mamadu, I., Atagbaza, A. O. (2026). *Size of the 2026 Ebola outbreak and
  risk of cross-border spillover from Bundibugyo virus in Ituri Province,
  DR Congo, and its implications for preparedness: a recalibrated stochastic
  modelling study.* The Lancet Infectious Diseases, 25 June 2026.
  DOI: [10.1016/S1473-3099(26)00320-8](https://doi.org/10.1016/S1473-3099(26)00320-8).
- **Onset-to-death delay reanalysis** that this work uses for
  delay priors — Funk, S. (2026). *bdbv-linelist-analysis:
  Bayesian reanalysis of the 2012 Isiro Bundibugyo line list.*
  <https://github.com/sbfnk/bdbv-linelist-analysis>.

## Funding

This work was funded by the National Institute for Health and Care Research (NIHR) Health Protection Research Unit in Health Analytics & Modelling, a partnership between the UK Health Security Agency, Imperial College London and the London School of Hygiene & Tropical Medicine (grant code NIHR207404).
The views expressed are those of the author(s) and not necessarily those of the NIHR, UK Health Security Agency or the Department of Health and Social Care.

## Further references

- Rosello et al., *Ebola virus disease in DRC, 1976–2014*, eLife 2015.
  Original Isiro 2012 onset-to-death gamma point estimate.
- Imai et al., *Estimating the potential total number of novel
  coronavirus cases in Wuhan City*, Imperial COVID-19 Response Team
  Report 1, 17 January 2020. Methodological template for Method 1.
- Charniga et al., *Best practices for estimating and reporting
  epidemiological delay distributions*, PLOS Computational Biology 2024.
  Followed for the delay-distribution reporting here.
