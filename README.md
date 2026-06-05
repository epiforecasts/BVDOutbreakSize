# Estimating the current size of the 2026 DRC Bundibugyo virus outbreak: a joint Bayesian re-analysis of the McCabe et al. report

**Authors:** Sam Abbott, Kath Sherratt, Samuel Brand and Sebastian Funk.

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://epiforecasts.io/BVDOutbreakSize/stable) [![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://epiforecasts.io/BVDOutbreakSize/dev) [![Tests](https://github.com/epiforecasts/BVDOutbreakSize/actions/workflows/test.yml/badge.svg)](https://github.com/epiforecasts/BVDOutbreakSize/actions/workflows/test.yml) [![codecov](https://codecov.io/gh/epiforecasts/BVDOutbreakSize/branch/main/graph/badge.svg)](https://codecov.io/gh/epiforecasts/BVDOutbreakSize) [![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl) [![Code Style: SciML](https://img.shields.io/static/v1?label=code%20style&message=SciML&color=9558b2&labelColor=389826)](https://github.com/SciML/SciMLStyle) [![DOI](https://zenodo.org/badge/1243778099.svg)](https://doi.org/10.5281/zenodo.20312758)

**Last updated:** 5 June 2026. This is a live report, re-run as new
data arrive, so the estimates change between updates.

**Data as of:** 28 May 2026. DRC counts come from the situation reports
of the Institut National de Santé Publique (INSP); Uganda imports come
from WHO. Estimates are reported as of this date; it can lag the update
date above.

**Abstract.** An outbreak of Ebola disease caused by Bundibugyo virus
(BVD) is ongoing in the Democratic Republic of the Congo (DRC),
with cases also detected across the border in Uganda. Estimating
the likely current size of the outbreak is useful for the response,
but most
cases are not yet reported and have to be inferred from the data
streams that are available. The Imperial College London report
(McCabe et al., [18 May 2026](https://www.imperial.ac.uk/mrc-global-infectious-disease-analysis/research-themes/preparedness-and-response-to-emerging-threats/report-ebola-18-05-2026/),
revised in a [20 May 2026 update](https://www.imperial.ac.uk/media/imperial-college/medicine/mrc-gida/Report-ebola-update-20-05-2026.pdf))
estimates the size with two analyses, geographic spread from the cases
exported to Uganda and back-calculation from suspected deaths in DRC.
Building on that work, we re-analyse the same problem as a single joint
Bayesian model over the latent infection count, fitting all streams
together with priors on the nuisance parameters that the report varies
in scenario sweeps. The model is a discrete-time renewal process on a
daily grid: infections drive a time-varying reproduction number,
symptom onsets follow after a sampled incubation period, and every
observed stream sits downstream of those onsets through its own delay.
Beyond the exported cases and DRC deaths the report uses, we condition
on two further streams, the reported cases in DRC (with an ascertainment
component) and the deaths among exported cases in Uganda. We fit the DRC
streams to the full situation-report trajectory, conditioning on the
between-vintage increments across successive sitreps rather than only
the latest cumulative total, which sharpens the reproduction number. We
also add a no-onward-transmission projected-deaths counterfactual, a
one-week-ahead forecast and an onset-to-death delay sensitivity
analysis. We report the joint posterior from current data, and to
separate newer data from the change in method we also fit the data as of
each report version in sequence (18 May, then the 20 May update),
comparing against both a joint reimplementation of the report's approach
and its published estimates. The cumulative infection count $C_T$ is the
headline latent quantity; we compare against the report's scenarios on
the expected confirmed cases, the ascertained quantity recovered
downstream.

**Scope.** This work is motivated by adding an external view of the
current situation, based on our understanding of real-time infectious
disease dynamics and the infection process that gives rise to observed
epidemic surveillance counts.
We are actively developing it and encourage feedback, so please get in
touch.
We fully support reuse and adaptation.
Find out more in the
[contributing guide](https://epiforecasts.io/BVDOutbreakSize/stable/contributing).

**Use of AI.** The model code and analysis were drafted by a language
model and reviewed and revised under human oversight; the named authors
are responsible for that oversight.

**Why our numbers differ from McCabe et al.** Two reasons. First, the
method: we fit all streams jointly in a single Bayesian model rather
than combining separate scenario analyses (see
[each way our method departs from the report](https://epiforecasts.io/BVDOutbreakSize/stable/analysis#What-we-do-differently-from-McCabe-et-al.)).
Second, the data: our cut-off is later than either McCabe et al.
report version (16 May for the 18 May report; 18 May for the 20 May
update), and we fit the DRC streams to the full run of situation
reports rather than a single total at the cut-off, so the reproduction
number is informed by the shape of the reported trajectory rather than
the prior alone. The joint posterior assumes a single common cut-off for
every data stream, so the deaths, exports and reported-case counts must
all be kept in sync to the same date.
<!-- SHARED:END -->

**See:**
[current outbreak size](https://epiforecasts.io/BVDOutbreakSize/stable/analysis#Summary) ·
[comparison with McCabe et al.](https://epiforecasts.io/BVDOutbreakSize/stable/analysis#Comparison-with-McCabe-et-al.) ·
[how the data streams compare](https://epiforecasts.io/BVDOutbreakSize/stable/analysis#How-the-data-streams-compare) ·
[limitations](https://epiforecasts.io/BVDOutbreakSize/stable/analysis#Limitations) ·
[full joint results](https://epiforecasts.io/BVDOutbreakSize/stable/analysis#Results).

## Installing the package

To use the model and the bundled outbreak data from your own Julia
environment, add the package from the repository:

```julia
using Pkg
Pkg.add(url = "https://github.com/epiforecasts/BVDOutbreakSize")
```

You can then load the model machinery and the data the report is
fitted to:

```julia
using BVDOutbreakSize
obs = load_observations()
```

This gives you the exported model components, constants and data
loaders, enough to build your own analysis on top of the package.
Reproducing the full report (fitting the models and writing the
result tables and plots) can be done in a few ways, described next.

## Running

There are a couple of ways to re-fit the model.

### Re-fit from a clone

```bash
git clone --recurse-submodules https://github.com/epiforecasts/BVDOutbreakSize
cd BVDOutbreakSize
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. scripts/run.jl
```

`scripts/run.jl` fits the models and writes the output CSVs (the
analysis literate is also run as part of the docs build). Running
`docs/examples/analysis.jl` directly instead steps through the full
narrative.

### Re-fit without cloning

`scripts/reproduce.jl` fetches the package, instantiates its
environment, and runs the fit:

```bash
curl -fsSL https://raw.githubusercontent.com/epiforecasts/BVDOutbreakSize/main/scripts/reproduce.jl | julia
```

Outputs land in `./bvd-output`; set `BVD_OUTPUT_DIR` to write them
elsewhere, or `BVD_REF` to a release tag to reproduce a specific
version. The script clones into a temporary directory and runs from
there, so it leaves your own Julia environments untouched.

### Render the docs page

Executes the literate and produces HTML at `docs/build/`:

```bash
julia --project=docs -e 'using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.instantiate()'
julia --project=docs docs/make.jl
```

Updating the observation counts for a new sitrep is a single-file
edit of `data/observations.toml`; the literate picks the new numbers
up automatically.

## Outputs and releases

Each push to `main` regenerates the model outputs as part of the
documentation build and publishes them as a GitHub Release. The
[latest release](https://github.com/epiforecasts/BVDOutbreakSize/releases/latest)
bundles the saved result tables and plots, thinned posterior draws, a
copy of the input `observations.toml` that produced them, a
`site.zip` snapshot of the rendered report site, and `analysis.html`,
a self-contained single-file copy of the report that opens offline
([download the latest](https://github.com/epiforecasts/BVDOutbreakSize/releases/latest/download/analysis.html));
the same artifacts are
written to the repository's `output/` directory on each build. Browse
[all releases](https://github.com/epiforecasts/BVDOutbreakSize/releases)
for earlier output bundles — major versions of the report are kept as
GitHub Releases.

The rendered report is published from the
[`gh-pages` branch](https://github.com/epiforecasts/BVDOutbreakSize/tree/gh-pages),
where past and development versions of the analysis page can be found.

## Submodules

- `external/bdbv-linelist-analysis` — Bayesian reanalysis of the 2012
  Isiro BDBV line list (Rosello et al. 2015). Source of the
  informative onset-to-death gamma shape and scale priors.

## Citation

If you use or build on this project, please cite the works this
repository depends on:

- **This project** — Abbott, S., Brand, S., Funk, S. (2026).
  *BVDOutbreakSize: joint forward-generative Turing model for the
  2026 DRC Bundibugyo outbreak.*
  <https://github.com/epiforecasts/BVDOutbreakSize>.
  DOI: [10.5281/zenodo.20312758](https://doi.org/10.5281/zenodo.20312758).
- **INSP situation reports** that supply the DRC suspected-case and
  suspected-death counts and the sitrep cumulative trajectory —
  Institut National de Santé Publique, Démocratique du Congo (2026).
  *Situation reports on the 17th Ebola epidemic.*
  <https://insp.cd/ebola-17eme-epidemie/>. Transcribed across health
  zones by INRB-UMIE/Ebola_DRC_2026,
  <https://github.com/INRB-UMIE/Ebola_DRC_2026>.
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
- **Imperial reports** that this work re-implements and compares
  against, in both released versions —
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
- **Onset-to-death delay reanalysis** that this work uses for
  delay priors — Funk, S. (2026). *bdbv-linelist-analysis:
  Bayesian reanalysis of the 2012 Isiro Bundibugyo line list.*
  <https://github.com/sbfnk/bdbv-linelist-analysis>.

## Further references

- Rosello et al., *Ebola virus disease in DRC, 1976–2014*, eLife 2015.
  Original Isiro 2012 onset-to-death gamma point estimate.
- Imai et al., *Estimating the potential total number of novel
  coronavirus cases in Wuhan City*, Imperial COVID-19 Response Team
  Report 1, 17 January 2020. Methodological template for Method 1.
- Charniga et al., *Best practices for estimating and reporting
  epidemiological delay distributions*, PLOS Computational Biology 2024.
  Followed for the delay-distribution reporting here.
