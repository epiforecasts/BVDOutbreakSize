# Observation data

This directory is the single source of truth for the observations the
analysis conditions on.
`load_observations()` (in `src/BVDOutbreakSize.jl`) reads it; nothing in the
model hardcodes counts.

## Files

| File | What it is |
|---|---|
| `observations.toml` | The manifest the model loads. Every stream is a value or a dated `dates`/`values` history plus a prose `source =` citation. Edit this to advance the analysis. |
| `insp_sitrep_scanned.csv` | Our own direct scan of the INSP SitRep PDFs, one row per report (`date de rapportage`), with a free-text `notes` column recording the headline tiles, laboratory section and table figures. The audit trail behind the PDF-sourced streams in `observations.toml`. |
| `onset_curve_scanned.csv` | Confirmed cases by symptom-onset date, digitised from the analytique-format SitReps' onset epidemic-curve figure (one block per vintage). Not fitted; see the section below. |
| `released_estimates.csv` | Published point estimates for comparison. |
| `report-snapshot*.toml` | Frozen Imperial report point estimates at fixed vintages. |

## Where the data comes from

DRC national figures are the INSP (Institut National de Santé Publique)
situation reports for the 17th Ebola epidemic.
Two sources feed the manifest:

- The cumulative confirmed-case, confirmed-death, recovered and isolation
  series are cross-checked against the INRB-UMIE mirror's clean `national_*`
  daily CSVs
  (<https://github.com/INRB-UMIE/BDBV2026-Data>, `data/insp_sitrep/processed`);
  `scripts/confirm_insp_data.jl` regenerates the confirmed streams from them.
- Every other stream (suspected daily inflow, the 24h analysed laboratory
  volume, isolation occupancy, bed capacity, recoveries and the CTE/CT/CI
  patient-movement rows) is read directly from the SitRep PDFs and recorded
  in `insp_sitrep_scanned.csv`.

INSP is the primary source: it publishes first, and since SitRep 059
(12 July 2026) its richer "analytique / édition quotidienne" PDFs carry
content the mirror does not transcribe (epidemic curve by symptom-onset
date, age/sex pyramids, five provinces).
The mirror usually lags INSP by a report or two.

## Symptom-onset epidemic curve (`onset_curve_scanned.csv`)

From SitRep 059 the analytique PDFs carry a figure of confirmed cases by
symptom-onset date (`courbe épidémique par date de début des symptômes`,
DHIS2 line list), split by outcome (Vivant / Décédé).
It is the only published source for the onset-date distribution, and it is a
raster bar chart with no accompanying data table.
`scripts/digitize_onset_curve.jl` (the dependency-free Julia reference)
recovers the daily counts from the figure pixels and writes one block per
vintage to `onset_curve_scanned.csv` (columns `sitrep`, `report_date`,
`onset_date`, `confirmed_alive`, `confirmed_dead`, `confirmed_total`).
A Python port (`scripts/digitize_onset_curve.py`) for the automated
data-updater produces a byte-identical file; see `scripts/README.md`.
Both self-calibrate each figure from its axis ticks; the only manual input is
each vintage's rightmost x-axis tick date (in the script `CONFIG`).

These counts are approximate.
The digitised per-vintage totals run about 2% below the printed figure `n`
(SitRep 064: 2018 vs printed n = 2 064), and individual daily bars carry
roughly ±1–2 cases of pixel noise.
SitReps 059 and 060 reuse one figure, as do 061 and 062, so the five scanned
vintages hold three distinct onset snapshots (report dates 12, 14 and
17 July).

This stream is **not fitted**: the model does not yet read
`onset_curve_scanned.csv`.
It is captured so an onset-based likelihood and a reporting-delay component
can be added later.

### Rough onset-to-report delay

Because each vintage redraws the same onset cohort at a later report date, the
scanned curves form a reporting triangle: old onset dates are stable across
vintages while recent ones fill in as more confirmations arrive (e.g. onset
10 July reads 4 → 9 → 27 across the 12, 14 and 17 July snapshots).
Taking the latest snapshot as the near-complete reference for onset dates at
least ~12 days old, the empirical proportion of eventually-reported confirmed
cases reported within `d` days of onset is roughly 65% by 7 days, 85% by
~11 days and 95%+ by ~2 weeks, near-complete by ~3 weeks (median ~5–6 days).
This is a coarse estimate: it rests on three digitised snapshots only, the
reference snapshot is itself right-truncated for its most recent onsets, and
the delay it measures is onset → confirmed-and-reported (it folds together
care-seeking, lab confirmation and reporting).

### Fetching a SitRep from INSP directly

INSP blocks some default user agents with HTTP 403, but a browser
User-Agent returns HTTP 200.
Each SitRep post embeds the real PDF URL in a `pdfemb-data` base64 blob.
To get the direct PDF for a post (find the post at
<https://insp.cd/ebola-17eme-epidemie/> or via `check_new_sitreps.jl`):

```sh
UA='Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/605.1.15'
# 1. list posts and pick the id/slug you want
curl -s -A "$UA" \
  'https://insp.cd/wp-json/wp/v2/posts?search=sitrep&per_page=100&_fields=id,slug,date'
# 2. decode the embedded PDF URL for that post id, then download it
curl -s -A "$UA" 'https://insp.cd/wp-json/wp/v2/posts/<POST_ID>' \
  | python3 -c 'import sys,json,re,base64; \
c=json.load(sys.stdin)["content"]["rendered"]; \
b=re.search(r"pdfemb-data=([A-Za-z0-9_-]+)",c).group(1); b+="="*(-len(b)%4); \
print(json.loads(base64.urlsafe_b64decode(b))["url"])'
```

The INRB-UMIE mirror's PDFs can also be pulled with
`scripts/download_sitreps.jl` (lags INSP, no `analytique` reports).

## Checking for new SitReps (do this before every refresh)

```sh
julia --project=scripts scripts/check_new_sitreps.jl
```

It lists the INSP-published SitReps, compares them with the latest report in
`insp_sitrep_scanned.csv`, and prints the gap (exit code 1 if behind, 0 if
current).
This is the guard against the manifest silently drifting weeks behind again.

## What to read from each SitRep

Read these from the PDF (the `date de rapportage`, not the filename or the
publication date, is the key each row is stored under) and record them in the
matching `observations.toml` history:

| From the SitRep | `observations.toml` stream |
|---|---|
| Headline `Cumul cas confirmés` | `confirmed_case_history` |
| Headline `Cumul décès confirmés` (+ lethality %) | `confirmed_death_history` |
| Headline `Guéris — cumul` | `recovered_history` |
| Headline `Patients en isolement` (+ occupancy %) | `isolation_history`; occupancy → `bed_capacity_history` (= occupancy ÷ rate) |
| Headline / Tableau 3 `Cas suspects du jour` | `suspected_daily_history` |
| Laboratory 24h `échantillons analysés` per province (§4.3 table, or §3.2 bullets in the analytique format) | `tests_analysed_daily_history` |
| Occupation table (Tableau 6/7, or Tableau 5/6 in the analytique format): `Total admissions (24h)`, `Sorties — décédés / non-cas / évadés`, `Patients au lit (J-1)`, `dont confirmés / suspects` | `treatment_admissions_history`, `treatment_deaths_history`, `treatment_ruleout_history`, `treatment_absconded_history`, `treatment_aulit_history`, `treatment_confirmed_incare_history`, `treatment_suspect_incare_history` |

Always cross-check the confirmed/death/recovered/isolation headline against
the INRB-UMIE `national_*` CSV for the same date; note any disagreement in
the `source =` string and prefer the auditable value.

## Inclusion rules and conventions

- **Cut-off 28 May 2026 (SitRep 014)**: the last vintage with a coherent
  national laboratory total. The suspected cumulative streams freeze at their
  last stable vintage (26 May). Confirmed cases/deaths and the post-cut-off
  daily streams are recorded and, for the streams the model reads
  (e.g. the 24h analysed denominators), fitted beyond the cut-off.
- **24h analysed denominator**: only include a day whose confirmed-case
  increment does not exceed that day's national `analysed` count (you cannot
  confirm more specimens than you analysed). A province that reports
  "reçus / en cours d'analyse" but no completed count contributes 0, like a
  non-reporting province. Sum only completed per-province analysed counts.
- **Unpublished reports**: some SitRep numbers were never issued (029/12 Jun,
  043/26 Jun, 045/28 Jun). The series simply step over the missing date; note
  it in the `source =` string.
- **Headline vs table disagreement**: prefer the auditable table sum and
  record the headline as the discrepancy (precedent: SitRep 009's 220 over
  the erroneous 119; SitRep 061's Fin J 736 over the page-1 typo 737).
- **Format change from SitRep 059**: the analytique format moved the 24h
  analysed counts to §3.2 Laboratoire bullets, dropped the page-1 suspect-
  death-of-day subtitle, and reports the daily suspect total in Tableau 3.
  The DHIS2 suspect-death-of-day measure is not continuous in this format, so
  `suspected_daily_deaths_history` is frozen at 11 July; see its `source =`
  string.

After editing, validate with the loader and its invariants:

```sh
julia --project=. -e 'using BVDOutbreakSize: load_observations; \
  load_observations(); println("ok")'
julia --project=. -e 'using Pkg; Pkg.test()'   # includes test_load_observations
```
