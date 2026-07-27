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

These counts are approximate, and the error is a few percent in either
direction per scan, independent between vintages.
Against the printed figure `n` it ranges from −3.0% (SitReps 069/070/071:
2260 against n = 2 329) to +1.6% (SitRep 068: 2344 against n = 2 308), with
SitRep 064 at −2.2% (2018 against n = 2 064) and SitRep 072 at +0.4% (2531
against n = 2 521).
Individual daily bars carry roughly ±1–2 cases of pixel noise.
Some of the shortfall is the faded bars inside the
`données potentiellement incomplètes` band at the right of each figure, whose
lightened fill falls outside the colour masks, but that mechanism can only
lose cases and so does not explain the overshoots; treat the sign as unknown.

One consequence deserves emphasis before anyone fits this stream.
Late reporting only ever adds cases, so an onset date's count must be
non-decreasing across vintages, and the scans do not respect that.
On onset dates more than three weeks before the earliest report date in the
file (12 July, so onsets before 21 June), the scanned totals move both ways
between consecutive distinct snapshots: 064 → 065 falls by a net 36 cases
across 34 of 54 such days, and every other consecutive pair falls somewhere
too.
A between-vintage increment of a few cases is therefore at or below the noise
floor, which bounds what a reporting-delay estimate built from those
increments can support.
See issue #488.

SitReps 059 and 060 reuse one figure, as do 061 and 062, and 069, 070 and
071, so the thirteen scanned vintages hold nine distinct onset snapshots
(report dates 12, 14, 17, 18, 19, 20, 21, 22 and 25 July).
SitRep 068 (21 July) is now included: its PDF was unreachable while the
INSP fetch was broken, and it does carry a figure (n = 2 308).
The embedded figure is re-rendered at whatever size the layout needs and moves
in both directions: it shrinks to 1009×583 at SitRep 069 (against 1257×698 for
SitRep 064) and then grows to 1277×799 at SitRep 072.
Both moves have broken detection once — the shrink by falling under an absolute
blue floor, the growth by pushing anti-aliasing over an absolute orange cut —
which is why the figure-detection and axis-tick thresholds in both scripts are
pixel fractions rather than counts.

This stream is **not fitted**: the model does not yet read
`onset_curve_scanned.csv`.
It is captured so an onset-based likelihood and a reporting-delay component
can be added later.

### Rough onset-to-report delay

Because each vintage redraws the same onset cohort at a later report date, the
scanned curves form a reporting triangle: old onset dates are stable across
vintages while recent ones fill in as more confirmations arrive (e.g. onset
10 July reads 4 → 9 → 27 → 33 → 37 → 39 → 40 across the 12, 14, 17, 18, 19,
20 and 22 July snapshots).
Taking the latest snapshot as the near-complete reference for onset dates at
least ~12 days old, the empirical proportion of eventually-reported confirmed
cases reported within `d` days of onset is roughly 60% by 7 days, 85% by
~10 days, 90%+ by ~12 days and 95% by ~2.5 weeks, near-complete (98–99%) by
~3 weeks (median ~5–6 days).
This supersedes the earlier three-snapshot estimate and shifts it slightly
later; the seven-snapshot curve is flatter through the second week than the
three-snapshot one suggested.
This is still a coarse estimate: it rests on seven digitised snapshots over a
ten-day window, the reference snapshot is itself right-truncated for its most
recent onsets, and the delay it measures is onset → confirmed-and-reported
(it folds together care-seeking, lab confirmation and reporting).

Two caveats have been added since these figures were computed, and neither has
been folded into them.
The 21 and 25 July snapshots now exist (nine distinct snapshots, not seven), so
the estimate is one report behind at both ends.
More importantly, the scan noise documented above is two-sided and of the same
order as the between-vintage increments this estimate is built from, so the
quoted percentages are almost certainly more precise than the data supports.
Treat them as indicative until the triangle is re-estimated with that noise
carried through.

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
matching `observations.toml` history.
**Advance every stream in the table below on each update, whenever the SitRep
prints the value** — not just the confirmed/suspected headlines. The
isolation, bed-capacity, recovered, laboratory and CTE patient-movement
streams move together with the headlines: the treatment/occupancy streams
share one Tableau, and the loader invariants tie the latest treatment-flow day
to the isolation window, so advancing the headlines while leaving these behind
both drops data and can break the loader test.

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
  non-reporting province. Sum only completed per-province analysed counts. A
  province worded "documentés dans le réseau de collecte" rather than
  "analysés" is a shade looser but usable (precedent: SitRep 065, 18 July,
  Ituri 167 + Nord-Kivu 127 = 294); note the wording in the `source =` string.
- **Occupation table (Tableau 6/7)**: use it only when it balance-closes
  (`Patients au lit (J-1) + Total admissions (24h) − Total sorties (24h) =
  Patients en isolement (Fin J)`) and the `dont confirmés + dont suspects`
  split sums to the Fin J total. Bed capacity = Fin J ÷ `taux d'occupation`,
  which should equal the printed `Nombre de lits`. A **carried-forward tile**
  (page-1 isolement count *and* occupancy rate byte-identical to the prior day
  with no fresh Tableau to corroborate) is a reporting gap: record the raw
  value in the CSV note but exclude it from the fitted series (precedent:
  30 June, 3–4 July).
- **Unpublished reports**: some SitRep numbers were never issued (029/12 Jun,
  043/26 Jun, 045/28 Jun). The series simply step over the missing date; note
  it in the `source =` string.
- **Headline vs table disagreement**: prefer the auditable table sum and
  record the headline as the discrepancy (precedent: SitRep 009's 220 over
  the erroneous 119; SitRep 061's Fin J 736 over the page-1 typo 737).
- **Retrospective harmonisation steps**: when a vintage's asterisked headline
  jumps far beyond its own 24h new-case count because a provincial base was
  integrated, record the harmonised headline (it is what the INRB-UMIE
  `national_*` CSV carries) and state the gross-vs-net split in the
  `source =` string, so the artefact is visible to whoever reads the
  increment the daily likelihood actually fits. Precedents: SitRep 065
  (18 July, Haut-Uélé/Nia-Nia reattachment, +83 gross vs +77 net) and
  SitRep 069 (22 July, Ituri/Nord-Kivu base integration, +97 gross vs +369
  net cases and +62 gross vs +236 net deaths — much the largest so far, and
  called out in the report's own footnote).
- **Format change from SitRep 059**: the analytique format moved the 24h
  analysed counts to §3.2 Laboratoire bullets, dropped the page-1 suspect-
  death-of-day subtitle, and reports the daily suspect total in Tableau 3.
  The DHIS2 suspect-death-of-day measure is not continuous in this format, so
  `suspected_daily_deaths_history` is frozen at 11 July; see its `source =`
  string.

## Other signals to scan for and flag (not yet fitted)

Each SitRep also prints indicators the model does not (yet) read. On **every**
update, scan the PDF for any such signal, and if it is present and plausibly
useful, surface it in the PR under a **"Data available but not fitted (for
@seabbs to consider tracking)"** heading with its current value, so the
maintainer can decide whether to add a stream. This keeps new signals from
being silently dropped as the report format evolves. Known candidates:

| Signal (SitRep location) | Note |
|---|---|
| `Taux de suivi des contacts` national % (+ Tableau 4 contacts sous suivi / vus) | Contact-tracing coverage; a surveillance-intensity signal. |
| `Alertes validées — décédées (comm.)` (Tableau 3) | Community suspect-death proxy — the candidate replacement basis for the frozen `suspected_daily_deaths_history` (needs a methodology decision, issue #431). |
| PoE/PoC screening (Tableau 5): travellers screened, alerts, corpses swabbed | Cross-border importation pressure. |
| SMSPS / PPL (Tableau 7/8): front-line-worker infections cumulative, psychosocial follow-up | Health-system-strain signals. |
| CTE bed-capacity strain (§ défis): per-province occupancy vs beds | Local saturation the single national `bed_capacity_history` cannot represent. |
| Symptom-onset epidemic curve (analytique figure) | Digitised separately to `onset_curve_scanned.csv`; not fitted (see above). |
| Alert-investigation throughput (Tableau 3): `Total alertes du jour`, `Alertes investiguées`, `Taux d'investigation (24 h)` | The denominator behind `Cas suspects du jour` — how much of the alert inflow was actually worked. A direct surveillance-effort covariate for suspect ascertainment; moves independently of the validated-suspect count (19–23 July: 82.6%, 84.1%, 79.5%, 79.8%). |
| Occupation table `Total admissions` (cumulative row, distinct from `Total admissions (24 h)`) | Running CTE/CT/CI admission total; a cumulative check on the fitted 24h admission inflow. |
| EDS throughput (§7, per province): death alerts, EDS investigations performed, corpses swabbed | The ascertainment funnel behind the community suspect-death count, rather than another count of it. Gives an observed denominator where the frozen stream's likelihood had to infer one, which reframes issue #431. Always §7 prose, never a table, and present from SitRep 059 — but intermittently: some vintages give both provinces, some only one, and 071 gives no numbers at all, so a missing province is not a zero. Only the both-provinces-together layout is new in 072. **The 059 boundary is deliberate, not the edge of the data.** 059 is where the §7 dashboard layout begins; 058 and earlier print the same quantities in narrative style with spelled-out numbers (058 gives Ituri `Trente-six (36) EDS ont été réalisés`). They are excluded because the earlier era decomposes its denominator differently — 058 reads `Soixante-huit (68) décès étaient à prendre en charge, dont 50 alertes du jour et 18 reports`, so it is unsettled whether `décès à prendre en charge` is the same field as the later `alertes de décès`, or whether `alertes du jour` is. To extend the series earlier, settle that mapping first: assuming equivalence would bury a definitional level shift inside the series, which is the failure mode #431 exists to avoid on the suspect-death stream. Not yet reconciled against Tableau 3's dead-alert total for the same day (98 against 107 on 25 July). |

If a genuinely new indicator appears that is not in this list or the fitted
table, add a row here in the same PR so the procedure stays current.
Before calling anything new, moved or dropped, grep the adjacent vintages for
it (`pdftotext -layout` over `sitrep_pdfs/`), because a two-reader check of one
report cannot see what the neighbouring reports do.
Three claims of that kind in PR #490 turned out to describe changes that had
happened ten or more vintages earlier, or not at all.

Values for these signals accumulate in `candidate_signals.csv`
(`signal`, `sitrep`, `report_date`, `value`, `unit`, `source_note`), one row
per signal per vintage, so a series builds from the day a signal first appears
rather than from the day someone decides to fit it.
Extend it on every update alongside the fitted streams, and open one issue per
signal proposing it for fitting — one per signal, not one per vintage.
Nothing in the model reads this file.

Every value is read twice by independent readers, as for the fitted streams.
The `eds_*` rows are backfilled to SitRep 059; see the EDS row above for why
the series stops there and what would have to be settled to extend it.
**These three series are not comparable across vintages**, because which
provinces print a count changes from report to report: Nord-Kivu appears
numerically almost throughout, Ituri only in 061, 062, 065, 069, 070 and 072,
Tshopo never, and 071 prints no EDS numbers at all.
The value is the sum over the provinces that printed the quantity and the
`source_note` names them, so read that before differencing consecutive rows —
a step can be pure coverage.
Two further cautions in the notes: `EDS réalisés` exceeds the day's alerts in
060, 066 and 067, which the same vintages' carried-over `reports` explain; and
the swab count is usually printed as the combined `corps swabés et sécurisés`
rather than swabs alone.

After editing, validate with the loader and its invariants:

```sh
julia --project=. -e 'using BVDOutbreakSize: load_observations; \
  load_observations(); println("ok")'
julia --project=. -e 'using Pkg; Pkg.test()'   # includes test_load_observations
```
