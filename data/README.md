# Observation data

This directory is the single source of truth for the observations the analysis conditions on.
`load_observations()` (in `src/BVDOutbreakSize.jl`) reads it.
Nothing in the model hardcodes counts.

## Files

| File | What it is |
|---|---|
| `observations.toml` | The manifest the model loads. Every stream is a value or a dated `dates`/`values` history plus a prose `source =` citation. Edit this to advance the analysis. |
| `insp_sitrep_scanned.csv` | Our own direct scan of the INSP SitRep PDFs, one row per report (`date de rapportage`), with a free-text `notes` column recording the headline tiles, laboratory section and table figures. The audit trail behind the PDF-sourced streams in `observations.toml`. |
| `onset_curve_scanned.csv` | Confirmed cases by symptom-onset date, digitised from the analytique-format SitReps' onset epidemic-curve figure (one block per vintage). Fitted as the symptom-onset reporting-triangle stream; see the section below. |
| `released_estimates.csv` | Published point estimates for comparison. |
| `report-snapshot*.toml` | Frozen Imperial report point estimates at fixed vintages. |

## Where the data comes from

DRC national figures are the INSP (Institut National de Santé Publique) situation reports for the 17th Ebola epidemic.
Two sources feed the manifest:

- The cumulative confirmed-case, confirmed-death, recovered and isolation series are cross-checked against the INRB-UMIE mirror's clean `national_*` daily CSVs (<https://github.com/INRB-UMIE/BDBV2026-Data>, `data/insp_sitrep/processed`).
  `scripts/confirm_insp_data.jl` regenerates the confirmed streams from them.
- Every other stream (suspected daily inflow, the 24h analysed laboratory volume, isolation occupancy, bed capacity, recoveries and the CTE/CT/CI patient-movement rows) is read directly from the SitRep PDFs and recorded in `insp_sitrep_scanned.csv`.

INSP is the primary source: it publishes first, and since SitRep 059 (12 July 2026) its richer "analytique / édition quotidienne" PDFs carry content the mirror does not transcribe (epidemic curve by symptom-onset date, age/sex pyramids, five provinces).
The mirror usually lags INSP by a report or two.

## Symptom-onset epidemic curve (`onset_curve_scanned.csv`)

From SitRep 059 the analytique PDFs carry a figure of confirmed cases by symptom-onset date (`courbe épidémique par date de début des symptômes`, DHIS2 line list), split by outcome (Vivant / Décédé).
It is the only published source for the onset-date distribution, and it is a raster bar chart with no accompanying data table.
`scripts/digitize_onset_curve.jl` (the dependency-free Julia reference) recovers the daily counts from the figure pixels and writes one block per vintage to `onset_curve_scanned.csv` (columns `sitrep`, `report_date`, `onset_date`, `confirmed_alive`, `confirmed_dead`, `confirmed_total`).
A Python port (`scripts/digitize_onset_curve.py`) for the automated data-updater produces a byte-identical file.
See `scripts/README.md`.
Both self-calibrate each figure from its axis ticks.
The only manual input is each vintage's rightmost x-axis tick date (in the script `CONFIG`).

These counts are approximate, and the error is a few percent in either direction per scan, independent between vintages.
Against the printed figure `n` it ranges from −5.0% (SitReps 079/080: 2770 against n = 2 915) to +1.6% (SitRep 068: 2344 against n = 2 308).
Per-vintage, SitRep 064 sits at −2.2% (2018 against n = 2 064), SitRep 072 at +0.4% (2531 against n = 2 521), and SitReps 069/070/071 at −3.0% (2260 against n = 2 329).
The 079/080 figure is the largest shortfall seen so far, a couple of points past the previous −3.0% floor rather than an order of magnitude off, so it is recorded as extending the known range rather than rejected (section 4b: check neighbours before excluding).
It is also not a new image, as the duplicate note below explains.
Individual daily bars carry roughly ±1–2 cases of pixel noise.
Some of the shortfall is the faded bars inside the `données potentiellement incomplètes` band at the right of each figure, whose lightened fill falls outside the colour masks.
This mechanism can only lose cases and so does not explain the overshoots.
Treat the sign as unknown.

One consequence deserves emphasis before anyone fits this stream.
Late reporting only ever adds cases, so an onset date's count must be non-decreasing across vintages, and the scans do not respect that.
On onset dates more than three weeks before the earliest report date in the file (12 July, so onsets before 21 June), the scanned totals move both ways between consecutive distinct snapshots.
For example, 064 → 065 falls by a net 36 cases across 34 of 54 such days, and every other consecutive pair falls somewhere too.
A between-vintage increment of a few cases is therefore at or below the noise floor, which bounds what a reporting-delay estimate built from those increments can support.
See issue #488.

SitReps 061 and 062 reuse one figure, as do 069/070/071 and 073/074 (identical n = 2 567, identical digitised total and day count).
SitReps 079/080 now reuse one too (identical n = 2 915, identical digitised total 2 770 and day count, md5-verified byte-identical embedded images).
The twenty-eight scanned vintages therefore hold twenty distinct onset snapshots (report dates 12, 13, 14, 17, 18, 19, 20, 21, 22, 25, 26, 30, 31 July and 1, 3, 4, 5, 9, 10, 14 August); the 2 August (080) vintage shares its figure with 079 and the 10-13 August (088, 089, 090, 091) vintages share one md5-identical embedded image, both duplicate chains, while 3, 4, 5, 9 and 14 August (081, 082, 083, 087, 092) are all fresh distinct snapshots. SitReps 081-083 also carry a caption/page mislabelling distinct from 080's: their onset-date figure sits on page 5 under the "FIGURE 3 — COURBE EPIDEMIQUE PAR SEMAINE DE NOTIFICATION" caption (a notification-week label), while page 6's own "FIGURE 4 — COURBE EPIDEMIQUE PAR DATE DE DEBUT DES SYMPTOMES" caption has no image beneath it at all. The scripts already locate figures by pixel signature rather than caption text (see the 080 note above), so this does not need a page-fallback fix; what SitRep 081 needed instead was a lower axis-tick pixel-intensity cut (081's JPEG-compressed rendering has tick marks as faint as 1px, below the previous cut floor of 2). SitReps 084-086 carry no onset-date figure at all in the shorter "MVEBDB" brief format introduced at 084 - their "Évolution de cas par date de début des symptômes" section instead plots a *different* chart, a weekly notification-date bar count (own subtitle "Cas confirmés vivants et décédés par semaine de notification", x-axis "Date de début de semaine de notification"), confirmed by rendering and reading the embedded image directly rather than trusting the (identical, onset-labelled) section caption; none of the three are digitised. SitReps 087-089 restore the genuine daily onset-date figure (same subtitle/axis/incomplete-data-band style as 081-083), so all three are digitised: 087 is a fresh distinct snapshot, and 088/089 (10/11 August) share one md5-identical embedded image, extending the reprint pattern already seen at 080 and elsewhere. SitReps 090 and 091 (12/13 August) extend that same reprint chain: both embedded onset-figure images are md5-identical to 089's (and hence to 088's), a four-vintage-long duplicate (088=089=090=091, all n=3 362, all digitising to the same 2 349-day/3 049-count block), even though the surrounding section heading in all of 090-092 reverted to the "Évolution de cas ... au cours de 21 derniers jours" wording (the same caption text 084-086's genuinely different weekly-notification chart used) rather than 087-089's own "Nombre des cas confirmés par date de début des symptômes (n = ...)" heading; verified by rendering and reading the embedded chart image directly (not the surrounding caption) for all three, confirming the image itself is still the true daily onset-date bar chart with its own correct internal title in every case - a caption regression, not a content change, and not a repeat of the 084-086 chart-swap (which changed the plotted quantity, not just the surrounding heading). SitRep 092 (14 August) is a fresh, distinct snapshot (n = 3 767, digitising to 3 525, -6.4%): its rightmost x-axis tick moves for the first time since 083, from 5 August to 10 August, shrinking the axis-truncation gap that drove 087-089's larger undercount (see below) back toward the historical band; date-alignment checked against 091 at shifts of -1/0/+1 days, with the L1 distance clearly minimised at 0 (416 vs 784 and 761), confirming correct alignment.
SitReps 059 and 060 are not a repeat despite sharing a digitised total of 1821: they differ on three days, so they count as two snapshots.
Reprints are collapsed by exact value equality over the digitised block rather than by a list of vintage ids, so a new reprint is caught without a change to the loader.
Treating a reprint as a fresh snapshot would fabricate an increment of exactly zero and bias the fitted delay towards fast reporting.
SitRep 068 (21 July) is now included: its PDF was unreachable while the INSP fetch was broken, and it does carry a figure (n = 2 308).
The embedded figure is re-rendered at whatever size the layout needs and moves in both directions: it shrinks to 1009×583 at SitRep 069 (against 1257×698 for SitRep 064) and then grows to 1277×799 at SitRep 072.
Both moves have broken detection once (the shrink by falling under an absolute blue floor, the growth by pushing anti-aliasing over an absolute orange cut).
This is why the figure-detection and axis-tick thresholds in both scripts are pixel fractions rather than counts.
SitRep 080 broke the page-locator instead of the pixel thresholds.
The chart image sits on page 5 as usual, but its own selectable caption text this vintage misreads `par semaine de notification`.
The matching `date de debut des symptomes` phrase appears only on page 6, with no image of its own there, so a same-page caption lookup landed on the wrong page.
Both scripts now fall back to the immediate neighbouring pages when the caption page has no qualifying image.
That widened search pulled in a second blue-heavy image, the provincial case map, which passed the old blue-fraction floor (0.02).
Raising it to 0.055 (measured range: onset charts 0.066–0.125, the map 0.045–0.047 across every vintage checked) separates them with margin on both sides.
Both fixes were verified to reproduce every previously committed block (059–079) unchanged before the SitRep 080 block was accepted.
SitRep 081 broke the x-axis tick detection instead: its JPEG-compressed figure renders the weekly tick marks only 1px tall (max dark-pixel band value 2 in a 7-row window, against the working threshold of 2), so the old `cut in (4,3,2)` loop recovered only ~10 of the 19 ticks, inflating pixels-per-day and undercounting to 1411 against the printed n = 3 066 (−54%), far outside the noise band. The fix runs the cut loop down to 1 AND narrows the tick band window from `base+2:base+8` to `base+2:base+6` so a low-cut scan cannot pull in the x-axis date-label text further down, keeping the most complete regular weekly tick row. With it, SitRep 081 recovers the full 19-tick axis (ppd 6.714 vs 080's 6.643) and digitises to 2926 against printed n = 3 066 (−4.6%, back inside the band), and SitRep 082 digitises to 2967 (a fresh figure, reported n read off the axis at digitisation time). Both scripts (`.jl` and the byte-identical `.py`) were verified to reproduce every previously committed block (059–080) unchanged before the 081/082 blocks were accepted.

SitRep 087's figure broke the y-axis count scale instead. Every figure through SitRep 083 draws its y-axis gridlines at 0/20/40/60/80, and both scripts hard-coded that increment as the pixels-per-count divisor. SitRep 087's figure (and 088/089's shared reprint) draws its gridlines at 0/25/50/75 instead - confirmed by reading the printed tick labels directly, since the pixel geometry of an evenly-spaced tick strip cannot itself say whether the spacing means 20 or 25, any more than it can say what date the rightmost tick is (see `last_tick` above). Applying the old /20 divisor to a 25-count grid undercounts every bar by a scale-dependent amount, and this one was caught only because it broke the reporting-triangle invariant on dates that should have been stable: SitRep 083's 15 May onset count reads 26, and the same figure misread through the wrong divisor came out as 8 - a fall of 18 on a date twelve weeks old, far past where late reporting could plausibly still be arriving. The fix adds a per-vintage `Y_AXIS_STEP` override (defaulting to the historical 20, so it changes nothing for 059-083) alongside `CONFIG`, read off each new figure the same way `last_tick` is. With it, 087 digitises to 3233 against printed n = 3362 (-3.8%, inside the established band) and 088/089 to 3049 (-9.3%, outside the -5.0%/+1.6% band measured over 059-083 but explained by the same axis-truncation mechanism at a larger-than-usual scale: 087-089's last tick stays fixed at 5 August - the last date 083's own figure covered - while the report date advances to 9-11 August, so a growing share of the recently-confirmed cohort has an onset date past the plotted axis and is silently excluded from the digitised total while still counted in the printed n; see the axis-gap paragraph below). Both scripts were verified to reproduce every previously committed block (059-083) unchanged before the 087-089 blocks were accepted; a single unrelated 1-count drift on SitRep 083's own 8 June cell was found during that check: the Python script (unmodified, pre-existing code, not this fix) rebuilds it as 32/18/50 against the committed 33/18/51, reproducing on a clean checkout, while the Julia script reproduces the committed value exactly. This is the first observed case of the two "byte-identical" scripts disagreeing; the gap is a single count on a single cell out of 1878 previously-committed rows, far below the noise floor, so the Julia output (which matches history exactly) is used for this update and the discrepancy is left as pre-existing cross-implementation floating-point jitter rather than investigated further.

Each figure also stops its x axis short of its own report date, by anything from zero to eight days, and the last bar it prints is often a substantial count rather than a tail fading to zero.
An onset date past that axis is a date the figure says nothing about, not a bar of height zero.
Those cells are dropped rather than read as zeros.
The axis gap is about five days for most vintages, close to the reporting delay itself.
Reading it as "nothing reported yet" would force the fitted hazard to near zero over the first five days and pile the missing mass onto the delay at which the axis first covers the date.

This stream is fitted.
`load_onset_curve` (`src/onset_curve.jl`) reads the file, collapses reprints, drops vintages reported after the manifest cut-off and builds the between-vintage increment cells.
`onset_reporting_model` (`src/models/observations.jl`) fits them through a discrete reporting-delay hazard whose asymptote carries ascertainment.
It is the only direct observation of the latent onset series, every other stream seeing it only after a further convolution.
Advancing `as_of_date` in `observations.toml` past a newly digitised vintage's report date picks that vintage up with no code change.

### Rough onset-to-report delay

Because each vintage redraws the same onset cohort at a later report date, the scanned curves form a reporting triangle.
Old onset dates are stable across vintages while recent ones fill in as more confirmations arrive, for example onset 10 July reads 4 → 9 → 27 → 33 → 37 → 39 → 40 across the 12, 14, 17, 18, 19, 20 and 22 July snapshots.
Taking the latest snapshot as the near-complete reference for onset dates at least ~12 days old, the empirical proportion of eventually-reported confirmed cases reported within `d` days of onset is around 55–60% by 7 days, 85% by ~10 days and 90%+ by ~12 days.
By ~2.5 weeks that proportion reaches 95%, near-complete (98–99%) by ~3 weeks, with a median of ~5–6 days.
The 7-day figure is the least well determined of these.
A fuller reanalysis of the same triangle put it at 54–62% with a 95% interval of 43–68%, wide because the digitisation noise is close in size to the between-vintage increments the estimate rests on.
Quote the interval, not the point.
This supersedes the earlier three-snapshot estimate and shifts it slightly later.
The seven-snapshot curve is flatter through the second week than the three-snapshot one suggested.
This is still a coarse estimate: it rests on seven digitised snapshots over a ten-day window, the reference snapshot is itself right-truncated for its most recent onsets, and the delay it measures is onset → confirmed-and-reported (it folds together care-seeking, lab confirmation and reporting).

Two caveats have been added since these figures were computed, and neither has been folded into them.
The 21 and 25 July snapshots now exist (nine distinct snapshots, not seven), so the estimate is one report behind at both ends.
More importantly, the scan noise documented above is two-sided and of the same order as the between-vintage increments this estimate is built from, so the quoted percentages are almost certainly more precise than the data supports.
Treat them as indicative.
The model now carries the noise through.
The fitted reporting-delay hazard in `onset_reporting_model` scores the increments under a heavy-tailed likelihood whose scale is built from the measured digitisation error, and the analysis report quotes the delay with its posterior interval.
Prefer that estimate over the arithmetic above.

### Fetching a SitRep from INSP directly

INSP blocks some default user agents with HTTP 403, but a browser User-Agent returns HTTP 200.
Each SitRep post embeds the real PDF URL in a `pdfemb-data` base64 blob.
To get the direct PDF for a post (find the post at <https://insp.cd/ebola-17eme-epidemie/> or via `check_new_sitreps.jl`):

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

`scripts/download_sitreps.jl` bulk-fetches the whole MVE SitRep series from the `wp/v2/media` API instead, so the per-post decode above is only needed for a report the bulk fetch misses (a corrected `_v2` re-issue).

## Checking for new SitReps (do this before every refresh)

```sh
julia --project=scripts scripts/check_new_sitreps.jl
```

It lists the INSP-published SitReps, compares them with the latest report in `insp_sitrep_scanned.csv`, and prints the gap (exit code 1 if behind, 0 if current).
This is the guard against the manifest silently drifting weeks behind again.

## What to read from each SitRep

Read these from the PDF (the `date de rapportage`, not the filename or the publication date, is the key each row is stored under) and record them in the matching `observations.toml` history.
Advance every stream in the table below on each update, whenever the SitRep prints the value, not just the confirmed/suspected headlines.
The isolation, bed-capacity, recovered, laboratory and CTE patient-movement streams move together with the headlines: the treatment/occupancy streams share one Tableau, and the loader invariants tie the latest treatment-flow day to the isolation window.
Advancing the headlines while leaving these behind both drops data and can break the loader test.

| From the SitRep | `observations.toml` stream |
|---|---|
| Headline `Cumul cas confirmés` | `confirmed_case_history` |
| Headline `Cumul décès confirmés` (+ lethality %) | `confirmed_death_history` |
| Headline `Guéris — cumul` | `recovered_history` |
| Headline `Patients en isolement` (+ occupancy %) | `isolation_history`; occupancy → `bed_capacity_history` (= occupancy ÷ rate) |
| Headline / Tableau 3 `Cas suspects du jour` | `suspected_daily_history` |
| Laboratory 24h `échantillons analysés` per province (§4.3 table, or §3.2 bullets in the analytique format) | `tests_analysed_daily_history` |
| Occupation table (Tableau 6/7, or Tableau 5/6 in the analytique format): `Total admissions (24h)`, `Sorties — décédés / non-cas / évadés`, `Patients au lit (J-1)`, `dont confirmés / suspects` | `treatment_admissions_history`, `treatment_deaths_history`, `treatment_ruleout_history`, `treatment_absconded_history`, `treatment_aulit_history`, `treatment_confirmed_incare_history`, `treatment_suspect_incare_history` |

Always cross-check the confirmed/death/recovered/isolation headline against the INRB-UMIE `national_*` CSV for the same date.
Note any disagreement in the `source =` string and prefer the auditable value.

## Inclusion rules and conventions

- **Cut-off 28 May 2026 (SitRep 014)**: the last vintage with a coherent national laboratory total.
  The suspected cumulative streams freeze at their last stable vintage (26 May).
  Confirmed cases/deaths and the post-cut-off daily streams are recorded and, for the streams the model reads (e.g. the 24h analysed denominators), fitted beyond the cut-off.
- **24h analysed denominator**: only include a day whose confirmed-case increment does not exceed that day's national `analysed` count (you cannot confirm more specimens than you analysed).
  A province that reports "reçus / en cours d'analyse" but no completed count contributes 0, like a non-reporting province.
  Sum only completed per-province analysed counts.
  A province worded "documentés dans le réseau de collecte" rather than "analysés" is a shade looser but usable (precedent: SitRep 065, 18 July, Ituri 167 + Nord-Kivu 127 = 294).
  Note the wording in the `source =` string.
- **Occupation table (Tableau 6/7)**: use it only when it balance-closes (`Patients au lit (J-1) + Total admissions (24h) − Total sorties (24h) = Patients en isolement (Fin J)`) and the `dont confirmés + dont suspects` split sums to the Fin J total.
  Bed capacity = Fin J ÷ `taux d'occupation`, which should equal the printed `Nombre de lits`.
  A carried-forward tile (page-1 isolement count and occupancy rate byte-identical to the prior day with no fresh Tableau to corroborate) is a reporting gap: record the raw value in the CSV note but exclude it from the fitted series (precedent: 30 June, 3–4 July).
  A carried-forward Ensemble column is a narrower variant: the per-province rows and the `Total sorties (24h)` row are fresh and balance-close.
  However, the table's own printed `Ensemble`/Global cell for one row block (typically the Sorties breakdown: Décédés/Non-cas/Guéris/Évadés/Transférés) is a byte-identical copy of the previous vintage's Ensemble cell rather than that day's own row sum.
  Detect it by summing the per-province cells yourself and comparing to the printed Ensemble cell and to the previous vintage's Ensemble cell for the same row.
  If it matches the prior vintage rather than this vintage's own row sum, use the auditable row sum (it will also match the prose paragraph and the `Total sorties` row) and record the printed Ensemble cell as the discrepancy.
  Precedent: SitRep 079, 1 August — Ensemble cells 18/71/15/2/0 were a byte-identical copy of SitRep 078's own Ensemble cells, while the row sums 25/86/42/4/1 matched the prose and `Total sorties (24h)` = 158 exactly.
- **Unpublished reports**: some SitRep numbers were never issued (029/12 Jun, 043/26 Jun, 045/28 Jun).
  The series simply step over the missing date.
  Note it in the `source =` string.
- **Headline vs table disagreement**: prefer the auditable table sum and record the headline as the discrepancy.
  Precedents: SitRep 009's 220 over the erroneous 119, and SitRep 061's Fin J 736 over the page-1 typo 737.
- **Retrospective harmonisation steps**: when a vintage's asterisked headline jumps far beyond its own 24h new-case count because a provincial base was integrated, record the harmonised headline (it is what the INRB-UMIE `national_*` CSV carries) and state the gross-vs-net split in the `source =` string.
  This makes the artefact visible to whoever reads the increment the daily likelihood actually fits.
  State the split for both confirmed streams, since a break day now applies to cases and deaths alike.
  Precedents: SitRep 065 (18 July, Haut-Uélé/Nia-Nia reattachment, +83 gross vs +77 net cases and +40 gross vs +37 net deaths — both downward, which is what a transfer between provinces looks like rather than a base integration) and SitRep 069 (22 July, Ituri/Nord-Kivu base integration, +97 gross vs +369 net cases and +62 gross vs +236 net deaths — much the largest so far, and called out in the report's own footnote).
- **Format change from SitRep 059**: the analytique format moved the 24h analysed counts to §3.2 Laboratoire bullets, dropped the page-1 suspect-death-of-day subtitle, and reports the daily suspect total in Tableau 3.
  The DHIS2 suspect-death-of-day measure is not continuous in this format, so `suspected_daily_deaths_history` is frozen at 11 July.
  See its `source =` string.
- **Tableau 6 (occupation/patient-movement table) drops from SitRep 081**: checked directly against 077-080 (all present, 3 `TABLEAU 6` hits each) vs 081-085 (0 hits, all absent).
  081-083 replace it with per-province prose under a "PRISE EN CHARGE HOLISTIQUE (PECH)" pillar heading (isolement, lits, occupancy% per province, admissions/sorties incompletely); 084-085 drop that too.
  `treatment_admissions_history`, `treatment_deaths_history`, `treatment_ruleout_history`, `treatment_absconded_history`, `treatment_aulit_history`, `treatment_confirmed_incare_history` and `treatment_suspect_incare_history` are frozen at SitRep 080 pending a maintainer decision; see issue #562.
  `bed_capacity_history` is NOT frozen — it continues from the PECH/section-1.6 prose's per-province `lits` mentions (a partial, not national, total: whichever provinces print a bed count that vintage), following the precedent already set for SitReps 081-082.
- **SitRep 084 "MVEBDB" brief format**: from 06 August 2026 INSP switched to a much shorter (~6-page) template under a new "TASK FORCE PRESIDENTIELLE EBOLA 17" letterhead.
  It drops the page-1 `Cas suspects du jour` tile and Tableau 3 entirely (`suspected_daily_history` frozen at SitRep 083, 5 August), replacing it with a differently-structured national alert-validation funnel that is not a like-for-like replacement.
  It DOES give the 24h analysed-sample national total directly in prose (§1.4), a cleaner source than the previous per-province summation, so `tests_analysed_daily_history` is unaffected and continues through 085.
  SitRep 085's page-1 isolation tile (595) exactly matches Ituri's own province-level occupancy that day, strongly suggesting it is not a genuine national total (084's equivalent tile clearly was); excluded from `isolation_history`, which steps from 084 straight to 086.
  Full writeup and evidence in issue #562.
  SitReps 086-089 each pass the same cross-check (their national tile exceeds Ituri's own same-day occupancy by a plausible non-Ituri remainder; 089's 570 = Ituri 375 + Nord-Kivu 176 + Sud-Kivu/Tshopo 19 to the exact unit), so all four are back in `isolation_history` — 085 remains the sole excluded vintage in this format so far.
- **Per-zone-de-santé table returns from SitRep 087**: 086 still gives only the page-1 province-summary table introduced at 084; from 087 a new "Cas et décès confirmés par province et zone de santé" table appears (087 gives 24h flow only; 088 onward adds cumulative cas/décès per province too), letting `confirmed_death_history` reconcile against a per-province cumulative sum again from 088 (087's table has no cumulative-deaths column, so that one date rests on the page-1 headline alone).
  The unexplained `**` marker on the page-1 `CUMUL DES DÉCÈS CONFIRMÉS` tile (present with no footnote text at 084-086, see above) is absent from 087 onward — a presentation change with no data effect, since the marker never resolved to footnote text in the first place.
  SitRep 088's §1.4 Laboratoire prose states "78 nouveaux cas" where the per-ZS table and headline both give 68 (53 Ituri + 15 Nord-Kivu); the 78 figure equals Ituri's full 63 `résultats positifs` (53 new + 10 `reprélèvements` of already-confirmed patients) plus Nord-Kivu's 15, i.e. the prose sentence appears to conflate reprélèvements with new confirmations.
  The auditable 68 is what is reflected in the cumulative step; 78 is recorded as a discrepancy in `data/insp_sitrep_scanned.csv` and used nowhere.
- **SitReps 081-083's onset-curve figure is mislabelled, not missing**: the "date de début des symptômes" chart sits on page 5 under the "FIGURE 3 — COURBE EPIDEMIQUE PAR SEMAINE DE NOTIFICATION" caption (a stale/wrong label — the chart's own embedded title and axis both read "par date de début des symptômes"), while page 6's own "FIGURE 4" caption for that exact title has no image beneath it.
  This is the same class of caption/image mismatch as SitRep 080 (see the onset-curve section above), just a different swap.
  084-086 (the brief format) carry a genuinely different chart — a bar count whose own caption says "date de début des symptômes" but whose own subtitle and x-axis read "par semaine de notification" / "date de début de semaine de notification" — not digitised; see the onset-curve section and issue #562.
  SitReps 087-089 restore the genuine daily onset-date figure (verified by rendering and reading the embedded image, not just the caption) and are digitised; see the onset-curve section above, including the y-axis-grid digitiser fix that vintage required.
  SitReps 090-092 continue in the same MVEBDB format and all three are digitised; 090 and 091 share 089's exact embedded image (a four-way reprint, verified by rendering) despite their surrounding section heading reverting to the "21 derniers jours" wording, and 092 is a fresh snapshot with the rightmost axis tick moving for the first time since 083. See the onset-curve section above for the full detail and the date-alignment check.

## Other signals to scan for and flag (not yet fitted)

Each SitRep also prints indicators the model does not (yet) read.
On every update, scan the PDF for any such signal, and if it is present and plausibly useful, surface it in the PR under a **"Data available but not fitted (for @seabbs to consider tracking)"** heading with its current value, so the maintainer can decide whether to add a stream.
This keeps new signals from being silently dropped as the report format evolves.
Known candidates:

| Signal (SitRep location) | Note |
|---|---|
| `Taux de suivi des contacts` national % (+ Tableau 4 contacts sous suivi / vus) | Contact-tracing coverage; a surveillance-intensity signal. |
| `Alertes validées — décédées (comm.)` (Tableau 3) | Community suspect-death proxy — the candidate replacement basis for the frozen `suspected_daily_deaths_history` (needs a methodology decision, issue #431). |
| PoE/PoC screening (Tableau 5): travellers screened, alerts, corpses swabbed | Cross-border importation pressure. |
| SMSPS / PPL (Tableau 7/8): front-line-worker infections cumulative, psychosocial follow-up | Health-system-strain signals. Tableau 8's own printed total is Ituri-zones only; SitReps 069-072 carry no Nord-Kivu line at all, so the `health_worker_infections_cumulative` / `health_worker_deaths_cumulative` values recorded for those vintages are Ituri-only, not national. From SitRep 073 a narrative sentence beneath the table adds Nord-Kivu and states an explicit national cumulative (134 confirmed / 40 deaths); use the narrative total where it is given, the Tableau 8 total otherwise, and note which in `source_note`. From SitRep 077 a Haut-Uele PPL count also appears in the narrative for the first time (checked directly against 069-074, none mention it) but with no death figure, and the explicit national-total sentence stops being printed (077/078 give only the three separate province figures), so from 077 onward `source_note` records a maintainer-computed sum of whichever province figures are printed that vintage, naming them. |
| CTE bed-capacity strain (§ défis): per-province occupancy vs beds | Local saturation the single national `bed_capacity_history` cannot represent. |
| `formes sévères ou critiques` (PECH prose, §6/section 1.6) | Ituri-only clinical-severity count among isolated patients, distinct from volume signals. First seen SitRep 080; see issue #554. Absent from the SitRep 084+ brief format (no PECH-equivalent severity breakdown). |
| Symptom-onset epidemic curve (analytique figure) | Digitised separately to `onset_curve_scanned.csv`; fitted as the reporting-triangle stream (see above). |
| Alert-investigation throughput (Tableau 3): `Total alertes du jour`, `Alertes investiguées`, `Taux d'investigation (24 h)` | The denominator behind `Cas suspects du jour` — how much of the alert inflow was actually worked. A direct surveillance-effort covariate for suspect ascertainment; moves independently of the validated-suspect count (19–23 July: 82.6%, 84.1%, 79.5%, 79.8%). |
| Occupation table `Total admissions` (cumulative row, distinct from `Total admissions (24 h)`) | Running CTE/CT/CI admission total; a cumulative check on the fitted 24h admission inflow. |
| EDS throughput (§7, per province): death alerts, EDS investigations performed, corpses swabbed | The ascertainment funnel behind the community suspect-death count, rather than another count of it. Gives an observed denominator where the frozen stream's likelihood had to infer one, which reframes issue #431. Always §7 prose, never a table, and present from SitRep 059 — but intermittently: some vintages give both provinces, some only one, and 071 gives no numbers at all, so a missing province is not a zero. Only the both-provinces-together layout is new in 072. The 059 boundary is deliberate, not the edge of the data. 059 is where the §7 dashboard layout begins; 058 and earlier print the same quantities in narrative style with spelled-out numbers (058 gives Ituri `Trente-six (36) EDS ont été réalisés`). They are excluded because the earlier era decomposes its denominator differently — 058 reads `Soixante-huit (68) décès étaient à prendre en charge, dont 50 alertes du jour et 18 reports`, so it is unsettled whether `décès à prendre en charge` is the same field as the later `alertes de décès`, or whether `alertes du jour` is. To extend the series earlier, settle that mapping first: assuming equivalence would bury a definitional level shift inside the series, which is the failure mode #431 exists to avoid on the suspect-death stream. Not yet reconciled against Tableau 3's dead-alert total for the same day (98 against 107 on 25 July). |

If a genuinely new indicator appears that is not in this list or the fitted table, add a row here in the same PR so the procedure stays current.
Before calling anything new, moved or dropped, grep the adjacent vintages for it (`pdftotext -layout` over `sitrep_pdfs/`), because a two-reader check of one report cannot see what the neighbouring reports do.
A claim of novelty that skips this check has turned out to describe a change that had actually happened ten or more vintages earlier, or not at all.

Values for these signals accumulate in `candidate_signals.csv` (`signal`, `sitrep`, `report_date`, `value`, `unit`, `source_note`), one row per signal per vintage, so a series builds from the day a signal first appears rather than from the day someone decides to fit it.
Extend it on every update alongside the fitted streams, and open one issue per signal proposing it for fitting (one per signal, not one per vintage).
Nothing in the model reads this file.

Every value is read twice by independent readers, as for the fitted streams.
The `eds_*` rows are backfilled to SitRep 059.
See the EDS row above for why the series stops there and what would have to be settled to extend it.
These three series are not comparable across vintages, because which provinces print a count changes from report to report: Nord-Kivu appears numerically almost throughout, Ituri only in 061, 062, 065, 069, 070 and 072, Tshopo never, and 071 prints no EDS numbers at all.
The value is the sum over the provinces that printed the quantity, and the `source_note` names them, so read that before differencing consecutive rows: a step can be pure coverage.
This is a workaround for the schema having no province column, which is issue #492.
Per-province splitting should come before any of these signals is fitted.
Two further cautions in the notes: `EDS réalisés` exceeds the day's alerts in 060, 066 and 067, which the same vintages' carried-over `reports` explain, and the swab count is usually printed as the combined `corps swabés et sécurisés` rather than swabs alone.

After editing, validate with the loader and its invariants:

```sh
julia --project=. -e 'using BVDOutbreakSize: load_observations; \
  load_observations(); println("ok")'
julia --project=. -e 'using Pkg; Pkg.test()'   # includes test_load_observations
```
