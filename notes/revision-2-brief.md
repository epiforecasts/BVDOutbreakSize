# Report revision 2 — full brief (do not lose)

Captures the large review batch. Obey `docs/src/contributing.md` prose rules.
Group A = model/science (needs refits + validation). B = methods prose/structure.
C = plots. D = results/summary. E = process/release. Cross-refs to issues #224/#225.

## A. Model / science changes (refit + validate each)

1. **Growth-rate prior, transform FORWARD.** Put the prior on the growth rate /
   doubling time (like main: `r ~ LogNormal(log(log2/20), 0.15)`), using the
   genetic estimate, and derive the first Rt forward via Euler-Lotka — NOT a
   prior on R0=1.6 transformed back. The genetic report (cuomodannenburg2026)
   gives R0 mean 1.31-1.55 with wide uncertainty; our 1.6 came from a different
   GT. We want the growth-rate prior slightly inflated in SD but unbiased vs the
   source. In the Rt model section show the assumed dependence for the initial
   growth rate; give the growth-rate prior itself in the follow-up
   (seeding/growth) section where it makes sense. Closes the #223 direction.
2. **Anchor = 14 days after TMRCA** (was 7), given TMRCA uncertainty; increase
   the doubling-count prior centre to m ≈ 4 (was 2) to match.
3. **Generation interval = Gamma shape/scale** from the source (bdbv line-list /
   the cited GT source) WITH the source's reported uncertainty — not mean/SD
   moment-matched LogNormal. Gamma shape+scale is AD-stable. Close #224 and
   remove the in-text issue flag; describe what we do and why it's a limitation
   (we inflate/assume the GT uncertainty ourselves — state that explicitly).
4. **Test sensitivity prior**: renewal currently `Beta(30,2)` (mean ~0.94, too
   sensitive); should be the LESS-sensitive `Beta(6,2)` (mean 0.75) as on main —
   the BDBV assay sensitivity is lower. Fix + correct the prose (it wrongly says
   "just below the GeneXpert whole-blood clinical sensitivity"; check main's
   parameterisation + citations).
5. **Naming**: stop calling the headline `C_T`. We report cumulative infections
   (the running sum of I_t). Refer to it as cumulative infections / I_t, not
   C(T) — it is confusing.
6. **Forecast**: let Rt continue to EVOLVE over the 7-day horizon (it currently
   carries the current growth rate forward as no-change). Use the same forecast
   settings statement; don't say "no further interventions and no saturation"
   in a way that implies Rt is frozen.

## B. Methods prose + structure

- **Epi-process section reorg**: infections->onsets (incubation convolution)
  belongs IN the epidemiological-process section, not its own thing — this
  removes the duplicated incubation-period definition (currently in onset
  incidence AND elsewhere). Structure: Epi processes = incubation period, then
  infections->onsets model, then the other delays + CFR. Open the epi section
  with a short para outlining what's in it (like the infections subsection does).
- **Infections subsection opener** reword: "The infection process is made up of
  several processes. There are ... . Each is described in a subsection below."
  (not "built in generative order: the reproduction number, the generation
  interval that drives ...").
- "parameters feed each observation" -> "inform", not "feed".
- Rt model: "the reproduction number follows" -> "is assumed to follow". Make
  clear Rt is INTERPOLATED between knots (piecewise linear) — not clear in text
  or maths. Remove "rather than switching at a single date". Remove "The seed
  magnitude does not depend on the growth rate, so r enters only the renewal
  that follows."
- Define the ANCHOR clearly somewhere ("The anchor is the grid day N days after
  the genetic TMRCA ...").
- **Define PMF on first use.** Cite **EpiNow2** (correct citation) early when we
  introduce the kind of model — it has many similarities (NOT the joint fitting).
- **Citations on first use, then refer back** — no forward references. Sources
  cited when a parameter is first defined; methods referred back to, not forward.
- **Delays**: the lognormal/AD/moment-match note (if true that not all are
  lognormal/moment-matched and others are AD-stable) goes in the FIRST delay
  section and is referenced back, not in the overview. Each delay must say its
  SOURCE and the uncertainty (cite on first use). State explicitly where the GT
  and other delay uncertainties are OUR assumption vs from the source.
- **Onset-to-report**: missing source/uncertainty — add.
- **Onset-to-death**: cite that McCabe et al. use t in their onset-to-death;
  use a Gamma shape/scale (as reported) not mean/SD.
- **Onset-to-detection** -> rename "Onset-to-hospitalisation"; state it's the
  assumption used in the export model.
- **Shared observation submodels** intro: shorten (it's repeated in the
  subsections). Say "several parameters are assumed shared across ...: these
  are ... . More detail in the subsections." Reword "ties the count likelihoods
  together; the ascertainment fractions scale ..." -> "we assume the ... datasets
  are overdispersed and share a common ...; the laboratory testing priors are
  shared between ...". Remove "An independent alternative drops the shared
  hyperprior and gives each system its own fraction." Fix "fraction is informed
  by essentially a single aggregate data point" (factually wrong now).
- **Lab priors**: change "explicitly" wording. Severity enrichment intro: "We
  assume more severe, more-likely-Ebola cases are preferentially tested, via an
  enrichment factor (δ0)." Remove "which falls as the background grows" (BVD
  cases also grow). Remove the confusing "How the per-window confirmation
  probability ties to the suspect-pool composition ... is set out ... below"
  sentence. The sensitivity-prior prose is out of date (see A4). The
  "early confirmed vintages (18-23 May) have no per-vintage ..." is out of date
  (we have more recent data too) — update. Define how the modelled laboratory
  volume is modelled (it's referenced as "against the modelled laboratory
  volume" without definition). We only need to define WHICH model we use — drop
  the composition-vs-positivity branch description. Say WHY the deaths are
  "enriched on the odds scale by m_death".
- **Observation likelihoods (suspected cases, confirmed cases, suspected deaths,
  confirmed deaths, exports)**: give the FULL maths of each model (equation
  blocks, not maths-in-text), using SUM notation for the convolution, not `conv`.
  Show the travel scaling in the export maths; show the zero-exports-before-
  first-observation term maths.
- **Confirmed deaths**: remove the #225 flag and the "cases-like lab process ...
  flagged in issue #225" sentence; close #225. Remove the "thinning ties ..."
  claim (not true). Say why m_death enrichment.
- **Exports**: state the one-way-travel assumption (we model travel out, not
  return, so likely overestimates infections on its own). Show the travel-scaling
  maths and the zero-before-first-detection maths.

## C. Plots

- **Rt over time**: higher trajectory alpha (currently too faint); plot 30/60/90
  credible ribbons, NOT a median line; weekly dates on x; upper-bound the y-axis
  at ~2x the 90% bound (rounded) so random trajectories don't force whitespace,
  while still showing the R=1 line; the horizontal grey dashed R=1 line is
  currently missing — restore it.
- **Surveillance posterior over time plot**: plot 30/60/90 credible intervals,
  not the median / current intervals.
- **Forecast**: Rt evolves (see A6).
- **Delay sensitivity**: do it over the DIFFERENT DELAYS as on main (read main's
  delay-sensitivity), NOT a made-up shorter/longer scenario.
- **Clock-rate sensitivity**: state what the clock rates ARE.
- **Cumulative infections by data stream**: scale the x-axis by a multiple of the
  joint-fit 90% bound (the confirmed-cases-only stream is ill-defined and useless
  on its own; the current cap hides the mass).
- **Estimate evolution**: released estimates are frozen values -> connect them
  with a line/ribbon; show 30/60/90 credible intervals for releases AND our
  refit estimates; add our CURRENT model+data estimate across the period as a
  ribbon. Colours: red = renewal refit, blue = released estimates, another colour
  = current-data+model estimate. Remove McCabe scenario range backdrop here.
- **McCabe comparison**: download the McCabe 18-May and 20-May reports, extract
  their scenario UNCERTAINTY, and show the scenarios WITH uncertainty. Compare to
  the date-1 and date-2 McCabe reports AND our earlier released estimates.
- **Delays pairplots/tables**: pairplots + tables are missing for many delays.
  Incubation period belongs in the infection model. Consider a delays table
  before the surveillance section.

## D. Results / summary

- **Restore the prior-vs-posterior "shift" reporting** that was on main (in the
  results summary) — we lost it.
- Summary wording: not statements — "X is estimated to ...". "confirmed cases
  capture only a small share" -> "... are estimated to ...". "the outbreak is
  estimated to have begun on ...". "the initial growth rate is estimated to have
  been ..." (same for Rt).
- **Frozen/predictive variants wording**: "The frozen and predictive variants
  reuse the same structure ..." — readers don't know these yet (no methods
  section). Reword to "other model variants reuse these models with different
  amounts of data". The "Each frozen re-fit is a reduced fit of 500 draws ..." —
  use the SAME settings statement consistently. Reword the weird "The same
  frozen-refit approach underlies the matched-in-time comparison ...".

## E. Process / release

- **Merge new main**: 6 June data (#229), provisional-estimate warning (#230),
  news tidy (#228). Keep renewal's loader field structure with main's data.
- **Make renewal the 1.4.0 release**: own news section showing diffs vs main's
  1.3.0 — many changes, keep relatively concise.
- **Close issues #224 and #225** (describe the limitations in-text instead).
- After all the above: PR up to date, rebuild docs, host locally.
