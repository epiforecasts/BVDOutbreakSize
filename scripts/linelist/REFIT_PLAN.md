# Refit plan: pinned breakpoint, onsets comparison

Written 2026-08-20. Delete this file once the refit is done and the results are read.

Everything below can be run from a cold session. Nothing in it depends on state held in memory.

## Why a refit is needed

The intervention breakpoint was moving with the data source.

`default_breakpoint(obs)` is `obs.n - obs.who_first_sitrep_days`, and `who_first_sitrep_days` is `n - reported_history.days[1] + 1` (`src/data.jl:359`), so the breakpoint is the grid day before the first `reported_case_history` vintage.
`reported_case_history` is one of the three streams the line list replaces.
So although `onsets` never conditions on that stream, substituting it moved the breakpoint:

| | first reported vintage | breakpoint |
|---|---|---|
| sitrep | 2026-05-18 | 94 |
| linelist_known | 2026-05-26 | 102 |

Eight days apart. A difference in R_t between the two arms was therefore not attributable to the onset data alone.

`fit_single.jl` now takes the breakpoint from the released manifest for every data mode, so both arms sit at 94.
The released manifest's first reported vintage is the first WHO joint situation report; the line list's is only the date its export begins.
`--breakpoint=N` overrides if a run needs a different one.

The triangle-size difference between the two sources is deliberately left alone.
The line-list triangle carries 904 cells over 35 loaded vintages against 419 over 19, and that is part of what the comparison is measuring, not a confound to remove.

## What this invalidates

The breakpoint is now part of the chain cache key (`..._bp<N>.jls`), because it changes the model and a chain fitted at one breakpoint must not be silently reused under another.

Every existing chain therefore has an old-style name and will refit.
For the `sitrep` arms that is avoidable: their breakpoint is unchanged at 94, so the cached chains are still valid and can simply be renamed. See the optional step below.

The `linelist_known` arms must genuinely refit: their breakpoint changes from 102 to 94.

## Environment

```bash
cd ~/Documents/Github/BVDOutbreakSize
export LINELIST_INPUT_DIR=$PWD/ignore/linelist/inputs
export LINELIST_RELEASED_MANIFEST=$PWD/data/observations.toml
export LINELIST_AS_OF=2026-08-10
export BVD_DELAY_DIR=$HOME/Documents/Github/bvd-internal-cmmid/results   # stage 2 only
```

`LINELIST_AS_OF` is required and must be the same for every run in the comparison.

## Optional: keep the sitrep chains

Only do this if the check passes. It saves about 30 minutes on stage 1 and about two hours on stage 2.

```bash
cd ignore/linelist/chains
# Confirm the sitrep breakpoint really is 94 before trusting the rename.
# A 10-sample pilot prints it and costs a few minutes.
cd ~/Documents/Github/BVDOutbreakSize
LINELIST_OUT_DIR=$(mktemp -d) julia -t 2 --project=docs \
  scripts/linelist/fit_single.jl onsets --pilot --data=sitrep --delays=repo 2>&1 \
  | grep -A3 "Info: breakpoint"
```

Expect `breakpoint = 94`, `native_to_this_manifest = 94`, `pinned = false`.
If that is what it prints:

```bash
cd ~/Documents/Github/BVDOutbreakSize/ignore/linelist/chains
for f in onsets_sitrep_*_2026-08-10.jls; do
    [ -e "$f" ] || continue
    cp "$f" "${f%.jls}_bp94.jls"
done
ls -la onsets_sitrep_*_bp94.jls
```

`cp` rather than `mv`, so the originals survive a mistake.

## Stage 1: the comparison

Both arms at `repo` delays. This is the result, not a sensitivity: `onsets` estimates its reporting delay from the triangle rather than taking a prior, so at `repo` it uses no cmmid estimate at all and takes no view on a generation interval fitted from line-list transmission pairs.

```bash
cd ~/Documents/Github/BVDOutbreakSize
scripts/linelist/run_grid.sh --stage=1 --dry-run    # confirm what will run
scripts/linelist/run_grid.sh --stage=1
```

Grid: `onsets × {sitrep, linelist_known} × repo`.

Time: about 30 minutes for sitrep and 50 for linelist_known, so roughly 1.5 hours, or 50 minutes if the sitrep chain was renamed above.

The script detaches itself with `nohup caffeinate -is` and writes its PID to `ignore/linelist/run_grid.pid`. Watch it with:

```bash
tail -f ignore/linelist/run_grid.log
```

## Stage 2: generation-interval sensitivity

Only run this once stage 1 is read and looks sane.

```bash
scripts/linelist/run_grid.sh --stage=2
```

Grid: `onsets × {sitrep, linelist_known} × {cmmid_gi_any, cmmid_gi_case, cmmid_gi_diag}` — six cells, all refitting because of the new cache key. Roughly 5 to 6 hours.

Only `cmmid_gi_*` appears because `onsets` estimates its reporting delay and `delay_config` refuses a report-delay override there rather than ignoring it.

This brackets what the cmmid generation interval does to the answer. It is the sensitivity around stage 1, not the headline.

## Then: build the report

```bash
julia --project=docs scripts/linelist/plot_rt.jl
```

Writes into `ignore/linelist/`:

| | |
|---|---|
| `comparison.csv` | read this first: one row per delay config and quantity, both sources side by side with difference and ratio |
| `diagnostics_all.csv` | R-hat, bulk ESS, divergences per run |
| `rt_onsets.png` | trajectories overlaid |
| `delay_comparison.csv` | each fit's estimated onset-to-report delay against cmmid's independent fit |
| `inputs_at_cutoff.csv` | what each observation set held at the cut-off |
| `rt_all.csv`, `rt_summary.csv`, `stream_estimates_all.csv` | the collected detail |

## What to check before believing anything

1. `diagnostics_all.csv` first. The script warns on stderr for R-hat above 1.01, bulk ESS below 400, or divergences above 1% of draws. The `onsets` runs so far sit at R-hat 1.01–1.03 with ESS 87–330, so expect warnings on ESS: these fits are usable but low-resolution, and the line-list arm has consistently sampled worse than the situation-report arm.
2. The breakpoint line in each run log. `grep -A3 "Info: breakpoint" ignore/linelist/run_onsets_*.log` — every run must print 94.
3. That `comparison.csv` has rows. A configuration run on one source but not the other is dropped rather than half-reported, so a missing row means a missing arm.
4. `plot_rt.jl` errors rather than continues if the series within a fit do not share a last fitted date. That guard exists because a stale file from an earlier cut-off would otherwise be paired against a current one.

## If ESS is still low

These are 500 draws × 2 chains. Raising `samples` would improve resolution but breaks comparability with the release, whose settings the fit registry supplies precisely so a refit stays comparable. If more resolution is wanted, raise it for both arms together and say so when the result is reported.

`--chains=N` exists for diagnosis and puts the chain count in the output tag, so a diagnostic run cannot be mistaken for a comparison arm.
