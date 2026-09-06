#!/usr/bin/env python3
# /// script
# requires-python = ">=3.9"
# dependencies = ["pillow", "numpy"]
# ///
#
# Digitise the "courbe epidemique par date de debut des symptomes (liste
# lineaire DHIS2)" figure that the INSP analytique-format SitReps carry from
# SitRep 059 onward. That figure is the only published source for confirmed
# cases by symptom-onset date; it is a raster bar chart with no accompanying
# data table, so the daily counts have to be recovered from the pixels.
#
# Successive vintages redraw the same onset cohort at later report dates, so
# the set of digitised curves forms a reporting triangle: recent onset dates
# fill in (backfill) as more confirmations arrive. That is the signal used to
# estimate a rough onset-to-report delay (see data/README.md).
#
# Method (per figure, all self-calibrated from the image):
#   * baseline (count 0) = the widest dark horizontal row in the lower panel;
#   * count scale = the y-axis tick marks (0/20/40/60), evenly spaced, giving
#     pixels-per-count = tick-spacing / 20;
#   * date scale = the weekly x-axis tick marks; anchored on the rightmost
#     tick (whose date is given per vintage in CONFIG, read off the axis)
#     stepping back 7 days per tick;
#   * each daily bar height = the 75th-percentile column in a one-day window,
#     flooded up from the baseline counting light-blue (Vivant) and crimson
#     (Decede) pixels, bridging the few-pixel anti-alias gap between the two
#     stacked segments but stopping at the wide white gap up to the floating
#     "premier resultat positif" label / dashed line above the bar.
#
# Accuracy: the error is a few percent in either direction, per scan, and
# it is independent between vintages. Against the printed n it ranges from
# -3.0% (SitRep 069/070/071: 2260 vs n=2 329) to +1.6% (SitRep 068: 2344 vs
# n=2 308), with SitRep 064 at -2.2% (2018 vs n=2 064) and 072 at +0.4%
# (2531 vs n=2 521). Individual daily bars carry roughly +/-1-2 cases of
# pixel noise. Part of the shortfall sits in the faded bars of the `donnees
# potentiellement incompletes` band, whose lightened fill falls outside the
# colour masks, but that mechanism is one-sided and does not explain the
# overshoots, so treat the sign as unknown.
#
# The consequence that matters: the scans do not preserve a property the
# underlying data has. Late reporting only ever adds cases, so an onset
# date's count must be non-decreasing across vintages, yet on onset dates
# more than three weeks before the earliest report date in the file (12
# July, so onsets before 21 June) the scanned totals move both ways between
# consecutive snapshots - 064 -> 065 falls by a net 36 cases across 34 of 54
# such days, and every other consecutive pair falls somewhere too. So a
# between-vintage increment of a few cases is at or
# below the noise floor, and anything built on those increments (a
# reporting-delay estimate, say) has to account for it.
#
# The values are approximate and are not fitted by the model; they are
# captured for later use. See #488.
#
# Dependencies: Pillow and numpy (image analysis) and poppler's pdfimages /
# pdftotext / pdfinfo (figure extraction). The script carries PEP 723 inline
# metadata, so with uv installed the Python deps are fetched automatically
# into a throwaway env:
#   uv run scripts/digitize_onset_curve.py
# Without uv, install into an isolated venv first:
#   python3 -m venv .venv && .venv/bin/pip install Pillow numpy
#   .venv/bin/python scripts/digitize_onset_curve.py
# poppler must be on PATH either way (apt install poppler-utils /
# brew install poppler). See scripts/README.md.
#
# Usage:
#   python3 scripts/digitize_onset_curve.py [pdf_dir] [out_csv]
# Defaults: pdf_dir = data/sitrep_pdfs, out_csv = data/onset_curve_scanned.csv
# Download the PDFs first with scripts/download_sitreps.jl.

import csv
import datetime as dt
import os
import subprocess
import sys
import tempfile

import numpy as np
from PIL import Image

# Per-vintage anchors. `report_date` is the SitRep rapportage date;
# `last_tick` is the date of the rightmost weekly x-axis tick, read off the
# figure (the axis range differs between vintages). To add a new vintage,
# append its SitRep number, rapportage date and last x-axis tick date.
CONFIG = {
    "059": ("2026-07-12", "2026-07-12"),
    "060": ("2026-07-13", "2026-07-12"),
    "061": ("2026-07-14", "2026-07-12"),
    "062": ("2026-07-15", "2026-07-12"),
    "064": ("2026-07-17", "2026-07-15"),
    "065": ("2026-07-18", "2026-07-15"),
    "066": ("2026-07-19", "2026-07-15"),
    "067": ("2026-07-20", "2026-07-15"),
    "068": ("2026-07-21", "2026-07-22"),
    "069": ("2026-07-22", "2026-07-22"),
    "070": ("2026-07-23", "2026-07-22"),
    "071": ("2026-07-24", "2026-07-22"),
    "072": ("2026-07-25", "2026-07-22"),
    "073": ("2026-07-26", "2026-07-22"),
    "074": ("2026-07-27", "2026-07-22"),
    "077": ("2026-07-30", "2026-07-29"),
    "078": ("2026-07-31", "2026-07-29"),
    "079": ("2026-08-01", "2026-07-29"),
    "080": ("2026-08-02", "2026-07-29"),
    "081": ("2026-08-03", "2026-07-29"),
    "082": ("2026-08-04", "2026-08-05"),
    "083": ("2026-08-05", "2026-08-05"),
    "087": ("2026-08-09", "2026-08-05"),
    "088": ("2026-08-10", "2026-08-05"),
    "089": ("2026-08-11", "2026-08-05"),
    "090": ("2026-08-12", "2026-08-05"),
    "091": ("2026-08-13", "2026-08-05"),
    "092": ("2026-08-14", "2026-08-10"),
    "093": ("2026-08-15", "2026-08-10"),
    "094": ("2026-08-16", "2026-08-17"),
    "095": ("2026-08-17", "2026-08-17"),
    "096": ("2026-08-18", "2026-08-17"),
    "097": ("2026-08-19", "2026-08-17"),
    # "098" is deliberately absent. It is the only vintage INSP embedded
    # losslessly rather than as JPEG, so the fixed colour thresholds below
    # keep a fringe of each bar that JPEG blur costs every other vintage,
    # and it reads about 7% high on the same underlying data. Excluding it
    # keeps a vintage on a different bias scale out of the between-vintage
    # increments this file feeds. The evidence, and the controls that rule
    # out the render size, are in data/README.md. Read them before adding
    # it back.
    "099": ("2026-08-21", "2026-08-17"),
    "100": ("2026-08-22", "2026-08-17"),
    "101": ("2026-08-23", "2026-08-24"),
    "102": ("2026-08-24", "2026-08-24"),
    "103": ("2026-08-25", "2026-08-24"),
    "104": ("2026-08-26", "2026-08-24"),
    "105": ("2026-08-27", "2026-08-24"),
    "106": ("2026-08-28", "2026-08-24"),
    "107": ("2026-08-29", "2026-08-24"),
    "108": ("2026-08-30", "2026-08-31"),
    "109": ("2026-08-31", "2026-08-31"),
    # "110" is deliberately absent: its epidemic-curve figure is plotted by
    # notification date, not symptom-onset date (see data/README.md and
    # issue #644), so it is not a snapshot of this stream at all.
    "111": ("2026-09-02", "2026-08-31"),
    "112": ("2026-09-03", "2026-08-31"),
    "113": ("2026-09-04", "2026-08-31"),
}

# Every figure through SitRep 083 draws its y-axis on a 0/20/40/60/80 grid,
# which `digitize` assumed as a hard-coded divisor. From SitRep 087 the
# brief-format figure switched to a 0/25/50/75 grid (confirmed by reading
# the printed tick labels directly - the pixel geometry is otherwise
# indistinguishable, so this cannot be self-calibrated any more than
# `last_tick` can). Applying the old /20 divisor to a 25-count grid
# undercounts every bar by a scale-dependent amount and was caught only
# because it made stable, weeks-old onset dates fall (SitRep 083's 15 May
# read 26; the same date misread through the old divisor came out as 8).
# Override per vintage here; anything absent keeps the historical 20.
Y_AXIS_STEP = {
    "087": 25,
    "088": 25,
    "089": 25,
    "090": 25,
    "091": 25,
    "092": 25,
    "093": 25,
    "094": 25,
    "095": 25,
    "096": 25,
    "097": 25,
    "099": 25,
    "100": 25,
    "101": 25,
    "102": 25,
    "103": 25,
    "104": 25,
    "105": 25,
    "106": 25,
    "107": 25,
    "108": 25,
    "109": 25,
    "111": 25,
    "112": 25,
    "113": 25,
}


def _masks(im):
    R, G, B = im[:, :, 0], im[:, :, 1], im[:, :, 2]
    blue = (B > 150) & (G > 150) & (R < 210) & (B >= R + 15)
    red = (R > 120) & (R >= G + 50) & (R >= B + 50)
    orange = (R > 200) & (G > 110) & (G < 195) & (B < 90)
    dark = (R < 120) & (G < 120) & (B < 120)
    return blue, red, orange, dark


def _is_onset_curve(im):
    # The onset figure is a blue-dominant daily bar chart with no orange
    # (the age/sex pyramids use orange; the notification-week chart uses a
    # darker steel blue and prints value labels). The tests are pixel
    # fractions, not counts, because INSP re-renders the figure at whatever
    # size the layout needs and an absolute threshold silently flips as the
    # size moves: the blue floor already had to be lowered once when the
    # figure shrank to 1009x583, and SitRep 072's larger 1277x799
    # rendering then pushed the crimson/pink-band anti-aliasing to 631
    # orange pixels, past a 500 cut.
    #
    # Measured over SitReps 059-080, on the caption page and its immediate
    # neighbours (the page-fallback in extract_onset_image widens the search
    # there, which brings the provincial case map into the candidate pool -
    # it is blue-heavy too, from the lake/river fill and legend swatches):
    #   blue fraction    onset 0.066-0.125   province map 0.045-0.047
    #   orange fraction  onset <= 0.0007     age/sex pyramids >= 0.053
    #   red fraction     onset >= 0.046      (a floor, not a discriminator:
    #                                        the notification-week chart is
    #                                        also crimson-heavy)
    # The map's blue fraction sits clear below every onset chart seen so far,
    # so 0.055 (roughly the midpoint of the two clusters) discriminates with
    # margin on both sides without needing a caption-text match.
    blue, red, orange, _ = _masks(im)
    npx = im.shape[0] * im.shape[1]
    return (blue.sum() / npx > 0.055 and orange.sum() / npx < 0.01
            and red.sum() / npx > 0.01)


def _onset_page(pdf):
    # The onset figure usually sits on the page whose text carries its
    # caption.
    npages = int(subprocess.run(
        ["pdfinfo", pdf], check=True, capture_output=True, text=True
    ).stdout.split("Pages:")[1].split()[0])
    for p in range(1, npages + 1):
        txt = subprocess.run(
            ["pdftotext", "-layout", "-f", str(p), "-l", str(p), pdf, "-"],
            check=True, capture_output=True, text=True,
        ).stdout.lower()
        if ("date de debut des symptom" in txt
                or "date de début des symptôm" in txt):
            return p, npages
    return None, npages


def _best_onset_image(pdf, page, workdir):
    subprocess.run(
        ["pdfimages", "-png", "-f", str(page), "-l", str(page), pdf,
         os.path.join(workdir, "p")],
        check=True, capture_output=True,
    )
    best = None
    for name in sorted(os.listdir(workdir)):
        if not name.endswith(".png"):
            continue
        im = np.asarray(Image.open(os.path.join(workdir, name)).convert("RGB"))
        im = im.astype(int)
        if _is_onset_curve(im) and (best is None or im.size > best.size):
            best = im
        os.remove(os.path.join(workdir, name))
    return best


def extract_onset_image(pdf, workdir):
    p, npages = _onset_page(pdf)
    if p is None:
        return None
    best = _best_onset_image(pdf, p, workdir)
    if best is not None:
        return best
    # SitRep 080 embeds the chart on page 5 under a mislabelled caption ("par
    # semaine de notification") while the matching "date de debut des
    # symptomes" caption text sits on page 6 with no image of its own, so the
    # caption-text page lookup lands one page short of the real figure. Widen
    # to the immediate neighbours only (not the whole document): the map and
    # other embedded figures elsewhere in the report are large enough, and
    # blue enough in places (lakes, legends), to satisfy _is_onset_curve too,
    # so a document-wide scan silently grabs the wrong image.
    for q in (p - 1, p + 1):
        if 1 <= q <= npages:
            best = _best_onset_image(pdf, q, workdir)
            if best is not None:
                return best
    return None


def _longest_run(colmask):
    best = cur = 0
    for v in colmask:
        cur = cur + 1 if v else 0
        best = max(best, cur)
    return best


def _cluster(idx, gap=3):
    out, cl = [], []
    for i in idx:
        if cl and i - cl[-1] <= gap:
            cl.append(i)
        else:
            if cl:
                out.append(int(np.mean(cl)))
            cl = [i]
    if cl:
        out.append(int(np.mean(cl)))
    return out


def _baseline_row(im, H):
    # The count-0 baseline is the plot's bottom border: a solid line running
    # almost the full chart width. Score rows by their longest contiguous run
    # under a near-gray threshold (<180); a run-length ranking under that
    # threshold correctly finds the border in every vintage, including
    # tighter-anti-aliased renders, unlike a per-row pixel sum.
    R, G, B = im[:, :, 0], im[:, :, 1], im[:, :, 2]
    line = (R < 180) & (G < 180) & (B < 180)
    line[: int(H * 0.4)] = False
    best_row, best_run = 0, 0
    for r in range(H):
        run = _longest_run(line[r])
        if run > best_run:
            best_run, best_row = run, r
    if best_run < 100:
        raise ValueError("no baseline row found")
    return best_row


def _y_axis_ticks(dark, base, W):
    # The 0/20/40/60 y-axis tick rows, read from the label strip just left
    # of the vertical axis line. Candidate strips are scored by the longest
    # dark vertical run (the axis line itself), but only among strips whose
    # rows form a plausible axis: at least three clusters, evenly spaced,
    # with the last one (the 0 tick) on the baseline. Taking the longest run
    # alone is not enough - in SitRep 067 a glyph stroke outruns the real
    # axis line and yields a scale that halves every count.
    # Slice bounds are the 0-based images of the reference's 1-based ranges:
    # candidate columns 30..floor(W*0.13), the label strip being the ten
    # columns immediately left of the candidate, and rows running down to
    # three past the baseline.
    best = None
    for x in range(29, int(W * 0.13)):
        seg = dark[: base + 4, max(0, x - 10):x].sum(axis=1)
        yt = _cluster([y for y in range(len(seg)) if seg[y] >= 3])
        if len(yt) < 3 or abs(yt[-1] - base) > 3:
            continue
        d = np.diff(yt)
        if d.min() <= 5 or d.max() > 1.15 * d.min():
            continue
        rank = (_longest_run(dark[: base + 1, x]), -x)  # tie-break leftmost
        if best is None or rank > best[0]:
            best = (rank, yt)
    return None if best is None else best[1]


def digitize(im, last_tick_date, y_step=20):
    H, W, _ = im.shape
    blue, red, _, dark = _masks(im)
    base = _baseline_row(im, H)  # count-0 baseline row
    # count scale from the y-axis ticks (0/20/40/60 through SitRep 083;
    # 0/25/50/75 from SitRep 087 - see Y_AXIS_STEP)
    yt = _y_axis_ticks(dark, base, W)
    if yt is None:
        # Only fall back to the <180 near-gray mask when the strict <120
        # mask finds nothing, so every already-committed vintage (059-111)
        # keeps digitising under the original threshold. SitRep 112 shrinks
        # the embedded render's height (771x433 against 775-802x475-479 for
        # its neighbours), anti-aliasing the tick-label strip past the
        # strict mask; the weekly x-axis ticks already carry this same
        # scoped fallback (added for SitRep 108's larger render).
        R, G, B = im[:, :, 0], im[:, :, 1], im[:, :, 2]
        line = (R < 180) & (G < 180) & (B < 180)
        yt = _y_axis_ticks(line, base, W)
    if yt is None:
        raise ValueError("no y-axis tick strip found")
    ppc = np.median(np.diff(yt)) / float(y_step)
    ytop, y0 = yt[0], yt[-1]
    # x scale from the weekly tick marks below the baseline. The tick marks
    # are only a few pixels tall and shrink with the embedded figure
    # resolution (5-6 dark rows in the 1257x698 SitRep 064 rendering, 4 in
    # SitRep 066's 1275x623, 3 in SitRep 069/070's 1009x583), so step the
    # cut down until a full weekly row of ticks resolves instead of fixing
    # it at 4 and losing the axis entirely on the smaller figures.
    # They sit just below the baseline (a few px) and, on the faint
    # JPEG-compressed figures (SitRep 081), can be only 1px tall, so cut must
    # come all the way down to 1 to resolve them; the window stops at base+7
    # so a wide low-cut scan cannot pick up the x-axis date labels further
    # down. Step down through the cuts and keep the most complete weekly tick
    # row (the true axis has a fixed number of weekly ticks, so a too-strict
    # cut silently drops every other tick rather than failing).
    def _weekly_ticks(mask):
        band = mask[base + 2:base + 7, :].sum(axis=0)
        best_n, best = 0, np.array([])
        for cut in (4, 3, 2, 1):
            cand = np.array(_cluster([x for x in range(W) if band[x] >= cut]))
            if len(cand) >= 8 and len(cand) > best_n:
                best_n = len(cand)
                best = cand
        return best

    best = _weekly_ticks(dark)
    if len(best) == 0:
        # Only fall back to the <180 near-gray mask when the strict mask
        # finds nothing, so every already-committed vintage (059-107) keeps
        # digitising under the original threshold.
        R, G, B = im[:, :, 0], im[:, :, 1], im[:, :, 2]
        line = (R < 180) & (G < 180) & (B < 180)
        best = _weekly_ticks(line)
    if len(best) == 0:
        raise ValueError("no x-axis weekly tick row found")
    xt = best
    ppd = np.median(np.diff(xt)) / 7.0  # pixels per day
    # The daily bar windows are laid out in the 1-based pixel frame the
    # Julia reference uses, and converted back to 0-based only at the point
    # of indexing. Both languages round half to even, but that rule is not
    # translation-invariant: a window edge landing exactly on .5 rounds down
    # in one frame and up in the other, so a 0-based layout reads windows
    # one column apart from the reference. A dropped column changes which
    # bar the 75th percentile lands on, so the error is a whole segment, not
    # a count or two.
    lastx = xt[-1] + 1                   # rightmost tick is always real
    lastdate = dt.date.fromisoformat(last_tick_date)
    # per-column stacked bar height, flooded up from the baseline
    bc = np.zeros(W)
    rc = np.zeros(W)
    for x in range(W):
        r, miss, b, rr = y0 - 1, 0, 0, 0
        while r > ytop - 2 and r >= 0:
            if blue[r, x]:
                b += 1
                miss = 0
            elif red[r, x]:
                rr += 1
                miss = 0
            else:
                miss += 1
                if miss > 6:
                    break
            r -= 1
        bc[x], rc[x] = b, rr
    tot = bc + rc
    nz = np.where(tot > 2)[0] + 1        # 1-based, matching lastx
    barmin, barmax = nz.min(), nz.max()
    rows = []
    for off in range(-105, 4):
        cx = lastx + off * ppd
        if cx < barmin - ppd or cx > barmax + ppd:
            continue
        lo, hi = int(round(cx - ppd * 0.45)), int(round(cx + ppd * 0.45))
        cols = range(max(1, lo) - 1, min(W, hi))
        bvals = [bc[c] for c in cols]
        rvals = [rc[c] for c in cols]
        if max(b + r for b, r in zip(bvals, rvals)) < 1:
            continue
        alive = round(float(np.percentile(bvals, 75)) / ppc)
        dead = round(float(np.percentile(rvals, 75)) / ppc)
        date = lastdate + dt.timedelta(days=off)
        rows.append((date.isoformat(), alive, dead))
    # drop trailing zero rows and isolated tiny strays past the curve tail
    while rows and (rows[-1][1] + rows[-1][2]) == 0:
        rows.pop()
    while len(rows) >= 2:
        gap = (dt.date.fromisoformat(rows[-1][0])
               - dt.date.fromisoformat(rows[-2][0])).days
        if gap > 1 and (rows[-1][1] + rows[-1][2]) <= 2:
            rows.pop()
        else:
            break
    return rows


def main():
    pdf_dir = sys.argv[1] if len(sys.argv) > 1 else "data/sitrep_pdfs"
    out_csv = (sys.argv[2] if len(sys.argv) > 2
               else "data/onset_curve_scanned.csv")
    out_rows = []
    for sr in sorted(CONFIG):
        report_date, last_tick = CONFIG[sr]
        pdf = os.path.join(pdf_dir, f"SitRep_MVE_{sr}_2026.pdf")
        if not os.path.isfile(pdf):
            print(f"skip {sr}: {pdf} not found", file=sys.stderr)
            continue
        with tempfile.TemporaryDirectory() as wd:
            im = extract_onset_image(pdf, wd)
        if im is None:
            print(f"skip {sr}: no onset curve found", file=sys.stderr)
            continue
        rows = digitize(im, last_tick, Y_AXIS_STEP.get(sr, 20))
        total = sum(a + d for _, a, d in rows)
        print(f"SitRep {sr} ({report_date}): {len(rows)} onset days, "
              f"total {total} confirmed")
        for onset, alive, dead in rows:
            out_rows.append((sr, report_date, onset, alive, dead, alive + dead))
    with open(out_csv, "w", newline="") as f:
        w = csv.writer(f, lineterminator="\n")  # LF, matching the Julia ref
        w.writerow(["sitrep", "report_date", "onset_date",
                    "confirmed_alive", "confirmed_dead", "confirmed_total"])
        w.writerows(out_rows)
    print(f"wrote {len(out_rows)} rows to {out_csv}")


if __name__ == "__main__":
    main()
