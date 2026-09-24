#!/usr/bin/env python3
# /// script
# requires-python = ">=3.9"
# dependencies = ["pillow", "numpy"]
# ///
#
# Digitise the "courbe epidemique par date de debut des symptomes (liste
# lineaire DHIS2)" figure that the INSP analytique-format SitReps carry from
# SitRep 059 onward. That figure is the only published source for confirmed
# cases by symptom-onset date; it is a raster bar chart with no data table,
# so the daily counts are recovered from the figure pixels.
#
# This is the Python port of scripts/digitize_onset_curve.jl, for the
# automated data-updater, which has Python (not Julia) access. The Julia
# script is the reference and data/onset_curve_scanned.csv is its output;
# this port must reproduce that file byte for byte, and
# test/test_onset_digitiser.jl checks that it does whenever the PDFs are
# present. Every function below is the reference's function of the same
# name, with the same pixel classes, thresholds, tie-breaks and rounding.
# The reference indexes pixels from 1, so the bar-window arithmetic is done
# in that frame and converted to 0-based only at the point of indexing:
# round-half-to-even is not translation-invariant, and a window edge on
# exactly .5 would otherwise land one column off.
#
# Method (per figure, all self-calibrated from the image):
#   * baseline (count 0) = the widest dark horizontal row in the lower panel;
#   * count scale = the y-axis tick marks (0/20/40/60 or 0/25/50/75), evenly
#     spaced, giving pixels-per-count = tick-spacing / y_step;
#   * date scale = the weekly x-axis tick marks. Candidate tick rows come
#     from a strict and a near-gray mask at several cuts, and the one whose
#     regular chain from the rightmost tick is longest wins. Pixels per day
#     is the least-squares slope over that chain and each day is anchored on
#     the nearest chain tick at or before it. The rightmost tick's date is
#     in CONFIG, read off the axis;
#   * each daily bar = the bar's own pixel columns, taken as the window
#     [cx - ppd/2, cx + ppd/2] clipped to the nearest outline column on each
#     side. Outline columns are those mostly dark over their run. Every
#     column is read as the run of non-page pixels up from the baseline,
#     bridging up to three page pixels when bar colour resumes and skipping
#     neutral pixels on the tick rows and tick columns (gridlines). The
#     run's top is its highest pixel darker than an anti-alias, which is the
#     bar's outline. The bar height is the height most of its interior
#     columns agree on, or the tallest interior column when fewer than three
#     agree. Half a pixel of outline is subtracted before dividing by
#     pixels-per-count. The dead segment is the count of crimson pixels in
#     the chosen column.
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
# Incremental by default. A run reuses the rows out_csv already carries and
# opens the PDF only for the CONFIG vintages missing from it. The rows are
# written back in CONFIG order either way, so an incremental run and a full
# one produce the same file. A change to the digitiser itself does not
# invalidate those reused rows, so re-run with --rebuild after touching the
# digitising code, which re-reads every vintage. Each run prints how many
# vintages it reused and how many it read.
#
# Usage:
#   python3 scripts/digitize_onset_curve.py [pdf_dir] [out_csv] [--rebuild]
# Defaults: pdf_dir = data/sitrep_pdfs, out_csv = data/onset_curve_scanned.csv
# Download the PDFs first with scripts/download_sitreps.jl.

import csv
import datetime as dt
import math
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
    # "110" is deliberately absent. Its page-4 figure carries the same
    # outer caption as every other vintage ("par date de début des
    # symptômes") but the embedded chart's own internal title and x-axis
    # read "par date de NOTIFICATION" (n = 5 710) - a genuine basis change,
    # confirmed by extracting and viewing the raw embedded image rather
    # than trusting the caption. Digitising it would silently inject a
    # different-basis series into the reporting-triangle stream. See
    # data/README.md and issue #644.
    "111": ("2026-09-02", "2026-08-31"),
    "112": ("2026-09-03", "2026-08-31"),
    "113": ("2026-09-04", "2026-08-31"),
    "114": ("2026-09-05", "2026-08-31"),
    "115": ("2026-09-06", "2026-09-07"),
    "116": ("2026-09-07", "2026-09-07"),
    # 117 and 118 keep 116's 07 September tick and plot a window that
    # retreats behind their own report dates (to 09-03 and 09-04). That
    # only shortens the pair's coverage intersection, which
    # `load_onset_curve` drops as unobserved.
    "117": ("2026-09-08", "2026-09-07"),
    "118": ("2026-09-09", "2026-09-07"),
    "119": ("2026-09-10", "2026-09-07"),
    "120": ("2026-09-11", "2026-09-07"),
    "121": ("2026-09-12", "2026-09-07"),
    "122": ("2026-09-13", "2026-09-14"),
    "123": ("2026-09-14", "2026-09-14"),
    "124": ("2026-09-15", "2026-09-14"),
    "125": ("2026-09-16", "2026-09-14"),
    "126": ("2026-09-17", "2026-09-14"),
    "127": ("2026-09-18", "2026-09-14"),
    "128": ("2026-09-19", "2026-09-14"),
    "129": ("2026-09-20", "2026-09-21"),
    "130": ("2026-09-21", "2026-09-21"),
}

# Every figure through SitRep 083 draws its y-axis on a 0/20/40/60/80 grid.
# From SitRep 087 the brief-format figure switched to a 0/25/50/75 grid
# (confirmed by reading the printed tick labels directly - the pixel
# geometry is otherwise indistinguishable, so this cannot be
# self-calibrated any more than `last_tick` can). Override per vintage
# here; anything absent keeps the historical 20.
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
    "114": 25,
    "115": 25,
    "116": 25,
    "117": 25,
    "118": 25,
    "119": 25,
    "120": 25,
    "121": 25,
    "122": 25,
    "123": 25,
    "124": 25,
    "125": 25,
    "126": 25,
    "127": 25,
    "128": 25,
    "129": 25,
    "130": 25,
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
    # size the layout needs. Measured over SitReps 059-080, on the caption
    # page and its immediate neighbours (which brings the provincial case
    # map into the candidate pool):
    #   blue fraction    onset 0.066-0.125   province map 0.045-0.047
    #   orange fraction  onset <= 0.0007     age/sex pyramids >= 0.053
    #   red fraction     onset >= 0.046      (a floor, not a discriminator)
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
    # pdfimages writes PPM (P6) for RGB images by default, which is what the
    # reference parses; Pillow reads the same file, so both see the same
    # bytes.
    subprocess.run(
        ["pdfimages", "-f", str(page), "-l", str(page), pdf,
         os.path.join(workdir, "p")],
        check=True, capture_output=True,
    )
    best = None
    for name in sorted(os.listdir(workdir)):
        if not name.endswith(".ppm"):
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
    # Cluster nearly-adjacent indices, returning the floored mean of each
    # cluster. The floor commutes with the 1-based shift of the reference.
    out, cl = [], []
    for i in idx:
        if cl and i - cl[-1] <= gap:
            cl.append(i)
        else:
            if cl:
                out.append(sum(cl) // len(cl))
            cl = [i]
    if cl:
        out.append(sum(cl) // len(cl))
    return out


def _baseline_row(im, H):
    # The count-0 baseline is the plot's bottom border: a solid line running
    # almost the full chart width. Score rows by their longest contiguous run
    # under a near-gray threshold (<180), searched over the lower 60% of the
    # image; the first row with the longest run wins.
    R, G, B = im[:, :, 0], im[:, :, 1], im[:, :, 2]
    line = (R < 180) & (G < 180) & (B < 180)
    best_row, best_run = int(H * 0.4), 0
    for r in range(int(H * 0.4), H):
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
    # with the last one (the 0 tick) on the baseline. Slice bounds are the
    # 0-based images of the reference's 1-based ranges: candidate columns
    # 30..floor(W*0.13), the label strip being the ten columns immediately
    # left of the candidate, and rows running down to three past the
    # baseline.
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
    if best is None:
        raise ValueError("no y-axis tick strip found")
    return best[1]


def _pixel_classes(im):
    # Page is white (with JPEG chroma noise) and the pink `donnees
    # potentiellement incompletes` band. Neutral is gridline gray. Light is
    # the pale, low-saturation pixel a bar's top edge leaves above its
    # outline; saturated fill and the outline are not light. Dark is the
    # outline.
    R, G, B = im[:, :, 0], im[:, :, 1], im[:, :, 2]
    lo = np.minimum(R, np.minimum(G, B))
    hi = np.maximum(R, np.maximum(G, B))
    spread = hi - lo
    page = (((lo >= 228) & (spread <= 25))
            | ((R >= 238) & (G >= 200) & (B >= 200) & (R - G >= 15)))
    neutral = (lo >= 170) & (spread <= 12)
    light = (hi >= 190) & (spread < 60)
    crimson = (R - np.maximum(G, B)) >= 25
    darkpx = (R < 150) & (G < 150) & (B < 150)
    return page, neutral, light, crimson, darkpx


def _column_runs(page, neutral, light, crimson, darkpx, y0, gridrows,
                 gridcols, gap=3):
    # Per-column run of non-page pixels up from the baseline. Gridlines lie
    # on the tick rows and tick columns, so a neutral pixel there is page.
    # Up to `gap` page pixels are bridged when a non-page pixel follows. The
    # run's top is the highest pixel darker than `light`. Returns the run
    # height, the crimson count and the dark count per column, each counted
    # over the non-page pixels from the baseline up to that top.
    #
    # The reference walks each column with a miss counter; here every
    # column is walked at once. Upward index i is image row y0 - 1 - i.
    H, W = page.shape
    grid = np.zeros((H, W), dtype=bool)
    for r in gridrows:
        grid[max(0, r - 1):r + 2, :] = True
    for c in gridcols:
        grid[:, max(0, c - 1):c + 2] = True
    ok = ~page & ~(grid & neutral)
    up = ok[:y0][::-1]
    n = up.shape[0]
    idx = np.arange(n)[:, None]
    # consecutive misses ending at each upward index; the walk stops at the
    # first index whose miss run exceeds `gap`, and reads nothing above it
    last_ok = np.maximum.accumulate(np.where(up, idx, -1), axis=0)
    brk = (idx - last_ok) > gap
    stop = np.where(brk.any(axis=0), brk.argmax(axis=0), n)
    seen = up & (idx < stop)
    cand = seen & ~light[:y0][::-1]
    top = np.where(cand.any(axis=0), n - 1 - cand[::-1].argmax(axis=0), -1)
    counted = seen & (idx <= top)
    h = top + 1
    nr = (counted & crimson[:y0][::-1]).sum(axis=0)
    nd = (counted & darkpx[:y0][::-1]).sum(axis=0)
    return h, nr, nd


def _tick_chain(xt):
    # The regular weekly chain ending on the rightmost tick, as (week index,
    # x) pairs. Walking left, a spacing of one or two weeks within 8% (at
    # least 2.5 px) of the median spacing continues the chain; anything else
    # ends it.
    s = float(np.median(np.diff(xt)))
    ks = [0]
    xs = [xt[-1]]
    for j in range(len(xt) - 2, -1, -1):
        d = xt[j + 1] - xt[j]
        k = round(d / s)
        if not (1 <= k <= 2 and abs(d - k * s) <= max(2.5, 0.08 * s)):
            break
        ks.insert(0, ks[0] - k)
        xs.insert(0, xt[j])
    return ks, xs


def _modal_height(h, cols, cx):
    # Most common value of `h` over `cols` and how many columns carry it;
    # ties go to the value nearest `cx` by column, then to the value seen
    # first. `h` and `cols` share one index frame.
    best = 0
    bestkey = (-1, -math.inf)
    for u in dict.fromkeys(h[x] for x in cols):
        c = sum(1 for x in cols if h[x] == u)
        d = min(abs(x - cx) for x in cols if h[x] == u)
        key = (c, -d)
        if key > bestkey:
            best = u
            bestkey = key
    return best, bestkey[0]


def calibrate(im, y_step=20):
    """Axis calibration of one figure: the count-0 baseline row `base`, the
    y-axis tick rows `yt`, pixels per count `ppc`, the weekly tick chain
    `ks` (week index) and `xs` (1-based column), and pixels per day `ppd`.
    All rows and columns are 0-based except `xs`."""
    H, W, _ = im.shape
    _, _, _, dark = _masks(im)
    base = _baseline_row(im, H)  # count-0 baseline row
    R, G, B = im[:, :, 0], im[:, :, 1], im[:, :, 2]
    line = (R < 180) & (G < 180) & (B < 180)
    # count scale from the y-axis ticks (0/20/40/60 through SitRep 083;
    # 0/25/50/75 from SitRep 087 - see Y_AXIS_STEP)
    try:
        yt = _y_axis_ticks(dark, base, W)
    except ValueError:
        # SitRep 112's smaller render (771x433) anti-aliases the tick marks
        # and the axis line into the 120-180 near-gray range, so the strict
        # <120 mask finds three of the four ticks but not the one sitting on
        # the baseline itself. Only tried when the strict mask finds nothing.
        yt = _y_axis_ticks(line, base, W)
    ppc = float(np.median(np.diff(yt))) / float(y_step)
    # x scale from the weekly tick marks 2-6 rows below the baseline. Both
    # masks are tried at every cut and the tick row whose regular weekly
    # chain from the rightmost tick is longest wins.
    best_n = 0
    xt = []
    for mask in (dark, line):
        band = mask[base + 2:base + 7, :].sum(axis=0)
        for cut in (4, 3, 2, 1):
            cand = _cluster([x for x in range(W) if band[x] >= cut])
            if len(cand) < 8:
                continue
            n = len(_tick_chain(cand)[0])
            if n > best_n:
                best_n = n
                xt = cand
    if not xt:
        raise ValueError("no x-axis weekly tick row found")
    ks, xs0 = _tick_chain(xt)
    n = len(xs0)
    if n < 2:
        raise ValueError("weekly tick chain too short")
    # pixels per day by least squares over the chain; every term is an
    # integer-valued float, so the sums are exact in any order
    xs = [x + 1 for x in xs0]  # 1-based, the frame the bar windows use
    w = [7.0 * k for k in ks]
    swx = sum(a * b for a, b in zip(w, xs))
    sw = sum(w)
    sx = sum(xs)
    sww = sum(a * a for a in w)
    ppd = (n * swx - sw * sx) / (n * sww - sw * sw)
    return {"H": H, "W": W, "base": base, "yt": yt, "ppc": ppc,
            "y0": yt[-1], "ks": ks, "xs": xs, "ppd": ppd}


def day_column(cal, off):
    """1-based x of the day `off` days from the rightmost tick, anchored on
    the nearest chain tick at or before it."""
    ks, xs = cal["ks"], cal["xs"]
    js = [i for i, k in enumerate(ks) if 7 * k <= off]
    j = js[-1] if js else 0
    return xs[j] + (off - 7 * ks[j]) * cal["ppd"]


def digitize(im, last_tick_date, y_step=20):
    cal = calibrate(im, y_step)
    W, yt, ppc, y0 = cal["W"], cal["yt"], cal["ppc"], cal["y0"]
    ks, xs, ppd = cal["ks"], cal["xs"], cal["ppd"]
    lastdate = dt.date.fromisoformat(last_tick_date)
    page, neutral, light, crimson, darkpx = _pixel_classes(im)
    h0, nr0, nd0 = _column_runs(
        page, neutral, light, crimson, darkpx, y0, yt, [x - 1 for x in xs]
    )
    # The bar windows are laid out in the reference's 1-based column frame,
    # so pad each per-column array with a leading dummy and index it with
    # the 1-based column directly.
    pad = np.zeros(1, dtype=h0.dtype)
    h = np.concatenate([pad, h0])
    nr = np.concatenate([pad, nr0])
    nd = np.concatenate([pad, nd0])
    # outline columns are mostly dark over their run (a one-count bar is
    # all outline, so the floor keeps it as interior); an outline drawn
    # across two columns leaves a softer second column that still carries
    # the neighbour's height, dropped when anything else is left
    isborder = (h > 4) & (nd >= 0.3 * h)
    soft = (h > 4) & (nd >= np.maximum(0.1 * h, 3))
    nz = np.flatnonzero((h > 2) & ~isborder)
    barmin, barmax = int(nz.min()), int(nz.max())
    rows = []
    for off in range(7 * ks[0] - 7, 4):
        # anchor on the nearest chain tick at or before the day
        cx = day_column(cal, off)
        if cx < barmin - ppd or cx > barmax + ppd:
            continue
        lo = max(1, math.ceil(cx - ppd / 2 + 0.5))
        hi = min(W, math.floor(cx + ppd / 2 - 0.5))
        # clip the window to the nearest outline column on each side: the
        # window starts after the left outline and ends on the right one,
        # as in the reference
        c = round(cx)
        reach = math.ceil(ppd)
        left = [x for x in range(max(1, c - reach), c) if isborder[x]]
        if left:
            lo = max(lo, left[-1] + 1)
        right = [x for x in range(c + 1, min(W, c + reach) + 1)
                 if isborder[x]]
        if right:
            hi = min(hi, right[0])
        cols = [x for x in range(lo, hi + 1) if not soft[x]]
        if not cols:
            cols = [x for x in range(lo, hi + 1) if not isborder[x]]
        if not cols:
            continue
        hb, support = _modal_height(h, cols, cx)
        if support < 3:
            hb = max(int(h[x]) for x in cols)
        if hb < 1:
            continue
        jb = next(x for x in cols if h[x] == hb)
        total = round(max(0.0, hb - 0.5) / ppc)
        dead = min(total, round(max(0.0, float(nr[jb]) - 0.5) / ppc))
        rows.append((lastdate + dt.timedelta(days=off), total - dead, dead))
    # drop trailing zero rows and isolated tiny strays past the curve tail
    while rows and (rows[-1][1] + rows[-1][2]) == 0:
        rows.pop()
    while len(rows) >= 2:
        gap = (rows[-1][0] - rows[-2][0]).days
        if gap > 1 and (rows[-1][1] + rows[-1][2]) <= 2:
            rows.pop()
        else:
            break
    return rows


HEADER = ["sitrep", "report_date", "onset_date",
          "confirmed_alive", "confirmed_dead", "confirmed_total"]


def digitised_rows(out_csv):
    """Rows out_csv already holds, keyed by SitRep, so a run can reuse a
    vintage it has already read rather than open the PDF again. An absent,
    empty or differently-headed file yields nothing, and every vintage is
    then read, which is what a first run does anyway."""
    rows = {}
    if not os.path.isfile(out_csv):
        return rows
    with open(out_csv, newline="") as f:
        r = csv.reader(f)
        try:
            if next(r) != HEADER:
                return {}
        except StopIteration:
            return {}
        for row in r:
            if row:
                rows.setdefault(row[0], []).append(tuple(row))
    return rows


def main():
    args = [a for a in sys.argv[1:] if a != "--rebuild"]
    rebuild = "--rebuild" in sys.argv[1:]
    pdf_dir = args[0] if len(args) > 0 else "data/sitrep_pdfs"
    out_csv = args[1] if len(args) > 1 else "data/onset_curve_scanned.csv"
    cached = {} if rebuild else digitised_rows(out_csv)
    reused = 0
    read_now = 0
    out_rows = []
    for sr, (report_date, last_tick) in CONFIG.items():
        # Already digitised, so its rows are carried through untouched. They
        # are written in CONFIG order like any other, so reusing them cannot
        # reorder the file.
        if sr in cached:
            out_rows.extend(cached[sr])
            reused += 1
            continue
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
        # An onset date can never sit later than the axis of the report
        # that draws it, and that axis runs at most a day past the
        # rapportage date (the date-de-publication lag). The window above
        # self-calibrates from pixel content and can read a few stray days
        # past the last labelled tick when the figure's own "donnees
        # potentiellement incompletes" band extends that far (SitRep 115);
        # drop those here rather than loosen the invariant
        # test/test_onset_digitiser.jl checks.
        cutoff = dt.date.fromisoformat(report_date) + dt.timedelta(days=1)
        rows = [r for r in rows if r[0] <= cutoff]
        total = sum(a + d for _, a, d in rows)
        print(f"SitRep {sr} ({report_date}): {len(rows)} onset days, "
              f"total {total} confirmed")
        for onset, alive, dead in rows:
            out_rows.append(
                (sr, report_date, onset.isoformat(), alive, dead,
                 alive + dead)
            )
        read_now += 1
    with open(out_csv, "w", newline="") as f:
        w = csv.writer(f, lineterminator="\n")  # LF, matching the Julia ref
        w.writerow(HEADER)
        w.writerows(out_rows)
    print(f"wrote {out_csv}: {read_now} vintages read, "
          f"{reused} reused from the existing file")
    if reused:
        print("re-run with --rebuild to re-read every vintage after "
              "changing the digitiser")


if __name__ == "__main__":
    main()
