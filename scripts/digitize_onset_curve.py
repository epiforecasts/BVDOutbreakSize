#!/usr/bin/env python3
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
#   * date scale = the weekly x-axis tick marks; anchored on the RIGHTMOST
#     tick (whose date is given per vintage in CONFIG, read off the axis)
#     stepping back 7 days per tick;
#   * each daily bar height = the 75th-percentile column in a one-day window,
#     flooded up from the baseline counting light-blue (Vivant) and crimson
#     (Decede) pixels, bridging the few-pixel anti-alias gap between the two
#     stacked segments but stopping at the wide white gap up to the floating
#     "premier resultat positif" label / dashed line above the bar.
#
# Accuracy: the digitised per-vintage totals run ~2% below the printed n
# (SitRep 064: 2019 vs printed n=2 064), and individual daily bars carry
# roughly +/-1-2 cases of pixel noise. The values are approximate and are
# NOT fitted by the model; they are captured for later use.
#
# Dependencies: Pillow and numpy (image analysis) and poppler's pdfimages
# (figure extraction). Create an isolated env if the system Python lacks
# them, e.g.:
#   python3 -m venv .venv && .venv/bin/pip install Pillow numpy
#   .venv/bin/python scripts/digitize_onset_curve.py
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
    # darker steel blue and prints value labels).
    blue, red, orange, _ = _masks(im)
    return blue.sum() > 60000 and orange.sum() < 500 and red.sum() > 20000


def _onset_page(pdf):
    # The onset figure sits on the page whose text carries its caption.
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
            return p
    return None


def extract_onset_image(pdf, workdir):
    p = _onset_page(pdf)
    if p is None:
        return None
    subprocess.run(
        ["pdfimages", "-png", "-f", str(p), "-l", str(p), pdf,
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
    return best


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


def digitize(im, last_tick_date):
    H, W, _ = im.shape
    blue, red, _, dark = _masks(im)
    drow = dark.sum(axis=1)
    drow[: int(H * 0.4)] = 0
    base = int(np.argmax(drow))  # count-0 baseline row
    # y scale from the 0/20/40/60 axis ticks left of the vertical axis line
    axc = max(range(30, int(W * 0.13)),
              key=lambda x: _longest_run(dark[:base, x]))
    seg = dark[: base + 3, max(0, axc - 10):axc - 1].sum(axis=1)
    yt = _cluster([y for y in range(len(seg)) if seg[y] >= 3])
    ppc = np.median(np.diff(yt)) / 20.0
    ytop, y0 = yt[0], yt[-1]
    # x scale from the weekly tick marks below the baseline
    band = dark[base + 2:base + 9, :].sum(axis=0)
    xt = np.array(_cluster([x for x in range(W) if band[x] >= 4]))
    ppd = np.median(np.diff(xt)) / 7.0  # pixels per day
    lastx = xt[-1]                       # rightmost tick is always real
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
    barmin, barmax = np.where(tot > 2)[0].min(), np.where(tot > 2)[0].max()
    rows = []
    for off in range(-105, 4):
        cx = lastx + off * ppd
        if cx < barmin - ppd or cx > barmax + ppd:
            continue
        lo, hi = int(round(cx - ppd * 0.45)), int(round(cx + ppd * 0.45))
        cols = range(max(0, lo), min(W, hi + 1))
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
        rows = digitize(im, last_tick)
        total = sum(a + d for _, a, d in rows)
        print(f"SitRep {sr} ({report_date}): {len(rows)} onset days, "
              f"total {total} confirmed")
        for onset, alive, dead in rows:
            out_rows.append((sr, report_date, onset, alive, dead, alive + dead))
    with open(out_csv, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["sitrep", "report_date", "onset_date",
                    "confirmed_alive", "confirmed_dead", "confirmed_total"])
        w.writerows(out_rows)
    print(f"wrote {len(out_rows)} rows to {out_csv}")


if __name__ == "__main__":
    main()
