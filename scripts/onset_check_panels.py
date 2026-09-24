#!/usr/bin/env python3
# /// script
# requires-python = ">=3.9"
# dependencies = ["pillow", "numpy"]
# ///
#
# Crop panels for the vision check of one digitised onset-curve vintage.
#
# Each panel is an upscaled crop of the embedded figure covering a run of
# days, with the day of month printed under every bar and a count ruler on
# the left drawn from the digitiser's own calibration. Two panels are
# written per run of days:
#   blind_<sr>_<k>.png  the crop and rulers only, for a blind read
#   check_<sr>_<k>.png  the same crop with the digitised total (green) and
#                       alive (magenta) tops drawn per day and the total
#                       printed above each bar
# The figure is extracted from the SitRep PDF and calibrated exactly as
# scripts/digitize_onset_curve.py does, so the rulers show where the
# digitiser placed each day and each count. The digitised block comes
# from data/onset_curve_scanned.csv.
#
# Usage:
#   uv run scripts/onset_check_panels.py SR [--days 12] [--scale 8]
#       [--pdf-dir data/sitrep_pdfs] [--csv data/onset_curve_scanned.csv]
#       [--out output/onset_panels]
# Prints one line per panel with its date range and the source columns.

import argparse
import csv
import datetime as dt
import os
import sys
import tempfile

import numpy as np
from PIL import Image, ImageDraw

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import digitize_onset_curve as d


def digitised_block(path, sr):
    rows = {}
    with open(path, newline="") as f:
        r = csv.reader(f)
        next(r)
        for row in r:
            if row and row[0] == sr:
                rows[row[2]] = (int(row[3]), int(row[4]))
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("sr")
    ap.add_argument("--days", type=int, default=12)
    ap.add_argument("--scale", type=int, default=8)
    ap.add_argument("--pdf-dir", default="data/sitrep_pdfs")
    ap.add_argument("--csv", default="data/onset_curve_scanned.csv")
    ap.add_argument("--out", default="output/onset_panels")
    a = ap.parse_args()
    sr, ndays, scale = a.sr, a.days, a.scale
    if sr not in d.CONFIG:
        sys.exit(f"{sr} is not in CONFIG; add its row to the digitiser first")
    pdf = os.path.join(a.pdf_dir, f"SitRep_MVE_{sr}_2026.pdf")
    with tempfile.TemporaryDirectory() as wd:
        im = d.extract_onset_image(pdf, wd)
    if im is None:
        sys.exit(f"no onset curve found in {pdf}")
    cal = d.calibrate(im, d.Y_AXIS_STEP.get(sr, 20))
    rows = digitised_block(a.csv, sr)
    if not rows:
        sys.exit(f"no digitised rows for {sr} in {a.csv}")
    os.makedirs(a.out, exist_ok=True)
    lastdate = dt.date.fromisoformat(d.CONFIG[sr][1])
    ppd, ppc, base = cal["ppd"], cal["ppc"], cal["y0"]
    ytop, H, W = cal["yt"][0], cal["H"], cal["W"]
    y0, y1 = max(0, ytop - 8), min(H, base + 4)
    ymax = int((base - ytop) / ppc)
    dates = sorted(rows)
    first, last_day = dt.date.fromisoformat(dates[0]), dt.date.fromisoformat(dates[-1])
    k = 0
    while first <= last_day:
        last = first + dt.timedelta(days=ndays - 1)
        off0, off1 = (first - lastdate).days, (last - lastdate).days
        # day_column is 1-based; the crop is indexed 0-based
        x0 = int(np.floor(d.day_column(cal, off0) - 1 - ppd)) - 1
        x1 = int(np.ceil(d.day_column(cal, off1) - 1 + ppd)) + 1
        x0, x1 = max(0, x0), min(W, x1)
        cw, ch = (x1 - x0) * scale, (y1 - y0) * scale
        pad_l, pad_b = 60, 40
        for mode in ("blind", "check"):
            canvas = Image.new("RGB", (cw + pad_l, ch + pad_b), (255, 255, 255))
            crop = Image.fromarray(im[y0:y1, x0:x1].astype(np.uint8))
            canvas.paste(crop.resize((cw, ch), Image.NEAREST), (pad_l, 0))
            dr = ImageDraw.Draw(canvas)
            for c in range(0, ymax + 1, 5):
                y = (base - c * ppc - y0) * scale
                col = (0, 0, 0) if c % 25 == 0 else (150, 150, 150)
                dr.line([(pad_l - 8, y), (pad_l + cw, y)], fill=col, width=1)
                dr.text((2, y - 5), str(c), fill=(0, 0, 0))
            for off in range(off0, off1 + 1):
                date = lastdate + dt.timedelta(days=off)
                cx = d.day_column(cal, off) - 1
                X = (cx - x0) * scale + pad_l
                dr.line([(X, ch), (X, ch + 6)], fill=(0, 0, 0))
                dr.text((X - 6, ch + 8), f"{date.day:02d}", fill=(0, 0, 0))
                if date.day == 1 or off == off0:
                    dr.text((X - 10, ch + 22), date.strftime("%b"), fill=(0, 0, 0))
                if mode == "check" and date.isoformat() in rows:
                    alive, dead = rows[date.isoformat()]
                    lo = (cx - ppd * 0.45 - x0) * scale + pad_l
                    hi = (cx + ppd * 0.45 - x0) * scale + pad_l
                    yt = (base - (alive + dead) * ppc - y0) * scale
                    ya = (base - alive * ppc - y0) * scale
                    dr.line([(lo, yt), (hi, yt)], fill=(0, 200, 0), width=2)
                    if alive > 0 and dead > 0:
                        dr.line([(lo, ya), (hi, ya)], fill=(255, 0, 255), width=2)
                    dr.text(((lo + hi) / 2 - 6, yt - 14), str(alive + dead),
                            fill=(0, 120, 0))
            canvas.save(os.path.join(a.out, f"{mode}_{sr}_{k}.png"))
        print(f"panel {k}: {first} .. {last}  cols {x0}-{x1}")
        first = last + dt.timedelta(days=1)
        k += 1


if __name__ == "__main__":
    main()
