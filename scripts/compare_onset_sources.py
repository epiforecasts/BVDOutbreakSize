#!/usr/bin/env python3
# /// script
# requires-python = ">=3.9"
# dependencies = []
# ///
#
# Compare each digitised SitRep onset block (data/onset_curve_scanned.csv)
# with the INRB-UMIE dashboard's national onset curve
# (data/onset_dashboard_history.csv) from the snapshot nearest its report
# date, within two days.
#
# Differences on settled dates (onset more than 56 days before the later
# of the report date and the snapshot date) measure whether the two sources
# carry the same records, since late reporting has finished there.
# Differences on the last 14 days measure reporting lag between the two
# line lists. The two are reported separately and only the settled-date
# one is a check: the run exits non-zero when any block's settled-date mean
# absolute difference exceeds THRESHOLD cases per day.
#
# Usage:
#   uv run scripts/compare_onset_sources.py [--latest N]
#       [--scan data/onset_curve_scanned.csv]
#       [--dashboard data/onset_dashboard_history.csv]
#   uv run scripts/compare_onset_sources.py --self-test
# One row per block: sitrep, report date, snapshot date, scan total,
# dashboard observed, observed plus imputed, common days, r on common days,
# mean absolute daily difference on settled dates and on the last 14 days,
# then the five onset dates with the largest absolute difference.

import argparse
import csv
import datetime as dt
import math
import os
import subprocess
import sys
import tempfile

THRESHOLD = 2.0
SETTLED_DAYS = 56
EDGE_DAYS = 14
WINDOW_DAYS = 2


def read_scan(path):
    blocks = {}
    report = {}
    with open(path, newline="") as f:
        for row in csv.DictReader(f):
            sr = row["sitrep"]
            report[sr] = dt.date.fromisoformat(row["report_date"])
            blocks.setdefault(sr, {})[
                dt.date.fromisoformat(row["onset_date"])
            ] = int(row["confirmed_total"])
    return blocks, report


def read_dashboard(path):
    snaps = {}
    with open(path, newline="") as f:
        for row in csv.DictReader(f):
            if row["level"] != "national":
                continue
            snap = dt.date.fromisoformat(row["snapshot_date"])
            snaps.setdefault(snap, {})[dt.date.fromisoformat(row["onset_date"])] = (
                int(row["observed"]),
                int(row["imputed"]),
            )
    return snaps


def nearest_snapshot(snaps, report_date):
    best = None
    for snap in snaps:
        d = abs((snap - report_date).days)
        if d > WINDOW_DAYS:
            continue
        key = (d, snap)
        if best is None or key < best:
            best = key
    return None if best is None else best[1]


def pearson(xs, ys):
    n = len(xs)
    if n < 2:
        return float("nan")
    mx, my = sum(xs) / n, sum(ys) / n
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    sxx = sum((x - mx) ** 2 for x in xs)
    syy = sum((y - my) ** 2 for y in ys)
    if sxx == 0 or syy == 0:
        return float("nan")
    return sxy / math.sqrt(sxx * syy)


def mean_abs(diffs):
    return sum(abs(d) for d in diffs) / len(diffs) if diffs else float("nan")


def compare(block, report_date, snap, snap_date):
    common = sorted(set(block) & set(snap))
    later = max(report_date, snap_date)
    diffs = {d: block[d] - snap[d][0] for d in common}
    settled = [diffs[d] for d in common if (later - d).days > SETTLED_DAYS]
    edge = [diffs[d] for d in common if (later - d).days <= EDGE_DAYS]
    worst = sorted(common, key=lambda d: (-abs(diffs[d]), d))[:5]
    return {
        "scan_total": sum(block.values()),
        "observed": sum(v[0] for v in snap.values()),
        "observed_imputed": sum(v[0] + v[1] for v in snap.values()),
        "common": len(common),
        "r": pearson([block[d] for d in common], [snap[d][0] for d in common]),
        "settled_mad": mean_abs(settled),
        "edge_mad": mean_abs(edge),
        "worst": [(d, block[d], snap[d][0]) for d in worst],
    }


def fmt(x):
    return "nan" if math.isnan(x) else f"{x:.2f}"


def self_test():
    """Exercise the pure functions and the threshold exit on synthetic
    inputs; raise AssertionError on the first failure."""
    assert abs(pearson([1, 2, 3], [2, 4, 6]) - 1.0) < 1e-12
    assert abs(pearson([1, 2, 3], [3, 2, 1]) + 1.0) < 1e-12
    assert math.isnan(pearson([1, 1, 1], [1, 2, 3]))
    assert math.isnan(pearson([1], [1]))
    assert mean_abs([1, -3]) == 2.0
    assert math.isnan(mean_abs([]))
    d = dt.date
    snaps = {d(2026, 8, 1), d(2026, 8, 4), d(2026, 8, 6)}
    assert nearest_snapshot(snaps, d(2026, 8, 5)) == d(2026, 8, 4)
    assert nearest_snapshot(snaps, d(2026, 8, 3)) == d(2026, 8, 4)
    assert nearest_snapshot(snaps, d(2026, 8, 2)) == d(2026, 8, 1)  # tie: earlier
    assert nearest_snapshot(snaps, d(2026, 8, 9)) is None
    report, snap_date = d(2026, 9, 1), d(2026, 8, 31)
    block = {
        report - dt.timedelta(days=60): 10,  # settled
        report - dt.timedelta(days=56): 10,  # neither settled nor edge
        report - dt.timedelta(days=3): 10,   # edge
        report - dt.timedelta(days=1): 4,    # edge
    }
    snap = {k: (v, 1) for k, v in block.items()}
    snap[report - dt.timedelta(days=60)] = (7, 0)
    snap[report - dt.timedelta(days=3)] = (5, 0)
    snap[report + dt.timedelta(days=1)] = (2, 0)  # not a common day
    c = compare(block, report, snap, snap_date)
    assert c["common"] == 4 and c["scan_total"] == 34
    assert c["observed"] == 28 and c["observed_imputed"] == 30
    assert c["settled_mad"] == 3.0 and c["edge_mad"] == 2.5
    assert c["worst"][0] == (report - dt.timedelta(days=3), 10, 5)
    # the threshold exit, end to end on two small files
    with tempfile.TemporaryDirectory() as wd:
        scan = os.path.join(wd, "scan.csv")
        dash = os.path.join(wd, "dash.csv")
        with open(scan, "w") as f:
            f.write("sitrep,report_date,onset_date,confirmed_alive,"
                    "confirmed_dead,confirmed_total\n")
            for day, (n, ok) in ((60, (10, 7)), (61, (10, 10))):
                on = (report - dt.timedelta(days=day)).isoformat()
                f.write(f"001,{report},{on},{n},0,{n}\n")
                f.write(f"002,{report},{on},{ok},0,{ok}\n")
        with open(dash, "w") as f:
            f.write("snapshot_date,commit_date,commit_sha,level,unit,"
                    "onset_date,observed,imputed\n")
            for day in (60, 61):
                on = (report - dt.timedelta(days=day)).isoformat()
                f.write(f"{snap_date},x,x,national,National,{on},7,0\n")
                f.write(f"{snap_date},x,x,province,P,{on},99,0\n")
        cmd = [sys.executable, __file__, "--scan", scan, "--dashboard", dash]
        r = subprocess.run(cmd, capture_output=True, text=True, check=False)
        assert r.returncode == 1, r
        assert "001 (3.00)" in r.stderr and "002" not in r.stderr, r.stderr
        assert r.stdout.count("largest:") == 2, r.stdout
        r = subprocess.run(cmd + ["--latest", "1"], capture_output=True,
                           text=True, check=False)
        assert r.returncode == 0 and "002 " in r.stdout, r
    print("self-test passed")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--latest", type=int, default=None)
    ap.add_argument("--scan", default="data/onset_curve_scanned.csv")
    ap.add_argument("--dashboard", default="data/onset_dashboard_history.csv")
    ap.add_argument("--self-test", action="store_true")
    a = ap.parse_args()
    if a.self_test:
        self_test()
        return
    blocks, report = read_scan(a.scan)
    snaps = read_dashboard(a.dashboard)
    order = sorted(blocks)
    if a.latest is not None:
        order = order[-a.latest:]
    header = (
        "sitrep report     snapshot   scan  dash  dash+imp common r      "
        "settled  edge"
    )
    print(header)
    failed = []
    for sr in order:
        snap_date = nearest_snapshot(snaps, report[sr])
        if snap_date is None:
            continue
        c = compare(blocks[sr], report[sr], snaps[snap_date], snap_date)
        print(
            f"{sr:6} {report[sr]} {snap_date} {c['scan_total']:5} "
            f"{c['observed']:5} {c['observed_imputed']:8} {c['common']:6} "
            f"{fmt(c['r']):6} {fmt(c['settled_mad']):7} {fmt(c['edge_mad'])}"
        )
        print(
            "       largest: "
            + ", ".join(f"{d} scan {s} dash {o}" for d, s, o in c["worst"])
        )
        if not math.isnan(c["settled_mad"]) and c["settled_mad"] > THRESHOLD:
            failed.append((sr, c["settled_mad"]))
    if failed:
        sys.stdout.flush()
        print(
            "settled-date mean absolute difference above "
            f"{THRESHOLD} cases per day: "
            + ", ".join(f"{sr} ({m:.2f})" for sr, m in failed),
            file=sys.stderr,
        )
        sys.exit(1)


if __name__ == "__main__":
    main()
