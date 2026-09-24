#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.11"
# dependencies = []
# ///
"""Extract the INRB-UMIE dashboard's symptom-onset curves from its git history.

The Trends page of https://inrb-umie.github.io/BDBV2026-Epidemic_Dashboard/
inlines one svglite bar chart of confirmed cases by symptom-onset date for
the country, each province and each health zone (stacked observed plus
imputed onset).  The page is rebuilt by CI and committed to the repository,
so every commit of ``trends.html`` is one vintage of those charts.

Bar heights are read back into counts with the chart's own axis: the
numeric y tick labels give cases per point, the major vertical gridlines
under the date labels give the day positions.

Usage::

    uv run scripts/extract_dashboard_onsets.py --repo path/to/dashboard-clone
    uv run scripts/extract_dashboard_onsets.py --html trends.html
    uv run scripts/extract_dashboard_onsets.py --self-test

National and province rows go to ``data/onset_dashboard_history.csv`` and
health-zone rows to ``data/onset_dashboard_history_zones.csv.gz``.  Consecutive
commits whose charts read back identically are kept once, under the first.

A blob-filtered clone is enough::

    git clone --filter=blob:none --no-checkout \\
        https://github.com/inrb-umie/BDBV2026-Epidemic_Dashboard.git

The script asks git for each commit's ``trends.html`` blob and git fetches
the missing blobs on demand.
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import gzip
import hashlib
import json
import re
import subprocess
import sys
from dataclasses import dataclass
from itertools import pairwise
from pathlib import Path

PAGE_FILE = "trends.html"
LEVELS = ("national", "province", "zone")
COLUMNS = (
    "snapshot_date",
    "commit_date",
    "commit_sha",
    "level",
    "unit",
    "onset_date",
    "observed",
    "imputed",
)
MONTHS = {
    m: i + 1
    for i, m in enumerate(
        ("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
    )
}
MAJOR_GRID = "stroke-width: 1.16; stroke: #EBEBEB"
PAYLOAD_RE = re.compile(r'id="payload" type="application/json">(.*?)</script>', re.DOTALL)
TEXT_RE = re.compile(r"<text([^>]*)>([^<]*)</text>")
RECT_RE = re.compile(
    r"<rect x='(-?[\d.]+)' y='(-?[\d.]+)' width='([\d.]+)' height='([\d.]+)'"
    r" style='([^']*)'"
)
CLIP_RE = re.compile(
    r"<clipPath id='[^']*'>\s*<rect x='([\d.]+)' y='([\d.]+)' width='([\d.]+)' height='([\d.]+)'"
)
LINE_RE = re.compile(r"<polyline points='([\d.]+),([\d.]+) ([\d.]+),([\d.]+) ' style='([^']*)'")
FILL_RE = re.compile(r"fill: (#[0-9A-Fa-f]{6})")
NUM_RE = re.compile(r"^-?\d+(\.\d+)?$")
MONTH_LABEL_RE = re.compile(r"^(" + "|".join(MONTHS) + r")(?: (\d{1,2}))?$")


@dataclass
class Chart:
    dates: list[dt.date]
    observed: list[int]
    imputed: list[int]
    max_frac_dev: float
    filled_gaps: int


class ChartError(ValueError):
    pass


def _texts(svg: str) -> list[tuple[float, float, str]]:
    out = []
    for tm in TEXT_RE.finditer(svg):
        attrs, content = tm.groups()
        m = re.search(r"translate\(([\d.]+),([\d.]+)\)", attrs)
        if m:
            x, y = float(m.group(1)), float(m.group(2))
        else:
            mx = re.search(r"\bx='(-?[\d.]+)'", attrs)
            my = re.search(r"\by='(-?[\d.]+)'", attrs)
            if not (mx and my):
                continue
            x, y = float(mx.group(1)), float(my.group(1))
        out.append((x, y, content.strip()))
    return out


def _panel(svg: str) -> tuple[float, float, float, float]:
    """The plot panel clip box (x0, y0, x1, y1): the largest clip that is not the full canvas."""
    clips = [tuple(map(float, m.groups())) for m in CLIP_RE.finditer(svg)]
    if not clips:
        raise ChartError("no clipPath")
    canvas = max(clips, key=lambda c: c[2] * c[3])
    inner = [c for c in clips if c != canvas]
    if not inner:
        raise ChartError("no panel clip")
    x, y, w, h = max(inner, key=lambda c: c[2] * c[3])
    return x, y, x + w, y + h


def _series_fills(
    svg: str, rects: list[tuple[float, float, float, float, str]], panel
) -> tuple[str, str]:
    """Map the legend keys to the observed and imputed fills by the label beside each key."""
    _, _, x1, y1 = panel
    keys = [(x, y, w, h, f) for x, y, w, h, f in rects if (y >= y1 or x >= x1) and w > 10]
    labels = [(x, y, t) for x, y, t in _texts(svg) if "Onset" in t or "onset" in t]
    fills: dict[str, str] = {}
    for x, y, w, h, f in keys:
        cy = y + h / 2
        near = [l for l in labels if abs(l[1] - cy) < h and l[0] > x]
        if not near:
            continue
        label = min(near, key=lambda l: l[0] - x)[2].lower()
        if "observed" in label:
            fills["observed"] = f
        elif "imputed" in label:
            fills["imputed"] = f
    if "observed" not in fills or "imputed" not in fills:
        raise ChartError(f"legend keys not resolved: {fills}")
    return fills["observed"], fills["imputed"]


def _fit_line(xs: list[float], ys: list[float]) -> tuple[float, float]:
    n = len(xs)
    mx, my = sum(xs) / n, sum(ys) / n
    sxx = sum((x - mx) ** 2 for x in xs)
    slope = sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / sxx
    return slope, my - slope * mx


def parse_chart(svg: str, snapshot: dt.date) -> Chart:
    panel = _panel(svg)
    x0, y0, x1, y1 = panel
    texts = _texts(svg)
    rects = []
    for rm in RECT_RE.finditer(svg):
        x, y, w, h, style = rm.groups()
        fill = FILL_RE.search(style)
        if fill:
            rects.append((float(x), float(y), float(w), float(h), fill.group(1).upper()))
    obs_fill, imp_fill = _series_fills(svg, rects, panel)

    # y scale: cases per point from the numeric tick labels left of the panel
    yt = [(y, float(t)) for x, y, t in texts if x < x0 and NUM_RE.match(t)]
    if len(yt) < 2:
        raise ChartError("fewer than two y tick labels")
    ppu, _ = _fit_line([v for _, v in yt], [y for y, _ in yt])
    ppu = -ppu
    if ppu <= 0:
        raise ChartError("y axis not increasing upwards")

    # x scale: major vertical gridlines sit under the date labels (a month, or a month and day)
    vgrid = sorted(
        {
            float(m.group(1))
            for m in LINE_RE.finditer(svg)
            if m.group(1) == m.group(3) and MAJOR_GRID in m.group(5)
        }
    )
    months = [
        (x, MONTH_LABEL_RE.match(t)) for x, y, t in texts if y > y1 and MONTH_LABEL_RE.match(t)
    ]
    if not vgrid or not months:
        raise ChartError("no month gridlines or labels")
    # the labels carry no year: the first break is the latest such date on or before the snapshot,
    # and later breaks that wrap past December move into the next year
    anchors = []
    prev = None
    for gx in vgrid:
        lab = min(months, key=lambda m: abs(m[0] - gx))[1]
        mon, day = MONTHS[lab.group(1)], int(lab.group(2) or 1)
        if prev is None:
            year = (
                snapshot.year if (mon, day) <= (snapshot.month, snapshot.day) else snapshot.year - 1
            )
        elif mon < prev:
            year += 1
        prev = mon
        anchors.append((gx, dt.date(year, mon, day)))

    # a one-day chart draws its bar wider than the panel, so only the vertical extent is checked
    bars = [
        (x, y, w, h, f)
        for x, y, w, h, f in rects
        if f in (obs_fill, imp_fill) and x + w > x0 and x < x1 and y0 - 1 <= y and y + h <= y1 + 1
    ]
    if not bars:
        raise ChartError("no bars")
    width = max(w for _, _, w, _, _ in bars)
    centres = sorted({round(x + w / 2, 3) for x, _, w, _, _ in bars})
    if len(anchors) >= 2:
        ppd, _ = _fit_line([d.toordinal() for _, d in anchors], [gx for gx, _ in anchors])
    else:
        steps = [b - a for a, b in pairwise(centres)]
        ppd = min(steps) if steps else width
    gx, gd = anchors[0]
    by_day: dict[dt.date, list[float]] = {}
    max_off = 0.0
    for x, y, w, h, f in bars:
        rel = (x + w / 2 - gx) / ppd
        k = round(rel)
        max_off = max(max_off, abs(rel - k))
        day = gd + dt.timedelta(days=k)
        cell = by_day.setdefault(day, [0.0, 0.0])
        cell[0 if f == obs_fill else 1] += h / ppu
    if max_off > 0.2:
        raise ChartError(f"bar centres are {max_off:.2f} days off the day grid")
    first, last = min(by_day), max(by_day)
    dates = [first + dt.timedelta(days=i) for i in range((last - first).days + 1)]
    obs, imp, dev, gaps = [], [], 0.0, 0
    for d in dates:
        if d not in by_day:
            gaps += 1
        o, i = by_day.get(d, (0.0, 0.0))
        dev = max(dev, abs(o - round(o)), abs(i - round(i)))
        obs.append(round(o))
        imp.append(round(i))
    if dev > 0.05:
        raise ChartError(f"bar height {dev:.3f} cases off an integer")
    return Chart(dates, obs, imp, dev, gaps)


def snapshot_date_of(entry: dict) -> dt.date:
    m = re.match(r"(\d{4}-\d{2}-\d{2})/", entry.get("file", ""))
    if not m:
        raise ChartError(f"no data date in file path {entry.get('file')!r}")
    return dt.date.fromisoformat(m.group(1))


def extract_page(html: str, levels: tuple[str, ...]) -> tuple[dt.date, list[tuple]]:
    """Rows (level, unit, onset_date, observed, imputed) for one page, plus its data date."""
    m = PAYLOAD_RE.search(html)
    if not m:
        raise ChartError("no payload script")
    trends = json.loads(m.group(1)).get("onset_trends") or {}
    units: list[tuple[str, str, dict]] = []
    if "national" in levels and trends.get("national", {}).get("svg"):
        units.append(("national", "National", trends["national"]))
    if "province" in levels:
        units += [
            ("province", k, v) for k, v in (trends.get("provinces") or {}).items() if v.get("svg")
        ]
    if "zone" in levels:
        units += [
            ("zone", k, v) for k, v in (trends.get("health_zones") or {}).items() if v.get("svg")
        ]
    if not units:
        raise ChartError("no onset charts on the page")
    snapshot = snapshot_date_of(units[0][2])
    rows = []
    for level, unit, entry in units:
        try:
            chart = parse_chart(entry["svg"], snapshot)
        except ChartError as err:
            print(f"  skip {level} {unit}: {err}", file=sys.stderr)
            continue
        if chart.filled_gaps:
            print(
                f"  {level} {unit}: {chart.filled_gaps} days without a bar filled with zero",
                file=sys.stderr,
            )
        rows += [
            (level, unit, d.isoformat(), o, i)
            for d, o, i in zip(chart.dates, chart.observed, chart.imputed)
        ]
    return snapshot, rows


def tag_repeats(snapshots):
    """Yield each (sha, commit_date, snapshot, rows) with a flag saying its
    rows equal the previous snapshot's, so consecutive rebuilds of the
    same data are kept once, under the first."""
    previous = None
    for item in snapshots:
        digest = hashlib.sha1(repr(item[3]).encode()).hexdigest()
        yield item, digest == previous
        previous = digest


def git(repo: Path, *args: str) -> str:
    return subprocess.run(
        ["git", "-C", str(repo), *args], check=True, capture_output=True, text=True
    ).stdout


def page_commits(repo: Path, commits: list[str] | None) -> list[tuple[str, str]]:
    """(sha, iso commit date) oldest first for every commit that touches the page."""
    if commits is None:
        log = git(repo, "log", "--reverse", "--format=%H %cI", "--", PAGE_FILE)
        return [tuple(line.split()) for line in log.splitlines()]
    out = []
    for c in commits:
        out.append(
            (git(repo, "rev-parse", c).strip(), git(repo, "log", "-1", "--format=%cI", c).strip())
        )
    return out


SELF_TEST_SVG = """
<svg>
<clipPath id='c0'><rect x='0' y='0' width='300' height='150'/></clipPath>
<clipPath id='c1'><rect x='50' y='10' width='200' height='100'/></clipPath>
<text x='30' y='110'>0</text>
<text x='30' y='60'>5</text>
<text x='30' y='10'>10</text>
<polyline points='60,10 60,110 ' style='stroke-width: 1.16; stroke: #EBEBEB;' />
<polyline points='130,10 130,110 ' style='stroke-width: 1.16; stroke: #EBEBEB;' />
<polyline points='50,60 250,60 ' style='stroke-width: 1.16; stroke: #EBEBEB;' />
<text x='60' y='125'>Aug 01</text>
<text x='130' y='125'>Aug 08</text>
<rect x='56' y='80' width='8' height='30' style='fill: #1F77B4;' />
<rect x='56' y='70' width='8' height='10' style='fill: #FF7F0E;' />
<rect x='66' y='90' width='8' height='20' style='fill: #1F77B4;' />
<rect x='86' y='100' width='8' height='10' style='fill: #1F77B4;' />
<rect x='60' y='130' width='15' height='10' style='fill: #1F77B4;' />
<text x='80' y='135'>Observed onset</text>
<rect x='150' y='130' width='15' height='10' style='fill: #FF7F0E;' />
<text x='170' y='135'>Imputed onset</text>
</svg>
"""


def self_test() -> int:
    """Parse a hand-written chart through the tick-label calibration and
    the bar reading, and check the repeat tagging; non-zero on failure."""
    chart = parse_chart(SELF_TEST_SVG, dt.date(2026, 8, 10))
    assert chart.dates[0] == dt.date(2026, 8, 1), chart.dates
    assert chart.dates[-1] == dt.date(2026, 8, 4), chart.dates
    assert chart.observed == [3, 2, 0, 1], chart.observed
    assert chart.imputed == [1, 0, 0, 0], chart.imputed
    assert chart.filled_gaps == 1
    # the labels carry no year: a break after the snapshot's month-day is
    # last year's
    chart = parse_chart(SELF_TEST_SVG, dt.date(2026, 7, 10))
    assert chart.dates[0] == dt.date(2025, 8, 1), chart.dates
    # a bar off the day grid is refused
    bad = SELF_TEST_SVG.replace("<rect x='86'", "<rect x='83'")
    try:
        parse_chart(bad, dt.date(2026, 8, 10))
    except ChartError:
        pass
    else:
        raise AssertionError("off-grid bar accepted")
    a = [("national", "National", "2026-08-01", 3, 1)]
    b = [("national", "National", "2026-08-01", 4, 1)]
    items = [("s1", "d1", None, a), ("s2", "d2", None, a), ("s3", "d3", None, b),
             ("s4", "d4", None, a)]
    flags = [repeat for _, repeat in tag_repeats(items)]
    assert flags == [False, True, False, False], flags
    print("self-test passed")
    return 0


def main() -> int:
    if sys.argv[1:] == ["--self-test"]:
        return self_test()
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--repo", type=Path, help="clone of INRB-UMIE/BDBV2026-Epidemic_Dashboard")
    src.add_argument(
        "--html", type=Path, help="one saved trends.html (written as a single snapshot)"
    )
    ap.add_argument(
        "--commits", type=Path, help="file of commit SHAs to read instead of the page's full log"
    )
    ap.add_argument(
        "--levels", default="national,province,zone", help="comma list from national,province,zone"
    )
    ap.add_argument("--self-test", action="store_true", help="check the parser on a built-in chart")
    ap.add_argument(
        "--out",
        type=Path,
        default=Path("data/onset_dashboard_history.csv"),
        help="national and province rows",
    )
    ap.add_argument(
        "--zones-out",
        type=Path,
        default=Path("data/onset_dashboard_history_zones.csv.gz"),
        help="health-zone rows (a much larger file, gzipped when the name ends in .gz)",
    )
    args = ap.parse_args()
    levels = tuple(l for l in args.levels.split(",") if l in LEVELS)

    snapshots: list[tuple[str, str, dt.date, list[tuple]]] = []
    if args.html:
        html = args.html.read_text()
        snapshot, rows = extract_page(html, levels)
        snapshots.append(("", "", snapshot, rows))
    else:
        commits = args.commits.read_text().split() if args.commits else None

        def parsed():
            for sha, cdate in page_commits(args.repo, commits):
                html = git(args.repo, "cat-file", "blob", f"{sha}:{PAGE_FILE}")
                try:
                    snapshot, rows = extract_page(html, levels)
                except ChartError as err:
                    print(f"{sha[:8]} {cdate}: {err}", file=sys.stderr)
                    continue
                yield sha, cdate, snapshot, rows

        for (sha, cdate, snapshot, rows), repeat in tag_repeats(parsed()):
            nat = [r for r in rows if r[0] == "national"]
            print(
                f"{sha[:8]} {cdate[:10]} data {snapshot} units {len({(r[0], r[1]) for r in rows})}"
                f" national observed {sum(r[3] for r in nat)} imputed {sum(r[4] for r in nat)}"
                + ("  (same as previous, dropped)" if repeat else ""),
                file=sys.stderr,
            )
            if not repeat:
                snapshots.append((sha, cdate, snapshot, rows))

    order = {l: i for i, l in enumerate(LEVELS)}
    outputs = [(args.out, ("national", "province")), (args.zones_out, ("zone",))]
    for path, keep in outputs:
        if not any(l in levels for l in keep):
            continue
        path.parent.mkdir(parents=True, exist_ok=True)
        n = 0
        opener = gzip.open if path.suffix == ".gz" else open
        with opener(path, "wt", newline="") as f:
            w = csv.writer(f, lineterminator="\n")
            w.writerow(COLUMNS)
            for sha, cdate, snapshot, rows in snapshots:
                for level, unit, day, o, i in sorted(rows, key=lambda r: (order[r[0]], r[1], r[2])):
                    if level in keep:
                        w.writerow([snapshot.isoformat(), cdate, sha, level, unit, day, o, i])
                        n += 1
        print(f"wrote {path}: {len(snapshots)} snapshots, {n} rows", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
