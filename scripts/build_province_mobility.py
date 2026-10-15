#!/usr/bin/env python3
"""Print the Flowminder mobility shares of each origin province.

Reads Flowminder's two cohort subscriber-day files, as published on the
Humanitarian Data Exchange and mirrored in
https://github.com/INRB-UMIE/BDBV2026-Data under
``data/flowminder_short_trips/raw``:

- ``drc-bvd_ituri-cohort_subscriber-days-2026_06_08-v1.0-external.csv``
- ``drc-bvd_nk-cohort_subscriber-days-2026_06_08-v1.0-external.csv``

Each gives, per health zone, the average presence days per member of an
anonymised subscriber cohort: subscribers present in the Ituri outbreak
zones (Bunia, Mongbwalu, Nyankunde, Rwampara) or the Nord-Kivu ones (Beni,
Butembo, Katwa) for at least two days during 4-17 May 2026. The origin
zones carry no value, so movement within them is not captured. This sums
the follow-up period (18 May - 8 June 2026) over the zones of each province
in ``PROVINCE_SOURCE_NAMES`` order and prints the
``PROVINCE_SOURCE_MOBILITY`` literal in src/constants.jl. Days in provinces
outside that list are dropped and their share is printed as a comment.

Needs only the Python standard library.

Usage:

    python3 scripts/build_province_mobility.py path/to/raw
"""

import csv
import sys
from pathlib import Path

# PROVINCE_SOURCE_NAMES order, keyed by the files' province name.
PROVINCES = [
    ("ituri", "Ituri"),
    ("nord_kivu", "Nord-Kivu"),
    ("sud_kivu", "Sud-Kivu"),
    ("haut_uele", "Haut-Uele"),
    ("tshopo", "Tshopo"),
    ("bas_uele", "Bas-Uele"),
    ("sud_ubangi", "Sud-Ubangi"),
]

# Origin province and its cohort file.
ORIGINS = [
    ("ituri", "drc-bvd_ituri-cohort_subscriber-days-2026_06_08-v1.0-external.csv"),
    ("nord_kivu", "drc-bvd_nk-cohort_subscriber-days-2026_06_08-v1.0-external.csv"),
]

COLUMN = "Avg days (follow-up)"


def main(raw):
    names = {name: key for key, name in PROVINCES}
    print("const PROVINCE_SOURCE_MOBILITY = Dict(")
    for origin, fname in ORIGINS:
        totals = {key: 0.0 for key, _ in PROVINCES}
        outside = 0.0
        with open(Path(raw) / fname, encoding="utf-8") as f:
            for row in csv.DictReader(f):
                if not row[COLUMN]:
                    continue
                value = float(row[COLUMN])
                key = names.get(row["Province"])
                if key is None:
                    outside += value
                else:
                    totals[key] += value
        away = sum(v for k, v in totals.items() if k != origin)
        dropped = outside / (outside + away)
        values = ", ".join(f"{totals[key]:.6g}" for key, _ in PROVINCES)
        print(
            f'    "{origin}" => [{values}],'
            f"  # {dropped:.1%} of away days outside these provinces"
        )
    print(")")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    main(sys.argv[1])
