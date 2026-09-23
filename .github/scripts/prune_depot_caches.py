#!/usr/bin/env python3
"""Prune the repository's Julia depot caches so the fit caches are not evicted.

GitHub evicts the least recently used caches once a repository holds more
than 10 GB. The depot caches written by `julia-actions/cache` are 0.5-1 GB
each and carry the run id in their key, so every run adds a new one. That
action deletes the superseded ones on a pull request's ref but never on the
default branch, so `main` keeps one per run and job until eviction removes
them, along with the fit caches the docs build takes hours to rebuild.

Only keys ending in `;run_id=<n>;run_attempt=<n>` are candidates, which is
the suffix `julia-actions/cache` appends. The `fit-*` caches, the
pre-commit cache and anything else never match, so they are never deleted.

Deletes, in order:

1. every depot cache that a newer one with the same ref and key prefix has
   superseded, since a restore takes the newest match;
2. depot caches on the default branch not used for `--stale-days`, which is
   what a bumped `cache-name` leaves behind;
3. depot caches on other refs, least recently used first, until all caches
   together fit in `--budget-gb`.

Lists what it would delete unless `--delete` is passed.
"""

import argparse
import json
import os
import re
import subprocess
import sys
from datetime import datetime, timedelta, timezone

DEPOT_KEY = re.compile(r"^(?P<prefix>.*;)run_id=\d+;run_attempt=\d+$")
GB = 1024**3


def gh(*args):
    return subprocess.run(
        ["gh", *args], check=True, capture_output=True, text=True
    ).stdout


def list_caches(repo):
    out = gh(
        "api",
        "--paginate",
        f"repos/{repo}/actions/caches?per_page=100",
        "--jq",
        ".actions_caches[]",
    )
    return [json.loads(line) for line in out.splitlines() if line.strip()]


def timestamp(value):
    return datetime.fromisoformat(value.replace("Z", "+00:00"))


def plan(caches, default_ref, budget, stale_days, now):
    """Return the caches to delete as (cache, reason) pairs."""
    depot = [c for c in caches if DEPOT_KEY.match(c["key"])]
    doomed = {}

    groups = {}
    for c in depot:
        prefix = DEPOT_KEY.match(c["key"])["prefix"]
        groups.setdefault((c["ref"], prefix), []).append(c)
    for group in groups.values():
        group.sort(key=lambda c: timestamp(c["created_at"]), reverse=True)
        for c in group[1:]:
            doomed[c["id"]] = (c, "superseded")

    cutoff = now - timedelta(days=stale_days)
    for c in depot:
        if (
            c["id"] not in doomed
            and c["ref"] == default_ref
            and timestamp(c["last_accessed_at"]) < cutoff
        ):
            doomed[c["id"]] = (c, f"unused for {stale_days} days")

    total = sum(c["size_in_bytes"] for c in caches if c["id"] not in doomed)
    others = sorted(
        (c for c in depot if c["id"] not in doomed and c["ref"] != default_ref),
        key=lambda c: timestamp(c["last_accessed_at"]),
    )
    for c in others:
        if total <= budget:
            break
        doomed[c["id"]] = (c, "over budget")
        total -= c["size_in_bytes"]

    return list(doomed.values())


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--repo", default=os.environ.get("GITHUB_REPOSITORY"))
    parser.add_argument("--default-branch", default="main")
    parser.add_argument("--budget-gb", type=float, default=8.0)
    parser.add_argument("--stale-days", type=int, default=7)
    parser.add_argument("--delete", action="store_true")
    args = parser.parse_args()
    if not args.repo:
        parser.error("--repo is required outside GitHub Actions")

    caches = list_caches(args.repo)
    doomed = plan(
        caches,
        f"refs/heads/{args.default_branch}",
        args.budget_gb * GB,
        args.stale_days,
        datetime.now(timezone.utc),
    )

    freed = sum(c["size_in_bytes"] for c, _ in doomed)
    total = sum(c["size_in_bytes"] for c in caches)
    kept = [c for c in caches if not DEPOT_KEY.match(c["key"])]
    lines = [
        f"{len(caches)} caches, {total / GB:.2f} GB; "
        f"{len(kept)} are not depot caches and are never deleted "
        f"({sum(c['size_in_bytes'] for c in kept) / GB:.2f} GB).",
        f"{'Deleting' if args.delete else 'Would delete'} {len(doomed)}, "
        f"{freed / GB:.2f} GB, leaving {(total - freed) / GB:.2f} GB.",
        "",
    ]
    for c, reason in sorted(doomed, key=lambda d: (d[1], d[0]["ref"])):
        lines.append(
            f"- {c['size_in_bytes'] / 2**20:7.0f} MB  {reason:<18} "
            f"{c['ref']}  `{c['key']}`"
        )
    report = "\n".join(lines)
    print(report)
    summary = os.environ.get("GITHUB_STEP_SUMMARY")
    if summary:
        with open(summary, "a") as f:
            f.write("### Depot cache prune\n\n" + report + "\n")

    if args.delete:
        failed = 0
        for c, _ in doomed:
            try:
                gh("cache", "delete", str(c["id"]), "-R", args.repo)
            except subprocess.CalledProcessError as e:
                # Another prune, or the PR cleanup, may have got there first.
                print(f"could not delete {c['id']}: {e.stderr.strip()}")
                failed += 1
        print(f"deleted {len(doomed) - failed} of {len(doomed)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
