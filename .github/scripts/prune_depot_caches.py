#!/usr/bin/env python3
"""Prune the repository's depot and fit caches so main's fits are not evicted.

GitHub evicts the least recently used caches once a repository holds more
than 10 GB. The depot caches written by `julia-actions/cache` are 0.5-1 GB
each and carry the run id in their key, so every run adds a new one. That
action deletes the superseded ones on a pull request's ref but never on the
default branch, so `main` keeps one per run and job until eviction removes
them, along with the fit caches the docs build takes hours to rebuild.

Depot caches are the keys ending in `;run_id=<n>;run_attempt=<n>`, the
suffix `julia-actions/cache` appends. Fit caches are `fit-<id>-<hash>`, one
fitted chain per fit id and content hash. A fit job restores only its own
hash, so an older fit of the same id on the same ref is never restored
again. The `warmstart-*`, `recovery-*` and pre-commit caches and anything
else match neither, so they are never deleted.

A `docs-env` cache whose run has not finished is never deleted either. The
docs build's later jobs restore it hours after `list` saves it, and load the
compiled package `list` built against it. If it is gone, the first job to
miss it compiles every dependency and saves a new depot under the same key.
Every later job then rejects that compiled package, whose dependencies have
different build ids, and compiles it again.

Deletes, in order:

1. every depot cache that a newer one with the same ref and key prefix has
   superseded, since a restore takes the newest match;
2. depot caches on the default branch not used for `--stale-days`, which is
   what a bumped `cache-name` leaves behind;
3. fit caches that a newer fit of the same id on the same ref has
   superseded, unless created or used within `--fit-grace-hours`, which
   covers a build still running and a re-run of its failed jobs;
4. depot caches on other refs, least recently used first, until all caches
   together fit in `--budget-gb`;
5. fit caches on other refs, least recently used first, until all caches
   together fit in `--budget-gb`, skipping any ref with a docs build queued
   or running and any fit inside `--fit-grace-hours`. A re-run of that pull
   request's build then refits them. The default branch's fits are kept,
   since every pull request's build restores from them.

Lists what it would delete unless `--delete` is passed.
"""

import argparse
import json
import os
import re
import subprocess
import sys
from datetime import datetime, timedelta, timezone

DEPOT_KEY = re.compile(r"^(?P<prefix>.*;)run_id=(?P<run_id>\d+);run_attempt=\d+$")
FIT_KEY = re.compile(r"^fit-(?P<id>.+)-(?P<hash>[0-9a-f]{16})$")
HANDOFF_PREFIX = "docs-env;"
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


def in_flight(repo, caches):
    """Return the run ids of `docs-env` caches whose run has not finished."""
    running = set()
    for c in caches:
        m = DEPOT_KEY.match(c["key"])
        if not (m and c["key"].startswith(HANDOFF_PREFIX)):
            continue
        try:
            status = gh(
                "api",
                f"repos/{repo}/actions/runs/{m['run_id']}",
                "--jq",
                ".status",
            ).strip()
        except subprocess.CalledProcessError:
            # A deleted run cannot restore anything.
            continue
        if status != "completed":
            running.add(m["run_id"])
    return running


def kept_for_run(cache, running):
    """Whether `cache` is a `docs-env` cache a running build may restore."""
    m = DEPOT_KEY.match(cache["key"])
    return bool(
        m and cache["key"].startswith(HANDOFF_PREFIX) and m["run_id"] in running
    )


def superseded_fits(caches, now, grace_hours):
    """Return the fit caches a newer fit of the same id and ref supersedes.

    The newest per id and ref is kept, and so is any created or used within
    `grace_hours`.
    """
    groups = {}
    for c in caches:
        m = FIT_KEY.match(c["key"])
        if m:
            groups.setdefault((c["ref"], m["id"]), []).append(c)
    recent = now - timedelta(hours=grace_hours)
    old = []
    for group in groups.values():
        group.sort(key=lambda c: timestamp(c["created_at"]), reverse=True)
        old.extend(
            c
            for c in group[1:]
            if timestamp(c["created_at"]) < recent
            and timestamp(c["last_accessed_at"]) < recent
        )
    return old


def docs_build_refs(repo):
    """Return the cache refs of docs builds queued or running, or `None`."""
    refs = set()
    for status in ("queued", "in_progress"):
        try:
            out = gh(
                "api",
                f"repos/{repo}/actions/workflows/docs.yml/runs"
                f"?status={status}&per_page=100",
                "--jq",
                ".workflow_runs[] | [.event, .head_branch, "
                '(.pull_requests | map(.number | tostring) | join(","))] | @tsv',
            )
        except subprocess.CalledProcessError:
            # Unknown builds could be restoring any ref's fits.
            return None
        for line in out.splitlines():
            event, branch, prs = (line.split("\t") + ["", "", ""])[:3]
            if event == "pull_request":
                refs |= {f"refs/pull/{n}/merge" for n in prs.split(",") if n}
            elif branch:
                refs.add(f"refs/heads/{branch}")
    return refs


def plan(
    caches,
    default_ref,
    budget,
    stale_days,
    now,
    running=frozenset(),
    fit_grace_hours=8.0,
    building_refs=None,
):
    """Return the caches to delete as (cache, reason) pairs.

    `running` holds the run ids whose `docs-env` caches are kept.
    `building_refs` holds the refs with a docs build queued or running,
    whose fits are kept; `None`, for unknown, keeps every fit over budget.
    """
    depot = [c for c in caches if DEPOT_KEY.match(c["key"])]
    doomed = {}

    def doom(c, reason):
        if not kept_for_run(c, running):
            doomed[c["id"]] = (c, reason)

    groups = {}
    for c in depot:
        prefix = DEPOT_KEY.match(c["key"])["prefix"]
        groups.setdefault((c["ref"], prefix), []).append(c)
    for group in groups.values():
        group.sort(key=lambda c: timestamp(c["created_at"]), reverse=True)
        for c in group[1:]:
            doom(c, "superseded")

    cutoff = now - timedelta(days=stale_days)
    for c in depot:
        if (
            c["id"] not in doomed
            and c["ref"] == default_ref
            and timestamp(c["last_accessed_at"]) < cutoff
        ):
            doom(c, f"unused for {stale_days} days")

    for c in superseded_fits(caches, now, fit_grace_hours):
        doom(c, "superseded fit")

    total = sum(c["size_in_bytes"] for c in caches if c["id"] not in doomed)
    others = sorted(
        (
            c
            for c in depot
            if c["id"] not in doomed
            and c["ref"] != default_ref
            and not kept_for_run(c, running)
        ),
        key=lambda c: timestamp(c["last_accessed_at"]),
    )
    for c in others:
        if total <= budget:
            break
        doom(c, "over budget")
        total -= c["size_in_bytes"]

    if building_refs is None:
        return list(doomed.values())
    recent = now - timedelta(hours=fit_grace_hours)
    fits = sorted(
        (
            c
            for c in caches
            if FIT_KEY.match(c["key"])
            and c["id"] not in doomed
            and c["ref"] != default_ref
            and c["ref"] not in building_refs
            and timestamp(c["created_at"]) < recent
            and timestamp(c["last_accessed_at"]) < recent
        ),
        key=lambda c: timestamp(c["last_accessed_at"]),
    )
    for c in fits:
        if total <= budget:
            break
        doom(c, "fit over budget")
        total -= c["size_in_bytes"]

    return list(doomed.values())


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--repo", default=os.environ.get("GITHUB_REPOSITORY"))
    parser.add_argument("--default-branch", default="main")
    parser.add_argument("--budget-gb", type=float, default=8.0)
    parser.add_argument("--stale-days", type=int, default=7)
    parser.add_argument("--fit-grace-hours", type=float, default=8.0)
    parser.add_argument("--delete", action="store_true")
    args = parser.parse_args()
    if not args.repo:
        parser.error("--repo is required outside GitHub Actions")

    caches = list_caches(args.repo)
    running = in_flight(args.repo, caches)
    doomed = plan(
        caches,
        f"refs/heads/{args.default_branch}",
        args.budget_gb * GB,
        args.stale_days,
        datetime.now(timezone.utc),
        running,
        args.fit_grace_hours,
        docs_build_refs(args.repo),
    )
    spared = [c for c in caches if kept_for_run(c, running)]

    freed = sum(c["size_in_bytes"] for c, _ in doomed)
    total = sum(c["size_in_bytes"] for c in caches)
    kept = [
        c for c in caches if not (DEPOT_KEY.match(c["key"]) or FIT_KEY.match(c["key"]))
    ]
    lines = [
        f"{len(caches)} caches, {total / GB:.2f} GB; "
        f"{len(kept)} are neither depot nor fit caches and are never deleted "
        f"({sum(c['size_in_bytes'] for c in kept) / GB:.2f} GB).",
        f"Keeping {len(spared)} `docs-env` caches of running builds "
        f"({sum(c['size_in_bytes'] for c in spared) / GB:.2f} GB).",
        f"{'Deleting' if args.delete else 'Would delete'} {len(doomed)}, "
        f"{freed / GB:.2f} GB, leaving {(total - freed) / GB:.2f} GB.",
        "",
    ]
    for c, reason in [(c, "running build") for c in spared] + sorted(
        doomed, key=lambda d: (d[1], d[0]["ref"])
    ):
        lines.append(
            f"- {c['size_in_bytes'] / 2**20:7.0f} MB  {reason:<18} "
            f"{c['ref']}  `{c['key']}`"
        )
    report = "\n".join(lines)
    print(report)
    summary = os.environ.get("GITHUB_STEP_SUMMARY")
    if summary:
        with open(summary, "a") as f:
            f.write("### Depot and fit cache prune\n\n" + report + "\n")

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
