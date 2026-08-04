# BVDOutbreakSize

See [`README.md`](README.md) for what this project estimates and
[`docs/src/contributing.md`](docs/src/contributing.md) for the layout, the
commands and the conventions.
See [`scripts/README.md`](scripts/README.md) before running a script, as
they differ in which Julia project they need.

This file carries only what those do not.

## Writing

Follow the "Analysis report prose" section of
[`docs/src/contributing.md`](docs/src/contributing.md).
It is the source of truth and it is not optional.

Most of this repository was drafted with LLM assistance and the prose has
drifted before.
Follow those rules rather than the surrounding text, which may predate
them.

Development history belongs in `docs/src/news.md` and nowhere else.

## Push early rather than waiting locally

Everything here is slow.
The full test suite takes a long time and a full docs build fits every
model, which has hit GitHub's six-hour ceiling.

Do not sit on a branch waiting for a local full run.
Open the pull request early, push updates to it, and let CI do the long
work while you carry on.

Run `task format` before every push, because a formatting failure wastes a
whole CI run.
Beyond that run only narrow checks: `task test-quick`, a single filtered
TestItem run, or `task docs-main` for the report page alone.
Scope a TestItem filter to the `test/` root, or stale copies in sibling
worktrees are collected too, and read the summary line rather than the
exit code.

Do not run the full suite or a full docs build to check a small change.

## Fits

Fits are cached under `logs/fit_cache`.
A cached fit does not survive a version change in Turing or its
dependencies.
Refit rather than debugging a `KeyError` on a stale chain.

Any change to the model, the priors or the data needs a refit before its
results mean anything.

## Pull requests

Upstream is `epiforecasts/BVDOutbreakSize`.
Pull requests go there, not to the `seabbs` origin remote.
