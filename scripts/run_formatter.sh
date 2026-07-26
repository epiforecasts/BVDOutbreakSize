#!/usr/bin/env bash
# Run JuliaFormatter (SciML style) over src/, test/, docs/, scripts/
# using the project's isolated test/formatter/ sub-environment. Invoked
# by the local pre-commit hook and re-usable from the command line.
set -euo pipefail

cd "$(dirname "$0")/.."
# The registry is refreshed before instantiating so this hook and
# test/package/CodeFormatting.jl always resolve the SAME JuliaFormatter: the
# sub-environment has no committed Manifest, and `Pkg.instantiate` resolves
# against whatever registry snapshot the depot already holds. It also keeps the
# exact pin in test/formatter/Project.toml resolvable on a depot last updated
# before that version was registered.
julia --project=test/formatter -e '
using Pkg
Pkg.Registry.update()
Pkg.instantiate()
using JuliaFormatter
dirs = ["src", "test", "docs", "scripts"]
# `map`, not `all`: `all` short-circuits, so an unformatted file in an early
# directory left every later one unformatted AND unreported, one round trip per
# directory. Every directory is rewritten in a single pass.
results = map(d -> JuliaFormatter.format(d; overwrite = true), dirs)
exit(all(results) ? 0 : 1)'
