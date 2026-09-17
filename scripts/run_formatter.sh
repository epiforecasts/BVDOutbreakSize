#!/usr/bin/env bash
# Run Runic over src/, test/, docs/, scripts/, benchmark/ and ext/ using the
# project's isolated test/formatter/ sub-environment. Called by `task format`
# and re-usable from the command line.
#
# The pre-commit hook does not come through here. It builds its own
# environment from the Runic version in .pre-commit-config.yaml, which
# test/package/CodeFormatting.jl checks against the pin in
# test/formatter/Project.toml.
set -euo pipefail

cd "$(dirname "$0")/.."
# The registry is refreshed before instantiating so this script and
# test/package/CodeFormatting.jl always resolve the same Runic: the
# sub-environment has no committed Manifest, and `Pkg.instantiate` resolves
# against whatever registry snapshot the depot already holds. It also keeps
# the exact pin resolvable on a depot last updated before that version was
# registered.
julia --project=test/formatter -e '
using Pkg
Pkg.Registry.update()
Pkg.instantiate()
using Runic
dirs = filter(isdir, ["src", "test", "docs", "scripts", "benchmark", "ext"])
exit(Runic.main(["--inplace", dirs...]))'
