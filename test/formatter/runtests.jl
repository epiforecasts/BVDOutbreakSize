#!/usr/bin/env julia
# Runic format check over the source trees, run from this isolated
# environment so the pin above is the only Runic in scope.
#
# `Runic.main` is called directly rather than through `julia -m Runic`,
# which needs Julia >= 1.12 while CI also runs the LTS floor.
#
#   julia --project=test/formatter test/formatter/runtests.jl

using Runic

project_root = dirname(dirname(@__DIR__))
dirs = ("src", "test", "docs", "scripts", "benchmark", "ext")
dirs_to_check = filter(isdir, [joinpath(project_root, d) for d in dirs])

if isempty(dirs_to_check)
    println(stderr, "no source directories found to format-check")
    exit(0)
end

# `--check` alone prints nothing, so a red run says only that something is
# wrong. `--verbose` names every file it visits and marks the offenders, and
# `--diff` prints what would change, which makes a red run actionable from
# the log alone.
rc = Runic.main(["--check", "--verbose", "--diff", dirs_to_check...])

if rc != 0
    println(stderr, "Run scripts/run_formatter.sh to reformat in place.")
end

exit(rc)
