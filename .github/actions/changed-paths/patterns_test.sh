#!/usr/bin/env bash
# Check the path patterns the workflows gate on against a table of realistic
# changes. The patterns are read out of the workflow files themselves, so
# this tests what CI actually runs rather than a second copy of it.
#
# A wrong pattern is silent in both directions and neither shows up as a
# failure: too narrow and a change skips the work that would have caught it,
# too wide and the gate saves nothing. Hence a table rather than a review.
set -euo pipefail
cd "$(dirname "$0")/../../.."

pattern_for() { # workflow file -> the `patterns:` its `changes` job passes
    python3 -c '
import sys, yaml
d = yaml.safe_load(open(sys.argv[1]))
for step in d["jobs"]["changes"]["steps"]:
    if step.get("id") == "filter":
        print(step["with"]["patterns"].strip())
        break
else:
    raise SystemExit(f"{sys.argv[1]}: no changed-paths filter step")
' "$1"
}

REPORT=$(pattern_for .github/workflows/docs.yml)
TESTS=$(pattern_for .github/workflows/test.yml)
COVERAGE=$(pattern_for .github/workflows/codecoverage.yaml)

if [ "$TESTS" != "$COVERAGE" ]; then
    echo "coverage gates on different paths from the test suite it measures:"
    echo "  tests:    $TESTS"
    echo "  coverage: $COVERAGE"
    exit 1
fi

failures=0
check() { # description, newline-separated paths, expect report, expect tests
    local name="$1" files="$2" want_report="$3" want_tests="$4" report tests
    report=$(grep -qE "$REPORT" <<<"$files" && echo true || echo false)
    tests=$(grep -qE "$TESTS" <<<"$files" && echo true || echo false)
    if [ "$report" = "$want_report" ] && [ "$tests" = "$want_tests" ]; then
        printf '  ok    %-34s report=%-5s tests=%s\n' "$name" "$report" "$tests"
    else
        printf '  FAIL  %-34s report=%-5s tests=%-5s (wanted %s / %s)\n' \
            "$name" "$report" "$tests" "$want_report" "$want_tests"
        failures=$((failures + 1))
    fi
}

# The report is built from the model, the data and the pages. The tests are
# the package plus the few build files test items include directly.
check "model change"           "src/models/joint.jl"                       true  true
check "new data vintage"       "data/observations.toml"                    true  true
check "package source"         "src/scoring.jl"                            true  true
# `ext/` is split. The Enzyme extension is never loaded by the docs build,
# since Enzyme is not in the docs environment, so it must not rebuild the
# report. The TensorBoardLogger one is loaded, so it must.
check "the Enzyme extension"   "ext/BVDOutbreakSizeEnzymeExt.jl"           false true
check "the logging extension"  "ext/BVDOutbreakSizeTensorBoardLoggerExt.jl" true true
check "Project.toml"           "Project.toml"                              true  true
check "test only"              "test/test_renewal.jl"                      false true
check "analysis page"          "docs/pages/estimates/national.jl"          true  false
check "news and contributing"  $'docs/src/news.md\ndocs/src/contributing.md' true false
check "README"                 "README.md"                                 true  false
# Test items include cache.jl, registry.jl, summary.jl and score_releases.jl
# directly, so these reach the tests as well as the report.
check "a fit-registry helper"  "docs/fits/registry.jl"                     true  true
check "a build script"         "scripts/score_releases.jl"                 true  true
# `scripts/` is split the same way: only the two the docs workflow runs
# rebuild the report. The SitRep downloader, the scanners and the backfill
# driver do not.
check "the SitRep downloader"  "scripts/download_sitreps.jl"               false true
check "benchmarks"             "benchmark/run.jl"                          false false
check "Taskfile"               "Taskfile.yml"                              false false
check "agent instructions"     "AGENTS.md"                                 false false
check "citation metadata"      "CITATION.cff"                              false false
check "an unrelated workflow"  ".github/workflows/release.yml"             false false
check "the docs workflow"      ".github/workflows/docs.yml"                true  false
check "the docs pkgimage action" ".github/actions/bvd-pkgimage/action.yaml" true false
check "the cache prune"        $'.github/workflows/cache-prune.yml\n.github/scripts/prune_depot_caches.py' false false
check "the test workflow"      ".github/workflows/test.yml"                false true
check "the coverage workflow"  ".github/workflows/codecoverage.yaml"       false true
check "this gate itself"       ".github/actions/changed-paths/action.yaml" true  true

if [ "$failures" -gt 0 ]; then
    echo "$failures path-gating case(s) wrong" >&2
    exit 1
fi
echo "all path-gating cases behave as intended"
