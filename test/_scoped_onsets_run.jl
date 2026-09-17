## Run only this checkout's onset-stream test items:
##
##   julia --project=. test/_scoped_onsets_run.jl
##
## `@run_package_tests` walks every file it can reach, which in a tree with
## sibling worktrees means it discovers test items belonging to other
## checkouts and reports their failures as this one's. Comparing the item's
## filename against an absolute path built from `pwd()` keeps the run inside
## the working copy it was launched from. This file defines no test items of
## its own, so `test/runtests.jl` never picks it up.
using TestItemRunner

target = joinpath(pwd(), "test", "test_onsets.jl")
@run_package_tests filter = ti -> string(ti.filename) == target
