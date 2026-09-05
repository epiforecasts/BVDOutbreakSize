## Run only this checkout's nowcast test items:
##
##   julia --project=. test/_scoped_nowcast_run.jl
##
## `@run_package_tests` walks every file it can reach, which in a tree with
## sibling worktrees means it discovers test items belonging to other
## checkouts and reports their failures as this one's. Comparing the item's
## filename against absolute paths built from `pwd()` keeps the run inside
## the working copy it was launched from. This file defines no test items of
## its own, so `test/runtests.jl` never picks it up.
using TestItemRunner

targets = [joinpath(pwd(), "test", f)
           for f in ("test_onsets.jl", "test_plots.jl")]
@run_package_tests filter = ti -> string(ti.filename) in targets &&
                                  occursin("nowcast", lowercase(ti.name))
