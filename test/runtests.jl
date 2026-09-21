using TestItemRunner

## The suite is split across CI jobs by tag so no single job carries every
## cost. `:quality` (Aqua/JET/format/doctest) and `:ad` (Mooncake gradients,
## ~19 min) do not vary by platform or Julia version, so each runs once in a
## job of its own rather than on all four matrix cells. See .github/workflows/
## test.yml.
if "downgrade" in ARGS
    # AD-gradient items exercise Mooncake against the downgraded dep
    # set; tolerances drift below the package's pinned versions. The
    # `:slow` items are full NUTS fits that likewise need working AD, so
    # skip them too.
    @run_package_tests filter = ti -> !(:quality in ti.tags) &&
        !(:ad in ti.tags) &&
        !(:slow in ti.tags)
elseif "fast" in ARGS
    # Platform-portability cell. The `:slow` NUTS fits do not vary by
    # platform and already run on the Linux cells, so a slower runner
    # re-running them buys no signal and spends hours doing it. What is
    # left still loads the package, the data and every model, which is what
    # a platform check is for.
    @run_package_tests filter = ti -> !(:quality in ti.tags) &&
        !(:ad in ti.tags) &&
        !(:slow in ti.tags)
elseif "skip_quality" in ARGS
    @run_package_tests filter = ti -> !(:quality in ti.tags) &&
        !(:ad in ti.tags)
elseif "quality_only" in ARGS
    @run_package_tests filter = ti -> :quality in ti.tags
elseif "ad_only" in ARGS
    @run_package_tests filter = ti -> :ad in ti.tags
else
    @run_package_tests
end
