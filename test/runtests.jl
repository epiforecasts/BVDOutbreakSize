using TestItemRunner

if "downgrade" in ARGS
    # AD-gradient items exercise Mooncake against the downgraded dep
    # set; tolerances drift below the package's pinned versions. The
    # `:slow` items are full NUTS fits that likewise need working AD, so
    # skip them too. Skip quality (Aqua/JET/format/doctest) — those are
    # infra checks that do not change with dep versions.
    @run_package_tests filter = ti -> !(:quality in ti.tags) &&
                                      !(:ad in ti.tags) &&
                                      !(:slow in ti.tags)
elseif "fast" in ARGS
    # Platform-portability cell. The `:slow` NUTS fits and the quality
    # (Aqua/JET/format/doctest) items do not vary by platform and already
    # run on the Linux cells, so a slower runner re-running them buys no
    # signal and spends hours doing it. What is left still loads the
    # package, the data and every model, which is what a platform check is
    # for.
    @run_package_tests filter = ti -> !(:quality in ti.tags) &&
                                      !(:slow in ti.tags)
elseif "skip_quality" in ARGS
    @run_package_tests filter = ti -> !(:quality in ti.tags)
elseif "quality_only" in ARGS
    @run_package_tests filter = ti -> :quality in ti.tags
else
    @run_package_tests
end
