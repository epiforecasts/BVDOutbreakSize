# Gradient benchmarks per component per AD backend.
#
# Each (component, backend) pair is smoke-tested first: the gradient is taken
# once and checked finite and non-trivial (`ADFixtures.gradient_is_finite`).
# A pair that throws or returns a degenerate gradient is omitted rather than
# aborting the run, so a known-broken combination — Enzyme on `bvd_joint`
# (epiforecasts/BVDOutbreakSize#445) is the standing example — leaves the
# rest of the suite reporting. Omitted pairs are listed on stderr so a
# newly-broken one is visible in the run log rather than silently missing.
#
# The gradient is benchmarked through the same `LogDensityFunction` the NUTS
# sampler calls, so the numbers are the sampler's per-leapfrog cost.

using LogDensityProblems: logdensity_and_gradient

SUITE["AD gradients"] = BenchmarkGroup()

for scen in SCENARIOS, entry in BACKENDS
    ## Printed before the smoke test, not after: the test compiles the
    ## gradient, which is the slow step, so a log with no line for a pair is
    ## a pair still compiling rather than a stalled run.
    println(stderr, "[benchmark] smoke test $(scen.name) / $(entry.name)")
    if !ADFixtures.gradient_is_finite(scen, entry.adtype)
        println(stderr,
            "[benchmark] skipping $(scen.name) / $(entry.name): ",
            "no finite gradient")
        continue
    end
    ldf, x = ADFixtures.log_density_function(scen, entry.adtype)
    ## `BenchmarkGroup` has no `get!`, so subgroups are created on first use.
    grp = SUITE["AD gradients"]
    haskey(grp, scen.group) || (grp[scen.group] = BenchmarkGroup())
    haskey(grp[scen.group], scen.name) ||
        (grp[scen.group][scen.name] = BenchmarkGroup())
    bench = @benchmarkable logdensity_and_gradient($ldf, $x)
    grp[scen.group][scen.name][entry.name] = bench
end
