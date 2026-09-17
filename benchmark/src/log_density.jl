# One unconstrained log-density evaluation per component, without AD.
#
# This is the denominator for the gradient numbers: a reverse-mode gradient
# costs a small multiple of the forward evaluation, so a component whose
# gradient is slow because its forward pass is slow needs a different fix
# from one whose reverse pass is disproportionate.

using LogDensityProblems: logdensity

SUITE["Log density"] = BenchmarkGroup()

for scen in SCENARIOS
    ## `BenchmarkGroup` has no `get!`, so subgroups are created on first use.
    haskey(SUITE["Log density"], scen.group) ||
        (SUITE["Log density"][scen.group] = BenchmarkGroup())
    ## No `adtype`, so the function evaluates the log-joint and nothing else.
    ldf, x = ADFixtures.log_density_function(scen, nothing)
    SUITE["Log density"][scen.group][scen.name] = @benchmarkable logdensity(
        $ldf, $x
    )
end
