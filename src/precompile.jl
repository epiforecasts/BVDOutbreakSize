# Precompile the expensive first-call work so a fresh process does not pay
# it on its first fit. The dominant cost is not the gradient itself but
# building the Mooncake reverse rule, which happens when the
# `LogDensityFunction` is constructed: measured per component on a 40-day
# grid, rule construction is roughly 87% of a cold build, against a ~38 s
# floor any model pays. Compiling one log-density gradient here bakes those
# rules into the package precompile cache, which CI persists, so the report
# build's fits skip it. A component measured at 126 s of rule construction
# cold builds in 14 s once it is in the cache.
#
# The workload compiles the headline fit itself, through
# [`production_joint`](@ref), rather than a synthetic model shaped like it.
# Mooncake caches a rule against a method signature, and a
# `DynamicPPL.Model`'s type carries the types of its arguments, so a
# workload that differs from the fit in the type of any argument compiles a
# rule the fit never reaches. Building both from `joint_fit_args` removes
# the chance to drift: an earlier hand-written workload passed
# `province_increments = missing` where the fit passes `Matrix{Int}`, and
# left the isolation and treatment histories empty, so those likelihood
# loops never ran and Mooncake, which derives a rule only for code that
# executes, cached nothing for the most expensive stream in the model.
#
# Off by default, since the workload makes package precompilation slow. A
# package preference switches it on, which the report build sets through
# `docs/LocalPreferences.toml`. Enable it elsewhere with
#   using Preferences
#   set_preferences!(BVDOutbreakSize, "precompile_workload" => true)
# Changing the preference triggers one recompilation.

using PrecompileTools: @setup_workload, @compile_workload
using Preferences: @load_preference
using LogDensityProblems: logdensity_and_gradient
using Turing.DynamicPPL: link, VarInfo, getlogjoint, LogDensityFunction

"""
$(TYPEDSIGNATURES)

The model the precompile workload compiles: the headline joint on the real
observations, which is the call `docs/fits/registry.jl` samples.

Separate from the workload block so a test can assert it is the same method
instance the fit builds. Gradient rules are keyed by type, not by data size,
so the grid length does not affect what is cached.
"""
function precompile_workload_joint()
    obs = load_observations()
    return production_joint(obs; breakpoint = default_breakpoint(obs))
end

@static if @load_preference("precompile_workload", false)
    @setup_workload begin
        @compile_workload begin
            ## A rule that cannot be built during precompilation must not
            ## break the package: the Mooncake rule for a
            ## `Distributions.censored` distribution `eval`s into the
            ## `Mooncake` module, which cannot run here, so those rules
            ## compile on the first fit instead. Cache what compiles and
            ## skip the rest. Missing data files are caught here too, so a
            ## checkout without them still precompiles.
            try
                m = precompile_workload_joint()
                vi = link(VarInfo(m), m)
                ldf = LogDensityFunction(
                    m, getlogjoint, vi; adtype = default_adtype()
                )
                logdensity_and_gradient(ldf, collect(vi[:]))
            catch
            end
        end
    end
end
