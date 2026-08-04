# Precompile the expensive first-call work so a fresh process does not pay it
# on its first fit. The dominant cost is Mooncake building the reverse rule
# for `bvd_joint` the first time it differentiates the model (minutes).
# Compiling a single log-density gradient of a synthetic full-stream joint
# here bakes those rules into the package precompile cache. CI persists that
# cache (julia-actions/cache), so a build skips the compile. The Mooncake
# rules are keyed by method signature, not data size, so a 40-day synthetic
# fit compiles what the full-size fits reuse. The synthetic joint mirrors
# the report's headline call, so every precompilable stream's rule is
# covered. The gradient path (not a NUTS fit) is the exact expensive
# operation and avoids sampler-adaptation fragility at tiny sizes.
#
# Opt-in: the workload makes package precompilation slow (it runs the joint
# Mooncake compile), so it is off by default, for quick local iteration.
# A package preference switches it on. The report build enables it through
# `docs/LocalPreferences.toml`. Enable it for another environment with
#   using Preferences
#   set_preferences!(BVDOutbreakSize, "precompile_workload" => true)
# Changing the preference triggers one recompilation.

using PrecompileTools: @setup_workload, @compile_workload
using Preferences: @load_preference
using LogDensityProblems: logdensity_and_gradient
using Turing.DynamicPPL: link, VarInfo, getlogjoint, LogDensityFunction

@static if @load_preference("precompile_workload", false)
    @setup_workload begin
        ## Small synthetic data for the streams the report's joint scores, so
        ## the reverse rule is compiled (days are 1-based into an n = 40 grid,
        ## ascending; counts positive). The isolation/bed stream and the
        ## genetic-seeding TMRCA are left out on purpose: both score a
        ## `Distributions.censored` distribution whose Mooncake rule `eval`s
        ## into the `Mooncake` module, and that cannot run during
        ## precompilation. Those two rules compile on the report's first fit
        ## instead. Everything else is cached here.
        dh = (; days = [13, 18, 40], counts = [10, 14, 18])
        rh = (; days = [13, 18, 40], counts = [340, 516, 905])
        ch = (; days = [13, 18, 40], counts = [9, 17, 27])
        cdh = (; days = [18, 40], counts = [2, 5])
        lh = (; days = [18, 40], counts = [30, 50])
        ldh = (; days = [18, 40], counts = [5, 8])
        sdh = (; days = [18, 40], counts = [20, 30])
        sddh = (; days = [18, 40], counts = [2, 3])
        rch = (; days = [18, 40], counts = [4, 7])
        ## Onset-reporting-triangle stream: a couple of synthetic increment
        ## cells so `increments[i] ~ safe_studentt(...)` in
        ## `onset_reporting_model` actually runs (an empty history makes
        ## that loop a no-op, and Mooncake never builds a reverse rule for a
        ## branch it never executes; see `onset_reporting_model`).
        och = (; onset_days = [10, 12], report_days = [18, 40],
            prev_report_days = [0, 18], increments = [4, 6])
        @compile_workload begin
            ## Mirror the headline `bvd_joint` call for the precompilable
            ## streams, differentiated once under the default (Mooncake)
            ## backend so the report's first joint fit reuses the cached rule.
            m = bvd_joint(40, 2, 18, 905, 0, 27, 50;
                confirmed_deaths = 5,
                recovered_cases = 12,
                deaths_history = dh,
                reported_history = rh,
                confirmed_history = ch,
                confirmed_deaths_history = cdh,
                lab_history = lh,
                lab_daily_history = ldh,
                suspected_daily_history = sdh,
                suspected_daily_deaths_history = sddh,
                recovered_history = rch,
                export_case_days = [20, 30],
                export_death_days = [35],
                onset_curve_history = och,
                breakpoint = 30,
                background_re = true,
                confirmed_positivity_link = :composition)
            ## A rule that cannot be built during precompilation (the
            ## `eval`-into-`Mooncake` barrier) must not break the package.
            ## Cache what compiles and skip the rest.
            try
                vi = link(VarInfo(m), m)
                ldf = LogDensityFunction(
                    m, getlogjoint, vi; adtype = default_adtype())
                logdensity_and_gradient(ldf, collect(vi[:]))
            catch
            end
        end
    end
end
