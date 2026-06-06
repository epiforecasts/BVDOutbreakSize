# Exploratory harness: can Enzyme reverse-mode differentiate the FULL
# renewal `bvd_joint` (all six observed streams + genetic seeding), and
# does its gradient match Mooncake?
#
# Run from this isolated Enzyme sub-environment:
#   julia --project=test/enzyme test/enzyme/joint_explore.jl
#
# The joint is built straight from `load_observations()`, mirroring the
# `chn_joint` call in docs/examples/analysis.jl, so the observed-data
# arguments are exactly those used in the real fit. We take a single
# unconstrained-space gradient under each backend at a fixed prior link
# point and compare elementwise. A single `logdensity_and_gradient` (not
# NUTS) isolates compile/differentiate viability from trajectory
# exploration into extreme regions.
#
# STATUS as of renewal HEAD 7dc5620 (composition-link DEFAULT positivity
# with the assay-specificity term `p = s·q + (1 − spec)(1 − q)`, the dated
# per-day Uganda export likelihood, and the genetic seeding bound; full
# real-data joint 62-dim): Enzyme DIFFERENTIATES the joint and matches
# Mooncake after the de-box fix on this branch. Three conditional/anonymous
# closures that Enzyme reverse mode could not handle were made
# type-stable (pure refactor, Mooncake bit-identical):
#   1. `bvd_joint`'s `if background_re … else … end` boxed the
#      `case_bg_re`/`death_bg_re` closure capture in a `Base.RefValue`;
#   2. the `:composition` positivity `map(do i)` in `confirmed_cases_model`
#      was a second boxed closure, extracted to `composition_positivity`
#      (which now also applies the sensitivity/specificity transform);
#   3. that helper's `map(do i)` itself tripped Enzyme's shadow handling,
#      so it was rewritten as an explicit loop (no closure).
# This harness asserts the Enzyme/Mooncake gradient match on the current
# joint: logp identical, single gradient maxabs ≈ 2.7e-12 / relerr ≈ 2e-15.
# See PR #201. The fix lives on this branch as a candidate for renewal.
#
# Note: a further renewal change (the m-induced-T seeding) is still pending
# and will need another sync of this branch once it lands.

using Test
using Enzyme
using Mooncake
using Turing: DynamicPPL
using LogDensityProblems: logdensity, logdensity_and_gradient
using Random: seed!
using BVDOutbreakSize: default_adtype, enzyme_adtype, bvd_joint,
                       load_observations, genetic_seeding_model

obs = load_observations()
breakpoint = obs.n - obs.who_first_sitrep_days

# Full joint with the real observed streams, identical to analysis.jl.
function build_joint()
    bvd_joint(
        obs.n, obs.exported_cases, obs.total_deaths,
        obs.reported_cases, obs.exports_deaths, obs.confirmed_cases,
        obs.tests_analysed;
        confirmed_deaths = obs.confirmed_deaths,
        deaths_history = obs.deaths_history,
        reported_history = obs.reported_history,
        confirmed_history = obs.confirmed_history,
        confirmed_deaths_history = obs.confirmed_deaths_history,
        lab_history = obs.lab_history,
        tests_received_history = obs.tests_received_history,
        breakpoint = breakpoint,
        background_re = true,
        confirmed_positivity_link = :composition,
        genetic = genetic_seeding_model,
        tmrca_days = obs.tmrca_days)
end

model = build_joint()
seed!(20260518)
vi = DynamicPPL.link(DynamicPPL.VarInfo(model), model)
x0 = collect(vi[:])
@info "joint dimension" dim = length(x0)

function ldf(adtype)
    DynamicPPL.LogDensityFunction(
        model, DynamicPPL.getlogjoint, vi; adtype = adtype)
end

# Wall-clock for a single gradient (compile excluded via a warm call).
function timed_grad(adtype)
    l = ldf(adtype)
    v, g = logdensity_and_gradient(l, x0)   # compile
    t = @elapsed (_, g) = logdensity_and_gradient(l, x0)
    return v, g, t
end

println("=== Mooncake ===")
v_mc, g_mc, t_mc = timed_grad(default_adtype())
@info "Mooncake" logp=v_mc seconds=t_mc

println("=== Enzyme ===")
enzyme_ok = true
v_en = NaN
g_en = similar(g_mc)
t_en = NaN
try
    global v_en, g_en, t_en = timed_grad(enzyme_adtype())
    @info "Enzyme" logp=v_en seconds=t_en
catch err
    global enzyme_ok = false
    println("ENZYME FAILED:")
    showerror(stdout, err)
    println()
    for (i, fr) in enumerate(stacktrace(catch_backtrace()))
        i > 25 && break
        println("  ", fr)
    end
end

# Load-robust timing: alternate backends so background drift hits both
# equally (a single before/after pair previously misread a load artifact
# as a 2x speedup). Report the median per-gradient wall-clock.
if enzyme_ok
    l_mc = ldf(default_adtype())
    l_en = ldf(enzyme_adtype())
    logdensity_and_gradient(l_mc, x0)   # warm
    logdensity_and_gradient(l_en, x0)
    ts_mc = Float64[]
    ts_en = Float64[]
    for _ in 1:25
        push!(ts_mc, @elapsed logdensity_and_gradient(l_mc, x0))
        push!(ts_en, @elapsed logdensity_and_gradient(l_en, x0))
    end
    sort!(ts_mc)
    sort!(ts_en)
    med_mc = ts_mc[13]
    med_en = ts_en[13]
    speedup = med_mc / med_en
    @info "interleaved median single-gradient" mooncake_s=med_mc enzyme_s=med_en
    @info "  speedup" speedup
end

@testset "Enzyme differentiates the full joint" begin
    @test enzyme_ok
    if enzyme_ok
        maxabs = maximum(abs.(g_en .- g_mc))
        relerr = maxabs / max(1.0, maximum(abs.(g_mc)))
        @info "gradient comparison" maxabs relerr
        @test isapprox(g_en, g_mc; rtol = 1e-6)
    end
end
