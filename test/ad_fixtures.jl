"""
    ADFixtures

Shared AD gradient scenarios and backend metadata for BVDOutbreakSize.

This file is included by `test/test_ad_gradients.jl` and by
`benchmark/benchmarks.jl`, so the surface the tests assert is
differentiable is the surface the benchmarks time. It is a plain file
rather than a path package because `[sources]` needs Julia 1.11 and the
LTS test cell runs 1.10, where an unregistered path dependency cannot
resolve. Its imports come from whichever environment includes it, and
both carry every package named below.

A scenario is one model plus a seeded unconstrained point. The units are
per component rather than the full joint: the observation submodels from
`src/models/observations.jl` on a fixed latent draw, the shared latent
process on its own, and the single-stream composers from
`src/models/joint.jl` that attach one observation submodel each to that
latent. The difference between a composer and the latent baseline is the
cost that stream's likelihood adds, which is what guides optimisation.
`bvd_joint` itself is a scenario only on request ([`scenarios`](@ref)): its
gradient is ~14 ms over 76 parameters behind a ~18 min cold compile, too
slow for an unattended run.

The latent inputs each submodel needs (daily onsets, the suspected-case
background and testing fraction) are prior draws taken once at a fixed
seed, so a submodel scenario carries no sampled renewal cost of its own.
"""
module ADFixtures

using BVDOutbreakSize: BVDOutbreakSize, default_adtype, enzyme_adtype,
                       infection_model, onset_incidence_model,
                       patch_infection_model,
                       reported_cases_model, confirmed_cases_model,
                       deaths_model, exports_model, treatment_flow_model,
                       onset_reporting_model, province_composition_model,
                       exports_only_model, deaths_only_model,
                       cases_only_model, confirmed_only_model,
                       treatment_only_model, onsets_only_model, bvd_joint
using LogDensityProblems: logdensity_and_gradient
using Random: seed!
using Turing: DynamicPPL, returned

export Scenario, scenarios, backends, linked_point, gradient_is_finite

"""
One AD scenario: a named model, the group it belongs to and the seed used
to draw its unconstrained evaluation point.

`group` is one of `"Submodel"` (an observation submodel on a fixed latent
draw), `"Latent"` (the shared infection and onset process), `"Composer"`
(a single-stream composer) or `"Joint"` (the full `bvd_joint`).
"""
struct Scenario
    name::String
    group::String
    model::DynamicPPL.Model
    seed::Int
end

## Daily grid every scenario runs on. Long enough for the delay kernels to
## fill and for several observation vintages, short enough that a gradient
## is seconds rather than minutes.
const N = 40

## Fixed seed for the latent prior draw the submodel scenarios share, and
## the default seed for a scenario's evaluation point.
const SEED = 20260518

"""
    linked_point(model; seed = SEED) -> (varinfo, x)

Unconstrained-space `VarInfo` and parameter vector at a seeded prior draw,
the point every gradient is taken at. Mirrors the pattern in
`test/enzyme/runtests.jl`.
"""
function linked_point(model::DynamicPPL.Model; seed::Integer = SEED)
    seed!(seed)
    vi = DynamicPPL.link(DynamicPPL.VarInfo(model), model)
    return vi, collect(vi[:])
end

"""
    log_density_function(scen, adtype) -> (ldf, x)

`LogDensityFunction` for `scen` under `adtype`, with the evaluation point.
"""
function log_density_function(scen::Scenario, adtype)
    vi, x = linked_point(scen.model; seed = scen.seed)
    ldf = DynamicPPL.LogDensityFunction(
        scen.model, DynamicPPL.getlogjoint, vi; adtype = adtype)
    return ldf, x
end

"""
    gradient_is_finite(scen, adtype) -> Bool

Smoke test for one (scenario, backend) pair: take the gradient once and
report whether it came back finite and non-trivial. A backend that throws
returns `false` rather than propagating, so the benchmark suite can omit a
known-broken pair and keep running unattended.
"""
function gradient_is_finite(scen::Scenario, adtype)
    return try
        ldf, x = log_density_function(scen, adtype)
        logp, grad = logdensity_and_gradient(ldf, x)
        isfinite(logp) && length(grad) == length(x) &&
            all(isfinite, grad) && any(!iszero, grad)
    catch
        false
    end
end

"""
    backends()

AD backends to time, as `(name, adtype)` pairs. Mooncake is the package
default and differentiates every model. Enzyme is the opt-in backend from
the package's Enzyme extension and is included only when Enzyme is loaded,
so this module carries no Enzyme dependency of its own and stays usable
from the main test environment, which deliberately excludes it.

Which Enzyme pairs actually register is left to
[`gradient_is_finite`](@ref): Enzyme differentiates the single-stream
composers but not `bvd_joint` (epiforecasts/BVDOutbreakSize#445).
"""
function backends()
    out = [(name = "Mooncake", adtype = default_adtype())]
    try
        push!(out, (name = "Enzyme", adtype = enzyme_adtype()))
    catch
        ## No Enzyme extension loaded: Mooncake alone.
    end
    return out
end

## --- Fixed latent inputs -----------------------------------------------
##
## The observation submodels take a daily onset series and, for the lab and
## occupancy streams, the suspected-case stream's background and testing
## fraction. Drawing these once at a fixed seed keeps a submodel scenario's
## gradient to that submodel's own cost.

"Seeded prior draw of the shared infection and onset process."
function latent_draw(; n::Integer = N, seed::Integer = SEED)
    m = BVDOutbreakSize._latent(n, missing, infection_model,
        onset_incidence_model)
    seed!(seed)
    return returned(m, rand(m))
end

"Seeded prior draw of the suspected-case stream on a given onset series."
function case_draw(onsets; seed::Integer = SEED + 1)
    m = reported_cases_model(
        (; days = Int[], counts = Int[]), missing, onsets, 5.0, 0.3)
    seed!(seed)
    return returned(m, rand(m))
end

## --- Fixed observation histories ---------------------------------------
##
## Cumulative per-vintage counts on the `N`-day grid, shaped like the real
## streams (a handful of vintages, rising counts) but synthetic, so the
## fixtures need no data files and never move when the data update.

const CONFIRMED_HISTORY = (; days = [20, 26, 32, 38], counts = [12, 34, 61, 96])
const LAB_HISTORY = (; days = [20, 26, 32, 38], counts = [60, 150, 280, 430])
const LAB_DAILY_HISTORY = (; days = [34, 36, 38], counts = [22, 27, 31])
const REPORTED_HISTORY = (;
    days = [18, 24, 30, 36], counts = [40, 95, 170, 260])
const SUSPECTED_DAILY_HISTORY = (; days = [30, 33, 36], counts = [14, 17, 21])
const DEATHS_HISTORY = (; days = [22, 28, 34, 40], counts = [5, 14, 27, 41])
const ISOLATION_HISTORY = (;
    days = [28, 30, 32, 34, 36], counts = [26, 33, 38, 42, 45])
const BED_CAPACITY_HISTORY = (;
    days = [28, 30, 32, 34, 36], counts = [60, 60, 80, 80, 80])
const EXPORT_CASE_DAYS = [21, 27, 33]
const ONSET_CURVE_HISTORY = (;
    onset_days = [10, 11, 12, 13, 10, 11, 12, 13, 14],
    report_days = [15, 15, 15, 15, 20, 20, 20, 20, 20],
    prev_report_days = [0, 0, 0, 0, 15, 15, 15, 15, 0],
    increments = [2, 3, 1, 0, 1, 2, 3, 4, 5])
## Three provinces by three confirmed vintages, as
## `province_composition_model` takes them.
const PROVINCE_OBS = [853 21 42; 77 2 5; 3 0 0]
const PROVINCE_MODELLED = [800.0 20.0 40.0; 70.0 2.5 4.0; 2.0 0.1 0.2]
## Three-patch importation kernels: the uncoupled default and a coupled one
## that puts the between-patch term and epsilon on the tape.
const PATCH_KERNEL_OFF = zeros(3, 3)
const PATCH_KERNEL_ON = [0.0 0.0 0.0; 1e-4 0.0 0.0; 1e-5 0.0 0.0]

"""
    scenarios(; n = N, joint = false)

Every scenario, in group order: the shared latent process, the observation
submodels on a fixed latent draw, then the single-stream composers. Pass
`joint = true` to append the full `bvd_joint`, which is excluded by default
for the compile cost given in the module docstring.
"""
function scenarios(; n::Integer = N, joint::Bool = false)
    lat = latent_draw(; n = n)
    onsets = lat.onsets
    infections = lat.infection_state.infections
    cases = case_draw(onsets)

    out = Scenario[]
    push!(out,
        Scenario("latent (infection + onset staging)", "Latent",
            BVDOutbreakSize._latent(n, missing, infection_model,
                onset_incidence_model), SEED))
    ## The multi-patch renewal is its own latent surface: a per-patch Rt walk
    ## and, when the importation kernel is non-zero, the between-patch term
    ## and the implied-national Rt inversion.
    push!(out,
        Scenario("patch_infection_model (uncoupled)", "Latent",
            patch_infection_model(60, 3;
                importation_kernel = PATCH_KERNEL_OFF), SEED))
    push!(out,
        Scenario("patch_infection_model (coupled)", "Latent",
            patch_infection_model(60, 3;
                importation_kernel = PATCH_KERNEL_ON), SEED))

    ## Observation submodels. `k` is the shared surveillance dispersion and
    ## `p_drc` the pooled ascertainment, both fixed here at plausible values
    ## so the submodel samples only its own parameters.
    k, p_drc = 5.0, 0.3
    push!(out,
        Scenario("reported_cases_model", "Submodel",
            reported_cases_model(REPORTED_HISTORY, missing, onsets, k,
                p_drc;
                suspected_daily_history = SUSPECTED_DAILY_HISTORY), SEED))
    push!(out,
        Scenario("confirmed_cases_model", "Submodel",
            confirmed_cases_model(CONFIRMED_HISTORY, missing, onsets, k,
                p_drc, cases.bg_daily, cases.τ_test, cases.bvd_reports_daily;
                lab_history = LAB_HISTORY,
                lab_daily_history = LAB_DAILY_HISTORY), SEED))
    push!(out,
        Scenario("deaths_model", "Submodel",
            deaths_model(DEATHS_HISTORY, missing, onsets, k), SEED))
    push!(out,
        Scenario("exports_model", "Submodel",
            exports_model(missing, infections, 0.02;
                export_case_days = EXPORT_CASE_DAYS,
                incubation_pmf = lat.incubation_pmf), SEED))
    push!(out,
        Scenario("treatment_flow_model", "Submodel",
            treatment_flow_model(ISOLATION_HISTORY, cases.bvd_reports_daily,
                cases.bg_daily, p_drc, 0.3;
                capacity_history = BED_CAPACITY_HISTORY), SEED))
    push!(out,
        Scenario("onset_reporting_model", "Submodel",
            onset_reporting_model(ONSET_CURVE_HISTORY, onsets), SEED))
    push!(out,
        Scenario("province_composition_model", "Submodel",
            province_composition_model(PROVINCE_OBS, PROVINCE_MODELLED), SEED))

    ## Single-stream composers: the latent above plus one observation
    ## submodel each, so a composer minus the latent baseline is that
    ## stream's marginal gradient cost.
    push!(out,
        Scenario("exports_only_model", "Composer",
            exports_only_model(n, missing;
                export_case_days = EXPORT_CASE_DAYS), SEED))
    push!(out,
        Scenario("deaths_only_model", "Composer",
            deaths_only_model(n, missing;
                deaths_history = DEATHS_HISTORY), SEED))
    push!(out,
        Scenario("cases_only_model", "Composer",
            cases_only_model(n, missing;
                reported_history = REPORTED_HISTORY,
                suspected_daily_history = SUSPECTED_DAILY_HISTORY), SEED))
    push!(out,
        Scenario("confirmed_only_model", "Composer",
            confirmed_only_model(n, missing;
                confirmed_history = CONFIRMED_HISTORY,
                lab_history = LAB_HISTORY,
                lab_daily_history = LAB_DAILY_HISTORY), SEED))
    push!(out,
        Scenario("treatment_only_model", "Composer",
            treatment_only_model(n;
                isolation_history = ISOLATION_HISTORY,
                bed_capacity_history = BED_CAPACITY_HISTORY), SEED))
    push!(out,
        Scenario("onsets_only_model", "Composer",
            onsets_only_model(n;
                onset_curve_history = ONSET_CURVE_HISTORY), SEED))

    joint && push!(out,
        Scenario("bvd_joint", "Joint",
            bvd_joint(n, 2, 3, 5, 1, 4, 10; breakpoint = 14), SEED))
    return out
end

end # module
