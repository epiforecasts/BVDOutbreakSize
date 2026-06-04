# De-risk the docs build: validate the predict() + per-vintage PPC panel
# code (the new tests-received panel and background_re-in-predict) against
# the saved reduced chain, without re-fitting. Mirrors the analysis.jl
# block. Run: julia --project=. scripts/validate_ppc.jl

using BVDOutbreakSize
using Turing: predict
using Dates: Day
using Serialization: deserialize

chn = deserialize("logs/joint_chain.jls")
obs = load_observations()
const _BREAKPOINT = obs.n - obs.who_first_sitrep_days
_days_only(h) = (; days = h.days, counts = Int[])

pp_joint = predict(
    bvd_joint(
        obs.n, missing, missing, missing, missing, missing, missing;
        confirmed_deaths = missing,
        deaths_history = _days_only(obs.deaths_history),
        reported_history = _days_only(obs.reported_history),
        confirmed_history = obs.confirmed_history,
        confirmed_deaths_history = obs.confirmed_deaths_history,
        lab_history = obs.lab_history,
        tests_received_history = _days_only(obs.tests_received_history),
        breakpoint = _BREAKPOINT,
        background_re = true),
    chn)
println("predict OK; ", length(keys(pp_joint)), " predicted variables")

function _vintage_replicates(pp, prefix)
    key = first(k for k in keys(pp)
    if occursin("$prefix.increments", string(k)))
    return collect(pp[key])
end
_vintage_dates(days) = string.(obs.seeding .+ Day.(days .- 1))

reported_panel = (; title = "Suspected cases",
    dates = _vintage_dates(obs.reported_history.days),
    replicates = _vintage_replicates(pp_joint, "cases_state.reported_increments"),
    observed = obs.reported_history.counts, colour = :steelblue)
deaths_panel = (; title = "Suspected deaths",
    dates = _vintage_dates(obs.deaths_history.days),
    replicates = _vintage_replicates(pp_joint, "deaths_state.death_increments"),
    observed = obs.deaths_history.counts, colour = :firebrick)
tests_received_panel = (; title = "Specimens received",
    dates = _vintage_dates(obs.tests_received_history.days),
    replicates = _vintage_replicates(
        pp_joint, "confirmed_state.received_increments"),
    observed = obs.tests_received_history.counts, colour = :seagreen)

fig = plot_vintage_conditional_ppc(
    [reported_panel, deaths_panel, tests_received_panel])
println("panel figure built OK: ", typeof(fig))
import CairoMakie
CairoMakie.save("logs/validate_ppc.png", fig)
println("saved logs/validate_ppc.png — all new PPC code paths validated")
