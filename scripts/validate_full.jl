# De-risk the full docs build: short fit of the NEW model (composition,
# late confirmed windows, per-vintage confirmed deaths, 3 June data), then
# run the predict + every report panel/forecast so any predict-key or
# panel error surfaces in ~10 min instead of after a ~40-min build.
# Run: JULIA_NUM_THREADS=2 julia --project=. scripts/validate_full.jl

using BVDOutbreakSize
using Turing: predict
using Dates: Day
using Statistics: quantile

obs = load_observations()
const BP = obs.n - obs.who_first_sitrep_days

function mk(; kw...)
    bvd_joint(obs.n, obs.exported_cases, obs.total_deaths,
        obs.reported_cases, obs.exports_deaths, obs.confirmed_cases,
        obs.tests_analysed; confirmed_deaths = obs.confirmed_deaths,
        deaths_history = obs.deaths_history, reported_history = obs.reported_history,
        confirmed_history = obs.confirmed_history,
        confirmed_deaths_history = obs.confirmed_deaths_history,
        lab_history = obs.lab_history,
        lab_daily_history = obs.lab_daily_history,
        export_case_days = obs.export_case_days,
        export_death_days = obs.export_death_days,
        breakpoint = BP, background_re = true,
        confirmed_positivity_link = :composition,
        genetic = genetic_seeding_model, tmrca_days = obs.tmrca_days, kw...)
end

chn = nuts_sample(mk(); samples = 200, chains = 2, progress = false)
println("fit OK")

_days_only(h) = (; days = h.days, counts = Int[])
pp = predict(
    bvd_joint(obs.n, missing, missing, missing, missing, missing, missing;
        confirmed_deaths = missing,
        deaths_history = _days_only(obs.deaths_history),
        reported_history = _days_only(obs.reported_history),
        confirmed_history = obs.confirmed_history,
        confirmed_deaths_history = _days_only(obs.confirmed_deaths_history),
        lab_history = obs.lab_history,
        lab_daily_history = obs.lab_daily_history,
        export_case_days = obs.export_case_days,
        export_death_days = obs.export_death_days,
        breakpoint = BP, background_re = true,
        confirmed_positivity_link = :composition,
        genetic = genetic_seeding_model, tmrca_days = obs.tmrca_days),
    chn)
println("predict OK")

function _vr(prefix)
    collect(first(k for k in keys(pp)
    if occursin("$prefix.increments", string(k))) |> k -> pp[k])
end
function _vrep(prefix)
    key = first(k for k in keys(pp) if occursin("$prefix.increments", string(k)))
    return collect(pp[key])
end
_dates(days) = string.(obs.seeding .+ Day.(days .- 1))

w = BVDOutbreakSize.confirmed_positivity_windows(obs.confirmed_history,
    obs.lab_history)
cwd = vcat(w.early_days, w.obs_days, w.late_days)
ce = _vrep("confirmed_state.early_increments")
co = collect(first(pp[k]
for k in keys(pp)
if occursin("confirmed_state.confirmed_positives.positives", string(k))))
cl = _vrep("confirmed_state.late_increments")
conf_panel = (; title = "Confirmed cases", dates = _dates(cwd),
    replicates = [vcat(collect(e), collect(p), collect(l))
                  for (e, p, l) in zip(vec(ce), vec(co), vec(cl))],
    observed = [(i = searchsortedlast(obs.confirmed_history.days, d);
                    i == 0 ? 0 : Int(obs.confirmed_history.counts[i]))
                for d in cwd],
    colour = :goldenrod)
panels = [
    (; title = "Suspected cases", dates = _dates(obs.reported_history.days),
        replicates = _vrep("cases_state.reported_increments"),
        observed = obs.reported_history.counts, colour = :steelblue),
    conf_panel,
    (; title = "Suspected deaths", dates = _dates(obs.deaths_history.days),
        replicates = _vrep("deaths_state.death_increments"),
        observed = obs.deaths_history.counts, colour = :firebrick),
    (; title = "Confirmed deaths",
        dates = _dates(obs.confirmed_deaths_history.days),
        replicates = _vrep("confirmed_deaths_state.cdeath_increments"),
        observed = obs.confirmed_deaths_history.counts, colour = :purple),
    (; title = "Specimens analysed",
        dates = _dates(obs.lab_history.days),
        replicates = _vrep("confirmed_state.analysed_increments"),
        observed = obs.lab_history.counts, colour = :seagreen)
]
fig = plot_vintage_conditional_ppc(panels)
println("5-panel vintage PPC built OK (incl. confirmed-case late windows + ",
    "per-vintage confirmed deaths)")

## Rt plot + forecast (the other new pieces).
rtf = plot_rt(chn; n = obs.n, breakpoint = BP,
    rt_start = clamp(obs.n - round(Int, obs.tmrca_days), 1, obs.n),
    rt_walk_start = clamp(BP - RT_WALK_LEAD, 1, obs.n),
    as_of_date = string(obs.cutoff), seeding = obs.seeding,
    ramp = RT_INTERVENTION_RAMP)
println("plot_rt OK")
fc = forecast_reported(chn; horizon = 7)
ft = forecast_table(fc)
println("forecast OK; columns: ", names(ft))
import CairoMakie
CairoMakie.save("logs/validate_full_ppc.png", fig)
println("\nALL REPORT CODE PATHS VALIDATED — safe to build")
