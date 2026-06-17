# Validate ALL new/changed report plots + forecast against the saved
# export-fixed chain, without a full docs build. Builds: vintage PPC with
# confirmed cases; confirmed-deaths + exports scalar PP; Rt trajectory
# (established window, per-draw mask); intervention-effect density;
# clipped joint C_T density; confirmed-inclusive forecast.
using BVDOutbreakSize
using BVDOutbreakSize: confirmed_positivity_windows, knot_days,
                       interpolate_knots, sigmoid_ramp
using Turing: predict
using Serialization: deserialize
using Dates: Day
using Statistics: median, quantile
import CairoMakie

chn = deserialize("logs/joint_chain.jls")
obs = load_observations()
const _BREAKPOINT = obs.n - obs.who_first_sitrep_days
_days_only(h) = (; days = h.days, counts = Int[])

pp = predict(
    bvd_joint(
        obs.n, missing, missing, missing, missing, missing, missing;
        confirmed_deaths = missing,
        deaths_history = _days_only(obs.deaths_history),
        reported_history = _days_only(obs.reported_history),
        suspected_daily_deaths_history =
        _days_only(obs.suspected_daily_deaths_history),
        confirmed_history = obs.confirmed_history,
        confirmed_deaths_history = obs.confirmed_deaths_history,
        lab_history = obs.lab_history,
        lab_daily_history = obs.lab_daily_history,
        export_case_days = obs.export_case_days,
        export_death_days = obs.export_death_days,
        breakpoint = _BREAKPOINT,
        background_re = true),
    chn)

function _vrep(prefix)
    collect(first(pp[k] for k in keys(pp)
    if occursin("$prefix.increments", string(k))))
end
_vdates(days) = string.(obs.seeding .+ Day.(days .- 1))
_srep(name) = vec(Array(first(pp[k] for k in keys(pp)
if occursin(name, string(k)))))

# --- Confirmed-case per-vintage cumulative ---
w = confirmed_positivity_windows(obs.confirmed_history, obs.lab_history,
    obs.lab_daily_history)
win_days = vcat(w.early_days, w.obs_days)
# observed cumulative confirmed at each window end-day
function _conf_at(day)
    i = searchsortedlast(obs.confirmed_history.days, day)
    return i == 0 ? 0 : Int(obs.confirmed_history.counts[i])
end
win_obs = [_conf_at(d) for d in win_days]
early = _vrep("confirmed_state.early_increments")
posv = collect(first(pp[k]
for k in keys(pp)
if occursin("confirmed_state.confirmed_positives.positives", string(k))))
# per-draw concatenated increment vector across windows (early then obs)
conf_reps = [vcat(collect(e), collect(p))
             for (e, p) in zip(vec(early), vec(posv))]
confirmed_panel = (; title = "Confirmed cases",
    dates = _vdates(win_days),
    replicates = conf_reps,
    observed = win_obs, colour = :goldenrod)

reported_panel = (; title = "Suspected cases",
    dates = _vdates(obs.reported_history.days),
    replicates = _vrep("cases_state.reported_increments"),
    observed = obs.reported_history.counts, colour = :steelblue)
deaths_panel = (; title = "Suspected deaths",
    dates = _vdates(obs.deaths_history.days),
    replicates = _vrep("deaths_state.death_increments"),
    observed = obs.deaths_history.counts, colour = :firebrick)
## Daily new suspected deaths (a per-day count, not cumulative): the deaths
## analogue of the new-suspects-per-day stream, picking up where the
## cumulative suspected-death panel freezes on 26 May.
suspected_daily_deaths_panel = (; title = "New suspected deaths/day",
    dates = _vdates(obs.suspected_daily_deaths_history.days),
    replicates = _vrep("deaths_state.suspected_daily_deaths"),
    observed = obs.suspected_daily_deaths_history.counts,
    colour = :indianred, cumulative = false)
tests_panel = (; title = "Specimens analysed",
    dates = _vdates(obs.lab_history.days),
    replicates = _vrep("confirmed_state.analysed_increments"),
    observed = obs.lab_history.counts, colour = :seagreen)

fig1 = plot_vintage_conditional_ppc(
    [reported_panel, confirmed_panel, deaths_panel,
    suspected_daily_deaths_panel, tests_panel])
CairoMakie.save("logs/val_vintage.png", fig1)
println("OK vintage PPC (with confirmed cases): ", "ok")

# --- scalar PP: exports, export-deaths, confirmed deaths ---
## The dated export likelihood scores one Poisson per detection/death day,
## stored as a single per-day count vector `<prefix>.counts`; the scalar
## PPC total sums each replicate's per-day vector across the dated days.
function _dated_total(prefix)
    k = first(k for k in keys(pp)
    if occursin("$prefix.counts", string(k)))
    [sum(v) for v in vec(Array(pp[k]))]
end
pp_exports = _dated_total("exports_state.export_obs")
pp_exports_deaths = _dated_total("exports_deaths_state.death_obs")
pp_conf_deaths = _srep("confirmed_deaths_state.confirmed_deaths")
fig2 = plot_posterior_predictive(
    pp_exports, nothing, obs.exported_cases, nothing;
    pp_exports_deaths = pp_exports_deaths,
    obs_exports_deaths = obs.exports_deaths,
    pp_confirmed_deaths = pp_conf_deaths,
    obs_confirmed_deaths = obs.confirmed_deaths)
CairoMakie.save("logs/val_scalar.png", fig2)
println("OK scalar PP (exports + confirmed deaths): ", "ok")

# --- Rt trajectory (established window) ---
fig3 = plot_rt(chn; n = obs.n, breakpoint = _BREAKPOINT,
    rt_walk_start = _BREAKPOINT,
    as_of_date = string(obs.cutoff), seeding = obs.seeding)
CairoMakie.save("logs/val_rt.png", fig3)
println("OK Rt trajectory: ", "ok")

# --- intervention effect ---
eff = vec(Array(chn[Symbol("rt_state.intervention_effect")]))
println("intervention_effect: median ", round(median(eff); digits = 3),
    " exp ", round(exp(median(eff)); digits = 3))

# --- forecast incl confirmed ---
fc = forecast_reported(chn; horizon = 7,
    obs_cases = obs.reported_cases, obs_deaths = obs.total_deaths,
    obs_confirmed = obs.confirmed_cases,
    obs_confirmed_deaths = obs.confirmed_deaths)
println("forecast cols: ", names(fc))
ft = forecast_table(fc)
println(ft)
fig4 = plot_forecast(fc)
CairoMakie.save("logs/val_forecast.png", fig4)
println("OK forecast plot: ", "ok")

println("ALL VALIDATIONS PASSED")
