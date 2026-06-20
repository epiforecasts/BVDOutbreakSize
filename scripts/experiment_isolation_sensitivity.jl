# Sensitivity experiment: does the isolation/treatment-bed stream pin the
# joint cumulative-infection estimate (C_T) low, and which single isolation
# assumption drives it?
#
# Background (see notes/isolation-pinning.md). In recent dev fits the joint
# C_T sits well below what the deaths / cases / confirmed streams imply on
# their own (per-stream table in the released analysis: isolation ~90% CI
# 1548-12869 with its lower bound anchoring the joint 1723-4112, while
# deaths/cases/confirmed all sit at 4800-28000). The isolation occupancy is
# fitted as a length-of-stay survival of the admitted suspect inflow:
#   bed demand ~= p_iso * admissions * (E[LOS] + 1)   (Little's law)
# so a HIGH admission fraction (p_iso ~ Beta(2,2), mean 0.5) times a LONG BVD
# stay (~12 d) makes a fixed observed bed count imply FEW suspects, hence a
# small outbreak. p_iso and LOS are confounded for the occupancy LEVEL, so
# the prior centres do the work and the implied outbreak size rides on them.
#
# This script quantifies that pull two ways:
#   (A) Leave-one-out on the JOINT: fit with isolation included vs excluded
#       (exclude by passing an empty isolation_history) and compare C_T. If
#       isolation pins C_T low, removing it should RAISE the joint C_T.
#   (B) Assumption variants on the ISOLATION SINGLE-STREAM
#       (treatment_only_model): baseline; a lower admission prior
#       (Beta(2,6), mean 0.25); and a shorter BVD length-of-stay. Whichever
#       variant moves the single-stream C_T the most is the dominant driver.
#
# Small NUTS (samples=200, chains=2) — enough to read the medians/intervals,
# not a production fit. Run: julia --project=. scripts/experiment_isolation_sensitivity.jl

using BVDOutbreakSize
using BVDOutbreakSize: treatment_admission_model, isolation_admission_model,
    censored_delay_model, cdf_nmax, lognormal_meansd
using Turing: truncated, Normal, Beta
using Statistics: median
using Printf: @printf

obs = load_observations()
const N = obs.n
const BP = obs.n - obs.who_first_sitrep_days
const EMPTY_HISTORY = (; days = Int[], counts = Int[])

const SAMPLES = 200
const CHAINS = 2

_C(chn) = vec(Array(chn[:C_T]))

function _summarise(label, draws)
    s = posterior_summary(draws)
    @printf("%-34s median=%8.0f   90%% CI [%8.0f, %8.0f]   60%% CI [%8.0f, %8.0f]\n",
        label, median(draws), s.lo90, s.hi90, s.lo60, s.hi60)
    return draws
end

# ---------------------------------------------------------------------------
# (A) Joint leave-one-out: isolation included vs excluded.
# ---------------------------------------------------------------------------
# Mirrors the analysis.jl joint call. Excluding isolation = passing an empty
# isolation_history (and empty bed_capacity_history): the treatment submodel
# then samples from the prior and contributes no likelihood, so its pull on
# C_T disappears while every other stream is unchanged.
function joint_model(; include_isolation::Bool)
    iso = include_isolation ? obs.isolation_history : EMPTY_HISTORY
    cap = include_isolation ? obs.bed_capacity_history : EMPTY_HISTORY
    return bvd_joint(obs.n, obs.exported_cases, obs.total_deaths,
        obs.reported_cases, obs.exports_deaths, obs.confirmed_cases,
        obs.tests_analysed;
        confirmed_deaths = obs.confirmed_deaths,
        recovered_cases = obs.recovered_cases,
        deaths_history = obs.deaths_history,
        reported_history = obs.reported_history,
        confirmed_history = obs.confirmed_history,
        confirmed_deaths_history = obs.confirmed_deaths_history,
        lab_history = obs.lab_history,
        lab_daily_history = obs.lab_daily_history,
        suspected_daily_history = obs.suspected_daily_history,
        suspected_daily_deaths_history = obs.suspected_daily_deaths_history,
        isolation_history = iso,
        bed_capacity_history = cap,
        recovered_history = obs.recovered_history,
        export_case_days = obs.export_case_days,
        export_death_days = obs.export_death_days,
        breakpoint = BP, background_re = true,
        confirmed_positivity_link = :composition,
        genetic = genetic_seeding_model, tmrca_days = obs.tmrca_days)
end

println("=== (A) Joint leave-one-out: isolation IN vs OUT ===")
chn_joint_with = nuts_sample(joint_model(; include_isolation = true);
    samples = SAMPLES, chains = CHAINS, progress = false)
chn_joint_without = nuts_sample(joint_model(; include_isolation = false);
    samples = SAMPLES, chains = CHAINS, progress = false)
C_joint_with = _C(chn_joint_with)
C_joint_without = _C(chn_joint_without)

# ---------------------------------------------------------------------------
# (B) Isolation single-stream under assumption variants.
# ---------------------------------------------------------------------------
# treatment_only_model takes a `treatment` keyword (default
# treatment_admission_model). We inject closures that forward to
# treatment_admission_model with alternate submodels/priors, keeping the call
# signature (isolation_history, bvd_reports_daily, bg_daily, p_drc; kwargs).
#
#  - baseline:        the shipped priors (p_iso ~ Beta(2,2), BVD LOS ~12 d).
#  - low admission:   admission = isolation_admission_model(p_prior=Beta(2,6))
#                     (mean 0.25 instead of 0.5) -> for a fixed bed count,
#                     fewer admitted per suspect implies MORE suspects -> C_T up.
#  - short BVD LOS:    BVD stay re-centred on ~6 d instead of ~12 d -> beds
#                     turn over faster, so a fixed occupancy implies a higher
#                     inflow -> C_T up.
#
# Whichever variant raises the single-stream C_T most is the assumption that
# most pins it low at baseline.

# Shorter BVD length-of-stay delay submodel (mean prior 6 d vs the shipped
# 12 d), matching the truncation style of the shipped default.
short_bvd_los() = censored_delay_model(
    cdf_nmax(lognormal_meansd(6.0, 5.0); q = 0.99);
    mean_prior = truncated(Normal(6.0, 3.0); lower = 1),
    sd_prior = truncated(Normal(5.0, 3.0); lower = 1))

treatment_baseline(args...; kw...) = treatment_admission_model(args...; kw...)

treatment_low_admission(args...; kw...) = treatment_admission_model(args...;
    admission = isolation_admission_model(; p_prior = Beta(2.0, 6.0)), kw...)

treatment_short_los(args...; kw...) = treatment_admission_model(args...;
    bvd_los = short_bvd_los(), kw...)

function treatment_only(; treatment = treatment_admission_model)
    return treatment_only_model(obs.n;
        isolation_history = obs.isolation_history,
        bed_capacity_history = obs.bed_capacity_history,
        breakpoint = BP, treatment = treatment)
end

println("\n=== (B) Isolation single-stream variants ===")
chn_iso_base = nuts_sample(treatment_only(; treatment = treatment_baseline);
    samples = SAMPLES, chains = CHAINS, progress = false)
chn_iso_lowadm = nuts_sample(treatment_only(; treatment = treatment_low_admission);
    samples = SAMPLES, chains = CHAINS, progress = false)
chn_iso_shortlos = nuts_sample(treatment_only(; treatment = treatment_short_los);
    samples = SAMPLES, chains = CHAINS, progress = false)

C_iso_base = _C(chn_iso_base)
C_iso_lowadm = _C(chn_iso_lowadm)
C_iso_shortlos = _C(chn_iso_shortlos)

# ---------------------------------------------------------------------------
# Comparison tables.
# ---------------------------------------------------------------------------
println("\n================  C_T comparison (cumulative infections at cut-off)  ================")
println("\n(A) Joint leave-one-out:")
_summarise("joint, isolation INCLUDED", C_joint_with)
_summarise("joint, isolation EXCLUDED", C_joint_without)
@printf("  -> excluding isolation moves the joint median by %+.0f infections\n",
    median(C_joint_without) - median(C_joint_with))

println("\n(B) Isolation single-stream assumption variants:")
_summarise("baseline (p_iso~Beta(2,2), LOS~12d)", C_iso_base)
_summarise("low admission (p_iso~Beta(2,6))", C_iso_lowadm)
_summarise("short BVD LOS (~6d)", C_iso_shortlos)
@printf("  -> low-admission moves median by %+.0f; short-LOS moves median by %+.0f\n",
    median(C_iso_lowadm) - median(C_iso_base),
    median(C_iso_shortlos) - median(C_iso_base))

# streams_table for a tidy DataFrame view (same helper the report uses).
tbl = streams_table(
    "joint (iso in)" => C_joint_with,
    "joint (iso out)" => C_joint_without,
    "iso baseline" => C_iso_base,
    "iso low-admission" => C_iso_lowadm,
    "iso short-LOS" => C_iso_shortlos)
println("\nstreams_table view:")
show(stdout, tbl)
println()

# Optional figure: overlay the C_T posteriors.
try
    import CairoMakie
    fig = CairoMakie.Figure(size = (900, 500))
    ax = CairoMakie.Axis(fig[1, 1]; xlabel = "cumulative infections C_T",
        ylabel = "density",
        title = "Isolation pull on the outbreak size (C_T)")
    for (label, draws, col) in [
            ("joint (iso in)", C_joint_with, :firebrick),
            ("joint (iso out)", C_joint_without, :seagreen),
            ("iso baseline", C_iso_base, :darkorange),
            ("iso low-admission", C_iso_lowadm, :steelblue),
            ("iso short-LOS", C_iso_shortlos, :purple)]
        CairoMakie.density!(ax, draws; label = label, color = (col, 0.25),
            strokecolor = col, strokewidth = 2)
    end
    CairoMakie.axislegend(ax)
    mkpath("logs")
    CairoMakie.save("logs/isolation_sensitivity_C_T.png", fig)
    println("saved logs/isolation_sensitivity_C_T.png")
catch err
    @warn "figure step skipped" err
end

println("\nDone. Read (A) to see whether isolation pins the joint low, and (B) ",
    "to see which assumption (admission fraction vs BVD length-of-stay) drives it.")
