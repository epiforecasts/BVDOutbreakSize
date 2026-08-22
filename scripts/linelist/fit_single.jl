# Refit one single-stream BVDOutbreakSize model, on line-list case counts or on
# the situation reports, and write its reproduction-number trajectory.
#
# The comparison this serves: the same model, the same cut-off, the same
# sampler, fitted to two constructions of the same outbreak. A difference in the
# R_t path is then a difference between the data and nothing else.
#
#   confirmed `confirmed_only_model`, on the cumulative laboratory-confirmed
#             history and the laboratory volumes. It runs the suspected-case
#             stream in predictive mode only (empty history, missing total), so
#             it conditions on laboratory confirmation alone.
#   onsets    `onsets_only_model`, on the onset-by-vintage reporting triangle,
#             which is also laboratory-confirmed cases, by symptom-onset date.
#   cases     `cases_only_model`, on the reported-case history and the daily
#             suspected series. Available but not part of the comparison: the
#             line-list and situation-report reported streams are different
#             case definitions (every alert raised against a triaged suspect
#             count), so a difference between them is not attributable.
#
# Both fits are single-stream, so ascertainment is not identified: level
# differences between the data sources land in `C_T` and `p_drc`, and only shape
# differences reach R_t. `C_T` from these fits is not a size estimate and should
# not be read as one.
#
# Neither fit forecasts. `forecast_reported` cannot read a single-stream chain,
# since it pulls the expected-deaths quantity these fits do not carry, and the
# release publishes no reported-cases forecast to compare against.
#
# Usage, from the repository root:
#
#   LINELIST_INPUT_DIR=<dir> LINELIST_AS_OF=2026-08-10 \
#     julia -t 2 --project=docs scripts/linelist/fit_single.jl cases \
#       --data=linelist_known --delays=repo
#
# See scripts/linelist/README.md for the inputs, the thread count and the
# sampler settings. Progress goes to logs/<fit>.log.
#
# Writes, into the output directory, every name tagged `<fit>_<data>_<delays>`
# so nothing in a grid overwrites anything else:
#
#   linelist_<tag>_stream_estimates.csv   the cut-off and basic reproduction
#       numbers and the final outbreak size, in the shape and column names the
#       release's own stream_estimates.csv uses
#   linelist_<tag>_rt.csv                 the daily R_t trajectory, median and
#       30/60/90% bounds, one row per established grid day
#   linelist_<tag>_rt_draws.csv           thinned per-draw daily R_t, so
#       downstream recomputes its own summaries rather than reusing these
#   linelist_<tag>_delay.csv              onset-to-report delay quantiles, for
#       fits carrying the onset hazard
#   linelist_<tag>_diagnostics.csv        worst R-hat, smallest bulk ESS and
#       the divergent-transition count, so a trajectory is not read without
#       the sampler quality behind it
#   chains/<tag>_<as_of>.jls              the fitted chain, reused on rerun
#
# Under `--pilot` every output name gains a `_pilot`, so a timing run can never
# be mistaken for the real one.

using Pkg: Pkg
Pkg.instantiate()

using BVDOutbreakSize
using CSV
using DataFrames
using Dates
using Serialization: serialize, deserialize
using Statistics: quantile
using TOML

const BOS = dirname(dirname(@__DIR__))
include(joinpath(@__DIR__, "manifest.jl"))
include(joinpath(@__DIR__, "delays.jl"))

const INPUT_DIR = linelist_input_dir()
const OUT_DIR = linelist_output_dir(BOS)

const FITS = ("confirmed", "onsets", "cases")

## Two observation sets, both dated the way the situation reports date theirs:
## a case enters a series when the reporting system first records it as that
## kind of case. An event-date construction (counting each case by the date its
## alert was notified) was tried and dropped: for the confirmed stream it is
## right-truncated by both the notification-to-confirmation interval and the
## lag to the record reaching the export, so its last retained day holds about
## a quarter of the cases that will eventually carry that date, and the model
## has no truncation correction to tell that from a fall in transmission.
const DATA_MODES = ("sitrep", "linelist_known")

flagval(name, default) =
    let a = filter(x -> startswith(x, "--$name="), ARGS)
        isempty(a) ? default : String(split(a[1], "=", limit = 2)[2])
    end

args = filter(a -> !startswith(a, "--"), ARGS)
const PILOT = "--pilot" in ARGS
const REFIT = "--refit" in ARGS
const FIT = isempty(args) ? "cases" : args[1]
const DATA = flagval("data", "linelist_known")
const DELAYS_NAME = flagval("delays", "repo")
const SAMPLES = PILOT ? 10 : 500

## Chain count. The registry's own setting is the default and is what every
## comparison run uses, since sampler settings coming from the registry rather
## than from here is what keeps a refit comparable with the release. `--chains`
## exists for diagnosis: two chains cannot distinguish "this posterior has two
## modes" from "one chain adapted badly and stopped exploring", and more chains
## can. A run that changes it is a diagnostic, not a comparison arm.
const REGISTRY_CHAINS = 2
const CHAINS = let v = flagval("chains", string(REGISTRY_CHAINS))
    n = tryparse(Int, v)
    (isnothing(n) || n < 1) &&
        error("--chains must be a positive integer, got `$v`")
    n
end

## A non-default chain count changes the answer, so it changes the name. Without
## this a four-chain run would find the two-chain chain in the cache and re-export
## it, reporting the old fit under the new settings.
const SUFFIX = (PILOT ? "_pilot" : "") *
               (CHAINS == REGISTRY_CHAINS ? "" : "_c$(CHAINS)")

FIT in FITS || error("unknown fit `$FIT`; expected one of " * join(FITS, ", "))
DATA in DATA_MODES ||
    error("unknown data mode `$DATA`; expected one of " *
          join(DATA_MODES, ", "))

## Every run in a comparison must sit on one grid. The two line-list stream
## constructions end on different days and the released manifest on a third, so
## without one pinned cut-off the R_t series are not on the same axis and the
## comparison is partly a comparison of dates.
const AS_OF = let v = strip(get(ENV, "LINELIST_AS_OF", ""))
    isempty(v) ? nothing : String(v)
end
DATA == "sitrep" && isnothing(AS_OF) &&
    error("set LINELIST_AS_OF for --data=sitrep. The baseline manifest has no " *
          "replacement streams to take a cut-off from, and it only serves as " *
          "a comparator at the same cut-off as the line-list fits.")

const TAG = "$(FIT)_$(DATA)_$(DELAYS_NAME)$(SUFFIX)"

## The delay configuration. `repo` is `nothing`, meaning no overrides, so the
## package defaults stand by construction. A report-delay override handed to the
## onsets fit is refused here rather than silently ignored: that fit estimates
## its reporting delay from the triangle.
const DELAYS = delay_config(DELAYS_NAME)
check_delay_config(DELAYS, FIT)
write_delay_provenance(DELAYS, OUT_DIR)

## The manifest, built here per run rather than shared, so a grid cannot have
## two runs reading each other's observations.
const RELEASED = released_manifest(OUT_DIR)
manifest_path = joinpath(OUT_DIR, "linelist_$(DATA)$(SUFFIX)_observations.toml")

manifest = if DATA == "sitrep"
    write_baseline_manifest(released = RELEASED, as_of = AS_OF,
        out = manifest_path)
else
    write_manifest(released = RELEASED,
        streams = joinpath(INPUT_DIR, "linelist_streams_known.csv"),
        as_of = AS_OF, out = manifest_path,
        source = "DHIS2 case line list, counted at the snapshot that first " *
                 "held each case rather than by notification date, so the " *
                 "series has no ragged edge")
end

## The triangle is read from a fixed filename beside the manifest, so it is
## placed there rather than left wherever the inputs live. Which triangle is
## named, not inferred: the situation-report one is digitised from the published
## epidemic-curve figures and the line-list one is built from recorded onsets,
## and a fit handed the wrong one would be a mixture that nothing downstream
## shows.
const TRIANGLE_SOURCE = DATA == "sitrep" ? :sitrep : :linelist
triangle = place_onset_curve(INPUT_DIR, OUT_DIR;
    source = TRIANGLE_SOURCE, root = BOS)

obs = load_observations(manifest)
och = obs.onset_curve_history

@info "observations" data=DATA delays=DELAYS_NAME cutoff=obs.cutoff n=obs.n triangle_source=triangle.source triangle_vintages=triangle.vintages
@info "onset triangle, as loaded" cells=length(och.increments) vintages=length(och.total_days) last_total=och.last_total
if FIT == "onsets" && isempty(och.increments)
    error("the onset triangle is empty, so `onsets` would fit no data. It is " *
          "read from $(joinpath(dirname(manifest), "onset_curve_scanned.csv")).")
end

## The spec exactly as the released report runs it, taken from this repository's
## registry rather than reconstructed here, so this cannot silently drift from
## the fit it is being compared against. `delays` is the only thing this run
## changes about the model, and `nothing` changes nothing.
include(joinpath(BOS, "docs", "fits", "registry.jl"))

## The intervention breakpoint, pinned across data modes.
##
## `default_breakpoint(obs)` is `obs.n - obs.who_first_sitrep_days`, and
## `who_first_sitrep_days` is `n - reported_history.days[1] + 1` (src/data.jl),
## so the breakpoint is the grid day before the first `reported_case_history`
## vintage. That stream is one of the three the line list replaces, so taking
## the breakpoint from each run's own manifest moves it with the data source
## even though `onsets` never conditions on that stream: on the August data the
## released manifest's first reported vintage is 2026-05-18 (breakpoint 94) and
## the line list's is 2026-05-26 (breakpoint 102), putting the two arms eight
## days apart.
##
## The breakpoint is a fact about the intervention timeline rather than about
## which series is being fitted, and the released manifest's first reported
## vintage is the first WHO joint situation report, where the line list's is
## only the date its export begins. So it is taken from the released manifest
## for every data mode. For `--data=sitrep` that is the value the fit would
## have used anyway, so nothing about that arm changes.
function released_breakpoint(path, obs)
    raw = TOML.parsefile(path)
    haskey(raw, "reported_case_history") || return nothing
    ds = [Date(string(x)) for x in raw["reported_case_history"]["dates"]]
    kept = filter(<=(obs.cutoff), ds)
    isempty(kept) && return nothing
    ## Grid day of a date is `n - (cutoff - date)`, and the breakpoint is the
    ## day before the first vintage, matching `default_breakpoint` exactly.
    return obs.n - Dates.value(obs.cutoff - minimum(kept)) - 1
end

const NATIVE_BREAKPOINT = default_breakpoint(obs)
const BREAKPOINT = let ovr = strip(flagval("breakpoint", ""))
    if !isempty(ovr)
        v = tryparse(Int, ovr)
        isnothing(v) && error("--breakpoint must be an integer, got `$ovr`")
        v
    else
        something(released_breakpoint(RELEASED, obs), NATIVE_BREAKPOINT)
    end
end
@info "breakpoint" breakpoint=BREAKPOINT native_to_this_manifest=NATIVE_BREAKPOINT pinned=(BREAKPOINT != NATIVE_BREAKPOINT)

specs = build_fit_specs(obs; run_sensitivity = false, samples = SAMPLES,
    chains = CHAINS, breakpoint = BREAKPOINT, delays = DELAYS)
spec = specs[findfirst(s -> s.id == FIT, specs)]

## Chain cache. These fits bypass the package fit cache (which is keyed on the
## released data and would miss on every line-list manifest anyway), so they get
## their own: a completed chain is written once and reused, which means the R_t
## export below can be corrected and re-run without paying for the fit again.
## The key carries the cut-off, so a rerun at new data refits rather than
## silently reusing an old grid.
mkpath(joinpath(OUT_DIR, "chains"))
## The breakpoint is part of the key: it changes the model, and a chain fitted
## at a different one must not be reused under the new value. That is the same
## silent-substitution the delay override and the chain count are keyed against.
chain_path = joinpath(OUT_DIR, "chains",
    "$(TAG)_$(isnothing(AS_OF) ? string(obs.cutoff) : AS_OF)_bp$(BREAKPOINT).jls")

chn = if isfile(chain_path) && !REFIT
    @info "chain cache hit, not refitting" chain_path
    repair_chain_keys(deserialize(chain_path))
else
    @info "fitting" id=spec.id cutoff=obs.cutoff n=obs.n samples=SAMPLES chains=CHAINS
    started = now()
    c = spec.thunk()
    @info "fitted" elapsed = canonicalize(round(now() - started, Minute))
    serialize(chain_path, c)
    @info "chain written" chain_path
    c
end

## Fit quality, written before anything derived from the chain. A comparison
## between two data sources is only a comparison if both sides sampled: the
## line-list arms have run consistently more divergent than the situation-report
## arms of the same pair, so a reader needs the R-hat, the bulk ESS and the
## divergence count next to the estimate rather than buried in a run log that is
## not kept beside the outputs. `fit_diagnostics` is the package's own, so these
## numbers mean what the release's diagnostics table means.
diag = fit_diagnostics(chn)
diag_out = joinpath(OUT_DIR, "linelist_$(TAG)_diagnostics.csv")
CSV.write(diag_out,
    DataFrame([(fit = FIT, data = DATA, delays = DELAYS_NAME,
        samples = SAMPLES, chains = CHAINS,
        max_rhat = diag.max_rhat, min_ess_bulk = diag.min_ess_bulk,
        divergences = diag.n_divergent)]))
@info "diagnostics written" diag_out max_rhat=diag.max_rhat min_ess_bulk=diag.min_ess_bulk divergences=diag.n_divergent

## The reproduction-number walk, rebuilt from the parameters every chain carries.
##
## `rt_start` and `rt_walk_start` are 1, as the release uses for every
## single-stream fit (`_fit_rt_draws` and the `stream_fits` table in
## docs/examples/analysis.jl): the alias `R_T` and the joint's later walk start
## live in `bvd_joint`, not in the shared latent submodel these fits compose.
## Matching that convention is the point — the numbers here are read against the
## numbers those produce.
rt = reconstruct_rt(chn; n = obs.n, breakpoint = BREAKPOINT,
    rt_start = 1, rt_walk_start = 1, ramp = RT_INTERVENTION_RAMP)

## Days before a draw's own established window are masked, and the mask differs
## between draws, so each day is summarised over the draws that reach it. This
## is `_rt_quantile` in src/plots.jl, which `plot_rt` uses and which is not
## exported; three lines are reproduced rather than reached for.
rt_at(d) = collect(skipmissing(@view rt[:, d]))
rt_draws = Float64[rt[i, obs.n] for i in axes(rt, 1)]

grid_date(day) = obs.cutoff - Day(obs.n - day)

function rt_trajectory(rt, obs)
    rows = NamedTuple[]
    for d in 1:(obs.n)
        col = rt_at(d)
        ## A day no draw reaches is dropped rather than written as a blank row:
        ## it is outside every draw's established window, not a missing
        ## estimate.
        isempty(col) && continue
        s = posterior_summary(col)
        push!(rows,
            (; fit = FIT, data = DATA, delays = DELAYS_NAME,
                day = d, date = grid_date(d),
                median = quantile(col, 0.5),
                lo30 = s.lo30, hi30 = s.hi30, lo60 = s.lo60, hi60 = s.hi60,
                lo90 = s.lo90, hi90 = s.hi90))
    end
    isempty(rows) &&
        error("no grid day had an established draw, so there is no R_t " *
              "trajectory to write")
    return DataFrame(rows)
end

traj = rt_trajectory(rt, obs)
rt_out = joinpath(OUT_DIR, "linelist_$(TAG)_rt.csv")
CSV.write(rt_out, traj)
@info "Rt trajectory written" rt_out days=nrow(traj) from=minimum(traj.date) to=maximum(traj.date)

## Thinned per-draw trajectory, so a downstream summary is computed from draws
## rather than from these intervals, matching what the release ships in
## `stream_draws.csv`.
const RT_DRAW_THIN = 10
draw_rows = NamedTuple[]
for (j, i) in enumerate(1:RT_DRAW_THIN:size(rt, 1)), d in 1:(obs.n)

    v = rt[i, d]
    ismissing(v) && continue
    push!(draw_rows,
        (; fit = FIT, data = DATA, delays = DELAYS_NAME,
            draw = j, day = d, date = grid_date(d), value = Float64(v)))
end
CSV.write(joinpath(OUT_DIR, "linelist_$(TAG)_rt_draws.csv"),
    DataFrame(draw_rows))

## Final outbreak size, carried by the shared latent submodel and so present in
## every fit. Not identified without a second stream: reported here because the
## release reports it in the same table, not because it is a size estimate.
c_draws = vec(Array(chn[:C_T]))

## Basic reproduction number, the renewal walk's starting value. Probed rather
## than assumed, as `r0_walk_draws` in docs/examples/_setup.jl does: an absent
## key throws rather than reading back empty, and a model built without its own
## walk base should drop the quantity instead of failing the run.
r0_draws = try
    exp.(vec(Array(chn[Symbol("rt_state.log_R0")])))
catch
    nothing
end

function estimate_row(quantity, draws)
    s = posterior_summary(draws)
    return (fit = FIT, data = DATA, delays = DELAYS_NAME, quantity = quantity,
        median = quantile(draws, 0.5),
        lo30 = s.lo30, hi30 = s.hi30, lo60 = s.lo60, hi60 = s.hi60,
        lo90 = s.lo90, hi90 = s.hi90)
end

estimates = DataFrame([estimate_row(q, d)
                       for (q, d) in (("R_T", rt_draws), ("C_T", c_draws),
                               ("R0", r0_draws))
                       if !isnothing(d)])

out = joinpath(OUT_DIR, "linelist_$(TAG)_stream_estimates.csv")
CSV.write(out, estimates)

## Onset-to-report delay, for the fits that estimate one.
##
## `onset_report_G` is the normalised delay distribution, so its quantiles are
## delays in days. The hazard drifts with calendar time, so a delay distribution
## belongs to an onset date rather than to the fit: this reads it at the last day
## of the fitted grid, the most recent position the data support, which is the
## part of the series a nowcast actually corrects.
##
## The triangle is laboratory-confirmed cases nationally with recorded onsets
## only, so a delay estimated from it is not a like-for-like of one estimated
## from a wider case definition or from region-resolved rows with imputed onsets.
## On the line-list triangle it is the same interval `bvd-internal-cmmid` fits
## independently, so the two are worth reading against each other.
function delay_summary(chn, och; probs = (0.5, 0.9))
    isempty(och.increments) && return nothing
    grid_start = minimum(och.onset_days)
    grid_end = maximum(och.report_days)
    haz = try
        reconstruct_onset_hazard(chn; grid_start = grid_start,
            grid_end = grid_end)
    catch
        return nothing
    end
    isempty(haz.logit_h0) && return nothing

    u = grid_end
    D = length(haz.logit_h0[1])
    rows = NamedTuple[]
    for p in probs
        draws = Float64[]
        for i in eachindex(haz.logit_h0)
            g = [onset_report_G(d, haz.logit_h0[i], haz.γ[i], u, grid_start)
                 for d in 0:(D - 1)]
            k = findfirst(>=(p), g)
            ## No delay in the hazard's support reaches `p`: the support's own
            ## edge is the honest answer, not an extrapolation past it.
            push!(draws, isnothing(k) ? Float64(D - 1) : Float64(k - 1))
        end
        s = posterior_summary(draws)
        push!(rows,
            (fit = FIT, data = DATA, delays = DELAYS_NAME,
                quantity = "delay_p$(Int(round(100 * p)))",
                median = quantile(draws, 0.5),
                lo30 = s.lo30, hi30 = s.hi30, lo60 = s.lo60, hi60 = s.hi60,
                lo90 = s.lo90, hi90 = s.hi90))
    end
    return DataFrame(rows)
end

delays_df = delay_summary(chn, och)
if !isnothing(delays_df)
    delay_out = joinpath(OUT_DIR, "linelist_$(TAG)_delay.csv")
    CSV.write(delay_out, delays_df)
    @info "delay written" delay_out
    println(delays_df)
end

@info "written" out samples = SAMPLES
println(estimates)
