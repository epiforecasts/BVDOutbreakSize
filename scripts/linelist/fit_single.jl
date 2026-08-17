# Refit one single-stream BVDOutbreakSize model on line-list case counts.
#
# The joint refit (scripts/linelist/fit_joint.jl) is partial: most of the
# manifest's streams have no line-list source, so its result mixes line-list
# cases with situation-report deaths, laboratory volumes, beds and exports. The
# single-stream fits have no such mixture. Each conditions on one stream, and for
# the two run here the line list supplies every observation that stream needs, so
# the comparison against the release's own fit of the same id is a difference
# between data sources and nothing else.
#
#   cases    `cases_only_model`, on the reported-case history and the daily
#            suspected series, both from linelist_streams.csv.
#   onsets   `onsets_only_model`, on the onset-by-vintage reporting triangle
#            alone, from onset_curve_scanned.csv.
#
# Neither fit forecasts. `forecast_reported` cannot read a single-stream chain,
# since it pulls the expected-deaths quantity these fits do not carry, and the
# release publishes no reported-cases forecast to compare against.
#
# Usage, from the repository root:
#
#   LINELIST_INPUT_DIR=<dir> julia -t 2 --project=docs \
#       scripts/linelist/fit_single.jl cases --pilot   # 10 samples
#   LINELIST_INPUT_DIR=<dir> julia -t 2 --project=docs \
#       scripts/linelist/fit_single.jl cases           # 500, 2 chains
#
# See scripts/linelist/README.md for the inputs, the thread count and the
# sampler settings. Progress goes to logs/<fit>.log.
#
# Writes, into the output directory:
#
#   linelist_<fit>_stream_estimates.csv   the cut-off and basic reproduction
#       numbers and the final outbreak size, in the shape and column names the
#       release's own stream_estimates.csv uses, so the two are read side by
#       side without translation
#   linelist_<fit>_delay.csv              onset-to-report delay quantiles, for
#       fits carrying the onset hazard
#
# Under `--pilot` every output name gains a `_pilot`, so a timing run can never
# be mistaken for the real one.

using Pkg: Pkg
Pkg.instantiate()

using BVDOutbreakSize
using CSV
using DataFrames
using Dates
using Statistics: quantile

const BOS = dirname(dirname(@__DIR__))
include(joinpath(@__DIR__, "manifest.jl"))

const INPUT_DIR = linelist_input_dir()
const OUT_DIR = linelist_output_dir(BOS)

const FITS = ("cases", "onsets")

args = filter(a -> !startswith(a, "--"), ARGS)
const PILOT = "--pilot" in ARGS
const FIT = isempty(args) ? "cases" : args[1]
const SAMPLES = PILOT ? 10 : 500
const SUFFIX = PILOT ? "_pilot" : ""

FIT in FITS || error("unknown fit `$FIT`; expected one of " * join(FITS, ", "))

## `--streams=PATH` fits an alternative construction of the case streams, which
## is how two indexings are compared. The tag in the streams filename carries
## through to the manifest and every output, so two constructions never overwrite
## each other's results.
##
##   linelist_streams.csv            no tag, counted by notification date
##   linelist_streams_known.csv      "known", counted at the snapshot that first
##                                   held each case
const STREAMS = let a = filter(x -> startswith(x, "--streams="), ARGS)
    isempty(a) ? nothing : abspath(split(a[1], "=", limit = 2)[2])
end
const TAG = if isnothing(STREAMS)
    ""
else
    m = match(r"^linelist_streams_(\w+)\.csv$", basename(STREAMS))
    isnothing(m) ? "" : "_" * m.captures[1]
end

## The manifest. Without `--streams` this is the one the joint refit wrote,
## reused rather than rebuilt so the fits share byte-identical observations.
## With it, the manifest is built here from the same shared builder the joint
## refit uses, under a name carrying both the tag and the pilot suffix so a
## ten-sample run cannot overwrite the manifest a real fit was run against.
manifest = if isnothing(STREAMS)
    joinpath(OUT_DIR, "linelist_observations.toml")
else
    write_manifest(
        released = released_manifest(OUT_DIR),
        streams = STREAMS,
        out = joinpath(OUT_DIR,
            "linelist$(TAG)$(SUFFIX)_observations.toml"),
        source = "DHIS2 case line list, counted at the snapshot that first " *
                 "held each case rather than by notification date, so the " *
                 "series has no ragged edge")
end
isfile(manifest) ||
    error("missing $manifest. Run scripts/linelist/fit_joint.jl first, or " *
          "pass --streams=PATH to build a manifest here.")

## The triangle is read from a fixed filename beside the manifest, so it is
## placed there rather than left wherever the inputs live.
place_onset_curve(INPUT_DIR, OUT_DIR)

obs = load_observations(manifest)
och = obs.onset_curve_history

@info "onset triangle" cells=length(och.increments) vintages=length(och.total_days) last_total=och.last_total
if FIT == "onsets" && isempty(och.increments)
    error("the onset triangle is empty, so `onsets` would fit no data. It is " *
          "read from $(joinpath(dirname(manifest), "onset_curve_scanned.csv")).")
end

## The spec exactly as the released report runs it, taken from this repository's
## registry rather than reconstructed here, so this cannot silently drift from
## the fit it is being compared against.
include(joinpath(BOS, "docs", "fits", "registry.jl"))
specs = build_fit_specs(obs; run_sensitivity = false, samples = SAMPLES)
spec = specs[findfirst(s -> s.id == FIT, specs)]

@info "fitting" id=spec.id cutoff=obs.cutoff n=obs.n samples=SAMPLES
started = now()
chn = spec.thunk()
@info "fitted" elapsed = canonicalize(round(now() - started, Minute))

## Cut-off reproduction number. The joint exposes `R_T` directly; the
## single-stream composers do not, since that alias lives in `bvd_joint` rather
## than the shared latent submodel, so it is rebuilt from the walk parameters
## every chain carries and read at the cut-off, the last day of the
## reconstructed path. `rt_start` and `rt_walk_start` are 1, as the release uses
## for every single-stream fit. This mirrors `_fit_rt_draws` in
## docs/examples/analysis.jl; matching it is the point, since the number this
## produces is compared against the number that produces.
rt = reconstruct_rt(chn; n = obs.n, breakpoint = default_breakpoint(obs),
    rt_start = 1, rt_walk_start = 1, ramp = RT_INTERVENTION_RAMP)
rt_draws = Float64[rt[i, obs.n] for i in axes(rt, 1)]

## Final outbreak size, carried by the shared latent submodel and so present in
## every fit.
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
    return (fit = FIT, quantity = quantity,
        median = quantile(draws, 0.5),
        lo30 = s.lo30, hi30 = s.hi30, lo60 = s.lo60, hi60 = s.hi60,
        lo90 = s.lo90, hi90 = s.hi90)
end

estimates = DataFrame([estimate_row(q, d)
                       for (q, d) in (("R_T", rt_draws), ("C_T", c_draws),
                               ("R0", r0_draws))
                       if !isnothing(d)])

out = joinpath(OUT_DIR,
    "linelist_$(FIT)$(TAG)$(SUFFIX)_stream_estimates.csv")
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
            (fit = FIT,
                quantity = "delay_p$(Int(round(100 * p)))",
                median = quantile(draws, 0.5),
                lo30 = s.lo30, hi30 = s.hi30, lo60 = s.lo60, hi60 = s.hi60,
                lo90 = s.lo90, hi90 = s.hi90))
    end
    return DataFrame(rows)
end

delays = delay_summary(chn, och)
if !isnothing(delays)
    delay_out = joinpath(OUT_DIR,
        "linelist_$(FIT)$(TAG)$(SUFFIX)_delay.csv")
    CSV.write(delay_out, delays)
    @info "delay written" delay_out
    println(delays)
end

@info "written" out samples = SAMPLES
println(estimates)
