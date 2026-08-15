# Refit BVDOutbreakSize on DHIS2 line-list case counts.
#
# The two models read different data, which is the deepest confound in
# docs/model-comparison.qmd: their outputs differ partly because the INSP
# situation reports and the DHIS2 line list differ. This refit replaces the DRC
# case streams in the BVDOutbreakSize observation manifest with counts derived
# from the line list, leaving the model itself untouched, so the remaining
# difference is a difference between models rather than between data sources.
#
# Only three streams have a line-list source (see
# bvd-analysis's stream builder). Deaths, laboratory volumes, beds,
# recoveries and the Uganda exports have none — the line list records no
# outcome — and are carried over from the released manifest unchanged. The
# refit is therefore partial, and its result is not a line-list-only fit.
#
# Run twice. On 2026-08-10, 6 hours 50 minutes and 135 divergent transitions,
# without the onset-reporting stream: `load_observations` reads the onset triangle
# from `onset_curve_scanned.csv` beside the manifest it is handed, and
# `load_onset_curve` returns an empty no-op rather than an error when that file is
# absent, so the stream vanished silently. On 2026-08-12, with the triangle in
# place, 5 hours 27 minutes and 94 divergences. The two agreed to within a
# rounding difference, so that stream is not what drives this fit's result.
#
# What does drive it, and why this script's output should not be read as the
# line-list answer: the streams it reads are counted by notification date, an
# event date, so their last fortnight is short by however long cases take to
# reach the export (median 4 days, 90th percentile 12). The model reads that
# ragged edge as a fall in transmission. bos_linelist_streams_known.csv counts
# each case at the snapshot that first held it instead, which has no such edge;
# scripts/linelist/fit_single.jl takes it via `--streams`. This script has
# not been re-run on that construction. Results and limitations are in
# docs/model-comparison.qmd.
#
# Usage, from the bvd-analysis repo root:
#
#   Rscript bvd-analysis's stream builder 20260803   # case streams
#   bvd-analysis's fetch_bos.sh                              # released assets
#   scripts/linelist/run_joint_fit.sh                   # the fit, detached
#
# Use the wrapper rather than calling Julia directly. Five to seven hours is long
# enough that a closed terminal or a sleeping laptop is the likeliest way the run
# ends, and the wrapper handles both. To run it in the foreground anyway:
#
#   julia -t 2 --project=docs \
#       scripts/linelist/fit_joint.jl
#
# `-t 2` matters: Julia defaults to one thread, and `nuts_sample` runs its
# chains with MCMCThreads, so without it the two chains run one after the other
# and the fit takes twice as long. Two is also the ceiling worth giving it,
# since the fit registry sets two chains and further threads sit idle. The first
# run instantiates and precompiles the BVDOutbreakSize docs environment, about
# 42 minutes, once.
#
# Sampler settings are the released report's own (500 samples, 2 chains), taken
# from that repo's fit registry rather than set here, so the refit stays
# comparable with the release it is being compared against. Progress goes to
# logs/joint.log.
#
# Reads bos_observations.toml and bos_linelist_streams.csv from INPUT_DIR;
# writes bos_linelist_forecast.csv and bos_linelist_posterior_summary.csv into
# OUT_DIR. Both default beside this repo and are overridable with
# LINELIST_INPUT_DIR and LINELIST_OUT_DIR.

using Pkg: Pkg
Pkg.instantiate()

using BVDOutbreakSize
using CSV
using DataFrames
using Dates
using Statistics: quantile
using TOML

# This repo is the model; the line-list inputs are produced elsewhere. INPUT_DIR
# points at whatever wrote them, defaulting to a sibling bvd-analysis clone.
const BOS = dirname(dirname(@__DIR__))
const INPUT_DIR = get(ENV, "LINELIST_INPUT_DIR",
    joinpath(dirname(BOS), "bvd-analysis", "ignore", "report"))
const OUT_DIR = get(ENV, "LINELIST_OUT_DIR", joinpath(BOS, "ignore", "linelist"))
mkpath(OUT_DIR)

## The manifest builder is shared with scripts/linelist/fit_single.jl,
## so a fit of one stream and a fit of all of them read the same substitution.
include(joinpath(@__DIR__, "manifest.jl"))

manifest = write_manifest(
    released = joinpath(INPUT_DIR, "bos_observations.toml"),
    streams = joinpath(INPUT_DIR, "bos_linelist_streams.csv"),
    out = joinpath(OUT_DIR, "bos_linelist_observations.toml"),
    source = "DHIS2 case line list (bvd-data-mirror), counted by " *
             "notification date; built by " *
             "bvd-analysis's stream builder")

obs = load_observations(manifest)

## The headline joint fit, exactly as the released report runs it. Taken from
## that repo's own registry rather than reconstructed here, so this cannot
## silently drift from the fit it is being compared against.
include(joinpath(BOS, "docs", "fits", "registry.jl"))
specs = build_fit_specs(obs; run_sensitivity = false)
joint = specs[findfirst(s -> s.id == "joint", specs)]

@info "fitting" cutoff = obs.cutoff n = obs.n
started = now()
chn = joint.thunk()
@info "fitted" elapsed = canonicalize(round(now() - started, Minute))

summary = summary_table(chn,
    [:r, :r0, :doubling_time, :T, :R_T, :CFR, :C_T,
        :p_drc, :p_uganda, :k, :tau_test, :lambda_bg]; digits = 2)
CSV.write(joinpath(OUT_DIR, "bos_linelist_posterior_summary.csv"), summary)

## `summary_table` publishes bounds without a median, which cannot fill a
## median-and-interval table. The release also ships `stream_estimates.csv`,
## which carries one, so the same shape is written here from the draws and the
## comparison reads both fits the same way.
function stream_estimates(chn, quantities)
    rows = NamedTuple[]
    for q in quantities
        draws = vec(Array(chn[q]))
        push!(rows, (; fit = "joint", quantity = string(q),
            median = quantile(draws, 0.5),
            lo30 = quantile(draws, 0.35), hi30 = quantile(draws, 0.65),
            lo60 = quantile(draws, 0.2), hi60 = quantile(draws, 0.8),
            lo90 = quantile(draws, 0.05), hi90 = quantile(draws, 0.95)))
    end
    return DataFrame(rows)
end

CSV.write(joinpath(OUT_DIR, "bos_linelist_stream_estimates.csv"),
    stream_estimates(chn, (:R_T, :C_T, :r, :doubling_time)))

runs = [(h, forecast_reported(chn; horizon = h,
             obs_cases = obs.reported_cases,
             obs_deaths = obs.total_deaths,
             obs_confirmed = obs.confirmed_cases,
             obs_confirmed_deaths = obs.confirmed_deaths,
             obs_recovered = obs.recovered_cases,
             grid_n = obs.n))
        for h in (7, 14, 21, 28)]

CSV.write(joinpath(OUT_DIR, "bos_linelist_forecast.csv"),
    forecast_archive(runs; made_date = obs.cutoff, thin = 5))

@info "written" dir = OUT_DIR
