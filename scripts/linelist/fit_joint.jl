# Refit BVDOutbreakSize on line-list case counts.
#
# This model reads the INSP situation reports. A model fitted to the same
# outbreak from a case line list reads different data, so its results differ for
# two reasons at once: the models differ and the data differ. This refit replaces
# the DRC case streams in the observation manifest with counts derived from the
# line list and leaves the model untouched, so the remaining difference is a
# difference between models.
#
# The refit is partial. Only three streams have a line-list source (cumulative
# confirmed cases, cumulative reported cases, daily new suspected cases). Deaths,
# laboratory volumes, beds, recoveries and the Uganda exports have none, since a
# case line list records no outcome, and are carried over from the released
# manifest unchanged. So this fit is not a line-list-only fit, and
# scripts/linelist/fit_single.jl is the comparison to trust.
#
# The streams are counted by notification date, an event date, so their last
# fortnight is short by however long cases take to reach the export (median 4
# days, 90th percentile 12). The model reads that ragged edge as a fall in
# transmission. A construction that counts each case at the snapshot which first
# held it has no such edge, and fit_single.jl takes one via
# `--data=linelist_known`.
#
# Usage, from the repository root. Five to seven hours, so run it detached with
# the wrapper rather than calling Julia directly:
#
#   LINELIST_INPUT_DIR=<dir> scripts/linelist/run_joint_fit.sh
#
# See scripts/linelist/README.md for the inputs, the thread count and the
# sampler settings. Progress goes to logs/joint.log.

using Pkg: Pkg
Pkg.instantiate()

using BVDOutbreakSize
using CSV
using DataFrames
using Dates
using Statistics: quantile
using TOML

## The manifest builder and the input plumbing are shared with
## scripts/linelist/fit_single.jl, so a fit of one stream and a fit of all of
## them read the same substitution.
const BOS = dirname(dirname(@__DIR__))
include(joinpath(@__DIR__, "manifest.jl"))

const INPUT_DIR = linelist_input_dir()
const OUT_DIR = linelist_output_dir(BOS)

manifest = write_manifest(
    released = released_manifest(OUT_DIR),
    streams = joinpath(INPUT_DIR, "linelist_streams.csv"),
    out = joinpath(OUT_DIR, "linelist_observations.toml"),
    source = "DHIS2 case line list, counted by notification date")

## The triangle has to sit beside the manifest to be read at all, and a fit that
## silently drops it costs seven hours to discover.
place_onset_curve(INPUT_DIR, OUT_DIR)

obs = load_observations(manifest)
och = obs.onset_curve_history
isempty(och.increments) &&
    error("the onset triangle loaded empty, so this fit would silently drop " *
          "the onset stream")
@info "onset triangle" cells=length(och.increments) vintages=length(och.total_days)

## The headline joint fit, exactly as the released report runs it. Taken from
## this repository's own registry rather than reconstructed here, so it cannot
## silently drift from the fit it is being compared against.
include(joinpath(BOS, "docs", "fits", "registry.jl"))
specs = build_fit_specs(obs; run_sensitivity = false)
joint = specs[findfirst(s -> s.id == "joint", specs)]

@info "fitting" cutoff=obs.cutoff n=obs.n
started = now()
chn = joint.thunk()
@info "fitted" elapsed = canonicalize(round(now() - started, Minute))

summary = summary_table(chn,
    [:r, :r0, :doubling_time, :T, :R_T, :CFR, :C_T,
        :p_drc, :p_uganda, :k, :tau_test, :lambda_bg]; digits = 2)
CSV.write(joinpath(OUT_DIR, "linelist_posterior_summary.csv"), summary)

## `summary_table` publishes bounds without a median, which cannot fill a
## median-and-interval table. The release also ships `stream_estimates.csv`,
## which carries one, so the same shape is written here from the draws and the
## comparison reads both fits the same way.
function stream_estimates(chn, quantities)
    rows = NamedTuple[]
    for q in quantities
        draws = vec(Array(chn[q]))
        push!(rows,
            (; fit = "joint", quantity = string(q),
                median = quantile(draws, 0.5),
                lo30 = quantile(draws, 0.35), hi30 = quantile(draws, 0.65),
                lo60 = quantile(draws, 0.2), hi60 = quantile(draws, 0.8),
                lo90 = quantile(draws, 0.05), hi90 = quantile(draws, 0.95)))
    end
    return DataFrame(rows)
end

CSV.write(joinpath(OUT_DIR, "linelist_stream_estimates.csv"),
    stream_estimates(chn, (:R_T, :C_T, :r, :doubling_time)))

runs = [(h,
            forecast_reported(chn; horizon = h,
                obs_cases = obs.reported_cases,
                obs_deaths = obs.total_deaths,
                obs_confirmed = obs.confirmed_cases,
                obs_confirmed_deaths = obs.confirmed_deaths,
                obs_recovered = obs.recovered_cases,
                grid_n = obs.n))
        for h in (7, 14, 21, 28)]

CSV.write(joinpath(OUT_DIR, "linelist_forecast.csv"),
    forecast_archive(runs; made_date = obs.cutoff, thin = 5))

@info "written" dir = OUT_DIR
