# Collect what the single-stream line-list refits (`fit_single.jl`) have
# written and reduce it to the comparison the grid exists to make: the same
# model at the same cut-off on two data sources, situation reports against
# the line list. Nothing here refits or reconstructs anything — it only reads
# what `fit_single.jl` already wrote, so a change to these numbers has to
# come from rerunning the grid, not from editing this script.
#
# `comparison.csv` is the file to read first. Everything else is either its
# inputs or the detail behind it.
#
# The onset-to-report delay comparison lives alongside the Rt figures for
# the same reason `delays.jl` fits it in the first place: `onsets` estimates
# that delay from the reporting triangle rather than taking it as a prior,
# so it is the one quantity in this grid with an independent, out-of-band
# check available (the `bvd-internal-cmmid` fit on the same interval).
# `delay_comparison.csv` puts the two side by side rather than leaving the
# comparison to be redone by hand from two separate files.
#
# Usage, from the repository root:
#
#   LINELIST_OUT_DIR=<dir> julia --project=docs scripts/linelist/plot_rt.jl
#
# `LINELIST_OUT_DIR` is resolved exactly as the fits resolve it (see
# `linelist_output_dir` in `manifest.jl`): it defaults to `ignore/linelist`
# and refuses any other path inside the repository, since this repository is
# public and these outputs derive from the line list.
#
# Writes, into the output directory:
#
#   comparison.csv          the headline: one row per (fit, delays,
#       quantity), the two data sources side by side with their difference
#       and ratio. Rt at the cut-off and at 7, 14 and 28 days before it,
#       then R_T, R0 and C_T
#   rt_all.csv              every `linelist_*_rt.csv` concatenated, columns
#       unchanged
#   stream_estimates_all.csv every `linelist_*_stream_estimates.csv`
#       concatenated, columns unchanged
#   rt_confirmed.png         band-and-line comparison, one series per
#   rt_onsets.png             (data, delays) combination — written only for
#       a fit that has at least one input file; a fit with none is skipped
#       rather than plotted empty
#   rt_summary.csv           median and 30/60/90% bounds at each series' own
#       last fitted date and at 7, 14 and 28 days before it
#   diagnostics_all.csv     worst R-hat, smallest bulk ESS and divergence
#       count per run, with a warning naming any run past the thresholds
#   delay_comparison.csv     every `linelist_*_delay.csv` row, plus the two
#       independently fitted `bvd-internal-cmmid` reference values
#
# Every input is validated against the schema `fit_single.jl` writes before
# it is used: a file with the wrong columns, or a value outside the fit,
# data or delays vocabulary those scripts use, is a loud error naming the
# file rather than a row silently dropped from a figure.

using Pkg: Pkg
Pkg.instantiate()

using BVDOutbreakSize
using CSV
using DataFrames
using Dates: Date, Day, date2epochdays, epochdays2date
using Statistics: quantile
using TOML
using CairoMakie

const BOS = dirname(dirname(@__DIR__))
include(joinpath(@__DIR__, "manifest.jl"))

const OUT_DIR = linelist_output_dir(BOS)

## The vocabulary `fit_single.jl` writes. Checked here rather than assumed,
## since a value outside it (a typo'd `--data=`, a delay config renamed
## upstream without this script following) would otherwise plot silently
## under whatever label it happened to carry.
## The comparison is `confirmed` and `onsets`, both laboratory-confirmed
## cases, against two observation sets. `cases` is still a valid fit and is
## listed so a stray file is not rejected, but it is not part of the
## comparison: the two reported-case streams are different case definitions.
const VALID_FITS = ("confirmed", "onsets", "cases")

## The fits the comparison is drawn and tabulated for. Files for the others
## are read if present, so a stray file is not rejected, but they are not
## compared and the script says so on stderr when it finds them.
##
## `cases`: the two reported-case streams are different case definitions, and
## its `cmmid_*` chains predate the delay overrides reaching that fit, so they
## carry package-default delays under a cmmid label.
##
## `confirmed`: on linelist_known it does not sample. Both chains freeze, one
## at the `inv_sqrt_k = 0` truncation boundary, giving a max R-hat of 2.63 and
## a bulk ESS of 2.4 -- the line-list confirmed history is too smooth for the
## negative-binomial observation model. Its sitrep arm samples, but a
## one-sided arm is not a comparison.
const COMPARISON_FITS = ("onsets",)

const VALID_DATA = ("sitrep", "linelist_known")
const VALID_QUANTITIES = ("delay_p50", "delay_p90")

## What `fit_single.jl` writes into stream_estimates.csv: the cut-off and
## basic reproduction numbers and the final size.
const VALID_STREAM_QUANTITIES = ("R_T", "C_T", "R0")

## The two observation sets, in the direction the contrast is written:
## `linelist_known` against `sitrep`, so a positive difference means the line
## list gives the larger value.
const CONTRAST_BASE = "sitrep"
const CONTRAST_ALT = "linelist_known"
is_valid_delays(x) = x == "repo" || startswith(x, "cmmid_")

const RT_COLS = ["fit", "data", "delays", "day", "date", "median",
    "lo30", "hi30", "lo60", "hi60", "lo90", "hi90"]
const DELAY_COLS = ["fit", "data", "delays", "quantity", "median",
    "lo30", "hi30", "lo60", "hi60", "lo90", "hi90"]
const STREAM_COLS = DELAY_COLS
const DIAG_COLS = ["fit", "data", "delays", "samples", "chains", "max_rhat",
    "min_ess_bulk", "divergences"]

## When a run is reported as not to be read rather than merely noted.
## R-hat and bulk ESS are the usual thresholds; the divergence one is a
## fraction of the run's own draw count rather than a fixed number, so it
## does not change meaning if the sample size does.
const RHAT_MAX = 1.01
const ESS_MIN = 400
const DIVERGENCE_FRACTION = 0.01

function check_columns(df, path, required)
    have = names(df)
    missing_cols = [c for c in required if !(c in have)]
    isempty(missing_cols) ||
        error("$path is missing required column(s): " *
              join(missing_cols, ", ") * "; found: " * join(have, ", "))
    return df
end

## Shared by both file kinds: both carry `fit`, `data` and `delays`, and a
## row outside that vocabulary is exactly the kind of silent-mismatch this
## script is meant to refuse rather than plot.
function check_domains(df, path)
    bad = unique(filter(!in(VALID_FITS), df.fit))
    isempty(bad) || error("$path has unexpected `fit` value(s): " *
          join(bad, ", ") * "; expected one of " * join(VALID_FITS, ", "))
    bad = unique(filter(!in(VALID_DATA), df.data))
    isempty(bad) || error("$path has unexpected `data` value(s): " *
          join(bad, ", ") * "; expected one of " * join(VALID_DATA, ", "))
    bad = unique(filter(!is_valid_delays, df.delays))
    isempty(bad) || error("$path has unexpected `delays` value(s): " *
          join(bad, ", ") *
          "; expected `repo` or a name starting `cmmid_`")
    return df
end

## `date` is read as an ISO string by CSV.jl's own type inference in every
## case seen so far, but that inference is not a guarantee: coerce
## explicitly and fail on the file rather than on the first date arithmetic
## three functions further down.
function ensure_dates!(df, path)
    eltype(df.date) <: Date && return df
    try
        df.date = Date.(string.(df.date))
    catch e
        error("$path: could not parse the `date` column as ISO dates ($e)")
    end
    return df
end

function read_rt(path)
    df = CSV.read(path, DataFrame)
    check_columns(df, path, RT_COLS)
    check_domains(df, path)
    ensure_dates!(df, path)
    return select(df, RT_COLS)
end

function read_delay(path)
    df = CSV.read(path, DataFrame)
    check_columns(df, path, DELAY_COLS)
    check_domains(df, path)
    bad = unique(filter(!in(VALID_QUANTITIES), df.quantity))
    isempty(bad) || error("$path has unexpected `quantity` value(s): " *
          join(bad, ", ") * "; expected one of " *
          join(VALID_QUANTITIES, ", "))
    return select(df, DELAY_COLS)
end

function read_diag(path)
    df = CSV.read(path, DataFrame)
    check_columns(df, path, DIAG_COLS)
    check_domains(df, path)
    return select(df, DIAG_COLS)
end

function read_stream(path)
    df = CSV.read(path, DataFrame)
    check_columns(df, path, STREAM_COLS)
    check_domains(df, path)
    bad = unique(filter(!in(VALID_STREAM_QUANTITIES), df.quantity))
    isempty(bad) || error("$path has unexpected `quantity` value(s): " *
          join(bad, ", ") * "; expected one of " *
          join(VALID_STREAM_QUANTITIES, ", "))
    return select(df, STREAM_COLS)
end

## `_pilot` outputs come from ten-sample timing runs. They carry the same
## schema and the same names as a real fit, so left in they would be drawn
## as a series on equal footing with the 500-sample fits and nothing on the
## figure would say which was which. `--pilot` exists so a timing run cannot
## be mistaken for the real one; honouring that here is the other half of it.
const PILOT_RE = r"_pilot_(rt|delay|stream_estimates)\.csv$"

rt_files = sort(filter(readdir(OUT_DIR)) do f
    occursin(r"^linelist_.*_rt\.csv$", f) && !occursin(PILOT_RE, f)
end)
isempty(rt_files) &&
    error("no linelist_*_rt.csv files found in $OUT_DIR (pilot runs are " *
          "ignored). Run `scripts/linelist/run_grid.sh --stage=1` first to " *
          "produce them.")

delay_files = sort(filter(readdir(OUT_DIR)) do f
    occursin(r"^linelist_.*_delay\.csv$", f) && !occursin(PILOT_RE, f)
end)

stream_files = sort(filter(readdir(OUT_DIR)) do f
    occursin(r"^linelist_.*_stream_estimates\.csv$", f) &&
        !occursin(PILOT_RE, f)
end)

diag_files = sort(filter(readdir(OUT_DIR)) do f
    occursin(r"^linelist_.*_diagnostics\.csv$", f) &&
        !occursin(r"_pilot_diagnostics\.csv$", f)
end)

rt_all = reduce(vcat,
    [read_rt(joinpath(OUT_DIR, f)) for f in rt_files])
delay_dfs = [read_delay(joinpath(OUT_DIR, f)) for f in delay_files]
stream_dfs = [read_stream(joinpath(OUT_DIR, f)) for f in stream_files]
stream_all = isempty(stream_dfs) ? DataFrame() : reduce(vcat, stream_dfs)
diag_dfs = [read_diag(joinpath(OUT_DIR, f)) for f in diag_files]
diag_all = isempty(diag_dfs) ? DataFrame() : reduce(vcat, diag_dfs)

## Every series within a fit has to sit on one grid. The comparison is the same
## model at the same cut-off, so a series ending on a different day is either a
## run given a different `LINELIST_AS_OF` or a file left over from an earlier
## cut-off; pairing it against the others would contrast two questions and the
## contrast table would not show it. This is the check that the one cut-off the
## grid is supposed to share was actually passed to every run.
function check_grid_alignment(rt)
    for g in groupby(rt, :fit; sort = true)
        ends = combine(groupby(DataFrame(g), [:data, :delays]),
            :date => maximum => :last_date)
        length(unique(ends.last_date)) <= 1 && continue
        listing = join(["  $(r.data)/$(r.delays) ends $(r.last_date)"
                        for r in eachrow(sort(ends, :last_date))], "\n")
        error("the `$(first(g.fit))` series do not share a last fitted " *
              "date, so they are not on one grid and must not be " *
              "compared:\n" * listing *
              "\nRerun with a single LINELIST_AS_OF, or clear the outputs " *
              "from the earlier cut-off out of $OUT_DIR.")
    end
    return nothing
end
check_grid_alignment(rt_all)

CSV.write(joinpath(OUT_DIR, "rt_all.csv"), rt_all)
if !isempty(stream_all)
    CSV.write(joinpath(OUT_DIR, "stream_estimates_all.csv"), stream_all)
end

## Sampler quality, and a loud line for any run that should not be read as a
## result. A difference between two data sources is only attributable to the
## data if both sides sampled; the line-list arms have run more divergent than
## their situation-report pairs, which is exactly the asymmetry that would
## otherwise be invisible in a table of medians.
if !isempty(diag_all)
    CSV.write(joinpath(OUT_DIR, "diagnostics_all.csv"), diag_all)
    for r in eachrow(diag_all)
        draws = r.samples * r.chains
        reasons = String[]
        ismissing(r.max_rhat) || isnan(r.max_rhat) || r.max_rhat <= RHAT_MAX ||
            push!(reasons, "max_rhat $(r.max_rhat) > $RHAT_MAX")
        ismissing(r.min_ess_bulk) || isnan(r.min_ess_bulk) ||
            r.min_ess_bulk >= ESS_MIN ||
            push!(reasons, "min_ess_bulk $(r.min_ess_bulk) < $ESS_MIN")
        r.divergences <= DIVERGENCE_FRACTION * draws ||
            push!(reasons,
            "$(r.divergences) divergent transitions in $draws draws")
        isempty(reasons) ||
            @warn "fit quality: $(r.fit)/$(r.data)/$(r.delays)" *
                  " -- " * join(reasons, "; ")
    end
end

## Nested 30/60/90% ribbons under a solid median line, one series per
## (data, delays) combination, in the alpha scheme `plot_rt` in
## src/plots.jl uses for the same shape of band. `sort = true` on the
## grouping is what gives a rerun the same colours: `wong_colors()` is
## assigned by position, so an unsorted group order would recolour the
## figure on every run for no reason but file-listing order.
function plot_fit_rt(df_fit, fit_name)
    ## The cap is taken from each series' typical 90% upper band (the median
    ## over its days) rather than the largest one anywhere, which is what
    ## `plot_rt_streams` does and for the reason it gives. The early weeks of a
    ## weakly-informed fit can carry a 90% band an order of magnitude above the
    ## rest: capping on the maximum lets that one spike squash every series flat
    ## against the axis and the comparison disappears. The median over days is
    ## robust to it while still scaling to whatever is actually being drawn.
    function series_top(g)
        v = collect(skipmissing(g.hi90))
        return isempty(v) ? 0.0 : quantile(v, 0.5)
    end
    groups = groupby(df_fit, [:data, :delays]; sort = true)
    tops = [series_top(g) for g in groups]
    ytop = max(1.2, 1.2 * (isempty(tops) ? 0.0 : maximum(tops)))

    fig = Figure(; size = (1100, 460))
    ax = Axis(fig[1, 1]; xlabel = "Date",
        ylabel = "Reproduction number Rt",
        title = "Estimated Rt: $(fit_name) fit, line-list comparison",
        limits = (nothing, (0.0, ytop)),
        xticklabelrotation = pi / 6)

    colours = CairoMakie.Makie.wong_colors()
    for (i, g) in enumerate(groups)
        s = sort(DataFrame(g), :date)
        x = Float64.(date2epochdays.(s.date))
        colour = colours[mod1(i, length(colours))]
        band!(ax, x, s.lo90, s.hi90; color = (colour, 0.15))
        band!(ax, x, s.lo60, s.hi60; color = (colour, 0.28))
        band!(ax, x, s.lo30, s.hi30; color = (colour, 0.42))
        lines!(ax, x, s.median; color = colour, linewidth = 2,
            label = "$(s.data[1]) / $(s.delays[1])")
    end
    hlines!(ax, [1.0]; color = (:grey, 0.8), linestyle = :dash,
        linewidth = 2)

    all_x = Float64.(date2epochdays.(df_fit.date))
    lo = floor(Int, minimum(all_x))
    hi = ceil(Int, maximum(all_x))
    CairoMakie.xlims!(ax, lo, hi)
    ## About a dozen dated ticks whatever the span. A fixed weekly step gives
    ## twenty-six labels over a six-month grid, which overlap into an unreadable
    ## band along the axis.
    step = max(7, 7 * cld(cld(hi - lo, 12), 7))
    ax.xticks = collect(lo:step:hi)
    ax.xtickformat = vals -> [string(epochdays2date(round(Int, v)))
                              for v in vals]

    Legend(fig[1, 2], ax, "data / delays"; framevisible = true)
    return fig
end

written = String[joinpath(OUT_DIR, "rt_all.csv")]
isempty(stream_all) || push!(written,
    joinpath(OUT_DIR, "stream_estimates_all.csv"))
for fit_name in COMPARISON_FITS
    df_fit = filter(:fit => ==(fit_name), rt_all)
    isempty(df_fit) && continue
    fig = plot_fit_rt(df_fit, fit_name)
    path = joinpath(OUT_DIR, "rt_$(fit_name).png")
    CairoMakie.save(path, fig)
    push!(written, path)
end

## Four rows per series where the series reaches that far back: the last
## fitted date, then 7, 14 and 28 days earlier. A day one of those offsets
## lands on but the series does not reach (a series shorter than 28 days,
## or a masked pre-established day within it) is omitted rather than
## interpolated or extrapolated, so every row in this file is a value the
## fit actually produced.
function summary_rows(rt_all)
    rows = NamedTuple[]
    for g in groupby(rt_all, [:fit, :data, :delays]; sort = true)
        s = sort(DataFrame(g), :date)
        last_date = maximum(s.date)
        for offset_days in (0, 7, 14, 28)
            idx = findfirst(==(last_date - Day(offset_days)), s.date)
            isnothing(idx) && continue
            r = s[idx, :]
            push!(rows,
                (; fit = r.fit, data = r.data, delays = r.delays,
                    offset_days, date = r.date, median = r.median,
                    lo30 = r.lo30, hi30 = r.hi30, lo60 = r.lo60,
                    hi60 = r.hi60, lo90 = r.lo90, hi90 = r.hi90))
        end
    end
    return DataFrame(rows)
end

rt_summary = summary_rows(rt_all)
CSV.write(joinpath(OUT_DIR, "rt_summary.csv"), rt_summary)
push!(written, joinpath(OUT_DIR, "rt_summary.csv"))

## The contrast the grid exists to make, in one file: for each fit and delay
## configuration, the two data sources side by side on one quantity per row.
##
## Only the medians are differenced. The two fits are not independent samples
## of one population — they are the same model on two constructions of the
## same outbreak — so a difference of medians is a description of how far
## apart the two answers sit, not a test, and no interval is put on it. The
## per-source 90% bounds are carried so a reader can see whether the gap is
## small or large next to the width of either fit.
##
## `C_T` is included because the fits write it, but neither single-stream fit
## identifies ascertainment, so it is not a size estimate and a ratio of two
## `C_T` values is a ratio of two quantities that absorb the level difference
## between the sources. Read the Rt rows; see scripts/linelist/README.md.
const QUANTITY_ORDER = ["Rt_cutoff", "Rt_cutoff_minus7", "Rt_cutoff_minus14",
    "Rt_cutoff_minus28", "R_T", "R0", "C_T"]
qorder(q) = something(findfirst(==(q), QUANTITY_ORDER),
    length(QUANTITY_ORDER) + 1)

function contrast_table(rt_summary, stream_all)
    long = NamedTuple[]
    for r in eachrow(rt_summary)
        q = r.offset_days == 0 ? "Rt_cutoff" :
            "Rt_cutoff_minus$(r.offset_days)"
        push!(long, (; fit = r.fit, delays = r.delays, quantity = q,
            data = r.data, date = r.date, median = r.median,
            lo90 = r.lo90, hi90 = r.hi90))
    end
    for r in eachrow(stream_all)
        push!(long, (; fit = r.fit, delays = r.delays, quantity = r.quantity,
            data = r.data, date = missing, median = r.median,
            lo90 = r.lo90, hi90 = r.hi90))
    end
    isempty(long) && return DataFrame()

    L = filter(:fit => in(COMPARISON_FITS), DataFrame(long))
    rows = NamedTuple[]
    for g in groupby(L, [:fit, :delays, :quantity]; sort = true)
        b = g[g.data .== CONTRAST_BASE, :]
        a = g[g.data .== CONTRAST_ALT, :]
        ## A configuration run on one source but not the other is not a
        ## contrast, so it is left out rather than half-reported.
        (nrow(b) == 1 && nrow(a) == 1) || continue
        bm, am = b.median[1], a.median[1]
        push!(rows,
            (; fit = first(g.fit), delays = first(g.delays),
                quantity = first(g.quantity), date = b.date[1],
                sitrep = bm, sitrep_lo90 = b.lo90[1],
                sitrep_hi90 = b.hi90[1],
                linelist = am, linelist_lo90 = a.lo90[1],
                linelist_hi90 = a.hi90[1],
                difference = am - bm,
                ratio = bm == 0 ? missing : am / bm))
    end
    isempty(rows) && return DataFrame()
    out = DataFrame(rows)
    sort!(out, [:fit, :delays, order(:quantity, by = qorder)])
    return out
end

comparison = contrast_table(rt_summary, stream_all)
if isempty(comparison)
    @info "no (fit, delays, quantity) has been run on both data sources " *
          "yet, so comparison.csv is not written"
else
    CSV.write(joinpath(OUT_DIR, "comparison.csv"), comparison)
    push!(written, joinpath(OUT_DIR, "comparison.csv"))
end

## The reference rows record `bvd-internal-cmmid`'s own fit of the onset to
## line-list-classification delay, the same interval `onsets` estimates from
## the reporting triangle: two independent fits of one quantity, worth
## reading side by side rather than as two files a reviewer has to open
## separately.
##
## Read from `BVD_DELAY_DIR` at run time, not written here. Those are fitted
## values from a private repository whose disclosure rules keep them there,
## and this file is tracked in a public one. Without the variable set the
## comparison is simply omitted, since a stage-1 grid has no reason to have
## it set and the fitted rows are still worth writing on their own.
##
## `Cas confirmé` is the stratum to compare against: the reporting triangle
## is laboratory-confirmed cases, so the pooled fit across all
## classifications is a different population.
## `delays = "NA"` rather than left blank, since blank already means "not
## queried" for the interval columns, and these are not a delay-config run.
const REFERENCE_CLASS = "Cas confirmé"
const REFERENCE_SOURCE = "bvd-internal-cmmid, $REFERENCE_CLASS, onset to " *
                         "line-list classification"

function reference_delays()
    dir = get(ENV, "BVD_DELAY_DIR", "")
    if isempty(dir)
        @info "BVD_DELAY_DIR unset, so delay_comparison.csv carries the " *
              "fitted rows only"
        return nothing
    end
    path = joinpath(dir, "bayes_delay_by_classification.csv")
    isfile(path) || error("missing $path, which BVD_DELAY_DIR should hold. " *
          "See scripts/linelist/delays.jl for what that variable points at.")
    df = CSV.read(path, DataFrame)
    check_columns(df, path, ["classification", "p", "median"])
    rows = NamedTuple[]
    for (p, q) in ((0.5, "delay_p50"), (0.9, "delay_p90"))
        sel = df[(df.classification .== REFERENCE_CLASS) .& (df.p .== p), :]
        nrow(sel) == 1 ||
            error("expected exactly one $REFERENCE_CLASS row at p = $p in " *
                  "$path, found $(nrow(sel))")
        push!(rows,
            (; fit = "reference", data = "cmmid", delays = "NA",
                quantity = q, median = Float64(sel.median[1]),
                lo30 = missing, hi30 = missing, lo60 = missing,
                hi60 = missing, lo90 = missing, hi90 = missing,
                source = REFERENCE_SOURCE))
    end
    return DataFrame(rows)
end

reference_rows = reference_delays()

fitted_delays = isempty(delay_dfs) ? DataFrame() : reduce(vcat, delay_dfs)
if !isempty(fitted_delays)
    fitted_delays.source = fill(missing, nrow(fitted_delays))
end
delay_comparison = if isempty(fitted_delays)
    reference_rows
elseif isnothing(reference_rows)
    fitted_delays
else
    vcat(fitted_delays, reference_rows)
end
if isnothing(delay_comparison) || isempty(delay_comparison)
    @info "no fitted delays and no reference, so delay_comparison.csv is " *
          "not written"
else
    CSV.write(joinpath(OUT_DIR, "delay_comparison.csv"), delay_comparison)
    push!(written, joinpath(OUT_DIR, "delay_comparison.csv"))
end

## What each observation set actually held at the cut-off, read back out of the
## manifests the fits were run against rather than recomputed from the inputs.
## The manifest is what the model saw, so it is the honest source for a table
## saying what the two data sources disagree about; recomputing from the stream
## files would re-derive it and could disagree with the fit.
function inputs_at_cutoff()
    rows = NamedTuple[]
    for f in sort(readdir(OUT_DIR))
        m = match(r"^linelist_(.+)_observations\.toml$", f)
        isnothing(m) && continue
        data = m.captures[1]
        occursin("_pilot", data) && continue
        raw = TOML.parsefile(joinpath(OUT_DIR, f))
        as_of = raw["as_of_date"]
        ## The last vintage at or before the cut-off, which is what
        ## `load_observations` reads: a block can run past it (the known-by
        ## streams do) and those later vintages are dropped.
        function at_cutoff(block)
            haskey(raw, block) || return (missing, missing)
            d, v = raw[block]["dates"], raw[block]["values"]
            i = findlast(x -> x <= as_of, d)
            isnothing(i) ? (missing, missing) : (d[i], v[i])
        end
        cd_, cv = at_cutoff("confirmed_case_history")
        rd_, rv = at_cutoff("reported_case_history")
        push!(rows,
            (; data, as_of, confirmed = cv, confirmed_date = cd_,
                reported = rv, reported_date = rd_,
                reported_scalar = get(get(raw, "reported_cases",
                    Dict{String, Any}()), "value", missing)))
    end
    return DataFrame(rows)
end

inputs = inputs_at_cutoff()
if !isempty(inputs)
    CSV.write(joinpath(OUT_DIR, "inputs_at_cutoff.csv"), inputs)
    push!(written, joinpath(OUT_DIR, "inputs_at_cutoff.csv"))
end

## A fit that was read but is not compared is said so out loud. Its rows are
## still in rt_all.csv and rt_summary.csv, and a reader who does not know it
## was excluded would take its absence from comparison.csv for "no difference".
noncomparison = unique(filter(!in(COMPARISON_FITS), rt_all.fit))
isempty(noncomparison) ||
    @warn "read files for fit(s) outside the comparison; they are in " *
          "rt_all.csv and rt_summary.csv but are not plotted and not in " *
          "comparison.csv" fits=join(noncomparison, ", ") compared=join(COMPARISON_FITS, ", ")

combos = unique(select(rt_all, [:fit, :data, :delays]))
sort!(combos, [:fit, :data, :delays])
println("Read $(length(rt_files)) rt file(s) and $(length(delay_files)) " *
        "delay file(s) from $OUT_DIR")
println("(fit, data, delays) combinations found:")
for r in eachrow(combos)
    println("  $(r.fit), $(r.data), $(r.delays)")
end
println("Written:")
for p in written
    println("  $p")
end
