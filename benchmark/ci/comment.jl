# Build the PR benchmark comment from AirspeedVelocity result JSON.
#
# Usage:
#   julia --project=benchmark/ci benchmark/ci/comment.jl \
#       <results-dir> <package> <base-rev> <head-rev> <out.md>
#
# `run_pair.jl` runs the suite; this turns what it wrote into the comment.
# AirspeedVelocity's own table is a flat list of every leaf with a median
# ratio, no neutral band, no threshold and no way to set one. Two things
# this repository needs are therefore not in it.
#
# The first is a band measured from the run rather than chosen. A band
# fixed below the harness's own resolution reports noise as regression.
#
# The second is the spread. Unrelated components share no cause, so all of
# them moving by one factor points at an environment difference rather than
# at the diff. A flat table cannot show it.
#
# Both are recoverable here because AirspeedVelocity writes the raw
# per-sample times into its JSON, not just a summary.
using JSON3

const COMMENT_MARKER = "<!-- benchmark-comparison -->"

## The load-time benchmark AirspeedVelocity injects. It is dropped rather
## than reported. BenchmarkTools runs a warmup evaluation before it samples,
## and that warmup performs the `using`, so the in-process sample times a
## warm re-import rather than a load. Further samples relaunch Julia and do
## measure a load, but this comment reports a minimum, so the warm sample
## always wins. `run_pair.jl` asks for one sample so the relaunches are not
## paid for either.
const LOAD_KEY = "time_to_load"

## Band floor and ceiling. The floor keeps a quiet run from claiming a
## precision no repeated timing has; the ceiling keeps the bucket edges
## ordered and says plainly that a run this noisy is not reporting.
const BAND_FLOOR = 0.02
const BAND_CEILING = 0.2
const BAND_UNMEASURED = 0.05

## Quantile of the per-benchmark sample spread that sets the band.
const BAND_QUANTILE = 0.9

## Summary row order: the AD-free evaluation first, then the backends.
const GROUP_ORDER = ["Log density", "Mooncake", "Enzyme"]

# ---- loading ---------------------------------------------------------------

## One benchmark on one revision: the samples, reduced to what the comment
## needs. `best` is the estimator the ratio uses, `spread` the dispersion
## across samples.
struct Stat
    best::Float64
    spread::Float64
    memory::Float64
end

function percentile(xs, p)
    isempty(xs) && return NaN
    sorted = sort(xs)
    return sorted[clamp(ceil(Int, p * length(sorted)), 1, length(sorted))]
end

## Minimum, not median: a benchmark's fastest sample is the one least
## disturbed by whatever else the runner was doing.
function stat_from(node)
    memory = haskey(node, "memory") ? Float64(node["memory"]) : NaN
    if haskey(node, "times") && !isempty(node["times"])
        times = collect(Float64, node["times"])
        lo, hi = percentile(times, 0.25), percentile(times, 0.75)
        med = percentile(times, 0.5)
        spread = med > 0 ? (hi - lo) / med : NaN
        return Stat(minimum(times), spread, memory)
    elseif haskey(node, "median")
        ## Fallback for a results file carrying only the summary.
        med = Float64(node["median"])
        lo = haskey(node, "25") ? Float64(node["25"]) : NaN
        hi = haskey(node, "75") ? Float64(node["75"]) : NaN
        spread = (med > 0 && !isnan(lo) && !isnan(hi)) ? (hi - lo) / med : NaN
        return Stat(med, spread, memory)
    end
    return Stat(NaN, NaN, memory)
end

## Walk the BenchmarkTools group. Leaves carry `times`, inner nodes `data`.
function flatten!(out, node, prefix)
    (node isa AbstractDict || node isa JSON3.Object) || return out
    if haskey(node, "times") || haskey(node, "median")
        isempty(prefix) || (out[prefix] = stat_from(node))
    elseif haskey(node, "data")
        for (k, v) in pairs(node["data"])
            key = String(k)
            flatten!(out, v, isempty(prefix) ? key : string(prefix, " / ", key))
        end
    end
    return out
end

## `run_pair.jl` names each arm's file after the label its caller passes,
## so the path is exact.
function find_results(dir, pkg, rev)
    path = joinpath(dir, "results_$(pkg)@$(rev).json")
    isfile(path) && return path
    available = filter(f -> endswith(f, ".json"), readdir(dir))
    return error("no results for $pkg@$rev in $dir; found $available")
end

function load_arm(dir, pkg, rev)
    data = JSON3.read(read(find_results(dir, pkg, rev), String))
    return flatten!(Dict{String, Stat}(), data, "")
end

# ---- formatting ------------------------------------------------------------

is_ad(name) = startswith(name, "AD gradients")
is_load(name) = occursin(LOAD_KEY, name)

function fmt_time(ns)
    isnan(ns) && return "—"
    if ns < 1.0e3
        return string(round(ns; digits = 1), " ns")
    elseif ns < 1.0e6
        return string(round(ns / 1.0e3; digits = 2), " μs")
    elseif ns < 1.0e9
        return string(round(ns / 1.0e6; digits = 2), " ms")
    else
        return string(round(ns / 1.0e9; digits = 2), " s")
    end
end

fmt_pct(x) = isnan(x) ? "—" : string(round(100x; digits = 1), "%")

function fmt_ratio(r, band)
    isnan(r) && return "—"
    marker = r > 1 + band ? "🔴" : r < 1 - band ? "🟢" : "⚪"
    return string(marker, " ", round(r; digits = 2), "×")
end

struct Row
    name::String
    main_time::Float64
    pr_time::Float64
    main_mem::Float64
    pr_mem::Float64
    noise::Float64
    time_ratio::Float64
    status::Symbol
end

function build_rows(pr, main)
    rows = Row[]
    for name in sort(collect(union(keys(pr), keys(main))))
        inpr, inmain = haskey(pr, name), haskey(main, name)
        pt = inpr ? pr[name].best : NaN
        mt = inmain ? main[name].best : NaN
        pm = inpr ? pr[name].memory : NaN
        mm = inmain ? main[name].memory : NaN
        spreads = filter(
            !isnan,
            [inpr ? pr[name].spread : NaN, inmain ? main[name].spread : NaN]
        )
        noise = isempty(spreads) ? NaN : maximum(spreads)
        ratio = (inpr && inmain && mt > 0) ? pt / mt : NaN
        status = inpr && inmain ? :both : inpr ? :new : :removed
        push!(rows, Row(name, mt, pt, mm, pm, noise, ratio, status))
    end
    return rows
end

sort_key(r) = isnan(r.time_ratio) ? Inf : abs(r.time_ratio - 1)

function status_note(r)
    r.status === :new && return " *(new)*"
    r.status === :removed && return " *(removed)*"
    return ""
end

function render_table(rows, band)
    isempty(rows) && return "_none_\n"
    io = IOBuffer()
    println(io, "| Benchmark | main | PR | time | spread | memory |")
    println(io, "|---|---|---|---|---|---|")
    for r in rows
        memratio = if isnan(r.pr_mem) || isnan(r.main_mem) || r.main_mem == 0
            NaN
        else
            r.pr_mem / r.main_mem
        end
        println(
            io, "| ", r.name, status_note(r),
            " | ", fmt_time(r.main_time),
            " | ", fmt_time(r.pr_time),
            " | ", fmt_ratio(r.time_ratio, band),
            " | ", fmt_pct(r.noise),
            " | ", fmt_ratio(memratio, band), " |"
        )
    end
    return String(take!(io))
end

function group_of(name)
    is_ad(name) && return String(split(name, " / ")[end])
    return "Log density"
end

function buckets(band)
    lo, hi = 100 * (1 - band), 100 * (1 + band)
    return [
        ("🟢 <50%", 50.0),
        ("🟢 50–75%", 75.0),
        ("🟢 75–$(round(Int, lo))%", lo),
        ("⚪ $(round(Int, lo))–$(round(Int, hi))%", hi),
        ("🔴 $(round(Int, hi))–125%", 125.0),
        ("🔴 125–150%", 150.0),
        ("🔴 >150%", Inf),
    ]
end

function bucket_index(pct, edges)
    for (i, (_, edge)) in enumerate(edges)
        pct < edge && return i
    end
    return length(edges)
end

function summary_table(rows, band)
    edges = buckets(band)
    counts = Dict{String, Vector{Int}}()
    for r in rows
        isnan(r.time_ratio) && continue
        v = get!(counts, group_of(r.name), zeros(Int, length(edges)))
        v[bucket_index(100 * r.time_ratio, edges)] += 1
    end
    groups = String[]
    for g in GROUP_ORDER
        haskey(counts, g) && push!(groups, g)
    end
    for g in sort(collect(keys(counts)))
        g in groups || push!(groups, g)
    end
    io = IOBuffer()
    println(io, "| Group | ", join(first.(edges), " | "), " |")
    println(io, "|---|", repeat("---|", length(edges)))
    for g in groups
        cells = [c == 0 ? "·" : string(c) for c in counts[g]]
        println(io, "| ", g, " | ", join(cells, " | "), " |")
    end
    return String(take!(io))
end

## The signal a flat table cannot show. Unrelated components share no cause,
## so one factor applied to all of them is the environment, not the diff.
function spread_note(rows, band)
    ratios = [r.time_ratio for r in rows if !isnan(r.time_ratio)]
    length(ratios) < 4 && return ""
    lo, hi = minimum(ratios), maximum(ratios)
    med = percentile(ratios, 0.5)
    io = IOBuffer()
    println(
        io, "\nAcross all ", length(ratios), " benchmarks the ratio runs ",
        round(lo; digits = 2), "× to ", round(hi; digits = 2),
        "×, median ", round(med; digits = 2), "×."
    )
    if hi / lo < 1.15 && abs(med - 1) > band
        println(io, "\n> [!WARNING]")
        println(
            io, "> Every benchmark moved by close to the same factor. ",
            "Unrelated components do not share a cause, so a uniform shift ",
            "points at the environment rather than at this diff. Read the ",
            "rows below with that in mind."
        )
    end
    return String(take!(io))
end

## A pair is benchmarked only when its gradient smoke test passes, so a pair
## present on one revision only is a differentiability change.
function differentiability_note(rows)
    changed = filter(r -> is_ad(r.name) && r.status !== :both, rows)
    isempty(changed) && return ""
    io = IOBuffer()
    println(
        io, "\n**Differentiability changed.** ",
        "A pair is benchmarked only when its gradient smoke test passes, ",
        "so these components started or stopped differentiating:\n"
    )
    for r in sort(changed; by = x -> x.name)
        println(
            io, "- `", r.name, "` — ",
            r.status === :new ? "now differentiates" : "no longer differentiates"
        )
    end
    return String(take!(io))
end

# ---- main ------------------------------------------------------------------

function main(args)
    length(args) == 5 ||
        error("usage: comment.jl <dir> <pkg> <base-rev> <head-rev> <out.md>")
    dir, pkg, base_rev, head_rev, out_file = args
    main_arm = load_arm(dir, pkg, base_rev)
    pr_arm = load_arm(dir, pkg, head_rev)
    rows = filter(r -> !is_load(r.name), build_rows(pr_arm, main_arm))
    noise = percentile(filter(!isnan, [r.noise for r in rows]), BAND_QUANTILE)
    band = isnan(noise) ? BAND_UNMEASURED : clamp(noise, BAND_FLOOR, BAND_CEILING)

    all_sorted = sort(rows; by = sort_key, rev = true)
    eval_rows = filter(r -> !is_ad(r.name), all_sorted)
    ad_rows = filter(r -> is_ad(r.name), all_sorted)

    io = IOBuffer()
    println(io, COMMENT_MARKER)
    println(io, "## Benchmark comparison vs `main`\n")
    println(
        io,
        "Minimum time per call, per model component. Both revisions were ",
        "benchmarked by AirspeedVelocity in one job on one machine, so a ",
        "ratio here compares two revisions and not two runners.\n"
    )
    println(
        io,
        "Buckets are **PR time as a % of `main`, so lower is faster** ",
        "(🟢 faster, ⚪ within ", round(Int, 100band), "%, 🔴 slower)."
    )
    if isnan(noise)
        println(
            io,
            "\nThat band is provisional. This run carried no sample spread ",
            "to measure, so it is a placeholder and not a calibration.\n"
        )
    elseif noise > BAND_CEILING
        println(
            io,
            "\n> [!WARNING]\n",
            "> The band is capped. ", round(Int, 100BAND_QUANTILE),
            "% of benchmarks varied by less than ", fmt_pct(noise),
            " across their own samples, which is wider than this comment ",
            "will call neutral. Treat everything below as unresolved.\n"
        )
    else
        println(
            io,
            "\nThe band is measured from this run: ",
            round(Int, 100BAND_QUANTILE), "% of benchmarks had a sample ",
            "interquartile range under ", fmt_pct(noise), " of their median.\n"
        )
    end
    println(
        io,
        "> [!NOTE]\n",
        "> That band is a lower bound on what this harness can resolve. ",
        "AirspeedVelocity times each revision once, so the spread is ",
        "dispersion within a revision's own samples, not drift between the ",
        "two revisions, which are separated by the time it takes to compile ",
        "the second one's gradients.\n"
    )
    println(io, "Counts of benchmarks per bucket:\n")
    print(io, summary_table(rows, band))
    print(io, spread_note(rows, band))
    print(io, differentiability_note(rows))

    println(
        io, "\n<details><summary><b>Log density</b> — ", length(eval_rows),
        " benchmarks (by time change)</summary>\n"
    )
    print(io, render_table(eval_rows, band))
    println(io, "\n</details>")

    println(
        io, "\n<details><summary><b>AD gradients</b> — ", length(ad_rows),
        " benchmarks (by time change)</summary>\n"
    )
    print(io, render_table(ad_rows, band))
    println(io, "\n</details>")

    println(
        io,
        "\n`spread` is the interquartile range of that benchmark's own ",
        "samples, as a fraction of its median."
    )
    write(out_file, String(take!(io)))
    return println("Wrote comparison comment to ", out_file)
end

main(ARGS)
