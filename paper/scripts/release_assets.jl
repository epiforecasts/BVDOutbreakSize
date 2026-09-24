# Shared helpers for the paper scripts: locating and downloading the assets
# of one `results-*` GitHub release of epiforecasts/BVDOutbreakSize, reading
# the tables the site snapshot carries, and formatting numbers for prose.
#
# Included by paper_numbers.jl and fig_current.jl; not a module.

using CSV
using DataFrames
using Dates
using Statistics: mean, median, quantile
using TOML

const REPO = "epiforecasts/BVDOutbreakSize"

## The release the paper quotes. Newest `results-*` release with the full
## asset set at the time of writing (24 September 2026), marked "Latest" on
## GitHub; the joint summaries it carries are identical to results-2339,
## published thirteen minutes earlier from the next commit on main.
const DEFAULT_TAG = "results-2338"

## Assets every script reads. `site.zip` is the rendered report site; the
## diagnostics and province tables are only published inside it.
const ASSETS = [
    "posterior_draws.csv", "posterior_summary.csv", "stream_draws.csv",
    "stream_estimates.csv", "forecast.csv", "forecast_validation.csv",
    "province_forecast.csv", "onsets_over_time.csv", "observations.toml",
    "site.zip",
]

paper_dir() = normpath(joinpath(@__DIR__, ".."))
repo_dir() = normpath(joinpath(paper_dir(), ".."))
release_dir() = joinpath(paper_dir(), "data", "release")
asset_path(name) = joinpath(release_dir(), name)

"""
    ensure_assets(tag)

Download every asset in `ASSETS` of release `tag` into `release_dir()`
unless it is already there. Needs the `gh` CLI and network access on the
first run only.
"""
function ensure_assets(tag::AbstractString)
    mkpath(release_dir())
    for name in ASSETS
        isfile(asset_path(name)) && continue
        @info "downloading $name from $tag"
        run(
            `gh release download $tag -R $REPO -p $name -D $(release_dir())`
        )
    end
    return nothing
end

"""
    release_meta(tag) -> NamedTuple

Tag, source commit and publication date of the release. Read from GitHub
through `gh release view` once and cached beside the assets as
`release.toml`, so later runs need no network.
"""
function release_meta(tag::AbstractString)
    cache = asset_path("release.toml")
    if isfile(cache)
        m = TOML.parsefile(cache)
        m["tag"] == tag && return (
            tag = m["tag"], commit = m["commit"], date = Date(m["date"]),
        )
    end
    json = read(
        `gh release view $tag -R $REPO --json tagName,publishedAt,body`,
        String
    )
    sha = match(r"for commit ([0-9a-f]{40})", json)
    sha === nothing && error("no source commit in the release body of $tag")
    published = match(r"\"publishedAt\":\"([0-9]{4}-[0-9]{2}-[0-9]{2})", json)
    published === nothing && error("no publishedAt in `gh release view`")
    meta = (
        tag = tag, commit = sha.captures[1],
        date = Date(published.captures[1]),
    )
    open(cache, "w") do io
        println(io, "tag = \"", meta.tag, "\"")
        println(io, "commit = \"", meta.commit, "\"")
        println(io, "date = \"", meta.date, "\"")
    end
    return meta
end

"""
    site_page(path) -> String

One rendered page of the site snapshot, e.g. `"estimates/province.html"`,
read straight out of `site.zip`.
"""
function site_page(path::AbstractString)
    inner = "BVDOutbreakSize/dev/" * path
    return read(`unzip -p $(asset_path("site.zip")) $inner`, String)
end

## Strip tags and unescape the few entities Documenter emits.
function _cell_text(s::AbstractString)
    t = replace(s, r"<[^>]+>" => "")
    t = replace(
        t, "&amp;" => "&", "&lt;" => "<", "&gt;" => ">", "&quot;" => "\"",
        "&#39;" => "'", "&nbsp;" => " ", "\u200b" => ""
    )
    return strip(t)
end

"""
    html_tables(html) -> Vector{Vector{Vector{String}}}

Every `<table>` on a page as rows of cell strings, in page order.
"""
function html_tables(html::AbstractString)
    tables = Vector{Vector{String}}[]
    for t in eachmatch(r"<table.*?</table>"s, html)
        rows = Vector{String}[]
        for r in eachmatch(r"<tr.*?</tr>"s, t.match)
            cells = [
                _cell_text(c.captures[1])
                    for c in eachmatch(r"<t[hd][^>]*>(.*?)</t[hd]>"s, r.match)
            ]
            isempty(cells) || push!(rows, cells)
        end
        push!(tables, rows)
    end
    return tables
end

"""
    table_row(tables, label; nth = 1) -> Vector{String}

The `nth` row across `tables` labelled `label`, from the label onwards.
"""
function table_row(tables, label::AbstractString; nth::Integer = 1)
    seen = 0
    for t in tables, r in t
        ## Documenter prefixes DataFrame tables with a `Row` index column,
        ## so the label may sit in the first or second cell; the cells
        ## from the label onwards are returned either way.
        i = findfirst(==(label), r[1:min(2, end)])
        i === nothing && continue
        seen += 1
        seen == nth && return r[i:end]
    end
    error("no row labelled `$label` (occurrence $nth) in the page tables")
end

"""
    headed_tables(html, level) -> Vector{Pair{String, Vector{Vector{String}}}}

Each heading of `level` (e.g. `"h4"`) paired with the first table that
follows it before the next heading of any level.
"""
function headed_tables(html::AbstractString, level::AbstractString)
    out = Pair{String, Vector{Vector{String}}}[]
    heads = collect(eachmatch(Regex("<$level[^>]*>(.*?)</$level>", "s"), html))
    anyhead = collect(eachmatch(r"<h[1-6][^>]*>"s, html))
    tables = collect(eachmatch(r"<table.*?</table>"s, html))
    for h in heads
        stop = findfirst(m -> m.offset > h.offset, anyhead)
        limit = stop === nothing ? length(html) : anyhead[stop].offset
        t = findfirst(m -> h.offset < m.offset < limit, tables)
        t === nothing && continue
        push!(out, _cell_text(h.captures[1]) => html_tables(tables[t].match)[1])
    end
    return out
end

## --- Formatting for prose ------------------------------------------------

"""
    fmt2(x) -> String

`x` to at most two significant figures, with thousands separators once it
reaches four digits: `20108.7` becomes `"20,000"`, `1.1137` `"1.1"`,
`0.4747` `"0.47"` and `0.5` `"0.5"` (trailing zeros are not padded, since
the summaries are themselves rounded to two decimals).
"""
function fmt2(x::Real; sig::Integer = 2)
    x == 0 && return "0"
    mag = floor(Int, log10(abs(x)))
    decimals = sig - 1 - mag
    if decimals <= 0
        return fmt_count(round(Int, round(x; digits = decimals)))
    end
    return string(round(x; digits = decimals))
end

"""
    fmt_count(n) -> String

An exact integer with thousands separators.
"""
function fmt_count(n::Integer)
    s = string(abs(n))
    parts = String[]
    while length(s) > 3
        pushfirst!(parts, s[(end - 2):end])
        s = s[1:(end - 3)]
    end
    pushfirst!(parts, s)
    return (n < 0 ? "-" : "") * join(parts, ",")
end

fmt_interval(lo, hi; f = fmt2) = string(f(lo), " to ", f(hi))

fmt_pct(x::Real) = fmt2(100x) * "%"

## --- Asset readers ---------------------------------------------------------

read_asset(name) = CSV.read(asset_path(name), DataFrame)

"""
    summary_row(df, quantity) -> DataFrameRow

The `posterior_summary.csv` row for `quantity`, whose columns are
`Quantity, Lower 90%, Lower 60%, Lower 30%, Upper 30%, Upper 60%, Upper 90%`.
"""
function summary_row(df::DataFrame, quantity::AbstractString)
    i = findfirst(==(quantity), df.Quantity)
    i === nothing && error("no `$quantity` row in posterior_summary.csv")
    return df[i, :]
end

"""
    stream_row(df, fit, quantity) -> DataFrameRow

The `stream_estimates.csv` row for one fit and quantity, with columns
`median, lo30, hi30, lo60, hi60, lo90, hi90`.
"""
function stream_row(
        df::DataFrame, fit::AbstractString, quantity::AbstractString
    )
    i = findfirst(
        i -> df.fit[i] == fit && df.quantity[i] == quantity, 1:nrow(df)
    )
    i === nothing && error("no `$fit`/`$quantity` row in stream_estimates.csv")
    return df[i, :]
end

"""
    history(obs, key) -> (dates, values)

A dated cumulative series from `observations.toml`, e.g.
`confirmed_case_history`.
"""
function history(obs::AbstractDict, key::AbstractString)
    h = obs[key]
    return Date.(h["dates"]), Int.(h["values"])
end

"""
    increment_over(dates, values, from, to) -> Union{Int, Missing}

`values[to] - values[from]` when both dates are in the series.
"""
function increment_over(dates, values, from::Date, to::Date)
    i = findfirst(==(from), dates)
    j = findfirst(==(to), dates)
    (i === nothing || j === nothing) && return missing
    return values[j] - values[i]
end

"""
    forecast_draws(df, stream; made, horizon) -> Vector{Float64}

The archived draws of one stream at one origin and horizon from a
`forecast.csv`-schema table.
"""
function forecast_draws(
        df::DataFrame, stream::AbstractString; made::Date, horizon::Integer
    )
    sel = (df.stream .== stream) .& (df.made_date .== made) .&
        (df.horizon .== horizon)
    any(sel) || error("no `$stream` forecast made $made at horizon $horizon")
    return Float64.(df.value[sel])
end

q90(v) = (quantile(v, 0.05), quantile(v, 0.95))
q60(v) = (quantile(v, 0.2), quantile(v, 0.8))
q30(v) = (quantile(v, 0.35), quantile(v, 0.65))

## --- Assets of other releases ---------------------------------------------

tag_asset_path(tag, name) = joinpath(release_dir(), tag, name)

"""
    ensure_tag_assets(tag, names)

Download `names` of release `tag` into `release_dir()/tag/` unless they
are already there. Used for the assets of earlier releases the paper
quotes beside `DEFAULT_TAG`, as fig_first.jl and release_estimates.jl do.
"""
function ensure_tag_assets(tag::AbstractString, names)
    dir = joinpath(release_dir(), tag)
    mkpath(dir)
    for name in names
        isfile(tag_asset_path(tag, name)) && continue
        @info "downloading $name from $tag"
        run(`gh release download $tag -R $REPO -p $name -D $dir`)
    end
    return nothing
end

## --- Cross-release forecast scores ----------------------------------------

"""
    forecast_skill(scores) -> DataFrame

One row per (release, stream) of `data/forecast_scores.csv` that carries
both a joint and a persistence-baseline score at horizon seven: the
joint's CRPS over the baseline's (`skill`), whether the joint's 90%
interval covered the observation (`covered`) and whether the release is
a reconstructed backfill rather than a forecast published at the time
(`backfill`). The rows fig_evaluation.jl draws in its panel A.
"""
function forecast_skill(scores::DataFrame)
    h7 = scores[scores.horizon .== 7, :]
    rows = DataFrame(
        release = String[], date = Date[], stream = String[],
        skill = Float64[], covered = Bool[], backfill = Bool[],
    )
    for r in eachrow(h7[h7.fit .== "joint", :])
        b = findfirst(
            i -> h7.fit[i] == "baseline" && h7.release[i] == r.release &&
                h7.stream[i] == r.stream,
            1:nrow(h7)
        )
        b === nothing && continue
        push!(
            rows,
            (
                r.release, Date(r.made_date), r.stream, r.crps / h7.crps[b],
                r.coverage_90, occursin("(backfill)", r.release),
            )
        )
    end
    return sort!(rows, :date)
end
