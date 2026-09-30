# Headline cumulative-infections estimate of every results release.
#
# Run from the worktree root:
#
#     julia --project=. paper/scripts/release_estimates.jl
#
# Uses only the Julia standard library and shells out to `gh`. For each
# release tag in paper/data/code_size.csv it downloads three assets of
# the matching `results-<tag>` GitHub release into
# paper/data/release/results-<tag>/ (git-ignored) and reads:
#
#   observations.toml      -> `as_of_date`, the data cut-off
#   posterior_summary.csv  -> the `Lower 90%` and `Upper 90%` columns of
#                             the headline row (`C_T`, or
#                             `cumulative_cases` before v1.4.0)
#   posterior_draws.csv    -> the median of the same column
#
# This is the definition scripts/refresh_releases.jl uses for
# data/released_estimates.csv, which the repository stopped refreshing
# after v1.6.0. posterior_summary.csv carries no median column, so the
# median is computed from the thinned draws the release itself
# publishes. A release that does not exist on GitHub keeps its row with
# empty numbers. Nothing is estimated.
#
# Output: paper/data/release_estimates.csv with columns tag, tag_date,
# cutoff, model_version, median, lower90, upper90, source_asset,
# quantity.

using Dates
using Statistics
using TOML

const REPO = normpath(joinpath(@__DIR__, "..", ".."))
const DATA = joinpath(REPO, "paper", "data")
const CACHE = joinpath(DATA, "release")
const GH_REPO = "epiforecasts/BVDOutbreakSize"
const ASSETS = (
    "observations.toml", "posterior_summary.csv",
    "posterior_draws.csv",
)
# Headline column name by era; the first present wins.
const HEADLINE = ("C_T", "cumulative_cases")

splitlines(s) = filter(!isempty, split(s, '\n'))
vernum(t) = VersionNumber(lstrip(t, ['v', 'V']))

# Model version by tag, as the paper names the three eras.
function model_version(t)
    v = vernum(t)
    v < v"1.4.0" && return "closed-form"
    v < v"2.0.0" && return "renewal"
    return "provincial"
end

# gh -R epiforecasts/BVDOutbreakSize release view results-<tag> --json assets
release_exists(tag) = success(
    pipeline(
        Cmd(
            [
                "gh", "-R", GH_REPO,
                "release", "view", tag, "--json", "assets",
            ]
        );
        stdout = devnull, stderr = devnull
    )
)

# gh -R epiforecasts/BVDOutbreakSize release download results-<tag>
#    -p <asset> -D paper/data/release/results-<tag>/
function fetch(tag, asset)
    dir = joinpath(CACHE, tag)
    dest = joinpath(dir, asset)
    isfile(dest) && return dest
    mkpath(dir)
    run(
        pipeline(
            Cmd(
                [
                    "gh", "-R", GH_REPO, "release", "download", tag,
                    "-p", asset, "-D", dir, "--clobber",
                ]
            ); stdout = devnull
        )
    )
    return isfile(dest) ? dest : nothing
end

# Minimal CSV reader: no quoted commas in these assets.
function read_csv(path)
    lines = splitlines(read(path, String))
    header = String.(split(lines[1], ','))
    rows = [String.(split(l, ',')) for l in lines[2:end]]
    return header, rows
end

function summary_row(path)
    header, rows = read_csv(path)
    lo = findfirst(==("Lower 90%"), header)
    hi = findfirst(==("Upper 90%"), header)
    for q in HEADLINE
        i = findfirst(r -> r[1] == q, rows)
        i === nothing && continue
        return (
            quantity = q, lower90 = parse(Float64, rows[i][lo]),
            upper90 = parse(Float64, rows[i][hi]),
        )
    end
    error("no headline row in $path; rows are $(first.(rows))")
end

function draws_median(path, quantity)
    header, rows = read_csv(path)
    j = findfirst(==(quantity), header)
    j === nothing && error("no $quantity column in $path")
    return median(parse(Float64, r[j]) for r in rows)
end

csvcell(x) = occursin(r"[,\"\n]", string(x)) ?
    "\"" * replace(string(x), "\"" => "\"\"") * "\"" : string(x)
fmt(x) = x === nothing ? "" : string(round(Int, x))

# Tags and tag dates come from code_size.csv so the two tables agree.
size_header, size_rows = read_csv(joinpath(DATA, "code_size.csv"))
tagcol = findfirst(==("tag"), size_header)
datecol = findfirst(==("date"), size_header)
tags = [(r[tagcol], r[datecol]) for r in size_rows]

out = map(tags) do (t, tag_date)
    rtag = "results-" * t
    if !release_exists(rtag)
        println(rpad(t, 8), "no release ", rtag)
        return (t, tag_date, "", model_version(t), "", "", "", "", "")
    end
    files = Dict(a => fetch(rtag, a) for a in ASSETS)
    missing_assets = [a for a in ASSETS if files[a] === nothing]
    if !isempty(missing_assets)
        println(rpad(t, 8), "missing ", join(missing_assets, ", "))
        return (t, tag_date, "", model_version(t), "", "", "", "", "")
    end
    cutoff = string(TOML.parsefile(files["observations.toml"])["as_of_date"])
    s = summary_row(files["posterior_summary.csv"])
    med = draws_median(files["posterior_draws.csv"], s.quantity)
    println(
        rpad(t, 8), cutoff, "  ", s.quantity, "  ", fmt(med),
        " (", fmt(s.lower90), ", ", fmt(s.upper90), ")"
    )
    (
        t, tag_date, cutoff, model_version(t), fmt(med), fmt(s.lower90),
        fmt(s.upper90), "posterior_summary.csv;posterior_draws.csv", s.quantity,
    )
end

open(joinpath(DATA, "release_estimates.csv"), "w") do io
    println(
        io, "tag,tag_date,cutoff,model_version,median,lower90,upper90,",
        "source_asset,quantity"
    )
    for r in out
        println(io, join(csvcell.(r), ","))
    end
end
println("release_estimates.csv: ", length(out), " rows")
