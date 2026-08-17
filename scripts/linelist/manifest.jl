# Inputs for a line-list refit: where they come from, where they land, and the
# observation manifest built by substituting the case streams into a released
# one.
#
# Shared by scripts/linelist/fit_joint.jl and scripts/linelist/fit_single.jl, so
# the two cannot drift: a fit of one stream and a fit of all of them must be
# reading the same substitution or the comparison between them means nothing.
#
# `include` this, do not run it.

using CSV
using DataFrames
using Dates
using TOML

## The blocks the line list can supply. Everything else in the manifest stays
## as the situation reports gave it.
const LINELIST_BLOCKS = ("confirmed_case_history", "reported_case_history",
    "suspected_daily_history")

## Where the case streams and the reporting triangle are read from. Required
## rather than defaulted: the files hold counts derived from individual patient
## records, so the operator says where they are instead of a default quietly
## finding a directory that happens to exist.
function linelist_input_dir()
    dir = get(ENV, "LINELIST_INPUT_DIR", "")
    isempty(dir) && error("set LINELIST_INPUT_DIR to a directory holding " *
          "linelist_streams.csv and onset_curve_scanned.csv. " *
          "scripts/linelist/README.md gives both schemas, and " *
          "test/fixtures/linelist holds a synthetic pair.")
    isdir(dir) || error("LINELIST_INPUT_DIR is not a directory: $dir")
    return abspath(dir)
end

## Where the manifest and every result are written. Defaults inside `ignore/`,
## which is git-ignored precisely so that line-list-derived output cannot be
## committed to a public repository, and refuses any other path inside the
## repository so a mistyped override cannot land in a tracked directory.
function linelist_output_dir(root)
    ## `rstrip`, because a root written with a trailing separator (or as
    ## `<dir>/test/..`) leaves one behind and a doubled separator makes every
    ## prefix test below silently false.
    base = rstrip(abspath(root), '/')
    dir = rstrip(
        abspath(get(ENV, "LINELIST_OUT_DIR",
            joinpath(base, "ignore", "linelist"))), '/')
    ignore = joinpath(base, "ignore")
    inside = startswith(dir, base * "/")
    if inside && !startswith(dir, ignore * "/") && dir != ignore
        error("LINELIST_OUT_DIR is inside the repository but outside " *
              "ignore/, which is the only ignored output root: $dir. " *
              "This repository is public and these outputs derive from the " *
              "line list.")
    end
    mkpath(dir)
    return dir
end

## The manifest the release being compared against was fitted to. Taken from
## that release rather than from the working tree, since the comparison is
## against that fit and reads its observations.
##
## `LINELIST_RELEASED_MANIFEST` overrides with a local path, which is how the
## synthetic fixture is run without network access.
function released_manifest(out; tag = get(ENV, "LINELIST_RELEASE", ""))
    local_path = get(ENV, "LINELIST_RELEASED_MANIFEST", "")
    if !isempty(local_path)
        isfile(local_path) ||
            error("LINELIST_RELEASED_MANIFEST is not a file: $local_path")
        @info "released manifest, local" path = local_path
        return abspath(local_path)
    end

    repo = "epiforecasts/BVDOutbreakSize"
    if isempty(tag)
        out_tags = read(`gh release list -R $repo -L 200 --json tagName
            --jq ".[].tagName"`, String)
        tags = filter(t -> occursin(r"^results-v[0-9]", t),
            split(strip(out_tags), '\n'))
        isempty(tags) && error("no results release found in $repo")
        tag = first(tags)
    end

    dest = joinpath(out, "released_observations.toml")
    run(`gh release download $tag -R $repo -p observations.toml
        -O $dest --clobber`)
    @info "released manifest, downloaded" tag dest
    return dest
end

## `load_observations` reads the reporting triangle from a fixed filename beside
## the manifest it is handed, and `load_onset_curve` returns an empty no-op
## rather than an error when that file is absent. A fit whose manifest and
## triangle sit in different directories therefore drops the onset stream and
## says nothing, which is why the triangle is copied next to the manifest here
## and its absence is fatal.
function place_onset_curve(input, out)
    src = joinpath(input, "onset_curve_scanned.csv")
    isfile(src) || error("missing $src. It is the onset-by-vintage reporting " *
          "triangle; scripts/linelist/README.md gives the " *
          "schema.")
    dest = joinpath(out, "onset_curve_scanned.csv")
    cp(src, dest; force = true)
    @info "onset triangle placed beside the manifest" dest
    return dest
end

## Read the case streams, failing on a file that is missing or has the wrong
## columns rather than on the first empty selection downstream.
function read_streams(path)
    isfile(path) || error("missing $path. It holds the case streams; " *
          "scripts/linelist/README.md gives the schema.")
    df = CSV.read(path, DataFrame)
    cols = ("stream", "date", "value")
    all(c -> c in names(df), cols) ||
        error("$path needs columns " * join(cols, ", ") * "; found " *
              join(names(df), ", "))
    return df
end

## Build the manifest: the released one with the case streams swapped out.
## Written to a file rather than passed in memory because `load_observations`
## reads a path, and because the manifest is then a recordable artefact of what
## was fitted.
##
## `source` records how the replacement streams were indexed, because that is
## the difference between a series with a ragged edge and one without, and it is
## not recoverable from the numbers alone.
function write_manifest(; released, streams, out,
        source = "DHIS2 case line list, counted by notification date")
    raw = TOML.parsefile(released)
    df = read_streams(streams)

    for block in LINELIST_BLOCKS
        rows = sort(df[df.stream .== block, :], :date)
        isempty(rows) && error("no rows for $block in $streams")
        raw[block] = Dict{String, Any}(
            "dates" => [string(d) for d in rows.date],
            "values" => Int.(rows.value),
            "source" => source
        )
    end

    ## The cut-off scalar has to move with the history it belongs to.
    ##
    ## The released manifest freezes `reported_case_history` at 26 May, with
    ## `[reported_cases] value = 1077` matching its last vintage, because INSP
    ## began revising the suspected count downward from SitRep 013 and the series
    ## stopped being trustworthy. `load_observations` reads that scalar as the
    ## observed cumulative reported total at the cut-off and the models condition
    ## on it directly.
    ##
    ## Replacing the history without replacing the scalar leaves the model told
    ## that the cut-off total is 1077 while the history it fits reaches thirteen
    ## thousand at the same date, and it can only reconcile the two by distorting
    ## ascertainment and the growth path. So the scalar is rewritten from the
    ## replacement history's own last value.
    ##
    ## `confirmed_cases` needs no equivalent: the released manifest carries no
    ## such scalar, so `load_observations` already derives it from whichever
    ## confirmed history it is given.
    let rows = sort(df[df.stream .== "reported_case_history", :], :date)
        raw["reported_cases"] = Dict{String, Any}(
            "value" => Int(rows.value[end]),
            "source" =>
                "Last vintage of the replacement " *
                "reported_case_history above, so the cut-off total and " *
                "the history it comes from cannot disagree."
        )
    end

    ## `confirmed_break_dates` marks days where INSP retrospectively harmonised
    ## its own confirmed series, with the gross counts printed in the situation
    ## reports. Those are artefacts of the situation-report series, not events
    ## in the outbreak: the line list carries its own revisions inside the case
    ## records and has no such steps. Kept, they would ask the model to absorb a
    ## harmonisation the data no longer contains, which is also why
    ## `load_observations` rejects them here. `occupancy_break_dates` stays,
    ## since the isolation stream is still the situation reports'.
    delete!(raw, "confirmed_break_dates")

    ## The cut-off moves to the last line-list day. `load_observations` drops
    ## vintages after the cut-off, so the situation-report streams are read to
    ## the same date and the two manifests stay comparable.
    last_day = maximum(df.date)
    raw["as_of_date"] = string(last_day)

    open(out, "w") do io
        TOML.print(io, raw; sorted = true)
    end
    @info "manifest written" out as_of = string(last_day)
    return out
end
