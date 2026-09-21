# Release notes and the version bump that follows a release, from the one
# place the changes are already written down. `docs/src/news.md` holds a
# section per version, newest first, and a release's GitHub notes have been
# that section copied across by hand. This does the copying, and then opens
# the next version's section so the following change has somewhere to go.
#
#   julia --project=scripts scripts/release_notes.jl notes
#   julia --project=scripts scripts/release_notes.jl next patch
#   julia --project=scripts scripts/release_notes.jl open patch
#
# `notes` prints the release body for the version in `Project.toml`. `next`
# prints the version a bump of that size would give. `open` writes both
# files, bumping `Project.toml` and inserting the new section above the one
# just released.
#
# The functions are pure and take the file contents, so `test/test_release_notes.jl`
# exercises them without touching the repository.

using Dates: today
using TOML: TOML

const NEWS_PATH = joinpath(dirname(@__DIR__), "docs", "src", "news.md")
const PROJECT_PATH = joinpath(dirname(@__DIR__), "Project.toml")
const CITATION_PATH = joinpath(dirname(@__DIR__), "CITATION.cff")

## A version heading and nothing else on the line, so a `###` subsection
## heading inside a release cannot be mistaken for the start of the next one.
const VERSION_HEADING = r"^##[ \t]+v([0-9]+\.[0-9]+\.[0-9]+)[ \t]*$"m

"""
    project_version(toml::AbstractString) -> String

The `version` field of a `Project.toml` given as text.
"""
function project_version(toml::AbstractString)
    v = get(TOML.parse(toml), "version", nothing)
    v === nothing && error("Project.toml carries no version field")
    return String(v)
end

"""
    news_sections(news::AbstractString) -> Vector{@NamedTuple{version::String, body::String}}

Every version section of `news.md`, newest first. `body` is everything under
the heading and above the next one, stripped of surrounding blank lines.
"""
function news_sections(news::AbstractString)
    starts = collect(eachmatch(VERSION_HEADING, news))
    isempty(starts) && return NamedTuple{(:version, :body), Tuple{String, String}}[]
    out = NamedTuple{(:version, :body), Tuple{String, String}}[]
    for (i, m) in pairs(starts)
        from = m.offset + ncodeunits(m.match)
        to = i < length(starts) ? prevind(news, starts[i + 1].offset) : lastindex(news)
        push!(
            out, (
                version = String(m.captures[1]),
                body = strip(news[from:to]),
            )
        )
    end
    return out
end

"""
    release_notes(news, version) -> String

The `news.md` section for `version`. Errors when the newest section is for
some other version, which means the release being cut and the notes written
for it have come apart.
"""
function release_notes(news::AbstractString, version::AbstractString)
    sections = news_sections(news)
    isempty(sections) && error("news.md carries no version sections")
    top = first(sections)
    top.version == version || error(
        "Project.toml is at $version but the newest news.md section is " *
            "v$(top.version). Write the notes for $version, or bump the " *
            "version to match, before releasing."
    )
    isempty(top.body) && error(
        "the v$version section of news.md is empty; there is nothing to " *
            "release"
    )
    return top.body
end

"""
    next_version(version, kind) -> String

The version a `"major"`, `"minor"` or `"patch"` bump of `version` gives.
"""
function next_version(version::AbstractString, kind::AbstractString)
    parts = tryparse.(Int, split(version, "."))
    (length(parts) == 3 && !any(isnothing, parts)) ||
        error("cannot parse version $version as major.minor.patch")
    major, minor, patch = parts
    kind == "major" && return "$(major + 1).0.0"
    kind == "minor" && return "$major.$(minor + 1).0"
    kind == "patch" && return "$major.$minor.$(patch + 1)"
    return error("unknown bump $kind; expected major, minor or patch")
end

"""
    bump_project(toml, old, new) -> String

`Project.toml` text with the `version` field moved from `old` to `new`.
"""
function bump_project(
        toml::AbstractString, old::AbstractString, new::AbstractString
    )
    ## A plain string rather than a pattern. The version is the only
    ## `version = "..."` line in the file; the compat entries are keyed by
    ## package name.
    target = "version = \"$old\""
    occursin(target, toml) ||
        error("Project.toml does not carry $target")
    return replace(toml, target => "version = \"$new\""; count = 1)
end

"""
    bump_citation(cff, released, date) -> String

`CITATION.cff` text with the version and release date moved to the release
just cut. It cites the released software, so it follows the release rather
than the version `Project.toml` moves on to.
"""
function bump_citation(
        cff::AbstractString, released::AbstractString,
        date::AbstractString
    )
    version = Regex("^version:[ \\t]*.*\$", "m")
    released_on = Regex("^date-released:[ \\t]*.*\$", "m")
    occursin(version, cff) ||
        error("CITATION.cff carries no version field")
    occursin(released_on, cff) ||
        error("CITATION.cff carries no date-released field")
    cff = replace(cff, version => "version: $released"; count = 1)
    return replace(
        cff, released_on => "date-released: \"$date\""; count = 1
    )
end

"""
    open_section(news, released, opening) -> String

`news.md` with an empty section for `opening` inserted above the section for
`released`, so the next merged change has a heading to write under.
"""
function open_section(
        news::AbstractString, released::AbstractString,
        opening::AbstractString
    )
    m = match(VERSION_HEADING, news)
    m === nothing && error("news.md carries no version sections")
    String(m.captures[1]) == released || error(
        "expected the newest news.md section to be v$released, found " *
            "v$(m.captures[1])"
    )
    header = "## v$opening\n\nChanges since v$released.\n\n"
    return news[1:prevind(news, m.offset)] * header * news[m.offset:end]
end

## --- Command line -------------------------------------------------------

function main(args)
    isempty(args) && error(
        "usage: release_notes.jl notes | next <bump> | open <bump>"
    )
    command = args[1]
    news = read(NEWS_PATH, String)
    project = read(PROJECT_PATH, String)
    version = project_version(project)
    if command == "notes"
        print(release_notes(news, version))
    elseif command in ("next", "open")
        length(args) == 2 ||
            error("$command needs a bump: major, minor or patch")
        upcoming = next_version(version, args[2])
        if command == "next"
            print(upcoming)
        else
            ## Read the notes first. They validate that the section and the
            ## version agree, and neither file should be written when they
            ## do not.
            release_notes(news, version)
            write(NEWS_PATH, open_section(news, version, upcoming))
            write(PROJECT_PATH, bump_project(project, version, upcoming))
            ## CITATION.cff cites the release just cut, not the version the
            ## repository moves on to, so it takes `version` rather than
            ## `upcoming`. Nothing else was updating it, so it would have
            ## drifted a version further behind with every release.
            write(
                CITATION_PATH,
                bump_citation(
                    read(CITATION_PATH, String), version, string(today())
                )
            )
            print(upcoming)
        end
    else
        error("unknown command $command; expected notes, next or open")
    end
    return nothing
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main(ARGS)
