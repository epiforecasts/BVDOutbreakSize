## Tests for the release helper (`scripts/release_notes.jl`): pulling a
## version's section out of news.md, refusing to release when the section and
## Project.toml disagree, and opening the next version's section.

@testitem "release_notes reads the newest section" begin
    include(joinpath(@__DIR__, "..", "scripts", "release_notes.jl"))

    news = """
    # News

    Release notes for BVDOutbreakSize.

    ## v2.1.0

    Changes since v2.0.0.

    ### Model

    - A model change (#756).

    ## v2.0.0

    Changes since v1.18.0

    - An older change.
    """

    sections = news_sections(news)
    @test [s.version for s in sections] == ["2.1.0", "2.0.0"]
    notes = release_notes(news, "2.1.0")
    @test startswith(notes, "Changes since v2.0.0.")
    @test occursin("- A model change (#756).", notes)
    ## The section stops at the next version heading, not at the `###`
    ## subsection heading inside it.
    @test occursin("### Model", notes)
    @test !occursin("An older change", notes)

    ## Releasing a version whose notes were never written is the mistake
    ## worth catching, so it is an error rather than an empty body.
    @test_throws ErrorException release_notes(news, "2.1.1")
    @test_throws ErrorException release_notes("# News\n\nNothing.\n", "2.1.0")
end

@testitem "version arithmetic and the Project.toml bump" begin
    include(joinpath(@__DIR__, "..", "scripts", "release_notes.jl"))

    @test next_version("2.1.0", "patch") == "2.1.1"
    @test next_version("2.1.0", "minor") == "2.2.0"
    @test next_version("2.1.0", "major") == "3.0.0"
    @test next_version("1.13.2", "patch") == "1.13.3"
    @test_throws ErrorException next_version("2.1.0", "huge")
    @test_throws ErrorException next_version("2.1", "patch")

    toml = """
    name = "BVDOutbreakSize"
    uuid = "7c2a1f8e-9d4b-4a16-bf3e-1f5c8c3a90e2"
    version = "2.1.0"

    [compat]
    CSV = "0.10, 1.0"
    """
    @test project_version(toml) == "2.1.0"
    bumped = bump_project(toml, "2.1.0", "2.1.1")
    @test project_version(bumped) == "2.1.1"
    ## The compat bounds are versions too; only the package's own moves.
    @test occursin("CSV = \"0.10, 1.0\"", bumped)
    @test_throws ErrorException bump_project(toml, "9.9.9", "10.0.0")
end

@testitem "open_section starts the next version above the released one" begin
    include(joinpath(@__DIR__, "..", "scripts", "release_notes.jl"))

    news = "# News\n\nPreamble.\n\n## v2.1.0\n\nChanges since v2.0.0.\n"
    opened = open_section(news, "2.1.0", "2.1.1")
    sections = news_sections(opened)
    @test [s.version for s in sections] == ["2.1.1", "2.1.0"]
    ## The new section is empty apart from the line saying what it follows,
    ## and the released one is untouched.
    @test first(sections).body == "Changes since v2.1.0."
    @test last(sections).body == release_notes(news, "2.1.0")
    @test startswith(opened, "# News\n\nPreamble.\n")

    ## Only the newest section is ever released, so once the next one is
    ## open the released version can no longer be cut again by mistake.
    @test_throws ErrorException release_notes(opened, "2.1.0")

    ## Opening against the wrong released version would bury a section.
    @test_throws ErrorException open_section(news, "2.0.0", "2.0.1")
end

@testitem "CITATION.cff follows the release, not the next version" begin
    include(joinpath(@__DIR__, "..", "scripts", "release_notes.jl"))

    cff = """
    cff-version: 1.2.0
    title: BVDOutbreakSize
    version: 2.1.0
    date-released: "2026-09-17"
    """
    ## It cites the released software, so after cutting v2.1.1 it says
    ## 2.1.1, not the 2.1.2 that Project.toml moves on to.
    bumped = bump_citation(cff, "2.1.1", "2026-09-21")
    @test occursin("version: 2.1.1", bumped)
    @test occursin("date-released: \"2026-09-21\"", bumped)
    @test !occursin("2026-09-17", bumped)
    @test occursin("cff-version: 1.2.0", bumped)

    @test_throws ErrorException bump_citation(
        "title: BVDOutbreakSize\n", "2.1.1", "2026-09-21"
    )
end

@testitem "CRLF files are read the same as LF ones" begin
    include(joinpath(@__DIR__, "..", "scripts", "release_notes.jl"))

    ## A Windows checkout has CRLF line endings, and the repository carries
    ## no `.gitattributes` to stop that. Every function here matches on line
    ## structure, so a stray `\r` before the line end silently found no
    ## version sections at all and the release read as a file with no notes.
    lf = """
    # News

    ## v2.1.1

    Changes since v2.1.0.

    ### Infrastructure

    - A change (#1).

    ## v2.1.0

    Changes since v2.0.0.
    """
    crlf = replace(lf, "\n" => "\r\n")

    @test [s.version for s in news_sections(crlf)] ==
        [s.version for s in news_sections(lf)]
    ## The body is a slice of the text it was given, so CRLF in means CRLF
    ## out. That is right: the notes go to the release verbatim. What has to
    ## hold is that the same content is found, not that the bytes match.
    @test replace(release_notes(crlf, "2.1.1"), "\r\n" => "\n") ==
        release_notes(lf, "2.1.1")
    @test occursin("- A change (#1).", release_notes(crlf, "2.1.1"))

    opened = open_section(crlf, "2.1.1", "2.1.2")
    @test [s.version for s in news_sections(opened)] ==
        ["2.1.2", "2.1.1", "2.1.0"]

    ## The bump leaves the line endings it found rather than rewriting the
    ## line as a bare newline and mixing the two.
    cff = replace(
        "cff-version: 1.2.0\nversion: 2.1.0\ndate-released: \"2026-09-17\"\n",
        "\n" => "\r\n"
    )
    bumped = bump_citation(cff, "2.1.1", "2026-09-21")
    @test occursin("version: 2.1.1\r\n", bumped)
    @test occursin("date-released: \"2026-09-21\"\r\n", bumped)
    @test !occursin("version: 2.1.1\n\n", bumped)

    toml = replace("name = \"X\"\nversion = \"2.1.0\"\n", "\n" => "\r\n")
    @test project_version(bump_project(toml, "2.1.0", "2.1.1")) == "2.1.1"
end

@testitem "the shipped news.md, Project.toml and CITATION.cff agree" begin
    include(joinpath(@__DIR__, "..", "scripts", "release_notes.jl"))

    ## The repository's own files, so a release section that drifts from the
    ## version is caught here rather than by the release workflow.
    version = project_version(read(PROJECT_PATH, String))
    news = read(NEWS_PATH, String)
    @test !isempty(release_notes(news, version))

    ## CITATION.cff cites the last release, so it sits at or below the
    ## version the repository is working towards, never above it.
    cff = read(CITATION_PATH, String)
    cited = match(r"^version:[ \t]*(\S+)"m, cff)
    @test cited !== nothing
    @test VersionNumber(cited.captures[1]) <= VersionNumber(version)
end
