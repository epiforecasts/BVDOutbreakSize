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

@testitem "the shipped news.md and Project.toml agree" begin
    include(joinpath(@__DIR__, "..", "scripts", "release_notes.jl"))

    ## The repository's own files, so a release section that drifts from the
    ## version is caught here rather than by the release workflow.
    version = project_version(read(PROJECT_PATH, String))
    news = read(NEWS_PATH, String)
    @test !isempty(release_notes(news, version))
end
