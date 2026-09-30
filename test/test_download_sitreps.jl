## Tests for the network-free parts of scripts/download_sitreps.jl. The
## script's `main` is guarded by a PROGRAM_FILE check, so including it here
## defines the helpers without making a request.

@testitem "download_sitreps parse_args reads only, limit and outdir" begin
    using BVDOutbreakSize: BVDOutbreakSize
    include(
        joinpath(pkgdir(BVDOutbreakSize), "scripts", "download_sitreps.jl")
    )

    d = parse_args(String[])
    @test d.only_numbers === nothing
    @test d.limit == typemax(Int)
    @test isempty(d.rest)

    a = parse_args(["--limit", "3", "out", "--only", "74,9"])
    @test a.only_numbers == ["074", "009"]
    @test a.limit == 3
    @test a.rest == ["out"]

    @test_throws ErrorException parse_args(["--limit"])
end

@testitem "download_sitreps missing_numbers takes the newest N, ascending" begin
    using BVDOutbreakSize: BVDOutbreakSize
    include(
        joinpath(pkgdir(BVDOutbreakSize), "scripts", "download_sitreps.jl")
    )

    mktempdir() do dir
        touch(joinpath(dir, "SitRep_MVE_133_2026.pdf"))
        numbers = ["130", "131", "132", "133", "134"]
        @test missing_numbers(dir, numbers, typemax(Int)) ==
            ["130", "131", "132", "134"]
        @test missing_numbers(dir, numbers, 2) == ["132", "134"]
        @test isempty(missing_numbers(dir, numbers, 0))
    end
end

@testitem "download_sitreps gaps_suspected skips a complete cache" begin
    using BVDOutbreakSize: BVDOutbreakSize
    include(
        joinpath(pkgdir(BVDOutbreakSize), "scripts", "download_sitreps.jl")
    )

    mktempdir() do dir
        ## An empty listing always needs the posts.
        @test gaps_suspected(dir, String[])
        for n in 1:5
            n == 3 && continue
            touch(joinpath(dir, "SitRep_MVE_$(lpad(n, 3, '0'))_2026.pdf"))
        end
        ## 003 is missing below the newest report on disk.
        @test gaps_suspected(dir, ["005"])
        touch(joinpath(dir, "SitRep_MVE_003_2026.pdf"))
        @test !gaps_suspected(dir, ["005"])
        ## A listed report beyond the disk is fetched by the media walk, and
        ## the reports never published are not gaps.
        @test !gaps_suspected(dir, ["004"])
    end
end

@testitem "download_sitreps fetch_mirror follows LFS and keeps only PDFs" begin
    using BVDOutbreakSize: BVDOutbreakSize
    include(
        joinpath(pkgdir(BVDOutbreakSize), "scripts", "download_sitreps.jl")
    )

    pointer = "version https://git-lfs.github.com/spec/v1\noid sha256:0\n"
    ## A fake `get` that serves `files[url]` with status 200, else 404.
    function fake_get(files)
        return (url, dest) -> begin
            haskey(files, url) || return 404
            write(dest, files[url])
            return 200
        end
    end
    name = "data/insp_sitrep/raw/SitRep_MVE_015-2026.pdf"
    mktempdir() do dir
        dest = joinpath(dir, "SitRep_MVE_015_2026.pdf")
        files = Dict(
            MIRROR_RAW * name => pointer, MIRROR_LFS * name => "%PDF-1.7 x"
        )
        @test fetch_mirror("015", dest; get = fake_get(files))
        @test read(dest, String) == "%PDF-1.7 x"

        ## A page that is not a PDF is discarded, not kept as the report.
        files = Dict(MIRROR_RAW * name => "<html>404</html>")
        @test !fetch_mirror("015", dest; get = fake_get(files))
        @test !isfile(dest)

        @test !fetch_mirror("015", dest; get = fake_get(Dict()))
        @test !isfile(dest)
    end
end
