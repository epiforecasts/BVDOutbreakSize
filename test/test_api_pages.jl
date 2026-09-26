## Every public docstring is on an API page. A docstring an `@autodocs` block
## does not list is left out of the report, and an `@ref` to it from a page
## that is rendered fails the docs build's combine step on a missing
## cross-reference (`production_joint` from `recovery.jl`).

@testitem "API pages: every public docstring's file is listed" tags = [
    :quality,
] begin
    using BVDOutbreakSize
    M = BVDOutbreakSize
    root = pkgdir(M)
    src = joinpath(root, "src")
    libdir = joinpath(root, "docs", "src", "lib")
    listed = String[]
    for f in readdir(libdir; join = true)
        endswith(f, ".md") || continue
        for block in eachmatch(r"```@autodocs(.*?)```"s, read(f, String))
            for m in eachmatch(r"\"([^\"]+\.jl)\"", block[1])
                push!(listed, m[1])
            end
        end
    end
    ispub = isdefined(Base, :ispublic) ? Base.ispublic : Base.isexported
    unlisted = Set{String}()
    for (binding, multidoc) in Base.Docs.meta(M)
        binding.mod === M && ispub(M, binding.var) || continue
        for doc in values(multidoc.docs)
            path = get(doc.data, :path, nothing)
            path === nothing && continue
            rel = relpath(path, src)
            any(p -> endswith(rel, p), listed) || push!(unlisted, rel)
        end
    end
    @test isempty(unlisted)
end
