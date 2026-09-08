## Julia boxes any local that a closure captures and something later
## reassigns. A boxed local is type-unstable at every use, which costs
## Mooncake a dictionary lookup per call site on every gradient. Results stay
## correct, so no other test can see it: the only symptom is roughly twice the
## gradient time.
##
## If this fails, look for a variable assigned inside an `if` (or rebound in a
## loop) that a nearby comprehension, `map`, or `do` block also reads, and
## assign the final name once instead. See the "Closures in model code"
## section of `docs/src/contributing.md`.
@testitem "no boxed captures in the model bodies" begin
    using BVDOutbreakSize

    M = BVDOutbreakSize
    models_dir = joinpath(pkgdir(M), "src", "models")

    ## `@model` puts the model body in a gensym-named evaluator method, not in
    ## the user-facing constructor, so `methods(confirmed_cases_model)` never
    ## reaches what AD differentiates. Walk every binding in the module,
    ## gensym'd ones included, and keep the methods defined under `src/models/`.
    offenders = Tuple{String, Int}[]
    for nm in names(M; all = true, imported = false)
        isdefined(M, nm) || continue
        f = getfield(M, nm)
        f isa Function || continue
        for m in methods(f)
            parentmodule(m) === M || continue
            startswith(string(m.file), models_dir) || continue
            ci = try
                Base.uncompressed_ast(m)
            catch
                nothing
            end
            ci === nothing && continue
            boxes = count(x -> occursin("Core.Box", string(x)), ci.code)
            boxes > 0 && push!(offenders, (string(nm), boxes))
        end
    end

    @test isempty(offenders)
    isempty(offenders) ||
        @info "boxed captures" offenders
end
