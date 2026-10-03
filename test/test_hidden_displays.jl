## Literate runs a code chunk as one block and shows only its last value, so
## a hidden display line followed by more code in the same chunk never
## reaches the page (#735). Each hidden display ends its chunk (`#-`).

@testitem "report pages: one hidden display per code chunk" tags = [
    :quality,
] begin
    using BVDOutbreakSize
    pages = joinpath(pkgdir(BVDOutbreakSize), "docs", "pages")
    ## A line Literate renders as markdown, or a chunk break. `##` is a code
    ## comment and `#src` lines are dropped before parsing.
    breaks(l) = occursin(r"^\h*#-", l) ||
        occursin(r"^\h*#md ", l) ||
        occursin(r"^\h*#( |$)", l)
    ## A hidden line that shows a value rather than assigning or loading one.
    shows(l) = occursin(r"#hide$", l) &&
        !occursin(r";\s*#hide$", l) &&
        !occursin(r"^\h*(using|import|include)\b", l) &&
        !occursin(r"^[^=(]*[^=!<>]=[^=]", l)
    shared = String[]
    for (dir, _, files) in walkdir(pages), f in files
        endswith(f, ".jl") || continue
        path = joinpath(dir, f)
        pending = nothing
        for (i, l) in enumerate(eachline(path))
            l = rstrip(l)
            (isempty(l) || endswith(l, "#src")) && continue
            if breaks(l)
                pending = nothing
                continue
            end
            pending === nothing ||
                push!(shared, "$(relpath(path, pages)):$pending")
            pending = shows(l) ? i : nothing
        end
    end
    @test isempty(shared)
end

## A `<details>` fold meant for construction code must not end up wrapping
## the display it was supposed to reveal: a fold whose only code line is a
## hidden display leaves that figure or table collapsed by default (#1056).

@testitem "forecast pages: a shown result is not left folded" tags = [
    :quality,
] begin
    using BVDOutbreakSize
    forecasts = joinpath(pkgdir(BVDOutbreakSize), "docs", "pages", "forecasts")
    opens(l) = occursin(r"^\h*#md # <details>", l)
    closes(l) = occursin(r"^\h*#md # </details>", l)
    is_code(l) = !isempty(l) && !occursin(r"^\h*#", l)
    ## A hidden line that shows a value rather than assigning or loading one.
    shows(l) = occursin(r"#hide$", l) &&
        !occursin(r";\s*#hide$", l) &&
        !occursin(r"^\h*(using|import|include)\b", l) &&
        !occursin(r"^[^=(]*[^=!<>]=[^=]", l)
    folded = String[]
    for (dir, _, files) in walkdir(forecasts), f in files
        endswith(f, ".jl") || continue
        path = joinpath(dir, f)
        open_at_depth1 = false
        saw_code = false
        only_shows = true
        for (i, l) in enumerate(eachline(path))
            l = rstrip(l)
            if opens(l)
                open_at_depth1 = true
                saw_code = false
                only_shows = true
                continue
            end
            if closes(l)
                if open_at_depth1 && saw_code && only_shows
                    push!(folded, "$(relpath(path, forecasts)):$i")
                end
                open_at_depth1 = false
                continue
            end
            (open_at_depth1 && is_code(l)) || continue
            saw_code = true
            only_shows &= shows(l)
        end
    end
    @test isempty(folded)
end
