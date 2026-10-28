# Cell rules for the blocks scripts/province_care_manifest.jl emits.
# Dependency-free so the test suite can include it.

"""
    is_misprinted(misprinted, sitrep, province, column) -> Bool

Whether the (SitRep, province, CSV column) cell is listed in `misprinted`, so
only that cell is left out of its block.
"""
is_misprinted(misprinted, sitrep, province, column) =
    haskey(misprinted, (sitrep, province, column))

"""
    drop_after!(blocks, last) -> Vector

Remove every `(date, value)` entry dated after `last` from each province's
pairs in `blocks`, and return the removed `(province, date, value)` entries
sorted by province.
"""
function drop_after!(blocks, last)
    dropped = Tuple{String, Any, Any}[]
    for prov in sort!(collect(keys(blocks)))
        pairs = blocks[prov]
        for (d, v) in pairs
            d > last && push!(dropped, (prov, d, v))
        end
        filter!(p -> first(p) <= last, pairs)
    end
    return dropped
end
