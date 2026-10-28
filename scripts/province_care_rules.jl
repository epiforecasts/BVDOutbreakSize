# Rules for the blocks scripts/province_care_manifest.jl emits; see
# "Province isolation occupancy and beds" in data/README.md. Dependency-free
# so the test suite can include it.

"""Last date of the national `isolation_history` in a parsed manifest."""
function last_isolation_date(manifest::AbstractDict)
    block = get(manifest, "isolation_history", nothing)
    block === nothing && error("the manifest has no [isolation_history] block")
    return maximum(Date.(block["dates"]))
end

"""Remove the `misprinted` cells of `column` from `cells`, keyed by (SitRep, province)."""
function drop_misprinted!(cells::AbstractDict, misprinted::AbstractDict, column::Symbol)
    for (sr, prov, col) in keys(misprinted)
        col == column && delete!(cells, (sr, prov))
    end
    return cells
end

"""Drop entries dated after `last`, returning them as (province, date, value)."""
function drop_after!(blocks::AbstractDict, last)
    dropped = Tuple{String, Any, Any}[]
    for prov in sort!(collect(keys(blocks)))
        kept = empty(blocks[prov])
        for (d, v) in blocks[prov]
            d > last ? push!(dropped, (prov, d, v)) : push!(kept, (d, v))
        end
        isempty(kept) ? delete!(blocks, prov) : (blocks[prov] = kept)
    end
    return dropped
end
