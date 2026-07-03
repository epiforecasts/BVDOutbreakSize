# Content-addressed caching for the expensive model fits in the analysis
# report. Each fit is serialised to `<cache_dir>/<key>.jls`, where `key`
# already encodes a content hash of the fit-relevant source and data, so a
# cached fit is reused only when the code and data that produced it are
# unchanged. This lets the docs build load precomputed chains — produced once,
# e.g. by a per-fit CI matrix or on the HPC — instead of refitting every model
# inline, and lets each fit be recomputed independently.

using SHA: SHA256_CTX, update!, digest!, sha256
using Serialization: serialize, deserialize

# A FlexiChain serialised by a fit job and deserialised by the render job can
# carry `VarName` keys whose type differs subtly from the ones the render
# environment constructs — a version-skew artefact of Julia's `Serialization`.
# The chain deserialises, but its `OrderedDict` then can't be indexed, so
# `chn[:C_T]` throws a `KeyError` even though the key is present. Rebuilding
# every key's `VarName` in THIS environment restores all access. `map_keys`
# iterates the dict (rather than indexing it) and rehashes into a fresh dict, so
# it works on the otherwise-unindexable chain; only the ~few-hundred keys are
# re-wrapped, the draw matrices are untouched.
function _looks_like_chain(x)
    hasfield(typeof(x), :_data) &&
        hasfield(typeof(x), :_metadata) && hasfield(typeof(x), :_structures)
end

function _rebuild_varname(vn)
    APL = parentmodule(typeof(vn))          # AbstractPPL
    return APL.VarName{APL.getsym(vn)}(APL.getoptic(vn))
end

"""
Re-key any deserialised `FlexiChain` (bare, or nested as the `chn` field of a
frozen fit's named tuple) so its parameter keys are constructed in the current
environment, making `chn[:name]` lookups work regardless of the version that
wrote it.
"""
function repair_chain_keys(x)
    if _looks_like_chain(x)
        FC = parentmodule(typeof(x))        # FlexiChains
        return FC.map_keys(
            k -> k isa FC.Parameter ? FC.Parameter(_rebuild_varname(FC.get_name(k))) : k,
            x)
    elseif x isa NamedTuple && haskey(x, :chn)
        return merge(x, (; chn = repair_chain_keys(x.chn)))
    else
        return x
    end
end

"SHA-256 hex digest of a file's bytes, or `\"absent\"` when the file is missing."
function file_sha256(path::AbstractString)
    isfile(path) || return "absent"
    return bytes2hex(sha256(read(path)))
end

"""
    tree_sha256(dir; exts = (".csv",)) -> String

SHA-256 hex digest over every file under `dir` (recursively) whose name ends in
one of `exts`, hashing the sorted `(relative path, contents)` pairs so the
result is order-independent and reproducible. Returns `"absent"` for a missing
directory.
"""
function tree_sha256(dir::AbstractString; exts = (".csv",))
    isdir(dir) || return "absent"
    files = String[]
    for (root, _, fs) in walkdir(dir), f in fs

        any(endswith(f, e) for e in exts) && push!(files, joinpath(root, f))
    end
    ctx = SHA256_CTX()
    for path in sort(files)
        update!(ctx, codeunits(relpath(path, dir)))
        update!(ctx, read(path))
    end
    return bytes2hex(digest!(ctx))
end

"""
    content_hash(source_files; data_dir = nothing, extra = "", len = 16) -> String

Short content hash combining the digests of each file in `source_files`, the
`tree_sha256` of `data_dir` (when given) and an `extra` string (sampler
settings, a schema version, ...). Used to build cache keys so any change to the
model source, data or settings yields a fresh key.
"""
function content_hash(source_files;
        data_dir = nothing, extra::AbstractString = "", len::Integer = 16)
    ctx = SHA256_CTX()
    for f in source_files
        update!(ctx, codeunits(file_sha256(f)))
    end
    data_dir === nothing || update!(ctx, codeunits(tree_sha256(data_dir)))
    update!(ctx, codeunits(extra))
    return bytes2hex(digest!(ctx))[1:len]
end

"""
    fit_or_load(key, thunk; cache_dir, refit = false) -> Any

Return the cached result at `<cache_dir>/<key>.jls` when it exists and `refit`
is false; otherwise run `thunk()`, serialise its result to that path and return
it. `key` should already encode a content hash (see [`content_hash`](@ref)) so
a stale cache is never silently reused. The result is written to a temporary
file and moved into place, so an interrupted fit does not leave a half-written
cache entry.
"""
function fit_or_load(key::AbstractString, thunk;
        cache_dir::AbstractString, refit::Bool = false)
    mkpath(cache_dir)
    path = joinpath(cache_dir, key * ".jls")
    if !refit && isfile(path)
        @info "fit cache hit" key
        return repair_chain_keys(deserialize(path))
    end
    @info "fit cache miss — fitting" key
    result = thunk()
    tmp = path * ".tmp"
    serialize(tmp, result)
    mv(tmp, path; force = true)
    return result
end
