# Content-addressed caching for the expensive model fits in the analysis
# report. Each fit is serialised to `<cache_dir>/<key>.jls`, where `key`
# already encodes a content hash of the fit-relevant source and data, so a
# cached fit is reused only when the code and data that produced it are
# unchanged. This lets the docs build load precomputed chains — produced once,
# e.g. by a per-fit CI matrix or on the HPC — instead of refitting every model
# inline, and lets each fit be recomputed independently.

using SHA: SHA256_CTX, update!, digest!, sha256
using Serialization: serialize, deserialize

# Rebuild a deserialised FlexiChain's parameter keys in the current environment
# so a FlexiChains version skew between the fit and render jobs doesn't break
# `chn[:name]` lookups (the chain loads but its keys can't be indexed).
function _is_chain(x)
    hasfield(typeof(x), :_data) && hasfield(typeof(x), :_metadata) &&
        hasfield(typeof(x), :_structures)
end
function _rebuild_varname(vn)
    (APL = parentmodule(typeof(vn));
        APL.VarName{APL.getsym(vn)}(APL.getoptic(vn)))
end
function repair_chain_keys(x)
    if _is_chain(x)
        FC = parentmodule(typeof(x))
        return FC.map_keys(
            k -> k isa FC.Parameter ? FC.Parameter(_rebuild_varname(FC.get_name(k))) :
                 k, x)
    elseif x isa NamedTuple && haskey(x, :chn)
        return merge(x, (; chn = repair_chain_keys(x.chn)))
    end
    return x
end

"SHA-256 hex digest of a file's bytes, or `\"absent\"` when the file is missing."
function file_sha256(path::AbstractString)
    isfile(path) || return "absent"
    return bytes2hex(sha256(read(path)))
end

"""
    tree_sha256(dir; exts = (".csv",), exclude = ()) -> String

SHA-256 hex digest over every file under `dir` (recursively) whose name ends in
one of `exts`, hashing the sorted `(relative path, contents)` pairs so the
result is order-independent and reproducible. Returns `"absent"` for a missing
directory. Files whose basename or path relative to `dir` matches an entry of
`exclude` are skipped, so a directory can hold non-input files (e.g. a
published-estimate overlay) without them contributing to the digest.
"""
function tree_sha256(dir::AbstractString; exts = (".csv",), exclude = ())
    isdir(dir) || return "absent"
    files = String[]
    for (root, _, fs) in walkdir(dir), f in fs

        any(endswith(f, e) for e in exts) || continue
        path = joinpath(root, f)
        (f in exclude || relpath(path, dir) in exclude) && continue
        push!(files, path)
    end
    ctx = SHA256_CTX()
    for path in sort(files)
        update!(ctx, codeunits(relpath(path, dir)))
        update!(ctx, read(path))
    end
    return bytes2hex(digest!(ctx))
end

"""
    content_hash(source_files; data_dir = nothing, data_exclude = (),
                 extra = "", len = 16) -> String

Short content hash combining the digests of each file in `source_files`, the
`tree_sha256` of `data_dir` (when given, skipping any file matched by
`data_exclude`) and an `extra` string (sampler settings, a schema version, ...).
Used to build cache keys so any change to the model source, data or settings
yields a fresh key.
"""
function content_hash(source_files;
        data_dir = nothing, data_exclude = (),
        extra::AbstractString = "", len::Integer = 16)
    ctx = SHA256_CTX()
    for f in source_files
        update!(ctx, codeunits(file_sha256(f)))
    end
    data_dir === nothing ||
        update!(ctx, codeunits(tree_sha256(data_dir; exclude = data_exclude)))
    update!(ctx, codeunits(extra))
    return bytes2hex(digest!(ctx))[1:len]
end

"""
    fit_or_load(key, thunk; cache_dir, refit = false, strict = false) -> Any

Return the cached result at `<cache_dir>/<key>.jls` when it exists and `refit`
is false; otherwise run `thunk()`, serialise its result to that path and return
it. `key` should already encode a content hash (see [`content_hash`](@ref)) so
a stale cache is never silently reused. The result is written to a temporary
file and moved into place, so an interrupted fit does not leave a half-written
cache entry.

Hits and misses are written straight to `stderr` (bypassing any output capture
by `Literate`/`Documenter`), so during a render the cache behaviour is visible
in the job log rather than silent.

`strict = true` turns a miss into an error instead of refitting. The render
loads every fit that the per-fit CI matrix already produced, so a miss there
means the render is looking in the wrong cache directory (or asking for a fit
the matrix never made) — a bug that should fail in seconds naming the key and
the directory, not silently refit the whole report for hours.
"""
function fit_or_load(key::AbstractString, thunk;
        cache_dir::AbstractString, refit::Bool = false, strict::Bool = false)
    mkpath(cache_dir)
    path = joinpath(cache_dir, key * ".jls")
    if !refit && isfile(path)
        println(stderr, "[fit cache] HIT  $key")
        return repair_chain_keys(deserialize(path))
    end
    if strict && !refit
        error("fit cache MISS for key $key in $cache_dir — expected this fit " *
              "to have been produced by the CI fit matrix and downloaded here. " *
              "Refusing to refit inline (strict mode). Check that the render's " *
              "BVD_FIT_CACHE points at the collected fits.")
    end
    println(stderr, "[fit cache] MISS — fitting $key")
    result = thunk()
    tmp = path * ".tmp"
    serialize(tmp, result)
    mv(tmp, path; force = true)
    return result
end
