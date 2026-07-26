using JuliaFormatter

project_root = dirname(dirname(@__DIR__))
dirs_to_check = filter(isdir,
    [joinpath(project_root, d) for d in ("src", "test", "docs", "scripts")])

# Every directory is checked. The previous `all(d -> format(d), dirs)`
# short-circuited on the first unformatted one, so a violation in `src` hid
# whatever `scripts` had, and `verbose = true` then printed every file it
# touched — which says nothing about which file was actually wrong. The verdict
# below still comes from `format(dir)`; the per-file pass underneath it only
# names the offenders so a red run is actionable from the log alone.
unformatted_dirs = String[]
for dir in dirs_to_check
    JuliaFormatter.format(dir; overwrite = false) || push!(unformatted_dirs, dir)
end

if !isempty(unformatted_dirs)
    println(stderr, "Unformatted Julia code under: ",
        join(relpath.(unformatted_dirs, project_root), ", "))
    for dir in unformatted_dirs, (root, _, files) in walkdir(dir), file in files
        endswith(file, ".jl") || continue
        path = joinpath(root, file)
        if !JuliaFormatter.format(path; overwrite = false)
            println(stderr, "  ", relpath(path, project_root))
        end
    end
    println(stderr, "Run scripts/run_formatter.sh to reformat in place.")
end

exit(isempty(unformatted_dirs) ? 0 : 1)
