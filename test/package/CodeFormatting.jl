@testitem "Code formatting" tags = [:quality] begin
    using Pkg
    formatter_env = joinpath(@__DIR__, "..", "formatter")
    if isdir(formatter_env) && isfile(joinpath(formatter_env, "Project.toml"))
        # Instantiate the formatter environment via a subprocess so the
        # active project of the test process is not mutated (otherwise
        # later @testitems lose access to BVDOutbreakSize).
        #
        # The registry is refreshed first. This job's Julia depot is restored
        # from a CI cache, so its registry snapshot is as old as that cache,
        # and `Pkg.instantiate` (unlike `Pkg.add`) never updates registries on
        # its own — it resolves against whatever the depot happens to hold.
        # Without this the exact pin in `test/formatter/Project.toml` is
        # unresolvable from a depot cached before that version was registered.
        resolve = "using Pkg; Pkg.Registry.update(); Pkg.instantiate()"
        run(`julia --startup-file=no --project=$formatter_env -e $resolve`)
        cmd = Cmd(
            `julia --startup-file=no --project=$formatter_env $(joinpath(formatter_env, "runtests.jl"))`;
            ignorestatus = true
        )
        result = run(pipeline(cmd, stdout = stdout, stderr = stderr); wait = true)
        @test result.exitcode == 0
    else
        @test_skip "Formatter environment not found"
    end
end

@testitem "Runic pins agree" tags = [:quality] begin
    using Pkg: TOML

    # The pre-commit hook builds its own environment and never reads
    # `test/formatter/`, so the Runic version is declared twice. A bump that
    # moves only one leaves the hook writing a style the check above
    # rejects, on a tree nobody changed. Compare them here so the failure
    # names both versions instead.
    root = joinpath(@__DIR__, "..", "..")
    env_toml = TOML.parsefile(joinpath(root, "test", "formatter", "Project.toml"))
    env_pin = get(get(env_toml, "compat", Dict()), "Runic", nothing)

    config = read(joinpath(root, ".pre-commit-config.yaml"), String)
    hook_match = match(r"'Runic@([0-9]+\.[0-9]+\.[0-9]+)'", config)

    @test env_pin !== nothing
    @test hook_match !== nothing
    if env_pin !== nothing && hook_match !== nothing
        hook_pin = "=" * hook_match[1]
        if env_pin != hook_pin
            @warn "Runic pins have drifted; set both to one version and " *
                "reformat the tree with `task format`." env_pin hook_pin
        end
        @test env_pin == hook_pin
    end
end
