@testitem "Code formatting" tags=[:quality] begin
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
        # its own — it resolves against whatever the depot happens to hold. The
        # pre-commit hook runs the same formatter over the same directories but
        # reaches it through `Pkg.add`, so without this the hook and this test
        # can run different JuliaFormatter versions and disagree about an
        # identical tree. It also makes the exact pin in
        # `test/formatter/Project.toml` resolvable from a depot cached before
        # that version was registered.
        resolve = "using Pkg; Pkg.Registry.update(); Pkg.instantiate()"
        run(`julia --project=$formatter_env -e $resolve`)
        cmd = Cmd(
            `julia --project=$formatter_env $(joinpath(formatter_env, "runtests.jl"))`;
            ignorestatus = true)
        result = run(pipeline(cmd, stdout = stdout, stderr = stderr); wait = true)
        @test result.exitcode == 0
    else
        @test_skip "Formatter environment not found"
    end
end
