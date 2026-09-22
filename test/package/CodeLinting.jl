@testitem "JET type-stability checks" tags = [:quality] begin
    using Pkg
    jet_env = joinpath(@__DIR__, "..", "jet")
    # Skip on experimental Julia (pre) where JET may not have a compatible
    # release yet. No CI job sets JULIA_CI_EXPERIMENTAL, so the guard is
    # inert there and this is a local escape hatch.
    if VERSION >= v"1.10" && get(ENV, "JULIA_CI_EXPERIMENTAL", "false") != "true" &&
            isdir(jet_env) && isfile(joinpath(jet_env, "Project.toml"))
        run(
            pipeline(
                `julia --startup-file=no --project=$jet_env -e "using Pkg; Pkg.instantiate()"`,
                stdout = stdout, stderr = stderr
            )
        )
        result = run(
            pipeline(
                Cmd(
                    `julia --startup-file=no --project=$jet_env $(joinpath(jet_env, "runtests.jl"))`;
                    ignorestatus = true
                ),
                stdout = stdout, stderr = stderr
            )
        )
        @test result.exitcode == 0
    else
        @test_skip "JET environment not run (pre Julia or missing env)"
    end
end
