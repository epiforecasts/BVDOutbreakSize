@testitem "Enzyme AD extension (isolated env)" tags = [:quality] begin
    using Pkg
    enzyme_env = joinpath(@__DIR__, "..", "enzyme")
    ## Enzyme reverse-mode is not viable on Windows for this model (the
    ## extension fails to precompile / segfaults), and may not have a
    ## compatible release on experimental Julia, so the Enzyme checks run
    ## only on Linux and macOS off the experimental matrix entry. They are
    ## isolated here so Enzyme never enters the main test environment.
    runnable = !Sys.iswindows() &&
        get(ENV, "JULIA_CI_EXPERIMENTAL", "false") != "true" &&
        isdir(enzyme_env) &&
        isfile(joinpath(enzyme_env, "Project.toml"))
    if runnable
        run(
            pipeline(
                `julia --startup-file=no --project=$enzyme_env -e "using Pkg; Pkg.instantiate()"`,
                stdout = stdout, stderr = stderr
            )
        )
        ## One wall-clock bound per part of the script. Enzyme can run its
        ## type analysis for hours on a model rather than fail on it, which
        ## on Julia 1.13.1 it does on the full joint. A part past its bound
        ## is killed and counts as the script failing to run, and the other
        ## part still reports. The joint has its own run so that its
        ## compile cannot take the per-component results with it.
        ##
        ## Measured on a 4-core machine: the components take about 18 min
        ## on Julia 1.13.1, and the joint throws after 15 min on 1.13.0.
        ## Each bound is several times that, and the two together fit the
        ## quality job's timeout on a cold depot.
        components_min = 60
        joint_min = 40
        function run_part(part, limit_min)
            proc = run(
                pipeline(
                    Cmd(`julia --startup-file=no --project=$enzyme_env $(joinpath(enzyme_env,
                    "runtests.jl")) $part`; ignorestatus = true),
                    stdout = stdout, stderr = stderr
                );
                wait = false
            )
            if timedwait(
                    () -> process_exited(proc), 60 * limit_min;
                    pollint = 10
                ) === :timed_out
                @warn "Enzyme $part run killed after $limit_min min"
                kill(proc, Base.SIGKILL)
                wait(proc)
                return nothing
            end
            return proc.exitcode
        end
        ## Enzyme's platform/version instability is tolerated: a crash or a
        ## killed run (an upstream Enzyme/LLVM failure, not a model issue)
        ## is recorded as broken rather than failing the suite. Mooncake,
        ## the default backend, is asserted to differentiate every model in
        ## the main suite.
        ##
        ## Exit code 2 is the script's own marker for a scenario that did
        ## not behave as `ADFixtures` declares, which is a regression
        ## rather than an upstream wobble and fails here. Tolerating it
        ## identically would leave the sweep unable to report anything it
        ## had not already been told to expect. Every other non-zero code
        ## is the script failing to run, which stays tolerated.
        for (part, limit_min) in (("components", components_min), ("joint", joint_min))
            code = run_part(part, limit_min)
            if code == 0
                @test true
            elseif code == 2
                @test code == 0
            else
                @test_broken code == 0
            end
        end
    else
        @test_skip "Enzyme environment not run (Windows / experimental / missing)"
    end
end
