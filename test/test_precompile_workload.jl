@testitem "precompile workload is the method instance the fit builds" begin
    using BVDOutbreakSize: precompile_workload_joint, production_joint,
        default_breakpoint, load_observations

    obs = load_observations()
    bp = default_breakpoint(obs)

    ## A `DynamicPPL.Model`'s type carries the types of its arguments, and
    ## that type is the key Mooncake caches a reverse rule against. Equal
    ## types mean the workload compiles the rule the fit reuses; unequal
    ## types mean the workload compiles a rule nothing reaches and every fit
    ## pays the cold build. This is the assertion that stops the two drifting
    ## apart: an argument added to one call alone fails here.
    @test typeof(precompile_workload_joint()) ==
        typeof(production_joint(obs; breakpoint = bp))
end

@testitem "joint fit arguments carry the streams the fit scores" begin
    using BVDOutbreakSize: joint_fit_args, patch_fit_args, default_breakpoint,
        load_observations

    obs = load_observations()
    args = joint_fit_args(obs; breakpoint = default_breakpoint(obs))

    ## The workload used to leave the isolation and treatment histories
    ## empty. Mooncake derives a rule only for code that executes, so a
    ## zero-iteration likelihood loop caches nothing, and `treatment_flow`
    ## is the most expensive component in the model. Non-empty here is what
    ## makes the precompile workload reach it.
    for stream in (
            :isolation_history, :bed_capacity_history,
            :treatment_admissions_history, :treatment_deaths_history,
            :confirmed_history, :reported_history, :deaths_history,
            :lab_history, :onset_curve_history,
        )
        @test haskey(args, stream)
    end
    @test !isempty(args.isolation_history.counts)
    @test !isempty(args.confirmed_history.counts)

    ## The province tables reach the model as integer matrices. Passing
    ## `missing` here, as an earlier workload did, is a different method
    ## instance from the one the headline fit builds.
    patch = patch_fit_args(obs)
    @test patch.province_increments isa AbstractMatrix{<:Integer}
    @test patch.province_death_increments isa AbstractMatrix{<:Integer}
    @test patch.n_patches > 1
end
