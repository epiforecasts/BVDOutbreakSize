## Tests for the relative-skill and per-release R0 additions in
## scripts/score_releases.jl: rel_skill_columns (log-CRPS skill to the
## persistence baseline and to a stream's individual fit) and r0_row
## (the established initial reproduction number, exp(rt_state.log_R0)).

@testitem "rel_skill_columns scores to baseline and individual fit" begin
    using Dates: Date
    using DataFrames: DataFrame

    ## Load the scorer's functions without its driver (guarded on
    ## PROGRAM_FILE), so rel_skill_columns can be exercised directly.
    include(joinpath(@__DIR__, "..", "scripts", "score_releases.jl"))

    ## One "confirmed cases" group with baseline, joint and an individual
    ## fit, plus a "recovered" group with joint and baseline only (no
    ## individual model). Rows are deliberately out of fit order.
    df = DataFrame(
        release = fill("results-vT", 5),
        made_date = fill(Date(2026, 7, 1), 5),
        stream = ["confirmed cases", "confirmed cases", "confirmed cases",
            "recovered", "recovered"],
        horizon = fill(7, 5),
        fit = ["baseline", "joint", "confirmed", "joint", "baseline"],
        log_crps = [0.8, 0.2, 0.4, 0.3, 0.6])

    to_base, to_indiv = rel_skill_columns(df)

    ## Every row of a group whose baseline was scored gets a baseline skill;
    ## the baseline row itself is 1.0.
    @test to_base == [1.0, 0.25, 0.5, 0.5, 1.0]  # 0.2/0.8, 0.4/0.8, 0.3/0.6

    ## The individual skill is defined only on the joint row of the group
    ## that carries an individual fit: confirmed cases' joint, 0.2/0.4.
    @test to_indiv[2] == 0.5
    ## Missing on the individual row, the baseline rows, and the recovered
    ## joint (no individual model fits it).
    @test all(ismissing, to_indiv[[1, 3, 4, 5]])
end

@testitem "rel_skill_columns leaves a frozen individual skill empty" begin
    using Dates: Date
    using DataFrames: DataFrame

    include(joinpath(@__DIR__, "..", "scripts", "score_releases.jl"))

    ## The frozen archive carries only the frozen fit and its baseline, no
    ## individual fits: baseline skill is defined, individual skill is not.
    df = DataFrame(
        release = fill("results-vT", 2),
        made_date = fill(Date(2026, 6, 20), 2),
        stream = ["confirmed cases", "confirmed cases"],
        horizon = fill(7, 2),
        fit = ["frozen", "baseline"],
        log_crps = [0.3, 0.6])

    to_base, to_indiv = rel_skill_columns(df)
    @test to_base == [0.5, 1.0]  # 0.3/0.6, then the baseline itself
    @test all(ismissing, to_indiv)
end

@testitem "rel_skill_cell empties missing, rounds a ratio" begin
    include(joinpath(@__DIR__, "..", "scripts", "score_releases.jl"))
    @test rel_skill_cell(missing) == ""
    @test rel_skill_cell(0.25) == "0.25"
    @test rel_skill_cell(1 / 3) == "0.3333"  # rounded to four places
end

@testitem "r0_row summarises exp(log_R0), skips when the column is absent" begin
    using Dates: Date

    include(joinpath(@__DIR__, "..", "scripts", "score_releases.jl"))

    dir = mktempdir()
    ## A draws asset carrying the raw walk base alongside R_T. exp of the
    ## three log values is [1, 2, 3], whose median is 2.
    path = joinpath(dir, "posterior_draws.csv")
    open(path, "w") do io
        println(io, "R_T,rt_state.log_R0")
        for v in (0.0, log(2.0), log(3.0))
            println(io, join((1.5, v), ','))
        end
    end
    row = r0_row("results-vT", path, Date(2026, 7, 1))
    @test row[1] == "results-vT"
    @test row[2] == "2026-07-01"
    @test row[3] ≈ 2.0  # median of exp.([0, log2, log3])

    ## A draws asset with no walk base is skipped, as is a missing asset.
    no_r0 = joinpath(dir, "no_r0.csv")
    open(no_r0, "w") do io
        println(io, "R_T")
        println(io, "1.5")
    end
    @test r0_row("results-vT", no_r0, Date(2026, 7, 1)) === nothing
    @test r0_row("results-vT", nothing, Date(2026, 7, 1)) === nothing
end
