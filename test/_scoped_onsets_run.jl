using TestItemRunner

target = joinpath(pwd(), "test", "test_onsets.jl")
@run_package_tests filter = ti -> string(ti.filename) == target
