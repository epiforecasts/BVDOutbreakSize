# The beds a province holds patients in on a recorded day, for the
# `[province_bed_capacity_history]` block that scripts/province_care_manifest.jl
# emits. Dependency-free so the test suite can include it.

"""
    effective_beds(printed, patients, rate_pct) -> Int

Effective beds of a province on a recorded day: the larger of the printed
beds and the patients held. When the patients exceed the printed beds, the
bed column is stale or the province is overcrowded, and the beds implied by
the printed occupancy rate (`patients / rate_pct` x 100, rounded, the rule
the national `bed_capacity_history` uses) count too. Below the printed beds
the rate only restates them, up to its rounding. `patients` and `rate_pct`
may be `nothing` when the report does not print them.
"""
function effective_beds(printed::Integer, patients, rate_pct)
    beds = Int(printed)
    (patients === nothing || patients <= beds) && return beds
    if rate_pct !== nothing && rate_pct > 0
        beds = max(beds, round(Int, 100 * patients / rate_pct))
    end
    return max(beds, Int(patients))
end
