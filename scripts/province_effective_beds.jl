# The beds a province holds patients in on a recorded day, for the
# `[province_bed_capacity_history]` block that scripts/province_care_manifest.jl
# emits. Dependency-free so the test suite can include it.

"""
    effective_beds(printed, patients, rate_pct) -> Int

Effective beds of a province on a recorded day: the largest of the printed
beds, the beds implied by the printed occupancy rate (`patients / rate_pct`
x 100, rounded, the rule the national `bed_capacity_history` uses) and the
patients held. A province holding more patients than its printed beds has
at least that many beds in practice. `patients` and `rate_pct` may be
`nothing` when the report does not print them.
"""
function effective_beds(printed::Integer, patients, rate_pct)
    beds = Int(printed)
    patients === nothing && return beds
    if rate_pct !== nothing && rate_pct > 0
        beds = max(beds, round(Int, 100 * patients / rate_pct))
    end
    return max(beds, Int(patients))
end
