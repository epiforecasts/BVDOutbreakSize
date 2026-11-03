# The beds a province holds patients in on a recorded day, for the
# `[province_bed_capacity_history]` block that scripts/province_care_manifest.jl
# emits. Dependency-free so the test suite can include it.

"""
    effective_beds(printed, patients, rate_pct) -> Int

Effective beds of a province on a recorded day: the printed beds, or the
beds implied by the printed occupancy rate (`patients / rate_pct` x 100,
rounded, the rule the national `bed_capacity_history` uses) when those are
more and the patients exceed the printed beds. The printed bed column can go
stale while the rate is computed against the current beds. Below the printed
beds the rate only restates them, up to its rounding. The patients held are
not a bed count, so a province can hold more patients than beds. `patients`
and `rate_pct` may be `nothing` when the report does not print them.
"""
function effective_beds(printed::Integer, patients, rate_pct)
    beds = Int(printed)
    (patients === nothing || patients <= beds) && return beds
    (rate_pct === nothing || rate_pct <= 0) && return beds
    return max(beds, round(Int, 100 * patients / rate_pct))
end
