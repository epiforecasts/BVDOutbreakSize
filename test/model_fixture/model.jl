## The model and the helpers shared by `generate.jl`, which writes
## `joint.toml`, and `test/test_model_fixture.jl`, which checks the package
## against it. Its imports come from whichever file includes it.

## A fixed past cut-off, so the daily data updates leave the model's data
## unchanged.
const FIXTURE_CUTOFF = "2026-09-10"
const FIXTURE_PATH = joinpath(@__DIR__, "joint.toml")

## The headline joint, as `production_joint` builds it, on the data as of
## `FIXTURE_CUTOFF`.
function fixture_joint()
    obs = freeze_observations(FIXTURE_CUTOFF)
    return production_joint(obs; breakpoint = default_breakpoint(obs))
end

## The log density NUTS samples, with the package's gradient backend.
function fixture_ldf(model, vi)
    return LogDensityFunction(
        model, getlogjoint_internal, vi; adtype = default_adtype()
    )
end

## Each sampled varname's slice of the unconstrained vector, by name.
function varname_ranges(ldf)
    return Dict(
        string(vn) => rt.range
            for (vn, rt) in pairs(get_all_ranges_and_transforms(ldf))
    )
end

## The recorded `:=` quantities at `x`: a scalar as itself and a vector as
## its sum. Sampled varnames are left out.
function reported_values(ldf, x, sampled)
    pws = ParamsWithStats(x, ldf)
    out = Dict{String, Float64}()
    for (vn, v) in pairs(pws.params)
        k = string(vn)
        k in sampled && continue
        if v isa Real
            out[k] = Float64(v)
        elseif v isa AbstractArray{<:Real}
            out[k] = Float64(sum(v))
        end
    end
    return out, pws.stats
end
