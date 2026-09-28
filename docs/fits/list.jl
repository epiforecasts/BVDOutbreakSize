# Print the fit ids as a JSON array. Set `BVD_RUN_SENSITIVITY=true` to include
# the sensitivity re-fits.
#
#   julia --project=docs docs/fits/list.jl   # ["joint","exports",...]
#
# `BVD_FIT_STAGE=base` or `dependent` narrows that listing to one stage.
#
# With `--groups` it prints one `<group>_ids=[...]` line per fit group
# instead, for the `fit_<group>` jobs in `.github/workflows/docs.yml`. Each
# render job waits only on the groups its page reads, so the groups are the
# unit a page can depend on. `fit_group` in `registry.jl` assigns them:
#
# - `joint`: the headline joint and its no-patches control
# - `streams`: the single-stream fits
# - `frozen`: every fit to data frozen at an earlier cut-off
# - `sensitivity`: the gated sensitivity re-fits (empty unless enabled)
# - `zone`: the health-zone fits melded from a cached parent chain
#
# The base groups print a JSON array of ids, which the job spreads over
# `matrix.id`. The `zone` group prints `matrix.include` entries instead,
# because each of its fits also needs the artifact pattern holding its parent
# chains, which the job downloads before fitting.
using Pkg: Pkg
Pkg.instantiate()

using BVDOutbreakSize
include(joinpath(@__DIR__, "registry.jl"))

json_string(s) = "\"" * s * "\""
json_array(items) = "[" * join(items, ",") * "]"
function json_object(pairs)
    return "{" *
        join((json_string(k) * ":" * json_string(v) for (k, v) in pairs), ",") *
        "}"
end

## The fit jobs upload their chains as `fit-<id>`, so a single parent is a
## plain artifact name and a longer list a brace glob `download-artifact`
## expands.
function parent_pattern(spec)
    needs = fit_needs(spec)
    length(needs) == 1 && return "fit-" * only(needs)
    return "fit-{" * join(needs, ",") * "}"
end

function matrix_entry(spec)
    return json_object(("id" => spec.id, "parents" => parent_pattern(spec)))
end

if "--groups" in ARGS
    specs = build_fit_specs(
        load_observations(); run_sensitivity = run_sensitivity_env()
    )
    for group in FIT_GROUPS
        in_group = [s for s in specs if fit_group(s) == group]
        items = group == "zone" ? map(matrix_entry, in_group) :
            [json_string(s.id) for s in in_group]
        println(group, "_ids=", json_array(items))
    end
else
    ## `BVD_FIT_STAGE` narrows the plain listing to one stage, the same
    ## variable `all.jl` runs a single pass with.
    ids = fit_ids(; stage = fit_stage_env(:all))
    println(json_array([json_string(id) for id in ids]))
end
