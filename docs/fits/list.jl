# Print the fit ids as a JSON array. Set `BVD_RUN_SENSITIVITY=true` to include
# the sensitivity re-fits.
#
#   julia --project=docs docs/fits/list.jl   # ["joint","exports",...]
#
# With `--groups` it prints one `<group>_ids=[...]` line per fit group
# instead, for the `fit_<group>` jobs in `.github/workflows/docs.yml`. Each
# render job waits only on the groups its page reads, so the groups are the
# unit a page can depend on:
#
# - `joint`: the headline joint and its no-patches control
# - `streams`: the single-stream fits
# - `frozen`: every fit to data frozen at an earlier cut-off
# - `sensitivity`: the gated sensitivity re-fits (empty unless enabled)
using Pkg: Pkg
Pkg.instantiate()

using BVDOutbreakSize
include(joinpath(@__DIR__, "registry.jl"))

json_array(ids) = "[" * join(("\"" * id * "\"" for id in ids), ",") * "]"

function fit_group(spec)
    spec.id in ("joint", "sens_no_patches") && return "joint"
    spec.kind === :frozen && return "frozen"
    startswith(spec.id, "sens_") && return "sensitivity"
    return "streams"
end

if "--groups" in ARGS
    specs = build_fit_specs(
        load_observations(); run_sensitivity = run_sensitivity_env()
    )
    for group in ("joint", "streams", "frozen", "sensitivity")
        ids = [s.id for s in specs if fit_group(s) == group]
        println(group, "_ids=", json_array(ids))
    end
else
    println(json_array(fit_ids()))
end
