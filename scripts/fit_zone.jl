# Fit the health-zone model from a parent joint fit, outside the docs build.
#
#   julia --project=docs scripts/fit_zone.jl \
#       --parent logs/fit_cache/joint__<hash>.parent.jls \
#       --out logs/zone_local [--samples 800] [--warmup 500] [--chains 2] \
#       [--depth 10] [--target-accept 0.8]
#
# `--parent` is a parent extract (`<key>.parent.jls`, the `fit-extras-joint`
# artefact of a CI run) or a full `joint__<hash>.jls` chain. The observations
# are frozen to the parent's cut-off. The chain is written to
# `<out>/zone_chain.jls` with the diagnostics bundle under
# `<out>/diagnostics/` and the fit summary printed and appended to
# `<out>/summary.md`.
using Pkg: Pkg
Pkg.instantiate()

using BVDOutbreakSize
using Serialization: serialize, deserialize
using Dates: Day
include(joinpath(@__DIR__, "..", "docs", "fits", "cache.jl"))
include(joinpath(@__DIR__, "..", "docs", "fits", "summary.jl"))

function parse_args(args)
    opts = Dict{String, String}(
        "out" => "logs/zone_local", "samples" => "800", "warmup" => "500",
        "chains" => "2", "depth" => "10", "target-accept" => "0.8"
    )
    i = 1
    while i <= length(args)
        a = args[i]
        startswith(a, "--") || error("unexpected argument $(a)")
        i + 1 <= length(args) || error("missing value for $(a)")
        opts[a[3:end]] = args[i + 1]
        i += 2
    end
    haskey(opts, "parent") || error("--parent <extract or chain .jls> is required")
    return opts
end

const OPTS = parse_args(ARGS)
const OUT = abspath(OPTS["out"])
mkpath(OUT)

parent = repair_chain_keys(deserialize(OPTS["parent"]))
n_parent = if parent isa NamedTuple && haskey(parent, :draws)
    parent.n
else
    length(first(vec(collect(parent[:infections_patch])))) ÷
        length(PROVINCE_NAMES)
end
obs = let o = load_observations()
    n_parent == o.n ? o : freeze_observations(o.seeding + Day(n_parent - 1))
end
@info "fit_zone.jl" parent = OPTS["parent"] cutoff = obs.cutoff n = obs.n out = OUT

chn = fit_zone(
    parent, obs;
    samples = parse(Int, OPTS["samples"]), chains = parse(Int, OPTS["chains"]),
    n_adapts = parse(Int, OPTS["warmup"]), max_depth = parse(Int, OPTS["depth"]),
    target_accept = parse(Float64, OPTS["target-accept"])
)
serialize(joinpath(OUT, "zone_chain.jls"), chn)
write_fit_bundle("local", chn, joinpath(OUT, "diagnostics"))
md = fit_summary_markdown("local", chn)
print(stdout, md)
open(io -> print(io, md), joinpath(OUT, "summary.md"), "a")
inputs = zone_fit_inputs(parent, obs)
d = zone_sampler_diagnostics(chn, inputs; max_depth = parse(Int, OPTS["depth"]))
@info "zone sampler" d.max_rhat d.min_ess_bulk d.n_divergent d.depth_cap_fraction d.step_size d.ebfmi
