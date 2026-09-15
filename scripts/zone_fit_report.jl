# Standalone fit report for the health-zone model, following the Bayesian
# workflow stage by stage: prior predictive check, simulation-based
# recovery, the real-data fit's sampler diagnostics, posterior predictive
# checks and the reporting table and map. Each stage writes its figures and
# an HTML fragment under the output directory, and `index.html` is
# reassembled from whatever fragments exist, so the report can be hosted
# while the long stages are still running.
#
#   julia --project=docs scripts/zone_fit_report.jl \
#       --parent logs/fit_cache/joint__<hash>.jls \
#       --chain logs/zone_report/zone_chain.jls \
#       --out logs/zone_report --stages a,c
#
# `--parent` is the joint chain the zone model melds from (required).
# `--chain` is a fitted zone chain for stages c, d and e; when absent those
# stages fit one at `--samples`/`--warmup` and save it there. `--stages`
# selects which of a, b, c, d, e to run (default all). The observations are
# frozen to the parent's cut-off. Stage a saves its prior draws next to the
# chain for stage d to reuse. `--recovery-samples`/`--recovery-warmup` set
# the recovery refit (default 200/200) and `--top` the zones shown per
# patch (default 10).

using BVDOutbreakSize
using BVDOutbreakSize: bvd_zone, zone_forward, zone_initial_shares,
                       zone_deviation_knots, zone_initial_params,
                       _zone_draws, _zone_stat, _draws, _has_key
using Serialization: serialize, deserialize
using Statistics: mean, median, quantile
using Random: Xoshiro
using Dates: Day, now
using DataFrames: DataFrame, names, eachrow
using Distributions: DirichletMultinomial
using Turing: sample, Prior
const FlexiChains = BVDOutbreakSize.FlexiChains
import CairoMakie
using CairoMakie: Figure, Axis, lines!, band!, scatter!, vlines!, save
using Base64: base64encode
include(joinpath(@__DIR__, "..", "docs", "fits", "cache.jl"))

## --- Arguments ------------------------------------------------------------

function parse_args(args)
    opts = Dict{String, String}("out" => "logs/zone_report",
        "stages" => "a,b,c,d,e", "samples" => "600", "warmup" => "400",
        "recovery-samples" => "200", "recovery-warmup" => "200",
        "prior-draws" => "200", "top" => "10")
    i = 1
    while i <= length(args)
        a = args[i]
        startswith(a, "--") || error("unexpected argument $(a)")
        i + 1 <= length(args) || error("missing value for $(a)")
        opts[a[3:end]] = args[i + 1]
        i += 2
    end
    haskey(opts, "parent") || error("--parent <joint chain .jls> is required")
    return opts
end

const OPTS = parse_args(ARGS)
const OUT = abspath(OPTS["out"])
mkpath(OUT)
const STAGES = Set(strip.(split(OPTS["stages"], ",")))
const TOP = parse(Int, OPTS["top"])
const NCOLS = 5
const CHAIN_PATH = get(OPTS, "chain", joinpath(OUT, "zone_chain.jls"))
const PRIOR_PATH = joinpath(OUT, "a_prior_chain.jls")

## --- Shared inputs --------------------------------------------------------

log_line(msg) = (println("[", now(), "] ", msg); flush(stdout))

log_line("loading the parent chain $(OPTS["parent"])")
const PARENT = repair_chain_keys(deserialize(OPTS["parent"]))
const OBS = let
    v = first(vec(collect(PARENT[:infections_patch])))
    n_parent = length(v) ÷ length(PROVINCE_NAMES)
    o = load_observations()
    n_parent == o.n ? o : freeze_observations(o.seeding + Day(n_parent - 1))
end
log_line("observations frozen to $(OBS.cutoff) (n = $(OBS.n))")
const INPUTS = zone_fit_inputs(PARENT, OBS)
const ZD = INPUTS.model_data
const NZ = length(INPUTS.zone_keys)
const K = length(INPUTS.knots)
log_line("$(NZ) zones, $(count(INPUTS.walking)) walking, $(K) knots, " *
         "$(length(ZD.cell_patch)) cells")

## --- HTML helpers ---------------------------------------------------------

function png_tag(path; width = "100%")
    data = base64encode(read(path))
    return "<img src=\"data:image/png;base64,$(data)\" " *
           "style=\"width:$(width);max-width:1500px\">"
end

function html_table(df::DataFrame; digits = 3)
    fmt(x) = x isa AbstractFloat ? (isnan(x) ? "" : string(round(x; digits))) :
             string(x)
    io = IOBuffer()
    print(io, "<table><thead><tr>")
    for c in names(df)
        print(io, "<th>", c, "</th>")
    end
    print(io, "</tr></thead><tbody>")
    for r in eachrow(df)
        print(io, "<tr>")
        for c in names(df)
            print(io, "<td>", fmt(r[c]), "</td>")
        end
        print(io, "</tr>")
    end
    print(io, "</tbody></table>")
    return String(take!(io))
end

function write_fragment(stage, title, body)
    open(joinpath(OUT, "fragment_$(stage).html"), "w") do io
        println(io, "<section id=\"stage-$(stage)\"><h2>$(title)</h2>")
        println(io, body)
        println(io, "<p class=\"stamp\">Written $(now()).</p></section>")
    end
    assemble_index()
end

const TITLES = (("a", "Prior predictive check"),
    ("b", "Simulation-based recovery"),
    ("c", "Real-data fit: sampler diagnostics"),
    ("d", "Posterior predictive checks"),
    ("e", "Zone ranking and map"))

function assemble_index()
    css = """
    body { font-family: system-ui, sans-serif; max-width: 1540px; margin: 2rem auto; padding: 0 1rem; color: #222; }
    table { border-collapse: collapse; font-size: 0.85rem; margin: 1rem 0; }
    th, td { border: 1px solid #ccc; padding: 0.2rem 0.5rem; text-align: right; }
    th { background: #f0f0f0; }
    td:first-child, th:first-child { text-align: left; }
    .stamp { color: #777; font-size: 0.8rem; }
    section { margin-bottom: 3rem; }
    """
    open(joinpath(OUT, "index.html"), "w") do io
        println(io, "<!doctype html><html><head><meta charset=\"utf-8\">")
        println(io, "<title>Health-zone model fit report</title>")
        println(io, "<style>$(css)</style></head><body>")
        println(io, "<h1>Health-zone model fit report</h1>")
        println(io,
            "<p>Parent chain <code>$(OPTS["parent"])</code>; " *
            "observations frozen to $(OBS.cutoff); $(NZ) zones, " *
            "$(count(INPUTS.walking)) walking, $(K) weekly knots, " *
            "$(length(ZD.cell_patch)) scored patch-vintage cells.</p>")
        for (stage, title) in TITLES
            path = joinpath(OUT, "fragment_$(stage).html")
            if isfile(path)
                print(io, read(path, String))
            else
                println(io,
                    "<section><h2>$(title)</h2><p>Not yet run.</p></section>")
            end
        end
        println(io, "</body></html>")
    end
    log_line("index.html assembled at $(joinpath(OUT, "index.html"))")
end

## --- Figure helpers -------------------------------------------------------

## The composition figures, one per patch, saved and embedded in order.
function composition_tags(chn, stem; prior_chain = nothing, label)
    figs = plot_zone_composition_ppc(chn, INPUTS; prior_chain, top = TOP,
        ncols = NCOLS, label)
    tags = String[]
    for (i, fig) in enumerate(figs)
        path = joinpath(OUT, "$(stem)_$(i).png")
        save(path, fig)
        push!(tags, png_tag(path))
    end
    return join(tags, "\n")
end

function trace_panels(keys_and_labels; path)
    n = length(keys_and_labels)
    fig = Figure(; size = (900, 160 * n + 20))
    for (k, (m, label)) in enumerate(keys_and_labels)
        ax = Axis(fig[k, 1]; ylabel = label,
            xlabel = k == n ? "Iteration" : "")
        for c in 1:size(m, 2)
            lines!(ax, 1:size(m, 1), m[:, c]; linewidth = 0.7)
        end
    end
    save(path, fig)
    return path
end

## Iteration-by-chain matrix of a scalar chain key, or of one element of a
## vector deterministic.
chain_matrix(chn, key) = Float64.(Array(chn[key]))
chain_matrix(chn, key, idx) =
    let m = chn[key]
        reshape(Float64[v[idx] for v in m], size(m))
    end

## Divergences per chain.
function divergences_per_chain(chn, nchains)
    d = _zone_stat(chn, :numerical_error)
    d === nothing && return fill(0, nchains)
    return [Int(sum(view(d, :, c))) for c in 1:nchains]
end

## `min – max` of the finite entries of a column, over the rows `keep`.
function range_string(df, col, keep = trues(size(df, 1)))
    v = filter(isfinite, df[keep, col])
    isempty(v) && return ""
    return "$(round(minimum(v); digits = 3)) – $(round(maximum(v); digits = 3))"
end

## The R-hat and ESS ranges of the cut-off quantities over the zones, the
## reproduction number over the walking zones only.
function range_table(chn, inputs)
    zdiag = zone_diagnostics_table(chn, inputs)
    walking = collect(inputs.walking)
    rows = [("R_T_zone (walking)", "R_T", walking),
        ("share_T_zone", "share_T", trues(NZ)),
        ("delta_T_zone", "delta_T", trues(NZ))]
    return DataFrame(quantity = first.(rows),
        rhat = [range_string(zdiag, "rhat_$(s)", k) for (_, s, k) in rows],
        ess_bulk = [range_string(zdiag, "ess_bulk_$(s)", k)
                    for (_, s, k) in rows],
        ess_tail = [range_string(zdiag, "ess_tail_$(s)", k)
                    for (_, s, k) in rows]),
    zdiag
end

## --- Stage a: prior predictive ---------------------------------------------

function prior_chain()
    if isfile(PRIOR_PATH)
        log_line("loading the prior draws $(PRIOR_PATH)")
        return deserialize(PRIOR_PATH)
    end
    ndraw = parse(Int, OPTS["prior-draws"])
    log_line("sampling $(ndraw) prior draws")
    chn = sample(bvd_zone(ZD), Prior(), ndraw;
        chain_type = FlexiChains.VNChain, progress = false)
    serialize(PRIOR_PATH, chn)
    return chn
end

function stage_a()
    chn = prior_chain()
    ndraw = size(chn, 1)
    body = """
    <p>$(ndraw) draws from the zone model prior pushed through the share
    renewal and the reporting delay, binned to the vintages, for the
    $(TOP) largest zones per patch. Points are the observed share of each
    patch's allocated confirmed cases at each vintage.</p>
    $(composition_tags(chn, "a_prior_predictive"; label = "prior predictive"))
    """
    write_fragment("a", "Prior predictive check", body)
end

## --- Stage b: simulation-based recovery -----------------------------------

function stage_b()
    log_line("stage b: recovery")
    rng = Xoshiro(20260915)
    ## One draw from a moderate part of the prior, rather than its tails, so
    ## the simulated dataset resembles the real one in dispersion.
    z_w = INPUTS.z_w_start .+ 0.2 .* randn(rng, NZ)
    σ_level = 0.25
    σ_δ = 0.06
    h = 42.0
    ρ = 0.03
    z_level = randn(rng, NZ)
    z_drift = randn(rng, max(ZD.n_walking * (K - 1), 1))
    φ = exp2(-7 / h)
    w0 = zone_initial_shares(z_w, ZD.patch_ranges, 2.0)
    δ_knots = zone_deviation_knots(z_level, z_drift, σ_level, σ_δ, φ,
        ZD.patch_ranges, ZD.walking, ZD.walk_index, ZD.n_walking, K)
    fw = zone_forward(ZD, δ_knots, w0, nothing)
    κ = (1 - ρ) / ρ
    counts = zeros(Int, NZ, length(INPUTS.days))
    for c in eachindex(ZD.cell_patch)
        zs = ZD.patch_ranges[ZD.cell_patch[c]]
        v = ZD.cell_vintage[c]
        π = max.(fw.increments[zs, v], 1e-12)
        π ./= sum(π)
        counts[zs, v] = rand(rng,
            DirichletMultinomial(ZD.cell_total[c], κ .* π))
    end
    ## The simulated observations in the manifest shape.
    hist = Dict{String, Dict{String, NamedTuple}}()
    for z in 1:NZ
        prov = INPUTS.zone_province[z]
        haskey(hist, prov) || (hist[prov] = Dict{String, NamedTuple}())
        hist[prov][INPUTS.zone_names[z]] = (; days = copy(INPUTS.days),
            counts = cumsum(counts[z, :]))
    end
    sim_obs = merge(OBS,
        (; zone_confirmed_history = hist,
            zone_death_history = Dict{String, Dict{String, NamedTuple}}()))
    ## The walking set is data-dependent, so the refit's own inputs say
    ## which zones carry a walk.
    sim_inputs = zone_fit_inputs(PARENT, sim_obs)
    nd = ZD.n - ZD.t0 + 1
    truth = (; share_T = fw.shares[nd, :], delta_T = δ_knots[:, K],
        R_T = [fw.infections[nd, z] / max(fw.forces[nd, z], floatmin())
               for z in 1:NZ],
        σ_δ, ρ, h, walking = sim_inputs.walking)
    serialize(joinpath(OUT, "b_truth.jls"), truth)
    t = time()
    chn = fit_zone(PARENT, sim_obs;
        samples = parse(Int, OPTS["recovery-samples"]),
        n_adapts = parse(Int, OPTS["recovery-warmup"]),
        callback = progress_callback(; path = joinpath(OUT, "b_fit.log")))
    minutes = round((time() - t) / 60; digits = 1)
    log_line("recovery fit took $(minutes) min")
    serialize(joinpath(OUT, "b_chain.jls"), chn)
    fig = Figure(; size = (1200, 800))
    cover = Float64[]
    shown = Int[]
    for (k, (key, tkey, label)) in enumerate((
        (:share_T_zone, :share_T, "share_T"),
        (:delta_T_zone, :delta_T, "δ_T"), (:R_T_zone, :R_T, "R_T")))
        draws = _zone_draws(chn, key, NZ)
        tr = getfield(truth, tkey)
        ax = Axis(fig[1, k]; title = label, xlabel = "true",
            ylabel = "posterior")
        lo = Float64[]
        hi = Float64[]
        md = Float64[]
        xs = Float64[]
        for z in 1:NZ
            d = filter(isfinite, draws[z])
            (isempty(d) || !isfinite(tr[z])) && continue
            ## A level-only zone's R_T is its patch's value times a
            ## prior-driven level, so only the walking zones are scored.
            key === :R_T_zone && !sim_inputs.walking[z] && continue
            push!(xs, tr[z])
            push!(lo, quantile(d, 0.05))
            push!(hi, quantile(d, 0.95))
            push!(md, median(d))
        end
        CairoMakie.rangebars!(ax, xs, lo, hi; color = (:steelblue, 0.6))
        scatter!(ax, xs, md; color = :steelblue)
        lines!(ax, [minimum(xs), maximum(xs)], [minimum(xs), maximum(xs)];
            color = :black, linestyle = :dash)
        push!(cover, mean(lo .<= xs .<= hi))
        push!(shown, length(xs))
    end
    scalars = ((:region_drift_sd_zone, :σ_δ, "σ_δ"),
        (:composition_rho_zone, :ρ, "ρ"), (:region_halflife_zone, :h, "h"))
    rows = NamedTuple[]
    for (k, (key, tkey, label)) in enumerate(scalars)
        d = _draws(chn, key)
        tr = getfield(truth, tkey)
        ax = Axis(fig[2, k]; title = label, xlabel = label)
        CairoMakie.hist!(ax, d; bins = 40, color = (:steelblue, 0.6))
        vlines!(ax, [tr]; color = :black, linestyle = :dash)
        push!(rows,
            (quantity = label, truth = tr, median = median(d),
                lower_90 = quantile(d, 0.05), upper_90 = quantile(d, 0.95),
                covered = quantile(d, 0.05) <= tr <= quantile(d, 0.95)))
    end
    path = joinpath(OUT, "b_recovery.png")
    save(path, fig)
    cov_df = DataFrame(
        quantity = ["share_T_zone", "delta_T_zone", "R_T_zone (walking)"],
        zones = shown, coverage_90 = cover)
    sd = zone_sampler_diagnostics(chn, sim_inputs; max_depth = 8)
    body = """
    <p>One dataset simulated from a known parameter draw on the real design
    (the real patch infections, vintages and allocated totals; counts drawn
    from the Dirichlet-multinomial at ρ = $(ρ)), refitted at
    $(OPTS["recovery-samples"])/$(OPTS["recovery-warmup"]) in $(minutes)
    min. Top row: true against posterior median with 90% intervals per
    zone. Bottom row: posterior of the scalars with the truth dashed.</p>
    $(png_tag(path))
    <h3>90% coverage over zones</h3>$(html_table(cov_df))
    <h3>Scalars</h3>$(html_table(DataFrame(rows); digits = 3))
    <p>Refit diagnostics: max R-hat $(round(sd.max_rhat; digits = 3))
    (walking zones' R_T $(round(sd.max_rhat_R_T_walking; digits = 3))),
    min bulk ESS $(round(sd.min_ess_bulk; digits = 0)), divergences
    $(sd.n_divergent), depth-cap fraction
    $(round.(sd.depth_cap_fraction; digits = 2)).</p>
    """
    write_fragment("b", "Simulation-based recovery", body)
end

## --- Stage c: real-data fit diagnostics -----------------------------------

function load_or_fit_chain()
    if isfile(CHAIN_PATH)
        log_line("loading the zone chain $(CHAIN_PATH)")
        return repair_chain_keys(deserialize(CHAIN_PATH))
    end
    log_line("fitting the zone model at $(OPTS["samples"])/$(OPTS["warmup"])")
    t = time()
    chn = fit_zone(PARENT, OBS; samples = parse(Int, OPTS["samples"]),
        n_adapts = parse(Int, OPTS["warmup"]),
        callback = progress_callback(; path = joinpath(OUT, "zone_fit.log")))
    minutes = round((time() - t) / 60; digits = 1)
    log_line("fit took $(minutes) min")
    serialize(CHAIN_PATH, chn)
    write(CHAIN_PATH * ".minutes", string(minutes))
    return chn
end

function stage_c(chn)
    log_line("stage c: diagnostics")
    sd = zone_sampler_diagnostics(chn, INPUTS; max_depth = 8)
    nchains = length(sd.step_size)
    init = zone_initial_params(bvd_zone(ZD), INPUTS; chains = nchains)
    per_chain = DataFrame(chain = 1:nchains,
        initial_logp = init.logp,
        step_size = sd.step_size,
        depth_cap_fraction = sd.depth_cap_fraction,
        ebfmi = sd.ebfmi,
        divergences = divergences_per_chain(chn, nchains))
    minutes_path = CHAIN_PATH * ".minutes"
    minutes = isfile(minutes_path) ? read(minutes_path, String) : "unknown"
    overall = DataFrame(
        quantity = ["max R-hat (all but R_T_zone)",
            "max R-hat R_T_zone (walking)", "min bulk ESS", "min tail ESS",
            "divergences", "draws x chains", "wall time (min)"],
        value = [sd.max_rhat, sd.max_rhat_R_T_walking, sd.min_ess_bulk,
            sd.min_ess_tail, sd.n_divergent, prod(size(chn)), minutes])
    ranges, zdiag = range_table(chn, INPUTS)
    big = sort(1:NZ; by = z -> -INPUTS.cumulative[z])[1:5]
    traces = [(chain_matrix(chn, :region_drift_sd_zone), "σ_δ"),
        (chain_matrix(chn, :composition_rho_zone), "ρ"),
        (chain_matrix(chn, :region_halflife_zone), "h (days)")]
    for z in big
        push!(traces, (chain_matrix(chn, :delta_T_zone, z),
            "δ_T " * INPUTS.zone_labels[z]))
    end
    path = trace_panels(traces; path = joinpath(OUT, "c_traces.png"))
    body = """
    <p>NUTS, Mooncake, diagonal metric, target acceptance 0.8, tree depth cap
    8, data-informed start with per-chain jitter. The initial log density is
    that of each chain's jittered start.</p>
    <h3>Per chain</h3>$(html_table(per_chain; digits = 4))
    <h3>Overall</h3>$(html_table(overall; digits = 3))
    <h3>Ranges over zones</h3>$(html_table(ranges))
    <h3>Traces</h3>$(png_tag(path))
    <h3>Per zone</h3>$(html_table(zdiag; digits = 2))
    """
    write_fragment("c", "Real-data fit: sampler diagnostics", body)
end

## --- Stage d: posterior predictive ----------------------------------------

function stage_d(chn)
    log_line("stage d: posterior predictive")
    zone_tags = composition_tags(chn, "d_ppc"; prior_chain = prior_chain(),
        label = "posterior predictive")
    ## Patch-level split of the allocated totals against the parent's
    ## modelled province shares at the parent's own vintages.
    np = length(INPUTS.patch_names)
    nv = length(INPUTS.days)
    obs_patch = fill(NaN, np, nv)
    for v in 1:nv
        tots = [sum(@view INPUTS.counts[zs, v]) for zs in INPUTS.patch_ranges]
        N = sum(tots)
        N > 0 || continue
        obs_patch[:, v] .= tots ./ N
    end
    x = float.(INPUTS.days)
    fig = Figure(; size = (300 * np, 260))
    have_parent = _has_key(PARENT, :province_shares)
    pdays = OBS.province_confirmed_history["ituri"].days
    for p in 1:np
        ax = Axis(fig[1, p]; title = INPUTS.patch_labels[p],
            xlabel = "Grid day", ylabel = "Share of allocated total")
        keep = .!isnan.(obs_patch[p, :])
        scatter!(ax, x[keep], obs_patch[p, keep]; color = :black,
            markersize = 5)
        have_parent || continue
        ms = [collect(v) for v in vec(collect(PARENT[:province_shares]))]
        nvp = min(size(first(ms), 2), length(pdays))
        m = [ms[i][p, v] for i in eachindex(ms), v in 1:nvp]
        xp = float.(pdays[1:nvp])
        band!(ax, xp, vec(mapslices(v -> quantile(v, 0.05), m; dims = 1)),
            vec(mapslices(v -> quantile(v, 0.95), m; dims = 1));
            color = (:firebrick, 0.25))
        lines!(ax, xp, vec(mapslices(median, m; dims = 1));
            color = :firebrick)
    end
    path = joinpath(OUT, "d_ppc_patches.png")
    save(path, fig)
    body = """
    <p>Observed against modelled share of each patch's allocated confirmed
    cases per vintage for the $(TOP) largest zones per patch: posterior
    5–95% and 25–75% bands in blue, the prior predictive bands of stage a in
    grey, observed as points.</p>
    $(zone_tags)
    <p>Patch-level split of the zone tables' allocated totals (points)
    against the parent joint's modelled province shares at the parent's
    own spatial-table vintages (red, 5–95%). The zone model conditions on
    the patch totals, so this is the parent's check, shown for context.</p>
    $(png_tag(path))
    """
    write_fragment("d", "Posterior predictive checks", body)
end

## --- Stage e: ranking and map ---------------------------------------------

function stage_e(chn)
    log_line("stage e: ranking")
    ov = zone_overview_table(chn, INPUTS)
    shown = ov[:, [:zone, :patch, :cases, :share, :R_T, :p_R_above_1,
        :delta_T, :walking]]
    body = "<h3>Zones ranked by P(R_T &gt; 1)</h3>" *
           html_table(shown; digits = 2)
    keep = findall(isfinite, ov.rt_median)
    if isdefined(BVDOutbreakSize, :plot_zone_ranking)
        ranking = DataFrame(label = ov.zone[keep],
            patch = ov.patch_index[keep], rt_median = ov.rt_median[keep],
            rt_lo90 = ov.rt_lo90[keep], rt_hi90 = ov.rt_hi90[keep],
            p_rt_above_one = ov.p_rt_above_one[keep],
            walking = ov.walking[keep])
        fig = BVDOutbreakSize.plot_zone_ranking(ranking;
            patch_labels = INPUTS.patch_labels)
        path = joinpath(OUT, "e_ranking.png")
        save(path, fig)
        body *= "<h3>Ranking</h3>" * png_tag(path; width = "60%")
    end
    if isdefined(BVDOutbreakSize, :plot_zone_map)
        ## The geojson keys a zone without the manifest's province prefix.
        zone_map_keys = [String(last(split(k, "."; limit = 2)))
                         for k in INPUTS.zone_keys]
        rt = [filter(isfinite, r) for r in _zone_draws(chn, :R_T_zone, NZ)]
        zs = findall(!isempty, rt)
        fig = BVDOutbreakSize.plot_zone_map([median(rt[z]) for z in zs],
            zone_map_keys[zs]; lower = [quantile(rt[z], 0.05) for z in zs],
            upper = [quantile(rt[z], 0.95) for z in zs], diverging_at = 1.0,
            scale = log10, title = "Reproduction number at the cut-off",
            colorbar_label = "R")
        path = joinpath(OUT, "e_map.png")
        save(path, fig)
        body *= "<h3>Current R_T by zone</h3>" * png_tag(path; width = "70%")
    end
    write_fragment("e", "Zone ranking and map", body)
end

## --- Run ------------------------------------------------------------------

assemble_index()
"a" in STAGES && stage_a()
if any(s -> s in STAGES, ("c", "d", "e"))
    chn = load_or_fit_chain()
    "c" in STAGES && stage_c(chn)
    "d" in STAGES && stage_d(chn)
    "e" in STAGES && stage_e(chn)
end
"b" in STAGES && stage_b()
log_line("done: $(joinpath(OUT, "index.html"))")
