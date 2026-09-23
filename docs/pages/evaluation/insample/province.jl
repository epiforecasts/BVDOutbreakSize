# # Province in-sample checks
#
# Whether the fitted joint model reproduces the per-province data it was fitted to.
# The national checks are on the [national in-sample checks](@ref "In-sample checks") page.
# The per-province confirmed cases and deaths enter the model as compositions of the national totals, so the in-sample checks here are on each province's share of those totals rather than its count.

#md # ```@raw html
#md # <details><summary>Load packages, data and fitted chains</summary>
#md # ```

## Shared setup: packages, observations, the fit registry and every model fit
## (loaded from the content-addressed cache). See `docs/pages/_setup.jl`.
using BVDOutbreakSize
include(joinpath(pkgdir(BVDOutbreakSize), "docs", "pages", "_setup.jl"))

#md # ```@raw html
#md # </details>
#md # ```

# ## Summary
#
# The overall bullets come first, then a short block per province, from the checks further down this page.
# Each block gives how well the case and death compositions reproduce that province's counts.
# Bias runs from −1 to 1 and is zero when the observed counts sit at the predictive median, negative when the model under-predicts.
# Coverage is nominally 0.9.

#md # ```@eval
#md # using Markdown, BVDOutbreakSize
#md # dir = joinpath(pkgdir(BVDOutbreakSize), "docs", "src", "summary_assets")
#md # Markdown.parse(read(joinpath(dir, "evaluation_insample_province.md"), String))
#md # ```

# ## Province prior predictive check
#
# Before any observation is taken into account, what does the prior imply about each province's share of the confirmed cases and deaths?
# The shared prior on the national page is drawn from a single population and carries no province quantities.
# The draws here come from the four-patch model's prior instead, with the province grids in place and every observation left out of the density.
# The prior shares should bracket the observed ones without pinning them.

#md # ```@raw html
#md # <details><summary>Draw from the four-patch prior</summary>
#md # ```

## The province compositions are passed their observed counts rather than
## `missing`. `Prior()` leaves every observation out of the density, so the
## counts only fix each vintage's total. With `missing` the predictive path
## draws the total from the modelled increments instead, which an extreme
## prior draw can push past the integer range.
prior_patch_chn = let
    m = bvd_joint(
        obs.n, missing, missing, missing, missing, missing;
        deaths_history = (; days = Int[], counts = Int[]),
        reported_history = (; days = Int[], counts = Int[]),
        confirmed_history = (; days = Int[], counts = Int[]),
        export_case_days = obs.export_case_days,
        export_death_days = obs.export_death_days,
        breakpoint = obs.n - obs.who_first_sitrep_days,
        background_pooling = background_pooling_model,
        genetic = genetic_seeding_model,
        tmrca_days = obs.tmrca_days,
        n_patches = N_PATCHES,
        province_increments = province_cases.increments,
        province_days = province_cases.days,
        province_testing_covariate = province_testing,
        province_death_increments = province_deaths.increments,
        province_death_days = province_deaths.days
    )
    sample(m, Prior(), 1_000; progress = false)
end;

prior_province_table = patch_overview_table(prior_patch_chn, N_PATCHES);

#md # ```@raw html
#md # </details>
#md # ```

#md # ```@raw html
#md # <details><summary>Show prior province summary table</summary>
#md # ```

MarkdownTable(prior_province_table) #hide

#md # ```@raw html
#md # </details>
#md # ```

# The prior composition checks draw the same bands as the [posterior ones](@ref province-compositions) below, from the prior shares.

#md # ```@raw html
#md # <details><summary>Province composition prior predictive checks</summary>
#md # ```

prior_case_ppc_fig = plot_province_composition_ppc(
    prior_patch_chn;
    share_key = :province_shares,
    obs_increments = province_cases.increments,
    days = province_cases.days, seeding = obs.seeding, n_patches = N_PATCHES,
    title = "Confirmed case share by province, prior"
);

prior_death_ppc_fig = plot_province_composition_ppc(
    prior_patch_chn;
    share_key = :province_death_shares,
    obs_increments = province_deaths.increments,
    days = province_deaths.days, seeding = obs.seeding, n_patches = N_PATCHES,
    title = "Confirmed death share by province, prior"
);

#md # ```@raw html
#md # </details>
#md # ```

prior_case_ppc_fig #hide

prior_death_ppc_fig #hide

# The pair plot sets the prior against the posterior for each province's reproduction number at the cut-off and its cumulative infections, the latter on the log scale.

#md # ```@raw html
#md # <details><summary>Province prior and posterior pair plot</summary>
#md # ```

## Per-province entries of a vector deterministic, one draw vector each,
## keyed `<name>_<patch>` so they can sit side by side in one table.
function _patch_columns(chn, key::Symbol, name::AbstractString; f = identity)
    draws = [collect(v) for v in vec(collect(chn[key]))]
    return [
        Symbol(name, "_", p) => Float64[f(d[p]) for d in draws]
            for p in 1:N_PATCHES
    ]
end
_province_pair_draws(chn) = NamedTuple(
    vcat(
        _patch_columns(chn, :R_T_patch, "R_T"),
        _patch_columns(chn, :C_T_patch, "log10_C_T"; f = x -> log10(x + 1))
    )
)
_province_pair_labels = Dict(
    Symbol(name, "_", p) => string(label, ", ", PROVINCE_LABELS[p])
        for p in 1:N_PATCHES
        for (name, label) in (("R_T", "R_T"), ("log10_C_T", "log10 C_T"))
)
province_pair_fig = plot_pair(
    _province_pair_draws(chn_joint);
    prior = _province_pair_draws(prior_patch_chn),
    labels = _province_pair_labels
);

#md # ```@raw html
#md # </details>
#md # ```

province_pair_fig #hide

# ## [Province compositions](@id province-compositions)
#
# The per-province confirmed cases and deaths are fitted as compositions conditional on the national total, so what the model predicts is each province's share rather than its count.
# The panels below show that modelled share at every spatial vintage against the observed one.
# Each panel carries two bands.
# The grey band is the posterior predictive interval on the observed share, built by pushing every posterior draw's expected shares back through the composition's own overdispersed allocation at that vintage's observed total.
# The overdispersion is what absorbs reporting lags between the provincial and national tables and the reassignment of cases between health zones.
# The observed points should fall inside it.
# The coloured ribbon inside the grey band is the expected share alone, which is the modelled centre the points scatter around.
# A point outside the grey band is a vintage the composition does not reproduce, and points consistently to one side of the coloured ribbon are a province the model splits wrongly on average.
# Each panel starts at zero and takes its own upper limit, because the shares differ by orders of magnitude.
# The vintages stop before the cut-off, so the panels end earlier than the [national posterior predictive checks](@ref "Posterior predictive checks").

#md # ```@raw html
#md # <details><summary>Province composition posterior predictive checks</summary>
#md # ```

province_case_ppc_fig = plot_province_composition_ppc(
    chn_joint;
    share_key = :province_shares,
    obs_increments = province_cases.increments,
    days = province_cases.days, seeding = obs.seeding, n_patches = N_PATCHES,
    title = "Confirmed case share by province"
);

province_death_ppc_fig = plot_province_composition_ppc(
    chn_joint;
    share_key = :province_death_shares,
    obs_increments = province_deaths.increments,
    days = province_deaths.days, seeding = obs.seeding, n_patches = N_PATCHES,
    title = "Confirmed death share by province"
);

#md # ```@raw html
#md # </details>
#md # ```

province_case_ppc_fig #hide

province_death_ppc_fig #hide

# ## Province stream calibration
#
# Each province's cases and deaths are scored vintage by vintage against the same composition predictive as the grey bands above, as counts at the observed national total.
# The columns are those of the national [stream calibration](@ref "Stream calibration").
# A vintage with no confirmed cases or deaths nationally has no split to score and is left out.
# The check is on the split alone, so a province can be well calibrated here while the national total it is a share of is not.

#md # ```@raw html
#md # <details><summary>Build the province calibration panels</summary>
#md # ```

province_case_panels = province_composition_panels(
    chn_joint;
    share_key = :province_shares,
    obs_increments = province_cases.increments,
    stream = "Confirmed cases", n_patches = N_PATCHES
);
province_death_panels = province_composition_panels(
    chn_joint;
    share_key = :province_death_shares,
    obs_increments = province_deaths.increments,
    stream = "Confirmed deaths", n_patches = N_PATCHES
);
province_panels = vcat(province_case_panels, province_death_panels);
province_calibration_table = stream_calibration(province_panels);
province_calibration_fig = plot_stream_calibration(
    province_calibration_table
);

#md # ```@raw html
#md # </details>
#md # ```

province_calibration_fig #hide

#md # ```@raw html
#md # <details><summary>Province calibration table</summary>
#md # ```

MarkdownTable(province_calibration_table) #hide

#md # ```@raw html
#md # </details>
#md # ```

# ## Province correlations and totals
#
# The heatmap is the posterior correlation between the national outbreak size ($C_T$) and, for each province, its reproduction number at the cut-off, its relative case ascertainment and its case-fatality ratio.
# The case composition identifies only the product of a province's ascertainment and its incidence, so a strong negative correlation between the two is expected, and the deaths are what separate them.
# Blue is positive, red negative.

#md # ```@raw html
#md # <details><summary>Province posterior correlation heatmap</summary>
#md # ```

## LaTeX label for one province's entry of a quantity. Spaces in a province
## name need escaping inside `\mathrm`.
_province_tex(sym, p) = string(
    sym, raw"\ (\mathrm{", replace(PROVINCE_LABELS[p], " " => raw"\ "), "})"
)
province_correlation_draws = NamedTuple(
    vcat(
        [:C_T => vec(Array(chn_joint[:C_T]))],
        _patch_columns(chn_joint, :R_T_patch, "R_T"),
        _patch_columns(chn_joint, :province_ascertainment, "asc"),
        _patch_columns(chn_joint, :CFR_patch, "CFR")
    )
);
province_correlation_labels = merge(
    Dict(:C_T => raw"C_T"),
    Dict(
        Symbol(name, "_", p) => _province_tex(tex, p)
            for p in 1:N_PATCHES
            for (name, tex) in (
                ("R_T", raw"R_T"), ("asc", raw"p_\mathrm{asc}"),
                ("CFR", raw"\mathrm{CFR}"),
            )
    )
);
province_correlation_fig = plot_correlation_heatmap(
    province_correlation_draws; labels = province_correlation_labels
);

#md # ```@raw html
#md # </details>
#md # ```

province_correlation_fig #hide

# The totals plot takes each posterior predictive draw from the calibration panels above, sums each province's cases and deaths over the spatial vintages, and marks the observed totals with a crosshair.
# The diagonal panels are the predictive spread of each total against the observed value.
# Within cases, and within deaths, the provinces split a national total that is held at its observed value, so their totals sum to that value in every draw.
# The off-diagonal panels between two provinces' cases, or two provinces' deaths, therefore trade off by construction.
# The panels pairing one province's cases with its deaths are the ones that show whether the two compositions agree.

#md # ```@raw html
#md # <details><summary>Province totals against observed</summary>
#md # ```

## Keyed by panel title, which names the stream and the province. A total
## that is the same in every draw (a patch with no deaths at any vintage) has
## no spread to plot and gives the corner plot a zero-width axis, so it is
## left out.
_varying_panels = filter(
    p -> length(unique(sum(r) for r in p.replicates)) > 1, province_panels
);
province_totals = NamedTuple(
    Symbol(p.title) => [sum(Float64.(r)) for r in p.replicates]
        for p in _varying_panels
);
province_observed = NamedTuple(
    Symbol(p.title) => Float64(sum(p.observed)) for p in _varying_panels
);
province_pairs_fig = plot_stream_pairs(province_totals, province_observed);

#md # ```@raw html
#md # </details>
#md # ```

province_pairs_fig #hide
# ## Saving province in-sample outputs

#md # ```@raw html
#md # <details><summary>Write the summary bullets</summary>
#md # ```

## The bullets under the summary heading at the top of the page. They read
## tables built further down, so they are written here and read back when
## the site is assembled.
evaluation_insample_province_summary = let
    fmt(x) = ismissing(x) || !isfinite(x) ? "n/a" :
        string(round(x; digits = 2))
    fitted = filter(r -> isfinite(r["Bias"]), province_calibration_table)
    worst = first(
        sort(fitted, "Bias"; by = abs, rev = true), min(3, size(fitted, 1))
    )
    low = first(sort(fitted, "90% coverage"), min(1, size(fitted, 1)))
    n_cov = count(>=(0.8), fitted[!, "90% coverage"])
    overall = [
        string(
            "- **Least well reproduced:** ",
            join(
                [
                    string(
                        r["Stream"], " (bias ", fmt(r["Bias"]),
                        ", 90% coverage ", fmt(r["90% coverage"]), ")"
                    )
                        for r in eachrow(worst)
                ], "; "
            ), "."
        ),
        string(
            "- **Coverage:** ", n_cov, " of ", size(fitted, 1),
            " province streams have 90% coverage of at least 0.8",
            size(low, 1) == 0 ? "." :
                string(
                    "; the lowest is ", low[1, "Stream"], " at ",
                    fmt(low[1, "90% coverage"]), "."
                )
        ),
    ]
    cal = Dict(r["Stream"] => r for r in eachrow(province_calibration_table))
    function calibration(kind, p)
        r = get(cal, string(kind, ", ", PROVINCE_LABELS[p]), nothing)
        r === nothing && return string("- ", kind, ": not scored.")
        return string(
            "- ", kind, ": bias ", fmt(r["Bias"]), " and 90% coverage ",
            fmt(r["90% coverage"]), " over ", r["Vintages"], " vintages."
        )
    end
    detail = [
        join(
            [
                string("**", PROVINCE_LABELS[p], "**"), "",
                calibration("Confirmed cases", p),
                calibration("Confirmed deaths", p),
            ], "\n"
        )
            for p in 1:N_PATCHES
    ]
    join(vcat([join(overall, "\n")], detail), "\n\n")
end
dashboard_dir = joinpath(
    pkgdir(BVDOutbreakSize), "docs", "src", "summary_assets"
)
mkpath(dashboard_dir)
write(joinpath(dashboard_dir, "evaluation_insample_province.md"), evaluation_insample_province_summary);

#md # ```@raw html
#md # </details>
#md # ```
