# Delay-prior configurations for the line-list refits: which distributions the
# fit is given for the generation interval and the onset-to-report delay, and
# where those numbers come from.
#
# Shared by scripts/linelist/fit_single.jl, which passes the result through
# `build_fit_specs(obs; delays = ...)`.
#
# `include` this, do not run it.
#
# ## Where the numbers live, and why not here
#
# The fresh estimates are fitted in `bvd-internal-cmmid`, which is private. Its
# DISCLOSURE.md permits fitted distribution parameters to be committed *there*,
# "without onward sharing". This repository is public, so the values are read at
# run time from `BVD_DELAY_DIR` and never written into a tracked file. What is
# tracked here is the arithmetic that turns those estimates into priors.
#
# ## Which delays are relevant
#
# Traced through the two composers this is used with:
#
#   cases_only_model  -> infection_model      generation interval
#                     -> onset_incidence_model  incubation
#                     -> reported_cases_model   onset_to_report
#
#   onsets_only_model -> infection_model      generation interval
#                     -> onset_incidence_model  incubation
#                     -> onset_reporting_model  a *fitted* hazard, no prior
#
# So the generation interval is swappable in both, `onset_to_report` in `cases`
# only, and the incubation period in neither: cmmid's own incubation fits are
# unidentified (means of 18-70 days on 17-64 contacts) and it recommends the
# MacNeil et al. (2010) estimate, which is already this repository's prior. The
# onsets fit *estimates* the delay cmmid fits, so it is a comparison rather than
# an input.
#
# A report-delay config handed to the onsets fit would therefore be a silent
# no-op. `check_delay_config` makes it an error instead.

using BVDOutbreakSize: cdf_nmax, lognormal_meansd
using CSV
using DataFrames
using Distributions: Gamma, Normal, truncated
using SHA: sha256

## The generation-interval row in the cmmid tables that each config name uses.
## Keys are the config suffix; values are the `source` column value shared by
## `gi_estimates.csv` and `si_estimates.csv`.
const GI_SOURCES = Dict(
    "any" => "any matched contact",
    "case" => "secondary a case",
    "diag" => "secondary diagnosed")

const DELAY_CONFIGS = ["repo",
    "cmmid_gi_any", "cmmid_gi_case", "cmmid_gi_diag",
    "cmmid_any", "cmmid_case", "cmmid_diag",
    "cmmid_rep"]

## Where the cmmid results tables are read from. Required rather than defaulted,
## for the disclosure reason above: the operator says where the private results
## are instead of a default quietly finding a directory.
function delay_dir()
    dir = get(ENV, "BVD_DELAY_DIR", "")
    isempty(dir) && error("set BVD_DELAY_DIR to the bvd-internal-cmmid " *
          "`results/` directory holding gi_estimates.csv, si_estimates.csv " *
          "and bayes_pooled_parameters.csv. See scripts/linelist/README.md.")
    isdir(dir) || error("BVD_DELAY_DIR is not a directory: $dir")
    return abspath(dir)
end

## Prior-width multiplier. The cmmid standard errors come from large samples
## (n = 3710 for the report delay), so the priors below are near-fixed and carry
## almost no delay uncertainty into R_t. This is the knob for a run that widens
## them, so "how much does the answer depend on the delay being exactly this"
## is one flag away rather than an edit.
delay_prior_inflate() = parse(Float64,
    get(ENV, "BVD_DELAY_PRIOR_INFLATE", "1"))

function _read_table(dir, name, cols)
    path = joinpath(dir, name)
    isfile(path) || error("missing $path. It is a bvd-internal-cmmid " *
          "results table; BVD_DELAY_DIR should point at that repository's " *
          "`results/` directory.")
    df = CSV.read(path, DataFrame)
    missing_cols = [c for c in cols if !(c in names(df))]
    isempty(missing_cols) ||
        error("$path is missing column(s) " * join(missing_cols, ", ") *
              "; found " * join(names(df), ", "))
    return df
end

function _one_row(df, col, value, path)
    rows = df[df[!, col] .== value, :]
    nrow(rows) == 1 ||
        error("expected exactly one row with $col = \"$value\" in $path, " *
              "found $(nrow(rows))")
    return rows[1, :]
end

## Generation interval, as a Gamma shape/scale prior pair.
##
## cmmid publishes the generation interval as a point moment-match: it
## deconvolves the fitted serial interval by the incubation variance and matches
## the result to a Gamma, so `shape` and `scale` come with no interval of their
## own. The uncertainty that *is* published is on the serial interval's
## location, as a 95% CI on `meanlog`.
##
## That is carried through here rather than replaced by a spread of our own,
## which is how this repository built its own generation-interval prior (priors
## centred on the source's point values, spreads from the source's reported
## uncertainty on the mean). For a lognormal the relative SE of the mean is the
## SE of `meanlog`, so `σ_m = gi_mean * se_meanlog`; the delta method then moves
## it onto the Gamma parameters with the SD held at its point value, since the
## deconvolution fixes the variance and only the location is being propagated:
##
##   α = (m/s)^2   ->  ∂α/∂m = 2m/s^2
##   θ = s^2/m     ->  ∂θ/∂m = -s^2/m^2
##
## The induced α-θ correlation is dropped, as this repository's own
## generation-interval prior also drops it.
function _gi_prior(dir, suffix; inflate)
    source = GI_SOURCES[suffix]
    gi_path = joinpath(dir, "gi_estimates.csv")
    si_path = joinpath(dir, "si_estimates.csv")
    gi = _read_table(dir, "gi_estimates.csv",
        ["source", "gi_mean", "gi_sd", "shape", "scale"])
    si = _read_table(dir, "si_estimates.csv",
        ["source", "meanlog_lower", "meanlog_upper"])

    g = _one_row(gi, "source", source, gi_path)
    s = _one_row(si, "source", source, si_path)

    m, sd = Float64(g.gi_mean), Float64(g.gi_sd)
    α, θ = Float64(g.shape), Float64(g.scale)
    se_meanlog = (Float64(s.meanlog_upper) - Float64(s.meanlog_lower)) /
                 (2 * 1.96)
    σ_m = m * se_meanlog
    σ_α = 2 * m * σ_m / sd^2
    σ_θ = sd^2 * σ_m / m^2

    ## `lower` matches this repository's own generation-interval prior, which
    ## truncates both parameters at 0.1 to keep the Gamma well defined.
    return (; alpha_prior = truncated(Normal(α, σ_α * inflate); lower = 0.1),
        theta_prior = truncated(Normal(θ, σ_θ * inflate); lower = 0.1),
        nmax = cdf_nmax(Gamma(α, θ)),
        source, alpha = α, theta = θ, sd_alpha = σ_α * inflate,
        sd_theta = σ_θ * inflate, mean = m, sd = sd)
end

## Onset-to-report, as a mean/SD prior pair on a LogNormal.
##
## cmmid fits onset -> line-list classification as a lognormal, selected over
## gamma and Weibull by both AIC and LOOIC, with double interval censoring and
## per-observation right truncation. `bayes_pooled_parameters.csv` carries it as
## a brms fit: `Intercept` is `meanlog`, and `sigma_Intercept` is `log(sdlog)`.
##
## `censored_delay_model` builds `lognormal_meansd(mean, sd)`, so centring the
## mean and SD priors on the moment-matched values reproduces the cmmid
## lognormal exactly at the prior centre — no re-parameterisation and no new
## submodel. The SEs move across by the delta method:
##
##   mean = exp(μ + σ²/2)        ∂/∂μ = mean          ∂/∂σ = mean·σ
##   sd   = mean·k, k = √(e^{σ²}−1)
##                               ∂/∂μ = sd
##                               ∂/∂σ = mean·σ·k + mean·σ·e^{σ²}/k
##
## This delay is matched to the *known-by* streams, which are indexed by
## `date_ll_classification`. It does not describe the notification-date streams:
## cmmid's onset -> notification is a different, shorter interval (national
## median 3 days) and is published empirically, without a truncation correction.
function _report_prior(dir; inflate)
    path = joinpath(dir, "bayes_pooled_parameters.csv")
    df = _read_table(dir, "bayes_pooled_parameters.csv",
        ["term", "Estimate", "Est.Error"])
    int = _one_row(df, "term", "Intercept", path)
    sig = _one_row(df, "term", "sigma_Intercept", path)

    μ = Float64(int.Estimate)
    se_μ = Float64(int[Symbol("Est.Error")])
    ## brms models `sigma` on the log scale, so `sdlog` is the exponential of
    ## the reported coefficient and its SE moves across by the same factor.
    σ = exp(Float64(sig.Estimate))
    se_σ = σ * Float64(sig[Symbol("Est.Error")])

    mean = exp(μ + σ^2 / 2)
    k = sqrt(exp(σ^2) - 1)
    sd = mean * k
    se_mean = hypot(mean * se_μ, mean * σ * se_σ)
    se_sd = hypot(sd * se_μ, (mean * σ * k + mean * σ * exp(σ^2) / k) * se_σ)

    ## `lower = 1` matches the other `censored_delay_model` priors in the
    ## package: a mean or SD below a day is outside what a daily grid resolves.
    return (; mean_prior = truncated(Normal(mean, se_mean * inflate);
            lower = 1),
        sd_prior = truncated(Normal(sd, se_sd * inflate); lower = 1),
        nmax = cdf_nmax(lognormal_meansd(mean, sd)),
        meanlog = μ, sdlog = σ, mean, sd,
        sd_mean = se_mean * inflate, sd_sd = se_sd * inflate)
end

## The configuration a fit is given. `repo` returns `nothing`, meaning "pass no
## overrides": the repository's own priors are then the code defaults by
## construction rather than by a transcription here that could drift from them.
function delay_config(name)
    name in DELAY_CONFIGS ||
        error("unknown delay config `$name`; expected one of " *
              join(DELAY_CONFIGS, ", "))
    name == "repo" && return nothing

    dir = delay_dir()
    inflate = delay_prior_inflate()

    ## `cmmid_rep` swaps the report delay and nothing else, so the generation
    ## interval stays at the package default. It exists because the cmmid
    ## generation interval is fitted from line-list transmission pairs, and a
    ## comparison that wants the report delay without taking a view on that
    ## generation interval has no other configuration to ask for: every
    ## `cmmid_*` and `cmmid_gi_*` name carries it. `gi` is `nothing` here, and
    ## the splat in docs/fits/registry.jl reads that as "no override" rather
    ## than restating the package default, which could drift from it.
    if name == "cmmid_rep"
        report = _report_prior(dir; inflate)
        return (; name, dir, inflate, has_report = true,
            gi_alpha_prior = nothing, gi_theta_prior = nothing,
            gi_nmax = nothing,
            report_mean_prior = report.mean_prior,
            report_sd_prior = report.sd_prior, report_nmax = report.nmax,
            gi = nothing, report)
    end

    suffix = String(last(split(name, "_")))
    haskey(GI_SOURCES, suffix) ||
        error("no generation-interval source for config `$name`")
    gi = _gi_prior(dir, suffix; inflate)

    ## `cmmid_gi_*` swaps the generation interval alone. `cmmid_*` swaps the
    ## report delay too. Running both is what separates the two effects: the
    ## report delay moves the onsets->reports mapping by about nine days, so a
    ## config that changes it and the generation interval together cannot say
    ## which did the work.
    has_report = !startswith(name, "cmmid_gi_")
    report = has_report ? _report_prior(dir; inflate) : nothing

    return (; name, dir, inflate, has_report,
        gi_alpha_prior = gi.alpha_prior, gi_theta_prior = gi.theta_prior,
        gi_nmax = gi.nmax,
        report_mean_prior = isnothing(report) ? nothing : report.mean_prior,
        report_sd_prior = isnothing(report) ? nothing : report.sd_prior,
        report_nmax = isnothing(report) ? nothing : report.nmax,
        gi, report)
end

## A report-delay override reaches `reported_cases_model`, which the onsets fit
## does not run: `onsets_only_model` composes the latent process with
## `onset_reporting_model`, whose delay is a fitted hazard rather than a prior.
## Passed there the override would be silently dropped and the run would look
## like a result. Refuse instead.
function check_delay_config(delays, fit)
    isnothing(delays) && return nothing
    if fit == "onsets" && delays.has_report
        error("delay config `$(delays.name)` overrides the onset-to-report " *
              "delay, but the `onsets` fit estimates that delay from the " *
              "reporting triangle rather than taking a prior, so the " *
              "override would be silently ignored. Use `cmmid_gi_*` for the " *
              "onsets fit, or compare its fitted delay against cmmid's.")
    end
    return nothing
end

## Every number the config carries, with the SHA-256 of each file it came from.
## Written beside the fit outputs (git-ignored) so a chain on disk can be traced
## back to the estimates that produced it, which the values themselves cannot
## be once they are inside a prior object.
function write_delay_provenance(delays, out)
    isnothing(delays) && return nothing
    rows = NamedTuple[]
    push!(rows,
        (; config = delays.name, item = "inflate",
            value = string(delays.inflate), source = "BVD_DELAY_PRIOR_INFLATE",
            sha256 = ""))

    files = String[]
    isnothing(delays.gi) || append!(files, ["gi_estimates.csv", "si_estimates.csv"])
    delays.has_report && push!(files, "bayes_pooled_parameters.csv")
    for f in files
        path = joinpath(delays.dir, f)
        push!(rows,
            (; config = delays.name, item = "file", value = f,
                source = path, sha256 = bytes2hex(sha256(read(path)))))
    end

    g = delays.gi
    isnothing(g) || for (k, v) in (("gi_source", g.source), ("gi_alpha", g.alpha),
        ("gi_alpha_sd", g.sd_alpha), ("gi_theta", g.theta),
        ("gi_theta_sd", g.sd_theta), ("gi_mean", g.mean), ("gi_sd", g.sd),
        ("gi_nmax", delays.gi_nmax))
        push!(rows,
            (; config = delays.name, item = k, value = string(v),
                source = "gi_estimates.csv + si_estimates.csv", sha256 = ""))
    end

    if delays.has_report
        r = delays.report
        for (k, v) in (("report_meanlog", r.meanlog), ("report_sdlog", r.sdlog),
            ("report_mean", r.mean), ("report_mean_sd", r.sd_mean),
            ("report_sd", r.sd), ("report_sd_sd", r.sd_sd),
            ("report_nmax", delays.report_nmax))
            push!(rows,
                (; config = delays.name, item = k, value = string(v),
                    source = "bayes_pooled_parameters.csv", sha256 = ""))
        end
    end

    path = joinpath(out, "delay_provenance_$(delays.name).csv")
    CSV.write(path, DataFrame(rows))
    @info "delay provenance written" path
    return path
end
