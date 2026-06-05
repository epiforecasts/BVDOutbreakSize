# Headline joint fit of ALL confirmed cases + confirmed deaths through
# 3 June via the missing-lab-data queue. The late (29 May-3 June) and
# early (18-22 May) dark windows have confirmed counts but no published
# analysed denominator, so they are fitted by the Poisson-thinned queue
# (`confirmed_queue = true`) rather than paired with the 28-May analysed
# total. Serializes the chain and writes a result summary.
#
#   julia --project=test scripts/headline_3june.jl

using BVDOutbreakSize
using BVDOutbreakSize: bvd_joint, load_observations, nuts_sample,
                       confirmed_q_re_model, genetic_seeding_model,
                       report_onset_offset, fit_diagnostics,
                       posterior_summary, tensorboard_callback
using Turing: predict, @varname
using TensorBoardLogger
using Serialization: serialize
using Statistics: median, quantile
import FlexiChains

include(joinpath(@__DIR__, "joint_setup.jl"))

const OUT = joinpath(@__DIR__, "..", "logs", "headline_3june_result.txt")
const CHN = joinpath(@__DIR__, "..", "logs", "chn_3june.jls")

ival90(x) = (round(quantile(x, 0.05); digits = 0),
    round(quantile(x, 0.95); digits = 0))

function main()
    obs = load_observations()
    fa = build_fit_args(obs)
    gs = T -> genetic_seeding_model(T, obs.genetic_tmrca_days;
        tmrca_days_sd = obs.genetic_tmrca_days_sd)
    g = growth_for(obs)

    model = bvd_joint(obs.exported_cases, fa.deaths, fa.reported,
        fa.export_deaths; fa.kw..., growth = g,
        first_export_detection_delta = obs.first_export_detection_delta,
        report_onset_offset = report_onset_offset(obs.as_of_date),
        confirmed_q_random_effect = confirmed_q_re_model, genetic = gs)

    chn = nuts_sample(model; samples = 500, chains = 4, seed = 1,
        target_accept = 0.99,
        callback = tensorboard_callback("logs/headline_3june"))
    serialize(CHN, chn)

    d = fit_diagnostics(chn)
    ess_tail = FlexiChains.ess(chn; kind = :tail)
    et = filter(isfinite,
        collect(Iterators.flatten(values(ess_tail))))
    min_ess_tail = isempty(et) ? NaN : minimum(et)

    C = vec(Array(chn[:cumulative_cases]))
    Cinf = vec(Array(chn[:cumulative_infections]))

    ## Posterior predictive: per-vintage confirmed cases and confirmed
    ## deaths. Rebuild the model with those streams set to `missing` so
    ## `predict` simulates replicate counts at every vintage.
    nconf = length(fa.kw.confirmed_cases)
    ncd = length(fa.kw.confirmed_deaths)
    pp_kw = (; fa.kw...,
        confirmed_cases = fill(missing, nconf),
        confirmed_deaths = fill(missing, ncd))
    pp_model = bvd_joint(missing, fa.deaths, fa.reported,
        fa.export_deaths; pp_kw..., growth = g,
        first_export_detection_delta = obs.first_export_detection_delta,
        report_onset_offset = report_onset_offset(obs.as_of_date),
        confirmed_q_random_effect = confirmed_q_re_model, genetic = gs)
    pp = predict(pp_model, chn)

    ## Per-vintage confirmed-case increments -> cumulative for comparison.
    conf_draws = pp[@varname(confirmed_cases)]   # vector of length-n vecs
    cd_draws = pp[@varname(confirmed_deaths)]
    cconf_obs = cumsum(Int.(fa.kw.confirmed_cases))
    ccd_obs = cumsum(Int.(fa.kw.confirmed_deaths))

    ## Build a draws x nvintage matrix of cumulative predicted counts.
    function cum_quantiles(draws, n)
        cum = [cumsum(collect(d)) for d in vec(draws)]
        med = [median(getindex.(cum, w)) for w in 1:n]
        lo = [quantile(getindex.(cum, w), 0.05) for w in 1:n]
        hi = [quantile(getindex.(cum, w), 0.95) for w in 1:n]
        return (med = med, lo = lo, hi = hi)
    end
    cq = cum_quantiles(conf_draws, nconf)
    cdq = cum_quantiles(cd_draws, ncd)

    open(OUT, "w") do io
        println(io, "Headline joint fit: ALL confirmed cases + confirmed ",
            "deaths through 3 June via the missing-lab-data queue")
        println(io, "as_of_date: ", obs.as_of_date,
            "  (confirmed_queue = true, e = 0 headline)")
        println(io, "samples = 500, chains = 4, target_accept = 0.99, ",
            "seed = 1")
        println(io, "")
        println(io, "=== Convergence ===")
        println(io, "max R-hat:     ", round(d.max_rhat; digits = 3))
        println(io, "min ESS bulk:  ", round(d.min_ess_bulk; digits = 1))
        println(io, "min ESS tail:  ", round(min_ess_tail; digits = 1))
        println(io, "divergences:   ", d.n_divergent)
        println(io, "")
        println(io, "=== Outbreak size C(T) ===")
        clo, chi = ival90(C)
        ilo, ihi = ival90(Cinf)
        println(io, "cumulative infections median: ",
            round(median(Cinf); digits = 0), "  90%: ", ilo, " to ", ihi)
        println(io, "cumulative cases median:      ",
            round(median(C); digits = 0), "  90%: ", clo, " to ", chi)
        println(io, "")
        println(io, "=== Confirmed cases per vintage (cumulative) ===")
        println(io, "offsets:  ", fa.kw.confirmed_offsets)
        println(io, "observed: ", cconf_obs)
        println(io, "PP med:   ", round.(cq.med; digits = 0))
        println(io, "PP 5%:    ", round.(cq.lo; digits = 0))
        println(io, "PP 95%:   ", round.(cq.hi; digits = 0))
        println(io, "(late dark windows 29 May-3 June are the last 6: ",
            "263, 282, 321, 344, 363, 381)")
        println(io, "")
        println(io, "=== Confirmed deaths per vintage (cumulative) ===")
        println(io, "offsets:  ", fa.kw.confirmed_death_offsets)
        println(io, "observed: ", ccd_obs)
        println(io, "PP med:   ", round.(cdq.med; digits = 0))
        println(io, "PP 5%:    ", round.(cdq.lo; digits = 0))
        println(io, "PP 95%:   ", round.(cdq.hi; digits = 0))
        println(io, "confirmed-deaths cum observed (final): ", ccd_obs[end],
            "  PP med: ", round(cdq.med[end]; digits = 0),
            "  90%: ", round(cdq.lo[end]; digits = 0), " to ",
            round(cdq.hi[end]; digits = 0))
        flush(io)
    end
    println("Wrote ", OUT)
    println("max R-hat=", round(d.max_rhat; digits = 3),
        " min ESS bulk=", round(d.min_ess_bulk; digits = 1),
        " div=", d.n_divergent)
end

main()
