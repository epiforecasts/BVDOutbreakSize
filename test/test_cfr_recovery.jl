## Parameter recovery for the case-fatality ratio at a truth near the low end
## of past outbreaks. An over-concentrated prior pulls such a truth toward
## its mean: the death series is built deterministically from a fixed CFR
## truth (no sampling noise), so a biased recovery can only come from the
## prior or the delay/ascertainment trade-off, not from chance in the
## simulated counts.

@testitem "deaths_model recovers a low case-fatality ratio" tags = [:slow] begin
    using BVDOutbreakSize: deaths_model, nuts_sample, bin_increments
    using Turing: fix, returned
    using Random: MersenneTwister
    using Statistics: median, quantile

    n = 90
    onsets = fill(30.0, n)
    k = 8.0
    true_cfr = 0.093
    history = (; days = [30, 60, 90], counts = Int[])

    truth_model = fix(
        deaths_model(history, missing, onsets, k);
        cfr_state = (; CFR = true_cfr)
    )
    vi = rand(MersenneTwister(1), truth_model)
    st = returned(truth_model, vi)
    increments = round.(Int, bin_increments(st.deaths_daily, history.days))
    history_obs = (; days = history.days, counts = cumsum(increments))

    model = deaths_model(history_obs, missing, onsets, k)
    chn = nuts_sample(
        model; samples = 500, n_adapts = 500, chains = 1, seed = 2,
        progress = false
    )
    cfr = vec(Array(chn[Symbol("cfr_state.CFR")]))
    lo, hi = quantile(cfr, [0.05, 0.95])
    @test lo <= true_cfr <= hi
    @test isapprox(median(cfr), true_cfr; atol = 0.01)
end
