## Tests for the isolation/treatment-bed occupancy stream ("Patients en
## isolement"), a prevalence (length-of-stay) stream: the suspect inflow
## (BVD treatment stay plus non-BVD rule-out stay) carried through a
## length-of-stay survival into a daily stock, scored against the modelled
## occupancy on each report day. Exercised through the `convolve_survival`
## helper, `treatment_only_model` and `bvd_joint`.

@testitem "convolve_survival: same-day discharge returns the inflow" begin
    using BVDOutbreakSize: convolve_survival
    x = [1.0, 2.0, 3.0, 4.0]
    ## A length-of-stay point mass at 0 (`pmf = [1]`) gives survival
    ## `S(0) = 1`, `S(τ>0) = 0`, so the admission day is the only occupancy
    ## day and the occupancy equals the inflow.
    @test convolve_survival(x, [1.0]) == x
end

@testitem "convolve_survival: fixed stay accumulates the right occupancy" begin
    using BVDOutbreakSize: convolve_survival
    ## A length-of-stay fixed at 2 days (`pmf = [0, 0, 1]`) means a patient
    ## occupies a bed on the admission day and the next two days, so the
    ## survival weights are `S(0)=S(1)=S(2)=1`, `S(τ≥3)=0`. With one
    ## admission per day the occupancy ramps 1, 2, 3 then holds at 3.
    los = [0.0, 0.0, 1.0]
    x = ones(5)
    occ = convolve_survival(x, los)
    @test occ == [1.0, 2.0, 3.0, 3.0, 3.0]
    ## Total occupancy equals the total inflow times `E[LOS] + 1` (here 3).
    @test sum(convolve_survival([0.0, 0.0, 1.0, 0.0, 0.0], los)) ≈ 3.0
end

@testitem "convolve_survival: survival weights are non-increasing" begin
    using BVDOutbreakSize: convolve_survival, discretise_censored,
                           lognormal_meansd
    ## A single unit admission on day 1 traces the survival curve directly:
    ## occupancy[t] = S(t-1), which must be non-increasing and start at 1.
    los = discretise_censored(lognormal_meansd(6.0, 4.0), 30)
    x = zeros(40)
    x[1] = 1.0
    occ = convolve_survival(x, los)
    @test occ[1] ≈ 1.0
    @test all(diff(occ) .<= 1e-10)
    @test all(occ .>= -1e-12)
end

@testitem "censoring_cap: bound from recorded capacity, never below obs" begin
    using BVDOutbreakSize: censoring_cap

    iso_days = [9, 10, 11]
    iso_obs = [260, 262, 315]
    capacity_history = (; days = [9, 10, 11], counts = [400, 410, 452])
    ## Each day takes its recorded capacity (all above the observed count), so
    ## the bound is the data, fixed, and never below the occupancy.
    cap = censoring_cap(iso_days, iso_obs, capacity_history)
    @test cap == [400.0, 410.0, 452.0]
    @test all(cap .>= iso_obs)

    ## A capacity below the count is floored at the count, so the censored NB
    ## never has zero probability on the observation (the -Inf wall is removed).
    low_cap = (; days = [9], counts = [200])
    @test censoring_cap([9], [260], low_cap) == [260.0]

    ## No recorded capacity gives a large finite no-op bound (not Inf, which
    ## safe_rate would map to eps), so the censoring is inactive.
    nocap = censoring_cap([9], [260], (; days = Int[], counts = Int[]))
    @test only(nocap) > 1.0e5

    ## The predictive generator passes the capacity history days-only (counts
    ## emptied) and the occupancy missing; the guard is on the counts, so this
    ## must not index the empty counts vector and returns the no-op cap.
    daysonly = (; days = [9, 10, 11], counts = Int[])
    @test only(censoring_cap([10], missing, daysonly)) > 1.0e5

    ## Nearest-day fill for a report day without its own capacity entry.
    @test censoring_cap([12], [300], capacity_history) == [452.0]
end

@testitem "two-day reporting gap: increments bin across the missing day" begin
    using BVDOutbreakSize: bin_increments, load_observations
    using Dates: Date, value

    ## A missing report day leaves a hole in the per-vintage day-indices (the
    ## latent daily grid is unbroken; only the observations skip a day). The
    ## vintage after the hole is two grid days on, so its increment bin must
    ## span both the missing day's and the report day's modelled values — a
    ## two-day bin — not assume a one-day step. With `daily[t] = t` the bin
    ## value reveals which grid days it summed.
    daily = Float64.(collect(1:10))
    days = [6, 7, 9]             # day 8 absent (no report); 9 is two days on
    binned = bin_increments(daily, days)
    @test binned[end] == daily[8] + daily[9]  # spans the missing + report day
    @test binned[2] == daily[7]               # the unbroken steps stay one-day
    ## Mass conservation: no grid day is double-counted or dropped across the
    ## hole.
    @test sum(binned) == sum(daily[1:days[end]])

    ## The same gap on the live merged manifest. SitRep 043 (26 June) was not
    ## published, so the histories step 25 June -> 27 June and 26 June is a
    ## latent grid day with no observation mapped to it.
    obs = load_observations()
    ## Only require the manifest to reach past the 26 June hole so the indices
    ## below are on-grid; the exact cut-off advances with each data update, so
    ## it is not pinned here.
    @test obs.cutoff >= Date("2026-06-27")
    ## Grid index of a calendar date: the cut-off is the last grid day, and a
    ## day `k` before it sits at index `n - k` (the latent grid is unbroken).
    idx_of(d) = obs.n - value(obs.cutoff - d)
    i26 = idx_of(Date("2026-06-26"))
    i27 = idx_of(Date("2026-06-27"))
    i25 = idx_of(Date("2026-06-25"))
    @test i27 - i25 == 2               # a genuine two-day jump over the hole
    ## No per-vintage history indexes the missing 26 June grid day.
    for nm in (:confirmed_history, :confirmed_deaths_history,
        :isolation_history, :treatment_admissions_history,
        :suspected_daily_history, :recovered_history)
        @test i26 ∉ getfield(obs, nm).days
    end
end

@testitem "bed_capacity_walk: positive capacity path over the grid" tags=[
    :slow
] begin
    using Turing: sample, Prior
    import FlexiChains
    using BVDOutbreakSize: bed_capacity_walk_model

    ## The walk returns a positive bed-capacity path; with a tight innovation
    ## SD it stays a gentle drift around the baseline rather than blowing up.
    chn = sample(bed_capacity_walk_model(30), Prior(), 100;
        chain_type = FlexiChains.VNChain, progress = false)
    ks = string.(collect(keys(chn)))
    @test any(k -> occursin("C0", k), ks)
    C0 = vec(Array(chn[:C0]))
    @test all(C0 .> 0)
end

@testitem "isolation occupancy: conditioned fit stays positive" tags=[
    :slow
] begin
    using Turing: sample, Prior
    import FlexiChains
    using BVDOutbreakSize: treatment_only_model

    ## A daily occupancy stock on later days, supplied as observed counts.
    isolation_history = (; days = [28, 29, 30, 31, 32, 33],
        counts = [206, 233, 258, 267, 283, 309])
    chn = sample(
        treatment_only_model(33; isolation_history),
        Prior(), 100;
        chain_type = FlexiChains.VNChain, progress = false
    )
    C_T = vec(Array(chn[:C_T]))
    @test length(C_T) == 100
    @test all(isfinite, C_T)
    @test all(C_T .> 0)
end

@testitem "isolation occupancy: predictive path samples the counts" tags=[
    :slow
] begin
    using Turing: sample, Prior
    import FlexiChains
    using BVDOutbreakSize: treatment_only_model

    ## Days but no counts: the occupancy is a predictive generator, so its
    ## per-day counts are sampled under the `treatment_state.isolation`
    ## submodel rather than conditioned.
    isolation_history = (; days = [28, 29, 30, 31, 32, 33], counts = Int[])
    chn = sample(
        treatment_only_model(33; isolation_history),
        Prior(), 50;
        chain_type = FlexiChains.VNChain, progress = false
    )
    ks = string.(collect(keys(chn)))
    @test any(k -> occursin("isolation", k), ks)
end

@testitem "isolation occupancy: empty history is a no-op" tags=[:slow] begin
    using Turing: sample, Prior
    import FlexiChains
    using BVDOutbreakSize: treatment_only_model

    ## With no isolation history the occupancy submodel scores nothing and
    ## adds no sampled occupancy keys.
    chn = sample(
        treatment_only_model(33),
        Prior(), 50;
        chain_type = FlexiChains.VNChain, progress = false
    )
    ks = string.(collect(keys(chn)))
    @test !any(k -> occursin("isolation.increments", k), ks)
    C_T = vec(Array(chn[:C_T]))
    @test all(isfinite, C_T)
    @test all(C_T .> 0)
end

@testitem "accumulate_occupancy: balance closes and stays non-negative" begin
    using BVDOutbreakSize: accumulate_occupancy, convolve_delay
    death_pmf = [0.0, 0.1, 0.2, 0.3, 0.2, 0.2]
    recover_pmf = vcat(zeros(6), [0.2, 0.3, 0.3, 0.2])
    ruleout_pmf = [0.0, 0.3, 0.4, 0.3]
    CFR_iso = 0.35

    ## A constant admission inflow with simple discharge series builds a
    ## non-negative running balance; the confirmed sub-stock is a subset of
    ## the occupied BVD stock and the suspect sub-stock is the complement, so
    ## the two sum to the total demand on every day.
    n = 30
    A_bvd = fill(4.0, n)
    A_bg = fill(6.0, n)
    deaths = convolve_delay(CFR_iso .* A_bvd, death_pmf)
    recover = convolve_delay((1 - CFR_iso) .* A_bvd, recover_pmf)
    ruleout = convolve_delay(A_bg, ruleout_pmf)
    conf_hazard = fill(0.55, n)
    acc = accumulate_occupancy(A_bvd, A_bg, deaths, recover, ruleout, 0.01,
        conf_hazard)
    @test all(acc.demand .>= -1e-9)
    @test all(acc.O_conf .>= -1e-9)
    @test all(acc.O_susp .>= -1e-9)
    ## Confirmed is a subset of the occupied BVD stock.
    @test all(acc.O_conf .<= acc.O_bvd .+ 1e-9)
    ## Suspect + confirmed = total demand.
    @test all(abs.(acc.O_conf .+ acc.O_susp .- acc.demand) .< 1e-6)
    ## A zero confirmation hazard leaves the confirmed sub-stock empty.
    acc0 = accumulate_occupancy(A_bvd, A_bg, deaths, recover, ruleout, 0.01,
        zeros(n))
    @test all(acc0.O_conf .== 0)
    @test all(acc0.O_susp .≈ acc0.demand)

    ## Declining-occupancy regime: a rising-then-falling admission series so the
    ## stock peaks and then drains, with a non-trivial abscond rate. A
    ## total-only abscond accounting can let `O_bvd` exceed `D` and silently
    ## lose mass; with the two sub-stocks each draining their own share of
    ## absconds the invariant `O_conf + O_susp == D` must hold across the
    ## whole grid, including the declining tail, with both sub-stocks
    ## non-negative and the confirmed subset never exceeding the demand.
    m = 60
    t = 1:m
    ramp = @. 30.0 * exp(-((t - 18.0)^2) / (2 * 8.0^2))  # rise then fall
    A_bvd_d = 0.6 .* ramp
    A_bg_d = 0.4 .* ramp
    deaths_d = convolve_delay(CFR_iso .* A_bvd_d, death_pmf)
    recover_d = convolve_delay((1 - CFR_iso) .* A_bvd_d, recover_pmf)
    ruleout_d = convolve_delay(A_bg_d, ruleout_pmf)
    conf_hazard_d = fill(0.4, m)
    acc_d = accumulate_occupancy(
        A_bvd_d, A_bg_d, deaths_d, recover_d, ruleout_d, 0.02, conf_hazard_d)
    @test all(acc_d.demand .>= -1e-9)
    @test all(acc_d.O_conf .>= -1e-9)
    @test all(acc_d.O_susp .>= -1e-9)
    @test all(acc_d.O_conf .<= acc_d.O_bvd .+ 1e-9)
    @test all(acc_d.O_bvd .<= acc_d.demand .+ 1e-9)
    @test all(abs.(acc_d.O_conf .+ acc_d.O_susp .- acc_d.demand) .< 1e-6)
    ## The series genuinely declines, so the tail exercises the regime where a
    ## total-only accounting would lose mass.
    peak = argmax(acc_d.demand)
    @test peak < m
    @test acc_d.demand[end] < acc_d.demand[peak]
end

@testitem "clinical_stay_survival: complement matches stock balance" begin
    using BVDOutbreakSize: clinical_stay_survival, accumulate_occupancy,
                           convolve_delay
    ## The clinical-stay survival is the complement of the mixed discharge CDF,
    ## `S_clin(d) = 1 − Σ_{j≤d}(CFR·death_pmf + (1−CFR)·recover_pmf)`.
    ## With no same-day discharge it starts at 1, non-increasing toward 0.
    death_pmf = [0.0, 0.30, 0.35, 0.20, 0.10, 0.05]
    recover_pmf = vcat(zeros(8), [0.15, 0.25, 0.30, 0.20, 0.10])
    CFR = 0.45
    S = clinical_stay_survival(death_pmf, recover_pmf, CFR)
    @test S[1] ≈ 1.0
    @test all(diff(S) .<= 1e-12)             # non-increasing
    @test all(S .>= -1e-12)                  # non-negative
    @test S[end] ≈ 0.0 atol = 1e-9           # everyone discharged by the tail
    @test length(S) == max(length(death_pmf), length(recover_pmf))

    ## Cohort identity: the abscond-free occupied BVD stock the running balance
    ## builds equals the survival-weighted admission convolution with exactly
    ## this `S_clin`, so the two-clock weight is the stock's own per-cohort
    ## weight (not the inclusive `P(stay ≥ d)` tail).
    n = 40
    t = 1:n
    A_bvd = 0.6 .* (@. 25.0 * exp(-((t - 16.0)^2) / (2 * 7.0^2)))
    A_bg = 0.4 .* (@. 25.0 * exp(-((t - 16.0)^2) / (2 * 7.0^2)))
    deaths = convolve_delay(CFR .* A_bvd, death_pmf)
    recover = convolve_delay((1 - CFR) .* A_bvd, recover_pmf)
    ruleout = convolve_delay(A_bg, [0.0, 0.3, 0.4, 0.3])
    ## No absconds → the BVD stock is exactly the cohort survival
    ## reconstruction.
    acc0 = accumulate_occupancy(A_bvd, A_bg, deaths, recover, ruleout, 0.0,
        fill(0.3, n))
    @test all(abs.(acc0.O_bvd .- convolve_delay(A_bvd, S)) .< 1e-8)
end

@testitem "two_clock_confirmed: cohort split is a subset of the BVD stock" begin
    using BVDOutbreakSize: two_clock_confirmed, clinical_stay_survival,
                           accumulate_occupancy, convolve_delay
    death_pmf = [0.0, 0.30, 0.35, 0.20, 0.10, 0.05]      # fast death
    recover_pmf = vcat(zeros(8), [0.15, 0.25, 0.30, 0.20, 0.10])  # slow recover
    CFR = 0.45
    S = clinical_stay_survival(death_pmf, recover_pmf, CFR)

    n = 60
    t = 1:n
    A_bvd = 0.6 .* (@. 25.0 * exp(-((t - 20.0)^2) / (2 * 9.0^2)))
    A_bg = 0.4 .* (@. 25.0 * exp(-((t - 20.0)^2) / (2 * 9.0^2)))
    deaths = convolve_delay(CFR .* A_bvd, death_pmf)
    recover = convolve_delay((1 - CFR) .* A_bvd, recover_pmf)
    ruleout = convolve_delay(A_bg, [0.0, 0.3, 0.4, 0.3])
    conf_hazard = [d < 8 ? 0.05 : 0.30 for d in 1:n]

    ## Zero abscond: the raw cohort confirmed stock is a strict subset of the
    ## occupied BVD stock everywhere, no clamp needed. The two clocks are
    ## both referenced to admission, so the confirmed-and-present cohort can
    ## never exceed the present cohort.
    acc0 = accumulate_occupancy(A_bvd, A_bg, deaths, recover, ruleout, 0.0,
        conf_hazard)
    raw0 = two_clock_confirmed(A_bvd, conf_hazard, S)
    @test all(raw0 .<= acc0.O_bvd .+ 1e-9)
    @test all(raw0 .>= -1e-9)

    ## A zero confirmation hazard leaves the confirmed sub-stock empty.
    @test all(two_clock_confirmed(A_bvd, zeros(n), S) .== 0)

    ## A just-admitted cohort is never already confirmed: with a single
    ## admission on day 1 only, the confirmed stock on day 1 is exactly 0 (the
    ## hazard applies from the next day), and it then rises as the cohort
    ## confirms while it survives.
    A1 = vcat([10.0], zeros(n - 1))
    o1 = two_clock_confirmed(A1, fill(0.4, n), S)
    @test o1[1] == 0
    @test o1[2] > 0

    ## Monotone in the hazard: a uniformly higher hazard confirms more of the
    ## surviving cohort, so the confirmed stock is pointwise at least as large.
    lo = two_clock_confirmed(A_bvd, fill(0.1, n), S)
    hi = two_clock_confirmed(A_bvd, fill(0.5, n), S)
    @test all(hi .>= lo .- 1e-9)

    ## Fast-death tail: against the proportional-share split, the two-clock
    ## confirmed share is lower early because cases that die before confirming
    ## are excluded from the confirmed pool (the comparator's whole point).
    ## The proportional split drains the confirmed pool at the pool-average
    ## discharge rate, over-attributing confirmed deaths in the early
    ## fast-death window.
    O_conf_prop = acc0.O_conf
    O_bvd = acc0.O_bvd
    early = 3:12
    share_prop = [O_bvd[i] > 0 ? O_conf_prop[i] / O_bvd[i] : 0.0 for i in early]
    share_tc = [O_bvd[i] > 0 ? raw0[i] / O_bvd[i] : 0.0 for i in early]
    @test all(share_tc .<= share_prop .+ 1e-9)
    @test sum(share_tc) < sum(share_prop)
end

@testitem "isolation abscond: driven by the two-clock suspect stock" begin
    using BVDOutbreakSize: treatment_flow_model
    using Turing: returned
    using Random: MersenneTwister
    ## The scored abscond flow is `κ · O_susp(t-1)` off the two-clock suspect
    ## stock `O_susp = demand − two-clock O_conf`, not the proportional-share
    ## suspect `accumulate_occupancy` carries. With no manual occupancy break
    ## the returned `suspect_incare` census is that two-clock `O_susp` (no
    ## offset), so the abscond flow must be a one-day-lagged multiple of it. A
    ## non-zero confirmation hazard makes the two-clock and proportional splits
    ## differ, so this pins the wiring onto the two-clock stock.
    n = 40
    t = 1:n
    bvd_reports_daily = @. 20.0 * exp(-((t - 18.0)^2) / (2 * 6.0^2))
    bg_daily = @. 12.0 * exp(-((t - 18.0)^2) / (2 * 6.0^2))
    conf_hazard = fill(0.25, n)
    isolation_history = (; days = collect(10:2:38),
        counts = fill(50, length(10:2:38)))
    model = treatment_flow_model(isolation_history, bvd_reports_daily, bg_daily,
        0.6, 0.4; conf_hazard_daily = conf_hazard)
    st = returned(model, rand(MersenneTwister(1), model))
    ab = st.abscond_daily
    κ = st.abscond_frac
    ## The suspect stock the abscond flow multiplies is the demand remainder
    ## after the two-clock confirmed sub-stock, `O_susp = D − O_conf`.
    ## Reconstruct it from the returned demand and confirmed census, not from
    ## `suspect_incare`, so the pin still holds under a non-zero occupancy-break
    ## offset (the census carries the offset; the abscond stock does not).
    O_susp = st.demand .- st.confirmed_incare
    ## Day 1 has no prior stock, so no abscond outflow.
    @test ab[1] == 0
    ## Every later day is exactly the two-clock suspect lagged one day, scaled
    ## by the abscond fraction — the definitional link to the two-clock split.
    @test all(abs(ab[i] - κ * O_susp[i - 1]) < 1e-9 for i in 2:n)
    ## With no occupancy break the returned suspect census is that same stock.
    @test st.suspect_incare ≈ O_susp
    ## The flow genuinely responds to the two-clock confirmed pool: with a
    ## non-zero hazard the confirmed sub-stock is populated, so the suspect
    ## stock sits strictly below the total demand and the abscond outflow is
    ## below `κ · demand(t-1)` — it consumes the two-clock, not the total,
    ## occupancy.
    @test any(>(0), st.confirmed_incare)
    @test any(O_susp[i] < st.demand[i] - 1e-9 for i in 1:n)
end

@testitem "admission_headroom: bound above obs, never on the boundary" begin
    using BVDOutbreakSize: admission_headroom
    ## Capacity 400 with previous-day occupancy 260 leaves 140 free beds; a
    ## slack admission count is censored at the headroom, a saturated count is
    ## kept strictly inside the support (obs + 0.5), so it never sits on the
    ## non-differentiable NegBin-CDF boundary.
    capacity_history = (; days = [9, 10], counts = [400, 410])
    isolation_history = (; days = [9, 10], counts = [260, 262])
    adm_days = [10, 11]
    adm_obs = [50, 300]
    head = admission_headroom(adm_days, adm_obs, capacity_history,
        isolation_history)
    @test head[1] ≈ 410 - 260         # day 10: capacity 410 less prev occ 260
    @test head[2] > adm_obs[2]        # day 11: saturated, strictly above obs
    @test all(head .> adm_obs .- 1e-9)
    ## No capacity record gives a large no-op headroom.
    nocap = admission_headroom([10], [50], (; days = Int[], counts = Int[]),
        isolation_history)
    @test only(nocap) > 1.0e5
end

@testitem "occupancy split: sub-stock parameters sampled, stays positive" tags=[
    :slow
] begin
    using Turing: sample, Prior
    import FlexiChains
    using BVDOutbreakSize: treatment_only_model

    ## Occupancy on days 28-33 with a published split on the last three days.
    ## On split days the two sub-stock censuses are scored instead of the total
    ## (a per-day total-or-split switch); the abscond fraction and the manual
    ## occupancy break-day step are sampled and stay in range. Day 32 is passed
    ## as a manual `occupancy_break_days`, so a break step is fitted there.
    ## The lab / confirmed stream is conditioned so the in-care confirmation
    ## hazard is non-zero and the split is identified (the coherent config).
    isolation_history = (; days = [28, 29, 30, 31, 32, 33],
        counts = [206, 233, 258, 267, 283, 309])
    confirmed_incare = (; days = [31, 32, 33], counts = [120, 130, 140])
    suspect_incare = (; days = [31, 32, 33], counts = [147, 153, 169])
    confirmed_history = (; days = [28, 30, 32], counts = [40, 70, 110])
    lab_history = (; days = [28, 30, 32], counts = [120, 200, 320])
    chn = sample(
        treatment_only_model(33; isolation_history,
            confirmed_history, confirmed_cases = 110, lab_history,
            treatment_confirmed_incare_history = confirmed_incare,
            treatment_suspect_incare_history = suspect_incare,
            occupancy_break_days = [32]),
        Prior(), 60;
        chain_type = FlexiChains.VNChain, progress = false)
    ks = string.(collect(keys(chn)))
    ## The manual occupancy break step (non-centred) is sampled, and the
    ## cut-off cumulative offset deterministic is exposed.
    @test any(k -> occursin("occupancy_step", k), ks)
    @test any(k -> occursin("occupancy_break", k), ks)
    @test !any(k -> occursin("total_break", k), ks)
    @test any(k -> occursin("abscond_frac", k), ks)
    ## With a non-zero borrowed hazard the in-care confirmation-rate modifier
    ## is sampled (its log parameter) and exposed as a finite positive
    ## deterministic ρ, the free lever the confirmed/suspected-in-care split
    ## identifies.
    @test any(k -> occursin("incare_confirm_log", k), ks)
    ρ_key = first(k for k in keys(chn)
    if occursin("incare_confirm_modifier", string(k)))
    ρ = vec(Array(chn[ρ_key]))
    @test all(isfinite, ρ)
    @test all(ρ .> 0)
    C_T = vec(Array(chn[:C_T]))
    @test all(isfinite, C_T)
    @test all(C_T .> 0)
    ## The abscond fraction stays a small non-negative loss-to-follow-up rate
    ## under the prior. Index by the chain key object itself (FlexiChains
    ## resolves the submodel-prefixed varname), not a reconstructed Symbol.
    ab_key = first(k for k in keys(chn) if occursin("abscond_frac", string(k)))
    ab = vec(Array(chn[ab_key]))
    @test all(ab .>= 0)
end

@testitem "occupancy split: predictive path samples sub-stock censuses" tags=[
    :slow
] begin
    using Turing: sample, Prior
    import FlexiChains
    using BVDOutbreakSize: treatment_only_model

    ## Days but no counts on the split histories: the two sub-stock censuses
    ## are predictive generators, so their per-day increments are sampled
    ## under the `confirmed_incare_obs` / `suspect_incare_obs` submodels. The
    ## split is only identified when the in-care confirmation overlay is
    ## non-zero, so the lab/confirmed stream is conditioned to supply the
    ## hazard (the coherent config); without it the split would be unscored
    ## (the guarded path is covered separately below).
    isolation_history = (; days = [28, 29, 30, 31, 32, 33],
        counts = [206, 233, 258, 267, 283, 309])
    confirmed_history = (; days = [28, 30, 32], counts = [40, 70, 110])
    lab_history = (; days = [28, 30, 32], counts = [120, 200, 320])
    chn = sample(
        treatment_only_model(33; isolation_history,
            confirmed_history, confirmed_cases = 110, lab_history,
            treatment_confirmed_incare_history = (; days = [31, 32, 33],
                counts = Int[]),
            treatment_suspect_incare_history = (; days = [31, 32, 33],
                counts = Int[])),
        Prior(), 40;
        chain_type = FlexiChains.VNChain, progress = false)
    ks = string.(collect(keys(chn)))
    @test any(k -> occursin("confirmed_incare_obs.increments", k), ks)
    @test any(k -> occursin("suspect_incare_obs.increments", k), ks)
end

@testitem "occupancy split: empty split history is a no-op" tags=[:slow] begin
    using Turing: sample, Prior
    import FlexiChains
    using BVDOutbreakSize: treatment_only_model

    ## With no split history the two sub-stock census likelihoods score
    ## nothing (no days → no scored or sampled increments), so no
    ## split-observation keys appear, while the total-occupancy backbone
    ## still runs.
    isolation_history = (; days = [28, 29, 30, 31, 32, 33],
        counts = [206, 233, 258, 267, 283, 309])
    chn = sample(
        treatment_only_model(33; isolation_history),
        Prior(), 50;
        chain_type = FlexiChains.VNChain, progress = false)
    ks = string.(collect(keys(chn)))
    @test !any(k -> occursin("confirmed_incare_obs.increments", k), ks)
    @test !any(k -> occursin("suspect_incare_obs.increments", k), ks)
    C_T = vec(Array(chn[:C_T]))
    @test all(isfinite, C_T)
    @test all(C_T .> 0)
end

@testitem "occupancy split: no lab data leaves it unscored and finite" tags=[
    :slow
] begin
    using Turing: sample, Prior
    import FlexiChains
    using BVDOutbreakSize: treatment_only_model

    ## A published split census but no lab/confirmed data: the borrowed
    ## in-care confirmation hazard `τ_test · p_pos` is then structurally
    ## zero, so the modelled confirmed sub-stock is empty by construction.
    ## The split-census likelihood must no-op (the split is identified only
    ## in the joint where the lab pipeline exists) rather than score the
    ## observed confirmed-in-care against a zero sub-stock and blow up. The
    ## fit must run and stay finite, with no split-observation keys scored,
    ## and the split days fall back to the total-occupancy backbone.
    isolation_history = (; days = [28, 29, 30, 31, 32, 33],
        counts = [206, 233, 258, 267, 283, 309])
    confirmed_incare = (; days = [31, 32, 33], counts = [120, 130, 140])
    suspect_incare = (; days = [31, 32, 33], counts = [147, 153, 169])
    chn = sample(
        treatment_only_model(33; isolation_history,
            treatment_confirmed_incare_history = confirmed_incare,
            treatment_suspect_incare_history = suspect_incare),
        Prior(), 60;
        chain_type = FlexiChains.VNChain, progress = false)
    ks = string.(collect(keys(chn)))
    ## The split census is unscored (no sampled or scored sub-stock increments).
    @test !any(k -> occursin("confirmed_incare_obs.increments", k), ks)
    @test !any(k -> occursin("suspect_incare_obs.increments", k), ks)
    ## With a structurally-zero borrowed hazard the in-care confirmation-rate
    ## modifier is unidentified, so it is not sampled (ρ = 1, a no-op) rather
    ## than adding a vestigial dimension that only follows its prior.
    @test !any(k -> occursin("incare_confirm_log", k), ks)
    ## The total-occupancy backbone still runs and stays finite (no 1e45 blow-up
    ## from scoring against a structurally-zero confirmed sub-stock).
    C_T = vec(Array(chn[:C_T]))
    @test all(isfinite, C_T)
    @test all(C_T .> 0)
    ## The cut-off bed demand stays a finite, bounded stock (an incoherent
    ## config would otherwise diverge to ~1e45). Index by the chain key
    ## object so FlexiChains resolves the submodel-prefixed varname.
    dem_key = first(k for k in keys(chn)
    if occursin("expected_bed_demand", string(k)))
    dem = vec(Array(chn[dem_key]))
    @test all(isfinite, dem)
    @test all(dem .< 1.0e6)
end

@testitem "occupancy split: occupancy outlasts the census (no-flow days)" begin
    using Turing: DynamicPPL
    using LogDensityProblems: logdensity
    using Random: seed!
    using BVDOutbreakSize: treatment_only_model

    ## A situation report that advances the occupancy/flow/case grid but leaves
    ## the suspect/confirmed in-care census frozen a day or two earlier (the
    ## 24-25 June gap: SitRep 042 extended occupancy/cases/lab but not the
    ## Tableau 6 split). The occupancy backbone and the daily flows run to the
    ## grid end (day 35) while the split census and the flow histories stop two
    ## days short (day 33). The per-day flow / split likelihoods must score only
    ## the days each history actually reports and skip the two trailing
    ## occupancy-only days cleanly: no day index reaches past the end of a
    ## shorter history, the split day-set is the census days alone (so the
    ## trailing occupancy days stay on the total-occupancy backbone), and the
    ## conditioned log density is finite.
    n = 35
    isolation_history = (; days = collect(28:35),
        counts = [206, 233, 258, 267, 283, 309, 301, 297])
    ## Flows reported only to day 33 (two days short of the occupancy grid).
    admissions = (; days = [31, 32, 33], counts = [60, 55, 61])
    incare_deaths = (; days = [31, 32, 33], counts = [9, 14, 9])
    ruleouts = (; days = [31, 32, 33], counts = [34, 26, 36])
    absconded = (; days = [31, 32, 33], counts = [1, 6, 1])
    ## Census also reported only to day 33.
    confirmed_incare = (; days = [31, 32, 33], counts = [120, 130, 140])
    suspect_incare = (; days = [31, 32, 33], counts = [147, 153, 169])
    ## A lab / confirmed stream so the in-care confirmation hazard is non-zero
    ## and the split census is identified (scored, not the guarded no-op path).
    confirmed_history = (; days = [28, 30, 32], counts = [40, 70, 110])
    lab_history = (; days = [28, 30, 32], counts = [120, 200, 320])
    ## A manual occupancy break day (33), two days short of the grid end. The
    ## per-day break offset must stay finite across the trailing no-flow days,
    ## not index past the grid.

    model = treatment_only_model(n; isolation_history,
        treatment_admissions_history = admissions,
        treatment_deaths_history = incare_deaths,
        treatment_ruleout_history = ruleouts,
        treatment_absconded_history = absconded,
        treatment_confirmed_incare_history = confirmed_incare,
        treatment_suspect_incare_history = suspect_incare,
        occupancy_break_days = [33],
        confirmed_history, confirmed_cases = 110, lab_history)

    ## A prior draw plus a conditioned log-density evaluation exercises every
    ## per-day flow / split likelihood against the shorter histories. The
    ## occupancy runs to day 35 with the census ending at day 33, so this is the
    ## exact reported-occupancy / no-flow gap; any history indexed past its end
    ## (or a mis-scored trailing day) would throw a bounds error or return a
    ## non-finite value here.
    seed!(424242)
    vi = DynamicPPL.VarInfo(model)
    ldf = DynamicPPL.LogDensityFunction(model, DynamicPPL.getlogjoint, vi)
    lp = logdensity(ldf, collect(vi[:]))
    @test isfinite(lp)

    ## Several further prior draws: none trips a bounds / dimension error from
    ## the trailing occupancy-only days.
    for s in 1:5
        seed!(1000 + s)
        vi_s = DynamicPPL.VarInfo(model)
        @test isfinite(logdensity(ldf, collect(vi_s[:])))
    end
end

@testitem "cumulative_occupancy_offset: accumulates from the break day" begin
    using BVDOutbreakSize: cumulative_occupancy_offset
    iso_days = [10, 11, 12, 13, 14]
    ## Two manual break days: +5 from day 12, −3 from day 14. Δ is zero before
    ## the first break, then cumulates the steps at or before each iso day.
    Δ = cumulative_occupancy_offset(iso_days, [12, 14], [5.0, -3.0])
    @test Δ == [0.0, 0.0, 5.0, 5.0, 2.0]
    ## No break days is a no-op (all zeros).
    @test cumulative_occupancy_offset(iso_days, Int[], Float64[]) ==
          zeros(5)
end

@testitem "isolation occupancy: no break days sample no offset step" tags=[
    :slow
] begin
    using Turing: sample, Prior
    import FlexiChains
    using BVDOutbreakSize: treatment_only_model

    isolation_history = (; days = [28, 29, 30, 31, 32, 33],
        counts = [206, 233, 258, 267, 283, 309])
    chn = sample(
        treatment_only_model(33; isolation_history),
        Prior(), 50;
        chain_type = FlexiChains.VNChain, progress = false
    )
    ks = collect(keys(chn))
    ## The opt-in offset is off by default: no sampled step, offset pinned zero.
    @test !any(k -> occursin("occupancy_step", string(k)), ks)
    occ = only(filter(k -> occursin("occupancy_break", string(k)), ks))
    brk = vec(Array(chn[occ]))
    @test all(==(0), brk)
end

@testitem "isolation occupancy: a manual break day fits an offset step" tags=[
    :slow
] begin
    using Turing: sample, Prior
    import FlexiChains
    using BVDOutbreakSize: treatment_only_model

    isolation_history = (; days = [28, 29, 30, 31, 32, 33],
        counts = [206, 233, 258, 267, 283, 309])
    ## Opt in to a single break on day 31; a step is sampled and the cut-off
    ## cumulative offset is finite (non-zero prior draws).
    chn = sample(
        treatment_only_model(33; isolation_history,
            occupancy_break_days = [31]),
        Prior(), 100;
        chain_type = FlexiChains.VNChain, progress = false
    )
    ks = collect(keys(chn))
    @test any(k -> occursin("occupancy_step", string(k)), ks)
    occ = only(filter(k -> occursin("occupancy_break", string(k)), ks))
    brk = vec(Array(chn[occ]))
    @test length(brk) == 100
    @test all(isfinite, brk)
    @test any(!=(0), brk)
end

@testitem "isolation occupancy: joint prior runs with the live data" tags=[
    :slow
] begin
    using Turing: sample, Prior
    import FlexiChains
    using BVDOutbreakSize: load_observations, bvd_joint, genetic_seeding_model

    obs = load_observations()
    @test !isempty(obs.isolation_history.counts)
    breakpoint = obs.n - obs.who_first_sitrep_days
    m = bvd_joint(obs.n, obs.exported_cases, obs.total_deaths,
        obs.reported_cases, obs.exports_deaths, obs.confirmed_cases,
        obs.tests_analysed;
        confirmed_deaths = obs.confirmed_deaths,
        deaths_history = obs.deaths_history,
        reported_history = obs.reported_history,
        confirmed_history = obs.confirmed_history,
        confirmed_deaths_history = obs.confirmed_deaths_history,
        lab_history = obs.lab_history,
        lab_daily_history = obs.lab_daily_history,
        suspected_daily_history = obs.suspected_daily_history,
        isolation_history = obs.isolation_history,
        bed_capacity_history = obs.bed_capacity_history,
        treatment_confirmed_incare_history =
        obs.treatment_confirmed_incare_history,
        treatment_suspect_incare_history =
        obs.treatment_suspect_incare_history,
        export_case_days = obs.export_case_days,
        export_death_days = obs.export_death_days,
        breakpoint = breakpoint,
        genetic = genetic_seeding_model,
        tmrca_days = obs.tmrca_days)
    chn = sample(m, Prior(), 30;
        chain_type = FlexiChains.VNChain, progress = false)
    C_T = vec(Array(chn[:C_T]))
    iso = vec(Array(chn[:expected_isolation_T]))
    dem = vec(Array(chn[:expected_bed_demand_T]))
    cap = vec(Array(chn[:bed_capacity]))
    @test length(C_T) == 30
    @test all(isfinite, C_T)
    @test all(C_T .> 0)
    @test all(isfinite, iso)
    @test all(iso .> 0)
    ## Occupancy is `min(demand, C)`, so it never exceeds the latent demand,
    ## and the supply-limited occupancy never exceeds the bed capacity.
    @test all(iso .<= dem .+ 1e-6)
    @test all(iso .<= cap .+ 1e-6)
    ## The severity skew is non-negative and admits BVD suspects at least as
    ## readily as the base (non-BVD rule-out) rate.
    skew = vec(Array(chn[:isolation_severity]))
    p_base = vec(Array(chn[:isolation_admission]))
    p_bvd = vec(Array(chn[:isolation_bvd_admission]))
    @test all(skew .>= -1e-9)
    @test all(p_bvd .>= p_base .- 1e-9)
    @test all(0 .<= p_bvd .<= 1)
end
