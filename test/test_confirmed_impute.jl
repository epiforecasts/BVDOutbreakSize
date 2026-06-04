## Tests for the imputed analysed-denominator extension. Confirmed
## vintages whose national analysed count is missing (the early 18-22 May
## and late 29-31 May lab windows) get a TIGHT log-random-walk denominator
## anchored to the observed 23-28 May increments, so they can be fitted
## without a free per-vintage denominator funnelling against positivity.

@testsnippet ImputeFixtures begin
    using BVDOutbreakSize: bvd_joint, load_observations, nuts_sample,
                           analysed_impute_model, analysed_curve_model,
                           confirmed_q_re_model, fit_diagnostics
    using Turing: sample, Prior
    import FlexiChains

    function _inc(values)
        out = similar(values, Int)
        prev = 0
        for i in eachindex(values)
            out[i] = values[i] - prev
            prev = values[i]
        end
        return out
    end

    ## Build confirmed vectors with the early 18-22 May vintages carrying a
    ## missing analysed denominator and the 23-28 May vintages observed.
    function early_extension(obs)
        ch = obs.confirmed_case_history
        sa = obs.tests_analysed_history
        sr = obs.tests_received_history
        keep = [i == 1 || sa.values[i] > sa.values[i - 1]
                for i in eachindex(sa.values)]
        aoff = collect(sa.offsets)[keep]
        analysed_base = _inc(sa.values[keep])
        ridx = [findfirst(==(off), sr.offsets) for off in aoff]
        recv_base = _inc([sr.values[i] for i in ridx])
        offs, ccum = Int[], Int[]
        a_vals = Union{Missing, Int}[]
        r_vals = Union{Missing, Int}[]
        for off in [10, 9, 8, 7, 6]
            j = findfirst(==(off), ch.offsets)
            push!(offs, off)
            push!(ccum, ch.values[j])
            push!(a_vals, missing)
            push!(r_vals, missing)
        end
        for (k, off) in enumerate(aoff)
            j = findfirst(==(off), ch.offsets)
            push!(offs, off)
            push!(ccum, ch.values[j])
            push!(a_vals, analysed_base[k])
            push!(r_vals, recv_base[k])
        end
        return (offsets = offs, confirmed = _inc(ccum),
            analysed = a_vals, received = r_vals)
    end

    ## Full-window confirmed vectors: the early 18-22 May AND late
    ## 29-31 May vintages carry a missing analysed denominator, with the
    ## observed 23-28 May block in between carrying the published
    ## cumulative analysed anchors. Used by the low-DOF curve method.
    function full_window(obs)
        ch = obs.confirmed_case_history
        sa = obs.tests_analysed_history
        late_off = [-1, -2, -3]
        late_cum = [263, 282, 321]
        offs = vcat(collect(ch.offsets), late_off)
        ccum = vcat(collect(ch.values), late_cum)
        confirmed = _inc(ccum)
        analysed_cum = Vector{Union{Missing, Float64}}(
            missing, length(offs))
        analysed = Vector{Union{Missing, Int}}(missing, length(offs))
        prev = 0.0
        for (k, off) in enumerate(sa.offsets)
            j = findfirst(==(off), offs)
            j === nothing && continue
            analysed_cum[j] = Float64(sa.values[k])
            ## Floor the analysed increment at the confirmed increment so the
            ## 24-25 May lab stall (ΔA = 0 with ΔC > 0) stays a valid
            ## Binomial denominator.
            analysed[j] = max(Int(sa.values[k] - prev), confirmed[j])
            prev = Float64(sa.values[k])
        end
        return (offsets = offs, confirmed = confirmed,
            analysed = analysed, analysed_cum = analysed_cum)
    end
end

@testitem "analysed_impute_model: tight log-RW around the anchor" tags=[:slow] setup=[ImputeFixtures] begin
    chn=sample(analysed_impute_model(5, log(120.0)), Prior(), 500;
        chain_type = FlexiChains.VNChain, progress = false)
    σ=vec(Array(chn[:σ_A]))
    @test all(σ .>= 0)
    ## Default walk SD is tight (truncated Normal(0, 0.3)).
    @test 0.0 < sum(σ) / length(σ) < 0.5
end

@testitem "impute extension: prior draws finite and valid" tags=[:slow] setup=[ImputeFixtures] begin
    obs=load_observations()
    c=early_extension(obs)
    rep=obs.reported_case_history
    dh=obs.death_history
    m=bvd_joint(obs.exported_cases, _inc(dh.values), _inc(rep.values),
        obs.export_deaths_daily;
        reported_offsets = rep.offsets,
        death_offsets = dh.offsets,
        confirmed_cases = c.confirmed,
        confirmed_offsets = c.offsets,
        samples_analysed = c.analysed,
        samples_received = c.received,
        tests_analysed = obs.cumulative_tests_analysed,
        tests_offset = 0,
        first_export_detection_delta = obs.first_export_detection_delta,
        confirmed_q_random_effect = confirmed_q_re_model,
        confirmed_analysed_impute = analysed_impute_model)
    chn=sample(m, Prior(), 100;
        chain_type = FlexiChains.VNChain, progress = false)
    C=vec(Array(chn[:cumulative_cases]))
    @test all(isfinite, C)
    @test all(C .> 0)
end

@testitem "impute extension: early window enters the gradient fit" tags=[:slow] setup=[ImputeFixtures] begin
    ## The imputed denominator lets the early 18-22 May confirmed vintages
    ## enter the gradient-based fit with finite, positive C(T). This only
    ## checks the Binomial + log-RW path runs: a full fit does NOT converge
    ## (the imputed ΔA funnels against the q random effect, R-hat ≈ 2),
    ## which is why the extension is off by default. A few draws suffice.
    obs=load_observations()
    c=early_extension(obs)
    rep=obs.reported_case_history
    dh=obs.death_history
    chn=nuts_sample(
        bvd_joint(obs.exported_cases, _inc(dh.values), _inc(rep.values),
            obs.export_deaths_daily;
            reported_offsets = rep.offsets,
            death_offsets = dh.offsets,
            confirmed_cases = c.confirmed,
            confirmed_offsets = c.offsets,
            samples_analysed = c.analysed,
            samples_received = c.received,
            tests_analysed = obs.cumulative_tests_analysed,
            tests_offset = 0,
            first_export_detection_delta =
            obs.first_export_detection_delta,
            confirmed_q_random_effect = confirmed_q_re_model,
            confirmed_analysed_impute = analysed_impute_model);
        samples = 5, chains = 1, seed = 1, progress = false)
    C=vec(Array(chn[:cumulative_cases]))
    @test all(isfinite, C)
    @test all(C .> 0)
end

@testitem "curve method: all confirmed vintages enter the fit" tags=[:slow] setup=[ImputeFixtures] begin
    ## Low-DOF parametric extrapolation: the cumulative-analysed curve
    ## (two shared parameters) supplies the missing 18-22 May and
    ## 29-31 May denominators as deterministic differences, so ALL
    ## confirmed vintages enter the gradient fit with finite, positive
    ## C(T) — without the per-vintage funnel of the free impute latent.
    obs=load_observations()
    c=full_window(obs)
    rep=obs.reported_case_history
    dh=obs.death_history
    ## Every confirmed vintage is present; only the early/late denominators
    ## are missing (supplied by the curve).
    n_missing=count(ismissing, c.analysed)
    @test n_missing == 8
    @test length(c.confirmed) == length(c.offsets) == 14
    chn=nuts_sample(
        bvd_joint(obs.exported_cases, _inc(dh.values), _inc(rep.values),
            obs.export_deaths_daily;
            reported_offsets = rep.offsets,
            death_offsets = dh.offsets,
            confirmed_cases = c.confirmed,
            confirmed_offsets = c.offsets,
            samples_analysed = c.analysed,
            tests_analysed = obs.cumulative_tests_analysed,
            tests_offset = 0,
            first_export_detection_delta =
            obs.first_export_detection_delta,
            report_onset_offset = 8,
            confirmed_q_random_effect = confirmed_q_re_model,
            confirmed_analysed_curve = analysed_curve_model,
            confirmed_analysed_cum = c.analysed_cum,
            confirmed_selection_clock = :volume);
        samples = 5, chains = 1, seed = 1, progress = false)
    C=vec(Array(chn[:cumulative_cases]))
    @test all(isfinite, C)
    @test all(C .> 0)
end

@testitem "analysed_curve_model: OLS extrapolation is monotone" tags=[:slow] setup=[ImputeFixtures] begin
    using BVDOutbreakSize: analysed_curve_model
    ## Six observed cumulative anchors (23-28 May), missing early (5) and
    ## late (3). The deterministic OLS fit must return a monotone
    ## cumulative passing near the anchors and extrapolating both ways.
    t=Float64.(0:13)
    oc=Vector{Union{Missing, Float64}}(missing, 14)
    oc[6:11]=[211.0, 295.0, 295.0, 403.0, 648.0, 755.0]
    chn=sample(analysed_curve_model(collect(t), oc), Prior(), 1;
        chain_type = FlexiChains.VNChain, progress = false)
    ## Reconstruct deterministically (T-invariant). Anchors fitted by OLS.
    obs_t=t[6:11]
    obs_y=log.(oc[6:11])
    tref=sum(obs_t)/length(obs_t)
    xc=obs_t .- tref
    ybar=sum(obs_y)/length(obs_y)
    b=sum(xc .* (obs_y .- ybar))/sum(abs2, xc)
    a=ybar
    A=[exp(a+b*(ti-tref)) for ti in t]
    @test issorted(A)                      # monotone cumulative
    @test b > 0                            # increasing
    @test A[1] < 211                       # backward extrapolation < first
    @test A[end] > 755                     # forward extrapolation > last
end
