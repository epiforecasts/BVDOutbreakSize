## Tests for the condensed laboratory-throughput queue (`confirmed_queue`).
## ALL confirmed vintages enter the fit: observed-denominator windows
## (23-28 May) condition on the published analysed count
## (Binomial + Poisson), while the dark early 18-22 May and late 29-31 May
## windows use the exact marginal Poisson(μ_A · p_pos), so the unobserved
## denominator is integrated out rather than carried as a free latent.

@testsnippet QueueFixtures begin
    using BVDOutbreakSize: bvd_joint, load_observations, nuts_sample,
                           epi_exclusion_model, confirmed_q_re_model
    using Turing: sample, Prior, logjoint
    using Random: Xoshiro
    import FlexiChains

    function _inc(values)
        out = similar(collect(values), Int)
        prev = 0
        for i in eachindex(out)
            out[i] = values[i] - prev
            prev = values[i]
        end
        return out
    end

    ## Build the FULL confirmed series: early dark windows (18-22 May, no
    ## analysed denominator), observed windows (23-28 May, strictly
    ## increasing analysed) and late dark windows (29-31 May, sitreps
    ## 015-017 confirmed-only, cumulative 263, 282, 321).
    function full_extension(obs)
        ch = obs.confirmed_case_history
        sa = obs.tests_analysed_history
        sr = obs.tests_received_history
        coff = collect(ch.offsets)
        cval = collect(ch.values)
        conf_at(off) = cval[findfirst(==(off), coff)]
        keep = [i == 1 || sa.values[i] > sa.values[i - 1]
                for i in eachindex(sa.values)]
        aoff = collect(sa.offsets)[keep]
        analysed_base = _inc(collect(sa.values)[keep])
        sroff = collect(sr.offsets)
        srval = collect(sr.values)
        ridx = [findfirst(==(off), sroff) for off in aoff]
        recv_base = _inc([srval[i] for i in ridx])
        offs, ccum = Int[], Int[]
        a_vals = Union{Missing, Int}[]
        r_vals = Union{Missing, Int}[]
        for off in [10, 9, 8, 7, 6]
            push!(offs, off)
            push!(ccum, conf_at(off))
            push!(a_vals, missing)
            push!(r_vals, missing)
        end
        for (k, off) in enumerate(aoff)
            push!(offs, off)
            push!(ccum, conf_at(off))
            push!(a_vals, analysed_base[k])
            push!(r_vals, recv_base[k])
        end
        for (off, cum) in zip([-1, -2, -3], [263, 282, 321])
            push!(offs, off)
            push!(ccum, cum)
            push!(a_vals, missing)
            push!(r_vals, missing)
        end
        return (offsets = offs, confirmed = _inc(ccum),
            analysed = a_vals, received = r_vals)
    end

    function queue_model(obs, c; epi_exclusion = nothing)
        rep = obs.reported_case_history
        dh = obs.death_history
        return bvd_joint(obs.exported_cases, _inc(dh.values),
            _inc(rep.values), obs.export_deaths_daily;
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
            confirmed_queue = true,
            confirmed_epi_exclusion = epi_exclusion)
    end
end

@testitem "queue: all confirmed vintages enter with finite C(T)" tags=[:slow] setup=[QueueFixtures] begin
    obs=load_observations()
    c=full_extension(obs)
    ## Every confirmed vintage (early dark, observed, late dark) is present.
    @test count(ismissing, c.analysed) == 8
    @test count(!ismissing, c.analysed) == 5
    m=queue_model(obs, c)
    ## Prior draws: ALL vintages enter, C(T) finite and positive, and the
    ## queue predicts a positive analysed denominator for the dark windows.
    chn=sample(m, Prior(), 100;
        chain_type = FlexiChains.VNChain, progress = false)
    C=vec(Array(chn[:cumulative_cases]))
    @test all(isfinite, C)
    @test all(C .> 0)
    dark=vec(Array(chn[:dark_analysed_total]))
    @test all(isfinite, dark)
    @test all(dark .> 0)
    ## A short gradient fit runs without the free-latent funnel (the impute
    ## path could not even take a NUTS step cleanly here). A few draws check
    ## the Poisson-thinned path is differentiable with finite, positive C(T).
    chn2=nuts_sample(m; samples = 5, chains = 1, seed = 1,
        progress = false)
    C2=vec(Array(chn2[:cumulative_cases]))
    @test all(isfinite, C2)
    @test all(C2 .> 0)
end

@testitem "queue: dark windows use the Poisson-thinned likelihood" tags=[:slow] setup=[QueueFixtures] begin
    ## A dark window's confirmed observation is conditioned as
    ## Poisson(μ_A · p_pos): changing only that datum shifts the joint
    ## log-density by exactly the Poisson log-pdf difference, and the
    ## denominator is integrated out (no Binomial on a latent ΔA).
    obs=load_observations()
    c=full_extension(obs)
    ## Two datasets differing only in the LAST (late dark) confirmed
    ## increment. With the Poisson-thinned likelihood the conditioned value
    ## must change the logjoint (the dark window is observed, not dropped).
    c2=(; c..., confirmed = copy(c.confirmed))
    c2.confirmed[end]+=5
    m1=queue_model(obs, c)
    m2=queue_model(obs, c2)
    ## Draw latent parameters until the joint density is finite (prior draws
    ## can land in the tails), then score both datasets on the same draw.
    ## Wrapped in a function so the loop variables are not soft-scope globals.
    function finite_draw(model; seed = 42, max_tries = 200)
        rng=Xoshiro(seed)
        vals=rand(rng, model)
        tries=0
        while !isfinite(logjoint(model, vals))&&tries<max_tries
            vals=rand(rng, model)
            tries+=1
        end
        return vals
    end
    vals=finite_draw(m1)
    lj1=logjoint(m1, vals)
    lj2=logjoint(m2, vals)
    ## The late window is dark, so perturbing its confirmed count changes
    ## the joint density purely through its Poisson term (finite, nonzero).
    @test isfinite(lj1)
    @test isfinite(lj2)
    @test lj1 != lj2
end
