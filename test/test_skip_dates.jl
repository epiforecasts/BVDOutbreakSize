@testitem "a skipped vintage merges into the next window" begin
    using BVDOutbreakSize
    using BVDOutbreakSize: load_observations
    using Dates: Date
    using TOML

    ## The fixture round-trips the manifest through the TOML parser and adds
    ## only the skip block, so it applies whatever the text layout. Dates are
    ## taken from the manifest's own blocks so a data update cannot
    ## invalidate them.
    path = joinpath(pkgdir(BVDOutbreakSize), "data", "observations.toml")
    raw = TOML.parsefile(path)
    write_manifest(r) = (
        p = joinpath(mktempdir(), "observations.toml");
        open(io -> TOML.print(io, r), p, "w"); p
    )
    ## The baseline carries no skip, whatever the manifest lists.
    delete!(raw, "skip_dates")
    base = load_observations(write_manifest(raw))

    deaths = raw["confirmed_death_history"]
    skip_date = String(deaths["dates"][end - 1])
    zone_date = String(raw["zone_confirmed_history"]["dates"][end - 1])
    raw["skip_dates"] = Dict(
        "confirmed_death_history" => [skip_date],
        "province_death_history" => [skip_date],
        "zone_confirmed_history" => [zone_date],
        "zone_death_history" => [zone_date],
    )
    obs = load_observations(write_manifest(raw))

    ## The skipped vintage is gone from the loaded history and every other
    ## vintage is unchanged, so the increment across the gap is the sum of
    ## the two it replaces.
    cd = base.confirmed_deaths_history
    sd = obs.confirmed_deaths_history
    @test length(sd.days) == length(cd.days) - 1
    @test sd.days == cd.days[[1:(end - 2); end]]
    @test sd.counts == cd.counts[[1:(end - 2); end]]
    @test sd.counts[end] - sd.counts[end - 1] ==
        cd.counts[end] - cd.counts[end - 2]

    ## Only the named block changes: cases keep every vintage.
    @test obs.confirmed_history == base.confirmed_history
    @test obs.province_confirmed_history == base.province_confirmed_history
    for (prov, h) in obs.province_death_history
        @test length(h.days) == length(base.province_death_history[prov].days) - 1
    end

    ## Zone blocks drop the date from every zone in both blocks.
    for blk in (:zone_confirmed_history, :zone_death_history)
        for (prov, zones) in getproperty(obs, blk), (zone, h) in zones
            b = getproperty(base, blk)[prov][zone]
            @test length(h.days) == length(b.days) - 1
        end
    end
end

@testitem "a skip date that matches no vintage is rejected" begin
    using BVDOutbreakSize
    using BVDOutbreakSize: load_observations
    using TOML

    ## A transposed digit would otherwise skip nothing while appearing to
    ## handle a revision, so the loader refuses it, as it does for an
    ## unknown block name.
    path = joinpath(pkgdir(BVDOutbreakSize), "data", "observations.toml")
    write_manifest(r) = (
        p = joinpath(mktempdir(), "observations.toml");
        open(io -> TOML.print(io, r), p, "w"); p
    )
    raw = TOML.parsefile(path)
    raw["skip_dates"] = Dict("confirmed_death_history" => ["2026-01-01"])
    @test_throws ErrorException load_observations(write_manifest(raw))

    raw = TOML.parsefile(path)
    raw["skip_dates"] = Dict(
        "no_such_history" => [String(raw["confirmed_death_history"]["dates"][end])],
    )
    @test_throws ErrorException load_observations(write_manifest(raw))
end

@testitem "skipping a zone vintage drops the revision figures on it" begin
    using BVDOutbreakSize
    using BVDOutbreakSize: load_observations
    using TOML

    ## A printed 24h figure must sit on a vintage the zone history carries,
    ## so the revision blocks follow their zone block's skip list. The date
    ## is the latest one carrying a revision figure in the manifest.
    path = joinpath(pkgdir(BVDOutbreakSize), "data", "observations.toml")
    raw = TOML.parsefile(path)
    delete!(raw, "skip_dates")
    write_manifest(r) = (
        p = joinpath(mktempdir(), "observations.toml");
        open(io -> TOML.print(io, r), p, "w"); p
    )
    base = load_observations(write_manifest(raw))
    rev_dates = String[
        String(d) for (prov, zones) in raw["zone_death_revision_24h"]
            if zones isa AbstractDict
            for sub in values(zones) for d in sub["dates"]
    ]
    date = maximum(rev_dates)
    raw["skip_dates"] = Dict(
        "zone_confirmed_history" => [date],
        "zone_death_history" => [date],
    )
    obs = load_observations(write_manifest(raw))
    day = only(
        setdiff(
            first(values(first(values(base.zone_death_history)))).days,
            first(values(first(values(obs.zone_death_history)))).days,
        )
    )
    for blk in (:zone_confirmed_revision_24h, :zone_death_revision_24h)
        for (prov, zones) in getproperty(obs, blk), (zone, h) in zones
            @test !(day in h.days)
        end
    end
    @test any(
        day in h.days for (prov, zones) in base.zone_death_revision_24h
            for (zone, h) in zones
    )
end
