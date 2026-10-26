## Tests for the per-health-zone data: the zone history blocks in the
## manifest, the health-zone metadata table and map, the increment helper
## that reshapes the histories for a zone composition model, and the
## cut-off freeze.

@testitem "zone histories partition the province and national cumulatives" begin
    using BVDOutbreakSize: load_observations, PROVINCE_SOURCE_NAMES

    obs = load_observations()
    for (zh, ph, nat) in (
            (
                obs.zone_confirmed_history, obs.province_confirmed_history,
                obs.confirmed_history,
            ),
            (
                obs.zone_death_history, obs.province_death_history,
                obs.confirmed_deaths_history,
            ),
        )
        @test !isempty(zh)
        @test issubset(keys(zh), PROVINCE_SOURCE_NAMES)
        days = first(values(first(values(zh)))).days
        @test issorted(days)
        @test allunique(days)
        ## One table per vintage, so every zone shares the vintage days.
        for zones in values(zh), z in values(zones)

            @test z.days == days
        end
        national = Dict(zip(nat.days, nat.counts))
        checked = 0
        for (i, day) in enumerate(days)
            total = 0
            for prov in keys(zh)
                @test haskey(zh[prov], "unallocated")
                s = sum(z.counts[i] for z in values(zh[prov]))
                if haskey(ph, prov)
                    j = findfirst(==(day), ph[prov].days)
                    if j !== nothing
                        @test s == ph[prov].counts[j]
                        checked += 1
                    end
                end
                total += s
            end
            @test haskey(national, day)
            @test total == national[day]
        end
        @test checked > 0
    end
    ## Both blocks describe the same zones on the same vintages.
    zc, zd = obs.zone_confirmed_history, obs.zone_death_history
    @test keys(zc) == keys(zd)
    for prov in keys(zc)
        @test keys(zc[prov]) == keys(zd[prov])
        for z in keys(zc[prov])
            @test zc[prov][z].days == zd[prov][z].days
            ## Deaths never exceed cases in a zone; the unallocated row
            ## can (A ventiler holds deaths with its cases cell NA).
            z == "unallocated" ||
                @test all(zd[prov][z].counts .<= zc[prov][z].counts)
        end
    end
end

@testitem "health-zone metadata matches the zone blocks" begin
    using BVDOutbreakSize: load_observations, load_health_zones,
        _default_health_zones, PROVINCE_SOURCE_NAMES

    obs = load_observations()
    rows = load_health_zones()
    in_blocks = Set(
        (p, z) for (p, zs) in obs.zone_confirmed_history
            for z in keys(zs) if z != "unallocated"
    )
    @test Set((r.province, r.zone) for r in rows) == in_blocks
    @test allunique(r.zone for r in rows)
    @test allunique(r.zscode for r in rows)
    for r in rows
        @test r.province in PROVINCE_SOURCE_NAMES
        @test r.population > 0
        @test !isempty(r.label)
        @test length(r.zscode) == 11
        ## Inside the DRC's bounding box.
        @test -14 <= r.lat <= 6
        @test 12 <= r.lon <= 32
    end
    ## Patch order, then alphabetical by key within a province.
    order = Dict(p => i for (i, p) in enumerate(PROVINCE_SOURCE_NAMES))
    @test issorted([(order[r.province], r.zone) for r in rows])
    ## The table the zone fit reads when no metadata is passed.
    @test _default_health_zones() == rows
end

@testitem "health-zone map carries every zone and stays small" begin
    using BVDOutbreakSize: load_health_zones

    path = joinpath(@__DIR__, "..", "src", "assets", "health_zones.geojson")
    @test filesize(path) < 1_500_000
    text = read(path, String)
    rows = load_health_zones()
    for r in rows
        @test occursin("\"zone\":\"$(r.zone)\"", text)
        @test occursin("\"zscode\":\"$(r.zscode)\"", text)
    end
    ## Every keyed feature is one of the metadata rows, once.
    keyed = [m.captures[1] for m in eachmatch(r"\"zone\":\"([a-z_]+)\"", text)]
    @test sort(keyed) == sort([r.zone for r in rows])
    ## The map covers the eight provinces' 178 zones, keyed or not.
    @test count("\"zone\":\"\"", text) > 0
    @test count("\"type\":\"Feature\"", text) ==
        length(keyed) + count("\"zone\":\"\"", text) == 178
    ## Nord-Ubangi's 11 zones are all drawn, Wapinda among the keyed ones.
    @test count("\"province\":\"nord_ubangi\"", text) == 11
    w = only(r for r in rows if r.zone == "wapinda")
    @test w.province == "nord_ubangi"
    @test w.zscode == "CSvyaLZJiEt"
end

@testitem "zone_increment_matrix totals equal the zone-sum increments" begin
    using BVDOutbreakSize: load_observations, zone_increment_matrix,
        zone_reattribution_days, PROVINCE_NAMES,
        PROVINCE_MEMBERS
    using Dates: Date, Day

    obs = load_observations()
    zh = obs.zone_confirmed_history
    res = zone_increment_matrix(zh, PROVINCE_NAMES)
    @test length(res) == length(PROVINCE_NAMES)
    @test [r.patch for r in res] == PROVINCE_NAMES
    days = first(values(first(values(zh)))).days
    for r in res
        @test r.days == days
        @test size(r.increments) == (length(r.zones), length(days))
        @test r.totals == vec(sum(r.increments; dims = 1))
        @test all(>=(0), r.increments)
        keep = [!(d in r.excluded) for d in days]
        for (i, (p, z)) in enumerate(r.zones)
            @test z != "unallocated"
            @test p in PROVINCE_MEMBERS[r.patch]
            c = zh[p][z].counts
            plain = max.(diff(vcat(0, c)), 0)
            @test r.increments[i, keep] == plain[keep]
            @test all(iszero, r.increments[i, .!keep])
        end
    end
    ## The reattribution vintages: Ituri's unallocated row falls on 18 June
    ## (cases and deaths) and on 22 July (cases), and no other province's
    ## ever does, so exactly those two Ituri columns are zeroed.
    dates(ds) = [obs.seeding + Day(d - 1) for d in ds]
    re = zone_reattribution_days(zh)
    @test collect(keys(re)) == ["ituri"]
    @test dates(re["ituri"]) == [Date("2026-06-18"), Date("2026-07-22")]
    @test dates(zone_reattribution_days(obs.zone_death_history)["ituri"]) ==
        [Date("2026-06-18")]
    ituri = res[findfirst(r -> r.patch == "ituri", res)]
    @test dates(ituri.excluded) == [Date("2026-06-18"), Date("2026-07-22")]
    @test all(r -> r.patch == "ituri" || isempty(r.excluded), res)
    ## Without a reattribution table every column is kept.
    plain = zone_increment_matrix(
        zh, PROVINCE_NAMES;
        reattribution = Dict{String, Vector{Int}}()
    )
    @test all(r -> isempty(r.excluded), plain)
    @test all(>=(0), plain[1].increments)
    ## Every zone lands in exactly one patch.
    all_zones = reduce(vcat, (r.zones for r in res))
    @test allunique(all_zones)
    @test length(all_zones) ==
        sum(count(!=("unallocated"), keys(zs)) for zs in values(zh))
    ## The unallocated rows are what separate the totals from the province.
    prov = obs.province_confirmed_history["ituri"]
    for (j, day) in enumerate(days)
        k = findfirst(==(day), prov.days)
        k === nothing && continue
        un = zh["ituri"]["unallocated"].counts[j]
        @test sum(zh["ituri"][z].counts[j] for (_, z) in ituri.zones) + un ==
            prov.counts[k]
    end
    ## No data, no rows; a patch with no zone data gets an empty matrix.
    @test isempty(zone_increment_matrix(Dict{String, Any}(), PROVINCE_NAMES))
    lone = zone_increment_matrix(
        zh, ["ituri", "nowhere"],
        Dict("ituri" => ["ituri"], "nowhere" => ["nowhere"])
    )
    @test length(lone) == 2
    @test isempty(lone[2].zones)
    @test size(lone[2].increments) == (0, 0)
    ## Zones on different vintages are refused.
    H = @NamedTuple{days::Vector{Int}, counts::Vector{Int}}
    bad = Dict(
        "ituri" => Dict{String, H}(
            "a" => (; days = [1, 2], counts = [1, 2]),
            "b" => (; days = [1, 3], counts = [1, 2])
        )
    )
    @test_throws ErrorException zone_increment_matrix(bad, ["ituri"])
end

@testitem "zone_cumulative_falls reports a revision out of a named zone" begin
    using BVDOutbreakSize: zone_cumulative_falls, zone_reattribution_days,
        zone_increment_matrix

    H = @NamedTuple{days::Vector{Int}, counts::Vector{Int}}
    ## `a` loses five cases at the third vintage and one at the fourth;
    ## `b` only rises, and the unallocated row never falls, so the
    ## unallocated rule alone sees nothing.
    hist = Dict(
        "ituri" => Dict{String, H}(
            "a" => (; days = [1, 4, 7, 9], counts = [10, 20, 15, 14]),
            "b" => (; days = [1, 4, 7, 9], counts = [1, 2, 8, 9]),
            "unallocated" => (; days = [1, 4, 7, 9], counts = [0, 0, 5, 6])
        )
    )
    falls = zone_cumulative_falls(hist)
    @test [(f.province, f.zone, f.day, f.from, f.to) for f in falls] ==
        [("ituri", "a", 7, 20, 15)]
    ## A one-unit correction is below the default guard and above `0`.
    @test [f.day for f in zone_cumulative_falls(hist; min_fall = 0)] == [7, 9]
    @test isempty(zone_cumulative_falls(hist; min_fall = 5))
    ## The unallocated rule alone leaves the province out; widening it to
    ## the named zones' falls adds that vintage, which the composition
    ## then drops.
    @test isempty(zone_reattribution_days(hist))
    widened = zone_reattribution_days(hist; include_zone_falls = true)
    @test widened == Dict("ituri" => [7])
    res = only(
        zone_increment_matrix(
            hist, ["ituri"], Dict("ituri" => ["ituri"]);
            reattribution = widened
        )
    )
    @test res.excluded == [7]
    @test all(iszero, res.increments[:, 3])
end

@testitem "a flagged zone fall takes the printed 24h figure" begin
    using BVDOutbreakSize: zone_increments, zone_increment_matrix,
        zone_cumulative_falls

    H = @NamedTuple{days::Vector{Int}, counts::Vector{Int}}
    days = [1, 4, 7, 9]
    ## `a` falls by five at day 7 (flagged) and by one at day 9 (below the
    ## guard); `b` rises at day 7. Each has a printed figure on its day.
    a = [10, 20, 15, 14]
    b = [1, 2, 8, 9]
    printed_a = (; days = [7, 9], counts = [3, 6])
    printed_b = (; days = [7], counts = [2])
    @test zone_increments(a, days) == [10, 10, 0, 0]
    ## The flagged fall takes its figure; the unflagged one keeps the clamp.
    @test zone_increments(a, days, printed_a) == [10, 10, 3, 0]
    ## A zone that does not fall keeps its cumulative difference.
    @test zone_increments(b, days, printed_b) == [1, 1, 6, 1]
    ## A printed negative is the revision itself and reads as no new count.
    @test zone_increments(a, days, (; days = [7], counts = [-2])) ==
        [10, 10, 0, 0]
    hist = Dict(
        "ituri" => Dict{String, H}(
            "a" => (; days, counts = a), "b" => (; days, counts = b)
        )
    )
    printed = Dict("ituri" => Dict{String, H}("a" => printed_a, "b" => printed_b))
    m = only(
        zone_increment_matrix(
            hist, ["ituri"], Dict("ituri" => ["ituri"]);
            reattribution = Dict{String, Vector{Int}}(), printed
        )
    )
    @test m.zones == [("ituri", "a"), ("ituri", "b")]
    @test m.increments == [10 10 3 0; 1 1 6 1]
    @test m.totals == [11, 11, 9, 1]
    plain = only(
        zone_increment_matrix(
            hist, ["ituri"], Dict("ituri" => ["ituri"]);
            reattribution = Dict{String, Vector{Int}}()
        )
    )
    @test plain.increments == [10 10 0 0; 1 1 6 1]
    @test [f.day for f in zone_cumulative_falls(hist)] == [7]
end

@testitem "the printed 24h revision figures sit on flagged falls" begin
    using BVDOutbreakSize: load_observations, zone_cumulative_falls,
        zone_increments, zone_increment_matrix, PROVINCE_NAMES
    using Dates: Date

    obs = load_observations()
    day(s) = (Date(s) - obs.seeding).value + 1
    ## Every figure sits on a flagged fall, and every flagged fall from
    ## SitRep 060 (13 July), the first table with 24h columns, has one.
    first_24h = day("2026-07-13")
    for (hist, printed) in (
            (obs.zone_confirmed_history, obs.zone_confirmed_revision_24h),
            (obs.zone_death_history, obs.zone_death_revision_24h),
        )
        flagged = Set(
            (f.province, f.zone, f.day) for f in zone_cumulative_falls(hist)
        )
        recorded = Set(
            (prov, zone, d)
                for (prov, zones) in printed for (zone, p) in zones
                for d in p.days
        )
        @test recorded ⊆ flagged
        @test Set(f for f in flagged if f[3] >= first_24h) == recorded
    end
    ## Bunia's deaths fall from 228 to 175 on 22 July; the printed 24h
    ## total of 7 replaces the clamped zero, and every other vintage keeps
    ## its cumulative difference.
    bunia = obs.zone_death_history["ituri"]["bunia"]
    v = findfirst(==(day("2026-07-22")), bunia.days)
    inc = zone_increments(
        bunia.counts, bunia.days, obs.zone_death_revision_24h["ituri"]["bunia"]
    )
    plain = zone_increments(bunia.counts, bunia.days)
    @test plain[v] == 0
    @test inc[v] == 7
    rest = setdiff(eachindex(inc), v)
    @test inc[rest] == plain[rest]
    ## Butembo's cases fall from 65 to 62 on 22 July with 4 printed. Nord-Kivu
    ## has no unallocated fall that day, so the cell reaches the composition,
    ## and it is the only case cell the figures change.
    res = zone_increment_matrix(
        obs.zone_confirmed_history, PROVINCE_NAMES;
        printed = obs.zone_confirmed_revision_24h
    )
    base = zone_increment_matrix(obs.zone_confirmed_history, PROVINCE_NAMES)
    nk = findfirst(r -> r.patch == "nord_kivu", res)
    i = findfirst(==(("nord_kivu", "butembo")), res[nk].zones)
    w = findfirst(==(day("2026-07-22")), res[nk].days)
    @test base[nk].increments[i, w] == 0
    @test res[nk].increments[i, w] == 4
    @test sum(count(r.increments .!= b.increments) for (r, b) in zip(res, base)) == 1
end

@testitem "freeze_observations truncates the zone blocks" begin
    using BVDOutbreakSize: load_observations, freeze_observations

    full = load_observations()
    frozen = freeze_observations("2026-07-15")
    zf, zfull = frozen.zone_confirmed_history, full.zone_confirmed_history
    @test keys(zf) == keys(zfull)
    for prov in keys(zf), z in keys(zf[prov])

        h, hf = zfull[prov][z], zf[prov][z]
        @test length(hf.days) < length(h.days)
        @test all(<=(frozen.n), hf.days)
        ## The frozen series is the head of the full one. Grid days count
        ## from the seeding date, so they do not move with the cut-off.
        @test hf.counts == h.counts[1:length(hf.counts)]
        @test hf.days == h.days[1:length(hf.days)]
    end
    ## 15 July is itself a vintage, so it is the last day kept.
    @test zf["ituri"]["bunia"].days[end] == frozen.n
    @test frozen.zone_death_history["ituri"]["bunia"].days ==
        zf["ituri"]["bunia"].days
end
