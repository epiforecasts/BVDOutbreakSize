## Tests for the per-health-zone data: the zone history blocks in the
## manifest, the health-zone metadata table and map, the increment helper
## that reshapes the histories for a zone composition model, and the
## cut-off freeze.

@testitem "zone histories partition the province and national cumulatives" begin
    using BVDOutbreakSize: load_observations, PROVINCE_SOURCE_NAMES

    obs = load_observations()
    for (zh, ph, nat) in (
        (obs.zone_confirmed_history, obs.province_confirmed_history,
        obs.confirmed_history),
        (obs.zone_death_history, obs.province_death_history,
        obs.confirmed_deaths_history))
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
                           PROVINCE_SOURCE_NAMES

    obs = load_observations()
    rows = load_health_zones()
    in_blocks = Set((p, z) for (p, zs) in obs.zone_confirmed_history
    for z in keys(zs) if z != "unallocated")
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
end

@testitem "health-zone map carries every zone and stays small" begin
    using BVDOutbreakSize: load_health_zones

    path = joinpath(@__DIR__, "..", "data", "health_zones.geojson")
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
    ## The map covers the seven provinces' zones, keyed or not.
    @test count("\"zone\":\"\"", text) > 0
    @test count("\"type\":\"Feature\"", text) ==
          length(keyed) + count("\"zone\":\"\"", text)
end

@testitem "zone_increment_matrix totals equal the zone-sum increments" begin
    using BVDOutbreakSize: load_observations, zone_increment_matrix,
                           PROVINCE_NAMES, PROVINCE_MEMBERS

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
        for (i, (p, z)) in enumerate(r.zones)
            @test z != "unallocated"
            @test p in PROVINCE_MEMBERS[r.patch]
            c = zh[p][z].counts
            @test r.increments[i, :] == max.(diff(vcat(0, c)), 0)
        end
    end
    ## Every zone lands in exactly one patch.
    all_zones = reduce(vcat, (r.zones for r in res))
    @test allunique(all_zones)
    @test length(all_zones) ==
          sum(count(!=("unallocated"), keys(zs)) for zs in values(zh))
    ## The unallocated rows are what separate the totals from the province.
    ituri = res[findfirst(r -> r.patch == "ituri", res)]
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
    lone = zone_increment_matrix(zh, ["ituri", "nowhere"],
        Dict("ituri" => ["ituri"], "nowhere" => ["nowhere"]))
    @test length(lone) == 2
    @test isempty(lone[2].zones)
    @test size(lone[2].increments) == (0, 0)
    ## Zones on different vintages are refused.
    H = @NamedTuple{days::Vector{Int}, counts::Vector{Int}}
    bad = Dict("ituri" => Dict{String, H}(
        "a" => (; days = [1, 2], counts = [1, 2]),
        "b" => (; days = [1, 3], counts = [1, 2])))
    @test_throws ErrorException zone_increment_matrix(bad, ["ituri"])
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
