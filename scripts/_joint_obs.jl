# Joint-observation packing extracted from docs/examples/analysis.jl so the
# investigation script reuses the exact headline configuration.

function _increments(v)
    d = similar(v, Int)
    prev = 0
    for i in eachindex(v)
        d[i] = v[i] - prev
        prev = v[i]
    end
    return d
end

function joint_obs(o; observe = true)
    _stream(h,
        s) = h === missing ?
             (Union{Missing, Int}[observe ? s : missing], [0]) :
             (observe ? _increments(h.values) :
              fill(missing, length(h.values)), h.offsets)
    rep, rep_off = _stream(o.reported_case_history, o.reported_cases)
    dth, dth_off = _stream(o.death_history, o.total_deaths)
    have_conf = o.confirmed_case_history !== missing ||
                o.confirmed_cases !== missing
    have_pervintage = have_conf &&
                      o.confirmed_case_history !== missing &&
                      o.tests_analysed_history !== missing
    if have_pervintage
        sa = o.tests_analysed_history
        sr = o.tests_received_history
        ch = o.confirmed_case_history
        conf_off = collect(ch.offsets)
        conf = observe ? Union{Missing, Int}[_increments(ch.values)...] :
               fill(missing, length(conf_off))
        keep = [i == 1 || sa.values[i] > sa.values[i - 1]
                for i in eachindex(sa.values)]
        aoff = collect(sa.offsets)[keep]
        a_inc = _increments(sa.values[keep])
        ridx = [findfirst(==(off), sr.offsets) for off in aoff]
        r_inc = _increments([sr.values[i] for i in ridx])
        aobs = Dict(off => a_inc[k] for (k, off) in enumerate(aoff))
        robs = Dict(off => r_inc[k] for (k, off) in enumerate(aoff))
        analysed = Union{Missing, Int}[get(aobs, off, missing)
                                       for off in conf_off]
        received = Union{Missing, Int}[get(robs, off, missing)
                                       for off in conf_off]
    else
        conf_total = o.confirmed_cases !== missing ? o.confirmed_cases :
                     o.confirmed_case_history === missing ? missing :
                     o.confirmed_case_history.values[end]
        conf,
        conf_off = have_conf ?
                   (Union{Missing, Int}[observe ? conf_total : missing],
            [0]) : (Union{Missing, Int}[], Int[])
        analysed = Union{Missing, Int}[]
        received = Union{Missing, Int}[]
    end
    ec_full = o.exported_cases_daily
    last_import = isempty(ec_full) ? nothing : findlast(!=(0), ec_full)
    export_last_offset = last_import === nothing ? 0 :
                         length(ec_full) - last_import
    _truncate(v) = v[1:max(length(v) - export_last_offset, 0)]
    ecases = isempty(ec_full) ? ec_full :
             (observe ? _truncate(ec_full) :
              fill(missing, length(_truncate(ec_full))))
    ed_trunc = _truncate(o.export_deaths_daily)
    edaily = observe ? ed_trunc : fill(missing, length(ed_trunc))
    if o.confirmed_death_history !== missing
        cdh = o.confirmed_death_history
        cdeath = observe ?
                 Union{Missing, Int}[_increments(cdh.values)...] :
                 fill(missing, length(cdh.values))
        cdeath_off = collect(cdh.offsets)
    else
        cdeath = Union{Missing, Int}[]
        cdeath_off = Int[]
    end
    return (deaths = dth, reported = rep, export_deaths = edaily,
        kw = (; reported_offsets = rep_off, death_offsets = dth_off,
            confirmed_cases = conf, confirmed_offsets = conf_off,
            samples_analysed = analysed,
            samples_received = received,
            confirmed_deaths = cdeath,
            confirmed_death_offsets = cdeath_off,
            exported_cases_daily = ecases,
            export_last_offset = export_last_offset,
            confirmed_epi_exclusion = nothing,
            tests_analysed = observe ? o.cumulative_tests_analysed :
                             missing, tests_offset = 0))
end
