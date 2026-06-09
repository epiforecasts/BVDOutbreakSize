using BVDOutbreakSize, Serialization, Statistics
chn = deserialize("logs/joint_chain.jls")
obs = load_observations()

function q(sym)
    v = vec(collect(chn[sym]))
    return (round(median(v); digits = 3),
        round(quantile(v, 0.05); digits = 3),
        round(quantile(v, 0.95); digits = 3))
end

for s in (Symbol("rt_state.log_R0"), Symbol("rt_state.sigma_rw"),
    Symbol("rt_state.intervention_effect"), :r0, :r, :R_T, :T,
    :doubling_time, :C_T)
    m, lo, hi = q(s)
    println(rpad(string(s), 30), rpad(m, 10), "[", lo, ", ", hi, "]")
end

lr0 = vec(collect(chn[Symbol("rt_state.log_R0")]))
println("\nimplied R0 = exp(log_R0) median: ", round(exp(median(lr0)); digits = 3))
println("tmrca_days (genetic lower bound on T): ", obs.tmrca_days)

c = obs.reported_history.counts
d = obs.deaths_history.counts
nday = obs.reported_history.days[end] - obs.reported_history.days[1]
println("observed case doubling (days 85-93): ",
    round(nday * log(2) / log(c[end] / c[1]); digits = 1), " d")
println("observed death doubling (days 85-93): ",
    round(nday * log(2) / log(d[end] / d[1]); digits = 1), " d")
