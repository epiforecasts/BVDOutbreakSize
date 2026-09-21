# Pieces both fit reports use: the shape a fit result comes back in, and how
# a diagnostic number is printed. `summary.jl` reports one fit to its job
# summary and `convergence.jl` decides whether the fits are good enough to
# publish, and the two had these in common. Neither file is in
# `FIT_SOURCE_FILES`, so moving them here costs no refit.

## The frozen fits return `(; cutoff, o, chn)` rather than the chain itself.
fit_chain(x) = x isa NamedTuple && haskey(x, :chn) ? x.chn : x

## A count, or `n/a` where the quantity is undefined. An effective sample
## size over a degenerate parameter has no value to round.
fmt_count(x) = isfinite(x) ? string(round(Int, x)) : "n/a"

## Three significant figures, or `n/a` on the same terms.
fmt_value(x) = isfinite(x) ? string(round(x; sigdigits = 3)) : "n/a"
