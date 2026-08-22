# Synthetic line-list fixture

Invented numbers, in the schema `scripts/linelist/` expects.
Nothing here came from a case line list or from any real outbreak, and no estimate produced from it means anything.

It exists so the line-list refit can be run and tested without the real inputs, which are individual-patient-record derived and live outside this repository.

`linelist_streams.csv` is a logistic epidemic over 1 May to 31 July 2026: cumulative confirmed cases, cumulative reported cases at about four times that, and the daily new suspected series differenced from the reported one.
It is a daily series, indexed by event date.

`linelist_streams_known.csv` is that same outbreak under the known-by indexing the fits actually use, which is what `fit_single.jl` reads: the daily series thinned to weekly snapshots, the two cumulative histories carrying the total known at each snapshot and the suspected series summed over the interval since the last one, so no case is dropped by the thinning.

`onset_curve_scanned.csv` is the same outbreak seen as a reporting triangle, three vintages of one bell-shaped onset curve, each vintage seeing more of the recent past than the last so the between-vintage increments are non-zero.

All three are small enough to read and to edit directly, which is what to do if the schema changes.

`test/test_linelist_manifest.jl` runs the manifest substitution on this fixture.
`scripts/linelist/README.md` shows how to run a ten-sample fit against it.
