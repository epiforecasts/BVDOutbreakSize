# Synthetic delay estimates

Invented numbers in the shape of the `bvd-internal-cmmid` results tables that
`scripts/linelist/delays.jl` reads. They exercise the arithmetic that turns
those tables into priors; they describe nothing.

The real estimates are not here and must not be. They are fitted in a private
repository whose disclosure rules permit fitted parameters to be committed there
and not shared onward, and this repository is public. `delays.jl` reads them at
run time from `BVD_DELAY_DIR`.

The values are chosen so the delta-method widths come out exact rather than
approximate, which is what makes the test a check on the arithmetic rather than
a restatement of it. For `any matched contact`: mean 20, SD 10 give shape 4 and
scale 5, and a `meanlog` 95% interval half-width of 0.098 gives a standard error
of 0.05, so the propagated widths are 0.4 on the shape and 0.25 on the scale.
