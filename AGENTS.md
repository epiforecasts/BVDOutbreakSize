# BVDOutbreakSize

Discrete-time renewal model of the 2026 BVD outbreak, fitted jointly to
multiple surveillance streams in Turing.jl.
The analysis report is `docs/examples/analysis.jl` (Literate + Documenter).

## Report prose standards

These apply to the prose in `docs/examples/analysis.jl` and any rendered
report text.
Use the existing report text as the template for tone.

- No code references in prose.
  Do not name functions, parameters, files, or `:symbols` in the narrative.
  Describe each quantity in words, and define a derived quantity in words
  the first time it appears near its figure or table.
- One sentence per line.
  Do not wrap at a column width; break on sentences.
- Concise and direct.
  Cut filler and adjectives.
  Avoid the LLM-indicator words: comprehensive, leverage, robust,
  framework (when vague), utilise, facilitate, novel, landscape, foster,
  harness, streamline, pivotal, nuanced, multifaceted, cornerstone,
  synergy, overarching.
- Report intervals as sentences, without a leading median.
  Write the credible interval as a phrase, not a "median (lower, upper)"
  construction.
- Minimise colons and dashes in prose; use them only when genuinely needed.
- UK English throughout.

## Code

- Max 80 characters per line; no trailing whitespace; no spurious blanks.
- Every model delay is sampled from a prior, never fixed.
- Run the SciML JuliaFormatter on changed files before committing
  (`bash scripts/run_formatter.sh`), and run the tests.
