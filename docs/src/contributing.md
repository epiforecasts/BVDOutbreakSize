# Contributing

Issues and pull requests are welcome at [epiforecasts/BVDOutbreakSize](https://github.com/epiforecasts/BVDOutbreakSize).
The [home page](index.md) covers installing the package, re-fitting the model, rendering the report and updating the data.
This page covers the conventions for changing the project.

## Tasks

`Taskfile.yml` wraps the common commands, and `task --list` describes each one.
The ones used while making a change:

- `task format` runs Runic over `src/`, `test/`, `docs/`, `scripts/`, `benchmark/` and `ext/`.
- `task test` runs the full test suite.
  `task test-quick` skips the quality and AD items.
- `task fetch-fits` downloads the fits from the latest successful docs build, and `task docs` renders the site from them.
- `BVD_FIT_ID=<id> task fit` fits and caches one model, and `task fit-all` fits them all.
  `julia --project=docs docs/fits/list.jl` lists the ids.
- `task check-convergence` runs the convergence gate CI applies before publishing.
- `task release-notes` prints the notes the next release would publish.

To render one report page from the cached fits, set `BVD_DOC_PAGE` to its path under `docs/pages/`:

```bash
BVD_DOC_PAGE=estimates/national julia --project=docs docs/execute.jl
```

The scripts under `scripts/` differ in which Julia project they need.
Read [`scripts/README.md`](https://github.com/epiforecasts/BVDOutbreakSize/blob/main/scripts/README.md) before running one.

## Model code

The model is built from Turing submodels in three layers under `src/models/`.
`priors.jl` holds the building blocks, one per parameter family, each owning its own priors.
`observations.jl` holds one observation submodel per data stream, each taking the growth state and adding its likelihood.
`joint.jl` holds the composers that stitch these into full models, from the single-stream fits up to `bvd_joint`.
The [Methods](methods.md) page describes the model and the [API reference](lib/api.md) lists every submodel.

- Submodels compose via `~ to_submodel(...)`.
- A composer includes only the likelihoods for the streams it carries.
  A single-stream composer never instantiates the other observation submodels, because a discrete stream left sampled trips Turing's model check.
- Pass a stream as `missing` to drop its likelihood.
  `bvd_joint` with every stream missing is the generator for the prior and posterior predictive checks.
- Pass an optional component in as a submodel rather than switching it on with a flag, as in `background_pooling = background_pooling_model`.
  A flag reaches the model as a value rather than a type, so both arms are inferred on every build.
- Keep the NaN and Inf safe clamps (`safe_nbinomial`, `eps`-flooring of expected counts) when editing a likelihood.
  They guard against extreme NUTS warmup proposals.
- The AD backend is Mooncake reverse mode.
  The hand-written rules in `src/ad_rules.jl` are timed against the backend's own by `task benchmark-rules`.
  Run it after an AD backend upgrade and delete any rule that no longer pays for itself.
- Code, code comments and docstrings keep to 80 characters per line.
- Runic formats all Julia code.
  Its version is pinned in both `test/formatter/Project.toml` and `.pre-commit-config.yaml`, and `test/package/CodeFormatting.jl` fails if the two drift.
- The DocStringExtensions templates in `src/docstrings.jl` add the signature to every docstring, so a docstring carries only its prose.
- Comments describe the code as it is.
  The history of a change belongs in [News](news.md).

### Closures in model code

A closure in a model body must not box its captures.

Julia boxes a local variable when a closure captures it and something later reassigns it.
The variable becomes a `Core.Box` and is type-unstable at every use.
Mooncake answers type instability with a `DynamicDerivedRule`, a dictionary lookup per call site on every gradient evaluation, and Enzyme's reverse mode cannot differentiate through a box at all.
Removing the boxes from the observation models roughly halved the gradient with the log-density bit-identical.

A box has no failure mode, which is why one can sit in the code for a long time.
Gradients, the log-density and the parameter count are all correct, so the only symptom is an unexplained ratio of gradient time to primal time.
A profiler does not localise it either, because the cost is charged to the frame that captures the variable rather than to the loop that pays it.

Three idioms produce a box.

- A widening guard that writes back in place, `if eltype(x) === Any; x = convert(...); end`, where a comprehension or `map` nearby reads `x`.
- An accumulator rebound inside a loop whose comprehensions capture it.
- A variable written on both arms of an `if`/`else` and then captured.

The mere possibility of reassignment creates the box, whether or not the branch ever runs.
A guard that almost never fires is therefore not free on the AD path.
The branch itself costs nothing.
The box it creates costs on every gradient.

The fix in every case is to compute into a separate binding and assign the captured name exactly once.
Write a branch as an expression whose value is assigned, `x = if cond ... else ... end`, or as a ternary, rather than assigning `x` inside it.

`test/test_boxed_captures.jl` enforces this for every method defined under `src/models/`.
It walks the module's bindings rather than calling `methods` on the model constructors, because `@model` puts the model body in a gensym-named evaluator method and that generated method, not the user-facing constructor, is what AD differentiates.
To check one method by hand:

```bash
julia --project=. -e '
using BVDOutbreakSize
m = first(methods(BVDOutbreakSize.gate_before))
println(any(x -> occursin("Core.Box", string(x)),
    Base.uncompressed_ast(m).code))'
```

`@code_warntype` reports the same thing, as a `Core.Box` in the variable list.

## Tests

The tests use TestItemRunner.
Each `test/test_*.jl` file holds the `@testitem`s for one feature.
`test/runtests.jl` collects every item with `@run_package_tests`, so a new file needs no registering.

Tags split the suite across CI jobs.

- `:quality` marks the Aqua, JET, formatting and doctest items.
- `:ad` marks the AD gradient checks.
- `:slow` marks the items that run full NUTS fits.

`runtests.jl` reads test arguments to choose among them.
`skip_quality` drops the quality and AD items, `quality_only` and `ad_only` run one tag, and `fast` and `downgrade` drop all three.

To run one file, point `TestItemRunner.run_tests` at the `test/` directory with a filter and run it with `--project=test`:

```julia
using TestItemRunner
root = joinpath(pwd(), "test")
TestItemRunner.run_tests(
    root; filter = ti -> ti.filename == joinpath(root, "test_renewal.jl")
)
```

Scope the filter to this `test/` directory, or copies of the test files in sibling worktrees are collected too.
Read the `Test Summary` line rather than the exit code.

## Report pages

Each rendered page is a Literate file under `docs/pages/`, in a folder per navigation group.
Every page includes `docs/pages/_setup.jl`, which loads the observations and every fit through the cache.
Anything two pages need lives there rather than on whichever page defined it first.

`docs/execute.jl` renders one page.
`docs/make.jl` renders the pages and assembles the Vitepress site, one stage at a time under `BVD_DOCS_STAGE`.
CI renders each page in its own job and then runs the combine stage.

A new page needs:

- its Literate file under `docs/pages/`
- an entry in `PAGES`, a render stage and a navigation entry in `docs/make.jl`
- an entry in the page list in `docs/execute.jl`
- an entry in the render matrix in `.github/workflows/docs.yml`
- its rendered markdown in `.gitignore`

Setup and table-construction code sits inside `<details>` dropdowns via `#md # @raw html` blocks.
The bare result object follows with `#hide`, so only the output renders.

The shared front matter (title, authors, abstract, scope) is single-sourced in `README.md`, up to the `<!-- SHARED:END -->` marker.
Edit it in `README.md` only.
`docs/front_matter.jl` reads it at build time and fills in the dates.
`docs/make.jl` copies the whole README to the home page, so do not duplicate it into a report page.

## Fits and the fit cache

Fits are cached under `logs/fit_cache`, keyed on a content hash.
`fit_content_hash` in `docs/fits/registry.jl` builds the hash from three inputs.

- The bytes of each file in `FIT_SOURCE_FILES`: the files under `src/models/`, `renewal.jl`, `sampling.jl`, `constants.jl`, `data.jl` and `onset_curve.jl` in `src/`, and the cache code itself.
- Every file under `data/` except those named in `FIT_DATA_EXCLUDE`.
- The cache schema version and the sampler settings, including `joint_sampler_args()` for the headline joint and its spatial control.

Any edit to one of those files, a comment included, changes the key for every fit.
In CI that is a cold refit of every model, which takes hours.
Plotting, summary and reporting code sits outside the key, so a report change reuses the cached fits.

Every file `scripts/score_releases.jl` writes into `data/` must be listed in `FIT_DATA_EXCLUDE`.
Otherwise the render's hash differs from the fit matrix's and every fit misses.

The render never fits.
A fit missing from the cache fails the build naming its key.
`BVD_FIT_STRICT=false` restores inline fitting, for a page run outside the cache entirely.

A cached fit does not survive a version change in Turing or its dependencies.
Refit rather than debugging a `KeyError` on a stale chain.
Any change to the model, the priors or the data needs a refit before its results mean anything.

Every NUTS fit writes a progress log to `logs/<fit>.log` and a TensorBoard run under `logs/tensorboard/<fit>/`.
`BVD_FIT_LOG` controls this (`all` when unset, or `progress`, `tensorboard`, `none`).
`task tensorboard` shows every fit in the worktree.

## Prose

- One sentence per line in markdown and write-up prose.
  Do not wrap prose at 80 characters.
- UK English throughout.
- Write in the present tense and describe the current state.
  Development history belongs in [News](news.md) and nowhere else.

### Analysis report prose

These apply to the narrative prose in the report pages under `docs/pages/`, and to write-up prose generally.
Use the existing report text as the template for tone.
The measured sentence- and paragraph-level rules below were reverse-engineered from a manuscript the maintainers are happy with.
The repo-specific rules that follow take precedence where the two disagree.
`AGENTS.md` carries a short version of this section for agent sessions, pointing back here for the full rule.
Keep the two in agreement: if this section changes in a way that affects the summary, update `AGENTS.md` too.

- No code references in the narrative.
  Do not name functions, parameters, files, or `:symbols` in the prose.
  Describe each quantity in words, and define a derived quantity in words the first time it appears, near its figure or table.
- Concise and direct.
  Cut filler and adjectives.
  Avoid the LLM-indicator words: comprehensive, leverage, robust, framework (when vague), utilise, facilitate, novel, landscape, foster, harness, streamline, pivotal, nuanced, multifaceted, cornerstone, synergy, overarching.
  Avoid filler and stance markers: very, significant, crucial, essential, clearly, obviously, "it is important to note that".
  Delete the LLM connective openers `Furthermore`, `Moreover`, `Additionally`, `In conclusion`, `Overall`.
  Restrict genuine connectives to `However,` `Whilst`, `Yet`, `Unlike`, `Similarly`, `Alternatively`, `Finally,`, `Instead`.
- Target 15-30 words per sentence and treat 40 as the ceiling.
  Exceed it only for a flat, colon-led enumeration, never for nested subordination.
  Use at most two commas per sentence for internal structure, a third only in a serial list.
- Never use a dash as punctuation.
  Replace it with a full stop, or with parentheses when the aside is a gloss or example.
  Use a colon only to introduce a list or an expansion, never to join two independent statements.
  Use a semicolon only to separate list items that themselves contain commas.
- Instead of a trailing qualifier (", which means that...", ", thereby enabling..."), cut the comma and start a new sentence with `This`, `These` or `It also`.
  Keep any surviving trailing `, which` clause to about one sentence in ten, and only where the clause carries the point rather than padding it.
- Write paragraphs of 2-6 sentences.
  Open each with its claim, not with context, and close it with a consequence or a limitation, never a restatement of the opening.
- Hedge with `can`, `could`, `may` or `potential`, one hedge per claim.
  Never stack hedges ("may potentially", "could possibly suggest").
- Report intervals as sentences, without a leading median.
  Write the credible interval as a phrase, not a "median (lower, upper)" construction.
- Report numbers with their provenance and units.
  Delete a vague quantity ("a large number of", "significantly more") rather than leave it unsourced.
- UK English throughout.
- Section and subsection titles are just the title.
  No descriptive suffix after a title (not "Reproduction number — weekly random walk with intervention ramp", just "Reproduction number"), and no detail-dump in the first sentence after a heading.
- Order the methods generatively, infections through to observation endpoints: the infection process first, then the epidemiological processes (delays, case-fatality ratio), then the observation models (surveillance streams before exports), then the joint model.
- Define every quantity before it is used.
  Define the reproduction number before the seeding that relies on it; introduce the initial infection count before describing how it arises; define every symbol and operator (including convolution) the first time it appears.
  Never use a symbol the reader has not met.
  Give an equation once in display maths and refer to its symbols inline afterwards; never narrate a displayed equation in words.
- Do not repeat.
  State a convention once (the credible-interval levels, the delay discretisation) and do not restate it per bullet or subsection.
  Cut sentences that duplicate earlier content.
- Cite the source of each prior and carry the uncertainty the source reports.
  When a source gives a distribution with uncertainty (a shape and scale with intervals), propagate that, not a self-assigned weakly-informative spread.
  Do not write "with an assumed weakly-informative spread" repeatedly.
  If a prior is our own choice, say so plainly ("we use a prior of ...").
- State assumptions as assumptions ("we assume a single seed case", "we assume the response scale-up takes about three weeks").
  Do not assert a false rationale for a modelling choice (not "a Poisson because the count is small").
- Do not editorialise or justify priors in the narrative (not "a diffuse prior would let the background absorb the whole stream").
  State what the model does.
- Methods belong in the methods.
  Do not leave model description (the intervention model, the counterfactual, the forecast, the evaluation) in the results.
  Move it to the methods and keep the results to findings.
- No project or issue history in the narrative: no PR or issue numbers, no "we previously used X", "this was changed to Y", "after refactoring".
  State the current design only.
  A difference from work being replicated, or a limitation of an external dependency, is the exception, written as a present-tense fact.
  An issue number may also stay where the issue is the record for something still unresolved, such as a workaround pointing at the open root cause it does not fix, or a data decision whose reasoning lives in the issue.
  The test is whether the number still tells the reader something they cannot get from the code.
  "Tracked by #495" earns its place; "changed in #495" does not.
- Label quantities accurately.
  Do not call suspected cases onsets, and prefer "current cumulative" over "final cumulative".
- For a latent quantity (infections, onsets, deaths) report the modelled estimate without overlaying observed data that sits downstream of unmodelled processes.
- Plots use the same credible-interval ribbons as the tables, not a bare median, and show only the period being estimated rather than greying out the rest.
- Model code shown in the report is clean.
  Strip working comments before it is displayed.
- Flag a future improvement as a GitHub issue, not a buried caveat in the prose.
- Bullet lists, bold and italics are advisory rather than enforced here, since this is a technical report rather than a journal manuscript.
  Do not add new bullet lists or bold to the narrative prose, and do not restructure an existing list that carries genuinely parallel content, but do strip bold or italics used for mid-sentence emphasis.

## Commits

Commit messages follow Conventional Commits, `type(scope): summary`, with a lower-case imperative summary.
The types in use are `feat`, `fix`, `docs`, `test`, `refactor`, `perf`, `ci`, `chore` and `style`, plus `data` for a change to the observations.
The scope names the area touched, as in `fix(scripts)`, `docs(news)` or `perf(ad)`.

## News

Every user-visible change adds an entry to [News](news.md) as part of the change itself.
The file holds one section per version, newest first.
The top section is the open one, and its heading matches the version in `Project.toml`.
Add to that section rather than starting a new one.

Within a version, entries sit under `### Model`, `### Data`, `### Report`, `### Performance`, `### Fixed`, `### Infrastructure` or `### Dependencies`.
An entry says what changed and why, in the present tense, with the pull request or issue number in brackets.
Say whether fitted values change.
Write one sentence per line.

## Pull requests and CI

`main` is branch-protected, so changes go through pull requests to `epiforecasts/BVDOutbreakSize`.
Everything here is slow.
The full test suite takes a long time and a full docs build fits every model.
Open the pull request early and let CI do the long work.
Run `task format` before every push, then only narrow checks locally: `task test-quick`, a single-file test run, or one rendered page.

CI runs the test suite (`.github/workflows/test.yml`) and builds the docs, publishing `output/` as a GitHub Release on each push to `main` (`.github/workflows/docs.yml`).

On a pull request each of those runs only when the change touches something it is built from.
The test suite and coverage need `src/`, `ext/`, `test/`, `data/`, `Project.toml`, and `docs/fits/` and `scripts/` because test items include files from both.
The report needs `src/`, `ext/`, `data/`, `docs/`, `scripts/`, `README.md` and `Project.toml`.
A workflow that skips says so in the summary of its `changes` job, so a skipped build is visible rather than being an absent check.
A push to `main`, a tag and a manual run are never gated.

The lists live in each workflow's `changes` job and are checked by `.github/actions/changed-paths/patterns_test.sh`, which pre-commit runs whenever one of them is edited.
Widen the list when something new feeds a build: a pattern that is too narrow skips the job that would have caught the change, and nothing reports that as a failure.

## Releases

A release is cut by commenting `@release` on any issue or pull request.
`.github/workflows/release.yml` tags `main`, publishes a GitHub release whose notes are the newest section of `docs/src/news.md`, and opens a pull request bumping the version and starting the next section.
`@release minor` and `@release major` choose the size of that bump; plain `@release` is a patch.
`task release-notes` prints what would be published, so the notes can be read before anything is cut.
