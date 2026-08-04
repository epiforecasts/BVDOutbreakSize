# Prose style guide

Reverse-engineered from `index.qmd` and `case-studies.qmd` of `EpiAware/ComposableProbabilisticIDModels` (439 prose sentences measured).
Every rule is grounded in a measured pattern in that source, with a quoted example from it.

Measured baseline: median 22 words per sentence (24 main text, 20 supplement), 90th percentile 36, a quarter of sentences 15 words or shorter, 3-4% over 45.
Per sentence: 1.2 commas, 0.03 semicolons, 0.09 colons, 0 dashes, 0.3 parenthetical pairs.
98% of source lines hold exactly one sentence. Median 4 sentences per paragraph.

## Part A: Rules

1. Target 15-30 words per sentence and treat 40 as the ceiling; two-thirds of sentences sit in that band.
   "Modelling infectious disease dynamics can be a valuable tool for synthesising information and developing quantitative evidence for understanding and control." (24 words)

2. Exceed 40 words only when the excess is a flat enumeration under a colon, never nested subordination.
   "We then replicate three published analyses ...: a COVID-19 analysis ...; an extension ...; and an ordinary differential equation model ..." (58 words, one clause plus three items)

3. One sentence per line in the source; every prose line ends at a full stop.
   Display maths, fenced code, figure captions and table rows are exempt and keep their own lines.

4. Use at most two commas per sentence for internal structure; a third only in a serial list.
   "Diseases affect populations heterogeneously across age, location, and risk group."

5. Never use a dash as punctuation; there are zero em or en dashes in either file.

6. Use a colon only to introduce a list or an expansion, never to join two independent statements.
   "An ARIMA process is composed of three components: an autoregressive (AR) process, a moving average (MA) process, and a differencing operation."

7. Use a semicolon only to separate items in a list that already contains commas; all 11 in the source do this.

8. Instead of appending a trailing subordinate clause, start a new sentence with `This`, `These`, `It also` or `Such`; this is the signature move, with `This` opening 40 sentences.
   "Second, structures contain other structures as fields ... / This pattern enables models to be assembled like building blocks."

9. Keep trailing `, which` and `, enabling ...` clauses under one per ten sentences, and only where the clause carries the point (measured 6 and 5 per 100).

10. Use parentheses freely, but only for cross-references, citations, glosses of jargon, or concrete examples (28-40 per 100 sentences).
    "multiple dispatch (a programming paradigm where function behaviour is determined by the types of all arguments rather than just the first)"

11. Write paragraphs of 2-6 sentences; go past 8 only for an enumeration of parallel points (median 4, range 2-13).

12. Open each paragraph with the claim, not with context.
    "A key strength of our proposed approach, demonstrated through our three case studies, is its modular design."

13. End a paragraph by adding a consequence or a limitation, never by restating the opening.
    "However, they also show that current approaches often struggle to meet all of these criteria, and even when they do, may not be adaptable to new contexts."

14. Restrict connective openers to `However,` `Whilst` `Yet` `Unlike` `Similarly` `Alternatively` `Finally,` `Instead`.
    `Furthermore`, `Moreover`, `Additionally` and `In conclusion` appear zero times.

15. Hedge with `can`, `could`, `may` and `potential`, and nothing else (12, 5, 3 and 3 per 100 sentences).
    "Composable frameworks may also be key for enabling robust large language model assisted model construction."

16. State limits as facts, without apology.
    "Since our prototype does not yet support uncertain delay distributions, we use the mean parameter values."

17. Use UK English and prefer `whilst` to `while` (31 against 3): `parameterise`, `specialised`, `discretised`, `modelling`.

18. Refer to code by bare backticked identifier, in a sentence whose subject is what the code does.
    "we use `AR` to generate the latent log $R_t$ trajectory, `Renewal` to compute expected infections ..."

19. Quote a full line of code in prose only when that exact line is the point of the paragraph.
    "The key line `@submodel ε_t = generate_latent(latent_model.ε_t, n - p)` enables composition by delegating to whatever error model was provided."

20. Give an equation once in display maths, then refer to its symbols inline; never narrate a displayed equation in words.
    "$g_t$ is the generation distribution probability mass function (pmf)."

21. Cross-reference with `@fig-`, `@sec-`, `@tbl-`, letting the reference be the grammatical subject or sit in parentheses.
    "@fig-mishra shows model checks for each component"; "(@fig-sir E)". Panels are "@fig-sir E" or "Panels A-C", never "Figure 5, panel E".

22. Forward-reference detail rather than duplicating it.
    "For further details about the DSL structure and backend implementation, see @sec-materials-methods."

23. Keep Results at the level of what was done and what was found; parameter values, priors and code belong in Methods or the supplement.
    Main text: "a negative binomial observation model for case counts". Supplement: `HalfNormal(0.1)`.

24. Describe each figure panel once, in one sentence, in panel order.
    "Panel D compares the continuous serial interval distribution with its discretised form."

25. Write self-contained figure captions naming every panel and every colour; captions are the one place where length and repetition of Results text are correct.

26. Report numbers with provenance and units.
    "the CIS tested over four million swabs from more than 150,000 households at a cost exceeding £500 million"

27. Contrast with prior work by naming the specific difference and the reason; never describe changes to your own project.
    "The key difference from _Mishra et al_ is in the initialization."

28. Use no bullet lists in prose (zero across both files); enumerations run inline under a colon, or become a table.

29. Use bold only as a label in a table column; all six occurrences are table theme names such as `**Model structure**`.

30. Use italics only for cited author names and Latin terms: `_Mishra et al_`, `*a priori*`.

31. Signpost internal structure with HTML comments rather than extra headings or bold.
    `<!-- Strengths and weaknesses of Turing.jl -->` marks each Discussion block without rendering.

32. Define jargon inline at first use, then use the short form bare.
    "domain-specific languages (programming languages tailored to specific application domains ...)" then "DSL" thereafter.

## Part B: Anti-patterns to strip

1. **Over-length sentences.** Anything past ~40 words that is not a colon-led list.
   The source's own worst line is a 120-word run-on starting "`EpiNow2` is built around three core modeling components ... These are: ...". Split at every sentence boundary and give each component its own sentence.

2. **Multiple sentences on one line.** Grep for `. ` mid-line and break at each boundary; the 2-3% of source lines that violate this read worse than their neighbours.

3. **Trailing comma-appended qualifiers.** `, which means that ...`, `, thereby enabling ...`, `, highlighting the importance of ...`.
   Cut the comma and start a new sentence with `This`.

4. **Dashes as punctuation.** Replace with a full stop, or parentheses when the aside is a gloss or example.

5. **Colon or semicolon joining independent clauses.** Replace with a full stop; keep the colon only where a list or expansion follows.

6. **Three or more structural commas in one sentence.** Split. Serial lists are the only exception.

7. **Restatement.** A paragraph ending by re-saying its first sentence, or a Results sentence repeating the Methods description of the same model. Delete the second instance and cross-reference.

8. **Methods content in Results.** Priors, tuning settings, package versions, solver choices. Move them and leave a `@sec-` reference.

9. **Project history.** Issue and PR numbers, "we previously used X", "this was changed to Y", "after refactoring", "in an earlier version". Zero of these appear in the source; state the current design only.
   Legitimate exceptions are differences from work being replicated and limitations of external dependencies, both written as present-tense facts.

10. **Pasted code or long expressions in prose.** Reduce to a bare backticked identifier, or move into a fenced chunk with one sentence above saying what it does.

11. **Equations narrated in words, or symbols redefined after the display block.**

12. **Hedge stacking.** "may potentially", "could possibly suggest", "it appears that ... might". One hedge per claim, from `can`, `could`, `may`, `potential`.

13. **Filler intensifiers and stance markers.** `very`, `significant`, `crucial`, `essential`, `clearly`, `obviously`, `it is important to note that`, `it is worth noting`. Near-zero in the source; delete without replacement.

14. **LLM connective openers.** `Furthermore`, `Moreover`, `Additionally`, `In conclusion`, `Overall`. Delete, or use `However,` where a genuine contrast exists.

15. **LLM vocabulary.** `comprehensive`, `practitioner`, `landscape`, `utilise`, `foster`, `harness`, `pivotal`, `nuanced`, `multifaceted`, `cornerstone`, `synergy`, `overarching`, `delve`, `realm`, `holistic`, `seamless`.
    `framework` is fine attached to a concrete noun ("Bayesian hierarchical framework") and must go when vague.

16. **Bold or italics for emphasis.** Strip; if a point needs weight, give it its own short sentence.

17. **Bullet lists in prose.** Convert to an inline colon-led list, or a table when there are more than about four structured items.

18. **Undefined jargon and bare acronyms.** Add a parenthetical gloss at first use, then use the short form.

19. **US spelling and `while`.** Convert to UK English and `whilst`.

20. **Vague numbers.** "a large number of", "significantly more". Replace with the figure and its source, or delete the claim.
