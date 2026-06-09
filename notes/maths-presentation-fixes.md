# Maths presentation fixes (sentence-maths style)

Roadmap from the maths-presentation audit. Apply as the FINAL presentation
pass, after the model changes land (delays natural params, lab-model bundle,
export-death fix, λ_bg_death switch-on), so the final equations are formatted
once. Line numbers are pre-model-change and will drift.

Style: one key quantity per display ```math block, prose sentence introducing
it, no inline clutter; split a mean/intensity from its likelihood so the
`~ NegBinomial(...)` sits on its own line (not as a trailing clause with a full
stop).

## Priority 1 (maintainer-flagged)
1. **Suspect-pool composition φ_v inline → display** (~1185-1188). Pull
   `φ_v = (p_drc·bvd)_v / ((p_drc·bvd)_v + λ_bg,v)` onto its own ```math line.
2. **Laboratory-volume likelihood split** (~1174-1179). `v_t = …` on its own
   line; then `Σ v_t ~ NegBinomial(…)` on a separate ```math line.
3. **Confirmed-deaths likelihood split** (~1264-1270). `p_cd = logistic(…)` on
   its own line; then `Σ p_cd·m_t ~ NegBinomial(…)` separately.

## Priority 2
4. **q_susp composition inline → display** (~1256-1258).
5. **Cases mean c_t split from likelihood** (~1107-1110). `c_t = p_drc·bvd_t +
   λ_bg` then the increment `~ NegBinomial`.
6. **Deaths mean m_t split from likelihood** (~1140-1144). NOTE: this equation
   gains a `+ λ_bg_death` background term once background deaths are switched
   on — format after that change.
7. **Convolution + bvd_t definitions inline → display** (~1095-1099). `(x*f)_t
   = Σ x_{t-s} f_s` and `bvd_t = Σ onsets_{t-s} f_rep,s` each on their own line.
8. **Export intensity: C_t, det_t inline → display** (~1305-1307). NOTE: the
   export-death model is changing (CFR on pre-ascertainment exported cases), so
   format the export equations after that.
9. **Export-death intensity μ_t / Λ_d split** (~1349-1351). Same caveat as 8.

All proposed rewrites follow: prose intro → display mean/intensity → prose →
display likelihood.
