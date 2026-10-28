# Vignette before/after audit — 2026-10-22

Baseline: `git archive 12e0145` (pre-correction HEAD) installed to an isolated
library. Current: working tree, R sources matching `../source-sha256.txt`.
Vignette sources identical at both revisions. Each of the 20 Rmd files was
knitted in a fresh R process per library (one BLAS thread), saving the .md,
figures and every non-function object of the knit environment.

## Result

- 19/20 knitted under both versions. `articles/class-enumeration.Rmd` FAILED
  under the current package: its hidden check matched the old thin-profile
  warning text (`grepl("members for", ...)`). Fixed to
  `"membership mass for"`; re-knit ok and matches the scratch-patched knit.
- Saved objects: after removing the additive `gaussian_diagnostics` field and
  the enumeration columns `min_standardized_eigenvalue`/`min_gaussian_mass`,
  every object is `identical()` except (a) reworded `latents_thin_profiles`
  warning strings and (b) stored calls built by `do.call()` that inline the
  changed `multilpa` body; their arguments are identical.
- Prose (non-output lines incl. inline `r`): 0 lines differ in all 20.
- Printed numbers: every baseline number appears in current output, except
  "35 columns" -> "37 columns" (workflow-lca/-lpa) and "19 tables" -> "20
  tables" (lpa summary).
- Figures: 0 of all PNGs differ (md5).
- No vignette uses zero sampling weights, full-covariance LTA at the floor,
  or staged bootstrap, so the three result-changing corrections are not
  exercised by any vignette.

## Not changed but worth a decision
- README.Rmd:177, enumeration.Rmd:129, articles/evaluation.Rmd:57 still say a
  thin profile "may be spurious" / "fewer effective members"; the warning
  now says "Gaussian membership mass ... not an identification test".
  Not wrong, but the tone now differs from the warning.
