# API selection and reporting review

The current 1.1.16 source passes its focused tests for the latest standard error, staged comparison, residual adjustment, wildcard and help page changes. Four older defects remain reproducible in selection and reporting. The help rewrite also breaks the registered documentation claims check, which runs in the equivalence workflow.

## Scope and evidence

Reviewed the current source against `8518dffe` (1.1.9, the CRAN version recorded in HANDOFF) and `a10a919b` (0.9.21, the older CRAN baseline used in NEWS). Local version and release status are repository evidence, not a fresh CRAN lookup. No production code, package tests, or git state was changed.

Read the relevant wildcard dispatch and table code, model comparison changes, categorical simulation, diagnostics boundary checks, cross-level results and the two CI workflows. Earlier Nov 4 findings were used as leads. Fresh synthetic cases with seeds 1106, 1116 and 1127 reproduce the four findings below; `probes.R` contains executable assertions, input structure inspection and evidence. `probes.log` retains output and session information. This is not a fresh mathematical audit of every estimator or a rerun of every external reference.

Environment: R 4.5.2, macOS arm64, one core, `NOT_CRAN=true`, single-thread BLAS settings.

## Findings still open

### High priority categorical bootstrap LRT loses joint item names

`R/enumeration.R:1102` converts the named list of categorical draws with `as.data.frame(drawn, stringsAsFactors = FALSE)`. Its default name repair changes `a:b` to `a.b`; the selection at line 1106 then requests the original `a:b` and fails. Categorical names containing spaces or punctuation have the same underlying exposure.

Fresh case: 225 rows in 45 groups, two categorical items joined with `local_dependence = "a:b"`, one continuous item, two profiles, one versus two group classes. Both fits converge. The observed LR statistic is 57.47968, but all three bootstrap replicates fail with `undefined columns selected`, leaving `p_value = NA` and `n_valid = 0`. Calling `.multilpa_simulate(pair_null)` directly gives the same error, so this is a simulation failure before any refit.

Smallest correction to evaluate: preserve column names with `check.names = FALSE`, then test joined items, non-syntactic single items and FIML partial pairs. Existing valid p-values would be unchanged; previously failed tests could return values.

History: the name-repair line already exists in 0.9.21. Joint local dependence was added after 0.9.21, exposing the older simulator bug through that feature. It is also present in 1.1.9.

### High priority categorical bounds receive clean diagnostics and no LRT warning

`R/diagnostics-checks.R:68` derives the check solely from `x$boundary` or `x$random_boundary`. `R/compare-lrt.R:48` uses the same stored flag for its boundary warning. Ordinary LCA's flag does not cover response probabilities at zero or one.

Fresh case: two-class LCA on 400 synthetic rows and five binary items. Two informative items are deterministic within class; the fitted response probabilities include `1e-10` and `1 - 1e-10`. The fit converges, `fit$boundary` is `FALSE`, and `get_checks()` returns `estimates at a bound: no`, status `ok`. `parameter_inference()` correctly refuses the same fit with `latents_boundary_fit` because its categorical probabilities are bound-active.

A nested comparison that frees a joint `c:d` response table reports LR 4.214299, df 2, p 0.121584 and emits no warning. This does not establish the correct reference distribution for that particular pair; it demonstrates that the package's promised boundary caution is absent where its own inference guard detects active bounds.

Correction to evaluate: a separate diagnostic predicate covering categorical and count boundaries, reused by the check, grid trust display and analytic LRT warning. Changing `fit$boundary` itself needs care because `bootstrap_lrt()` uses that flag to refuse inputs.

History: boundary diagnostics were introduced after 0.9.21, with the current predicate present by `480354e7` (0.97.5 work). Analytic LRT warning logic arrived in `c979dce9` (1.1.8). Both predate 1.1.9 and remain in 1.1.16.

### High priority cross level prevalence changes its meaning by form

`R/cross-level.R:452` stores `sum_h omega_h pi(k|h)` as the full model's `profile_prevalences`; `R/cross-level.R:530` describes that quantity as a profile's share of individuals, and line 544 reports it as `prevalence`. In the full form the formula represents the average group's composition. With unequal group sizes this differs from the share of individuals. Restricted fits report the individual share because their composition is common across group classes.

Fresh case: 60 groups, alternating sizes 4 and 35, with different generating compositions. Full model table: profile prevalences 0.579963 and 0.420037. Average individual posteriors: 0.325551 and 0.674449, a 25 percentage point difference. Restricted table: 0.325731 and 0.674269, matching its average individual posteriors within `3e-7`.

The group-based quantity is meaningful, but the table and comment do not explain its difference from the restricted model. Choose either a common individual-based column or explicit documentation of the full form's group-based quantity. The former changes reported prevalences; the latter changes their stated interpretation.

History: the column was introduced in `be277875` (0.9.32, Oct 16); it is new since 0.9.21 and already present in 1.1.9.

### Medium priority cross level sample counts include zero weight units

`R/cross-level.R:588` uses `x$n_groups` and `x$n_obs` in the fit table; the printer at lines 604 onward uses those totals too. The fitted object already has `n_used_groups` and `n_used_obs`, and `nobs()` reports the used groups correctly.

Fresh case: six of 60 groups have weight zero. The fit table reports 60 groups and 1170 rows, while `nobs()` reports 54 and the positive-weight subset table reports 54 groups and 1053 rows. Log likelihood and BIC agree between the weighted and subset fits to about `2e-12`, so this finding concerns displayed counts, not likelihood or criterion arithmetic.

Correction to evaluate: use the stored used counts, retaining a fallback for older saved fits. The raw-count table code predates 0.9.21; the later zero-weight contract leaves its reporting inconsistent with current `nobs()`.

## Latest changes verified

`test-review-1105-fixes.R`: 23 expectations pass. This includes an independent Hessian for held categorical measurement, charging a staged model its measurement parameters against a joint model, retaining conditional parameter counts when both fits hold the same measurement, leaving joint-only comparisons unchanged, BH defaults matching between the residual verbs, raw `adjust = "none"`, invalid adjustment refusal, and the new cross-level criterion warning.

The fresh cross-level fit also emits `latents_cross_level_criteria`, explicitly warning against choosing class count or form from its criteria. The Nov 4 finding about missing caution is therefore addressed; the prevalence inconsistency above remains separate.

`test-wildcards.R`: 207 expectations pass, including delegation to existing views, weighted posterior-derived displays, two-level output, all supported fitted families in its fixtures, failed regression candidate filtering, table stacking, availability refusals and overview return behavior.

`test-plot-docs.R`: seven expectations pass, covering individual help files and aliases and length constraints. Passing these checks does not establish that every argument explanation or external attribution is substantively correct.

Combined: 237 passed, zero failed, zero errors, zero warnings and zero skips. Saved structured test results: `test-results.rds`.

## Documentation claim validation fails

Fresh command: `Rscript equivalence/run.R '^claims$'` fails in `equivalence/test-claims.R:43`, which requires registered phrases still to appear in their documented source. The recent help rewrite removes or changes phrases registered in `equivalence/CLAIMS.csv`, including inference references and the empty-row convention. Evidence: `claims.log`. Root's independent count identifies 23 missing phrases out of 44 registered rows.

`.github/workflows/equivalence.yaml` invokes this check under its documentation claims step. This is a concrete CI regression in the uncommitted help rewrite. It is not evidence that estimates disagree with external software. Reconcile the registry with surviving documentation and decide where the backed external comparisons should be documented; retain substantive reference tests rather than removing them to make the registry green.

## Limits

No external software fitting, live CI result retrieval, real fork/socket recovery, bootstrap size calibration, full vignette rendering or visual inspection was performed in this slice. The root review coordinates the full package suite and the other mathematical areas. Earlier Nov 4 minor findings outside the four reproduced cases above should retain their prior-review status until rerun.
