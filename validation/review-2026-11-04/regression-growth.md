# Accuracy review 2026-11-04: `regression-growth`

Mixture regressions and growth mixtures, version 1.1.15, branch `features`. Probes are in `validation/review-2026-11-04/probes/regression-growth-*.R`; every number below came from one of them, run in this session.

**Result in one line.** The likelihoods, parameter counts, scores, standard errors, predictions and empirical Bayes effects all checked out. The defects are in what surrounds them: models that are not identified but report success, unused factor levels counted as parameters, a column name that corrupts tables, and a few labels.

## 1. What was reviewed

| File (R/) | Depth |
|---|---|
| mixture-regression-engine.R, -starts.R, -inference.R, -methods.R, -ordinal.R, mixture-regression.R | Read in full; checked by probes |
| mixture-regression-predict.R | Prediction and simulation read in full; plot helpers read only for the numbers they compute |
| mixture-regression-enumerate.R | Read in full. Grid columns checked through the fit tables. **BLRT was read, not run** |
| growth-mixture.R, -engine.R, -correlated.R, -inference.R, -methods.R, -predict.R, -simulate.R, trajectory-views.R | Read in full; checked by probes (three-step input helpers in trajectory-views.R read only) |
| growth-mixture-parameter-inference.R | Wald path run. **Bootstrap path read, not run** |
| growth-mixture-plot.R | Skimmed for computed quantities (observed means, spread, ellipse). Lines 1-75 not read; nothing rendered |
| kernel-structures.R (regression structures, multinomial M-step, scores, packing) | Read in full, because the scope depends on it; transition and profile structures not checked |

**Not reviewed:** `.latents_em()`, `.latents_quasi_newton()`, `.latents_cluster_bootstrap()`, the entropy helpers beyond one hand ICL, and BLRT calibration. Not run: robust covariance of a multilevel growth fit, OPG, and external comparisons the verification matrix already marks strong (flexmix, Latent GOLD, Mplus, lcmm).

## 2. Findings

Counts: CRITICAL 0, MAJOR 3, MINOR 6, EDGE 3. All are CONFIRMED except where a part is marked PLAUSIBLE.

### MAJOR 1. Models that are not identified report success, with no warning (CONFIRMED)

**Where:** R/mixture-regression.R:318-324 (the only growth identification refusal, a random intercept when every unit has one row), R/mixture-regression.R:655 and R/mixture-regression-inference.R:382 (the singularity rule, smallest eigenvalue at most 1e-10 of the largest).

**(a) Random intercept and slope with two occasions per person.** A pre/post design has 3 distinct covariance moments and 4 covariance parameters, so there is an exact ridge: `g00 - d, g01 + d, g11 - 2d, sigma^2 + d`. On 150 simulated persons (truth `g00 = 1, g01 = 0, g11 = 0, sigma^2 = 0.49`):

| d | LL | g00 | g01 | g11 | sigma^2 |
|---|---|---|---|---|---|
| -0.10 | -454.6319 | 1.097 | -0.209 | 0.451 | 0.390 |
| 0 (reported) | -454.6319 | 0.997 | -0.109 | 0.251 | 0.490 |
| 0.03 | -454.6319 | 0.967 | -0.079 | 0.191 | 0.520 |

- The fit reports `converged = TRUE`, `random_boundary = FALSE`, `n_parameters = 6`, and raises no warning.
- The numerically differenced information has a smallest correlation-form eigenvalue of 7.5e-8. That is differencing noise, but it clears the 1e-10 rule, so the standard errors come out finite.
- The warning depends on the path: with `max_iter = 200, tol = 1e-5` on the same data scaled by 1.0000001, `latents_singular_information` did fire.
- Two classes on 200 persons (`equal` covariance, `equal` variance, truth slope SD 0, residual SD 0.7): slope SD 0.54 (SE 476), residual SD 0.534 (SE 241), `converged = TRUE`, `random_boundary = FALSE`, `n_best_replicated = 4`. The class coefficients and shares have usable standard errors.
- The BIC counts the unidentified parameter, which adds `log(n)` per ridge (5.0 at 150 persons).

**(b) Binary outcomes with two rows per unit and no predictors.** `mixture_regression(b ~ 1, n_classes = 2, family = "binomial", id = , class_level = "group")` has 3 parameters and 2 free probabilities. Result on 150 units: no warning, standard errors 19.45 and 13.81, prior SE 3.58. All 6 starts report LL -193.0514 and `n_best_replicated = 6`, so the replication count reads as reassurance on a flat ridge. The row-level refusal does not apply because the nesting is `"group"`.

**Who is affected:** anyone fitting `random = "time"` to two-wave data, the commonest place this arises; and small binary repeated-measures designs. A reader sees `converged = TRUE` and may report variance components the data cannot determine; only the very large standard errors signal it, and print does not show them.

**Smallest fix:** a rank check before fitting. For a growth fit, stack one row per distinct pattern and pair `(j, l)`: `z_j kronecker z_l` plus the `delta_jl` column for `sigma^2`. If the rank is below `q(q+1)/2 + 1`, raise `latents_not_identified` or warn. Do the same for a group-level binary or ordinal mixture without predictors, comparing parameters with distinct response patterns. Tightening the 1e-10 rule alone would also catch (a), but risks flagging fits that are merely ill-conditioned.

**Numbers changed:** only unidentified fits, which would be refused or warned. Owner to decide between a refusal and a warning.

**Probes:** regression-growth-07-more.R (section 4), regression-growth-08-identify.R (A, B).

### MAJOR 2. Unused factor levels in `membership`, `group_membership` and `random` become parameters (CONFIRMED)

**Where:** R/mixture-regression.R:891 (`.mixture_membership_design`, a model frame without `drop.unused.levels`) and R/mixture-regression.R:690 (`random_design`). The regression formula does drop them (R/mixture-regression.R:606-607), so the behaviour is inconsistent within one call.

A factor keeps its level after `subset()` or after `missing = "omit"`. With `membership = ~ grp` and an unused level `zz`, on 60 persons:

| | k | LL | BIC | SEs |
|---|---|---|---|---|
| Unused level kept | 9 | -330.3227 | 697.4945 | all NA, `latents_singular_information` |
| `factor(grp)` first | 8 | -330.3227 | 693.4001 | present |

- The likelihood is identical, while k is one too many per free class and the BIC is `log(60) = 4.09` too high. Grid comparisons shift accordingly.
- The warning blames "an empty class, a boundary solution or separation", which sends the user in the wrong direction.
- `random = ~ 0 + grp` with an unused level gave k = 13, a reported variance of 0.12 for a level no row has, and all intervals NA.

**The zero-weight guard has the same hole.** R/mixture-regression.R:598-604 warns that rows carrying a level only they have are left out because "the model has no parameter for it". For a membership factor the design still has that column:

- k is 9 against 8 for the same model without those rows, and BIC is 633.3644 against 629.5358, with equal LL (-299.4533). This contradicts the owner's zero-weight rule.
- With the default robust covariance the fit stops with `latents_too_few_groups`, "The unit scores are rank deficient".

**Who is affected:** anyone with a subsetted factor in these formulas, a common situation.

**Smallest fix:** `drop.unused.levels = TRUE` in both model frames. Keep the stored `xlevels` for prediction.

**Numbers changed:** yes for affected fits only (k, BIC, and standard errors restored).

**Probes:** regression-growth-05-invariants.R, regression-growth-06-misc.R (E).

### MAJOR 3. A column named `class` as `id` or `cluster` corrupts the assignment tables without a word (CONFIRMED)

**Where:** R/mixture-regression-methods.R:283-284 and :311, R/growth-mixture-methods.R:586-592, :610 and :698-701. Students in classes makes this a natural name in the package's own field.

| Fit | Table | What happens |
|---|---|---|
| LCGA, `id = "class"` | `groups`, `assignments` | The id column is overwritten by the modal class: ids 501, 502, 503 all read 1 |
| Growth, `id = "class"` | `assignments` | Two columns named `class`; the printed "Class" column shows 501, 502, 503 |
| Growth, `id = "class"` | `individual` | Selecting by name takes the id twice (`class`, `class.1`); the modal class is gone |
| Growth, `cluster = "class"` | `assignments` | `a$class` returns the cluster number, and the printed "Class" column shows it |

A column named `row` or `posterior` (R/mixture-regression-methods.R:282-285) or `n_rows` (:309) collides the same way (read, not run). PLAUSIBLE: trajectory plots read `class` from the individual table, so they would mislabel too; settle it by drawing `plot_trajectories()` on such a fit.

**Who is affected:** any user with such a column name. Estimates are right; assignments are wrong or lost.

**Smallest fix:** refuse reserved names at fit time with `latents_bad_argument`, or write the id under a non-colliding name.

**Numbers changed:** no estimate; the affected tables only.

**Probe:** regression-growth-06-misc.R (D).

### MINOR 4. Correlated-residual fits report the replication of a different model

R/growth-mixture-correlated.R:318 (`fit <- uncorrelated`) keeps the starts table and `n_best_replicated` of the model without the process. The process model is fitted once, from the best uncorrelated solution.

- Example: reported LL -581.43, `selected` start -594.47 (exactly the no-process LL), `n_best_replicated = 6`.
- On `growth_scores`, seeds 1-4 all gave LL -2721.343 with `n_best = 5`, selected start -2743.118.
- So the estimates were stable here, but the count says nothing about the reported maximum. PLAUSIBLE: a better correlated maximum may exist from another class placement; settle it by refitting the process from each completed start.

**Fix:** report the count as `NA`, or refit the process from several starts. **Numbers changed:** that column; a multi-start refit could also change estimates. Probes 07 (section 3), 08 (C).

### MINOR 5. Two-level "Groups assigned" counts rows

R/mixture-regression-methods.R:462 labels everything outside `observation` nesting "Groups". In a two-level fit on 40 groups and 480 rows the Classes table prints "Groups assigned" 279 and 201. The recovery table's unit label, from the same line, reads "Groups" for row counts under `by = "class"` (the label follows from the line, not run). **Fix:** make the label depend on nesting. **Numbers changed:** none. Probe 02.

### MINOR 6. `exp_estimate` of an intercept is printed as a ratio

R/mixture-regression-methods.R:71-74 and :470. An intercept's `exp_estimate` is a baseline odds or rate, printed under "Odds ratio" or "Rate ratio" (0.38 and 2.58 in probe 06; 1.29 in probe 01). The help says the column appears "where the link makes it an odds or rate ratio". Likewise `odds_ratio` on membership intercepts (:164; 0.76 in probe 02). **Fix:** blank them, as thresholds already are, or relabel. **Numbers changed:** that column's intercept rows, if blanked.

### MINOR 7. `recovery` accepts the data in another row order and returns a wrong table

R/mixture-regression-methods.R:385-390 checks only the row count, although `.trajectory_check_data()` exists. On separated row-level classes: correct table 154 / 0 / 1 / 145; with `data[order(data$x), ]`, 76 / 78 / 79 / 67, no message. Group-level fits are usually caught by the within-group check. **Fix:** compare the model columns on `kept_rows`. **Numbers changed:** none for correct calls. Probe 06 (F).

### MINOR 8. `fitted` table's `observed` is not the observed value for numeric ordinal outcomes

R/mixture-regression-methods.R:324. Data 0, 1, 2 show as `observed` 1, 2, 3. `fitted` is on the same 1..C scale and the help documents that scoring for class means, so `residual` is right; only the column name misleads. Trajectory views share this (`.trajectory_observed`, read). **Fix:** say so on the page, or rename. **Numbers changed:** none. Probe 07 (section 5).

### MINOR 9. "`bic` on the independent units" is untrue for `id` without group classes

man/get_results.Rd:246 against `.mixture_weighted_sample()`. With `id` and `n_group_classes = 1`, `bic` equals `bic_rows` (480 rows, 40 ids), while the `weights` and `vcov_type` documentation calls the `id` groups the independent units. The audit of 2026-11-01 logged the numeric consequence as an owner decision; this finding is only the wording. **Fix:** the sentence. **Numbers changed:** none. Probe 02.

### EDGE

- A one-class growth fit gives `share_conf_low` and `share_conf_high` NaN (R/growth-mixture-methods.R:237-238, logit of 1). Probe 07.
- Classes are ordered by rows, not units, under `"group"` (R/mixture-regression-starts.R:356): unbalanced data gave `class_1` share 0.33 (20 persons) and `class_2` 0.67 (40). The reference class is then not the largest. No page claims otherwise. Probe 06 (G).
- Single-row groups pass the binary and ordinal refusals under `"group"`. Both warned singular with NA standard errors; the binary one also warned separation and estimated -734. Probe 08 (B).

## 3. Verified sound

| Item | Check | Result |
|---|---|---|
| Gaussian, one class | `lm` ML; hand SE | LL diff 0; coef 4e-16; SE rel 4e-13; `sigma` and SE(`sigma`) equal; AIC, BIC, SABIC equal hand |
| Poisson with offset; binomial counts; Bernoulli | `glm` | LL within 3.4e-13; SE rel at most 2.7e-7; k correct; `exp_estimate` equals exp(coef) |
| Negative binomial, one class | `MASS::glm.nb`; hand NB2 | LL diff 0; alpha = 1/theta to 8 digits. glm.nb SEs differ up to 2.7% because it conditions on theta and uses expected information, which I take as expected and did not verify separately |
| Ordinal, one class | `MASS::polr` | LL equal to 8 decimals; estimates and SEs agree to the printed digits; sign convention right |
| No-intercept Gaussian | `lm` | coef 2e-16 |
| Two-level with both covariates | Hand triple sum | LL 1.1e-13; k = 11; `bic` on groups, `bic_rows`, SABIC, ICL equal hand; SE vs numDeriv 7e-9 |
| Group Poisson with membership | Hand | LL 0; BIC on groups; logit reference class 1 (first column exactly 0); SE 1.7e-8 |
| Two classes: binomial counts, NB | Hand | Equal to all 9 decimals printed |
| Growth, one class | `nlme`, `lme4` (ML) | LL equal to 8 decimals with all three; fixed effects 4e-10; SE ratio 1.00006 at most (I attribute this to observed against conditional information; not separately verified); k = lme4 df; `bic_rows` = lme4 BIC |
| Growth diagonal, slope only | `nlme` | LL equal to 8 decimals |
| Growth two classes: varying, equal, proportional, diagonal, with membership | Hand multivariate normal | LL diff at most 8e-13; k = length(theta) |
| Same four | numDeriv Hessian | SE within 7e-4. LEARNINGS records numDeriv as weak at this scale, so this is an upper bound |
| Random-table delta SEs (variance, SD, correlation) | numDeriv Jacobian with reported vcov | 7e-11 |
| Delta intervals (log, Fisher) | `nlme::intervals` | Agree to about 1e-3 or better; nlme uses its own approximate covariance |
| Common + offset + two classes | Hand | LL equal; SE 1.8e-5; `class_response` exact |
| Multilevel growth | Hand triple sum; Richardson second difference | LL 4.6e-13; SE 5e-6; `bic` on clusters, `bic_persons`, nobs = clusters |
| AR and Brownian residuals | Hand; numDeriv | LL equal; score vs gradient at most 1.4e-7; SE at most 6e-4; EB by hand equal |
| Robust growth covariance | Hand sandwich from per-person gradients, CR1 | 1.1e-10 |
| `predict` response, conditional, multilevel | Hand marginalisation over classes and group classes | Within 2e-14; row alignment correct under shuffled rows |
| Empirical Bayes | `nlme::ranef`; fitted level 1 | 1.5e-7; 1e-7 |
| Trajectory table and band | Hand | 0 |
| Row order, rescaling | y × 1000 + 5e4, t × 12 + 2000 | Coefficients 1e-15; LL shift exactly `-n log 1000`; slope and SE ratio 83.33 |
| Duplicates against integer weights | Growth, group, two-level, observation | Coefficients within 1.2e-6; LL scales by the weight normalisation |
| BH set in the coefficient table | Hand | Non-intercepts only, jointly over classes and common; intercepts and thresholds NA |
| ICL, two-level | Hand joint entropy | Equal to 6 decimals |
| Recovery shares | Hand | Share of assigned class, as documented |

**Quadrature:** there is none. Gaussian random effects integrate in closed form, and other families are refused with `random`, so no quadrature error exists.

**A cross-family difference, not a defect:** `parameter_inference()` on a growth fit adjusts over every tested row including intercepts and membership logits (probe 03), while the regression table adjusts over non-intercepts. Each matches its own help page, but the same model family gets different sets.

## 4. Test gaps

Method: formula mutations applied in memory, running 11 regression and growth test files from a scratchpad copy so no snapshot was touched (probes 09, 09b, 09c). All mutations were active together.

**Limits of this harness.** The copy has no DESCRIPTION, so the files ran under testthat edition 2 and the snapshot test could not run. The baseline also has two "relabelling" failures in test-review-1026-mixreg.R (dimnames), which I attribute to this harness without having verified it under edition 3. A surviving mutation means these 11 files did not catch it; the other test files were only searched by text.

**Caught:** recovery share; growth standard errors (the three "time in days/hours/calendar" tests caught at least one of the random-table, share or residual-SD doublings; which ones was not separated).

**Survived, all active together:**

| Formula | Mutation | Note |
|---|---|---|
| Regression `p_adjusted` set | BH over all rows including intercepts | No regression test asserts it beyond thresholds being NA |
| `exp_estimate` | Set to the estimate | No regression test mentions it |
| Regression `bic_rows` | Set to `bic` | Not mentioned in the regression and growth test files |
| Growth `bic_rows`, `bic_persons` | Set to `bic` | `bic_persons` not mentioned in those files |
| Growth `sabic` | +1 | |
| Regression `sabic` | +1 | Survived here, but the print snapshot pins `SABIC 6313.40`, so the real suite would catch it. Regression pin, not a formula check |
| `group_share`, group-class SE | Inverted, doubled | |
| `sigma_std_error`, `dispersion_std_error`, `prior_std_error` | Doubled | |
| Share and random-table lower limits | Halved, widened | No test mentions `share_conf` |

`equivalence/engine-golden/` fingerprints outputs and would flag a change, but it records what the code gave, so it cannot show a formula is right.

**Also missing:** any test that an unidentified growth design is flagged, any test with an unused membership level, and any test with `id` or `cluster` named `class`.
