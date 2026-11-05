# Review 2026-11-04, area `surface-docs`

Tables, plots, print methods, descriptives, data sets and conditions against what they say, and vignettes, articles, README, NEWS 1.1.15 and DESCRIPTION against what the output shows and what the statistics allow.

Probes are in `validation/review-2026-11-04/probes/surface-docs-*.R`. All were run through the two-slot wrapper. No file outside this report and those probes was changed.

Result: 0 CRITICAL, 2 MAJOR, 14 MINOR, 5 EDGE. The tables and plots of one fit agree everywhere I compared them; the descriptives and all four simulated data sets are correct. Most findings are prose that the output does not support.

## 1. What was reviewed

| File | Depth |
|---|---|
| R/accessors.R, R/methods.R, R/latents-table.R, R/descriptives.R | Read line by line; tables probed against hand computation |
| R/get-data.R | Code of the catalogue read (lines 480-1171); roxygen prose only for `residuals`, `standardized` and recovery, as the help pages were audited twice this week |
| R/wildcards.R | Read in full; every classification verb probed against its table in eight fit types |
| R/plot-gg.R, R/plot.R | Read in full; layer data compared with tables |
| R/plot-named.R | Code read to the start of `plot_diagnostics()`; the rest only through probes |
| R/conditions.R, R/condition-calls.R | Registry compared by grep with raising sites. R/condition-calls.R was not read closely. The inheritance claims were not re-derived (audit 2026-11-01 did) |
| R/data.R, data/, data-raw/ | Each generator re-run with its save step removed and compared with the shipped data at tolerance 0; documented values checked |
| vignettes/lpa.Rmd, lca.Rmd | Read in full |
| vignettes/articles/*.Rmd (18), README.Rmd | Read in full. Hard-coded numbers re-computed by running the relevant calls, not by rendering |
| NEWS.md 1.1.15, DESCRIPTION | Read; claims about the verbs and the grid summary checked by running them |

Also read outside the list because a finding led there: R/compare-models.R (the penalty only), R/simulation-service.R:210-250, R/selection.R:1-200, R/enumeration.R:660-800, R/transitions.R:208-235, R/growth-mixture-plot.R:106-128, R/mixture-regression-methods.R:445-492.

**Not reviewed or not re-run:**

- No article was rendered, and figures were not looked at. Plots were checked through their data and labels.
- The p-values quoted in growth-mixture.Rmd (p = 0.01, p = 0.40, 99 replicates, about an hour) were not re-run.
- `validation/ADDITIVE_SIMULATION.md` and the coverage statements in additive.Rmd taken from it.
- The statement in workflow-lca.Rmd:141 that the categorical residual test rejected 0.25% and 0.08% at nominal 5%. I found no record of it in LEARNINGS.md and did not simulate it, as the brief rules out simulations. It is unverified.
- `get_tna()` and `get_group_tna()` pooling, as described in case-engagement-transitions.Rmd.
- The internal arithmetic of LMR, three-step and inference tables, beyond matching the prose to the output. Those belong to other areas.

## 2. Findings

### MAJOR

#### M1. `compare_models()` compares a staged fit with a joint fit on the wrong parameter count, and gives the staged fit all the BIC weight

CONFIRMED. R/compare-models.R:418 (`n_parameters = as.integer(fit$n_parameters)`, with the criteria taken from the fit's stored values). The file has no handling of `fixed` or `n_parameters_with_measurement`.

**What is wrong.** `?fit_staged` (man/fit_staged.Rd:126-127) says a comparison with a jointly estimated fit "uses `n_parameters_with_measurement`". `compare_models()` uses the 3 parameters of the second stage. The staged fit has the lower likelihood and the same model, yet its BIC is penalised for 3 parameters against the joint fit's 23.

**Evidence** (`surface-docs-wild.R`, `surface-docs-arith.R`), 2 profiles and 2 group classes on `course_engagement`:

| Model | k | LogLik | BIC | ΔBIC | Weight |
|---|---|---|---|---|---|
| joint | 23 | -8439.82 | 16986.89 | 92.74 | 0.000 |
| staged | 3 | -8440.08 | 16894.15 | 0.00 | 1.000 |

With 23 parameters the staged BIC is 16987.42, which is 0.53 above the joint fit, so the ranking reverses. No warning is raised. `compare_models(test = "lrt")` does refuse the same pair (`latents_bad_nesting`), so only the criteria path is affected.

**Who is affected.** Anyone comparing `fit_staged()` with `multilpa()`. case-engagement-covariates.Rmd:326 invites this kind of comparison.

**Smallest fix.** Compare on `n_parameters_with_measurement` when the fits hold different blocks, or refuse with `latents_incomparable_models`. **Changes reported numbers** in that comparison only.

#### M2. `p_adjusted` from `get_results(fit, "residuals")` is not adjusted, and `get_residuals(fit)` returns different numbers under the same column name

CONFIRMED. R/get-data.R:859-865 with R/inference.R:1750 (the default `adjust` resolves to `"none"`); R/report.R:150 and R/wildcards.R:1106-1108 use BH.

**What is wrong.** By default `p_adjusted` equals `p_value`. This is documented in `?get_results` ("`none` by default") and recorded in LEARNINGS as by design, so the code does what its page says. But the column name says otherwise, the two verbs that return this table disagree, and an article uses the unadjusted column as an adjusted filter.

**Evidence** (`surface-docs-misc2.R`), three-class LCA of `student_esm`, 84 tests, and workflow-lca.Rmd:145-149 `subset(get_results(fit, "residuals"), p_adjusted < 0.01)`:

| Source | Pairs with `p_adjusted < 0.01` | `listened_music` x `walking` |
|---|---|---|
| `get_results(fit, "residuals")` | 6 | 2.1986e-05 |
| `get_residuals(fit)` | 1 | 0.0018468 |

On a two-level LPA the first row reads 6.37e-12 against 1.27e-10 (`surface-docs-readme-misc.R`). The article then calls a pair passing this filter "strong evidence of dependence".

**Who is affected.** Every reader of the default residual table. case-engagement-evaluation.Rmd prints the same column as well, though its prose rightly treats the p-values as descriptive.

**Smallest fix.** Either make `adjust = "BH"` the default of the table, or leave `p_adjusted` as `NA` under `adjust = "none"`. Correct workflow-lca.Rmd either way. **Changes a reported column.** The owner should decide, given the 2026-10 decision that the default is unadjusted.

### MINOR

#### m1. case-esm-lca.Rmd has the two student classes the wrong way round

CONFIRMED. vignettes/articles/case-esm-lca.Rmd:171-173 and 183-185.

The prose says class 1 holds 58.4% of students with 0.974 of prompts in the low-media profile, and class 2 holds 41.6% with 0.565 media-heavy. The same call gives group class 1 = 0.416 (0.565 media-heavy) and group class 2 = 0.584 (0.974 low-media) (`surface-docs-esm.R`). NEWS 1.1.15 warns that group classes may now be numbered the other way round; this article's prose was not updated. It is the only hard-coded article with this problem: the other numbers in it, and in the case-engagement articles, workflows, lta, growth-mixture, trajectory-classes, mixture-regression and README, reproduce (section 3). The article has no `stopifnot()` guard. Fix: swap the labels and add a guard. No number changes.

#### m2. Three documents say the bootstrap p-value is withheld if any replicate fails; the code withholds it only below half

CONFIRMED. vignettes/articles/model-selection.Rmd:669-671, vignettes/articles/class-enumeration.Rmd:528-529, README.Rmd:320, against R/simulation-service.R:227-236.

The code requires `2 * n_valid >= iter` and otherwise computes `(1 + b) / (n_valid + 1)` on the valid replicates. With 15 valid of 19, the helper returns p = 0.125 rather than `NA` (`surface-docs-misc2.R`). The code comment explains this as a deliberate choice, and print() and the warning state it correctly, so the prose is what is wrong. Its denominator `B + 1` should also say valid replicates. No number changes.

#### m3. trajectory-classes.Rmd calls the bootstrap test "the stronger evidence" against over-extraction, and it over-extracts here too

CONFIRMED. vignettes/articles/trajectory-classes.Rmd:223-225.

The article's point is that latent class growth analysis without random effects needs extra classes to absorb within-class spread. The bootstrap null has the same misspecification, so the test cannot correct for it. On `growth_scores` (three generating kinds, no random effects) every step is rejected at the smallest p attainable (`surface-docs-lcga-blrt.R`, 19 replicates, 3 starts):

| Comparison | LR | p |
|---|---|---|
| 3 v 4 | 224.03 | 0.05 |
| 4 v 5 | 218.15 | 0.05 |
| 5 v 6 | 227.40 | 0.05 |

Fix: say the test shares this misspecification and that the remedy is to compare with the random-effects model. No number changes.

#### m4. Enumeration plots report candidates never estimated as "did not converge"

CONFIRMED. R/plot-gg.R:1642 and 1683-1684; R/plot-gg.R:1814 and 1827-1828.

`failures <- sum(!grid$converged)` counts candidates refused as not identified. A grid with one such candidate is subtitled "6 candidates, 1 did not converge (x)", while its table says "Multiple group classes are not identifiable" and `plot_trust()` says "not fitted" (`surface-docs-wild.R`). The grid in model-selection.Rmd has three such candidates. `print.multilpa_enumeration` already separates the two cases. Label only.

#### m5. The class table of a two-level regression labels rows as groups

CONFIRMED. R/mixture-regression-methods.R:462 and 476.

`units` is "Rows" only when `nesting` is `"observation"`. For `mixture_regression(score ~ hours, id = "student", n_group_classes = 2)` the class table prints "Groups assigned 557, 343", which sums to the 900 rows; there are 150 students (`surface-docs-wild.R`). The data frame is correct. Print only.

#### m6. A printed coefficient rounds to 0.00 beside an interval in four decimals

CONFIRMED. R/latents-table.R:96-97.

A column whose largest value is at least 100 prints two decimals whatever its smallest value is. With an intercept near 200 and a slope of 0.0040 the table shows (`surface-docs-misc2.R`):

```
Class 1  x   0.00  [ 0.0040,  0.0041]  <.001
Class 2  x  -0.00  [-0.0030, -0.0030]  <.001
```

A predictor in large units (income, for instance) gives this on realistic data. The data frame is correct. Print only.

#### m7. Two wrong statements about AIC's penalty and path

CONFIRMED.

- vignettes/articles/model-selection.Rmd:322 says "AIC, with the lightest penalty". At the 106 students of this example SABIC's per-parameter penalty is 1.504, below AIC's 2, and CLC has none. The article itself gives the 175.34 crossover at lines 244-246. `sabic_groups` also picks the largest candidate.
- vignettes/articles/class-enumeration.Rmd:233-234 says AIC's weak penalty "leaves it decreasing over the entire range". It rises from 5 to 6 profiles at one group class (17077 to 17093) and at two (16735 to 16737) (`surface-docs-grid.R`).

#### m8. The article's account of AvePP against entropy has the mechanism backwards

CONFIRMED by arithmetic. vignettes/articles/model-selection.Rmd:524-527.

It says a level with mostly well-classified units and a tail of uncertain ones has "moderate AvePP together with low entropy" (`surface-docs-arith.R`):

| Posteriors | Relative entropy | Mean largest posterior |
|---|---|---|
| K = 2, 90% at .99, 10% flat | 0.827 | 0.941 |
| K = 4, 80% at .97, 20% at .70 | 0.768 | 0.916 |
| K = 2, everyone at .90 | 0.531 | 0.900 |
| K = 2, everyone at .80 | 0.278 | 0.800 |

A tail leaves relative entropy high. Low entropy beside good AvePP comes from moderate uncertainty spread over all units.

#### m9. "typically far better separated" is not what the article's own grid shows

CONFIRMED. vignettes/articles/class-enumeration.Rmd:499-501.

For the chosen `profiles4_groups3`, relative entropy is 0.823 for groups and 0.739 for profiles. For `profiles3_groups3` it is 0.797 and 0.792. In 3 of the 15 two-level candidates the group level is worse: `profiles2_groups3` (0.816 against 0.940), `profiles2_groups4` (0.798 against 0.942) and `profiles3_groups4` (0.772 against 0.793) (`surface-docs-grid.R`). model-selection.Rmd:533-534 words the same point correctly.

#### m10. `drop_empty = FALSE` does not restore the earlier individual-level BIC as the article says

CONFIRMED. vignettes/articles/model-selection.Rmd:189-190, against `?multilpa` (R/multilpa.R:991-997, which is right).

With 150 empty rows, `bic_individual` is 14967 under both settings, the value for n = 1272; n = 1422 would give 14971 (`surface-docs-srl.R`). Entropy and counts do change (0.76988 against 0.70763). Whether this reproduces 0.9.37 was not tested; the claim as written cannot hold for the criteria the preceding sentences list.

#### m11. Standardized means and descriptives ignore sampling weights

CONFIRMED, previously logged (audit 2026-11-01, "owner decision"); new evidence here. R/accessors.R:620-628, also used by R/plot-gg.R:263 and the heatmap.

The centre and spread are unweighted and include zero-weight rows. With weights of 3 and 1 by student type, the model-implied overall standardized mean, which the heatmap presents as "SDs from observed mean", is -0.13 to -0.22 rather than about 0 (`surface-docs-readme-misc.R`). `descriptives()` of a weighted fit returns the unweighted table, and neither page says so. Fixing the centre **changes reported numbers** on weighted fits; stating it on the page does not.

#### m12. vignettes/lca.Rmd shows a call that fails on its own fit

CONFIRMED. vignettes/lca.Rmd:106-112.

The `eval = FALSE` block runs `plot_profiles(lca)` on a fit with "no"/"yes" categories, which raises `latents_bad_argument`: "Categorical profile summaries need distinct finite numeric category labels" (`surface-docs-readme-misc.R`).

#### m13. "the grid above raises once"

CONFIRMED. vignettes/articles/model-selection.Rmd:126. The grid raises `latents_thin_profiles` twice, for profile masses of 4.9 and 9.1 (`surface-docs-grid.R`).

#### m14. "a residual above 5" names the wrong column

CONFIRMED. vignettes/articles/workflow-lpa.Rmd:226. The filter is `statistic > 5`, a z statistic. The residuals in those rows are 0.38 to 0.69 (`surface-docs-workflow.R`).

### EDGE

- `get_results(fit, "assignments", truth = )` with a group-level truth on a fit with `n_group_classes = 1` cross-tabulates against the single group class (62 and 44 in class 1, proportion 1) rather than the profiles. Not wrong, but says nothing (`surface-docs-misc2.R`).
- case-engagement-evaluation.Rmd says 0.226 matches "to six decimal places"; the values are 0.2256397 and 0.2256388, which differ in the sixth (`surface-docs-case2.R`).
- enumeration.Rmd says the result "holds 24 fitted candidates" and then that three could not be estimated.
- R/latents-table.R:66: a predictor named `all` or `common` would print as "All classes". From reading the code; not run.
- `latents_partial_pair` stays on the conditions page as "Not raised" (already logged).

## 3. Verified sound

**Descriptives.**
- ICC equals `(MSB - MSW) / (MSB + (n0 - 1) MSW)` with `n0 = (n - sum(n_i^2)/n) / (k - 1)`: 0.52413 by `anova(lm())` and by the package, on unbalanced groups.
- Stratified output, including the `NA` stratum, adds up.
- The formula matches the prose in lpa.Rmd and case-engagement-profiles.Rmd.

**Tables.**
- Profile table standard errors match `parameter_inference()` for all 14 covariance structures, with none missing for means or variances.
- Standardized means and their standard errors equal a hand computation from `get_results(fit, "data")`.
- `delta_profile` equals a hand computation from the candidates table.
- `get_tree()` masses, `share_of_from` and `share_of_to` equal `crossprod(p, q) / n` by hand.
- The `odds_correct_classification` formula checks by hand on one row (11.45).
- LTA `prevalence` is the posterior expected share of a class's enrolments, as case-engagement-transitions.Rmd says (read in R/transitions.R:220).

**Table against plot.**
- `plot_bars()`, `plot_heatmap()`, `plot_prevalence()` and `plot_transitions()` match their tables.
- The minimum ring matches `get_best()`, and `plot_deltas()` matches `get_deltas()`.
- `plot_sizes()`, `plot_posteriors()`, `plot_entropy()` and `plot_avepp()` match `get_classification()`, `get_entropy()` and `get_avepp()`, at each level, on a two-level LPA, single-level and two-level regression, multilevel growth mixture, additive, full cross-level, general LTA and a weighted two-level regression.
- The 80% band in growth plots is `qnorm(0.9)` on the random-effect variance alone, as growth-mixture.Rmd says.

**Formulas in the articles** (derived, not run unless noted):
- AIC, KIC, BIC, SABIC, CAIC, ICL, AWE and CLC as tabled in model-selection.Rmd; the 175.34 SABIC crossover; relative entropy; AvePP; OCC; BLRT `(b + 1) / (B + 1)`.
- Knee geometry in R/selection.R against the numbered steps of the article.
- Median and tied-mode rules in lca.Rmd against R/plot-gg.R:405-408.

**Hard-coded numbers reproduced by running the call:**
- lpa-related: README (83%, 79%, entropy 0.94).
- case-engagement-profiles, -evaluation, -covariates and -transitions: every number listed in those articles that I compared, including class numbering.
- case-esm-lca: all except m1.
- workflow-lpa (BIC pattern, 4.7, p = 0.005, 69% and 31%, 0.74, 0.88) and workflow-lca (13.2, 1.3%, 58%, 4%, 39%, 0.66, 5.30, 4.63, 4.53).
- class-enumeration SRL section (seven profiles on the boundary, masses 3 and 10, eight profiles with two nearly empty profiles).
- growth-mixture (weights, class labels, 4, 0.7, -0.25, programme sign, 40 units).
- trajectory-classes (BIC falling to five, curvature not rewarded).
- mixture-regression (sleep and motivation signs, classes).
- additive (38 of 49 effective groups, recovery prose, boundary fit).
- lta (odds ratios near one).

**Claims I tried to refute and could not** (`surface-docs-misc3.R`):
- model-selection.Rmd:217-218, deleting empty rows gives identical results: at `tol = 1e-10` the log likelihood (-7358.3), `bic_individual` (14966.8) and relative entropies (0.76992, 0.92621) agree. At the default `tol` the entropies differ in the fourth decimal (0.76988 against 0.76997), which is convergence, not a defect.
- covariates.Rmd:181-183, "the profiles do not owe their separation to previous grades": the attendance effect is 0.036 (0.0014 to 0.071), forum posting 0.004, BIC 16986.89 without against 16992.07 with, and the profile means barely move (attendance 0.649 to 0.642, -0.921 to -0.910). The claim covers only the two indicators tested, but the output supports it.
- tree-plot.Rmd's "offline" class: in the two-class solution it has elevated friends (0.214 against 0.144 overall) and walking (0.195 against 0.127) and almost no media use, so the label holds.

**Data sets.** All four simulated generators reproduce the shipped data exactly. The generating values in R/data.R match the scripts, and stationary shares follow from the transition matrices (0.8 and 0.25). Dimensions, ranges, 49.7% `worried == 1`, 15 to 51 prompts per student (median 24) and the `srl` correlations match. `student_esm` and `srl` were checked against their documentation only, since their sources are not in this repository.

**Conditions.** Every class documented on latents-conditions has a raising site by grep, apart from `latents_partial_pair` (documented as not raised). The only `latents_` class in the code that is absent from the page and is not an object or attribute class is `latents_task_error`, which is internal to the parallel collector. No class name is built dynamically.

**NEWS and DESCRIPTION.**
- Grid summary warnings, per-family classification plots, refusal of `T1 ... T10` labels and the emEM default for EII, EEI and EEE (printed "emEM: 50 screened") all hold.
- The version is 1.1.15 in both files.
- `%||%` is defined in R/utils.R, so `R (>= 4.1.0)` is consistent with its use.
- Sampling-weight statements in lpa.Rmd hold: Kish n is printed and `vcov_type = "observed"` is refused.

## 4. Test gaps

- **M1:** no test compares a staged fit with a joint fit through `compare_models()`. The grep for `staged` or `fit_staged` in tests/testthat/test-compare\*.R found nothing.
- **M2:** tests/testthat/test-wildcards.R:130 checks `get_residuals()` against `diagnostics()` only, so it would not notice the disagreement with `get_results()`.
- **Article prose:** `stopifnot()` guards exist in several articles but not in case-esm-lca.Rmd, the workflows' hard-coded numbers or the case-engagement articles, so a change in class numbering does not stop a build. That is how m1 survived. A guard on each named class would have caught it.
- **Print layer:** I found no test of the "Groups assigned" label on two-level fits (m5), of decimals in mixed-magnitude columns (m6), or of the failure wording in enumeration plots (m4).
- **Weights:** the standardized table has no test on a weighted fit (m11).
- **Sound:** the ICC has a hand-computed test (tests/testthat/test-audit-reporting-repairs.R:10), which would fail if the formula broke.
