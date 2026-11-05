# Review 2026-11-04: `additive-crosslevel`

Package at branch `features`, version 1.1.15. R 4.5.2, one core, single-threaded BLAS, no parallel workers started. Probes are in `validation/review-2026-11-04/probes/additive-crosslevel-p*.R`. No file in R/, man/, tests/ or vignettes/ was changed.

Summary: the additive likelihood, variance components, conditional intercepts and intervals are correct, and so is the cross-level likelihood as a computation. The problems are in what the cross-level numbers are taken to mean. Findings: 0 CRITICAL, 2 MAJOR, 6 MINOR, 3 EDGE.

## 1. What was reviewed

| File | Depth |
|---|---|
| R/additive.R | Every function read. Density, E-step and fit checked by quadrature, dense covariance, `nlme` and BFGS. |
| R/additive-inference.R | Every function read. Effective groups and intervals recomputed by hand. Scores and information not redone: the 2026-11-01 audit checked them against numDeriv on all five specifications. |
| R/additive-methods.R | Tables and accessors read. Parameters, group classes and conditional intercepts recomputed. Plot bodies read only for the numbers they use. |
| R/additive-compare.R | Read. Simulator tested only by mutation. Grid run with a mocked lost worker. |
| R/cross-level.R | Every function read. Likelihood, counts, criteria and compositions recomputed. Stationarity by finite differences, weights against duplicated groups, zero weights against a subset. |
| R/parallel.R, R/parallel-report.R | Read line by line. Multi-core paths judged by reading. Retry and loss paths run by mocking `.latents_fork_pass` and `.latents_socket_pass` in one process. |
| R/utils.R | Read. `.multilpa_select_start()` run on constructed ties. |

Also read: HANDOFF.md, LEARNINGS.md, validation/ADDITIVE_MODEL_DESIGN.md, validation/ADDITIVE_SIMULATION.md, validation/additive-diagnostic-protocol.md, validation/audit-2026-11-01/additive-compare.md, equivalence/latentgold-families/REPORT.md and refine.csv.

Not reviewed:
- No worker process was started, so real fork and socket behaviour rests on reading the code and on the package's own tests.
- The Windows `tasklist` path and the socket-worker retirement logic were read but not exercised.
- `search = "emEM"` was not probed on these families.
- The vignette was not rendered, and the 50 threshold was not re-estimated (settled study, not re-run).
- The Latent GOLD gaps were not re-run. refine.csv still shows 0.054 in log likelihood on g08 after refining the quadrature, against a declared 1e-4. The report and LEARNINGS attribute this to Latent GOLD's quadrature. I did not check it.
- Cross-level families have no external reference, so correctness here means agreement with the formula in R/cross-level.R:7, computed independently.

## 2. Findings

### F1. MAJOR, CONFIRMED: cross-level criteria select group classes and the full form on data with no group structure

Where: R/cross-level.R:437-440 (criteria), R/multilpa.R:1136-1137 and R/get-data.R (the "comparable" statement), vignettes/articles/additive.Rmd:280-284 (the reading of the likelihood gain).

The statement "comparable only between cross-level fits of the same data" is true of the scale (see section 3). It does not make the comparison informative, and the pages and the vignette read as if it did. `fit$bic` and `bic_individual` apply ordinary penalties to a likelihood that counts the ratings twice.

Evidence: every member draws her profile independently of her group, so the true number of group classes is 1 and composition does not differ between groups. There are 2 true profiles. Three seeds per row, `n_starts = 3`, probe p5. Negative values mean the larger model wins.

| Scenario | restricted: BIC(H=2) - BIC(H=1) | full: BIC(H=2) - BIC(H=1) | full - restricted at H=2: log L (1 extra parameter) | full - restricted at H=2: BIC | profile family (ordinary likelihood): BIC(H=2) - BIC(H=1) |
|---|---|---|---|---|---|
| d=2, n=6, J=100, equal/equal | -46.5, -60.4, -41.0 | -104.2, -112.0, -84.4 | +31.2, +28.1, +24.0 | -57.7, -51.5, -43.3 | +8.6, +9.2, +9.2 |
| d=1, n=6, J=100, equal/equal | +7.8, +1.4, +4.3 | -47.0, -45.7, -36.6 | +29.7, +25.8, +22.8 | -54.8, -47.0, -40.9 | +8.4, +9.2, +9.2 |
| d=2, n=6, J=100, varying/varying | -37.3, -52.9, -34.5 | -95.3, -104.6, -78.4 | +31.3, +28.2, +24.3 | -58.0, -51.7, -44.0 | +8.6, +9.2, +9.2 |
| d=2, n=20, J=60, equal/equal | -24.0, -38.7, -34.2 | -50.7, -76.8, -55.6 | +15.4, +21.1, +12.7 | -26.7, -38.1, -21.4 | +7.8, +8.1, +8.2 |

Probe p4 gave the same result on six more seeds of the first scenario: full -66.6 to -113.1, restricted -26.5 to -48.7, profile family +0.9 to +9.2. A spurious full fit reports compositions of about 0.33 and 0.62 (seed 1) when the true share is 0.5 in every group.

Mechanism:
- A group mean is a function of its members' realised profiles. A split of the groups by their means therefore also separates them by composition, and the full model is credited for the same split twice. This holds in one dimension too (d=1 row).
- With d >= 2, composition makes group means correlated across indicators while the group-level density is diagonal. The restricted model then also gains, which the d=1 row does not show.

Who is affected: anyone choosing H, or choosing full over restricted, from the `fit` table. The vignette does this: "letting the composition differ by class raises the likelihood far more than its one extra parameter costs". The table shows a gain of 13 to 31 log-likelihood units with one extra parameter when composition does not differ at all, so that sentence is not evidence. These are 12 + 6 datasets, not a study, so rates are not estimated. All 18 agree in direction for the full form. ADDITIVE_SIMULATION.md lists simulation evidence for these families as not established.

The arithmetic is correct and this follows from the manifest specification, as I understand it from the repo's notes on Houle et al. (I did not re-read the paper). So this is a defect in what the numbers are said to mean, not in the estimator.

Smallest fix: state on the `multilpa` and `get_results` pages and in the vignette that these criteria do not calibrate H or full against restricted, and remove the vignette's inference. Reported numbers do not change. Whether to withhold or rename the criteria is the owner's decision, and that would change the `fit` table.

### F2. MAJOR, CONFIRMED: `prevalence` in the cross-level `profiles` table measures different things in the two forms

Where: R/cross-level.R:452-454 (computed), :520 (comment: "a profile's share of the individuals"), :534 (reported).

Restricted reports the share of individuals. Full reports `sum_h omega_h pi_{k|h}`, the share in the average group. They agree only when group size is unrelated to class. No help page defines the column for this family; the only statement of its meaning is the comment at :520, which is false for the full form.

Evidence, same profile on the same data (probes p2, p4):

| Data | full `prevalence` | restricted `prevalence` | full: mean posterior over individuals |
|---|---|---|---|
| class-1 groups n=5, class-2 groups n=15, J=80 | 0.504 / 0.496 | 0.379 / 0.621 | 0.379 / 0.621 |
| class-1 groups n=4, class-2 groups n=40, J=60 | 0.479 / 0.521 | 0.700 / 0.300 | 0.699 / 0.301 |

Who is affected: full cross-level users whose group size is related to composition (large and small schools, for example). The full value is a correct estimate of its own quantity, but read as labelled it is off by 0.12 to 0.22 here. The test at tests/testthat/test-review-1015-families.R:171 pins this formula on equal-sized groups, where the two meanings cannot be told apart.

Smallest fix, for the owner to choose: (a) keep the number and document the column as "average group's share" for full, and correct the comment; no number changes. (b) Report the weighted mean of `subject_posteriors`, as restricted effectively does; this changes the full form's column when sizes are unequal, and the test.

### F3. MINOR, CONFIRMED: the help page misdefines `effective_groups`

Where: R/get-data.R:252-253, man/get_results.Rd:267.

The page says "the number of perfectly classified groups that would estimate the class weight as precisely". That number is `w(1-w)/Var(w)`. The statistic is `w^2(1-w)/Var(w)`, which is `w` times it: the groups expected in this class among those perfectly classified groups. On one fit (probe p1), `effective_groups` is 57.1 and 50.9 while the page's quantity is 108.06 for both classes. The internal comment at R/additive-inference.R:480-490 and the protocol have it right, and the statistic itself is correct.

Fix: reword to "the number of groups in this class, if every group's class were known, that would estimate the weight as precisely". No number changes.

### F4. MINOR, CONFIRMED: cross-level `fit` table and print count zero-weight units

Where: R/cross-level.R:578 (table), :600 (print).

This breaks the owner's zero-weight rule. With 5 of 60 groups at weight zero (probe p2), the table and print give `n_groups = 60, n_obs = 323`. The same fit on the subset gives 55 and 303, with identical log likelihood (-799.6237) and BIC (1675.387). The additive family's table gives 55 and 303 on these data (R/additive-methods.R uses `n_used_*`), and `nobs()` already gives 55. Both cross-level forms are affected. Criteria are correct.

Fix: use `x$n_used_groups` and `x$n_used_obs`, with `%||%` for older fits. Changes these two columns and the print line, only on fits with zero weights. R/cross-level.R:217 has the same gap; see E1.

### F5. MINOR, CONFIRMED by reading plus a one-core probe: a family grid with `seed = NULL` and `cores > 1` loses `set.seed()` reproducibility

Where: R/additive-compare.R:239-242 and :276, with the claim at R/parallel.R:3-6 and :118.

`.additive_enumerate()` runs all of `multilpa(seed = NULL)` on the worker, including the random starts (`.additive_initialize()`, `sample.int()`, R/additive.R).
- Socket worker: a fresh process with no `.Random.seed`, so it seeds itself.
- Fork: `mclapply()` with its default `mc.set.seed = TRUE`. From my memory of `parallel`, under Mersenne-Twister this removes the child's seed, so the child reseeds from time and PID. With L'Ecuyer it gives a stream that is repeatable but differs from one core. I did not run this.

In both cases `set.seed(1); enumerate_classes(family = ..., cores = 2)` differs from one core and, except possibly under L'Ecuyer, from run to run, with no warning. The profile and transition grids avoid this by preparing unseeded starts in the parent (R/start-pool.R:63-92).

The help page limits its promise to "With `seed`", so the page is accurate. The code comment "changes the time an enumeration takes and nothing it returns" is false for this path.

One-core evidence that the RNG state matters (probe p3, H=3 and H=4, `n_starts = 2`): `set.seed(1)` gives -951.6748, -950.7681 and repeats exactly; `set.seed(2)` gives -951.6748, -950.2929; `set.seed(3)` gives -951.6246, -950.7681; `set.seed(4)` gives -951.6748, -950.4966. Differences are small here and would shrink with more starts.

Who is affected: users who set one seed at the top of a script and pass `cores` to a family grid. No test covers this: every family-grid parallel test passes `seed`.

Fix: draw starts in the parent, or draw one seed per candidate from the caller's stream before dispatch. The second changes unseeded one-core results, so it needs the owner's decision.

### F6. MINOR, CONFIRMED by mock: a sensitivity seed whose worker is lost twice crashes with an unrelated error

Where: R/sensitivity.R:171, R/parallel.R:487.

Without `on_lost` the mapper's terminal value is `list(error = ...)`, which is not a condition, so `inherits(fit, "condition")` is false and `vapply()` fails. Probe p4 gives `simpleError: values must be length 1, but FUN(X[[2]]) result is length 0`, instead of the documented NA row and `latents_sensitivity_dropped`. The numbers are not changed silently: the call fails, after `latents_worker_lost`.

Fix: pass `on_lost = function(seed, message) errorCondition(message, class = "latents_candidate_failed", call = NULL)`, as R/start-pool.R does. Changes no successful result.

### F7. MINOR, CONFIRMED: `latents_weak_class` is misdescribed

Where: R/conditions.R:368-372, man/latents-conditions.Rd:389.

The page says `multilpa(family = "additive")` only, and fewer than 50 effective groups. The code also raises it for dispersion fits, and raises it at fit time for collapsed classes (R/additive.R, `.additive_fit()`), with no effective-groups count involved. Probe p4 raises it from a dispersion fit.

Related: the warning from `summary(fit, vcov_type = "robust")` sends the reader to `get_results(x, "group_classes")`, which uses the default covariance (R/additive-methods.R:57). In probe p3 the warning says 37.4 and the table says 36.5. Near 50 the two can disagree about whether the class is flagged. The code comment acknowledges this for the default table only.

Fix: reword the condition's entry and say which covariance the table uses. No number changes.

### F8. MINOR, CONFIRMED: nothing supports or guards the only comparison the page allows

Where: R/diagnostics-checks.R:417, R/compare-models.R, R/cross-level.R:620-623.

- `diagnostics()` on a cross-level fit says "meaningful only against other models: compare_models()", but `compare_models()` refuses cross-level fits (`latents_bad_argument`). `bootstrap_lrt()` and `enumerate_classes()` refuse them too.
- The `compare_models()` refusal says it accepts fits from `multilpa()`, which a cross-level fit is. The additive family got its own message for this reason (R/compare-models.R:171-172); cross-level did not.
- The base generics, which the package does not guard, mix families silently. `AIC(full, additive, profiles)` prints 1766.9, 2113.3 and 1579.2 with no warning (probe p2), though the page calls these incomparable.

Fix: correct the diagnostic detail and the refusal, and consider a marker on the `logLik` object. No fitted number changes.

### Edge cases (not chased)

- E1. R/cross-level.R:217 checks `n_group_classes` against all groups, not those of positive weight. The fit still refuses later, with a message about indicator rows (probe p3). The additive check is correct.
- E2. PLAUSIBLE. In `.additive_polish()`, releasing a coordinate under `additive_dispersion` with equal between variances sets `0.01 * within` per class, so an "equal" block briefly differs across classes before the first M-step pools it. A first likelihood outside the constraint set could read as a decrease in `.latents_em()`. Settle by running that path with a released coordinate. Not run.
- E3. PLAUSIBLE. `.latents_check_cores(1e10)`: `as.integer()` gives NA, so the check fails with a bare missing-value error, not a classed one. Read, not run.

## 3. Verified sound

**Additive likelihood and moments**
- Tiny case (3 groups of sizes 1, 2, 4; two indicators; class-specific W and T): `integrate()` over the intercept against the closed form differs by at most 1.8e-15 in log density. Log likelihood -18.69404 both ways. Conditional mean 2.117647 and variance 0.08823529 match quadrature (p1).
- Fitted log likelihood -2673.204 matches a dense Cholesky evaluation from raw ratings (p1).
- Conditional-intercept table row (t7, class 2, y2) matches the formula after the canonical reorder (mean 0.2929523, variance 0.2211846).

**Variance components**
- One class on unbalanced data against `nlme::lme(method = "ML")`, independently per indicator: log likelihood -1191.739 both; means, T (0.6608181, 0.2974423) and W (1.045177, 0.4257031) agree to the printed digits. SE of the mean 0.10493 / 0.06962 against nlme's 0.10493 / 0.06959. I did not investigate the small gap on the second indicator (p3).
- All five specifications, two classes: an independent base-R likelihood equals the reported value at the estimate, and BFGS from there finds nothing higher to the printed digits (p3).

**Class weights, parameter counts, sample sizes**
- `count = J * weight` at the optimum (79.286 = 150 * 0.5285733).
- Additive counts 9, 11, 9, 13, 11 equal `length(coef)` on the free scale. Cross-level counts 19 (full) and 18 (restricted) match hand counts.
- BIC uses `log(J)` and `bic_individual` uses `log(N)` in both families, by hand. The weighted additive penalty uses `log(150)`, the unweighted count of groups.

**Effective groups**
- `w^2(1-w)/Var(w)` by hand: 57.11601, 50.94092 against the table.
- Unit-free: identical after multiplying the ratings by 10, while the log likelihood shifts by exactly `N d log 10` (3785.45).

**Intervals**
- Logit and log bounds recomputed by hand match to the printed digits.
- Log-scale SE 0.262173 equals delta SE over estimate.

**Weights and invariants, additive**
- Weight 2 against duplicated groups: parameters agree within 4.1e-11; log likelihood scales by exactly J+/sum(w) (-2609.382 * 200/150 = -3479.176).
- Label permutation leaves the likelihood unchanged (difference 0).
- `bootstrap_lrt()` refuses a weighted family fit (`latents_unsupported_weights`), so unweighted replicates cannot meet a weighted statistic.

**Cross-level**
- Independent raw-rating likelihood equals reported: -864.4304 (full), -920.2713 (restricted).
- Full is stationary: largest finite-difference gradient 8.1e-4 over 19 free parameters.
- H=1 gives the same fit in both forms (-959.4913, 13 parameters).
- Restricted composition equals the posterior-weighted formula by hand (difference 0); full composition is the model's `pi_{k|h}`; rows sum to 1. The plot subtitle and help page describe the two correctly.
- Weights against duplicates, both forms: profile and class parameters within 3.6e-5, log likelihood scales exactly (-865.1979 * 80/60 = -1153.597).
- Zero weights equal the subset fit: full -799.6237, restricted -853.2859.

**The page's scale statement**
- Multiplying units by 10 shifts the working log likelihood by `(N + J) d log 10` (1763.78 observed and predicted), against `N d log 10` for additive. So differences among cross-level fits of the same data are unit-free, while differences against other families depend on the units. That part of the statement holds; F1 is about what the comparison can support.

**Seeds, workers and ties** (by reading unless noted)
- Seeded candidates give the same RNG path on any process: the parent's RNG kinds, including `sample.kind`, are carried into the task, and `set.seed()` runs inside it.
- Start values are drawn in the parent in start order, and EM draws no random numbers. Results return in `indices` order whatever the dispatch or completion order (`order(positions)`, R/parallel.R:493).
- A lost seeded task retried by mock gives values identical to serial (1.372198, 2.571629, 3.215942, 4.364452).
- A failed worker cannot silently change a result. A loss raises `latents_worker_lost` and leaves a recorded failure: an `error` row in the family grid (run by mock, and `get_best()` then leaves that candidate out) or a failed start with `latents_failed_starts`. The reporter sums losses across phases. A candidate's own R error is caught separately and propagates. F6 is a loud failure, not a silent one.
- `.multilpa_select_start()` is deterministic, independent of arrival order: highest within `tie`, then converged, then earliest. Constructed ties gave 1, 4 and 4 as derived; `-Inf` and `NaN` are never candidates.

## 4. Test gaps

Each mutation below was applied inside a probe in memory, against the listed files only (p6). The rest of the suite was not run against them, so "passes" means these files.

| Mutation | Files | Result |
|---|---|---|
| `.cross_level_count`: full counts prevalences once instead of H times | test-cross-level, test-review-1015-families, test-review-followup, test-review-families | **All pass.** No test hand-counts cross-level parameters, so the BIC penalty is unpinned. |
| `.additive_simulate`: between SD set to the variance (missing `sqrt`) | test-additive-compare, test-audit-1101-additive-compare, test-review-1026-additive | **All pass.** The simulator's test checks only the mean (test-additive-compare.R:83-87), so the bootstrap reference distribution could be wrong unnoticed. |
| `.multilpa_wald_bounds`: log-interval spread doubled | test-additive-methods, test-review-1015-families | **All pass.** Only `conf_low > 0` and equality with `confint()` are checked. No test file calls the helper by name. |
| `.multilpa_wald_bounds`: logit lower bound replaced by Wald | test-additive-methods, test-review-1015-families | **All pass.** |
| Restricted composition replaced by a pooled share | the cross-level files above | Caught (test-review-followup, 1 failure). |
| Full prevalence replaced by an unweighted class average | the cross-level files above | Caught (test-review-1015-families, 1 failure). |

Other gaps, by reading:
- The header of tests/testthat/test-cross-level.R says the reference likelihood is computed from the published formula, but that file contains no such test. The real one is in test-review-families.R (`review_cross_loglik`), so coverage exists and the header is misplaced.
- No test covers a family grid with `seed = NULL` and `cores > 1` (F5).
- No test covers terminal worker loss in `sensitivity()` (F6).
- No test covers `n_groups` and `n_obs` in the cross-level `fit` table under zero weights (F4).
- No test covers `prevalence` on unequal group sizes (F2).
- Nothing tests cross-level criteria on data without group structure (F1).
