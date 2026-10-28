# Astra — changes since release and pre-release review

**Review date:** 27 October 2026  
**Release baseline:** `24773e6954bb7723aa75a34484a8e4638c32f1f6`, version 0.97.1  
**Reviewed branch/HEAD:** `features`, `5596a15c932a306bf7d2e04a8aca1b365809c13f`, version 0.97.4  
**Candidate:** the working tree at review start, including the uncommitted diagnostics, validation and presentation changes described in the existing handoff as intended for 0.97.5. `DESCRIPTION` still says 0.97.4.

## Repair follow-up — 27 October 2026

**Current status: all 17 findings are repaired and locally verified.** The user subsequently authorized repairs to **all 17 findings**. The implementation changes below preserve the pre-existing working changes and the package's statistical conventions. `DESCRIPTION` remains 0.97.4; no commit, push or release has been made.

The numbered review beginning below is retained as a **historical account of the pre-fix candidate**. Its hold recommendation and failing results describe that snapshot. Current repair evidence is retained separately in [validation/astra-fixes-2026-10-27](validation/astra-fixes-2026-10-27/).

| Finding | Repair |
|---|---|
| A01 | Updated enumeration vignette/article winner assertions and explanations to the existing complete-data entropy convention; regenerated README output. Assertions remain active. |
| A02 | LTA covariate centering/scaling now uses the same contributing initial/transition rows as likelihood and identification checks, in fitting and inference. |
| A03 | Gaussian degeneracy checks account for offsets and intercepts. Valid constant-response models match independent `lm()` estimates and likelihoods; genuinely degenerate models remain rejected. |
| A04 | Forwarded arguments receive collision-free bindings, including collisions with caller symbols and the fitting function. |
| A05 | `multilca()` forwards its already-evaluated arguments; data generators and argument expressions execute once. |
| A06 | LTA probability-boundary checks exclude inactive initial/transition probabilities while retaining checks on active boundaries. |
| A07 | Ordinary categorical support is derived from positive-weight observations. Original excluded labels are preserved for data retrieval; declared joined-item grids retain their full combinations. |
| A08 | Bootstrap covariance rows and columns are NA for coefficients marked unidentified, consistently with the coefficient table. |
| A09 | Classification readers share eligible-unit masks and weights, including cross-level groups and general-LTA empty/excluded observations. Cross-level assignment counts also exclude zero-weight groups. |
| A10 | Regression nesting is read from its specification; group-level regressions and hierarchical growth report the levels they actually classify. |
| A11 | Sensitivity plots distinguish converged, unfinished and failed starts; an unfinished high likelihood is no longer presented as a replicated maximum. The deterministic Ward start policy is unchanged. |
| A12 | Replication log likelihood is checked against independent reference likelihoods and transformed under affine/duplication invariants. Eight search-bookkeeping fields have explicit, reasoned exemptions. |
| A13 | The shared gate rejects partial/absent coverage, unexplained disagreements, failed suites and stale known-open entries; strict mode rejects skipped suites. All 27 formerly partial fields are now fully compared using independent mixed-residual and excluded-observation references. |
| A14 | CI installs the reference dependencies and current package, runs harness self-tests and every registered suite in strict mode, and retains newly generated comparison evidence. |
| A15 | Sensitivity-method dots arguments are documented and Rd documentation regenerated. |
| A16 | Newly supported diagnostic families now supply entropy, classification and average-posterior tables; data-frame conversion respects the requested table. |
| A17 | Regression-grid classification quality uses eligible units and sampling weights consistently with each fitted model. |

The additional extreme-growth starting-value limitation reported in the review is also repaired: uncorrelated growth initialization and optimization use contributing persons, then restore expectations for every original person. Appending zero-weight outcomes of magnitude 1e100 now leaves estimates, likelihood and the start table unchanged.

Two packaging issues found during integration were also corrected: the Markdown NEWS version heading now parses correctly, and the existing `lme4` reference test was moved intact to the external equivalence suite. Its four assertions pass there; the package gains no dependency.

### Repair validation

| Check | Result and scope |
|---|---|
| Package unit suite | **213 files, 1,846 tests, 12,664 passing assertions; zero failures/errors.** Two model-diagnostic warnings from fixtures; 16 fork tests skipped because Accelerate BLAS cannot safely fork. Socket-worker tests ran. This complete run preceded the final empty-data ordering correction described below. |
| Final installed-package regressions | **33 tests, 333 assertions; zero failures/errors/warnings.** Reran the core regressions, mixture regression and error-contract tests against the final built package after the empty-data correction. |
| Full-output equivalence | **All 15 suites represented; 9,068/9,068 quantities agree over 3,143,431 values. Zero partial or uncovered outputs; zero skipped suites.** The final matrix combines the full run with four affected-suite reruns, as explained below. Four conservative Gaussian-support warnings arise in the hard-mclust fixtures. |
| Harness self-tests | **31 tests, 138 assertions pass**, including partial-coverage rejection, independent replication-likelihood checks and 21 new mixed-residual oracle assertions. |
| Additional CI equivalence tests | **30 files, 65 tests, 409 assertions pass; zero failures/errors/skips.** Four tutorial-model advisories: one small Gaussian class and three residual-score warnings of 1.6% of a standard error. These remain visible in the log. |
| Claims and retained references | **394 assertions pass**, no warnings or skips. Covers documentation registry, missing-pattern and retained Mplus inference/robust/weight/LMR checks. |
| Normal source build | **Pass**, with all four package vignettes built and their assertions active. The class-enumeration article also renders separately. |
| Final package check | `R CMD check --as-cran --no-manual --no-tests`: **zero errors, zero warnings, one note** because the offline environment cannot verify the current time. Remote incoming checks were disabled; tests were run separately as above. Examples, including `--run-donttest`, and vignette rebuilding pass. |
| Lint and whitespace | **Zero lints** and `git diff --check` clean for the changed package, harness and workflow files. |

The initial complete equivalence run correctly failed the new strict gate: it exposed **27 existing partial comparisons** and an empty-data error-class regression introduced by the Gaussian reference-scale guard. The guard now runs after the public sample-size validation; a new assertion verifies `latents_bad_data` for an empty regression input. The final installed-package regression run verifies this correction together with valid constant-response models.

The partial comparisons were repaired substantively. Mixed Gaussian/categorical residuals now use an independent posterior-weighted correlation reference, tested against ordinary Pearson correlations on explicitly replicated observations. Added zero-weight persons and groups now have their complete posterior/assignment outputs checked with a separate Gaussian Bayes calculation under the baseline parameters; exhaustive enumeration of all class paths agrees with that calculation to floating-point precision. No output was exempted to hide missing comparisons, and no numerical tolerance was loosened.

The [final combined harness report](validation/astra-fixes-2026-10-27/harness/final/REPORT.md) retains the original measured results for 11 unaffected suites and replaces the complete results of `depmixs4`, `hard-depmixs4`, `invariance` and `refusals` with their passing reruns. Every CSV row identifies its run of origin. This is explicitly a combined validation matrix, **not a claim that the entire 15-suite campaign was run a second time**. The strict gate passes on the combined result, with an additional assertion that every registered suite is present.

The original seven core regression blocks fail or error on the pre-fix installation and pass after repair. The permanent core file now contains 79 passing assertions, including the final empty-data guard; diagnostics has 80 and surface behavior has 22. The files are `tests/testthat/test-astra-core.R`, `test-astra-diagnostics.R` and `test-astra-surface.R`.

Evidence is retained under [validation/astra-fixes-2026-10-27](validation/astra-fixes-2026-10-27/): `unit.log` and result CSV/RDS; `final-installed-regression.log`; `harness/final/` and its original/rerun folders; `harness-self-tests-final.log`; `ci-equivalence.log`; `retained-references.log`; `build.log`; `00check.log`; `lint.log`; and the built source archive with its SHA-256 hash. `unit-source-manifest.json` identifies the complete unit-run snapshot; `source-manifest.json` identifies the final package source; `validation-source-manifest.json` records the final harness and workflow. The final package differs from the complete unit-run snapshot only in the ordering of the Gaussian guard and its added empty-data assertion.

### Remaining review boundaries

The fixes have not been committed, pushed or released. Version remains 0.97.4. The updated Linux CI workflow has been verified through its local test commands but has not been run remotely; the Mac fork skips therefore remain a platform-specific validation boundary. No new Monte Carlo calibration or fresh commercial-software campaign is claimed.

The synthetic HTML previews (sensitivity, diagnostics, README, class enumeration) were deleted at the user's request on 27 October 2026; the user does not want toy preview pages. These remaining boundaries do not reopen any of A01–A17, but they should remain visible when deciding to publish a release.

### Historical review

## 1. Release recommendation

**Hold the release.** The changes contain substantial, well-supported improvements, but the candidate is not ready to ship. A normal source build fails in the enumeration vignette, a newly introduced LTA covariate transformation can change estimates when only zero-weight people are added, and several public API and diagnostic regressions are reproducible. The focused equivalence rerun also fails its acceptance conditions.

The important distinction is between evidence that many established calculations remain correct and evidence that the complete release candidate is ready. The former is strong; the latter is presently negative. Passing unit tests or citing an older equivalence report does not resolve the concrete counterexamples below.

This review **did not modify production code or existing tests**, commit changes, publish a release, or change the project's statistical conventions. It records defects and proposed repairs for a subsequent implementation pass. Synthetic reproductions and validation logs are retained in [validation/astra-review-2026-10-27](validation/astra-review-2026-10-27/).

## 2. Scope, baseline and method

### 2.1 What “since release” means here

The checkout has no local release tags. The latest commit explicitly titled “Release” and the tip of local `main`/`origin/main` is 0.97.1 at `24773e6`, dated 25 October 2026. That is the review baseline. No remote fetch, GitHub release lookup or CRAN publication verification was performed; this is a repository-based release comparison.

Four commits follow that baseline:

| Commit | Date, local commit timestamp | Declared version/change | Review assessment |
|---|---|---|---|
| `8ac3576` | 26 Oct, 10:16 | Audit repairs; version advances to 0.97.2 | Broad numerical, inference, adapter, prediction and validation repairs; most added lines are retained evidence rather than runtime code. |
| `d5d4ae8` | 26 Oct, 23:36 | Full-output equivalence harness and two review rounds | Largest behavioral change: weighting, complete-data entropy, relative bounds, stopping, centering inference, growth coordinates, LTA identification and expanded validation. Several findings below originate here. |
| `285eca6` | 27 Oct, 09:27 | 0.97.3; interaction repairs and equivalence-kit conventions | Repairs interactions involving partial pairs, extreme positive weights, weighted summaries, Gaussian floors, growth identification and reference fixtures. Important remaining zero-weight consumers are still inconsistent. |
| `5596a15` | 27 Oct, 10:32 | 0.97.4; bounded condition calls | Addresses oversized error calls, but its new binding mechanism introduces name collisions and repeated expression evaluation. |

The committed range changes **399 files**, with **552,112 inserted and 3,281 deleted lines**. Within `R/`, **80 files** change, with **5,707 insertions and 960 deletions**. The broader inventory is:

| Area | Changed files | Nature of additions |
|---|---:|---|
| Runtime R source | 80 | Estimators, inference, diagnostics, adapters, search and public API |
| Package tests | 62 | New regression suites, changed fixtures/assertions and test bookkeeping |
| Rd documentation | 15 | Argument defaults, numerical/statistical conventions and output contracts |
| Equivalence infrastructure/evidence | 74 | Full-output harness, suites, fixtures, comparison reports and provenance |
| Validation archives | 155 | Previous audits, scripts, tables and retained execution evidence |
| Vignette article sources | 6 | Updates to workflows, model selection, regression and tree documentation |
| Deferred/future code | 4 | Supporting deferred-feature updates; not the released runtime surface |
| Other | 3 | DESCRIPTION, NEWS and a tutorial document |

The [complete committed file inventory](validation/astra-review-2026-10-27/committed-file-inventory.csv) lists every changed path and line count. Large generated evidence files were reviewed through their provenance, summaries, harness and representative checks, not by claiming a line-by-line statistical verification of every retained value.

At review start there were **42 tracked paths with uncommitted changes**, plus new runtime files `R/input-checks.R` and `R/diagnostics-checks.R`, three new review test files, and many pre-existing untracked evidence artifacts. Those package changes are included in the candidate assessment. Findings explicitly distinguish new committed regressions, current uncommitted additions and older gaps exposed by the newly adopted conventions.

### 2.2 How the review was conducted

- Read the release-to-HEAD diffs, current overlapping changes, package metadata, CI workflows, prior audit reports, handoff and learning log.
- Partitioned independent source review into model engines, inference/statistics and public API/reporting, as required by the applicable repository instructions; reconciled overlapping findings in this report.
- Copied the current package source to a separate scratch directory, recorded per-file SHA-256 hashes and installed it into an isolated library. Socket workers therefore see the same candidate package as the test process.
- Ran the package suite with `NOT_CRAN=true`; ran package checks, lint, harness self-tests and three selected equivalence suites; attempted the normal vignette-building source build.
- Used small synthetic examples to verify findings, including a real release-baseline comparison for the constant-response regression and a controlled in-memory chart correction to isolate the LTA defect.
- Preserved the user's working changes. Review output is this file and supporting review evidence/project notes.

The [source manifest](validation/astra-review-2026-10-27/source-manifest.json) identifies the snapshot used by the installed-package checks. A final hash comparison confirmed that none of its 410 source files changed during this review. Source locations below refer to the reviewed working tree; line numbers may move after repairs.

## 3. Findings at a glance

**P1:** release-blocking build failure or materially incorrect fitted result. **P2:** reproducible correctness/API defect or a validation gap that should be repaired before release. **P3:** contract/documentation or lower-impact consistency issue. Severity is separate from whether a defect was newly introduced.

| ID | Priority | Finding | Provenance |
|---|---|---|---|
| A01 | P1 | Normal package build fails on stale enumeration-vignette winners | Committed criterion changes; vignette not updated |
| A02 | P1 | Zero-weight people distort LTA covariate scaling and fitted estimates | Introduced in `d5d4ae8` |
| A03 | P2 | Constant-response guard rejects valid offset/no-intercept regressions | `d5d4ae8`, weighted counterpart `285eca6` |
| A04 | P2 | Condition-call binding overwrites arguments with colliding names | Introduced in `5596a15` |
| A05 | P2 | `multilca()` evaluates argument expressions twice | Committed forwarding; current early data check broadens it |
| A06 | P2 | Excluded people trigger LTA boundary flags and inference refusal | Older behavior retained under the new zero-weight convention |
| A07 | P2 | Categories seen only on excluded rows increase parameter count/BIC | Older category-discovery behavior retained |
| A08 | P2 | Unidentified bootstrap coefficients retain finite covariance rows | Incomplete repair in `d5d4ae8` |
| A09 | P2 | New diagnostics reintroduce excluded groups and empty LTA occasions | Uncommitted diagnostics helpers |
| A10 | P2 | Regression diagnostics omit or mislabel classification levels | Uncommitted family dispatch |
| A11 | P2 | Sensitivity plot calls unconverged starting states replicated maxima | Uncommitted plot addition |
| A12 | P2 | New replication output breaks equivalence coverage/invariance checks | Uncommitted output; harness not adapted |
| A13 | P2 | Harness acceptance gate permits partially compared columns | New full-output harness |
| A14 | P2 | CI does not run the new full-output harness or its self-tests | Integration missing for new infrastructure |
| A15 | P2 | Sensitivity documentation causes an Rd usage warning | Uncommitted method documentation |
| A16 | P3 | Newly supported diagnostics families violate advertised table/coercion contract | Uncommitted support expansion |
| A17 | P2 | Regression enumeration quality views discard sampling weights | Committed adapter inconsistent with round-three repairs |

Additional lower-priority observations are described below. The table does not count duplicate confirmations by multiple reviewers as separate defects.

## 4. Detailed findings and repair guidance

### A01 — The standard source build cannot complete

**Location:** [vignettes/enumeration.Rmd](vignettes/enumeration.Rmd), lines 145–163, especially 160–162. The same assumptions also occur in `vignettes/articles/class-enumeration.Rmd`, lines 169–171. Existing rendered README output still lists the old ICL winner.

A fresh `R CMD build --no-manual source`, with vignette building enabled, fails with:

```text
Quitting from enumeration.Rmd:145-163 [models-values]
icl_groups_pick == "profiles4_groups2" is not TRUE
Error: Vignette re-building failed.
```

The candidate's actual results using the vignette's unchanged data, grid, starts and seed are:

| Criterion, groups convention | Vignette asserts | Current result | Current value |
|---|---|---|---:|
| BIC | profiles4_groups3 | profiles4_groups3 | 16852.14 |
| ICL | profiles4_groups2 | profiles2_groups3 | 17114.68 |
| AWE | profiles3_groups2 | profiles2_groups2 | 17292.24 |

The changed complete-data entropy convention intentionally affects ICL/AWE rankings. The documentation was not fully carried through to the executable vignette. This is a hard source-build failure, not merely stale illustrative prose. The other three shipped vignettes (`lca`, `lpa`, `lta`) rebuilt successfully in the same attempt.

**Required repair:** verify and update the selected candidates, associated explanations, article counterpart and rendered README values. Preserve the assertions as executable documentation checks. Rebuild all shipped vignettes and review affected article/tutorial claims rather than merely deleting the failed assertion.

**Evidence:** [build-vignettes.log](validation/astra-review-2026-10-27/build-vignettes.log), [enumeration-values.R](validation/astra-review-2026-10-27/enumeration-values.R), [enumeration-values.log](validation/astra-review-2026-10-27/enumeration-values.log).

### A02 — Excluded people change LTA covariate coordinates and the solution

**Location:** [R/lta-engine.R](R/lta-engine.R), `.lta_design_chart()`, lines 125–139, and its fitting caller around 981–990; also inspect the corresponding inference caller.

The new rank checks use active initial/move masks, but the new coordinate chart does not. Means and standard deviations include all initial rows and all `layout$within` transitions, including zero-weight people and potentially moves beyond the last informative occasion.

A synthetic 35-person, three-occasion fit with two profiles and initial/transition covariates gives log likelihood **−226.9712124143479**, converged, without a boundary flag. Adding five people whose weights are zero and whose covariates are `1e8` gives approximately **−228.4333**, a loss of **1.4621**, while still reporting convergence. The largest changed transition coefficient moves **1.577333 robust standard errors** relative to the original fit.

Computing the chart from the original 35 contributing people, while retaining the complete data in the actual fitting process, restores the original likelihood exactly. That controlled diagnostic replacement isolates the chart as the cause. It was performed in memory only.

**Required repair:** estimate both fitting and inference charts on the same active initial/move support already used by identification checks. Continue retaining full designs for prediction of excluded people. Test likelihood, coefficients, convergence, boundary flags and inference after zero-weight padding and after trailing empty occasions.

**Evidence:** [core review](validation/astra-review-2026-10-27/probes/astra-core-review.md), `astra-core-more.R`, `astra-lta-chart-proof.R` and `astra-lta-magnitude.R` with matching logs in the [probes directory](validation/astra-review-2026-10-27/probes/).

### A03 — Constant raw responses do not imply a degenerate regression

**Location:** [R/mixture-regression.R](R/mixture-regression.R), lines 628–638 and weighted counterpart 507–515.

The new guard refuses every constant Gaussian response. A constant response can have positive residual variance under a no-intercept design or a varying fixed offset. Both are supported model specifications.

```r
set.seed(1027)
d <- data.frame(y = rep(2, 60), x = seq(-2, 3, length.out = 60), o = rnorm(60))
# Each with n_classes = 1, n_starts = 0, vcov_type = "none":
mixture_regression(y ~ 0 + x, data = d, n_classes = 1,
                   n_starts = 0, vcov_type = "none")
mixture_regression(y ~ 1 + offset(o), data = d, n_classes = 1,
                   n_starts = 0, vcov_type = "none")
```

The release at `24773e6` fits both and agrees with independent `lm()` maximum-likelihood calculations: residual variances **3.5840188** and **0.9624063**, log likelihoods **−123.4308542** and **−83.9867552**. The candidate instead raises `latents_bad_data` for both.

**Required repair:** assess degeneracy conditional on the actual design and offset. Cover weighted and unweighted cases and keep the residual floor meaningful. Retain a helpful refusal for genuinely degenerate intercept-only constant-response fits.

**Evidence:** `astra-core-probes.R/.log`, `astra-constant-release.R/.log` in the probes directory. This finding includes execution against the actual release source, not just a hypothetical reference calculation.

### A04 — Safe-call forwarding is not collision-safe

**Location:** [R/condition-calls.R](R/condition-calls.R), lines 30–35.

Caller symbols retain their own names, while non-symbol arguments are bound under generated names such as `.vars`. Those naming schemes overlap:

```r
.vars <- data.frame(x = rnorm(100), y = rnorm(100))
multilpa(.vars, vars = c("x", "y"), id = NULL,
         n_profiles = 1, n_starts = 1, seed = 1)
```

The generated binding for the `vars` expression overwrites the data-frame binding. The function then reports that `.vars` is a character vector. Direct `multilca()` has the same failure.

**Required repair:** allocate unique internal binding names independently of user symbols, keeping original labels separately for error messages. Test collisions with generated names and internal function names while retaining bounded condition calls.

### A05 — Forwarded expressions can execute twice

**Location:** [R/multilca.R](R/multilca.R), lines 47–65, and [R/condition-calls.R](R/condition-calls.R), line 30.

`multilca()` forces `vars` and `list(...)`, then re-evaluates their original syntax in the caller frame during binding. The uncommitted early data check additionally forces `data` before that re-evaluation. A counter-based data factory runs **twice** under direct `multilca()` and once under `lca()`. A stateful `seed = next_seed()` expression also executes twice.

Consequences include repeated expensive data construction and fitting with different random/stateful inputs from those initially inspected. This violates normal once-only R promise evaluation.

**Required repair:** forward already evaluated formal promises/dots; use the original syntax solely for call metadata. Test evaluation counts for data, variable selection, seed and other forwarded controls. A04 and A05 need related but distinct tests.

**Evidence for A04–A05:** [surface review](validation/astra-review-2026-10-27/probes/astra-surface-review.md), `astra-surface-probe2`, `probe3` and `probe5` scripts/logs.

### A06 — Zero-weight covariates can falsely disable LTA inference

**Location:** [R/lta-inference.R](R/lta-inference.R), lines 662–685; call in `R/lta-engine.R` around 1407–1408.

The boundary checker does not receive the active likelihood masks or sampling weights. In the A02 fixture, using excluded covariates of **10000** preserves the log likelihood to **1.99e−13**. The original fit has `boundary = FALSE` and a valid **14 × 14** covariance. The padded fit has `boundary = TRUE`, and `vcov()` raises `latents_boundary_fit` solely because probabilities approach a bound for excluded people.

The chart-only correction used for A02 does not remove this false flag, so this is an independent consumer of inactive rows. It predates the new chart but contradicts the newly documented absent-unit convention.

**Required repair:** apply boundary checks only where the observed weighted likelihood uses the probability. Test that adding excluded units preserves covariance availability, not just the likelihood.

### A07 — Excluded-only categories still change model dimension

**Location:** [R/categorical.R](R/categorical.R), lines 33–63.

The contributing-row mask checks variation but does not determine the ordinary fitted category set. In a one-class LCA with 60 active rows, a binary item and a ternary item, appending one zero-weight row containing a new category changes:

| Quantity | Original | Zero-weight row added |
|---|---:|---:|
| Parameters | 3 | 4 |
| Log likelihood | −107.505568153683 | −107.505568159683 |
| BIC | 227.294169994033 | 231.388514568255 |

The meaningful change is the added parameter and **4.09434457** BIC penalty; the tiny likelihood difference comes from the extra category's probability floor. This is retained older behavior, not a newly introduced regression, but it prevents the new “zero-weight units are absent everywhere” claim from being true.

**Required repair:** specify ordinary category discovery on active support and how excluded rows with unseen categories are represented. Treat deliberately declared joined-category grids separately. Test character, factor and numeric categories, not only numeric continuous indicators.

**Evidence for A06–A07:** `astra-core-final-checks.R/.log` and the core report.

### A08 — Bootstrap table masking does not reach the covariance attribute

**Location:** [R/lta-bootstrap.R](R/lta-bootstrap.R), lines 65–71.

The bootstrap marks unidentified trailing-occasion logits as `estimated = FALSE` and sets their reported standard errors/intervals to `NA`, then returns `stats::cov(kept)` without masking corresponding rows and columns.

An occasion-specific LTA with 80 people, four observed occasions and an all-missing fifth occasion for 20 people produced four usable bootstrap draws. Two unidentified logits have `NA` standard errors but covariance diagonal entries **0.04989773** and **0.06113789**, plus finite cross-covariances. Downstream contrasts can therefore consume apparent uncertainty estimates for parameters the table correctly declares unidentified.

**Required repair:** mask the covariance rows and columns consistently. Keep raw replicate values only as clearly identified raw evidence. Test covariance/table consistency, including interior and trailing empty occasions. Four draws are sufficient to expose this contract error; this experiment is not a bootstrap-calibration study.

**Evidence:** [inference review](validation/astra-review-2026-10-27/probes/astra-inference-review.md), `astra-lta-bootstrap-probe` scripts/logs and saved fixture.

### A09 — New posterior-based diagnostics use the wrong eligible units

**Location:** [R/diagnostics-checks.R](R/diagnostics-checks.R), lines 171–197 and 315–325; family dispatch in [R/report.R](R/report.R), lines 123–146.

There are two independently reproduced manifestations of the same missing eligibility policy:

1. **Zero-weight cross-level groups.** A 60-group restricted cross-level fit has ten zero-weight groups and `nobs(fit) = 50`. Existing fitted expected group counts are **30/20**. New diagnostics report **36/24** and print the smallest class as `40.0% (n = 24)`. The helper normalizes to all 60 posterior rows and also counts all modal assignments. A second independently generated fixture confirmed the same inflation.
2. **Empty general-LTA occasions.** A `multilpa_lta` fit with 340 rows, 320 informative, and `drop_empty = TRUE` falls through to the same posterior-only path. Empty occasions still have finite smoothed posteriors, so checking posterior completeness does not exclude them. Expected counts sum to **340 instead of 320**; entropy is **0.8543787 instead of 0.8888741** computed directly over informative rows.

**Required repair:** define one family-aware eligibility/weight adapter and use it for effective counts, modal counts, AvePP, OCC and entropy. Respect stored informative flags, `drop_empty`, positive weights and row-to-group mapping. Test both LTA engines and both empty-row policies. The current “every family” fixture exercises homogeneous LTA, which does not cover the general-LTA fallback.

### A10 — Regression diagnostics lose the model's nesting

**Location:** [R/diagnostics-checks.R](R/diagnostics-checks.R), `.latents_classifying_levels()`, lines 203–217.

Mixture regression stores class counts/nesting under `x$spec`; the new dispatcher reads top-level `x$n_classes` and `x$n_group_classes`. A fitted two-level, two-row-class/two-group-class model consequently receives only individual-level checks, omitting available group classification. A `class_level = "group"` fit labels its checks `individuals`. A one-class regression bypasses the intended “nothing to classify” branch.

**Required repair:** derive each family's levels from its actual metadata and correct posterior matrices. Cover ordinary, group-classified and two-level regression with one and several classes. The analogous omission of second-level growth checks was identified statically but was not independently fitted here; do not treat it as an equally verified counterexample.

### A11 — Sensitivity plots claim maxima without convergence

**Location:** [R/sensitivity.R](R/sensitivity.R), lines 197–206 and 220–221.

```r
s <- sensitivity(fit, seeds = 1:3, max_iter = 0, n_starts = 1)
p <- plot(s)
```

A real probe returns three rows with `converged = FALSE`, `iterations = 0`, and `best = TRUE`. The new plot titles them **“3 of 3 seeds reached the best maximum”**. It checks finite likelihood but ignores convergence. No EM step was performed.

**Required repair:** distinguish unfinished fits and count only converged fits as maxima, or describe them as observed likelihood values without a maximum claim. Keep failed and unfinished seeds visible. This is separate from the pre-existing unresolved design choice about deterministic Ward starts; changing the start strategy is not necessary to correct the false title.

**Evidence for A09–A11:** inference and surface reports and their retained synthetic probe scripts/logs. Plot labels/data were inspected programmatically; no visual implementation was changed.

### A12 — The replication table is not integrated into the equivalence harness

**Locations:** [R/get-data.R](R/get-data.R), new replication catalogue entry around 1212; [equivalence/harness/outputs.R](equivalence/harness/outputs.R), exemptions around 51–61; [equivalence/harness/suites/invariance.R](equivalence/harness/suites/invariance.R), expected-table transformations.

A fresh selected-suite run reports:

| Suite | Compared quantities | Agree | Values | Uncovered output columns |
|---|---:|---:|---:|---:|
| mclust | 465 | 465 | 22,402 | **45** |
| invariance | 898 | 894 | 38,912 | 0 |
| refusals | 50 | 50 | 50 | N/A |
| Total | **1,413** | **1,409** | **61,364** | **45** |

The 45 uncovered outputs are nine numeric replication fields across five fits. The four disagreements are the new replication table's `tolerance` and `log_likelihood` under outcome rescaling and data duplication. The likelihood itself transforms correctly; the harness leaves its copied expected replication fields untransformed.

These are validation integration defects, not evidence of four new likelihood-estimation errors. Nevertheless, the normal harness acceptance rules reject this candidate. The custom review runner saved results without overwriting the repository's existing generated reports; its successful process exit is **not** reported as a passing equivalence gate.

**Required repair:** explicitly classify each new output as independently checked or exempt for a justified search-bookkeeping reason, and apply the correct mathematical transformations where comparisons remain. Do not solve this by broad exemptions for actual estimates or by loosening tolerances.

**Evidence:** [harness-suites.csv](validation/astra-review-2026-10-27/harness-suites.csv), [harness-coverage.csv](validation/astra-review-2026-10-27/harness-coverage.csv), [harness-comparisons.csv](validation/astra-review-2026-10-27/harness-comparisons.csv).

### A13 — Partial output coverage is not a failing condition

**Location:** [equivalence/harness/outputs.R](equivalence/harness/outputs.R), lines 323–357; [equivalence/harness/run.R](equivalence/harness/run.R), uncovered-output filtering and final exit condition.

`output_ledger()` reports `partial` if some but not all values in a column were compared. The runner fails only on `status == "not covered"`. A synthetic inventory with three produced values and one comparison yields `partial`, yet the runner's uncovered gate detects zero rows. The remaining two values have no required exemption/reason.

**Required repair:** fail unaccounted partial coverage, or explicitly represent and justify every remaining value. Add a harness self-test that exercises the actual acceptance predicate, not merely ledger status creation. This does not prove the retained numerical comparisons are wrong; it disproves the stronger claim that the gate always enforces complete value coverage.

**Evidence:** [ledger-probe.R](validation/astra-review-2026-10-27/ledger-probe.R) and [ledger-probe.log](validation/astra-review-2026-10-27/ledger-probe.log).

### A14 — CI does not exercise the new harness

**Location:** [.github/workflows/equivalence.yaml](.github/workflows/equivalence.yaml), lines 41–51; other workflow files checked for harness invocations.

The workflow runs `tests/equivalence` plus a selected set of top-level equivalence files for claims and retained inference references. It does **not** invoke `equivalence/harness/test-harness.R` or `equivalence/harness/run.R`. No other checked workflow invokes them either.

Thus the newly added full-output coverage machinery, its stale API scan, invariance/refusal suites and its acceptance rules can break while the configured equivalence workflow remains green. The current A12 failure illustrates that gap. This is an integration omission for new functionality; the workflow itself was not changed in the reviewed commit range.

**Required repair:** add an explicit CI/release job for harness self-tests and the intended required suites; install their reference dependencies and define which skips are acceptable. Preserve source/version provenance and archive reports. A passing workflow must say which validation layer it actually ran.

### A15 — Sensitivity method documentation produces a check warning

**Location:** [man/sensitivity.Rd](man/sensitivity.Rd), usage lines 20–22 and argument list; roxygen source in `R/sensitivity.R`.

The new `as.data.frame.latents_sensitivity()` and `plot.latents_sensitivity()` signatures include `...`, but the shared Rd argument list does not document it. Fresh `R CMD check` reports:

```text
Undocumented arguments in Rd file 'sensitivity.Rd'
  ‘...’
```

**Required repair:** document the actual dots behavior in the roxygen source, regenerate Rd, and rerun usage/documentation checks. This is a genuine candidate warning, distinct from the temporary no-vignette archive warnings described in the validation section.

### A16 — New diagnostics support does not match the advertised table contract

**Location:** [R/report.R](R/report.R), lines 123–146 and return documentation; diagnostics catalogue in `R/get-data.R`.

For general LTA and newly accepted additive/cross-level/regression/growth fits, diagnostics may populate only `replication` and `checks`. On the executed general-LTA example, `get_results(dg, "entropy")` errors and `as.data.frame(dg)` returns the replication table, despite the documented promise that coercion returns entropy.

These families previously lacked diagnostics support, so this is a newly advertised contract mismatch rather than a removed formerly successful method. Older supported diagnostics objects were also emulated and checked: printing still works, and their coercion remains the entropy table.

**Required repair:** either supply consistent entropy tables with explicit unsupported rows, or clearly define and document a different primary table for these families. Do not let coercion meaning depend accidentally on the first table available.

### A17 — Regression enumeration plots mix weighted and unweighted quality measures

**Location:** [R/grid-views.R](R/grid-views.R), lines 264–267 and 285–288; public plotting consumer in [R/plot-gg.R](R/plot-gg.R), around line 1682.

The quality-view adapter explicitly returns no weights for regression grids. Round-three repairs made the candidate fits' classification tables weighted, but this consumer still reads unweighted posterior matrices. The primary individual entropy is copied from the grid and remains weighted, while other panels use a different convention.

Two installed-library reproductions show the disagreement:

| Quantity | Enumeration quality view | Candidate's weighted result |
|---|---:|---:|
| Minimum AvePP, weighted single-level regression | 0.6876206 | 0.6902239 |
| Group entropy, weighted two-level regression | 0.6751874 | 0.6669745 |

The second value was also independently calculated from the group posteriors and weights. This is a reporting defect, not a new estimate of the model's likelihood. The same adapter reads full posterior matrices, so eligibility must be handled together with weights.

**Required repair:** use the candidate's own classification summaries or a common correct row/group weight and eligibility adapter. Assert that every plotted quality quantity matches the corresponding candidate table under unequal weights and zero-weight padding.

**Evidence:** `astra-regression-grid-probe.R/.log` and `astra-regression-group-grid-probe.R/.log`; follow-up section in the inference review. The numeric data consumed by the plot were checked without changing plot implementation.

## 5. Review of the additions by functional area

### 5.1 Core Gaussian and categorical estimation

The shared EM updates now stop relative to observed-value counts rather than absolute log-likelihood magnitude. Default Gaussian floors use a relative data-scale bound. VEI/EVE/VVE shape updates no longer compare dimensionless shape quantities to a dimensional variance floor. SQUAREM's reported step count tracks steps retained rather than attempted acceleration steps. These changes address real unit-dependence and improve reproducibility under scale changes.

Categorical changes preserve category labels and encode pair provenance explicitly. Ordinary labels containing `::` are no longer mistaken for partially observed joined items. Under FIML, a half-observed declared pair contributes the marginal probability of its observed part. Negative-binomial shared-dispersion updates avoid a premature floor decision. These are statistical-result changes and need to remain prominent in release notes.

Coverage includes newly added kernel/categorical round-one and round-two regression files, ordinary-label interaction tests and reference-harness fixtures. Remaining issues are A07 and incomplete propagation of the active-support convention. No claim is made that every covariance structure, missingness pattern and optional feature interaction has been exhaustively refitted in this review.

### 5.2 Weights, multilevel models and information criteria

Weights now normalize to raw-positive support; extremely small positive ratios remain present rather than silently becoming absent. Zero-weight sample sizes, likelihood normalization and weighted classification summaries have been repaired across many families. Complete-data classification entropy includes group uncertainty plus conditional profile uncertainty for nontransition models, under both information-criterion sample-size conventions. The two conventions differ in sample-size penalties, not by discarding individual uncertainty from one of them.

These changes intentionally alter likelihood scale, ICL/AWE/CLC, model rankings and some reported counts relative to 0.97.1. Modal assignments remain actual counts where documented, while estimated shares and posterior summaries use sampling weights. The residual gaps in A02, A06, A07 and A09 demonstrate why the convention needs end-to-end tests, not just likelihood tests for ordinary Gaussian examples.

### 5.3 Inference and uncertainty

Grand-centering/scaling uncertainty is incorporated through stacked estimating equations; person transformations default to robust covariance and refuse unsupported observed/OPG paths. Covariate inference and imputation pooling carry that choice through. Boundary response probabilities now hold only constrained coordinates, and upper complements receive unavailable standard errors consistently. Numerical information calculations use better-conditioned coordinates for offsets and small-scale parameters. Label encoding protects lookup of terms containing punctuation.

The previous review records substantial calibration and independent numerical-derivative evidence. That evidence is useful but historical; this review did not rerun the entire Monte Carlo calibration campaign. The fresh LTA bootstrap counterexample A08 shows an output consistency gap even when the main standard-error table looks correct.

### 5.4 Mixture regression and growth mixtures

Regression changes include dispersion updates, conditioned inference, prediction alignment, weighted classification and stronger data/refusal checks. Growth fitting standardizes outcome and predictors, then transforms estimates, covariance and predictions back; Brownian bounds respect time units; singleton random-intercept identification is checked. These are meaningful improvements for time in years/days/calendar origins and outcomes in small units.

A03 is a concrete overbroad guard introduced while improving validation. Growth's absent-person starting path still fails at an extreme numeric scale (see Section 8), but moderate-scale controls and eight two-class synthetic fixtures held their deletion-reference likelihoods. Do not generalize that extreme failure into a claim that ordinary growth fitting is broken.

### 5.5 Transition models

Changes include active-design rank checks, conditioned covariate fitting/inference, explicit handling of unidentified logits/empty blocks, better bootstrap metadata, unit-aware bounds and improved sequence/data alignment. A02, A06 and A08 are specific remaining problems in those pathways. General LTA is a separate result class from homogeneous transitions; both must be included in public-API coverage, as A09 shows.

The earlier nonconverged AR multilevel-growth comparison remains a limitation, not a newly established defect. Agreement of two nonconverged transformed fits does not establish an independently verified maximum.

### 5.6 Three-step methods, resampling and selection

Empty rows are removed from classification-error matrices, BCH and step-three analyses using retained original positions for outcome alignment. ML distal search uses a finite deterministic start/optimization strategy and checks stationary, conditioned solutions. Group bootstrap draws are generated before refitting so refit RNG consumption cannot shift later resamples. Class matching incorporates dispersion information where equal means do not identify labels.

BLRT now uses valid replicates, reports invalid counts and returns NA below half valid. This is an explicit project decision, not an arithmetic accident. Its 50% rule is operational: it does not establish nominal calibration if refit failures depend on simulated outcomes. Historical calibration under other failure conditions should not be generalized automatically.

Comparison repairs preserve named fits whose labels overlap filtering arguments and improve data equivalence checks. Selection plots distinguish unavailable replication from actual replication and avoid claiming that repeated starts prove a global maximum.

### 5.7 Accessors, condition handling, diagnostics and presentation

Committed work improves objective-data identity checks, reconstruction after CSV round trips, punctuation-safe labels, prediction/refit alignment, TNA identifiers and sequence eligibility. Bounded condition calls address a real host failure in which data values could be deparsed into enormous messages. The security/resource intent is sound; A04/A05 identify defects in its implementation.

Uncommitted additions provide classed input errors, fitted-object type checks, default-method refusals, readable prints, compact `str()`, a replication table, sensitivity plot methods and structured diagnostic checks. The no-cutoff design follows the recorded user preference: report entropy/AvePP/OCC/class shares as values and flag factual convergence/boundary/empty-assignment facts. Repairs should preserve that choice rather than reintroduce rejected rule-of-thumb thresholds.

The main remaining work is correct family metadata, consistent contributing-unit eligibility and accurate plot/table contracts, not changing the chosen presentation philosophy.

### 5.8 Parallelism and resource behavior

Core requests above detected hardware count now produce a warning and are capped. Worker errors are bounded before storage/propagation. Parallel seed sensitivity fixes each refit's seed and uses one core per refit. The installed snapshot matters because PSOCK workers load installed namespaces rather than an arbitrary source-loaded parent namespace.

The new cap is not a memory-based admission controller, and this review does not claim memory use is bounded for every data size. Cross-platform/fork coverage is reported separately from current macOS tests. No new stress benchmark or live Windows/Linux execution was performed here.

### 5.9 Validation and documentation infrastructure

The new output inventory, per-value comparisons, explicit exemptions, provenance hashes, known-open register, hard reference suites and invariance/refusal checks are a major improvement over checking a few selected likelihoods. They also need their own acceptance tests and CI integration. A12–A14 are actionable gaps, not a reason to discard the harness.

The retained report says it was generated against **0.97.2**, with **8,916 quantities / 3,142,823 values** agreeing. Those are historical results for its recorded source, not fresh results for this dirty 0.97.4 candidate. The current three-suite rerun must not be replaced by that older headline. An empty known-open register only means no registered exception exists; it does not imply no defect remains outside those checks.

## 6. Validation performed during this review

Environment: **R 4.5.2**, macOS arm64, Apple Accelerate BLAS; `testthat` 3.3.2, `pkgload` 1.4.1. Heavy R commands used one BLAS thread. No new package dependencies were installed. All model reproductions used synthetic/public package data.

| Check | Fresh result | Interpretation |
|---|---|---|
| Install frozen source into isolated library | Passed | Package loads and workers can target the candidate install. |
| Full package unit/invariant suite, `NOT_CRAN=true` | **210 files / 1,827 tests; 12,488 passing expectations, 0 failures/errors, 2 warnings, 16 skips** | Combined full run and eight affected-file reruns, replacing those files without double counting; remaining skips are fork/Accelerate restrictions. |
| Harness self-tests | 28 tests, 103 passing expectations | Does not cover the partial-coverage acceptance hole in A13. |
| Selected harness: mclust, invariance, refusals | 1,413 quantities, 61,364 values; 4 disagreements and 45 uncovered outputs | Candidate fails normal harness acceptance conditions; details in A12. |
| Synthetic partial-ledger test | Reproduced acceptance gap | A13. |
| Lint with repository configuration | 0 lints | Source style checks passed; not evidence of statistical correctness. |
| Source build without vignettes/manual | Passed | Used to inspect remaining package checks independently of documentation rebuilding. |
| Check of that archive, examples enabled, tests/vignette rebuilding disabled | 0 errors, 3 warnings, 1 note | One real Rd warning; two missing-vignette warnings caused by this intentionally reduced build; NEWS parsing note. |
| Normal source build with vignette rebuilding | Failed in enumeration; other three vignettes passed | A01 is a real standard-build blocker. |
| Check examples | Passed | From the reduced-build check, not a full release check. |
| Focused code/probe reviews | Confirmed findings above | Includes fresh installed-library repeats for key statistical and reporting defects. |

The initial isolated full run had six local-socket errors and extra skips because external validation fixtures were absent from the staging copy. Eight affected test files were rerun against the same installed package with local socket access and the existing fixtures made available. All six errors cleared; the rerun produced 1,076 passing expectations and no failures/errors. The final totals above replace those eight files in the original results rather than adding duplicate tests. The remaining **16 skips** explicitly require a fork-safe BLAS; they were not forced to run under Accelerate. The **two warnings** are small effective membership in `test-review-1015-transitions.R` and collapsed group classes in `test-search.R`; they are visible fixture warnings, not hidden passes. See [final-unit-summary.txt](validation/astra-review-2026-10-27/final-unit-summary.txt), [final-unit-results.csv](validation/astra-review-2026-10-27/final-unit-results.csv) and [rerun.log](validation/astra-review-2026-10-27/rerun.log).

The reduced-build check also notes NEWS headings it cannot parse as version sections. This includes historical headings; it is not presented as a newly introduced estimator defect. Its network dependency-index lookups could not reach CRAN/Bioconductor in the restricted environment, while installed-dependency checks completed. The two missing-`inst/doc` warnings are consequences of building with `--no-build-vignettes`, not additional independently discovered source failures.

The principal commands were `R CMD INSTALL -l /private/tmp/astra-review/lib /private/tmp/astra-review/source`, the archived `run-suite.R` using `testthat::test_dir(..., package = "latents", load_package = "installed")`, `Rscript equivalence/harness/test-harness.R`, and the archived `run-harness.R` calling `run_equivalence(c("mclust", "invariance", "refusals"))`. The reduced archive was checked with `R CMD check --no-manual --no-tests --no-build-vignettes --no-vignettes`; the separate normal build was `R CMD build --no-manual source`. `NOT_CRAN=true` was set for the unit suite, and `R_LIBS` pointed at the isolated install. The scripts and logs preserve the exact invocation details and session information.

A complete clean `R CMD check --as-cran` on a successfully built final source archive has **not** passed. It cannot be claimed while A01 remains. No fresh all-15-suite external harness, full top-level external-kit rerun, live Mplus/Latent GOLD session, Windows/Linux check or complete article/tutorial rendering campaign was performed. Prior evidence for those activities remains historical.

## 7. Compatibility and release communication

The release notes should explicitly tell users which results require refitting or reinterpretation:

| Change | User-visible consequence | Release communication needed |
|---|---|---|
| Relative default variance floor | Fits/boundary status may differ from old absolute-floor results, especially under small units | Document `NULL` semantics and stored effective floor; an explicit numeric bound retains its dimensional meaning. |
| Observed-value EM stopping | Iteration counts and near-optimum estimates may differ | Do not promise bitwise compatibility with 0.97.1. |
| Complete-data entropy for ICL/AWE/CLC | Model rankings can change materially | Update executable examples, README and tutorials; A01 is an actual example. |
| Zero-weight normalization/support | Sample size, likelihood scale, criteria and classification summaries change | State the positive-support convention and repair remaining inconsistencies before making universal claims. |
| Grand/person transform inference | SEs change; some requested covariance types are now refused | Explain why robust defaults and refusals differ from saved analyses. |
| Half-observed declared categorical pairs | FIML uses information previously discarded | Existing estimates may change; this is a likelihood correction. |
| Bound response probabilities | Previously finite/zero SEs can become NA | Explain that NA reflects unavailable regular inference rather than a new fitting failure. |
| BLRT valid-replicate policy | Some formerly NA p-values become numeric | Report invalid counts and the conditional-valid-replicate policy prominently. |
| Saved multigroup fits | Missing stored conditional entropy can make corrected criteria unavailable | Tell users to refit; do not silently reconstruct an unsupported quantity. |
| Classed errors and compact object display | Condition classes, call payloads and printed output change | Preserve numerical/data access through explicit tables and raw `unclass()` inspection. |
| Sensitivity result class | Remains a data frame but gains a specialized plot method | Document plot semantics, especially failures and nonconvergence. |
| Expanded diagnostics | New checks and supported families | Resolve A09/A10/A16; avoid implying all requested tables exist for every family. |

`numDeriv` was added to Suggests. There are no new Imports in the committed DESCRIPTION diff. Exact dependency versions from the local run are recorded in the test/session logs; no claim about current online package versions is made.

## 8. Lower-priority observations and limits

- **Extreme growth starting values:** adding four zero-weight people with outcomes around ±1e100 to a one-class random-intercept fit raises `latents_em_decrease` although active data alone converge. At ±1e5, ±1e8, ±1e10 and ±1e20, the same padding preserves likelihood within 1.9e−13. Eight two-class moderate-scale fixtures also agree within 1.2e−13. This is a confirmed extreme-scale active-support limitation, not evidence of a common multi-class failure.
- **Deterministic first starts:** the prior handoff leaves a user decision open about random-only sensitivity refits. A deterministic Ward first start can produce identical fits across seeds. Do not interpret such replication as independent evidence of a global maximum, and do not silently change the fitting strategy during repairs.
- **Bivariate residual calibration:** the project already documents conservative categorical chi-square behavior; no parametric-bootstrap calibration was added here. Statistical test availability and unsupported pairs should remain explicit.
- **SQUAREM scale dependence:** the prior low-priority step-length observation remains; the review did not establish a new practical failure beyond the specific cases reported.
- **Cross-level mapping:** the new diagnostic fallback says weighted individual-level checks lack a row-to-group map, although a reviewed fit retains `sufficient_statistics$index`. A correct adapter may recover that capability. Returning n/a is preferable to inventing unweighted values, but the stated limitation needs checking.
- **Modal-versus-estimated shares:** under unequal weights, an unweighted modal share and weighted estimated share describe different populations. Their difference should not be labelled solely as classification uncertainty.
- **Whitespace:** committed generated session/evidence logs contain trailing whitespace under `git diff --check`; package-source/current-change checks examined here did not show the same issue. This is housekeeping, not a numerical release blocker.
- **New review file packaging:** `Astra.md` is a requested repository report and is not currently excluded by `.Rbuildignore`. Before building the eventual release from the repository, exclude the report if it should not ship, or otherwise handle that top-level artifact deliberately. The frozen package tested here predates creation of this report.

## 9. Required work before release

1. **Repair A01 and A02 first.** Restore a normal source build and active-support LTA chart invariance. Keep a reproduction proving the old code fails and the corrected code passes.
2. **Repair the public API regressions A03–A05.** Exercise supported offset/no-intercept regressions, collision-prone caller names and exactly-once expression evaluation.
3. **Complete the support/identification rules A06–A10 and A17.** Test likelihood, criteria, effective counts, covariance availability and all classification levels, including plotted enumeration summaries; do not stop after testing one scalar output.
4. **Correct sensitivity and diagnostics contracts A11/A15/A16.** Preserve the user's no-rule-of-thumb-cutoff design and the pending deterministic-start choice.
5. **Restore validation gates A12–A14.** Adapt replication coverage/transforms, fail unaccounted partial coverage and run the harness in CI with explicit dependency/skip policy.
6. **Update result-bearing documentation.** Rebuild shipped vignettes, review article counterparts and refresh rendered README results. Audit claims tied to old criterion/weight conventions.
7. **Freeze the exact final candidate.** Ensure new runtime files/tests are included, choose the release version, and record a source manifest. Do not rely on a version string shared by multiple dirty source states.
8. **Run final checks on that artifact.** Full unit suite, affected and full required equivalence layers, harness self-tests, normal build and clean `R CMD check --as-cran`, plus supported-platform CI. Record all skips and reasons. Re-run after any repair affecting the tested behavior.

The implementation effort should preserve the useful architecture already added: explicit conventions, source-provenance evidence, classified errors, full-output inventories and transparent diagnostics. The remaining failures are sufficiently concrete to support targeted repairs and regression tests.

## 10. Evidence index

- [All committed changed files and line counts](validation/astra-review-2026-10-27/committed-file-inventory.csv)
- [Frozen candidate source SHA-256 manifest](validation/astra-review-2026-10-27/source-manifest.json)
- [Core-model review and detailed coverage](validation/astra-review-2026-10-27/probes/astra-core-review.md)
- [Inference/statistical review and detailed coverage](validation/astra-review-2026-10-27/probes/astra-inference-review.md)
- [Public API/reporting review and detailed coverage](validation/astra-review-2026-10-27/probes/astra-surface-review.md)
- [Initial package unit-suite log](validation/astra-review-2026-10-27/unit.log), [affected-file rerun](validation/astra-review-2026-10-27/rerun.log), [final combined summary](validation/astra-review-2026-10-27/final-unit-summary.txt)
- [Package check log](validation/astra-review-2026-10-27/check.log)
- [Normal build/vignette failure log](validation/astra-review-2026-10-27/build-vignettes.log)
- [Harness self-tests](validation/astra-review-2026-10-27/harness-self-tests.log)
- [Selected-suite results](validation/astra-review-2026-10-27/harness-suites.csv), [coverage](validation/astra-review-2026-10-27/harness-coverage.csv), [comparisons](validation/astra-review-2026-10-27/harness-comparisons.csv)
- [Lint result](validation/astra-review-2026-10-27/lint.log)

Reproduction scripts are archived with their logs under `validation/astra-review-2026-10-27/probes/`. Some preserve their original `/tmp/astra-*` fixture/library paths; regenerate their fixtures or adjust those paths when replaying elsewhere. Controlled in-memory modifications used to isolate a cause are explicitly described and are not production fixes. Historical project reports are cited as historical evidence and have not been overwritten.
