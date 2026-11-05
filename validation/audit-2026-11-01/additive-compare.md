# Audit 2026-11-01 — slice `additive-compare`

Files: R/additive.R, R/additive-compare.R, R/additive-inference.R,
R/additive-methods.R, R/tree.R, R/compare-models.R, R/compare-lrt.R
(116 top-level objects: 114 functions and 2 constants).

Checking scripts (synthetic data; `OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
NOT_CRAN=true`, `devtools::load_all()`, one R process at a time, no workers)
are in the session scratchpad
`/private/tmp/claude-503/-Users-mohammedsaqr-Documents-Github-latents/f561d1fc-3262-4e98-94dd-fd3bab536edf/scratchpad/addcmp/`:

- `decomp.R`: a two-class additive truth (120 groups of 3-8 ratings, 2
  indicators) drawn with the test helper `additive_draw()`; the fitted log
  likelihood against a dense Cholesky evaluation of every group's
  `N(mu 1, W I + T 11')` mixture (difference 0); class means against the
  M-step fixed point (posterior average of the conditional intercepts, max
  difference 3e-12); a conditional intercept by the formula
  `mu + nT/(W + nT)(a - mu)` (equal); intercept means as the posterior mix of
  the conditional means (difference 0); parameter counts of the five
  family/between-variance specifications at H = 3, d = 2 by hand
  (16, 12, 12, 20, 16; all equal); the `.additive_nested()` 5 x 5 table.
- `inf.R`: analytic group scores against `numDeriv::grad`, observed
  covariance against `solve(-numDeriv::hessian)`, natural-scale SEs against a
  numDeriv delta method, robust covariance against `V S'S V G/(G-1)` by hand,
  for all five specifications (largest relative difference 1e-4, numDeriv's
  own accuracy on a log likelihood of -2400); BH adjustment over the four
  mean rows only (checked by hand).
- `tree.R`: mass, share_of_from and share_of_to sums and the flow margins on
  an LPA grid (individuals, groups, `from`/`to`, `candidates`, `"all"`), a
  weighted LPA grid, a family grid and a family grid with a zero-weight
  group: every sum within 7e-15 of 1, margins within 7e-16, no NA.
- `ic.R`: `compare_models()` against `get_results(fit, "information_criteria",
  format = "wide")` (16 columns, difference 0) and the grid's candidates
  table (29 columns, difference 0); delta BIC and weights by hand; an
  EEI/VVI/VVV LRT by hand (19.74732 on 2 df, p = 5.151e-05, equal); the
  family grid table against each candidate's `get_results(fit, "fit")`.
- `cmp2.R`, `meth.R`, `enum.R`: argument and data refusals of every verb in
  the slice, weighted and zero-weight additive fits (zero-weight groups give
  the same log likelihood and BIC as dropping them), the recovery table.
- `misc.R`, `sim2.R`: exact class matching against brute force over 200
  random overlap matrices; greedy matching; colour de-duplication;
  `.additive_simulate()` moments over 3000 draws (z = 1.43 and -1.18 for the
  two grand means; pooled within variance 2.058 against 2.064); a family
  `bootstrap_lrt()` with `iter = 2`; every additive and tree plot returns a
  ggplot or a `latents_plots` list (visual output unchanged, not rendered);
  `.latents_lrt_spans()` on an affine recoding (TRUE) and an unrelated
  predictor (FALSE).

`.additive_nested()` and the analytic LRT (`.latents_lrt_*`,
`.multilpa_lrt_nesting()`) have no pair in common: `compare_models()` refuses
group-class family fits (class `latents_bad_argument`, now with a message
that says where they are compared), so `.additive_nested()` serves
`bootstrap_lrt()` alone. It was checked on all 25 ordered pairs of the five
specifications against the block-rank rule (a block shared by every class is
a special case of a class-specific one); the analytic LRT was checked by hand
on the profile pairs above and by its own regression files.

## 1. Function table

| function | file:line | checked how | verdict |
|---|---|---|---|
| `.additive_identity_tolerance` (constant) | R/additive.R:14 | read; its roxygen block also carried a stray copy of `.additive_prepare()`'s description | fixed (doc) |
| `.additive_used_groups` | R/additive.R:25 | ran (`meth.R`, zero-weight groups) | ok |
| `.additive_prepare` | R/additive.R:62 | ran (`decomp.R`, `meth.R`); test-additive-engine.R | ok |
| `.additive_structure` | R/additive.R:158 | ran (`decomp.R`, all five specifications) | ok |
| `.additive_check_parameters` | R/additive.R:176 | read; test-additive-engine.R | ok |
| `.additive_row` | R/additive.R:209 | ran via `.additive_log_density()` | ok |
| `.additive_log_density` | R/additive.R:226 | ran (`decomp.R`: dense likelihood difference 0) | ok |
| `.additive_expectation` | R/additive.R:262 | ran (`decomp.R`: conditional intercept and intercept means by hand) | ok |
| `.additive_slice` | R/additive.R:312 | ran via the M-step | ok |
| `.additive_between_score` | R/additive.R:330 | read; test-additive-engine.R (KKT) | ok |
| `.additive_pool` | R/additive.R:350 | ran via the M-step of shared blocks (`decomp.R`) | ok |
| `.additive_initialize` | R/additive.R:369 | read; ran via every fit | ok |
| `.additive_maximization` | R/additive.R:474 | ran (`decomp.R`: fixed point 3e-12) | ok |
| `.additive_accelerator` | R/additive.R:549 | read; ran via every fit | ok |
| `.additive_observed` | R/additive.R:582 | read | ok |
| `.additive_reference_scales` | R/additive.R:595 | read | ok |
| `.additive_em` | R/additive.R:619 | ran via every fit | ok |
| `.additive_polish` | R/additive.R:661 | ran (family grid: boundary candidate `additive_groups3`); test-review-1026-additive-r2.R | ok |
| `.additive_collapsed_classes` | R/additive.R:754 | read; test-review-1026-additive-r2.R (ADD2-5) | ok |
| `.additive_arguments` | R/additive.R:792 | read | ok |
| `.additive_rows` | R/additive.R:800 | ran via the counts | ok |
| `.additive_count_parameters` | R/additive.R:810 | ran (`decomp.R`: five hand counts); new test | ok |
| `.additive_family_label` | R/additive.R:819 | ran via print | ok |
| `.additive_fit` | R/additive.R:837 | ran (`decomp.R`, `meth.R`, weighted and zero-weight fits) | ok |
| `print.multilpa_additive` | R/additive.R:1058 | ran (`decomp.R`); _snaps/additive-fit.md | ok |
| `.additive_model_code` | R/additive-compare.R:9 | ran via the family grid | ok |
| `.additive_simulate` | R/additive-compare.R:22 | ran (`sim2.R`: moments over 3000 draws; ids and sizes kept) | ok |
| `.additive_nested` | R/additive-compare.R:49 | ran (`decomp.R`: 25 ordered pairs against the block ranks); new test | ok |
| `.additive_bootstrap_lrt` | R/additive-compare.R:63 | ran (`misc.R`, iter = 2); test-additive-compare.R | ok |
| `.additive_enumerate` | R/additive-compare.R:189 | ran (`enum.R`: 0, 2.5, NA, "2", integer() made grids of failed candidates or an unclassed error) | fixed |
| `get_results.latents_family_enumeration` | R/additive-compare.R:288 | ran (`tree.R`: unknown `what`/`criterion` gave simpleError) | fixed |
| `as.data.frame.latents_family_enumeration` | R/additive-compare.R:307 | ran (`ic.R`) | ok |
| `print.latents_family_enumeration` | R/additive-compare.R:314 | ran (`tree.R`) | ok |
| `summary.latents_family_enumeration` | R/additive-compare.R:328 | read; test-additive-families.R | ok |
| `plot.latents_family_enumeration` | R/additive-compare.R:336 | ran (`misc.R`, returns ggplot) | ok |
| `.additive_block_names` | R/additive-inference.R:13 | ran via `coef()` names | ok |
| `.additive_block_pack` | R/additive-inference.R:20 | ran (`inf.R`) | ok |
| `.additive_block_unpack` | R/additive-inference.R:26 | ran (`inf.R`) | ok |
| `.additive_block_combine` | R/additive-inference.R:34 | ran via the scores of shared blocks (`inf.R`) | ok |
| `.additive_names` | R/additive-inference.R:48 | ran (`inf.R`, confint row names) | ok |
| `.additive_pack` | R/additive-inference.R:66 | ran (`inf.R`) | ok |
| `.additive_unpack` | R/additive-inference.R:79 | ran (`inf.R`: pack/unpack round trip in the numDeriv likelihood) | ok |
| `.additive_group_scores` | R/additive-inference.R:108 | ran (`inf.R`: against numDeriv gradient, five specifications) | ok |
| `.additive_information` | R/additive-inference.R:149 | ran (`inf.R`: against numDeriv Hessian) | ok |
| `.additive_inference` | R/additive-inference.R:186 | ran (`inf.R`: observed and robust by hand; `step.R`: six bad `step` values gave simpleError) | fixed |
| `.additive_sample_groups` | R/additive-inference.R:269 | read | ok |
| `.additive_weak_class_threshold` (constant) | R/additive-inference.R:280 | read | ok |
| `.additive_warn_weak` | R/additive-inference.R:285 | ran (warning seen in `inf.R`) | ok |
| `.additive_natural` | R/additive-inference.R:311 | ran (`inf.R`: Jacobian against numDeriv) | ok |
| `.additive_block_labels` | R/additive-inference.R:335 | ran via the parameter table | ok |
| `.additive_parameter_table` | R/additive-inference.R:353 | ran (`inf.R`: SEs, BH over the four mean rows by hand) | ok |
| `.additive_newton` | R/additive-inference.R:422 | ran via every converged interior fit; test-additive-fit.R | ok |
| `.additive_effective_groups` | R/additive-inference.R:493 | ran (`inf.R`, group_classes table) | ok |
| `.additive_tables` | R/additive-methods.R:7 | read | ok |
| `.additive_parameters_or_notice` | R/additive-methods.R:15 | read; test-additive-methods.R | ok |
| `.additive_default_vcov` | R/additive-methods.R:32 | read | ok |
| `.additive_effective_or_na` | R/additive-methods.R:43 | ran (`ic.R`, family grid NA for boundary fits) | ok |
| `.additive_entropy` | R/additive-methods.R:62 | ran via the fit table | ok |
| `.additive_recovery` | R/additive-methods.R:68 | ran (`meth.R`: table, missing args, non-constant truth) | ok |
| `.additive_table` | R/additive-methods.R:107 | ran (`meth.R`: every table, weighted and zero-weight) | ok |
| `get_results.multilpa_additive` | R/additive-methods.R:214 | ran (`meth.R`: unknown what/adjust/level classed) | ok |
| `as.data.frame.multilpa_additive` | R/additive-methods.R:237 | ran (`meth.R`) | ok |
| `coef.multilpa_additive` | R/additive-methods.R:248 | ran (`inf.R`, `meth.R`) | ok |
| `vcov.multilpa_additive` | R/additive-methods.R:265 | ran (`inf.R`: robust by hand) | ok |
| `confint.multilpa_additive` | R/additive-methods.R:286 | ran (`inf.R`, `meth.R`) | ok |
| `logLik.multilpa_additive` | R/additive-methods.R:310 | ran (`meth.R`) | ok |
| `nobs.multilpa_additive` | R/additive-methods.R:319 | ran (`meth.R`: 118 with two zero-weight groups of 120) | ok |
| `parameter_inference.multilpa_additive` | R/additive-methods.R:323 | ran (`meth.R`: within-group reordering accepted, changed data refused, bootstrap refused) | ok |
| `summary.multilpa_additive` | R/additive-methods.R:375 | ran (`meth.R`) | ok |
| `print.summary_multilpa_additive` | R/additive-methods.R:391 | read; test-additive-methods.R | ok |
| `as.data.frame.summary_multilpa_additive` | R/additive-methods.R:405 | ran (`meth.R`) | ok |
| `get_results.summary_multilpa_additive` | R/additive-methods.R:412 | ran (`meth.R`) | ok |
| `.additive_key` | R/additive-methods.R:418 | ran via plots | ok |
| `plot.multilpa_additive` | R/additive-methods.R:436 | ran (`misc.R`: three views return ggplot; unknown `what` classed) | ok |
| `.latents_tree_specs` | R/tree.R:28 | ran (`tree.R`, all paths) ; test-tree.R | ok |
| `.latents_tree_candidates` | R/tree.R:199 | ran (`from = c(1, 3)` on a family grid named the one-group-class candidate) | fixed |
| `.latents_check_tree_arguments` | R/tree.R:255 | read; test-tree.R | ok |
| `.latents_tree_word` | R/tree.R:282 | read | ok |
| `.latents_tree_table_of` | R/tree.R:296 | ran (`tree.R`) | ok |
| `.latents_tree_fit` | R/tree.R:320 | read; test-tree.R | ok |
| `.latents_singular` | R/tree.R:335 | read | ok |
| `.latents_check_tree_filter` | R/tree.R:345 | read; test-tree.R | ok |
| `.latents_tree_posteriors` | R/tree.R:364 | ran (`tree.R`: weighted, zero-weight group) | ok |
| `.latents_permutations` | R/tree.R:411 | ran (`misc.R`: 120 unique rows for n = 5) | ok |
| `.latents_match_classes` | R/tree.R:428 | ran (`misc.R`: optimal against brute force, 200 matrices) | ok |
| `.latents_greedy_match` | R/tree.R:446 | ran (`misc.R`, 3 x 2) | ok |
| `.latents_distinct_colours` | R/tree.R:465 | ran (`misc.R`) | ok |
| `.latents_tree_layout` | R/tree.R:487 | ran (`tree.R`: shares, flows) | ok |
| `.latents_tree_anchor` | R/tree.R:555 | ran via `plot(grid, what = "tree")`; test-tree.R | ok |
| `.latents_tree_results` | R/tree.R:578 | ran (`tree.R`) | ok |
| `.latents_tree_table` | R/tree.R:599 | ran (`tree.R`: invariants); new test | ok |
| `.gg_tree_frames` | R/tree.R:635 | ran via tree plots (`misc.R`) | ok |
| `.gg_tree_plot` | R/tree.R:690 | ran via tree plots | ok |
| `.gg_view_tree` | R/tree.R:763 | ran via tree plots | ok |
| `compare_models` | R/compare-models.R:185 | ran (`cmp2.R`, `ic.R`: `cutoff = 2` gave simpleError; additive fits got a misleading refusal; a regression grid was told "takes an enumeration grid") | fixed |
| `plot.latents_comparison` | R/compare-models.R:323 | ran (`cmp2.R`, ggplot) | ok |
| `.multilpa_compare_candidates` | R/compare-models.R:365 | ran (`cmp2.R`: names, counts, selection, failed candidates) | ok |
| `.multilpa_compare_fits` | R/compare-models.R:428 | ran (`ic.R`: every criterion against the fit's own table); `cutoff` check moved to the entry | fixed |
| `.latents_compared_model` | R/compare-models.R:526 | ran (`cmp2.R`) | ok |
| `.latents_comparison_classification` | R/compare-models.R:558 | ran (`cmp2.R`) | ok |
| `.latents_comparison_measurement_fit` | R/compare-models.R:595 | read | ok |
| `.latents_comparison_diagnostic_view` | R/compare-models.R:605 | read; test-lta-compare.R | ok |
| `.latents_comparison_structure` | R/compare-models.R:616 | read | ok |
| `.latents_same_comparison_data` | R/compare-models.R:626 | read; test-review-1026-codex.R | ok |
| `.latents_comparison_signature` | R/compare-models.R:652 | read; test-features-joint-items.R | ok |
| `.latents_comparison_responses` | R/compare-models.R:694 | read; test-review-1101-compare-lrt.R | ok |
| `.latents_comparison_lrt` | R/compare-lrt.R:17 | ran (`ic.R`: three-model LRT by hand) | ok |
| `.latents_lrt_pair` | R/compare-lrt.R:76 | ran (`cmp2.R`: class counts, df <= 0 refused) | ok |
| `.multilpa_lrt_nesting` | R/compare-lrt.R:161 | read; test-compare-lrt.R, test-review-1101-compare-lrt.R | ok |
| `.latents_lrt_option_rank` | R/compare-lrt.R:253 | read; test-review-1101-compare-lrt.R (F6) | ok |
| `.latents_lrt_spans` | R/compare-lrt.R:270 | ran (`misc.R`: affine TRUE, unrelated FALSE) | ok |
| `.multilpa_lrt_predictors_nested` | R/compare-lrt.R:297 | read; test-review-1101-compare-lrt.R (F1) | ok |
| `.lta_lrt_designs_nested` | R/compare-lrt.R:344 | read; test-review-1101-compare-lrt.R | ok |
| `.latents_lrt_structure` | R/compare-lrt.R:364 | read | ok |
| `.latents_lrt_singular_covariance` | R/compare-lrt.R:395 | read (EII/VII within EVE, VVE, EEV, VEV with >= 2 profiles and >= 2 continuous); test-review-1101-compare-lrt.R (F4) | ok |
| `.multilpa_lrt_covariance_nested` | R/compare-lrt.R:421 | read; test-compare-lrt.R | ok |

## 2. Fixed

**A1. A family grid's tables refused an unknown name with an unclassed
error.** `get_results(grid, "bogus")` and `get_results(grid, "candidates",
criterion = "xx")` on a `latents_family_enumeration` raised `simpleError`
(`'arg' should be one of ...`), where the transition grid and every fit raise
`latents_bad_argument`. Cause: bare `match.arg()`. Change: R/additive-compare.R:295
and :297 use `.latents_match_formal()` (partial matching kept:
`get_results(grid, "cand")` is unchanged). Test: "a family grid refuses an
unknown table or criterion with a class".

**A2. Invalid group-class counts became a grid of failed candidates.**
`enumerate_classes(..., family = "additive", n_group_classes = c(0, 2))`
returned a grid with rows `additive_groups0` (NA everywhere, converged FALSE);
`2.5`, `NA` and `"2"` gave grids where every candidate failed;
`integer()` raised `simpleError` ("'names' attribute [1] must be the same
length as the vector [0]"). The profile grid refuses these with
`latents_bad_argument`, but enumerate_classes() branches to the family path
before its check. Change: R/additive-compare.R:207-213 applies the profile
grid's check, with its message. Test: "a family grid refuses invalid
group-class counts before fitting".

**A3. `compare_models(cutoff = )` refused out-of-range values with an
unclassed error.** `compare_models(a = f, b = g, cutoff = 2)` raised
`simpleError` from a `stopifnot()` in `.multilpa_compare_fits()`. Change: the
check moved to the entry, R/compare-models.R:217-219, as `.latents_require()`
(`latents_bad_argument`); the helper's `stopifnot()` was removed. The
`@return` text now names the class. Test: "compare_models() refuses a bad
`cutoff` with a class" (with the LRT by hand).

**A4. Group-class family fits and regression grids got a misleading
refusal.** `compare_models(a = additive_equal, b = additive)` (two
`multilpa(family = ...)` fits) raised `latents_bad_argument` with "takes ...
two or more fits of one kind: from multilpa(), lpa(), lca() ...", which these
are; a regression grid was told it takes "an enumeration grid". Change:
R/compare-models.R:223-237 gives regression grids the family grid's "compared
in their own table" refusal and family fits a refusal that names
`enumerate_classes(family = )`, `as.data.frame(grid)` and `bootstrap_lrt()`.
Class unchanged (`latents_bad_argument`); the test asserts the class and the
word `bootstrap_lrt`, since only the message changed. Test: "group-class
family fits get their own refusal from compare_models()".

**A5. A tree's `from`/`to` dropped a second count on a family grid.**
`get_results(family_grid, "tree", from = c(1, 3), to = c(2, 3), model =
"additive")` returned the flows from `additive_groups1` to `additive_groups2`:
`resolve()` reads only `value[1L]` on a family grid (its candidates have no
profile count). Change: R/tree.R:218-221 accepts one count on a family grid,
as the refusal message already said. Test: "a tree's `from` and `to` take one
count on a family grid".

**A7. A bad `step` for additive inference raised an unclassed error.**
Reported by the coordinator from another agent's reading; reproduced
(`step.R`): `parameter_inference(additive_fit, step = )` with `-1`, `0`,
`NA`, `c(1e-4, 1e-3)`, `"a"` or `Inf` each raised `simpleError` ("`step` must
be finite and positive") from a bare `stopifnot()` in `.additive_inference()`.
Change: R/additive-inference.R:189-192 uses `.latents_require()`
(`latents_bad_argument`). Test: "additive inference refuses a bad `step`
with a class".

**A6 (documentation).** The roxygen block of `.additive_identity_tolerance`
(R/additive.R:7) began with a copy of `.additive_prepare()`'s description.
Removed; `@noRd`, no Rd change.

Each fix was broken once on purpose (A1 back to `match.arg()`, A2 check
short-circuited, A3 check disabled, A4 branch disabled, A5 back to `1:2`, A7 back to `stopifnot()`);
the matching test failed every time and the files were restored.

New invariant and calibration tests (no fix behind them, requested by the
slice brief): tree flows are a joint distribution (mass sums to 1 per pair,
share_of_from per from-class, share_of_to per to-class, margins equal the
shares) on a family grid and on an LPA grid at both levels and in
`get_results(grid, "all")`; the additive log likelihood against a dense
evaluation, the M-step fixed point of the class means, parameter counts by
hand; the `.additive_nested()` table; `compare_models()` against
`get_results(fit, "information_criteria")` and the grid's candidates table.

## 3. Logged, not fixed

- `compare_models(test = "wald")` raises `match.arg()`'s `simpleError`;
  tests/testthat/test-compare-lrt.R:148-149 asserts `class = "simpleError"`,
  so the classed fix (`.latents_match_arg(test, c("none", "lrt"), "test")`)
  waits for that test (see Needs another slice). Tried and reverted.
- `parameter_inference(additive_fit, data = )` and the family
  `bootstrap_lrt(data = )` compare group ids in first-occurrence order: the
  same data with rows reordered across groups is refused
  (`latents_bad_inference_data`); reordering within groups is accepted.
- `.additive_recovery()` sorts by class label as text: with 10 or more group
  classes `group_class_10` sorts before `group_class_2`.
- Tree table: a class with exactly zero posterior mass in a row gives NaN
  `share_of_from`/`share_of_to` and a NaN centre in the layout ordering; not
  reachable with converged fits on realistic data.
- `subset()` of a `latents_comparison` drops its `models` attribute, so
  `get_results(<subset>, "classification")` refuses with
  `latents_incomplete_fit` (classed, message says rows were taken out); the
  comment in `.latents_comparison_classification()` saying subsetting keeps
  attributes holds for `[` with some arguments only.
- `.additive_enumerate()`: a count above the number of groups still becomes a
  failed candidate row (`latents_unidentified`), as a too-large count does in
  the profile grid.

## 4. Needs another slice

- tests/testthat/test-compare-lrt.R:148-149 encodes `compare_models(test =
  "wald")` as `simpleError`; change it to `latents_bad_argument` and switch
  R/compare-models.R:216 to `.latents_match_arg(test, c("none", "lrt"),
  "test")` (owner of that test file).

## 5. Test runs

Run with `OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NOT_CRAN=true`,
`testthat::test_file()` per file after `devtools::load_all()` (expectation
counts summed per file):

```
test-audit-1101-additive-compare.R: PASS 55 FAIL 0 SKIP 0 ERROR 0
test-additive-compare.R: PASS 29 FAIL 0 SKIP 0 ERROR 0
test-additive-engine.R: PASS 29 FAIL 0 SKIP 0 ERROR 0
test-additive-families.R: PASS 24 FAIL 0 SKIP 0 ERROR 0
test-additive-fit.R: PASS 43 FAIL 0 SKIP 0 ERROR 0
test-additive-methods.R: PASS 82 FAIL 0 SKIP 0 ERROR 0
test-review-1026-additive-r2.R: PASS 32 FAIL 0 SKIP 0 ERROR 0
test-review-1026-additive.R: PASS 9 FAIL 0 SKIP 0 ERROR 0
test-tree.R: PASS 131 FAIL 0 SKIP 0 ERROR 0
test-review-pairs-trees.R: PASS 58 FAIL 0 SKIP 0 ERROR 0
test-compare-lrt.R: PASS 34 FAIL 0 SKIP 0 ERROR 0
test-review-1101-compare-lrt.R: PASS 29 FAIL 0 SKIP 0 ERROR 0
test-bootstrap-lrt-extended.R: PASS 24 FAIL 0 SKIP 0 ERROR 0
test-lta-compare.R: PASS 31 FAIL 0 SKIP 0 ERROR 0
test-api-compat-0921.R: PASS 44 FAIL 0 SKIP 0 ERROR 0
test-compat-0921.R: PASS 10 FAIL 0 SKIP 0 ERROR 0
test-release-audit-comparison.R: PASS 13 FAIL 0 SKIP 0 ERROR 0
test-audit-recent-api.R: PASS 44 FAIL 0 SKIP 0 ERROR 0
test-post-submission-review.R: PASS 60 FAIL 0 SKIP 0 ERROR 0
test-review-1026-codex.R: PASS 56 FAIL 0 SKIP 0 ERROR 0
test-grid-contract.R: PASS 237 FAIL 0 SKIP 0 ERROR 0
test-review-1026-selection.R: PASS 40 FAIL 0 SKIP 0 ERROR 0
test-review-1026-selection-r2.R: PASS 41 FAIL 0 SKIP 0 ERROR 0
test-informative-rows.R: PASS 57 FAIL 0 SKIP 0 ERROR 0
test-review-1015-lead.R: PASS 19 FAIL 0 SKIP 0 ERROR 0
test-audit-inference-repairs.R: PASS 86 FAIL 0 SKIP 0 ERROR 0
test-features-pair-covariance.R: PASS 69 FAIL 0 SKIP 0 ERROR 0
test-features-joint-items.R: PASS 61 FAIL 0 SKIP 0 ERROR 0
test-local-dependence.R: PASS 50 FAIL 0 SKIP 0 ERROR 0
test-enumeration.R: PASS 315 FAIL 0 SKIP 0 ERROR 0
```

(test-compare-lrt.R had one ERROR while `test` was classed: its line 148
asserts `simpleError`; that change was reverted, then 34/0/0/0.) The audit
file's 55 includes the `step` test added after the coordinator's request;
the three additive files were rerun after it (82, 43, 32 passes).

LRT equivalence kit, `run_package_kit_tests("lrt", output = <scratchpad>)`:
`lrt kit: PASS 2 FAIL 0 SKIP 0 ERROR 0` (27 comparisons with mclust, glca
and Mplus, 0 disagreeing).

`Rscript equivalence/engine-golden/compare.R` (full catalogue):
`identical: 242; within_tolerance: 0; DIFFERENT: 18; ERROR: 0`. Thirteen are
the bootstrap attribute rename (`replicate_rows` -> `replicate_identity`).
The other five are not from this slice: `cov_weights`,
`family_additive_weights` (probe `inference_observed` message) and
`reg_gaussian_weights_observed_refused` (refusal message) follow the
R/weights.R change that names `method = "bootstrap"` only where a family has
one (the additive family has none, so its message rightly lost it);
`lta_refusal_no_time` and `svc_fit_staged_one_class` changed `error_class`
in other slices' files. No difference in a number, table or print.

`devtools::document()`: only man/compare_models.Rd rewritten (the
`@return` sentence on `cutoff` and family fits; the file already carried
today's uncommitted LRT text); NAMESPACE unchanged.
