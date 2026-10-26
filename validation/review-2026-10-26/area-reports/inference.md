# Inference partition review (slug: inference), 2026-10-26

## Scope

Files: `R/inference.R`, `R/inference-service.R`, `R/bootstrap-inference.R`,
`R/pool-imputations.R`, `R/sensitivity.R`, `R/alignment-service.R`,
`R/simulation-service.R`, `R/accessors.R`. Codex commit 8ac3576 changed four of
them: accessors.R, alignment-service.R, inference.R, simulation-service.R.
The other four are unchanged since HEAD~1; I reviewed them as written.

Method: small synthetic fits under `devtools::load_all()`, with
OMP/VECLIB threads set to 1. Oracles were independent hand-written
likelihoods differentiated by `numDeriv`, `sandwich::vcovCL`,
`mice::pool.scalar`, analytic Hessians, and the pre-Codex tree (`git archive
HEAD~1`, read-only) to tell regressions from older behaviour. Scratch scripts
are in `.../scratchpad/inference/` (c1.R to c22.R).

**No production file was edited.** Both confirmed regressions either need a
fix outside this partition or need a tolerance decision. Details below.

## Verdict on each Codex change

| Change | Verdict | Evidence |
|---|---|---|
| `.latents_parameter_component()` escaping of `% : , [ ] \r` in categorical response labels (inference.R:5, 146, 1087). The accessor key (accessors.R:672) and plot-gg.R:750 changed to match. | **Correct for `multilpa` fits, but regression for covariate fits (INF-1), and incomplete (INF-4).** | `multilpa` category `"x, y"`: before the change, `get_results(fit,"responses")` SE was NA; now it is 0.0513. `multilpa_covariates` categories `"p:q"` and `"x, y"`: before, SEs were 0.0500 and 0.0775; now both are NA. Under `boundary="fix"` in a covariate fit, 0 of 2 bound logits are matched. |
| `.multilpa_same_values()` now compares numbers exactly (`tolerance = 0`), and refuses mixed numeric/text IDs whose text forms collide. | Correct for its purpose. | Origin 1e6 with one value shifted by 1e-6: `.multilpa_check_alignment` now refuses it for all 5 centring/scaling combinations; the pre-Codex tree accepted all 5. |
| `.multilpa_alignment_values()` recognises a write.csv round trip; the continuous check now reapplies the forward transform. | **Regression for person-centred fits (INF-2).** | CSV round trip of the original data. Before: accepted for all 5 transform combinations. Now: refused for `person/none` and `person/person`, and accepted in 0 of 10 random person-centred datasets. `sensitivity(fit, data = read.csv(...))` now errors (`latents_bad_inference_data`, "`a`, `b` disagree"); pre-Codex it ran. |
| `.multilpa_inference_matrix()`: reconstruction branch, plus exact comparison (`tolerance = 0` instead of `1e-12`). | Correct, with a behaviour change (INF-3). | `data = NULL` and `data = original` give "ok" for all 5 transform combinations. A frame changed by 1e-6 at origin 1e6 is refused. A CSV round trip is now refused for every fit (before, untransformed fits accepted it). |
| alignment-service.R: Gaussian log-variance-ratio and correlation features; NB `log1p` spread features; the dispersion replication fix `rep(1, nrow(count_means))`. | Correct. | Variances are bounded below by `min_variance > 0` (enforced), so the log is finite. The NB variance is >= 0, so `log1p` is finite. When means are equal, the profiles are now separated by their spread: recovered order `2 1`, signature finite. Known permutation `c(3,1,2)` is recovered as `order(perm)` for diagonal and full covariance; coefficients identical. `block$variances` is present in every caller: profile, covariate, transition (`blocks_of`) and lta measurement blocks (checked on a diagonal and a full transitions fit). |
| `.multilpa_start_variances()` builds the matrix with `byrow` instead of `t(vapply())`. | Correct (the old code was wrong for 1 indicator). | Test: a 1-indicator, 2-profile start gives a 2x1 matrix `c(2,5)`; a 2x2 start gives `matrix(c(1,4,3,6),2)`. |
| simulation-service.R: documentation of the Box-Muller limitation. | Correct (documentation only). | The caller's stream is restored after a seeded bootstrap and after an error. A fresh session stays fresh. The Box-Muller cached variate is not preserved (FALSE), exactly as documented. |

## Findings

| id | severity | file:line | problem | reproduction (actual vs expected) | status |
|---|---|---|---|---|---|
| INF-1 | medium | Producer: inference.R:144-147 and 1077-1095 (escaped). Consumer: accessors.R:672 (escaped). Out-of-partition producer: covariate-inference.R:958, 967 (raw `sprintf("%s:%s")`). Consumers: covariate-inference.R:395, 575. | Codex escaped the `multilpa` labels but not the covariate-family labels. The accessor key and `.multilpa_response_bound()` are now escaped, while covariate inference terms are raw. Result: covariate fits lose response standard errors, and covariate `boundary = "fix"` holds nothing whenever a category or indicator name contains `: , [ ] %`. | c4.R / c5.R: covariate fit with categories `"p:q"`, `"x, y"`, `"z"`. `get_results(cf,"responses")$probability_standard_error` is `NA NA NA NA 0.0617 0.0566`; expected (pre-Codex) `0.0500 0.0524 0.0775 0.0665 0.0617 0.0566`. With `min_probability` set to 0.2, `.multilpa_response_bound()` returns `c:p%3Aq`, `c:x%2C y`, and 0 of 2 match the estimation names `c:p:q`, `c:x, y`. Skipped regression test in the new file: it fails on HEAD and passes on HEAD~1. | reported, not fixed. The fix belongs in covariate-inference.R: escape both components in the covariate label builder so every family spells terms one way. plot-gg.R:749 probably has the same mismatch for covariate fits (by inspection only, not run). |
| INF-2 | medium | accessors.R:353-365 and 283-300 | For person centring, `.multilpa_uncentered_indicators()` differs from the original data by 1 ulp in many rows. `.multilpa_alignment_values()` serialises those reconstructed values, so the 15-digit text can differ from the original's CSV text, and a legitimate round trip is refused. | c7.R: column `b`, 30 rows differ, e.g. `0.35483730983513062` vs `0.35483730983513057`. `alignment_values(csv, reconstructed)` is FALSE; `alignment_values(csv, original)` is TRUE. c8.R: 0 of 10 person-centred fits accepted. c9.R: `sensitivity()` errors on HEAD and runs on HEAD~1. Skipped regression test in the new file: fails on HEAD, passes on HEAD~1. | reported, not fixed (needs a tolerance decision). Suggested fix: accept when `abs(supplied - reconstructed)` is at most half a unit in the 15th significant digit plus the reconstruction rounding bound (a few ulps of `max(abs(value), abs(offset))`). That is the precision CSV itself carries, not an origin-proportional band. |
| INF-3 | low | inference.R:1313-1330 | Behaviour change by design: `parameter_inference(fit, data = read.csv(...))` is now refused for every fit. Untransformed fits accepted it before. | c6.R table: `csv` column is BAD for all 5 combinations on HEAD; HEAD~1 gave "ok" for `none/none`. | reported (documented as "checked exactly"; the lead decides). |
| INF-4 | low | extra-indicators.R:542 (outside partition); inference.R:80-85 and 203-206 | Escaping is incomplete and the problem predates Codex. Ordinal category labels and continuous indicator names are not escaped, so a comma or colon still corrupts the tidy columns. | c22.R: ordered level `"mid, ok"` gives term `o:mid`. c17.R: continuous `"x,1"` gives `measurement.mean.profile_1.x` and covariance term `x:1`; `get_results(fit,"profiles")` gives NA mean SE for `x,1` and NA variance SE for both. Same on HEAD~1. | reported, not fixed (spans several files). |
| INF-5 | low | inference.R:1429-1465 | Escaped text reaches the user. `coef()`, `parameter_inference()$term` and `confint(parm=)` show `c:x%2C y` and `c:p%3Aq`, and the covariate family shows raw `c:x, y`. | c2.R output rows 11-12. | reported (cosmetic / consistency). |
| INF-6 | info | lta-bootstrap.R:201-205 (outside) | The transition signature leaves out `covariances`, so the new correlation feature never applies to transition fits. Both sides are built the same way, so this is not a mismatch. | c11.R: full transitions signature has 4 columns (2 means and 2 log-variance ratios). | note |
| INF-7 | info | alignment-service.R:72-80 | Above 7 classes the match is a greedy nearest match, not an optimal assignment (Hungarian). This is as documented. | (reading only) | note |

## Algorithms verified correct, with the check used

| Algorithm | Check | Max discrepancy |
|---|---|---|
| Observed information, VVI, 1 group class (optimHess of the scaled score) | Independent hand likelihood, `numDeriv::hessian` | Hessian relative 1.07e-8; unconstrained covariance relative 1.8e-8 |
| Delta method to natural scale (means, variances, simplex) | `numDeriv::jacobian` of an independent map, `J V J'` | SE absolute 1.64e-9 |
| Two-level VVV + categorical (log-Cholesky, response logits, group logits) | Independent group-class likelihood, numDeriv | Hessian relative 1.50e-7; natural SE relative 2.17e-7; log-likelihood difference 0 |
| Robust CR1 sandwich `G/(G-1) A^-1 B A^-1` (White 1982; Cameron and Miller 2015) | numDeriv per-group scores of the independent likelihood | relative 3.0e-8 (VVI), 7.6e-9 (VVV) |
| OPG `B^-1` (BHHH) | Same per-group scores | relative 5.3e-10 |
| MLR scaling correction `tr(A^-1 B)/q` | Same | 1.017552 vs 1.017552 |
| `.inference_sandwich` / `.inference_cluster_factor` | `sandwich::vcovCL(type="HC0", cadjust=TRUE)` on a clustered glm | absolute 1.4e-17 |
| Covariance charts EEI, VVI, EEE, VEI, EVE, VVE | numDeriv Hessian of the package likelihood in chart coordinates | Hessian relative 1.7e-9 to 6.2e-6; SE relative up to 4.1e-6 (EVE/VVE orientation charts) |
| `boundary = "fix"` conditional information | numDeriv Hessian over the free coordinates with the held logit fixed | covariance relative 4.8e-7; natural SE relative 6.8e-7 |
| `.inference_score_information` (central and five-point; plus1, max1 and absolute step rules) | Analytic Hessian of a normal log-likelihood | all within 1e-7 (test) |
| `.inference_invert` / positive-definiteness refusal | Indefinite and near-singular matrices raise `latents_singular_information` | test |
| Wald intervals (symmetric, log-scale, logit-scale) | Hand formulas | 0 |
| Rubin's rules: total variance, riv, df `(m-1)(1+1/r)^2`, fmi `(r+2/(df+3))/(r+1)` (Rubin 1987; van Buuren 2018 section 2.3) | `mice::pool.scalar(n = Inf)` and hand computation | estimate 0, total 0, df 1.8e-15, riv 2.2e-16, fmi 1.1e-16 |
| `pool_imputations()` end to end | Pooled from `get_results(.,"imputations")` by hand; vcov diagonal | estimate 2.2e-16; SE 0; vcov diagonal 1.4e-17 |
| Label alignment and permutation | Known permutation recovered; relabelled fit's SE table matched to the original's | coefficients identical; SE relative 2.7e-10 (diagonal), 1.9e-10 (full) |
| Translation and scale equivariance (`a -> 1000 a + 1e6`) | Ratios of SEs | mean 1000.001 / 999.997; variance ratio / 1e6 = 1.000002; robust and OPG 1000.00 |
| Seed scoping (`.latents_local_seed`) | Stream identical after a seeded bootstrap and after an error; no `.Random.seed` created in a fresh session; seeded replicates reproducible | exact |
| Ordinal rare category in the bootstrap | 12/12 replicates valid, no column misalignment (levels are kept) | n/a |
| Constancy pattern (`sd <= 0`) in this partition | Only `.multilpa_standardization` (accessors.R:534-537) uses it. `stats::sd()` of constant vectors (n = 1..200, values 0.1 and 1/3) is exactly 0, and `multilpa()` refuses constant indicators earlier (`latents_bad_data`). Not affected. | 0 nonzero out of 200 |

## Tests added

File: `tests/testthat/test-review-1026-inference.R`, 11 tests:
1. observed, robust (CR1) and OPG covariances match independent numDeriv computations;
2. `boundary = 'fix'` inverts the Hessian of the free coordinates only;
3. the CR1 sandwich composition equals `sandwich::vcovCL(type = 'HC0')`;
4. score-difference information recovers an analytic Hessian under both schemes (plus the invert refusals);
5. Rubin's rules match hand computation and `mice::pool.scalar`;
6. alignment recovers a known permutation and relabelling leaves the SE table invariant (including equal means separated by spread);
7. standard errors are translation invariant and scale equivariant;
8. a seeded bootstrap restores the caller's stream and is reproducible;
9. starting variances keep profiles in rows with one indicator;
10. (skipped, INF-1) covariate fits keep response standard errors for reserved category names;
11. (skipped, INF-2) a CSV round trip of a person-centred fit's data still aligns.

Suggests packages are guarded with `skip_if_not_installed()`; seeds are local.

Result (final run): `tests=11 expectations=64 failed=0 skipped=2 errors=0`, about 3.4 s wall clock.

The two skipped tests were checked with the `skip()` lines removed. On HEAD
both fail: SE `anyNA` is TRUE, and the alignment check raises
`latents_bad_inference_data`. On the HEAD~1 tree both pass (2/2). Remove the
`skip()` once each fix lands.

## Related existing test files run (all on the current tree, NOT_CRAN=true)

| file | tests | expectations | failed | errors |
|---|---|---|---|---|
| test-inference.R | 7 | 54 | 0 | 0 |
| test-tidy-inference.R | 10 | 58 | 0 | 0 |
| test-inference-cache.R | 4 | 21 | 0 | 0 |
| test-categorical-inference.R | 6 | 27 | 0 | 0 |
| test-alignment-service.R | 12 | 63 | 0 | 0 |
| test-assignment-alignment.R | 8 | 26 | 0 | 0 |
| test-audit-nb-alignment-repairs.R | 6 | 20 | 0 | 0 |
| test-audit-inference-repairs.R | 17 | 85 | 0 | 0 |
| test-math-inference-audit.R | 10 | 41 | 0 | 0 |
| test-bootstrap-inference.R | 15 | 66 | 0 | 0 |
| test-pool-imputations.R | 8 | 78 | 0 | 0 |
| test-sensitivity.R | 13 | 65 | 0 | 0 |
| test-review-1015-accessors.R | 12 | 34 | 0 | 0 |
| test-review-1015-inference.R | 10 | 39 | 0 | 0 |
| test-review-1020-inference.R | 6 | 77 | 0 | 0 |
| test-tidy-accessors.R | 19 | 248 | 0 | 0 |
| test-staged-identity.R | 5 | 54 | 0 | 0 |

No regressions. Nothing in the partition was edited, and other reviewers'
concurrent edits to other R files were loaded as they stood.

## Not checked, and why

- The claim in the `parameter_inference()` docs that Mplus MLR applies no
  `G/(G-1)` factor: Mplus is not installed, so I could not verify it.
- Barnard-Rubin small-sample df: not implemented by design (the docs say
  "large complete-data sample"), so there is nothing to verify.
- The plot-gg.R covariate-fit interval mismatch implied by INF-1: inferred
  from the code, not run (outside the partition).
- Greedy alignment for more than 7 classes: not stress-tested.
- Bootstrap percentile coverage: not studied (that would be a Monte Carlo
  study, out of time budget).
