# Additive and cross-level family review, 1–2 October 2026

Reviewed the current sources in `R/additive.R`, `R/additive-inference.R`,
`R/additive-methods.R`, `R/additive-compare.R`, and `R/cross-level.R`, their
`multilpa()`/enumeration/bootstrap dispatchers, and the existing additive,
cross-level and sampling-weight tests. No git mutations were performed.

## Confirmed issues and fixes

| Issue | Correction | Regression evidence |
| --- | --- | --- |
| Additive bootstrap LRT compared means/scatter but omitted group sizes; supplied data were checked against means only. | Both fitted and supplied data must reproduce group IDs, sizes, means and within scatter. | Changing two ratings by opposite amounts preserves a group mean but is refused; adding a row equal to a group's mean preserves means/scatter but changes its size and is refused. Both previously passed the checks. |
| One-class additive and dispersion fits were accepted as different nested models despite being the same model. | Refuse comparisons between any two one-class group-class fits. | The actual additive/dispersion fits formerly produced an LRT and now raise `latents_bad_argument`. |
| Invalid family confidence levels generated unusable bounds or recycled vectors; invalid differentiation steps reached numerical linear algebra. | Validate finite scalar `0 < level < 1` centrally in the parameter table and finite positive scalar `step` centrally in inference. | Invalid levels through `parameter_inference()`, `confint()` and `summary()`, and invalid steps through `parameter_inference()`, are refused. These tests failed before the fix. |
| Additive bootstrap controls could reach allocations/refits with nonfinite or invalid values. | Validate integer iteration/start counts and finite positive tolerance before simulation. | Sixteen invalid-control cases cover `iter`, `n_starts`, `max_iter` and `tol`. |
| A weighted restricted cross-level fit overwrote an indicator called `.sampling_weight` with temporary row weights. | Pick a unique internal weight-column name using `make.unique()`. | Renaming `y1` to `.sampling_weight` preserves fitted means and likelihood, which also agrees with an independent raw-data formula. |
| Cross-level counts exceeding R's integer range reached `as.integer()`. | Reject overflowing profile/class/start/iteration counts before conversion. | All four count arguments reject `.Machine$integer.max + 1`. |
| Missing entries in enumerated between-variance restrictions caused an unclassed missing-condition error. | Explicit `anyNA()` validation. | `between_variance = c("equal", NA_character_)` raises `latents_bad_argument`. |
| Cross-level inference refusal said group means were used twice. | State that the likelihood includes ratings and group means computed from those same ratings. | The inference refusal is verified for every cross-level restriction. |

## Numerical coverage

`tests/testthat/test-review-families.R` provides a reproducible synthetic fixture:
90 groups, 630 observations, group sizes 4/7/10, two numeric indicators,
unequal class proportions and unequal group weights. The fixture's structure,
head, summary and actual fitted object were inspected before validation.

| Model | Between variance | Within variance | New numerical checks |
| --- | --- | --- | --- |
| Additive | Varying | Shared | Dense Gaussian likelihood, gradient, information; weighted EM stationarity before Newton |
| Additive | Shared | Shared | Same |
| Dispersion | Shared | Varying | Same; shared-location synthetic data for EM check |
| Additive-dispersion | Varying | Varying | Same |
| Additive-dispersion | Shared | Varying | Same |
| Restricted cross-level | Shared and varying | Shared and varying | All four combinations: direct raw-data working likelihood and full numerical gradient |
| Full cross-level | Shared and varying | Shared and varying | All four combinations: direct raw-data working likelihood and full numerical gradient |

The additive reference builds each group's full stacked covariance
`kron(diag(W), I) + kron(diag(T), 11')` and evaluates its Cholesky density.
It does not use the engine's sufficient-statistic density. At perturbed
parameters, weighted total scores agree with numerical derivatives to
`1e-6`, observed information with an independent numerical Hessian to `1e-4`,
and log likelihoods to `1e-10`, for every restriction.

The weighted EM-only fits converge in 17–34 iterations, remain interior
(smallest between variance 0.154–0.374), increase likelihood monotonically,
and have maximum absolute analytic scores 5.11e-5–1.23e-4. Those scores were
independently validated against the dense likelihood above. Checking the
unpolished EM point prevents Newton from concealing M-step errors.

The cross-level reference independently evaluates
`sum_j w_j log(sum_h omega_h f(mean_j | h) prod_i sum_k pi_kh f(y_ij | k))`
from raw ratings and returned parameters. Restricted prevalences are the
weighted overall profile proportions, rather than the descriptive composition
table. Likelihoods agree to `1e-8`; gradients across every free block
(profile/group means, log variances, profile-prevalence logits and group-class
logits) are below `2e-3` at `tol = 1e-12`. Reported group and individual
posterior rows sum to one to `1e-12` even with weights. Shared-variance blocks
are checked directly.

A weighted one-class additive inference check reconstructs the CR0 sandwich
from the returned weighted group scores and observed Hessian. It agrees with
the inference covariance and `vcov(scale = "unconstrained")` to `1e-10`.
Observed covariance is explicitly refused for weighted fits.

The equal-between cross-level M-step divides by the number of groups, which
is correct here because sampling weights are normalized to sum to that
number. No change was made to this denominator.

## Checks and reproducibility

Baseline and post-fix original suites passed without failures:

```r
pkgload::load_all(quiet = TRUE)
testthat::test_local(filter = "^(additive|cross-level)", reporter = "summary")
```

The added review tests passed without failures or warnings. Their final
output is retained at `validation/review-2026-10-01/test-families.log`:

```r
testthat::test_local(filter = "review-families", reporter = "summary")
```

`lintr::lint()` reported no lints in any of the five reviewed production
files. The root agent is responsible for the current-source full-suite and
package-check evidence and consolidated `HANDOFF.md`, `LEARNINGS.md` and
`CHANGES.md` updates.

## Intentional limitations and next-agent double checks

Cross-level models use a working likelihood because the between indicators
are means of the same ratings used at the individual level. These likelihoods
should be compared only among cross-level fits of the same data. Ordinary
likelihood Wald inference remains explicitly unsupported; this review does
not establish sandwich inference for that working likelihood.

Additive families retain their complete continuous raw-rating contract.
Boundary fits refuse Wald intervals; the existing boundary/KKT, small-class,
classification, parameter-table, natural/unconstrained covariance, start,
enumeration, simulation, bootstrap, print and plot tests remain in force.
Additive bootstrap parameter inference remains explicitly unsupported, and
weighted parametric bootstrap LRTs remain refused by the shared dispatcher.

No new external-program equivalence or repeated-sampling coverage claim was
made. Retained Latent GOLD and simulation evidence should be assessed through
the existing equivalence/validation reports. The independent checks here
establish formulas, optimization and the targeted surface corrections.

The next agent should rerun the two commands above on the complete integrated
worktree, check that the root full-suite/package checks cover these files,
and inspect bootstrap nesting/data identity and the temporary weight-column
fix. No unresolved confirmed defect remains in this review's assigned scope.
