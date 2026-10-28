# Regression inference repairs — 2026-10-28

Fixed REL-1 and REL-3 in `R/mixture-regression-inference.R`. No version, NEWS, git history, or unrelated source files were changed by this agent.

- Shared Gaussian/NB dispersion packing extracts a scalar with `[[1L]]`, so incidental class names cannot concatenate into `log_sigma.all` or `log_dispersion.all`. Public coefficient/covariance names and class SD/dispersion estimates and errors remain available.
- `.mixture_inference()` now uses `.mixture_carrying_spec(fit$spec)` consistently. Its chart, score rows, rank checks, cluster indexing and CR1 unit count all see the same positive-weight support as fitting. The original fitted object keeps its complete data and predictions for excluded units.

New tests: `tests/testthat/test-release-audit-regression.R` (58 passing expectations, no warnings/skips). Covers one and two classes for both shared-dispersion families, shared coefficients, excluded individual rows, observation-level classes with clustered inference, group-level classes with membership predictors, and two-level classes with both membership designs. Group tests use unequal positive weights and an excluded group that sorts before the retained groups, exercising reindexing. They compare full public covariances, coefficient tables, likelihoods, retained posterior rows and cluster metadata.

The final tests were also run in a separate process against the original HEAD packer and inference implementation: **25 failures and four singular-information warnings**, confirming that the tests detect the original defects (`before-final-tests.log`). No original implementation was restored on disk.

Validation:

- New file: **58/58** expectations (`after.log`, `timing.log`, `timing.csv`).
- Six existing files: **257/257** expectations, no warnings/failures/skips (`existing-tests.log`): mixture-regression, mixture-regression-labels, mixture-regression-prediction-design, review-1020-zero-weights, review-1026-mixreg and review-1026-mixreg-r2.
- Normal unmodified harness command `Rscript equivalence/harness/run.R flexmix`: **433/433 quantities, 29,933 values**, zero warnings/skips/uncovered outputs, source unchanged through the run (`flexmix/REPORT.md` and CSVs). No in-memory patch was used for this pass.
- `git diff --check` passes for the owned source and test files.

Commands used `NOT_CRAN=true` for package tests and `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1`; each process loaded current source. Root coordinates the final complete suite and release checks after all agents finish.
