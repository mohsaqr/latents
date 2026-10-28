# Core release fixes — 2026-10-28

## Changes

- REL-4 (`R/additive.R`): compute additive posterior probabilities with row-max-shifted exponentials normalized by their row sum. Keep the existing stable log likelihood calculation. A large negative group log density can no longer round away the posterior normalizer when its log-sum-exp is subtracted. This also protects the final collapsed-class check, which was the failing call in the regression fixture.
- REL-6 (`R/diagnostics-checks.R`): count actual classification posterior columns for profile/class/transition fits, including the noise alternative. One Gaussian profile plus noise now receives the entropy, class-size and average-posterior checks that its tables support. An ordinary one-profile fit still reports classification as unavailable.
- `tests/testthat/test-release-audit-core.R`: two regression cases, 20 expectations. Excluded sentinel groups must preserve likelihood, means, both variance components, class weights and contributing posteriors; all reported posterior rows remain normalized. Noise diagnostics must agree with the entropy/classification tables, and their residuals remain unavailable.

No changes to the weighted modal-proportion comparison, sensitivity behavior, package version, NEWS or git history.

## Validation

`before.log`: run with the new regression tests and old production source. One additive error (`posterior rows must sum to one`) and four failed noise-diagnostic expectations.

`after.log`: the same regression file with the fixes; all 18 expectations passed.

Warning-quality followup: `inspect-warnings.R` / `warnings.log` record the actual warnings from each fixture. The baseline and augmented additive fits each emit only `latents_boundary` (one between-group variance estimated at zero); the noise and ordinary one-profile fits emit no warnings. The new tests now explicitly expect the two `latents_boundary` warnings and use no blanket warning suppression. `after-warning-check.log`: final version **20 passed, no unexpected warnings**, elapsed **1.372 s** (1.280 s user, 0.034 s system), excluding package loading. Production source was unchanged by this followup.

`targeted.log`: existing additive-engine, additive-families, additive-fit and review-1027-diagnostics tests, together with the new regression file. All five files passed, with no failures or skips. Run with `NOT_CRAN=true` and BLAS/OpenMP limited to one thread.

Commands were run from the repository root using `testthat::test_file()` after `pkgload::load_all()` for the before/after proof, and `testthat::test_local(filter="^(additive-(engine|fit|families)|review-1027-diagnostics|release-audit-core)$", reporter="summary", stop_on_failure=TRUE)` for the existing checks. `git diff --check` passed for the changed production files and new regression file.
