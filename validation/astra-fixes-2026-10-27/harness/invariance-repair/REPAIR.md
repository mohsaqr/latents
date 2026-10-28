# Invariance coverage repair — A13 follow-up

Only `equivalence/harness/suites/invariance.R` was changed. Existing replication likelihood mappings were preserved. No package runtime code, comparison tolerance, exemption or ledger rule changed.

The zero-weight padding case fitted 40 additional persons in four groups but supplied expected posterior/assignment rows only for the original 150 persons and 15 groups. Eleven fields were partially compared. The expected tables now include every row returned by the padded fit.

New helper `.expect_excluded_gaussian_groups()` derives the additional outputs from the original fit's means, diagonal variances, profile shares and group-class shares, using only base normal densities and Bayes' rule. For each group class it integrates profile densities, multiplies those likelihood contributions across the group, normalizes group-class probabilities, then marginalizes conditional profile probabilities. The new rows include person posteriors/modal flags; group size, raw log density, posterior and modal flag; and assignment classes, uncertainty and profile posteriors. The zero weights continue to leave the fitted objective and aggregate classification summaries unchanged. No output from the padded fit or package E-step/prediction helper is used to compute these expectations.

## Verification

Strict selected suite command:

```sh
OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NOT_CRAN=true LATENTS_EQUIVALENCE_REQUIRE_ALL=true LATENTS_EQUIVALENCE_OUTPUT=validation/astra-fixes-2026-10-27/harness/invariance-repair Rscript equivalence/harness/run.R invariance
```

Result: **826/826 quantities agree, 39,232 values compared, zero partial or uncovered outputs, zero warnings, exit 0**, 7.93 seconds. `run.log`, `report.csv`, `coverage.csv`, `suites.csv`, `SESSION.txt` and `REPORT.md` retain the results.

`verify-oracle.R` provides an additional independent exhaustive check: enumerate every one of the `2^10` row-profile paths under each of two group classes, separately for each of the four added groups. Against this direct joint enumeration, the factored helper's largest log-density discrepancy is **1.136868e-13**, and its largest person/group posterior discrepancy is **1.110223e-16**. The same script verifies the old expected table has exactly 11 partial fields, while the corrected expected table has none and every comparison agrees. Run it from the repository root with one BLAS thread; `oracle-fixture.rds` contains the synthetic baseline/padded data and fits.

The affected case's ledger changes from **65 compared / 28 exempt / 11 partial** to **76 compared / 28 exempt / zero partial**. Existing exemptions are unchanged. `git diff --check -- equivalence/harness/suites/invariance.R` passes.
