# Supplemental public API audit — 2026-10-28

Review-only follow-up of recent changes in `R/get-data.R`, `R/accessors.R`, `R/compare-models.R`, `R/report.R`, `R/grid-views.R`, `R/input-checks.R`, `R/condition-calls.R` and `R/pool-imputations.R`. Current source loaded with `pkgload::load_all('.')`; no production files changed.

## Confirmed P2: comparison wrappers silently omit a requested model

`R/compare-models.R:115-120` reconstructs the position of fit-valued arguments named `n_profiles`, `n_group_classes`, or `cutoff`, to preserve their historical use as model labels. It counts the arguments in `sys.call()`, which leaves a forwarding `...` as **one** expression, although `fits <- list(...)` has already expanded it into several models. Consequently `order(positions)` is shorter than `c(fits, restored)`, and subsetting discards models without warning.

Reproduction with three successfully fitted models on the same simulated outcome (one-profile VVI, two-profile VVI, and two-profile EEI):

```r
compare_models(first = a, second = b, n_profiles = c) # all three models
wrapper <- function(...) compare_models(..., n_profiles = c)
wrapper(first = a, second = b)                       # silently only a and b
```

The restored model `c` is omitted from both the comparison rows and the attached `models` list. The BIC weights are consequently normalized over the wrong candidate set. A generic wrapper `function(...) compare_models(...)` also fails to preserve direct-call behavior when passed the formerly valid `n_profiles` model label.

Evidence: `probe-public-api.R` and `probe-public-api.log`. This tests the real public functions with synthetic data and requires three rows directly and two through the faulty wrapper. The issue is narrower than ordinary `compare_models(a, b)`: it needs a forwarding `...` plus a model supplied under one of those three restored labels.

Recommended repair: compute ordering from a call representation with forwarded dots expanded, or restore the models without a position vector whose length can differ from the fit list. Add direct-versus-forwarded parity tests for all three restored labels, multiple forwarded fits, and labels supplied through the wrapper itself.

## Other inspected paths and successful probes

- `get_results(fit, 'all')` for an interior one-profile fit returns data-frame tables; the new replication table is available and the trivial group classification is excluded from the classification catalogue as intended.
- `diagnostics(fit, plots = FALSE)` and direct `get_results(fit, 'classification')` agree. The entropy table still contains the trivial group row by the documented separate rule.
- New pooled `coef()`/`vcov()` names align, and `confint()` reproduces the stored pooled intervals at the stored confidence level for two independently perturbed completed data frames. The earlier fresh independent pool-imputations equivalence suite also passed 92 comparisons.
- Read the changed input-validation, bounded-condition, data reconstruction, grid classification and reporting code. No further reproducible numerical or data-loss defect was confirmed within this targeted review. Existing condition-call tests cover name collisions and data redaction; they were not redundantly rerun here while the parent runs the full package suite.

Limitations: targeted review of recent changes and small public API probes, not exhaustive coverage of historical serialized fits, every unusual column label, plotting output, all mixed-family grids, or every valid wrapper. The already-known sensitivity selection issue was not counted as new. No rendering changes were made.
