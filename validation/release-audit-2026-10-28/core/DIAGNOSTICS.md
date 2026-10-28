# Focused diagnostics review — 2026-10-28

Read `R/diagnostics-checks.R`, diagnostics assembly and replication reporting in `R/report.R`, classification/entropy adapters, and bivariate residual handling. Exercised current working-tree code using `pkgload::load_all()`. No production edits.

## DIAG-AUDIT-1 — P2: one profile plus noise is reported as nothing to classify

**Location:** `R/diagnostics-checks.R:220`; consequent early return at lines 264–269.

The adapter uses `x$n_profiles` as the number of classified alternatives. A noise fit has another alternative, represented by class 0 in its classification table. For `n_profiles = 1, noise = TRUE`, the adapter reports one class and the check function suppresses the class-size checks and returns `n/a` for entropy and average posterior with the incorrect explanation **“one class: nothing to classify”**.

Reproduced with a synthetic two-dimensional sample of 250 normal observations and 50 uniform background observations (seed 821). `lpa(..., n_profiles=1, noise=TRUE, n_starts=2, seed=1)` yields:

- Entropy table: **2** alternatives and relative entropy **0.9442635**.
- Classification table: **252** assigned to profile 1 and **48** assigned to noise/class 0.
- Average posteriors: **0.9904360** and **0.9877783**.
- Checks: entropy and lowest average posterior both **n/a**, claiming nothing is classified; class-size checks omitted.

Use the number of actual classification posterior columns, including noise, when deciding whether the level classifies. The existing three-profile noise test only checks that residuals are unavailable and misses this case.

## DIAG-AUDIT-2 — P3: weighted modal-versus-estimated gap mixes sample and weighted population shares

**Location:** `R/diagnostics-checks.R:287-290`.

The underlying classification table deliberately retains raw `n_modal` and `proportion_modal` and correctly weights `estimated_proportion`. The new check compares those differently weighted quantities without explaining that distinction. Its displayed gap can be entirely due to the sampling weights, even when classification is certain.

Reproduction (seed 229): 150 well-separated observations in each of two classes; weights 1 in one class and 9 in the other. Every maximum posterior is exactly 1 to printed/machine precision. Raw modal proportions are **[0.5, 0.5]**; weighted estimated proportions are **[0.1, 0.9]**; independently computed weighted modal proportions are also **[0.1, 0.9]**. Check key `modal vs estimated proportions` reports **“largest gap 0.400”**, while the consistently weighted modal-versus-soft gap is **0**.

This check is **status `info`, with no threshold or warning**. It does not alter estimation and is not a false statistical rejection. Either compare consistently weighted proportions for this diagnostic, keeping raw counts in the underlying table, or explicitly label the comparison as unweighted sample versus weighted estimated shares. The table itself is not identified as defective.

## Checks passed

- Weighted categorical fit: 240 observations including 80 zero-weight rows; classification totals and entropy `n_units` both correctly count **160** contributing observations.
- Every categorical pair is returned for both profiles; under weights `statistic`, `df`, `p_value` and `p_adjusted` are unavailable. The residual check correctly returns `info`, criterion `not tested`, with an explicit sampling-weight explanation.
- Reviewed the mixed/count residual path: it constructs cross-kind pairs separately from continuous and discrete pairs and reports residual correlations without unsupported p-values. A weighted continuous/categorical/count synthetic probe additionally exercises the six expected pair/profile rows.
- Reviewed finished-versus-screened start counting: replication uses final-stage rows and records screened rows separately; no new defect confirmed there.
- Noise residuals are correctly unavailable; the confirmed noise error is in classification checks only.

## Reproduce

Run from the repository root:

```sh
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 Rscript validation/release-audit-2026-10-28/core/probe-diagnostics.R
```

Output is saved in `probe-diagnostics.log`. These are focused synthetic checks, not an exhaustive validation of every family/empty-cell configuration or the asymptotic calibration of residual tests.
