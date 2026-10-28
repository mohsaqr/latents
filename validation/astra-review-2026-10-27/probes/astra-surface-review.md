# Surface, validation and diagnostic review

Review date: 2026-10-27. Scope: committed changes from `24773e6` to `5596a15`, plus the pre-existing uncommitted candidate, especially condition forwarding, input validation, accessors, tables/printing, diagnostics, sensitivity and TNA. No package source or package test files were changed by this reviewer. Synthetic probes are `/tmp/astra-surface-probe{,2,3,4,5}.R` and matching `.log` files. Probes 4–5 used the parent's fresh installed package `/private/tmp/astra-review/lib`; earlier probes used `pkgload::load_all('.')`. All fitting probes used one BLAS thread.

## Findings requiring correction

### SURF-1 — P2: argument binding can overwrite a valid argument with another argument

**Introduced in committed 0.97.4; remains in current tree.** `R/condition-calls.R:30–35` uses the caller's symbol verbatim, but uses `.argument_name` for expressions. These namespaces overlap. A valid data frame named `.vars` is overwritten when the `vars` expression is bound to `.vars`; both data and vars then reference the latter value.

Reproduction (probe 2):

```r
set.seed(1027)
d <- data.frame(id = rep(1:20, each = 5), x = rnorm(100), y = rnorm(100))
.vars <- d
multilpa(.vars, vars = c("x", "y"), id = NULL,
         n_profiles = 1, n_starts = 1, seed = 1)
```

Actual: `latents_bad_data`, “`data` must be a data frame, but `.vars` is a character vector of length 2.” Direct `multilca()` fails identically. The `lpa()` wrapper happens to avoid the collision because it forwards its own formal symbols. Any user symbol that collides with generated argument names is at risk; this is not a restriction documented on data frame names.

**Repair:** allocate unique internal binding names independently of caller names, while retaining caller labels separately for messages. Test collisions against generated names and the reserved function binding. Keep condition calls bounded.

### SURF-2 — P2: `multilca()` evaluates argument expressions twice

**Committed regression for `vars` and `...`; the uncommitted early data-frame check additionally exposes `data`.** `R/multilca.R:47–65` forces `vars` and `list(...)`, then `.latents_bind_values()` evaluates the original syntax again in the caller frame (`R/condition-calls.R:30`). Under current code the early data check also forces `data` before the second evaluation.

Probe 3:

```r
calls <- 0L
make_data <- function() { calls <<- calls + 1L; d }
multilca(make_data(), vars = c("x", "y"), id = NULL,
         n_profiles = 1, n_starts = 1, seed = 1)
calls  # 2; lca(make_data(), ...) gives 1
```

Probe 5, using the installed candidate, independently counted two executions of `seed = next_seed()` where `next_seed()` increments a counter. That path already exists in committed 0.97.4 because `list(...)` was present before the early data check. Expensive data generation is repeated; random or stateful expressions can yield a different model input/control from the value initially validated. Ordinary R argument promises are expected to evaluate once.

**Repair:** bind values obtained from the existing formal promises and evaluated dots, preserving missing-argument behavior and the original call solely as metadata. Regression tests should count evaluations for data, vars, seed and other forwarded controls, not only condition-call size.

### SURF-3 — P2: new cross-level diagnostics reintroduce zero-weight groups into class counts

**Uncommitted diagnostic redesign.** `R/diagnostics-checks.R:179–196` removes missing posterior rows but does not remove zero-weight rows; it normalizes weights to *all* posterior rows and tabulates all modal assignments. This contradicts the candidate's explicit convention that zero-weight units are absent, and produces wrong reported effective class sizes.

Probe 4 fitted a restricted cross-level model to 25 groups, 4 observations per group, with weight 1 for 20 groups and 0 for 5 groups:

```r
set.seed(1027)
d <- data.frame(id = rep(1:25, each = 4), x = rnorm(100), y = rnorm(100),
                w = rep(c(rep(1, 20), rep(0, 5)), each = 4))
f <- multilpa(d, c("x", "y"), id = "id", n_profiles = 2,
              n_group_classes = 2, family = "restricted_cross_level",
              weights = "w", n_starts = 2, seed = 1)
latents:::.latents_level_classification(f, "groups")
```

Actual effective counts: 4.068015 and 20.931985, summing to 25; modal counts 5 and 20, also summing to 25. Positive support is 20 groups. Correct effective counts under the established normalization sum to 20 (approximately 3.254412 and 16.745588). The fit likelihood itself is not implicated by this probe. A class occupied only by excluded groups could also evade the diagnostic's “no assigned group” flag.

**Repair:** filter positive weights before posterior classification, compute normalizers/counts on positive support, and test that appending zero-weight groups leaves diagnostics invariant.

### SURF-4 — P2: mixture-regression diagnostics lose nesting and class-count metadata

**Uncommitted diagnostic redesign.** `R/diagnostics-checks.R:203–217` reads class counts from `x$n_classes`/`x$n_group_classes`, but `latents_mixture_regression` stores them under `x$spec`. It therefore treats the number of individual classes as unknown and group classes as one.

Probe 2 fitted `mixture_regression(y ~ x, d, id = "id", n_classes = 2, n_group_classes = 2, n_starts = 2, seed = 1)`. The fit table correctly reports two group classes and a two-level model, but diagnostic checks contain only `model` and `individuals`, omitting the group level entirely although group entropy and posteriors are available (`x$expectation$rho`; the fit table's `group_entropy`).

The same defect means `class_level = "group"` reports all classification checks as `individuals`. Probe 3's one-class regression reports average posterior 1.000 as classification information rather than the intended `n/a: one class: nothing to classify` behavior, because the one-class branch cannot detect its class count.

**Repair:** derive levels/counts from each family's actual metadata; use the correct group-posterior table for two-level regression; test ordinary, group-classified and two-level regression, with one- and multiple-class cases. Growth's early return also omits any second-level class checks for hierarchical growth; this was identified statically, not independently fitted in this review.

### SURF-5 — P2: new sensitivity plot describes unconverged initial states as replicated maxima

**Uncommitted plotting addition.** `R/sensitivity.R:197–206,220–221` checks only whether log likelihood is nonmissing. It does not inspect `converged`, despite that column being available. Every finite unfinished run is labelled a “best maximum” or “lower maximum”.

Probe 4:

```r
f <- lpa(d, c("x", "y"), n_profiles = 2, n_starts = 1, seed = 1)
s <- sensitivity(f, seeds = 1:3, max_iter = 0, n_starts = 1)
p <- plot(s)
s[, c("converged", "iterations", "best")]
p$labels$title
```

Actual: all three rows are `converged = FALSE`, `iterations = 0`, `best = TRUE`; plot title is “3 of 3 seeds reached the best maximum”, and every point is “best maximum”. No EM step was taken. The table exposes nonconvergence, but the new visual drops that distinction. This is separate from the intentionally unresolved decision about deterministic Ward starts.

**Repair:** show unfinished refits explicitly and count only converged refits as maxima; alternatively describe finite log-likelihood values without claiming maxima. Preserve failed and unfinished seeds on the axis. Test all-unconverged and mixed converged/unconverged tables.

## Other observations and compatibility notes

- **Diagnostic coercion documentation:** newly supported additive/regression/growth/cross-level diagnostics have no entropy table in their catalogue. Consequently `as.data.frame(diagnostics(additive_fit, plots = FALSE))` returns the replication table, because it becomes the first available table. This does not break an old supported family, but contradicts the unqualified “returns the entropy table” documentation in `R/report.R`/`man/diagnostics.Rd`. Either document the family-specific primary table or create a consistent table for all families.
- **Unnecessary cross-level limitation:** the new helper states that weighted cross-level fits store no map from rows to groups, but the reviewed fit stores `sufficient_statistics$index`, which is the row-to-group map used during fitting. Individual-level diagnostics currently return n/a rather than wrong weighted values, so this is a capability/documentation improvement rather than a release blocker.
- **Modal versus estimated share:** new diagnostics subtract an unweighted modal proportion from a weighted estimated proportion. This is literally what its fields contain, but under sampling weights the gap mixes sample composition and population classification. Label the distinction or compute both on the same weighting convention; do not present the gap as solely assignment uncertainty.
- **Deterministic Ward sensitivity:** HANDOFF explicitly records the unresolved design question. `n_starts = 1` can return identical fits across seeds because the start is deterministic. Preserve the owner's pending decision; do not silently switch all starts to random. Explain that seed stability with a deterministic first start is not independent evidence of global optimality. This review's zero-iteration sensitivity reproduction also demonstrates identical seed outputs.
- **Missing unsupported quantities:** documentation promises n/a with a reason for family-unsupported quantities; several such diagnostics omit OCC or other unsupported classification rows entirely. This is a consistency issue, distinct from the concrete loss of an available regression group level above.

## Positively verified behavior and coverage

- Checked committed accessor changes: exact objective-data identity, CSV-compatible alignment, centered/scaled reconstruction bounds, typed and half-observed pair reconstruction, delimiter-safe inference term matching, and one-indicator covariance-start reshaping. No new concrete defect was established in these paths in this pass; numerical inference is delegated to the inference review.
- Read TNA sequence changes: classified rows are obtained through the shared sequence table; numeric identifiers are matched by native value rather than default printed strings. General transition fits still explicitly refuse unsupported aggregation. A plain LPA without ordering raises the intended ordering error; it does not accidentally return NULL through the new default method.
- The invalid-column `multilca()` probe returns `latents_bad_data` with NULL `conditionCall`, confirming the original oversized-call repair on this representative path.
- Emulated a saved diagnostics object from before checks existed by removing `$checks` and `$replication`: it prints and its default coercion equals its entropy table. Replication is reconstructed through the fallback. The legacy one-class print still contains NA entropy, but that is preserved old behavior, not a newly introduced failure.
- Existing diagnostic tests cover standard families, screened/final starts, positive sampling weights, missing/noise residuals, unsupported plotting, an indicator named `profile`, and failed sensitivity seeds. They do not cover zero sampling weights in the posterior-based diagnostic family, nested regression diagnostics, argument-name collisions, repeated expression evaluation, or unconverged sensitivity plotting. This explains why the reproduced defects can coexist with a passing general suite.
- The new no-cutoff status design follows the recorded owner preference: entropy, AvePP, OCC, class share and replication count are informational; only concrete convergence, bounds, empty assignment and residual-test facts are flagged. Do not restore discarded rule-of-thumb thresholds as a repair.

## Validation limits

This is an independent focused code/probe review, not a full visual inspection or a fresh external statistical-equivalence run. Plot metadata was inspected programmatically; no plot implementation was edited. Full package, build, harness and CI validation are handled by the parent reviewer. No assertion of defect-free behavior is made for unexecuted combinations.
