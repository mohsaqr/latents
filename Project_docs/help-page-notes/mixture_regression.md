# Implementation notes removed from `?mixture_regression` (2026-11-03)

- **seed**: `NULL` advances the caller's random stream; an integer seed
  restores `.Random.seed` afterwards, but R does not store a cached
  Box-Muller normal draw in `.Random.seed`, so that hidden draw cannot be
  preserved.
- **min_variance**: in a growth mixture the floor is applied in the
  standardized coordinates the model is fitted in (outcome and predictors
  centred and scaled), so it does not depend on the units of time.
- **vcov_type = "robust"**: CR1 factor `G / (G - 1)`, clustered on the
  top-level unit, or on `id` for a single-level fit given `id`. Weighted fits
  reject `"observed"` and `"opg"`.
- **weights**: the information criteria use the weighted log likelihood and
  the entropy term of ICL is weighted. A zero-weight row whose factor level
  (or ordinal outcome category) no positive-weight row has is dropped with a
  `latents_rows_dropped` warning, since the model has no parameter for it.
- **residual_correlation**: the model without the process is fitted first
  from every start; the process is then estimated jointly with the other
  parameters. The bound warning is `latents_random_boundary`; without
  `random` the argument raises `latents_bad_argument`. Brownian motion with
  negative times raises `latents_bad_data`.
- **ordinal identification**: a single-level mixture is rejected when it has
  more free parameters than `P (C - 1)`, the number of free category
  probabilities in `P` distinct covariate patterns (for example two classes
  and one binary factor). The thresholds replace the formula's intercept, so
  the intercept must not be removed. An ordinal threshold gap at zero (a
  class giving no row some middle category) is a boundary; that gap and the
  thresholds above it have no standard error.
- **Conditions** (formerly a section): `latents_bad_argument` (inconsistent
  specification, `common` naming no term); `latents_bad_data` (outcome
  outside the family's support, covariate varying within a group, factor
  with one value among positive-weight rows); `latents_not_identified`;
  `latents_missing_data` (`missing = "error"`); `latents_rows_dropped`
  (warning); `latents_no_valid_start` (every start degenerated);
  `latents_unconverged` (warning); `latents_degenerate_start` (warning);
  `latents_separation` (warning, extreme logits or log means);
  `latents_failed_starts` (warning, growth-mixture starts marked "failed");
  `latents_boundary` (warning); `latents_random_boundary` (warning);
  `latents_singular_information` (warning, some standard errors missing).
- **n_starts**: growth-mixture random starts first run the model without
  random effects for up to 200 iterations, and a second extra start
  ("trajectories") is added beside the residuals start (`kind` in the starts
  table). Growth EM runs at `max(tol, 1e-6)`; the quasi-Newton finish uses
  `tol`. Random effects are integrated analytically; estimation is EM then
  quasi-Newton with analytic scores.
- **weights**: scaled to sum to the number of units of positive weight.
- **random intercept with single-row units**: refused as
  `latents_not_identified`.
