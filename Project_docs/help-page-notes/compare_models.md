# Implementation notes removed from `?compare_models` (2026-11-03)

- **Conditions** (formerly in Value): `latents_bad_argument` (incomparable
  fits, fewer than two models, duplicate labels, family grid, candidates
  named together with `n_profiles`/`n_group_classes`, those arguments without
  a grid, `cutoff` outside (0, 1], family fits, regression grids, mixture
  regressions under `test = "lrt"`); `latents_unknown_candidate`;
  `latents_failed_candidate`; `latents_incomparable_models`;
  `latents_bad_nesting` (also `latents_singular_restriction` for a singular
  covariance restriction); `latents_reversed_likelihood`;
  `latents_no_converge`; `latents_unsupported_weights`,
  `latents_unsupported_prior`, `latents_unsupported_noise`; a variance at
  `min_variance` warns `latents_boundary`.
- **smallest_share**: relative to all observation mass including a noise
  component. The values are those of the grid's candidates table.
- **Nesting detail**: "covariances of chosen pairs (`local_dependence`)
  within more pairs or a full covariance, fewer joint item pairs"; the first
  model's covariate design must be recovered by setting the second's added
  coefficients to zero; the restrictions listed are regular at a generic
  null.
