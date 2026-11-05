# Implementation notes removed from `?r3step` (2026-11-03)

- **x**: refusals are `latents_unsupported_three_step`; `level = "groups"` on
  a trajectory model, `by_group_class` outside its case, and the Markov
  covariates on a non-transition fit are `latents_bad_argument`; a covariate
  inside the trajectory regression is `latents_bad_covariate`.
- **data**: too little identifying information warns that alignment cannot
  be verified.
- **vcov_type = "robust"**: CR1 factor `G / (G - 1)`; the group score
  contributions sum to zero at the estimate, hence the requirement of more
  groups than coefficients (`latents_too_few_groups`, "observed" remains).
- **Value**: a coefficient with zero standard error has `NA_real_` statistic
  and p-value.
