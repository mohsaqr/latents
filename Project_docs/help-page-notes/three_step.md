# Implementation notes removed from `?three_step` (2026-11-03)

- **x**: refusals are `latents_unsupported_three_step`; `level = "groups"`
  on a trajectory model is `latents_bad_argument`; `variance` with another
  method or a categorical outcome is `latents_bad_argument`.
- **data**: a frame with too little identifying information warns that
  alignment cannot be verified.
- **Standard errors**: one outer product per independent unit, CR1 factor
  `G / (G - 1)`; contributions sum to zero at the estimate so the rank is at
  most `G - 1`; one group or as many groups as classes raises
  `latents_too_few_groups`; a class carried by one unit (Kish effective
  number below 1.5) warns `latents_small_classes`.
- **ML**: log likelihood sum_i log sum_t pi_t D[t, W_i] phi(z_i; mu_t,
  sigma_t), maximised by EM then Newton steps; equal-variance fits also use
  bounded BFGS from quantile-based initial values; accepted solutions must
  converge, have a small score and a positive-definite information matrix;
  twenty random starts with a fixed internal seed; the caller's
  `.Random.seed` is restored (Box-Muller cache excepted). Categorical ML
  with BCH shares out of range: maximum on the bound, `latents_boundary`
  warning; BCH out-of-range share warns `latents_bch_out_of_range`.
- **step_one = "estimated"**: adds G V1 G', V1 the step-one covariance on the
  unconstrained scale, G the derivative of the step-three estimates through
  the posteriors.
- **Value**: a contrast with zero or missing standard error has `NA_real_`
  statistic and p-values; the overall test is on the baseline-category
  logits.
