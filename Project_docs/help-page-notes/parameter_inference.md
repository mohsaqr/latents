# Implementation notes removed from `?parameter_inference` (2026-11-03)

- **vcov_type**: the sandwich is A^-1 B A^-1 with B the outer product of the
  per-group scores times the CR1 factor G / (G - 1); these errors are
  sqrt(G / (G - 1)) times Mplus's MLR. OPG and observed information are
  asymptotically equivalent under a correctly specified regular model; OPG
  also needs group scores of full rank. Weighted fits raise
  `latents_unsupported_weights` for any type but "robust"; person-centred or
  person-scaled fits raise `latents_unsupported_inference` for "observed" or
  "opg". Grand centring/scaling treats the sample mean and sd as estimating
  equations.
- **method = "wald" coordinates**: EEI, VVI, EEE, VVV in log variances or
  log-Cholesky coordinates; the ten constrained structures in log volumes,
  determinant-one log shapes and Cayley-transform orientations (VEE and EVV:
  determinant-one log-Cholesky factors). Constraints can make the
  natural-scale covariance singular; the sum constraints on class
  probabilities do too.
- **seed**: R does not expose the Box-Muller cached normal variate, so that
  cache cannot be restored (see base::Random).
- **Covariate fits**: the boundary refusal is for a membership logit above 20
  in coordinates independent of the covariates' units and origin (a slope
  per standard deviation, an intercept at the covariates' means); the fit
  warns `latents_extreme_coefficients`. Bootstrap raises
  `latents_unsupported_inference`.
- **Transition fits**: scores follow the Fisher identity; the observed
  information is the numerical Jacobian of the analytic score. The bootstrap
  without `data` raises `latents_bad_argument`. A no-information transition
  row raises `latents_boundary_fit`.
- **Growth mixtures**: `step` is not used (the information is the derivative
  of the analytic score); the bootstrap refits each resample from `n_starts`
  starts; multilevel growth mixture bootstrap raises
  `latents_unsupported_inference`; `data` must be the fitted data.
- **Value details**: a categorical response probability at or below 1e-6 is
  treated as on its bound; a mixing probability is refused only at zero. The
  table for a `fixed` fit has `n_parameters` rows plus one per set of
  probabilities whose reference category the estimation scale drops;
  `fit_staged()` bootstrap sets attribute `stages_resampled` to 2. Bootstrap
  attributes: `covariance`, `level`, `method`, `iter`, `n_valid`,
  `replicates`, `replicate_number`, `messages`, `structure`, `fixed`;
  `covariance_unconstrained` is NULL. Wald attributes also hold `level`,
  `step`, `method`, `vcov_type`; `scaling_correction` is tr(A^-1 B) / q with B
  without the CR1 factor, NA unless robust/opg.
- **Computation**: two score evaluations per free parameter; the information
  matrix is cached within a session for the same fit, data and `step`, and
  reused by plotting and `get_results()`.
- **Conditions** (formerly a section): `latents_no_converge`,
  `latents_boundary_fit`, `latents_singular_information`,
  `latents_no_free_parameters`, `latents_bad_inference_data`,
  `latents_too_few_groups`, `latents_unsupported_inference` (prior fit under
  Wald; bootstrap for additive, covariate or `fixed` fits),
  `latents_unsupported_noise`, `latents_unconverged` (warning),
  `latents_bad_argument` (invalid level/step/iter/n_starts/max_iter/tol;
  mixture regression), `latents_bad_data`, `latents_not_a_fit`,
  `latents_bootstrap_failed` (fewer than two usable resamples),
  `latents_bootstrap_dropped` (warning).
