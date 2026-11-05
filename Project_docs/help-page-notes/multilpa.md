# Implementation notes removed from `?multilpa` (2026-11-03)

Material cut from the help page when it was rewritten for readers. It is
still true of the code; it is kept here for maintainers.

- **tol**: the count of observed indicator values is each non-missing cell in
  a row of positive weight, with a minimum count of one.
- **select_start**: "likelihood" selects among starts within
  `tol * max(1, n_observed)` of the highest log likelihood, preferring
  converged starts in that set; "converged" falls back to all starts when
  none converged. Remaining ties resolve by start order. The replication
  tolerance `1e-6 * max(1, n_observed)` is stored as `replication_tolerance`
  and does not select a start.
- **acceleration**: SQUAREM extrapolates along the last two EM steps in
  unit-free coordinates (means divided by reference standard deviations,
  variances and covariances by the corresponding products; probabilities
  unchanged); reference standard deviations use observed values in rows of
  positive weight. Each cycle ends with an ordinary EM step; convergence is
  assessed over a complete cycle; `iterations` and `max_iter` count EM steps
  under either setting. For EVE and VVE each EM iteration takes a fixed number
  of warm-started orientation steps (a generalized EM step), which keeps the
  EM map a fixed function as SQUAREM requires.
- **volume/shape/orientation**: EVE and VVE use the minorize-maximize step of
  Browne and McNicholas (2014). Structures other than EEI, VVI, EEE and VVV
  are maximized across all profiles jointly and cannot be combined with a held
  `variances` block. Wald standard errors for all fourteen are computed in the
  structure's free coordinates (log volumes, determinant-one log shapes,
  orientations) and delta-method transformed.
- **scaling**: SQUAREM's search path is equivariant under grand scaling in
  exact arithmetic when variance bounds are inactive or transformed with the
  model; rounding near convergence can change iteration counts and nearly
  tied starts can change the selected labelling. Grand scaling adds the
  density Jacobian to each start's log likelihood. A scalar variance or
  eigenvalue floor is not preserved by unequal rescaling. Standard deviations
  use denominator n - 1. Person-scaled groups observed exactly twice raise
  `latents_atomic_groups`; when every group is like that the fit is rejected
  with `latents_bad_data`.
- **families**: between variances estimated at zero are reached by a boundary
  maximization with a Karush-Kuhn-Tucker check (`latents_boundary`); interior
  fits are finished by Newton steps.
- **weights**: under `missing = "error"` a missing categorical, ordinal or
  count value on a zero-weight unit is not rejected, a missing continuous
  value is. An indicator constant among positive-weight units is rejected
  regardless of zero-weight values. `parameter_inference()` rejects
  `vcov_type = "observed"` and `"opg"` (`latents_unsupported_weights`).
- **cores**: workers are forked on macOS and Linux, a socket cluster on
  Windows and on macOS with Accelerate BLAS. Lost worker results are retried
  once with half as many workers (rounded up), with one `latents_worker_lost`
  warning; starts still lost are recorded as failed. Socket workers retain
  completed chunks, but tasks within a failed chunk may be retried. Under emEM
  both screening and final runs are distributed. `enumerate_classes()` shares
  its pool across candidates and starts only for models without covariates or
  joined categorical `local_dependence` pairs; other engines use one task per
  candidate.
- **noise**: the hypervolume is the smaller of the axis-aligned and
  principal-component-aligned bounding boxes (mclust `hypvol()`); the fit adds
  `noise_probability`, `noise_posteriors`, `n_noise` and `hypervolume`.
- **prior**: the decrease check made under maximum likelihood is not applied
  because a MAP step can lower the unpenalized likelihood. The resolved
  hyperparameters are stored as `prior_parameters` in centred indicator units.
- **direct_effects**: estimated by a conditional maximization step in EM
  (Meng & Rubin, 1993); parameter `direct_effect`, term
  `"indicator:covariate"` in `parameter_inference()`. A term that is not
  `"indicator ~ covariate"` or an indicator that is not continuous raises
  `latents_bad_argument`; a covariate absent, non-numeric, missing or constant
  raises `latents_bad_data`.
- **local_dependence**: a pair mixing categorical and continuous, naming an
  indicator twice or sharing one with another pair raises
  `latents_bad_argument`; a categorical item whose categories contain ":"
  raises `latents_bad_data`. `predict()`, `parameter_inference(data = )` and
  the other functions that take `data` join the pair in the caller's data.
- **search**: with a supplied start and `n_starts = 1`, screening is skipped;
  `n_starts = 1` under emEM finalizes Ward's start whenever it screens
  successfully even if another start screened higher. Low effective Gaussian
  membership warns `latents_thin_profiles`. The selected search is recorded on
  the fit.
- **time**: zero-padded numbers and ISO dates are accepted; `latents_bad_time`
  is also of class `latents_bad_data`.
- **id = NULL**: the internal unit column is `.observation`; a `data` with a
  column of that name raises `latents_bad_data`; the fit raises a
  `latents_single_level` message.
- **Value**: `n_informative` is the number of rows with at least one observed
  indicator; `effective_profile_counts` and `effective_group_counts` are
  weighted under `weights`; `small_classes` warns `latents_small_classes`.
