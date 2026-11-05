# Implementation notes removed from `?lta` (2026-11-03)

- **select_start**: "likelihood" selects among starts within
  `tol * max(1, n_observed)` of the highest log likelihood (observed
  indicator values in positive-weight rows), preferring converged starts;
  "converged" falls back to all starts when none converged; ties resolve by
  start order.
- **acceleration**: an extrapolated point that leaves a probability
  non-positive or a row not summing to one is replaced by the plain EM point;
  convergence is assessed over a complete cycle; `iterations` and `max_iter`
  count EM steps under either setting.
- **cores**: forked workers on macOS and Linux, a socket cluster on Windows
  and on macOS with Accelerate BLAS; lost results retried once with half as
  many workers, one `latents_worker_lost` warning; a start still lost is
  recorded as failed; all starts failing raises `latents_all_starts_failed`;
  socket workers use bounded chunks and completed tasks in a failed chunk may
  be retried.
- **transitions = "occasion"**: an unidentified move's logits are reported
  with `estimated = FALSE` and no standard error. A profile never occupied
  before a later occasion raises `latents_empty_transition_row`.
- **transition_covariates with occasions = "grid"**: the skipped-position
  refusal is `latents_bad_data`.
- **time**: `latents_bad_time` is also of class `latents_bad_data`;
  zero-padded numbers and ISO dates are accepted.
- **n_profiles < 2**: `latents_bad_transition`.
- **model**: the homogeneous engine rejects Wald inference for the ten
  structures other than EEI, VVI, EEE, VVV with
  `latents_unsupported_inference`.
- **multilpa_lta tables**: Wald standard errors are observed, robust or OPG
  from analytic scores; a move only zero-weight units make is left out of the
  implied transition table; the tables take `data` and `truth`; any other
  argument raises `latents_bad_argument`; `summary()` takes `level` and
  `vcov_type`.
- **engine**: a diagonal `multilpa_lta` fit without an explicit `model` is
  finished by a quasi-Newton search on the exact likelihood; other fits by
  EM. The general engine warns `latents_unconverged` even at `max_iter = 0`;
  the homogeneous engine does not.
- **first-order occasion-specific transitions**: the move into the first
  observed occasion is also flagged unidentified when earlier occasions are
  empty.
