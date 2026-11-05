# Implementation notes removed from `?bootstrap_lrt` (2026-11-03)

- **iter**: transition comparisons accept any positive integer; the others
  need at least two.
- **data**: for group-class families the stored sufficient statistics are
  used when `data` is omitted; for covariate models the data are rebuilt
  from the alternative.
- **Conditions** (formerly a section): `latents_bad_nesting` (class counts,
  fixed blocks, covariate specifications, parameter counts);
  `latents_incomparable_models` (data or fitting specification);
  `latents_bad_argument` (unsupported nesting for transition and family
  comparisons; transition models without `data`; bad seed; incompatible
  families); `latents_unsupported_bootstrap` (person centring or scaling);
  `latents_unsupported_prior`, `latents_unsupported_weights`,
  `latents_unsupported_noise`; `latents_bad_data` (unsupported missing
  indicators); `latents_bad_inference_data`; `latents_no_converge` (also
  active variance bounds for profile models);
  `latents_reversed_likelihood`; `latents_failed_replicates` (warning).
- **Monte Carlo standard error**: summarises uncertainty from the finite
  number of replicates.
