# Implementation notes removed from `?fit_staged` (2026-11-03)

- **measurement**: an incompatible model raises `latents_bad_stage`; a
  non-`multilpa` object, `n_group_classes` below two, or an unknown
  `variance_model`, `covariance_model` or `missing` raises
  `latents_bad_argument`.
- **n_screen**: a value with `search = "kmeans"` or below `n_starts` raises
  `latents_bad_argument`.
