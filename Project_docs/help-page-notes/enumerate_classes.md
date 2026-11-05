# Implementation notes removed from `?enumerate_classes` (2026-11-03)

- **Parallel estimation**: socket workers load the installed package, which
  must match the session's version; each worker needs memory for its tasks;
  tasks are individual starts or whole candidates depending on the family;
  warnings and messages are reported in candidate order; without `seed`,
  starting values differ between runs; lost results retried once with fewer
  workers, `latents_worker_lost` warning; unrecovered tasks recorded as
  failed starts or candidates.
- **Conditions** (formerly a section): `latents_bad_argument` for an unknown
  `model`, a structure given twice, `model` with no continuous indicator, a
  bad `seed`, an argument `multilpa()` does not take, `n_group_classes` other
  than 1 with `id = NULL`, missing `id`, `between_variance` with the profile
  family, `n_profiles` or `model` with a group-class family or `id = NULL`
  there, `time` with a group-class family or `between_variance`, `cores`
  below 1, an invalid or conflicting `type`.
- **Value**: the print header reports the resolved type;
  `get_results(x, "best")`, `"classification"`, `"deltas"`, `"knees"` are the
  `what` values behind the getters.
