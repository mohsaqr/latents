# Implementation notes removed from `?get_results` (2026-11-03)

- **Alignment**: before combining supplied `data` with results, row alignment
  is checked on identifiers, occasions and indicators the fit retains;
  values serialised by `utils::write.csv()` are accepted when distinct values
  survive; columns absent from the fit cannot establish alignment.
- **"responses"**: a probability on its bound has an `NA` standard error;
  zero rows with the same columns when no indicator is categorical
  (likewise "profiles" when every indicator is categorical).
- **"posteriors"**: in wide format excluded rows have `NA` profile and
  posterior columns.
- **"classification" weighting**: `estimated_n`/`estimated_proportion` are
  weighted posterior totals, `average_posterior` a weighted mean;
  `n_modal`/`proportion_modal` stay counts.
- **"classification_errors"**: sum_i w_i p_it 1(W_i = s) / sum_i w_i p_it;
  noise fits raise `latents_unsupported_noise` and "all"/summary omit it.
- **"bch_weights"**: `level = "both"` should not feed a weighted regression.
- **"entropy" weighting**: `entropy_sum = sum_i w_i e_i`,
  `relative_entropy = 1 - sum_i w_i e_i / (sum_i w_i log K)`.
- **"sequences"**: rows with nothing observed are left out (kept under
  `drop_empty = FALSE`); `latents_no_time` without `time`.
- **"gaussian_diagnostics"**: `positive_rows` is the same for every
  component; `min_variance_ratio` uses the sampling-weighted marginal
  variance on the fitted scale including between-component variation;
  `min_standardized_eigenvalue` standardises by the weighted marginal sd;
  missing reference variation gives `NA`.
- **"candidates"**: a class with no assigned unit is excluded from the AvePP
  min/max/mean; a level with one class has `NA` summaries; `what =
  "criteria"` is a deprecated name for "best".
- **"knees"**: the knee is measured perpendicular to the chord; a line needs
  four candidates for `knee_without_largest` and `stable`, else `NA` with
  `note`.
- **"local_dependence"**: unobserved combinations sit on the
  `min_probability` floor and count as parameters; `latents_no_local_dependence`
  / `latents_no_direct_effects` without them.
- **"lmr"**: p-values by Imhof's method to about 1e-9; warnings raised while
  scoring a candidate are recorded in `note`; `latents_no_lmr` with no
  eligible pair; Vermunt (2024) cited for the adjusted test.
- **"tree"**: every pair is listed including zero mass; units are the rows
  and groups the fit counts under `drop_empty`, weighted at both levels;
  `latents_no_tree`, `latents_unsupported_noise`; the summary of a grid
  carries the default table only.
- **Other families**: `confint()` of a regression uses the stored covariance
  and refuses other arguments (`latents_bad_argument`); `vcov()` warns
  `latents_unused_argument` for `vcov_type =` in place of `type =`;
  growth mixtures have no "groups", "fitted" or "classification" table;
  additive `effective_groups = w^2 (1 - w) / Var(w)` and the
  `latents_weak_class` warning below 50; transition/family grid "best"
  reports the lowest value and does not select; with one group class every
  family is the same model.
- **"all"**: only the conditions meaning the table does not exist are
  skipped; any other error is raised; "all" with extra arguments raises
  `latents_bad_argument`.
- **Information criteria**: `deviance`, `aic`, `kic` carry
  `convention = NA`; `clc` is the same in both rows; Cavanaugh's general
  form -2L + 3q is 3 smaller for every model; a two-level fit saved by
  0.97.1 or earlier lacks the rows' within-class entropies (`awe`, `icl`,
  `clc` NA); `n_informative_groups`; transition "groups" criteria use
  group-class entropy, "individuals" occasion-level profile entropy; under
  weights `L` is the pseudo log likelihood and `EN` weighted.
- **Classification quality**: `average_posterior` is `NA` for a class with
  no modal members; `n_modal` and `estimated_n` need not coincide;
  `level = "groups"` on a fit without group classes raises
  `latents_no_group_classes`.
- **Residuals**: pooled over profiles under `by = "overall"`; a nominal
  indicator's category order matters for the mixed/count correlation; a
  half-observed joined pair enters only where both items are observed;
  design-based inference not implemented for weighted fits.
- **Recovery**: several truth columns are stacked; group-level recovery
  counts each group once; unused factor levels are omitted.
