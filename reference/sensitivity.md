# Refit a model from several random seeds and check whether it reaches the same solution

`sensitivity()` refits a fitted mixture model from several random seeds
and reports whether the refits reach the same solution. Mixture-model
estimation can converge to different local maxima, yielding different
class assignments for the same data and model specification, so a
solution that only one seed finds should not be interpreted.

Each seed gives one refit of the same model to the same data. The
function records the refit's log likelihood, whether it converged and
how many iterations it took, and compares its class assignments with
those of the reference fit, after aligning the class labels. Log
likelihoods within `tolerance` of each other are grouped as one optimum,
so the table shows how many seeds reached the best one and whether those
that did assign the same observations to the same classes.

## Usage

``` r
sensitivity(
  x,
  data = NULL,
  seeds = 1:10,
  n_starts = NULL,
  max_iter = NULL,
  tol = NULL,
  tolerance = NULL,
  cores = 1L
)
```

## Arguments

- x:

  The reference fit. Accepted:

  - a profile or class fit from
    [`lpa()`](https://pak.dynasite.org/latents/reference/lpa.md),
    [`lca()`](https://pak.dynasite.org/latents/reference/lca.md),
    [`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md)
    or
    [`multilca()`](https://pak.dynasite.org/latents/reference/multilca.md)
    without membership covariates. A
    [`fit_staged()`](https://pak.dynasite.org/latents/reference/fit_staged.md)
    result varies the second-stage starts while holding its first-stage
    measurement solution fixed;

  - a transition fit from
    [`lta()`](https://pak.dynasite.org/latents/reference/lta.md), which
    also needs `data`;

  - a growth mixture
    ([`mixture_regression()`](https://pak.dynasite.org/latents/reference/mixture_regression.md)
    with `random`).

  Refused with `latents_unsupported_sensitivity`: a membership-covariate
  fit, and a directly fitted model with measurement blocks held fixed
  (`fixed`), because its original free starting values are not retained.
  Refused with `latents_not_a_fit`: the other
  [`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md)
  families (additive, dispersion and cross-level) and a mixture
  regression without `random`.

- data:

  The data the model was fitted to, in its original row order. Optional
  for a profile or class fit and for a growth mixture, which carry the
  columns they were built from; required for an
  [`lta()`](https://pak.dynasite.org/latents/reference/lta.md) fit
  (`latents_bad_argument` otherwise). Shared identifier and indicator
  columns are checked against the fit (`latents_bad_inference_data`).

- seeds:

  At least two distinct integer seeds accepted by
  [`set.seed()`](https://rdrr.io/r/base/Random.html). Default `1:10`.
  Invalid values raise `latents_bad_argument`. Reusing the reference
  seed reproduces its initialisation when the other fitting controls are
  unchanged.

- n_starts:

  Starts per refit. Defaults to the number the reference fit used.

- max_iter, tol:

  Passed to each refit. Default to the reference fit's (for a fit that
  does not record them, to the defaults of the function that fitted it).

- tolerance:

  Maximum log-likelihood difference for grouping refits. Each group is
  anchored at its highest log likelihood and spans no more than
  `tolerance`. Only converged refits receive an optimum index or a best
  flag. Defaults to the reference fit's replication tolerance, 1e-6 per
  observed value.

- cores:

  Number of processor cores (default 1). The seeds are refitted in
  parallel, one seed per worker at a time, each refit on one core. Every
  refit sets its own seed, so the result is identical to `cores = 1`.
  Workers are forked processes on macOS and Linux, and a socket cluster
  on Windows and on macOS with Apple's Accelerate BLAS, as in
  [`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md);
  socket workers load the installed copy of the package. More cores than
  seeds are not used.

## Value

A `data.frame` of class `latents_sensitivity` (which
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) draws), one row
per seed, in the order of `seeds`, with the columns

- `seed`:

  integer: the seed of the refit.

- `log_likelihood`:

  numeric: the refit's maximised log likelihood.

- `converged`:

  logical: whether the refit met its convergence criterion.

- `iterations`:

  integer: the number of iterations the refit reports.

- `optimum`:

  integer: 1 for the group of converged refits holding the highest log
  likelihood, 2 for the next lower group, and so on; `NA` for a refit
  that did not converge.

- `best`:

  logical: whether the refit converged to optimum 1.

- `agreement`:

  numeric, 0 to 1: the proportion of units given the same class as in
  the reference fit, after label matching. The units are the
  observations of a profile or class fit, the occasions of a transition
  fit (profiles matched on the measurement) and the persons of a growth
  mixture (classes matched by their posterior memberships of the fitted
  persons). For profile and transition models, only observations
  included in both fits under their `drop_empty` settings are compared.

A failed refit has `NA` in every column except `seed`, and the call
warns with `latents_sensitivity_dropped`. The reference fit itself is
not a row of the table. Invalid arguments raise `latents_bad_argument`;
see `x` and `data` for the other refusals.

[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) returns a
plain `data.frame`.
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) accepts `main`
and `subtitle` and returns a ggplot object showing each refit's
log-likelihood difference from the best converged refit and its
classification agreement with the reference fit. Unfinished refits are
distinguished from converged maxima; if none converged, the highest
observed log likelihood is the reference instead.

## Details

Replication across seeds complements replication across starts within
one fit. Each refit uses the reference fit's estimation controls unless
alternatives are supplied.

A seed whose refit fails gives a row with `NA` estimates, and a
`latents_sensitivity_dropped` warning reports how many seeds failed.

## Interpretation

Consistent results across converged refits support numerical
reproducibility, but do not guarantee a global maximum or establish the
validity of the model. Similar likelihood values alone do not imply
identical parameter estimates. The analysis assesses sensitivity to
initialisation, rather than uncertainty arising from sampling or model
selection.

## See also

[`enumerate_classes()`](https://pak.dynasite.org/latents/reference/enumerate_classes.md)
to vary the number of classes, and
[`parameter_inference()`](https://pak.dynasite.org/latents/reference/parameter_inference.md)
for uncertainty about the parameters of a single solution.

## Examples

``` r
set.seed(7)
example_data <- data.frame(
  school = rep(seq_len(12), each = 10),
  score_a = c(rnorm(60, -1.5), rnorm(60, 1.5)),
  score_b = c(rnorm(60, -1), rnorm(60, 1))
)
fit <- multilpa(example_data, c("score_a", "score_b"), "school",
                n_profiles = 2, n_group_classes = 1, n_starts = 3, seed = 1)
sensitivity(fit, seeds = 1:3, n_starts = 3)
#>   seed log_likelihood converged iterations optimum best agreement
#> 1    1      -396.1633      TRUE          9       1 TRUE         1
#> 2    2      -396.1633      TRUE          9       1 TRUE         1
#> 3    3      -396.1633      TRUE          9       1 TRUE         1

# Growth mixtures are known for local solutions.
growth <- mixture_regression(score ~ wave, data = subset(growth_scores, student <= 60),
                             n_classes = 2, id = "student", class_level = "group",
                             random = ~ 1, n_starts = 1, seed = 1)
sensitivity(growth, seeds = 1:2)
#>   seed log_likelihood converged iterations optimum best agreement
#> 1    1      -1175.788      TRUE         81       1 TRUE         1
#> 2    2      -1175.788      TRUE         81       1 TRUE         1
```
