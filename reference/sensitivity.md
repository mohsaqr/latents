# Assess sensitivity of a fitted mixture model to initialisation

`sensitivity()` evaluates how consistently a fitted mixture model is
recovered from different random initialisations. Mixture-model
estimation can converge to different local maxima, yielding different
class assignments for the same data and model specification. Examining
this dependence on initialisation helps assess the numerical reliability
of the solution being interpreted.

The function repeats estimation under several random seeds and compares
the resulting likelihood values and class assignments with the reference
fit. Class labels are aligned before assignments are compared. By
considering likelihood and classification agreement together, the
analysis shows whether repeated searches recover a comparable fit and
whether they assign the same observations to the same classes.

Results include convergence status, log likelihood, and classification
agreement for each refit. Similar likelihood values are grouped using a
specified tolerance. The data and model specification remain fixed
throughout the analysis.

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

  The reference fit: a
  [`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md)
  model, an [`lta()`](https://pak.dynasite.org/latents/reference/lta.md)
  model, or a growth mixture
  ([`mixture_regression()`](https://pak.dynasite.org/latents/reference/mixture_regression.md)
  with `random`). A
  [`fit_staged()`](https://pak.dynasite.org/latents/reference/fit_staged.md)
  result varies the second-stage starts while holding its first-stage
  measurement solution fixed. A directly fitted model with measurement
  blocks held fixed is refused (`latents_unsupported_sensitivity`),
  because its original free starting values are not retained. A
  membership-covariate fit raises `latents_unsupported_sensitivity`. For
  an [`lta()`](https://pak.dynasite.org/latents/reference/lta.md) fit,
  `agreement` is the share of occasions given the same profile, after
  matching the refit's profiles on the measurement; for a growth mixture
  it is the share of persons given the same class, after matching the
  refit's classes by their posterior memberships of the fitted persons.

- data:

  Optional for a
  [`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md)
  fit or a growth mixture, which carry the columns they were built from.
  Required for an
  [`lta()`](https://pak.dynasite.org/latents/reference/lta.md) fit
  (`latents_bad_argument` otherwise). It must be the data the model was
  fitted to, in its original row order; shared identifier and indicator
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
  `tolerance`. Defaults to the reference fit's replication tolerance,
  1e-6 per observed value.

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
per seed, with columns `seed`, `log_likelihood`, `converged`,
`iterations`, `optimum` (1 for the group containing the highest attained
log likelihood, 2 for the next lower group, and so on), `best` (whether
this seed reached optimum 1) and `agreement` (the proportion of
classified observations given the same profile as the reference fit,
after label matching). For profile and transition models, this uses
observations included in both fits under their `drop_empty` settings. A
failed refit has `NA` in every column except `seed`. The `optimum` and
`best` columns group attained likelihood values without filtering on
convergence; inspect `converged` before interpreting them as evidence of
replicated maxima.

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

## References

Hipp, J. R., & Bauer, D. J. (2006). Local solutions in the estimation of
growth mixture models. *Psychological Methods, 11*, 36–53.

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

# Growth mixtures are known for local solutions (Hipp and Bauer 2006).
growth <- mixture_regression(score ~ wave, data = subset(growth_scores, student <= 60),
                             n_classes = 2, id = "student", class_level = "group",
                             random = ~ 1, n_starts = 1, seed = 1)
sensitivity(growth, seeds = 1:2)
#>   seed log_likelihood converged iterations optimum best agreement
#> 1    1      -1175.788      TRUE         81       1 TRUE         1
#> 2    2      -1175.788      TRUE         81       1 TRUE         1
```
