# Refit a model under several seeds

Refits the model under each of several random seeds and returns one row
per seed with the maximised log likelihood, convergence, which distinct
optimum the refit reached, and the proportion of observations assigned
to the same profile as in the reference fit. Profile labels of each
refit are matched to the reference fit's before agreement is computed.

## Usage

``` r
sensitivity(
  x,
  data = NULL,
  seeds = 1:10,
  n_starts = NULL,
  max_iter = NULL,
  tol = NULL,
  tolerance = NULL
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
  membership-covariate fit is refused with
  `latents_unsupported_sensitivity`. For an
  [`lta()`](https://pak.dynasite.org/latents/reference/lta.md) fit,
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

  Seeds to refit under: at least two distinct whole numbers that
  [`set.seed()`](https://rdrr.io/r/base/Random.html) accepts. Default
  `1:10`. Anything else raises `latents_bad_argument`. If the reference
  fit's own seed is included, that row reproduces the reference fit.

- n_starts:

  Starts per refit. Defaults to the number the reference fit used.

- max_iter, tol:

  Passed to each refit. Default to the reference fit's (for an
  [`lta()`](https://pak.dynasite.org/latents/reference/lta.md) fit, to
  the defaults of
  [`lta()`](https://pak.dynasite.org/latents/reference/lta.md)).

- tolerance:

  Two log likelihoods within this of each other count as the same
  optimum. Defaults to the reference fit's replication tolerance (for a
  growth mixture, `1e-6 * (1 + abs(logLik))`).

## Value

A base `data.frame`, one row per seed, with columns `seed`,
`log_likelihood`, `converged`, `iterations`, `optimum` (1 for the best
maximum found across the seeds, 2 for the next distinct one, and so on),
`best` (whether this seed reached optimum 1) and `agreement` (the
proportion of observations given the same profile as the reference fit,
after label matching). A failed refit has `NA` in every column except
`seed`.

## Details

A fit's own replication count describes only the starts within one
seed's random-number stream; refitting under other seeds shows whether
other streams reach the same maximum and the same classification.

A seed whose refit fails gives a row with `NA` estimates, and a
`latents_sensitivity_dropped` warning reports how many seeds failed.

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
#> 1    1      -1175.788      TRUE         70       1 TRUE         1
#> 2    2      -1175.788      TRUE         70       1 TRUE         1
```
