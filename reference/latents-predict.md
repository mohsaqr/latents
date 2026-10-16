# Predict from and simulate a fitted model

[`predict()`](https://rdrr.io/r/stats/predict.html) evaluates a fit on
new rows. For
[`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md)
and
[`multilca()`](https://pak.dynasite.org/latents/reference/multilca.md)
fits without membership covariates it is the E-step at the estimates,
applied through the fit's own preparation (categorical levels,
centering, missing-data handling), so predicting the training rows
reproduces the fit's posteriors exactly. In a two-level model the groups
in `newdata` are new groups, each classified from its own rows. For
[`mixture_regression()`](https://pak.dynasite.org/latents/reference/mixture_regression.md)
fits it returns fitted outcomes or class probabilities. Growth mixtures
refuse with `latents_unsupported_prediction`; their fitted curves are
the `"trajectories"` and `"individual"` tables of
[`get_results()`](https://pak.dynasite.org/latents/reference/get_results.md).
Membership- covariate and latent transition fits refuse with the same
class; their fitted rows' posteriors are `get_results(x, "posteriors")`.

## Usage

``` r
# S3 method for class 'latents_growth_mixture'
simulate(object, nsim = 1, seed = NULL, ...)

# S3 method for class 'latents_growth_mixture'
predict(object, ...)

# S3 method for class 'latents_mixture_regression'
predict(
  object,
  newdata = NULL,
  type = c("response", "class_response", "posterior", "probabilities"),
  ...,
  trials = NULL
)

# S3 method for class 'latents_mixture_regression'
simulate(object, nsim = 1, seed = NULL, ...)

# S3 method for class 'multilpa'
predict(object, newdata = NULL, type = c("class", "posterior", "density"), ...)

# S3 method for class 'multilpa_covariates'
predict(object, ...)

# S3 method for class 'multilpa_transitions'
predict(object, ...)
```

## Arguments

- object:

  A fitted model.

- nsim:

  Number of simulated outcome vectors.

- seed:

  `NULL` uses and advances the caller's random stream. An integer makes
  simulation reproducible and restores the caller's random state.

- ...:

  Unused.

- newdata:

  `NULL` for the fitted rows, or a data frame with the indicator or
  predictor columns (and `id` and membership covariates where the model
  uses them).

- type:

  For a profile or class model: `"class"` (modal profile and group class
  with their posteriors), `"posterior"` (every profile's probability as
  well) or `"density"` (each row's log density, marginal over group
  classes, for scoring or held-out likelihood). For a mixture
  regression: `"response"` (the mean outcome over classes; expected
  successes for a binomial count response), `"class_response"` (one row
  per row and class), `"posterior"` (class probabilities given the
  outcome, which `newdata` must contain) or `"probabilities"` (each
  category under each class; ordinal family only).

- trials:

  For binomial mixture-regression predictions, the number of trials: one
  non-negative whole number, a vector with one value per row, or a
  column name in the prediction data. `NULL` uses the fitted trials when
  `newdata = NULL`, or infers them from the response columns when
  present. A count response otherwise requires `trials`; a Bernoulli fit
  defaults to one trial. For `type = "posterior"`, any explicit value
  must agree with the observed response. Predictions remain on the
  expected success-count scale; `trials = 1` requests probabilities for
  each trial.

## Value

A base `data.frame` with one row per row of `newdata` (per row and class
for `"class_response"`, per row, class and category for
`"probabilities"`), starting with `row`, its position.
[`simulate()`](https://rdrr.io/r/stats/simulate.html) returns one row
per fitted row and one column per simulation, `sim_1`, `sim_2`, ...

## Details

[`simulate()`](https://rdrr.io/r/stats/simulate.html) draws new outcomes
on the fitted design of a mixture regression or growth mixture: class
memberships from the fitted mixing model (group classes first in a
two-level model), random effects from each class's covariance, then
outcomes from each row's class regression.

## Conditions

`latents_bad_data` when `newdata` lacks an indicator or the `id` column,
has a category the fit never saw, or a missing value the fit cannot
handle, or a regression design cannot be evaluated. Binomial count
predictions without trial counts also raise `latents_bad_data`;
`latents_bad_argument` for `"probabilities"` outside the ordinal family;
`latents_unsupported_prediction` for
[`predict()`](https://rdrr.io/r/stats/predict.html) on a growth mixture,
a membership-covariate fit or a latent transition fit.

## Examples

``` r
fit <- mixture_regression(score ~ hours, data = study_hours, n_classes = 2,
                          n_starts = 2, seed = 1)
predict(fit, data.frame(hours = c(2, 10)), type = "class_response")
#>   row   class     prior   fitted
#> 1   1 class_1 0.5647883 43.49612
#> 2   2 class_1 0.5647883 80.00444
#> 3   1 class_2 0.4352117 56.52624
#> 4   2 class_2 0.4352117 62.96269
head(simulate(fit, nsim = 2, seed = 1))
#>      sim_1    sim_2
#> 1 85.96028 84.66181
#> 2 43.06464 61.01232
#> 3 57.85893 62.03282
#> 4 51.20502 72.29366
#> 5 78.31283 84.47117
#> 6 59.89099 71.01467
```
