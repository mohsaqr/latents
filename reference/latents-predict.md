# Predict from and simulate a fitted model

[`predict()`](https://rdrr.io/r/stats/predict.html) applies a fitted
model to the original sample or new observations to obtain class
membership probabilities or predicted outcomes.

For latent profile and class models, it estimates each observation's
probability of belonging to each class, given its observed indicators,
and identifies the most probable class. These probabilities express
uncertainty in classification; the assigned class provides a discrete
summary of that uncertainty.

For mixture regression models, it predicts outcomes from the
relationships estimated within each class, either separately by class or
averaged over class membership probabilities. Growth mixture models also
provide individual trajectories informed by the person's observed
outcomes.

[`simulate()`](https://rdrr.io/r/stats/simulate.html) generates
replicate outcomes under a fitted mixture regression or growth mixture
model, preserving the original predictors and grouping structure. The
simulations reflect the model's assumed variation in class membership
and outcomes while holding the estimated parameters fixed.

## Usage

``` r
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
```

## Arguments

- object:

  A fitted model.

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
  per observation and class), `"posterior"` (class probabilities given
  the outcome, which `newdata` must contain) or `"probabilities"` (each
  category under each class; ordinal family only). For a growth mixture:
  `"response"` (the class trajectories averaged with each person's prior
  class probabilities), `"class_response"` (one row per observation and
  class, with random effects set to zero), and, using outcomes in
  `newdata`, `"posterior"` (one row per person), `"class_conditional"`
  (one row per observation and class: the class trajectory plus the
  person's random effects predicted under that class, E(b \| y, class)),
  `"conditional"` (those averaged with the posterior class
  probabilities) and `"random_effects"` (one row per person and class
  with its posterior probability; their posterior-weighted mean is
  [`lcmm::predictRE()`](https://cecileproust-lima.github.io/lcmm/reference/predictRE.html)).

- ...:

  Unused.

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

- nsim:

  Number of simulated outcome vectors.

- seed:

  `NULL` uses and advances the caller's random stream. An integer makes
  simulation reproducible and restores the caller's `.Random.seed`. R
  does not store the cached Box-Muller normal draw in `.Random.seed`, so
  the next normal draw is not preserved when that generator has a cached
  draw.

## Value

A base `data.frame` with one row per row of `newdata` (per row and class
for `"class_response"` and `"class_conditional"`, per row, class and
category for `"probabilities"`), starting with `row`, its position;
`"posterior"` and `"random_effects"` of a growth mixture have one row
per person (and class), starting with the `id` column.
[`simulate()`](https://rdrr.io/r/stats/simulate.html) returns one row
per fitted row and one column per simulation, `sim_1`, `sim_2`, ...

## Details

For
[`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md)
and
[`multilca()`](https://pak.dynasite.org/latents/reference/multilca.md)
models without membership covariates, it evaluates the E-step at the
estimated parameters using the fitted categorical levels, centering, and
missing-data rules. Applying it to the estimation data reproduces the
fitted posterior probabilities for classified observations. In a
two-level model, groups in `newdata` are treated as new groups and
classified using their supplied observations. With the default
`drop_empty = TRUE`, rows with every indicator missing retain their
position but have `NA` for the predicted profile and its probabilities.
A group is unclassified only when none of its rows has an observed
indicator. Fitted rows with zero sampling weight are unclassified under
either setting of `drop_empty`. New-data prediction does not apply
sampling weights. These exclusions apply to classification output;
density predictions retain each row's evaluated log density. For a model
with a noise component, prediction supports complete rows and, under
`missing = "fiml"`, rows with every indicator missing. Partially
observed rows are unsupported because the fitted noise hypervolume does
not specify their marginal density. An entirely missing row has log
density zero, including under the noise component.

For growth mixtures (`mixture_regression(random = )`), individuals and
clusters in `newdata` are treated as new units. Class mean trajectories
require only predictors; posterior class probabilities, individual
trajectories, and empirical Bayes random effects also require observed
outcomes. Prediction is not supported for membership-covariate or latent
transition fits and raises `latents_unsupported_prediction`. For these
fits, use `get_results(x, "posteriors")` to extract posterior
probabilities for the estimation data.

[`simulate()`](https://rdrr.io/r/stats/simulate.html) draws new outcomes
on the fitted design of a mixture regression or growth mixture: class
memberships from the fitted mixing model (group classes first in a
two-level model), random effects from each class's covariance, then
outcomes from each row's class regression.

## Conditions

`latents_bad_data` when `newdata` lacks an indicator or the `id` column,
contains a category absent from the estimation data, contains
unsupported missing values, or does not allow evaluation of the
regression design. Binomial count predictions without trial counts also
raise `latents_bad_data`; `latents_bad_argument` for `"probabilities"`
outside the ordinal family; `latents_unsupported_prediction` for
[`predict()`](https://rdrr.io/r/stats/predict.html) on a
membership-covariate fit or a latent transition fit.
`latents_unsupported_noise` for partially observed prediction rows in a
model with a noise component. A growth mixture prediction raises
`latents_bad_data` when `newdata` lacks `id`, the cluster column or (for
the outcome-based types) the outcome, and warns with
`latents_unused_argument` for arguments it ignores.

## Examples

``` r
fit <- mixture_regression(score ~ hours, data = study_hours, n_classes = 2,
                          n_starts = 2, seed = 1)
predict(fit, data.frame(hours = c(2, 10)), type = "class_response")
#>   row   class     prior   fitted
#> 1   1 class_1 0.5645888 43.49393
#> 2   2 class_1 0.5645888 80.00640
#> 3   1 class_2 0.4354112 56.52305
#> 4   2 class_2 0.4354112 62.96629
head(simulate(fit, nsim = 2, seed = 1))
#>      sim_1    sim_2
#> 1 85.96107 84.66292
#> 2 43.06356 61.01025
#> 3 57.86180 62.03635
#> 4 51.20585 72.29776
#> 5 78.31739 84.47418
#> 6 59.89100 71.01171

# Predict an individual trajectory conditional on observed scores.
growth <- mixture_regression(score ~ wave, data = subset(growth_scores, student <= 60),
                             n_classes = 2, id = "student", class_level = "group",
                             random = ~ 1, n_starts = 1, seed = 1)
predict(growth, subset(growth_scores, student == 61), type = "conditional")
#>   row   fitted
#> 1   1 58.88460
#> 2   2 65.58452
#> 3   3 72.28445
#> 4   4 75.63441
```
