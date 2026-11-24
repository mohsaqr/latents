# Predict class membership or outcomes from a fitted model, and simulate data from it

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

  A fitted model with a prediction method:

  - a profile or class model from
    [`lpa()`](https://pak.dynasite.org/latents/reference/lpa.md),
    [`lca()`](https://pak.dynasite.org/latents/reference/lca.md),
    [`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md)
    or
    [`multilca()`](https://pak.dynasite.org/latents/reference/multilca.md)
    without membership covariates;

  - a mixture regression from
    [`mixture_regression()`](https://pak.dynasite.org/latents/reference/mixture_regression.md);

  - a growth mixture (`mixture_regression(random = )`).

  A membership-covariate fit and an
  [`lta()`](https://pak.dynasite.org/latents/reference/lta.md) fit raise
  `latents_unsupported_prediction`; the additive, dispersion and
  cross-level families of
  [`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md)
  have no [`predict()`](https://rdrr.io/r/stats/predict.html) method.
  [`simulate()`](https://rdrr.io/r/stats/simulate.html) takes a mixture
  regression or a growth mixture only.

- newdata:

  `NULL` (default) for the fitted rows, or a data frame with the
  indicator or predictor columns (and `id` and membership covariates
  where the model uses them).

- type:

  What to predict. The first value listed for each kind of model is its
  default.

  For a profile or class model:

  - `"class"` (default): each row's most probable profile and its
    posterior probability, and in a two-level model the group's most
    probable group class and its posterior probability;

  - `"posterior"`: the same, plus the posterior probability of every
    profile;

  - `"density"`: each row's log density under the fitted model, marginal
    over profiles and group classes, for scoring or a held-out
    likelihood.

  For a mixture regression:

  - `"response"` (default): the mean outcome averaged over the classes
    with their prior probabilities (expected successes for a binomial
    count response);

  - `"class_response"`: the mean outcome under each class, one row per
    observation and class;

  - `"posterior"`: each row's class probabilities given its outcome,
    which `newdata` must contain, and the most probable class;

  - `"probabilities"`: the probability of each outcome category under
    each class; ordinal family only.

  For a growth mixture:

  - `"response"` (default): the class mean trajectories averaged with
    each person's prior class probabilities;

  - `"class_response"`: each class's mean trajectory, with random
    effects set to zero, one row per observation and class;

  - `"posterior"`: each person's class probabilities given the observed
    outcomes, one row per person;

  - `"class_conditional"`: one row per observation and class: the class
    trajectory plus the person's random effects predicted under that
    class, E(b \| y, class);

  - `"conditional"`: the `"class_conditional"` trajectories averaged
    with the posterior class probabilities, one row per observation;

  - `"random_effects"`: the predicted random effects, one row per person
    and class with its posterior probability; their posterior-weighted
    mean is `lcmm::predictRE()`.

  The last four use the outcomes, which `newdata` must then contain.

- ...:

  Unused. A profile, class or growth mixture
  [`predict()`](https://rdrr.io/r/stats/predict.html) warns with
  `latents_unused_argument` for an argument it does not take.

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

  Number of simulated outcome vectors (default 1).

- seed:

  `NULL` (default) uses and advances the caller's random stream. An
  integer makes the simulation reproducible and restores the caller's
  `.Random.seed`. R does not store the cached Box-Muller normal draw in
  `.Random.seed`, so the next normal draw is not preserved when that
  generator has a cached draw.

## Value

A base `data.frame`. Its rows and columns depend on the model and
`type`; `row` is the row's position in `newdata` (or in the fitted
data).

Profile or class model, one row per row of the data:

- `"class"`: `row`, the group identifier (two-level models), `profile`
  (the most probable profile; 0 for the noise component), `posterior`
  (its probability) and, with more than one group class, `group_class`
  and `group_posterior`;

- `"posterior"`: the same, plus `probability_profile_1`,
  `probability_profile_2` and so on;

- `"density"`: `row` and `log_density`.

Mixture regression:

- `"response"`: `row` and `fitted`, one row per observation;

- `"class_response"`: `row`, `class`, `prior` (the class's prior
  probability for the row) and `fitted`, one row per observation and
  class;

- `"posterior"`: `row`, the `id` column (when classes are assigned to
  groups), `class` (the most probable class), `posterior` (its
  probability) and `probability_class_1`, `probability_class_2` and so
  on, one row per observation;

- `"probabilities"`: `row`, `class`, `category` and `probability`, one
  row per observation, class and category.

Growth mixture:

- `"response"` and `"conditional"`: `row` and `fitted`, one row per
  observation;

- `"class_response"`: `row`, `class`, `prior` and `fitted`, one row per
  observation and class;

- `"class_conditional"`: `row`, `class`, `posterior` (the person's
  posterior probability of the class) and `fitted`, one row per
  observation and class;

- `"posterior"`: the `id` column, `class` (the most probable class),
  `probability` (its probability) and `posterior_class_1`,
  `posterior_class_2` and so on, one row per person;

- `"random_effects"`: the `id` column, `class`, `posterior` and one
  column per random effect, one row per person and class.

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
indicator. Fitted rows with zero sampling weight keep their predicted
profile and probabilities under either setting of `drop_empty`, as in
`get_results(x, "assignments")`. New-data prediction does not apply
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
model with a noise component. A profile or class model's
[`predict()`](https://rdrr.io/r/stats/predict.html) raises
`latents_bad_argument` for a `type` it does not offer and warns with
`latents_unused_argument` for arguments it ignores, such as `data =` in
place of `newdata =`. A growth mixture prediction raises
`latents_bad_data` when `newdata` lacks `id`, the cluster column or (for
the outcome-based types) the outcome, and warns with
`latents_unused_argument` for arguments it ignores.

## Examples

``` r
# Classify new observations with a fitted profile model.
profiles <- lpa(faithful, c("eruptions", "waiting"), n_profiles = 2,
                n_starts = 2, seed = 1)
predict(profiles, data.frame(eruptions = c(2, 4.5), waiting = c(55, 80)),
        type = "posterior")
#>   row profile posterior probability_profile_1 probability_profile_2
#> 1   1       2         1          3.110134e-11          1.000000e+00
#> 2   2       1         1          1.000000e+00          1.264747e-23

# The mean outcome under each class of a mixture regression.
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
