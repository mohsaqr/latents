#' Predict from and simulate a fitted model
#'
#' `predict()` evaluates a fit on new rows. For [multilpa()] and [multilca()]
#' fits without membership covariates it is the E-step at the estimates,
#' applied through the fit's own preparation (categorical levels, centering,
#' missing-data handling), so predicting the training rows reproduces the
#' fit's posteriors exactly. In a two-level model the groups in `newdata` are
#' new groups, each classified from its own rows. For [mixture_regression()]
#' fits it returns fitted outcomes or class probabilities. Growth mixtures
#' refuse with `latents_unsupported_prediction`; their fitted curves are the
#' `"trajectories"` and `"individual"` tables of [get_results()]. Membership-
#' covariate and latent transition fits refuse with the same class; their
#' fitted rows' posteriors are `get_results(x, "posteriors")`.
#'
#' `simulate()` draws new outcomes on the fitted design of a mixture
#' regression or growth mixture: class memberships from the fitted mixing
#' model (group classes first in a two-level model), random effects from
#' each class's covariance, then outcomes from each row's class regression.
#'
#' @param object A fitted model.
#' @param newdata `NULL` for the fitted rows, or a data frame with the
#'   indicator or predictor columns (and `id` and membership covariates where
#'   the model uses them).
#' @param type For a profile or class model: `"class"` (modal profile and
#'   group class with their posteriors), `"posterior"` (every profile's
#'   probability as well) or `"density"` (each row's log density, marginal
#'   over group classes, for scoring or held-out likelihood). For a mixture
#'   regression: `"response"` (the mean outcome over classes; expected
#'   successes for a binomial count response),
#'   `"class_response"` (one row per row and class), `"posterior"` (class
#'   probabilities given the outcome, which `newdata` must contain) or
#'   `"probabilities"` (each category under each class; ordinal family only).
#' @param trials For binomial mixture-regression predictions, the number of
#'   trials: one non-negative whole number, a vector with one value per row,
#'   or a column name in the prediction data. `NULL` uses the fitted trials
#'   when `newdata = NULL`, or infers them from the response columns when
#'   present. A count response otherwise requires `trials`; a Bernoulli fit
#'   defaults to one trial. For `type = "posterior"`, any explicit value must
#'   agree with the observed response. Predictions remain on the expected
#'   success-count scale; `trials = 1` requests probabilities for each trial.
#' @param nsim Number of simulated outcome vectors.
#' @param seed `NULL` uses and advances the caller's random stream. An integer
#'   makes simulation reproducible and restores the caller's random state.
#' @param ... Unused.
#' @return A base `data.frame` with one row per row of `newdata` (per row and
#'   class for `"class_response"`, per row, class and category for
#'   `"probabilities"`), starting with `row`, its position. `simulate()`
#'   returns one row per fitted row and one column per simulation, `sim_1`,
#'   `sim_2`, ...
#' @section Conditions: `latents_bad_data` when `newdata` lacks an indicator
#'   or the `id` column, has a category the fit never saw, or a missing value
#'   the fit cannot handle, or a regression design cannot be evaluated.
#'   Binomial count predictions without trial counts also raise
#'   `latents_bad_data`; `latents_bad_argument` for `"probabilities"`
#'   outside the ordinal family; `latents_unsupported_prediction` for
#'   `predict()` on a growth mixture, a membership-covariate fit or a latent
#'   transition fit.
#' @examples
#' fit <- mixture_regression(score ~ hours, data = study_hours, n_classes = 2,
#'                           n_starts = 2, seed = 1)
#' predict(fit, data.frame(hours = c(2, 10)), type = "class_response")
#' head(simulate(fit, nsim = 2, seed = 1))
#' @name latents-predict
NULL
