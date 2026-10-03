# Simulation from a growth mixture model, and the verbs that refuse growth
# fits clearly where they are not implemented.

#' One simulated outcome vector from a growth mixture fit
#'
#' Each person's class is drawn from their membership probabilities (in a
#' multilevel growth mixture, given their cluster's group class, drawn first),
#' their
#' random effects from that class's covariance, and their outcomes around
#' the class trajectory plus random effects with the class's residual
#' variance, on the fitted design (persons, times and covariates kept).
#' @param object A `latents_growth_mixture` fit.
#' @return A numeric vector, one value per fitted row, in the fitted order.
#' @noRd
.growth_draw <- function(object) {
  spec <- object$spec
  params <- object$params
  structure <- .growth_structure(spec)
  prior <- if (identical(structure$nesting, "observation")) {
    exp(.mixture_log_softmax(spec$w, params$gamma))
  } else {
    # Each cluster's group class first, then each person's class given it.
    group_prior <- exp(.mixture_log_softmax(structure$v, params$delta))
    group_class <- vapply(seq_len(structure$n_groups), function(j) {
      sample.int(structure$n_group_classes, 1L, prob = group_prior[j, ])
    }, integer(1))
    by_group_class <- lapply(seq_len(structure$n_group_classes), function(h) {
      exp(.mixture_log_softmax(.mixture_two_level_design(structure, h),
                               params$class_logits))
    })
    person_group_class <- group_class[structure$group_index]
    Reduce(`+`, lapply(seq_len(structure$n_group_classes), function(h) {
      by_group_class[[h]] * (person_group_class == h)
    }))
  }
  n_persons <- nrow(prior)
  classes <- vapply(seq_len(n_persons), function(i) {
    sample.int(spec$n_classes, 1L, prob = prior[i, ])
  }, integer(1))
  covariances <- .growth_covariances(spec, params)
  n_random <- ncol(spec$random_design)
  effects <- t(vapply(seq_len(n_persons), function(i) {
    as.vector(.growth_factor(covariances[[classes[i]]]) %*% stats::rnorm(n_random))
  }, numeric(n_random))) |> matrix(n_persons, n_random)
  design <- cbind(spec$x, spec$z)
  person <- spec$group_index
  row_class <- classes[person]
  mean <- vapply(seq_len(spec$n), function(r) {
    sum(design[r, ] * .growth_coefficients(params, row_class[r]))
  }, numeric(1)) + spec$offset + rowSums(spec$random_design * effects[person, , drop = FALSE])
  stats::rnorm(spec$n, mean, sqrt(params$sigma2[row_class]))
}

#' Simulate outcomes from a growth mixture model
#'
#' Draws new outcomes on the fitted design: each person's class from their
#' membership probabilities, their random effects from that class's
#' covariance, and their outcomes around their own trajectory with the
#' class's residual variance.
#'
#' @param object A `latents_growth_mixture` fit.
#' @param nsim Number of simulated outcome vectors.
#' @param seed `NULL` or an integer seed; the caller's random state is
#'   restored afterwards.
#' @param ... Unused.
#' @return A base `data.frame` with one row per fitted row and one column per
#'   simulation, `sim_1`, `sim_2`, ...
#' @examples
#' \donttest{
#' fit <- mixture_regression(score ~ wave, growth_scores, n_classes = 3,
#'                           id = "student", class_level = "group",
#'                           random = "wave", random_covariance = "equal",
#'                           n_starts = 2, seed = 1)
#' head(simulate(fit, nsim = 2, seed = 1))
#' }
#' @export
simulate.latents_growth_mixture <- function(object, nsim = 1, seed = NULL, ...) {
  stopifnot("`nsim` must be a single positive integer" = .mixture_is_count(nsim))
  .mixture_with_seed(seed, {
    draws <- lapply(seq_len(nsim), function(i) .growth_draw(object))
    stats::setNames(as.data.frame(draws), paste0("sim_", seq_len(nsim)))
  })
}

#' Prediction for a growth mixture model (not implemented)
#'
#' Predictions for new persons are not implemented yet. The fitted class
#' trajectories and each fitted person's own predicted curve are tables of
#' [get_results.latents_growth_mixture()]: `"trajectories"` and
#' `"individual"`.
#'
#' @param object A `latents_growth_mixture` fit.
#' @param ... Unused.
#' @return Never returns; raises `latents_unsupported_prediction`.
#' @export
predict.latents_growth_mixture <- function(object, ...) {
  stop(errorCondition(paste(
    "predict() is not implemented for growth mixture models yet. The class",
    "trajectories are get_results(x, \"trajectories\"), and each fitted",
    "person's class trajectory and own predicted curve are",
    "get_results(x, \"individual\")."),
    class = "latents_unsupported_prediction", call = NULL))
}
