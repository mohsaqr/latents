# Simulation from a growth mixture model.

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
  if (!.growth_is_correlated(spec)) {
    return(stats::rnorm(spec$n, mean, sqrt(params$sigma2[row_class])))
  }
  # Correlated residuals: each person's errors jointly, process plus noise.
  residual <- numeric(spec$n)
  rows <- split(seq_len(spec$n), person)
  invisible(lapply(seq_along(rows), function(i) {
    r <- rows[[i]]
    covariance <- .growth_process_covariance(spec, params$correlation, spec$residual_times[r]) +
      diag(params$sigma2[classes[i]], length(r))
    residual[r] <<- as.vector(crossprod(chol(covariance), stats::rnorm(length(r))))
  }))
  mean + residual
}

#' @rdname latents-predict
#' @export
simulate.latents_growth_mixture <- function(object, nsim = 1, seed = NULL, ...) {
  .latents_require("`nsim` must be a single positive integer" = .mixture_is_count(nsim))
  .mixture_with_seed(seed, {
    draws <- lapply(seq_len(nsim), function(i) .growth_draw(object))
    stats::setNames(as.data.frame(draws), paste0("sim_", seq_len(nsim)))
  })
}
