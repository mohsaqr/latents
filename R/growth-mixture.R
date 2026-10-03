# Fitting the growth mixture model: argument checks, starting values from the
# latent class growth fit, start screening and selection, class ordering and
# the parameter count. The engine is in growth-mixture-engine.R.

#' Refuse growth mixture requests the model does not define
#' @return `NULL`, invisibly; raises `latents_bad_argument` otherwise.
#' @noRd
.growth_check_arguments <- function(random, family, id, class_level,
                                    n_group_classes, random_diagonal,
                                    cluster = NULL) {
  bad_argument <- function(message) {
    stop(errorCondition(message, class = "latents_bad_argument", call = NULL))
  }
  if (is.null(random)) {
    bad_argument(paste(
      "`cluster` nests persons in clusters for the multilevel growth mixture",
      "model, which needs random effects within classes: give `random`, for",
      "example `random = \"intercept\"`."))
  }
  if (!inherits(random, "formula") || length(random) != 2L) {
    bad_argument("`random` must be a one-sided formula such as `~ 1 + time`.")
  }
  if (!identical(family, "gaussian")) {
    bad_argument(paste(
      "Random effects within classes are implemented for the gaussian family;",
      "binomial and Poisson growth mixtures need numerical integration and",
      "are not supported yet."))
  }
  if (is.null(id) || !identical(class_level, "group")) {
    bad_argument(paste(
      "A growth mixture model classifies persons: pass `id` naming the",
      "persons and `class_level = \"group\"`."))
  }
  if (n_group_classes > 1L && is.null(cluster)) {
    bad_argument(paste(
      "Group classes of a growth mixture are classes of clusters of persons:",
      "pass `cluster` naming the clusters (schools, say) the persons are in."))
  }
  if (!is.null(cluster) && identical(cluster, id)) {
    bad_argument("`cluster` must name a different column from `id`.")
  }
  if (!is.logical(random_diagonal) || length(random_diagonal) != 1L ||
      is.na(random_diagonal)) {
    bad_argument("`random_diagonal` must be TRUE or FALSE.")
  }
  invisible(NULL)
}

#' The cluster level of a multilevel growth mixture
#'
#' Every person must sit in one cluster. The group-class logits read
#' cluster-level covariates (`group_membership`, constant within clusters);
#' the trajectory-class logits have one intercept per group class and the
#' person covariates' slopes (`membership`, constant within persons), shared
#' across group classes, as in the two-level regression mixture.
#'
#' @param cluster_values The cluster of every kept row.
#' @param group_index The person of every kept row.
#' @param row_membership The row-level person membership design (with its
#'   intercept).
#' @param row_group_membership The row-level cluster membership design.
#' @param n_group_classes Number of group classes.
#' @return A list of cluster fields for the specification: `cluster_index`
#'   (each person's cluster), `cluster_levels`, `n_clusters`,
#'   `n_cluster_classes`, `cluster_v`, `cluster_w` and
#'   `cluster_intercept_only`.
#' @noRd
.growth_cluster_spec <- function(cluster_values, group_index, row_membership,
                                 row_group_membership, n_group_classes) {
  first <- match(seq_len(max(group_index)), group_index)
  person_cluster <- cluster_values[first]
  if (any(cluster_values != person_cluster[group_index])) {
    stop(errorCondition(paste(
      "Every person must belong to one cluster: `cluster` varies within some",
      "persons (`id`)."), class = "latents_bad_data", call = NULL))
  }
  cluster_levels <- sort(unique(person_cluster))
  cluster_index <- match(person_cluster, cluster_levels)
  if (length(cluster_levels) <= n_group_classes) {
    stop(errorCondition(sprintf("%d clusters cannot support %d group classes.",
                                length(cluster_levels), n_group_classes),
                        class = "latents_bad_data", call = NULL))
  }
  person_membership <- .mixture_group_design(row_membership, group_index, "membership")
  cluster_v <- .mixture_group_design(row_group_membership, cluster_index[group_index],
                                     "group_membership")
  list(cluster_index = cluster_index, cluster_levels = cluster_levels,
       n_clusters = length(cluster_levels), n_cluster_classes = n_group_classes,
       cluster_v = cluster_v,
       cluster_w = person_membership[, colnames(person_membership) != "(Intercept)",
                                     drop = FALSE],
       cluster_intercept_only = ncol(cluster_v) == 1L &&
         identical(colnames(cluster_v), "(Intercept)"))
}

#' Starting values of the growth mixture model from a latent class growth fit
#'
#' The class trajectories, class shares and memberships come from the fit
#' without random effects. The random-effect covariance starts from the spread
#' of per-person least-squares deviations from the class trajectory, halved,
#' and the residual variance from half the class residual variance, because
#' without random effects the residual absorbs all individual variation.
#'
#' @param spec The model specification.
#' @param stats Per-person cross-products.
#' @param result One EM result of the model without random effects.
#' @param first Whether this is the first start, which splits the clusters of
#'   a multilevel growth mixture by their composition rather than at random.
#' @return A growth-mixture parameter list.
#' @noRd
.growth_start <- function(spec, stats, result, first = FALSE) {
  params <- result$params
  posterior <- result$expectation$group_tau %||% result$expectation$tau
  modal <- max.col(posterior, ties.method = "first")
  n_random <- ncol(spec$random_design)
  rows <- split(seq_len(spec$n), spec$group_index)
  design <- cbind(spec$x, spec$z)
  y <- spec$y - spec$offset
  deviations <- t(vapply(seq_along(rows), function(i) {
    at <- rows[[i]]
    if (length(at) <= n_random) return(rep(NA_real_, n_random))
    coefficients <- .growth_coefficients(params, modal[i])
    residual <- y[at] - as.vector(design[at, , drop = FALSE] %*% coefficients)
    random <- spec$random_design[at, , drop = FALSE]
    fit <- tryCatch(qr.coef(qr(random), residual), error = function(e) NULL)
    if (is.null(fit) || anyNA(fit)) rep(NA_real_, n_random) else fit
  }, numeric(n_random))) |> matrix(length(rows), n_random)
  fallback <- diag(max(stats::var(y) * 0.1, spec$min_variance), n_random)
  start_covariance <- function(members) {
    usable <- deviations[members, , drop = FALSE]
    usable <- usable[stats::complete.cases(usable), , drop = FALSE]
    if (nrow(usable) <= n_random + 1L) return(fallback)
    covariance <- stats::cov(usable) / 2
    eigen_values <- eigen(covariance, symmetric = TRUE, only.values = TRUE)$values
    if (min(eigen_values) <= 1e-8 * max(1, max(abs(eigen_values)))) {
      covariance <- covariance + diag(max(1e-4 * max(abs(diag(covariance))),
                                          spec$min_variance), n_random)
    }
    if (isTRUE(spec$random_diagonal)) covariance <- diag(diag(covariance), n_random)
    covariance
  }
  everyone <- seq_along(rows)
  params$random_covariance <- if (identical(spec$random_covariance, "varying")) {
    lapply(seq_len(spec$n_classes), function(k) start_covariance(which(modal == k)))
  } else list(start_covariance(everyone))
  params$random_scale <- rep(1, spec$n_classes)
  params$sigma2 <- pmax(params$sigma2 / 2, spec$min_variance)
  if (is.null(spec$cluster)) params else
    .growth_two_level_start(spec, params, posterior, first)
}

#' Starting membership logits of a multilevel growth mixture
#'
#' The first start splits the clusters into group classes by how much of the
#' first class their persons hold (quantile groups), since group classes
#' exist to differ in that composition; every other start splits them at
#' random. Each group class's trajectory-class logits start from the class
#' shares of its clusters' persons under the start's posteriors, and the
#' covariate slopes at zero. A start with identical group classes would be a
#' fixed point of EM, which the split avoids.
#'
#' @param spec The growth specification, with its cluster fields.
#' @param params Starting parameters with person-level `gamma`.
#' @param posterior Persons-by-classes starting posteriors.
#' @param first Whether to split by composition.
#' @return `params` with `delta` and `class_logits` in place of `gamma`.
#' @noRd
.growth_two_level_start <- function(spec, params, posterior, first = FALSE) {
  n_classes <- spec$n_classes
  n_group_classes <- spec$n_cluster_classes
  labels <- if (first) {
    composition <- as.vector(tapply(posterior[, 1L], spec$cluster_index, mean))
    ranks <- rank(composition, ties.method = "first")
    as.integer(ceiling(ranks * n_group_classes / spec$n_clusters))
  } else sample(rep_len(seq_len(n_group_classes), spec$n_clusters))
  person_label <- labels[spec$cluster_index]
  logit <- function(counts) {
    shares <- pmax(counts, .Machine$double.xmin)
    log(shares) - log(shares[1L])
  }
  intercepts <- t(vapply(seq_len(n_group_classes), function(h) {
    logit(colSums(posterior[person_label == h, , drop = FALSE]) + 0.5)
  }, numeric(n_classes))) |> matrix(n_group_classes, n_classes)
  params$class_logits <- rbind(intercepts,
                               matrix(0, ncol(spec$cluster_w), n_classes))
  params$delta <- rbind(logit(tabulate(labels, n_group_classes)),
                        matrix(0, ncol(spec$cluster_v) - 1L, n_group_classes))
  params$gamma <- NULL
  params
}

#' A start from the persons' own trajectories
#'
#' Each person's least-squares coefficients on the class-specific design
#' (the pooled ones for a person with too few observations), standardized
#' and grouped by Ward's method -- deterministic, as the first start of the
#' profile models is -- give the classes; one M-step of the model without
#' random effects from those classes gives the trajectories.
#'
#' @param spec The growth specification.
#' @return A seed as `.mixture_run_starts()` returns them: `params` and an
#'   `expectation` with the persons' (softened) class weights `group_tau`.
#' @noRd
.growth_trajectory_seed <- function(spec) {
  # Shared (`common`) terms are left out of a start; EM estimates them.
  y <- spec$y - spec$offset
  pooled <- qr.coef(qr(spec$x), y)
  pooled[is.na(pooled)] <- 0
  rows <- split(seq_len(spec$n), spec$group_index)
  coefficients <- t(vapply(rows, function(at) {
    if (length(at) <= ncol(spec$x)) return(pooled)
    fitted <- qr.coef(qr(spec$x[at, , drop = FALSE]), y[at])
    ifelse(is.na(fitted), pooled, fitted)
  }, numeric(ncol(spec$x)))) |> matrix(length(rows), ncol(spec$x))
  spread <- apply(coefficients, 2L, stats::sd)
  usable <- is.finite(spread) & spread > 0
  standardized <- if (any(usable)) {
    scale(coefficients[, usable, drop = FALSE], scale = spread[usable])
  } else matrix(0, length(rows), 1L)
  labels <- .multilpa_ward_assignments(unname(as.matrix(standardized)),
                                       spec$n_classes)
  expectation <- .mixture_seed_expectation(spec, labels)
  list(params = .mixture_maximization(spec, .mixture_null_params(spec), expectation),
       expectation = expectation)
}

#' Fit the growth mixture model from several starts
#'
#' Every start of the model without random effects becomes a start here;
#' starts are screened for 50 iterations and the better half run to
#' convergence, as for `mixture_regression()`.
#'
#' @return The fitted `latents_growth_mixture` object (without inference).
#' @noRd
.growth_fit <- function(spec, n_starts, max_iter, tol, select_start) {
  stats <- .growth_statistics(spec)
  # EM to a loose tolerance, then the quasi-Newton finish to `tol`: EM for a
  # mixed model crawls near its maximum.
  em_tol <- max(tol, 1e-6)
  # The starts of the model without random effects, then one start from the
  # persons' own trajectories: without random effects the classes often
  # split persons by their level, which a random intercept absorbs, so every
  # seed of that model can sit in the same wrong basin of this one.
  seeds <- c(.mixture_run_starts(spec, n_starts, min(max_iter, 200L), 1e-6),
             list(.growth_trajectory_seed(spec)))
  screening <- min(max_iter, 50L)
  screened <- lapply(seq_along(seeds), function(index) {
    start <- tryCatch(.growth_start(spec, stats, seeds[[index]], first = index == 1L),
                      error = function(e) e)
    if (inherits(start, "error")) {
      return(list(failed = conditionMessage(start)))
    }
    result <- .growth_em(spec, stats, start, screening, em_tol)
    result$stage <- "screened"
    result
  })
  failed <- vapply(screened, function(r) !is.null(r$failed), logical(1))
  if (all(failed)) {
    stop(errorCondition(sprintf("Every start failed: %s", screened[[1L]]$failed),
                        class = "latents_no_valid_start", call = NULL))
  }
  value <- vapply(screened, function(r) if (is.null(r$failed))
    r$expectation$log_likelihood else -Inf, numeric(1))
  finalists <- utils::head(order(-value), max(2L, ceiling(length(screened) / 2)))
  results <- lapply(seq_along(screened), function(start) {
    result <- screened[[start]]
    if (!is.null(result$failed)) return(result)
    if (!start %in% finalists || result$converged || max_iter <= screening) {
      if (start %in% finalists || result$converged) result$stage <- "completed"
      return(result)
    }
    continued <- .growth_em(spec, stats, result$params, max_iter - screening, em_tol)
    continued$iterations <- continued$iterations + result$iterations
    continued$history <- c(result$history, continued$history[-1L])
    continued$stage <- "completed"
    continued
  })
  results <- lapply(results, function(result) {
    if (!is.null(result$failed) || !identical(result$stage, "completed") ||
        isTRUE(result$degenerate)) return(result)
    finish <- .growth_quasi_newton(spec, stats, result$params, tol)
    result$params <- finish$params
    result$expectation <- .growth_expectation(spec, stats, finish$params)
    result$converged <- isTRUE(result$converged) && finish$converged
    result$history <- c(result$history, result$expectation$log_likelihood)
    result
  })
  .growth_select(spec, stats, results, select_start)
}

#' Choose among growth-mixture starts and assemble the fit
#' @noRd
.growth_select <- function(spec, stats, results, select_start) {
  usable_result <- function(r) is.null(r$failed)
  log_likelihood <- vapply(results, function(r) if (usable_result(r))
    r$expectation$log_likelihood else -Inf, numeric(1))
  converged <- vapply(results, function(r) usable_result(r) && isTRUE(r$converged),
                      logical(1))
  degenerate <- vapply(results, function(r) !usable_result(r) || isTRUE(r$degenerate),
                       logical(1))
  usable <- is.finite(log_likelihood) & !degenerate
  if (!any(usable)) {
    stop(errorCondition(sprintf(paste(
      "Every one of the %d starts degenerated: a class emptied or a residual",
      "variance collapsed. Fewer classes or more starts may help."), length(results)),
      class = "latents_no_valid_start", call = NULL))
  }
  candidates <- which(usable)
  best <- candidates[.multilpa_select_start(log_likelihood[candidates],
                                            converged[candidates], select_start)]
  chosen <- results[[best]]
  if (!chosen$converged) {
    warning(warningCondition(sprintf(paste(
      "The selected start did not converge in %d iterations; the estimates",
      "are not a maximum. Raise `max_iter`."), chosen$iterations),
      class = "latents_unconverged", call = NULL))
  }
  params <- .growth_order_classes(spec, chosen$params, chosen$expectation)
  expectation <- .growth_expectation(spec, stats, params)
  degenerate_classes <- .growth_degenerate_covariances(spec, params)
  if (length(degenerate_classes) > 0L) {
    intercept_only <- identical(colnames(spec$random_design), "(Intercept)")
    advice <- if (intercept_only) {
      paste("The random intercept adds nothing there: drop `random`, or keep it",
            "only if another class needs it.")
    } else {
      paste("Simplify the random effects, for example `random = \"intercept\"`",
            if (identical(spec$random_covariance, "varying"))
              "or `random_covariance = \"equal\"`." else ".")
    }
    warning(warningCondition(sprintf(paste(
      "The random-effect covariance of %s is degenerate (a variance at its floor",
      "or random effects perfectly correlated): the data cannot separate those",
      "random effects. %s Its standard errors are not reported."),
      paste(degenerate_classes, collapse = " and "), advice),
      class = "latents_random_boundary", call = NULL))
  }
  best_value <- log_likelihood[best]
  completed <- vapply(results, function(r) identical(r$stage, "completed"), logical(1))
  starts <- data.frame(
    start = seq_along(results),
    kind = c("residuals", rep("random", length(results) - 2L), "trajectories"),
    log_likelihood = log_likelihood, converged = converged,
    iterations = vapply(results, function(r) r$iterations %||% 0L, integer(1)),
    stage = vapply(results, function(r) r$stage %||% "failed", character(1)),
    degenerate = degenerate, selected = seq_along(results) == best)
  structure(list(
    spec = spec, stats = stats, params = params, expectation = expectation,
    log_likelihood = expectation$log_likelihood,
    n_parameters = .growth_count_parameters(spec),
    random_boundary = length(degenerate_classes) > 0L,
    degenerate_classes = degenerate_classes,
    converged = chosen$converged, iterations = chosen$iterations,
    history = chosen$history, starts = starts,
    n_best_replicated = sum(usable & completed &
                              abs(log_likelihood - best_value) <=
                              1e-6 * (1 + abs(best_value)))),
    class = "latents_growth_mixture")
}

#' Classes whose random-effect covariance is on its boundary
#'
#' A variance at its floor, a correlation within 0.005 of +-1, or a smallest
#' eigenvalue below 1e-6 of the largest: the likelihood is flat there and
#' Wald errors of that covariance are not meaningful.
#'
#' @return The names of the degenerate classes (`"all classes"` when one
#'   shared matrix is degenerate); empty when none is.
#' @noRd
.growth_degenerate_covariances <- function(spec, params) {
  covariances <- .growth_covariances(spec, params)
  degenerate <- vapply(covariances, function(covariance) {
    variances <- diag(covariance)
    at_floor <- any(variances <= spec$min_variance * (1 + 1e-6))
    correlation <- covariance / sqrt(outer(variances, variances))
    collinear <- nrow(covariance) > 1L &&
      any(abs(correlation[upper.tri(correlation)]) > 0.995)
    values <- eigen(covariance, symmetric = TRUE, only.values = TRUE)$values
    at_floor || collinear || min(values) <= 1e-6 * max(values)
  }, logical(1))
  if (!any(degenerate)) return(character())
  if (identical(spec$random_covariance, "varying")) {
    sprintf("class %d", which(degenerate))
  } else "all classes"
}

#' Order classes by size (then first coefficient), and group classes by
#' composition, and relabel every block
#' @noRd
.growth_order_classes <- function(spec, params, expectation) {
  structure <- .growth_structure(spec)
  orders <- .mixture_class_orders(structure, params, expectation$structure)
  order <- orders$class_order
  covariances <- .growth_covariances(spec, params)
  params <- .mixture_permute(structure, params, order, orders$group_order)
  if (identical(spec$random_covariance, "varying")) {
    params$random_covariance <- params$random_covariance[order]
  }
  if (identical(spec$random_covariance, "proportional")) {
    # Re-express with the new last class as the unit of scale.
    scale <- params$random_scale[order]
    reference <- scale[length(scale)]
    params$random_scale <- scale / reference
    params$random_covariance <- list(params$random_covariance[[1L]] * reference^2)
  }
  stopifnot("relabelling must leave every class covariance unchanged" =
              isTRUE(all.equal(.growth_covariances(spec, params), covariances[order],
                               tolerance = 1e-10)))
  params
}

#' Number of free parameters of the growth mixture model
#' @noRd
.growth_count_parameters <- function(spec) {
  n_classes <- spec$n_classes
  n_random <- ncol(spec$random_design)
  per_matrix <- if (isTRUE(spec$random_diagonal)) n_random else
    n_random * (n_random + 1L) / 2L
  random <- switch(spec$random_covariance,
                   varying = n_classes * per_matrix,
                   equal = per_matrix,
                   proportional = per_matrix + n_classes - 1L)
  as.integer(n_classes * ncol(spec$x) + ncol(spec$z) + random +
               (if (identical(spec$variance, "equal")) 1L else n_classes) +
               length(.growth_structure_names(.growth_structure(spec))))
}
