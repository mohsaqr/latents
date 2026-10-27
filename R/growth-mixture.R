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
#' @param check_count Whether to require more clusters than group classes,
#'   as fitting does; prediction for a few new clusters does not.
#' @return A list of cluster fields for the specification: `cluster_index`
#'   (each person's cluster), `cluster_levels`, `n_clusters`,
#'   `n_cluster_classes`, `cluster_v`, `cluster_w` and
#'   `cluster_intercept_only`.
#' @noRd
.growth_cluster_spec <- function(cluster_values, group_index, row_membership,
                                 row_group_membership, n_group_classes,
                                 check_count = TRUE) {
  first <- match(seq_len(max(group_index)), group_index)
  person_cluster <- cluster_values[first]
  if (any(cluster_values != person_cluster[group_index])) {
    stop(errorCondition(paste(
      "Every person must belong to one cluster: `cluster` varies within some",
      "persons (`id`)."), class = "latents_bad_data", call = NULL))
  }
  cluster_levels <- sort(unique(person_cluster))
  cluster_index <- match(person_cluster, cluster_levels)
  if (check_count && length(cluster_levels) <= n_group_classes) {
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
    # A rank-deficient person design leaves NA coefficients: no deviation.
    fit <- qr.coef(qr(random), residual)
    if (anyNA(fit)) rep(NA_real_, n_random) else fit
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
  # Each person's level is read at the first value of every predictor (the
  # baseline), not at its mean: levels at the mean can coincide for classes
  # that differ in slope, and then carry only noise into the grouping. The
  # origin at the minimum keeps the start free of the predictors' units and
  # origin (review round 2: centring at the mean lost a slope-defined class
  # split, 11 log-likelihood units below the truth-started fit).
  design <- spec$x
  if ("(Intercept)" %in% colnames(design)) {
    shifted <- colnames(design) != "(Intercept)"
    design[, shifted] <- sweep(design[, shifted, drop = FALSE], 2L,
                               apply(design[, shifted, drop = FALSE], 2L, min))
  }
  pooled <- qr.coef(qr(design), y)
  pooled[is.na(pooled)] <- 0
  rows <- split(seq_len(spec$n), spec$group_index)
  coefficients <- t(vapply(rows, function(at) {
    if (length(at) <= ncol(design)) return(pooled)
    fitted <- qr.coef(qr(design[at, , drop = FALSE]), y[at])
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
    # A start that cannot be built is kept, marked failed with its message, and
    # reported by a `latents_failed_starts` warning once the fit is selected.
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
  weights <- spec$sampling_weights %||% 1
  results <- lapply(results, function(result) {
    if (!is.null(result$failed) || !identical(result$stage, "completed") ||
        isTRUE(result$degenerate)) return(result)
    finish <- .growth_quasi_newton(spec, stats, result$params, tol)
    # Newton steps on the analytic score finish what L-BFGS-B leaves, and
    # convergence is judged in standard errors (as for correlated residuals):
    # L-BFGS-B's relative stop left fits 0.1-0.35 standard errors short while
    # reporting convergence (review round 2, GR2-4).
    bounds <- .growth_correlated_bounds(spec, .growth_pack(spec, finish$params), 1)
    polished <- .growth_newton_finish(spec, stats, finish$params, weights, bounds)
    result$params <- polished$params
    result$expectation <- .growth_expectation(spec, stats, polished$params)
    # Where the information gives no standard-error scale (it is singular,
    # or the random-effect covariance is degenerate and has no Wald errors),
    # convergence is the quasi-Newton finish's own.
    judged <- is.finite(polished$displacement) &&
      length(.growth_degenerate_covariances(spec, polished$params)) == 0L
    result$converged <- isTRUE(result$converged) &&
      (if (judged) polished$converged else finish$converged)
    result$history <- c(result$history, result$expectation$log_likelihood)
    result
  })
  .growth_select(spec, stats, results, select_start)
}

#' Put a growth mixture on standard coordinates for fitting
#'
#' The fit is not free of the units of the data on its own. EM stops on a
#' change in the log likelihood, which shifts by `-N log s` when the outcome's
#' units change; the quasi-Newton finish works in raw coefficients and
#' covariance factors; and the variance floor, the Cholesky jitter and the
#' degeneracy rule compare random-effect variances whose units are those of
#' the random design. At `y * 1e-5` a fit stopped 0.06 standard errors short
#' (review 2026-10-26, G2); with time in days instead of years a random-slope
#' variance was pinned at the floor, 11.3 log-likelihood units low and
#' reported converged; with time as a calendar year EM stopped with "likelihood
#' decreased" (round 2, GR2-1 to GR2-4).
#'
#' The model is therefore fitted on standard coordinates and mapped back
#' exactly by [.growth_unstandardize_fit()]: the outcome becomes
#' `(y - offset - center) / scale`, and every column of the class-specific,
#' common and random designs except an intercept is centred and divided by its
#' standard deviation. A column is centred only where an intercept absorbs the
#' centre: class-specific columns when the class design has an intercept,
#' common columns when either design has one, random columns when the random
#' design has an intercept and the covariance is not restricted to be
#' diagonal (centring would mix a diagonal covariance). An uncentred column is
#' divided by its root mean square. The variance floor is carried into the
#' outcome's standard units, so it applies to the residual variance and to
#' each random effect's conditional variance in these coordinates.
#'
#' @param spec The growth-mixture specification.
#' @return A list with the standardized `spec` and the transform: `center`,
#'   `scale` and `center_on` (`"beta"`, `"common"` or `"none"`) of the
#'   outcome, and `x_center`, `x_spread`, `z_center`, `z_spread`,
#'   `random_center`, `random_spread` of the design columns.
#' @noRd
.growth_standardize <- function(spec) {
  carrying <- (spec$row_weights %||% rep(1, spec$n)) > 0
  response <- (spec$y - spec$offset)[carrying]
  x_intercept <- "(Intercept)" %in% colnames(spec$x)
  z_intercept <- "(Intercept)" %in% colnames(spec$z)
  center_on <- if (x_intercept) "beta" else if (z_intercept) "common" else "none"
  center <- if (identical(center_on, "none")) 0 else mean(response)
  scale <- stats::sd(response)
  # A constant outcome is refused before this; a single row has no spread.
  if (!is.finite(scale) || scale <= 0) scale <- 1
  columns <- function(design, centre) {
    if (ncol(design) == 0L) return(list(center = numeric(), spread = numeric()))
    pieces <- lapply(seq_len(ncol(design)), function(j) {
      values <- design[carrying, j]
      constant <- max(values) == min(values)
      if (identical(colnames(design)[j], "(Intercept)") || constant) return(c(0, 1))
      if (centre) c(mean(values), stats::sd(values)) else c(0, sqrt(mean(values^2)))
    })
    list(center = vapply(pieces, `[`, numeric(1), 1L),
         spread = vapply(pieces, `[`, numeric(1), 2L))
  }
  random_intercept <- "(Intercept)" %in% colnames(spec$random_design)
  x_columns <- columns(spec$x, x_intercept)
  z_columns <- columns(spec$z, x_intercept || z_intercept)
  random_columns <- columns(spec$random_design,
                            random_intercept && !isTRUE(spec$random_diagonal))
  standard <- list(center = center, scale = scale, center_on = center_on,
                   x_center = x_columns$center, x_spread = x_columns$spread,
                   z_center = z_columns$center, z_spread = z_columns$spread,
                   random_center = random_columns$center,
                   random_spread = random_columns$spread)
  c(list(spec = .growth_apply_standard(spec, standard)), standard)
}

#' Apply a stored standard transform to a growth specification
#'
#' Used for the fitted data and for new persons (prediction), so both are in
#' the coordinates of the fitted standard parameters.
#'
#' @param spec A growth specification in the data's units.
#' @param standard The transform from [.growth_standardize()].
#' @return The specification in standard coordinates.
#' @noRd
.growth_apply_standard <- function(spec, standard) {
  recode <- function(design, center, spread) {
    if (ncol(design) == 0L) return(design)
    design[] <- (design - rep(center, each = nrow(design))) /
      rep(spread, each = nrow(design))
    design
  }
  standardized <- spec
  standardized$y <- (spec$y - standard$center) / standard$scale
  standardized$offset <- spec$offset / standard$scale
  standardized$min_variance <- spec$min_variance / standard$scale^2
  standardized$x <- recode(spec$x, standard$x_center %||% numeric(ncol(spec$x)),
                           standard$x_spread %||% rep(1, ncol(spec$x)))
  standardized$z <- recode(spec$z, standard$z_center %||% numeric(ncol(spec$z)),
                           standard$z_spread %||% rep(1, ncol(spec$z)))
  standardized$random_design <- recode(
    spec$random_design,
    standard$random_center %||% numeric(ncol(spec$random_design)),
    standard$random_spread %||% rep(1, ncol(spec$random_design)))
  standardized
}

#' The map of the random effects from standard to original coordinates
#'
#' With `Z = Z* M` (M carries the random design's centres and spreads) and the
#' outcome divided by `scale`, `Z b = scale Z* b*`, so `b = N b*` with
#' `N = scale M^-1` and `G = N G* N'`.
#'
#' @param names The random design's column names.
#' @param standard The transform from [.growth_standardize()].
#' @return The q x q matrix `N`.
#' @noRd
.growth_random_map <- function(names, standard) {
  n_random <- length(names)
  spread <- standard$random_spread %||% rep(1, n_random)
  center <- standard$random_center %||% numeric(n_random)
  m <- diag(spread, n_random)
  intercept <- match("(Intercept)", names)
  if (!is.na(intercept)) m[intercept, -intercept] <- center[-intercept]
  standard$scale * solve(m)
}

#' Map class-specific and common coefficients to original coordinates
#'
#' A centred column's centre times its coefficient moves into the intercept
#' that absorbs it (each class's own, or the common one).
#'
#' @param beta Class-specific coefficients (columns of x by classes).
#' @param common Common coefficients.
#' @param x_names,z_names The design column names.
#' @param standard The transform from [.growth_standardize()].
#' @return A list with `beta` and `common` in the data's units.
#' @noRd
.growth_unstandardize_coefficients <- function(beta, common, x_names, z_names, standard) {
  scale <- standard$scale
  x_center <- standard$x_center %||% numeric(length(x_names))
  x_spread <- standard$x_spread %||% rep(1, length(x_names))
  z_center <- standard$z_center %||% numeric(length(z_names))
  z_spread <- standard$z_spread %||% rep(1, length(z_names))
  slopes <- beta / x_spread
  common_slopes <- common / z_spread
  # What the centred common columns contribute at the origin of the data.
  common_shift <- sum(z_center * common_slopes)
  out_beta <- scale * slopes
  out_common <- scale * common_slopes
  x_intercept <- match("(Intercept)", x_names)
  z_intercept <- match("(Intercept)", z_names)
  if (!is.na(x_intercept)) {
    shift <- colSums(x_center * slopes) + common_shift
    out_beta[x_intercept, ] <- scale * (beta[x_intercept, ] - shift) +
      (if (identical(standard$center_on, "beta")) standard$center else 0)
  } else if (!is.na(z_intercept)) {
    out_common[z_intercept] <- scale * (common[z_intercept] - common_shift) +
      (if (identical(standard$center_on, "common")) standard$center else 0)
  }
  list(beta = out_beta, common = out_common)
}

#' Map growth-mixture parameters from standard coordinates to the data's
#'
#' Coefficients by [.growth_unstandardize_coefficients()], each random-effect
#' covariance as `N G* N'` ([.growth_random_map()]), residual variances by
#' `scale^2` and a residual process's standard deviation by `scale`; class
#' shares, rates and the proportional-covariance ratios are unchanged.
#'
#' @param params A parameter list in standard coordinates.
#' @param standard The list from [.growth_standardize()].
#' @param spec The specification in the data's units (for the column names).
#' @return The parameter list in the data's units.
#' @noRd
.growth_unstandardize_params <- function(params, standard, spec) {
  coefficients <- .growth_unstandardize_coefficients(
    params$beta, params$common, colnames(spec$x), colnames(spec$z), standard)
  params$beta[] <- coefficients$beta
  params$common[] <- coefficients$common
  params$sigma2 <- params$sigma2 * standard$scale^2
  map <- .growth_random_map(colnames(spec$random_design), standard)
  params$random_covariance <- lapply(params$random_covariance, function(covariance) {
    mapped <- map %*% covariance %*% t(map)
    (mapped + t(mapped)) / 2
  })
  if (!is.null(params$correlation) && "sd" %in% names(params$correlation)) {
    params$correlation[["sd"]] <- params$correlation[["sd"]] * standard$scale
  }
  params
}

#' Each person's pattern of random-effect posterior covariance
#' @param spec A growth specification.
#' @param stats Its statistics.
#' @return One pattern index per person, indexing `expectation$covariances`.
#' @noRd
.growth_person_patterns <- function(spec, stats) {
  if (!.growth_is_correlated(spec)) return(stats$pattern)
  pattern <- integer(length(stats$n))
  invisible(lapply(seq_along(stats$dense), function(g) {
    pattern[stats$dense[[g]]$members] <<- g
  }))
  pattern
}

#' Carry an E-step from standard coordinates to the data's
#'
#' Each person's class log densities move by `-n_i log scale`, which leaves
#' the posteriors unchanged; the random-effect posterior means and
#' covariances map by `N` ([.growth_random_map()]). Computing the E-step in
#' standard coordinates keeps it accurate when the data's coordinates are
#' ill-conditioned (a calendar-year random slope).
#'
#' @param spec,stats The specification and statistics in the data's units.
#' @param params The parameters in the data's units (membership logits).
#' @param standard_spec,standard_stats The same in standard coordinates.
#' @param expectation The E-step in standard coordinates.
#' @param standard The transform.
#' @return The E-step in the data's units (without the dense terms).
#' @noRd
.growth_map_expectation <- function(spec, stats, params, standard_spec, standard_stats,
                                    expectation, standard) {
  log_density <- expectation$log_density - stats$n * log(standard$scale)
  out <- .growth_structure_terms(spec, params, log_density)
  map <- .growth_random_map(colnames(spec$random_design), standard)
  out$means <- lapply(expectation$means, function(means) means %*% t(map))
  # Persons sharing a pattern in the data's units share it in standard ones
  # (each is a function of Z'Z or of the times and design rows).
  original <- .growth_person_patterns(spec, stats)
  standard_pattern <- .growth_person_patterns(standard_spec, standard_stats)
  representative <- standard_pattern[match(seq_len(max(original)), original)]
  out$covariances <- lapply(expectation$covariances, function(covariances) {
    array(vapply(representative, function(g) {
      map %*% matrix(covariances[, , g], dim(covariances)[1L]) %*% t(map)
    }, matrix(0, nrow(map), nrow(map))), c(nrow(map), nrow(map), length(representative)))
  })
  out$dense <- NULL
  out
}

#' Return a growth mixture fitted in standard coordinates to the data's units
#'
#' The parameters are mapped back exactly and the E-step is carried over
#' ([.growth_map_expectation()]). The log likelihood differs from the
#' standard-scale one by a constant, `-(sum of row weights) log scale`, so each
#' start's value and the iteration history move by the same shift and keep
#' their order. The standard parameters and transform are kept in
#' `fit$standard`, where inference and prediction are computed.
#'
#' @param fit The fit made on `standard$spec`.
#' @param spec The original specification (with any residual-process fields).
#' @param standard The list from [.growth_standardize()].
#' @return The fit in the data's units.
#' @noRd
.growth_unstandardize_fit <- function(fit, spec, standard) {
  params <- .growth_unstandardize_params(fit$params, standard, spec)
  stats <- .growth_statistics(spec)
  expectation <- .growth_map_expectation(spec, stats, params, fit$spec, fit$stats,
                                         fit$expectation, standard)
  shift <- expectation$log_likelihood - fit$expectation$log_likelihood
  fit$standard <- c(standard[setdiff(names(standard), "spec")],
                    list(params = fit$params))
  fit$spec <- spec
  fit$stats <- stats
  fit$params <- params
  fit$expectation <- expectation
  fit$log_likelihood <- expectation$log_likelihood
  fit$history <- fit$history + shift
  fit$starts$log_likelihood <- fit$starts$log_likelihood + shift
  fit
}

#' Map growth-mixture parameters from the data's coordinates to standard ones
#'
#' The exact inverse of [.growth_unstandardize_params()], for parameters that
#' were set in the data's units (by a caller, or a fit made before the
#' standard coordinates existed).
#'
#' @param params A parameter list in the data's units.
#' @param standard The transform.
#' @param spec The specification in the data's units.
#' @return The parameter list in standard coordinates.
#' @noRd
.growth_standardize_params <- function(params, standard, spec) {
  n_x <- ncol(spec$x)
  n_classes <- ncol(params$beta)
  n_common <- length(params$common)
  size <- n_x * n_classes + n_common
  image <- function(v) {
    mapped <- .growth_unstandardize_coefficients(
      matrix(v[seq_len(n_x * n_classes)], n_x, n_classes),
      v[n_x * n_classes + seq_len(n_common)], colnames(spec$x), colnames(spec$z),
      standard)
    c(mapped$beta, mapped$common)
  }
  origin <- image(numeric(size))
  linear <- vapply(seq_len(size), function(j) image(replace(numeric(size), j, 1)) - origin,
                   numeric(size))
  solved <- solve(matrix(linear, size, size), c(params$beta, params$common) - origin)
  params$beta[] <- solved[seq_len(n_x * n_classes)]
  params$common[] <- solved[n_x * n_classes + seq_len(n_common)]
  params$sigma2 <- params$sigma2 / standard$scale^2
  inverse <- solve(.growth_random_map(colnames(spec$random_design), standard))
  params$random_covariance <- lapply(params$random_covariance, function(covariance) {
    mapped <- inverse %*% covariance %*% t(inverse)
    (mapped + t(mapped)) / 2
  })
  if (!is.null(params$correlation) && "sd" %in% names(params$correlation)) {
    params$correlation[["sd"]] <- params$correlation[["sd"]] / standard$scale
  }
  params
}

#' A fit's parameters in its standard coordinates
#'
#' The stored standard parameters when the fit's own parameters are their
#' image (as fitted); otherwise, when `params` was set afterwards, the exact
#' inverse map of `params`, so the two never disagree.
#'
#' @param fit A `latents_growth_mixture` fit with `standard`.
#' @return A parameter list in standard coordinates.
#' @noRd
.growth_standard_params <- function(fit) {
  standard <- fit$standard
  stored <- standard$params
  image <- .growth_pack(fit$spec, .growth_unstandardize_params(stored, standard, fit$spec))
  current <- .growth_pack(fit$spec, fit$params)
  same <- all(abs(image - current) <= 1e-12 * pmax(abs(image), abs(current), 1))
  if (same) stored else .growth_standardize_params(fit$params, standard, fit$spec)
}

#' The E-step of a fitted growth mixture on a specification in data units
#'
#' Computed in the fit's standard coordinates and carried back, so an
#' ill-conditioned design (a calendar-year random slope) loses no accuracy.
#'
#' @param object A `latents_growth_mixture` fit.
#' @param spec A specification in the data's units (the fit's, or new persons').
#' @param stats Its statistics.
#' @return The E-step in the data's units.
#' @noRd
.growth_fitted_expectation <- function(object, spec, stats = .growth_statistics(spec)) {
  standard <- object$standard
  if (is.null(standard)) return(.growth_expectation(spec, stats, object$params))
  standard_spec <- .growth_apply_standard(spec, standard)
  standard_stats <- .growth_statistics(standard_spec)
  expectation <- .growth_expectation(standard_spec, standard_stats,
                                     .growth_standard_params(object))
  .growth_map_expectation(spec, stats, object$params, standard_spec, standard_stats,
                          expectation, standard)
}

#' The fit's standard coordinates: specification, statistics and parameters
#' @param fit A `latents_growth_mixture` fit with `standard`.
#' @return A list with `spec`, `stats` and `params`.
#' @noRd
.growth_standard_view <- function(fit) {
  spec <- .growth_apply_standard(fit$spec, fit$standard)
  list(spec = spec, stats = .growth_statistics(spec), params = .growth_standard_params(fit))
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
  failed <- vapply(results, function(r) !usable_result(r), logical(1))
  if (any(failed)) {
    # As in every other family: the fit stands on the surviving starts, and
    # the failed ones are named, not only marked in the starts table.
    warning(warningCondition(sprintf(paste(
      "%d of %d starts failed (%s); see get_results(fit, \"starts\")."),
      sum(failed), length(results),
      paste(unique(vapply(results[failed], `[[`, character(1), "failed")),
            collapse = "; ")), class = "latents_failed_starts", call = NULL))
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
#' Judged on the free matrices (one per class, or the shared one, which for a
#' proportional covariance is the last class's) in the standard coordinates
#' the model is fitted in ([.growth_standardize()]), where the random design's
#' columns have unit spread: a random effect's conditional variance (its
#' squared Cholesky diagonal) at the floor, a correlation within 0.005 of
#' +-1, or a smallest eigenvalue below 1e-6 of the largest, of the covariance
#' or of its correlation matrix. The likelihood is flat there and Wald errors
#' of that covariance are not meaningful. On a covariance in the data's own
#' units the eigenvalue rule compared variances in different units: a random
#' slope with time in hours was flagged (review round 2, GR2-3).
#'
#' @param spec The specification in standard coordinates.
#' @param params The parameters in standard coordinates.
#' @return The names of the degenerate classes (`"all classes"` when one
#'   shared matrix is degenerate); empty when none is.
#' @noRd
.growth_degenerate_covariances <- function(spec, params) {
  degenerate <- vapply(params$random_covariance, function(covariance) {
    variances <- diag(covariance)
    conditional <- diag(.growth_factor(covariance))^2
    at_floor <- any(conditional <= spec$min_variance * (1 + 1e-6))
    correlation <- covariance / sqrt(outer(variances, variances))
    collinear <- nrow(covariance) > 1L &&
      any(abs(correlation[upper.tri(correlation)]) > 0.995)
    ratio <- function(m) {
      values <- eigen(m, symmetric = TRUE, only.values = TRUE)$values
      min(values) <= 1e-6 * max(values)
    }
    at_floor || collinear || ratio(covariance) || ratio(correlation)
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
               length(.growth_correlation_names(spec)) +
               (if (identical(spec$variance, "equal")) 1L else n_classes) +
               length(.growth_structure_names(.growth_structure(spec))))
}
