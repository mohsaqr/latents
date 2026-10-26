# Adversarial Gaussian mixtures against mclust (Scrucca, Fraley, Murphy &
# Raftery 2023, Model-Based Clustering, Classification, and Density Estimation
# Using mclust in R. Chapman & Hall/CRC).
#
# The `mclust` suite compares well-behaved published data. This one simulates
# data built to break an implementation -- a component of nine rows, heavy
# overlap, near-singular covariances, a hundredfold variance ratio, a large
# offset, indicators eleven orders of magnitude apart, more parameters than
# rows, ties, and a near-degenerate univariate component -- and makes two
# separate claims on each, so that a computation defect and a search defect
# cannot be mistaken for each other:
#
#   (A) computation: started from mclust's solution, carried by mclust's own EM
#       to its fixed point (`.hard_polish()`), the package stays at that mode
#       and reports every output as the definitions applied to mclust's
#       estimates give it (`compare_fit()`, full ledger), plus scale-free
#       restatements of the estimates (`.hard_scale_free()`);
#   (B) search: from its own 50 random starts the package reaches a maximum no
#       lower than mclust's. A higher maximum than mclust's is recorded, after
#       checking with mclust's own EM that it is a genuine mode. A higher
#       maximum that sits on the package's variance floor is a different
#       problem from mclust's (which has no floor) and is recorded as such; the
#       search claim is then made over the starts off the floor.
#
# The fixtures are simulated here in base R from fixed seeds; no package code
# is used to build them.

#' Draw from a multivariate normal with base R
#' @param n Rows.
#' @param mean Mean vector.
#' @param sigma Covariance matrix.
#' @return An `n` x `length(mean)` matrix.
.hard_normal <- function(n, mean, sigma) {
  draws <- matrix(stats::rnorm(n * length(mean)), n) %*% chol(sigma)
  sweep(draws, 2L, mean, "+")
}

#' The adversarial fixtures
#'
#' Each element carries the data, the mclust models to fit, the number of
#' components, the variance floor to pass to the package, and why the fixture
#' is hard.
#' @return A list of fixtures.
.hard_fixtures <- function() {
  with_seed <- function(seed, code) {
    set.seed(seed)
    code
  }
  list(
    # 9 of 300 rows (3%) in a compact third component: its covariance is
    # estimated from nine points, and its weight sits near the empty boundary.
    list(name = "near-empty component", models = c("VVI", "VVV"), g = 3L,
         data = with_seed(101L, rbind(
           .hard_normal(150L, c(0, 0), diag(2)),
           .hard_normal(141L, c(4, 0), diag(c(1, 2))),
           .hard_normal(9L, c(2, 7), diag(c(0.3, 0.3)))))),
    # Two components 1.3 within-SD apart: relative entropy 0.24 (VVI) and 0.31
    # (VVV) at mclust's solutions, so most rows are split between components
    # and the likelihood surface is flat along the mixing direction.
    list(name = "heavy overlap", models = c("VVI", "VVV"), g = 2L,
         data = with_seed(202L, rbind(
           .hard_normal(150L, c(0, 0), diag(2)),
           .hard_normal(150L, c(1.2, 0.6), diag(c(1.5, 0.7)))))),
    # Within-component correlation 0.99: the covariance has condition number
    # 199, and its inverse amplifies any rounding in the residuals by as much.
    list(name = "correlation 0.99", models = c("VVV", "EEE"), g = 2L,
         data = with_seed(303L, {
           sigma <- matrix(c(1, 0.99, 0.99, 1), 2L)
           rbind(.hard_normal(150L, c(0, 0), sigma), .hard_normal(150L, c(2, -1), sigma))
         })),
    # Component variances 1 and 100.
    list(name = "variance ratio 100", models = c("VVI"), g = 2L,
         data = with_seed(404L, rbind(
           .hard_normal(150L, c(0, 0), diag(2)),
           .hard_normal(150L, c(3, 3), diag(c(100, 100)))))),
    # Means near 1e6 with SD 1: a sum of squares formed without centring would
    # lose about 1e12 * 2.2e-16 = 2e-4 of a unit variance to cancellation, and
    # a mean is represented only to 1e6 * 2.2e-16 = 2.2e-10. mclust's M-step
    # centres: its estimates here and on the data less 1e6 agree to 5.8e-10 in
    # the means and 1.1e-10 in the variances. The floor stays at its default:
    # it is absolute, the within-component variances are about 1, so 1e-6
    # sits six orders below them and is inactive (the compared
    # `gaussian_diagnostics$boundary` column states it per component).
    list(name = "offset 1e6", models = c("VVI", "VVV"), g = 2L,
         data = with_seed(505L, rbind(
           .hard_normal(150L, c(1e6, 1e6), diag(2)),
           .hard_normal(150L, c(1e6 + 4, 1e6 - 3), diag(2))))),
    # One indicator in units of 1e-4, the other of 1e5. The floor is absolute
    # (squared input units, 1e-6 by default), and the first indicator's
    # within-component variance is about 1e-8, under it. Left at its default
    # the floor would bind and the package would report a constrained maximum
    # that mclust, which has no such floor, does not; carried into the new
    # units with the smallest squared unit, as the invariance suite does, the
    # floor keeps its meaning (1e-6 of a unit variance) and stays inactive.
    #
    # mclust fails here under VVV with its defaults: Mclust() returns NULL,
    # because EM's singularity guard compares the reciprocal condition number
    # of each covariance with `eps` (machine epsilon, 2.2e-16), and these units
    # make it about 1 / 1.7e18 although the covariance is well conditioned on
    # its own scale. The guard is not scale-invariant; `eps = 1e-30` lets
    # mclust proceed, and its maximum then equals that of the standardised
    # data less n * sum(log(scale)) (-1762.303 both ways).
    list(name = "units 1e-4 and 1e5", models = c("VVI", "VVV"), g = 2L,
         min_variance = 1e-6 * (1e-4)^2, mclust_eps = 1e-30,
         data = with_seed(606L, {
           x <- rbind(.hard_normal(150L, c(0, 0), diag(2)),
                      .hard_normal(150L, c(3, 3), diag(c(1, 1.5))))
           cbind(x[, 1L] * 1e-4, x[, 2L] * 1e5)
         })),
    # 40 rows, 5 indicators, 3 components under VVV: 62 free parameters for
    # 40 rows, and each covariance is estimated from 13 or 14 points in five
    # dimensions.
    list(name = "40 rows, 62 parameters", models = c("VVV"), g = 3L,
         data = with_seed(707L, rbind(
           .hard_normal(14L, rep(0, 5L), diag(5L)),
           .hard_normal(13L, c(4, 0, 4, 0, 4), diag(5L)),
           .hard_normal(13L, c(0, 4, 0, 4, -4), diag(5L))))),
    # Values rounded to 0.5: many exact ties, so posteriors of tied rows are
    # identical and modal assignment depends on how ties are broken.
    list(name = "ties at 0.5", models = c("VVI", "VVV"), g = 2L,
         data = with_seed(808L, round(2 * rbind(
           .hard_normal(150L, c(0, 0), diag(2)),
           .hard_normal(150L, c(3, 2), diag(c(1, 2))))) / 2)),
    # Four components of unequal size and orientation, two of them overlapping
    # with opposite correlations: a search problem more than a computational
    # one. mclust's single hierarchical start reaches -1266.694 under VVV,
    # below its own VVI maximum (-1266.481) although VVV nests VVI, so its VVV
    # start misses the top mode.
    list(name = "four components, two overlapping", models = c("VVI", "VVV"), g = 4L,
         data = with_seed(1012L, rbind(
           .hard_normal(120L, c(0, 0), matrix(c(1, 0.5, 0.5, 1), 2L)),
           .hard_normal(80L, c(1, -0.6), matrix(c(0.6, -0.3, -0.3, 0.8), 2L)),
           .hard_normal(100L, c(5, 4), diag(c(2, 0.5))),
           .hard_normal(40L, c(-3, 5), diag(c(0.4, 2)))))),
    # Univariate, a component of SD 0.01 inside the range of two of SD 1: its
    # density is about 40 at its centre, against 0.4 for the others.
    list(name = "univariate SD 0.01 component", models = c("V"), g = 3L,
         data = with_seed(909L, matrix(c(stats::rnorm(140L, 0, 1),
                                         stats::rnorm(100L, 3, 1),
                                         stats::rnorm(60L, 1.5, 0.01)))))
  )
}

suite_hard_mclust <- function() {
  require_suite_packages("mclust", "numDeriv",
                         reason = "the Gaussian-mixture reference and its observed information")
  # Mclust() resolves its model-search function by name, so the namespace has to
  # be attached rather than only loaded.
  if (!"package:mclust" %in% search()) {
    suppressMessages(library(mclust))
    on.exit(detach("package:mclust", unload = FALSE), add = TRUE)
  }
  helpers <- new.env(parent = globalenv())
  sys.source(file.path("equivalence", "harness", "suites", "mclust.R"), envir = helpers)
  cases <- unlist(lapply(.hard_fixtures(), function(fixture) {
    lapply(fixture$models, function(model) {
      .hard_mclust_case(fixture$name, fixture$data, model, fixture$g,
                        fixture$min_variance %||% 1e-6,
                        fixture$mclust_eps %||% .Machine$double.eps, helpers)
    })
  }), recursive = FALSE)
  stack_cases(cases)
}

#' The package's starting list from an mclust solution
#' @param parameters mclust's `parameters`.
#' @param g,j Components and indicators.
#' @param full Whether the covariance model is full.
#' @return A list for `multilpa(start = )`.
.hard_start <- function(parameters, g, j, full) {
  means <- matrix(as.numeric(parameters$mean), g, j, byrow = TRUE)
  sigma <- if (j == 1L) array(rep_len(parameters$variance$sigmasq, g), c(1L, 1L, g)) else
    array(parameters$variance$sigma, c(j, j, g))
  start <- list(means = means, profile_probabilities = matrix(parameters$pro, 1L),
                group_probabilities = 1)
  if (full) {
    start$covariances <- sigma
  } else {
    start$variances <- t(matrix(vapply(seq_len(g), function(k) diag(matrix(sigma[, , k], j, j)),
                                       numeric(j)), j, g))
  }
  start
}

#' Scale-free comparisons of the estimates
#'
#' `compare_block()` measures a gap relative to `max(|reference|, 1)`. That is
#' strict for a unit-scale quantity and meaningless for one far from unit
#' scale: a mean of 1e6 passes a 1e-6 tolerance with a gap of one standard
#' deviation, and a variance of 1e-8 passes with a gap a hundred times its
#' size. These rows restate the estimates on their own scale -- mean gaps in
#' within-component standard deviations, variance, covariance and standard
#' error gaps relative to the reference value -- so the hard cases are held to
#' the same standard as unit-scale data.
#' @return A comparison frame (`output` is `NA`: these restate compared columns).
.hard_scale_free <- function(label, fit, reference) {
  expected <- reference_tables(align_reference(reference, fit), fit)
  obtained <- fit_tables(fit)
  profiles <- merge(expected$profiles, obtained$profiles, by = c("profile", "indicator"),
                    suffixes = c(".reference", ".obtained"))
  stopifnot("every profile row must be matched" = nrow(profiles) == nrow(expected$profiles))
  covariances <- merge(expected$covariances, obtained$covariances,
                       by = c("profile", "indicator", "indicator_2"),
                       suffixes = c(".reference", ".obtained"))
  scale_of <- function(profile, indicator) {
    profiles$variance.reference[profiles$profile == profile & profiles$indicator == indicator]
  }
  covariance_scale <- sqrt(vapply(seq_len(nrow(covariances)), function(i) {
    scale_of(covariances$profile[[i]], covariances$indicator[[i]]) *
      scale_of(covariances$profile[[i]], covariances$indicator_2[[i]])
  }, numeric(1)))
  gap <- function(quantity, standardised, tolerance) {
    compare_block(sprintf("%s: %s", label, quantity), rep(0, length(standardised)),
                  standardised, tolerance = tolerance)
  }
  relative <- function(column) {
    (profiles[[paste0(column, ".obtained")]] - profiles[[paste0(column, ".reference")]]) /
      profiles[[paste0(column, ".reference")]]
  }
  rbind(
    gap("mean gap in within-component SDs",
        (profiles$mean.obtained - profiles$mean.reference) /
          sqrt(profiles$variance.reference), 1e-6),
    gap("variance gap relative to the variance", relative("variance"), 1e-6),
    gap("covariance gap relative to sqrt(var_i var_j)",
        (covariances$covariance.obtained - covariances$covariance.reference) /
          covariance_scale, 1e-6),
    gap("mean standard error, relative gap", relative("mean_standard_error"), 1e-4),
    gap("variance standard error, relative gap", relative("variance_standard_error"), 1e-4))
}

#' One adversarial fixture and one mclust model against the package
#' @param name Fixture name.
#' @param x Data matrix.
#' @param model mclust model name.
#' @param g Number of components.
#' @param min_variance Variance floor passed to the package.
#' @param eps mclust's EM singularity tolerance (`emControl(eps = )`).
#' @param helpers Environment holding the mclust suite's helpers.
#' @return A `harness_case`.
.hard_mclust_case <- function(name, x, model, g, min_variance, eps, helpers) {
  j <- ncol(x)
  indicators <- sprintf("y%d", seq_len(j))
  colnames(x) <- indicators
  frame <- as.data.frame(x)
  frame$unit <- seq_len(nrow(frame))
  # Mclust() names columns from deparse(substitute(data)), so the argument is a
  # plain symbol; a univariate fixture goes in as a vector, as in the mclust suite.
  model_data <- if (j == 1L) as.numeric(x) else x
  label <- sprintf("%s, %s, %d components", name, model, g)
  reference <- mclust::Mclust(model_data, G = g, modelNames = model, verbose = FALSE,
                              control = mclust::emControl(tol = c(1e-15, 1e-15),
                                                          itmax = c(1e5, 1e5), eps = eps))
  if (is.null(reference)) {
    stop(errorCondition(sprintf("%s: mclust returned no solution", label),
                        class = "multilpa_bad_suite", call = NULL))
  }
  reference <- .hard_polish(.hard_refine(reference, model_data, eps), model_data)
  covariance <- helpers$.mclust_covariance(model)
  full <- identical(covariance$covariance_model, "full")
  fit_with <- function(...) {
    multilpa(frame, vars = indicators, id = "unit", n_profiles = g, n_group_classes = 1L,
             tol = 1e-15, max_iter = 20000L, min_variance = min_variance,
             covariance_model = covariance$covariance_model,
             variance_model = covariance$variance_model, ...)
  }

  # (A) Computation at mclust's mode.
  at_mode <- fit_with(start = .hard_start(reference$parameters, g, j, full), n_starts = 1L)
  reference_frame <- .hard_reference(reference, model, frame[indicators])
  moved <- as.numeric(logLik(at_mode)) - reference$loglik
  case <- compare_fit(
    sprintf("%s (A, at mclust's mode)", label), at_mode, reference_frame,
    # As in the mclust suite: both sets of standard errors are second
    # derivatives, the reference's numerical (Richardson extrapolation).
    tolerance = c(default = 1e-6, `profiles$mean_standard_error` = 1e-4,
                  `profiles$variance_standard_error` = 1e-4),
    extra = rbind(
      # How close the reference itself is to mclust's fixed point, after
      # .hard_polish(): the last EM step and the distance it implies, both on
      # the scale-free parameter scale. A comparison at 1e-6 needs the
      # reference well inside that.
      compare_values(sprintf("%s (A): mclust reference's remaining distance from its EM fixed point (%d polishing steps, %s)",
                             label, reference$polish$steps,
                             if (reference$polish$at_floor) "stopped at the rounding floor" else
                               if (is.na(reference$polish$rate)) "already there" else
                                 sprintf("rate %.4f", reference$polish$rate)),
                     # A hundredth of the 1e-6 the comparisons are made at.
                     0, reference$polish$distance, tolerance = 1e-8),
      compare_values(sprintf("%s (A): log-likelihood change from mclust's fixed point", label),
                     0, moved, tolerance = 1e-8),
      .hard_scale_free(sprintf("%s (A)", label), at_mode, reference_frame)))

  # (B) Search from the package's own starts.
  searched <- fit_with(n_starts = 50L, seed = 1L)
  search_rows <- NULL
  if (isTRUE(searched$boundary) && as.numeric(logLik(searched)) > reference$loglik) {
    # The package's floor makes a component that collapses onto a few rows a
    # finite, bound-active maximum; mclust has no floor and refuses such a
    # solution as singular, so the two maxima are of different problems. This
    # is the documented `boundary = "allow"` behaviour (with its warning), not
    # a search defect: it is recorded, and the search claim is restated over
    # the starts that keep every eigenvalue off the floor.
    smallest_eigenvalue <- min(get_results(searched, "gaussian_diagnostics")$min_standardized_eigenvalue)
    search_rows <- compare_values(
      sprintf(paste("%s (B): the default search (boundary = \"allow\") returns a bound-active",
                    "solution, log-likelihood %.4f against mclust's %.4f, smallest",
                    "standardized eigenvalue %.3g (recorded; mclust cannot represent it)"),
              label, as.numeric(logLik(searched)), reference$loglik, smallest_eigenvalue),
      1, as.numeric(searched$boundary), tolerance = 0.5)
    searched <- fit_with(n_starts = 50L, seed = 1L, boundary = "exclude")
  }
  package_ll <- as.numeric(logLik(searched))
  shortfall <- reference$loglik - package_ll
  search_rows <- rbind(search_rows, compare_values(
    sprintf("%s (B): shortfall of the package's own search against mclust", label),
    0, max(shortfall, 0), tolerance = 1e-6))
  if (-shortfall > 1e-6) {
    # mclust starts once, from model-based hierarchical clustering, and can stop
    # at a lower mode. That is not a package defect; the higher maximum is
    # checked to be genuine by continuing mclust's own EM from the package's
    # estimates and evaluating mclust's own density there.
    search_rows <- rbind(search_rows, .hard_confirm_higher(label, searched, model,
                                                           model_data, reference$loglik, eps))
  }
  case$comparisons <- rbind(case$comparisons, .with_output_columns(search_rows))
  case$comparisons <- transform(case$comparisons,
                                source = sprintf("mclust %s", utils::packageVersion("mclust")),
                                dataset = name, precision = "machine")
  case
}

#' mclust's solution in the harness's reference shape, with well-scaled SEs
#'
#' As `.mclust_reference()` in the mclust suite, except for the standard
#' errors, which come from [.hard_standard_errors()]: the mclust suite's
#' numerical Hessian starts from steps of 10% of each parameter and fails on
#' these fixtures: a step of 10% of a mean near 1e6 is 1e5 standard deviations
#' (its mean standard error came out 101.5 where the package reports 0.0761
#' and sd / sqrt(n) is about 0.08); on correlation 0.99 and on the units
#' fixture the steps leave the positive-definite cone and `chol()` fails; on
#' the near-empty fixture a step drives the dependent 3% proportion below
#' zero and every standard error is NaN.
#' @return A `harness_reference`.
.hard_reference <- function(reference, model, data) {
  indicators <- names(data)
  g <- reference$G
  j <- length(indicators)
  means <- matrix(as.numeric(reference$parameters$mean), g, j, byrow = TRUE,
                  dimnames = list(NULL, indicators))
  variance <- reference$parameters$variance
  sigma <- if (j == 1L) array(rep_len(variance$sigmasq, g), c(1L, 1L, g)) else
    array(variance[["sigma"]], c(j, j, g))
  dimnames(sigma) <- list(indicators, indicators, NULL)
  harness_reference(
    log_likelihood = reference$loglik, n_parameters = reference$df,
    posteriors = reference$z, means = means, covariances = sigma,
    standard_errors = .hard_standard_errors(means, sigma, reference$parameters$pro,
                                            model, as.matrix(data)),
    profile_probabilities = matrix(reference$parameters$pro, 1L),
    group_probabilities = 1,
    group_log_likelihood = stats::setNames(
      mclust::dens(as.matrix(data), modelName = reference$modelName,
                   parameters = reference$parameters, logarithm = TRUE),
      seq_len(nrow(data))),
    data = data)
}

#' Observed-information standard errors, computed in whitened coordinates
#'
#' The from-scratch likelihood of the mclust suite's
#' `.mclust_standard_errors()` (Cholesky log densities, log-sum-exp, the free
#' natural parameters of the mclust model), evaluated on data moved by an
#' affine map `z = A (x - c)` chosen so that every parameter is of order one:
#' `c` is the data mean, and `A` scales each indicator by its pooled
#' within-component SD, or, for a full-covariance model, whitens by the pooled
#' within-component covariance (a non-diagonal `A` would break a diagonal
#' model's structure). The natural parameters in `x` are an affine function of
#' those in `z` (`mu = A^-1 mu_z + c`, `Sigma = A^-1 Sigma_z A^-T`, proportions
#' unchanged), and the log-likelihood differs only by the constant
#' `n log|det A|`, so the inverse observed information maps exactly:
#' `V_x = J V_z J'`, with `J` read off the map column by column (exact, since
#' the map is affine). The Hessian in `z` is Richardson-extrapolated
#' (numDeriv) at mclust's estimates, as in the mclust suite.
#'
#' @return A list with `mean` and `variance`, components-by-indicators.
.hard_standard_errors <- function(means, sigma, pro, model, data) {
  require_suite_packages("numDeriv", reason = "the reference observed information")
  g <- nrow(means)
  j <- ncol(means)
  shared <- substr(model, 1L, 1L) == "E"
  diagonal <- nchar(model) == 1L || substr(model, 3L, 3L) == "I"
  cells <- if (diagonal) cbind(seq_len(j), seq_len(j)) else
    which(lower.tri(diag(j), diag = TRUE), arr.ind = TRUE)
  blocks <- if (shared) 1L else g
  n_mean <- g * j
  n_cov <- nrow(cells) * blocks

  # The proportions enter through the G - 1 smallest, the largest being one
  # minus their sum. numDeriv's first Richardson step is 10% of each
  # parameter, which on the near-empty fixture (0.508, 0.462, 0.030) drives a
  # dependent 0.030 below zero when it is the one left out; with the largest
  # left out no step can. Which proportion is dependent does not change the
  # information about the means and covariances.
  free <- order(pro)[seq_len(g - 1L)]
  dependent <- setdiff(seq_len(g), free)
  pooled <- Reduce(`+`, lapply(seq_len(g), function(k) pro[[k]] * matrix(sigma[, , k], j, j)))
  centre <- colMeans(data)
  to_z <- if (diagonal) diag(1 / sqrt(diag(pooled)), j) else solve(t(chol(pooled)))
  from_z <- solve(to_z)

  pack <- function(mu, covariances, first) {
    c(as.vector(mu),
      as.vector(vapply(seq_len(blocks), function(k) covariances[[k]][cells],
                       numeric(nrow(cells)))),
      first)
  }
  unpack <- function(theta) {
    mu <- matrix(theta[seq_len(n_mean)], g, j)
    values <- matrix(theta[n_mean + seq_len(n_cov)], nrow(cells), blocks)
    covariances <- lapply(seq_len(blocks), function(k) {
      m <- matrix(0, j, j)
      m[cells] <- values[, k]
      m[cells[, 2:1, drop = FALSE]] <- values[, k]
      m
    })
    list(mu = mu, covariances = covariances, first = theta[n_mean + n_cov + seq_len(g - 1L)])
  }
  # The map from z-space parameters to x-space ones.
  to_x <- function(theta) {
    parts <- unpack(theta)
    pack(t(from_z %*% t(parts$mu) + centre),
         lapply(parts$covariances, function(m) from_z %*% m %*% t(from_z)), parts$first)
  }

  z_data <- t(to_z %*% (t(data) - centre))
  z_theta <- pack(t(to_z %*% (t(means) - centre)),
                  lapply(seq_len(blocks), function(k) to_z %*% matrix(sigma[, , k], j, j) %*% t(to_z)),
                  pro[free])
  log_likelihood <- function(theta) {
    parts <- unpack(theta)
    proportions <- numeric(g)
    proportions[free] <- parts$first
    proportions[dependent] <- 1 - sum(parts$first)
    if (any(proportions <= 0)) return(-Inf)
    by_component <- vapply(seq_len(g), function(k) {
      root <- chol(parts$covariances[[if (shared) 1L else k]])
      standardised <- backsolve(root, t(z_data) - parts$mu[k, ], transpose = TRUE)
      log(proportions[[k]]) - 0.5 * colSums(standardised^2) -
        sum(log(diag(root))) - 0.5 * j * log(2 * pi)
    }, numeric(nrow(z_data)))
    by_component <- matrix(by_component, nrow(z_data))
    peak <- apply(by_component, 1L, max)
    sum(peak + log(rowSums(exp(by_component - peak))))
  }
  information <- -numDeriv::hessian(log_likelihood, z_theta)
  z_covariance <- solve(information)
  origin <- to_x(numeric(length(z_theta)))
  jacobian <- vapply(seq_along(z_theta), function(i) {
    to_x(replace(numeric(length(z_theta)), i, 1)) - origin
  }, numeric(length(z_theta)))
  errors <- sqrt(diag(jacobian %*% z_covariance %*% t(jacobian)))

  variance_rows <- matrix(errors[n_mean + seq_len(n_cov)], nrow(cells), blocks)
  on_diagonal <- cells[, 1L] == cells[, 2L]
  variance_errors <- t(vapply(seq_len(g), function(k) {
    values <- variance_rows[on_diagonal, if (shared) 1L else k]
    values[order(cells[on_diagonal, 1L])]
  }, numeric(j)))
  list(mean = matrix(errors[seq_len(n_mean)], g, j, dimnames = dimnames(means)),
       variance = matrix(variance_errors, g, j, dimnames = dimnames(means)))
}

#' Continue mclust's EM from its own solution
#'
#' The mclust suite's `.mclust_refine()`, with mclust's singularity tolerance
#' passed through (see the units fixture), so that a fixture that needs
#' `eps` below its default can be refined at all.
#' @return `reference` with refined `parameters`, `z` and `loglik`.
.hard_refine <- function(reference, data, eps) {
  refined <- mclust::em(data = data, modelName = reference$modelName,
                        parameters = reference$parameters,
                        control = mclust::emControl(tol = 0, itmax = 20000L, eps = eps))
  stopifnot("mclust's refining EM must finish without an error code" =
              identical(as.numeric(attr(refined, "returnCode")), 0))
  log_likelihood <- sum(mclust::dens(data, modelName = reference$modelName,
                                     parameters = refined$parameters, logarithm = TRUE))
  stopifnot("refining must not lose likelihood" = log_likelihood >= reference$loglik - 1e-9)
  reference$parameters <- refined$parameters
  reference$z <- refined$z
  reference$loglik <- log_likelihood
  reference
}

#' Carry mclust's EM the rest of the way to its fixed point
#'
#' `mclust::em()` stops when the log-likelihood stops changing in floating
#' point, which `tol = 0` does not prevent: on a flat surface the likelihood
#' gain of an EM step falls under its rounding while the parameters still move.
#' On the heavy-overlap VVI fixture, after `.mclust_refine()` returned, one more
#' of mclust's own EM steps still moved a parameter by 3.3e-8 with contraction
#' rate 0.9858 per step, so the estimate sat about 3.3e-8 / (1 - 0.9858) = 2.3e-6
#' from the fixed point -- the size of the posterior gaps a 1e-6 comparison
#' then reported. This continues mclust's own E- and M-steps (`estep()`,
#' `mstep()`) until a step moves no parameter by more than 1e-13 on its own
#' scale (proportions as they are, means in within-component SDs, covariances
#' relative to `sqrt(s_ii s_jj)`). The estimate is still mclust's: the map
#' iterated is mclust's, only the stopping rule is changed. The final step size
#' and the estimated remaining distance, `step * rate / (1 - rate)`, are kept so
#' the case can report how far from the fixed point the reference still is.
#'
#' @param reference A refined `Mclust` fit.
#' @param data The data it was fitted to.
#' @param max_steps Bound on the number of EM steps.
#' @return `reference` with polished `parameters`, `z`, `loglik`, and a
#'   `polish` list (`steps`, `last_step`, `rate`, `distance`).
.hard_polish <- function(reference, data, max_steps = 50000L) {
  model <- reference$modelName
  em_step <- function(parameters) {
    z <- mclust::estep(data = data, modelName = model, parameters = parameters)$z
    mclust::mstep(data = data, modelName = model, z = z)$parameters
  }
  # Parameters on their own scale, with the scale taken from `at` so that the
  # two ends of a step are measured in one set of units: dividing a mean near
  # 1e6 by each end's own SD would turn a 1e-10 wobble in the SD into a 1e-4
  # change of the ratio.
  scale_free <- function(parameters, at) {
    sigma_of <- function(p) {
      # `[[` and not `$`: a univariate fit has `sigmasq` and no `sigma`, and
      # `$` would partially match the one for the other.
      p$variance[["sigma"]] %||%
        array(rep_len(p$variance$sigmasq, length(p$pro)), c(1L, 1L, length(p$pro)))
    }
    sigma <- sigma_of(parameters)
    scale <- sigma_of(at)
    j <- dim(sigma)[1L]
    means <- matrix(parameters$mean, j)
    c(parameters$pro, unlist(lapply(seq_along(parameters$pro), function(k) {
      sd <- sqrt(diag(matrix(scale[, , k], j)))
      c(means[, k] / sd, matrix(sigma[, , k], j) / outer(sd, sd))
    })))
  }
  state <- list(parameters = reference$parameters, steps = 0L, last_step = Inf,
                smallest = Inf, since_smallest = 0L, rate = NA_real_)
  # EM is a sequential fixed-point map: each step needs the previous one, so
  # the iteration cannot be vectorised. It stops at a step under 1e-13, or when
  # twenty steps in a row fail to shrink the step: the step has then reached
  # the rounding floor of the data (on the offset-1e6 fixture a mean is only
  # represented to 1e6 * 2.2e-16 = 2.2e-10, so 1e-13 cannot be reached).
  while (state$last_step > 1e-13 && state$since_smallest < 20L && state$steps < max_steps) {
    following <- em_step(state$parameters)
    size <- max(abs(scale_free(following, state$parameters) -
                        scale_free(state$parameters, state$parameters)))
    shrank <- size < state$smallest
    state <- list(parameters = following, steps = state$steps + 1L, last_step = size,
                  smallest = min(size, state$smallest),
                  since_smallest = if (shrank) 0L else state$since_smallest + 1L,
                  rate = if (is.finite(state$last_step) && state$last_step > 0)
                    size / state$last_step else NA_real_)
  }
  log_likelihood <- sum(mclust::dens(data, modelName = model,
                                     parameters = state$parameters, logarithm = TRUE))
  stopifnot("polishing must not lose likelihood" = log_likelihood >= reference$loglik - 1e-9)
  reference$parameters <- state$parameters
  reference$z <- mclust::estep(data = data, modelName = model,
                               parameters = state$parameters)$z
  reference$loglik <- log_likelihood
  # In the linear regime the remaining distance is step * rate / (1 - rate).
  # At the rounding floor the step no longer contracts and the distance is the
  # floor itself, taken as the smallest step seen.
  at_floor <- state$since_smallest >= 20L
  reference$polish <- list(steps = state$steps, last_step = state$last_step,
                           rate = state$rate, at_floor = at_floor,
                           distance = if (at_floor || is.na(state$rate) || state$rate >= 1)
                             state$smallest else state$last_step * state$rate / (1 - state$rate))
  reference
}

#' Confirm with mclust that the package's higher maximum is a genuine mode
#'
#' mclust's own M-step is applied to the package's posteriors (at a fixed point
#' that reproduces the package's estimates), mclust's EM is continued from there
#' without a stopping tolerance, and mclust's own density gives the maximum.
#' @return A comparison frame with `output = NA`.
.hard_confirm_higher <- function(label, searched, model, model_data, mclust_ll, eps) {
  posteriors <- unname(as.matrix(searched$subject_posteriors))
  stepped <- mclust::mstep(data = model_data, modelName = model, z = posteriors)
  continued <- mclust::em(data = model_data, modelName = model,
                          parameters = stepped$parameters,
                          control = mclust::emControl(tol = 0, itmax = 20000L, eps = eps))
  mclust_there <- sum(mclust::dens(model_data, modelName = model,
                                   parameters = continued$parameters, logarithm = TRUE))
  package_ll <- as.numeric(logLik(searched))
  compare_values(
    sprintf(paste("%s (B): mclust's own EM, continued from the package's higher mode,",
                  "keeps its log-likelihood (mclust's single start reached %.6f,",
                  "lower by %.3g; recorded, not a package defect)"),
            label, mclust_ll, package_ll - mclust_ll),
    mclust_there, package_ll, tolerance = 1e-6, scale = "relative")
}
