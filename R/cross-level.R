# Cross-level families of multilpa(family = ): individual profiles (L1) and
# group classes (L2) estimated from the same indicators, following the
# manifest-aggregation specification of Houle, Morin & Harvey (2026,
# supplement): each group's manifest means of the indicators are the
# between-level indicators of its group class.
#
#   L_j = sum_h omega_h f(m_j | h) prod_i sum_k pi_{k|h} f(x_ij | k)
#
# Restricted cross-level: pi_{k|h} = pi_k, so the likelihood factorizes into
# a profile model of the ratings and a profile model of the group means.
# Full cross-level: pi_{k|h} varies with the group class (the profile model of
# multilpa()), and the group means also inform the group classes.
# The group means are computed from the same ratings, so this is the working
# likelihood of that specification, not a likelihood of the ratings alone.

utils::globalVariables(c("group_class", "share", "profile", "level"))

#' Arguments `multilpa()` accepts for a cross-level family
#' @noRd
.cross_level_arguments <- function() {
  c("data", "vars", "id", "n_profiles", "n_group_classes", "family",
    "variance_model", "between_variance", "n_starts", "max_iter", "tol",
    "min_variance", "seed", "weights")
}

#' Class-conditional log densities of the group means, groups x classes
#' @noRd
.cross_level_group_density <- function(group_means, means, variances) {
  matrix(vapply(seq_len(nrow(means)), function(h) {
    rowSums(stats::dnorm(group_means,
                         matrix(means[h, ], nrow(group_means), ncol(group_means),
                                byrow = TRUE),
                         matrix(sqrt(variances[h, ]), nrow(group_means),
                                ncol(group_means), byrow = TRUE), log = TRUE))
  }, numeric(nrow(group_means))), nrow = nrow(group_means))
}

#' Means and variances of the group means given class weights (M-step)
#' @noRd
.cross_level_between_update <- function(group_means, posterior, between_variance,
                                        min_variance) {
  counts <- colSums(posterior)
  means <- crossprod(posterior, group_means) / counts
  spread <- matrix(vapply(seq_len(ncol(posterior)), function(h) {
    colSums(posterior[, h] * sweep(group_means, 2L, means[h, ])^2)
  }, numeric(ncol(group_means))), nrow = ncol(posterior), byrow = TRUE)
  variances <- if (identical(between_variance, "equal")) {
    matrix(colSums(spread) / nrow(group_means), ncol(posterior),
           ncol(group_means), byrow = TRUE)
  } else spread / counts
  list(means = unname(means), variances = unname(pmax(variances, min_variance)))
}

#' Starting values for the full cross-level model from one partition of the
#' groups (by their means), so both parts share the group-class labels
#' @noRd
.cross_level_initialize <- function(x, group_index, group_means, n_profiles,
                                    n_group_classes, variance_model,
                                    between_variance, min_variance, start_index) {
  profile <- .multilpa_initialize(x, group_index, n_profiles, n_group_classes,
                                  variance_model, min_variance, start_index,
                                  hierarchical = start_index == 1L)
  spread <- apply(group_means, 2L, stats::sd)
  spread[!is.finite(spread) | spread <= 0] <- 1
  scaled <- sweep(group_means, 2L, spread, "/")
  hard <- if (n_group_classes == 1L) {
    rep(1L, nrow(group_means))
  } else if (start_index == 1L) {
    .multilpa_ward_assignments(scaled, n_group_classes)
  } else {
    centres <- scaled[sample.int(nrow(scaled), n_group_classes), , drop = FALSE]
    distance <- vapply(seq_len(n_group_classes), function(h) {
      colSums((t(scaled) - centres[h, ])^2)
    }, numeric(nrow(scaled)))
    max.col(-matrix(distance, nrow(scaled)), ties.method = "first")
  }
  group_posterior <- 0.9 * outer(hard, seq_len(n_group_classes), "==") +
    0.1 / n_group_classes
  between <- .cross_level_between_update(group_means, group_posterior,
                                         between_variance, min_variance)
  # Profile prevalences per group class: the members' profile posteriors under
  # the starting profiles, averaged within each class of the partition.
  log_density <- vapply(seq_len(n_profiles), function(k) {
    rowSums(stats::dnorm(x, matrix(profile$means[k, ], nrow(x), ncol(x), byrow = TRUE),
                         matrix(sqrt(profile$variances[k, ]), nrow(x), ncol(x),
                                byrow = TRUE), log = TRUE))
  }, numeric(nrow(x)))
  log_density <- matrix(log_density, nrow(x))
  member <- exp(log_density - .multilpa_row_max(log_density))
  member <- member / rowSums(member)
  weight <- group_posterior[group_index, , drop = FALSE]
  composition <- crossprod(weight, member) + 0.1 / n_profiles
  profile$profile_probabilities <- composition / rowSums(composition)
  profile$group_probabilities <- colSums(group_posterior) / nrow(group_means)
  c(profile, list(between_means = between$means,
                  between_variances = between$variances))
}

#' EM for the full cross-level model from one start
#' @noRd
.cross_level_em <- function(x, group_index, group_means, parameters,
                            variance_model, between_variance, min_variance,
                            max_iter, tol, sampling_weights = NULL) {
  evaluate <- function(point) {
    .multilpa_expectation(x, group_index, c(point, list(
      group_log_density = .cross_level_group_density(
        group_means, point$between_means, point$between_variances))),
      weights = sampling_weights)
  }
  expectation <- evaluate(parameters)
  history <- expectation$log_likelihood
  converged <- FALSE
  iteration <- 0L
  # Sequential by nature: each step starts from the last point.
  while (iteration < max_iter && !converged) {
    updated <- .multilpa_maximization(x, expectation, variance_model, min_variance)
    between <- .cross_level_between_update(group_means, expectation$group_posteriors,
                                           between_variance, min_variance)
    updated$between_means <- between$means
    updated$between_variances <- between$variances
    updated_expectation <- evaluate(updated)
    gain <- updated_expectation$log_likelihood - expectation$log_likelihood
    if (gain < -1e-10 * (1 + abs(expectation$log_likelihood))) {
      stop("EM likelihood decreased beyond numerical roundoff.")
    }
    converged <- abs(gain) <= tol * (1 + abs(expectation$log_likelihood))
    iteration <- iteration + 1L
    history <- c(history, updated_expectation$log_likelihood)
    parameters <- updated
    expectation <- updated_expectation
  }
  list(parameters = parameters, expectation = expectation, converged = converged,
       iterations = iteration, history = history)
}

#' Free parameters of a cross-level model
#' @noRd
.cross_level_count <- function(family, n_profiles, n_group_classes, d,
                               variance_model, between_variance) {
  rows <- function(mode, n) if (identical(mode, "equal")) 1L else n
  prevalences <- if (identical(family, "full_cross_level")) {
    n_group_classes * (n_profiles - 1L)
  } else n_profiles - 1L
  as.integer(n_profiles * d + rows(variance_model, n_profiles) * d + prevalences +
               n_group_classes * d + rows(between_variance, n_group_classes) * d +
               n_group_classes - 1L)
}

#' Fit a cross-level family
#'
#' Called by `multilpa(family = "restricted_cross_level" |
#' "full_cross_level")`.
#' @return An object of class `multilpa_cross_level`.
#' @noRd
.cross_level_fit <- function(data, vars, id, n_profiles, n_group_classes, family,
                             variance_model, between_variance, n_starts,
                             max_iter, tol, min_variance, seed, call,
                             weights = NULL) {
  counts <- list(n_profiles = n_profiles, n_group_classes = n_group_classes,
                 n_starts = n_starts)
  invisible(lapply(names(counts), function(field) {
    value <- counts[[field]]
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
        value < 1 || value != floor(value)) {
      stop(errorCondition(sprintf("%s must be a positive integer.", field),
                          class = "latents_bad_argument", call = NULL))
    }
  }))
  if (!is.numeric(max_iter) || length(max_iter) != 1L || !is.finite(max_iter) ||
      max_iter < 0 || max_iter != floor(max_iter) ||
      !is.numeric(tol) || length(tol) != 1L || !is.finite(tol) || tol <= 0 ||
      !is.numeric(min_variance) || length(min_variance) != 1L ||
      !is.finite(min_variance) || min_variance <= 0) {
    stop(errorCondition("max_iter, tol and min_variance must be valid numbers.",
                        class = "latents_bad_argument", call = NULL))
  }
  .multilpa_check_seed(seed)
  n_profiles <- as.integer(n_profiles)
  n_group_classes <- as.integer(n_group_classes)
  stats <- .additive_prepare(data, vars, id, weights)
  sampling_weights <- stats$sampling_weights
  if (n_group_classes > stats$n_groups) {
    stop(errorCondition(sprintf("%d group classes cannot be estimated from %d groups.",
                                n_group_classes, stats$n_groups),
                        class = "latents_unidentified", call = NULL))
  }
  x <- matrix(as.numeric(as.matrix(data[, vars, drop = FALSE])), nrow(data),
              length(vars), dimnames = list(NULL, vars))
  group_index <- stats$index
  group_means <- unname(stats$averages)
  if (!is.null(seed)) {
    had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    old_seed <- if (had_seed) get(".Random.seed", envir = .GlobalEnv) else NULL
    on.exit({
      if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)  # nolint: object_name_linter. R's name for the RNG state.
      else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    }, add = TRUE)
    set.seed(seed)
  }
  fitted <- if (identical(family, "full_cross_level")) {
    attempts <- lapply(seq_len(n_starts), function(start_index) {
      tryCatch({
        start <- .cross_level_initialize(x, group_index, group_means, n_profiles,
                                         n_group_classes, variance_model,
                                         between_variance, min_variance,
                                         start_index)
        .cross_level_em(x, group_index, group_means, start, variance_model,
                        between_variance, min_variance, max_iter, tol,
                        sampling_weights)
      }, error = function(error) list(error = conditionMessage(error)))
    })
    best <- .cross_level_select(attempts, n_starts)
    # Weighted posteriors are counts; each unit is reported by its own.
    if (!is.null(sampling_weights)) {
      reported <- .multilpa_expectation(x, group_index, c(best$parameters, list(
        group_log_density = .cross_level_group_density(
          group_means, best$parameters$between_means,
          best$parameters$between_variances))))
      reported$log_likelihood <- best$expectation$log_likelihood
      best$expectation <- reported
    }
    best
  } else {
    .cross_level_restricted(x, group_index, group_means, n_profiles,
                            n_group_classes, variance_model, between_variance,
                            n_starts, max_iter, tol, min_variance,
                            sampling_weights)
  }
  result <- .cross_level_result(fitted, family, stats, x, group_means, vars, id,
                                n_profiles, n_group_classes, variance_model,
                                between_variance, min_variance, call)
  if (!is.null(sampling_weights)) {
    result$weights <- weights
    result$sampling_weights <- stats::setNames(sampling_weights, stats$ids)
  }
  result
}

#' Keep the best start of the full cross-level EM, with every start's record
#' @noRd
.cross_level_select <- function(attempts, n_starts) {
  valid <- vapply(attempts, function(a) is.null(a$error), logical(1))
  if (!any(valid)) {
    stop(errorCondition(sprintf("All %d starts failed: %s", n_starts,
      paste(unique(vapply(attempts, `[[`, character(1), "error")), collapse = "; ")),
      class = "latents_all_starts_failed", call = NULL))
  }
  scores <- vapply(attempts, function(a) {
    if (is.null(a$error)) a$expectation$log_likelihood else -Inf
  }, numeric(1))
  converged <- vapply(attempts, function(a) isTRUE(a$converged), logical(1))
  best_start <- .multilpa_select_start(scores, converged)
  best <- attempts[[best_start]]
  best$starts <- data.frame(
    start = seq_len(n_starts), log_likelihood = ifelse(valid, scores, NA_real_),
    converged = converged,
    error = vapply(attempts, function(a) a$error %||% NA_character_, character(1)),
    selected = seq_len(n_starts) == best_start)
  best
}

#' The restricted cross-level model through its factorization
#'
#' With profile prevalences equal across group classes the likelihood is the
#' product of a profile model of every rating and a profile model of the group
#' means, so each is fitted on its own by the package's single-level EM.
#' @noRd
.cross_level_restricted <- function(x, group_index, group_means, n_profiles,
                                    n_group_classes, variance_model,
                                    between_variance, n_starts, max_iter, tol,
                                    min_variance, sampling_weights = NULL) {
  # Under sampling weights the factorization still holds: a group's weight
  # multiplies its members' ratings and its own mean, so the rating part is
  # a row-weighted fit and the group part a group-weighted one.
  single <- function(values, k, variance, unit_weights) {
    frame <- as.data.frame(values)
    indicators <- names(frame)
    if (!is.null(unit_weights)) frame$.sampling_weight <- unit_weights
    withCallingHandlers(
      multilpa(frame, indicators, NULL, n_profiles = k, variance_model = variance,
               n_starts = n_starts, max_iter = max_iter, tol = tol,
               min_variance = min_variance, acceleration = "none",
               weights = if (is.null(unit_weights)) NULL else ".sampling_weight"),
      latents_single_level = function(notice) invokeRestart("muffleMessage"))
  }
  row_weights <- if (is.null(sampling_weights)) NULL else sampling_weights[group_index]
  individual <- single(x, n_profiles, variance_model, row_weights)
  groups <- single(group_means, n_group_classes, between_variance, sampling_weights)
  # Each part's pseudo likelihood is scaled to its own unit count; put the
  # rating part back on the groups' weight scale before adding the two.
  individual_scale <- if (is.null(row_weights)) 1 else sum(row_weights) / length(row_weights)
  group_posterior <- unname(groups$subject_posteriors)
  subject <- unname(individual$subject_posteriors)
  parameters <- list(
    means = unname(individual$means), variances = unname(individual$variances),
    profile_probabilities = matrix(
      if (is.null(row_weights)) colMeans(subject) else
        colSums(subject * row_weights) / sum(row_weights),
      n_group_classes, n_profiles,
                                   byrow = TRUE),
    group_probabilities = unname(if (is.null(sampling_weights)) colMeans(group_posterior) else
      colSums(group_posterior * sampling_weights) / sum(sampling_weights)),
    between_means = unname(groups$means),
    between_variances = unname(groups$variances))
  list(parameters = parameters,
       expectation = list(log_likelihood = individual_scale *
                            individual$log_likelihood + groups$log_likelihood,
                          group_posteriors = group_posterior,
                          subject_posteriors = subject),
       converged = isTRUE(individual$converged) && isTRUE(groups$converged),
       iterations = NA_integer_, history = NA_real_,
       starts = data.frame(part = c("individual profiles", "group classes"),
                           log_likelihood = c(individual$log_likelihood,
                                              groups$log_likelihood),
                           converged = c(individual$converged, groups$converged),
                           n_best_replicated = c(individual$n_best_replicated,
                                                 groups$n_best_replicated)))
}

#' Assemble the cross-level result object
#' @noRd
.cross_level_result <- function(fitted, family, stats, x, group_means, vars, id,
                                n_profiles, n_group_classes, variance_model,
                                between_variance, min_variance, call) {
  profile_names <- paste0("profile_", seq_len(n_profiles))
  class_names <- paste0("group_class_", seq_len(n_group_classes))
  p <- fitted$parameters
  dimnames(p$means) <- dimnames(p$variances) <- list(profile_names, vars)
  dimnames(p$between_means) <- dimnames(p$between_variances) <- list(class_names, vars)
  group_posteriors <- matrix(fitted$expectation$group_posteriors, stats$n_groups,
                             dimnames = list(stats$ids, class_names))
  subject_posteriors <- matrix(fitted$expectation$subject_posteriors, nrow(x),
                               dimnames = list(NULL, profile_names))
  # Profile composition of each group class: model-implied for the full
  # family; for the restricted family, where prevalences do not depend on the
  # class, the posterior-weighted composition of its members (descriptive).
  composition <- if (identical(family, "full_cross_level")) p$profile_probabilities else {
    weight <- group_posteriors[stats$index, , drop = FALSE]
    crossprod(weight, subject_posteriors) / colSums(weight)
  }
  dimnames(composition) <- list(class_names, profile_names)
  n_parameters <- .cross_level_count(family, n_profiles, n_group_classes,
                                     length(vars), variance_model, between_variance)
  log_likelihood <- fitted$expectation$log_likelihood
  boundary <- any(p$variances <= min_variance * (1 + 1e-8)) ||
    any(p$between_variances <= min_variance * (1 + 1e-8))
  result <- list(
    call = call, family = family, vars = vars, id = id,
    variance_model = variance_model, between_variance = between_variance,
    profile_means = p$means, profile_variances = p$variances,
    composition = composition,
    group_means = p$between_means, group_mean_variances = p$between_variances,
    group_probabilities = stats::setNames(p$group_probabilities, class_names),
    group_posteriors = group_posteriors, subject_posteriors = subject_posteriors,
    log_likelihood = log_likelihood, n_parameters = n_parameters,
    aic = -2 * log_likelihood + 2 * n_parameters,
    bic = -2 * log_likelihood + n_parameters * log(stats$n_groups),
    bic_individual = -2 * log_likelihood + n_parameters * log(nrow(x)),
    n_groups = stats$n_groups, n_obs = nrow(x),
    group_sizes = stats::setNames(stats$sizes, stats$ids),
    converged = isTRUE(fitted$converged), iterations = fitted$iterations,
    log_likelihood_history = fitted$history, boundary = boundary,
    starts = fitted$starts, sufficient_statistics = stats,
    min_variance = min_variance)
  class(result) <- "multilpa_cross_level"
  if (!result$converged) {
    warning(warningCondition("The cross-level fit did not converge; increase max_iter.",
                             class = "latents_unconverged", call = NULL))
  }
  if (boundary) {
    warning(warningCondition("A variance reached min_variance; this is a bound-active fit.",
                             class = "latents_boundary", call = NULL))
  }
  result
}

.cross_level_tables <- function() {
  c("profiles", "group_classes", "composition", "groups", "assignments", "fit",
    "starts")
}

#' Tables of a cross-level fit
#'
#' Tidy tables of a fit from `multilpa(family = "restricted_cross_level")` or
#' `multilpa(family = "full_cross_level")`.
#'
#' @param x A `multilpa_cross_level` fit.
#' @param what Which table; `"all"` returns a named list of every table.
#' @param ... Unused.
#' @section Tables:
#' \describe{
#'   \item{`profiles`}{One row per individual profile and indicator: `level =
#'     "individual"`, `profile`, `indicator`, `mean`, `variance`; then one row
#'     per group class and indicator for the group means: `level = "group"`,
#'     the class in `profile`, and the class's `mean` and `variance` of the
#'     group means.}
#'   \item{`group_classes`}{One row per group class: `weight`, `count`
#'     (summed posterior) and `n_assigned`.}
#'   \item{`composition`}{One row per group class and profile: the share of
#'     the class's members in each profile. Model-implied for the full
#'     family; for the restricted family a posterior-weighted description,
#'     since its prevalences do not depend on the class.}
#'   \item{`groups`}{One row per group: size, modal group class, its
#'     posterior and every class probability.}
#'   \item{`assignments`}{One row per data row: its modal `profile` and
#'     posterior, and its group's modal `group_class` (inherited).}
#'   \item{`fit`}{One row: family, log likelihood, `n_parameters`, `aic`,
#'     `bic` (by groups), `bic_individual` (by rows), sizes, convergence and
#'     `boundary`. The likelihood is the working likelihood of the manifest
#'     specification (the group means are computed from the same ratings), so
#'     it is comparable only between cross-level fits of the same data.}
#'   \item{`starts`}{The start records (full family) or the two parts'
#'     likelihoods (restricted family).}
#' }
#' @return A base `data.frame`, or a named list of them for `what = "all"`.
#' @examples
#' set.seed(1)
#' ratings <- data.frame(
#'   team = rep(seq_len(40), each = 6),
#'   climate = rep(rnorm(40, rep(c(-1, 1), each = 20), 0.5), each = 6) +
#'     rnorm(240, sample(c(-1, 1), 240, replace = TRUE), 0.6))
#' fit <- multilpa(ratings, "climate", "team", n_profiles = 2,
#'                 n_group_classes = 2, family = "full_cross_level",
#'                 n_starts = 3, seed = 1)
#' get_results(fit)
#' get_results(fit, "composition")
#' @export
get_results.multilpa_cross_level <- function(x, what = "profiles", ...) {
  what <- match.arg(what, c(.cross_level_tables(), "all"))
  if (identical(what, "all")) {
    return(stats::setNames(lapply(.cross_level_tables(), function(name) {
      .cross_level_table(x, name)
    }), .cross_level_tables()))
  }
  .cross_level_table(x, what)
}

#' Build one table of a cross-level fit
#' @noRd
.cross_level_table <- function(x, what) {
  long <- function(means, variances, level) {
    cells <- expand.grid(row = seq_len(nrow(means)), indicator = seq_along(x$vars))
    data.frame(level = level, profile = rownames(means)[cells$row],
               indicator = x$vars[cells$indicator], mean = as.vector(means),
               variance = as.vector(variances), stringsAsFactors = FALSE)
  }
  class_names <- names(x$group_probabilities)
  assigned <- max.col(x$group_posteriors, ties.method = "first")
  switch(what,
    profiles = rbind(long(x$profile_means, x$profile_variances, "individual"),
                     long(x$group_means, x$group_mean_variances, "group")),
    group_classes = data.frame(group_class = class_names,
                               weight = unname(x$group_probabilities),
                               count = unname(colSums(x$group_posteriors)),
                               n_assigned = tabulate(assigned, length(class_names)),
                               stringsAsFactors = FALSE),
    composition = {
      cells <- expand.grid(group_class = seq_len(nrow(x$composition)),
                           profile = seq_len(ncol(x$composition)))
      data.frame(group_class = rownames(x$composition)[cells$group_class],
                 profile = colnames(x$composition)[cells$profile],
                 share = as.vector(x$composition), stringsAsFactors = FALSE)
    },
    groups = cbind(stats::setNames(data.frame(rownames(x$group_posteriors),
                                              stringsAsFactors = FALSE), x$id),
                   data.frame(n = unname(x$group_sizes),
                              group_class = class_names[assigned],
                              posterior = x$group_posteriors[cbind(seq_along(assigned),
                                                                   assigned)],
                              stats::setNames(as.data.frame(unname(x$group_posteriors)),
                                              paste0("probability_", class_names)))),
    assignments = {
      profile <- max.col(x$subject_posteriors, ties.method = "first")
      index <- x$sufficient_statistics$index
      cbind(data.frame(row = seq_along(index)),
            stats::setNames(data.frame(rownames(x$group_posteriors)[index],
                                       stringsAsFactors = FALSE), x$id),
            data.frame(profile = colnames(x$subject_posteriors)[profile],
                       posterior = x$subject_posteriors[cbind(seq_along(profile), profile)],
                       group_class = class_names[assigned[index]],
                       stringsAsFactors = FALSE))
    },
    fit = data.frame(family = x$family, n_profiles = nrow(x$profile_means),
                     n_group_classes = length(class_names),
                     log_likelihood = x$log_likelihood, n_parameters = x$n_parameters,
                     aic = x$aic, bic = x$bic, bic_individual = x$bic_individual,
                     n_groups = x$n_groups, n_obs = x$n_obs, converged = x$converged,
                     boundary = x$boundary, stringsAsFactors = FALSE),
    starts = x$starts)
}

#' @rdname get_results.multilpa_cross_level
#' @param row.names,optional Unused; part of the generic.
#' @export
as.data.frame.multilpa_cross_level <- function(x, row.names = NULL, optional = FALSE,
                                               what = "profiles", ...) {
  get_results.multilpa_cross_level(x, what = what)
}

#' @rdname get_results.multilpa_cross_level
#' @export
print.multilpa_cross_level <- function(x, ...) {
  label <- if (identical(x$family, "full_cross_level")) "Full" else "Restricted"
  cat(sprintf(paste0(
    "%s cross-level model: %d profiles, %d group classes\n",
    "%d groups, %d observations; working log-likelihood %.3f, %d parameters\n",
    "converged: %s\n"),
    label, nrow(x$profile_means), length(x$group_probabilities), x$n_groups,
    x$n_obs, x$log_likelihood, x$n_parameters, if (x$converged) "yes" else "no"))
  .latents_print_weights(x)
  cat(paste("Tables: get_results(x, what = ), e.g. \"profiles\", \"composition\",",
            "\"group_classes\", \"groups\"; summary(x).\n"))
  invisible(x)
}

#' @rdname get_results.multilpa_cross_level
#' @param object A `multilpa_cross_level` fit.
#' @export
summary.multilpa_cross_level <- function(object, ...) {
  get_results.multilpa_cross_level(object, "all")
}

#' @rdname get_results.multilpa_cross_level
#' @export
logLik.multilpa_cross_level <- function(object, ...) {
  structure(object$log_likelihood, df = object$n_parameters,
            nobs = object$n_groups, class = "logLik")
}

#' @rdname get_results.multilpa_cross_level
#' @export
nobs.multilpa_cross_level <- function(object, ...) object$n_groups

#' @rdname parameter_inference
#' @export
parameter_inference.multilpa_cross_level <- function(x, ...) {
  stop(errorCondition(paste(
    "Standard errors are not available for the cross-level families: their",
    "working likelihood uses the group means twice."),
    class = "latents_unsupported_inference", call = NULL))
}

#' Plot a cross-level fit
#'
#' `"profiles"` shows the individual profiles' means per indicator;
#' `"group_means"` the group classes' means of the group means;
#' `"composition"` the profile shares within each group class. Profiles and
#' classes are told apart by colour and point shape.
#' @param x A `multilpa_cross_level` fit.
#' @param what Which view.
#' @param main,subtitle Optional title and subtitle.
#' @param ... Unused.
#' @return A ggplot object. Raises `latents_missing_package` without ggplot2.
#' @export
plot.multilpa_cross_level <- function(x, what = c("profiles", "group_means",
                                                  "composition"),
                                      main = NULL, subtitle = NULL, ...) {
  what <- match.arg(what)
  .gg_require()
  if (identical(what, "composition")) {
    frame <- .cross_level_table(x, "composition")
    profiles <- unique(frame$profile)
    colours <- stats::setNames(rep_len(.gg_okabe_ito, length(profiles)), profiles)
    return(ggplot2::ggplot(frame, ggplot2::aes(x = group_class, y = share,
                                               fill = profile)) +
      ggplot2::geom_col(colour = "grey20", linewidth = 0.3) +
      ggplot2::geom_text(ggplot2::aes(label = profile),
                         position = ggplot2::position_stack(vjust = 0.5), size = 3) +
      ggplot2::scale_fill_manual(values = colours) +
      .gg_theme() + .gg_top_legend() +
      ggplot2::labs(x = NULL, y = "Share of members") +
      .gg_titles(main, subtitle, "Profile composition of each group class",
                 if (identical(x$family, "full_cross_level")) "Model-implied" else
                   "Posterior-weighted description", NULL))
  }
  chosen <- if (identical(what, "profiles")) "individual" else "group"
  frame <- subset(.cross_level_table(x, "profiles"), level == chosen)
  frame$indicator <- factor(frame$indicator, levels = x$vars)
  series <- unique(frame$profile)
  colours <- stats::setNames(rep_len(.gg_okabe_ito, length(series)), series)
  shapes <- stats::setNames(rep_len(.gg_shapes, length(series)), series)
  ggplot2::ggplot(frame, ggplot2::aes(x = indicator, y = mean, colour = profile,
                                      fill = profile, shape = profile,
                                      group = profile)) +
    ggplot2::geom_line(linewidth = 0.6, alpha = 0.7) +
    ggplot2::geom_point(size = 3, colour = "grey20") +
    ggplot2::scale_colour_manual(values = colours) +
    ggplot2::scale_fill_manual(values = colours) +
    ggplot2::scale_shape_manual(values = shapes) +
    .gg_theme() + .gg_top_legend() +
    ggplot2::labs(x = NULL, y = "Mean") +
    .gg_titles(main, subtitle,
               if (chosen == "individual") "Individual profiles" else
                 "Group classes (means of the group means)", NULL, NULL)
}
