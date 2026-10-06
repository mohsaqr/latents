# Comparing group-class families: enumeration over families and numbers of
# group classes, and the parametric bootstrap likelihood-ratio test.

utils::globalVariables(c("criterion", "series", "value", "n_group_classes"))

#' Model code of a group-class fit: the family, with `_equal` when its
#' between-group variances are shared and the family also allows them to vary
#' @noRd
.additive_model_code <- function(family, between_variance) {
  ifelse(family != "dispersion" & between_variance == "equal",
         paste0(family, "_equal"), family)
}

#' Simulate raw ratings from a fitted group-class model
#'
#' Same groups, sizes and identifiers as the fitted data; each group's class is
#' drawn from the class weights, its intercepts from the class's between
#' distribution, and its ratings from the class's within distribution.
#' @param x A `multilpa_additive` fit.
#' @return A data frame with the fit's `id` column and indicators.
#' @noRd
.additive_simulate <- function(x) {
  stats <- x$sufficient_statistics
  n_classes <- length(x$group_probabilities)
  d <- length(x$vars)
  classes <- sample.int(n_classes, stats$n_groups, replace = TRUE,
                        prob = x$group_probabilities)
  intercepts <- matrix(stats::rnorm(stats$n_groups * d,
                                    x$means[classes, , drop = FALSE],
                                    sqrt(x$between_variances[classes, , drop = FALSE])),
                       stats$n_groups, d)
  row_class <- classes[stats$index]
  noise <- matrix(stats::rnorm(stats$n_obs * d), stats$n_obs, d) *
    sqrt(x$within_variances[row_class, , drop = FALSE])
  ratings <- intercepts[stats$index, , drop = FALSE] + noise
  colnames(ratings) <- x$vars
  cbind(stats::setNames(data.frame(stats$ids[stats$index],
                                   stringsAsFactors = FALSE), x$id),
        as.data.frame(ratings))
}

#' Is one group-class model nested in another?
#'
#' Nested when the alternative has at least as many classes and each block
#' (means, between and within variances) is at least as free: a block shared
#' by every class is a special case of a class-specific one. One class nests
#' in every family, since all families coincide there.
#' @noRd
.additive_nested <- function(null_model, alternative_model) {
  if (length(null_model$group_probabilities) == 1L) return(TRUE)
  rank <- c(equal = 0L, varying = 1L)
  blocks <- c("means", "between", "within")
  length(alternative_model$group_probabilities) >=
    length(null_model$group_probabilities) &&
    all(rank[unlist(alternative_model$structure[blocks])] >=
          rank[unlist(null_model$structure[blocks])])
}

#' Parametric bootstrap likelihood-ratio test for group-class fits
#' @noRd
.additive_bootstrap_lrt <- function(null_model, alternative_model, data, iter,
                                    n_starts, max_iter, tol, seed, call) {
  if (!inherits(null_model, "multilpa_additive") ||
      !inherits(alternative_model, "multilpa_additive")) {
    stop(errorCondition(paste(
      "bootstrap_lrt() compares two group-class fits (multilpa(family = ...))",
      "or two profile fits, not one of each."),
      class = "latents_bad_argument", call = NULL))
  }
  same_statistics <- function(left, right) {
    identical(left$ids, right$ids) && identical(left$sizes, right$sizes) &&
      isTRUE(all.equal(left$averages, right$averages)) &&
      isTRUE(all.equal(left$scatter, right$scatter))
  }
  same_data <- identical(null_model$vars, alternative_model$vars) &&
    identical(null_model$id, alternative_model$id) &&
    same_statistics(null_model$sufficient_statistics,
                     alternative_model$sufficient_statistics)
  if (!same_data) {
    stop(errorCondition("The two fits were not made from the same data and indicators.",
                        class = "latents_bad_inference_data", call = NULL))
  }
  if (!is.null(data)) {
    supplied <- .additive_prepare(data, null_model$vars, null_model$id)
    if (!same_statistics(supplied, null_model$sufficient_statistics)) {
      stop(errorCondition("`data` does not reproduce the data these fits were made from.",
                          class = "latents_bad_inference_data", call = NULL))
    }
  }
  equivalent_one_class <- length(null_model$group_probabilities) == 1L &&
    length(alternative_model$group_probabilities) == 1L
  if (!.additive_nested(null_model, alternative_model) || equivalent_one_class ||
      identical(c(length(null_model$group_probabilities), null_model$structure),
                c(length(alternative_model$group_probabilities),
                  alternative_model$structure))) {
    stop(errorCondition(paste(
      "The null model must be nested in, and differ from, the alternative: at",
      "most as many group classes, and every block (means, between and within",
      "variances) shared wherever the alternative shares it."),
      class = "latents_bad_argument", call = NULL))
  }
  if (!isTRUE(null_model$converged) || !isTRUE(alternative_model$converged)) {
    stop(errorCondition("Both fits must have converged.",
                        class = "latents_no_converge", call = NULL))
  }
  stopifnot("`iter` must be a single positive integer" =
              is.numeric(iter) && length(iter) == 1L && is.finite(iter) && iter >= 1 &&
              iter == floor(iter) && iter <= .Machine$integer.max,
            "`n_starts` must be a single positive integer" =
              is.numeric(n_starts) && length(n_starts) == 1L && is.finite(n_starts) &&
              n_starts >= 1 && n_starts == floor(n_starts) &&
              n_starts <= .Machine$integer.max,
            "`max_iter` must be a single positive integer" =
              is.numeric(max_iter) && length(max_iter) == 1L && is.finite(max_iter) &&
              max_iter >= 1 && max_iter == floor(max_iter) &&
              max_iter <= .Machine$integer.max,
            "`tol` must be finite and positive" =
              is.numeric(tol) && length(tol) == 1L && is.finite(tol) && tol > 0)
  .multilpa_check_seed(seed)
  .latents_local_seed(seed)
  # Convergence noise allowed below zero, as in the profile models.
  reversal_window <- 100 * tol * (1 + abs(null_model$log_likelihood))
  observed <- 2 * (alternative_model$log_likelihood - null_model$log_likelihood)
  if (observed < -reversal_window) {
    stop(errorCondition(sprintf(paste(
      "Alternative has lower likelihood by %.3g, beyond the %.3g that",
      "`tol = %.3g` can explain; improve its optimization first."),
      -observed, reversal_window, tol),
      class = "latents_reversed_likelihood", call = NULL))
  }
  observed <- max(0, observed)
  refit <- function(model, simulated) {
    multilpa(simulated, model$vars, model$id,
             n_group_classes = length(model$group_probabilities),
             family = model$family, between_variance = model$between_variance,
             n_starts = n_starts, max_iter = max_iter, tol = tol,
             min_variance = model$min_variance)
  }
  replicates <- .latents_blrt_replicates(iter, function(i) {
    simulated <- .additive_simulate(null_model)
    models <- lapply(list(null_model, alternative_model), refit, simulated)
    statistic <- 2 * (models[[2L]]$log_likelihood - models[[1L]]$log_likelihood)
    list(statistic = statistic,
         valid = all(vapply(models, function(m) isTRUE(m$converged), logical(1))) &&
           statistic >= -reversal_window,
         boundary = any(vapply(models, function(m) {
           any(m$between_zero) || any(m$within_floor)
         }, logical(1))),
         logits_settled = NA,
         null_replications = models[[1L]]$n_best_replicated,
         alternative_replications = models[[2L]]$n_best_replicated)
  }, muffle_warnings = TRUE)
  tally <- .latents_blrt_p_value(replicates, observed, iter)
  result <- list(statistic = observed, p_value = tally$p_value,
       monte_carlo_se = tally$monte_carlo_se,
       iter = iter, n_valid = tally$n_valid, replicates = replicates,
       fixed = character(),
       null_profiles = NA_integer_,
       null_group_classes = length(null_model$group_probabilities),
       alternative_profiles = NA_integer_,
       alternative_group_classes = length(alternative_model$group_probabilities),
       null_family = .additive_model_code(null_model$family, null_model$between_variance),
       alternative_family = .additive_model_code(alternative_model$family,
                                                 alternative_model$between_variance),
       call = call)
  class(result) <- "multilpa_bootstrap_lrt"
  result
}

#' Enumerate group-class families and numbers of group classes
#' @noRd
.additive_enumerate <- function(data, vars, id, family, n_group_classes,
                                between_variance, seed, extra, call) {
  families <- c("additive", "dispersion", "additive_dispersion")
  if (!is.character(family) || length(family) < 1L || anyNA(family) ||
      !all(family %in% families)) {
    stop(errorCondition(sprintf(
      "`family` must name group-class families from %s, or be \"profiles\" alone.",
      paste(sprintf("\"%s\"", families), collapse = ", ")),
      class = "latents_bad_argument", call = NULL))
  }
  if (!is.character(between_variance) || length(between_variance) < 1L ||
      anyNA(between_variance) ||
      !all(between_variance %in% c("varying", "equal"))) {
    stop(errorCondition("`between_variance` must be \"varying\", \"equal\" or both.",
                        class = "latents_bad_argument", call = NULL))
  }
  unknown <- setdiff(names(extra), c("n_starts", "max_iter", "tol", "min_variance",
                                     "weights"))
  if (length(unknown) > 0L) {
    stop(errorCondition(sprintf(
      "%s cannot be used when enumerating group-class families.",
      paste(sprintf("`%s`", unknown), collapse = ", ")),
      class = "latents_bad_argument", call = NULL))
  }
  grid <- unique(do.call(rbind, lapply(unique(family), function(f) {
    expand.grid(family = f, n_group_classes = unique(n_group_classes),
                between_variance = if (f == "dispersion") "equal" else
                  unique(between_variance),
                stringsAsFactors = FALSE)
  })))
  grid <- grid[order(grid$n_group_classes, match(grid$family, families),
                     grid$between_variance != "varying"), , drop = FALSE]
  runs <- lapply(seq_len(nrow(grid)), function(i) {
    warnings <- character()
    error_text <- NA_character_
    fit <- tryCatch(withCallingHandlers(do.call(multilpa, c(
      list(data = data, vars = vars, id = id,
           n_group_classes = grid$n_group_classes[i], family = grid$family[i],
           between_variance = grid$between_variance[i], seed = seed), extra)),
      warning = function(warning) {
        warnings <<- c(warnings, conditionMessage(warning))
        invokeRestart("muffleWarning")
      }), error = function(error) {
        error_text <<- conditionMessage(error)
        NULL
      })
    fit_row <- if (is.null(fit)) NULL else get_results(fit, "fit")
    row <- data.frame(
      model = .additive_model_code(grid$family[i], grid$between_variance[i]),
      family = grid$family[i], between_variance = grid$between_variance[i],
      n_group_classes = grid$n_group_classes[i],
      log_likelihood = fit_row$log_likelihood %||% NA_real_,
      n_parameters = fit_row$n_parameters %||% NA_integer_,
      aic = fit_row$aic %||% NA_real_, bic = fit_row$bic %||% NA_real_,
      bic_individual = fit_row$bic_individual %||% NA_real_,
      entropy = fit_row$entropy %||% NA_real_,
      smallest_class = fit_row$smallest_class %||% NA_real_,
      min_effective_groups = fit_row$min_effective_groups %||% NA_real_,
      converged = fit_row$converged %||% FALSE,
      boundary = fit_row$boundary %||% NA,
      n_best_replicated = fit_row$n_best_replicated %||% NA_integer_,
      warnings = paste(unique(warnings), collapse = "; "), error = error_text,
      stringsAsFactors = FALSE)
    list(fit = fit, row = row)
  })
  table <- do.call(rbind, lapply(runs, `[[`, "row"))
  rownames(table) <- NULL
  structure(list(table = table, fits = lapply(runs, `[[`, "fit"), call = call),
            class = "latents_family_enumeration")
}

#' Tables of a group-class family enumeration
#'
#' The result of `enumerate_classes(family = )` with group-class families:
#' one fitted candidate per family, between-variance restriction and number
#' of group classes.
#'
#' @param x A `latents_family_enumeration` result.
#' @param what `"candidates"` (one row per candidate) or `"best"` (the
#'   candidate with the lowest value of `criterion`).
#' @param criterion `"bic"` (penalized by the number of groups, the default),
#'   `"aic"` or `"bic_individual"`.
#' @param ... Unused.
#' @return A base `data.frame`, one row per candidate: `model` (the code
#'   [candidate_fit()] takes), `family`, `between_variance`,
#'   `n_group_classes`, `log_likelihood`, `n_parameters`, `aic`, `bic`,
#'   `bic_individual`, relative `entropy`, `smallest_class`,
#'   `min_effective_groups`, `converged`, `boundary`, `n_best_replicated`,
#'   and the candidate's `warnings` and `error`. With one group class every
#'   family is the same model, so those rows share a likelihood.
#' @examples
#' set.seed(1)
#' ratings <- data.frame(
#'   team = rep(seq_len(40), each = 6),
#'   climate = rep(rnorm(40, rep(c(-1, 1), each = 20), 0.5), each = 6) +
#'     rnorm(240))
#' grid <- enumerate_classes(ratings, "climate", "team",
#'                           family = c("additive", "dispersion"),
#'                           n_group_classes = 1:2, n_starts = 3, seed = 1)
#' get_results(grid)
#' get_results(grid, "best")
#' @export
get_results.latents_family_enumeration <- function(x, what = c("candidates", "best"),
                                                   criterion = c("bic", "aic",
                                                                 "bic_individual"),
                                                   ...) {
  what <- match.arg(what)
  criterion <- match.arg(criterion)
  table <- x$table
  if (identical(what, "candidates")) return(table)
  usable <- which(is.finite(table[[criterion]]))
  if (length(usable) == 0L) {
    stop(errorCondition("No candidate was fitted.",
                        class = "latents_failed_candidate", call = NULL))
  }
  best <- usable[which.min(table[[criterion]][usable])]
  table[best, , drop = FALSE]
}

#' @rdname get_results.latents_family_enumeration
#' @param row.names,optional Unused; part of the generic.
#' @export
as.data.frame.latents_family_enumeration <- function(x, row.names = NULL,
                                                     optional = FALSE, ...) {
  x$table
}

#' @rdname get_results.latents_family_enumeration
#' @export
print.latents_family_enumeration <- function(x, ...) {
  cat(sprintf("Group-class family enumeration: %d candidates\n", nrow(x$table)))
  print(x$table[, c("model", "n_group_classes", "log_likelihood", "n_parameters",
                    "bic", "aic", "entropy", "min_effective_groups", "converged",
                    "boundary")], digits = 4, row.names = FALSE)
  cat("Pick one with candidate_fit(x, n_group_classes = , model = ).\n")
  invisible(x)
}

#' @rdname get_results.latents_family_enumeration
#' @param object A `latents_family_enumeration` result.
#' @export
summary.latents_family_enumeration <- function(object, ...) {
  object$table
}

#' Plot information criteria across group-class families
#'
#' One line per model code, criterion against the number of group classes;
#' models are told apart by colour, point shape and line type.
#' @param x A `latents_family_enumeration` result.
#' @param criterion `"bic"`, `"aic"` or `"bic_individual"`.
#' @param main,subtitle Optional title and subtitle.
#' @param ... Unused.
#' @return A ggplot object. Raises `latents_missing_package` without ggplot2.
#' @examples
#' set.seed(1)
#' ratings <- data.frame(
#'   team = rep(seq_len(40), each = 6),
#'   climate = rep(rnorm(40, rep(c(-1, 1), each = 20), 0.5), each = 6) +
#'     rnorm(240))
#' grid <- enumerate_classes(ratings, "climate", "team",
#'                           family = c("additive", "dispersion"),
#'                           n_group_classes = 1:2, n_starts = 3, seed = 1)
#' if (requireNamespace("ggplot2", quietly = TRUE)) plot(grid)
#' @export
plot.latents_family_enumeration <- function(x, criterion = c("bic", "aic",
                                                             "bic_individual"),
                                            main = NULL, subtitle = NULL, ...) {
  criterion <- match.arg(criterion)
  .gg_require()
  frame <- x$table
  frame$value <- frame[[criterion]]
  frame <- subset(frame, is.finite(value))
  models <- unique(frame$model)
  frame$series <- factor(frame$model, levels = models)
  colours <- stats::setNames(rep_len(.gg_okabe_ito, length(models)), models)
  shapes <- stats::setNames(rep_len(.gg_shapes, length(models)), models)
  linetypes <- stats::setNames(rep_len(c("solid", "dashed", "dotdash", "dotted",
                                         "longdash"), length(models)), models)
  ggplot2::ggplot(frame, ggplot2::aes(x = n_group_classes, y = value,
                                      colour = series, fill = series,
                                      shape = series, linetype = series)) +
    ggplot2::geom_line(linewidth = 0.6) +
    ggplot2::geom_point(size = 2.8, colour = "grey20") +
    ggplot2::scale_colour_manual(values = colours) +
    ggplot2::scale_fill_manual(values = colours) +
    ggplot2::scale_shape_manual(values = shapes) +
    ggplot2::scale_linetype_manual(values = linetypes) +
    ggplot2::scale_x_continuous(breaks = sort(unique(frame$n_group_classes))) +
    .gg_theme() + .gg_top_legend() +
    ggplot2::labs(x = "Group classes", y = toupper(criterion)) +
    .gg_titles(main, subtitle, "Group-class families compared",
               "Lower is better", NULL)
}
