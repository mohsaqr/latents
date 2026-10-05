# predict() for fitted latent profile and class models.

#' Classify new observations with a fitted model
#'
#' Evaluates the fitted mixture on new rows: the E-step of the fit, at its
#' estimates, applied to data it has not seen. Rows go through the fit's own
#' preparation -- the same categorical levels, the same centering, the same
#' missing-data handling -- so a prediction on the training data reproduces
#' the fit's own posteriors exactly.
#'
#' @section Groups:
#' A two-level model classifies a row using its group: the group class is
#' inferred from all of that group's rows, and it shifts the row's profile
#' probabilities. The groups in `newdata` (the fit's `id` column) are treated
#' as *new* groups, each classified from its own rows alone, which is the case
#' of applying a model to a new cohort. A single-level fit, or a fit given
#' `id = NULL`, treats every row as its own unit.
#'
#' @param object A fit from [multilpa()] or [multilca()] without membership
#'   covariates.
#' @param newdata `NULL` for the fitted rows, or a data frame with the fit's
#'   indicator columns (and its `id` column for a two-level fit).
#' @param type `"class"`: the modal profile and its probability, and the
#'   modal group class. `"posterior"`: the same with every profile's
#'   probability. `"density"`: each row's log density under the fitted
#'   mixture, marginal over the group classes (the density of a row from a
#'   new group), for scoring, outlier detection or held-out likelihood.
#' @param ... Unused.
#' @return A base `data.frame` with one row per row of `newdata`: `row` (its
#'   position in `newdata`), the `id` column for a two-level fit, and
#'   - for `"class"`: `profile` (`0` is the noise component of a
#'     `noise = TRUE` fit) and its `posterior`, and for two-level fits
#'     `group_class` and `group_posterior`;
#'   - for `"posterior"`: those columns and one `probability_profile_k`
#'     column per profile (and `probability_noise`);
#'   - for `"density"`: `log_density`.
#' @section Conditions:
#'   `latents_bad_data` when `newdata` lacks an indicator or the `id` column,
#'   has a missing value the fit cannot handle (`missing = "error"`), a
#'   category the fit never saw, or a non-numeric continuous indicator.
#' @examples
#' activity <- c("browse", "lectures", "forum_read")
#' training <- subset(course_engagement, student <= 60)
#' fit <- multilpa(training, activity, id = "student", n_profiles = 2,
#'                 n_group_classes = 2, n_starts = 2, seed = 1)
#' new_students <- subset(course_engagement, student > 100)
#' head(predict(fit, new_students))
#' head(predict(fit, new_students, type = "density"))
#' @importFrom stats predict
#' @export
predict.multilpa <- function(object, newdata = NULL,
                             type = c("class", "posterior", "density"), ...) {
  type <- match.arg(type)
  prepared <- if (is.null(newdata)) .multilpa_training_rows(object) else
    .multilpa_prediction_rows(object, newdata)
  parameters <- object[intersect(c("means", "variances", "covariances",
                                   "profile_probabilities",
                                   "group_probabilities",
                                   "response_probabilities", "ordinal_intercepts",
                                   "ordinal_locations", "count_means",
                                   "count_dispersion"),
                                 names(object))]
  if (isTRUE(object$noise)) {
    # The fit reports the noise share apart from the Gaussian profiles' shares
    # (together they sum to one); the E-step takes it as a last column.
    parameters$profile_probabilities <- cbind(parameters$profile_probabilities,
                                              noise = object$noise_probability)
    parameters$noise_log_density <- rep(-log(object$hypervolume),
                                        nrow(prepared$x))
  }
  if (identical(type, "density")) {
    return(data.frame(row = prepared$rows,
                      log_density = .multilpa_row_log_density(
                        prepared, parameters)))
  }
  expectation <- .multilpa_expectation(prepared$x, prepared$group_index,
                                       parameters, prepared$codes,
                                       extra = prepared$extra)
  .multilpa_prediction_table(object, prepared, expectation,
                             posteriors = identical(type, "posterior"))
}

#' The fitted rows, prepared as the fit prepared them
#' @noRd
.multilpa_training_rows <- function(object) {
  .multilpa_prediction_rows(object, .multilpa_fit_data(object))
}

#' The data frame a fit was made from, as it was handed in
#' @noRd
.multilpa_fit_data <- function(object) {
  data <- get_results(object, "data")
  if (!is.null(object$id) && !isTRUE(object$single_level) &&
      !object$id %in% names(data)) {
    data[[object$id]] <- object$group_values
  }
  data
}

#' Prepare new rows as the fit prepared its own
#'
#' @param object A `multilpa` fit.
#' @param newdata A data frame.
#' @return A list of the continuous matrix `x`, categorical `codes` (or
#'   `NULL`), `group_index`, the group labels `groups`, and `rows`.
#' @noRd
.multilpa_prediction_rows <- function(object, newdata) {
  bad_data <- function(message) {
    stop(errorCondition(message, class = "latents_bad_data", call = NULL))
  }
  if (!is.data.frame(newdata)) bad_data("`newdata` must be a data frame.")
  categorical <- object$categorical %||% character()
  continuous <- .multilpa_continuous_names(object)
  two_level <- !isTRUE(object$single_level) && !is.null(object$id)
  needed <- c(object$vars, if (two_level) object$id)
  absent <- setdiff(needed, names(newdata))
  if (length(absent) > 0L) {
    bad_data(sprintf("`newdata` lacks the column%s %s the fit uses.",
                     if (length(absent) > 1L) "s" else "",
                     paste(sprintf("`%s`", absent), collapse = ", ")))
  }
  n <- nrow(newdata)
  if (n == 0L) bad_data("`newdata` has no rows.")
  if (!all(vapply(newdata[continuous], function(value) {
    is.numeric(value) && is.null(dim(value))
  }, logical(1)))) {
    bad_data("Every continuous indicator in `newdata` must be a numeric vector.")
  }
  # A zero-column frame converts to a logical matrix; the E-step needs double.
  x <- matrix(as.numeric(as.matrix(newdata[continuous])), nrow = n,
              ncol = length(continuous), dimnames = list(NULL, continuous))
  codes <- if (length(categorical) == 0L) NULL else {
    matrix(vapply(categorical, function(indicator) {
      value <- newdata[[indicator]]
      code <- match(as.character(value), object$categorical_levels[[indicator]])
      unseen <- !is.na(value) & is.na(code)
      if (any(unseen)) {
        bad_data(sprintf(paste(
          "`%s` has the categor%s %s, which the fit never saw, so it has no",
          "estimated probability."), indicator,
          if (length(unique(value[unseen])) > 1L) "ies" else "y",
          paste(sprintf("`%s`", unique(as.character(value[unseen]))),
                collapse = ", ")))
      }
      code
    }, integer(n)), nrow = n, dimnames = list(NULL, categorical))
  }
  extra <- .latents_prepare_extra_like(object, newdata)
  if (identical(object$missing %||% "error", "error") &&
      (anyNA(x) || (!is.null(codes) && anyNA(codes)) ||
       anyNA(.latents_extra_matrix(extra)))) {
    bad_data(paste(
      "`newdata` has missing indicator values, and the fit was made with",
      "`missing = \"error\"`. Refit with `missing = \"fiml\"` to classify",
      "incomplete rows."))
  }
  if (any(is.infinite(x))) bad_data("`newdata` has non-finite indicator values.")
  groups <- if (two_level) newdata[[object$id]] else seq_len(n)
  if (two_level && anyNA(groups)) bad_data("The `id` column has missing values.")
  group_levels <- unique(groups)
  group_index <- match(groups, group_levels)
  centering <- object$centering %||% "none"
  if (identical(centering, "grand")) {
    x <- sweep(x, 2L, object$centering_offsets[1L, continuous], "-")
  } else if (identical(centering, "person")) {
    x <- .multilpa_center_indicators(x, group_index, length(group_levels),
                                     "person")$x
  }
  list(x = x, codes = codes, extra = extra, group_index = group_index,
       groups = groups, rows = seq_len(n), two_level = two_level)
}

#' Log density of each row under the fitted mixture
#'
#' Marginal over the group classes: the profile shares are the group classes'
#' shares averaged with the group-class probabilities, which is the density
#' of a row from a group not yet seen.
#' @noRd
.multilpa_row_log_density <- function(prepared, parameters) {
  shares <- drop(parameters$group_probabilities %*%
                   parameters$profile_probabilities)
  single <- parameters
  single$profile_probabilities <- matrix(shares, nrow = 1L)
  single$group_probabilities <- 1
  expectation <- .multilpa_expectation(prepared$x,
                                       seq_len(nrow(prepared$x)), single,
                                       prepared$codes, extra = prepared$extra)
  expectation$group_log_likelihood
}

#' Tidy prediction table
#' @noRd
.multilpa_prediction_table <- function(object, prepared, expectation,
                                       posteriors) {
  posterior <- expectation$subject_posteriors
  noisy <- isTRUE(object$noise)
  modal <- max.col(posterior, ties.method = "first")
  n_profiles <- object$n_profiles
  table <- data.frame(row = prepared$rows)
  if (prepared$two_level) table[[object$id]] <- prepared$groups
  table$profile <- ifelse(noisy & modal == ncol(posterior), 0L, modal)
  table$posterior <- posterior[cbind(seq_len(nrow(posterior)), modal)]
  if (prepared$two_level && object$n_group_classes > 1L) {
    group_modal <- max.col(expectation$group_posteriors, ties.method = "first")
    table$group_class <- group_modal[prepared$group_index]
    table$group_posterior <- expectation$group_posteriors[
      cbind(prepared$group_index, table$group_class)]
  }
  if (posteriors) {
    columns <- paste0("probability_profile_", seq_len(n_profiles))
    if (noisy) columns <- c(columns, "probability_noise")
    probability <- as.data.frame(posterior)
    names(probability) <- columns
    table <- cbind(table, probability)
  }
  table
}
