#' Classification entropy sum
#' @param probabilities Posterior probability matrix, rows summing to one.
#' @return Sum of `-p log p` over every cell, the classification entropy.
#' @noRd
.multilpa_entropy_sum <- function(probabilities) {
  stopifnot("`probabilities` must be a numeric matrix" =
              is.matrix(probabilities) && is.numeric(probabilities))
  positive <- probabilities > 0
  -sum(probabilities[positive] * log(probabilities[positive]))
}

#' Relative classification entropy
#' @param probabilities Posterior probability matrix, rows summing to one.
#' @return Entropy on the zero-to-one scale, or `NA_real_` for a single class.
#' @noRd
.multilpa_relative_entropy <- function(probabilities) {
  stopifnot("`probabilities` must be a numeric matrix" =
              is.matrix(probabilities) && is.numeric(probabilities))
  if (ncol(probabilities) == 1L) return(NA_real_)
  1 - .multilpa_entropy_sum(probabilities) /
    (nrow(probabilities) * log(ncol(probabilities)))
}

#' The formula behind each criterion, keyed by its name
#' @param criterion Character vector of criterion names.
#' @return Character vector of formulas, the same length as `criterion`.
#' @noRd
.multilpa_criterion_definitions <- function(criterion) {
  definitions <- c(deviance = "-2L",
                   aic = "-2L + 2q",
                   kic = "-2L + 3(q + 1)",
                   bic = "-2L + q log(n)",
                   sabic = "-2L + q log((n + 2) / 24)",
                   caic = "-2L + q (log(n) + 1)",
                   awe = "-2(L - EN) + 2q (1.5 + log(n))",
                   icl = "-2L + q log(n) + 2 EN",
                   clc = "-2L + 2 EN")
  stopifnot("Every criterion must have a definition" =
              all(criterion %in% names(definitions)))
  unname(definitions[criterion])
}

#' Information criteria for a fitted multilevel latent profile model
#'
#' Returns every information criterion the package computes, in one tidy table,
#' with the sample-size convention stated on each row. Multilevel mixtures admit
#' two defensible sample sizes: the number of independent groups and the number
#' of individuals. Criteria that depend on a sample size are therefore reported
#' once per convention rather than silently fixing one.
#'
#'   The verb is named for what it returns. Classification sharpness and class
#'   separation are reported by [classification_table()] and [entropy_table()]
#'   rather than folded into this table, so `information_criteria()` stays a
#'   pure likelihood-penalty summary.
#'
#' @param object A fitted `multilpa` model.
#' @param definitions `FALSE`, the default, returns the numbers alone. `TRUE`
#'   appends a `definition` column carrying each criterion's formula, which is
#'   the same text as the Details section below.
#' @seealso [classification_table()] and [entropy_table()] for the
#'   classification diagnostics that accompany these criteria. [logLik()] for
#'   the maximized log likelihood itself.
#' @return A base `data.frame` of class `data.frame`, one row per criterion and
#'   sample-size convention, with the columns
#'   \describe{
#'     \item{`criterion`}{character: `"deviance"`, `"aic"`, `"kic"`, `"bic"`,
#'       `"sabic"`, `"caic"`, `"awe"`, `"icl"` or `"clc"`.}
#'     \item{`convention`}{character: `"groups"` or `"individuals"`, and
#'       `NA_character_` for a criterion that uses no sample size and no
#'       level-specific entropy.}
#'     \item{`n`}{integer: the sample size that convention supplies, and
#'       `NA_integer_` where no sample size enters.}
#'     \item{`value`}{numeric: the criterion. **Lower is better on every row.**}
#'     \item{`penalty`}{numeric: what the criterion adds to the deviance, so
#'       `value == deviance + penalty` holds on every row, and `0` on the
#'       `"deviance"` row itself.}
#'     \item{`definition`}{character, present only when
#'       `definitions = TRUE`: the criterion's formula.}
#'   }
#'   The table reports the deviance `-2L` rather than the log likelihood `L`, so
#'   that the `value` column points one way throughout; `logLik(object)` returns
#'   the log likelihood. `deviance`, `aic` and `kic` do not depend on a sample
#'   size and carry `convention = NA_character_` with `n = NA_integer_`. `clc`
#'   does not depend on one either, but its convention selects which level's
#'   classification uncertainty it penalizes, so it is reported once per
#'   convention. `value` is `NA_real_` where the entropy a criterion needs is
#'   undefined, which is the group level of a continuous random-intercept fit.
#' @details Let `q` be the number of free parameters, `n` the chosen sample
#'   size, and `EN` the classification entropy of the level matching that
#'   convention. The criteria are `deviance = -2L`, `aic = -2L + 2q`,
#'   `bic = -2L + q log(n)`,
#'   `sabic = -2L + q log((n + 2) / 24)`, `caic = -2L + q (log(n) + 1)`,
#'   `awe = -2(L - EN) + 2q (1.5 + log(n))`, `icl = -2L + q log(n) + 2 EN`,
#'   `kic = -2L + 3(q + 1)`, and `clc = -2L + 2 EN`. Note that `clc` uses the
#'   entropy sum, as `icl` and `awe` here do; `tidyLPA` reports a `CLC` built
#'   from relative entropy instead, which is bounded by one and so penalizes
#'   almost nothing, and the two numbers are not comparable.
#'   The entropy-based criteria use group-level posteriors under the `"groups"`
#'   convention and individual-level posteriors under the `"individuals"`
#'   convention. This is a stated per-level choice, not a unique multilevel
#'   definition; use one convention consistently across compared candidates.
#'   Individual sample sizes exclude rows with no observed indicators. For a
#'   continuous random-intercept fit, group classification entropy is undefined,
#'   so group-level `awe` and `icl` are `NA`.
#' @references Schwarz, G. (1978). Estimating the dimension of a model.
#'   Annals of Statistics, 6, 461--464. Sclove, S. L. (1987). Application of
#'   model-selection criteria to some problems in multivariate analysis.
#'   Psychometrika, 52, 333--343. Bozdogan, H. (1987). Model selection and
#'   Akaike's information criterion. Psychometrika, 52, 345--370.
#'   Banfield, J. D., & Raftery, A. E. (1993). Model-based Gaussian and
#'   non-Gaussian clustering. Biometrics, 49, 803--821. Biernacki, C., Celeux,
#'   G., & Govaert, G. (2000). Assessing a mixture model for clustering with the
#'   integrated completed likelihood. IEEE Transactions on Pattern Analysis and
#'   Machine Intelligence, 22, 719--725. Biernacki, C., & Govaert, G. (1997).
#'   Using the classification likelihood to choose the number of clusters.
#'   Computing Science and Statistics, 29, 451--457. Cavanaugh, J. E. (1999).
#'   A large-sample model selection criterion based on Kullback's symmetric
#'   divergence. Statistics and Probability Letters, 42, 333--343.
#' @examples
#' set.seed(7)
#' example_data <- data.frame(
#'   school = rep(seq_len(12), each = 10),
#'   score_a = rnorm(120), score_b = rnorm(120)
#' )
#' fit <- multilpa(example_data, c("score_a", "score_b"), "school",
#'                   n_profiles = 2, n_group_classes = 1, n_starts = 2, seed = 1)
#' information_criteria(fit)
#' information_criteria(fit, definitions = TRUE)
#' @export
information_criteria <- function(object, definitions = FALSE) {
  stopifnot("`object` must be a fitted model of this package" = .multilpa_any_fit(object),
            "`definitions` must be TRUE or FALSE" =
              isTRUE(definitions) || isFALSE(definitions))
  q <- object$n_parameters
  deviance <- -2 * object$log_likelihood
  conventions <- data.frame(
    convention = c("groups", "individuals"),
    n = c(object$n_groups, object$n_informative %||% object$n_observations))
  entropies <- c(groups = if (is.null(object$group_posteriors)) NA_real_ else
                   .multilpa_entropy_sum(object$group_posteriors),
                 individuals = .multilpa_entropy_sum(object$subject_posteriors))
  scale_free <- data.frame(
    criterion = c("deviance", "aic", "kic"),
    convention = NA_character_, n = NA_integer_,
    penalty = c(0, 2 * q, 3 * (q + 1)))
  scaled <- do.call(rbind, lapply(seq_len(nrow(conventions)), function(row) {
    convention <- conventions$convention[row]
    n <- conventions$n[row]
    entropy_sum <- entropies[[convention]]
    penalties <- c(bic = q * log(n),
                   sabic = q * log((n + 2) / 24),
                   caic = q * (log(n) + 1),
                   awe = 2 * q * (1.5 + log(n)) + 2 * entropy_sum,
                   icl = q * log(n) + 2 * entropy_sum,
                   # The only criterion here that does not use the sample size.
                   # Its convention still matters, because it decides which
                   # level's classification uncertainty is being penalized.
                   clc = 2 * entropy_sum)
    data.frame(criterion = names(penalties), convention = convention, n = n,
               penalty = unname(penalties))
  }))
  result <- rbind(scale_free, scaled)
  # Every row is the same deviance plus its own penalty, so the `value` column
  # means one thing throughout: smaller is a better-supported model.
  result$value <- deviance + result$penalty
  result <- result[c("criterion", "convention", "n", "value", "penalty")]
  if (isTRUE(definitions)) {
    result$definition <- .multilpa_criterion_definitions(result$criterion)
  }
  row.names(result) <- NULL
  result
}

#' The posterior matrices of the levels a diagnostic was asked for
#' @param object A fitted model of this package.
#' @param level One of `"individuals"`, `"groups"` or `"both"`.
#' @return A named list of posterior matrices, in reporting order.
#' @noRd
.multilpa_posterior_levels <- function(object, level) {
  levels_wanted <- if (identical(level, "both")) c("individuals", "groups") else level
  if (is.null(object$group_posteriors)) {
    if (identical(level, "groups")) {
      stop(errorCondition("This model has no discrete group classes.",
                          class = "multilpa_no_group_classes", call = NULL))
    }
    levels_wanted <- "individuals"
  }
  posteriors <- list(individuals = object$subject_posteriors,
                     groups = object$group_posteriors)
  posteriors[levels_wanted]
}

#' Classification quality for a fitted multilevel latent profile model
#'
#' Summarizes how sharply the posterior probabilities separate classes, at the
#' individual level, the group level, or both. Reports the diagnostics normally
#' expected alongside a mixture solution: modal counts, model-estimated class
#' sizes, average posterior probability in the assigned class, and the odds of
#' correct classification.
#'
#' @param object A fitted `multilpa` model.
#' @param level `"individuals"` for latent profiles, `"groups"` for latent group
#'   classes, or `"both"` to stack them in one table.
#' @param detail Removed. It used to change the columns this verb returns, which
#'   made one verb answer with two incompatible shapes. The cross-tabulation it
#'   produced is now the separate verb [average_posteriors()]. Supplying it
#'   raises an error of class `multilpa_removed_argument`.
#' @return A base `data.frame`, one row per level and class, with the columns
#'   \describe{
#'     \item{`level`}{character: `"individuals"` or `"groups"`.}
#'     \item{`class`}{integer: the class index within that level.}
#'     \item{`n_modal`}{integer: units whose modal class is this one.}
#'     \item{`proportion_modal`}{numeric: `n_modal` over the units at that level.}
#'     \item{`estimated_n`}{numeric: the model-estimated class size, the column
#'       sum of the posteriors.}
#'     \item{`estimated_proportion`}{numeric: `estimated_n` over the units.}
#'     \item{`average_posterior`}{numeric: mean posterior probability of this
#'       class among the units assigned to it, so higher is better.}
#'     \item{`odds_correct_classification`}{numeric: see Details, so higher is
#'       better.}
#'   }
#'   `average_posterior` is `NA_real_` for a class with no modal members, where
#'   the mean is taken over nothing and is undefined rather than zero.
#'   `odds_correct_classification` is `NA_real_` wherever it is undefined, that
#'   is whenever `average_posterior` is missing, or it or
#'   `estimated_proportion` is zero or one; a single-class solution therefore
#'   always reports `NA_real_` odds.
#' @details The odds of correct classification for class `k` is
#'   `(p / (1 - p)) / (r / (1 - r))`, where `p` is the average posterior in the
#'   assigned class and `r` is the model-estimated class proportion. Values near
#'   one indicate classification no better than the class proportion alone;
#'   values of five or more are conventionally read as adequate separation.
#'   Modal assignment discards classification uncertainty, so `n_modal` and
#'   `estimated_n` differ whenever entropy is below one.
#'   Continuous random-intercept fits have individual profiles but no discrete
#'   group classes; `"both"` returns individuals only and `"groups"` raises an
#'   error of class `multilpa_no_group_classes`.
#' @seealso [average_posteriors()] for the full average-posterior
#'   cross-tabulation these diagonals come from, [classification_errors()] for
#'   the complementary table conditioned on the true class, and
#'   [entropy_table()] for entropy at each level.
#' @references Nagin, D. S. (2005). Group-Based Modeling of Development.
#'   Harvard University Press.
#' @examples
#' set.seed(7)
#' example_data <- data.frame(
#'   school = rep(seq_len(12), each = 10),
#'   score_a = rnorm(120), score_b = rnorm(120)
#' )
#' fit <- multilpa(example_data, c("score_a", "score_b"), "school",
#'                   n_profiles = 2, n_group_classes = 1, n_starts = 2, seed = 1)
#' classification_table(fit)
#' @export
classification_table <- function(object, level = c("individuals", "groups", "both"),
                                 detail) {
  stopifnot("`object` must be a fitted model of this package" = .multilpa_any_fit(object))
  if (!missing(detail)) {
    stop(errorCondition(paste(
      "`detail` was removed from classification_table().",
      "The cross-tabulation it returned is now average_posteriors()."),
      class = "multilpa_removed_argument", call = NULL))
  }
  level <- match.arg(level)
  posteriors <- .multilpa_posterior_levels(object, level)
  result <- do.call(rbind, lapply(names(posteriors), function(which_level) {
    .multilpa_classification_summary(posteriors[[which_level]], which_level)
  }))
  row.names(result) <- NULL
  result
}

#' Average posterior probabilities by modal class
#'
#' The cross-tabulation mixture software reports beside a solution: for the
#' units assigned to each class, their mean posterior probability of belonging
#' to every class. The diagonal is the `average_posterior` column of
#' [classification_table()]; the off-diagonal entries say which classes a unit
#' is confused with, and each assigned class's row is a probability
#' distribution summing to one.
#'
#' This table conditions on the *assigned* class. [classification_errors()]
#' conditions on the *true* class instead, and the two are different numbers,
#' not two spellings of one table.
#'
#' @param object A fitted `multilpa` model.
#' @param level `"individuals"` for latent profiles, `"groups"` for latent group
#'   classes, or `"both"` to stack them in one table.
#' @return A base `data.frame`, one row per level, assigned class and class,
#'   with the columns
#'   \describe{
#'     \item{`level`}{character: `"individuals"` or `"groups"`.}
#'     \item{`assigned_class`}{integer: the modal class the units were assigned
#'       to, which is what the row conditions on.}
#'     \item{`class`}{integer: the class whose posterior probability is averaged.}
#'     \item{`n_assigned`}{integer: how many units carry that `assigned_class`,
#'       repeated across the row's classes.}
#'     \item{`average_posterior`}{numeric: the mean posterior probability of
#'       `class` among those units, `NA_real_` when `n_assigned` is zero and the
#'       mean is therefore undefined.}
#'   }
#' @seealso [classification_table()] for the one-row-per-class summary and
#'   [classification_errors()] for the same square conditioned the other way.
#' @examples
#' set.seed(7)
#' example_data <- data.frame(
#'   school = rep(seq_len(12), each = 10),
#'   score_a = rnorm(120), score_b = rnorm(120)
#' )
#' fit <- multilpa(example_data, c("score_a", "score_b"), "school",
#'                   n_profiles = 2, n_group_classes = 1, n_starts = 2, seed = 1)
#' average_posteriors(fit)
#' @export
average_posteriors <- function(object, level = c("individuals", "groups", "both")) {
  stopifnot("`object` must be a fitted model of this package" = .multilpa_any_fit(object))
  level <- match.arg(level)
  posteriors <- .multilpa_posterior_levels(object, level)
  result <- do.call(rbind, lapply(names(posteriors), function(which_level) {
    .multilpa_average_posterior_rows(posteriors[[which_level]], which_level)
  }))
  row.names(result) <- NULL
  result
}

#' Mean posterior of every class within each modal assignment group
#' @param probabilities Posterior probability matrix.
#' @return A square matrix, rows indexed by assigned class, columns by class.
#' @noRd
.multilpa_average_posterior_matrix <- function(probabilities) {
  stopifnot("`probabilities` must be a numeric matrix" =
              is.matrix(probabilities) && is.numeric(probabilities))
  n_classes <- ncol(probabilities)
  assigned <- max.col(probabilities, ties.method = "first")
  t(vapply(seq_len(n_classes), function(class) {
    rows <- assigned == class
    if (!any(rows)) rep(NA_real_, n_classes) else
      colMeans(probabilities[rows, , drop = FALSE])
  }, numeric(n_classes)))
}

#' One row per class of classification diagnostics for one level
#' @param probabilities Posterior probability matrix.
#' @param level Label recorded in the `level` column.
#' @return A tidy `data.frame`, one row per class.
#' @noRd
.multilpa_classification_summary <- function(probabilities, level) {
  average <- .multilpa_average_posterior_matrix(probabilities)
  n_classes <- ncol(probabilities)
  classes <- seq_len(n_classes)
  assigned <- max.col(probabilities, ties.method = "first")
  n_modal <- tabulate(assigned, nbins = n_classes)
  estimated_n <- colSums(probabilities)
  estimated_proportion <- estimated_n / nrow(probabilities)
  diagonal <- average[cbind(classes, classes)]
  data.frame(level = level, class = classes, n_modal = n_modal,
             proportion_modal = n_modal / nrow(probabilities),
             estimated_n = estimated_n,
             estimated_proportion = estimated_proportion,
             average_posterior = diagonal,
             odds_correct_classification =
               .multilpa_odds_correct(diagonal, estimated_proportion))
}

#' One row per assigned class and class of average posteriors for one level
#' @param probabilities Posterior probability matrix.
#' @param level Label recorded in the `level` column.
#' @return A tidy `data.frame`, one row per assigned class and class.
#' @noRd
.multilpa_average_posterior_rows <- function(probabilities, level) {
  average <- .multilpa_average_posterior_matrix(probabilities)
  n_classes <- ncol(probabilities)
  classes <- seq_len(n_classes)
  assigned <- max.col(probabilities, ties.method = "first")
  data.frame(level = level,
             assigned_class = rep(classes, each = n_classes),
             class = rep(classes, times = n_classes),
             n_assigned = rep(tabulate(assigned, nbins = n_classes),
                              each = n_classes),
             average_posterior = as.vector(t(average)))
}

#' Odds of correct classification
#' @param average Average posterior probability in the assigned class.
#' @param proportion Model-estimated class proportion.
#' @return Odds ratio, or `NA_real_` where the ratio is undefined.
#' @noRd
.multilpa_odds_correct <- function(average, proportion) {
  stopifnot("`average` and `proportion` must have equal length" =
              length(average) == length(proportion))
  defined <- !is.na(average) & average > 0 & average < 1 &
    proportion > 0 & proportion < 1
  odds <- rep(NA_real_, length(average))
  odds[defined] <- (average[defined] / (1 - average[defined])) /
    (proportion[defined] / (1 - proportion[defined]))
  odds
}

#' Relative entropy of a fitted multilevel latent profile model
#'
#' Reports the zero-to-one relative entropy at each level, alongside the raw
#' classification entropy the entropy-penalized information criteria use.
#'
#' @param object A fitted `multilpa` model.
#' @return A base `data.frame`, one row per level the fit has, with the columns
#'   \describe{
#'     \item{`level`}{character: `"individuals"`, and `"groups"` when the fit
#'       has discrete group classes.}
#'     \item{`n_classes`}{integer: classes at that level.}
#'     \item{`n_units`}{integer: units at that level.}
#'     \item{`entropy_sum`}{numeric: the classification entropy `EN`, the
#'       quantity `awe`, `icl` and `clc` penalize, so lower is sharper.}
#'     \item{`relative_entropy`}{numeric: `1 - EN / (n log K)`, on the zero-to-one
#'       scale, so higher is sharper.}
#'   }
#'   `relative_entropy` is `NA_real_` when a level has a single class, where it
#'   is undefined rather than perfect. Continuous random-intercept fits have no
#'   discrete group classes and return the individual level only.
#' @seealso [classification_table()] for per-class classification quality and
#'   [information_criteria()] for the criteria that use `entropy_sum`.
#' @examples
#' set.seed(7)
#' example_data <- data.frame(
#'   school = rep(seq_len(12), each = 10),
#'   score_a = rnorm(120), score_b = rnorm(120)
#' )
#' fit <- multilpa(example_data, c("score_a", "score_b"), "school",
#'                   n_profiles = 2, n_group_classes = 1, n_starts = 2, seed = 1)
#' entropy_table(fit)
#' @export
entropy_table <- function(object) {
  stopifnot("`object` must be a fitted model of this package" = .multilpa_any_fit(object))
  posteriors <- list(individuals = object$subject_posteriors,
                     groups = object$group_posteriors)
  posteriors <- Filter(Negate(is.null), posteriors)
  result <- do.call(rbind, lapply(names(posteriors), function(level) {
    probabilities <- posteriors[[level]]
    data.frame(level = level, n_classes = ncol(probabilities),
               n_units = nrow(probabilities),
               entropy_sum = .multilpa_entropy_sum(probabilities),
               relative_entropy = .multilpa_relative_entropy(probabilities))
  }))
  row.names(result) <- NULL
  result
}
