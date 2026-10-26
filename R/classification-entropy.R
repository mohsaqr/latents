# The entropy of the complete-data classification, which ICL, AWE and CLC
# charge (review 2026-10-26 round 2, S-1).

#' Each row's expected classification entropy within its group class
#'
#' `E_C[H(X_i | C, y)] = -sum_h sum_k J_hik log(J_hik / omega_jh)`, with
#' `J_hik = P(X_i = k, C_j = h | y_j)` the E-step's `joint` and `omega_jh`
#' the posterior of the row's group. Added over the rows of a group and to the
#' group-class entropy of the group, it gives the entropy of the group's
#' complete-data classification (its group class and its rows' profiles).
#' With one group class it is the row's own posterior entropy.
#'
#' @param joint The E-step's list of rows-by-profiles matrices, one per group
#'   class, summing over the classes to the row posteriors.
#' @param group_posteriors Groups by group classes.
#' @param group_index The group of each row.
#' @return One non-negative number per row.
#' @noRd
.latents_row_conditional_entropy <- function(joint, group_posteriors, group_index) {
  stopifnot("`joint` must hold one matrix per group class" =
              is.list(joint) && length(joint) == ncol(group_posteriors),
            "`group_index` must give the group of every row" =
              length(group_index) == nrow(joint[[1L]]))
  Reduce(`+`, lapply(seq_along(joint), function(class) {
    piece <- joint[[class]]
    share <- group_posteriors[group_index, class]
    positive <- piece > 0
    terms <- matrix(0, nrow(piece), ncol(piece))
    # A positive joint probability has a positive group share, so both logs
    # are finite.
    terms[positive] <- piece[positive] *
      (log(piece[positive]) - log(share[row(piece)[positive]]))
    -rowSums(terms)
  }))
}

#' The entropy of a fit's complete-data classification
#'
#' The sum over the classified groups of their group-class entropy and over
#' the classified rows of their entropy within the group class,
#' `H(C | y) + sum_h P(C = h | y) sum_i H(X_i | C = h, y)`, each unit
#' weighted by its sampling weight. This is the entropy term of the
#' classification likelihood `log L_c = log L - EN` in its soft (posterior
#' expectation) form (Biernacki, Celeux & Govaert 2000; McLachlan & Peel 2000,
#' section 6.10.3). A fit with one group class has no group-class entropy, and
#' its rows' entropy within the class is their posterior entropy, which needs
#' no stored quantity. A fit with several group classes needs the rows'
#' within-class entropy stored when it was fitted
#' (`classification_entropy_rows`); a fit saved before that has none, and its
#' entropy is `NA`.
#'
#' @param x A profile or class model (not a transition model, whose profiles
#'   are a Markov chain within a person).
#' @return A single number, or `NA_real_`.
#' @noRd
.multilpa_complete_entropy <- function(x) {
  rows <- .multilpa_informative_rows(x)
  row_weights <- .multilpa_level_weights(x, "individuals", sum(rows)) %||% 1
  groups <- x$group_posteriors
  several <- !is.null(groups) && ncol(groups) > 1L
  within <- x$classification_entropy_rows
  if (is.null(within)) {
    if (several) return(NA_real_)
    posteriors <- .multilpa_class_posteriors(x)
    within <- -rowSums(ifelse(posteriors > 0, posteriors * log(posteriors), 0))
  }
  stopifnot("one within-class entropy per row" = length(within) == length(rows))
  row_part <- sum(row_weights * within[rows])
  if (!several) return(row_part)
  counted <- .multilpa_classified_group_posteriors(x)
  row_part + .multilpa_entropy_sum(
    counted, .multilpa_level_weights(x, "groups", nrow(counted)))
}

#' Relative entropy of one level of a fit, weighted as the fit is
#'
#' @param fit A fitted profile, class or transition model.
#' @param level `"individuals"` or `"groups"`.
#' @return The level's relative entropy (sampling-weighted for a weighted
#'   fit), or `NA_real_` for a level the fit lacks or with one class.
#' @noRd
.multilpa_level_relative_entropy <- function(fit, level) {
  posteriors <- .multilpa_level_posteriors(fit, level)
  if (is.null(posteriors)) return(NA_real_)
  .multilpa_relative_entropy(posteriors,
                             .multilpa_level_weights(fit, level, nrow(posteriors)))
}

#' The AvePP summary of one level of a fit, weighted as the fit is
#'
#' @param fit A fitted profile, class or transition model.
#' @param level `"individuals"` or `"groups"`.
#' @param cutoff The AvePP below which a class is counted.
#' @return The vector of `.multilpa_avepp_summary()`.
#' @noRd
.multilpa_level_avepp <- function(fit, level, cutoff) {
  posteriors <- .multilpa_level_posteriors(fit, level)
  weights <- if (is.null(posteriors)) NULL else
    .multilpa_level_weights(fit, level, nrow(posteriors))
  .multilpa_avepp_summary(posteriors, cutoff, weights)
}

#' The classified posteriors of one level
#' @param fit A fitted profile, class or transition model.
#' @param level `"individuals"` or `"groups"`.
#' @return A posterior matrix, or `NULL` for a level the fit lacks.
#' @noRd
.multilpa_level_posteriors <- function(fit, level) {
  if (identical(level, "groups")) .multilpa_classified_group_posteriors(fit) else
    .multilpa_classified_posteriors(fit)
}
