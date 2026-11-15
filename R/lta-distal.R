# Transition-aware three-step analysis of a distal outcome: a person-level
# outcome related to the profile an lta() chain occupies at one occasion, or
# to the pair of profiles at two occasions (a transition). Step one is the
# fitted lta(); step two classifies each occasion from its own indicators
# (.lta_occasion_posteriors(), the step two of the three-step latent Markov
# model of r3step()) and estimates one pooled classification-error matrix D;
# step three relates the outcome to the latent profile(s).
#
# - BCH: an occasion's assignment errs independently of the other occasions'
#   given the profiles, so a person's weight for the pair (k, l) is the
#   product of the rows of D^-1 for the two assignments, and for one occasion
#   the row for that occasion's assignment (Bolck, Croon and Hagenaars 2004;
#   Vermunt 2010).
# - ML: the assignments at every occasion are indicators of the chain with
#   the error rates D, the chain's initial and transition probabilities are
#   held at the lta() estimates (as the profile-level ML step holds the class
#   sizes), and the outcome is normal (or multinomial) given the profile at
#   the chosen occasion(s). A person's fixed weight on profile k is
#   P(s_t = k | assignments) from a forward-backward pass of that chain, so
#   the person's other occasions inform the profile at occasion t.

#' The step-three units, assignments and weights of a distal outcome of a chain
#'
#' @param x A covariate-free homogeneous `lta()` fit (class
#'   `multilpa_transitions`) with one group class.
#' @param data The data it was fitted to, with the (coded) outcome.
#' @param outcome The outcome column, constant within each person.
#' @param occasion One or two increasing positions in the chain.
#' @param method The three-step method.
#' @return A list with `pieces` (as `.multilpa_level_assignments()`, plus
#'   `step_weights` for BCH or `ml_weights` for ML) and `values` (one per
#'   person entering step three).
#' @noRd
.lta_distal_inputs <- function(x, data, outcome, occasion, method) {
  unsupported <- function(...) {
    stop(errorCondition(paste(...), class = "latents_unsupported_three_step", call = NULL))
  }
  bad <- function(...) {
    stop(errorCondition(paste(...), class = "latents_bad_argument", call = NULL))
  }
  if (!inherits(x, "multilpa_transitions")) {
    if (inherits(x, "multilpa_lta")) .multilpa_check_three_step_fit(x)
    bad("`occasion` relates a distal outcome to the profile of an lta() chain at",
        "that occasion; it applies to lta() fits only.")
  }
  if (outcome %in% .multilpa_measured_items(x)) {
    stop(errorCondition(sprintf(
      "`%s` helped define the fitted classes, so it is not a distal outcome.",
      outcome), class = c("multilpa_indicator_reused", "latents_bad_outcome"),
      call = NULL))
  }
  if (x$n_group_classes > 1L) {
    unsupported(
      "A distal outcome of a transition chain is implemented for one group class:",
      "the step-three chain would also need the group classes.")
  }
  if (!is.numeric(occasion) || !length(occasion) %in% 1:2 || anyNA(occasion) ||
      any(occasion != round(occasion)) || any(occasion < 1)) {
    bad("`occasion` must be one or two positions in the sequence (whole numbers,",
        "the first occasion being 1).")
  }
  occasion <- as.integer(occasion)
  if (length(occasion) == 2L && occasion[1L] >= occasion[2L]) {
    bad("Two occasions must be given in increasing order, the earlier first:",
        "the pair is the profile at the first and the profile at the second.")
  }
  .lta_check_data(x, data)
  a <- .lta_arguments(x)
  if (outcome %in% c(a$id, a$time)) {
    stop(errorCondition(sprintf(
      "`%s` is the fit's %s column, not an outcome.", outcome,
      if (identical(outcome, a$id)) "id" else "time"),
      class = "latents_bad_outcome", call = NULL))
  }
  groups <- .multilpa_prepare_groups(data[[a$id]])
  time_values <- .multilpa_time_values(data, a$time, a$id, a$vars)
  layout <- .multilpa_sequence_layout(groups$index, time_values, groups$n, a$occasions)
  if (max(occasion) > layout$n_occasions) {
    bad(sprintf("`occasion` is beyond the longest sequence, which has %d occasions.",
                layout$n_occasions))
  }
  by_person <- split(data[[outcome]], groups$index)
  if (!all(vapply(by_person, function(v) length(unique(v)) == 1L, logical(1)))) {
    stop(errorCondition(paste(
      "A distal outcome of a transition chain is measured once per person, so",
      "it must be constant over each person's rows."),
      class = "latents_bad_outcome", call = NULL))
  }
  person_values <- vapply(by_person, function(v) v[[1L]], numeric(1))

  k_profiles <- x$n_profiles
  posteriors <- .lta_occasion_posteriors(x)
  observed <- x$informative %||% (rowSums(!is.na(data[a$vars])) > 0L)
  modal <- max.col(posteriors, ties.method = "first")
  errors <- .multilpa_error_matrix(list(posteriors = posteriors[observed, , drop = FALSE],
                                        modal = modal[observed], n_classes = k_profiles))
  # A slot holds the row of that occasion, NA where the person has none or
  # nothing was measured there.
  slot <- layout$slot
  slot[!is.na(slot) & !observed[pmax(slot, 1L)]] <- NA_integer_
  pair <- length(occasion) == 2L
  n_cells <- if (pair) k_profiles^2 else k_profiles
  # Cell (k, l) of a pair is column (k - 1) K + l: the later profile varies
  # fastest.
  combine <- function(first, second) {
    if (!pair) return(first)
    do.call(cbind, lapply(seq_len(k_profiles), function(k) first[, k] * second))
  }
  if (identical(method, "ml")) {
    # Every person with a measured occasion whose sequence reaches the
    # chosen occasions; a chosen occasion with nothing measured is inferred
    # from the person's other occasions through the chain.
    kept <- which(rowSums(!is.na(slot)) > 0L & layout$within[, max(occasion)])
    chain <- .lta_distal_chain(x, errors, modal, slot[kept, , drop = FALSE], occasion)
    pieces <- list(posteriors = chain, modal = max.col(chain, ties.method = "first"),
                   ml_weights = chain)
  } else {
    # The weighting methods need an assignment at every chosen occasion.
    rows <- slot[, occasion, drop = FALSE]
    kept <- which(rowSums(is.na(rows)) == 0L)
    rows <- rows[kept, , drop = FALSE]
    first <- rows[, 1L]
    second <- rows[, length(occasion)]
    pieces <- switch(method,
      bch = {
        inverse <- .multilpa_invert_errors(errors)
        list(step_weights = combine(inverse[modal[first], , drop = FALSE],
                                    inverse[modal[second], , drop = FALSE]))
      },
      proportional = list(posteriors = combine(posteriors[first, , drop = FALSE],
                                               posteriors[second, , drop = FALSE])),
      modal = list(modal = if (pair) (modal[first] - 1L) * k_profiles + modal[second] else
        modal[first]))
    # The proportional shares stand in as `posteriors` for every method.
    pieces$posteriors <- pieces$posteriors %||%
      combine(posteriors[first, , drop = FALSE], posteriors[second, , drop = FALSE])
    pieces$modal <- pieces$modal %||% max.col(pieces$posteriors, ties.method = "first")
  }
  if (length(kept) == 0L) {
    stop(errorCondition(sprintf(
      "No person has a measured response at occasion %s.", paste(occasion, collapse = " and ")),
      class = "latents_bad_outcome", call = NULL))
  }
  pieces$n_classes <- n_cells
  pieces$units <- kept
  pieces$group_index <- seq_along(kept)
  pieces$kept <- kept
  pieces$weights <- NULL
  pieces$classification_errors <- errors
  list(pieces = pieces, values = unname(person_values[kept]), occasion = occasion,
       n_profiles = k_profiles)
}

#' The step-three chain's posterior of the profile(s) at chosen occasions
#'
#' A forward-backward pass of the homogeneous chain whose single indicator at
#' each occasion is the assigned profile, observed with the error rates `D`,
#' with the initial and transition probabilities held at the lta() estimates.
#'
#' @param x The lta() fit.
#' @param errors The classification-error matrix, rows the true profile.
#' @param modal Each row's assigned profile.
#' @param slot Persons-by-occasions rows (NA where nothing was assigned).
#' @param occasion One or two positions.
#' @return A persons-by-cells matrix whose rows sum to one: P(s_t = k | m) or
#'   P(s_t1 = k, s_t2 = l | m), cell (k - 1) K + l.
#' @noRd
.lta_distal_chain <- function(x, errors, modal, slot, occasion) {
  k_profiles <- ncol(errors)
  n <- nrow(slot)
  last <- max(occasion)
  initial <- as.vector(x$initial_probabilities[1L, ])
  transition <- matrix(x$transition_probabilities[, , 1L], k_profiles, k_profiles)
  emission <- lapply(seq_len(last), function(t) {
    block <- matrix(1, n, k_profiles)
    here <- !is.na(slot[, t])
    block[here, ] <- t(errors[, modal[slot[here, t]], drop = FALSE])
    block
  })
  # Occasions after the last chosen one: their backward message.
  after <- if (ncol(slot) > last) seq.int(last + 1L, ncol(slot)) else integer()
  backward_tail <- Reduce(function(beta, t) {
    here <- !is.na(slot[, t])
    block <- matrix(1, n, k_profiles)
    block[here, ] <- t(errors[, modal[slot[here, t]], drop = FALSE])
    message <- (block * beta) %*% t(transition)
    message / rowSums(message)
  }, rev(after), init = matrix(1, n, k_profiles))
  normalize <- function(m) m / rowSums(m)
  # Forward messages, normalized each occasion (their scale cancels).
  alpha <- vector("list", last)
  alpha[[1L]] <- normalize(emission[[1L]] * matrix(initial, n, k_profiles, byrow = TRUE))
  for (t in seq_len(last)[-1L]) {  # recursion over occasions
    alpha[[t]] <- normalize((alpha[[t - 1L]] %*% transition) * emission[[t]])
  }
  if (length(occasion) == 1L) {
    return(normalize(alpha[[last]] * backward_tail))
  }
  # P(s_t1 = k, s_t2 = l, m) up to a constant: alpha_t1(k) times the paths
  # from k at t1 to l at t2 through the emissions between, times beta_t2(l).
  start <- occasion[1L]
  joint <- do.call(cbind, lapply(seq_len(k_profiles), function(k) {
    path <- matrix(0, n, k_profiles)
    path[, k] <- alpha[[start]][, k]
    for (u in seq.int(start + 1L, last)) {  # recursion over occasions
      path <- (path %*% transition) * emission[[u]]
    }
    path * backward_tail
  }))
  normalize(joint)
}

#' Label a step-three table with the occasion(s) and the profile pairs
#' @param result The table from the step-three path.
#' @param occasion One or two positions.
#' @param k_profiles The number of profiles.
#' @return The table with attribute `occasion`, and for a pair the columns
#'   `from` and `to` (and `reference_from`, `reference_to` for contrasts).
#' @noRd
.lta_distal_label <- function(result, occasion, k_profiles) {
  pair <- length(occasion) == 2L
  title <- attr(result, "title")
  where <- if (pair) sprintf("profiles at occasions %d and %d", occasion[1L], occasion[2L]) else
    sprintf("profile at occasion %d", occasion)
  if (!is.null(title)) {
    attr(result, "title") <- sub("by class", paste("by", where),
                                 sub("class differences", paste(where, "differences"), title))
  }
  if (pair) {
    from <- function(cell) (cell - 1L) %/% k_profiles + 1L
    to <- function(cell) (cell - 1L) %% k_profiles + 1L
    attrs <- attributes(result)
    if ("class" %in% names(result)) {
      position <- match("class", names(result))
      added <- data.frame(from = from(result$class), to = to(result$class))
      if ("reference_class" %in% names(result)) {
        added$reference_from <- from(result$reference_class)
        added$reference_to <- to(result$reference_class)
      }
      kept_attrs <- attrs[setdiff(names(attrs), c("names", "row.names"))]
      result <- cbind(as.data.frame(unclass(result)[seq_len(position)],
                                    stringsAsFactors = FALSE),
                      added,
                      as.data.frame(unclass(result)[-seq_len(position)],
                                    stringsAsFactors = FALSE))
      attributes(result)[names(kept_attrs)] <- kept_attrs
      display <- attr(result, "display")
      if (!is.null(display)) {
        column <- .latents_column
        display <- unlist(lapply(display, function(entry) {
          if (identical(entry$column, "class")) {
            list(column("From", "from", "integer"), column("To", "to", "integer"))
          } else if (identical(entry$column, "reference_class")) {
            list(column("vs from", "reference_from", "integer"),
                 column("vs to", "reference_to", "integer"))
          } else list(entry)
        }), recursive = FALSE)
        attr(result, "display") <- display
      }
    }
  }
  attr(result, "occasion") <- occasion
  result
}
