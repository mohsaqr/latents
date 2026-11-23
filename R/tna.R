# Handing a fitted transition model to the tna package.
#
# `lta()` estimates exactly what a transition network is: a square
# matrix of probabilities over the same states, and where the sequences start.
# tna builds its models from those two things, so the handover is the estimates
# themselves rather than a re-analysis of the assignments.
#
# That is the difference worth keeping. tna's own method for a fitted mixture
# model, `group_model.mhmm()`, takes the model's *grouping* and counts
# transitions between modal assignments, discarding the model's estimated
# matrices. Here the estimates are what is handed over: they are maximised over
# the full posteriors, so an observation split 0.6/0.4 between two profiles
# contributes to both, where a modal count gives it entirely to one.

#' Get the transitions between a fitted model's profiles as one transition
#' network for the whole sample
#'
#' @description
#' `get_tna()` returns the transitions between the profiles (classes) of a
#' fitted model as one transition network for the whole sample, an object of the
#' tna package. The profiles are the nodes, named `profile_1`, `profile_2`, and
#' so on. The weight of the edge from one profile to another is the probability
#' of moving to it at the next occasion, and the weight of a loop is the
#' probability of staying.
#'
#' Where the probabilities come from depends on the model:
#'
#' * for a transition fit from [lta()], they are the transition
#'   probabilities the model estimated, pooled over its group classes, and
#'   the network starts from the estimated initial probabilities;
#' * for a profile or class fit made with `time`, they are counted by
#'   [tna::tna()] from each group's sequence of assigned (most likely)
#'   profiles.
#'
#' @details
#' For a transition fit, the pooled matrix is built from the expected
#' transition counts, summed over the group classes and divided by their
#' row totals. Each group class therefore counts in proportion to the
#' transitions it is expected to make, not to its class probability. If a
#' profile has no expected moves out of it in any class, its row is the
#' class rows averaged by class probability, and a warning of class
#' `latents_empty_transition_row` says that the row is not estimated from
#' transitions; a missing transition is never counted as staying.
#'
#' @section Which fits are supported:
#' * An [lta()] fit of class `multilpa_transitions` (the basic transition
#'   model: continuous or categorical indicators and one first-order
#'   transition matrix for all occasions, without covariates, weights or a
#'   stayer class): always.
#' * An [lta()] fit of class `multilpa_lta` (ordinal or count indicators,
#'   weights, covariates of the initial profiles, occasion-specific
#'   measurement): when its transitions are first order, the same at every
#'   occasion and without transition covariates. A fit whose transitions
#'   vary by occasion, by covariate or by the previous profile (`order = 2`)
#'   has no single matrix and raises `latents_unsupported_tna`; read its
#'   transitions with `get_results(x, "transitions")`.
#' * A mover-stayer fit (`mover_stayer = TRUE`): not by `get_tna()`, which
#'   raises `latents_unsupported_tna` because the stayer class records no
#'   expected moves to pool. [get_group_tna()] returns one network per
#'   mover class and one for the stayers.
#' * A profile or class fit from [lpa()], [lca()], [multilpa()] or
#'   [multilca()] (with or without membership covariates): when it was
#'   made with `time`; without it the fit raises `latents_no_time`.
#' * Any other model (a mixture regression, a growth mixture, an additive
#'   or cross-level fit, pooled imputations) raises `latents_not_a_fit`.
#'
#' @section Interpretation:
#' A network counted from assigned profiles treats each assignment as
#' known and does not carry the uncertainty of the classification. The
#' network of a transition fit uses probabilities estimated from the
#' posterior probabilities, so an observation that may belong to two
#' profiles contributes to both. A network pooled over the group classes
#' can hide differences between them; [get_group_tna()] returns their
#' separate networks.
#'
#' @param x A transition fit from [lta()], or a profile or class fit made
#'   with `time`. See *Which fits are supported*.
#' @param ... Further arguments passed to [tna::tna()].
#' @return An object of class `tna`, as [tna::tna()] returns it: a list
#'   with `weights` (the square matrix of transition probabilities, rows
#'   the current profile and columns the next, each row summing to one),
#'   `inits` (the probability of starting in each profile), `labels` (the
#'   profile names) and `data` (the sequences, for a profile or class
#'   fit; `NULL` for a transition fit). The functions of the tna package
#'   apply to it, among them `centralities()`, `communities()`,
#'   `cliques()` and `plot()`.
#'
#'   Conditions raised: `latents_no_time`, `latents_unsupported_tna` and
#'   `latents_not_a_fit` as described under *Which fits are supported*;
#'   `latents_missing_package` when the tna package is not installed; and
#'   the warning `latents_empty_transition_row`.
#' @seealso [get_group_tna()] for one network per latent group class;
#'   [plot_transitions()] for the transition matrix as a plot;
#'   `get_results(x, "transitions")` for it as a table.
#' @examples
#' if (requireNamespace("tna", quietly = TRUE)) {
#'   activity <- c("browse", "lectures", "forum_read")
#'   moves <- lta(course_engagement, vars = activity, id = "student",
#'                time = "sequence", n_profiles = 2, n_group_classes = 2,
#'                n_starts = 2, seed = 1)
#'   get_tna(moves)
#' }
#' @export
get_tna <- function(x, ...) {
  UseMethod("get_tna")
}

#' @noRd
#' @export
get_tna.default <-function(x, ...) {
  .latents_require_fit(x, substitute(x))
}

#' @noRd
#' @export
get_tna.multilpa_transitions <-function(x, ...) {
  .multilpa_require_tna()
  aggregate <- .multilpa_aggregate_transitions(x)
  if (any(!aggregate$estimated)) {
    .multilpa_warn_tna_empty_rows()
  }
  tna::tna(aggregate$probabilities, inits = aggregate$initial, ...)
}

#' Get one transition network per latent group class, to compare how their
#' members move between profiles
#'
#' @description
#' `get_group_tna()` returns one transition network per latent group
#' class, as a grouped object of the tna package, so that the group classes
#' can be compared on how their members stay in and move between profiles.
#' Every network has the same nodes, the profiles (classes) `profile_1`,
#' `profile_2`, and so on, and its edge weights are the probabilities of
#' moving from one profile to another at the next occasion.
#'
#' Where the probabilities come from depends on the model:
#'
#' * for a transition fit from [lta()], each network is the transition
#'   matrix the model estimated for that group class, and starts from the
#'   class's estimated initial probabilities;
#' * for a profile or class fit made with `time`, the groups are split by
#'   their assigned (most likely) group class, and [tna::group_tna()]
#'   counts each network from the sequences of assigned profiles of the
#'   groups in that class.
#'
#' @details
#' A transition-matrix row with no expected moves out of its profile keeps
#' its fitted probabilities, and a warning of class
#' `latents_empty_transition_row` says that the row is not estimated from
#' transitions.
#'
#' @section Which fits are supported:
#' The same fits as [get_tna()], and mover-stayer fits as well:
#'
#' * An [lta()] fit of class `multilpa_transitions`: always.
#' * An [lta()] fit of class `multilpa_lta`: when its transitions are
#'   first order, the same at every occasion and without transition
#'   covariates; otherwise it raises `latents_unsupported_tna`.
#' * A mover-stayer fit (`mover_stayer = TRUE`): one network per mover
#'   class and one named `"Stayers"`, whose matrix is the identity.
#' * A profile or class fit from [lpa()], [lca()], [multilpa()] or
#'   [multilca()] (with or without membership covariates): when it was
#'   made with `time`; without it the fit raises `latents_no_time`.
#' * Any other model raises `latents_not_a_fit`.
#'
#' @section Interpretation:
#' The networks describe latent group classes, not the individual observed
#' groups. Networks counted from assigned profiles treat both the profile
#' and the group-class assignments as known and carry the uncertainty of
#' neither.
#'
#' @param x A transition fit from [lta()], or a profile or class fit made
#'   with `time`. See *Which fits are supported*.
#' @param label One string: what the classes are called in the names of
#'   the networks and in tna's output (default `"Group class"`, which
#'   gives `"Group class 1"`, `"Group class 2"`, and so on).
#' @param ... Further arguments passed to [tna::tna()] for each class of a
#'   transition fit, or to [tna::group_tna()] for a profile or class fit.
#' @return An object of class `group_tna`: a named list of `tna` networks,
#'   one per latent group class, each as [get_tna()] describes. The
#'   grouped functions of the tna package apply to it: `centralities()`
#'   returns one table with a `group` column, `compare()` contrasts two
#'   classes, and `plot()` draws every network.
#'
#'   Conditions raised: `latents_no_time`, `latents_unsupported_tna` and
#'   `latents_not_a_fit` as described under *Which fits are supported*;
#'   `latents_missing_package` when the tna package is not installed;
#'   `latents_bad_argument` when `label` is not one string; and the warning
#'   `latents_empty_transition_row`.
#' @seealso [get_tna()] for one network over the whole sample;
#'   [plot_transitions()] for the transition matrices as a plot.
#' @examples
#' if (requireNamespace("tna", quietly = TRUE)) {
#'   activity <- c("browse", "lectures", "forum_read")
#'   moves <- lta(course_engagement, vars = activity, id = "student",
#'                time = "sequence", n_profiles = 2, n_group_classes = 2,
#'                n_starts = 2, seed = 1)
#'   get_group_tna(moves)
#' }
#' @export
get_group_tna <- function(x, ...) {
  UseMethod("get_group_tna")
}

#' @noRd
#' @export
get_group_tna.default <-function(x, ...) {
  .latents_require_fit(x, substitute(x))
}

#' @rdname get_group_tna
#' @export
get_group_tna.multilpa_transitions <- function(x, label = "Group class", ...) {
  .multilpa_require_tna()
  .latents_require("`label` must be a single string" =
                     is.character(label) && length(label) == 1L && !is.na(label))
  if (any(x$empty_transition_rows %||% FALSE)) {
    .multilpa_warn_tna_empty_rows()
  }
  states <- .multilpa_state_labels(x)
  names <- sprintf("%s %d", label, seq_len(x$n_group_classes))
  models <- stats::setNames(lapply(seq_len(x$n_group_classes), function(class) {
    probabilities <- matrix(x$transition_probabilities[, , class],
                            length(states), length(states),
                            dimnames = list(states, states))
    tna::tna(probabilities, inits = as.numeric(x$initial_probabilities[class, ]),
             ...)
  }), names)
  .multilpa_as_group_tna(models, label = label)
}

#' @noRd
#' @export
get_tna.multilpa_lta <-function(x, ...) {
  .multilpa_require_tna()
  if (isTRUE(x$mover_stayer)) {
    stop(errorCondition(paste(
      "get_tna() pools the classes by their expected moves, which the stayer",
      "class does not record. Use get_group_tna(x): one network per mover",
      "class and one for the stayers."), class = "latents_unsupported_tna",
      call = NULL))
  }
  networks <- .lta_class_networks(x, "get_tna()")
  states <- rownames(networks[[1L]]$probabilities)
  # Expected moves, not probabilities, are summed over the classes, as for a
  # homogeneous transition fit (.multilpa_aggregate_transitions()).
  counts <- Reduce(`+`, lapply(networks, function(network) {
    network$probabilities * network$origin_count
  }))
  totals <- rowSums(counts)
  if (any(totals <= 0)) {
    .multilpa_warn_tna_empty_rows(
      "the rows of get_results(fit, \"transitions\") whose `estimated` is FALSE")
  }
  averaged <- Reduce(`+`, Map(function(network, share) network$probabilities * share,
                              networks, x$group_probabilities))
  probabilities <- counts / ifelse(totals > 0, totals, 1)
  probabilities[totals <= 0, ] <- averaged[totals <= 0, , drop = FALSE]
  initial <- Reduce(`+`, Map(function(network, share) network$initial * share,
                             networks, x$group_probabilities))
  tna::tna(probabilities, inits = stats::setNames(initial / sum(initial), states), ...)
}

#' @noRd
#' @export
get_group_tna.multilpa_lta <-function(x, label = "Group class", ...) {
  .multilpa_require_tna()
  .latents_require("`label` must be a single string" =
                     is.character(label) && length(label) == 1L && !is.na(label))
  networks <- .lta_class_networks(x, "get_group_tna()")
  if (any(x$empty_transition_rows %||% FALSE)) {
    .multilpa_warn_tna_empty_rows(
      "the rows of get_results(fit, \"transitions\") whose `estimated` is FALSE")
  }
  classes <- names(x$group_probabilities)
  names <- ifelse(classes == "stayers", "Stayers",
                  sprintf("%s %d", label, seq_along(classes)))
  models <- stats::setNames(lapply(networks, function(network) {
    tna::tna(network$probabilities, inits = network$initial, ...)
  }), names)
  .multilpa_as_group_tna(models, label = label)
}

#' Each group class's transition matrix, initial distribution and support
#'
#' Read from the fit's own transition and initial tables, which average the
#' model-implied probabilities over the groups; under homogeneous first-order
#' transitions without transition covariates every group has the same
#' matrix, so the table is that matrix.
#' @param x A `multilpa_lta` fit.
#' @param verb The calling verb, for the refusal.
#' @return A list over group classes (stayers included), each with
#'   `probabilities` (states by states, named), `initial` (named) and
#'   `origin_count` (expected moves out of each state; `NA` for stayers).
#' @noRd
.lta_class_networks <- function(x, verb) {
  if (!identical(x$transitions, "homogeneous") ||
      length(x$transition_covariates) > 0L || (x$order %||% 1L) != 1L) {
    .lta_refuse_tna(verb)
  }
  transitions <- .lta_table(x, "transitions")
  initial <- .lta_table(x, "initial")
  states <- sprintf("profile_%d", seq_len(x$n_profiles))
  lapply(names(x$group_probabilities), function(class) {
    rows <- transitions[transitions$group_class == class, , drop = FALSE]
    starts <- initial[initial$group_class == class, , drop = FALSE]
    probabilities <- matrix(0, length(states), length(states),
                            dimnames = list(states, states))
    probabilities[cbind(match(rows$from, states), match(rows$to, states))] <-
      rows$probability
    origin <- rows[!duplicated(rows$from), , drop = FALSE]
    list(probabilities = probabilities,
         initial = stats::setNames(starts$probability[match(states, starts$profile)],
                                   states),
         origin_count = origin$origin_expected_count[match(states, origin$from)])
  })
}

#' Refuse a transition network for a general transition fit
#' @param verb The verb, for the message.
#' @return Never returns; raises `latents_unsupported_tna`.
#' @noRd
.lta_refuse_tna <- function(verb) {
  stop(errorCondition(sprintf(paste(
    "%s needs one transition matrix per group class, and this general",
    "transition fit has none: its transitions vary by occasion, covariate or",
    "previous profile. Its transitions are get_results(x,",
    "\"transitions\")."), verb), class = "latents_unsupported_tna", call = NULL))
}

#' Warn when a network contains fitted transition rows without move data
#' @param where Where the user can see which rows those are, for the message.
#' @return `NULL`, invisibly.
#' @noRd
.multilpa_warn_tna_empty_rows <- function(
    where = "get_results(fit, \"transitions\", estimated = FALSE)") {
  warning(warningCondition(paste0(
    "A transition row has no expected outgoing moves, so its probabilities ",
    "are not estimated from transitions. The network retains the fitted ",
    "row; inspect ", where, "."),
    class = "latents_empty_transition_row", call = NULL))
  invisible(NULL)
}

#' Is the tna package available?
#' @return `NULL`, invisibly; raises `latents_missing_package` when it is not.
#' @noRd
.multilpa_require_tna <- function() {
  if (!requireNamespace("tna", quietly = TRUE)) {
    stop(errorCondition(paste(
      "The tna package is needed to build a transition network.",
      "Install it with install.packages(\"tna\")."),
      class = "latents_missing_package", call = NULL))
  }
  invisible(NULL)
}

#' State labels for a transition network
#' @param x A fitted model from [lta()].
#' @return A character vector, one per profile.
#' @noRd
.multilpa_state_labels <- function(x) {
  labels <- dimnames(x$transition_probabilities)[[1L]]
  if (is.null(labels)) sprintf("profile_%d", seq_len(x$n_profiles)) else labels
}

#' The transition model marginal over latent group classes
#'
#' @param x A fitted model from [lta()].
#' @return A list with `probabilities`, a square matrix; `initial`, the state
#'   distribution; and `estimated`, marking rows with positive expected counts.
#' @noRd
.multilpa_aggregate_transitions <- function(x) {
  states <- .multilpa_state_labels(x)
  ## Expected counts, not probabilities: summing the counts weights each class
  ## by the transition mass it carries, which is the marginal probability;
  ## averaging the matrices would weight it by class size alone.
  counts <- apply(x$transition_counts, seq_len(2L), sum)
  totals <- rowSums(counts)
  ## A state with no outgoing moves has no count-based transition estimate.
  ## Keep the model's own row as a labelled fallback. A self-transition of one
  ## would falsely turn missing transition information into evidence of staying.
  probabilities <- counts / ifelse(totals > 0, totals, 1)
  ## The fallback for every row at once: flattening the from x to x class
  ## array to (from, to) x class turns the class-probability average into one
  ## matrix product. Only the empty rows keep it.
  n_states <- x$n_profiles
  averaged <- matrix(matrix(x$transition_probabilities,
                            n_states * n_states, x$n_group_classes) %*%
                       x$group_probabilities, n_states, n_states)
  empty <- totals <= 0
  probabilities[empty, ] <- averaged[empty, , drop = FALSE] /
    rowSums(averaged[empty, , drop = FALSE])
  dimnames(probabilities) <- list(states, states)
  initial <- as.numeric(x$group_probabilities %*% x$initial_probabilities)
  list(probabilities = probabilities,
       initial = stats::setNames(initial / sum(initial), states),
       estimated = totals > 0)
}

#' Collect per-class models into tna's grouped object
#'
#' tna builds a `group_tna` from sequence data, so there is no public
#' constructor for one made from estimated matrices. The object is a named list
#' of `tna` models carrying the attributes tna's grouped verbs read; `groups`
#' and `cols` describe the sequence data a model built from matrices does not
#' have, and are empty rather than invented.
#'
#' @param models A named list of `tna` models, one per class.
#' @param label What the classes are called.
#' @return An object of class `group_tna`.
#' @noRd
.multilpa_as_group_tna <- function(models, label) {
  structure(models,
            groups = vector("list", length(models)),
            label = label,
            levels = names(models),
            na.rm = TRUE,
            cols = character(0L),
            groupwise = FALSE,
            type = "relative",
            scaling = character(0L),
            class = "group_tna")
}

#' @noRd
#' @export
get_tna.multilpa <-function(x, ...) {
  .multilpa_require_tna()
  tna::tna(.multilpa_modal_sequences(x)$states, ...)
}

#' @noRd
#' @export
get_tna.multilpa_covariates <-get_tna.multilpa

#' @noRd
#' @export
get_group_tna.multilpa <-function(x, label = "Group class", ...) {
  .multilpa_require_tna()
  .latents_require("`label` must be a single string" =
                     is.character(label) && length(label) == 1L && !is.na(label))
  sequences <- .multilpa_modal_sequences(x)
  grouped <- cbind(sequences$states,
                   .group = sprintf("%s %d", label, sequences$group_class))
  tna::group_tna(grouped, group = ".group", ...)
}

#' @noRd
#' @export
get_group_tna.multilpa_covariates <-get_group_tna.multilpa

#' Each group's modal profiles in time order, one row per group
#'
#' Only classified rows enter (`get_results(x, "sequences")` leaves the others
#' out): a row with no observed indicator has the class
#' prior for its posterior, so its modal profile is not an observed state. It
#' is left out, as an unclassified row is left out of the classification
#' statistics and the three-step assignments
#' (`.multilpa_informative_rows()`), rather than marked `NA`: a sequence here
#' already runs from one observed row to the next, so an empty row is treated
#' exactly as an absent one, and the network does not depend on which of the
#' two the data used. A group with no classified row has no sequence.
#' `drop_empty = FALSE` keeps every row, as before.
#'
#' @param x A fitted model made with `time =`.
#' @return A list with `states`, a data frame with one row per group with a
#'   classified row and one column per position (`T1`, `T2`, ...), holding
#'   `profile_k` labels and `NA` after a group's last classified row; and
#'   `group_class`, each group's modal group class, in the same row order.
#' @noRd
.multilpa_modal_sequences <- function(x) {
  # The sequence table already leaves the unclassified rows out.
  sequences <- get_results(x, "sequences")
  sequences <- sequences[order(sequences$group, sequences$time), , drop = FALSE]
  ids <- unique(sequences$group)
  groups <- factor(match(sequences$group, ids), levels = seq_along(ids))
  by_group <- split(sprintf("profile_%d", sequences$profile), groups)
  width <- max(lengths(by_group))
  states <- as.data.frame(
    do.call(rbind, lapply(by_group, function(profiles) {
      c(profiles, rep(NA_character_, width - length(profiles)))
    })), stringsAsFactors = FALSE)
  names(states) <- paste0("T", seq_len(width))
  rownames(states) <- NULL
  first <- match(ids, sequences$group)
  stopifnot("one row per group" = nrow(states) == nlevels(groups))
  list(states = states, group_class = sequences$group_class[first])
}
