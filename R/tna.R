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

#' A transition network for the whole sample
#'
#' Builds one [tna::tna()] model from a fitted transition model, aggregated
#' over every latent group class.
#'
#' The aggregate is formed from the expected transition counts, summed across
#' classes and then normalised, rather than by averaging the class transition
#' matrices. The two differ: averaging weights each class by how probable it is,
#' and summing counts weights it by how much transition mass it actually
#' contributes, which is what a marginal transition probability means. A class
#' holding a tenth of the students but a fifth of the observed moves counts for
#' the latter.
#'
#' @param x A fitted model from [lta()].
#' @param ... Passed to [tna::tna()].
#' @return An object of class `tna`, as [tna::tna()] returns: every verb of that
#'   package applies to it, including `centralities()`, `communities()`,
#'   `cliques()` and `plot()`.
#' @seealso [get_group_tna()] for one network per latent class.
#' @examples
#' if (requireNamespace("tna", quietly = TRUE)) {
#'   activity <- c("browse", "lectures", "forum_read")
#'   moves <- lta(course_engagement, vars = activity, id = "student",
#'                            time = "sequence", n_profiles = 2,
#'                            n_group_classes = 2, n_starts = 2, seed = 1)
#'   get_tna(moves)
#' }
#' @export
get_tna <- function(x, ...) {
  UseMethod("get_tna")
}

#' @rdname get_tna
#' @export
get_tna.multilpa_transitions <- function(x, ...) {
  .multilpa_require_tna()
  aggregate <- .multilpa_aggregate_transitions(x)
  tna::tna(aggregate$probabilities, inits = aggregate$initial, ...)
}

#' A transition network for each latent group class
#'
#' Builds one [tna::tna()] model per latent group class of a fitted transition
#' model, collected into the `group_tna` object tna's grouped verbs expect.
#'
#' @param x A fitted model from [lta()].
#' @param label What the classes are called in tna's output.
#' @param ... Passed to [tna::tna()] for each class.
#' @return An object of class `group_tna`, one `tna` model per latent class:
#'   `centralities()` returns one tidy table with a `group` column,
#'   `compare()` contrasts two classes, and `plot()` draws them together.
#' @seealso [get_tna()] for one network over the whole sample.
#' @examples
#' if (requireNamespace("tna", quietly = TRUE)) {
#'   activity <- c("browse", "lectures", "forum_read")
#'   moves <- lta(course_engagement, vars = activity, id = "student",
#'                            time = "sequence", n_profiles = 2,
#'                            n_group_classes = 2, n_starts = 2, seed = 1)
#'   get_group_tna(moves)
#' }
#' @export
get_group_tna <- function(x, ...) {
  UseMethod("get_group_tna")
}

#' @rdname get_group_tna
#' @export
get_group_tna.multilpa_transitions <- function(x, label = "Group class", ...) {
  .multilpa_require_tna()
  stopifnot("`label` must be a single string" =
              is.character(label) && length(label) == 1L)
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

#' Is the tna package available?
#' @return `NULL`, invisibly; raises `multilpa_missing_package` when it is not.
#' @noRd
.multilpa_require_tna <- function() {
  if (!requireNamespace("tna", quietly = TRUE)) {
    stop(errorCondition(paste(
      "The tna package is needed to build a transition network.",
      "Install it with install.packages(\"tna\")."),
      class = "multilpa_missing_package", call = NULL))
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
#' @return A list with `probabilities`, a square matrix, and `initial`.
#' @noRd
.multilpa_aggregate_transitions <- function(x) {
  states <- .multilpa_state_labels(x)
  ## Expected counts, not probabilities: summing the counts weights each class
  ## by the transition mass it carries, which is the marginal probability;
  ## averaging the matrices would weight it by class size alone.
  counts <- apply(x$transition_counts, seq_len(2L), sum)
  totals <- rowSums(counts)
  ## A state no observation ever leaves contributes no row to normalise. Its row
  ## is reported as a self-transition rather than as NaN, which is what "never
  ## left" means and what every downstream centrality can read.
  probabilities <- counts / ifelse(totals > 0, totals, 1)
  empty <- which(totals <= 0)
  if (length(empty) > 0L) probabilities[cbind(empty, empty)] <- 1
  dimnames(probabilities) <- list(states, states)
  initial <- as.numeric(x$group_probabilities %*% x$initial_probabilities)
  list(probabilities = probabilities,
       initial = stats::setNames(initial / sum(initial), states))
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
