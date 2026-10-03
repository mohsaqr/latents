# The simulation service: seed scoping, categorical draws, the cluster
# bootstrap and the bootstrap likelihood-ratio replicate loop.
#
# Every verb that draws random numbers takes a `seed`, sets it, and gives the
# caller's random state back on exit; every engine that simulates draws
# classes and categories by inverse CDF; the bootstrap verbs resample groups
# and refit, or simulate and refit. Each engine had its own copy of these.
# They live here once, with the arithmetic and the order of the draws exactly
# as the engines had them, so a seeded result is the same number it was.
#
# What stays with the engines is what genuinely differs: how each model
# draws its indicators (`rnorm(n, mu, sd)` and `mu + sd * rnorm(n)` differ in
# the last bit, and `sample.int(prob = )` and the inverse CDF below give
# different categories from the same stream), the order in which classes and
# occasions are visited, and each verb's policy for a failed replicate.

#' Set a seed for the calling function and restore the random state on its
#' exit
#'
#' Seeded verbs must not leave the caller's random stream changed: a user who
#' sets a seed, fits a model with its own `seed`, and then draws again should
#' get the draw they would have got without the fit. The previous
#' `.Random.seed` is restored when the *calling* function exits, normally or
#' by an error; when there was none, the one `set.seed()` created is removed,
#' so a fresh session stays fresh. Registering the restore on the caller's
#' frame (rather than wrapping its body) lets a call site be one line placed
#' exactly where its own block used to set the seed.
#'
#' The seed is not validated here: each verb validates it with
#' `.multilpa_check_seed()` where it always has, so the refusal comes at the
#' same point as before.
#'
#' @param seed `NULL` (nothing is done) or a seed for [set.seed()].
#' @param frame The function frame whose exit restores the state.
#' @param after Passed to [on.exit()]: whether the restore runs after (the
#'   default) or before the caller's other exit handlers.
#' @return `NULL`, invisibly.
#' @noRd
.latents_local_seed <- function(seed, frame = parent.frame(), after = TRUE) {
  if (is.null(seed)) return(invisible(NULL))
  had_seed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  old_seed <- if (had_seed) get(".Random.seed", envir = globalenv())
  restore <- function() {
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = globalenv())  # nolint: object_name_linter. R's name for the RNG state.
    } else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      rm(".Random.seed", envir = globalenv())
    }
  }
  # `on.exit()` registers on the frame it is evaluated in, so it is called in
  # the caller's frame with a call to `restore` as its expression.
  do.call(base::on.exit, list(as.call(list(restore)), add = TRUE, after = after),
          envir = frame)
  set.seed(seed)
  invisible(NULL)
}

#' Draw one category per row from a matrix of row-wise probabilities
#'
#' Inverse-CDF sampling in one pass, one uniform per row: the category is the
#' first whose cumulative probability reaches the uniform. The last cumulative
#' probability is set to exactly one, so that a rounding shortfall (a row
#' summing to slightly less than one, as exponentiated log probabilities do)
#' cannot leave a uniform above every cumulative total; it falls in the last
#' category, not past it, and not back to the first.
#'
#' A one-column matrix still consumes one uniform per row, so a model with a
#' single class draws the same stream as one with several.
#'
#' @param probabilities Matrix whose rows each sum to one.
#' @return An integer vector, one column index per row.
#' @noRd
.latents_draw_rows <- function(probabilities) {
  stopifnot("`probabilities` must be a matrix" = is.matrix(probabilities),
            "`probabilities` must have at least one column" = ncol(probabilities) >= 1L)
  cumulative <- t(apply(probabilities, 1L, cumsum))
  if (ncol(probabilities) == 1L) cumulative <- matrix(1, nrow(probabilities), 1L)
  cumulative[, ncol(cumulative)] <- 1
  max.col(cumulative >= stats::runif(nrow(probabilities)), ties.method = "first")
}

#' One resample of the groups, as a data frame
#'
#' A group drawn twice has to become two groups, or the refit would pool the two
#' copies into one unit and the resample would carry fewer groups than it drew.
#'
#' @param data The fitting data.
#' @param rows_by_group Row indices of each group, in group order.
#' @param id The group column's name.
#' @param drawn Which groups were drawn, with replacement.
#' @return A data frame with one relabelled group per draw.
#' @noRd
.multilpa_resample_groups <- function(data, rows_by_group, id, drawn) {
  stopifnot(is.data.frame(data), is.list(rows_by_group), is.character(id),
            length(id) == 1L)
  taken <- rows_by_group[drawn]
  resampled <- data[unlist(taken, use.names = FALSE), , drop = FALSE]
  resampled[[id]] <- rep(seq_along(drawn), lengths(taken))
  row.names(resampled) <- NULL
  resampled
}

#' Cluster bootstrap replicates of a fit's estimates
#'
#' Each replicate draws the groups with replacement, refits, and reads the
#' estimates in the original fit's labels; a refit that fails, does not
#' converge or cannot be matched is a row of `NA` with its reason. One
#' replicate's draw, refit and alignment happen together, so replicate `i`
#' sees the same stream whatever the refit consumes.
#'
#' @param data The fitting data.
#' @param group_index The group of each row (integers `1..n_groups`).
#' @param n_groups Number of groups.
#' @param id The group column's name.
#' @param iter How many resamples to draw.
#' @param refit A function of the resampled data returning a fit; it may raise
#'   (the message is recorded) and handles its own warnings.
#' @param extract A function of a converged refit returning its estimates in
#'   the reference labels, or a string saying why it could not.
#' @param reference_names The estimates' names, which also fix their number.
#' @return A list with `estimates` (`iter` rows, `NA` for a failed resample)
#'   and `messages` (why each failed, `NA` otherwise).
#' @noRd
.latents_cluster_bootstrap <- function(data, group_index, n_groups, id, iter,
                                       refit, extract, reference_names) {
  rows_by_group <- split(seq_len(nrow(data)), group_index)
  replicates <- lapply(seq_len(iter), function(i) {
    drawn <- sample.int(n_groups, n_groups, replace = TRUE)
    resampled <- .multilpa_resample_groups(data, rows_by_group, id, drawn)
    failure <- function(message) {
      list(estimate = rep(NA_real_, length(reference_names)), message = message)
    }
    fit <- tryCatch(refit(resampled), error = function(error) conditionMessage(error))
    if (is.character(fit)) return(failure(fit))
    if (!isTRUE(fit$converged)) return(failure("the replicate did not converge"))
    estimate <- extract(fit)
    if (is.character(estimate)) return(failure(estimate))
    list(estimate = estimate, message = NA_character_)
  })
  estimates <- do.call(rbind, lapply(replicates, `[[`, "estimate"))
  colnames(estimates) <- reference_names
  list(estimates = estimates,
       messages = vapply(replicates, `[[`, character(1), "message"))
}

#' The replicate loop of a bootstrap likelihood-ratio test
#'
#' Each replicate simulates from the null and refits both models inside one
#' call of `replicate`, so the draws interleave with the refits as they always
#' have. A replicate that raises is recorded with its message; the warnings it
#' raised are recorded too, and either muffled (the transition and additive
#' verbs, whose refits warn routinely) or passed on to the caller (the profile
#' verb, which has always shown them).
#'
#' @param iter Number of replicates.
#' @param replicate A function of the replicate number returning a list with
#'   `statistic`, `valid`, `boundary`, `logits_settled`, `null_replications`
#'   and `alternative_replications`.
#' @param muffle_warnings Whether a replicate's warnings stop at the record.
#' @return A data frame, one row per replicate.
#' @noRd
.latents_blrt_replicates <- function(iter, replicate, muffle_warnings) {
  do.call(rbind, lapply(seq_len(iter), function(i) {
    warning_text <- character()
    tryCatch(withCallingHandlers({
      outcome <- replicate(i)
      valid <- outcome$valid
      data.frame(replicate = i,
                 statistic = if (valid) max(0, outcome$statistic) else NA_real_,
                 valid = valid, boundary = outcome$boundary,
                 logits_settled = outcome$logits_settled,
                 null_replications = outcome$null_replications,
                 alternative_replications = outcome$alternative_replications,
                 warnings = if (length(warning_text) == 0L) NA_character_ else
                   paste(unique(warning_text), collapse = "; "),
                 error = if (valid) NA_character_ else
                   "Nonconvergence or reversed likelihood")
    }, warning = function(warning) {
      warning_text <<- c(warning_text, conditionMessage(warning))
      if (muffle_warnings) invokeRestart("muffleWarning")
    }), error = function(error) {
      data.frame(replicate = i, statistic = NA_real_, valid = FALSE, boundary = NA,
                 logits_settled = NA, null_replications = NA_integer_,
                 alternative_replications = NA_integer_,
                 warnings = paste(unique(warning_text), collapse = "; "),
                 error = conditionMessage(error))
    })
  }))
}

#' The bootstrap likelihood-ratio p-value from its replicates
#'
#' The p-value counts the observed statistic among the replicates,
#' `(1 + #{T* >= T}) / (iter + 1)`, and is reported only when every replicate
#' is valid: dropping the failed ones would condition the reference
#' distribution on the fits that happened to succeed.
#'
#' @param replicates The replicate table.
#' @param observed The observed statistic.
#' @param iter Number of replicates.
#' @return A list with `p_value`, `monte_carlo_se` and `n_valid`; warns
#'   (`latents_failed_replicates`) when the p-value is withheld.
#' @noRd
.latents_blrt_p_value <- function(replicates, observed, iter) {
  valid <- all(replicates$valid)
  p_value <- if (valid) (1 + sum(replicates$statistic >= observed)) / (iter + 1) else NA_real_
  if (!valid) {
    warning(warningCondition(paste(
      "Some bootstrap fits failed validation, so p_value is NA.",
      "summary() reports every replicate; improve fitting and rerun."),
      class = "latents_failed_replicates"))
  }
  list(p_value = p_value,
       monte_carlo_se = if (valid) sqrt(p_value * (1 - p_value) / (iter + 1)) else NA_real_,
       n_valid = sum(replicates$valid))
}
