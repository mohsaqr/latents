# Nonparametric bootstrap inference for transition models: resample persons,
# refit the same specification, put each replicate's labels back on the
# original's, and read standard errors and percentile intervals off the
# replicates (Efron & Tibshirani, 1993; Davison & Hinkley, 1997, chapter 3).

#' Percentile bootstrap inference for a transition fit
#'
#' Called by [parameter_inference()] on a `multilpa_transitions` or
#' `multilpa_lta` fit when `method = "bootstrap"`. Persons (the `id` groups)
#' are the resampled units, so dependence across a person's occasions is kept.
#' The table has the Wald path's rows and columns; `standard_error` is the
#' standard deviation of the replicates and `conf_low`/`conf_high` their
#' percentile interval, while `statistic`, `p_value` and `p_adjusted` are `NA`.
#'
#' @param x A transition fit.
#' @param data The data it was fitted to.
#' @param level Interval level.
#' @param iter,n_starts,max_iter,tol Resamples, and the controls of each refit.
#' @param adjust Kept for a frame of the same shape as the Wald path's.
#' @return A base `data.frame` with attributes `covariance`, `level`,
#'   `method`, `iter`, `n_valid`, `replicates` and `messages`.
#' @noRd
.lta_bootstrap_inference <- function(x, data, level, iter, n_starts, max_iter, tol,
                                     adjust) {
  if (is.null(data)) {
    stop(errorCondition(paste(
      "The bootstrap for a transition model needs `data`: each resample keeps",
      "its persons' occasions and covariates."),
      class = "latents_bad_argument", call = NULL))
  }
  .lta_check_data(x, data)
  if (!isTRUE(x$converged)) {
    stop(errorCondition(
      "The original fit did not converge; bootstrap it only once it has.",
      class = "latents_no_converge", call = NULL))
  }
  table <- .lta_bootstrap_frame(x)
  drawn <- .lta_bootstrap_replicates(x, data, iter, n_starts, max_iter, tol,
                                     table$name)
  valid <- stats::complete.cases(drawn$estimates)
  if (sum(valid) < 2L) {
    reasons <- drawn$messages[!is.na(drawn$messages)]
    stop(errorCondition(sprintf(
      "Only %d of %d resamples produced a usable fit; the first reason was: %s",
      sum(valid), iter,
      if (length(reasons) > 0L) reasons[1L] else "no reason was recorded"),
      class = "latents_bootstrap_failed", call = NULL))
  }
  if (sum(valid) < iter) {
    warning(warningCondition(sprintf(
      "%d of %d resamples did not produce a usable fit and were dropped.",
      iter - sum(valid), iter), class = "latents_bootstrap_dropped", call = NULL))
  }
  kept <- drawn$estimates[valid, , drop = FALSE]
  quantiles <- apply(kept, 2L, stats::quantile,
                     probs = c((1 - level) / 2, (1 + level) / 2), names = FALSE)
  table$standard_error <- apply(kept, 2L, stats::sd)
  # The interval is the inference: a Wald statistic would put back the normal
  # approximation the percentile interval avoids.
  table$statistic <- NA_real_
  table$p_value <- NA_real_
  table$p_adjusted <- NA_real_
  table$conf_low <- quantiles[1L, ]
  table$conf_high <- quantiles[2L, ]
  covariance <- stats::cov(kept)
  dimnames(covariance) <- list(table$name, table$name)
  columns <- setdiff(names(table), "name")
  result <- table[, columns, drop = FALSE]
  row.names(result) <- NULL
  structure(result, covariance = covariance, covariance_unconstrained = NULL,
            level = level, method = "bootstrap", iter = iter, n_valid = sum(valid),
            replicates = kept, messages = drawn$messages, adjust = adjust)
}

#' The rows a transition bootstrap reports, with their estimates
#'
#' The same rows, labels and column order as the Wald table of the fit's
#' class, plus `name`, the coefficient each row is matched on.
#'
#' @param x A transition fit.
#' @return A base `data.frame`.
#' @noRd
.lta_bootstrap_frame <- function(x) {
  estimates <- .lta_bootstrap_estimates(x)
  if (inherits(x, "multilpa_transitions")) {
    labels <- .multilpa_transition_labels(names(estimates), .multilpa_transition_view(x))
    return(data.frame(labels, estimate = unname(estimates), standard_error = NA_real_,
                      statistic = NA_real_, p_value = NA_real_, p_adjusted = NA_real_,
                      conf_low = NA_real_, conf_high = NA_real_, name = names(estimates),
                      row.names = NULL, stringsAsFactors = FALSE))
  }
  pieces <- strsplit(names(estimates), ".", fixed = TRUE)
  block <- vapply(pieces, `[`, character(1), 1L)
  moves <- strsplit(vapply(pieces, `[`, character(1), 3L), "->", fixed = TRUE)
  second <- block == "transition2"
  frame <- data.frame(
    block = block,
    group_class = vapply(pieces, `[`, character(1), 2L),
    # An initial logit names one profile, a move names its origin before it.
    from = vapply(moves, function(m) if (length(m) < 2L) NA_character_ else
      m[length(m) - 1L], character(1)),
    to = vapply(moves, function(m) m[length(m)], character(1)),
    term = vapply(pieces, function(p) paste(p[-(1:3)], collapse = "."), character(1)),
    estimate = unname(estimates), standard_error = NA_real_, statistic = NA_real_,
    p_value = NA_real_, p_adjusted = NA_real_, conf_low = NA_real_,
    conf_high = NA_real_, name = names(estimates), stringsAsFactors = FALSE)
  if (any(second)) {
    frame$previous <- ifelse(second, vapply(moves, `[`, character(1), 1L), NA_character_)
  }
  # The Wald table's block order: transitions, initial, then second order.
  frame[order(match(frame$block, c("transition", "initial", "transition2")),
              seq_len(nrow(frame))), , drop = FALSE]
}

#' The estimates a transition bootstrap reports, named
#'
#' A homogeneous fit reports its natural parameters (measurement, group shares,
#' initial and transition probabilities); a general fit its initial,
#' transition and second-order logits.
#'
#' @param x A transition fit.
#' @return A named numeric vector.
#' @noRd
.lta_bootstrap_estimates <- function(x) {
  if (inherits(x, "multilpa_transitions")) return(.multilpa_transition_encode(x, "natural"))
  theta <- .lta_pack(.lta_parameters(x), x$variance_model, x$extra_data)
  names(theta) <- .lta_names(x)
  theta[sub("\\..*$", "", names(theta)) %in% c("transition", "initial", "transition2")]
}

#' Bootstrap replicates of a transition fit's reported estimates
#'
#' @param x A transition fit.
#' @param data The data it was fitted to.
#' @param iter,n_starts,max_iter,tol Resamples, and the controls of each refit.
#' @param reference_names The coefficient names each replicate is matched to.
#' @return A list with `estimates` (`iter` rows, `NA` for a failed resample)
#'   and `messages` (why each failed, `NA` otherwise).
#' @noRd
.lta_bootstrap_replicates <- function(x, data, iter, n_starts, max_iter, tol,
                                      reference_names) {
  .latents_cluster_bootstrap(
    data, x$group_index, x$n_groups, x$id, iter,
    # A replicate's own warnings (an unconverged start, a boundary) are judged
    # from the fit itself, so they are not repeated once per resample.
    refit = function(resampled) {
      withCallingHandlers(
        .lta_refit(x, resampled, n_starts, max_iter, tol),
        warning = function(w) invokeRestart("muffleWarning"),
        message = function(m) invokeRestart("muffleMessage"))
    },
    extract = function(fit) {
      aligned <- .lta_align_estimates(fit, x)
      if (!setequal(names(aligned), reference_names)) {
        return("the replicate's coefficients could not be matched")
      }
      unname(aligned[reference_names])
    },
    reference_names = reference_names)
}

#' A replicate's reported estimates, in the original fit's labels
#'
#' Profiles are matched on the measurement, then group classes on what the
#' matched profiles make of them.
#'
#' @param replicate A transition fit from a resample.
#' @param reference The original fit.
#' @return The replicate's estimates, named as the original's.
#' @noRd
.lta_align_estimates <- function(replicate, reference) {
  profile_order <- .multilpa_match_order(.lta_profile_signature(reference, reference),
                                         .lta_profile_signature(replicate, reference))
  if (inherits(reference, "multilpa_transitions")) {
    permuted <- .lta_permute_transitions(replicate, profile_order,
                                         seq_len(replicate$n_group_classes))
    class_order <- .multilpa_match_order(.lta_class_signature(reference),
                                         .lta_class_signature(permuted))
    return(.multilpa_transition_encode(
      .lta_permute_transitions(permuted, seq_len(permuted$n_profiles), class_order),
      "natural"))
  }
  relabelled <- .lta_relabel_profiles(.lta_bootstrap_estimates(replicate), profile_order)
  classes <- names(reference$group_probabilities)
  class_order <- .multilpa_match_order(
    .lta_class_signature(reference, .lta_bootstrap_estimates(reference)),
    .lta_class_signature(replicate, relabelled))
  .lta_relabel_classes(relabelled, class_order, classes)
}

#' What a transition fit's profiles look like, for matching
#'
#' The measurement signature of every block, in the units
#' `.latents_measurement_signature()` uses, with no mixing column: a
#' transition fit's profile shares change from occasion to occasion.
#'
#' @param x A transition fit.
#' @param reference The fit whose scale both signatures are measured in.
#' @return A profiles-by-features matrix.
#' @noRd
.lta_profile_signature <- function(x, reference) {
  blocks_of <- function(fit) {
    blocks <- if (inherits(fit, "multilpa_lta")) fit$measurement else
      list(list(means = fit$means, variances = fit$variances,
                response_probabilities = fit$response_probabilities))
    lapply(blocks, function(block) c(block, list(n_profiles = fit$n_profiles)))
  }
  signature <- .latents_measurement_signature(blocks_of(x), blocks_of(reference))
  stopifnot("a transition fit must have measurement to match profiles on" =
              !is.null(signature))
  signature
}

#' What a transition fit's group classes look like, for matching
#'
#' A class's share and its initial and transition parameters, once the
#' profiles agree. A stayer class can only match a stayer class.
#'
#' @param x A transition fit, its profiles already matched.
#' @param estimates For a general fit, its logits with matched profile labels.
#' @return A classes-by-features matrix.
#' @noRd
.lta_class_signature <- function(x, estimates = NULL) {
  shares <- unname(x$group_probabilities)
  if (inherits(x, "multilpa_transitions")) {
    n_classes <- x$n_group_classes
    return(cbind(shares, unname(x$initial_probabilities),
                 t(vapply(seq_len(n_classes), function(h) {
                   as.vector(x$transition_probabilities[, , h])
                 }, numeric(x$n_profiles^2)))))
  }
  classes <- names(x$group_probabilities)
  pieces <- strsplit(names(estimates), ".", fixed = TRUE)
  owner <- vapply(pieces, `[`, character(1), 2L)
  key <- vapply(pieces, function(p) paste(p[-2L], collapse = "."), character(1))
  keys <- sort(unique(key))
  values <- t(vapply(classes, function(h) {
    row <- numeric(length(keys))
    mine <- owner == h
    row[match(key[mine], keys)] <- estimates[mine]
    row
  }, numeric(length(keys))))
  stayer <- as.numeric(.lta_parameters(x)$stayer %||% rep(FALSE, length(classes)))
  # A stayer class has no transition parameters to compare; the flag, scaled
  # far beyond any logit, keeps it from being matched to a mover class.
  cbind(shares, values, 1e6 * stayer)
}

#' Rename profiles in general-fit logits, rebasing the initial logits
#'
#' Transitions and second-order transitions are logits against staying, which
#' a relabelling maps to staying, so their names move and their values do not.
#' Initial logits are against the last profile: when another profile becomes
#' last, every initial logit is shifted by the new reference's.
#'
#' @param estimates Named logits of one fit.
#' @param order Integer permutation: new profile `j` is old profile `order[j]`.
#' @return The logits, named in the new labels.
#' @noRd
.lta_relabel_profiles <- function(estimates, order) {
  n_profiles <- length(order)
  new_of_old <- match(seq_len(n_profiles), order)
  pieces <- strsplit(names(estimates), ".", fixed = TRUE)
  block <- vapply(pieces, `[`, character(1), 1L)
  rename <- function(profiles) {
    numbers <- as.integer(sub("^profile_", "", strsplit(profiles, "->", fixed = TRUE)[[1L]]))
    paste(paste0("profile_", new_of_old[numbers]), collapse = "->")
  }
  moved <- vapply(seq_along(pieces), function(i) {
    p <- pieces[[i]]
    paste(c(p[1L], p[2L], rename(p[3L]), p[-(1:3)]), collapse = ".")
  }, character(1))
  result <- stats::setNames(unname(estimates), moved)
  initial <- block == "initial"
  if (!any(initial)) return(result)
  # One baseline-category logit vector per group class and term.
  owner <- vapply(pieces[initial], function(p) paste(c(p[2L], p[-(1:3)]), collapse = "\r"),
                  character(1))
  old_profile <- as.integer(sub("^profile_", "", vapply(pieces[initial], `[`,
                                                        character(1), 3L)))
  rebased <- unlist(unname(lapply(split(seq_along(owner), owner), function(rows) {
    eta <- numeric(n_profiles)
    eta[old_profile[rows]] <- estimates[initial][rows]
    shifted <- .latents_rebase_logits(matrix(eta, 1L), order)[1L, ]
    p <- pieces[initial][[rows[1L]]]
    stats::setNames(shifted[-n_profiles],
                    paste(p[1L], p[2L], paste0("profile_", seq_len(n_profiles - 1L)),
                          paste(p[-(1:3)], collapse = "."), sep = "."))
  })))
  c(result[!initial], rebased)
}

#' Rename group classes in general-fit logits
#' @param estimates Named logits, profiles already matched.
#' @param order Integer permutation: new class `i` is old class `order[i]`.
#' @param classes The class names, in position order.
#' @return The logits, named in the new class labels.
#' @noRd
.lta_relabel_classes <- function(estimates, order, classes) {
  new_name <- stats::setNames(classes[match(seq_along(classes), order)], classes)
  pieces <- strsplit(names(estimates), ".", fixed = TRUE)
  names(estimates) <- vapply(pieces, function(p) {
    paste(c(p[1L], new_name[[p[2L]]], p[-(1:2)]), collapse = ".")
  }, character(1))
  estimates
}

#' Reorder a homogeneous transition fit's profiles and group classes
#' @param x A `multilpa_transitions` fit.
#' @param profiles,classes Integer permutations: new `j` is old `order[j]`.
#' @return The fit with every profile- and class-indexed parameter reordered.
#' @noRd
.lta_permute_transitions <- function(x, profiles, classes) {
  stopifnot(inherits(x, "multilpa_transitions"),
            identical(sort(as.integer(profiles)), seq_len(x$n_profiles)),
            identical(sort(as.integer(classes)), seq_len(x$n_group_classes)))
  keep_names <- function(new, old) {
    dimnames(new) <- dimnames(old)
    new
  }
  x$means <- keep_names(x$means[profiles, , drop = FALSE], x$means)
  x$variances <- keep_names(x$variances[profiles, , drop = FALSE], x$variances)
  if (!is.null(x$covariances)) {
    x$covariances <- keep_names(x$covariances[, , profiles, drop = FALSE], x$covariances)
  }
  if (length(x$response_probabilities) > 0L) {
    x$response_probabilities <- lapply(x$response_probabilities, function(block) {
      keep_names(block[profiles, , drop = FALSE], block)
    })
  }
  x$initial_probabilities <- keep_names(
    x$initial_probabilities[classes, profiles, drop = FALSE], x$initial_probabilities)
  x$transition_probabilities <- keep_names(
    x$transition_probabilities[profiles, profiles, classes, drop = FALSE],
    x$transition_probabilities)
  x$group_probabilities <- stats::setNames(x$group_probabilities[classes],
                                           names(x$group_probabilities))
  x
}

#' Validate bootstrap controls, scope the seed, and run the transition bootstrap
#' @param x A transition fit.
#' @param data,level,iter,n_starts,max_iter,tol,seed,adjust As in
#'   [parameter_inference()].
#' @return The bootstrap inference table.
#' @noRd
.lta_bootstrap_dispatch <- function(x, data, level, iter, n_starts, max_iter, tol,
                                    seed, adjust) {
  whole <- function(value, minimum) is.numeric(value) && length(value) == 1L &&
    is.finite(value) && value >= minimum && value == floor(value)
  stopifnot(
    "`iter` must be a single integer of at least two" = whole(iter, 2),
    "`n_starts` must be a single positive integer" = whole(n_starts, 1),
    "`max_iter` must be a single positive integer" = whole(max_iter, 1),
    "`tol` must be a single positive number" =
      is.numeric(tol) && length(tol) == 1L && is.finite(tol) && tol > 0)
  .multilpa_check_seed(seed)
  .latents_local_seed(seed, after = FALSE)
  .lta_bootstrap_inference(x, data, level, as.integer(iter), as.integer(n_starts),
                           as.integer(max_iter), tol, adjust)
}
