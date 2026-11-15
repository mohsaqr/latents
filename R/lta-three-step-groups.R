# Group classes in the three-step latent Markov model. With an lta() fit of
# H > 1 group classes, step three is a mixture of H hidden Markov chains for
# the modal assignments (the mixed latent Markov model of Vermunt, Tran and
# Magidson 2008, in the bias-adjusted three-step form of Di Mari, Oberski and
# Vermunt 2016): a person belongs to group class h with probability pi_h
# (logits against the last group class) and, given h, follows the step-three
# chain of `.lta_three_step_fit()` with its own initial and transition
# coefficients, intercepts and slopes, as lta(n_group_classes =,
# initial_covariates =, transition_covariates =) has them in one step. The
# classification errors are the pooled step-two matrix, shared by every
# group class because the measurement is. The group classes are latent in
# step three as in step one; their only indicators are the assignments.

#' Step-one group-class posteriors, in the layout's person order
#' @param x The step-one transition fit.
#' @param groups From `.multilpa_prepare_groups()` on the same data.
#' @return A persons-by-group-classes matrix whose rows sum to one.
#' @noRd
.lta_three_step_group_posteriors <- function(x, groups) {
  posteriors <- x$group_posteriors
  rows <- if (!is.null(rownames(posteriors))) match(groups$ids, rownames(posteriors)) else
    seq_len(groups$n)
  stopifnot("step-one group posteriors must cover every person" =
              !anyNA(rows) && nrow(posteriors) == groups$n)
  unname(posteriors[rows, , drop = FALSE])
}

#' The step-one chain of each group class, when the fit reports one
#'
#' A homogeneous first-order fit reports per group class an initial
#' distribution and a transition matrix of the latent profiles; their logits
#' (and the group-class shares') are the natural start of step three, whose
#' chain is the same one seen through the classification errors.
#'
#' @param x The step-one transition fit.
#' @return A list with `shares`, `initial` (classes by profiles) and
#'   `transition` (profiles by profiles by classes), or `NULL`.
#' @noRd
.lta_three_step_step_one_chain <- function(x) {
  k <- x$n_profiles
  h <- x$n_group_classes
  initial <- x$initial_probabilities
  transition <- x$transition_probabilities
  shares <- x$group_probabilities
  usable <- is.matrix(initial) && identical(dim(initial), c(h, k)) &&
    is.array(transition) && identical(dim(transition), c(k, k, h)) &&
    length(shares) == h && all(is.finite(c(initial, transition, shares)))
  if (!usable) return(NULL)
  list(shares = as.numeric(shares), initial = unname(initial),
       transition = unname(transition))
}

#' Intercept-only start of every group class's chain from the step-one chain
#' @param chain From `.lta_three_step_step_one_chain()`.
#' @param p0,p Numbers of initial and transition design columns.
#' @return The standardized-chart start, as `.lta_three_step_group_start()`.
#' @noRd
.lta_three_step_chain_start <- function(chain, p0, p) {
  floor_log <- function(v) log(pmax(v, 1e-4))
  n_classes <- length(chain$shares)
  k_profiles <- ncol(chain$initial)
  per_class <- lapply(seq_len(n_classes), function(h) {
    initial <- matrix(0, p0, k_profiles - 1L)
    initial[1L, ] <- floor_log(chain$initial[h, -k_profiles]) -
      floor_log(chain$initial[h, k_profiles])
    transition <- array(0, c(p, k_profiles - 1L, k_profiles))
    transition[1L, , ] <- vapply(seq_len(k_profiles), function(k) {
      others <- setdiff(seq_len(k_profiles), k)
      floor_log(chain$transition[k, others, h]) - floor_log(chain$transition[k, k, h])
    }, numeric(k_profiles - 1L))
    c(as.vector(initial), as.vector(transition))
  })
  c(floor_log(chain$shares[-n_classes]) - floor_log(chain$shares[n_classes]),
    unlist(per_class))
}

#' Intercept-only starts of each group class's chain from weighted counts
#'
#' The single-class rule of `.lta_three_step_fit()`, with every person's
#' assignments weighted by its posterior for the group class: initial and
#' move counts of the modal profiles, half a count added to every cell,
#' covariate effects at zero.
#'
#' @param inputs From `.lta_three_step_inputs()`.
#' @param weights Persons by group classes.
#' @param p0,p Numbers of initial and transition design columns.
#' @return The standardized-chart start: group logits, then each group
#'   class's initial and transition coefficients.
#' @noRd
.lta_three_step_group_start <- function(inputs, weights, p0, p) {
  layout <- inputs$layout
  k_profiles <- inputs$n_profiles
  n_classes <- ncol(weights)
  first <- inputs$modal[layout$slot[, 1L]]
  moves <- do.call(rbind, lapply(seq_len(layout$n_occasions)[-1L], function(t) {
    observed <- which(!is.na(layout$slot[, t - 1L]) & !is.na(layout$slot[, t]))
    cbind(person = observed, from = inputs$modal[layout$slot[observed, t - 1L]],
          to = inputs$modal[layout$slot[observed, t]])
  }))
  moves <- moves[!is.na(moves[, "from"]) & !is.na(moves[, "to"]), , drop = FALSE]
  shares <- colSums(weights) + 0.5
  per_class <- lapply(seq_len(n_classes), function(h) {
    seen <- !is.na(first)
    first_counts <- vapply(seq_len(k_profiles), function(k) {
      sum(weights[seen, h][first[seen] == k])
    }, numeric(1)) + 0.5
    initial <- matrix(0, p0, k_profiles - 1L)
    initial[1L, ] <- log(first_counts[-k_profiles] / first_counts[k_profiles])
    cell <- (moves[, "from"] - 1L) * k_profiles + moves[, "to"]
    move_counts <- matrix(vapply(seq_len(k_profiles^2), function(index) {
      sum(weights[moves[cell == index, "person"], h])
    }, numeric(1)), k_profiles, byrow = TRUE) + 0.5
    transition <- array(0, c(p, k_profiles - 1L, k_profiles))
    transition[1L, , ] <- vapply(seq_len(k_profiles), function(k) {
      others <- setdiff(seq_len(k_profiles), k)
      log(move_counts[k, others] / move_counts[k, k])
    }, numeric(k_profiles - 1L))
    c(as.vector(initial), as.vector(transition))
  })
  c(log(shares[-n_classes] / shares[n_classes]), unlist(per_class))
}

#' Log likelihood, scores and group posteriors of the mixture of chains
#'
#' Each group class's chain is evaluated by `.lta_three_step_evaluate()`;
#' a person's likelihood is sum_h pi_h L_h. By Fisher's identity the score of
#' group class h's coefficients is its chain score weighted by the person's
#' posterior r_h, and the score of a group logit is r_g - pi_g.
#'
#' @param theta Group logits, then each class's standardized coefficients.
#' @param inputs From `.lta_three_step_inputs()`, with standardized designs.
#' @param n_classes Number of group classes.
#' @return A list with `log_likelihood` (per person), `scores` and
#'   `posteriors` (persons by group classes).
#' @noRd
.lta_three_step_group_evaluate <- function(theta, inputs, n_classes) {
  gamma <- theta[seq_len(n_classes - 1L)]
  chain <- matrix(theta[-seq_len(n_classes - 1L)], ncol = n_classes)
  per_class <- lapply(seq_len(n_classes), function(h) {
    .lta_three_step_evaluate(chain[, h], inputs)
  })
  logits <- c(gamma, 0)
  log_share <- logits - max(logits)
  log_share <- log_share - log(sum(exp(log_share)))
  joint <- sweep(vapply(per_class, `[[`, numeric(length(per_class[[1L]]$log_likelihood)),
                        "log_likelihood"), 2L, log_share, "+")
  if (!is.matrix(joint)) joint <- matrix(joint, ncol = n_classes)
  top <- .multilpa_row_max(joint)
  log_likelihood <- top + log(rowSums(exp(joint - top)))
  posteriors <- exp(joint - log_likelihood)
  group_scores <- sweep(posteriors[, -n_classes, drop = FALSE], 2L,
                        exp(log_share[-n_classes]))
  chain_scores <- do.call(cbind, lapply(seq_len(n_classes), function(h) {
    per_class[[h]]$scores * posteriors[, h]
  }))
  list(log_likelihood = log_likelihood, scores = cbind(group_scores, chain_scores),
       posteriors = posteriors)
}

#' Fit the mixture of step-three chains and its covariance
#'
#' Up to three starts, all from step one: its per-class chains (when the fit
#' reports them), and counts of the modal assignments weighted by its group
#' posteriors and by its modal group classes. The fitted group classes are
#' then matched to the
#' step-one classes by the permutation that maximizes the agreement of the
#' two sets of posteriors, so `group_class_h` means the same in both steps.
#'
#' @param inputs From `.lta_three_step_inputs()`, with `group_posteriors`.
#' @param vcov_type `"robust"` or `"observed"`.
#' @return As `.lta_three_step_fit()`, with `group_class` in the labels and
#'   `posteriors`, the step-three group-class posteriors.
#' @noRd
.lta_three_step_group_fit <- function(inputs, vcov_type) {
  designs <- inputs$designs
  layout <- inputs$layout
  k_profiles <- inputs$n_profiles
  step_one <- inputs$group_posteriors
  n_classes <- ncol(step_one)
  contributing <- inputs$informative_groups %||% rep(TRUE, nrow(layout$slot))
  inside <- lapply(seq_along(designs$transition), function(t) {
    layout$within[, t + 1L] & contributing
  })
  stacked <- do.call(rbind, designs$transition)
  initial_basis <- .lta_three_step_basis(designs$initial, which(contributing))
  transition_basis <- .lta_three_step_basis(stacked, unlist(inside))
  standardized <- inputs
  standardized$designs <- list(
    initial = designs$initial %*% initial_basis,
    transition = lapply(designs$transition, function(d) d %*% transition_basis))
  p0 <- ncol(designs$initial)
  p <- ncol(stacked)
  n_chain <- p0 * (k_profiles - 1L) + p * (k_profiles - 1L) * k_profiles
  # The quasi-Newton search asks for the value and the gradient at the same
  # point: one evaluation serves both.
  cache <- new.env(parent = emptyenv())
  evaluate <- function(theta) {
    if (!identical(cache$theta, theta)) {
      cache$value <- .lta_three_step_group_evaluate(theta, standardized, n_classes)
      cache$theta <- theta
    }
    cache$value
  }
  objective <- function(theta) -sum(evaluate(theta)$log_likelihood[contributing])
  total <- function(theta) colSums(evaluate(theta)$scores[contributing, , drop = FALSE])
  hard <- diag(n_classes)[max.col(step_one, ties.method = "first"), , drop = FALSE]
  starts <- c(
    if (!is.null(inputs$step_one_chain))
      list(.lta_three_step_chain_start(inputs$step_one_chain, p0, p)),
    list(.lta_three_step_group_start(inputs, step_one, p0, p),
         .lta_three_step_group_start(inputs, hard, p0, p)))
  search <- function(theta, iterations) tryCatch(stats::optim(
    theta, objective, function(theta) -total(theta), method = "BFGS",
    control = list(maxit = iterations, reltol = 1e-12)), error = function(e) NULL)
  # A short quasi-Newton search from every start; the best is carried on by
  # Newton, with a longer search first only when Newton cannot finish and
  # the coefficients are not running off to a separated boundary.
  searched <- Filter(function(f) !is.null(f) && is.finite(f$value),
                     lapply(starts, search, iterations = 100L))
  if (!length(searched)) {
    stop(errorCondition(
      "The step-three model with group classes could not be evaluated from any start.",
      class = "latents_no_converge", call = NULL))
  }
  best <- searched[[which.min(vapply(searched, `[[`, numeric(1), "value"))]]
  fitted <- .lta_three_step_newton(best$par, objective, total)
  if (max(abs(total(fitted$par))) > 1e-6 && max(abs(fitted$par)) < 20) {
    longer <- search(fitted$par, 2000L)
    if (!is.null(longer)) fitted <- .lta_three_step_newton(longer$par, objective, total)
  }
  # Match the step-three group classes to step one's.
  posteriors <- evaluate(fitted$par)$posteriors
  orders <- .lta_three_step_permutations(n_classes)
  agreement <- vapply(orders, function(o) {
    sum(step_one[contributing, , drop = FALSE] * posteriors[contributing, o, drop = FALSE])
  }, numeric(1))
  order <- orders[[which.max(agreement)]]
  theta <- .lta_three_step_group_permute(fitted$par, order, n_classes, n_chain)
  evaluated <- evaluate(theta)
  scores <- evaluated$scores[contributing, , drop = FALSE]
  if (max(abs(colSums(scores))) > 1e-3 * max(1, sqrt(nrow(scores)))) {
    stop(errorCondition(sprintf(paste(
      "The step-three transition model with group classes did not converge",
      "(largest score %.3g)."), max(abs(colSums(scores)))),
      class = "latents_no_converge", call = NULL))
  }
  gamma <- theta[seq_len(n_classes - 1L)]
  .r3step_require_interior(matrix(gamma, 1L, n_classes - 1L),
                           matrix(1, sum(contributing), 1L))
  transition_design <- do.call(rbind, Map(function(design, keep) {
    design[keep, , drop = FALSE]
  }, standardized$designs$transition, inside))
  chain <- matrix(theta[-seq_len(n_classes - 1L)], ncol = n_classes)
  n_initial <- p0 * (k_profiles - 1L)
  invisible(lapply(seq_len(n_classes), function(h) {
    .r3step_require_interior(
      matrix(chain[seq_len(n_initial), h], p0, k_profiles - 1L),
      standardized$designs$initial[contributing, , drop = FALSE])
    coefficients <- array(chain[-seq_len(n_initial), h], c(p, k_profiles - 1L, k_profiles))
    lapply(seq_len(k_profiles), function(k) {
      .r3step_require_interior(matrix(coefficients[, , k], p, k_profiles - 1L),
                               transition_design)
    })
  }))
  information <- .inference_score_information(total, theta, 1e-5, scheme = "central",
                                              h_rule = "plus1")
  inverse <- .mixture_invert(information)
  covariance <- if (identical(vcov_type, "robust")) {
    .inference_sandwich(inverse, crossprod(scores), nrow(scores))
  } else inverse
  initial_block <- kronecker(diag(k_profiles - 1L), initial_basis)
  transition_block <- kronecker(diag((k_profiles - 1L) * k_profiles), transition_basis)
  chain_basis <- rbind(
    cbind(initial_block, matrix(0, nrow(initial_block), ncol(transition_block))),
    cbind(matrix(0, nrow(transition_block), ncol(initial_block)), transition_block))
  natural_basis <- as.matrix(.lta_block_diagonal(
    c(list(diag(n_classes - 1L)), rep(list(chain_basis), n_classes))))
  profiles <- paste0("profile_", seq_len(k_profiles))
  classes <- paste0("group_class_", seq_len(n_classes))
  chain_labels <- function(h) {
    initial_labels <- expand.grid(term = colnames(designs$initial),
                                  to = profiles[-k_profiles], stringsAsFactors = FALSE)
    transition_labels <- do.call(rbind, lapply(seq_len(k_profiles), function(k) {
      expand.grid(term = colnames(stacked), to = profiles[-k], from = profiles[k],
                  stringsAsFactors = FALSE)
    }))
    rbind(
      data.frame(block = "initial", group_class = classes[h], from = NA_character_,
                 to = initial_labels$to, term = initial_labels$term,
                 stringsAsFactors = FALSE),
      data.frame(block = "transition", group_class = classes[h],
                 from = transition_labels$from, to = transition_labels$to,
                 term = transition_labels$term, stringsAsFactors = FALSE))
  }
  labels <- rbind(
    data.frame(block = "group", group_class = classes[-n_classes], from = NA_character_,
               to = NA_character_, term = "(Intercept)", stringsAsFactors = FALSE),
    do.call(rbind, lapply(seq_len(n_classes), chain_labels)))
  posteriors <- evaluated$posteriors
  dimnames(posteriors) <- list(NULL, classes)
  list(estimates = as.vector(natural_basis %*% theta),
       covariance = natural_basis %*% covariance %*% t(natural_basis),
       log_likelihood = -objective(theta), labels = labels, posteriors = posteriors)
}

#' Newton-Raphson polish of a quasi-Newton search
#'
#' The mixture likelihood is flat along the group-class directions, where
#' BFGS creeps; a few Newton steps on the numerically differentiated score
#' reach the maximum to the precision the standard errors need. A step that
#' lowers the likelihood is halved; the polish stops when the largest score
#' is below 1e-8 or no step improves.
#'
#' @param theta The search's end point.
#' @param objective Negative log likelihood.
#' @param total Total score (gradient of the log likelihood).
#' @param max_steps Most Newton steps.
#' @return A list with `par` and `value`, as `stats::optim()`.
#' @noRd
.lta_three_step_newton <- function(theta, objective, total, max_steps = 30L) {
  state <- list(par = theta, value = objective(theta), score = total(theta))
  step <- function(state) {
    information <- .inference_score_information(total, state$par, 1e-5, scheme = "central",
                                                h_rule = "plus1")
    direction <- tryCatch(solve(information, state$score), error = function(e) NULL)
    if (is.null(direction) || any(!is.finite(direction))) return(NULL)
    # Halve the step until the likelihood does not fall (at most 14 halvings).
    try_length <- function(length) {
      if (length < 1e-4) return(NULL)
      candidate <- state$par + length * direction
      value <- objective(candidate)
      if (is.finite(value) && value <= state$value + 1e-12) {
        list(par = candidate, value = value)
      } else try_length(length / 2)
    }
    tried <- try_length(1)
    if (is.null(tried)) return(NULL)
    c(tried, list(score = total(tried$par)))
  }
  Reduce(function(current, index) {
    # Done, stuck, or running off to a separated boundary (which the
    # interior check then refuses).
    if (max(abs(current$score)) <= 1e-8 || isTRUE(current$stopped) ||
        max(abs(current$par)) > 50) return(current)
    following <- step(current)
    if (is.null(following)) c(current[c("par", "value", "score")], list(stopped = TRUE)) else
      following
  }, seq_len(max_steps), state)[c("par", "value")]
}

#' Block-diagonal matrix of a list of matrices (base R)
#' @noRd
.lta_block_diagonal <- function(blocks) {
  rows <- vapply(blocks, nrow, integer(1))
  cols <- vapply(blocks, ncol, integer(1))
  out <- matrix(0, sum(rows), sum(cols))
  row_start <- cumsum(c(0L, rows))
  col_start <- cumsum(c(0L, cols))
  # Fills disjoint diagonal blocks of one matrix in place.
  invisible(lapply(seq_along(blocks)[rows > 0L & cols > 0L], function(b) {
    out[row_start[b] + seq_len(rows[b]), col_start[b] + seq_len(cols[b])] <<- blocks[[b]]
  }))
  out
}

#' Every ordering of n group classes
#' @noRd
.lta_three_step_permutations <- function(n) {
  if (n == 1L) return(list(1L))
  do.call(c, lapply(seq_len(n), function(first) {
    lapply(.lta_three_step_permutations(n - 1L), function(rest) {
      c(first, setdiff(seq_len(n), first)[rest])
    })
  }))
}

#' Relabel the group classes of a mixture-of-chains parameter vector
#'
#' New class h is old class `order[h]`: its chain block moves, and the group
#' logits are re-expressed against the new last class.
#' @noRd
.lta_three_step_group_permute <- function(theta, order, n_classes, n_chain) {
  logits <- c(theta[seq_len(n_classes - 1L)], 0)[order]
  chain <- matrix(theta[-seq_len(n_classes - 1L)], n_chain, n_classes)[, order, drop = FALSE]
  c(logits[-n_classes] - logits[n_classes], as.vector(chain))
}
