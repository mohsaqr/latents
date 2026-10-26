# An independent reference for latent class models with jointly modelled
# pairs of categorical items, written from the formula. Nothing here calls a
# latents function or internal (checked in the test with codetools).
#
# Model: P(x_i) = sum_j pi_j prod_v P_j(x_iv), where a pair (a, b) enters as
# one indicator whose categories are the observed combinations of a and b, so
# P_j(a, b) is an unrestricted joint table per class. An indicator missing in
# a row is left out of that row's product; a pair is missing when either of
# its items is (the convention of latents). `ref_marginal_*` implement the
# other convention, Latent GOLD's with `missing includeall`: a row observing
# one item of a pair contributes that item's marginal P_j(a) = sum_b P_j(a, b).

#' Code the indicators after joining the pairs
#'
#' @param data A data frame.
#' @param vars The item names.
#' @param pairs Pairs written "a:b".
#' @return A list: `codes`, an integer matrix with one column per indicator
#'   (pairs joined, named "a:b", in the place of their first item), and
#'   `levels`, the labels each indicator's codes stand for.
ref_code <- function(data, vars, pairs = character()) {
  split <- strsplit(pairs, ":", fixed = TRUE)
  firsts <- vapply(split, `[`, "", 1L)
  seconds <- vapply(split, `[`, "", 2L)
  indicators <- setdiff(vars, seconds)
  columns <- lapply(indicators, function(v) {
    at <- match(v, firsts)
    if (is.na(at)) return(as.character(data[[v]]))
    a <- as.character(data[[v]])
    b <- as.character(data[[seconds[at]]])
    ifelse(is.na(a) | is.na(b), NA_character_, paste(a, b, sep = "\r"))
  })
  names(columns) <- ifelse(indicators %in% firsts,
                           paste(indicators, seconds[match(indicators, firsts)], sep = ":"),
                           indicators)
  levels <- lapply(columns, \(x) sort(unique(x[!is.na(x)])))
  codes <- vapply(seq_along(columns), \(v) match(columns[[v]], levels[[v]]),
                  integer(nrow(data)))
  codes <- matrix(codes, nrow = nrow(data), dimnames = list(NULL, names(columns)))
  # The joined labels are written with ":" for reading; "\r" kept them apart.
  list(codes = codes, levels = lapply(levels, \(l) gsub("\r", ":", l, fixed = TRUE)))
}

#' Log density of every row under every class
#' @param codes From `ref_code()`.
#' @param priors Class proportions.
#' @param probabilities List per indicator of k x C matrices, columns in the
#'   order of `levels`.
#' @return An n x k matrix of log P(x_i, class j).
ref_log_joint <- function(codes, priors, probabilities) {
  k <- length(priors)
  per_item <- lapply(seq_len(ncol(codes)), function(v) {
    observed <- !is.na(codes[, v])
    out <- matrix(0, nrow(codes), k)
    out[observed, ] <- t(log(probabilities[[v]][, codes[observed, v], drop = FALSE]))
    out
  })
  Reduce(`+`, per_item) + matrix(log(priors), nrow(codes), k, byrow = TRUE)
}

ref_log_sum_exp <- function(m) {
  top <- apply(m, 1L, max)
  top + log(rowSums(exp(m - top)))
}

#' The log likelihood: sum_i w_i log sum_j pi_j prod_v P_j(x_iv)
ref_loglik <- function(codes, priors, probabilities, weights = NULL) {
  weights <- weights %||% rep(1, nrow(codes))
  sum(weights * ref_log_sum_exp(ref_log_joint(codes, priors, probabilities)))
}

#' Number of free parameters: k - 1 class logits and C - 1 per indicator and class
ref_n_parameters <- function(levels, k) {
  (k - 1L) + k * sum(lengths(levels) - 1L)
}

#' EM from random starts; the best maximum
#'
#' @param codes,levels From `ref_code()`.
#' @param k Number of classes.
#' @param n_starts Random starts (random posteriors, Dirichlet(1)).
#' @param seed Seed of the starts; the global stream is restored.
#' @param tol Stop when a nonnegative log likelihood gain is less than `tol`
#'   (absolute). Negative numerical steps do not establish convergence.
#' @param weights Optional case weights (used as given).
#' @param start Optional list(priors, probabilities) used as the only start.
#' @return The best complete state visited, including the initial parameters:
#'   log_likelihood, priors, probabilities, posteriors, best_iteration;
#'   iterations/converged and termination_log_likelihood/termination_gain
#'   describe the observed stopping point. `all` holds each start's best value.
ref_em <- function(codes, levels, k, n_starts = 20L, seed = 1L, tol = 1e-12,
                   max_iter = 20000L, weights = NULL, start = NULL) {
  weights <- weights %||% rep(1, nrow(codes))
  if (exists(".Random.seed", globalenv())) {
    saved <- get(".Random.seed", globalenv())
    on.exit(assign(".Random.seed", saved, globalenv()), add = TRUE)
  }
  set.seed(seed)
  m_step <- function(tau) {
    wt <- tau * weights
    priors <- colSums(wt) / sum(wt)
    probabilities <- lapply(seq_len(ncol(codes)), function(v) {
      observed <- !is.na(codes[, v])
      indicator <- outer(codes[observed, v], seq_along(levels[[v]]), `==`) * 1
      counts <- crossprod(wt[observed, , drop = FALSE], indicator)
      counts / rowSums(counts)
    })
    list(priors = priors, probabilities = probabilities)
  }
  run <- function(parameters) {
    history <- -Inf
    iteration <- 1L
    best <- NULL
    # EM is iterative by nature; a bounded loop over its steps.
    repeat {
      joint <- ref_log_joint(codes, parameters$priors, parameters$probabilities)
      marginal <- ref_log_sum_exp(joint)
      ll <- sum(weights * marginal)
      tau <- exp(joint - marginal)
      if (!is.finite(ll) || any(!is.finite(tau))) {
        stop("Independent EM produced non-finite likelihood or posteriors.", call. = FALSE)
      }
      # Finite arithmetic can lower the objective at a stationary point.
      # Retain parameters, objective and posteriors together, including start.
      if (is.null(best) || ll > best$log_likelihood) {
        best <- c(parameters, list(log_likelihood = ll, posteriors = tau,
                                   best_iteration = iteration))
      }
      done <- iteration > 1L && ll >= history && (ll - history) < tol
      if (done || iteration >= max_iter) {
        return(c(best, list(iterations = iteration, converged = done,
                            termination_log_likelihood = ll,
                            termination_gain = ll - history)))
      }
      history <- ll
      parameters <- m_step(tau)
      iteration <- iteration + 1L
    }
  }
  starts <- if (!is.null(start)) list(start) else lapply(seq_len(n_starts), function(s) {
    draws <- matrix(stats::rexp(nrow(codes) * k), nrow(codes), k)
    m_step(draws / rowSums(draws))
  })
  fits <- lapply(starts, run)
  values <- vapply(fits, `[[`, numeric(1), "log_likelihood")
  best <- fits[[which.max(values)]]
  best$all <- values
  best
}

#' Log likelihood under Latent GOLD's partial-pair convention
#'
#' @param data A data frame with the raw items.
#' @param pair c(first, second).
#' @param others The other items.
#' @param priors Class proportions.
#' @param joint A k x C1 x C2 array of the pair's joint probabilities, with
#'   dimnames the categories.
#' @param probabilities Named list of k x C matrices for `others` (columns
#'   named by category).
ref_marginal_loglik <- function(data, pair, others, priors, joint, probabilities) {
  ref_marginal_parts(data, pair, others, priors, joint, probabilities)$log_likelihood
}

ref_marginal_parts <- function(data, pair, others, priors, joint, probabilities) {
  k <- length(priors)
  n <- nrow(data)
  a <- match(as.character(data[[pair[1]]]), dimnames(joint)[[2]])
  b <- match(as.character(data[[pair[2]]]), dimnames(joint)[[3]])
  pair_density <- vapply(seq_len(k), function(j) {
    table <- matrix(joint[j, , ], dim(joint)[2])
    ifelse(!is.na(a) & !is.na(b), table[cbind(pmax(a, 1L), pmax(b, 1L))],
           ifelse(!is.na(a), rowSums(table)[pmax(a, 1L)],
                  ifelse(!is.na(b), colSums(table)[pmax(b, 1L)], 1)))
  }, numeric(n))
  others_density <- Reduce(`+`, lapply(others, function(v) {
    x <- match(as.character(data[[v]]), colnames(probabilities[[v]]))
    out <- matrix(0, n, k)
    out[!is.na(x), ] <- t(log(probabilities[[v]][, x[!is.na(x)], drop = FALSE]))
    out
  }), matrix(0, n, k))
  log_joint <- log(matrix(pair_density, n, k)) + others_density +
    matrix(log(priors), n, k, byrow = TRUE)
  marginal <- ref_log_sum_exp(log_joint)
  list(log_likelihood = sum(marginal), posteriors = exp(log_joint - marginal), a = a, b = b)
}

#' EM under Latent GOLD's partial-pair convention
#'
#' The pair's table runs over every combination of the two items' categories;
#' a row with one item observed spreads its expected count over the other item
#' in proportion to P_j(b | a).
ref_marginal_em <- function(data, pair, others, k, n_starts = 20L, seed = 1L,
                            tol = 1e-12, max_iter = 20000L) {
  if (exists(".Random.seed", globalenv())) {
    saved <- get(".Random.seed", globalenv())
    on.exit(assign(".Random.seed", saved, globalenv()), add = TRUE)
  }
  set.seed(seed)
  n <- nrow(data)
  categories <- lapply(c(pair, others), \(v) sort(unique(as.character(stats::na.omit(data[[v]])))))
  names(categories) <- c(pair, others)
  c1 <- length(categories[[1]])
  c2 <- length(categories[[2]])
  m_step <- function(tau, joint) {
    priors <- colMeans(tau)
    a <- match(as.character(data[[pair[1]]]), categories[[1]])
    b <- match(as.character(data[[pair[2]]]), categories[[2]])
    counts <- array(0, dim = c(k, c1, c2), dimnames = list(NULL, categories[[1]], categories[[2]]))
    invisible(lapply(seq_len(k), function(j) {
      table <- matrix(joint[j, , ], c1)
      both <- !is.na(a) & !is.na(b)
      full <- matrix(0, c1, c2)
      cells <- tapply(tau[both, j], list(factor(a[both], seq_len(c1)),
                                         factor(b[both], seq_len(c2))), sum)
      full <- full + ifelse(is.na(cells), 0, unclass(cells))
      only_a <- !is.na(a) & is.na(b)
      by_a <- tapply(tau[only_a, j], factor(a[only_a], seq_len(c1)), sum)
      full <- full + ifelse(is.na(as.vector(by_a)), 0, as.vector(by_a)) * table / rowSums(table)
      only_b <- is.na(a) & !is.na(b)
      by_b <- tapply(tau[only_b, j], factor(b[only_b], seq_len(c2)), sum)
      full <- full + t(ifelse(is.na(as.vector(by_b)), 0, as.vector(by_b)) * t(table) /
                         colSums(table))
      counts[j, , ] <<- full / sum(full)
    }))
    probabilities <- lapply(stats::setNames(others, others), function(v) {
      x <- match(as.character(data[[v]]), categories[[v]])
      observed <- !is.na(x)
      indicator <- outer(x[observed], seq_along(categories[[v]]), `==`) * 1
      out <- crossprod(tau[observed, , drop = FALSE], indicator)
      out <- out / rowSums(out)
      colnames(out) <- categories[[v]]
      out
    })
    list(priors = priors, joint = counts, probabilities = probabilities)
  }
  run <- function(tau) {
    uniform <- array(1 / (c1 * c2), dim = c(k, c1, c2))
    parameters <- m_step(tau, uniform)
    history <- -Inf
    iteration <- 1L
    # EM is iterative by nature; a bounded loop over its steps.
    repeat {
      parts <- ref_marginal_parts(data, pair, others, parameters$priors, parameters$joint,
                                  parameters$probabilities)
      if (iteration > 1L && (parts$log_likelihood - history) < tol || iteration >= max_iter) {
        return(c(parameters, list(log_likelihood = parts$log_likelihood,
                                  posteriors = parts$posteriors,
                                  converged = iteration < max_iter)))
      }
      history <- parts$log_likelihood
      parameters <- m_step(parts$posteriors, parameters$joint)
      iteration <- iteration + 1L
    }
  }
  fits <- lapply(seq_len(n_starts), function(s) {
    draws <- matrix(stats::rexp(n * k), n, k)
    run(draws / rowSums(draws))
  })
  values <- vapply(fits, `[[`, numeric(1), "log_likelihood")
  best <- fits[[which.max(values)]]
  best$all <- values
  best
}

#' The joint and independence table of one pair from a k x K joint matrix
#'
#' @param joint k x K matrix of the pair's joint probabilities; colnames are
#'   the joined labels "a:b".
#' @return A data frame: profile, category_1, category_2, probability,
#'   independent (product of the two marginals within the class), ratio.
ref_dependence_table <- function(joint) {
  labels <- strsplit(colnames(joint), ":", fixed = TRUE)
  first <- vapply(labels, `[`, "", 1L)
  second <- vapply(labels, `[`, "", 2L)
  do.call(rbind, lapply(seq_len(nrow(joint)), function(j) {
    p <- joint[j, ]
    m1 <- vapply(first, \(c) sum(p[first == c]), numeric(1))
    m2 <- vapply(second, \(c) sum(p[second == c]), numeric(1))
    data.frame(profile = j, category_1 = first, category_2 = second,
               probability = unname(p), independent = unname(m1 * m2),
               ratio = unname(p / (m1 * m2)), stringsAsFactors = FALSE)
  }))
}
