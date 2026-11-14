# Where PROC LCA's stopping rule leaves EM. Run from the project root:
#   Rscript equivalence/sas-published/mad-stop.R
#
# PROC LCA stops EM when no gamma or rho moved more than 1e-6 in the last
# iteration (MAD, users' guide 1.3.2, section 3). On a flat likelihood, or one
# whose maximum lies on a bound, EM still has far to go when its steps are
# that small. This script runs a plain EM with that rule from 50 random rho
# starts (PROC LCA's own start: random rhos, equal gammas) for each PROC LCA
# case of the kit and reports, over the starts that reach the top log
# likelihood: how far the stopped estimates lie from latents' maximum, and
# how close the nearest stopped point comes to PROC LCA's printed values. For
# the drinking example it also computes, at that point, the standard errors
# from the inverted Hessian of the log likelihood in the probabilities (PROC
# LCA's errors) and compares them with the printed ones; at latents' maximum
# the same computation must give latents' observed-information errors.
#
# Writes mad-stop.csv beside this script.

suppressMessages(pkgload::load_all(".", quiet = TRUE))
directory <- file.path("equivalence", "sas-published")
source(file.path(directory, "compare-functions.R"))

#' Log likelihood and posteriors of a categorical LCA
#' @param codes Integer matrix, persons by items (codes 1..C_j).
#' @param rho List over items of C_j x K probability matrices.
.lca_e_step <- function(codes, gamma, rho) {
  log_joint <- Reduce(`+`, lapply(seq_along(rho), function(j)
    log(pmax(rho[[j]][codes[, j], , drop = FALSE], 1e-300))), init = 0) +
    matrix(log(gamma), nrow(codes), length(gamma), byrow = TRUE)
  top <- apply(log_joint, 1L, max)
  scaled <- exp(log_joint - top)
  list(log_likelihood = sum(top + log(rowSums(scaled))), posterior = scaled / rowSums(scaled))
}

#' EM stopped by PROC LCA's rule: largest absolute change of a gamma or rho
#' below `criterion`
.mad_em <- function(codes, categories, k, seed, criterion = 1e-6, max_iter = 5000L) {
  set.seed(seed)
  rho <- lapply(categories, function(c) {
    draws <- matrix(stats::runif(c * k), c, k)
    sweep(draws, 2L, colSums(draws), "/")
  })
  gamma <- rep(1 / k, k)
  # EM is iterative: each step needs the previous estimates.
  for (iteration in seq_len(max_iter)) {
    posterior <- .lca_e_step(codes, gamma, rho)$posterior
    gamma_new <- colMeans(posterior)
    rho_new <- lapply(seq_along(categories), function(j) {
      counts <- rowsum(posterior, factor(codes[, j], levels = seq_len(categories[[j]])),
                       reorder = TRUE)
      sweep(counts, 2L, colSums(counts), "/")
    })
    change <- max(abs(c(gamma_new - gamma, unlist(rho_new) - unlist(rho))))
    gamma <- gamma_new
    rho <- rho_new
    if (is.finite(change) && change < criterion) break
  }
  list(gamma = gamma, rho = rho, iterations = iteration,
       log_likelihood = .lca_e_step(codes, gamma, rho)$log_likelihood)
}

#' The printed parameters of a case, as estimates in PROC LCA's class order
.printed_vector <- function(printed) {
  stats::setNames(as.numeric(printed$estimate),
                  ifelse(printed$parameter == "gamma", sprintf("gamma %d", printed$class),
                         sprintf("rho %d %s %d", printed$class, printed$item, printed$category)))
}

#' A stopped point as the same named vector, its classes in the order `perm`
.point_vector <- function(point, items, perm) {
  gamma <- stats::setNames(point$gamma[perm], sprintf("gamma %d", seq_along(perm)))
  rho <- unlist(lapply(seq_along(items), function(j) {
    m <- point$rho[[j]][, perm, drop = FALSE]
    stats::setNames(as.vector(m), sprintf("rho %d %s %d", rep(seq_along(perm), each = nrow(m)),
                                          items[[j]], rep(seq_len(nrow(m)), length(perm))))
  }))
  c(gamma, rho)
}

#' The class order of a point closest to the printed values
.closest <- function(point, items, printed_values) {
  orders <- permutations(length(point$gamma))
  gaps <- vapply(orders, function(perm) {
    v <- .point_vector(point, items, perm)
    max(abs(v[names(printed_values)] - printed_values))
  }, numeric(1))
  list(perm = orders[[which.min(gaps)]], gap = min(gaps))
}

#' Standard errors from the inverted Hessian in the free probabilities
#' (gamma_1..K-1, rho_1..C-1 of each item and class), at a point
#'
#' The Hessian is the central difference of the analytic score, accurate to
#' about 1e-8 here (a second difference of the log likelihood itself would
#' lose its digits to cancellation).
.hessian_errors <- function(codes, point, items, perm, step = 1e-7) {
  k <- length(point$gamma)
  categories <- vapply(point$rho, nrow, integer(1))
  free <- c(point$gamma[-k], unlist(lapply(point$rho, function(m) m[-nrow(m), , drop = FALSE])))
  unpack <- function(x) {
    gamma <- c(x[seq_len(k - 1L)], 1 - sum(x[seq_len(k - 1L)]))
    ends <- (k - 1L) + cumsum((categories - 1L) * k)
    rho <- lapply(seq_along(categories), function(j) {
      m <- matrix(x[(ends[[j]] - (categories[[j]] - 1L) * k) + seq_len((categories[[j]] - 1L) * k)],
                  categories[[j]] - 1L, k)
      rbind(m, 1 - colSums(m))
    })
    list(gamma = gamma, rho = rho)
  }
  score <- function(x) {
    p <- unpack(x)
    posterior <- .lca_e_step(codes, p$gamma, p$rho)$posterior
    gamma_score <- colSums(posterior[, -k, drop = FALSE]) / p$gamma[-k] -
      sum(posterior[, k]) / p$gamma[[k]]
    rho_score <- unlist(lapply(seq_along(categories), function(j) {
      c_j <- categories[[j]]
      counts <- rowsum(posterior, factor(codes[, j], levels = seq_len(c_j)), reorder = TRUE)
      counts[-c_j, , drop = FALSE] / p$rho[[j]][-c_j, , drop = FALSE] -
        matrix(counts[c_j, ] / p$rho[[j]][c_j, ], c_j - 1L, k, byrow = TRUE)
    }))
    c(gamma_score, rho_score)
  }
  hessian <- vapply(seq_along(free), function(i) {
    up <- free
    down <- free
    up[i] <- up[i] + step
    down[i] <- down[i] - step
    (score(up) - score(down)) / (2 * step)
  }, numeric(length(free)))
  covariance <- solve(-(hessian + t(hessian)) / 2)
  # Errors of every gamma and rho, the last of each set by the delta method
  # (its derivative is exact: minus the sum of the others).
  values <- function(x) {
    p <- unpack(x)
    .point_vector(list(gamma = p$gamma, rho = p$rho), items, perm)
  }
  base <- values(free)
  derivative <- vapply(seq_along(free), function(i) {
    up <- free
    up[i] <- up[i] + 1
    values(up) - base
  }, numeric(length(base)))
  stats::setNames(sqrt(rowSums((derivative %*% covariance) * derivative)), names(base))
}

#' One case: codes, categories, latents' fit and the printed blocks
.mad_case <- function(case, persons, k, gamma, rho, with_errors = FALSE) {
  items <- names(persons)
  codes <- vapply(persons, as.integer, integer(nrow(persons)))
  categories <- vapply(persons, nlevels, integer(1))
  printed <- .proc_lca_estimates(gamma, rho)
  printed_values <- .printed_vector(printed)
  fit <- .proc_lca_fit(persons, k)
  latents <- .proc_lca_latents(fit)
  latents_point <- list(gamma = latents$estimate[latents$parameter == "gamma"],
                        rho = lapply(items, function(item) {
                          rows <- latents[latents$parameter == "rho" & latents$item == item, ]
                          matrix(rows$estimate[order(rows$profile, rows$category)],
                                 categories[[item]], k)
                        }))
  runs <- lapply(seq_len(50L), function(seed) .mad_em(codes, categories, k, seed))
  log_likelihoods <- vapply(runs, `[[`, numeric(1), "log_likelihood")
  top <- runs[log_likelihoods > as.numeric(stats::logLik(fit)) - 0.005]
  latents_order <- .closest(latents_point, items, printed_values)
  latents_values <- .point_vector(latents_point, items, latents_order$perm)
  distance <- vapply(top, function(point) {
    order <- .closest(point, items, latents_values)
    max(abs(.point_vector(point, items, order$perm)[names(latents_values)] - latents_values))
  }, numeric(1))
  nearest <- vapply(top, function(point) .closest(point, items, printed_values)$gap, numeric(1))
  best <- top[[which.min(nearest)]]
  summary <- data.frame(
    case = case, starts_at_top = length(top),
    latents_to_printed = latents_order$gap,
    stopped_to_latents_min = min(distance), stopped_to_latents_max = max(distance),
    nearest_stopped_to_printed = min(nearest),
    iterations = paste(range(vapply(top, `[[`, integer(1), "iterations")), collapse = "-"),
    se_gap_latents = NA_real_, se_gap_stopped = NA_real_, se_check_at_maximum = NA_real_)
  if (with_errors) {
    errors <- stats::setNames(as.numeric(printed$standard_error), names(printed_values))
    order <- .closest(best, items, printed_values)
    at_stop <- .hessian_errors(codes, best, items, order$perm)
    latents_errors <- latents$standard_error[match(
      paste(printed$parameter, latents_order$perm[printed$class], printed$item, printed$category),
      paste(latents$parameter, latents$profile, latents$item, latents$category))]
    summary$se_gap_latents <- max(abs(latents_errors - errors))
    summary$se_gap_stopped <- max(abs(at_stop[names(errors)] - errors))
    # The script's own check: at latents' maximum its Hessian errors are
    # latents' observed-information errors.
    at_maximum <- .hessian_errors(codes, latents_point, items, latents_order$perm)
    summary$se_check_at_maximum <- max(abs(at_maximum[names(errors)] - latents_errors))
  }
  summary
}

data <- proc_lca_data(directory)
cases <- .proc_lca_cases[.proc_lca_cases$k > 1L, ]
result <- do.call(rbind, lapply(seq_len(nrow(cases)), function(i) {
  printed <- .proc_lca_printed[[cases$case[[i]]]]
  .mad_case(cases$case[[i]], data[[cases$data[[i]]]], cases$k[[i]], printed$gamma, printed$rho,
            with_errors = cases$case[[i]] == "proclca_drinking_k3")
}))
print(result, row.names = FALSE)
utils::write.csv(result, file.path(directory, "mad-stop.csv"), row.names = FALSE)
