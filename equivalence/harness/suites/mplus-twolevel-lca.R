# Mplus User's Guide example 10.7: two-level latent class analysis with
# categorical latent class indicators and a between-level categorical latent
# variable. This is the nonparametric two-level mixture of Vermunt (2003) in
# Mplus's own notation -- cw(4) profiles measured by ten three-category items,
# cb(5) group classes that shift profile prevalence through `cw#1-cw#3 ON cb`,
# and a measurement model shared across group classes.
#
# Artifacts in equivalence/mplus/twolevel-lca/ are the genuine Mplus 8.8 files
# from https://www.statmodel.com/usersguide/chap10/chapter10.zip
# (sha256 ae29930be639669311d8652af790c531a861d8aa5414783766cd175f7c32b1d0).
# Reference values are parsed from ex10.7.out rather than transcribed.
#
# What Mplus retained, and where each piece is compared
# ------------------------------------------------------
# ex10.7.out holds: the 80 thresholds and the 19 latent-class logits, each with
# its MLR standard error; the response probabilities ("RESULTS IN PROBABILITY
# SCALE", no standard errors); for CW and CB the class counts and proportions
# based on the estimated posteriors and on the most likely class, the entropy,
# the average posteriors by most likely class (AvePP) and the classification
# probabilities; the log-likelihood, AIC, BIC, sample-size adjusted BIC and the
# MLR scaling correction; TECH1 (parameter numbering and starting values, which
# are inputs, not results) and TECH8 (the iteration history). No SAVEDATA file
# was retained, so Mplus's own row posteriors are not available.
#
# Mplus prints three outputs this package has no counterpart for, so they are
# not compared: the counts of the joint (CB, CW) patterns, the overall entropy
# of those patterns, and the logits of the classification probabilities.
#
# Mplus did not run to its maximum. Its EM stops on an absolute change of 1e-3
# in the log-likelihood, and TECH8 shows the winning start (set 24) still
# gaining 0.00096 per iteration, at a ratio of 0.84 per step, when it stopped.
# Climbing from Mplus's published estimates gains a further 0.0065 and moves
# parameters by up to 0.022 -- forty times their printed rounding -- and moves
# the CB posterior counts by 29 observations. Mplus's printed results are
# therefore compared where they were computed, at Mplus's estimates, and the
# package's converged fit is compared with Mplus's estimates carried on to the
# maximum. Two cases:
#
#   1. At Mplus's published estimates. The package is evaluated there
#      (`max_iter = 0`, started at the printed thresholds and logits). Every
#      output it reports is compared, to machine precision, with a reference
#      computed at the same parameters by a two-level likelihood written here in
#      base R. Every result Mplus printed is compared with the same package fit
#      at Mplus's printed precision plus the most that rounding the estimates to
#      three decimals can move it (see `.lca2_rounding_allowance()`), because
#      Mplus computed its results from unrounded estimates.
#   2. Converged. The package's converged fit is compared with the reference
#      carried from Mplus's published estimates to the maximum by Newton steps
#      on the from-scratch likelihood, standard errors included, which come from
#      the from-scratch observed information. The claims this suite made before
#      the ledger existed are kept here as extra rows.

# The rows' profile entropies within their group class, which the
# complete-data entropy of ICL, AWE and CLC needs, come from the shared
# two-level reference (`two_level_row_conditional_entropy()`).
source(file.path("equivalence", "harness", "two-level-lca-reference.R"), local = environment())

#' Read the Mplus output shipped with example 10.7
#' @param path Path to the gzipped Mplus output file.
#' @return Character vector of output lines.
.mplus_output_lines <- function(path) {
  connection <- gzfile(path, open = "rt")
  on.exit(close(connection), add = TRUE)
  readLines(connection)
}

#' Pull a labelled scalar from an Mplus fit-statistics block
#' @param lines Mplus output lines.
#' @param label Text preceding the value.
#' @return The single numeric value reported on that line.
.mplus_scalar <- function(lines, label) {
  matched <- grep(label, lines, fixed = TRUE, value = TRUE)
  if (length(matched) < 1L) stop(sprintf("Mplus output has no line matching '%s'.", label))
  value <- utils::tail(strsplit(trimws(matched[[1L]]), " +")[[1L]], 1L)
  as.numeric(value)
}

#' Parse the RESULTS IN PROBABILITY SCALE block
#'
#' @param lines Mplus output lines.
#' @param items Character vector of item names, upper case as Mplus prints them.
#' @param n_categories Number of categories per item.
#' @return A named list of `pattern` (a two-column matrix of cb and cw indices)
#'   and `probabilities` (a patterns-by-item-by-category array).
.mplus_probability_scale <- function(lines, items, n_categories) {
  start <- grep("^RESULTS IN PROBABILITY SCALE", lines)
  if (length(start) != 1L) stop("Expected exactly one probability-scale block.")
  block <- lines[seq(start, length(lines))]
  headers <- grep("^Latent Class Pattern ", block)
  pattern <- do.call(rbind, lapply(headers, function(index) {
    as.integer(strsplit(trimws(sub("^Latent Class Pattern ", "", block[[index]])), " +")[[1L]])
  }))
  probabilities <- vapply(headers, function(index) {
    window <- block[seq(index, min(index + 6L * length(items), length(block)))]
    vapply(items, function(item) {
      at <- grep(sprintf("^ %s$", item), window)[[1L]]
      rows <- window[at + seq_len(n_categories)]
      as.numeric(sub("^ +Category [0-9]+ +", "", rows))
    }, numeric(n_categories))
  }, matrix(0, n_categories, length(items)))
  list(pattern = pattern,
       probabilities = aperm(probabilities, c(3L, 2L, 1L)))
}

#' Parse the categorical latent variable logits into prevalence probabilities
#' @param lines Mplus output lines.
#' @param n_profiles Number of within-level classes.
#' @param n_group_classes Number of between-level classes.
#' @return A list with `group_probabilities` and `profile_probabilities`.
.mplus_prevalence <- function(lines, n_profiles, n_group_classes) {
  start <- grep("^Categorical Latent Variables", lines)[[1L]]
  block <- lines[seq(start, length(lines))]
  value_of <- function(name, within) {
    hits <- grep(sprintf("^ +%s +-?[0-9]", name), block)
    if (length(hits) < 1L) stop(sprintf("no estimate printed for %s", name))
    chosen <- if (within) hits[[1L]] else hits[[length(hits)]]
    as.numeric(strsplit(trimws(block[[chosen]]), " +")[[1L]][[2L]])
  }
  intercepts <- vapply(seq_len(n_profiles - 1L),
                       \(s) value_of(sprintf("CW#%d", s), within = TRUE), numeric(1))
  slope_start <- grep("^ CW#1 +ON", block)[[1L]]
  slopes <- vapply(seq_len(n_profiles - 1L), function(s) {
    at <- grep(sprintf("^ CW#%d +ON", s), block)[[1L]]
    window <- block[at + seq_len(n_group_classes)]
    vapply(seq_len(n_group_classes - 1L), function(j) {
      as.numeric(strsplit(trimws(window[[j]]), " +")[[1L]][[2L]])
    }, numeric(1))
  }, numeric(n_group_classes - 1L))
  means_at <- grep("^ Means", block)
  means_at <- means_at[means_at > slope_start][[1L]]
  group_logits <- vapply(seq_len(n_group_classes - 1L), function(j) {
    as.numeric(strsplit(trimws(block[[means_at + j]]), " +")[[1L]][[2L]])
  }, numeric(1))
  softmax <- function(logits) exp(logits) / sum(exp(logits))
  profile_probabilities <- t(vapply(seq_len(n_group_classes), function(j) {
    dummy <- as.numeric(seq_len(n_group_classes - 1L) == j)
    softmax(c(intercepts + as.numeric(crossprod(slopes, dummy)), 0))
  }, numeric(n_profiles)))
  list(group_probabilities = softmax(c(group_logits, 0)),
       profile_probabilities = profile_probabilities)
}

#' Renormalize probabilities that a program printed to fixed precision
#'
#' Published output rounds each cell, so rows sum to one only up to the printed
#' precision. This checks the discrepancy is no larger than that rounding can
#' explain, then renormalizes, so that rounded published values can be fed back
#' in as an exact parameter set.
#'
#' @param block A matrix of probabilities with rows summing to one before rounding.
#' @param precision The printed resolution, for example 1e-3 for three decimals.
#' @return The matrix with rows renormalized to sum to one.
.renormalize_printed <- function(block, precision) {
  drift <- abs(rowSums(block) - 1)
  if (any(drift > ncol(block) * precision / 2)) {
    stop("Row sums depart from one by more than the printed rounding can explain.")
  }
  block / rowSums(block)
}

# ---------------------------------------------------------------------------
# More of the Mplus output: the estimates with their standard errors, and the
# classification results of each level
# ---------------------------------------------------------------------------

#' Parse the thresholds of MODEL RESULTS, with their standard errors
#' @param lines Mplus output lines.
#' @return A `data.frame`, one row per printed threshold: `group_class`,
#'   `profile`, `item`, `cut` (Mplus's `$1`, `$2`), `estimate`,
#'   `standard_error`.
.mplus_thresholds <- function(lines) {
  from <- grep("^MODEL RESULTS$", lines)
  to <- grep("^QUALITY OF NUMERICAL RESULTS", lines)
  stopifnot("Mplus output has one MODEL RESULTS block" =
              length(from) == 1L && length(to) == 1L && to > from)
  block <- lines[seq(from, to)]
  headers <- grep("^Latent Class Pattern ", block)
  ends <- c(headers[-1L], length(block))
  do.call(rbind, lapply(seq_along(headers), function(k) {
    pattern <- as.integer(strsplit(trimws(sub("^Latent Class Pattern ", "",
                                               block[[headers[[k]]]])), " +")[[1L]])
    rows <- grep("^ +U[0-9]+\\$[0-9]+ ", block[seq(headers[[k]], ends[[k]])], value = TRUE)
    fields <- strsplit(trimws(rows), " +")
    labels <- vapply(fields, `[[`, character(1), 1L)
    data.frame(group_class = pattern[[1L]], profile = pattern[[2L]],
               item = as.integer(sub("^U([0-9]+)\\$.*$", "\\1", labels)),
               cut = as.integer(sub("^.*\\$", "", labels)),
               estimate = as.numeric(vapply(fields, `[[`, character(1), 2L)),
               standard_error = as.numeric(vapply(fields, `[[`, character(1), 3L)))
  }))
}

#' Parse the latent-class logits of MODEL RESULTS, with their standard errors
#'
#' In Mplus's parameter order: the CW intercepts, the `CW#s ON CB#b` slopes
#' (CW by CW, CB fastest), the CB means. CW's last class and CB's last class
#' are the reference categories.
#'
#' @param lines Mplus output lines.
#' @param n_profiles,n_group_classes Numbers of CW and CB classes.
#' @return A `data.frame` with `label`, `estimate`, `standard_error`.
.mplus_class_logits <- function(lines, n_profiles, n_group_classes) {
  start <- grep("^Categorical Latent Variables", lines)
  end <- grep("^QUALITY OF NUMERICAL RESULTS", lines)
  stopifnot("Mplus output has one block of latent-class logits" =
              length(start) == 1L && length(end) == 1L)
  block <- lines[seq(start, end)]
  read_after <- function(anchor, n, prefix = "") {
    at <- grep(anchor, block)
    stopifnot("each logit block is printed once" = length(at) == 1L)
    fields <- strsplit(trimws(block[at + seq_len(n)]), " +")
    data.frame(label = paste0(prefix, vapply(fields, `[[`, character(1), 1L)),
               estimate = as.numeric(vapply(fields, `[[`, character(1), 2L)),
               standard_error = as.numeric(vapply(fields, `[[`, character(1), 3L)))
  }
  out <- rbind(
    read_after("^ Intercepts", n_profiles - 1L),
    do.call(rbind, lapply(seq_len(n_profiles - 1L), function(s) {
      read_after(sprintf("^ CW#%d +ON", s), n_group_classes - 1L, sprintf("CW#%d ON ", s))
    })),
    read_after("^ Means", n_group_classes - 1L))
  expected <- c(sprintf("CW#%d", seq_len(n_profiles - 1L)),
                sprintf("CW#%d ON CB#%d", rep(seq_len(n_profiles - 1L), each = n_group_classes - 1L),
                        rep(seq_len(n_group_classes - 1L), n_profiles - 1L)),
                sprintf("CB#%d", seq_len(n_group_classes - 1L)))
  stopifnot("the logits are in Mplus's parameter order" = identical(out$label, expected))
  out
}

#' Counts and proportions of one latent class variable
#'
#' Mplus prints these per variable twice, based on the estimated posteriors and
#' on the most likely class. It counts both variables over observations, so a
#' CB count is the sum of each cluster's posterior times its size.
#'
#' @param lines Mplus output lines.
#' @param basis `"ESTIMATED POSTERIOR"` or `"MOST LIKELY"`.
#' @param variable `"CW"` or `"CB"`.
#' @param k Number of classes of the variable.
#' @return A list with `count` and `proportion`, in Mplus's class labels.
.mplus_variable_classes <- function(lines, basis, variable, k) {
  headers <- grep("^FINAL CLASS COUNTS AND PROPORTIONS FOR EACH LATENT CLASS VARIABLE", lines)
  at <- headers[grepl(basis, lines[headers + 1L], fixed = TRUE)]
  stopifnot("one per-variable count block on this basis" = length(at) == 1L)
  first <- grep(sprintf("^ +%s +1 +[0-9.]+ +[0-9.]+ *$", variable), lines)
  first <- first[first > at][[1L]]
  fields <- strsplit(trimws(lines[first + seq_len(k) - 1L]), " +")
  values <- vapply(fields, function(f) as.numeric(utils::tail(f, 2L)), numeric(2))
  list(count = values[1L, ], proportion = values[2L, ])
}

#' A k-by-k matrix printed under the first matching header after line `from`
#'
#' For the classification probabilities, which the AvePP parser of
#' `helper-mplus-classification.R` does not read (it looks for its own header).
#' @param lines Mplus output lines.
#' @param header Pattern of the header line.
#' @param from Line index to search after.
#' @param k Number of classes.
#' @return The matrix, in Mplus's labels.
.mplus_square_after <- function(lines, header, from, k) {
  at <- grep(header, lines)
  at <- at[at > from][1L]
  stopifnot("the matrix is printed after the anchor" = !is.na(at))
  rows <- grep("^\\s+[0-9]+(\\s+-?[0-9.]+)+\\s*$", lines)
  rows <- rows[rows > at][seq_len(k + 1L)][-1L]
  values <- do.call(rbind, lapply(strsplit(trimws(lines[rows]), "\\s+"),
                                  function(fields) as.numeric(fields)[-1L]))
  stopifnot("the matrix is k by k" = identical(dim(values), as.integer(c(k, k))))
  values
}

#' Mplus's classification results of each level, in Mplus's labels
#' @param lines Mplus output lines.
#' @param n_profiles,n_group_classes Numbers of CW and CB classes.
#' @return A list with `individuals` (CW) and `groups` (CB), each holding
#'   `count`, `proportion`, `n_modal`, `proportion_modal`, `entropy`,
#'   `average_posteriors` (rows: most likely class) and `classification`
#'   (rows: latent class, columns: most likely class).
.mplus_classification <- function(lines, n_profiles, n_group_classes) {
  level <- function(variable, k) {
    anchor <- grep(sprintf("Classification Quality for %s\\s*$", variable), lines)
    stopifnot("each level has one classification block" = length(anchor) == 1L)
    posterior <- .mplus_variable_classes(lines, "ESTIMATED POSTERIOR", variable, k)
    modal <- .mplus_variable_classes(lines, "MOST LIKELY", variable, k)
    list(count = posterior$count, proportion = posterior$proportion,
         n_modal = modal$count, proportion_modal = modal$proportion,
         entropy = mplus_entropy_after(lines, anchor),
         average_posteriors = mplus_avepp_matrix_after(lines, anchor, k),
         classification = .mplus_square_after(
           lines, "^Classification Probabilities for the Most Likely", anchor, k))
  }
  list(individuals = level("CW", n_profiles), groups = level("CB", n_group_classes))
}

# ---------------------------------------------------------------------------
# The two-level latent class model, from scratch
#
# Parameters are Mplus's, in Mplus's order: for each CW class s and item j the
# thresholds tau_sj1 < tau_sj2 with P(u <= c | s) = plogis(tau_sjc); the CW
# intercepts a_s and slopes b_sg, so that log P(cw = s | cb = g) / P(cw = S |
# cb = g) = a_s + b_sg (b_sG = 0); the CB means m_g = log P(cb = g) / P(cb = G).
# A cluster's likelihood is sum_g P(cb = g) prod_i sum_s P(cw = s | g) f_s(y_i).
# Nothing here calls the package.
# ---------------------------------------------------------------------------

#' Dimensions of the model
#' @param n_profiles,n_items,n_group_classes Its sizes (three categories per item).
#' @return A list with `S`, `J`, `B`, `n_thresholds` and `q`.
.lca2_dimensions <- function(n_profiles, n_items, n_group_classes) {
  n_thresholds <- 2L * n_items * n_profiles
  list(S = n_profiles, J = n_items, B = n_group_classes, n_thresholds = n_thresholds,
       q = n_thresholds + (n_profiles - 1L) * n_group_classes + n_group_classes - 1L)
}

#' Mplus's parameters to probabilities
#' @param theta Parameter vector in Mplus's order.
#' @param dims From [.lca2_dimensions()].
#' @return A list with `thresholds` (S x J x 2), `gamma` (B x S, P(cw | cb)),
#'   `delta` (P(cb)) and `prob` (S x J x 3, P(u = c | cw)).
.lca2_unpack <- function(theta, dims) {
  S <- dims$S
  B <- dims$B
  nt <- dims$n_thresholds
  thresholds <- aperm(array(theta[seq_len(nt)], c(2L, dims$J, S)), c(3L, 2L, 1L))
  intercepts <- theta[nt + seq_len(S - 1L)]
  slopes <- matrix(theta[nt + S - 1L + seq_len((S - 1L) * (B - 1L))], S - 1L, B - 1L,
                   byrow = TRUE)
  means <- theta[nt + (S - 1L) * B + seq_len(B - 1L)]
  softmax <- function(v) exp(v - max(v)) / sum(exp(v - max(v)))
  logits <- cbind(t(intercepts + cbind(slopes, 0)), 0)
  cumulative <- stats::plogis(thresholds)
  list(thresholds = thresholds,
       gamma = t(apply(logits, 1L, softmax)),
       delta = softmax(c(means, 0)),
       prob = array(c(cumulative[, , 1L], cumulative[, , 2L] - cumulative[, , 1L],
                      1 - cumulative[, , 2L]), c(S, dims$J, 3L)))
}

#' Probabilities to Mplus's parameters, the inverse of [.lca2_unpack()]
#' @param gamma,delta,prob As returned by [.lca2_unpack()].
#' @return The parameter vector in Mplus's order.
.lca2_pack <- function(gamma, delta, prob) {
  S <- ncol(gamma)
  B <- nrow(gamma)
  first <- stats::qlogis(prob[, , 1L, drop = FALSE])
  second <- stats::qlogis(prob[, , 1L, drop = FALSE] + prob[, , 2L, drop = FALSE])
  thresholds <- as.vector(aperm(array(c(first, second), c(dim(prob)[1:2], 2L)),
                                c(3L, 2L, 1L)))
  logits <- log(gamma[, -S, drop = FALSE] / gamma[, S])
  intercepts <- logits[B, ]
  slopes <- t(logits[-B, , drop = FALSE]) - intercepts
  c(thresholds, intercepts, as.vector(t(slopes)), log(delta[-B] / delta[B]))
}

#' The E step: likelihood and posteriors of both levels
#'
#' Within a cluster the joint posterior of (cb = g, cw = s) for row i is
#' `z_jg gamma_gs f_s(y_i) / L_ig`, with `L_ig = sum_s gamma_gs f_s(y_i)`; the
#' marginal CW posterior sums it over g. Densities are scaled by their row
#' maximum, so nothing underflows.
#'
#' @param parameters From [.lca2_unpack()].
#' @param y Rows-by-items matrix of category codes 1 to 3.
#' @param group Cluster index of each row, 1 to the number of clusters.
#' @return A list with `log_likelihood`, `group_log_likelihood`,
#'   `group_posteriors` (clusters x B), `posteriors` (rows x S), and the
#'   `scaled` densities and `ratio` (z_jg / L_ig, scaled) the scores reuse.
.lca2_estep <- function(parameters, y, group) {
  log_density <- Reduce(`+`, lapply(seq_len(ncol(y)), function(j) {
    t(log(parameters$prob[, j, ])[, y[, j], drop = FALSE])
  }))
  peak <- do.call(pmax, as.data.frame(log_density))
  scaled <- exp(log_density - peak)
  mix <- scaled %*% t(parameters$gamma)
  by_group <- rowsum(log(mix) + peak, group, reorder = TRUE) +
    rep(log(parameters$delta), each = max(group))
  top <- do.call(pmax, as.data.frame(by_group))
  group_ll <- top + log(rowSums(exp(by_group - top)))
  group_posteriors <- exp(by_group - group_ll)
  ratio <- group_posteriors[group, , drop = FALSE] / mix
  list(log_likelihood = sum(group_ll), group_log_likelihood = group_ll,
       group_posteriors = group_posteriors,
       posteriors = scaled * (ratio %*% parameters$gamma), scaled = scaled, ratio = ratio)
}

#' Each cluster's score: the derivative of its log-likelihood in Mplus's
#' parameters
#'
#' By Fisher's identity the score is the posterior expectation of the
#' complete-data score: for a threshold, the CW-posterior-weighted derivative
#' of the log category probability; for the logit of P(cw = s | cb = g), the
#' cluster's sum of `w_igs - z_jg gamma_gs`; for a CB mean, `z_jg - delta_g`.
#'
#' @inheritParams .lca2_estep
#' @param theta Parameter vector in Mplus's order.
#' @param dims From [.lca2_dimensions()].
#' @return A clusters-by-parameters matrix.
.lca2_scores <- function(theta, y, group, dims) {
  parameters <- .lca2_unpack(theta, dims)
  e <- .lca2_estep(parameters, y, group)
  cumulative <- stats::plogis(parameters$thresholds)
  density <- cumulative * (1 - cumulative)
  grid <- expand.grid(cut = 1:2, item = seq_len(dims$J), profile = seq_len(dims$S))
  thresholds <- vapply(seq_len(nrow(grid)), function(k) {
    s <- grid$profile[[k]]
    j <- grid$item[[k]]
    observed <- y[, j]
    lower <- cumulative[s, j, 1L]
    upper <- cumulative[s, j, 2L]
    derivative <- if (grid$cut[[k]] == 1L) {
      ifelse(observed == 1L, 1 - lower,
             ifelse(observed == 2L, -density[s, j, 1L] / (upper - lower), 0))
    } else {
      ifelse(observed == 2L, density[s, j, 2L] / (upper - lower),
             ifelse(observed == 3L, -upper, 0))
    }
    e$posteriors[, s] * derivative
  }, numeric(nrow(y)))
  sizes <- as.vector(table(group))
  logit_scores <- lapply(seq_len(dims$B), function(g) {
    joint <- rowsum(e$ratio[, g] * e$scaled, group, reorder = TRUE)
    sweep(joint, 2L, parameters$gamma[g, ], "*") -
      outer(sizes * e$group_posteriors[, g], parameters$gamma[g, ])
  })
  intercepts <- Reduce(`+`, logit_scores)[, seq_len(dims$S - 1L), drop = FALSE]
  slopes <- do.call(cbind, lapply(seq_len(dims$S - 1L), function(s) {
    vapply(seq_len(dims$B - 1L), function(g) logit_scores[[g]][, s], numeric(max(group)))
  }))
  means <- sweep(e$group_posteriors[, seq_len(dims$B - 1L), drop = FALSE], 2L,
                 parameters$delta[seq_len(dims$B - 1L)])
  cbind(rowsum(thresholds, group, reorder = TRUE), intercepts, slopes, means)
}

#' Observed information: central differences of the analytic score
#' @inheritParams .lca2_scores
#' @param step Difference step on the logit scale.
#' @return The symmetric Hessian of the log-likelihood.
.lca2_hessian <- function(theta, y, group, dims, step = 1e-4) {
  columns <- vapply(seq_along(theta), function(k) {
    shift <- replace(numeric(length(theta)), k, step)
    (colSums(.lca2_scores(theta + shift, y, group, dims)) -
       colSums(.lca2_scores(theta - shift, y, group, dims))) / (2 * step)
  }, numeric(length(theta)))
  (columns + t(columns)) / 2
}

#' Carry estimates to the maximum by Newton steps with a fixed Hessian
#'
#' The Hessian at the starting point is reused, so the iteration converges
#' linearly at the rate by which it differs from the Hessian on the way; each
#' step must not lose likelihood, so the climb stays on the starting mode.
#'
#' @inheritParams .lca2_scores
#' @param hessian Hessian at `theta`.
#' @param tolerance Largest absolute score accepted as a maximum.
#' @param max_steps Bound on the number of steps.
#' @return A list with `theta`, `log_likelihood`, `steps` and `score` (the
#'   largest absolute score at the result).
.lca2_refine <- function(theta, hessian, y, group, dims, tolerance = 1e-8, max_steps = 50L) {
  log_likelihood_at <- function(t) .lca2_estep(.lca2_unpack(t, dims), y, group)$log_likelihood
  climb <- function(t, log_likelihood, steps) {
    score <- colSums(.lca2_scores(t, y, group, dims))
    if (max(abs(score)) < tolerance) {
      return(list(theta = t, log_likelihood = log_likelihood, steps = steps,
                  score = max(abs(score))))
    }
    if (steps >= max_steps) {
      stop(sprintf("Newton refinement did not reach a score below %g in %d steps",
                   tolerance, max_steps), call. = FALSE)
    }
    proposal <- t - solve(hessian, score)
    proposed <- log_likelihood_at(proposal)
    stopifnot("a refining step must not lose likelihood" = proposed >= log_likelihood - 1e-9)
    climb(proposal, proposed, steps + 1L)
  }
  climb(theta, log_likelihood_at(theta), 0L)
}

#' Every probability of the model as one vector
#' @inheritParams .lca2_scores
#' @return `prob` (S x J x 3), then `gamma` (B x S), then `delta`, as vectors.
.lca2_natural <- function(theta, dims) {
  parameters <- .lca2_unpack(theta, dims)
  c(as.vector(parameters$prob), as.vector(parameters$gamma), parameters$delta)
}

#' The covariances Mplus and the package report, from scratch
#'
#' `observed` is the inverse observed information; `robust` the sandwich
#' `A B A` with `B` the sum of the cluster score cross-products, which is
#' Mplus's MLR (it reproduces Mplus's printed standard errors without a
#' G / (G - 1) factor); `scaling` is `tr(A B) / q`, Mplus's MLR scaling
#' correction factor. `natural_jacobian` carries any of them to the
#' probabilities by the delta method.
#'
#' @inheritParams .lca2_scores
#' @param hessian Hessian at `theta`.
#' @return A list with `observed`, `robust`, `scaling` and `natural_jacobian`.
.lca2_covariances <- function(theta, y, group, dims, hessian) {
  scores <- .lca2_scores(theta, y, group, dims)
  bread <- solve(-hessian)
  meat <- crossprod(scores)
  step <- 1e-6
  jacobian <- vapply(seq_along(theta), function(k) {
    shift <- replace(numeric(length(theta)), k, step)
    (.lca2_natural(theta + shift, dims) - .lca2_natural(theta - shift, dims)) / (2 * step)
  }, numeric(length(.lca2_natural(theta, dims))))
  list(observed = bread, robust = bread %*% meat %*% bread,
       scaling = sum(diag(bread %*% meat)) / length(theta), natural_jacobian = jacobian)
}

#' The from-scratch solution at `theta` in the harness's reference shape
#'
#' @inheritParams .lca2_scores
#' @param ids Cluster identifier of each cluster index.
#' @param data The indicator data, in the fit's row order.
#' @param levels Named list of each item's category labels, as the fit has them.
#' @param covariances From [.lca2_covariances()], or `NULL` for no standard
#'   errors.
#' @return A `harness_reference`, in Mplus's class order.
.lca2_reference <- function(theta, y, group, dims, ids, data, levels, covariances = NULL) {
  parameters <- .lca2_unpack(theta, dims)
  e <- .lca2_estep(parameters, y, group)
  items <- names(levels)
  responses <- stats::setNames(lapply(seq_along(items), function(j) {
    matrix(parameters$prob[, j, ], dims$S, 3L, dimnames = list(NULL, levels[[j]]))
  }), items)
  standard_errors <- if (!is.null(covariances)) {
    natural <- covariances$natural_jacobian
    errors <- sqrt(diag(natural %*% covariances$observed %*% t(natural)))
    by_cell <- array(errors[seq_len(dims$S * dims$J * 3L)], c(dims$S, dims$J, 3L))
    list(probability = stats::setNames(lapply(seq_along(items), function(j) {
      matrix(by_cell[, j, ], dims$S, 3L, dimnames = list(NULL, levels[[j]]))
    }), items))
  }
  group_posteriors <- e$group_posteriors
  rownames(group_posteriors) <- as.character(ids)
  harness_reference(
    log_likelihood = e$log_likelihood, n_parameters = dims$q, posteriors = e$posteriors,
    groups = ids[group], responses = responses, standard_errors = standard_errors,
    profile_probabilities = parameters$gamma, group_probabilities = parameters$delta,
    group_posteriors = group_posteriors,
    group_log_likelihood = stats::setNames(e$group_log_likelihood, as.character(ids)),
    data = data,
    row_conditional_entropy = two_level_row_conditional_entropy(
      e$scaled, parameters$gamma, e$group_posteriors[group, , drop = FALSE]))
}

# ---------------------------------------------------------------------------
# How far rounding the published estimates can move a result
# ---------------------------------------------------------------------------

#' The most each classification result of one level can move under rounding
#'
#' The estimates are printed to three decimals, so Mplus's own estimates lie
#' within `half_width` of each printed one. To first order a smooth result `R`
#' moves by at most `half_width * sum_k |dR / d theta_k|`, the bound used for the
#' counts, proportions, entropy and, with the most likely classes held, the
#' average posteriors and classification probabilities. A most likely class
#' can change only for a unit whose leading posterior some other class can
#' overtake within those bounds; such units are counted as moved, and each one
#' is charged its full effect on the results that depend on the assignment.
#'
#' @param p Posterior matrix of the level, units by classes, in Mplus's labels.
#' @param jacobian Units x classes x parameters array, `d p / d theta`.
#' @param half_width Half the printed resolution of the estimates.
#' @param weight Observations per unit (Mplus counts CB over observations).
#' @return A list shaped like one level of [.mplus_classification()], holding
#'   the allowances.
.lca2_level_allowance <- function(p, jacobian, half_width, weight) {
  n <- nrow(p)
  k <- ncol(p)
  q <- dim(jacobian)[3L]
  flat <- matrix(jacobian, n * k, q)
  reach <- half_width * matrix(rowSums(abs(flat)), n, k)
  modal <- max.col(p, ties.method = "first")
  n_modal <- tabulate(modal, k)
  stopifnot("every class is some unit's most likely class" = all(n_modal > 0L))
  leading <- cbind(seq_len(n), modal)
  can_take <- (p + reach >= p[leading] - reach[leading]) & col(p) != modal
  can_leave <- rowSums(can_take) > 0L
  moved <- tabulate(modal[can_leave], k) + colSums(can_take)
  stopifnot("rounding cannot empty a class" = all(moved < n_modal))
  count_gradient <- apply(jacobian, c(2L, 3L), sum)
  count <- half_width * rowSums(abs(count_gradient))
  entropy_gradient <- as.vector(crossprod(as.vector(ifelse(p > 0, log(p) + 1, 0)), flat))
  assigned <- array(rowsum(matrix(jacobian, n, k * q), modal, reorder = TRUE), c(k, k, q))
  assigned_mass <- rowsum(p, modal, reorder = TRUE)
  totals <- colSums(p)
  average_posteriors <- half_width * apply(abs(assigned), c(1L, 2L), sum) / n_modal +
    moved / (n_modal - moved)
  involved <- (outer(modal, seq_len(k), `==`) & can_leave) | can_take
  classification <- t(vapply(seq_len(k), function(true_class) {
    vapply(seq_len(k), function(assigned_class) {
      gradient <- (assigned[assigned_class, true_class, ] * totals[[true_class]] -
                     assigned_mass[assigned_class, true_class] * count_gradient[true_class, ]) /
        totals[[true_class]]^2
      half_width * sum(abs(gradient)) +
        sum(p[involved[, assigned_class], true_class]) / totals[[true_class]]
    }, numeric(1))
  }, numeric(k)))
  list(count = weight * count, proportion = count / n,
       n_modal = weight * moved, proportion_modal = moved / n,
       entropy = half_width * sum(abs(entropy_gradient)) / (n * log(k)),
       average_posteriors = average_posteriors, classification = classification)
}

#' Rounding allowances of every classification result, both levels
#' @inheritParams .lca2_scores
#' @param half_width Half the printed resolution of the estimates.
#' @param weight Observations per cluster.
#' @param step Difference step for the posterior Jacobian.
#' @return A list with `individuals` and `groups`, from
#'   [.lca2_level_allowance()].
.lca2_rounding_allowance <- function(theta, y, group, dims, half_width, weight,
                                     step = 1e-5) {
  posteriors_at <- function(t) {
    e <- .lca2_estep(.lca2_unpack(t, dims), y, group)
    c(as.vector(e$posteriors), as.vector(e$group_posteriors))
  }
  jacobian <- vapply(seq_along(theta), function(k) {
    shift <- replace(numeric(length(theta)), k, step)
    (posteriors_at(theta + shift) - posteriors_at(theta - shift)) / (2 * step)
  }, numeric(nrow(y) * dims$S + max(group) * dims$B))
  e <- .lca2_estep(.lca2_unpack(theta, dims), y, group)
  within <- seq_len(nrow(y) * dims$S)
  list(individuals = .lca2_level_allowance(
         e$posteriors, array(jacobian[within, ], c(nrow(y), dims$S, length(theta))),
         half_width, 1),
       groups = .lca2_level_allowance(
         e$group_posteriors,
         array(jacobian[-within, ], c(max(group), dims$B, length(theta))),
         half_width, weight))
}

#' How far rounding the estimates moves the MLR standard errors
#'
#' A standard error is a function of the Hessian, so a first-order bound would
#' need a third derivative per parameter. The shift is measured instead:
#' printed rounding errors are modelled as independent and uniform within
#' `half_width`, the MLR standard errors are recomputed at `n_draws` such
#' points, and the allowance is four times the root-mean-square shift (about
#' four standard deviations of the shift the actual rounding causes). The seed
#' is fixed and the caller's random-number state restored.
#'
#' @inheritParams .lca2_scores
#' @param half_width Half the printed resolution of the estimates.
#' @param at The standard errors and scaling factor at `theta`.
#' @param n_draws Number of points.
#' @return One allowance per parameter, then one for the scaling factor.
.lca2_standard_error_allowance <- function(theta, y, group, dims, half_width, at,
                                           n_draws = 8L) {
  had_seed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  if (had_seed) saved <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
  on.exit(if (had_seed) assign(".Random.seed", saved, envir = globalenv()) else
    rm(".Random.seed", envir = globalenv()), add = TRUE)
  set.seed(20261026L)
  shifts <- vapply(seq_len(n_draws), function(draw) {
    moved <- theta + stats::runif(length(theta), -half_width, half_width)
    covariances <- .lca2_covariances(moved, y, group, dims, .lca2_hessian(moved, y, group, dims))
    c(sqrt(diag(covariances$robust)), covariances$scaling) - at
  }, numeric(length(theta) + 1L))
  4 * sqrt(rowMeans(shifts^2))
}

# ---------------------------------------------------------------------------
# The package's results in Mplus's labels
# ---------------------------------------------------------------------------

#' One level's classification results of a fit, in Mplus's class labels
#' @param tables From [fit_tables()].
#' @param which `"individuals"` or `"groups"`.
#' @param map The fit's class for each Mplus class.
#' @param weight Observations per unit, to count as Mplus does.
#' @return A list shaped like one level of [.mplus_classification()].
.package_level <- function(tables, which, map, weight) {
  k <- length(map)
  in_order <- function(table) {
    rows <- subset(table, level == which)
    rows[order(rows$class), , drop = FALSE]
  }
  counts <- in_order(tables$counts)
  classes <- in_order(tables$classification)
  entropy <- subset(tables$entropy, level == which)
  square <- function(table, row, column, value) {
    rows <- subset(table, level == which)
    out <- matrix(NA_real_, k, k)
    out[cbind(rows[[row]], rows[[column]])] <- rows[[value]]
    out[map, map, drop = FALSE]
  }
  list(count = weight * counts$effective_count[map],
       proportion = counts$effective_proportion[map],
       n_modal = weight * classes$n_modal[map],
       proportion_modal = classes$proportion_modal[map],
       entropy = entropy$relative_entropy,
       average_posteriors = square(tables$average_posteriors, "assigned_class", "class",
                                   "average_posterior"),
       classification = square(tables$classification_errors, "true_class", "assigned_class",
                               "probability"))
}

#' Compare Mplus's printed classification results with a fit's, one level
#' @param label Row label prefix.
#' @param printed,obtained,allowance Lists shaped like one level of
#'   [.mplus_classification()].
#' @param noun `"CW"` or `"CB"`.
#' @return A comparison `data.frame`.
.compare_classification <- function(label, printed, obtained, allowance, noun) {
  k <- length(printed$count)
  cells <- sprintf("%d,%d", rep(seq_len(k), times = k), rep(seq_len(k), each = k))
  # Mplus's printed decimals: counts and proportions five, the rest three.
  printed_rows <- function(what, digits, output, names) {
    compare_values(sprintf("%s: Mplus's %s %s %s", label, noun, what, names),
                   as.numeric(printed[[what]]), as.numeric(obtained[[what]]),
                   tolerance = 0.5 * 10^-digits + as.numeric(allowance[[what]]),
                   output = output)
  }
  rbind(
    printed_rows("count", 5, "counts$effective_count", seq_len(k)),
    printed_rows("proportion", 5, "counts$effective_proportion", seq_len(k)),
    printed_rows("n_modal", 0, "classification$n_modal", seq_len(k)),
    printed_rows("proportion_modal", 5, "classification$proportion_modal", seq_len(k)),
    printed_rows("entropy", 3, "entropy$relative_entropy", ""),
    printed_rows("average_posteriors", 3, "average_posteriors$average_posterior",
                 paste("(most likely, class)", cells)),
    printed_rows("classification", 3, "classification_errors$probability",
                 paste("(class, most likely)", cells)))
}

suite_mplus_twolevel_lca <- function() {
  artifacts <- file.path("equivalence", "mplus", "twolevel-lca")
  require_suite_files(
    file.path(artifacts, c("ex10.7.out.gz", "ex10.7.dat.gz")),
    reason = "the retained Mplus 8.8 run of User's Guide example 10.7 and its data")
  # The classification parsers shared with the Mplus classification tests, into
  # the environment this suite was sourced into, where its parsers look them up.
  sys.source(file.path("equivalence", "helper-mplus-classification.R"),
             envir = environment(sys.function()))
  items <- paste0("u", seq_len(10))
  n_profiles <- 4L
  n_group_classes <- 5L
  n_categories <- 3L
  dims <- .lca2_dimensions(n_profiles, length(items), n_group_classes)
  # Mplus prints every estimate to three decimals.
  half_width <- 5e-4

  lines <- .mplus_output_lines(file.path(artifacts, "ex10.7.out.gz"))
  stopifnot(
    "Mplus did not terminate normally" =
      any(grepl("THE MODEL ESTIMATION TERMINATED NORMALLY", lines)),
    "Mplus did not replicate the best loglikelihood" =
      any(grepl("THE BEST LOGLIKELIHOOD VALUE HAS BEEN REPLICATED", lines)),
    "Mplus reported an error" = !any(grepl("*** ERROR", lines, fixed = TRUE))
  )

  observed <- utils::read.table(gzfile(file.path(artifacts, "ex10.7.dat.gz")),
                                col.names = c(items, "dumb", "dumw", "clus"))
  ids <- unique(observed$clus)
  group <- match(observed$clus, ids)
  cluster_size <- unique(as.vector(table(group)))
  # Mplus counts CB over observations; with clusters of one size that is the
  # per-cluster count times the size, and its averages are per-cluster averages.
  stopifnot("every cluster has the same size" = length(cluster_size) == 1L)
  y <- as.matrix(observed[items]) + 1L
  stopifnot("items are coded 0, 1, 2" = all(y %in% seq_len(n_categories)))

  # The package stops when the log-likelihood changes by less than
  # tol * (1 + |LL|) in an iteration. At tol = 1e-12 that is 1e-7 here, and on
  # this flat likelihood it stops 2.4e-8 below the maximum with a largest score
  # of 1.2e-3 and estimates 5.8e-5 from the maximum (measured 2026-10-26), too
  # far to compare posteriors to 1e-6. At 1e-15, as in the mclust suite, its
  # log-likelihood matches the from-scratch maximum to 1e-10.
  fit <- multilpa(observed, vars = items, id = "clus",
                  n_profiles = n_profiles, n_group_classes = n_group_classes,
                  categorical = items, n_starts = 20L, seed = 20260919L,
                  tol = 1e-15, max_iter = 20000L)
  levels <- fit$categorical_levels[items]
  stopifnot("the package labels the categories 0, 1, 2" =
              all(vapply(levels, identical, logical(1), c("0", "1", "2"))))

  # ---- Mplus's published estimates --------------------------------------
  thresholds <- .mplus_thresholds(lines)
  # The measurement model is shared across group classes, so any one group
  # class carries it. Verify that before relying on it.
  spread <- aggregate(cbind(estimate, standard_error) ~ profile + item + cut,
                      data = thresholds, FUN = function(v) diff(range(v)))
  stopifnot("Mplus did not hold the measurement model fixed across group classes" =
              all(spread$estimate < 1e-9 & spread$standard_error < 1e-9))
  thresholds <- subset(thresholds, group_class == 1L)
  thresholds <- thresholds[order(thresholds$profile, thresholds$item, thresholds$cut), ]
  logits <- .mplus_class_logits(lines, n_profiles, n_group_classes)
  theta_published <- c(thresholds$estimate, logits$estimate)
  printed_errors <- c(thresholds$standard_error, logits$standard_error)
  stopifnot("Mplus printed one estimate per free parameter" =
              length(theta_published) == dims$q &&
              identical(as.integer(.mplus_scalar(lines, "Number of Free Parameters")), dims$q))
  published <- .lca2_unpack(theta_published, dims)

  # The package evaluated at Mplus's published estimates: started there and not
  # iterated. The thresholds and logits are carried to probabilities exactly,
  # so the package holds the printed estimates themselves.
  at_mplus <- multilpa(observed, vars = items, id = "clus",
                       n_profiles = n_profiles, n_group_classes = n_group_classes,
                       categorical = items, n_starts = 1L, max_iter = 0L,
                       start = starting_values(list(
                         profile_probabilities = published$gamma,
                         group_probabilities = published$delta,
                         response_probabilities = lapply(seq_along(items), function(j) {
                           matrix(published$prob[, j, ], n_profiles, n_categories)
                         }))))

  # ---- The from-scratch reference at the published estimates and at the max
  hessian_published <- .lca2_hessian(theta_published, y, group, dims)
  covariances_published <- .lca2_covariances(theta_published, y, group, dims,
                                             hessian_published)
  refined <- .lca2_refine(theta_published, hessian_published, y, group, dims)
  hessian_max <- .lca2_hessian(refined$theta, y, group, dims)
  covariances_max <- .lca2_covariances(refined$theta, y, group, dims, hessian_max)
  data_items <- observed[items]
  reference_published <- .lca2_reference(theta_published, y, group, dims, ids,
                                         data_items, levels)
  reference_max <- .lca2_reference(refined$theta, y, group, dims, ids, data_items, levels,
                                   covariances_max)

  # ---- Case 1: at Mplus's published estimates -----------------------------
  label_published <- "ex10.7 at Mplus's published estimates"
  exempt_published <- c(
    `responses$probability_standard_error` = paste(
      "the package withholds standard errors at a point it did not converge to",
      "(max_iter = 0; 'Inference requires a converged fit'); they are compared",
      "at the maximum in the converged case"))
  # Both sides compute every table from the same parameters, so they agree to
  # rounding error (observed: 9e-14 at worst); 1e-10 is the claim made.
  exact <- c(default = 1e-10, `classification$odds_correct_classification` = 1e-10,
             `responses$threshold` = 1e-10)
  matched_published <- compare_fit(label_published, at_mplus, reference_published,
                                   tolerance = exact, exempt = exempt_published)

  tables <- fit_tables(at_mplus)
  # The package's class for each Mplus class.
  profile_map <- align_classes(
    flatten_blocks(lapply(seq_along(items), function(j) published$prob[, j, ])),
    flatten_blocks(at_mplus$response_probabilities))
  group_map <- align_classes(published$gamma,
                             unname(at_mplus$profile_probabilities)[, profile_map, drop = FALSE])

  # Thresholds: the package's `threshold` is logit P(u <= c), the category
  # labels 0, 1, 2 are Mplus's categories 1, 2, 3, so Mplus's `U$1` and `U$2`
  # are the package's thresholds of categories "0" and "1", with the same sign.
  # Evaluated at the printed estimates, the package must give back the printed
  # digits; a reversed sign or a P(u > c) convention would not.
  responses <- tables$responses
  response_key <- paste(responses$profile, responses$indicator, responses$category)
  threshold_key <- paste(profile_map[thresholds$profile], items[thresholds$item],
                         c("0", "1")[thresholds$cut])
  probability_scale <- .mplus_probability_scale(lines, toupper(items), n_categories)
  first_group <- which(probability_scale$pattern[, 1L] == 1L)
  printed_probabilities <- probability_scale$probabilities[
    first_group[order(probability_scale$pattern[first_group, 2L])], , , drop = FALSE]
  probability_grid <- expand.grid(profile = seq_len(n_profiles), item = seq_len(length(items)),
                                  category = seq_len(n_categories))
  probability_key <- paste(profile_map[probability_grid$profile], items[probability_grid$item],
                           c("0", "1", "2")[probability_grid$category])
  # Mplus printed P(u = c) from its unrounded thresholds: P(u = 1) = F(tau_1)
  # moves by at most f(tau_1) times the rounding, P(u = 3) by f(tau_2), and
  # P(u = 2) by their sum (f is the logistic density).
  logistic_density <- stats::dlogis(published$thresholds)
  category_reach <- half_width * array(
    c(logistic_density[, , 1L], logistic_density[, , 1L] + logistic_density[, , 2L],
      logistic_density[, , 2L]), c(n_profiles, length(items), n_categories))

  # The logits, from the package's prevalence table: CW's last class and CB's
  # last class are Mplus's reference categories.
  prevalence <- tables$profile_probabilities
  gamma <- matrix(NA_real_, n_group_classes, n_profiles)
  gamma[cbind(prevalence$group_class, prevalence$profile)] <- prevalence$probability
  gamma <- gamma[group_map, profile_map, drop = FALSE]
  delta <- vapply(group_map, function(g) {
    unique(prevalence$group_class_probability[prevalence$group_class == g])
  }, numeric(1))
  package_logits <- .lca2_pack(gamma, delta, published$prob)[-seq_len(dims$n_thresholds)]
  intercepts_and_slopes <- seq_len((n_profiles - 1L) * n_group_classes)

  # The log-likelihood at rounded estimates: to second order it moves by at
  # most h sum|g| + h^2 sum|H| / 2 with g the score and H the Hessian there.
  score_published <- colSums(.lca2_scores(theta_published, y, group, dims))
  likelihood_reach <- half_width * sum(abs(score_published)) +
    half_width^2 * sum(abs(hessian_published)) / 2
  criteria <- subset(tables$information_criteria, select = c(log_likelihood, aic, bic_individual,
                                                             sabic_individual, n_parameters))
  allowance <- .lca2_rounding_allowance(theta_published, y, group, dims, half_width,
                                        weight = cluster_size)
  printed_classification <- .mplus_classification(lines, n_profiles, n_group_classes)

  errors_published <- sqrt(diag(covariances_published$robust))
  error_reach <- .lca2_standard_error_allowance(
    theta_published, y, group, dims, half_width,
    at = c(errors_published, covariances_published$scaling))
  parameter_labels <- c(sprintf("U%d$%d in CW class %d", thresholds$item, thresholds$cut,
                                thresholds$profile), logits$label)

  printed_published <- rbind(
    compare_printed(sprintf("%s: Mplus's threshold %s", label_published,
                            parameter_labels[seq_len(dims$n_thresholds)]),
                    thresholds$estimate,
                    responses$threshold[match(threshold_key, response_key)],
                    digits = 3, output = "responses$threshold"),
    compare_values(sprintf("%s: Mplus's probability-scale P(%s = %d | CW class %d)",
                           label_published, toupper(items[probability_grid$item]),
                           probability_grid$category, probability_grid$profile),
                   as.vector(printed_probabilities),
                   responses$probability[match(probability_key, response_key)],
                   tolerance = half_width + as.vector(category_reach),
                   output = "responses$probability"),
    compare_printed(sprintf("%s: Mplus's logit %s", label_published,
                            logits$label[intercepts_and_slopes]),
                    logits$estimate[intercepts_and_slopes],
                    package_logits[intercepts_and_slopes], digits = 3,
                    output = "profile_probabilities$probability"),
    compare_printed(sprintf("%s: Mplus's logit %s", label_published,
                            logits$label[-intercepts_and_slopes]),
                    logits$estimate[-intercepts_and_slopes],
                    package_logits[-intercepts_and_slopes], digits = 3,
                    output = "profile_probabilities$group_class_probability"),
    compare_values(sprintf("%s: Mplus's %s", label_published,
                           c("H0 log-likelihood", "AIC", "BIC", "sample-size adjusted BIC")),
                   c(.mplus_scalar(lines, "H0 Value"), .mplus_scalar(lines, "Akaike (AIC)"),
                     .mplus_scalar(lines, "Bayesian (BIC)"),
                     .mplus_scalar(lines, "Sample-Size Adjusted BIC")),
                   c(criteria$log_likelihood, criteria$aic, criteria$bic_individual,
                     criteria$sabic_individual),
                   tolerance = half_width + c(1, 2, 2, 2) * likelihood_reach,
                   output = c("information_criteria$log_likelihood", "information_criteria$aic",
                              "information_criteria$bic_individual",
                              "information_criteria$sabic_individual")),
    compare_printed(sprintf("%s: Mplus's %s", label_published,
                            c("number of free parameters", "number of observations",
                              "number of clusters")),
                    c(.mplus_scalar(lines, "Number of Free Parameters"),
                      .mplus_scalar(lines, "Number of observations"),
                      .mplus_scalar(lines, "Number of clusters")),
                    c(criteria$n_parameters, tables$model$n_observations, tables$model$n_groups),
                    digits = 0, output = c("information_criteria$n_parameters",
                                           "model$n_observations", "model$n_groups")),
    .compare_classification(label_published, printed_classification$individuals,
                            .package_level(tables, "individuals", profile_map, 1),
                            allowance$individuals, "CW"),
    .compare_classification(label_published, printed_classification$groups,
                            .package_level(tables, "groups", group_map, cluster_size),
                            allowance$groups, "CB"),
    # Kept from before the ledger: the likelihood at Mplus's published
    # estimates, now the exact printed thresholds and logits rather than the
    # three-decimal probability scale.
    compare_values("log-likelihood at Mplus's published estimates",
                   .mplus_scalar(lines, "H0 Value"), as.numeric(logLik(at_mplus)),
                   tolerance = 1e-6, scale = "relative",
                   output = "information_criteria$log_likelihood")
  ) |> transform(precision = "printed")

  # Mplus's standard errors have no column in the package's tables (they are
  # MLR, on the threshold and logit scale); they check the from-scratch
  # information that supplies the reference standard errors below.
  reference_checks <- rbind(
    compare_values(sprintf("%s: Mplus's MLR S.E. of %s, from-scratch sandwich", label_published,
                           parameter_labels),
                   printed_errors, errors_published,
                   tolerance = half_width + error_reach[seq_len(dims$q)]),
    compare_values(sprintf("%s: Mplus's MLR scaling correction factor, from scratch",
                           label_published),
                   .mplus_scalar(lines, "H0 Scaling Correction Factor"),
                   covariances_published$scaling,
                   tolerance = 5e-5 + error_reach[[dims$q + 1L]])
  ) |> transform(precision = "printed")

  comparisons_published <- rbind(
    transform(matched_published$comparisons, precision = "machine"),
    printed_published, reference_checks)
  case_published <- harness_case(
    comparisons_published,
    output_ledger(label_published, fit_outputs(at_mplus), comparisons_published,
                  exempt_published))

  # ---- Case 2: converged --------------------------------------------------
  label_converged <- "ex10.7 converged"
  # The fit's class for each Mplus class.
  fit_profiles <- align_classes(
    flatten_blocks(lapply(seq_along(items), function(j) reference_max$responses[[j]])),
    flatten_blocks(fit$response_probabilities))
  fit_groups <- align_classes(
    reference_max$profile_probabilities,
    unname(fit$profile_probabilities)[, fit_profiles, drop = FALSE])

  # Where the package stopped. Even at tol = 1e-15 its stopping rule ends on a
  # change of 1e-10 in a log-likelihood of 1e5, which leaves its estimates up to
  # about 1e-5 from the maximum along the flattest directions (a cluster's
  # posterior, which sums fifty rows' evidence, moves by 1.8e-6 for it). That
  # is measured here, not assumed: every table is computed from scratch at the
  # package's own estimates and at the maximum, and the gap between the two,
  # column by column, is the part of any disagreement the stopping point
  # explains. The tolerance at the maximum is the harness's own for the column
  # plus that measured gap. The package's cluster posteriors are also checked
  # against the from-scratch evaluation at its own estimates, to machine
  # precision, which shows the gap is the stopping point and not a computation.
  theta_package <- .lca2_pack(
    unname(fit$profile_probabilities), unname(fit$group_probabilities),
    aperm(vapply(fit$response_probabilities[items], function(block) unname(as.matrix(block)),
                 matrix(0, n_profiles, n_categories)), c(1L, 3L, 2L)))
  reference_package <- .lca2_reference(theta_package, y, group, dims, ids, data_items, levels)
  stopping <- compare_tables("stopping point",
                             reference_tables(align_reference(reference_package, fit), fit),
                             reference_tables(align_reference(reference_max, fit), fit),
                             tolerance = c(default = 1))
  stopping <- subset(stopping, !is.na(output))
  harness_tolerance <- vapply(strsplit(stopping$output, "$", fixed = TRUE), function(key) {
    .output_tolerance(.output_default_tolerances, key[[1L]], key[[2L]])
  }, numeric(1))
  tolerance_converged <- c(
    stats::setNames(harness_tolerance + stopping$difference, stopping$output),
    # Both sides invert an observed information at the same point; they differ
    # by the numerical second derivatives, observed at 5e-10.
    `responses$probability_standard_error` = 1e-6)
  at_package <- .lca2_estep(.lca2_unpack(theta_package, dims), y, group)
  package_groups <- fit_tables(fit)$group_posteriors
  package_groups <- package_groups[order(package_groups$group_class,
                                         match(package_groups$group, as.character(ids))), ]
  stopping_rows <- compare_block(
    sprintf("%s: cluster posteriors, against the from-scratch likelihood at the package's own estimates",
            label_converged),
    as.vector(at_package$group_posteriors), package_groups$posterior,
    tolerance = 1e-10, output = "group_posteriors$posterior")
  matched_converged <- compare_fit(label_converged, fit, reference_max,
                                   tolerance = tolerance_converged)

  # The package's MLR ("robust") standard errors against the from-scratch
  # sandwich at the same maximum. The package multiplies the cluster
  # cross-product by G / (G - 1) (CR1); Mplus does not, which is why the
  # from-scratch sandwich without it reproduces Mplus's printed S.E. above.
  robust <- parameter_inference(fit, vcov_type = "robust")
  n_clusters <- length(ids)
  natural <- covariances_max$natural_jacobian
  sandwich <- sqrt(diag(natural %*% covariances_max$robust %*% t(natural)) *
                     n_clusters / (n_clusters - 1))
  natural_grid <- rbind(
    data.frame(parameter = "response", outcome = sprintf("profile_%d",
                                                         fit_profiles[probability_grid$profile]),
               term = paste(items[probability_grid$item],
                            c("0", "1", "2")[probability_grid$category], sep = ":")),
    data.frame(parameter = "probability",
               outcome = sprintf("profile_%d", rep(fit_profiles, each = n_group_classes)),
               term = sprintf("group_class_%d", rep(fit_groups, times = n_profiles))),
    data.frame(parameter = "probability", outcome = sprintf("group_class_%d", fit_groups),
               term = NA_character_))
  robust_key <- paste(robust$parameter, robust$outcome, robust$term)
  robust_errors <- robust$standard_error[match(paste(natural_grid$parameter, natural_grid$outcome,
                                                     natural_grid$term), robust_key)]
  blocks <- list(response = seq_len(dims$S * dims$J * n_categories),
                 profile = dims$S * dims$J * n_categories + seq_len(n_group_classes * n_profiles),
                 group = length(sandwich) - rev(seq_len(n_group_classes)) + 1L)
  # Both sides are at the same maximum and both use exact scores, so they
  # differ only by the numerical second derivatives on each side.
  robust_rows <- do.call(rbind, lapply(names(blocks), function(block) {
    compare_block(sprintf("%s: robust (MLR, CR1) S.E. of the %s probabilities, against the from-scratch sandwich",
                          label_converged, block),
                  sandwich[blocks[[block]]], robust_errors[blocks[[block]]], tolerance = 1e-6)
  }))

  prevalence_mplus <- .mplus_prevalence(lines, n_profiles, n_group_classes)
  mplus_responses <- lapply(seq_along(items), function(item) {
    .renormalize_printed(matrix(printed_probabilities[, item, ], n_profiles, n_categories), 1e-3)
  })
  profile_order <- align_classes(flatten_blocks(mplus_responses),
                                 flatten_blocks(fit$response_probabilities))
  fit_responses <- lapply(fit$response_probabilities,
                          function(block) unname(as.matrix(block))[profile_order, , drop = FALSE])
  group_order <- align_classes(prevalence_mplus$profile_probabilities,
                               unname(fit$profile_probabilities)[, profile_order, drop = FALSE])
  profile_probabilities <- unname(fit$profile_probabilities)[group_order, profile_order,
                                                             drop = FALSE]
  group_probabilities <- unname(fit$group_probabilities)[group_order]
  cell_names <- sprintf("item-response P(%s = %d | profile %d)",
                        rep(items, each = n_profiles * n_categories),
                        rep(rep(seq_len(n_categories), each = n_profiles), length(items)),
                        rep(seq_len(n_profiles), n_categories * length(items)))

  # Kept from before the ledger: the converged fit against Mplus's printed
  # estimates. These four are compared on a relative scale. Mplus prints them to
  # three decimals, but multilpa does not reproduce those digits: it converges a
  # little past Mplus's reported maximum, to -98736.121 against -98736.127 (see
  # the header: Mplus stopped short). On a log-likelihood of magnitude 98,736 an
  # absolute gap of 0.006 is the wrong thing to report, and "reproduces the
  # printed digits" is the wrong question; what matters is that the two agree
  # to within 1e-7 of their own size.
  kept <- rbind(
    compare_values("maximised log-likelihood",
                   .mplus_scalar(lines, "H0 Value"), as.numeric(logLik(fit)),
                   tolerance = 1e-6, scale = "relative",
                   output = "information_criteria$log_likelihood"),
    compare_values("free parameters",
                   .mplus_scalar(lines, "Number of Free Parameters"),
                   fit$n_parameters, tolerance = 0.5, output = "information_criteria$n_parameters"),
    compare_values("Akaike (AIC)", .mplus_scalar(lines, "Akaike (AIC)"),
                   fit$aic, tolerance = 1e-6, scale = "relative",
                   output = "information_criteria$aic"),
    compare_values("Bayesian (BIC)", .mplus_scalar(lines, "Bayesian (BIC)"),
                   fit$bic_individual, tolerance = 1e-6, scale = "relative",
                   output = "information_criteria$bic_individual"),
    compare_values("sample-size adjusted BIC",
                   .mplus_scalar(lines, "Sample-Size Adjusted BIC"),
                   subset(fit_tables(fit)$information_criteria, select = sabic_individual)[[1L]],
                   tolerance = 1e-6, scale = "relative",
                   output = "information_criteria$sabic_individual"),
    # Prevalence is recovered from logits printed to three decimals whose
    # reported standard errors are 0.22 to 0.41, so this tolerance sits about
    # two orders of magnitude below the sampling uncertainty Mplus reports.
    compare_values(sprintf("group-class probability %d", seq_len(n_group_classes)),
                   prevalence_mplus$group_probabilities, group_probabilities,
                   tolerance = 5e-3, output = "profile_probabilities$group_class_probability"),
    compare_values(sprintf("profile prevalence P(profile %d | group class %d)",
                           rep(seq_len(n_profiles), each = n_group_classes),
                           rep(seq_len(n_group_classes), n_profiles)),
                   as.numeric(prevalence_mplus$profile_probabilities),
                   as.numeric(profile_probabilities), tolerance = 5e-3,
                   output = "profile_probabilities$probability"),
    compare_values(cell_names,
                   as.numeric(vapply(mplus_responses, identity,
                                     matrix(0, n_profiles, n_categories))),
                   as.numeric(vapply(fit_responses, identity,
                                     matrix(0, n_profiles, n_categories))),
                   tolerance = 1.5e-3, output = "responses$probability")
  ) |> transform(precision = "printed")

  converged_checks <- rbind(
    stopping_rows,
    robust_rows,
    compare_values(sprintf("%s: MLR scaling correction factor, against the from-scratch tr(A B) / q",
                           label_converged),
                   covariances_max$scaling, attr(robust, "scaling_correction"),
                   tolerance = 1e-6),
    compare_values(sprintf("%s: log-likelihood of Mplus's estimates refined to the maximum", label_converged),
                   refined$log_likelihood, as.numeric(logLik(fit)),
                   tolerance = 1e-10, scale = "relative")
  ) |> transform(precision = "machine")

  comparisons_converged <- rbind(
    transform(matched_converged$comparisons, precision = "machine"),
    kept, .with_output_columns(converged_checks))
  case_converged <- harness_case(
    comparisons_converged,
    output_ledger(label_converged, fit_outputs(fit), comparisons_converged))

  case <- stack_cases(list(case_published, case_converged))
  case$comparisons <- transform(case$comparisons, source = "Mplus 8.8 (UG example 10.7)",
                                dataset = "ex10.7")
  case
}
