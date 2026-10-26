# The two-level latent class model, written from its definition in base R.
#
# Shared by the glca and multilevLCA suites (each suite file sources it).
# Neither program reports every quantity the package does: glca gives no
# per-group log-likelihood and drops rows with no observed item, and neither
# reports observed-information standard errors or a usable per-group
# outer-product (BHHH) estimate (see the suites' headers for why theirs are
# not used). What they do not report is computed here, at their estimates,
# from the model's definition (Vermunt 2003):
#
#   f(y_g) = sum_w delta_w prod_{i in g} sum_c gamma_wc prod_m rho_mc(y_im),
#
# with a missing item contributing a factor of one. Every program-reported
# quantity this code can also produce is checked against the program's own
# value before it is trusted (see `.two_level_check()`), so a disagreement in
# the harness cannot come from a slip here passing silently.
#
# Vermunt, J. K. (2003). Multilevel latent class models. Sociological
# Methodology, 33(1), 213-239.

#' Category codes of the indicators, in the fit's category order
#' @param data A data frame holding the indicators.
#' @param levels Named list of category labels, as `fit$categorical_levels`.
#' @return An integer matrix, rows by indicators, `NA` where missing.
.two_level_codes <- function(data, levels) {
  codes <- vapply(names(levels), function(item) {
    match(as.character(data[[item]]), levels[[item]])
  }, integer(nrow(data)))
  codes <- matrix(codes, nrow(data), dimnames = list(NULL, names(levels)))
  stopifnot("every observed response must be one of the fit's categories" =
              identical(unname(is.na(codes)), unname(is.na(as.matrix(data[names(levels)])))))
  codes
}

#' Evaluate the two-level model at given parameters
#'
#' @param delta Group-class probabilities.
#' @param gamma Group-classes-by-profiles matrix of P(profile | group class).
#' @param rho Named list, one profiles-by-categories matrix per indicator.
#' @param codes Category codes from [.two_level_codes()].
#' @param groups Group identifier of each row.
#' @return A list: `log_likelihood`; `group_log_likelihood`, named by group in
#'   first-occurrence order; `group_posteriors`, groups by group classes with
#'   those row names; `posteriors`, the rows' marginal profile posteriors
#'   `P(c | y_g) = sum_w P(w | y_g) P(c | y_i, w)`; and
#'   `row_conditional_entropy`, each row's profile entropy within its group
#'   class, `sum_w P(w | y_g) H(c | y_i, w)` (see
#'   [two_level_row_conditional_entropy()]).
two_level_evaluate <- function(delta, gamma, rho, codes, groups) {
  group_ids <- unique(as.character(groups))
  group_index <- match(as.character(groups), group_ids)
  n_profiles <- ncol(gamma)
  # log f(y_i | c), the observed items only: a missing item integrates to one.
  log_conditional <- Reduce(`+`, lapply(colnames(codes), function(item) {
    observed <- !is.na(codes[, item])
    out <- matrix(0, nrow(codes), n_profiles)
    out[observed, ] <- t(log(rho[[item]]))[codes[observed, item], , drop = FALSE]
    out
  }))
  peak <- apply(log_conditional, 1L, max)
  scaled <- exp(log_conditional - peak)
  mixed <- scaled %*% t(gamma)
  per_group <- rowsum(log(mixed) + peak, group_index, reorder = TRUE)
  joint <- sweep(per_group, 2L, log(delta), "+")
  top <- apply(joint, 1L, max)
  group_log_likelihood <- top + log(rowSums(exp(joint - top)))
  group_posteriors <- exp(joint - group_log_likelihood)
  posteriors <- Reduce(`+`, lapply(seq_along(delta), function(w) {
    sweep(scaled, 2L, gamma[w, ], "*") * (group_posteriors[group_index, w] / mixed[, w])
  }))
  row_conditional_entropy <- two_level_row_conditional_entropy(
    scaled, gamma, group_posteriors[group_index, , drop = FALSE])
  dimnames(group_posteriors) <- list(group_ids, NULL)
  list(log_likelihood = sum(group_log_likelihood),
       group_log_likelihood = stats::setNames(as.vector(group_log_likelihood), group_ids),
       group_posteriors = group_posteriors, posteriors = posteriors,
       row_conditional_entropy = row_conditional_entropy)
}

#' Each row's profile entropy within its group class
#'
#' `E_C[H(c | y_i, C)] = -sum_w P(w | y_g) sum_c q_icw log q_icw`, with
#' `q_icw = gamma_wc f_c(y_i) / sum_c' gamma_wc' f_c'(y_i)` the profile
#' posterior of row `i` given group class `w`. Added over the rows and to the
#' group-class entropy it is the entropy of the complete-data classification
#' (group classes and profiles) that ICL, AWE and CLC charge. Written from the
#' definition, for `harness_reference(row_conditional_entropy = )`.
#'
#' @param densities Rows-by-profiles `f_c(y_i)`, each row scaled by any
#'   positive constant (the posteriors given a class do not depend on it).
#' @param gamma Group-classes-by-profiles matrix of P(profile | group class).
#' @param row_group_posteriors Rows-by-group-classes matrix: each row's group's
#'   posterior `P(w | y_g)`.
#' @return One non-negative number per row.
two_level_row_conditional_entropy <- function(densities, gamma, row_group_posteriors) {
  stopifnot("one group posterior row per data row" =
              nrow(row_group_posteriors) == nrow(densities),
            "gamma must have one row per group class" =
              nrow(gamma) == ncol(row_group_posteriors))
  Reduce(`+`, lapply(seq_len(nrow(gamma)), function(w) {
    joint <- sweep(densities, 2L, gamma[w, ], "*")
    within <- joint / rowSums(joint)
    row_group_posteriors[, w] * -rowSums(ifelse(within > 0, within * log(within), 0))
  }))
}

#' Standard errors of every probability, from the observed information and
#' from the outer product of the group scores
#'
#' The log-likelihood of [two_level_evaluate()] is parameterised by baseline
#' logits (last category the reference) so that numerical steps never leave the
#' simplex. `observed` inverts minus its Hessian, taken by Richardson
#' extrapolation (numDeriv); `opg` inverts the cross-product of the per-group
#' scores, the numDeriv Jacobian of the groups' log-likelihood contributions
#' (the BHHH estimator). At a maximum either inverse transforms with the
#' Jacobian of the reparameterisation, so `J V J'`, with `J` the softmax
#' Jacobian, is the covariance of the probabilities themselves, including each
#' row's last category.
#'
#' A response row with a probability at or below `bound` is held at its
#' estimate and the others are reported conditionally on it, as
#' `parameter_inference(boundary = "fix")` documents; every probability of a
#' held row is on a bound (one at the lower, its complement at the upper) and
#' has an `NA` standard error, as `?get_results` documents for `"responses"`.
#'
#' @inheritParams two_level_evaluate
#' @param bound The probability at or below which a row is held.
#' @return A list with `observed` and, when there are more groups than free
#'   parameters, `opg`; each a list: `delta`, like `delta`; `gamma`, like
#'   `gamma`; `probability`, like `rho`.
two_level_standard_errors <- function(delta, gamma, rho, codes, groups, bound = 1e-6) {
  require_suite_packages("numDeriv", reason = "the reference observed information")
  blocks <- c(list(matrix(delta, 1L)),
              lapply(seq_len(nrow(gamma)), function(w) gamma[w, , drop = FALSE]), rho)
  n_gamma <- nrow(gamma)
  # Every block is a stack of probability rows; each row has its own softmax.
  rows <- do.call(rbind, lapply(seq_along(blocks), function(b) {
    data.frame(block = b, row = seq_len(nrow(blocks[[b]])), width = ncol(blocks[[b]]))
  }))
  values <- do.call(c, lapply(blocks, function(block) {
    lapply(seq_len(nrow(block)), function(r) unname(block[r, ]))
  }))
  rows$held <- vapply(values, function(p) any(p <= bound), logical(1))
  stopifnot(
    "a class probability on its bound is not handled here" =
      !any(rows$held & rows$block <= 1L + n_gamma),
    "only two-category response rows are held at a bound here" =
      all(rows$width[rows$held] == 2L))
  free_rows <- which(!rows$held)
  free_width <- rows$width[free_rows] - 1L
  free_end <- cumsum(free_width)
  theta <- unlist(lapply(values[free_rows], function(p) log(p[-length(p)] / p[length(p)])))
  to_values <- function(theta) {
    current <- values
    current[free_rows] <- lapply(seq_along(free_rows), function(i) {
      logits <- c(theta[free_end[[i]] - free_width[[i]] + seq_len(free_width[[i]])], 0)
      exp(logits - max(logits)) / sum(exp(logits - max(logits)))
    })
    current
  }
  evaluate <- function(theta) {
    current <- to_values(theta)
    by_block <- lapply(seq_along(blocks), function(b) do.call(rbind, current[rows$block == b]))
    two_level_evaluate(as.vector(by_block[[1L]]), do.call(rbind, by_block[1L + seq_len(n_gamma)]),
                       stats::setNames(by_block[-seq_len(1L + n_gamma)], names(rho)),
                       codes, groups)
  }
  # The cross-product of G group scores has rank at most G, so the OPG exists
  # only with more groups than free parameters (the package's
  # `latents_too_few_groups` rule).
  opg_possible <- length(unique(groups)) > length(theta)
  covariances <- Filter(Negate(is.null), list(
    observed = solve(-numDeriv::hessian(function(theta) evaluate(theta)$log_likelihood, theta)),
    opg = if (opg_possible) solve(crossprod(numDeriv::jacobian(function(theta) {
      evaluate(theta)$group_log_likelihood
    }, theta)))))
  # Softmax Jacobian, block diagonal over the free rows: d p_k / d a_j =
  # p_k (1[k = j] - p_j) for the free logits j of the same row.
  natural_end <- cumsum(rows$width)
  cells <- do.call(rbind, lapply(seq_along(free_rows), function(i) {
    width <- rows$width[[free_rows[[i]]]]
    p <- values[[free_rows[[i]]]]
    grid <- expand.grid(k = seq_len(width), j = seq_len(width - 1L))
    data.frame(natural = natural_end[[free_rows[[i]]]] - width + grid$k,
               free = free_end[[i]] - free_width[[i]] + grid$j,
               value = p[grid$k] * ((grid$k == grid$j) - p[grid$j]))
  }))
  jacobian <- matrix(0, sum(rows$width), length(theta))
  jacobian[cbind(cells$natural, cells$free)] <- cells$value
  held_cells <- rep(rows$held, rows$width)
  cell_block <- rep(rows$block, rows$width)
  lapply(covariances, function(covariance) {
    errors <- sqrt(diag(jacobian %*% covariance %*% t(jacobian)))
    errors[held_cells] <- NA_real_
    shaped <- lapply(seq_along(blocks), function(b) {
      out <- matrix(errors[cell_block == b], nrow(blocks[[b]]), ncol(blocks[[b]]),
                    byrow = TRUE)
      dimnames(out) <- dimnames(blocks[[b]])
      out
    })
    list(delta = as.vector(shaped[[1L]]),
         gamma = do.call(rbind, shaped[1L + seq_len(n_gamma)]),
         probability = stats::setNames(shaped[-seq_len(1L + n_gamma)], names(rho)))
  })
}

#' Check this file's evaluation against a program's own reported value
#'
#' Stops, rather than reporting a comparison, because a failure here is a
#' defect in the reference code or a misreading of the program's output, not
#' a finding about the package.
#'
#' @param what What is being checked, for the message.
#' @param ours,theirs Numeric vectors or matrices of the same size.
#' @param tolerance Largest allowed absolute difference.
#' @return `TRUE`, invisibly.
.two_level_check <- function(what, ours, theirs, tolerance) {
  gap <- max(abs(as.numeric(ours) - as.numeric(theirs)))
  if (!(length(ours) == length(theirs) && is.finite(gap) && gap <= tolerance)) {
    stop(sprintf("the two-level reference does not reproduce %s (largest gap %.3g, allowed %.3g)",
                 what, gap, tolerance), call. = FALSE)
  }
  invisible(TRUE)
}

#' Project response probabilities onto the package's parameter space
#'
#' `multilpa()` bounds every response probability below by `min_probability`
#' (1e-10 by default), so where another program's maximum puts a probability
#' at zero, the package's maximum puts it at the bound. The projection is the
#' comparable reference: it changes the log-likelihood by about the bound
#' times the number of responses, and it makes the logit `threshold` of such a
#' row a finite number (-23.03 at 1e-10) rather than `-Inf`.
#'
#' @param rho Named list of profiles-by-categories matrices.
#' @param floor The package's `min_probability`.
#' @return `rho`, each row raised to `floor` and renormalised.
two_level_floor <- function(rho, floor) {
  lapply(rho, function(block) {
    raised <- pmax(block, floor)
    raised / rowSums(raised)
  })
}

#' The class permutations `align_reference()` applies
#'
#' The same exhaustive matches, profiles on the response probabilities and
#' group classes on the profile mix, for the tables `compare_fit()` does not
#' carry (the class-probability standard errors).
#'
#' @param rho Named list of the reference's profiles-by-categories matrices.
#' @param gamma The reference's group-classes-by-profiles matrix.
#' @param fit The fitted model.
#' @return A list with `profiles` and `group_classes`: `rho[[m]][profiles, ]`
#'   and `gamma[group_classes, profiles]` are in the fit's order.
two_level_alignment <- function(rho, gamma, fit) {
  ours <- flatten_blocks(lapply(fit$response_probabilities[names(rho)],
                                function(block) unname(as.matrix(block))))
  theirs <- flatten_blocks(lapply(rho, function(block) unname(as.matrix(block))))
  profiles <- align_classes(ours, theirs)
  group_classes <- if (nrow(gamma) > 1L) {
    align_classes(unname(fit$profile_probabilities), unname(gamma[, profiles, drop = FALSE]))
  } else 1L
  list(profiles = profiles, group_classes = group_classes)
}

#' Compare `parameter_inference()` standard errors with a reference's
#'
#' One comparison row per parameter level, over every parameter of it. The
#' rows carry `output = NA`: the class-probability standard errors are not a
#' column of any `get_results()` table, and the OPG ones are not what the
#' `responses` table reports.
#'
#' @param label Label of the fit.
#' @param fit The fitted model.
#' @param reference_errors One estimator's list from
#'   [two_level_standard_errors()], in the reference's class order, with the
#'   fit's category labels.
#' @param alignment From [two_level_alignment()].
#' @param vcov_type The package's matching `vcov_type`.
#' @param levels The parameter levels to compare.
#' @return A comparison `data.frame`.
two_level_error_rows <- function(label, fit, reference_errors, alignment, vcov_type,
                                 levels = c("measurement", "profile", "group")) {
  profiles <- alignment$profiles
  group_classes <- alignment$group_classes
  measurement <- do.call(rbind, lapply(names(reference_errors$probability), function(item) {
    block <- reference_errors$probability[[item]][profiles, , drop = FALSE]
    grid <- expand.grid(profile = seq_len(nrow(block)), category = seq_len(ncol(block)))
    data.frame(level = "measurement", outcome = sprintf("profile_%d", grid$profile),
               term = paste0(item, ":", colnames(block)[grid$category]),
               reference = block[cbind(grid$profile, grid$category)])
  }))
  gamma <- reference_errors$gamma[group_classes, profiles, drop = FALSE]
  grid <- expand.grid(group_class = seq_len(nrow(gamma)), profile = seq_len(ncol(gamma)))
  expected <- rbind(
    measurement,
    data.frame(level = "profile", outcome = sprintf("profile_%d", grid$profile),
               term = sprintf("group_class_%d", grid$group_class),
               reference = gamma[cbind(grid$group_class, grid$profile)]),
    data.frame(level = "group", outcome = sprintf("group_class_%d", seq_along(group_classes)),
               term = "", reference = reference_errors$delta[group_classes]))
  # `boundary = "fix"` holds response probabilities on their bound, as the
  # reference does; on a fit with none it changes nothing.
  table <- parameter_inference(fit, vcov_type = vcov_type, boundary = "fix")
  obtained <- data.frame(level = table$level, outcome = table$outcome,
                         term = ifelse(is.na(table$term), "", table$term),
                         obtained = table$standard_error)
  joined <- merge(expected, obtained, by = c("level", "outcome", "term"), all.x = TRUE,
                  sort = TRUE)
  do.call(rbind, lapply(levels, function(level) {
    at <- joined$level == level
    compare_block(sprintf("%s: %s-level standard errors, %s, against the reference",
                          label, level, vcov_type),
                  joined$reference[at], joined$obtained[at],
                  tolerance = .output_default_tolerances[["responses$probability_standard_error"]])
  }))
}
