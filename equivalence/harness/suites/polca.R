# Single-level limit against poLCA (Linzer & Lewis 2011, J Stat Soft 42(10)).
#
# With one group class the two-level likelihood factorizes over individuals, so
# multilpa must reduce exactly to ordinary latent class analysis. Every output a
# fit reports is compared through `compare_fit()` (see outputs.R), and four
# claims that are not outputs of the fit are kept beside them:
#   1. optimum      - both programs find the same maximised log-likelihood;
#   2. likelihood   - multilpa scores poLCA's own estimates identically, which
#                     tests the likelihood function rather than the optimizer;
#   3. formula      - a closed-form latent class log-likelihood written here,
#                     independent of both programs, agrees at those estimates;
#   4. errors       - poLCA's own standard errors agree with the package's
#                     estimator of the same kind, the outer product of the
#                     scores (see `.polca_inference_rows()`); the
#                     observed-information column the fit reports is compared
#                     against errors computed here from the textbook
#                     likelihood (see `.polca_standard_errors()`).

#' Closed-form latent class log-likelihood, written from the definition
#'
#' @param codes Integer matrix of category codes, observations by items.
#' @param shares Numeric vector of class shares.
#' @param probabilities List of class-by-category matrices, one per item.
#' @return The observed-data log-likelihood, a single numeric value.
.textbook_lca_log_likelihood <- function(codes, shares, probabilities) {
  per_class <- vapply(seq_along(shares), function(class) {
    contributions <- vapply(seq_len(ncol(codes)), function(item) {
      code <- codes[, item]
      value <- log(probabilities[[item]][class, ][code])
      value[is.na(code)] <- 0          # observed-data likelihood: skip missing
      value
    }, numeric(nrow(codes)))
    log(shares[[class]]) + rowSums(contributions)
  }, numeric(nrow(codes)))
  peak <- apply(per_class, 1L, max)
  sum(peak + log(rowSums(exp(per_class - peak))))
}

suite_polca <- function() {
  require_suite_packages(
    "poLCA", "numDeriv",
    reason = paste("the single-level categorical reference, its five bundled datasets",
                   "and the reference observed information"))
  cases <- list(
    list(dataset = "carcinoma", items = c("A", "B", "C", "D", "E", "F", "G"),
         classes = 2:4),
    list(dataset = "cheating",
         items = c("LIEEXAM", "LIEPAPER", "FRAUD", "COPYEXAM"), classes = 2:3),
    list(dataset = "gss82",
         items = c("PURPOSE", "ACCURACY", "UNDERSTA", "COOPERAT"), classes = 2:3),
    list(dataset = "values", items = c("A", "B", "C", "D"), classes = 2L),
    list(dataset = "election",
         items = c("MORALG", "CARESG", "KNOWG", "LEADG", "DISHONG", "INTELG"),
         classes = 2:3)
  )
  stack_cases(unlist(lapply(cases, function(case) {
    lapply(case$classes, function(k) .polca_case(case$dataset, case$items, k))
  }), recursive = FALSE))
}

#' Run one dataset-by-class-count poLCA comparison
#' @param dataset Name of a poLCA dataset.
#' @param items Character vector of manifest item names.
#' @param k Number of latent classes.
#' @return A `harness_case`.
.polca_case <- function(dataset, items, k) {
  raw <- get(utils::data(list = dataset, package = "poLCA",
                         envir = environment()), envir = environment())
  # Compared on complete cases only. poLCA's `na.rm = FALSE` path does not
  # report the observed-data log-likelihood: on `election` it returns
  # -10848.839695 where both multilpa and a textbook likelihood written from
  # the definition give -10834.784505 at poLCA's own estimates, and its class
  # shares sum to 0.99216 rather than one. Equivalence under missing data is
  # established against glca instead, which does agree to 6e-08.
  observed <- raw[stats::complete.cases(raw[, items, drop = FALSE]), items, drop = FALSE]
  row.names(observed) <- NULL
  frame <- observed
  frame$unit <- factor(seq_len(nrow(frame)))

  formula <- stats::as.formula(paste("cbind(", paste(items, collapse = ","), ") ~ 1"))
  set.seed(20260919L)
  searched <- poLCA::poLCA(formula, data = observed, nclass = k, nrep = 30L,
                           maxiter = 10000L, tol = 1e-10, verbose = FALSE,
                           na.rm = TRUE, calc.se = FALSE)
  # poLCA stops on an absolute change in the log-likelihood. Where the surface
  # is flat that leaves its estimates short of the fixed point: on gss82 with
  # two classes a run stopped at tol = 1e-12 still has a response probability
  # 2e-6 from it. Restarting poLCA() from its own probabilities does not help,
  # because a restart resets the class shares to uniform. The winning solution
  # is therefore continued with poLCA's own E and M steps, without a stopping
  # tolerance, as the mclust suite continues mclust's, so the reference is
  # poLCA's estimate and not where its stopping rule happened to fire.
  reference <- .polca_refine(searched)

  fit <- multilpa(frame, vars = items, id = "unit", n_profiles = k,
                  n_group_classes = 1L, categorical = items, n_starts = 30L,
                  seed = 20260919L, tol = 1e-16, max_iter = 20000L)
  # multilpa stops when a step gains less than `tol * (1 + |log L|)`. At the
  # 1e-14 this suite used before, that is a gain of 2.8e-11 on gss82 with two
  # classes, and the fit stopped 2.0e-11 below poLCA's fixed point with a
  # response probability 3.2e-7 from it (EM continued from multilpa's
  # estimates reaches poLCA's to 3.7e-14) and posteriors 2.9e-6 off: the
  # surface is flat there (that probability's standard error is 0.062), so a
  # likelihood gap at the rounding level of |log L| leaves parameters well
  # short. At 1e-16 the bound is below that rounding, so the rule becomes
  # "stop when the likelihood stops rising", the rule the reference is refined
  # by, and the fit must still report convergence. Even then a likelihood rule
  # cannot resolve the last digits of a flat surface: on carcinoma with four
  # classes multilpa stops 1.4e-12 below the bounded fixed point with a
  # probability 3.6e-7 from it (EM on the bounded space continued from its
  # estimates reaches the reference's to 7e-13), and its posteriors agree to
  # 5.8e-7, inside the default 1e-6.
  stopifnot("multilpa must converge at the tolerance asked of it" = isTRUE(fit$converged))
  floor <- fit$min_probability %||% 1e-10

  # poLCA's estimates on the package's parameter space. multilpa bounds every
  # response probability below by `min_probability` (1e-10) and holds a
  # category whose unconstrained share falls below it exactly there (the KKT
  # solution in `.multilpa_bound_probabilities()`); poLCA lets the same
  # probability fall to exactly zero (carcinoma with three and four classes).
  # Both are the maximum of their own space, and the two differ by no more
  # than 1e-10 in any probability, but a
  # threshold (a logit of a cumulative probability) of a zero is infinite where
  # the bounded one is -23.03. The reference is moved onto the bounded space
  # exactly as the KKT solution defines it, rather than any tolerance loosened.
  # The projection is the bounded maximum to within rounding of the floor: EM
  # on the bounded space (poLCA's E and M steps, then the bound) started from
  # it moves no probability by more than 1.1e-9 on carcinoma with four
  # classes, which has 32 cells on the bound.
  # The posteriors, the likelihood and each row's contribution are poLCA's own
  # (its E step and its likelihood expression) at the bounded estimates, so
  # that every reference table describes one point. The posteriors differ from
  # those at the unbounded estimates only where a floored probability enters,
  # by 2.6e-8 at most (election with three classes), but a residual table whose
  # chi-square is of order 1e-10 has a Cramer's V of order 1e-5 (V is a square
  # root), and taking the posteriors at the unbounded point moved V on
  # carcinoma with four classes from multilpa's 1.5e-5 to exactly 0.
  bounded <- lapply(reference$probs, .polca_bound, floor = floor)
  point <- .polca_evaluate(bounded, reference$P, reference$codes)
  # Raising a zero to the floor costs about (class size) x floor of likelihood
  # per held cell, 2.4e-8 on carcinoma with two classes; it must stay well
  # inside the 1e-6 at which the likelihoods are compared below.
  stopifnot("bounding must not move the maximum by a comparable amount" =
              abs(point$llik - reference$llik) <= 1e-7)
  labelled <- .polca_label(bounded, reference$codes, observed, fit$categorical_levels)

  # Score poLCA's own estimates under the multilpa likelihood.
  at_reference <- multilpa(frame, vars = items, id = "unit", n_profiles = k,
                           n_group_classes = 1L, categorical = items, n_starts = 1L,
                           max_iter = 0L,
                           start = starting_values(list(
                             profile_probabilities = matrix(reference$P, 1L, k),
                             group_probabilities = 1,
                             response_probabilities = labelled)))
  textbook <- .textbook_lca_log_likelihood(reference$codes, reference$P, bounded)

  label <- sprintf("%s, %d classes", dataset, k)
  n_categories <- vapply(fit$response_probabilities, ncol, integer(1))

  # Parameters are only worth comparing where they are identified. The
  # contingency table has prod(categories) - 1 degrees of freedom; when that
  # barely exceeds the number of free parameters the likelihood has a flat
  # ridge, and two programs reaching the same maximum can still report quite
  # different estimates. `cheating` with three classes is such a case: four
  # binary items give 15 degrees of freedom against 14 free parameters, and
  # poLCA and multilpa agree on the log-likelihood to 2e-12, and on every
  # unit's contribution to it to 9e-13, while reporting class shares of
  # 0.035/0.111/0.854 and 0.045/0.078/0.877.
  residual_df <- prod(n_categories) - 1L - fit$n_parameters
  identified <- residual_df >= 3L

  extra <- rbind(
    # The absolute statement the relative default (1e-6 of |log L|) is looser
    # than: on election that default would admit a gap of 0.01.
    compare_values(sprintf("%s: maximised log-likelihood", label),
                   point$llik, as.numeric(logLik(fit)), tolerance = 1e-6),
    compare_values(sprintf("%s: log-likelihood at poLCA's estimates", label),
                   point$llik, as.numeric(logLik(at_reference)), tolerance = 1e-6),
    compare_values(sprintf("%s: textbook log-likelihood at poLCA's estimates", label),
                   point$llik, textbook, tolerance = 1e-6),
    # poLCA's own criteria, as it prints them (the outputs are compared
    # against criteria written from their definitions by `compare_fit()`).
    compare_values(sprintf("%s: poLCA's BIC", label),
                   searched$bic, fit$bic_individual, tolerance = 1e-5),
    compare_values(sprintf("%s: poLCA's AIC", label),
                   searched$aic, fit$aic, tolerance = 1e-5),
    # poLCA defines this as min(N, prod(categories) - 1) - npar: the response
    # pattern table bounds the degrees of freedom once it is smaller than the
    # sample.
    compare_values(sprintf("%s: residual degrees of freedom", label),
                   searched$resid.df,
                   min(fit$n_informative, prod(n_categories) - 1L) - fit$n_parameters,
                   tolerance = 0.5))

  # The reference observed information needs an identified maximum; on a flat
  # ridge it is singular, and nothing parameter-level is compared there.
  errors <- if (identified) .polca_standard_errors(reference$codes, reference$P, bounded)
  harness <- harness_reference(
    log_likelihood = point$llik, n_parameters = searched$npar,
    posteriors = point$posterior,
    responses = labelled,
    standard_errors = if (identified) list(probability = .polca_label(
      errors$probability, reference$codes, observed, fit$categorical_levels)),
    profile_probabilities = matrix(reference$P, 1L),
    group_probabilities = 1,
    group_log_likelihood = stats::setNames(point$by_row, seq_len(nrow(observed))),
    data = observed)

  case <- if (identified) {
    compare_fit(label, fit, harness,
                # A p-value is the chi-square tail of its statistic, so it moves
                # by the chi-square density times the statistic's change, and
                # for one degree of freedom that density is unbounded as the
                # statistic goes to zero. On carcinoma with four classes the
                # statistics agree to 5.6e-7 (multilpa is 3.6e-7 short of the
                # fixed point, see above) and the worst p-value cell, at a
                # statistic of 0.021 where the density is 2.7, moves by
                # 1.1e-6. 1e-5 admits a statistic gap of 1e-6 where the
                # density is up to 10 (a statistic above 0.0016).
                tolerance = c(`residuals$p_value` = 1e-5, `residuals$p_adjusted` = 1e-5),
                extra = rbind(extra, .polca_inference_rows(label, fit, reference$codes,
                                                           reference$P, bounded, point,
                                                           errors, labelled, observed)))
  } else {
    .polca_unidentified_case(label, fit, harness, extra)
  }
  case$comparisons <- transform(case$comparisons,
                                source = sprintf("poLCA %s", utils::packageVersion("poLCA")),
                                dataset = dataset, precision = "machine")
  case
}

#' Continue poLCA's EM from its own solution to the fixed point
#'
#' Uses poLCA's own E step (`poLCA.postClass.C`) and M step (`poLCA.probHat.C`,
#' shares as mean posteriors), the two calls its `poLCA()` loop makes, until no
#' parameter moves by more than 1e-15.
#'
#' @param searched A `poLCA` fit.
#' @param max_steps Bound on the number of EM steps; reaching it is an error.
#' @param tolerance Largest parameter change at the fixed point.
#' @return A list with `probs`, `P`, `llik` (by [.polca_evaluate()]), `codes`
#'   and `steps`.
.polca_refine <- function(searched, max_steps = 500000L, tolerance = 1e-15) {
  posterior_of <- utils::getFromNamespace("poLCA.postClass.C", "poLCA")
  m_step <- utils::getFromNamespace("poLCA.probHat.C", "poLCA")
  vectorize <- utils::getFromNamespace("poLCA.vectorize", "poLCA")
  unvectorize <- utils::getFromNamespace("poLCA.unvectorize", "poLCA")
  # poLCA returns the data it fitted with factors restored; internally it works
  # on `cbind()` of the items, which is each factor's level code and each
  # numeric item's value.
  n <- nrow(searched$y)
  codes <- vapply(searched$y, as.integer, integer(n))
  codes <- matrix(codes, n, dimnames = list(NULL, names(searched$y)))
  k <- length(searched$P)
  vp <- vectorize(searched$probs)
  shares <- searched$P
  steps <- 0L
  moved <- Inf
  # EM is a fixed-point iteration: each step starts from the last one, so the
  # sequence has no vectorised form. It is bounded and the bound is checked.
  while (moved > tolerance && steps < max_steps) {
    posterior <- posterior_of(matrix(shares, n, k, byrow = TRUE), vp, codes)
    probabilities <- m_step(posterior, codes, vp)
    updated <- colMeans(posterior)
    moved <- max(abs(c(probabilities - vp$vecprobs, updated - shares)))
    vp$vecprobs <- probabilities
    shares <- updated
    steps <- steps + 1L
  }
  stopifnot("poLCA's EM must reach its fixed point within the step bound" =
              moved <= tolerance)
  probs <- stats::setNames(unvectorize(vp), names(searched$probs))
  llik <- .polca_evaluate(probs, shares, codes)$llik
  stopifnot("refining must not lose likelihood" = llik >= searched$llik - 1e-9)
  list(probs = probs, P = shares, llik = llik, codes = codes, steps = steps)
}

#' poLCA's posteriors and likelihood at given estimates
#'
#' poLCA's own E step, and the log-likelihood by the expression its `poLCA()`
#' loop uses (`poLCA.ylik.C` returns each class's likelihood scaled by
#' `.Machine$double.xmax`, which the expression takes back out).
#'
#' @param probs Classes-by-categories matrices, in poLCA's code order.
#' @param shares Class shares.
#' @param codes poLCA's code matrix.
#' @return A list with `posterior`, `by_row` (each row's log-likelihood
#'   contribution) and `llik`.
.polca_evaluate <- function(probs, shares, codes) {
  posterior_of <- utils::getFromNamespace("poLCA.postClass.C", "poLCA")
  likelihoods <- utils::getFromNamespace("poLCA.ylik.C", "poLCA")
  vectorize <- utils::getFromNamespace("poLCA.vectorize", "poLCA")
  vp <- vectorize(lapply(probs, unname))
  prior <- matrix(shares, nrow(codes), length(shares), byrow = TRUE)
  by_row <- log(rowSums(prior * likelihoods(vp, codes))) - log(.Machine$double.xmax)
  list(posterior = posterior_of(prior, vp, codes), by_row = by_row, llik = sum(by_row))
}

#' Move one item's probabilities onto multilpa's bounded simplex
#'
#' The Karush-Kuhn-Tucker solution of the bounded multinomial maximum: every
#' category below the floor is held at it, and the others keep their
#' proportions, rescaled to the probability left over.
#'
#' @param block Classes-by-categories probability matrix.
#' @param floor The lower bound.
#' @return The bounded matrix.
.polca_bound <- function(block, floor) {
  bounded <- t(apply(block, 1L, function(p) {
    held <- p < floor
    p[!held] <- p[!held] * (1 - sum(held) * floor) / sum(p[!held])
    p[held] <- floor
    p
  }))
  stopifnot("rescaling must not push another category below the floor" =
              all(bounded >= floor * (1 - 1e-12)))
  matrix(bounded, nrow(block), dimnames = dimnames(block))
}

#' Give poLCA's category columns the fit's category labels and order
#'
#' poLCA codes each item 1..C (a factor by its level codes, a number by its
#' value); the label of code `c` is read off the data rather than assumed.
#'
#' @param blocks Classes-by-categories matrices in poLCA's code order.
#' @param codes poLCA's code matrix.
#' @param observed The data poLCA was fitted to.
#' @param levels The fit's `categorical_levels`.
#' @return The same matrices with columns named and ordered as the fit's.
.polca_label <- function(blocks, codes, observed, levels) {
  stats::setNames(lapply(seq_along(blocks), function(j) {
    item <- colnames(codes)[[j]]
    pairs <- unique(data.frame(code = codes[, j], label = as.character(observed[[item]])))
    stopifnot("each poLCA code must stand for exactly one category" =
                !anyDuplicated(pairs$code) && !anyDuplicated(pairs$label),
              "poLCA and the fit must see the same categories" =
                setequal(pairs$label, levels[[item]]))
    block <- blocks[[j]]
    colnames(block) <- pairs$label[match(seq_len(ncol(block)), pairs$code)]
    rownames(block) <- NULL
    block[, levels[[item]], drop = FALSE]
  }), colnames(codes))
}

#' Observed-information standard errors of a latent class model, from scratch
#'
#' The log-likelihood is [.textbook_lca_log_likelihood()]. Its Hessian at
#' poLCA's estimates is taken by Richardson extrapolation (numDeriv) in
#' baseline-category logits: in each class and item, the log-odds of every
#' category off the boundary against the largest one, and the log-odds of
#' every class share against the largest share. The inverse is carried to the
#' natural parameters, every probability and share, by the delta method with
#' the softmax Jacobian written out below. At a maximum the inverse observed
#' information transforms with the Jacobian of any reparameterisation, so these
#' are the standard errors of the natural parameters whatever scale either
#' program works on.
#'
#' The Hessian is not taken in the natural parameters themselves, as the
#' mclust suite takes it for means and variances. A step there has to be a
#' small fraction of the smallest probability to stay inside the simplex, and
#' at such steps the second differences are rounding noise: on gss82 with two
#' classes the error of P(ACCURACY = Mostly true) in the class where it is 0.03
#' came out 0.062247, 0.062236 and 0.062331 for numDeriv steps of 1e-2, 3e-3
#' and 1e-3 of each coordinate, and the information was not positive definite
#' at 1e-4. Logits are unbounded, so a step of 0.1 is always admissible; every
#' error on gss82 then moves by at most 3e-7 over first steps from 0.4 to
#' 0.02, and on carcinoma with four classes by at most 2e-6 down to 0.05.
#'
#' A probability on its boundary (at or below 1e-6, where its expected count
#' is below 0.002 of an observation in the largest dataset here, so the data
#' cannot tell it from zero) is not an interior maximum: the likelihood's
#' curvature there does not describe its sampling distribution (Self & Liang
#' 1987), so it has no Wald standard error and the reference says `NA`. It is
#' held at its value, and the other parameters' errors are conditional on it,
#' which is what the package's `boundary = "fix"` documents. A category whose
#' siblings are all held sits at one minus their sum, on the upper boundary,
#' and is `NA` for the same reason.
#'
#' @param codes Integer code matrix.
#' @param shares Class shares.
#' @param probabilities Classes-by-categories matrices, in poLCA's code order.
#' @param boundary The probability at or below which a cell is on its bound.
#' @param step The first step of the Richardson sequence, on the logit scale.
#' @return A list with `probability` (matrices like `probabilities`) and
#'   `share` (one per class).
.polca_standard_errors <- function(codes, shares, probabilities, boundary = 1e-6,
                                   step = 1e-1) {
  k <- length(shares)
  # One softmax per class and item with at least two categories off the
  # boundary: the held categories keep their values, and the rest share what
  # is left, `mass`, through logits against the largest of them.
  rows <- unlist(lapply(seq_along(probabilities), function(item) {
    block <- probabilities[[item]]
    lapply(seq_len(k), function(class) {
      off_bound <- which(block[class, ] > boundary)
      if (length(off_bound) < 2L) return(NULL)
      baseline <- off_bound[[which.max(block[class, off_bound])]]
      list(item = item, class = class, free = setdiff(off_bound, baseline),
           baseline = baseline, mass = sum(block[class, off_bound]))
    })
  }), recursive = FALSE)
  rows <- Filter(Negate(is.null), rows)
  widths <- vapply(rows, function(row) length(row$free), integer(1))
  ends <- cumsum(widths)
  rows <- Map(function(row, end, width) {
    row$positions <- end - width + seq_len(width)
    row
  }, rows, ends, widths)
  n_free <- sum(widths)
  share_baseline <- which.max(shares)
  free_shares <- setdiff(seq_len(k), share_baseline)
  theta <- c(unlist(lapply(rows, function(row) {
    block <- probabilities[[row$item]]
    log(block[row$class, row$free] / block[row$class, row$baseline])
  })), log(shares[free_shares] / shares[[share_baseline]]))
  softmax <- function(logits) {
    exponentiated <- exp(c(logits, 0) - max(c(logits, 0)))
    exponentiated / sum(exponentiated)
  }
  unpack <- function(theta) {
    blocks <- Reduce(function(blocks, row) {
      blocks[[row$item]][row$class, c(row$free, row$baseline)] <-
        row$mass * softmax(theta[row$positions])
      blocks
    }, rows, probabilities)
    unpacked_shares <- numeric(k)
    unpacked_shares[c(free_shares, share_baseline)] <-
      softmax(theta[n_free + seq_along(free_shares)])
    list(probabilities = blocks, shares = unpacked_shares)
  }
  log_likelihood <- function(theta) {
    parameters <- unpack(theta)
    .textbook_lca_log_likelihood(codes, parameters$shares, parameters$probabilities)
  }
  # One absolute first step for every logit (numDeriv's own steps are relative,
  # which shrinks them to rounding noise for a logit near zero).
  information <- -numDeriv::hessian(log_likelihood, theta,
                                    method.args = list(d = 0, eps = step, zero.tol = Inf))
  covariance <- chol2inv(chol(information))
  # The softmax Jacobian: d p_c / d logit_d = p_c (1[c = d] - p_d / mass), for
  # every category c of the softmax (the baseline included) and every logit d.
  softmax_rows <- function(values, categories, free, positions, mass) {
    t(vapply(categories, function(category) {
      out <- numeric(length(theta))
      out[positions] <- values[[category]] * ((free == category) - values[free] / mass)
      out
    }, numeric(length(theta))))
  }
  cell_errors <- lapply(probabilities, function(block) {
    matrix(NA_real_, nrow(block), ncol(block), dimnames = dimnames(block))
  })
  cell_errors <- Reduce(function(errors, row) {
    values <- probabilities[[row$item]][row$class, ]
    categories <- c(row$free, row$baseline)
    jacobian <- softmax_rows(values, categories, row$free, row$positions, row$mass)
    errors[[row$item]][row$class, categories] <-
      sqrt(rowSums((jacobian %*% covariance) * jacobian))
    errors
  }, rows, cell_errors)
  share_jacobian <- softmax_rows(shares, seq_len(k), free_shares,
                                 n_free + seq_along(free_shares), 1)
  list(probability = cell_errors,
       share = sqrt(rowSums((share_jacobian %*% covariance) * share_jacobian)))
}

#' poLCA's own standard errors, and the share errors, against the package
#'
#' `poLCA.se()` builds the information as the cross-product of each
#' observation's score (`info <- t(s) %*% s`, scores on the baseline-logit
#' scale of the probabilities and the shares), inverts it with `MASS::ginv()`
#' and maps it to the probabilities by the delta method. That is the
#' outer-product-of-gradients (BHHH) estimator, not the observed information,
#' so poLCA's errors are compared with `parameter_inference(vcov_type =
#' "opg")` and not with the observed-information column of the responses
#' table, which is compared against [.polca_standard_errors()] instead. They
#' are evaluated by poLCA's own function at the reference point (the refined
#' estimates on the bounded space, and the posteriors there).
#'
#' Cells on the boundary are left out: their score column is zero to rounding,
#' `ginv()` drops that direction and the delta-method factor `p (1 - p)`
#' vanishes, so poLCA reports zero where the package holds the cell and
#' reports `NA`. Dropping the direction is holding the cell, so the interior
#' errors are conditional on the same thing in both.
#'
#' @param label Label for the rows.
#' @param fit The fitted model.
#' @param codes poLCA's code matrix.
#' @param shares,bounded The reference class shares and response
#'   probabilities, in poLCA's class and code order.
#' @param point The list [.polca_evaluate()] returned at them.
#' @param errors The list [.polca_standard_errors()] returned at them.
#' @param labelled `bounded` with the fit's category labels.
#' @param observed The data poLCA was fitted to.
#' @return Comparison rows, none of them an output of the fit.
.polca_inference_rows <- function(label, fit, codes, shares, bounded, point, errors,
                                  labelled, observed) {
  n <- nrow(codes)
  k <- length(shares)
  polca_errors <- utils::getFromNamespace("poLCA.se", "poLCA")(
    y = codes, x = matrix(1, n, 1L), probs = lapply(bounded, unname),
    prior = matrix(shares, n, k, byrow = TRUE), rgivy = point$posterior)
  as_fit <- function(blocks) .polca_label(blocks, codes, observed, fit$categorical_levels)
  order_classes <- align_classes(flatten_blocks(fit$response_probabilities),
                                 flatten_blocks(labelled))
  in_fit_order <- function(blocks) {
    lapply(blocks, function(block) block[order_classes, , drop = FALSE])
  }
  polca_cells <- in_fit_order(as_fit(polca_errors$probs))
  observed_cells <- in_fit_order(as_fit(errors$probability))
  cells <- do.call(rbind, lapply(names(polca_cells), function(item) {
    block <- polca_cells[[item]]
    data.frame(outcome = sprintf("profile_%d", rep(seq_len(k), ncol(block))),
               term = paste0(item, ":", rep(colnames(block), each = k)),
               polca = as.vector(block),
               interior = !is.na(as.vector(observed_cells[[item]])),
               stringsAsFactors = FALSE)
  }))
  shares <- data.frame(outcome = sprintf("profile_%d", seq_len(k)), term = "group_class_1",
                       polca = polca_errors$P[order_classes],
                       observed = errors$share[order_classes], stringsAsFactors = FALSE)
  # A parameter the package's table does not report gets `NA`, which
  # `compare_block()` counts as a disagreement against a reference value.
  lookup <- function(table, parameter, wanted) {
    rows <- table[table$parameter == parameter, , drop = FALSE]
    rows$standard_error[match(paste(wanted$outcome, wanted$term),
                              paste(rows$outcome, rows$term))]
  }
  opg <- as.data.frame(parameter_inference(fit, vcov_type = "opg", boundary = "fix"))
  information <- as.data.frame(parameter_inference(fit, boundary = "fix"))
  interior <- cells[cells$interior, , drop = FALSE]
  rbind(
    compare_block(sprintf(paste("%s: poLCA's standard errors (OPG) of the response",
                                "probabilities, against the package's OPG"), label),
                  interior$polca, lookup(opg, "response", interior), tolerance = 1e-4),
    compare_block(sprintf(paste("%s: poLCA's standard errors (OPG) of the class shares,",
                                "against the package's OPG"), label),
                  shares$polca, lookup(opg, "probability", shares), tolerance = 1e-4),
    compare_block(sprintf("%s: observed-information standard errors of the class shares", label),
                  shares$observed, lookup(information, "probability", shares),
                  tolerance = 1e-4))
}

#' Compare what a non-identified fit can be compared on
#'
#' On a flat ridge the two programs reach the same maximum at different
#' estimates, so the estimates, the posteriors and everything computed from
#' them describe different points of the ridge. What does not depend on where
#' on the ridge a program stopped is compared: the maximum, the parameter
#' count, every criterion that does not charge the individual classification
#' entropy, each unit's log-likelihood contribution (the fitted probability of
#' its response pattern, the same at every point of the ridge), the
#' group-level tables (one group class), the residual degrees of freedom and
#' the model-implied residual association, which is zero by construction.
#'
#' @param label Label for the rows.
#' @param fit The fitted model.
#' @param harness The `harness_reference`.
#' @param extra Further comparison rows.
#' @return A `harness_case`.
.polca_unidentified_case <- function(label, fit, harness, extra) {
  ridge <- paste("cheating with three classes is not identified: four binary items give",
                 "15 degrees of freedom against 14 free parameters, and poLCA and multilpa",
                 "reach the same maximum, every unit's log-likelihood contribution agreeing",
                 "to 9e-13, at different estimates (class shares 0.035/0.111/0.854 against",
                 "0.045/0.078/0.877)")
  from_posteriors <- sprintf(paste("computed from the posteriors, which on this flat ridge",
                                   "describe a different point for each program (%s)"), ridge)
  entropy_charge <- sprintf(paste("charges the individual classification entropy, computed from",
                                  "posteriors that are not identified on this ridge (%s)"), ridge)
  exempt <- c(
    `responses$probability` = ridge, `responses$threshold` = ridge,
    `responses$probability_standard_error` = sprintf(paste(
      "the observed information is singular along a ridge of maxima, so the",
      "errors either program computes are rounding (multilpa reports 894 for a",
      "probability, the reference gives class-share errors of 71, 49 and 22) at",
      "different points of the ridge (%s)"), ridge),
    `profile_probabilities$probability` = ridge,
    counts = from_posteriors, posteriors = from_posteriors,
    `assignments$profile` = from_posteriors, `assignments$uncertainty` = from_posteriors,
    classification = from_posteriors, average_posteriors = from_posteriors,
    classification_errors = from_posteriors, bch_weights = from_posteriors,
    `entropy$entropy_sum` = from_posteriors, `entropy$relative_entropy` = from_posteriors,
    `residuals$observed` = from_posteriors,
    `residuals$residual` = from_posteriors, `residuals$effective_n` = from_posteriors,
    `residuals$statistic` = from_posteriors, `residuals$p_value` = from_posteriors,
    `residuals$p_adjusted` = from_posteriors,
    `information_criteria$awe_individual` = entropy_charge,
    `information_criteria$icl_individual` = entropy_charge,
    `information_criteria$clc_individual` = entropy_charge,
    # Both conventions charge the same complete-data entropy (round 2), so
    # the groups columns carry the ridge's posteriors too.
    `information_criteria$awe_groups` = entropy_charge,
    `information_criteria$icl_groups` = entropy_charge,
    `information_criteria$clc_groups` = entropy_charge)
  inventory <- fit_outputs(fit)
  # The posterior columns of the assignments are named by profile.
  posterior_columns <- grep("^assignments[$]posterior_profile_", inventory$output, value = TRUE)
  exempt[posterior_columns] <- from_posteriors
  expected <- reference_tables(align_reference(harness, fit), fit)
  obtained <- fit_tables(fit)
  obtained$model$.row <- 1L
  # Every table keeps its keys and the value columns not exempted above.
  kept <- stats::setNames(lapply(names(expected), function(name) {
    if (name %in% names(exempt)) return(NULL)
    table <- expected[[name]]
    dropped <- sub("^[^$]*[$]", "", names(exempt)[startsWith(names(exempt), paste0(name, "$"))])
    table[setdiff(names(table), dropped)]
  }), names(expected))
  kept <- Filter(Negate(is.null), kept)
  comparisons <- rbind(compare_tables(label, obtained, kept, .output_default_tolerances),
                       .with_output_columns(extra))
  harness_case(comparisons, output_ledger(label, inventory, comparisons, exempt))
}
